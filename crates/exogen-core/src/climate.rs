// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! The energy balance (PLAN.md 4.6, Phase 5A, F235): the year's mean surface temperature of
//! every cell of a coarse cubesphere, where the sunlight it absorbs, the heat it radiates
//! to space and the heat its atmosphere carries to and from its neighbours balance (North
//! 1975):
//!
//! ```text
//! Q (1 - a) - OLR(T) + D div grad T = 0
//! ```
//!
//! Q is the year's sunlight on the cell (F234); a the planet's albedo there; OLR the grey
//! atmosphere's outgoing radiation, sigma T^4 / (1 + 3 tau / 4), through the greenhouse
//! optical depth of its gases at the cell's temperature (F21, its water vapour's now
//! saturating, F235); D the transport the atmosphere's mass and the body's rotation allow
//! (Williams and Kasting 1997, F33), zero on an airless body, whose every cell balances
//! alone. Under a scattering atmosphere the albedo rises as the sun's mean height falls (the
//! low sun and the clouds of North and Coakley 1979), and where water is the working
//! liquid, ice and snow brighten the ground the year leaves below freezing (North 1975), less
//! under a redder star, whose near infrared they reflect little of (F250).
//!
//! The balance is solved by Newton's method from the zero-dimensional mean temperature, or
//! where ice and snow brighten the ground from its balance without them (F249), the warm
//! branch, each step's linear system by conjugate gradients; the local slope is floored
//! at a share of the radiation's own so that the ice's feedback, which can make it negative,
//! leaves the system positive definite.
//!
//! Through the days and seasons (F236) the places of a coarser layout run in step under
//! their air (F245): the air over each place holds its column's heat and exchanges heat with
//! its land and its sea by their shares of the place, and the transport carries heat from
//! the air over one place to the air over the next, a share of it through the sea's water
//! (the two layers of Rose and Marshall 2009), so that the transport answers the seasons
//! and reaches the ground only as its air touches it.

use std::f64::consts::{PI, TAU};

use rayon::prelude::*;

use crate::constants::*;
use crate::fields::Field;
use crate::grid::projection::{PaddedField, sample_bilinear};
use crate::grid::{CubeSphere, direction_from_lon_lat, lon_lat_from_direction};
use crate::insolation::{PathStep, Sky, SunTable, Sunlight, diurnal_share, read_grid, read_row};
use crate::params::ParamGraph;
use crate::params::derive::is_icy;
use crate::params::derive::physics::{
    blackbody_share_below, mean_heat_capacity, mean_molecular_mass, regolith_heat,
    saturation_temperature, species_columns, sunlight_pattern, tau_dry, vapour_pressure, volatile,
    water_column, water_optical_depth,
};

/// What the energy balance reads.
#[derive(Clone, Copy, Debug)]
pub struct EnergyBalance {
    pub sky: Sky,
    /// The planet's Bond albedo with the sun at its mean height and no ice.
    pub albedo: f64,
    /// How far the albedo rises per unit fall of the cosine of the sun's mean height.
    pub low_sun: f64,
    /// Whether water, the working liquid, freezes where the year leaves it cold.
    pub ice: bool,
    /// How strongly ice and snow brighten the ground under the star's light, against under
    /// the Sun's, under which their albedos are given (F250).
    pub brightening: f64,
    /// The atmosphere's heat transport, W/m²/K.
    pub transport: f64,
    /// The greenhouse optical depth of the dry gases.
    pub tau_dry: f64,
    /// The water a warm atmosphere can hold vapour from, metres of global layer.
    pub water: f64,
    /// The zero-dimensional mean temperature the solution starts from, K.
    pub start: f64,
    /// The body's own heat flowing out through its surface, W/m², which every cell radiates
    /// with its sunlight: a cell the star never warms stands where it radiates this alone.
    pub heat_flux: f64,
    /// The year, the orbit about the star, s.
    pub year: f64,
    /// The frost of the atmosphere's bulk gas (F239).
    pub frost: Option<Frost>,
    /// The mean anomaly of the snapshot's season, rad from perihelion.
    pub snapshot: f64,
    /// The surface pressure, Pa: its mean over the year.
    pub pressure: f64,
    /// The heat the air exchanges with the ground per kelvin by its turbulence, rho c_p C_H U
    /// (F237), W/m²/K, none without an atmosphere; and, where water is the working liquid,
    /// the vapour its sea gives the air per pascal of saturation pressure the sea stands
    /// above the air's, rho C_H U (M_w / M) / p, kg/m²/s/Pa (F245).
    pub exchange: f64,
    pub evaporation: f64,
    /// The working liquid's vapour, whose latent heat the sea and the wet land give the air
    /// (F269): water's, methane's, any liquid's under an atmosphere.
    pub liquid: Option<&'static crate::constants::Volatile>,
    /// Whether a fast spinner's places stand on its latitudes' rows alone (F259): the years
    /// of the orbit's cycles, run cheaply.
    pub rows: bool,
    /// What the wind reads of the body and its air (F267).
    pub circulation: crate::wind::Air,
}

/// The year's mean temperature of a surface at its own albedo, facing where its slope faces
/// (F256): the climate's year-mean balance there, the year's sunlight on the surface at its
/// albedo and what else the climate brings the ground, its body's heat and its air's
/// transport, radiated through the greenhouse over it. The stages that need the ground's
/// mean temperature read it here: on an airless body the sunlight and the heat alone, its
/// year's sunlight in the direction it faces exact for any obliquity, orbit and resonance.
#[derive(Clone, Debug)]
pub struct YearMean {
    sun: SunTable,
    heat: f64,
    /// Under an atmosphere, on the climate's grid: the heat the air brings each cell, W/m²,
    /// and the share of a black body's radiation the greenhouse over it lets out.
    air: Option<(CubeSphere, PaddedField, PaddedField)>,
}

impl YearMean {
    pub fn new(balance: &EnergyBalance, climate: &Climate) -> Self {
        let grid = climate.temperature.grid();
        let airless = balance.transport <= 0.0 && balance.tau_dry <= 0.0 && balance.water <= 0.0;
        let air = (!airless).then(|| {
            let emission: Vec<f32> = climate
                .temperature
                .as_slice()
                .iter()
                .map(|&t| {
                    let t = t as f64;
                    (balance.olr(t) / (SIGMA_SB * t.powi(4)).max(1e-30)) as f32
                })
                .collect();
            (
                grid,
                PaddedField::new(&grid, climate.transport.as_slice()),
                PaddedField::new(&grid, &emission),
            )
        });
        Self {
            sun: balance.sky.annual_table(),
            heat: balance.heat_flux,
            air,
        }
    }

    /// An airless body's, under its sky, its own heat flowing out at `heat` W/m².
    pub fn airless(sky: &Sky, heat: f64) -> Self {
        Self {
            sun: sky.annual_table(),
            heat: heat.max(0.0),
            air: None,
        }
    }

    /// The year's mean temperature, K, of a surface of `albedo` facing `d`.
    pub fn temperature(&self, d: [f64; 3], albedo: f64) -> f64 {
        let absorbed = self.sun.flux_at(d) * (1.0 - albedo);
        let (brought, emission) = match &self.air {
            Some((grid, brought, emission)) => (
                sample_bilinear(grid, brought, d) as f64,
                (sample_bilinear(grid, emission, d) as f64).max(1e-6),
            ),
            None => (0.0, 1.0),
        };
        radiative_temperature(absorbed + brought + self.heat, emission)
    }

    /// The warmest year-mean temperature anywhere on the body of a surface of `albedo`, K.
    pub fn warmest(&self, albedo: f64) -> f64 {
        self.sun
            .directions()
            .map(|d| self.temperature(d, albedo))
            .fold(0.0, f64::max)
    }
}

/// The temperature, K, at which a surface radiates what it takes, `taken` W/m², a black body's
/// radiation's `emission` share of it getting out; never below the cosmic background.
pub fn radiative_temperature(taken: f64, emission: f64) -> f64 {
    (taken.max(0.0) / (emission * SIGMA_SB))
        .powf(0.25)
        .max(COSMIC_BACKGROUND_K)
}

/// The year's mean climate on the climate grid.
#[derive(Clone, Debug)]
pub struct Climate {
    /// Surface temperature, K.
    pub temperature: Field<f32>,
    /// The sunlight absorbed, W/m².
    pub absorbed: Field<f32>,
    /// The albedo at the year's mean temperature.
    pub albedo: Field<f32>,
    /// The albedo without its ice and snow (F244).
    pub ground_albedo: Field<f32>,
    /// The albedo of ice and snow under the star's light over the ground without its caps
    /// (F250), which they brighten only where it is the brighter.
    pub snow_albedo: Field<f32>,
    /// The heat the atmosphere brings the cell, W/m² (negative where it carries heat away).
    pub transport: Field<f32>,
    /// Newton's steps to converge.
    pub iterations: usize,
}

/// The ground under the climate (F236): what its days and years read.
#[derive(Clone, Copy, Debug)]
pub struct Ground {
    /// Thermal inertia, J m⁻² K⁻¹ s^-1/2.
    pub inertia: f64,
    /// Volumetric heat capacity, J m⁻³ K⁻¹.
    pub heat: f64,
    /// The heat capacity of the air over a square metre, J m⁻² K⁻¹: its column's, which
    /// the cycle's air holds (F245).
    pub air: f64,
}

impl Ground {
    pub fn from_graph(graph: &ParamGraph) -> Self {
        let empty = std::collections::BTreeMap::new();
        let composition = graph
            .get("atmosphere.composition")
            .as_map()
            .unwrap_or(&empty);
        let air = if graph.f64("atmosphere.surface_pressure") > 0.0 {
            mean_heat_capacity(composition) * graph.f64("atmosphere.column_mass")
        } else {
            0.0
        };
        Self {
            inertia: graph.f64("surface.thermal_inertia"),
            heat: regolith_heat(graph.f64("cratering.target_ice")),
            air,
        }
    }
}

/// The frost of the atmosphere's bulk gas (F239): where the ground would cool below its frost
/// point at its partial pressure, the gas condenses on it, and its latent heat holds the
/// ground at the frost point while it lies (Leighton and Murray 1966).
#[derive(Clone, Copy, Debug)]
pub struct Frost {
    pub species: &'static str,
    /// The frost point at the gas's partial pressure, K.
    pub point: f64,
    /// The latent heat of sublimation, J/kg.
    pub latent: f64,
    /// The frost's albedo in the energy balance (F257) and its colour.
    pub albedo: f64,
    pub colour: [f32; 3],
    /// The gas's mean partial pressure, Pa, and its column over the body, kg/m².
    pub partial: f64,
    pub column: f64,
    /// The gas's vapour curve.
    pub volatile: &'static crate::constants::Volatile,
}

impl Frost {
    /// The condensable gas with the largest partial pressure below its triple point, which
    /// alone has the mass to cover the ground; none without an atmosphere, or on a body whose
    /// frost its own stage lays (Io's).
    pub fn from_graph(graph: &ParamGraph) -> Option<Self> {
        let pressure = graph.f64("atmosphere.surface_pressure");
        if pressure <= 0.0 || crate::stages::io::is_io_regime(graph) {
            return None;
        }
        let empty = std::collections::BTreeMap::new();
        let composition = graph
            .get("atmosphere.composition")
            .as_map()
            .unwrap_or(&empty);
        let total: f64 = composition.values().sum::<f64>().max(1e-12);
        composition
            .iter()
            .filter_map(|(name, share)| {
                let v = volatile(name)?;
                let partial = pressure * share / total;
                let colour = frost_colour(name)?;
                (partial > 0.0 && partial < v.p_triple).then_some((partial, v, colour))
            })
            .max_by(|a, b| a.0.total_cmp(&b.0))
            .map(|(partial, v, colour)| {
                let columns = species_columns(composition, graph.f64("atmosphere.column_mass"));
                Self {
                    species: v.name,
                    point: saturation_temperature(v, partial),
                    latent: v.l_sub,
                    albedo: frost_albedo(v.name).unwrap_or(0.0),
                    colour,
                    partial,
                    column: columns.get(v.name).copied().unwrap_or(0.0),
                    volatile: v,
                }
            })
    }

    /// The frost's emissivity at a temperature and its slope per kelvin (F247): nitrogen's,
    /// its alpha phase's below its transition and its beta phase's above, the step blended
    /// smoothly across it; every other frost's a black body's.
    pub fn emissivity(&self, t: f64) -> (f64, f64) {
        if self.species != "N2" {
            return (1.0, 0.0);
        }
        let (alpha, beta) = N2_EMISSIVITY;
        let x = ((t - N2_PHASE_TRANSITION_K) / N2_PHASE_HALF_WIDTH_K).clamp(-1.0, 1.0);
        let step = (2.0 + 3.0 * x - x * x * x) / 4.0;
        let slope = 0.75 * (1.0 - x * x) / N2_PHASE_HALF_WIDTH_K;
        (alpha + (beta - alpha) * step, (beta - alpha) * slope)
    }

    /// The frost point at a partial pressure of the gas, K; none at none.
    pub fn point_at(&self, partial: f64) -> f64 {
        if partial <= 0.0 {
            0.0
        } else {
            saturation_temperature(self.volatile, partial)
        }
    }
}

/// A sea's mixed layer and the ice that forms on it (F244): the layer's heat capacity,
/// J/m²/K; the sea's freezing point, K; its ice's latent heat of fusion, J/kg, density,
/// kg/m³, and conductivity, W/m/K; the heat the sea brings the ice's base, W/m²; and the
/// ice's colour, drawn over the sea.
#[derive(Clone, Copy, Debug)]
pub struct Ocean {
    pub capacity: f64,
    pub freezing: f64,
    pub latent: f64,
    pub density: f64,
    pub conductivity: f64,
    pub base_flux: f64,
    pub colour: [f32; 3],
}

impl Ocean {
    /// The body's sea where its working liquid is water (F244): its mixed layer
    /// OCEAN_MIXED_LAYER_M deep, freezing SEA_SALT_DEPRESSION_K under pure water's point as
    /// Earth's salt water does; the ices of other seas wait on their own constants.
    pub fn from_graph(graph: &ParamGraph) -> Option<Self> {
        if graph.str("volatiles.working_liquid") != "water" {
            return None;
        }
        let water = volatile("H2O")?;
        Some(Self {
            capacity: water.liquid_density * SEA_WATER_HEAT_CAPACITY * OCEAN_MIXED_LAYER_M,
            freezing: water.t_triple - SEA_SALT_DEPRESSION_K,
            latent: WATER_FUSION_HEAT,
            density: ICE_DENSITIES[0].1,
            conductivity: SEA_ICE_CONDUCTIVITY,
            base_flux: SEA_ICE_BASE_FLUX,
            colour: SEA_ICE_COLOUR,
        })
    }
}

/// A sea's year (F244, F245), on the climate's grid: its surface's mean through the year,
/// its coldest and warmest, K, its ice at the snapshot's season, kg/m², and the share of the
/// year ice lies on it; the land's temperatures where a place of the cycle has no sea.
#[derive(Clone, Debug)]
pub struct Seas {
    pub mean: Field<f32>,
    pub min: Field<f32>,
    pub max: Field<f32>,
    pub ice_now: Field<f32>,
    pub ice_share: Field<f32>,
    pub ocean: Ocean,
    /// The surface's mean, coldest and warmest at each place of the cycle's layout, each
    /// place without a sea given the mean of the places with one beside it, spreading from the
    /// sea (F273); the ice at the snapshot's season and the share of the sea it covers at
    /// each place; and the layout's latitudes and longitudes, a fast spinner's one.
    mean_rows: Vec<f64>,
    min_rows: Vec<f64>,
    max_rows: Vec<f64>,
    ice_rows: Vec<f64>,
    cover_rows: Vec<f64>,
    layout: (usize, usize),
    /// Where the ice at the snapshot's season ends between the places (F251).
    guide: Guide,
    margin: Margin,
}

impl Seas {
    /// The sea's ice at the snapshot's season in a direction, kg/m², read along the cycle's
    /// latitudes, a resonant body's longitudes too; between places under ice and open ones,
    /// thinning to its edge where the margin sets it (F251).
    pub fn ice_at(&self, d: [f64; 3]) -> f64 {
        match self.edge(d) {
            Some((depth, corners)) => depth * self.among(&self.ice_rows, &corners, true),
            None => self.read(&self.ice_rows, d),
        }
    }

    /// The sea's surface's mean through the year in a direction, K, read smoothly between the
    /// places (F273): a sea under ice part of its year stands between its ice's top and its
    /// open water through it, with no edge drawn through its year's mean.
    pub fn mean_at(&self, d: [f64; 3]) -> f64 {
        self.smooth(&self.mean_rows, d)
    }

    /// The sea's surface's coldest and warmest through the year in a direction, K, read as
    /// its mean is (F253, F273).
    pub fn min_at(&self, d: [f64; 3]) -> f64 {
        self.smooth(&self.min_rows, d)
    }

    pub fn max_at(&self, d: [f64; 3]) -> f64 {
        self.smooth(&self.max_rows, d)
    }

    /// `rows` read smoothly between the places in a direction (F273).
    fn smooth(&self, rows: &[f64], d: [f64; 3]) -> f64 {
        let (lon, lat) = lon_lat_from_direction(d);
        crate::insolation::smooth_grid(rows, self.layout.1, lat, lon)
    }

    /// The share of the sea its ice covers at the snapshot's season in a direction (F245):
    /// each place's, 1 - exp(-m / SEA_ICE_COVER_MASS), read between the places; between
    /// places under ice and open ones, from none at the ice's edge (F251) to the share at the
    /// places under ice, as a marginal ice zone spans the distance.
    pub fn cover_at(&self, d: [f64; 3]) -> f64 {
        match self.edge(d) {
            Some((depth, corners)) => depth * self.among(&self.cover_rows, &corners, true),
            None => self.read(&self.cover_rows, d),
        }
    }

    /// `rows` read between the places about a direction that have a sea (F266).
    fn read(&self, rows: &[f64], d: [f64; 3]) -> f64 {
        self.corners(d).iter().map(|&(p, w)| w * rows[p]).sum()
    }

    /// The places about a direction with their weights, those with a sea among them (F266).
    fn corners(&self, d: [f64; 3]) -> [(usize, f64); 4] {
        let (lon, lat) = lon_lat_from_direction(d);
        self.margin
            .ground(&reader(self.layout.0, self.layout.1, lat, lon))
    }

    /// How far into the ice a direction stands where the places about it differ (F251), and
    /// those places with their weights.
    fn edge(&self, d: [f64; 3]) -> Option<(f64, [(usize, f64); 4])> {
        let corners = self.corners(d);
        let depth = self.margin.depth(&corners, self.guide.at(d))?;
        Some((depth, corners))
    }

    /// The weighted mean of `rows` over the places a direction reads the ice holds, `held`,
    /// or leaves open (`Margin::kind`).
    fn among(&self, rows: &[f64], corners: &[(usize, f64); 4], held: bool) -> f64 {
        let (sum, weight) = self
            .margin
            .kind(corners, held)
            .iter()
            .fold((0.0, 0.0), |(s, t), &(p, w)| (s + w * rows[p], t + w));
        if weight > 0.0 { sum / weight } else { 0.0 }
    }
}

/// `rows`, a value at each place of a layout of `n_lat` latitudes and `n_lon` longitudes, held
/// only where `own`, the rest given the mean of their neighbours that hold or were given one,
/// spreading from the places that hold it until all hold one (F273): so that read smoothly
/// about a coast a sea's field reads the sea's.
pub(crate) fn spread_from(rows: &[f64], own: &[bool], (n_lat, n_lon): (usize, usize)) -> Vec<f64> {
    let n = n_lat * n_lon;
    let mut out = rows.to_vec();
    let mut known = own.to_vec();
    if !known.iter().any(|&k| k) {
        return out;
    }
    let neighbours = |p: usize| {
        let (i, j) = (p / n_lon, p % n_lon);
        let mut list = Vec::with_capacity(4);
        if i > 0 {
            list.push(p - n_lon);
        }
        if i + 1 < n_lat {
            list.push(p + n_lon);
        }
        if n_lon > 1 {
            list.push(i * n_lon + (j + n_lon - 1) % n_lon);
            list.push(i * n_lon + (j + 1) % n_lon);
        }
        list
    };
    while known.iter().any(|&k| !k) {
        let next: Vec<Option<f64>> = (0..n)
            .map(|p| {
                if known[p] {
                    return None;
                }
                let (sum, count) = neighbours(p)
                    .into_iter()
                    .filter(|&q| known[q])
                    .fold((0.0, 0), |(s, c), q| (s + out[q], c + 1));
                (count > 0).then(|| sum / count as f64)
            })
            .collect();
        if next.iter().all(Option::is_none) {
            break;
        }
        for (p, v) in next.into_iter().enumerate() {
            if let Some(v) = v {
                (out[p], known[p]) = (v, true);
            }
        }
    }
    out
}

/// The pattern the year's sunlight sets over the body (F251): the diffusive balance's
/// temperature at each cosine of the angle from its axis (`sunlight_pattern`), the point a
/// body locked to its star centres its sunlight on, any other's spin axis: smooth everywhere
/// and as round about its axis as the sunlight. It guides where an ice's edge falls between
/// the cycle's places.
#[derive(Clone, Debug)]
struct Guide {
    profile: Vec<f64>,
    axis: [f64; 3],
}

impl Guide {
    fn of(balance: &EnergyBalance, climate: &Climate) -> Guide {
        let locked = balance.sky.resonance == Some((1, 1));
        let axis = if locked {
            let grid = climate.absorbed.grid();
            let sum = grid.cells().fold([0.0; 3], |acc, c| {
                let (d, q) = (grid.center(c), climate.absorbed.get(c) as f64);
                let w = q * grid.cell_area(c);
                [acc[0] + w * d[0], acc[1] + w * d[1], acc[2] + w * d[2]]
            });
            let norm = (sum[0] * sum[0] + sum[1] * sum[1] + sum[2] * sum[2]).sqrt();
            if norm > 0.0 {
                [sum[0] / norm, sum[1] / norm, sum[2] / norm]
            } else {
                [1.0, 0.0, 0.0]
            }
        } else {
            [0.0, 0.0, 1.0]
        };
        let e = balance.sky.eccentricity;
        let flux = balance.sky.flux / (1.0 - e * e).max(1e-9).sqrt();
        let tilt = balance.sky.obliquity;
        let absorbed = |mu: f64| {
            let q = if locked {
                flux * mu.max(0.0)
            } else {
                flux * crate::insolation::zonal_annual_share(mu.clamp(-1.0, 1.0).asin(), tilt)
            };
            (1.0 - balance.albedo) * q
        };
        let (_, slope) = balance.olr_and_slope(balance.start.max(1.0));
        let profile = sunlight_pattern((slope, balance.transport), absorbed);
        Guide { profile, axis }
    }

    fn at(&self, d: [f64; 3]) -> f64 {
        let mu = (d[0] * self.axis[0] + d[1] * self.axis[1] + d[2] * self.axis[2]).clamp(-1.0, 1.0);
        let n = self.profile.len();
        let x = ((mu + 1.0) / 2.0 * n as f64 - 0.5).clamp(0.0, (n - 1) as f64);
        let i = (x.floor() as usize).min(n - 2);
        let f = x - i as f64;
        self.profile[i] * (1.0 - f) + self.profile[i + 1] * f
    }
}

/// Where an ice's edge falls between the cycle's places (F251). The cycle tells, place by
/// place, whether the ice holds; between a place it holds and one beside it it does not, its
/// edge stands where the guide crosses the midpoint of the two places' guides, smoothed along
/// the edge over MARGIN_SMOOTHING spacings of the places, so that it runs along the guide's
/// isoline round its regions rather than along the layout's steps, as the caps' edges follow
/// their traps read smoothly between the climate's cells (F241). The ice thins to its edge
/// over half the guide's step between the two places. F266: each edge counts by the share of
/// its ground, the sea's or the land's, the less of its two places hold, in the smoothing and
/// in where it stands: a place stands at the edge by the most such share with a place of the
/// other kind beside it, and a point reads the guided edge by how much of the places about it
/// stand at it, the places' own ice for the rest; so that a coastal place's sliver of sea,
/// frozen by the land about it, neither pins the edge to the places' steps nor ends the
/// guide's reach along a line between them.
#[derive(Clone, Debug)]
struct Margin {
    /// Whether the ice holds at each place, none where the place has none of its ground.
    held: Vec<Option<bool>>,
    /// How far each place stands at the ice's edge, 0 to 1: the most share of its ground it
    /// and a place of the other kind beside it both hold.
    beside: Vec<f64>,
    /// At each place near the edge: the guide the edge stands at and the guide over which the
    /// ice thins to it, smoothed along it, and the nearest place under ice and open place;
    /// NaN and none elsewhere.
    edge: Vec<f64>,
    thinning: Vec<f64>,
    near: Vec<Option<(usize, usize)>>,
}

impl Margin {
    fn new(
        (n_lat, n_lon): (usize, usize),
        held: Vec<Option<bool>>,
        shares: &[f64],
        guide: &Guide,
    ) -> Margin {
        let n = n_lat * n_lon;
        let lat_of = |p: usize| -PI / 2.0 + PI * (p / n_lon) as f64 / (n_lat - 1) as f64;
        let lon_of = |p: usize| -PI + TAU * (p % n_lon) as f64 / n_lon as f64;
        let direction = |p: usize| direction_from_lon_lat(lon_of(p), lat_of(p));
        // A fast spinner's place is its latitude's: its guide the guide's mean along it.
        let values: Vec<f64> = (0..n)
            .map(|p| match n_lon {
                1 => {
                    let k = SEA_SHARE_LONGITUDES;
                    (0..k)
                        .map(|m| {
                            let lon = -PI + TAU * m as f64 / k as f64;
                            guide.at(direction_from_lon_lat(lon, lat_of(p)))
                        })
                        .sum::<f64>()
                        / k as f64
                }
                _ => guide.at(direction(p)),
            })
            .collect();
        let neighbours = |p: usize| {
            let (i, j) = (p / n_lon, p % n_lon);
            let mut out = Vec::with_capacity(4);
            if i > 0 {
                out.push(p - n_lon);
            }
            if i + 1 < n_lat {
                out.push(p + n_lon);
            }
            if n_lon > 1 {
                out.push(i * n_lon + (j + n_lon - 1) % n_lon);
                out.push(i * n_lon + (j + 1) % n_lon);
            }
            out
        };
        // The edge between each place under ice and each open place beside it whose guide is
        // the warmer.
        let mut links: Vec<Link> = Vec::new();
        let mut beside = vec![0.0f64; n];
        for p in 0..n {
            let Some(here) = held[p] else { continue };
            for q in neighbours(p) {
                if held[q] != Some(!here) {
                    continue;
                }
                beside[p] = beside[p].max(shares[p].min(shares[q]).max(0.0));
                if !here {
                    continue;
                }
                let (cold, warm) = (p, q);
                let step = values[warm] - values[cold];
                if step <= 0.0 {
                    continue;
                }
                let (a, b) = (direction(cold), direction(warm));
                let mid = [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
                let norm = (mid[0] * mid[0] + mid[1] * mid[1] + mid[2] * mid[2]).sqrt();
                let lat = 0.5 * (lat_of(cold) + lat_of(warm));
                let at = [mid[0] / norm, mid[1] / norm, mid[2] / norm];
                links.push(Link {
                    guide: values[cold] + 0.5 * step,
                    half: 0.5 * step,
                    at,
                    lat,
                    places: (cold, warm),
                    share: shares[cold].min(shares[warm]).max(0.0),
                });
            }
        }
        // Smoothed along the edge: a fast spinner's by latitude alone.
        let width = MARGIN_SMOOTHING * PI / (n_lat - 1) as f64;
        let mut edge = vec![f64::NAN; n];
        let mut thinning = vec![f64::NAN; n];
        let mut near = vec![None; n];
        for p in 0..n {
            let here = direction(p);
            let (mut sum, mut thin, mut total, mut best) = (0.0, 0.0, 0.0, (0.0, None));
            for l in &links {
                let angle = match n_lon {
                    1 => (lat_of(p) - l.lat).abs(),
                    _ => {
                        let at = l.at;
                        let dot = here[0] * at[0] + here[1] * at[1] + here[2] * at[2];
                        dot.clamp(-1.0, 1.0).acos()
                    }
                };
                let w = l.share * (-(angle / width).powi(2)).exp();
                (sum, thin, total) = (sum + w * l.guide, thin + w * l.half, total + w);
                if w > best.0 {
                    best = (w, Some(l.places));
                }
            }
            if total > MARGIN_REACH {
                (edge[p], thinning[p], near[p]) = (sum / total, thin / total, best.1);
            }
        }
        Margin {
            held,
            beside,
            edge,
            thinning,
            near,
        }
    }

    /// The places about a point on the margin's ground, the sea's or the land's, with their
    /// weights over their sum (F266): by a coast a point reads its own ground's places, not
    /// the other's, which hold none of its ice; all of them where none is on its ground.
    fn ground(&self, corners: &[(usize, f64); 4]) -> [(usize, f64); 4] {
        let total: f64 = corners
            .iter()
            .filter(|&&(p, w)| w > 0.0 && self.held[p].is_some())
            .map(|c| c.1)
            .sum();
        if total <= 0.0 {
            return *corners;
        }
        corners.map(|(p, w)| {
            (
                p,
                if self.held[p].is_some() {
                    w / total
                } else {
                    0.0
                },
            )
        })
    }

    /// How far into the ice a point stands, its guide `at`, among `corners`, places with
    /// their weights: by the guided edge, 0 where its guide stands past the edge's, rising to
    /// 1 half the guide's step into the ice, as far as the places stand at the edge (F266),
    /// and by the share of them the ice holds at for the rest; none where they are all of a
    /// kind and none stands at the edge.
    fn depth(&self, corners: &[(usize, f64); 4], at: f64) -> Option<f64> {
        let live = || {
            corners
                .iter()
                .filter(|&&(p, w)| w > 0.0 && self.held[p].is_some())
        };
        let total: f64 = live().map(|c| c.1).sum();
        if total <= 0.0 {
            return None;
        }
        let standing = live().map(|&(p, w)| w * self.beside[p]).sum::<f64>() / total;
        let holds = live()
            .filter(|&&(p, _)| self.held[p] == Some(true))
            .map(|c| c.1)
            .sum::<f64>()
            / total;
        if standing <= 0.0 && (holds <= 0.0 || holds >= 1.0) {
            return None;
        }
        let (mut edge, mut thin, mut weight) = (0.0, 0.0, 0.0);
        for &(p, w) in live().filter(|&&(p, _)| self.edge[p].is_finite()) {
            (edge, thin, weight) = (
                edge + w * self.edge[p],
                thin + w * self.thinning[p],
                weight + w,
            );
        }
        let guided = if weight > 0.0 {
            ((edge - at * weight) / thin.max(1e-300)).clamp(0.0, 1.0)
        } else {
            holds
        };
        Some(standing * guided + (1.0 - standing) * holds)
    }

    /// The places of a kind a point between places reads, `held` or open, with their
    /// weights: the corners of that kind, or, where none is, the nearest such places to the
    /// corners along the edge.
    fn kind(&self, corners: &[(usize, f64); 4], held: bool) -> Vec<(usize, f64)> {
        let own: Vec<(usize, f64)> = corners
            .iter()
            .copied()
            .filter(|&(p, w)| w > 0.0 && self.held[p] == Some(held))
            .collect();
        if !own.is_empty() {
            return own;
        }
        corners
            .iter()
            .filter(|c| c.1 > 0.0)
            .filter_map(|&(p, w)| {
                self.near[p].map(|(cold, warm)| (if held { cold } else { warm }, w))
            })
            .collect()
    }
}

/// An edge between a place under ice and an open place beside it (F251): its guide, half the
/// guide's step between them, where it stands and its latitude, the two places, cold and
/// warm, and the share of its ground the less of them holds (F266).
#[derive(Clone, Copy, Debug)]
struct Link {
    guide: f64,
    half: f64,
    at: [f64; 3],
    lat: f64,
    places: (usize, usize),
    share: f64,
}

/// A sea's mixed layer at a place (F244): the layer's temperature, K, its ice, kg/m², and
/// the surface's temperature, the ice's top where it lies.
#[derive(Clone, Copy, Debug)]
struct Sea {
    t: f64,
    ice: f64,
    surface: f64,
}

/// A volatile's ice that can lie on the ground (F240): its gas, its colour and the albedo
/// that is its colour's mean.
#[derive(Clone, Copy, Debug)]
pub struct Ice {
    pub species: &'static str,
    pub colour: [f32; 3],
    pub volatile: &'static crate::constants::Volatile,
}

impl Ice {
    /// The inventory's ices the model draws (F240): carbon dioxide's, nitrogen's and
    /// methane's in their frosts' colours (F239), water's in `water`, clean water ice's,
    /// except where the water is an icy body's own shell.
    pub fn of_inventory(graph: &ParamGraph, water: [f32; 3]) -> Vec<Ice> {
        let empty = std::collections::BTreeMap::new();
        let inventory = graph.get("volatiles.inventory").as_map().unwrap_or(&empty);
        let shell = is_icy(graph.str("body.composition_class"));
        inventory
            .iter()
            .filter(|&(name, &m)| m > 0.0 && !(shell && name == "H2O"))
            .filter_map(|(name, _)| {
                let v = volatile(name)?;
                let colour = if name == "H2O" {
                    water
                } else {
                    frost_colour(name)?
                };
                Some(Ice {
                    species: v.name,
                    colour,
                    volatile: v,
                })
            })
            .collect()
    }

    /// The ice's albedo in the energy balance: a frost's apart from its colour (F257), water
    /// ice's its colour's mean.
    pub fn albedo(&self) -> f64 {
        frost_albedo(self.species).unwrap_or((self.colour.iter().sum::<f32>() / 3.0) as f64)
    }
}

/// Where each ice of the inventory would last (F240), on the climate's grid.
#[derive(Clone, Debug)]
pub struct Traps {
    /// The frost's own gas's ice, lying all year from an unlimited store: the frost it gains
    /// in a year, kg/m², positive where a cap would grow (Leighton and Murray 1966).
    pub growth: Option<(Ice, Field<f32>)>,
    /// Every other ice: the year-mean saturation pressure over it, brought to the datum's air
    /// by the place's height, Pa; the lower, the colder the trap, and the coldest's holds the
    /// air's vapour (Ingersoll 1970).
    pub vapour: Vec<(Ice, Field<f32>)>,
    /// Each of `vapour`'s ices' warmest through the year, K (F241).
    pub warmest: Vec<Field<f32>>,
}

impl Traps {
    /// The climate's grid the traps stand on.
    pub fn grid(&self) -> CubeSphere {
        self.growth
            .as_ref()
            .map(|g| g.1.grid())
            .or_else(|| self.vapour.first().map(|v| v.1.grid()))
            .unwrap_or_else(|| CubeSphere::new(CLIMATE_FACE_SIZE))
    }

    /// An ice's coldest trap: its latitude and longitude, degrees, and the saturation
    /// pressure it holds the datum's air at, Pa.
    pub fn coldest(&self, species: &str) -> Option<(f64, f64, f64)> {
        let (_, field) = self.vapour.iter().find(|(ice, _)| ice.species == species)?;
        let grid = field.grid();
        let (k, &v) = field
            .as_slice()
            .iter()
            .enumerate()
            .filter(|(_, v)| v.is_finite())
            .min_by(|a, b| a.1.total_cmp(b.1))?;
        let (lon, lat) = lon_lat_from_direction(grid.center(grid.cell(k)));
        Some((lat.to_degrees(), lon.to_degrees(), v as f64))
    }
}

/// The frost's own gas's caps on the climate's grid, the reservoir its air draws on (F258):
/// the share of each cell they cover, the ice they hold over that share, kg/m², and the albedo
/// of the ground about them, under the other ices' caps.
#[derive(Clone, Debug)]
pub struct Reservoir {
    pub share: Field<f32>,
    pub mass: Field<f32>,
    pub ground: Field<f32>,
}

/// A cell of the climate's grid in the traps (F240): its latitude and longitude, rad, the
/// share of the datum's pressure its height keeps, its year's mean temperature, K, the
/// places of the cycle's layout its air is read from, with their weights, and how much
/// colder than the datum's its air stands at its height, K (F245).
#[derive(Clone, Copy, Debug)]
struct TrapCell {
    lat: f64,
    lon: f64,
    thin: f64,
    steady: f64,
    air: [(usize, f64); 4],
    cooling: f64,
}

impl TrapCell {
    /// The frost at a step of its cycle at the cell's partial pressure, `partials` the
    /// datum's.
    fn laid(&self, frost: Option<&Frost>, partials: &[f64], k: usize) -> Option<Frost> {
        frost.map(|f| Frost {
            point: f.point_at(partials.get(k).copied().unwrap_or(f.partial) * self.thin),
            ..*f
        })
    }

    /// The air over the cell at a step of its cycle, K, read from the cycle's air at a tilt
    /// over its layout's `places`, at the cell's height; the cell's own year's temperature on
    /// an airless body.
    fn air_at(&self, air: &[f64], places: usize, k: usize) -> f64 {
        if air.is_empty() {
            return self.steady;
        }
        self.air
            .iter()
            .map(|&(p, w)| w * air[k * places + p])
            .sum::<f64>()
            - self.cooling
    }
}

/// The colour of a gas's frost, linear reflectance (F239); none for a gas whose frost the
/// model does not draw.
fn frost_colour(species: &str) -> Option<[f32; 3]> {
    match species {
        "CO2" => Some(CO2_FROST_COLOUR),
        "N2" => Some(N2_FROST_COLOUR),
        "CH4" => Some(CH4_FROST_COLOUR),
        _ => None,
    }
}

/// A gas's frost's albedo in the energy balance (F257), apart from its colour; none for a gas
/// whose frost the model does not draw.
fn frost_albedo(species: &str) -> Option<f64> {
    match species {
        "CO2" => Some(CO2_FROST_ALBEDO),
        "N2" => Some(N2_FROST_ALBEDO),
        "CH4" => Some(CH4_FROST_ALBEDO),
        _ => None,
    }
}

/// The year's surface temperatures through the days and seasons (F236), K, the frost (F239)
/// and the sea (F244), under the air that carries the transport through them (F245).
#[derive(Clone, Debug)]
pub struct Cycle {
    /// The mean over the year: the land's, where a place has land and sea.
    pub mean: Field<f32>,
    /// The warmest the surface stands.
    pub max: Field<f32>,
    /// The coldest.
    pub min: Field<f32>,
    /// The share of the year frost lies on the ground, 1 under a residual cap.
    pub frost_share: Field<f32>,
    /// The frost on the ground at the snapshot's season, kg/m².
    pub frost_now: Field<f32>,
    /// The frost's gas, if the atmosphere has one.
    pub frost: Option<Frost>,
    /// The surface pressure's lowest and highest through the year as the frost takes and
    /// gives the gas back, Pa, where it lays frost.
    pub pressure_range: Option<(f64, f64)>,
    /// The frost that outlasts the year, kg/m² a year over the body: its caps' growth, which
    /// the caps keep (F241), where it lays frost.
    pub cap_growth: Option<f64>,
    /// What the frost's own gas's caps gain through the year where they hold its air (F258),
    /// kg/m² a year, over each place's area on the climate's grid and over the body: more
    /// than nothing where the year's pressure stands above the one they would hold it at.
    pub reservoir: Option<(Field<f32>, f64)>,
    /// The tilts the year was run at, and at each the frost's gas's partial pressure at the
    /// datum through the year's steps, Pa, where it lays frost (F240).
    pub tilts: Vec<f64>,
    pub partials: Vec<Vec<f64>>,
    /// The sea's year, where the body has one (F244).
    pub seas: Option<Seas>,
    /// The ices that come and go through the year, under an atmosphere (F246).
    pub cryosphere: Option<Cryosphere>,
    /// The air over each place of the cycle's layout through its last cycle at each tilt,
    /// step after step, K, none on an airless body; and the layout's latitudes and
    /// longitudes (F245).
    air: Vec<Vec<f64>>,
    layout: (usize, usize),
    /// The cycle's year step by step and its clock, under an atmosphere (F253).
    pub series: Option<Series>,
    pub clock: Option<Clock>,
    /// The wind through the year at the cycle's places, under an atmosphere (F267), and the
    /// seas' currents it drives (F272).
    pub wind: Option<crate::wind::Wind>,
    pub currents: Option<crate::currents::Currents>,
    /// The cycle's places under an atmosphere, as the wind and the moisture read them, and
    /// each place's frost, its sea's ice and its land's snow (F271), kg/m², and its land's
    /// surface and its sea's, K, at the end of each step, place after place, the tilts' mean
    /// (F268); and the snow its land's melt took through each step, kg/m², NaN where its year's
    /// snow is unknown (F277).
    pub sites: Vec<crate::wind::Site>,
    pub masses: [Vec<f32>; 4],
    pub surfaces: [Vec<f32>; 2],
    /// Each place's surface temperatures through the cycle, as quantiles, and the albedo its
    /// ground stood at (F256).
    spread: Vec<[f32; CYCLE_SPREAD]>,
    spread_albedo: Vec<f32>,
}

/// The cycle's year step by step (F253): each step's end's phase of the year from the
/// northern spring equinox, the star's longitude seen from the body, its declination and
/// its flux there, the surface pressure, and the frost's and the sea's ice's mean over each
/// hemisphere, north then south, kg/m² of the hemisphere.
#[derive(Clone, Debug, Default)]
pub struct Series {
    pub phase: Vec<f64>,
    pub ls_deg: Vec<f64>,
    pub declination_deg: Vec<f64>,
    pub flux: Vec<f64>,
    pub pressure: Vec<f64>,
    pub frost: [Vec<f64>; 2],
    pub sea_ice: [Vec<f64>; 2],
    pub snow: [Vec<f64>; 2],
}

/// The clock the cryosphere's phases run on (F253): the cycle they count through, s, the
/// year or a resonant body's whole repeat, its orbits and its steps; the solar day, s, none
/// on a body locked to its star; the orbit's tilt, eccentricity and the star's longitude at
/// perihelion; and the phases of the northern spring equinox, summer solstice, autumn
/// equinox and winter solstice in the cycle's first orbit.
#[derive(Clone, Copy, Debug)]
pub struct Clock {
    pub cycle_s: f64,
    pub orbits: u32,
    pub steps: usize,
    pub solar_day_s: Option<f64>,
    pub obliquity_deg: f64,
    pub eccentricity: f64,
    pub perihelion_ls_deg: f64,
    pub seasons: [f64; 4],
}

/// An ice's window through the year at a place (F246): the phases of the year, its share
/// from the northern spring equinox, it freezes at and melts out at, 0 and 1 where it lies
/// all year, its peak, kg/m² (m once brought to a grid), and the phase of its peak; all none
/// where it never lies.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Window {
    pub freeze: f64,
    pub melt: f64,
    pub peak: f64,
    pub at: f64,
}

impl Window {
    /// Whether the ice lies all year.
    pub fn perennial(&self) -> bool {
        self.peak > 0.0 && self.freeze == 0.0 && self.melt == 1.0
    }

    /// The ice's peak at a phase of the year as an engine draws it from the cryosphere's
    /// maps (F253): none outside its window, from its freezing round to its melting; within
    /// it its peak times (t / rise)^a, t the time since it froze, rise the time from its
    /// freezing to its peak's phase, then 1 - (g / fall)^b, g the time since its peak, fall
    /// the time from its peak to its melting, `curve` (a, b); its peak all year where it lies
    /// all year. In the units of its peak.
    pub fn at_phase(&self, phase: f64, (rise, fall): (f64, f64)) -> f64 {
        if self.peak <= 0.0 {
            return 0.0;
        }
        if self.perennial() {
            return self.peak;
        }
        let span = (self.melt - self.freeze).rem_euclid(1.0);
        let since = (phase - self.freeze).rem_euclid(1.0);
        if since >= span {
            return 0.0;
        }
        let up = (self.at - self.freeze).rem_euclid(1.0).min(span);
        if since <= up && up > 0.0 {
            self.peak * (since / up).powf(rise)
        } else {
            let down = (span - up).max(f64::MIN_POSITIVE);
            self.peak * (1.0 - ((since - up) / down).clamp(0.0, 1.0).powf(fall))
        }
    }

    /// The window an ice lies in through a cycle of `masses`, one at the end of each step
    /// (F246), a mass past `cover` counting as lying, as half the ground's colour hidden:
    /// the longest run of steps it lies through, which it freezes into and melts out of
    /// where its mass crosses `cover`, between the steps' ends (F255); its peak the most it
    /// holds in that run. `phase` takes a time in steps from the cycle's start to the phase
    /// of the year.
    fn of(masses: &[f32], cover: f64, phase: impl Fn(f64) -> f64) -> Window {
        match Lying::longest(masses, cover) {
            None => Window::default(),
            Some(run) if run.len == masses.len() => Window {
                freeze: 0.0,
                melt: 1.0,
                peak: masses[run.peak] as f64,
                at: phase(run.peak as f64 + 1.0),
            },
            Some(run) => {
                let (opens, closes) = run.crossings(masses, cover);
                Window {
                    freeze: phase(run.start as f64 + opens),
                    melt: phase(run.start as f64 + closes),
                    peak: masses[run.peak] as f64,
                    at: phase(run.peak as f64 + 1.0),
                }
            }
        }
    }
}

/// An ice that comes and goes as the cryosphere's maps draw it (F253): its name and its
/// species, whether it lies on the sea, its colour, linear reflectance, its density, kg/m³,
/// the mass that hides 1 - 1/e of what lies under it, kg/m², and its year curve.
#[derive(Clone, Debug)]
pub struct SeasonalIce {
    pub name: String,
    pub species: &'static str,
    pub sea: bool,
    pub colour: [f32; 3],
    pub density: f64,
    pub cover_mass: f64,
    pub curve: Option<(f64, f64)>,
}

/// The longest run of a cycle's steps an ice lies through (F246), its mass past `cover` at
/// each step's end: its first step, its length and the step of its peak, from a step it does
/// not lie in round the cycle; all the cycle where it never leaves.
#[derive(Clone, Copy, Debug)]
struct Lying {
    start: usize,
    len: usize,
    peak: usize,
}

impl Lying {
    /// When the ice's mass crosses `cover` on its way into the run and out of it, in steps
    /// from the run's first step's start (F255): the masses stand at their steps' ends, so
    /// that it opens within its first step and closes within the step after its last, where
    /// the mass between them crosses. A step that ends with none of the ice tells nothing of
    /// when in it the ice came or went, a step taking no more than lies: there the ice comes
    /// or goes no slower than through the run's step beside it, so that the places whose ice
    /// runs out within one step do not all close at its end.
    fn crossings(&self, masses: &[f32], cover: f64) -> (f64, f64) {
        let n = masses.len();
        let m = |k: usize| masses[k % n] as f64;
        let (first, last) = (self.start, self.start + self.len - 1);
        let (before, after) = (m(first + n - 1), m(last + 1));
        let (gain, loss) = (m(first) - before, m(last) - after);
        let steady = |by: f64, beside: f64, none: bool| match none && self.len > 1 {
            true => by.max(beside),
            false => by,
        };
        let gain_in = steady(gain, m(first + 1) - m(first), before <= 0.0);
        let loss_out = steady(loss, m(last + n - 1) - m(last), after <= 0.0);
        let opens = match gain_in > gain {
            true => 1.0 - (m(first) - cover) / gain_in,
            false => (cover - before) / gain.max(f64::MIN_POSITIVE),
        }
        .clamp(0.0, 1.0);
        let closes =
            self.len as f64 + ((m(last) - cover) / loss_out.max(f64::MIN_POSITIVE)).clamp(0.0, 1.0);
        (opens, closes)
    }

    fn longest(masses: &[f32], cover: f64) -> Option<Lying> {
        let n = masses.len();
        let lying: Vec<bool> = masses.iter().map(|&m| m as f64 > cover).collect();
        let count = lying.iter().filter(|&&l| l).count();
        let peak_in = |steps: &mut dyn Iterator<Item = usize>| {
            steps
                .map(|k| (k, masses[k] as f64))
                .max_by(|a, b| a.1.total_cmp(&b.1))
                .map_or(0, |(k, _)| k)
        };
        if count == 0 {
            return None;
        }
        if count == n {
            let peak = peak_in(&mut (0..n));
            return Some(Lying {
                start: 0,
                len: n,
                peak,
            });
        }
        let open = lying.iter().position(|&l| !l).unwrap_or(0);
        let (mut best, mut longest, mut start, mut run) = (0, 0, 0, 0);
        for j in 1..=n {
            let k = (open + j) % n;
            if lying[k] {
                if run == 0 {
                    start = k;
                }
                run += 1;
                if run > longest {
                    (best, longest) = (start, run);
                }
            } else {
                run = 0;
            }
        }
        let peak = peak_in(&mut (0..longest).map(|j| (best + j) % n));
        Some(Lying {
            start: best,
            len: longest,
            peak,
        })
    }
}

/// An ice's year curve (F253): over the places its ice comes and goes at, the exponents of
/// its rise from its window's opening to its peak, its peak times (t / rise)^a, and of its
/// fall to its window's close, its peak times 1 - (g / fall)^b, fitted by least squares to
/// its masses through the cycle, `masses` place after place, `steps` a place, each place's
/// over its peak; none where none comes and goes with a peak past YEAR_CURVE_PEAK times
/// `cover`, the mass it lies at.
fn year_curve(masses: &[f32], steps: usize, cover: f64) -> Option<(f64, f64)> {
    let mut rising: Vec<(f64, f64)> = Vec::new();
    let mut falling: Vec<(f64, f64)> = Vec::new();
    for place in masses.chunks(steps) {
        let Some(run) = Lying::longest(place, cover).filter(|r| r.len < steps) else {
            continue;
        };
        let peak = place[run.peak] as f64;
        if peak < YEAR_CURVE_PEAK * cover {
            continue;
        }
        // Times in steps from the run's first step's start: it opens and closes where its
        // mass crosses the mass it lies at (F255) and peaks at its peak step's end.
        let at = ((run.peak + steps - run.start) % steps) as f64 + 1.0;
        let (open, close) = run.crossings(place, cover);
        for j in 0..run.len {
            let t = j as f64 + 1.0;
            let m = place[(run.start + j) % steps] as f64;
            let y = (m / peak).clamp(0.0, 1.0);
            if t <= at {
                rising.push(((t - open) / (at - open), y));
            }
            if t >= at {
                falling.push(((t - at) / (close - at), y));
            }
        }
    }
    if rising.is_empty() || falling.is_empty() {
        return None;
    }
    let best = |error: &dyn Fn(f64) -> f64| -> f64 {
        let (lo, hi) = YEAR_CURVE_RANGE;
        let per = (1.0 / YEAR_CURVE_STEP).round();
        ((lo * per).round() as i64..=(hi * per).round() as i64)
            .map(|i| i as f64 / per)
            .map(|a| (a, error(a)))
            .min_by(|x, y| x.1.total_cmp(&y.1))
            .map_or(1.0, |(a, _)| a)
    };
    let rise = best(&|a| rising.iter().map(|&(t, y)| (y - t.powf(a)).powi(2)).sum());
    let fall = best(&|b| {
        falling
            .iter()
            .map(|&(g, y)| (y - (1.0 - g.powf(b))).powi(2))
            .sum()
    });
    Some((rise, fall))
}

/// The ices that come and go through the year (F246): at each place of the cycle's layout
/// its land's frost's window and its sea's ice's, and the phase of the snapshot's season.
#[derive(Clone, Debug)]
pub struct Cryosphere {
    land: Vec<Window>,
    sea: Vec<Window>,
    /// The land's snow's windows (F271).
    snow: Vec<Window>,
    layout: (usize, usize),
    /// The phase of the year of the snapshot's season.
    pub snapshot: f64,
    /// Where the sea's ice and the land's frost end between the places through the year
    /// (F251, F252).
    guide: Guide,
    sea_margin: Margin,
    land_margin: Margin,
    /// Each ice's year curve, the land's frost's, the sea's ice's and the land's snow's (F253,
    /// F271).
    pub curves: [Option<(f64, f64)>; 3],
    /// The land's snow's year at each place, run again on finer ground (F273), and the
    /// phase of the year of the northern spring equinox, as the cycle's steps count it.
    snow_year: Option<SnowYear>,
    equinox: f64,
}

/// The land's snow's year at the cycle's places (F273): each place's precipitation through
/// the cycle's steps, kg/m²/s, and its air at the datum, K, none where its land holds no
/// snow's year; the body's melt for each kelvin of its air past the melting point, kg/m²/s/K,
/// what melted over the air's degrees past the point where the snow lay short of its most,
/// Earth's snow's where none both lay and melted;
/// the step, s, the point, K, and the steps.
#[derive(Clone, Debug)]
struct SnowYear {
    fall: Vec<Vec<f32>>,
    air: Vec<Vec<f32>>,
    melt: f64,
    dt: f64,
    point: f64,
    steps: usize,
}

impl Cryosphere {
    /// Whether an ice, its land's frost or, `sea`, its sea's ice, lies anywhere in the year.
    pub fn holds(&self, sea: bool) -> bool {
        let windows = if sea { &self.sea } else { &self.land };
        windows.iter().any(|w| w.peak > 0.0)
    }

    /// Whether the land's snow lies anywhere in the year (F271).
    pub fn holds_snow(&self) -> bool {
        self.snow.iter().any(|w| w.peak > 0.0)
    }

    /// The window of the land's snow in a direction on ground `above` m over the datum
    /// (F273): its year run on the ground there, the precipitation of the places about it
    /// falling as snow by the share the air's temperature at its height gives (F268), the
    /// air the places' lapsed to it by `lapse` K/m, and melting by the body's melt for each
    /// kelvin of the air past the point; so that the snow lies the longer and the deeper on
    /// higher ground, along the terrain's heights rather than the places' cells.
    pub fn snow_over(&self, d: [f64; 3], above: f64, lapse: f64) -> Window {
        let Some(year) = &self.snow_year else {
            return Window::default();
        };
        let (lon, lat) = lon_lat_from_direction(d);
        let corners = reader(self.layout.0, self.layout.1, lat, lon);
        let live: Vec<(usize, f64)> = corners
            .iter()
            .copied()
            .filter(|&(p, w)| w > 0.0 && year.fall[p].len() == year.steps)
            .collect();
        let total: f64 = live.iter().map(|c| c.1).sum();
        if total <= 0.0 {
            return Window::default();
        }
        let steps = year.steps;
        let (mut fall, mut air) = (vec![0.0f64; steps], vec![0.0f64; steps]);
        for &(p, w) in &live {
            for k in 0..steps {
                fall[k] += w / total * year.fall[p][k] as f64;
                air[k] += w / total * year.air[p][k] as f64;
            }
        }
        let cold = lapse * above.max(0.0);
        let share = |t: f64| {
            let x = (t - year.point - RAIN_SNOW_OFFSET) / RAIN_SNOW_WIDTH;
            1.0 / (1.0 + x.clamp(-50.0, 50.0).exp())
        };
        if (0..steps).all(|k| fall[k] * share(air[k] - cold) <= 0.0) {
            return Window::default();
        }
        let (mut m, mut series) = (0.0f64, vec![0.0f32; steps]);
        for round in 0..SNOW_FINE_YEARS {
            for k in 0..steps {
                let t = air[k] - cold;
                m = (m + year.dt * fall[k] * share(t)).min(SNOW_MOST_MASS);
                m = (m - year.dt * year.melt * (t - year.point).max(0.0)).max(0.0);
                if round + 1 == SNOW_FINE_YEARS {
                    series[k] = m as f32;
                }
            }
        }
        let cover = SNOW_COVER_MASS * std::f64::consts::LN_2;
        Window::of(&series, cover, |t| {
            (t / steps as f64 - self.equinox).rem_euclid(1.0)
        })
    }

    /// The window of an ice in a direction (F246), its land's frost's or, `sea`, its sea's
    /// ice's: its peak read between the cycle's places, as its mass is, and its phases those
    /// of the places it lies at, by their weights, on the circle of the year, its window
    /// lengthening toward the places it lies at all year by their share (F255). Between places a sea's ice lies at and places
    /// it never lies at, it lies only within its margin's edge (F251), its peak thinning to it.
    pub fn at(&self, d: [f64; 3], sea: bool) -> Window {
        if sea {
            self.window(d, (&self.sea, &self.sea_margin, SEA_ICE_COVER_MASS))
        } else {
            self.window(d, (&self.land, &self.land_margin, FROST_COVER_MASS))
        }
    }

    /// The window of an ice of `windows` at the places, `margin` its edge between them and
    /// `lies` the mass it lies at, in a direction, as `at` reads it.
    fn window(&self, d: [f64; 3], (windows, margin, lies): (&[Window], &Margin, f64)) -> Window {
        let (lon, lat) = lon_lat_from_direction(d);
        let corners = margin.ground(&reader(self.layout.0, self.layout.1, lat, lon));
        let depth = margin.depth(&corners, self.guide.at(d));
        if depth.is_some_and(|depth| depth <= 0.0) {
            return Window::default();
        }
        // Within the margin, the windows of the places the ice lies at that it reads.
        let lying: Vec<(f64, Window)> = match depth {
            Some(_) => margin
                .kind(&corners, true)
                .into_iter()
                .map(|(p, w)| (w, windows[p]))
                .collect(),
            None => corners
                .iter()
                .filter(|&&(p, w)| w > 0.0 && windows[p].peak > 0.0)
                .map(|&(p, w)| (w, windows[p]))
                .collect(),
        };
        let peak: f64 = match depth {
            Some(depth) => {
                let (sum, weight) = lying.iter().fold((0.0, 0.0), |(s, t), (w, window)| {
                    (s + w * window.peak, t + w)
                });
                depth * sum / weight.max(1e-300)
            }
            None => corners.iter().map(|&(p, w)| w * windows[p].peak).sum(),
        };
        if lying.is_empty() {
            return Window::default();
        }
        // Within the margin, ice thinner than the mass it lies at never lies (F252).
        if depth.is_some() && peak <= lies * std::f64::consts::LN_2 {
            return Window::default();
        }
        let circular = |among: &[(f64, Window)], pick: fn(&Window) -> f64| {
            let (s, c) = among.iter().fold((0.0, 0.0), |(s, c), (w, window)| {
                let angle = TAU * pick(window);
                (s + w * angle.sin(), c + w * angle.cos())
            });
            (s.atan2(c) / TAU).rem_euclid(1.0)
        };
        let at = circular(&lying, |w| w.at);
        let seasonal: Vec<(f64, Window)> = lying
            .iter()
            .copied()
            .filter(|(_, w)| !w.perennial())
            .collect();
        let weight = |among: &[(f64, Window)]| among.iter().map(|(w, _)| w).sum::<f64>();
        // The share of the places it reads that it lies at all year.
        let all_year = 1.0 - weight(&seasonal) / weight(&lying);
        let (freeze, melt) = if seasonal.is_empty() {
            (0.0, 1.0)
        } else {
            let (freeze, melt) = (
                circular(&seasonal, |w| w.freeze),
                circular(&seasonal, |w| w.melt),
            );
            // Toward the places it lies at all year, it freezes the earlier and melts the
            // later, by their share, until it lies all year at them (F255).
            let span = (melt - freeze).rem_euclid(1.0);
            let longer = span + all_year * (1.0 - span);
            let grow = 0.5 * (longer - span);
            let (freeze, melt) = (
                (freeze - grow).rem_euclid(1.0),
                (melt + grow).rem_euclid(1.0),
            );
            // Toward its edge the ice lies the more briefly, its window closing on its peak
            // with its depth (F252), as the frost at a cap's margin comes late and goes early.
            let (before, after) = ((at - freeze).rem_euclid(1.0), (melt - at).rem_euclid(1.0));
            match depth {
                _ if longer >= 1.0 - 1e-9 => (0.0, 1.0),
                Some(depth) if before + after <= (melt - freeze).rem_euclid(1.0) + 1e-9 => (
                    (at - depth * before).rem_euclid(1.0),
                    (at + depth * after).rem_euclid(1.0),
                ),
                _ => (freeze, melt),
            }
        };
        Window {
            freeze,
            melt,
            peak,
            at,
        }
    }
}

/// A place's year: its land's mean, warmest and coldest temperature, the share of the year
/// frost lies on it and its frost at the snapshot's season; and its sea's (F244): the
/// surface's mean, coldest and warmest, its ice at the snapshot's season and the share of
/// the year ice lies on it, the land's temperatures where the place has no sea.
#[derive(Clone, Copy, Debug)]
struct Year {
    mean: f64,
    max: f64,
    min: f64,
    frost_share: f64,
    frost_now: f64,
    sea: [f64; 5],
    /// The surface's temperatures through the cycle, as quantiles, and the albedo its ground
    /// stood at (F256).
    spread: [f32; CYCLE_SPREAD],
    albedo: f32,
}

/// The cycle at a tilt (F245): the places' years; where the atmosphere lays frost, the
/// surface pressure's range, Pa, the caps' growth in a year over the body, kg/m², and the
/// frost's gas's partial pressure at each step; and the air over each place at each step,
/// K.
struct TiltYear {
    years: Vec<Year>,
    range: Option<(f64, f64)>,
    growth: Option<f64>,
    partials: Vec<f64>,
    air: Vec<f64>,
    /// Each place's frost, its sea's ice and its land's snow at the end of each step, kg/m²,
    /// place after place (F246, F271), and its land's surface and its sea's, K (F268); and the
    /// snow its land's melt took through each step, kg/m² (F277).
    masses: [Vec<f32>; 4],
    surfaces: [Vec<f32>; 2],
    /// What the frost's own gas's caps gained over each place's area in a year, kg/m²
    /// (F258).
    caps: Vec<f64>,
}

/// The lockstep through its last cycle (F239, F245): the places' runs; where the atmosphere
/// lays frost, the surface pressure's range, Pa, the frost the caps grow by in the cycle over
/// the body, kg/m², and the frost's gas's partial pressure after each step; and the air over
/// each place after each step, K.
struct Lockstep {
    runs: Vec<Run>,
    range: Option<(f64, f64)>,
    growth: Option<f64>,
    partials: Vec<f64>,
    air: Vec<f64>,
    /// Each place's frost, its sea's ice and its land's snow at the end of each step, kg/m²,
    /// place after place (F246, F271), and its land's surface and its sea's, K (F268); and the
    /// snow its land's melt took through each step, kg/m² (F277).
    masses: [Vec<f32>; 4],
    surfaces: [Vec<f32>; 2],
    /// What the frost's own gas's caps gained over each place's area in the cycle, kg/m²
    /// (F258).
    caps: Vec<f64>,
    /// The gas in the air and its frost at the end, kg/m² over the body, apart from what the
    /// caps gave or took (F258).
    active: f64,
}

/// A lockstep's end on the latitudes' rows, to start the places from (F258): the gas in the
/// air and its frost, kg/m² over the body, and each row's ground.
struct Warm {
    active: f64,
    rows: Vec<WarmRow>,
}

/// A row's ground at a lockstep's end (F258): its land's column and frost, its caps' column
/// and its sea.
struct WarmRow {
    t: Vec<f64>,
    mass: f64,
    cap: Option<Vec<f64>>,
    sea: Option<Sea>,
}

impl Lockstep {
    /// Its state at its end on the latitudes' rows, to start the places from (F258).
    fn warm(self) -> Warm {
        Warm {
            active: self.active,
            rows: self
                .runs
                .into_iter()
                .map(|r| WarmRow {
                    t: r.t,
                    mass: r.mass,
                    cap: r.cap.map(|c| c.t),
                    sea: r.sea,
                })
                .collect(),
        }
    }
}

impl Year {
    /// A year without frost or sea, its surface's temperatures `spread` at `albedo`.
    fn bare(mean: f64, max: f64, min: f64, (spread, albedo): ([f32; CYCLE_SPREAD], f64)) -> Self {
        Self {
            mean,
            max,
            min,
            frost_share: 0.0,
            frost_now: 0.0,
            sea: [mean, min, max, 0.0, 0.0],
            spread,
            albedo: albedo as f32,
        }
    }
}

/// The quantiles of a place's surface temperatures through its cycle (F256): CYCLE_SPREAD of
/// them, each the middle of an equal share of the sampled times.
fn spread_of(samples: &mut [f64]) -> [f32; CYCLE_SPREAD] {
    let mut out = [0.0f32; CYCLE_SPREAD];
    if samples.is_empty() {
        return out;
    }
    samples.sort_by(|a, b| a.total_cmp(b));
    let n = samples.len();
    for (i, q) in out.iter_mut().enumerate() {
        let k = ((i as f64 + 0.5) / CYCLE_SPREAD as f64 * n as f64) as usize;
        *q = samples[k.min(n - 1)] as f32;
    }
    out
}

/// The mean of a rate through a place's cycle (F256), its temperatures `spread` brought from
/// the albedo they stood at to another, `albedos` (from, to), as the radiation they balance
/// scales: what a surface brighter or darker than the climate's ground takes from the same
/// sun.
pub fn mean_rate(spread: &[f64], (from, to): (f64, f64), rate: impl Fn(f64) -> f64) -> f64 {
    let scale = ((1.0 - to).max(0.0) / (1.0 - from).max(1e-6)).powf(0.25);
    spread.iter().map(|&t| rate(t * scale)).sum::<f64>() / spread.len().max(1) as f64
}

/// F250: how strongly ice and snow brighten the ground under a star of effective temperature
/// `teff` K, against under the Sun: by its light's share in the visible and its weighted share
/// in the near band, which they reflect, against the Sun's.
pub fn ice_brightening(teff: f64) -> f64 {
    let reflected = |t: f64| {
        let (visible, near) = ICE_BANDS_UM;
        let below = blackbody_share_below(visible, t);
        below + ICE_NEAR_BAND_WEIGHT * (blackbody_share_below(near, t) - below)
    };
    if teff > 0.0 {
        reflected(teff) / reflected(T_SUN)
    } else {
        1.0
    }
}

impl EnergyBalance {
    pub fn from_graph(graph: &ParamGraph) -> Self {
        let empty = std::collections::BTreeMap::new();
        let composition = graph
            .get("atmosphere.composition")
            .as_map()
            .unwrap_or(&empty);
        let columns = species_columns(composition, graph.f64("atmosphere.column_mass"));
        let inventory = graph
            .get("volatiles.inventory")
            .as_map()
            .and_then(|m| m.get("H2O").copied())
            .unwrap_or(0.0);
        // As the zero-dimensional climate: a water-bearing atmosphere holds vapour.
        let water = inventory
            + if composition.get("H2O").copied().unwrap_or(0.0) > 0.0 {
                1.0
            } else {
                0.0
            };
        let airless = graph.f64("atmosphere.surface_pressure") <= 0.0;
        let (ventilation, heat) = ventilation(graph);
        let liquid = crate::stages::caps::liquid_species(graph.str("volatiles.working_liquid"))
            .first()
            .and_then(|s| volatile(s));
        let evaporation = match liquid {
            Some(w) if ventilation > 0.0 => {
                ventilation * GAS_CONSTANT
                    / w.r_specific
                    / mean_molecular_mass(composition)
                    / graph.f64("atmosphere.surface_pressure")
            }
            _ => 0.0,
        };
        Self {
            sky: Sky::from_graph(graph),
            albedo: graph.f64("atmosphere.bond_albedo"),
            low_sun: LOW_SUN_ALBEDO
                * (graph.f64("atmosphere.column_mass") / LOW_SUN_COLUMN).clamp(0.0, 1.0),
            ice: graph.bool("climate.ice_albedo_feedback")
                && graph.str("volatiles.working_liquid") == "water",
            brightening: ice_brightening(graph.f64("star.teff")),
            transport: if airless {
                0.0
            } else {
                graph.f64("climate.transport_coefficient")
            },
            tau_dry: if airless { 0.0 } else { tau_dry(&columns) },
            water: if airless { 0.0 } else { water },
            start: graph.f64("atmosphere.surface_temperature_mean"),
            heat_flux: graph.f64("interior.surface_heat_flux").max(0.0),
            year: graph.f64("orbit.year"),
            frost: Frost::from_graph(graph),
            pressure: graph.f64("atmosphere.surface_pressure").max(0.0),
            snapshot: snapshot_anomaly(
                graph.f64("climate.season"),
                graph.f64("orbit.perihelion_season").to_radians(),
                graph.f64("orbit.star_eccentricity"),
            ),
            exchange: ventilation * heat,
            evaporation,
            liquid: if evaporation > 0.0 { liquid } else { None },
            rows: false,
            circulation: crate::wind::Air::from_graph(graph),
        }
    }

    /// The year's balance at another state of the orbit's cycles (F259): at an obliquity, rad,
    /// and a season at perihelion, rad from the northern spring equinox, its places on its
    /// latitudes' rows.
    pub fn at_state(&self, obliquity: f64, perihelion: f64) -> Self {
        let mut state = *self;
        state.sky.obliquity = obliquity;
        state.sky.perihelion_season = perihelion;
        state.rows = true;
        state
    }

    /// The heat the air over a place exchanges with its ground per kelvin (F245), W/m²/K, at
    /// the ground's temperature: by its turbulence (F237), and by the radiation its
    /// greenhouse takes from the ground and gives back, the share 3 tau / 4 / (1 + 3 tau / 4)
    /// of the ground's own radiation's slope the greenhouse holds; none without an
    /// atmosphere.
    fn coupling(&self, t: f64) -> f64 {
        if self.exchange <= 0.0 {
            return 0.0;
        }
        let g = 0.75 * (self.tau_dry + water_optical_depth(t, self.water));
        self.exchange + 4.0 * SIGMA_SB * t.powi(3) * g / (1.0 + g)
    }

    /// The latent heat a wet surface at `t` gives the air per kelvin it stands above it, the
    /// air held saturated, at the stated wind (F245), W/m²/K; none without a working liquid.
    fn latent_slope(&self, t: f64) -> f64 {
        match self.liquid {
            Some(w) if self.evaporation > 0.0 => {
                let l = if t >= w.t_triple { w.l_vap } else { w.l_sub };
                l * self.evaporation * vapour_pressure(w, t) * l / (w.r_specific * t * t)
            }
            _ => 0.0,
        }
    }

    /// The heat the air over a place exchanges with its ground per kelvin at the place's
    /// contact (F269): its turbulence at the place's own wind, the radiation as before.
    fn coupling_at(&self, t: f64, contact: &Contact) -> f64 {
        let turbulence = self.exchange * (contact.breeze - 1.0);
        (self.coupling(t) + turbulence).max(0.0)
    }

    /// The heat the air exchanges with a place's land per kelvin at its contact (F269): its
    /// ground's, and the latent heat of what its wet ground gives the air.
    fn land_coupling(&self, t: f64, contact: &Contact) -> f64 {
        self.coupling_at(t, contact) + self.latent_slope(t) * contact.breeze * contact.wet
    }

    /// The heat the air exchanges with a place's sea per kelvin at its contact (F269).
    fn sea_coupling_at(&self, t: f64, contact: &Contact) -> f64 {
        self.coupling_at(t, contact) + self.latent_slope(t) * contact.breeze
    }

    /// The land's radiation and the heat it gives air at `air` K, its latent heat's among
    /// it, at the place's contact (F269), with their slope, at the ground's temperature.
    fn land_shed(&self, air: f64, contact: Contact) -> impl Fn(f64) -> (f64, f64) + Copy + '_ {
        move |t: f64| {
            let (r, slope) = self.olr_and_slope(t);
            let h = self.land_coupling(t, &contact);
            (r + h * (t - air), slope + h)
        }
    }

    /// The frost's radiation and the heat it gives air at `air` K at the place's contact
    /// (F247, F269), with their slope.
    fn frost_shed_at(
        &self,
        air: f64,
        frost: &Frost,
        contact: Contact,
    ) -> impl Fn(f64) -> (f64, f64) + Copy + '_ {
        let f = *frost;
        move |t: f64| {
            let (r, slope) = self.olr_and_slope(t);
            let (e, de) = f.emissivity(t);
            let h = self.coupling_at(t, &contact);
            (e * r + h * (t - air), e * slope + de * r + h)
        }
    }

    /// The ground's radiation and the heat it gives air at `air` K (F245), with their slope,
    /// at the ground's temperature.
    fn shed(&self, air: f64) -> impl Fn(f64) -> (f64, f64) + Copy + '_ {
        move |t: f64| {
            let (r, slope) = self.olr_and_slope(t);
            let h = self.coupling(t);
            (r + h * (t - air), slope + h)
        }
    }

    /// The frost's radiation at its emissivity and the heat it gives air at `air` K (F247),
    /// with their slope, at the frost's temperature.
    fn frost_shed(&self, air: f64, frost: &Frost) -> impl Fn(f64) -> (f64, f64) + Copy + '_ {
        let f = *frost;
        move |t: f64| {
            let (r, slope) = self.olr_and_slope(t);
            let (e, de) = f.emissivity(t);
            let h = self.coupling(t);
            (e * r + h * (t - air), e * slope + de * r + h)
        }
    }

    /// A place's land's albedo at its temperature (F245): its ground's, brightening under ice
    /// and snow of albedo `snow` where water is the working liquid by the ramp the year's
    /// balance reads (F235), on the land's own temperature through the year rather than the
    /// year's mean (North and Coakley 1979).
    fn land_albedo(&self, (ground, snow): (f64, f64), t: f64) -> f64 {
        let (share, _) = self.ice_share(t);
        ground + (snow - ground).max(0.0) * share
    }

    /// A place's land's albedo under its snow (F271) at its temperature: its ground's,
    /// brightening under snow of albedo `snow` by the snow's cover, and the less as the snow
    /// nears its melting, by the ramp the energy balance's ice and snow brighten the ground by
    /// (F235, F245): melting snow is the darker and the patchier.
    fn snow_albedo(&self, (ground, snow): (f64, f64), mass: f64, t: f64) -> f64 {
        let (share, _) = self.ice_share(t);
        ground + (snow - ground).max(0.0) * snow_cover(mass) * share
    }

    /// The albedo of ice or snow whose albedo under the Sun's light is `ice`, over ground of
    /// albedo `ground`, under the star's light (F250).
    fn ice_over(&self, ice: f64, ground: f64) -> f64 {
        ground + self.brightening * (ice - ground)
    }

    /// An ice's albedo under the star's light (F250): water's, its colour's under the Sun's,
    /// brightening the planet's ground as the energy balance's ice and snow do; the other
    /// ices', their colours'.
    pub fn ice_albedo(&self, ice: &Ice) -> f64 {
        if ice.species == "H2O" {
            self.ice_over(ice.albedo(), self.albedo)
        } else {
            ice.albedo()
        }
    }

    /// The outgoing radiation and its slope at a surface temperature.
    fn olr_and_slope(&self, t: f64) -> (f64, f64) {
        let h = 0.01;
        (self.olr(t), (self.olr(t + h) - self.olr(t - h)) / (2.0 * h))
    }

    /// The outgoing radiation at a surface temperature, W/m².
    pub fn olr(&self, t: f64) -> f64 {
        let tau = self.tau_dry + water_optical_depth(t, self.water);
        SIGMA_SB * t.powi(4) / (1.0 + 0.75 * tau)
    }

    /// The share of a cell ice and snow cover at a temperature, from none at freezing to
    /// all at the ice line, and its slope.
    fn ice_share(&self, t: f64) -> (f64, f64) {
        if !self.ice {
            return (0.0, 0.0);
        }
        let (lo, hi) = ICE_LINE_RANGE;
        if t >= hi {
            (0.0, 0.0)
        } else if t <= lo {
            (1.0, 0.0)
        } else {
            ((hi - t) / (hi - lo), -1.0 / (hi - lo))
        }
    }

    /// The year's mean climate on a grid.
    pub fn solve(&self, grid: CubeSphere) -> Climate {
        self.solve_over(grid, None)
    }

    /// The ground's albedo on a grid without ice and snow (F235): the body's, with the
    /// low-sun term about the sunlight-weighted mean height, so that without ice the planet
    /// reflects its Bond albedo; and under `caps` of ice (F241), by the share of each cell
    /// they cover.
    fn albedos(
        &self,
        grid: CubeSphere,
        sun: &Sunlight,
        caps: Option<&(Field<f32>, Field<f32>)>,
    ) -> (Vec<f64>, Vec<f64>) {
        let n = grid.cell_count();
        let area: Vec<f64> = (0..n).map(|k| grid.cell_area(grid.cell(k))).collect();
        let flux: Vec<f64> = sun.flux.as_slice().iter().map(|&q| q as f64).collect();
        let height: Vec<f64> = sun.height.as_slice().iter().map(|&m| m as f64).collect();
        let weight: f64 = (0..n).map(|k| flux[k] * area[k]).sum();
        let mean_height = if weight > 0.0 {
            (0..n).map(|k| flux[k] * area[k] * height[k]).sum::<f64>() / weight
        } else {
            0.0
        };
        let bare: Vec<f64> = (0..n)
            .map(|k| (self.albedo + self.low_sun * (mean_height - height[k])).clamp(0.0, 1.0))
            .collect();
        let base: Vec<f64> = (0..n)
            .map(|k| {
                let a = bare[k];
                match caps {
                    Some((share, albedo)) => {
                        let c = grid.cell(k);
                        let s = share.get(c) as f64;
                        a + s * (albedo.get(c) as f64 - a)
                    }
                    None => a,
                }
            })
            .collect();
        (bare, base)
    }

    /// The ground's albedo on a grid under `caps` of ice (F241), without ice and snow: the
    /// albedo the year's places run their land at (F258).
    pub fn ground_albedo(
        &self,
        grid: CubeSphere,
        caps: Option<&(Field<f32>, Field<f32>)>,
    ) -> Field<f32> {
        let sun = self.sky.annual_sunlight(grid);
        let (_, base) = self.albedos(grid, &sun, caps);
        Field::from_fn(grid, |c, _| base[grid.index(c)] as f32)
    }

    /// The year's mean climate on a grid over caps of ice (F241): `caps` the share of each
    /// cell they cover and their albedo there, which each cell reflects at by its share.
    pub fn solve_over(&self, grid: CubeSphere, caps: Option<&(Field<f32>, Field<f32>)>) -> Climate {
        let sun = self.sky.annual_sunlight(grid);
        let n = grid.cell_count();
        let area: Vec<f64> = (0..n).map(|k| grid.cell_area(grid.cell(k))).collect();
        let links = links(grid);
        let flux: Vec<f64> = sun.flux.as_slice().iter().map(|&q| q as f64).collect();
        let (bare, base) = self.albedos(grid, &sun, caps);
        // Ice and snow under the star's light over the ground without its caps (F250).
        let snow: Vec<f64> = bare.iter().map(|&a| self.ice_over(ICE_ALBEDO, a)).collect();
        let albedo = |k: usize, t: f64| -> (f64, f64) {
            let (share, slope) = self.ice_share(t);
            let brighter = (snow[k] - base[k]).max(0.0);
            (base[k] + brighter * share, brighter * slope)
        };
        let d = self.transport;
        let mut t = vec![self.start.max(1.0); n];
        // The warm branch (F249): where ice and snow brighten the cold ground, the balance
        // starts from its climate without them, so that a body whose mean stands below
        // freezing keeps the warmth its sunlight's pattern gives it where that is a balance,
        // as a body locked to its star keeps its sea under its star.
        if self.ice {
            let bare = EnergyBalance {
                ice: false,
                ..*self
            }
            .solve_over(grid, caps);
            t = bare
                .temperature
                .as_slice()
                .iter()
                .map(|&v| v as f64)
                .collect();
        }
        // An airless cell's start: its own balance with the year's sunlight and its heat.
        if d <= 0.0 {
            for k in 0..n {
                let (a, _) = albedo(k, t[k]);
                t[k] = ((flux[k] * (1.0 - a)).max(0.0) / SIGMA_SB + self.heat_flux / SIGMA_SB)
                    .powf(0.25)
                    .max(COSMIC_BACKGROUND_K);
            }
        }
        let mut iterations = 0;
        for it in 0..CLIMATE_MAX_ITERATIONS {
            iterations = it + 1;
            let (residual, slope): (Vec<f64>, Vec<f64>) = (0..n)
                .into_par_iter()
                .map(|k| {
                    let (a, da) = albedo(k, t[k]);
                    let h = 0.01;
                    let olr = self.olr(t[k]);
                    let dolr = (self.olr(t[k] + h) - self.olr(t[k] - h)) / (2.0 * h);
                    let mut r = area[k] * (flux[k] * (1.0 - a) + self.heat_flux - olr);
                    for &(j, w) in &links[k] {
                        r += d * w * (t[j] - t[k]);
                    }
                    // d(-F)/dT, the ice's feedback floored below the radiation's slope.
                    let local = (dolr + flux[k] * da).max(CLIMATE_SLOPE_FLOOR * dolr);
                    (r, area[k] * local)
                })
                .unzip();
            let step = conjugate_gradients(&slope, &links, d, &residual);
            let mut largest: f64 = 0.0;
            for k in 0..n {
                let s = step[k].clamp(-CLIMATE_MAX_STEP_K, CLIMATE_MAX_STEP_K);
                t[k] = (t[k] + s).max(COSMIC_BACKGROUND_K);
                largest = largest.max(s.abs());
            }
            if largest < CLIMATE_TOLERANCE_K {
                break;
            }
        }
        let reflected: Vec<f64> = (0..n).map(|k| albedo(k, t[k]).0).collect();
        let brought: Vec<f64> = (0..n)
            .map(|k| {
                links[k]
                    .iter()
                    .map(|&(j, w)| d * w * (t[j] - t[k]))
                    .sum::<f64>()
                    / area[k]
            })
            .collect();
        Climate {
            temperature: Field::from_fn(grid, |c, _| t[grid.index(c)] as f32),
            absorbed: Field::from_fn(grid, |c, _| {
                let k = grid.index(c);
                (flux[k] * (1.0 - reflected[k])) as f32
            }),
            albedo: Field::from_fn(grid, |c, _| reflected[grid.index(c)] as f32),
            ground_albedo: Field::from_fn(grid, |c, _| base[grid.index(c)] as f32),
            snow_albedo: Field::from_fn(grid, |c, _| snow[grid.index(c)] as f32),
            transport: Field::from_fn(grid, |c, _| brought[grid.index(c)] as f32),
            iterations,
        }
    }

    /// The year's surface temperatures through the days and seasons on the climate's grid
    /// (F236): at each place a column of ground under the star's path, its sunlight absorbed
    /// at its albedo. A body that spins many times a year is run through the year on the
    /// day's mean sunlight until its seasons repeat, then through a day at each season; a
    /// body in a spin-orbit resonance through the orbits it repeats in, until they repeat. An
    /// airless body's places stand alone, each at the albedo its year's balance found; under
    /// an atmosphere they run in step under the air that carries the transport between them
    /// (`lockstep`), each over its land and, `sea` its share of each cell of the climate's
    /// grid, its sea (F245), and its share of the frost's own gas's caps, the `reservoir`
    /// its air draws on (F258).
    pub fn cycle(
        &self,
        climate: &Climate,
        ground: &Ground,
        sea: Option<(&Ocean, &Field<f32>)>,
        reservoir: Option<&Reservoir>,
    ) -> Cycle {
        self.cycle_with(climate, ground, sea, reservoir, None)
    }

    /// The year as `cycle` runs it, its air's exchange with the ground at each place's
    /// contact from a year's wind and moisture, `coupling` (F269).
    pub fn cycle_with(
        &self,
        climate: &Climate,
        ground: &Ground,
        sea: Option<(&Ocean, &Field<f32>)>,
        reservoir: Option<&Reservoir>,
        coupling: Option<&Coupling>,
    ) -> Cycle {
        let grid = climate.temperature.grid();
        let all = self.sky.obliquities();
        let tilts: Vec<f64> = if all.len() > CYCLE_SWING_SAMPLES {
            (0..CYCLE_SWING_SAMPLES)
                .map(|j| all[j * all.len() / CYCLE_SWING_SAMPLES])
                .collect()
        } else {
            all
        };
        let (n_lat, n_lon) = self.cycle_grid();
        let paths: Vec<Vec<PathStep>> = match self.sky.resonance {
            Some(_) => tilts.iter().map(|&eps| self.resonant_steps(eps)).collect(),
            None => Vec::new(),
        };
        let ocean = sea.map(|s| s.0);
        let layout = (self.exchange > 0.0)
            .then(|| self.layout(climate, sea.map(|s| s.1), reservoir, coupling));
        // A fast spinner whose caps hold its air spins its seasons up on its latitudes' rows
        // first, its places started from them (F258).
        let rows = layout
            .as_ref()
            .filter(|l| {
                l.n_lon > 1 && self.sky.resonance.is_none() && l.places.iter().any(|p| p.cap > 0.0)
            })
            .map(|_| {
                let rows = (CYCLE_LATITUDES, 1);
                self.layout_on(climate, sea.map(|s| s.1), reservoir, rows, coupling)
            });
        let each: Vec<TiltYear> = if let Some(layout) = &layout {
            tilts
                .iter()
                .enumerate()
                .map(|(j, &eps)| match self.sky.resonance {
                    None => self.seasonal_years(ground, layout, eps, ocean, rows.as_ref()),
                    Some((_, orbits)) => {
                        self.seasonal_orbits(ground, layout, (orbits, &paths[j]), ocean)
                    }
                })
                .collect()
        } else {
            let lat_of = |i: usize| -PI / 2.0 + PI * i as f64 / (n_lat - 1) as f64;
            let lon_of = |j: usize| -PI + TAU * j as f64 / n_lon as f64;
            let at = |lat: f64, lon: f64| -> (f64, f64, f64) {
                let (c, _, _) = grid.locate(direction_from_lon_lat(lon, lat));
                (
                    climate.temperature.get(c) as f64,
                    climate.albedo.get(c) as f64,
                    climate.transport.get(c) as f64,
                )
            };
            tilts
                .iter()
                .enumerate()
                .map(|(j, &eps)| TiltYear {
                    years: (0..n_lat * n_lon)
                        .into_par_iter()
                        .map(|k| {
                            let place = (lat_of(k / n_lon), lon_of(k % n_lon));
                            let start = at(place.0, place.1);
                            match self.sky.resonance {
                                None => {
                                    let (steady, albedo, brought) = start;
                                    let lat = place.0;
                                    self.zonal_cycle(ground, lat, eps, steady, albedo, brought)
                                }
                                Some((_, orbits)) => {
                                    self.resonant_cycle(ground, orbits, &paths[j], place, start)
                                }
                            }
                        })
                        .collect(),
                    range: None,
                    growth: None,
                    partials: Vec::new(),
                    air: Vec::new(),
                    masses: [Vec::new(), Vec::new(), Vec::new(), Vec::new()],
                    surfaces: [Vec::new(), Vec::new()],
                    caps: Vec::new(),
                })
                .collect()
        };
        // The tilts' years combined: each place's means their mean and its extremes the
        // widest, the pressure's range the widest and the caps' growth the mean.
        let n = each.len().max(1) as f64;
        let rows: Vec<Year> = (0..n_lat * n_lon)
            .map(|k| {
                each.iter().map(|e| e.years[k]).fold(
                    Year {
                        mean: 0.0,
                        max: f64::MIN,
                        min: f64::MAX,
                        frost_share: 0.0,
                        frost_now: 0.0,
                        sea: [0.0, f64::MAX, f64::MIN, 0.0, 0.0],
                        spread: [0.0; CYCLE_SPREAD],
                        albedo: 0.0,
                    },
                    |acc, e| Year {
                        mean: acc.mean + e.mean / n,
                        max: acc.max.max(e.max),
                        min: acc.min.min(e.min),
                        frost_share: acc.frost_share + e.frost_share / n,
                        frost_now: acc.frost_now + e.frost_now / n,
                        sea: [
                            acc.sea[0] + e.sea[0] / n,
                            acc.sea[1].min(e.sea[1]),
                            acc.sea[2].max(e.sea[2]),
                            acc.sea[3] + e.sea[3] / n,
                            acc.sea[4] + e.sea[4] / n,
                        ],
                        spread: std::array::from_fn(|i| acc.spread[i] + e.spread[i] / n as f32),
                        albedo: acc.albedo + e.albedo / n as f32,
                    },
                )
            })
            .collect();
        let pressure_range = each
            .iter()
            .filter_map(|e| e.range)
            .reduce(|a, b| (a.0.min(b.0), a.1.max(b.1)));
        let growths: Vec<f64> = each.iter().filter_map(|e| e.growth).collect();
        let cap_growth =
            (!growths.is_empty()).then(|| growths.iter().sum::<f64>() / growths.len() as f64);
        // What the caps that hold the air gained through the year at each place, the tilts'
        // mean, and over the body (F258).
        let reservoir = layout
            .as_ref()
            .filter(|l| l.places.iter().any(|p| p.cap > 0.0))
            .map(|l| {
                let gains: Vec<f64> = (0..n_lat * n_lon)
                    .map(|k| each.iter().map(|e| e.caps[k]).sum::<f64>() / n)
                    .collect();
                let whole: f64 = l.places.iter().map(|p| p.weight).sum();
                let total = gains
                    .iter()
                    .zip(&l.places)
                    .map(|(g, p)| g * p.weight)
                    .sum::<f64>()
                    / whole;
                let field = Field::from_fn(grid, |_, d| {
                    let (lon, lat) = lon_lat_from_direction(d);
                    let value = match n_lon {
                        1 => read_row(&gains, lat),
                        _ => read_grid(&gains, n_lon, lat, lon),
                    };
                    value as f32
                });
                (field, total)
            });
        let field = |pick: fn(&Year) -> f64| {
            let values: Vec<f64> = rows.iter().map(pick).collect();
            Field::from_fn(grid, |_, d| {
                let (lon, lat) = lon_lat_from_direction(d);
                let value = match n_lon {
                    1 => read_row(&values, lat),
                    _ => read_grid(&values, n_lon, lat, lon),
                };
                value as f32
            })
        };
        // Where the sea's ice ends between the places, guided by the steady balance (F251):
        // whether it holds at each place with a sea.
        let guide = Guide::of(self, climate);
        let wet = |p: usize| layout.as_ref().is_some_and(|l| l.places[p].sea > 0.0);
        // The share of each place's ground the sea holds and the land does (F266).
        let sea_shares: Vec<f64> = (0..n_lat * n_lon)
            .map(|p| layout.as_ref().map_or(0.0, |l| l.places[p].sea))
            .collect();
        let land_shares: Vec<f64> = sea_shares.iter().map(|s| 1.0 - s).collect();
        let held = |holds: &dyn Fn(usize) -> bool| -> Vec<Option<bool>> {
            (0..n_lat * n_lon)
                .map(|p| wet(p).then(|| holds(p)))
                .collect()
        };
        // The ices' windows: each place's masses through its cycle, the tilts' mean step by
        // step, the phase of the year counted from the northern spring equinox.
        let places = n_lat * n_lon;
        let steps = each
            .first()
            .map_or(0, |e| e.masses[0].len() / places.max(1));
        let orbits = self.sky.resonance.map_or(1, |(_, q)| q);
        let e = self.sky.eccentricity;
        let equinox = snapshot_anomaly(0.0, self.sky.perihelion_season, e) / TAU / orbits as f64;
        let phase = |t: f64| (t / steps as f64 - equinox).rem_euclid(1.0);
        // The tilts' mean masses, the land's frost's, the sea's ice's and the land's snow's,
        // and the snow its melt took, place after place, step after step.
        let masses: [Vec<f32>; 4] = std::array::from_fn(|j| {
            (0..if steps > 0 { places * steps } else { 0 })
                .map(|i| each.iter().map(|t| t.masses[j][i]).sum::<f32>() / each.len() as f32)
                .collect()
        });
        let surfaces: [Vec<f32>; 2] = std::array::from_fn(|j| {
            (0..if steps > 0 { places * steps } else { 0 })
                .map(|i| each.iter().map(|t| t.surfaces[j][i]).sum::<f32>() / each.len() as f32)
                .collect()
        });
        let covers = [
            FROST_COVER_MASS * std::f64::consts::LN_2,
            SEA_ICE_COVER_MASS * std::f64::consts::LN_2,
            SNOW_COVER_MASS * std::f64::consts::LN_2,
        ];
        let cryosphere = (steps > 0).then(|| {
            let windows = |j: usize| -> Vec<Window> {
                masses[j]
                    .chunks(steps)
                    .map(|m| Window::of(m, covers[j], phase))
                    .collect()
            };
            let (land, sea, snow) = (windows(0), windows(1), windows(2));
            let sea_margin = Margin::new(
                (n_lat, n_lon),
                held(&|p| sea[p].peak > 0.0),
                &sea_shares,
                &guide,
            );
            let dry = |p: usize| layout.as_ref().is_some_and(|l| l.places[p].sea < 1.0);
            let land_held = (0..places)
                .map(|p| dry(p).then(|| land[p].peak > 0.0))
                .collect();
            // The land's snow's year at each place, to be run again on finer ground (F273):
            // the body's melt for each kelvin of the air past the point, what melted over the
            // air's degrees past it by the share of the ground the snow covered, where it lay
            // short of its most.
            let snow_year = layout.as_ref().filter(|_| self.ice).and_then(|l| {
                let point = self.liquid?.t_triple;
                let dt = self.year / steps as f64;
                let air: Vec<Vec<f32>> = (0..places)
                    .map(|p| {
                        (0..steps)
                            .map(|k| {
                                let sum: f64 = each.iter().map(|t| t.air[k * places + p]).sum();
                                (sum / each.len().max(1) as f64) as f32
                            })
                            .collect()
                    })
                    .collect();
                let (mut melted, mut degrees) = (0.0, 0.0);
                for (p, place) in l.places.iter().enumerate() {
                    let fall = &l.snowfall[p];
                    let lay = &masses[2][p * steps..(p + 1) * steps];
                    if fall.len() != steps || lay.iter().any(|&m| m as f64 >= SNOW_MOST_MASS) {
                        continue;
                    }
                    let w = place.weight * place.land();
                    for k in 0..steps {
                        let before = lay[(k + steps - 1) % steps] as f64;
                        let melt = (before + dt * fall[k] as f64 - lay[k] as f64).max(0.0);
                        let past = (air[p][k] as f64 - point).max(0.0) * dt;
                        let covered = w * snow_cover(before);
                        (melted, degrees) = (melted + w * melt, degrees + covered * past);
                    }
                }
                let melt = if degrees > 0.0 {
                    melted / degrees
                } else {
                    SNOW_MELT_PER_DEGREE
                };
                Some(SnowYear {
                    fall: l.falling.clone(),
                    air,
                    melt,
                    dt,
                    point,
                    steps,
                })
            });
            Cryosphere {
                land,
                sea,
                snow,
                layout: (n_lat, n_lon),
                snapshot: snapshot_phase(self.snapshot, steps, orbits as usize, equinox),
                guide: guide.clone(),
                sea_margin,
                land_margin: Margin::new((n_lat, n_lon), land_held, &land_shares, &guide),
                curves: [0, 1, 2].map(|j| year_curve(&masses[j], steps, covers[j])),
                snow_year,
                equinox,
            }
        });
        // The year step by step (F253): the star's place in the sky at each step's end, the
        // pressure, and each ice's mean over each hemisphere, its land's or its sea's share.
        let series = layout.as_ref().filter(|_| steps > 0).map(|l| {
            let sky = &self.sky;
            let mut series = Series::default();
            let partial = |k: usize| -> f64 {
                each.iter().filter_map(|t| t.partials.get(k)).sum::<f64>() / each.len() as f64
            };
            let others = self
                .frost
                .map_or(self.pressure, |f| (self.pressure - f.partial).max(0.0));
            for k in 0..steps {
                let t = k as f64 + 1.0;
                let mean = (TAU * orbits as f64 * t / steps as f64).rem_euclid(TAU);
                let big_e = crate::insolation::eccentric_anomaly(mean, e);
                let nu = 2.0
                    * ((1.0 + e).sqrt() * (big_e / 2.0).sin())
                        .atan2((1.0 - e).sqrt() * (big_e / 2.0).cos());
                let ls = (nu + sky.perihelion_season).rem_euclid(TAU);
                series.phase.push(phase(t));
                series.ls_deg.push(ls.to_degrees());
                series
                    .declination_deg
                    .push((sky.obliquity.sin() * ls.sin()).asin().to_degrees());
                series.flux.push(sky.flux / (1.0 - e * big_e.cos()).powi(2));
                series.pressure.push(match self.frost {
                    Some(_) => others + partial(k),
                    None => self.pressure,
                });
                let kinds = [&mut series.frost, &mut series.sea_ice, &mut series.snow];
                for (j, kind) in kinds.into_iter().enumerate() {
                    let (mut sums, mut areas) = ([0.0; 2], [0.0; 2]);
                    for (p, place) in l.places.iter().enumerate() {
                        let share = if j == 1 { place.sea } else { place.land() };
                        let north = if place.lat > 1e-9 {
                            1.0
                        } else if place.lat < -1e-9 {
                            0.0
                        } else {
                            0.5
                        };
                        let m = masses[j][p * steps + k] as f64;
                        for (h, part) in [(0, north), (1, 1.0 - north)] {
                            sums[h] += place.weight * part * share * m;
                            areas[h] += place.weight * part;
                        }
                    }
                    for h in 0..2 {
                        kind[h].push(sums[h] / areas[h].max(1e-30));
                    }
                }
            }
            series
        });
        // The clock: the cycle, its solar day, and its seasons' phases in its first orbit.
        let clock = (steps > 0).then(|| {
            let sky = &self.sky;
            let season = |ls: f64| {
                let t = snapshot_anomaly(ls, sky.perihelion_season, e) / TAU / orbits as f64;
                (t - equinox).rem_euclid(1.0)
            };
            // A resonant body's spin is its resonance's; one locked to its star has no day.
            let spins = sky
                .resonance
                .map_or(sky.spins, |(p, q)| p as f64 / q as f64);
            let solar = (spins - 1.0).abs();
            Clock {
                cycle_s: self.year * orbits as f64,
                orbits,
                steps,
                solar_day_s: (solar > 1e-9).then(|| self.year / solar),
                obliquity_deg: sky.obliquity.to_degrees(),
                eccentricity: e,
                perihelion_ls_deg: sky.perihelion_season.to_degrees(),
                seasons: [0.0, 0.25, 0.5, 0.75].map(season),
            }
        });
        // The wind through the year at the places (F267), from the air over each at each step,
        // the tilts' mean; not on the latitudes' rows of the orbit's cycles.
        let sites: Vec<crate::wind::Site> = layout.as_ref().map_or_else(Vec::new, |l| {
            l.places
                .iter()
                .map(|p| crate::wind::Site {
                    lat: p.lat,
                    lon: p.lon,
                    weight: p.weight,
                    sea: p.sea,
                    ice: p.cap,
                })
                .collect()
        });
        let wind = layout.as_ref().filter(|_| !self.rows).and_then(|l| {
            let water = ocean.map_or(0.0, |_| OCEAN_TRANSPORT_SHARE);
            let links = l.air_links(self.transport, water);
            let air: Vec<&[f64]> = each.iter().map(|e| &e.air[..]).collect();
            crate::wind::blow(&self.circulation, &sites, (n_lat, n_lon), &links, &air)
        });
        // The places with a sea, whose fields the seas spread to the others (F273).
        let wet_places: Vec<bool> = (0..n_lat * n_lon).map(wet).collect();
        // The seas' currents the year's wind drives (F272).
        let currents = match (&wind, ocean) {
            (Some(w), Some(_)) => self.currents(w, &sites, &each),
            _ => None,
        };
        let seas = ocean.map(|ocean| Seas {
            mean: field(|y| y.sea[0]),
            min: field(|y| y.sea[1]),
            max: field(|y| y.sea[2]),
            ice_now: field(|y| y.sea[3]),
            ice_share: field(|y| y.sea[4]),
            ocean: *ocean,
            mean_rows: spread_from(
                &rows.iter().map(|y| y.sea[0]).collect::<Vec<_>>(),
                &wet_places,
                (n_lat, n_lon),
            ),
            min_rows: spread_from(
                &rows.iter().map(|y| y.sea[1]).collect::<Vec<_>>(),
                &wet_places,
                (n_lat, n_lon),
            ),
            max_rows: spread_from(
                &rows.iter().map(|y| y.sea[2]).collect::<Vec<_>>(),
                &wet_places,
                (n_lat, n_lon),
            ),
            ice_rows: rows.iter().map(|y| y.sea[3]).collect(),
            cover_rows: rows
                .iter()
                .map(|y| 1.0 - (-y.sea[3].max(0.0) / SEA_ICE_COVER_MASS).exp())
                .collect(),
            margin: Margin::new(
                (n_lat, n_lon),
                held(&|p| rows[p].sea[3] > SEA_ICE_COVER_MASS * std::f64::consts::LN_2),
                &sea_shares,
                &guide,
            ),
            guide: guide.clone(),
            layout: (n_lat, n_lon),
        });
        Cycle {
            mean: field(|y| y.mean),
            max: field(|y| y.max),
            min: field(|y| y.min),
            frost_share: field(|y| y.frost_share),
            frost_now: field(|y| y.frost_now),
            frost: self.frost,
            pressure_range,
            cap_growth,
            reservoir,
            tilts,
            partials: each.iter().map(|e| e.partials.clone()).collect(),
            seas,
            cryosphere,
            air: each.into_iter().map(|e| e.air).collect(),
            layout: (n_lat, n_lon),
            series,
            clock,
            wind,
            currents,
            sites,
            masses,
            surfaces,
            spread: rows.iter().map(|y| y.spread).collect(),
            spread_albedo: rows.iter().map(|y| y.albedo).collect(),
        }
    }

    /// The seas' currents the year's `wind` drives at the cycle's places (F272): its stress on
    /// the sea through the year, rho_a C_D |V| V, the air's density at the surface pressure and
    /// the tilts' mean air over each place at each step, the sea's drag the law of the wall's
    /// over its roughness.
    fn currents(
        &self,
        wind: &crate::wind::Wind,
        sites: &[crate::wind::Site],
        each: &[TiltYear],
    ) -> Option<crate::currents::Currents> {
        let (n_lat, n_lon) = wind.layout;
        let n = n_lat * n_lon;
        let steps = wind.steps.max(1);
        if each.is_empty() || each.iter().any(|e| e.air.len() != steps * n) {
            return None;
        }
        let drag = (VON_KARMAN / (WIND_REFERENCE_HEIGHT / SEA_ROUGHNESS).ln()).powi(2);
        let air = &self.circulation;
        let stress: Vec<[f64; 2]> = (0..n)
            .map(|p| {
                let mut s = [0.0; 2];
                for k in 0..steps {
                    let i = k * n + p;
                    let t = each.iter().map(|e| e.air[i]).sum::<f64>() / each.len() as f64;
                    let rho = air.pressure * air.molar_mass / (GAS_CONSTANT * t.max(1.0));
                    let (u, v) = (wind.east[i] as f64, wind.north[i] as f64);
                    let pull = rho * drag * u.hypot(v) / steps as f64;
                    s[0] += pull * u;
                    s[1] += pull * v;
                }
                s
            })
            .collect();
        let density = self.liquid.map_or(1_000.0, |l| l.liquid_density);
        let spin = (air.omega, air.radius, density);
        let links = place_links((n_lat, n_lon));
        let currents = crate::currents::circulate(sites, (n_lat, n_lon), &stress, spin, &links);
        if let (Some(c), true) = (&currents, std::env::var_os("EXOGEN_DEBUG").is_some()) {
            eprintln!("currents: the strongest gyre {:.1} Sv", c.strongest() / 1e6);
        }
        currents
    }

    /// The star's path through a resonant body's sky as its cycle steps it (F245): an airless
    /// body's at the insolation's samples, which resolve its noons; under an atmosphere, whose
    /// air is solved at every step, at the steps of a fast spinner's year an orbit, or of its
    /// day a solar day where its days are the shorter.
    fn resonant_steps(&self, eps: f64) -> Vec<PathStep> {
        let Some((spins, orbits)) = self.sky.resonance else {
            return Vec::new();
        };
        let samples = if self.exchange > 0.0 {
            let days = (spins as f64 / orbits as f64 - 1.0).abs();
            CYCLE_YEAR_STEPS.max((CYCLE_DAY_STEPS as f64 * days).ceil() as usize)
        } else {
            INSOLATION_ORBIT_SAMPLES
        };
        self.sky.resonant_path_sampled(eps, orbits, samples)
    }

    /// The places the cycle runs at under an atmosphere (F245): a fast spinner's latitudes, a
    /// resonant body's latitudes and longitudes, each with its share of the sphere's area, its
    /// year's balance, its ground's albedo without ice and snow, and the share of it the sea
    /// covers, `sea` on the climate's grid, a latitude's its mean along it; and each place's
    /// links to the places beside it, the transport's conductance to each per unit of the
    /// place's area at unit transport: on the unit sphere, the length of the edge they share
    /// over the distance between them, over the place's area.
    /// The cycle's places (F245, F254): under an atmosphere by latitude and longitude, a fast
    /// spinner's latitudes 2° apart, a resonant body's 4°, so that a fast spinner's land, its
    /// seas and its heights shape its seasons along its latitudes, as they do a resonant
    /// body's, through the air between them; an airless fast spinner's by its latitudes
    /// alone, its places apart, with no air between them.
    fn cycle_grid(&self) -> (usize, usize) {
        match self.sky.resonance {
            None if self.exchange > 0.0 && !self.rows => CYCLE_FAST_GRID,
            None => (CYCLE_LATITUDES, 1),
            Some(_) => CYCLE_RESONANT_GRID,
        }
    }

    fn layout(
        &self,
        climate: &Climate,
        sea: Option<&Field<f32>>,
        reservoir: Option<&Reservoir>,
        coupling: Option<&Coupling>,
    ) -> Layout {
        self.layout_on(climate, sea, reservoir, self.cycle_grid(), coupling)
    }

    /// The cycle's places on a grid of `n_lat` latitudes and `n_lon` longitudes, as `layout`.
    fn layout_on(
        &self,
        climate: &Climate,
        sea: Option<&Field<f32>>,
        reservoir: Option<&Reservoir>,
        (n_lat, n_lon): (usize, usize),
        coupling: Option<&Coupling>,
    ) -> Layout {
        let grid = climate.temperature.grid();
        let (dphi, dlon) = (PI / (n_lat - 1) as f64, TAU / n_lon as f64);
        let lat_of = |i: usize| -PI / 2.0 + dphi * i as f64;
        let edges = |i: usize| {
            let lat = lat_of(i);
            (
                (lat - dphi / 2.0).max(-PI / 2.0),
                (lat + dphi / 2.0).min(PI / 2.0),
            )
        };
        let value_at = |f: &Field<f32>, lat: f64, lon: f64| {
            f.get(f.grid().locate(direction_from_lon_lat(lon, lat)).0) as f64
        };
        // A pair of values summed at a place: along its latitude a fast spinner's, by
        // SEA_SHARE_LONGITUDES samples of `sample` at each longitude, and a place's own on a
        // grid by longitude.
        let along = |sample: &dyn Fn(f64) -> (f64, f64), lon: f64| match n_lon {
            1 => (0..SEA_SHARE_LONGITUDES)
                .map(|m| sample(TAU * m as f64 / SEA_SHARE_LONGITUDES as f64))
                .fold((0.0, 0.0), |a, b| (a.0 + b.0, a.1 + b.1)),
            _ => sample(lon),
        };
        let samples = match n_lon {
            1 => SEA_SHARE_LONGITUDES as f64,
            _ => 1.0,
        };
        let places = (0..n_lat * n_lon)
            .map(|k| {
                let (i, j) = (k / n_lon, k % n_lon);
                let (lat, lon) = (lat_of(i), -PI + dlon * j as f64);
                let (lo, hi) = edges(i);
                let (c, _, _) = grid.locate(direction_from_lon_lat(lon, lat));
                let (sea, _) = along(
                    &|lon| (sea.map_or(0.0, |f| value_at(f, lat, lon)), 0.0),
                    lon,
                );
                let sea = (sea / samples).clamp(0.0, 1.0);
                // The share of the place the caps that hold the air cover and the ice they
                // hold over it (F258), the ground about them at its own albedo.
                let (cap, cap_mass, albedo) = match reservoir {
                    Some(r) => {
                        let (share, held) = along(
                            &|lon| {
                                let s = value_at(&r.share, lat, lon);
                                (s, s * value_at(&r.mass, lat, lon))
                            },
                            lon,
                        );
                        let mass = if share > 0.0 { held / share } else { 0.0 };
                        let cap = (share / samples).clamp(0.0, 1.0 - sea);
                        (cap, mass, r.ground.get(c) as f64)
                    }
                    None => (0.0, 0.0, climate.ground_albedo.get(c) as f64),
                };
                Place {
                    lat,
                    lon,
                    weight: (hi.sin() - lo.sin()) / n_lon as f64,
                    steady: climate.temperature.get(c) as f64,
                    albedo,
                    snow: climate.snow_albedo.get(c) as f64,
                    sea,
                    cap,
                    cap_mass,
                    contact: coupling.map_or_else(Contact::default, |c| c.at((n_lat, n_lon), k)),
                }
            })
            .collect();
        let snowfall = (0..n_lat * n_lon)
            .map(|k| match coupling.filter(|_| self.ice) {
                Some(c) => c.snowfall_at((n_lat, n_lon), CYCLE_YEAR_STEPS, k),
                None => Vec::new(),
            })
            .collect();
        let falling = (0..n_lat * n_lon)
            .map(|k| match coupling.filter(|_| self.ice) {
                Some(c) => c.falling_at((n_lat, n_lon), CYCLE_YEAR_STEPS, k),
                None => Vec::new(),
            })
            .collect();
        let currents = coupling
            .and_then(|c| c.currents.clone())
            .filter(|c| c.layout == (n_lat, n_lon));
        let wind = coupling.map_or_else(Vec::new, |c| c.wind_on((n_lat, n_lon)));
        Layout {
            places,
            links: place_links((n_lat, n_lon)),
            n_lon,
            snowfall,
            falling,
            currents,
            wind,
            modes: place_modes(n_lon),
        }
    }

    /// The years of a fast spinner under its atmosphere (F239, F245): its latitudes run in
    /// step through the year on the day's mean sunlight, then each through a day at each
    /// season from the ground, its frost, the frost point, the air over it and its land's
    /// albedo as the season left them, the air holding through the day.
    fn seasonal_years(
        &self,
        ground: &Ground,
        layout: &Layout,
        eps: f64,
        ocean: Option<&Ocean>,
        rows: Option<&Layout>,
    ) -> TiltYear {
        let solar_days = (self.sky.spins - eps.cos()).abs().max(1e-3);
        let day = self.year / solar_days;
        let column = Column::new(ground, day, self.year);
        let steps = CYCLE_YEAR_STEPS;
        let sun: Vec<(f64, f64)> = (0..steps)
            .map(|m| {
                let anomaly = TAU * (m as f64 + 0.5) / steps as f64;
                self.sky.at_mean_anomaly(anomaly, eps)
            })
            .collect();
        let light = |m: usize, p: &Place| sun[m].0 * diurnal_share(p.lat, sun[m].1);
        let timing = (
            steps,
            self.year / steps as f64,
            snapshot_step(self.snapshot, steps),
        );
        let every = steps / CYCLE_SEASONS;
        // The seasons spun up on the latitudes' rows where caps hold the air (F258): the frost
        // their air lays takes years to come to its own, cheaply on the rows.
        let warm = rows.map(|rows| {
            let spin = (0, None, CYCLE_SPIN_UP_YEARS);
            self.lockstep(ground, &column, rows, timing, light, ocean, spin)
                .warm()
        });
        // On the rows alone, caps that hold the air take the spin-up's years (F258).
        let years = match layout.n_lon == 1 && layout.places.iter().any(|p| p.cap > 0.0) {
            true => CYCLE_SPIN_UP_YEARS,
            false => CYCLE_MAX_YEARS,
        };
        let lockstep = self.lockstep(
            ground,
            &column,
            layout,
            timing,
            light,
            ocean,
            (every, warm.as_ref(), years),
        );
        let dt = day / CYCLE_DAY_STEPS as f64;
        let frost = self.frost.as_ref();
        let years = lockstep
            .runs
            .par_iter()
            .zip(layout.places.par_iter())
            .map(|(r, p)| {
                let mut work = Work::new(column.len());
                let mut trial = Trial::new(column.len());
                let (sp, cp) = p.lat.sin_cos();
                let (mut mean, mut hi, mut lo) = (0.0, f64::MIN, f64::MAX);
                let mut samples = Vec::new();
                for record in &r.records {
                    let f = frost.map(|f| Frost {
                        point: record.point,
                        ..*f
                    });
                    let (flux, decl) = sun[record.k];
                    let (sd, cd) = decl.sin_cos();
                    let shed = self.land_shed(record.air, p.contact);
                    let mut t = record.t.clone();
                    let mut mass = record.mass;
                    let mut sum = 0.0;
                    for d in 0..CYCLE_DAYS {
                        for h in 0..CYCLE_DAY_STEPS {
                            let hour = TAU * (h as f64 + 0.5) / CYCLE_DAY_STEPS as f64 - PI;
                            let lit = flux * (sp * sd + cp * cd * hour.cos()).max(0.0);
                            let bare = (1.0 - record.albedo) * lit;
                            match &f {
                                Some(f) => {
                                    let gains = (bare, (1.0 - f.albedo) * lit);
                                    let rows = (&mut work, &mut trial);
                                    let heat = self.heat_flux;
                                    let sheds =
                                        (shed, self.frost_shed_at(record.air, f, p.contact));
                                    column.frosty_step(
                                        &mut t, rows, dt, gains, heat, sheds, f, &mut mass,
                                    );
                                }
                                None => {
                                    column.step(&mut t, &mut work, dt, bare, self.heat_flux, shed)
                                }
                            }
                            if d + 1 == CYCLE_DAYS {
                                sum += t[0];
                                hi = hi.max(t[0]);
                                lo = lo.min(t[0]);
                                samples.push(t[0]);
                            }
                        }
                    }
                    mean += sum / CYCLE_DAY_STEPS as f64;
                }
                let mean = mean / r.records.len().max(1) as f64;
                let albedo =
                    r.records.iter().map(|x| x.albedo).sum::<f64>() / r.records.len().max(1) as f64;
                Year {
                    mean,
                    max: hi,
                    min: lo,
                    frost_share: r.frosted as f64 / steps as f64,
                    frost_now: r.now,
                    sea: r.sea_year(steps, (mean, hi, lo)),
                    spread: spread_of(&mut samples),
                    albedo: albedo as f32,
                }
            })
            .collect();
        TiltYear {
            years,
            range: lockstep.range,
            growth: lockstep.growth,
            partials: lockstep.partials,
            air: lockstep.air,
            masses: lockstep.masses,
            surfaces: lockstep.surfaces,
            caps: lockstep.caps,
        }
    }

    /// The cycles of a resonant body under its atmosphere (F239, F245): its places run in
    /// step through the orbits it repeats in, along the star's `path`.
    fn seasonal_orbits(
        &self,
        ground: &Ground,
        layout: &Layout,
        (orbits, path): (u32, &[PathStep]),
        ocean: Option<&Ocean>,
    ) -> TiltYear {
        let cycle = self.year * orbits as f64;
        let column = Column::new(ground, cycle, cycle);
        let e = self.sky.eccentricity;
        let scale = self.sky.flux / (1.0 - e * e).sqrt();
        let light = |k: usize, p: &Place| {
            let (s, (sp, cp)) = (&path[k], p.lat.sin_cos());
            scale * s.flux * (sp * s.sin_decl + cp * s.cos_decl * (s.hour + p.lon).cos()).max(0.0)
        };
        let steps = path.len();
        let snapshot = snapshot_step(self.snapshot, steps / orbits.max(1) as usize);
        let timing = (steps, cycle / steps as f64, snapshot);
        let year = (0, None, CYCLE_MAX_YEARS);
        let lockstep = self.lockstep(ground, &column, layout, timing, light, ocean, year);
        let years = lockstep
            .runs
            .iter()
            .map(|r| {
                let mean = r.sums[0] / steps as f64;
                let mut samples: Vec<f64> = r.surface.iter().map(|&t| t as f64).collect();
                Year {
                    mean,
                    max: r.hi,
                    min: r.lo,
                    frost_share: r.frosted as f64 / steps as f64,
                    frost_now: r.now,
                    sea: r.sea_year(steps, (mean, r.hi, r.lo)),
                    spread: spread_of(&mut samples),
                    albedo: r.albedo as f32,
                }
            })
            .collect();
        let orbits = orbits.max(1) as f64;
        TiltYear {
            years,
            range: lockstep.range,
            growth: lockstep.growth.map(|g| g / orbits),
            partials: lockstep.partials,
            air: lockstep.air,
            masses: lockstep.masses,
            surfaces: lockstep.surfaces,
            caps: lockstep.caps.iter().map(|g| g / orbits).collect(),
        }
    }

    /// Places run in step through their cycles under one atmosphere (F239, F245).
    ///
    /// The air over each place holds its column's heat (F236) and exchanges heat with the
    /// ground under it, its land's and its sea's by their shares of the place
    /// (`land_coupling`, `sea_coupling_at`, at the place's contact with its ground, F269),
    /// and with the air beside it by the transport, which the cycle answers
    /// through its seasons as the year's balance does through its year (F235): at each step
    /// the air is solved implicitly from the ground as the last step left it, then the ground
    /// under it, each place's land, its column's surface exchanging heat with the air, and its
    /// sea (`sea_step`). The transport reaches the ground only through the air, so that a
    /// thin atmosphere, whose air barely touches its ground, carries little of it there.
    ///
    /// The frost of the atmosphere's gas draws on that one atmosphere: the frost the cold
    /// ground takes is the atmosphere's own, so that its pressure, and the frost point with
    /// it, falls as the caps grow and rises as they sublime (Leighton and Murray 1966). At
    /// each step every place tries its step, and the point is the one at which the frost they
    /// take together leaves the gas at the pressure whose point it is: a thin atmosphere's
    /// pressure follows its frost's balance through the step. The gas's whole mass is set so
    /// that its pressure over the cycle averages the body's where the frost allows it. Frost
    /// that lies all through the cycle and grows is a cap's growth, which the caps keep
    /// (F241) where the gas is given back: it is counted, and taken off in the share given.
    ///
    /// The frost's own gas's caps are the reservoir its air draws on (F258): each place's
    /// share of them stands on its own column at the frost point through the cycle, gaining
    /// what the air leaves it and giving the air what it sublimes, so that a cap in its summer
    /// holds the air up against the frost its winter lays elsewhere. What they gain or lose
    /// through the last cycle is counted, the gas they give the air held in the gas's whole
    /// mass as their frost is.
    ///
    /// The cycle is `steps` steps of `dt` with the snapshot's step, each place lit at a step
    /// as `light` gives, the ground recorded every `every` steps (never at 0), the places
    /// started from the latitudes' rows `warm` left where it is given, through at most
    /// `years` cycles.
    #[allow(clippy::too_many_arguments)]
    fn lockstep(
        &self,
        ground: &Ground,
        column: &Column,
        layout: &Layout,
        (steps, dt, snapshot): (usize, f64, usize),
        light: impl Fn(usize, &Place) -> f64 + Sync,
        ocean: Option<&Ocean>,
        (every, warm, years): (usize, Option<&Warm>, usize),
    ) -> Lockstep {
        let places = &layout.places;
        let n = places.len();
        let frost = self.frost.as_ref();
        let whole = places.iter().map(|p| p.weight).sum::<f64>().max(1e-30);
        let over = |runs: &[Run], pick: fn(&Run) -> f64| {
            runs.iter()
                .zip(places)
                .map(|(r, p)| pick(r) * p.weight)
                .sum::<f64>()
                / whole
        };
        let column_gas = frost.map_or(1.0, |f| f.column.max(1e-30));
        let partial = |air: f64| frost.map_or(0.0, |f| f.partial * air.max(0.0) / column_gas);
        let mut runs: Vec<Run> = places
            .iter()
            .zip(&layout.snowfall)
            .map(|(p, fall)| Run::new((column.len(), steps), p, ocean, frost, fall))
            .collect();
        let capped = runs.iter().any(|r| r.cap.is_some());
        let mut air: Vec<f64> = places.iter().map(|p| p.steady).collect();
        let mut series = vec![0.0; steps * n];
        let (mut diag, mut rhs) = (vec![0.0; n], vec![0.0; n]);
        let hold = ground.air / dt;
        let water = ocean.map_or(0.0, |_| OCEAN_TRANSPORT_SHARE);
        let air_links = layout.air_links(self.transport, water);
        let mut inflow = vec![0.0; n];
        // The sea ice's drift, where the year's wind is known on these places (F276).
        let mut drift = ocean
            .filter(|_| n > 0 && layout.n_lon > 0 && layout.wind.len() == steps * n)
            .map(|o| {
                let sea: Vec<f64> = places.iter().map(|p| p.sea).collect();
                let weight: Vec<f64> = places.iter().map(|p| p.weight).collect();
                let grid = (n / layout.n_lon, layout.n_lon);
                let spin = (self.circulation.omega, self.circulation.radius);
                let water = self.liquid.map_or(1_000.0, |l| l.liquid_density);
                crate::drift::Drift::new((&sea, &weight), grid, spin, (water, o.density))
            });
        let (mut held, mut load) = (frost.map_or(0.0, |f| f.column), 0.0);
        if let Some(w) = warm {
            for (r, p) in runs.iter_mut().zip(places) {
                let at = (p.lat + PI / 2.0) / PI * (w.rows.len() - 1) as f64;
                let row = &w.rows[(at.round() as usize).min(w.rows.len() - 1)];
                r.t.clone_from(&row.t);
                r.mass = row.mass;
                if let (Some(c), Some(t)) = (r.cap.as_mut(), &row.cap) {
                    c.t.clone_from(t);
                }
                if let (Some(sea), Some(from)) = (r.sea.as_mut(), row.sea) {
                    *sea = from;
                }
            }
            held = w.active;
            load = over(&runs, Run::load);
        }
        let mut range = (f64::MAX, f64::MIN);
        let mut growth = 0.0;
        let mut partials = Vec::with_capacity(steps);
        let mut last = vec![f64::NAN; 2 * n];
        for _ in 0..years {
            range = (f64::MAX, f64::MIN);
            partials.clear();
            runs.iter_mut().for_each(Run::begin);
            for k in 0..steps {
                // The air, from the ground under it as the last step left it: its land's, its
                // sea's and its caps'.
                for (q, (r, p)) in runs.iter().zip(places).enumerate() {
                    let contact = &p.contact;
                    let land = self.land_coupling(r.t[0], contact) * p.land();
                    let (wet, surface) = match &r.sea {
                        Some(s) => (self.sea_coupling_at(s.surface, contact) * p.sea, s.surface),
                        None => (0.0, 0.0),
                    };
                    let (iced, top) = match &r.cap {
                        Some(c) => (self.coupling_at(c.t[0], contact) * c.share, c.t[0]),
                        None => (0.0, 0.0),
                    };
                    // The snow that falls on its land froze in it, its latent heat of fusion
                    // the air's, as the ground's melting of it takes it back (F271).
                    let frozen = match (self.liquid.filter(|_| self.ice), r.snowfall.get(k)) {
                        (Some(l), Some(&fall)) => l.l_fus * fall as f64 * p.land(),
                        _ => 0.0,
                    };
                    diag[q] = hold + land + wet + iced;
                    rhs[q] = hold * air[q] + land * r.t[0] + wet * surface + iced * top + frozen;
                }
                transport_solve(layout, &air_links, (&diag, &rhs), &mut air, true);
                series[k * n..(k + 1) * n].copy_from_slice(&air);
                if let Some(ocean) = ocean.filter(|_| water > 0.0) {
                    // Where the seas' currents carry their heat (F272), the diffusion carries
                    // the water's share of the transport the wind does not drive.
                    let undriven = match layout.currents {
                        Some(_) => OCEAN_UNDRIVEN_SHARE,
                        None => 1.0,
                    };
                    let step = (self.transport * water * undriven, dt);
                    sea_transport(layout, &mut runs, step, ocean, &mut inflow);
                    if let Some(currents) = &layout.currents {
                        sea_currents(layout, &mut runs, (currents, dt), ocean, &mut inflow);
                    }
                }
                let highest = frost.map(|f| f.point_at(partial(held)));
                // Bare ground stays clear of frost above every point the air can reach, the
                // point of all its gas; with caps that hold the air, which may give it more,
                // above the point of the air as it stands, the step tried again past it if
                // the point settles above (F258).
                let ceiling = match (frost, highest) {
                    (Some(f), Some(h)) if capped => Some(h.max(f.point_at(partial(held - load)))),
                    _ => highest,
                };
                runs.par_iter_mut()
                    .zip(places.par_iter())
                    .zip(air.par_iter().zip(inflow.par_iter()))
                    .for_each(|((r, p), (&a, &brought))| {
                        let lit = light(k, p);
                        // The snow that falls on the land through the step lies on it, and
                        // brightens it by its cover where the year's snow is known (F271).
                        r.albedo = match r.snowfall.get(k) {
                            Some(&fall) => {
                                r.snow = (r.snow + dt * fall as f64).min(SNOW_MOST_MASS);
                                self.snow_albedo((p.albedo, p.snow), r.snow, r.t[0])
                            }
                            None => self.land_albedo((p.albedo, p.snow), r.t[0]),
                        };
                        r.start.copy_from_slice(&r.t);
                        let bare = (1.0 - r.albedo) * lit;
                        let shed = self.land_shed(a, p.contact);
                        let heat = self.heat_flux;
                        match (frost, ceiling) {
                            (Some(f), Some(top)) => {
                                let gains = (bare, (1.0 - f.albedo) * lit);
                                let rows = (&mut r.work, &mut r.trial);
                                column.trial(&r.t, rows, dt, gains, heat, shed, r.mass, top);
                                if let Some(c) = r.cap.as_mut() {
                                    let rows = (&mut c.work, &mut c.trial);
                                    column.trial(&c.t, rows, dt, gains, heat, shed, c.mass, top);
                                }
                            }
                            _ => column.step(&mut r.t, &mut r.work, dt, bare, heat, shed),
                        }
                        if let (Some(sea), Some(ocean)) = (r.sea.as_mut(), ocean) {
                            let exchange = self.sea_coupling_at(sea.surface, &p.contact);
                            let water = (ocean, brought);
                            self.sea_step(sea, dt, (lit, p.albedo), (a, exchange), water);
                        }
                    });
                if let (Some(drift), Some(ocean)) = (drift.as_mut(), ocean) {
                    let wind = &layout.wind[k * n..(k + 1) * n];
                    self.drift_ice(drift, &mut runs, (wind, &air), dt, ocean);
                }
                let point = match (frost, highest, ceiling) {
                    (Some(f), Some(highest), Some(top)) => {
                        let gas = (held - load, whole);
                        let mut point =
                            self.frost_point(column, (&runs, places, &air), gas, dt, f, highest);
                        if point > top {
                            runs.par_iter_mut()
                                .zip(places.par_iter())
                                .zip(air.par_iter())
                                .for_each(|((r, p), &a)| {
                                    if r.mass <= 0.0 && r.trial.clear && r.trial.free[0] < point {
                                        let lit = light(k, p);
                                        let gains =
                                            ((1.0 - r.albedo) * lit, (1.0 - f.albedo) * lit);
                                        let shed = self.land_shed(a, p.contact);
                                        let heat = self.heat_flux;
                                        let rows = (&mut r.work, &mut r.trial);
                                        let open = f64::INFINITY;
                                        column.trial(&r.t, rows, dt, gains, heat, shed, 0.0, open);
                                    }
                                });
                            let runs = (&runs[..], &places[..], &air[..]);
                            point = self.frost_point(column, runs, gas, dt, f, highest);
                        }
                        let (olr, slope) = self.olr_and_slope(point);
                        let (e, de) = f.emissivity(point);
                        let each = runs.par_iter_mut().zip(places.par_iter());
                        each.zip(air.par_iter()).for_each(|((r, p), &a)| {
                            let exchange = self.coupling_at(point, &p.contact);
                            let shed = (
                                e * olr + exchange * (point - a),
                                e * slope + de * olr + exchange,
                            );
                            let lasted = column.finish(
                                &mut r.t,
                                &r.trial,
                                dt,
                                point,
                                shed,
                                f.latent,
                                &mut r.mass,
                            );
                            // The frost that ran out within the step leaves its ground the
                            // rest of it bare, to frost again (F275).
                            if let Some(lasts) = lasted {
                                let gains = ((1.0 - r.albedo) * light(k, p), self.heat_flux);
                                let bare = self.land_shed(a, p.contact);
                                let rows = (&mut r.work, &mut r.trial);
                                let ends = (point, shed.0, f.latent);
                                let (t, start) = (&mut r.t, &r.start);
                                let step = (dt, lasts);
                                column.rest_bare(
                                    t,
                                    start,
                                    rows,
                                    step,
                                    gains,
                                    bare,
                                    ends,
                                    &mut r.mass,
                                );
                            }
                            if let Some(c) = r.cap.as_mut() {
                                let (t, mass) = (&mut c.t, &mut c.mass);
                                column.finish(t, &c.trial, dt, point, shed, f.latent, mass);
                            }
                        });
                        load = over(&runs, Run::load);
                        let now = partial(held - load);
                        range = (range.0.min(now), range.1.max(now));
                        partials.push(now);
                        point
                    }
                    _ => f64::NAN,
                };
                // The land's snow (F271): where it lies and the step left the ground above its
                // melting point, the ground stood at the point through the step instead and the
                // step's warmth melted the snow.
                if let Some(liquid) = self.liquid.filter(|_| self.ice) {
                    let melt = (liquid.t_triple, liquid.l_fus);
                    runs.par_iter_mut()
                        .zip(places.par_iter())
                        .zip(air.par_iter())
                        .for_each(|((r, p), &a)| {
                            r.melted = 0.0;
                            if r.snow > 0.0 && r.t[0] > melt.0 {
                                let gains = ((1.0 - r.albedo) * light(k, p), self.heat_flux);
                                let shed = self.land_shed(a, p.contact)(melt.0);
                                let rows = (&mut r.work, &mut r.trial);
                                let (t, snow) = (&mut r.t, &mut r.snow);
                                let before = *snow;
                                column.thaw(t, &r.start, rows, dt, gains, shed, melt, snow);
                                r.melted = before - r.snow;
                            }
                        });
                }
                runs.par_iter_mut()
                    .zip(air.par_iter())
                    .for_each(|(r, &a)| r.record(k, snapshot, every, point, a));
            }
            // The frost that lay all cycle and grew through it: a cap's growth; and what the
            // caps that hold the air gained through it (F258).
            growth = over(&runs, |r| r.growth() * r.bare_share());
            for r in runs.iter_mut() {
                if let Some(c) = r.cap.as_mut() {
                    c.gained = c.mass - c.start;
                }
            }
            let means: Vec<f64> = runs
                .iter()
                .flat_map(|r| [r.sums[0] / steps as f64, r.sea_sum / steps as f64])
                .collect();
            let settled = means
                .iter()
                .zip(&last)
                .all(|(m, l)| (m - l).abs() < CYCLE_SETTLED_K);
            last = means;
            if settled {
                break;
            }
            // The gas's whole mass: the atmosphere's, and the frost's that comes and goes in
            // the share of the body's pressure the air keeps through the cycle: all of it
            // where the frost leaves the air its pressure, so that the air's mean comes to the
            // body's; little where the frost's own balance holds the air below it, as no store
            // of gas would raise it. A cap's growth comes off in the same share, as the gas it
            // took is given back: the cap a closed atmosphere's frost lays stays.
            //
            // Caps that hold the air are such a store (F258): the air's mean comes to the
            // body's, the caps going on giving or taking through the next cycle what they did
            // through this one, their kilometres barely touched.
            if let Some(f) = frost {
                let mean_load = over(&runs, |r| r.frost_sum) / steps as f64;
                if capped {
                    let drift = over(&runs, |r| {
                        r.cap.as_ref().map_or(0.0, |c| c.share * c.gained)
                    });
                    for r in runs.iter_mut() {
                        r.mass -= r.growth();
                    }
                    load = over(&runs, Run::load);
                    held = (f.column + mean_load + drift - growth).max(load);
                } else {
                    let kept = ((held - mean_load) / column_gas).clamp(0.0, 1.0);
                    for r in runs.iter_mut() {
                        r.mass -= kept * r.growth();
                    }
                    load = over(&runs, Run::load);
                    held = (f.column + (mean_load - growth) * kept).max(load);
                }
            }
            for r in runs.iter_mut() {
                let sums = std::mem::take(&mut r.sums);
                settle(&mut r.t, &sums, steps);
                r.sums = sums;
                if let Some(c) = r.cap.as_mut() {
                    settle(&mut c.t, &c.sums, steps);
                }
            }
        }
        if let (Some(d), true) = (&drift, std::env::var_os("EXOGEN_DEBUG").is_some()) {
            let (speed, spent, subs, sweeps) = d.summary();
            eprintln!(
                "drift: the sea ice's mean speed {speed:.3} m/s; {spent:.2} s, {subs:.1} sub-steps and {sweeps:.1} sweeps a step"
            );
        }
        // The surface pressure: the other gases' and the frost's gas's.
        let range = frost.map(|f| {
            let others = (self.pressure - f.partial).max(0.0);
            (others + range.0, others + range.1)
        });
        let masses = std::array::from_fn(|j| {
            runs.iter()
                .flat_map(|r| r.masses[j].iter().copied())
                .collect()
        });
        let surfaces = [
            runs.iter()
                .flat_map(|r| r.surface.iter().copied())
                .collect(),
            runs.iter()
                .flat_map(|r| r.sea_surface.iter().copied())
                .collect(),
        ];
        let caps = runs
            .iter()
            .map(|r| r.cap.as_ref().map_or(0.0, |c| c.share * c.gained))
            .collect();
        let active = held - over(&runs, Run::taken);
        Lockstep {
            runs,
            range,
            growth: frost.map(|_| growth),
            partials,
            air: series,
            masses,
            surfaces,
            caps,
            active,
        }
    }

    /// The frost point a step of the lockstep settles at (F239): the point at which the frost
    /// the places take through the step leaves the gas, `gas` kg/m² over the body before it,
    /// at the partial pressure whose point it is, each place's ground shedding its heat to
    /// space and to the air over it (F245). The frost taken rises with the point and the
    /// point of what is left falls, so the two cross once, between none and `highest`, the
    /// point of all the gas; found by false position, halving the side that stays, the low
    /// side, which leaves the air the gas the frost does not take. The caps that hold the air
    /// may give it more than all the gas (F258), and the point then stands above `highest`.
    fn frost_point(
        &self,
        column: &Column,
        (runs, places, air): (&[Run], &[Place], &[f64]),
        (gas, whole): (f64, f64),
        dt: f64,
        frost: &Frost,
        highest: f64,
    ) -> f64 {
        // Each place's ground bare of the caps where it takes or gives frost, and its caps:
        // its trial, its frost, its share of the body and the air over it.
        let mut taking: Vec<(&Trial, f64, f64, f64, Contact)> = Vec::new();
        for ((r, p), &a) in runs.iter().zip(places).zip(air) {
            let share = p.weight * r.bare_share() / whole;
            if r.mass > 0.0 || !r.trial.clear {
                taking.push((&r.trial, r.mass, share, a, p.contact));
            }
            if let Some(c) = &r.cap {
                taking.push((&c.trial, c.mass, p.weight * c.share / whole, a, p.contact));
            }
        }
        let gap = |point: f64| {
            let olr = self.olr(point) * frost.emissivity(point).0;
            let taken: f64 = taking
                .iter()
                .map(|&(trial, mass, w, a, contact)| {
                    let exchange = self.coupling_at(point, &contact);
                    let radiated = olr + exchange * (point - a);
                    column
                        .deposit(trial, dt, point, radiated, frost.latent, mass)
                        .1
                        * w
                })
                .sum();
            point - frost.point_at(frost.partial * (gas - taken).max(0.0) / frost.column.max(1e-30))
        };
        let (mut lo, mut hi) = (0.0, highest);
        let (mut g_lo, mut g_hi) = (gap(lo), gap(hi));
        if g_hi <= 0.0 {
            if runs.iter().all(|r| r.cap.is_none()) {
                return hi;
            }
            // The point of what the air holds at `highest`, past the crossing.
            (lo, g_lo) = (hi, g_hi);
            hi -= g_hi;
            g_hi = gap(hi);
            if g_hi <= 0.0 {
                return hi;
            }
        }
        let mut side = 0;
        for _ in 0..FROST_POINT_ITERATIONS {
            let mid = if g_hi > g_lo {
                (lo * g_hi - hi * g_lo) / (g_hi - g_lo)
            } else {
                0.5 * (lo + hi)
            };
            let g = gap(mid);
            if g <= 0.0 && g > -FROST_POINT_TOLERANCE_K {
                return mid;
            }
            if g > 0.0 {
                (hi, g_hi) = (mid, g);
                if side > 0 {
                    g_lo *= 0.5;
                }
                side = 1;
            } else {
                (lo, g_lo) = (mid, g);
                if side < 0 {
                    g_hi *= 0.5;
                }
                side = -1;
            }
            if hi - lo < FROST_POINT_TOLERANCE_K {
                break;
            }
        }
        lo
    }

    /// Where each ice would last (F240): every cell of the climate's grid runs a column of
    /// the ice through its year at its height's pressure, `heights` m above the datum's air
    /// thinning it by the scale height, under the frost's seasons (F239), at each tilt the
    /// year was run at, under the cycle's air over it, as much colder at its height as the
    /// ground is (F237), which exchanges heat with the ice as with the ground about it
    /// (F245). The frost's own gas's ice lies all year from an
    /// unlimited store, held at its point at the cell's partial pressure, and gains or loses what the surface's
    /// balance leaves (Leighton and Murray 1966). Any other ice lies under the frost where
    /// the frost lays it, its albedo its surface's, and its trap is the year-mean saturation
    /// pressure over it through the days and seasons, divided by the share of the datum's
    /// pressure its height keeps, so that it is the vapour pressure the datum's air must hold
    /// to keep it (Ingersoll 1970).
    pub fn traps(
        &self,
        climate: &Climate,
        ground: &Ground,
        cycle: &Cycle,
        ices: &[Ice],
        heights: &Field<f32>,
        downscale: &Downscale,
    ) -> Traps {
        let grid = climate.temperature.grid();
        let thins: Vec<f64> = heights
            .as_slice()
            .iter()
            .map(|&h| {
                let h = if h.is_finite() { h as f64 } else { 0.0 };
                if downscale.scale_height > 0.0 {
                    (-h / downscale.scale_height).exp()
                } else {
                    1.0
                }
            })
            .collect();
        let (_, n_lon) = cycle.layout;
        let cell_at = |lat: f64, lon: f64, thin: f64| {
            self.trap_cell((climate, cycle, downscale), lat, lon, thin)
        };
        let frost = cycle.frost.as_ref();
        // A fast spinner's traps are its latitudes' at a few heights, the frost's point the
        // only thing a height changes, read at each cell between them; a resonant body's are
        // each cell's own.
        let levels: Option<Vec<f64>> = match n_lon {
            1 => {
                let (lo, hi) = thins
                    .iter()
                    .fold((f64::MAX, f64::MIN), |(a, b), &t| (a.min(t), b.max(t)));
                let k = if frost.is_some() && hi > lo * (1.0 + 1e-6) {
                    TRAP_HEIGHT_LEVELS
                } else {
                    1
                };
                Some(
                    (0..k)
                        .map(|j| match k {
                            1 => 1.0,
                            _ => (lo.ln() + (hi / lo).ln() * j as f64 / (k - 1) as f64).exp(),
                        })
                        .collect(),
                )
            }
            _ => None,
        };
        let places: Vec<TrapCell> = match &levels {
            Some(levels) => (0..CYCLE_LATITUDES)
                .flat_map(|i| {
                    let lat = -PI / 2.0 + PI * i as f64 / (CYCLE_LATITUDES - 1) as f64;
                    levels.iter().map(move |&thin| (lat, thin))
                })
                .map(|(lat, thin)| cell_at(lat, 0.0, thin))
                .collect(),
            None => (0..grid.cell_count())
                .map(|k| {
                    let (lon, lat) = lon_lat_from_direction(grid.center(grid.cell(k)));
                    cell_at(lat, lon, thins[k])
                })
                .collect(),
        };
        let mut traps = Traps {
            growth: None,
            vapour: Vec::new(),
            warmest: Vec::new(),
        };
        for ice in ices {
            let own = frost.filter(|f| f.species == ice.species);
            let runs: Vec<(f64, f64)> = places
                .par_iter()
                .map(|cell| self.trap_value(ground, cycle, cell, ice))
                .collect();
            let values: Vec<f64> = runs.iter().map(|r| r.0).collect();
            let hottest: Vec<f64> = runs.iter().map(|r| r.1).collect();
            let read = |values: &[f64], c: crate::grid::Cell, d: [f64; 3], log: bool| {
                let k = grid.index(c);
                match &levels {
                    Some(levels) => {
                        let (_, lat) = lon_lat_from_direction(d);
                        read_table(values, levels, lat, thins[k], log)
                    }
                    None => values[k],
                }
            };
            match own {
                Some(_) => {
                    let field = Field::from_fn(grid, |c, d| read(&values, c, d, false) as f32);
                    traps.growth = Some((*ice, field));
                }
                None => {
                    // The vapour at the datum's air: the saturation pressure over the ice over
                    // the share of that air's pressure its height keeps.
                    let field = Field::from_fn(grid, |c, d| {
                        (read(&values, c, d, true) / thins[grid.index(c)]) as f32
                    });
                    let warmest = Field::from_fn(grid, |c, d| read(&hottest, c, d, false) as f32);
                    traps.vapour.push((*ice, field));
                    traps.warmest.push(warmest);
                }
            }
        }
        traps
    }

    /// A place of the traps (F240) at a latitude and longitude, rad, and a share of the
    /// datum's pressure its height keeps. Its air is the cycle's over its place, as much colder
    /// the higher it stands as its ground is (F237): near its ground the air follows the
    /// ground, a thin atmosphere's barely cooling with height; its height the one its
    /// pressure's share gives.
    fn trap_cell(
        &self,
        (climate, cycle, downscale): (&Climate, &Cycle, &Downscale),
        lat: f64,
        lon: f64,
        thin: f64,
    ) -> TrapCell {
        let grid = climate.temperature.grid();
        let (n_lat, n_lon) = cycle.layout;
        let (c, _, _) = grid.locate(direction_from_lon_lat(lon, lat));
        let height = if thin > 0.0 && thin != 1.0 {
            -downscale.scale_height * thin.ln()
        } else {
            0.0
        };
        let steady = climate.temperature.get(c) as f64;
        TrapCell {
            lat,
            lon,
            thin,
            steady,
            air: reader(n_lat, n_lon, lat, lon),
            cooling: steady - downscale.at(steady, steady, height),
        }
    }

    /// An ice's trap at a place through the tilts' mean (F240), and its warmest at any: the
    /// frost's own gas's ice's growth, kg/m² a year, every other ice's year-mean saturation
    /// pressure over it, Pa, at the place's own pressure.
    fn trap_value(&self, ground: &Ground, cycle: &Cycle, cell: &TrapCell, ice: &Ice) -> (f64, f64) {
        let frost = cycle.frost.as_ref();
        let own = frost.filter(|f| f.species == ice.species);
        let n = cycle.tilts.len().max(1) as f64;
        cycle
            .tilts
            .iter()
            .enumerate()
            .fold((0.0, f64::MIN), |(sum, hot), (j, &eps)| {
                let partials = cycle.partials.get(j).map_or(&[][..], |p| &p[..]);
                let air = cycle.air.get(j).map_or(&[][..], |a| &a[..]);
                let along = (partials, air);
                let (value, warmest) = match own {
                    Some(f) => (self.ice_growth(ground, cell, eps, along, f), f64::NAN),
                    None => self.ice_vapour(ground, cell, eps, along, frost, ice),
                };
                (sum + value / n, hot.max(warmest))
            })
    }

    /// The cold traps on the cycle's latitudes' rows at `levels`, shares of the datum's
    /// pressure (F259), the year run on its rows: each ice's trap at each latitude of the rows
    /// and each level, latitude after latitude, as `trap_value` finds it there.
    pub fn trap_table(
        &self,
        (climate, ground, cycle): (&Climate, &Ground, &Cycle),
        ices: &[Ice],
        levels: &[f64],
        downscale: &Downscale,
    ) -> Vec<(Ice, Vec<f64>)> {
        let places: Vec<TrapCell> = (0..CYCLE_LATITUDES)
            .flat_map(|i| {
                let lat = -PI / 2.0 + PI * i as f64 / (CYCLE_LATITUDES - 1) as f64;
                levels.iter().map(move |&thin| (lat, thin))
            })
            .map(|(lat, thin)| self.trap_cell((climate, cycle, downscale), lat, 0.0, thin))
            .collect();
        ices.iter()
            .map(|ice| {
                let values = places
                    .par_iter()
                    .map(|cell| self.trap_value(ground, cycle, cell, ice).0)
                    .collect();
                (*ice, values)
            })
            .collect()
    }

    /// The sea ice's drift through a step (F276): the wind's stress on the ice at each place,
    /// rho_a C_a |V| V at the density of the air over it, its drag the law of the wall's over
    /// the ice's roughness; the ice the drift brings to open water melts in the water's warmth
    /// over its freezing point, and a sea the drift leaves bare stands open at it.
    fn drift_ice(
        &self,
        drift: &mut crate::drift::Drift,
        runs: &mut [Run],
        (wind, air): (&[[f32; 2]], &[f64]),
        dt: f64,
        ocean: &Ocean,
    ) {
        let before: Vec<f64> = runs
            .iter()
            .map(|r| r.sea.as_ref().map_or(0.0, |s| s.ice.max(0.0)))
            .collect();
        if before.iter().all(|&m| m <= 0.0) {
            return;
        }
        let drag = (VON_KARMAN / (WIND_REFERENCE_HEIGHT / ICE_ROUGHNESS).ln()).powi(2);
        let c = &self.circulation;
        let stress: Vec<[f64; 2]> = wind
            .iter()
            .zip(air)
            .map(|(v, &t)| {
                let rho = c.pressure * c.molar_mass / (GAS_CONSTANT * t.max(1.0));
                let (u, w) = (v[0] as f64, v[1] as f64);
                let pull = rho * drag * u.hypot(w);
                [pull * u, pull * w]
            })
            .collect();
        let mut ice = before.clone();
        drift.step(&mut ice, &stress, dt);
        for ((r, &was), &now) in runs.iter_mut().zip(&before).zip(&ice) {
            let Some(sea) = r.sea.as_mut() else {
                continue;
            };
            if now == was {
                continue;
            }
            if was > 0.0 {
                sea.ice = now;
                if now <= 0.0 {
                    (sea.t, sea.surface) = (ocean.freezing, ocean.freezing);
                }
                continue;
            }
            // Ice brought to open water melts in its warmth over the freezing point.
            let melts = ocean.capacity * (sea.t - ocean.freezing).max(0.0) / ocean.latent;
            if now <= melts {
                sea.t -= now * ocean.latent / ocean.capacity;
                sea.surface = sea.t;
            } else {
                sea.ice = now - melts;
                (sea.t, sea.surface) = (ocean.freezing, ocean.freezing);
            }
        }
    }

    /// One step of a sea's mixed layer and its ice (F244, F245), lit by `lit` W/m², its open
    /// water at its ground's albedo `open`, under air at `air` K exchanging `exchange` W/m²/K
    /// with it. Open, the layer takes its sunlight, radiates and exchanges heat with the air,
    /// implicitly; cooled past its freezing point, the cold it is left with freezes ice at the
    /// ice's latent heat, and the layer stays at the point. Under ice, the ice's top stands
    /// where its sunlight, its radiation, the air's heat and the heat conducted up through the
    /// ice balance (Semtner 1976's zero layers), above the point no further, its excess there
    /// melting the ice from the top; it grows at its base by what it conducts, less what the
    /// sea brings it, its own heat from below and `brought` W/m², the water's share of the
    /// transport (`sea_transport`); ice melted through gives the layer the rest. The ice
    /// reflects as the year's balance's ice and snow do (F235), as melting ice once it melts,
    /// under the star's light (F250), over the share of the sea it covers: thin ice is dark.
    fn sea_step(
        &self,
        sea: &mut Sea,
        dt: f64,
        (lit, open): (f64, f64),
        (air, exchange): (f64, f64),
        (ocean, brought): (&Ocean, f64),
    ) {
        let melting = sea.surface >= ocean.freezing - SEA_ICE_MELT_MARGIN_K;
        let ice = self
            .ice_over(
                if melting {
                    SEA_ICE_MELT_ALBEDO
                } else {
                    ICE_ALBEDO
                },
                open,
            )
            .max(open);
        let cover = 1.0 - (-sea.ice.max(0.0) / SEA_ICE_COVER_MASS).exp();
        let frozen = (1.0 - open - (ice - open) * cover) * lit;
        if sea.ice <= 0.0 {
            let (r, slope) = self.olr_and_slope(sea.t);
            let gain = (1.0 - open) * lit - r + exchange * (air - sea.t);
            let t = sea.t + dt * gain / (ocean.capacity + dt * (slope + exchange));
            if t >= ocean.freezing {
                (sea.t, sea.surface) = (t, t);
                return;
            }
            sea.ice = ocean.capacity * (ocean.freezing - t) / ocean.latent;
            (sea.t, sea.surface) = (ocean.freezing, ocean.freezing);
            return;
        }
        let through = ocean.conductivity / (sea.ice / ocean.density).max(SEA_ICE_THINNEST_M);
        let mut top = sea.surface.min(ocean.freezing);
        for _ in 0..SEA_ICE_NEWTON_STEPS {
            let (r, slope) = self.olr_and_slope(top);
            let gap = frozen - r + exchange * (air - top) + through * (ocean.freezing - top);
            top = (top + gap / (slope + exchange + through)).max(COSMIC_BACKGROUND_K);
        }
        let mut melted = 0.0;
        if top > ocean.freezing {
            top = ocean.freezing;
            let left = frozen - self.olr(top) + exchange * (air - top);
            melted = left.max(0.0) * dt / ocean.latent;
        }
        let below = ocean.base_flux + brought;
        let grown = (through * (ocean.freezing - top) - below) * dt / ocean.latent;
        sea.ice += grown - melted;
        sea.surface = top;
        if sea.ice < 0.0 {
            sea.t = ocean.freezing - sea.ice * ocean.latent / ocean.capacity;
            (sea.ice, sea.surface) = (0.0, sea.t);
        }
    }

    /// A cell's sunlight through its cycle at a tilt (F240): a fast spinner's year in
    /// CYCLE_YEAR_STEPS steps on the day's mean, with each step's flux and declination for
    /// its days; a resonant body's orbits along its path. Returns the steps' sunlight, their
    /// flux and declination where the days are run, the step's length and the cycle's in
    /// years.
    fn trap_sun(&self, cell: &TrapCell, eps: f64) -> (Vec<f64>, Vec<(f64, f64)>, f64, f64) {
        match self.sky.resonance {
            None => {
                let sun: Vec<(f64, f64)> = (0..CYCLE_YEAR_STEPS)
                    .map(|m| {
                        let anomaly = TAU * (m as f64 + 0.5) / CYCLE_YEAR_STEPS as f64;
                        self.sky.at_mean_anomaly(anomaly, eps)
                    })
                    .collect();
                let light = sun
                    .iter()
                    .map(|&(flux, decl)| flux * diurnal_share(cell.lat, decl))
                    .collect();
                (light, sun, self.year / CYCLE_YEAR_STEPS as f64, 1.0)
            }
            Some((_, orbits)) => {
                let path = self.resonant_steps(eps);
                let e = self.sky.eccentricity;
                let scale = self.sky.flux / (1.0 - e * e).sqrt();
                let (sp, cp) = cell.lat.sin_cos();
                let light = path
                    .iter()
                    .map(|s| {
                        let cz = sp * s.sin_decl + cp * s.cos_decl * (s.hour + cell.lon).cos();
                        scale * s.flux * cz.max(0.0)
                    })
                    .collect();
                let cycle = self.year * orbits as f64;
                (light, Vec::new(), cycle / path.len() as f64, orbits as f64)
            }
        }
    }

    /// The column a cell's ice runs on at a tilt, and its day: a fast spinner's resolves its
    /// day, a resonant body's its whole cycle.
    fn trap_column(&self, ground: &Ground, eps: f64) -> (Column, f64) {
        match self.sky.resonance {
            None => {
                let solar_days = (self.sky.spins - eps.cos()).abs().max(1e-3);
                let day = self.year / solar_days;
                (Column::new(ground, day, self.year), day)
            }
            Some((_, orbits)) => {
                let cycle = self.year * orbits as f64;
                (Column::new(ground, cycle, cycle), cycle)
            }
        }
    }

    /// The frost a cell's ice of the frost's own gas gains in a year at a tilt (F240), kg/m²:
    /// lying all year from an unlimited store, held at its point at the cell's partial
    /// pressure, it gains or loses what the surface's balance leaves (F239), under the cycle's
    /// air over the cell (F245), `along` the cycle's partial pressures and air at the tilt,
    /// the deep ground settled under it.
    fn ice_growth(
        &self,
        ground: &Ground,
        cell: &TrapCell,
        eps: f64,
        (partials, air): (&[f64], &[f64]),
        frost: &Frost,
    ) -> f64 {
        let (column, _) = self.trap_column(ground, eps);
        let (light, _, dt, years) = self.trap_sun(cell, eps);
        let n = column.len();
        let places = air.len() / light.len().max(1);
        let mut t = vec![cell.steady; n];
        let (mut work, mut trial) = (Work::new(n), Trial::new(n));
        let gain = |lit: f64| (1.0 - frost.albedo) * lit;
        let mut gained = 0.0;
        for _ in 0..TRAP_YEARS {
            let mut mass = TRAP_STORE;
            let mut sums = vec![0.0; n];
            for (k, &lit) in light.iter().enumerate() {
                let f = Frost {
                    point: frost
                        .point_at(partials.get(k).copied().unwrap_or(frost.partial) * cell.thin),
                    ..*frost
                };
                let a = cell.air_at(air, places, k);
                let sheds = (self.shed(a), self.frost_shed(a, &f));
                let rows = (&mut work, &mut trial);
                let gains = (gain(lit), gain(lit));
                column.frosty_step(
                    &mut t,
                    rows,
                    dt,
                    gains,
                    self.heat_flux,
                    sheds,
                    &f,
                    &mut mass,
                );
                for (s, v) in sums.iter_mut().zip(&t) {
                    *s += v;
                }
            }
            gained = (mass - TRAP_STORE) / years;
            settle(&mut t, &sums, light.len());
        }
        gained
    }

    /// The year-mean saturation pressure over a cell's ice of another gas at a tilt (F240),
    /// Pa, and the warmest it stands, K: the frost's gas laying frost on it where the cell's
    /// partial pressure's point allows (F239), under the cycle's air over the cell (F245),
    /// `along` the cycle's partial pressures and air at the tilt; run through its year until
    /// it repeats and then, a fast spinner's, through a day at each season, the air holding
    /// through the day.
    fn ice_vapour(
        &self,
        ground: &Ground,
        cell: &TrapCell,
        eps: f64,
        (partials, air): (&[f64], &[f64]),
        frost: Option<&Frost>,
        ice: &Ice,
    ) -> (f64, f64) {
        let (column, day) = self.trap_column(ground, eps);
        let (light, sun, dt, _) = self.trap_sun(cell, eps);
        let albedo = self.ice_albedo(ice);
        let n = column.len();
        let places = air.len() / light.len().max(1);
        let (mut t, mut mass) = (vec![cell.steady; n], 0.0);
        let (mut work, mut trial) = (Work::new(n), Trial::new(n));
        let mut last = f64::NAN;
        let mut seasons: Vec<(usize, Vec<f64>, f64)> = Vec::new();
        let (mut vapour, mut hot) = (0.0, f64::MIN);
        let every = (light.len() / CYCLE_SEASONS).max(1);
        for _ in 0..CYCLE_MAX_YEARS {
            seasons.clear();
            (vapour, hot) = (0.0, f64::MIN);
            let mut sums = vec![0.0; n];
            for (k, &lit) in light.iter().enumerate() {
                let f = cell.laid(frost, partials, k);
                let state = (&mut t, &mut mass, &mut work, &mut trial);
                let a = cell.air_at(air, places, k);
                self.ice_step(&column, state, dt, (lit, albedo), f.as_ref(), a);
                for (s, v) in sums.iter_mut().zip(&t) {
                    *s += v;
                }
                vapour += vapour_pressure(ice.volatile, t[0]);
                hot = hot.max(t[0]);
                if !sun.is_empty() && k.is_multiple_of(every) {
                    seasons.push((k, t.clone(), mass));
                }
            }
            let surface = sums[0] / light.len() as f64;
            if (surface - last).abs() < CYCLE_SETTLED_K {
                break;
            }
            last = surface;
            settle(&mut t, &sums, light.len());
        }
        if sun.is_empty() {
            return (vapour / light.len() as f64, hot);
        }
        // A fast spinner's days at each season, from the ice and its frost as the season left
        // them: the saturation pressure through the last of them.
        let (sp, cp) = cell.lat.sin_cos();
        let dt = day / CYCLE_DAY_STEPS as f64;
        let mut sum = 0.0;
        for (k, state, frost_mass) in &seasons {
            let (flux, decl) = sun[*k];
            let (sd, cd) = decl.sin_cos();
            let f = cell.laid(frost, partials, *k);
            let a = cell.air_at(air, places, *k);
            let (mut t, mut mass) = (state.clone(), *frost_mass);
            for d in 0..CYCLE_DAYS {
                for h in 0..CYCLE_DAY_STEPS {
                    let hour = TAU * (h as f64 + 0.5) / CYCLE_DAY_STEPS as f64 - PI;
                    let lit = flux * (sp * sd + cp * cd * hour.cos()).max(0.0);
                    let state = (&mut t, &mut mass, &mut work, &mut trial);
                    self.ice_step(&column, state, dt, (lit, albedo), f.as_ref(), a);
                    if d + 1 == CYCLE_DAYS {
                        sum += vapour_pressure(ice.volatile, t[0]);
                        hot = hot.max(t[0]);
                    }
                }
            }
        }
        (sum / (seasons.len().max(1) * CYCLE_DAY_STEPS) as f64, hot)
    }

    /// A cell's ice through a step (F240, F245): lit by `lit` W/m² at its `albedo`, the
    /// frost laying its frost on it where `frost` allows (F239), under air at `air` K.
    fn ice_step(
        &self,
        column: &Column,
        (t, mass, work, trial): (&mut Vec<f64>, &mut f64, &mut Work, &mut Trial),
        dt: f64,
        (lit, albedo): (f64, f64),
        frost: Option<&Frost>,
        air: f64,
    ) {
        let shed = self.shed(air);
        let own = (1.0 - albedo) * lit;
        match frost {
            Some(f) => {
                let gains = (own, (1.0 - f.albedo) * lit);
                let sheds = (shed, self.frost_shed(air, f));
                column.frosty_step(t, (work, trial), dt, gains, self.heat_flux, sheds, f, mass);
            }
            None => column.step(t, work, dt, own, self.heat_flux, shed),
        }
    }

    /// A fast spinner's year at a latitude and tilt.
    fn zonal_cycle(
        &self,
        ground: &Ground,
        lat: f64,
        eps: f64,
        steady: f64,
        albedo: f64,
        brought: f64,
    ) -> Year {
        let olr = |t: f64| self.olr_and_slope(t);
        let solar_days = (self.sky.spins - eps.cos()).abs().max(1e-3);
        let day = self.year / solar_days;
        let column = Column::new(ground, day, self.year);
        let mut work = Work::new(column.len());
        let mut t = vec![steady; column.len()];
        let dt = self.year / CYCLE_YEAR_STEPS as f64;
        let mut seasons: Vec<(Vec<f64>, f64, f64)> = Vec::new();
        let mut last = f64::NAN;
        for _ in 0..CYCLE_MAX_YEARS {
            seasons.clear();
            let mut sums = vec![0.0; t.len()];
            for m in 0..CYCLE_YEAR_STEPS {
                let anomaly = TAU * (m as f64 + 0.5) / CYCLE_YEAR_STEPS as f64;
                let (flux, decl) = self.sky.at_mean_anomaly(anomaly, eps);
                let gain = (1.0 - albedo) * flux * diurnal_share(lat, decl) + brought;
                column.step(&mut t, &mut work, dt, gain, self.heat_flux, olr);
                for (s, v) in sums.iter_mut().zip(&t) {
                    *s += v;
                }
                if m % (CYCLE_YEAR_STEPS / CYCLE_SEASONS) == 0 {
                    seasons.push((t.clone(), flux, decl));
                }
            }
            let surface = sums[0] / CYCLE_YEAR_STEPS as f64;
            if (surface - last).abs() < CYCLE_SETTLED_K {
                break;
            }
            last = surface;
            settle(&mut t, &sums, CYCLE_YEAR_STEPS);
        }
        // A day at each season, from the ground as the season left it.
        let (sp, cp) = lat.sin_cos();
        let dt = day / CYCLE_DAY_STEPS as f64;
        let (mut mean, mut hi, mut lo) = (0.0, f64::MIN, f64::MAX);
        let mut samples = Vec::new();
        for (state, flux, decl) in &seasons {
            let (sd, cd) = decl.sin_cos();
            let mut t = state.clone();
            let mut sum = 0.0;
            for d in 0..CYCLE_DAYS {
                for h in 0..CYCLE_DAY_STEPS {
                    let hour = TAU * (h as f64 + 0.5) / CYCLE_DAY_STEPS as f64 - PI;
                    let cz = (sp * sd + cp * cd * hour.cos()).max(0.0);
                    let gain = (1.0 - albedo) * flux * cz + brought;
                    column.step(&mut t, &mut work, dt, gain, self.heat_flux, olr);
                    if d + 1 == CYCLE_DAYS {
                        sum += t[0];
                        hi = hi.max(t[0]);
                        lo = lo.min(t[0]);
                        samples.push(t[0]);
                    }
                }
            }
            mean += sum / CYCLE_DAY_STEPS as f64;
        }
        let spread = (spread_of(&mut samples), albedo);
        Year::bare(mean / seasons.len().max(1) as f64, hi, lo, spread)
    }

    /// A resonant body's cycle at a place (latitude, longitude), started from its year's
    /// balance (temperature, albedo, transport).
    fn resonant_cycle(
        &self,
        ground: &Ground,
        orbits: u32,
        path: &[PathStep],
        (lat, lon): (f64, f64),
        (steady, albedo, brought): (f64, f64, f64),
    ) -> Year {
        let olr = |t: f64| self.olr_and_slope(t);
        let cycle = self.year * orbits as f64;
        let column = Column::new(ground, cycle, cycle);
        let mut work = Work::new(column.len());
        let mut t = vec![steady; column.len()];
        let dt = cycle / path.len() as f64;
        let e = self.sky.eccentricity;
        let scale = self.sky.flux / (1.0 - e * e).sqrt();
        let (sp, cp) = lat.sin_cos();
        let (mut mean, mut hi, mut lo) = (0.0, f64::MIN, f64::MAX);
        let mut last = f64::NAN;
        let mut samples = Vec::with_capacity(path.len());
        for _ in 0..CYCLE_MAX_YEARS {
            let mut sums = vec![0.0; t.len()];
            (hi, lo) = (f64::MIN, f64::MAX);
            samples.clear();
            for s in path {
                let cz = (sp * s.sin_decl + cp * s.cos_decl * (s.hour + lon).cos()).max(0.0);
                let gain = (1.0 - albedo) * scale * s.flux * cz + brought;
                column.step(&mut t, &mut work, dt, gain, self.heat_flux, olr);
                for (sum, v) in sums.iter_mut().zip(&t) {
                    *sum += v;
                }
                hi = hi.max(t[0]);
                lo = lo.min(t[0]);
                samples.push(t[0]);
            }
            mean = sums[0] / path.len() as f64;
            if (mean - last).abs() < CYCLE_SETTLED_K {
                break;
            }
            last = mean;
            settle(&mut t, &sums, path.len());
        }
        Year::bare(mean, hi, lo, (spread_of(&mut samples), albedo))
    }
}

/// The mean anomaly of a season, an orbit fraction from the northern spring equinox, on an
/// orbit whose perihelion falls at a solar longitude, rad.
fn snapshot_anomaly(season: f64, perihelion: f64, e: f64) -> f64 {
    let nu = (TAU * season - perihelion).rem_euclid(TAU);
    let big_e = 2.0 * (((1.0 - e) / (1.0 + e)).sqrt() * (nu / 2.0).tan()).atan();
    (big_e - e * big_e.sin()).rem_euclid(TAU)
}

/// The step of an orbit sampled in `steps` that holds a mean anomaly.
fn snapshot_step(anomaly: f64, steps: usize) -> usize {
    ((anomaly / TAU * steps as f64) as usize).min(steps.saturating_sub(1))
}

/// The phase of a cycle of `steps` over `orbits` its snapshot is taken at (F284), counted
/// from its equinox, `equinox` of the cycle after perihelion: the end of the step of its first
/// orbit that holds its season's mean `anomaly`, as each step's record is its end and the
/// windows' phases are their steps' ends; not the season's own, up to a step before it.
fn snapshot_phase(anomaly: f64, steps: usize, orbits: usize, equinox: f64) -> f64 {
    let k = snapshot_step(anomaly, steps / orbits.max(1));
    ((k as f64 + 1.0) / steps as f64 - equinox).rem_euclid(1.0)
}

/// Moves each layer of the ground below the surface so that its mean over the last cycle
/// is the surface's: the deep ground settles at the surface's mean, as it does once the
/// cycles repeat, without waiting the many cycles its depth takes.
fn settle(t: &mut [f64], sums: &[f64], steps: usize) {
    let surface = sums[0] / steps as f64;
    for (v, s) in t.iter_mut().zip(sums).skip(1) {
        *v += surface - s / steps as f64;
    }
}

/// What the downscaling reads (F237): how the greenhouse above a place thins with its height,
/// and how the air holds a slope to its surroundings.
#[derive(Clone, Copy, Debug)]
pub struct Downscale {
    /// The dry gases' optical depths at the datum, each with the power of its column it
    /// grows by, tau = a u^b: carbon dioxide, methane, nitrogen's collisions.
    pub gases: [(f64, f64); 3],
    /// The water a warm atmosphere can hold vapour from, metres of global layer.
    pub water: f64,
    /// The pressure's scale height, m.
    pub scale_height: f64,
    /// The lapse rate, K/m: no surface cools with height faster than the air convecting
    /// over it.
    pub lapse: f64,
    /// The heat the air exchanges with the ground per kelvin between them, W/m²/K.
    pub exchange: f64,
    /// The body's radius, m.
    pub radius: f64,
}

impl Downscale {
    pub fn from_graph(graph: &ParamGraph) -> Self {
        let empty = std::collections::BTreeMap::new();
        let composition = graph
            .get("atmosphere.composition")
            .as_map()
            .unwrap_or(&empty);
        let pressure = graph.f64("atmosphere.surface_pressure");
        if pressure <= 0.0 {
            return Self {
                gases: [(0.0, 1.0); 3],
                water: 0.0,
                scale_height: 0.0,
                lapse: 0.0,
                exchange: 0.0,
                radius: graph.f64("body.radius"),
            };
        }
        let columns = species_columns(composition, graph.f64("atmosphere.column_mass"));
        let tau = |species: &str, (a, b): (f64, f64)| {
            (a * columns.get(species).copied().unwrap_or(0.0).powf(b), b)
        };
        let water = EnergyBalance::from_graph(graph).water;
        let (ventilation, heat) = ventilation(graph);
        Self {
            gases: [
                tau("CO2", TAU_CO2),
                tau("CH4", TAU_CH4),
                tau("N2", TAU_N2_CIA),
            ],
            water,
            scale_height: graph.f64("atmosphere.scale_height"),
            lapse: graph.f64("climate.lapse_rate"),
            exchange: ventilation * heat,
            radius: graph.f64("body.radius"),
        }
    }

    /// The surface temperature at a place, from the climate's temperature there (`flat`) and
    /// at the place its slope faces (`facing`), at a height above the datum, K. The ground
    /// sheds its heat to space and to the air; the share it sheds by radiation follows its
    /// own sunlight and greenhouse, the rest the air over it. So its slope moves it toward
    /// the facing place's temperature by the radiation's share, and its height cools it by
    /// the air's share of the lapse rate and the radiation's share of the thinner greenhouse
    /// above it, the dry gases thinning with the pressure and the water vapour with its own
    /// scale height; never faster than the lapse rate. A height that is not finite is taken
    /// at the datum.
    pub fn at(&self, flat: f64, facing: f64, height: f64) -> f64 {
        let height = if height.is_finite() { height } else { 0.0 };
        let keep = self.keep(flat);
        let t = flat + keep * (facing - flat);
        if self.scale_height <= 0.0 || height == 0.0 {
            return t;
        }
        let ratio = (-height / self.scale_height).exp();
        let dry_datum: f64 = self.gases.iter().map(|g| g.0).sum();
        let dry_here: f64 = self.gases.iter().map(|&(tau, b)| tau * ratio.powf(b)).sum();
        let water_datum = water_optical_depth(flat, self.water);
        let water_here = if self.water > 0.0 {
            let scale = flat * flat / (WATER_L_OVER_RV * self.lapse.max(1e-9));
            let column = water_column(flat, self.water) * (-height / scale).exp();
            TAU_H2O.0 * (column / TAU_H2O.1).ln_1p()
        } else {
            0.0
        };
        let thinning =
            (1.0 + 0.75 * (dry_here + water_here)) / (1.0 + 0.75 * (dry_datum + water_datum));
        let change = keep * t * (thinning.powf(0.25) - 1.0) - (1.0 - keep) * self.lapse * height;
        let cap = self.lapse * height.abs();
        t + change.clamp(-cap, cap)
    }

    /// The year's mean surface temperature on the surface's own grid from the climate's: each
    /// place at its height above the sea's surface (or the mean radius without a sea), its
    /// slope from its neighbours' heights, the sea's surface flat.
    /// The share of the sunlight a slope faces that its temperature keeps against its
    /// surroundings', at a temperature: what the ground radiates per kelvin against what the
    /// air exchanges with it.
    fn keep(&self, t: f64) -> f64 {
        let radiative = 4.0 * SIGMA_SB * t.powi(3);
        radiative / (radiative + self.exchange).max(1e-30)
    }

    /// The frost on the ground at the snapshot's season brought to a grid (F239), kg/m²: a
    /// place takes the frost of the place its slope faces, as its temperature does, held
    /// toward its surroundings' by the air as its temperature is, at the climate's `mean`
    /// temperatures; a sea's surface is flat. The frost does not yet see a height's thinner
    /// air and its lower frost point.
    pub fn frost(
        &self,
        frost: &Field<f32>,
        mean: &Field<f32>,
        height: &Field<f32>,
        sea_level: f64,
    ) -> Field<f32> {
        let src = frost.grid();
        let lying = PaddedField::new(&src, frost.as_slice());
        let warmth = PaddedField::new(&src, mean.as_slice());
        let sea = sea_level.is_finite();
        Field::from_fn(height.grid(), |c, d| {
            let flat = sample_bilinear(&src, &lying, d) as f64;
            if sea && height.get(c) as f64 <= sea_level {
                return flat.max(0.0) as f32;
            }
            let normal = surface_normal(height, c, d, self.radius);
            let facing = sample_bilinear(&src, &lying, normal) as f64;
            let facing = if facing.is_finite() { facing } else { flat };
            let keep = self.keep(sample_bilinear(&src, &warmth, d) as f64);
            (flat + keep * (facing - flat)).max(0.0) as f32
        })
    }

    /// The ices' windows on a grid (F246): a sea's ice's where the grid's liquid lies, its
    /// peak's thickness, m, at `densities.1`; a land's frost's where `frost` reads its frost,
    /// the place its slope faces where its own sunlight rules it, its peak's thickness brought
    /// there as the frost is, at `densities.0`. Each its freeze and melt phases, its peak's
    /// thickness and its peak's phase, a field apiece, and which ice lies (F253): 1 the
    /// land's frost, 2 the sea's ice, 3 the land's snow, its year run on the ground's own
    /// height over the sea at `sea_level` (F271, F273), 0 none.
    pub fn seasons(
        &self,
        cryosphere: &Cryosphere,
        mean: &Field<f32>,
        (height, liquid, sea_level): (&Field<f32>, &Field<f32>, f64),
        densities: (f64, f64),
    ) -> [Field<f32>; 5] {
        let src = mean.grid();
        let warmth = PaddedField::new(&src, mean.as_slice());
        let grid = height.grid();
        let windows: Vec<(Window, f32)> = (0..grid.cell_count())
            .into_par_iter()
            .map(|k| {
                let c = grid.cell(k);
                let d = grid.center(c);
                if liquid.get(c) > 0.0 {
                    let w = cryosphere.at(d, true);
                    let kind = if w.peak > 0.0 { 2.0 } else { 0.0 };
                    return (
                        Window {
                            peak: w.peak / densities.1,
                            ..w
                        },
                        kind,
                    );
                }
                let flat = cryosphere.at(d, false);
                let facing = cryosphere.at(surface_normal(height, c, d, self.radius), false);
                let keep = self.keep(sample_bilinear(&src, &warmth, d) as f64);
                let peak = (flat.peak + keep * (facing.peak - flat.peak)).max(0.0);
                let read = if (keep > 0.5 && facing.peak > 0.0) || flat.peak <= 0.0 {
                    facing
                } else {
                    flat
                };
                // The land's snow where it covers more of the ground than the frost (F271).
                let frosted = peak > 0.0 && read.peak > 0.0;
                let h = height.get(c) as f64;
                let above = if sea_level.is_finite() {
                    (h - sea_level).max(0.0)
                } else {
                    h
                };
                let snow = cryosphere.snow_over(d, above, self.lapse);
                let snowier = !frosted || snow.peak / SNOW_COVER_MASS > peak / FROST_COVER_MASS;
                if snow.peak > 0.0 && snowier {
                    let peak = snow.peak / SNOW_DENSITY;
                    return (Window { peak, ..snow }, 3.0);
                }
                if !frosted {
                    return (Window::default(), 0.0);
                }
                (
                    Window {
                        peak: peak / densities.0,
                        ..read
                    },
                    1.0,
                )
            })
            .collect();
        let field = |pick: fn(&Window) -> f64| {
            Field::from_fn(grid, |c, _| pick(&windows[grid.index(c)].0) as f32)
        };
        [
            field(|w| w.freeze),
            field(|w| w.melt),
            field(|w| w.peak),
            field(|w| w.at),
            Field::from_fn(grid, |c, _| windows[grid.index(c)].1),
        ]
    }

    pub fn apply(&self, climate: &Field<f32>, height: &Field<f32>, sea_level: f64) -> Field<f32> {
        let grid = height.grid();
        let src = climate.grid();
        let padded = PaddedField::new(&src, climate.as_slice());
        if std::env::var_os("EXOGEN_DEBUG").is_some() {
            let bad = height.as_slice().iter().filter(|h| !h.is_finite()).count();
            if bad > 0 {
                eprintln!("climate: {bad} heights not finite, taken at the datum");
            }
        }
        let sea = sea_level.is_finite();
        Field::from_fn(grid, |c, d| {
            let flat = sample_bilinear(&src, &padded, d) as f64;
            let h = height.get(c) as f64;
            let under_sea = sea && h <= sea_level;
            let facing = if under_sea {
                flat
            } else {
                sample_bilinear(&src, &padded, surface_normal(height, c, d, self.radius)) as f64
            };
            let above = if sea { (h - sea_level).max(0.0) } else { h };
            self.at(flat, facing, above) as f32
        })
    }

    /// The poles' shadows (F243): ground whose horizon stands, in every direction, above the
    /// highest the star's disc climbs there through the year never sees it (`Shade`). It
    /// stands at what the terrain it sees and the body's own heat hold it at: the terrain's
    /// share of its sky, the mean of sin² of its horizon's elevation, at what that terrain
    /// radiates through the year, `emitted` W/m² on the climate's grid, its own heat with
    /// it; and the air holds it toward its surroundings as it does a slope. The horizon is
    /// traced in SHADOW_AZIMUTHS directions out to SHADOW_STEPS cells, over the body's
    /// curve, where the star never climbs past SHADOW_SUN_LIMIT_DEG. Returns the share of
    /// the body in shadow.
    pub fn shade(
        &self,
        temperature: &mut Field<f32>,
        height: &Field<f32>,
        emitted: &Field<f32>,
        shade: &Shade,
    ) -> f64 {
        let grid = height.grid();
        let src = emitted.grid();
        let padded = PaddedField::new(&src, emitted.as_slice());
        let step = (4.0 * PI / grid.cell_count() as f64).sqrt();
        let limit = SHADOW_SUN_LIMIT_DEG.to_radians();
        let level = |d: [f64; 3]| {
            let h = height.get(grid.locate(d).0) as f64;
            if h.is_finite() { h } else { 0.0 }
        };
        // The horizon in a direction: the highest the ground stands over a place, over the
        // body's curve; traced until it passes `enough`.
        let horizon =
            |d: [f64; 3], h0: f64, (east, north): ([f64; 3], [f64; 3]), phi: f64, enough: f64| {
                let mut most = f64::MIN;
                for s in 1..=SHADOW_STEPS {
                    let theta = s as f64 * step;
                    let p = crate::geom::offset(d, east, north, theta, phi);
                    let distance = theta * self.radius;
                    let drop = distance * distance / (2.0 * self.radius);
                    most = most.max(((level(p) - h0 - drop) / distance).atan());
                    if most > enough {
                        break;
                    }
                }
                most
            };
        let shadowed: Vec<(usize, f32)> = (0..grid.cell_count())
            .into_par_iter()
            .filter_map(|k| {
                let c = grid.cell(k);
                let d = grid.center(c);
                let (_, lat) = lon_lat_from_direction(d);
                if shade.highest(lat, 0.0) > limit {
                    return None;
                }
                let h0 = level(d);
                let basis = crate::geom::local_basis(d);
                let equator = if lat > 0.0 { -PI / 2.0 } else { PI / 2.0 };
                let alphas = (0..SHADOW_AZIMUTHS).map(|j| TAU * j as f64 / SHADOW_AZIMUTHS as f64);
                // Lit if the star clears the horizon in any direction, the equator's first.
                for alpha in alphas.clone() {
                    let sun = shade.highest(lat, alpha);
                    if horizon(d, h0, basis, equator + alpha, sun) <= sun {
                        return None;
                    }
                }
                let terrain = alphas
                    .map(|alpha| {
                        horizon(d, h0, basis, equator + alpha, f64::MAX)
                            .max(0.0)
                            .sin()
                            .powi(2)
                    })
                    .sum::<f64>()
                    / SHADOW_AZIMUTHS as f64;
                let walls = sample_bilinear(&src, &padded, d) as f64 + shade.heat;
                let own = ((terrain * walls + shade.heat) / SIGMA_SB).powf(0.25);
                let flat = temperature.get(c) as f64;
                Some((k, (flat + self.keep(flat) * (own - flat)) as f32))
            })
            .collect();
        let mut area = 0.0;
        for &(k, t) in &shadowed {
            temperature.as_mut_slice()[k] = t;
            area += grid.cell_area(grid.cell(k));
        }
        area / (4.0 * PI)
    }
}

/// What the poles' shadows read (F243): the body's largest tilt to its star's orbit through
/// its precession, which sets the highest the star climbs at a pole, rad; the star's
/// angular radius, rad; and the body's own heat flowing out, W/m².
#[derive(Clone, Copy, Debug)]
pub struct Shade {
    pub tilt: f64,
    pub sun: f64,
    pub heat: f64,
}

impl Shade {
    pub fn from_graph(graph: &ParamGraph, balance: &EnergyBalance) -> Self {
        let tilt = balance
            .sky
            .obliquities()
            .iter()
            .fold(0.0f64, |a, &e| a.max(e.min(PI - e)));
        let distance = graph.f64("orbit.star_distance").max(1.0);
        Self {
            tilt,
            sun: (graph.f64("star.radius") * R_SUN / distance).atan(),
            heat: balance.heat_flux,
        }
    }

    /// The highest the star's disc climbs over the year at a latitude, rad, seen in a
    /// direction `alpha` from the equator's: near a pole, the star's azimuth is its hour
    /// angle, and at the solstice it stands at sin(phi) sin(tilt) + cos(phi) cos(tilt)
    /// cos(alpha).
    fn highest(&self, lat: f64, alpha: f64) -> f64 {
        let (sp, cp) = lat.abs().sin_cos();
        let (st, ct) = self.tilt.sin_cos();
        (sp * st + cp * ct * alpha.cos()).clamp(-1.0, 1.0).asin() + self.sun
    }
}

/// The most an ice of the frost's own gas would grow by anywhere, kg/m² a year, and the
/// share of the body's area on which it would grow (F240).
pub fn growth_reach(growth: &Field<f32>) -> (f64, f64) {
    let grid = growth.grid();
    let (mut most, mut grown, mut whole) = (f64::MIN, 0.0, 0.0);
    for (k, &g) in growth.as_slice().iter().enumerate() {
        let a = grid.cell_area(grid.cell(k));
        most = most.max(g as f64);
        whole += a;
        if g > 0.0 {
            grown += a;
        }
    }
    (most, grown / whole.max(1e-30))
}

/// A value of a fast spinner's traps at a latitude and a share of the datum's pressure
/// (F240): `values` its CYCLE_LATITUDES latitudes' at each of `levels`, read linearly between
/// them, in the logarithm where `log`, and in the logarithm of the share.
pub(crate) fn read_table(values: &[f64], levels: &[f64], lat: f64, thin: f64, log: bool) -> f64 {
    let k = levels.len();
    let n = values.len() / k;
    let x = ((lat + PI / 2.0) / PI * (n - 1) as f64).clamp(0.0, (n - 1) as f64);
    let i = (x.floor() as usize).min(n - 2);
    let wx = x - i as f64;
    let (j, wy) = if k == 1 {
        (0, 0.0)
    } else {
        let span = (levels[k - 1] / levels[0]).ln();
        let y = ((thin / levels[0]).ln() / span * (k - 1) as f64).clamp(0.0, (k - 1) as f64);
        let j = (y.floor() as usize).min(k - 2);
        (j, y - j as f64)
    };
    let at = |i: usize, j: usize| {
        let v = values[i * k + j];
        if log { v.max(1e-300).ln() } else { v }
    };
    let row = |i: usize| match k {
        1 => at(i, 0),
        _ => at(i, j) * (1.0 - wy) + at(i, j + 1) * wy,
    };
    let v = row(i) * (1.0 - wx) + row(i + 1) * wx;
    if log { v.exp() } else { v }
}

/// A field's area-weighted mean over each cell of a coarser grid (F240), its non-finite
/// values left out.
pub fn coarsen(fine: &Field<f32>, coarse: CubeSphere) -> Field<f32> {
    let grid = fine.grid();
    let n = coarse.cell_count();
    let (mut sum, mut area) = (vec![0.0f64; n], vec![0.0f64; n]);
    for (k, &v) in fine.as_slice().iter().enumerate() {
        if !v.is_finite() {
            continue;
        }
        let c = grid.cell(k);
        let (to, _, _) = coarse.locate(grid.center(c));
        let j = coarse.index(to);
        let a = grid.cell_area(c);
        sum[j] += v as f64 * a;
        area[j] += a;
    }
    Field::from_fn(coarse, |c, _| {
        let j = coarse.index(c);
        if area[j] > 0.0 {
            (sum[j] / area[j]) as f32
        } else {
            0.0
        }
    })
}

/// The ground's colour under its frost (F239): the frost hides it as 1 - exp(-m / M), m the
/// frost and M FROST_COVER_MASS.
pub fn frosted_albedo(
    albedo: &Field<[f32; 3]>,
    frost: &Field<f32>,
    colour: [f32; 3],
) -> Field<[f32; 3]> {
    Field::from_fn(albedo.grid(), |c, _| {
        let cover = frost_cover(frost.get(c));
        let a = albedo.get(c);
        std::array::from_fn(|i| (a[i] as f64 * (1.0 - cover) + colour[i] as f64 * cover) as f32)
    })
}

/// The share of the land snow of `m` kg/m² of water covers (F271).
pub fn snow_cover(m: f64) -> f64 {
    1.0 - (-m.max(0.0) / SNOW_COVER_MASS).exp()
}

/// The share of the ground's colour a frost of `m` kg/m² hides (F239).
pub fn frost_cover(m: f32) -> f64 {
    1.0 - (-(m.max(0.0) as f64) / FROST_COVER_MASS).exp()
}

/// The unit normal of the ground at a cell, from its four neighbours' heights by least
/// squares in the plane of the cell.
fn surface_normal(height: &Field<f32>, c: crate::grid::Cell, d: [f64; 3], radius: f64) -> [f64; 3] {
    let grid = height.grid();
    let (east, north) = crate::geom::local_basis(d);
    let h0 = height.get(c) as f64;
    let (mut sxx, mut sxy, mut syy, mut sxh, mut syh) = (0.0, 0.0, 0.0, 0.0, 0.0);
    for (di, dj) in [(1, 0), (-1, 0), (0, 1), (0, -1)] {
        let nb = grid.neighbor(c, di, dj);
        let p = grid.center(nb);
        let dp = [
            (p[0] - d[0]) * radius,
            (p[1] - d[1]) * radius,
            (p[2] - d[2]) * radius,
        ];
        let (x, y) = (crate::geom::dot(dp, east), crate::geom::dot(dp, north));
        let dh = height.get(nb) as f64 - h0;
        sxx += x * x;
        sxy += x * y;
        syy += y * y;
        sxh += x * dh;
        syh += y * dh;
    }
    let det = sxx * syy - sxy * sxy;
    if det.abs() < 1e-30 {
        return d;
    }
    let gx = (syy * sxh - sxy * syh) / det;
    let gy = (sxx * syh - sxy * sxh) / det;
    crate::geom::normalize([
        d[0] - gx * east[0] - gy * north[0],
        d[1] - gx * east[1] - gy * north[1],
        d[2] - gx * east[2] - gy * north[2],
    ])
}

impl Cycle {
    /// What its air exchanges with its ground reads of the year (F269): at each place its
    /// wind's speed through the weather about it, sqrt(|V|² + U²), U the body's mean wind, as
    /// the moisture reads it, over the stated SURFACE_WIND_SPEED, and, with `moisture`, the
    /// air's mean saturation over it and the share of its land kept wet; none without the
    /// year's wind.
    pub fn coupling(&self, moisture: Option<&crate::moisture::Moisture>) -> Option<Coupling> {
        let wind = self.wind.as_ref()?;
        let n = wind.layout.0 * wind.layout.1;
        let steps = wind.steps.max(1);
        let contacts = (0..n)
            .map(|p| {
                let breeze = (0..steps)
                    .map(|k| {
                        let (u, v) = (wind.east[k * n + p] as f64, wind.north[k * n + p] as f64);
                        u.hypot(v).hypot(wind.speed[k])
                    })
                    .sum::<f64>()
                    / steps as f64
                    / SURFACE_WIND_SPEED;
                let (humidity, wet) =
                    moisture
                        .filter(|m| m.saturation.len() == n)
                        .map_or((1.0, 0.0), |m| {
                            (
                                m.saturation[p].clamp(0.0, 1.0),
                                m.wetness[p].clamp(0.0, 1.0),
                            )
                        });
                Contact {
                    breeze,
                    humidity,
                    wet,
                }
            })
            .collect();
        // The snow that falls on each place's land through the year (F271).
        let series = |share: bool| match moisture.filter(|m| m.falling.len() == n * m.steps) {
            Some(m) => (0..n * m.steps)
                .map(|i| {
                    let (p, k) = (i / m.steps, i % m.steps);
                    let falls = m.falling[k * n + p];
                    if share {
                        falls * m.snow[k * n + p]
                    } else {
                        falls
                    }
                })
                .collect(),
            None => Vec::new(),
        };
        let (snowfall, falling) = (series(true), series(false));
        Some(Coupling {
            layout: wind.layout,
            contacts,
            snowfall,
            falling,
            steps: moisture.map_or(0, |m| m.steps),
            currents: self.currents.clone(),
            wind: wind
                .east
                .iter()
                .zip(&wind.north)
                .map(|(&u, &v)| [u, v])
                .collect(),
            wind_steps: wind.steps,
        })
    }

    /// The moisture through the cycle (F268) of the air `air`, each place's ground standing
    /// `heights` m over the sea: none without the cycle's wind.
    pub fn moisture(
        &self,
        air: &crate::moisture::Air,
        heights: &[f64],
    ) -> Option<crate::moisture::Moisture> {
        let wind = self.wind.as_ref()?;
        let (n_lat, n_lon) = self.layout;
        let n = n_lat * n_lon;
        let steps = wind.steps;
        let year = self.clock?.cycle_s;
        if self.sites.len() != n
            || self.air.is_empty()
            || self.surfaces[0].len() != n * steps
            || self.masses[1].len() != n * steps
        {
            return None;
        }
        let tilts = self.air.len() as f64;
        let mean: Vec<f64> = (0..steps * n)
            .map(|i| self.air.iter().map(|a| a[i]).sum::<f64>() / tilts)
            .collect();
        let at = |field: &[f32], k: usize| -> Vec<f64> {
            (0..n).map(|p| field[p * steps + k] as f64).collect()
        };
        let weather = |k: usize| crate::moisture::Weather {
            air: &mean[k * n..(k + 1) * n],
            east: &wind.east[k * n..(k + 1) * n],
            north: &wind.north[k * n..(k + 1) * n],
            land: at(&self.surfaces[0], k),
            sea: at(&self.surfaces[1], k),
            ice: at(&self.masses[1], k),
            snow: (self.masses[3].len() == n * steps)
                .then(|| (at(&self.masses[2], k), at(&self.masses[3], k))),
            gust: wind.speed[k],
        };
        let layout = Layout::of_sites(&self.sites, (n_lat, n_lon));
        let solve = |links: &[Vec<(usize, f64)>], diag: &[f64], rhs: &[f64], x: &mut [f64]| {
            transport_solve(&layout, links, (diag, rhs), x, true);
        };
        crate::moisture::rain(
            air,
            &self.sites,
            (n_lat, n_lon),
            heights,
            (year, steps),
            &weather,
            &solve,
        )
    }

    /// The surface's temperatures through its cycle in a direction (F256): the quantiles of
    /// the places about it, read between them, and the albedo their ground stood at.
    pub fn spread_at(&self, d: [f64; 3]) -> ([f64; CYCLE_SPREAD], f64) {
        let (lon, lat) = lon_lat_from_direction(d);
        let corners = reader(self.layout.0, self.layout.1, lat, lon);
        let (mut out, mut albedo) = ([0.0; CYCLE_SPREAD], 0.0);
        for &(p, w) in corners.iter().filter(|c| c.1 > 0.0) {
            for (o, &q) in out.iter_mut().zip(&self.spread[p]) {
                *o += w * q as f64;
            }
            albedo += w * self.spread_albedo[p] as f64;
        }
        (out, albedo)
    }

    /// The frost point through the year where the air keeps `thin` of the datum's pressure,
    /// K (F258): its mean, its lowest and its highest through the year's steps at each tilt,
    /// where the atmosphere lays frost. The frost's own gas's caps stand at it.
    pub fn frost_point_year(&self, thin: f64) -> Option<(f64, f64, f64)> {
        let f = self.frost?;
        let (mut sum, mut lo, mut hi, mut n) = (0.0, f64::MAX, f64::MIN, 0usize);
        for &p in self.partials.iter().flatten() {
            let t = f.point_at(p * thin);
            sum += t;
            (lo, hi) = (lo.min(t), hi.max(t));
            n += 1;
        }
        (n > 0).then(|| (sum / n as f64, lo, hi))
    }

    /// The coldest and the warmest the surface stands anywhere through the year, K.
    pub fn extremes(&self) -> (f64, f64) {
        let coldest = self
            .min
            .as_slice()
            .iter()
            .fold(f64::MAX, |a, &v| a.min(v as f64));
        let warmest = self
            .max
            .as_slice()
            .iter()
            .fold(f64::MIN, |a, &v| a.max(v as f64));
        (coldest, warmest)
    }
}

/// A column of ground: each layer's heat capacity per square metre and the conductance from
/// each layer to the one below.
#[derive(Clone, Debug)]
struct Column {
    capacity: Vec<f64>,
    conductance: Vec<f64>,
}

/// The tridiagonal solve's working rows.
#[derive(Clone)]
struct Work {
    a: Vec<f64>,
    b: Vec<f64>,
    c: Vec<f64>,
    d: Vec<f64>,
    /// The column as a step began, kept to step again where its frost runs out (F275).
    start: Vec<f64>,
}

impl Work {
    fn new(n: usize) -> Self {
        Self {
            a: vec![0.0; n],
            b: vec![0.0; n],
            c: vec![0.0; n],
            d: vec![0.0; n],
            start: Vec::with_capacity(n),
        }
    }
}

/// A step's outcomes before its frost point is settled (F239): bare ground's step left to
/// itself, and the column below a surface held at a point, each layer linear in the point.
#[derive(Clone)]
struct Trial {
    /// Bare ground's free step.
    free: Vec<f64>,
    /// The held column: each layer at base + point × lift.
    base: Vec<f64>,
    lift: Vec<f64>,
    /// The surface as the step found it, and what it gains, W/m².
    start: f64,
    gain: f64,
    /// Bare ground the step leaves above every point it could settle at.
    clear: bool,
}

impl Trial {
    fn new(n: usize) -> Self {
        Self {
            free: vec![0.0; n],
            base: vec![0.0; n],
            lift: vec![0.0; n],
            start: 0.0,
            gain: 0.0,
            clear: true,
        }
    }
}

/// What the air over a place exchanges with its ground reads of the year's wind and moisture
/// (F269): its wind's speed through the weather about it over the stated SURFACE_WIND_SPEED
/// the exchange was given at, the air's saturation over it, and the share of its land the
/// bucket keeps wet, whose latent heat ties it to its air as the sea's ties the sea. By
/// default the stated wind, saturated air and dry land, the cycle as it ran before the
/// moisture.
#[derive(Clone, Copy, Debug)]
pub struct Contact {
    pub breeze: f64,
    pub humidity: f64,
    pub wet: f64,
}

impl Default for Contact {
    fn default() -> Self {
        Self {
            breeze: 1.0,
            humidity: 1.0,
            wet: 0.0,
        }
    }
}

/// Each place's contact from a year's wind and moisture (F269), on the layout they ran at, and
/// the snow that falls on its land through the year's `steps`, kg/m²/s, place after place,
/// step after step within each (F271).
#[derive(Clone, Debug)]
pub struct Coupling {
    pub layout: (usize, usize),
    pub contacts: Vec<Contact>,
    pub snowfall: Vec<f32>,
    /// The precipitation on each place through the year's steps, kg/m²/s, laid out as the
    /// snow is (F273).
    pub falling: Vec<f32>,
    pub steps: usize,
    /// The seas' currents the year's wind drives (F272).
    pub currents: Option<crate::currents::Currents>,
    /// The wind at each place through the year's `wind_steps`, east and north, m/s, step
    /// after step, place after place within each, that drifts the sea ice (F276).
    pub wind: Vec<[f32; 2]>,
    pub wind_steps: usize,
}

impl Coupling {
    /// The year's wind on a layout of `n_lat` latitudes and `n_lon` longitudes (F276), as
    /// `series_at` reads the snow: each place's own on the coupling's layout, its row's mean on
    /// a layout by latitude alone, none elsewhere; step after step, place after place within
    /// each.
    fn wind_on(&self, (n_lat, n_lon): (usize, usize)) -> Vec<[f32; 2]> {
        let (m_lat, m_lon) = self.layout;
        let (m, steps) = (m_lat * m_lon, self.wind_steps);
        if steps == 0 || self.wind.len() != m * steps {
            return Vec::new();
        }
        if (m_lat, m_lon) == (n_lat, n_lon) {
            return self.wind.clone();
        }
        if n_lon == 1 && m_lat == n_lat && m_lon > 0 {
            return (0..steps * n_lat)
                .map(|i| {
                    let (k, row) = (i / n_lat, i % n_lat);
                    let mut mean = [0.0f32; 2];
                    for v in &self.wind[k * m + row * m_lon..k * m + (row + 1) * m_lon] {
                        mean[0] += v[0] / m_lon as f32;
                        mean[1] += v[1] / m_lon as f32;
                    }
                    mean
                })
                .collect();
        }
        Vec::new()
    }

    /// The contact at row `i`, place `k`, of a layout of `n_lat` latitudes and `n_lon`
    /// longitudes: the place's own on the coupling's layout, its row's mean on a layout by
    /// latitude alone, the default elsewhere.
    fn at(&self, (n_lat, n_lon): (usize, usize), k: usize) -> Contact {
        let (m_lat, m_lon) = self.layout;
        if (m_lat, m_lon) == (n_lat, n_lon) {
            return self.contacts[k];
        }
        if n_lon == 1 && m_lat == n_lat && m_lon > 0 {
            let row = &self.contacts[k * m_lon..(k + 1) * m_lon];
            let mean = |f: fn(&Contact) -> f64| row.iter().map(f).sum::<f64>() / m_lon as f64;
            return Contact {
                breeze: mean(|c| c.breeze),
                humidity: mean(|c| c.humidity),
                wet: mean(|c| c.wet),
            };
        }
        Contact::default()
    }

    /// The snow that falls on place `k`'s land through a year of `steps` on a layout of
    /// `n_lat` latitudes and `n_lon` longitudes (F271), as `at` reads its contact: the place's
    /// own on the coupling's layout, its row's mean on a layout by latitude alone; none where
    /// the year runs in other steps.
    fn snowfall_at(&self, layout: (usize, usize), steps: usize, k: usize) -> Vec<f32> {
        self.series_at(&self.snowfall, layout, steps, k)
    }

    /// The precipitation on place `k` through a year of `steps`, read as the snow is (F273).
    fn falling_at(&self, layout: (usize, usize), steps: usize, k: usize) -> Vec<f32> {
        self.series_at(&self.falling, layout, steps, k)
    }

    /// A series of `values` at each place through the year's steps, place `k`'s on a layout
    /// of `n_lat` latitudes and `n_lon` longitudes, as `snowfall_at` reads the snow.
    fn series_at(
        &self,
        values: &[f32],
        (n_lat, n_lon): (usize, usize),
        steps: usize,
        k: usize,
    ) -> Vec<f32> {
        let (m_lat, m_lon) = self.layout;
        if steps != self.steps || values.len() != m_lat * m_lon * steps {
            return Vec::new();
        }
        let series = |p: usize| &values[p * steps..(p + 1) * steps];
        if (m_lat, m_lon) == (n_lat, n_lon) {
            return series(k).to_vec();
        }
        if n_lon == 1 && m_lat == n_lat && m_lon > 0 {
            let mut mean = vec![0.0f32; steps];
            for j in 0..m_lon {
                for (m, &v) in mean.iter_mut().zip(series(k * m_lon + j)) {
                    *m += v / m_lon as f32;
                }
            }
            return mean;
        }
        Vec::new()
    }
}

/// A place of the cycle's layout (F239, F245): its latitude and longitude, rad, its share of
/// the sphere's area, its year's mean temperature, its ground's albedo without ice and snow,
/// the share of it the sea covers, and the share the frost's own gas's caps cover with the
/// ice they hold there, kg/m² (F258).
#[derive(Clone, Copy, Debug)]
struct Place {
    lat: f64,
    lon: f64,
    weight: f64,
    steady: f64,
    albedo: f64,
    snow: f64,
    sea: f64,
    cap: f64,
    cap_mass: f64,
    /// What its air exchanges with its ground reads of the year's wind and moisture (F269).
    contact: Contact,
}

impl Place {
    /// The share of the place its land covers bare of the caps.
    fn land(&self) -> f64 {
        (1.0 - self.sea - self.cap).max(0.0)
    }
}

/// The places the cycle runs at under an atmosphere (F245): each place's links to its
/// neighbours, the transport's conductance to each per unit of the place's area at unit
/// transport, W/m²/K; and the layout's longitudes, a fast spinner's one.
struct Layout {
    places: Vec<Place>,
    links: Vec<Vec<(usize, f64)>>,
    n_lon: usize,
    /// The snow that falls on each place's land through the year (F271), kg/m²/s, step after
    /// step, where water is the working liquid and the year's snow is known.
    snowfall: Vec<Vec<f32>>,
    /// The precipitation on each place through the year, kg/m²/s, where its snow is (F273).
    falling: Vec<Vec<f32>>,
    /// The seas' currents, where they ran on these places (F272).
    currents: Option<crate::currents::Currents>,
    /// The year's wind at each place, east and north, m/s, step after step, place after place
    /// within each, where it is known on these places (F276).
    wind: Vec<[f32; 2]>,
    /// A resonant body's longitudes' Fourier modes: the cos and sin of each mode's phase at
    /// each longitude, mode after mode, and longitude after longitude.
    modes: [Vec<f64>; 4],
}

/// The links of a layout of `n_lat` latitudes and `n_lon` longitudes (F245): each place's
/// neighbours and its conductance to each per unit of its area at unit transport, on the unit
/// sphere the length of the edge they share over the distance between them, over the place's
/// area.
pub(crate) fn place_links((n_lat, n_lon): (usize, usize)) -> Vec<Vec<(usize, f64)>> {
    let (dphi, dlon) = (PI / (n_lat - 1) as f64, TAU / n_lon as f64);
    let edges = |i: usize| {
        let lat = -PI / 2.0 + dphi * i as f64;
        (
            (lat - dphi / 2.0).max(-PI / 2.0),
            (lat + dphi / 2.0).min(PI / 2.0),
        )
    };
    (0..n_lat * n_lon)
        .map(|k| {
            let (i, j) = (k / n_lon, k % n_lon);
            let (lo, hi) = edges(i);
            let band = hi.sin() - lo.sin();
            let mut out = Vec::with_capacity(4);
            if i > 0 {
                out.push((k - n_lon, lo.cos() / (dphi * band)));
            }
            if i + 1 < n_lat {
                out.push((k + n_lon, hi.cos() / (dphi * band)));
            }
            if n_lon > 1 {
                let across = (hi - lo) / (dlon * (0.5 * (lo + hi)).cos() * dlon * band);
                let row = i * n_lon;
                out.push((row + (j + n_lon - 1) % n_lon, across));
                out.push((row + (j + 1) % n_lon, across));
            }
            out
        })
        .collect()
}

/// A resonant body's longitudes' Fourier modes over `n_lon` longitudes (F245): the cos and sin
/// of each mode's phase at each longitude, mode after mode, and longitude after longitude.
fn place_modes(n_lon: usize) -> [Vec<f64>; 4] {
    let count = if n_lon > 1 { n_lon / 2 + 1 } else { 0 };
    let phase = |m: usize, j: usize| TAU * (m * j % n_lon) as f64 / n_lon as f64;
    let by_mode = |f: fn(f64) -> f64| -> Vec<f64> {
        (0..count * n_lon)
            .map(|k| f(phase(k / n_lon, k % n_lon)))
            .collect()
    };
    let by_longitude = |f: fn(f64) -> f64| -> Vec<f64> {
        (0..count * n_lon)
            .map(|k| f(phase(k % count, k / count)))
            .collect()
    };
    [
        by_mode(f64::cos),
        by_mode(f64::sin),
        by_longitude(f64::cos),
        by_longitude(f64::sin),
    ]
}

impl Layout {
    /// A layout of the cycle's places as the wind and the moisture keep them (F268), for the
    /// transport's solver: their latitudes, longitudes, areas and shares of sea and caps.
    fn of_sites(sites: &[crate::wind::Site], (n_lat, n_lon): (usize, usize)) -> Self {
        Layout {
            places: sites
                .iter()
                .map(|s| Place {
                    lat: s.lat,
                    lon: s.lon,
                    weight: s.weight,
                    steady: 0.0,
                    albedo: 0.0,
                    snow: 0.0,
                    sea: s.sea,
                    cap: s.ice,
                    cap_mass: 0.0,
                    contact: Contact::default(),
                })
                .collect(),
            links: place_links((n_lat, n_lon)),
            n_lon,
            snowfall: Vec::new(),
            falling: Vec::new(),
            currents: None,
            wind: Vec::new(),
            modes: place_modes(n_lon),
        }
    }

    /// The transport through the air at `transport` W/m²/K (F245): each place's conductance to
    /// each of its neighbours, less the share `water` a sea's water carries across the sea's
    /// share of their boundary.
    fn air_links(&self, transport: f64, water: f64) -> Vec<Vec<(usize, f64)>> {
        self.links
            .iter()
            .enumerate()
            .map(|(p, links)| {
                links
                    .iter()
                    .map(|&(q, g)| {
                        let wet = self.places[p].sea.min(self.places[q].sea);
                        (q, transport * g * (1.0 - water * wet))
                    })
                    .collect()
            })
            .collect()
    }
}

/// A step of a place's ground recorded through the lockstep (F239, F245): the step, the
/// ground, its frost, the frost point, the air over it and the land's albedo.
struct Record {
    k: usize,
    t: Vec<f64>,
    mass: f64,
    point: f64,
    air: f64,
    albedo: f64,
}

/// A place's ground, its frost and its sea through the lockstep (F239, F245), and its last
/// cycle.
struct Run {
    t: Vec<f64>,
    mass: f64,
    work: Work,
    trial: Trial,
    sea: Option<Sea>,
    /// The land's albedo at the step.
    albedo: f64,
    /// Each layer's sum through the cycle, the surface's warmest and coldest, the steps frost
    /// lay and the frost at the snapshot's step.
    sums: Vec<f64>,
    hi: f64,
    lo: f64,
    frosted: usize,
    now: f64,
    /// The frost the cycle began with, its sum through the cycle and its least.
    frost_start: f64,
    frost_sum: f64,
    frost_least: f64,
    /// The sea's surface's sum through the cycle, its coldest and warmest, the steps ice lay
    /// on it and its ice at the snapshot's step.
    sea_sum: f64,
    sea_lo: f64,
    sea_hi: f64,
    iced: usize,
    ice_now: f64,
    /// The steps recorded.
    records: Vec<Record>,
    /// The frost, the sea's ice and the land's snow at the end of each step of the cycle,
    /// kg/m² (F246, F271), and the snow the land's melt took through each step, kg/m², NaN
    /// where its year's snow is unknown (F277); and the snow the step's melt took, kg/m².
    masses: [Vec<f32>; 4],
    melted: f64,
    /// The land's snow, kg/m² of its water, the snow that falls on it through each step of the
    /// cycle, kg/m²/s, none where its year's snow is unknown, and the column as the step began
    /// (F271).
    snow: f64,
    snowfall: Vec<f32>,
    start: Vec<f64>,
    /// The surface's temperature at the end of each step of the cycle (F256), and its sea's
    /// (F268).
    surface: Vec<f32>,
    sea_surface: Vec<f32>,
    /// The place's share of the frost's own gas's caps (F258).
    cap: Option<Patch>,
}

/// A place's share of the caps that hold the air through the lockstep (F258): the share of
/// the place they cover, their column and each layer's sum through the cycle, the ice they
/// hold, kg/m², what they were laid with and held as the cycle began, and what the last
/// cycle left them.
struct Patch {
    share: f64,
    t: Vec<f64>,
    sums: Vec<f64>,
    mass: f64,
    laid: f64,
    start: f64,
    gained: f64,
    work: Work,
    trial: Trial,
}

impl Run {
    fn new(
        (n, steps): (usize, usize),
        place: &Place,
        ocean: Option<&Ocean>,
        frost: Option<&Frost>,
        snowfall: &[f32],
    ) -> Self {
        // The snow the land's melt takes through each step, NaN where its year's snow is not
        // known (F277).
        let melt_unknown = if snowfall.len() == steps {
            0.0
        } else {
            f32::NAN
        };
        let sea = ocean.filter(|_| place.sea > 0.0).map(|o| {
            let t = place.steady.max(o.freezing);
            Sea {
                t,
                ice: 0.0,
                surface: t,
            }
        });
        let cap = frost
            .filter(|_| place.cap > 0.0 && place.cap_mass > 0.0)
            .map(|f| Patch {
                share: place.cap,
                t: vec![place.steady.min(f.point); n],
                sums: vec![0.0; n],
                mass: place.cap_mass,
                laid: place.cap_mass,
                start: place.cap_mass,
                gained: 0.0,
                work: Work::new(n),
                trial: Trial::new(n),
            });
        Self {
            t: vec![place.steady; n],
            mass: 0.0,
            work: Work::new(n),
            trial: Trial::new(n),
            sea,
            cap,
            albedo: place.albedo,
            sums: vec![0.0; n],
            hi: f64::MIN,
            lo: f64::MAX,
            frosted: 0,
            now: 0.0,
            frost_start: 0.0,
            frost_sum: 0.0,
            frost_least: f64::MAX,
            sea_sum: 0.0,
            sea_lo: f64::MAX,
            sea_hi: f64::MIN,
            iced: 0,
            ice_now: 0.0,
            records: Vec::new(),
            masses: [
                vec![0.0; steps],
                vec![0.0; steps],
                vec![0.0; steps],
                vec![melt_unknown; steps],
            ],
            melted: 0.0,
            snow: 0.0,
            snowfall: if snowfall.len() == steps {
                snowfall.to_vec()
            } else {
                Vec::new()
            },
            start: vec![place.steady; n],
            surface: vec![0.0; steps],
            sea_surface: vec![0.0; steps],
        }
    }

    /// The frost the cycle added to frost that lay all through it, kg/m².
    fn growth(&self) -> f64 {
        if self.frost_least > 0.0 {
            (self.mass - self.frost_start).max(0.0)
        } else {
            0.0
        }
    }

    /// The share of the place its frost lies on: its land bare of the caps (F258).
    fn bare_share(&self) -> f64 {
        self.cap.as_ref().map_or(1.0, |c| 1.0 - c.share)
    }

    /// The frost's gas the place holds out of the air over its area, kg/m²: its frost on the
    /// ground bare of the caps, and what its caps took since they were laid (F258).
    fn load(&self) -> f64 {
        self.mass * self.bare_share() + self.taken()
    }

    /// What the place's caps took out of the air since they were laid, kg/m² of the place.
    fn taken(&self) -> f64 {
        self.cap
            .as_ref()
            .map_or(0.0, |c| c.share * (c.mass - c.laid))
    }

    /// A new cycle.
    fn begin(&mut self) {
        self.sums.fill(0.0);
        (self.hi, self.lo) = (f64::MIN, f64::MAX);
        self.frosted = 0;
        self.frost_start = self.mass;
        (self.frost_sum, self.frost_least) = (0.0, f64::MAX);
        (self.sea_sum, self.sea_lo, self.sea_hi) = (0.0, f64::MAX, f64::MIN);
        self.iced = 0;
        self.records.clear();
        if let Some(c) = self.cap.as_mut() {
            c.sums.fill(0.0);
            c.start = c.mass;
        }
    }

    /// The record of the cycle's `k`th step, just taken at `point` under air at `air` K.
    fn record(&mut self, k: usize, snapshot: usize, every: usize, point: f64, air: f64) {
        for (s, v) in self.sums.iter_mut().zip(&self.t) {
            *s += v;
        }
        (self.hi, self.lo) = (self.hi.max(self.t[0]), self.lo.min(self.t[0]));
        if self.mass > 0.0 {
            self.frosted += 1;
        }
        self.frost_sum += self.load();
        self.frost_least = self.frost_least.min(self.mass);
        if let Some(c) = self.cap.as_mut() {
            for (s, v) in c.sums.iter_mut().zip(&c.t) {
                *s += v;
            }
        }
        if let Some(sea) = &self.sea {
            self.sea_sum += sea.surface;
            (self.sea_lo, self.sea_hi) =
                (self.sea_lo.min(sea.surface), self.sea_hi.max(sea.surface));
            if sea.ice > 0.0 {
                self.iced += 1;
            }
            if k == snapshot {
                self.ice_now = sea.ice;
            }
        }
        if k == snapshot {
            self.now = self.mass;
        }
        self.masses[0][k] = self.mass as f32;
        self.masses[1][k] = self.sea.map_or(0.0, |s| s.ice) as f32;
        self.masses[2][k] = self.snow as f32;
        if !self.snowfall.is_empty() {
            self.masses[3][k] = self.melted as f32;
        }
        self.surface[k] = self.t[0] as f32;
        self.sea_surface[k] = self.sea.map_or(0.0, |s| s.surface) as f32;
        if every > 0 && k.is_multiple_of(every) {
            self.records.push(Record {
                k,
                t: self.t.clone(),
                mass: self.mass,
                point,
                air,
                albedo: self.albedo,
            });
        }
    }

    /// The place's sea's year through a cycle of `steps` (F244): its surface's mean, coldest
    /// and warmest, its ice at the snapshot's step and the share of the cycle ice lay; the
    /// land's mean, warmest and coldest, `land`, where the place has no sea.
    fn sea_year(&self, steps: usize, (mean, hi, lo): (f64, f64, f64)) -> [f64; 5] {
        match self.sea {
            Some(_) => [
                self.sea_sum / steps as f64,
                self.sea_lo,
                self.sea_hi,
                self.ice_now,
                self.iced as f64 / steps as f64,
            ],
            None => [mean, lo, hi, 0.0, 0.0],
        }
    }
}

impl Column {
    /// Layers from COLUMN_TOP_SHARE of the skin depth of `day` down to three skin depths of
    /// `deep_period`, each COLUMN_GROWTH times thicker than the one above.
    fn new(ground: &Ground, day: f64, deep_period: f64) -> Self {
        let conductivity = ground.inertia * ground.inertia / ground.heat;
        let diffusivity = conductivity / ground.heat;
        let skin = |period: f64| (diffusivity * period / PI).sqrt();
        let bottom = 3.0 * skin(deep_period.max(day));
        let mut dz = Vec::new();
        let (mut z, mut h) = (0.0, COLUMN_TOP_SHARE * skin(day));
        while (z < bottom || dz.len() < 8) && dz.len() < COLUMN_MAX_LAYERS {
            dz.push(h);
            z += h;
            h *= COLUMN_GROWTH;
        }
        let capacity: Vec<f64> = dz.iter().map(|d| ground.heat * d).collect();
        let conductance = dz
            .windows(2)
            .map(|w| conductivity / (0.5 * (w[0] + w[1])))
            .collect();
        Self {
            capacity,
            conductance,
        }
    }

    fn len(&self) -> usize {
        self.capacity.len()
    }

    /// One step with the frost at its point (F239): the step's trial, ended at the point,
    /// `sheds` the bare ground's radiation and the frosted ground's (F247).
    #[allow(clippy::too_many_arguments)]
    fn frosty_step(
        &self,
        t: &mut [f64],
        (w, tr): (&mut Work, &mut Trial),
        dt: f64,
        gains: (f64, f64),
        bottom: f64,
        (bare, frosted): (
            impl Fn(f64) -> (f64, f64) + Copy,
            impl Fn(f64) -> (f64, f64) + Copy,
        ),
        f: &Frost,
        mass: &mut f64,
    ) {
        let mut start = std::mem::take(&mut w.start);
        start.clear();
        start.extend_from_slice(t);
        self.trial(t, (w, tr), dt, gains, bottom, bare, *mass, f.point);
        let at = frosted(f.point);
        if let Some(lasts) = self.finish(t, tr, dt, f.point, at, f.latent, mass) {
            let ends = (f.point, at.0, f.latent);
            let step = (dt, lasts);
            self.rest_bare(
                t,
                &start,
                (w, tr),
                step,
                (gains.0, bottom),
                bare,
                ends,
                mass,
            );
        }
        w.start = start;
    }

    /// A step's trial under the frost (F239). While frost lies on the ground, or the ground
    /// would cool below the frost point, the ground stands at the point and the step's energy
    /// sublimes or condenses frost instead, at its latent heat (Leighton and Murray 1966); the
    /// gains are the bare ground's and the frosted ground's, the ground under frost gaining
    /// between them by the frost's cover (F257), as its colour shows it. The trial takes bare
    /// ground's step left to itself and, unless that leaves it above `highest`, the column
    /// below a surface held at a point, which stands linear in the point, so that the point
    /// may be settled after.
    #[allow(clippy::too_many_arguments)]
    fn trial(
        &self,
        t: &[f64],
        (w, tr): (&mut Work, &mut Trial),
        dt: f64,
        (bare, frosted): (f64, f64),
        bottom: f64,
        olr: impl Fn(f64) -> (f64, f64) + Copy,
        mass: f64,
        highest: f64,
    ) {
        tr.start = t[0];
        tr.clear = false;
        tr.gain = if mass > 0.0 {
            bare + (frosted - bare) * frost_cover(mass as f32)
        } else {
            bare
        };
        if mass <= 0.0 {
            tr.free.copy_from_slice(t);
            self.step(&mut tr.free, w, dt, bare, bottom, olr);
            if tr.free[0] >= highest {
                tr.clear = true;
                return;
            }
        }
        self.held(t, w, tr, dt, bottom);
    }

    /// The column below a surface held at a point through a step (F239), as the point's
    /// linear response: each layer stands at base + point × lift.
    fn held(&self, t: &[f64], w: &mut Work, tr: &mut Trial, dt: f64, bottom: f64) {
        let n = t.len();
        tr.base[0] = 0.0;
        tr.lift[0] = 1.0;
        if n < 2 {
            return;
        }
        for (i, &ti) in t.iter().enumerate().skip(1) {
            w.b[i] = self.capacity[i] / dt;
            tr.base[i] = self.capacity[i] / dt * ti;
            tr.lift[i] = 0.0;
            w.a[i] = 0.0;
            w.c[i] = 0.0;
            let g = self.conductance[i - 1];
            if i == 1 {
                tr.lift[i] = g;
            } else {
                w.a[i] = -g;
            }
            w.b[i] += g;
            if i + 1 < n {
                let g = self.conductance[i];
                w.c[i] = -g;
                w.b[i] += g;
            }
        }
        tr.base[n - 1] += bottom;
        for i in 2..n {
            let m = w.a[i] / w.b[i - 1];
            w.b[i] -= m * w.c[i - 1];
            tr.base[i] -= m * tr.base[i - 1];
            tr.lift[i] -= m * tr.lift[i - 1];
        }
        tr.base[n - 1] /= w.b[n - 1];
        tr.lift[n - 1] /= w.b[n - 1];
        for i in (1..n - 1).rev() {
            tr.base[i] = (tr.base[i] - w.c[i] * tr.base[i + 1]) / w.b[i];
            tr.lift[i] = (tr.lift[i] - w.c[i] * tr.lift[i + 1]) / w.b[i];
        }
    }

    /// The frost a trial step takes at a point, kg/m² (F239), from the surface layer's balance
    /// through the step: what it needs to reach the point, less what it gains and what its
    /// column below brings it, more what it radiates, the balance left over latent heat. As
    /// found, and bounded: bare ground only takes frost, and frost gives back no more than lies.
    fn deposit(
        &self,
        tr: &Trial,
        dt: f64,
        point: f64,
        radiated: f64,
        latent: f64,
        mass: f64,
    ) -> (f64, f64) {
        if mass <= 0.0 && (tr.clear || tr.free[0] >= point) {
            return (0.0, 0.0);
        }
        let t1 = tr.base.get(1).map_or(point, |b| b + point * tr.lift[1]);
        let below = self.conductance.first().map_or(0.0, |g| g * (t1 - point));
        let short = self.capacity[0] * (point - tr.start) / dt - tr.gain - below + radiated;
        let found = short * dt / latent;
        (found, found.max(-mass.max(0.0)))
    }

    /// Ends a trial step at the point its air settled at (F239): bare ground the step leaves
    /// above the point keeps its free step; frosted ground, or ground that would cool below
    /// the point, stands at it and takes or gives frost, and where the frost runs out what was
    /// left to sublime warms the ground past the point. Returns, where the frost ran out
    /// within the step, the share of the step it lasted, its mass over what the step would
    /// have sublimed (F275).
    #[allow(clippy::too_many_arguments)]
    fn finish(
        &self,
        t: &mut [f64],
        tr: &Trial,
        dt: f64,
        point: f64,
        (radiated, slope): (f64, f64),
        latent: f64,
        mass: &mut f64,
    ) -> Option<f64> {
        if *mass <= 0.0 && (tr.clear || tr.free[0] >= point) {
            t.copy_from_slice(&tr.free);
            return None;
        }
        let (found, taken) = self.deposit(tr, dt, point, radiated, latent, *mass);
        let lay = *mass;
        for (v, (b, l)) in t.iter_mut().zip(tr.base.iter().zip(&tr.lift)) {
            *v = b + point * l;
        }
        *mass = (*mass + taken).max(0.0);
        if taken > found {
            let g = self.conductance.first().copied().unwrap_or(0.0);
            let respond = self.capacity[0] + dt * (slope + g);
            t[0] = point + (taken - found) * latent / respond;
            return Some((lay / -found).clamp(0.0, 1.0));
        }
        None
    }

    /// The step again where its frost ran out within it (F275): from the column `start` the
    /// step began with, the ground stands at the point through the share `lasts` of the step
    /// the frost lasted, the column below it as under a surface held there; through the rest
    /// the bare ground goes free at its own gain and radiation, `gain` and `bare`, and where it
    /// would cool past the point there it stands at the point again and takes frost, at the
    /// frosted ground's radiation there, `radiated`: so that a frost darker than its ground,
    /// which gains enough under the sun to sublime all it holds, leaves its ground the rest of
    /// the step to frost again, rather than a step bare.
    #[allow(clippy::too_many_arguments)]
    fn rest_bare(
        &self,
        t: &mut [f64],
        start: &[f64],
        (w, tr): (&mut Work, &mut Trial),
        (dt, lasts): (f64, f64),
        (gain, bottom): (f64, f64),
        bare: impl Fn(f64) -> (f64, f64) + Copy,
        (point, radiated, latent): (f64, f64, f64),
        mass: &mut f64,
    ) {
        let (held, rest) = (lasts * dt, (1.0 - lasts) * dt);
        t.copy_from_slice(start);
        if held > 0.0 {
            tr.start = start[0];
            self.held(start, w, tr, held, bottom);
            for (v, (b, l)) in t.iter_mut().zip(tr.base.iter().zip(&tr.lift)) {
                *v = b + point * l;
            }
        }
        *mass = 0.0;
        if rest <= 0.0 {
            return;
        }
        tr.free.copy_from_slice(t);
        self.step(&mut tr.free, w, rest, gain, bottom, bare);
        if tr.free[0] >= point {
            t.copy_from_slice(&tr.free);
            return;
        }
        // The bare ground would cool past the point: it stands there and takes frost.
        let from = t.to_vec();
        (tr.start, tr.gain, tr.clear) = (from[0], gain, false);
        self.held(&from, w, tr, rest, bottom);
        let (found, _) = self.deposit(tr, rest, point, radiated, latent, 0.0);
        for (v, (b, l)) in t.iter_mut().zip(tr.base.iter().zip(&tr.lift)) {
            *v = b + point * l;
        }
        *mass = found.max(0.0);
    }

    /// Ends a step on snow (F271): where snow lies and the step left the surface above the
    /// snow's melting point `point`, the surface stands at the point through the share of the
    /// step the snow lasts instead, from the column `start` the step began with, the column
    /// below it as it stands under a surface held there (F239), and the step's warmth melts
    /// snow at its latent heat of fusion `latent`; past it the ground warms as the free step
    /// left it.
    #[allow(clippy::too_many_arguments)]
    fn thaw(
        &self,
        t: &mut [f64],
        start: &[f64],
        (w, tr): (&mut Work, &mut Trial),
        dt: f64,
        (gain, bottom): (f64, f64),
        (radiated, _): (f64, f64),
        (point, latent): (f64, f64),
        snow: &mut f64,
    ) {
        tr.start = start[0];
        tr.gain = gain;
        tr.clear = false;
        self.held(start, w, tr, dt, bottom);
        let (found, _) = self.deposit(tr, dt, point, radiated, latent, *snow);
        if found >= 0.0 {
            return;
        }
        // The snow lasts the share of the step it can take the step's warmth for: through it
        // the ground stands at the point, past it the ground warms as the free step left it.
        let melted = (-found).min(*snow);
        let lasts = melted / -found;
        *snow -= melted;
        for (v, (b, l)) in t.iter_mut().zip(tr.base.iter().zip(&tr.lift)) {
            *v = lasts * (b + point * l) + (1.0 - lasts) * *v;
        }
    }

    /// One implicit step of `dt`: the surface gains `gain` W/m² and radiates as `olr` gives
    /// (its value and slope), linearised about its temperature; `bottom` W/m² flows into the
    /// lowest layer from below.
    fn step(
        &self,
        t: &mut [f64],
        w: &mut Work,
        dt: f64,
        gain: f64,
        bottom: f64,
        olr: impl Fn(f64) -> (f64, f64) + Copy,
    ) {
        let n = t.len();
        for (i, &ti) in t.iter().enumerate() {
            w.b[i] = self.capacity[i] / dt;
            w.d[i] = self.capacity[i] / dt * ti;
            w.a[i] = 0.0;
            w.c[i] = 0.0;
            if i > 0 {
                let g = self.conductance[i - 1];
                w.a[i] = -g;
                w.b[i] += g;
            }
            if i + 1 < n {
                let g = self.conductance[i];
                w.c[i] = -g;
                w.b[i] += g;
            }
        }
        let (r, dr) = olr(t[0]);
        w.b[0] += dr;
        w.d[0] += gain - r + dr * t[0];
        w.d[n - 1] += bottom;
        for i in 1..n {
            let m = w.a[i] / w.b[i - 1];
            w.b[i] -= m * w.c[i - 1];
            w.d[i] -= m * w.d[i - 1];
        }
        t[n - 1] = w.d[n - 1] / w.b[n - 1];
        for i in (0..n - 1).rev() {
            t[i] = (w.d[i] - w.c[i] * t[i + 1]) / w.b[i];
        }
        for v in t.iter_mut() {
            *v = v.max(COSMIC_BACKGROUND_K);
        }
    }
}

/// A field's area-weighted mean over the body, and over the cells within ten degrees of the
/// equator and beyond eighty of latitude.
pub fn bands(f: &Field<f32>) -> (f64, f64, f64) {
    let g = f.grid();
    let mut acc = [(0.0, 0.0); 3];
    for c in g.cells() {
        let a = g.cell_area(c);
        let v = f.get(c) as f64;
        let (_, lat) = lon_lat_from_direction(g.center(c));
        let lat = lat.to_degrees();
        acc[0].0 += v * a;
        acc[0].1 += a;
        if lat.abs() < 10.0 {
            acc[1].0 += v * a;
            acc[1].1 += a;
        }
        if lat.abs() > 80.0 {
            acc[2].0 += v * a;
            acc[2].1 += a;
        }
    }
    let mean = |(sum, area): (f64, f64)| if area > 0.0 { sum / area } else { f64::NAN };
    (mean(acc[0]), mean(acc[1]), mean(acc[2]))
}

/// Each cell's four neighbours with the weight of the edge between them, its length over
/// the distance between the centres, so that the sum of weight times difference is the
/// cell's integral of the Laplacian on the unit sphere.
fn links(grid: CubeSphere) -> Vec<[(usize, f64); 4]> {
    let delta = grid.delta();
    let angle = |a: [f64; 3], b: [f64; 3]| {
        let c = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
        let x = [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ];
        (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt().atan2(c)
    };
    (0..grid.cell_count())
        .into_par_iter()
        .map(|k| {
            let c = grid.cell(k);
            let (a, b) = grid.angles(c);
            let centre = grid.center(c);
            let mut out = [(0usize, 0.0f64); 4];
            for (slot, (di, dj)) in [(1, 0), (-1, 0), (0, 1), (0, -1)].into_iter().enumerate() {
                let nb = grid.neighbor(c, di, dj);
                let (ea, eb) = (a + di as f64 * delta / 2.0, b + dj as f64 * delta / 2.0);
                // The shared edge's two corners.
                let (p, q) = if di != 0 {
                    (
                        CubeSphere::direction(c.face, ea, eb - delta / 2.0),
                        CubeSphere::direction(c.face, ea, eb + delta / 2.0),
                    )
                } else {
                    (
                        CubeSphere::direction(c.face, ea - delta / 2.0, eb),
                        CubeSphere::direction(c.face, ea + delta / 2.0, eb),
                    )
                };
                let w = angle(p, q) / angle(centre, grid.center(nb)).max(1e-12);
                out[slot] = (grid.index(nb), w);
            }
            out
        })
        .collect()
}

/// The air over a body's ground (F237): the mass of it the surface wind brings a square
/// metre of ground through the bulk transfer coefficient, rho C_H U, kg/m²/s, none without an
/// atmosphere, and its heat capacity, J/kg/K.
fn ventilation(graph: &ParamGraph) -> (f64, f64) {
    let empty = std::collections::BTreeMap::new();
    let composition = graph
        .get("atmosphere.composition")
        .as_map()
        .unwrap_or(&empty);
    let heat = mean_heat_capacity(composition);
    let pressure = graph.f64("atmosphere.surface_pressure");
    if pressure <= 0.0 {
        return (0.0, heat);
    }
    let mu = mean_molecular_mass(composition);
    let t = graph.f64("atmosphere.surface_temperature_mean").max(1.0);
    let density = pressure * mu / (GAS_CONSTANT * t);
    (
        density * AIR_EXCHANGE_COEFFICIENT * SURFACE_WIND_SPEED,
        heat,
    )
}

/// The places of a layout of `n_lat` latitudes and `n_lon` longitudes a point is read from,
/// with their weights, as `read_row` and `read_grid` read them.
pub(crate) fn reader(n_lat: usize, n_lon: usize, lat: f64, lon: f64) -> [(usize, f64); 4] {
    let y = ((lat + PI / 2.0) / PI * (n_lat - 1) as f64).clamp(0.0, (n_lat - 1) as f64);
    let i = (y.floor() as usize).min(n_lat - 2);
    let fy = y - i as f64;
    if n_lon == 1 {
        return [(i, 1.0 - fy), (i + 1, fy), (0, 0.0), (0, 0.0)];
    }
    let x = (lon + PI).rem_euclid(TAU) / TAU * n_lon as f64;
    let j = (x.floor() as usize) % n_lon;
    let fx = x - x.floor();
    let j1 = (j + 1) % n_lon;
    [
        (i * n_lon + j, (1.0 - fx) * (1.0 - fy)),
        (i * n_lon + j1, fx * (1.0 - fy)),
        ((i + 1) * n_lon + j, (1.0 - fx) * fy),
        ((i + 1) * n_lon + j1, fx * fy),
    ]
}

/// Solves a transport's implicit step over a layout (F245): (diag_p + sum_q k_pq) x_p -
/// sum_q k_pq x_q = rhs_p at every place p, `links` each place's neighbours q with the
/// conductance k_pq per unit of its area, `x` holding the start. A fast spinner's latitudes
/// directly; a resonant body's latitudes and longitudes by conjugate gradients on the system
/// weighted by the places' areas, which makes it symmetric, preconditioned, where the system
/// is `zonal`, nearly the same along each latitude, by its mean along them solved exactly
/// (`Zonal`), and otherwise by each latitude's row solved exactly, since the rows' places
/// crowd together toward the poles. Returns the conjugate gradients' steps.
fn transport_solve(
    layout: &Layout,
    links: &[Vec<(usize, f64)>],
    (diag, rhs): (&[f64], &[f64]),
    x: &mut [f64],
    zonal: bool,
) -> usize {
    let n = diag.len();
    let total: Vec<f64> = (0..n)
        .map(|p| diag[p] + links[p].iter().map(|l| l.1).sum::<f64>())
        .collect();
    let link = |p: usize, q: usize| links[p].iter().find(|l| l.0 == q).map_or(0.0, |l| l.1);
    if layout.n_lon == 1 {
        let (mut b, mut c, mut r) = (total.clone(), vec![0.0; n], rhs.to_vec());
        for p in 0..n {
            if p + 1 < n {
                c[p] = -link(p, p + 1);
            }
            if p > 0 {
                let m = -link(p, p - 1) / b[p - 1];
                b[p] -= m * c[p - 1];
                r[p] -= m * r[p - 1];
            }
        }
        x[n - 1] = r[n - 1] / b[n - 1];
        for p in (0..n - 1).rev() {
            x[p] = (r[p] - c[p] * x[p + 1]) / b[p];
        }
        return 0;
    }
    let n_lon = layout.n_lon;
    let w: Vec<f64> = layout.places.iter().map(|p| p.weight).collect();
    let rows: Vec<Cyclic> = if zonal {
        Vec::new()
    } else {
        (0..n / n_lon)
            .map(|i| {
                let row = i * n_lon;
                let along: Vec<f64> = (0..n_lon)
                    .map(|j| link(row + j, row + (j + 1) % n_lon))
                    .collect();
                Cyclic::new(&total[row..row + n_lon], &along)
            })
            .collect()
    };
    let mean = zonal.then(|| Zonal::new(layout, &link, &total));
    let mut spectrum = vec![0.0; 2 * (n_lon / 2 + 1) * (n / n_lon)];
    let apply = |v: &[f64], out: &mut [f64]| {
        for p in 0..n {
            let mut s = total[p] * v[p];
            for &(q, k) in &links[p] {
                s -= k * v[q];
            }
            out[p] = w[p] * s;
        }
    };
    let mut scaled = vec![0.0; n_lon];
    let mut precondition = |r: &[f64], z: &mut [f64]| {
        if let Some(mean) = &mean {
            for p in 0..n {
                z[p] = r[p] / w[p];
            }
            mean.solve(layout, z, &mut spectrum);
            return;
        }
        for (i, row) in rows.iter().enumerate() {
            let span = i * n_lon..(i + 1) * n_lon;
            for (v, &x) in scaled.iter_mut().zip(&r[span.clone()]) {
                *v = x / w[span.start];
            }
            row.solve(&scaled, &mut z[span]);
        }
    };
    let dot = |a: &[f64], b: &[f64]| a.iter().zip(b).map(|(x, y)| x * y).sum::<f64>();
    let target: Vec<f64> = (0..n).map(|p| w[p] * rhs[p]).collect();
    let norm = dot(&target, &target).sqrt();
    if norm == 0.0 {
        x.fill(0.0);
        return 0;
    }
    let mut q = vec![0.0; n];
    apply(x, &mut q);
    let mut r: Vec<f64> = (0..n).map(|p| target[p] - q[p]).collect();
    let mut z = vec![0.0; n];
    precondition(&r, &mut z);
    let mut s = z.clone();
    let mut rz = dot(&r, &z);
    let mut steps = 0;
    while steps < CLIMATE_CG_ITERATIONS {
        if dot(&r, &r).sqrt() < CLIMATE_CG_TOLERANCE * norm {
            break;
        }
        steps += 1;
        apply(&s, &mut q);
        let alpha = rz / dot(&s, &q);
        for p in 0..n {
            x[p] += alpha * s[p];
            r[p] -= alpha * q[p];
        }
        precondition(&r, &mut z);
        let next = dot(&r, &z);
        let beta = next / rz;
        rz = next;
        for p in 0..n {
            s[p] = z[p] + beta * s[p];
        }
    }
    steps
}

/// A transport's system over a resonant body's layout averaged along its latitudes (F245):
/// each latitude's mean diagonal and its mean couplings along it and to the latitudes
/// beside it, solved exactly by the longitudes' Fourier modes, each a tridiagonal system over
/// the latitudes, factored once for the right-hand sides it solves.
struct Zonal {
    n_lat: usize,
    modes: usize,
    /// Each mode's tridiagonal, factored: the coupling to the latitude below over the pivot
    /// below it, the coupling to the latitude above and the pivots, mode after mode.
    lower: Vec<f64>,
    upper: Vec<f64>,
    pivots: Vec<f64>,
}

impl Zonal {
    fn new(layout: &Layout, link: &dyn Fn(usize, usize) -> f64, total: &[f64]) -> Self {
        let n_lon = layout.n_lon;
        let n_lat = total.len() / n_lon;
        let modes = n_lon / 2 + 1;
        let mean = |i: usize, f: &dyn Fn(usize, usize) -> f64| {
            (0..n_lon).map(|j| f(i * n_lon + j, j)).sum::<f64>() / n_lon as f64
        };
        let diag: Vec<f64> = (0..n_lat).map(|i| mean(i, &|p, _| total[p])).collect();
        let along: Vec<f64> = (0..n_lat)
            .map(|i| mean(i, &|p, j| link(p, i * n_lon + (j + 1) % n_lon)))
            .collect();
        let below: Vec<f64> = (0..n_lat)
            .map(|i| {
                if i == 0 {
                    0.0
                } else {
                    mean(i, &|p, _| link(p, p - n_lon))
                }
            })
            .collect();
        let above: Vec<f64> = (0..n_lat)
            .map(|i| {
                if i + 1 == n_lat {
                    0.0
                } else {
                    mean(i, &|p, _| link(p, p + n_lon))
                }
            })
            .collect();
        let size = modes * n_lat;
        let (mut lower, mut upper, mut pivots) =
            (vec![0.0; size], vec![0.0; size], vec![0.0; size]);
        for m in 0..modes {
            let shift = 2.0 * (TAU * m as f64 / n_lon as f64).cos();
            for i in 0..n_lat {
                let k = m * n_lat + i;
                let b = diag[i] - shift * along[i];
                upper[k] = -above[i];
                if i == 0 {
                    pivots[k] = b;
                } else {
                    lower[k] = -below[i] / pivots[k - 1];
                    pivots[k] = b - lower[k] * upper[k - 1];
                }
            }
        }
        Self {
            n_lat,
            modes,
            lower,
            upper,
            pivots,
        }
    }

    /// Solves the averaged system in place: `z` the right-hand side per unit area, then the
    /// solution; `spectrum` room for the modes.
    fn solve(&self, layout: &Layout, z: &mut [f64], spectrum: &mut [f64]) {
        let (n_lat, n_lon, modes) = (self.n_lat, layout.n_lon, self.modes);
        let [cos, sin, cos_t, sin_t] = &layout.modes;
        let (a, b) = spectrum.split_at_mut(modes * n_lat);
        // Each latitude's modes, latitude after latitude.
        a.fill(0.0);
        b.fill(0.0);
        for i in 0..n_lat {
            let (am, bm) = (
                &mut a[i * modes..(i + 1) * modes],
                &mut b[i * modes..(i + 1) * modes],
            );
            for (j, &x) in z[i * n_lon..(i + 1) * n_lon].iter().enumerate() {
                let (c, s) = (
                    &cos_t[j * modes..(j + 1) * modes],
                    &sin_t[j * modes..(j + 1) * modes],
                );
                for ((u, v), (c, s)) in am.iter_mut().zip(bm.iter_mut()).zip(c.iter().zip(s)) {
                    *u += x * c;
                    *v += x * s;
                }
            }
        }
        // Each mode's tridiagonal over the latitudes.
        for m in 0..modes {
            let k = |i: usize| m * n_lat + i;
            for v in [&mut *a, &mut *b] {
                for i in 1..n_lat {
                    v[i * modes + m] -= self.lower[k(i)] * v[(i - 1) * modes + m];
                }
                v[(n_lat - 1) * modes + m] /= self.pivots[k(n_lat - 1)];
                for i in (0..n_lat - 1).rev() {
                    v[i * modes + m] = (v[i * modes + m]
                        - self.upper[k(i)] * v[(i + 1) * modes + m])
                        / self.pivots[k(i)];
                }
            }
        }
        z.fill(0.0);
        for i in 0..n_lat {
            let row = &mut z[i * n_lon..(i + 1) * n_lon];
            for m in 0..modes {
                let weight = if m == 0 || 2 * m == n_lon { 1.0 } else { 2.0 } / n_lon as f64;
                let (am, bm) = (weight * a[i * modes + m], weight * b[i * modes + m]);
                let (c, s) = (
                    &cos[m * n_lon..(m + 1) * n_lon],
                    &sin[m * n_lon..(m + 1) * n_lon],
                );
                for ((v, c), s) in row.iter_mut().zip(c).zip(s) {
                    *v += am * c + bm * s;
                }
            }
        }
    }
}

/// A cyclic tridiagonal system of a layout's row (F245): its diagonal, and `along` each
/// place's coupling to the next, -along off the diagonal, the last's to the first; factored
/// once for the right-hand sides it solves, by Sherman and Morrison's correction of the
/// tridiagonal system without its corners.
struct Cyclic {
    along: Vec<f64>,
    gamma: f64,
    pivots: Vec<f64>,
    z: Vec<f64>,
    denominator: f64,
}

impl Cyclic {
    fn new(diag: &[f64], along: &[f64]) -> Self {
        let n = diag.len();
        let corner = along[n - 1];
        let gamma = -diag[0];
        let mut b = diag.to_vec();
        b[0] -= gamma;
        b[n - 1] -= corner * corner / gamma;
        let mut pivots = vec![0.0; n];
        pivots[0] = b[0];
        for j in 1..n {
            pivots[j] = b[j] - along[j - 1] * along[j - 1] / pivots[j - 1];
        }
        let mut row = Self {
            along: along.to_vec(),
            gamma,
            pivots,
            z: Vec::new(),
            denominator: 1.0,
        };
        let mut u = vec![0.0; n];
        (u[0], u[n - 1]) = (gamma, -corner);
        let mut z = vec![0.0; n];
        row.tridiagonal(&u, &mut z);
        row.denominator = 1.0 + z[0] - corner * z[n - 1] / gamma;
        row.z = z;
        row
    }

    /// The tridiagonal system without the corners.
    fn tridiagonal(&self, r: &[f64], x: &mut [f64]) {
        let n = r.len();
        x[0] = r[0];
        for j in 1..n {
            x[j] = r[j] + self.along[j - 1] / self.pivots[j - 1] * x[j - 1];
        }
        x[n - 1] /= self.pivots[n - 1];
        for j in (0..n - 1).rev() {
            x[j] = (x[j] + self.along[j] * x[j + 1]) / self.pivots[j];
        }
    }

    fn solve(&self, r: &[f64], x: &mut [f64]) {
        let n = r.len();
        let corner = self.along[n - 1];
        self.tridiagonal(r, x);
        let fact = (x[0] - corner * x[n - 1] / self.gamma) / self.denominator;
        for (v, z) in x.iter_mut().zip(&self.z) {
            *v -= fact * z;
        }
    }
}

/// The sea's water's share of the transport through a step (F245): the open seas' mixed
/// layers exchange heat implicitly across the sea's share of each boundary between their
/// places, the lesser of the two, at `conductance` per unit of the layout's links; a sea under
/// ice stands at its freezing point, and the heat the water brings it melts its ice from
/// below: each place's `inflow`, W/m² of its sea.
fn sea_transport(
    layout: &Layout,
    runs: &mut [Run],
    (conductance, dt): (f64, f64),
    ocean: &Ocean,
    inflow: &mut [f64],
) {
    let places = &layout.places;
    let n = places.len();
    let state = |r: &Run| r.sea.as_ref().map(|s| (s.ice > 0.0, s.t));
    let states: Vec<Option<(bool, f64)>> = runs.iter().map(state).collect();
    let hold = ocean.capacity / dt;
    let (mut diag, mut rhs, mut x) = (vec![1.0; n], vec![0.0; n], vec![0.0; n]);
    let mut links: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
    for p in 0..n {
        match states[p] {
            Some((false, t)) => {
                diag[p] = places[p].sea * hold;
                rhs[p] = diag[p] * t;
                x[p] = t;
                for &(q, g) in &layout.links[p] {
                    let k = conductance * g * places[p].sea.min(places[q].sea);
                    match states[q] {
                        Some((false, _)) if k > 0.0 => links[p].push((q, k)),
                        Some((true, _)) if k > 0.0 => {
                            diag[p] += k;
                            rhs[p] += k * ocean.freezing;
                        }
                        _ => {}
                    }
                }
            }
            Some((true, _)) => (rhs[p], x[p]) = (ocean.freezing, ocean.freezing),
            None => {}
        }
    }
    transport_solve(layout, &links, (&diag, &rhs), &mut x, false);
    for p in 0..n {
        inflow[p] = 0.0;
        match states[p] {
            Some((false, _)) => {
                if let Some(sea) = runs[p].sea.as_mut() {
                    (sea.t, sea.surface) = (x[p], x[p]);
                }
            }
            Some((true, _)) => {
                let brought: f64 = layout.links[p]
                    .iter()
                    .filter(|&&(q, _)| matches!(states[q], Some((false, _))))
                    .map(|&(q, g)| {
                        conductance * g * places[p].sea.min(places[q].sea) * (x[q] - ocean.freezing)
                    })
                    .sum();
                inflow[p] = brought / places[p].sea.max(1e-12);
            }
            None => {}
        }
    }
}

/// The seas' currents through a step (F272): the water each boundary between two places with
/// seas passes brings its upwind place's mixed layer's heat to the place downwind of it, in
/// steps that move no more than OCEAN_CFL of a place's water, each place's water changing by
/// what it brings over what it holds, so that where the currents meet the land a sea of one
/// temperature stays at it; a sea under ice stands at its freezing point, and the heat the
/// water brings it melts its ice from below, adding to each place's `inflow`, W/m² of its sea.
fn sea_currents(
    layout: &Layout,
    runs: &mut [Run],
    (currents, dt): (&crate::currents::Currents, f64),
    ocean: &Ocean,
    inflow: &mut [f64],
) {
    let places = &layout.places;
    let (n, n_lon) = (places.len(), layout.n_lon);
    if currents.east.len() != n || n_lon == 0 {
        return;
    }
    let depth = OCEAN_MIXED_LAYER_M;
    let area: Vec<f64> = places
        .iter()
        .map(|p| TAU * currents.radius * currents.radius * p.weight)
        .collect();
    let water: Vec<f64> = (0..n).map(|p| depth * area[p] * places[p].sea).collect();
    let iced: Vec<Option<bool>> = runs
        .iter()
        .map(|r| r.sea.as_ref().map(|s| s.ice > 0.0))
        .collect();
    let mut temps: Vec<f64> = runs
        .iter()
        .zip(&iced)
        .map(|(r, i)| match (r.sea.as_ref(), i) {
            (Some(_), Some(true)) => ocean.freezing,
            (Some(s), _) => s.t,
            (None, _) => 0.0,
        })
        .collect();
    // Each boundary between two places with seas, and the water it passes from the first to
    // the second, m³/s.
    let mut faces = Vec::with_capacity(2 * n);
    for p in 0..n {
        let (i, j) = (p / n_lon, p % n_lon);
        let east = i * n_lon + (j + 1) % n_lon;
        let north = p + n_lon;
        for (q, f) in [(east, currents.east[p]), (north, currents.north[p])] {
            if q < n && q != p && iced[p].is_some() && iced[q].is_some() && f != 0.0 {
                faces.push(if f > 0.0 { (p, q, f) } else { (q, p, -f) });
            }
        }
    }
    let mut given = vec![0.0; n];
    for &(from, _, f) in &faces {
        given[from] += f;
    }
    let fastest = (0..n)
        .filter(|&p| water[p] > 0.0)
        .map(|p| given[p] / water[p])
        .fold(0.0, f64::max);
    let subs = ((dt * fastest / OCEAN_CFL).ceil() as usize).clamp(1, OCEAN_MOST_SUBSTEPS);
    let h = dt / subs as f64;
    let mut brought = vec![0.0; n];
    let mut change = vec![0.0; n];
    for _ in 0..subs {
        change.fill(0.0);
        for &(from, to, f) in &faces {
            change[to] += f * h * (temps[from] - temps[to]);
        }
        for p in 0..n {
            match iced[p] {
                Some(true) => brought[p] += change[p],
                Some(false) if water[p] > 0.0 => temps[p] += change[p] / water[p],
                _ => {}
            }
        }
    }
    let heat = ocean.capacity / depth;
    for (p, r) in runs.iter_mut().enumerate() {
        match (r.sea.as_mut(), iced[p]) {
            (Some(sea), Some(false)) => (sea.t, sea.surface) = (temps[p], temps[p]),
            (Some(_), Some(true)) if area[p] * places[p].sea > 0.0 => {
                inflow[p] += brought[p] * heat / (dt * area[p] * places[p].sea);
            }
            _ => {}
        }
    }
}

/// Solves (diag(slope) + D L) x = r by conjugate gradients with the diagonal as the
/// preconditioner, L the weighted graph Laplacian of the links.
fn conjugate_gradients(slope: &[f64], links: &[[(usize, f64); 4]], d: f64, r: &[f64]) -> Vec<f64> {
    let n = slope.len();
    let apply = |v: &[f64]| -> Vec<f64> {
        (0..n)
            .into_par_iter()
            .map(|k| {
                let mut s = slope[k] * v[k];
                for &(j, w) in &links[k] {
                    s += d * w * (v[k] - v[j]);
                }
                s
            })
            .collect()
    };
    let diag: Vec<f64> = (0..n)
        .map(|k| slope[k] + d * links[k].iter().map(|l| l.1).sum::<f64>())
        .collect();
    let dot = |a: &[f64], b: &[f64]| -> f64 { a.par_iter().zip(b).map(|(x, y)| x * y).sum() };
    let mut x = vec![0.0; n];
    let mut res = r.to_vec();
    let mut z: Vec<f64> = (0..n).map(|k| res[k] / diag[k]).collect();
    let mut p = z.clone();
    let mut rz = dot(&res, &z);
    let norm0 = dot(r, r).sqrt();
    if norm0 == 0.0 {
        return x;
    }
    for _ in 0..CLIMATE_CG_ITERATIONS {
        let ap = apply(&p);
        let alpha = rz / dot(&p, &ap);
        for k in 0..n {
            x[k] += alpha * p[k];
            res[k] -= alpha * ap[k];
        }
        if dot(&res, &res).sqrt() < CLIMATE_CG_TOLERANCE * norm0 {
            break;
        }
        z = (0..n).map(|k| res[k] / diag[k]).collect();
        let rz_next = dot(&res, &z);
        let beta = rz_next / rz;
        rz = rz_next;
        for k in 0..n {
            p[k] = z[k] + beta * p[k];
        }
    }
    x
}

#[cfg(test)]
mod tests {
    use super::*;

    fn climate_of(name: &str) -> (Climate, EnergyBalance) {
        let graph = crate::presets::preset(name).unwrap().to_graph().0;
        let eb = EnergyBalance::from_graph(&graph);
        (eb.solve(CubeSphere::new(CLIMATE_FACE_SIZE)), eb)
    }

    /// Prints each preset's year of climate: the mean, the tropics and the poles, and the
    /// steps it took.
    #[test]
    #[ignore = "diagnostic: the presets' climates"]
    fn print_climate() {
        for name in crate::presets::PRESET_NAMES {
            let (climate, eb) = climate_of(name);
            let (mean, equator, poles) = bands(&climate.temperature);
            println!(
                "{name:10} mean {mean:6.1} K  equator {equator:6.1}  poles {poles:6.1}  (D {:.3}, albedo {:.2}, low sun {:.3}, ice {}, {} steps)",
                eb.transport, eb.albedo, eb.low_sun, eb.ice, climate.iterations
            );
        }
    }

    /// A body's sea, where water is its working liquid, over Earth's share of every cell of
    /// the climate's grid, 0.71, or EXOGEN_SEA_SHARE's.
    fn sea_of(graph: &ParamGraph) -> Option<(Ocean, Field<f32>)> {
        let ocean = Ocean::from_graph(graph)?;
        let share = std::env::var("EXOGEN_SEA_SHARE").map_or(0.71, |v| v.parse().unwrap());
        let grid = CubeSphere::new(CLIMATE_FACE_SIZE);
        Some((ocean, Field::filled(grid, share as f32)))
    }

    /// A body's year under its air, over its sea (F245).
    fn graph_cycle(graph: &ParamGraph) -> (Cycle, Climate) {
        let eb = EnergyBalance::from_graph(graph);
        let climate = eb.solve(CubeSphere::new(CLIMATE_FACE_SIZE));
        let sea = sea_of(graph);
        let sea = sea.as_ref().map(|(o, f)| (o, f));
        (
            eb.cycle(&climate, &Ground::from_graph(graph), sea, None),
            climate,
        )
    }

    fn cycle_of(name: &str) -> (Cycle, Climate) {
        graph_cycle(&crate::presets::preset(name).unwrap().to_graph().0)
    }

    /// A year's mean over the body, its land's and its sea's by their shares (F245).
    fn surface_mean(cycle: &Cycle, graph: &ParamGraph) -> f64 {
        let (land, ..) = bands(&cycle.mean);
        match (&cycle.seas, sea_of(graph)) {
            (Some(seas), Some((_, share))) => {
                let f = share.as_slice()[0] as f64;
                (1.0 - f) * land + f * bands(&seas.mean).0
            }
            _ => land,
        }
    }

    fn at(f: &Field<f32>, lat_deg: f64, lon_deg: f64) -> f64 {
        let g = f.grid();
        let d = direction_from_lon_lat(lon_deg.to_radians(), lat_deg.to_radians());
        f.get(g.locate(d).0) as f64
    }

    /// Prints each preset's year through its days and seasons against its steady balance.
    #[test]
    #[ignore = "diagnostic: the presets' days and seasons"]
    fn print_cycle() {
        let names = std::env::var("EXOGEN_BODY").unwrap_or_default();
        for name in crate::presets::PRESET_NAMES {
            if !names.is_empty() && !names.split(',').any(|n| n == *name) {
                continue;
            }
            let start = std::time::Instant::now();
            let graph = crate::presets::preset(name).unwrap().to_graph().0;
            let (cycle, climate) = graph_cycle(&graph);
            let took = start.elapsed().as_secs_f64();
            let (_, equator, poles) = bands(&cycle.mean);
            let mean = surface_mean(&cycle, &graph);
            let (steady, ..) = bands(&climate.temperature);
            println!(
                "{name:10} mean {mean:6.1} K (steady {steady:6.1})  equator {equator:6.1} from {:6.1} to {:6.1}  poles {poles:6.1} from {:6.1} to {:6.1}  90° east at noon {:6.1}  ({took:.1} s)",
                at(&cycle.min, 0.0, 0.0),
                at(&cycle.max, 0.0, 0.0),
                at(&cycle.min, 89.0, 0.0),
                at(&cycle.max, 89.0, 0.0),
                at(&cycle.max, 0.0, 90.0),
            );
        }
    }

    /// F236: the Moon's equator through its days and nights, its noons near 390 K and its
    /// nights near 95 K (Diviner; Vasavada et al. 2012), its mean far under the radiative
    /// balance with the year's mean sunlight.
    #[test]
    fn the_moon_s_equator_takes_its_days_and_nights() {
        let (cycle, climate) = cycle_of("Moon");
        let max = at(&cycle.max, 0.0, 0.0);
        let min = at(&cycle.min, 0.0, 0.0);
        let mean = at(&cycle.mean, 0.0, 0.0);
        let steady = at(&climate.temperature, 0.0, 0.0);
        assert!((375.0..400.0).contains(&max), "noon {max}");
        assert!((85.0..110.0).contains(&min), "night {min}");
        assert!(
            (195.0..240.0).contains(&mean) && mean < steady - 40.0,
            "mean {mean} steady {steady}"
        );
        // F256: its temperatures through its days, their quantiles from its night's to its
        // noon's, their mean its mean.
        let (spread, _) = cycle.spread_at(direction_from_lon_lat(0.0, 0.0));
        let (coldest, warmest) = (spread[0], spread[CYCLE_SPREAD - 1]);
        let middle = spread.iter().sum::<f64>() / CYCLE_SPREAD as f64;
        assert!(
            coldest > min - 1.0 && coldest < min + 15.0,
            "coldest {coldest}"
        );
        assert!(
            warmest > max - 5.0 && warmest < max + 1.0,
            "warmest {warmest}"
        );
        assert!((middle - mean).abs() < 3.0, "mean {middle} against {mean}");
    }

    /// F236: Mercury's hot longitudes reach some 700 K at perihelion's noon and its warm
    /// longitudes some 570 K at aphelion's, its nights near 100 K (Vasavada, Paige and Wood
    /// 1999).
    #[test]
    fn mercury_s_hot_longitudes_reach_seven_hundred_kelvin() {
        let (cycle, _) = cycle_of("Mercury");
        let hot = at(&cycle.max, 0.0, 0.0);
        let warm = at(&cycle.max, 0.0, 90.0);
        let night = at(&cycle.min, 0.0, 0.0);
        assert!((680.0..720.0).contains(&hot), "hot {hot}");
        assert!((545.0..590.0).contains(&warm), "warm {warm}");
        assert!((80.0..130.0).contains(&night), "night {night}");
    }

    /// F236: Mars's equator swings through its days, the dust-free ground of its tropics
    /// near 290 K at summer's noon and 180 K before dawn (Curiosity's ground temperatures,
    /// Martín-Torres et al. 2015), about its 210 K mean.
    #[test]
    fn mars_s_equator_swings_through_its_days() {
        let (cycle, _) = cycle_of("Mars");
        let max = at(&cycle.max, 0.0, 0.0);
        let min = at(&cycle.min, 0.0, 0.0);
        let mean = at(&cycle.mean, 0.0, 0.0);
        assert!((260.0..310.0).contains(&max), "noon {max}");
        assert!((150.0..200.0).contains(&min), "night {min}");
        assert!((195.0..235.0).contains(&mean), "mean {mean}");
    }

    /// F236, F245: a thick atmosphere holds the day and the seasons: Venus's and Titan's means
    /// through them stand at their steady balances; the Earth-like's, its land's and its
    /// sea's by their shares, within a few kelvin of its steady balance and the gate's 5 K of
    /// Earth's.
    #[test]
    fn a_thick_atmosphere_keeps_its_steady_mean() {
        for (name, within) in [("Venus", 0.5), ("Titan", 0.5), ("Earth", 3.0)] {
            let graph = crate::presets::preset(name).unwrap().to_graph().0;
            let (cycle, climate) = graph_cycle(&graph);
            let mean = surface_mean(&cycle, &graph);
            let (steady, ..) = bands(&climate.temperature);
            assert!(
                (mean - steady).abs() < within,
                "{name}: {mean} against {steady}"
            );
            if name == "Earth" {
                assert!((mean - 288.0).abs() < 5.0, "Earth's mean {mean}");
            }
        }
    }

    /// F236: a regolith's thermal inertia by the gas in its pores: the Moon's 55 in vacuum,
    /// Mars's near the 200 its global maps give for fine-grained ground (Putzig and Mellon
    /// 2007), dry sand's some 500 under a bar.
    #[test]
    fn the_regolith_conducts_by_the_gas_in_its_pores() {
        let inertia = |name: &str| {
            let graph = crate::presets::preset(name).unwrap().to_graph().0;
            graph.f64("surface.thermal_inertia")
        };
        assert!((inertia("Moon") - 55.0).abs() < 1.0);
        assert!(
            (150.0..230.0).contains(&inertia("Mars")),
            "Mars {}",
            inertia("Mars")
        );
        assert!(
            (400.0..560.0).contains(&inertia("Earth")),
            "Earth {}",
            inertia("Earth")
        );
    }

    /// Prints the frost of each body whose atmosphere lays one: its gas and frost point, and
    /// by latitude the coldest the ground stands, the share of the year frost lies on it and
    /// its frost at the snapshot's season.
    #[test]
    #[ignore = "diagnostic: the frost's seasons"]
    fn print_frost() {
        let names = std::env::var("EXOGEN_BODY").unwrap_or_else(|_| "Mars,Triton,Pluto".into());
        for name in names.split(',') {
            let (cycle, _) = cycle_of(name);
            let Some(f) = cycle.frost else {
                println!("{name}: no frost");
                continue;
            };
            println!(
                "{name}: {} frost at {:.1} K, albedo {:.2}, partial {:.4e} Pa, column {:.4} kg/m², pressure {:?}, caps grow {:?} kg/m² a year",
                f.species,
                f.point,
                f.albedo,
                f.partial,
                f.column,
                cycle.pressure_range,
                cycle.cap_growth
            );
            for lat in [
                -90.0, -80.0, -70.0, -60.0, -50.0, -30.0, 0.0, 30.0, 50.0, 60.0, 70.0, 80.0, 90.0,
            ] {
                println!(
                    "  {lat:6.1}°  coldest {:6.1} K  mean {:6.1} K  frost {:4.2} of the year  now {:9.2} kg/m²",
                    at(&cycle.min, lat, 0.0),
                    at(&cycle.mean, lat, 0.0),
                    at(&cycle.frost_share, lat, 0.0),
                    at(&cycle.frost_now, lat, 0.0)
                );
            }
        }
    }

    /// F239, the frost's seasons on Mars (Leighton and Murray 1966): carbon dioxide condenses
    /// on the winter poles and holds them at its frost point, near 145 K; its seasonal caps
    /// reach past 60° of latitude in both hemispheres (James et al. 1992), the northern
    /// laying a metre of frost at the pole by the end of its winter, the snapshot's season
    /// (MOLA's 1 to 1.5 m, Smith et al. 2001), the southern only beginning; and the pressure
    /// swings by a quarter to a third through the year as the caps take the gas and give it
    /// back (the Viking landers, Hess et al. 1980).
    #[test]
    fn mars_s_winter_poles_hold_their_frost() {
        let (cycle, _) = cycle_of("Mars");
        let f = cycle.frost.expect("Mars lays frost");
        assert_eq!(f.species, "CO2");
        for lat in [-90.0, -80.0, 80.0, 90.0] {
            let coldest = at(&cycle.min, lat, 0.0);
            assert!(
                coldest > f.point - 3.0 && coldest < f.point + 1.0,
                "{lat}: {coldest} against {}",
                f.point
            );
        }
        // At 60° a third of the year or more: at 64° north from Ls 225 to Ls 12, 0.39 of the
        // year (TES, Titus 2005).
        for lat in [-60.0, 60.0] {
            assert!(at(&cycle.frost_share, lat, 0.0) > 0.3, "{lat}");
        }
        for lat in [-30.0, 0.0, 30.0] {
            assert_eq!(at(&cycle.frost_share, lat, 0.0), 0.0, "{lat}");
        }
        let north = at(&cycle.frost_now, 90.0, 0.0);
        let south = at(&cycle.frost_now, -90.0, 0.0);
        assert!((700.0..2000.0).contains(&north), "north {north}");
        assert!(south < 100.0, "south {south}");
        let (lo, hi) = cycle.pressure_range.expect("its pressure swings");
        let swing = (hi - lo) / hi;
        assert!((0.2..0.4).contains(&swing), "{lo} to {hi}");
    }

    /// F239: ground that stays above its gas's frost point lays none: the Earth-like's poles
    /// far above carbon dioxide's 129 K, Titan's uplands above methane's 87 K; Venus's gases
    /// stand past their triple points, and the Moon has none.
    #[test]
    fn no_frost_where_the_ground_stays_above_its_point() {
        for name in ["Earth", "Titan", "Venus", "Moon"] {
            let (cycle, _) = cycle_of(name);
            let most = cycle
                .frost_share
                .as_slice()
                .iter()
                .fold(0.0f32, |a, &v| a.max(v));
            assert_eq!(most, 0.0, "{name}");
            assert_eq!(
                cycle.frost.is_some(),
                matches!(name, "Earth" | "Titan"),
                "{name}"
            );
        }
    }

    /// F239, F247: a thin atmosphere's frost draws on its own gas: Pluto's frost holds no
    /// more than the 1.6 kg/m² its atmosphere weighs, and the pressure follows the frost's
    /// balance through the year. Its nitrogen frost radiating at its phases' emissivities, the
    /// pressure climbs through the year to near the 1 Pa measured (Stern et al. 2015), and
    /// Triton's to near the 1.4 to 1.6 Pa Voyager 2 found (Tyler et al. 1989), where black-body
    /// frost held them at a hundredth and a tenth of it.
    #[test]
    fn a_thin_atmosphere_s_frost_draws_on_its_own_gas() {
        let (cycle, _) = cycle_of("Pluto");
        let f = cycle.frost.expect("Pluto lays frost");
        assert_eq!(f.species, "N2");
        let (mean, ..) = bands(&cycle.frost_now);
        assert!(
            mean > 0.0 && mean < 1.1 * f.column,
            "{mean} against {}",
            f.column
        );
        for (name, measured) in [("Pluto", 1.0), ("Triton", 1.4)] {
            let (cycle, _) = cycle_of(name);
            let (lo, hi) = cycle
                .pressure_range
                .expect("its pressure follows the frost");
            println!("{name}: {lo:.3e} to {hi:.3e} Pa");
            assert!(
                lo < hi && hi > 0.5 * measured && hi < 1.5 * measured,
                "{name}: {lo} to {hi}"
            );
        }
    }

    /// F247: nitrogen frost radiates as its phase does, 0.3 in its alpha phase and 0.8 in
    /// its beta phase above 35.61 K, smoothly across the transition; carbon dioxide frost as a
    /// black body.
    #[test]
    fn nitrogen_frost_radiates_as_its_phase_does() {
        let pluto = Frost::from_graph(&crate::presets::preset("Pluto").unwrap().to_graph().0)
            .expect("Pluto's frost");
        let near = |(e, slope): (f64, f64), want: f64| (e - want).abs() < 1e-12 && slope == 0.0;
        assert!(near(pluto.emissivity(30.0), 0.3) && near(pluto.emissivity(40.0), 0.8));
        let (mid, slope) = pluto.emissivity(N2_PHASE_TRANSITION_K);
        assert!((mid - 0.55).abs() < 1e-12 && slope > 0.0, "{mid} {slope}");
        let mars = Frost::from_graph(&crate::presets::preset("Mars").unwrap().to_graph().0)
            .expect("Mars's frost");
        assert_eq!(mars.emissivity(140.0), (1.0, 0.0));
    }

    /// F239: a slope takes the frost of the place it faces, as its temperature does: on Mars
    /// at the end of the northern winter, ground at 60° north tilted 20° toward the equator
    /// holds a small part of the frost of the flat ground about it.
    #[test]
    fn a_slope_facing_the_sun_holds_less_frost() {
        let (cycle, _) = cycle_of("Mars");
        let mars = downscale_of("Mars");
        let grid = CubeSphere::new(256);
        let rise = 20f64.to_radians().tan() * mars.radius;
        let height = Field::from_fn(grid, |_, d| {
            let (_, lat) = lon_lat_from_direction(d);
            (rise * lat.max(0.0)) as f32
        });
        let frost = mars.frost(&cycle.frost_now, &cycle.mean, &height, f64::NAN);
        let here = at(&frost, 60.0, 0.0);
        let flat = at(&cycle.frost_now, 60.0, 0.0);
        assert!(flat > 100.0 && here < 0.3 * flat, "{here} flat {flat}");
    }

    /// F239: frost hides the ground's colour as it thickens, a millimetre's worth most of it.
    #[test]
    fn frost_hides_the_ground_as_it_thickens() {
        let grid = CubeSphere::new(4);
        let ground = Field::filled(grid, [0.1f32, 0.2, 0.3]);
        let colour = CO2_FROST_COLOUR;
        for (m, share) in [(0.0f32, 0.0), (1.0, 1.0 - (-1.0f64).exp()), (20.0, 1.0)] {
            let a = frosted_albedo(&ground, &Field::filled(grid, m), colour).as_slice()[0];
            for i in 0..3 {
                let expect = 0.1 * (i + 1) as f64 * (1.0 - share) + colour[i] as f64 * share;
                assert!((a[i] as f64 - expect).abs() < 1e-5, "{m}: {a:?}");
            }
        }
    }

    /// Prints the frost of a Mars locked to its star: its night side's cold trap.
    #[test]
    #[ignore = "diagnostic: a locked body's frost"]
    fn print_locked_frost() {
        let mut graph = crate::presets::preset("Mars").unwrap().to_graph().0;
        let au: f64 = std::env::var("EXOGEN_AU").map_or(1.524, |v| v.parse().unwrap());
        graph.set("orbit.semi_major_axis", au * 1.496e11).unwrap();
        let year = graph.f64("orbit.year");
        graph.set("rotation.period", year).unwrap();
        let eb = EnergyBalance::from_graph(&graph);
        println!("resonance {:?}", eb.sky.resonance);
        let climate = eb.solve(CubeSphere::new(CLIMATE_FACE_SIZE));
        let started = std::time::Instant::now();
        let cycle = eb.cycle(&climate, &Ground::from_graph(&graph), None, None);
        let f = cycle.frost.unwrap();
        let (mean, ..) = bands(&cycle.frost_now);
        println!(
            "{:.1} s: pressure {:?}, caps grow {:?} kg/m² a year, frost {mean:.2} kg/m² over the body, column {:.1}",
            started.elapsed().as_secs_f64(),
            cycle.pressure_range,
            cycle.cap_growth,
            f.column
        );
        for lon in [0.0, 90.0, 180.0, -90.0] {
            println!(
                "  lon {lon:6.1}: coldest {:6.1} K, mean {:6.1} K, frost {:4.2} of the year, now {:9.2} kg/m²",
                at(&cycle.min, 0.0, lon),
                at(&cycle.mean, 0.0, lon),
                at(&cycle.frost_share, 0.0, lon),
                at(&cycle.frost_now, 0.0, lon)
            );
        }
    }

    /// Prints the transport's local response, D sum(w) / area, W/m²/K, against the
    /// radiation's slope, at the equator and the pole of each body.
    #[test]
    #[ignore = "diagnostic: the transport's local response"]
    fn print_transport_response() {
        for name in ["Earth", "Mars", "Venus", "Titan", "Pluto", "Triton"] {
            let graph = crate::presets::preset(name).unwrap().to_graph().0;
            let eb = EnergyBalance::from_graph(&graph);
            let grid = CubeSphere::new(CLIMATE_FACE_SIZE);
            let climate = eb.solve(grid);
            let l = links(grid);
            for lat in [0.0f64, 60.0, 89.0] {
                let (c, _, _) = grid.locate(direction_from_lon_lat(0.0, lat.to_radians()));
                let k = grid.index(c);
                let w: f64 = l[k].iter().map(|x| x.1).sum();
                let response = eb.transport * w / grid.cell_area(c);
                let t = climate.temperature.get(c) as f64;
                let (_, slope) = eb.olr_and_slope(t);
                println!(
                    "{name} {lat:4.0}°: response {response:.3} W/m²/K, radiation {slope:.3}, brought {:.2} W/m², steady {t:.1} K",
                    climate.transport.get(c)
                );
            }
        }
    }

    /// Prints each body's cold traps (F240), on ground at the datum: the frost's own gas's
    /// growth and every other ice's year-mean saturation pressure by latitude.
    #[test]
    #[ignore = "diagnostic: the cold traps"]
    fn print_traps() {
        let names =
            std::env::var("EXOGEN_BODY").unwrap_or_else(|_| "Mars,Earth,Pluto,Triton".into());
        for name in names.split(',') {
            let graph = crate::presets::preset(name).unwrap().to_graph().0;
            let eb = EnergyBalance::from_graph(&graph);
            let (cycle, climate) = graph_cycle(&graph);
            let ground = Ground::from_graph(&graph);
            let ices = Ice::of_inventory(&graph, [0.80, 0.82, 0.86]);
            let heights = Field::new(climate.temperature.grid());
            let downscale = Downscale::from_graph(&graph);
            let started = std::time::Instant::now();
            let traps = eb.traps(&climate, &ground, &cycle, &ices, &heights, &downscale);
            println!(
                "{name}: {:?} in {:.1} s",
                ices.iter().map(|i| i.species).collect::<Vec<_>>(),
                started.elapsed().as_secs_f64()
            );
            let lats = [
                -90.0, -85.0, -80.0, -70.0, -60.0, -40.0, 0.0, 40.0, 60.0, 70.0, 80.0, 85.0, 90.0,
            ];
            if let Some((ice, g)) = &traps.growth {
                let row: Vec<String> = lats
                    .iter()
                    .map(|&l| format!("{l}:{:.0}", at(g, l, 0.0)))
                    .collect();
                println!("  {} grows, kg/m² a year: {}", ice.species, row.join(" "));
            }
            for (ice, v) in &traps.vapour {
                let row: Vec<String> = lats
                    .iter()
                    .map(|&l| format!("{l}:{:.2e}", at(v, l, 0.0)))
                    .collect();
                println!("  {} year-mean vapour, Pa: {}", ice.species, row.join(" "));
                println!(
                    "  {} coldest trap {:?}",
                    ice.species,
                    traps.coldest(ice.species)
                );
            }
        }
    }

    /// A body's cold traps (F240) on ground at `height`, m, by longitude.
    fn traps_of(name: &str, height: impl Fn(f64) -> f64 + Sync) -> Traps {
        let graph = crate::presets::preset(name).unwrap().to_graph().0;
        let eb = EnergyBalance::from_graph(&graph);
        let (cycle, climate) = graph_cycle(&graph);
        let ground = Ground::from_graph(&graph);
        let ices = Ice::of_inventory(&graph, [0.80, 0.82, 0.86]);
        let heights = Field::from_fn(climate.temperature.grid(), |_, d| {
            height(lon_lat_from_direction(d).0) as f32
        });
        eb.traps(
            &climate,
            &ground,
            &cycle,
            &ices,
            &heights,
            &Downscale::from_graph(&graph),
        )
    }

    /// F240: each ice finds its coldest trap: Mars's water its north pole, its summer under
    /// the aphelion milder than the south's under the perihelion (the north residual cap's
    /// water ice, Kieffer et al. 1976); Earth's a pole; Pluto's methane its low latitudes,
    /// which its tilt of 120° leaves the least of the year's sunlight (Sputnik Planitia's
    /// 25° north, Hamilton et al. 2016).
    #[test]
    fn each_ice_finds_its_coldest_trap() {
        for (name, species, (lo, hi)) in [
            ("Mars", "H2O", (75.0, 90.0)),
            ("Earth", "H2O", (-90.0, -75.0)),
            ("Pluto", "CH4", (-35.0, 35.0)),
        ] {
            let traps = traps_of(name, |_| 0.0);
            let (lat, _, vapour) = traps.coldest(species).unwrap();
            let lat = if name == "Earth" { -lat.abs() } else { lat };
            assert!((lo..=hi).contains(&lat), "{name} {species} at {lat}");
            assert!(vapour > 0.0, "{name} {species} {vapour}");
        }
    }

    /// F258: the frost's own gas's caps hold its air: Pluto's nitrogen laid in its traps and
    /// its year run again with them as the reservoir its air draws on, its pressure no longer
    /// falling a hundredfold as its winter's frost takes the air, its year's mean the body's;
    /// the caps standing at the frost point, the warmer the lower they lie.
    #[test]
    fn the_caps_hold_the_air() {
        let graph = crate::presets::preset("Pluto").unwrap().to_graph().0;
        let eb = EnergyBalance::from_graph(&graph);
        let grid = CubeSphere::new(CLIMATE_FACE_SIZE);
        let climate = eb.solve(grid);
        let ground = Ground::from_graph(&graph);
        let bare = eb.cycle(&climate, &ground, None, None);
        let ices = Ice::of_inventory(&graph, [0.80, 0.82, 0.86]);
        let downscale = Downscale::from_graph(&graph);
        let traps = eb.traps(
            &climate,
            &ground,
            &bare,
            &ices,
            &Field::new(grid),
            &downscale,
        );
        let mut layers = crate::fields::LayerSet::new(CubeSphere::new(64));
        let laid = crate::stages::caps::lay(&graph, &mut layers, &traps, &[], None);
        let f = eb.frost.expect("Pluto's nitrogen frost");
        let reservoir =
            crate::stages::caps::reservoir(&laid, grid, &eb, f.species).expect("nitrogen caps");
        let over = crate::stages::caps::on_climate(&laid, grid, &eb);
        let held = eb.cycle(
            &eb.solve_over(grid, Some(&over)),
            &ground,
            None,
            Some(&reservoir),
        );
        let (lowest, _) = bare.pressure_range.unwrap();
        let (lo, hi) = held.pressure_range.unwrap();
        let series = held.series.as_ref().unwrap();
        let mean = series.pressure.iter().sum::<f64>() / series.pressure.len() as f64;
        let gain = held
            .reservoir
            .as_ref()
            .map(|r| r.1)
            .expect("the caps' gain");
        println!(
            "Pluto: pressure {lowest:.4} Pa at its lowest without its caps, {lo:.4} to {hi:.4} with them, mean {mean:.4}; the caps gain {gain:.3} kg/m² a year"
        );
        assert!(lo > 5.0 * lowest, "{lo} against {lowest}");
        let body = graph.f64("atmosphere.surface_pressure");
        assert!((mean / body - 1.0).abs() < 0.05, "{mean} against {body}");
        let (high, ..) = held.frost_point_year(1.0).unwrap();
        let (low, ..) = held.frost_point_year(1.2).unwrap();
        assert!(low > high, "{low} lower down, {high} at the datum");
    }

    /// F240: a basin traps better than a height: the higher pressure of its air raises the
    /// frost's point, so that Mars's carbon dioxide ice 5 km down grows more than 5 km up,
    /// and its water there needs less vapour in the air (the trapping that gathers Pluto's
    /// nitrogen in Sputnik Planitia, Bertrand and Forget 2016).
    #[test]
    fn a_basin_traps_better_than_a_height() {
        let traps = traps_of("Mars", |lon| 5000.0 * lon.sin());
        let (_, growth) = traps.growth.as_ref().expect("Mars's CO2 ice");
        let (low, high) = (at(growth, 80.0, -90.0), at(growth, 80.0, 90.0));
        println!("CO2 at 80° north grows {low:.0} kg/m² a year 5 km down, {high:.0} 5 km up");
        assert!(low > high + 10.0, "CO2 {low} down, {high} up");
        let (_, vapour) = traps.vapour.iter().find(|v| v.0.species == "H2O").unwrap();
        let (low, high) = (at(vapour, 80.0, -90.0), at(vapour, 80.0, 90.0));
        println!("H2O at 80° north needs {low:.2e} Pa 5 km down, {high:.2e} 5 km up");
        assert!(low < high, "H2O {low} down, {high} up");
    }

    /// A tidally locked ocean world (F244, F249, F250): the Earth-like as TRAPPIST-1e in the
    /// first habitable case of the THAI comparison (Fauchez et al. 2020; Sergeev et al. 2022),
    /// about a star of 0.0898 solar masses (Van Grootel et al. 2018) at 900 W/m², or
    /// EXOGEN_LOCKED_FLUX's, on a circular orbit, without its moon or its tilt, its spin
    /// locked to its orbit.
    fn locked_ocean_world() -> ParamGraph {
        let mut g = crate::presets::preset("Earth").unwrap().to_graph().0;
        g.set("star.mass", 0.0898).unwrap();
        let luminosity = g.f64("star.luminosity") * L_SUN;
        let flux: f64 = std::env::var("EXOGEN_LOCKED_FLUX").map_or(900.0, |v| v.parse().unwrap());
        g.set(
            "orbit.semi_major_axis",
            (luminosity / (4.0 * PI * flux)).sqrt(),
        )
        .unwrap();
        g.set("orbit.eccentricity", 0.0).unwrap();
        g.set("satellite.mass", 0.0).unwrap();
        g.set("rotation.obliquity", 0.0).unwrap();
        let period = g.f64("orbit.period");
        g.set("rotation.period", period).unwrap();
        g
    }

    /// Prints each preset's warmest place's year-mean temperature as the parameters estimate
    /// it (F249) against the energy balance's warmest cell, and the liquids at each.
    #[test]
    #[ignore = "diagnostic: the warmest place"]
    fn print_warmest() {
        let mut graphs: Vec<(String, ParamGraph)> = crate::presets::PRESET_NAMES
            .iter()
            .map(|n| {
                (
                    n.to_string(),
                    crate::presets::preset(n).unwrap().to_graph().0,
                )
            })
            .collect();
        graphs.push(("locked".to_string(), locked_ocean_world()));
        let locked = &graphs.last().unwrap().1;
        println!(
            "star teff: the Sun's {:.1} K, the locked world's {:.1} K at {:.6} L_sun; ice brightening {:.4} and {:.3}",
            graphs[0].1.f64("star.teff"),
            locked.f64("star.teff"),
            locked.f64("star.luminosity"),
            ice_brightening(graphs[0].1.f64("star.teff")),
            ice_brightening(locked.f64("star.teff"))
        );
        for (name, graph) in &graphs {
            let eb = EnergyBalance::from_graph(graph);
            let climate = eb.solve(CubeSphere::new(CLIMATE_FACE_SIZE));
            let warmest = climate
                .temperature
                .as_slice()
                .iter()
                .fold(f32::MIN, |a, &v| a.max(v));
            println!(
                "{name:10} mean {:6.1} K, warmest estimated {:6.1}, balance {:6.1} ({:.0} W/m²); liquid {} (at the mean {})",
                graph.f64("atmosphere.surface_temperature_mean"),
                graph.f64("atmosphere.surface_temperature_warmest"),
                warmest,
                graph.f64("orbit.insolation"),
                graph.str("volatiles.working_liquid"),
                graph.str("volatiles.mean_liquid"),
            );
        }
    }

    /// A body's sea through its year (F244, F245).
    fn seas_of(graph: &ParamGraph) -> (Seas, Cycle, Climate) {
        let (cycle, climate) = graph_cycle(graph);
        (cycle.seas.clone().expect("a water sea"), cycle, climate)
    }

    /// The locked ocean world's year over a sea at every place, THAI's aquaplanet (F250).
    fn locked_seas(world: &ParamGraph) -> (EnergyBalance, Seas, Climate) {
        let eb = EnergyBalance::from_graph(world);
        let grid = CubeSphere::new(CLIMATE_FACE_SIZE);
        let climate = eb.solve(grid);
        let ocean = Ocean::from_graph(world).expect("a water sea");
        let share = Field::filled(grid, 1.0f32);
        let cycle = eb.cycle(
            &climate,
            &Ground::from_graph(world),
            Some((&ocean, &share)),
            None,
        );
        (eb, cycle.seas.expect("a water sea"), climate)
    }

    /// A locked world's measures as THAI gives them (Sergeev et al. 2022): the body's and its
    /// dayside's mean surface temperatures, K, its planetary albedo, and the share of its
    /// dayside its sea's ice leaves open, the substellar point at longitude 0.
    fn thai_measures(eb: &EnergyBalance, seas: &Seas, climate: &Climate) -> [f64; 4] {
        let g = seas.mean.grid();
        let freezing = seas.ocean.freezing;
        let (mut all, mut day) = ((0.0, 0.0), (0.0, 0.0));
        let (mut lit, mut reflected, mut open) = (0.0, 0.0, 0.0);
        for c in g.cells() {
            let (a, d) = (g.cell_area(c), g.center(c));
            let t = seas.mean_at(d);
            all = (all.0 + t * a, all.1 + a);
            let (lon, lat) = lon_lat_from_direction(d);
            if lat.cos() * lon.cos() <= 0.0 {
                continue;
            }
            day = (day.0 + t * a, day.1 + a);
            let cover = seas.cover_at(d);
            if cover < 0.5 {
                open += a;
            }
            let base = climate.ground_albedo.get(c) as f64;
            let ice = if t >= freezing - SEA_ICE_MELT_MARGIN_K {
                SEA_ICE_MELT_ALBEDO
            } else {
                ICE_ALBEDO
            };
            let albedo = base + cover * (eb.ice_over(ice, base).max(base) - base);
            let q = climate.absorbed.get(c) as f64 / (1.0 - climate.albedo.get(c) as f64);
            lit += q * a;
            reflected += q * a * albedo;
        }
        [all.0 / all.1, day.0 / day.1, reflected / lit, open / day.1]
    }

    /// F250: ice and snow brighten the ground less under a redder star, by the share of its
    /// light in the bands they reflect: within 0.07 of the contrasts Shields et al. 2013 found
    /// over snow for an F2V star, a K2V star and AD Leo against the Sun, and under
    /// TRAPPIST-1's light some 0.3 of the Sun's.
    #[test]
    fn the_ice_dims_under_a_redder_star() {
        let visible = blackbody_share_below(0.7, T_SUN);
        assert!((visible - 0.488).abs() < 0.002, "{visible}");
        assert!(blackbody_share_below(1e3, T_SUN) > 0.999);
        assert!((ice_brightening(T_SUN) - 1.0).abs() < 1e-12);
        for (teff, shields) in [(6750.0, 1.040), (5084.0, 0.923), (3400.0, 0.503)] {
            let b = ice_brightening(teff);
            println!("{teff} K: {b:.3} against Shields et al.'s {shields}");
            assert!(
                (b - shields).abs() < 0.07,
                "{teff} K: {b} against {shields}"
            );
        }
        let trappist = ice_brightening(2566.0);
        assert!((0.25..0.35).contains(&trappist), "{trappist}");
    }

    /// F249, F250: TRAPPIST-1e as THAI's first habitable case, an aquaplanet locked to its
    /// star at 900 W/m², keeps open sea under its star, its working liquid water at its
    /// warmest place though its mean stands below freezing, and its night under ice: its sea
    /// open over 0.3 to 0.7 of its dayside, where THAI's four models leave almost half open
    /// (Sergeev et al. 2022). F251: its ice ends round its star along the sunlight's pattern
    /// between the cycle's places rather than in their steps: where it covers half its sea,
    /// its distance from the substellar point steps by under 0.2° from one bearing to the
    /// next, and ranges less than the places' reading does.
    #[test]
    fn a_locked_ocean_world_keeps_its_sea_under_its_star() {
        let world = locked_ocean_world();
        assert_eq!(world.str("rotation.mode"), "locked");
        assert!(world.f64("atmosphere.surface_temperature_mean") < 273.15);
        assert_eq!(world.str("volatiles.working_liquid"), "water");
        let (eb, seas, climate) = locked_seas(&world);
        let [mean, day, albedo, open] = thai_measures(&eb, &seas, &climate);
        println!(
            "TRAPPIST-1e: mean {mean:.1} K, dayside {day:.1} K, planetary albedo {albedo:.3}, open {open:.2} of the dayside (THAI: 232 to 246 K, 261 to 268 K, 0.22 to 0.28, almost half)"
        );
        let substellar = direction_from_lon_lat(0.0, 0.0);
        let night = direction_from_lon_lat(PI, 0.0);
        assert!(seas.mean_at(substellar) > seas.ocean.freezing + 5.0);
        assert!(seas.ice_at(substellar) < SEA_ICE_COVER_MASS);
        assert!(seas.ice_at(night) > 10.0 * SEA_ICE_COVER_MASS);
        assert!((0.3..0.7).contains(&open), "open {open}");
        // The edge on each bearing from the substellar point, a degree apart: guided, and
        // read between the places as before; its widest step from one bearing to the next and
        // its range.
        let edges = |cover: &dyn Fn([f64; 3]) -> f64| -> (f64, f64) {
            let found: Vec<f64> = (0..360)
                .map(|b| {
                    let (sb, cb) = (b as f64).to_radians().sin_cos();
                    (0..=2400)
                        .map(|k| 30.0 + k as f64 * 0.025)
                        .find(|&deg: &f64| {
                            let (sd, cd) = deg.to_radians().sin_cos();
                            cover([cd, sd * sb, sd * cb]) > 0.5
                        })
                        .unwrap_or(f64::NAN)
                })
                .collect();
            let step = (0..360)
                .map(|b| (found[(b + 1) % 360] - found[b]).abs())
                .fold(0.0, f64::max);
            let lo = found.iter().cloned().fold(f64::MAX, f64::min);
            let hi = found.iter().cloned().fold(f64::MIN, f64::max);
            (step, hi - lo)
        };
        let guided = edges(&|d| seas.cover_at(d));
        let stepped = edges(&|d| seas.read(&seas.cover_rows, d));
        println!(
            "the eye's edge: guided, steps up to {:.2}° between bearings, ranging over {:.2}°; read between the places, {:.2}° and {:.2}°",
            guided.0, guided.1, stepped.0, stepped.1
        );
        assert!(guided.0 < 0.2, "{guided:?}");
        assert!(guided.1 < stepped.1, "{guided:?} {stepped:?}");
    }

    /// F244: the sea holds its seasons: at the Earth-like's 45° its surface swings by under
    /// 10 K where its land swings by tens, and stays open the year round; by the end of the
    /// northern winter its sea ice near the pole stands a metre or a few thick (the Arctic's
    /// perennial ice some 3 m before its thinning, Rothrock et al. 1999).
    #[test]
    fn the_sea_holds_its_seasons_and_freezes_at_the_pole() {
        let earth = crate::presets::preset("Earth").unwrap().to_graph().0;
        let (seas, cycle, _) = seas_of(&earth);
        let swing = at(&seas.max, 45.0, 0.0) - at(&seas.min, 45.0, 0.0);
        let land = at(&cycle.max, 45.0, 0.0) - at(&cycle.min, 45.0, 0.0);
        assert!(swing < 10.0 && land > 40.0, "sea {swing}, land {land}");
        assert_eq!(at(&seas.ice_share, 45.0, 0.0), 0.0);
        let pole = at(&seas.ice_now, 85.0, 0.0);
        assert!((500.0..4000.0).contains(&pole), "ice {pole}");
    }

    /// F245: the air carries the seasons' transport: the Earth-like's mean, its land's and
    /// its sea's, within a kelvin of its steady balance; the air over 80° north swings by some
    /// 30 K through the year, as the Arctic's does from some 243 K in its winter to the
    /// melting point in its summer, and over 45° by some 10 K over a sea that is most of the
    /// latitude; the sea's ice lies all year near the poles and comes and goes at 60°, as
    /// Earth's seasonal ice zone does (Parkinson and Cavalieri 2008).
    #[test]
    fn the_air_carries_the_seasons() {
        let earth = crate::presets::preset("Earth").unwrap().to_graph().0;
        let (cycle, climate) = graph_cycle(&earth);
        let (mean, steady) = (surface_mean(&cycle, &earth), bands(&climate.temperature).0);
        assert!((mean - steady).abs() < 1.0, "{mean} against {steady}");
        let (n_lat, n_lon) = cycle.layout;
        let swing = |lat: f64| {
            // The place at the latitude and the first longitude, its air step by step.
            let i = ((lat + 90.0) / 180.0 * (n_lat - 1) as f64).round() as usize * n_lon;
            let air = &cycle.air[0];
            let places = n_lat * n_lon;
            let (lo, hi) = (0..air.len() / places)
                .map(|k| air[k * places + i])
                .fold((f64::MAX, f64::MIN), |(a, b), t| (a.min(t), b.max(t)));
            hi - lo
        };
        let (arctic, temperate) = (swing(80.0), swing(45.0));
        assert!((20.0..40.0).contains(&arctic), "80°: {arctic} K");
        assert!((5.0..15.0).contains(&temperate), "45°: {temperate} K");
        let seas = cycle.seas.as_ref().expect("a water sea");
        for lat in [-85.0, 85.0] {
            assert_eq!(at(&seas.ice_share, lat, 0.0), 1.0, "{lat}");
        }
        for lat in [-60.0, 60.0] {
            let share = at(&seas.ice_share, lat, 0.0);
            assert!((0.2..0.8).contains(&share), "{lat}: {share}");
        }
        assert_eq!(at(&seas.ice_share, 45.0, 0.0), 0.0);
    }

    /// F246: an ice's window through the year: the longest run of steps it lies through,
    /// round the cycle's end, its phases from the northern spring equinox, opening and
    /// closing where its mass crosses the mass it lies at (F255); all year as 0 and 1; none
    /// where it never lies.
    #[test]
    fn an_ice_s_window_is_its_longest_run() {
        let phase = |t: f64| (t / 10.0 - 0.25).rem_euclid(1.0);
        // Lying from step 8 round to step 2, its most at step 0; a shorter run at 5. None lies
        // before it or after it, so it comes in at 3 a step, as through its next step, its
        // mass crossing 0.5 at time 8 + 5/6, and goes at 3 a step, as through step 2, the 2 it
        // has left crossing 0.5 at time 3.5.
        let masses = [9.0f32, 5.0, 2.0, 0.0, 0.0, 3.0, 0.0, 0.0, 1.0, 4.0];
        let w = Window::of(&masses, 0.5, phase);
        let near = |a: f64, b: f64| (a - b).abs() < 1e-9;
        assert!(
            near(w.freeze, 0.6 + 1.0 / 30.0) && near(w.melt, 0.1),
            "{w:?}"
        );
        assert!(near(w.peak, 9.0) && near(w.at, 0.85), "{w:?}");
        let all = Window::of(
            &[2.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            0.5,
            phase,
        );
        assert!(all.perennial() && near(all.at, 0.85), "{all:?}");
        assert_eq!(Window::of(&[0.1; 10], 0.5, phase), Window::default());
    }

    /// F246: the ices come and go with the seasons: Mars's northern carbon dioxide frost at
    /// 70° north freezes in its autumn and melts out in its spring, its peak near the end of
    /// its winter, its southern frost half a year off (James et al. 1992); the Earth-like's
    /// sea ice at 65° north freezes in its autumn and melts out in its summer, its polar ice
    /// lying all year. F252: toward its edge between the places Mars's frost lies the more
    /// briefly, its window at its edge under half the window 3° in. F253: its year curve is
    /// fitted, and its year step by step holds its pressure's range, the star at the
    /// northern summer solstice's phase near 90°.
    #[test]
    fn the_ices_come_and_go_with_the_seasons() {
        let (cycle, _) = cycle_of("Mars");
        let cryo = cycle.cryosphere.as_ref().expect("Mars's frost");
        let at =
            |lat: f64, sea: bool| cryo.at(direction_from_lon_lat(0.0, f64::to_radians(lat)), sea);
        let (north, south) = (at(70.0, false), at(-70.0, false));
        println!("Mars 70° north {north:?}, south {south:?}");
        assert!(north.peak > 100.0 && !north.perennial(), "{north:?}");
        assert!((0.45..0.85).contains(&north.freeze), "{north:?}");
        assert!((0.05..0.4).contains(&north.melt), "{north:?}");
        let gap = (south.at - north.at).rem_euclid(1.0);
        assert!((0.3..0.7).contains(&gap), "{north:?} {south:?}");
        let length = |w: &Window| (w.melt - w.freeze).rem_euclid(1.0);
        let first = (0..=400)
            .map(|i| 30.0 + i as f64 * 0.1)
            .find(|&lat| at(lat, false).peak > 0.0)
            .expect("the northern frost's edge");
        let (edge, inside) = (at(first, false), at(first + 3.0, false));
        println!("Mars's northern frost's edge at {first:.1}°: {edge:?}; 3° in {inside:?}");
        assert!(length(&edge) < 0.5 * length(&inside), "{edge:?} {inside:?}");
        let curve = cryo.curves[0].expect("Mars's frost's year curve");
        println!(
            "Mars's frost's year curve: rise {:.2}, fall {:.2}",
            curve.0, curve.1
        );
        assert!((0.3..3.0).contains(&curve.0) && (0.5..5.0).contains(&curve.1));
        let series = cycle.series.as_ref().expect("Mars's year");
        let clock = cycle.clock.expect("Mars's clock");
        let (lo, hi) = cycle.pressure_range.expect("Mars's pressure");
        let low = series.pressure.iter().cloned().fold(f64::MAX, f64::min);
        let high = series.pressure.iter().cloned().fold(f64::MIN, f64::max);
        assert!(
            (low - lo).abs() < 1.0 && (high - hi).abs() < 1.0,
            "{low} {high}"
        );
        let summer = (0..series.phase.len())
            .min_by(|&a, &b| {
                let off = |k: usize| {
                    let d = (series.phase[k] - clock.seasons[1]).rem_euclid(1.0);
                    d.min(1.0 - d)
                };
                off(a).total_cmp(&off(b))
            })
            .unwrap();
        println!(
            "Mars's clock {clock:?}; Ls {:.1}° at the summer solstice's phase",
            series.ls_deg[summer]
        );
        assert!((series.ls_deg[summer] - 90.0).abs() < 5.0);
        let (cycle, _) = cycle_of("Earth");
        let cryo = cycle.cryosphere.as_ref().expect("the Earth-like's sea ice");
        let at = |lat: f64| cryo.at(direction_from_lon_lat(0.0, f64::to_radians(lat)), true);
        let (seasonal, polar) = (at(65.0), at(88.0));
        println!("Earth's sea 65° north {seasonal:?}, 88° {polar:?}");
        assert!(
            !seasonal.perennial() && seasonal.peak > 100.0,
            "{seasonal:?}"
        );
        assert!((0.5..0.95).contains(&seasonal.freeze), "{seasonal:?}");
        assert!((0.2..0.6).contains(&seasonal.melt), "{seasonal:?}");
        assert!(polar.perennial(), "{polar:?}");
    }

    /// F255: an ice's window opens and closes where its mass crosses the mass it lies at,
    /// within the step: its masses stand at their steps' ends.
    #[test]
    fn a_window_opens_where_its_mass_crosses() {
        let masses = [0.0f32, 0.5, 1.5, 2.0, 1.0, 0.2, 0.0, 0.0];
        let phase = |t: f64| t / masses.len() as f64;
        let w = Window::of(&masses, 1.0, phase);
        // Up through 1.0 halfway from 0.5 at time 2 to 1.5 at time 3; down to it at time 5.
        assert!((w.freeze - 2.5 / 8.0).abs() < 1e-12, "{w:?}");
        assert!((w.melt - 5.0 / 8.0).abs() < 1e-12, "{w:?}");
        assert!((w.peak - 2.0).abs() < 1e-12 && (w.at - 4.0 / 8.0).abs() < 1e-12);
    }

    /// F255: an ice that runs out within a step, the step's end holding none, goes through
    /// it as it went through the step before: two places losing 10 a step, one left with 2
    /// and one with 8, cross 0.5 at 0.15 and 0.75 of the way into the step after, where a
    /// line to the step's end would take both near it.
    #[test]
    fn an_ice_that_runs_out_within_a_step_goes_as_it_went() {
        let phase = |t: f64| t / 8.0;
        let melt = |left: f32| {
            let masses = [0.0, 0.0, 40.0, left + 20.0, left + 10.0, left, 0.0, 0.0];
            Window::of(&masses, 0.5, phase).melt
        };
        assert!((melt(2.0) - 6.15 / 8.0).abs() < 1e-9, "{}", melt(2.0));
        assert!((melt(8.0) - 6.75 / 8.0).abs() < 1e-9, "{}", melt(8.0));
    }

    /// F254: a fast spinner's cycle runs by longitude: at 60° north the Earth-like's land,
    /// where a continent stands, swings through its year far more widely than its sea beside
    /// it at the same latitude, as Siberia does against the Norwegian Sea.
    #[test]
    fn a_continent_s_year_swings_wider_than_its_sea_s() {
        let earth = crate::presets::preset("Earth").unwrap().to_graph().0;
        let eb = EnergyBalance::from_graph(&earth);
        let grid = CubeSphere::new(CLIMATE_FACE_SIZE);
        let climate = eb.solve(grid);
        let ocean = Ocean::from_graph(&earth).expect("a water sea");
        // Sea over the western half of the longitudes, land over the eastern.
        let share = Field::from_fn(grid, |_, d| {
            let (lon, _) = lon_lat_from_direction(d);
            if lon < 0.0 { 1.0f32 } else { 0.0 }
        });
        let cycle = eb.cycle(
            &climate,
            &Ground::from_graph(&earth),
            Some((&ocean, &share)),
            None,
        );
        let seas = cycle.seas.as_ref().expect("the sea");
        let land = direction_from_lon_lat(90f64.to_radians(), 60f64.to_radians());
        let sea = direction_from_lon_lat(-90f64.to_radians(), 60f64.to_radians());
        let read = |f: &Field<f32>, d: [f64; 3]| f.get(f.grid().locate(d).0) as f64;
        let land_swing = read(&cycle.max, land) - read(&cycle.min, land);
        let sea_swing = seas.max_at(sea) - seas.min_at(sea);
        println!("60° north: the land swings {land_swing:.1} K, the sea {sea_swing:.1} K");
        assert!(land_swing > 2.0 * sea_swing, "{land_swing} {sea_swing}");
    }

    /// F253: an ice read at a phase of the year as an engine reads the cryosphere's maps,
    /// against the vectors the docs give ports to check themselves by.
    #[test]
    fn the_ices_read_as_their_vectors_say() {
        let path = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../docs/cryosphere-vectors.json"
        );
        let text = std::fs::read_to_string(path).expect("the vectors");
        let vectors: serde_json::Value = serde_json::from_str(&text).unwrap();
        let list = vectors["vectors"].as_array().unwrap();
        assert!(list.len() >= 20);
        for v in list {
            let f = |key: &str| v[key].as_f64().unwrap();
            let window = Window {
                freeze: f("freeze"),
                melt: f("melt"),
                peak: f("peak"),
                at: f("at"),
            };
            let got = window.at_phase(f("phase"), (f("rise"), f("fall")));
            assert!((got - f("mass")).abs() < 1e-9, "{v}: {got}");
        }
    }

    /// Writes the cryosphere's vectors (F253): windows that wrap the year's end and do not,
    /// one that lies all year, one that peaks as it freezes, a brief one and none, each read
    /// at phases through its year.
    #[test]
    #[ignore = "writes docs/cryosphere-vectors.json"]
    fn write_cryosphere_vectors() {
        let windows = [
            ((0.6, 0.1, 0.9, 150.0), (0.9, 2.0)),
            ((0.2, 0.6, 0.35, 3.0), (0.7, 2.5)),
            ((0.0, 1.0, 0.5, 7.0), (1.0, 1.0)),
            ((0.3, 0.5, 0.3, 40.0), (1.0, 2.0)),
            ((0.95, 0.02, 0.98, 0.9), (1.0, 1.0)),
            ((0.0, 0.0, 0.0, 0.0), (1.0, 1.0)),
        ];
        let phases = [
            0.0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.35, 0.45, 0.5, 0.6, 0.75, 0.9, 0.97,
        ];
        let vectors: Vec<serde_json::Value> = windows
            .iter()
            .flat_map(|&((freeze, melt, at, peak), (rise, fall))| {
                let window = Window {
                    freeze,
                    melt,
                    peak,
                    at,
                };
                phases.iter().map(move |&phase| {
                    serde_json::json!({
                        "freeze": freeze, "melt": melt, "at": at, "peak": peak,
                        "rise": rise, "fall": fall, "phase": phase,
                        "mass": window.at_phase(phase, (rise, fall)),
                    })
                })
            })
            .collect();
        let doc = serde_json::json!({
            "about": "ExoGen's cryosphere maps read at a phase of the year (docs/adr/0006): an ice's window from freeze round to melt, its peak's phase at, its peak (any unit) and its curve's exponents rise and fall; mass is what the reading gives at the phase, in the peak's unit.",
            "vectors": vectors,
        });
        let path = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../docs/cryosphere-vectors.json"
        );
        std::fs::write(path, serde_json::to_string_pretty(&doc).unwrap()).unwrap();
    }

    /// Prints the Earth-like's sea ice through its year near its north pole, kg/m², every
    /// fourth step from perihelion, over its sea at every latitude (EXOGEN_SEA_SHARE).
    #[test]
    #[ignore = "diagnostic: the sea ice's year"]
    fn print_sea_series() {
        let earth = crate::presets::preset("Earth").unwrap().to_graph().0;
        let eb = EnergyBalance::from_graph(&earth);
        let climate = eb.solve(CubeSphere::new(CLIMATE_FACE_SIZE));
        let (ocean, share) = sea_of(&earth).unwrap();
        let layout = eb.layout(&climate, Some(&share), None, None);
        let ground = Ground::from_graph(&earth);
        let year = eb.seasonal_years(&ground, &layout, eb.sky.obliquity, Some(&ocean), None);
        let steps = year.masses[1].len() / layout.places.len();
        for lat in [89.0f64, 85.0, 80.0, 75.0, 70.0] {
            let i = ((lat + 90.0) / 2.0).round() as usize;
            let row: Vec<String> = (0..steps)
                .step_by(4)
                .map(|k| format!("{:.0}", year.masses[1][i * steps + k]))
                .collect();
            println!("{lat}: {}", row.join(" "));
        }
    }

    /// Prints the Earth-like's sea by latitude against its land, and the locked ocean
    /// world's along its equator from its substellar point.
    #[test]
    #[ignore = "diagnostic: the seas"]
    fn print_seas() {
        let earth = crate::presets::preset("Earth").unwrap().to_graph().0;
        let (seas, cycle, _) = seas_of(&earth);
        let (_, climate) = graph_cycle(&earth);
        println!(
            "Earth: the body's mean {:.1} K, its steady {:.1} K",
            surface_mean(&cycle, &earth),
            bands(&climate.temperature).0
        );
        for lat in [
            -90.0, -80.0, -75.0, -70.0, -65.0, -60.0, -55.0, -45.0, -30.0, 0.0, 30.0, 45.0, 55.0,
            60.0, 65.0, 70.0, 75.0, 80.0, 90.0,
        ] {
            println!(
                "Earth {lat:5.0}°: sea {:6.1} K ({:6.1} to {:6.1}), ice now {:8.1} kg/m², {:4.2} of the year; land {:6.1} ({:6.1} to {:6.1}); steady {:6.1}",
                at(&seas.mean, lat, 0.0),
                at(&seas.min, lat, 0.0),
                at(&seas.max, lat, 0.0),
                at(&seas.ice_now, lat, 0.0),
                at(&seas.ice_share, lat, 0.0),
                at(&cycle.mean, lat, 0.0),
                at(&cycle.min, lat, 0.0),
                at(&cycle.max, lat, 0.0),
                at(&climate.temperature, lat, 0.0)
            );
        }
        let (n_lat, n_lon) = cycle.layout;
        for lat in [-60.0f64, 0.0, 45.0, 60.0, 80.0] {
            let i = ((lat + 90.0) / 180.0 * (n_lat - 1) as f64).round() as usize * n_lon;
            let air = &cycle.air[0];
            let steps = air.len() / (n_lat * n_lon);
            let row: Vec<String> = (0..steps)
                .step_by(steps / 12)
                .map(|k| format!("{:.0}", air[k * n_lat * n_lon + i]))
                .collect();
            println!("Earth {lat:5.0}° air through the year: {}", row.join(" "));
        }
        // The sea ice's reach through the snapshot's season over the body, by its cover.
        let g = seas.ice_now.grid();
        let share = sea_of(&earth).map_or(0.0, |s| s.1.as_slice()[0] as f64);
        let (mut north, mut south) = (0.0, 0.0);
        for c in g.cells() {
            let m = seas.ice_now.get(c) as f64;
            if m > SEA_ICE_COVER_MASS * std::f64::consts::LN_2 {
                let (_, lat) = lon_lat_from_direction(g.center(c));
                let a = g.cell_area(c) * share / (4.0 * PI);
                if lat > 0.0 {
                    north += a;
                } else {
                    south += a;
                }
            }
        }
        println!(
            "Earth: sea ice at the snapshot over {:.2}% of the body, {:.2}% north, {:.2}% south",
            100.0 * (north + south),
            100.0 * north,
            100.0 * south
        );
        let world = locked_ocean_world();
        println!(
            "locked world: star {:.4} L_sun, insolation {:.1} W/m², mode {}, period {:.2} d, liquid {}, mean {:.1} K",
            world.f64("star.luminosity"),
            world.f64("orbit.insolation"),
            world.str("rotation.mode"),
            world.f64("orbit.period") / 86400.0,
            world.str("volatiles.working_liquid"),
            world.f64("atmosphere.surface_temperature_mean")
        );
        if Ocean::from_graph(&world).is_none() {
            let started = std::time::Instant::now();
            let (cycle, climate) = graph_cycle(&world);
            println!(
                "locked world: no sea at its zero-dimensional mean; its year under its air in {:.1} s, mean {:.1} K (steady {:.1})",
                started.elapsed().as_secs_f64(),
                bands(&cycle.mean).0,
                bands(&climate.temperature).0
            );
            for lon in [0.0, 60.0, 90.0, 120.0, 180.0] {
                println!(
                    "locked {lon:5.0}°: steady {:6.1} K, mean {:6.1} K",
                    at(&climate.temperature, 0.0, lon),
                    at(&cycle.mean, 0.0, lon)
                );
            }
            return;
        }
        let started = std::time::Instant::now();
        let (eb, seas, climate) = locked_seas(&world);
        let [mean, day, albedo, open] = thai_measures(&eb, &seas, &climate);
        println!(
            "locked world in {:.1} s: ice brightening {:.3}; mean {mean:.1} K, dayside {day:.1} K, planetary albedo {albedo:.3}, open {open:.2} of the dayside",
            started.elapsed().as_secs_f64(),
            eb.brightening
        );
        for lon in [0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 150.0, 180.0] {
            println!(
                "locked {lon:5.0}° from the substellar point: steady {:6.1} K, sea {:6.1} K, ice {:9.1} kg/m²",
                at(&climate.temperature, 0.0, lon),
                at(&seas.mean, 0.0, lon),
                at(&seas.ice_now, 0.0, lon)
            );
        }
        for lat in [30.0, 60.0, 89.0] {
            println!(
                "locked lat {lat}° at the substellar longitude: sea {:6.1} K, ice {:9.1}",
                at(&seas.mean, lat, 0.0),
                at(&seas.ice_now, lat, 0.0)
            );
        }
        // The eye's shape: the sea at one distance from the substellar point along its
        // meridian, its equator and the diagonal between, read from the cycle's places.
        for distance in [50.0f64, 55.0, 60.0, 65.0, 70.0] {
            let (cos, diagonal) = (distance.to_radians().cos(), distance.to_radians().cos());
            let x = diagonal.sqrt().acos().to_degrees();
            let read = |lat: f64, lon: f64| {
                let d = direction_from_lon_lat(lon.to_radians(), lat.to_radians());
                (seas.mean_at(d), seas.cover_at(d))
            };
            let (meridian, equator, across) =
                (read(distance, 0.0), read(0.0, distance), read(x, x));
            println!(
                "locked {distance}° from the substellar point (cos {cos:.3}): meridian {:.1} K cover {:.2}, equator {:.1} K cover {:.2}, diagonal {:.1} K cover {:.2}",
                meridian.0, meridian.1, equator.0, equator.1, across.0, across.1
            );
        }
    }

    /// F245: the zonal mean of a system the same along each latitude is the system: the
    /// preconditioned conjugate gradients solve it in a step or two.
    #[test]
    fn a_uniform_air_is_its_own_zonal_mean() {
        let mut graph = crate::presets::preset("Mars").unwrap().to_graph().0;
        let year = graph.f64("orbit.year");
        graph.set("rotation.period", year).unwrap();
        let eb = EnergyBalance::from_graph(&graph);
        let climate = eb.solve(CubeSphere::new(CLIMATE_FACE_SIZE));
        let layout = eb.layout(&climate, None, None, None);
        let n = layout.places.len();
        let links: Vec<Vec<(usize, f64)>> = layout
            .links
            .iter()
            .map(|l| l.iter().map(|&(q, g)| (q, 0.2 * g)).collect())
            .collect();
        let diag = vec![0.5; n];
        let rhs: Vec<f64> = (0..n)
            .map(|p| 200.0 + 50.0 * (p as f64 * 0.37).sin())
            .collect();
        let mut solutions = Vec::new();
        for zonal in [true, false] {
            let mut x = vec![200.0; n];
            let steps = transport_solve(&layout, &links, (&diag, &rhs), &mut x, zonal);
            let worst = (0..n)
                .map(|p| {
                    let total = diag[p] + links[p].iter().map(|l| l.1).sum::<f64>();
                    let off: f64 = links[p].iter().map(|&(q, k)| k * x[q]).sum();
                    (total * x[p] - off - rhs[p]).abs()
                })
                .fold(0.0, f64::max);
            println!("zonal {zonal}: {steps} steps, worst residual {worst:.3e}");
            assert!(worst < 1e-5, "zonal {zonal}: residual {worst}");
            assert!(!zonal || steps <= 2, "{steps} steps");
            solutions.push(x);
        }
        let apart = solutions[0]
            .iter()
            .zip(&solutions[1])
            .fold(0.0f64, |a, (x, y)| a.max((x - y).abs()));
        assert!(apart < 1e-6, "{apart} K apart");
    }

    fn downscale_of(name: &str) -> Downscale {
        Downscale::from_graph(&crate::presets::preset(name).unwrap().to_graph().0)
    }

    /// F243: a crater at the Moon's pole keeps its floor from the sun: its rim stands over
    /// the 1.8° the sun's disc climbs there at most, and its floor, seeing only its walls,
    /// stands near the 25 to 70 K Diviner finds on the floors of the Moon's polar craters
    /// (Paige et al. 2010); the same crater at 60° south sees the sun at noon.
    #[test]
    fn a_polar_crater_s_floor_never_sees_the_sun() {
        let (cycle, climate) = cycle_of("Moon");
        let moon = downscale_of("Moon");
        let graph = crate::presets::preset("Moon").unwrap().to_graph().0;
        let shade = Shade::from_graph(&graph, &EnergyBalance::from_graph(&graph));
        let grid = CubeSphere::new(256);
        let pole = direction_from_lon_lat(0.0, -PI / 2.0);
        let sixty = direction_from_lon_lat(0.0, -PI / 3.0);
        let (radius, depth) = (30_000.0, 3_000.0);
        let bowl = |d: [f64; 3], centre: [f64; 3]| {
            let r = crate::geom::angle_between(d, centre) * moon.radius / radius;
            if r < 1.0 { -depth * (1.0 - r * r) } else { 0.0 }
        };
        let height = Field::from_fn(grid, |_, d| (bowl(d, pole) + bowl(d, sixty)) as f32);
        let mut t = moon.apply(&cycle.mean, &height, f64::NAN);
        let share = moon.shade(&mut t, &height, &climate.absorbed, &shade);
        let floor = at(&t, -90.0, 0.0);
        let lit = at(&t, -60.0, 0.0);
        println!(
            "shadowed {share:.2e} of the body, the polar floor {floor:.1} K, at 60° {lit:.1} K"
        );
        assert!(share > 0.0, "no shadow");
        assert!((25.0..75.0).contains(&floor), "polar floor {floor}");
        assert!(lit > 150.0, "60° floor {lit}");
    }

    /// F237: a place's height cools it by the air's share of the lapse rate and the
    /// radiation's share of the thinner greenhouse: an Earth-like's 3 km summit some 15 K
    /// under its sea level, 5 K/km against the air's 6.3; Venus's 11 km summit near 650 K, as
    /// Maxwell Montes stands (Pettengill et al. 1992); Mars's thin air cools its heights by
    /// less than half a kelvin a kilometre; an airless body's not at all.
    #[test]
    fn a_summit_stands_under_less_greenhouse() {
        let earth = downscale_of("Earth");
        let t = earth.at(288.0, 288.0, 3000.0);
        assert!((270.0..278.0).contains(&t), "Earth {t}");
        let t = downscale_of("Venus").at(755.0, 755.0, 11_000.0);
        assert!((630.0..680.0).contains(&t), "Venus {t}");
        let t = downscale_of("Mars").at(210.0, 210.0, 20_000.0);
        assert!(t > 200.0 && t < 210.0, "Mars {t}");
        assert_eq!(downscale_of("Moon").at(250.0, 250.0, 5000.0), 250.0);
    }

    /// F237: a slope takes the sunlight of the place it faces: on the airless Moon ground at
    /// 60° north tilted 20° toward the equator stands at the temperature of 40° north; Titan's
    /// thick air holds a slope at its surroundings'.
    #[test]
    fn a_slope_faces_the_sunlight_of_another_latitude() {
        let (cycle, _) = cycle_of("Moon");
        let moon = downscale_of("Moon");
        let grid = CubeSphere::new(256);
        let rise = 20f64.to_radians().tan() * moon.radius;
        let height = Field::from_fn(grid, |_, d| {
            let (_, lat) = lon_lat_from_direction(d);
            (rise * lat.max(0.0)) as f32
        });
        let t = moon.apply(&cycle.mean, &height, f64::NAN);
        let here = at(&t, 60.0, 0.0);
        let flat = at(&cycle.mean, 60.0, 0.0);
        let faced = at(&cycle.mean, 40.0, 0.0);
        assert!(
            (here - faced).abs() < 3.0 && here > flat + 5.0,
            "{here} flat {flat} faced {faced}"
        );
        assert!((downscale_of("Titan").at(91.0, 95.0, 0.0) - 91.0).abs() < 0.1);
    }

    /// F235, the gate of 5A: an Earth-like at 1 AU around a G2V star within 5 K of Earth's
    /// global mean of 288 K, its tropics near Earth's 299 K and its poles below freezing.
    #[test]
    fn an_earth_like_world_is_within_five_kelvin_of_earth() {
        let (climate, _) = climate_of("Earth");
        let (mean, equator, poles) = bands(&climate.temperature);
        assert!((mean - 288.0).abs() < 5.0, "mean {mean}");
        assert!((equator - 299.0).abs() < 6.0, "equator {equator}");
        assert!(poles < 268.0 && poles > 230.0, "poles {poles}");
    }

    /// F235: Mars, Venus and Titan at their measured means, Venus and Titan nearly the same
    /// everywhere under their thick, slowly turning atmospheres (the Pioneer Venus probes a few
    /// kelvin apart from 60°N to 30°S near the ground, Seiff et al. 1980; Titan's poles 3 K
    /// under its equator, Jennings et al. 2016).
    #[test]
    fn mars_venus_and_titan_take_their_measured_means() {
        let (mars, _) = climate_of("Mars");
        let (mean, equator, poles) = bands(&mars.temperature);
        assert!((mean - 210.0).abs() < 8.0, "Mars mean {mean}");
        assert!(
            equator - poles > 25.0,
            "Mars equator {equator} poles {poles}"
        );
        let (venus, _) = climate_of("Venus");
        let (mean, equator, poles) = bands(&venus.temperature);
        assert!((mean / 737.0 - 1.0).abs() < 0.04, "Venus mean {mean}");
        assert!(
            equator - poles < 10.0,
            "Venus equator {equator} poles {poles}"
        );
        let (titan, _) = climate_of("Titan");
        let (mean, equator, poles) = bands(&titan.temperature);
        assert!((mean - 93.7).abs() < 4.0, "Titan mean {mean}");
        assert!(
            equator - poles < 5.0,
            "Titan equator {equator} poles {poles}"
        );
    }

    /// The zonal mean of a wind's year at a latitude, deg: its eastward and northward wind.
    fn zonal_wind(wind: &crate::wind::Wind, lat: f64) -> (f64, f64) {
        let (n_lat, n_lon) = wind.layout;
        let i = ((lat + 90.0) / 180.0 * (n_lat - 1) as f64).round() as usize;
        let (east, north, _) = wind.year_mean();
        let row = i * n_lon..(i + 1) * n_lon;
        (
            east[row.clone()].iter().sum::<f64>() / n_lon as f64,
            north[row].iter().sum::<f64>() / n_lon as f64,
        )
    }

    /// F267: the Earth-like's winds through its year blow in Earth's belts: easterly trades
    /// toward the equator within some 30°, westerlies between some 35 and 55°, easterlies
    /// poleward of them, as Earth's polar easterlies blow at some 60 to 70°; at some 3 to 9
    /// m/s over the body, its year's mean wind without the weather's gusts, whose speeds
    /// over Earth's seas average some 7 m/s (Wentz et al. 2007); the thermal equator
    /// following the sun to either side of the equator through the year.
    #[test]
    fn the_earth_like_s_winds_blow_in_its_belts() {
        let earth = crate::presets::preset("Earth").unwrap().to_graph().0;
        let (cycle, _) = graph_cycle(&earth);
        let wind = cycle.wind.as_ref().expect("a wind under the air");
        for side in [1.0, -1.0] {
            let (trades, toward) = zonal_wind(wind, 14.0 * side);
            let (westerlies, away) = zonal_wind(wind, 46.0 * side);
            let (polar, _) = zonal_wind(wind, 66.0 * side);
            assert!(
                trades < -1.0 && toward * side < 0.0,
                "{side}: {trades}, {toward}"
            );
            assert!(
                westerlies > 1.0 && away * side > 0.0,
                "{side}: {westerlies}, {away}"
            );
            assert!(polar < 0.0, "{side}: {polar}");
        }
        let brief = wind.summary();
        println!("Earth-like: {brief:?}");
        assert!((3.0..9.0).contains(&brief.speed), "{} m/s", brief.speed);
        assert!(
            brief.equator.0 < -2.0 && brief.equator.1 > 2.0,
            "{:?}",
            brief.equator
        );
    }

    /// F268: the Earth-like's water through its year: some 1,000 kg/m² a year falls over the
    /// body, as Earth's 1,000 to 1,100 do (Trenberth et al. 2007), most of it where its air
    /// converges, at its thermal equator and under its westerlies' storms, its subtropics
    /// evaporating more than falls on them; its land sheds some 300 a year as runoff, Earth's
    /// some 300 to 400.
    #[test]
    fn the_earth_like_s_water_falls_where_its_air_converges() {
        let earth = crate::presets::preset("Earth").unwrap().to_graph().0;
        let eb = EnergyBalance::from_graph(&earth);
        let (cycle, _) = graph_cycle(&earth);
        let air = crate::moisture::Air::from_graph(&earth, eb.transport, OCEAN_TRANSPORT_SHARE)
            .expect("water");
        let heights = vec![0.0; cycle.sites.len()];
        let water = cycle.moisture(&air, &heights).expect("a moisture");
        let brief = water.summary();
        println!("Earth-like: {brief:?}");
        assert!((700.0..1500.0).contains(&brief.precipitation), "{brief:?}");
        assert!((150.0..600.0).contains(&brief.runoff), "{brief:?}");
        let (n_lat, n_lon) = water.layout;
        let year = SECONDS_PER_DAY * 365.25;
        let zonal = |v: &[f64], lat: f64| {
            let i = ((lat + 90.0) / 180.0 * (n_lat - 1) as f64).round() as usize;
            v[i * n_lon..(i + 1) * n_lon].iter().sum::<f64>() / n_lon as f64 * year
        };
        let (p, e) = (&water.precipitation, &water.evaporation);
        for side in [1.0, -1.0] {
            let (rising, sinking, storms) = (0.0, 26.0 * side, 48.0 * side);
            assert!(zonal(p, rising) > 3.0 * zonal(p, sinking), "{side}");
            assert!(zonal(p, storms) > 2.0 * zonal(p, sinking), "{side}");
            assert!(zonal(e, sinking) > zonal(p, sinking), "{side}");
            assert!(zonal(p, storms) > zonal(e, storms), "{side}");
        }
        assert!(zonal(p, 0.0) > zonal(e, 0.0));
    }

    /// F269: the year again at each place's own wind and its wet ground's latent heat: the
    /// Earth-like's mean stays within 5 K of Earth's 288 K, 5A's gate, and its water near
    /// Earth's some 1,000 kg/m² a year; Titan's ground, its exchange with its air at its
    /// own slow wind, stands farther from its air, its contrast from its equator to its poles
    /// some four times its stated wind's.
    #[test]
    fn the_year_runs_again_at_its_own_wind() {
        let mut contrasts = Vec::new();
        for name in ["Earth", "Titan"] {
            let graph = crate::presets::preset(name).unwrap().to_graph().0;
            let eb = EnergyBalance::from_graph(&graph);
            let ground = Ground::from_graph(&graph);
            let climate = eb.solve(CubeSphere::new(CLIMATE_FACE_SIZE));
            let ocean = sea_of(&graph);
            let sea = ocean.as_ref().map(|(o, f)| (o, f));
            let water = if sea.is_some() {
                OCEAN_TRANSPORT_SHARE
            } else {
                0.0
            };
            let air = crate::moisture::Air::from_graph(&graph, eb.transport, water);
            let first = eb.cycle(&climate, &ground, sea, None);
            let heights = vec![0.0; first.sites.len()];
            let wet = air.and_then(|a| first.moisture(&a, &heights));
            let coupling = first.coupling(wet.as_ref()).expect("a coupling");
            let second = eb.cycle_with(&climate, &ground, sea, None, Some(&coupling));
            let again = air.and_then(|a| second.moisture(&a, &heights));
            let ((mean, equator, poles), (_, e1, p1)) = (bands(&second.mean), bands(&first.mean));
            println!(
                "{name}: mean {mean:.2} K, contrast {:.2} K where {:.2} at the stated wind; {:?}",
                equator - poles,
                e1 - p1,
                again.as_ref().map(|m| m.summary())
            );
            contrasts.push((equator - poles, e1 - p1));
            if name == "Earth" {
                assert!((mean - 288.0).abs() < 5.0, "{mean}");
                let rain = again.expect("water").summary().precipitation;
                assert!((800.0..1400.0).contains(&rain), "{rain}");
            }
        }
        let (titan, stated) = contrasts[1];
        assert!(titan > 3.0 * stated, "{titan} {stated}");
    }

    /// F269: exploration, THAI's measures of the locked ocean world, coupled.
    #[test]
    #[ignore]
    fn explore_thai_coupled() {
        let world = locked_ocean_world();
        let eb = EnergyBalance::from_graph(&world);
        let grid = CubeSphere::new(CLIMATE_FACE_SIZE);
        let climate = eb.solve(grid);
        let ocean = Ocean::from_graph(&world).expect("a water sea");
        let share = Field::filled(grid, 1.0f32);
        let sea = Some((&ocean, &share));
        let ground = Ground::from_graph(&world);
        let first = eb.cycle(&climate, &ground, sea, None);
        let air = crate::moisture::Air::from_graph(&world, eb.transport, OCEAN_TRANSPORT_SHARE);
        let heights = vec![0.0; first.sites.len()];
        let wet = air.and_then(|a| first.moisture(&a, &heights));
        let coupling = first.coupling(wet.as_ref()).expect("a coupling");
        let second = eb.cycle_with(&climate, &ground, sea, None, Some(&coupling));
        for (label, cycle) in [("uncoupled", &first), ("coupled", &second)] {
            let seas = cycle.seas.clone().expect("a water sea");
            let [mean, day, albedo, open] = thai_measures(&eb, &seas, &climate);
            println!(
                "TRAPPIST-1e {label}: mean {mean:.1} K, dayside {day:.1} K, planetary albedo {albedo:.3}, open {open:.2} of the dayside (THAI: 232 to 246 K, 261 to 268 K, 0.22 to 0.28, almost half)"
            );
        }
    }

    /// F275: frost darker than its ground, gaining enough under the sun to sublime all it
    /// holds early in a step, leaves the ground the rest of the step: where the bare ground's
    /// own gain would cool it past the point it frosts again, the step ending with thin frost
    /// at the point rather than bare; where its gain keeps it warm the step ends bare past it.
    #[test]
    fn frost_that_runs_out_leaves_its_ground_the_rest_of_its_step() {
        let graph = crate::presets::preset("Mars").unwrap().to_graph().0;
        let ground = Ground::from_graph(&graph);
        let column = Column::new(&ground, 88_775.0, 5.94e7);
        let n = column.len();
        let dt = 5.94e7 / CYCLE_YEAR_STEPS as f64;
        let (point, latent) = (148.0, 5.9e5);
        let shed = |t: f64| (0.95 * SIGMA_SB * t.powi(4), 3.8 * SIGMA_SB * t.powi(3));
        for (bare, frosts) in [(10.0, true), (200.0, false)] {
            let start = vec![point; n];
            let mut t = start.clone();
            let (mut w, mut tr) = (Work::new(n), Trial::new(n));
            let mut mass = 5.0;
            column.trial(
                &t,
                (&mut w, &mut tr),
                dt,
                (bare, 150.0),
                0.0,
                shed,
                mass,
                point,
            );
            let lasted = column.finish(&mut t, &tr, dt, point, shed(point), latent, &mut mass);
            let lasts = lasted.expect("the frost runs out within the step");
            let ends = (point, shed(point).0, latent);
            let rows = (&mut w, &mut tr);
            column.rest_bare(
                &mut t,
                &start,
                rows,
                (dt, lasts),
                (bare, 0.0),
                shed,
                ends,
                &mut mass,
            );
            println!(
                "bare gain {bare} W/m²: lasted {lasts:.3} of the step, ends {mass:.2} kg/m² at {:.2} K",
                t[0]
            );
            assert!(lasts > 0.0 && lasts < 0.2, "{lasts}");
            if frosts {
                assert!(mass > 0.0 && (t[0] - point).abs() < 1e-6, "{mass} {}", t[0]);
            } else {
                assert!(mass == 0.0 && t[0] > point, "{mass} {}", t[0]);
            }
        }
    }

    /// F271: snow under a warm step holds the ground at its melting point and melts by the
    /// step's warmth over its latent heat of fusion; a thin snow melts away and the ground
    /// warms past the point.
    #[test]
    fn snow_melts_at_its_point() {
        let graph = crate::presets::preset("Earth").unwrap().to_graph().0;
        let ground = Ground::from_graph(&graph);
        let column = Column::new(&ground, 86_400.0, 3.156e7);
        let n = column.len();
        let dt = 3.156e7 / CYCLE_YEAR_STEPS as f64;
        let shed = |t: f64| (0.6 * SIGMA_SB * t.powi(4), 2.4 * SIGMA_SB * t.powi(3));
        let (point, latent) = (273.16, 3.34e5);
        let mut melted = Vec::new();
        for lying in [500.0, 0.01] {
            let start = vec![272.0; n];
            let mut t = start.clone();
            let (mut w, mut tr) = (Work::new(n), Trial::new(n));
            column.step(&mut t, &mut w, dt, 400.0, 0.06, shed);
            assert!(t[0] > point, "{}", t[0]);
            let mut snow = lying;
            let rows = (&mut w, &mut tr);
            let gains = (400.0, 0.06);
            column.thaw(
                &mut t,
                &start,
                rows,
                dt,
                gains,
                shed(point),
                (point, latent),
                &mut snow,
            );
            println!("snow {lying} -> {snow:.3} kg/m², the ground {:.2} K", t[0]);
            melted.push((lying - snow, t[0]));
            if lying > 1.0 {
                assert!(
                    (t[0] - point).abs() < 1e-9 && snow < lying && snow > 0.0,
                    "{snow}"
                );
            } else {
                assert!(snow == 0.0 && t[0] > point, "{snow} {}", t[0]);
            }
        }
        // The step's warmth past the point: what the ground would gain over what it sheds and
        // takes below, some 400 - 0.6 σ T⁴ W/m² through the step, over the latent heat.
        let most = (400.0 - shed(point).0) * dt / latent;
        assert!(
            melted[0].0 > 0.2 * most && melted[0].0 < 1.2 * most,
            "{:?} {most}",
            melted
        );
    }

    /// F271: the Earth-like's year at its own wind lays snow on its northern land through its
    /// winter and melts it through its summer: Earth's northern snow covers some 47 million
    /// km² in January, nearly half of its land, and some 2 million in August (Robinson et al.
    /// 1993). F277: its snow lying until it melts, no place's land sheds more than falls on it.
    #[test]
    fn the_earth_like_s_snow_comes_and_goes() {
        let graph = crate::presets::preset("Earth").unwrap().to_graph().0;
        let eb = EnergyBalance::from_graph(&graph);
        let ground = Ground::from_graph(&graph);
        let climate = eb.solve(CubeSphere::new(CLIMATE_FACE_SIZE));
        let ocean = sea_of(&graph);
        let sea = ocean.as_ref().map(|(o, f)| (o, f));
        let air = crate::moisture::Air::from_graph(&graph, eb.transport, OCEAN_TRANSPORT_SHARE);
        let first = eb.cycle(&climate, &ground, sea, None);
        let heights = vec![0.0; first.sites.len()];
        let wet = air.and_then(|a| first.moisture(&a, &heights));
        let coupling = first.coupling(wet.as_ref()).expect("a coupling");
        let year = eb.cycle_with(&climate, &ground, sea, None, Some(&coupling));
        let n = year.sites.len();
        let steps = year.masses[2].len() / n.max(1);
        assert!(steps > 0, "no snow's year");
        // The share of the northern land the snow covers at each step.
        let covered: Vec<f64> = (0..steps)
            .map(|k| {
                let (mut sum, mut land) = (0.0, 0.0);
                for (p, site) in year.sites.iter().enumerate().filter(|(_, s)| s.lat > 0.0) {
                    let l = site.weight * (1.0 - site.sea - site.ice).max(0.0);
                    sum += l * snow_cover(year.masses[2][p * steps + k] as f64);
                    land += l;
                }
                sum / land.max(1e-30)
            })
            .collect();
        let (most, least) = covered
            .iter()
            .fold((0.0f64, 1.0f64), |(a, b), &c| (a.max(c), b.min(c)));
        let (mean, _, _) = bands(&year.mean);
        println!(
            "the northern land's snow covers from {least:.3} to {most:.3}; the mean {mean:.2} K"
        );
        assert!(most > 0.2 && least < 0.6 * most, "{least} {most}");
        assert!((mean - 288.0).abs() < 5.0, "{mean}");
        // F277: under the year's snowpack the land's bucket takes its rain and its snow's melt,
        // its snow lying until it melts: no place's land sheds more than falls on it.
        let water = air
            .and_then(|a| year.moisture(&a, &heights))
            .expect("the year's water");
        let most_shed = (0..n)
            .filter(|&p| water.land[p] > 0.0 && water.precipitation[p] > 0.0)
            .map(|p| water.runoff[p] / water.precipitation[p])
            .fold(0.0f64, f64::max);
        println!("the land's runoff at most {most_shed:.3} of its precipitation");
        assert!(most_shed <= 1.02, "{most_shed}");
    }

    /// F268: the locked ocean world's water falls under its star, where its air converges,
    /// and hardly at all on its night side, which its air leaves cold and dry.
    #[test]
    fn a_locked_world_s_water_falls_under_its_star() {
        let world = locked_ocean_world();
        let eb = EnergyBalance::from_graph(&world);
        let grid = CubeSphere::new(CLIMATE_FACE_SIZE);
        let climate = eb.solve(grid);
        let ocean = Ocean::from_graph(&world).expect("a water sea");
        let share = Field::filled(grid, 1.0f32);
        let sea = Some((&ocean, &share));
        let cycle = eb.cycle(&climate, &Ground::from_graph(&world), sea, None);
        let air = crate::moisture::Air::from_graph(&world, eb.transport, OCEAN_TRANSPORT_SHARE)
            .expect("water");
        let water = cycle
            .moisture(&air, &vec![0.0; cycle.sites.len()])
            .expect("a moisture");
        let (n_lat, n_lon) = water.layout;
        let year = SECONDS_PER_DAY * 365.25;
        let at = |lon: f64| {
            let (i, j) = (
                n_lat / 2,
                ((lon + 180.0) / 360.0 * n_lon as f64).round() as usize,
            );
            water.precipitation[i * n_lon + j % n_lon] * year
        };
        let brief = water.summary();
        println!(
            "TRAPPIST-1e: {brief:?}; {:.0} kg/m² a year under its star, {:.0} opposite",
            at(0.0),
            at(180.0)
        );
        assert!(at(0.0) > 20.0 * at(180.0), "{} {}", at(0.0), at(180.0));
        assert!((200.0..3000.0).contains(&brief.precipitation), "{brief:?}");
    }

    /// F267: the locked ocean world's air blows toward the point under its star: along its
    /// equator the wind runs toward the substellar longitude from either side, its night
    /// side's air flowing out to its day.
    #[test]
    fn a_locked_world_s_air_blows_toward_its_star() {
        let world = locked_ocean_world();
        let eb = EnergyBalance::from_graph(&world);
        let grid = CubeSphere::new(CLIMATE_FACE_SIZE);
        let climate = eb.solve(grid);
        let ocean = Ocean::from_graph(&world).expect("a water sea");
        let share = Field::filled(grid, 1.0f32);
        let sea = Some((&ocean, &share));
        let cycle = eb.cycle(&climate, &Ground::from_graph(&world), sea, None);
        let wind = cycle.wind.as_ref().expect("a wind under the air");
        let (n_lat, n_lon) = wind.layout;
        let (east, ..) = wind.year_mean();
        let row = (n_lat / 2) * n_lon;
        let at = |lon: f64| {
            let j = ((lon + 180.0) / 360.0 * n_lon as f64).round() as usize % n_lon;
            east[row + j]
        };
        let brief = wind.summary();
        println!(
            "TRAPPIST-1e: {brief:?}; at its equator {:.2} m/s at 40° W, {:.2} at 40° E",
            at(-40.0),
            at(40.0)
        );
        assert!(at(-40.0) > 0.5 && at(40.0) < -0.5);
        assert!((1.0..15.0).contains(&brief.speed), "{} m/s", brief.speed);
    }

    /// F235: an airless body's cells balance alone, each at the temperature that radiates
    /// the year's sunlight it absorbs.
    #[test]
    fn an_airless_cell_balances_its_own_sunlight() {
        let (moon, eb) = climate_of("Moon");
        assert_eq!(eb.transport, 0.0);
        let g = moon.temperature.grid();
        for c in g.cells().step_by(97) {
            let a = moon.absorbed.get(c) as f64;
            let t = moon.temperature.get(c) as f64;
            let expect = ((a + eb.heat_flux) / SIGMA_SB).powf(0.25);
            assert!((t - expect).abs() < 0.05, "{t} against {expect}");
        }
    }

    /// F284: the snapshot's phase is the moment its maps were taken: the end of the step that
    /// holds its season, as the windows' phases are their steps' ends, after its season's
    /// own phase by less than a step, on a cycle of one orbit or of several.
    #[test]
    fn the_snapshot_s_phase_is_its_step_s_end() {
        let steps = 360;
        for (season, orbits, equinox) in [(0.0, 1, 0.0), (0.37, 1, 0.12), (0.5, 3, 0.05)] {
            let anomaly = TAU * season;
            let phase = snapshot_phase(anomaly, steps, orbits, equinox);
            let k = snapshot_step(anomaly, steps / orbits);
            let end = ((k as f64 + 1.0) / steps as f64 - equinox).rem_euclid(1.0);
            assert!((phase - end).abs() < 1e-12, "{phase} {end}");
            let asked = (anomaly / TAU / orbits as f64 - equinox).rem_euclid(1.0);
            let after = (phase - asked).rem_euclid(1.0);
            assert!(
                after > 0.0 && after <= 1.0 / steps as f64 + 1e-12,
                "{after}"
            );
        }
    }
}
