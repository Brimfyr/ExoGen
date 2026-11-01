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
use crate::insolation::{PathStep, Sky, diurnal_share, read_grid, read_row};
use crate::params::ParamGraph;
use crate::params::derive::is_icy;
use crate::params::derive::physics::{
    blackbody_share_below, mean_heat_capacity, mean_molecular_mass, regolith_heat,
    saturation_temperature, species_columns, tau_dry, vapour_pressure, volatile, water_column,
    water_optical_depth,
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
    /// The frost's albedo and colour.
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
                    albedo: (colour.iter().sum::<f32>() / 3.0) as f64,
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
    /// The surface's mean, the ice at the snapshot's season and the share of the sea it
    /// covers at each place of the cycle's layout, and the layout's longitudes, a fast
    /// spinner's one.
    mean_rows: Vec<f64>,
    ice_rows: Vec<f64>,
    cover_rows: Vec<f64>,
    n_lon: usize,
}

impl Seas {
    /// The sea's ice at the snapshot's season in a direction, kg/m², read along the cycle's
    /// latitudes, a resonant body's longitudes too.
    pub fn ice_at(&self, d: [f64; 3]) -> f64 {
        self.read(&self.ice_rows, d)
    }

    /// The sea's surface's mean through the year in a direction, K.
    pub fn mean_at(&self, d: [f64; 3]) -> f64 {
        self.read(&self.mean_rows, d)
    }

    /// The share of the sea its ice covers at the snapshot's season in a direction (F245):
    /// each place's, 1 - exp(-m / SEA_ICE_COVER_MASS), read between the places, so that the
    /// ice's margin spans the distance between a place under ice and an open one, as a
    /// marginal ice zone does, where its mass would step from one to the other.
    pub fn cover_at(&self, d: [f64; 3]) -> f64 {
        self.read(&self.cover_rows, d)
    }

    fn read(&self, rows: &[f64], d: [f64; 3]) -> f64 {
        let (lon, lat) = lon_lat_from_direction(d);
        match self.n_lon {
            1 => read_row(rows, lat),
            n => read_grid(rows, n, lat, lon),
        }
    }
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

    pub fn albedo(&self) -> f64 {
        (self.colour.iter().sum::<f32>() / 3.0) as f64
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

    /// The window an ice lies in through a cycle of `masses`, one at the end of each step
    /// (F246), a mass past `cover` counting as lying, as half the ground's colour hidden:
    /// the longest run of steps it lies through, which it freezes in the middle of its first
    /// and melts in the middle of the step after its last; its peak the most it holds in
    /// that run. `phase` takes a time in steps from the cycle's start to
    /// the phase of the year.
    fn of(masses: &[f32], cover: f64, phase: impl Fn(f64) -> f64) -> Window {
        let n = masses.len();
        let lying: Vec<bool> = masses.iter().map(|&m| m as f64 > cover).collect();
        let count = lying.iter().filter(|&&l| l).count();
        let peak_in = |steps: &mut dyn Iterator<Item = usize>| {
            steps
                .map(|k| (k, masses[k] as f64))
                .max_by(|a, b| a.1.total_cmp(&b.1))
                .unwrap_or((0, 0.0))
        };
        if count == 0 {
            return Window::default();
        }
        if count == n {
            let (k, peak) = peak_in(&mut (0..n));
            return Window {
                freeze: 0.0,
                melt: 1.0,
                peak,
                at: phase(k as f64 + 1.0),
            };
        }
        // The longest run, from a step it does not lie in round the cycle.
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
        let (k, peak) = peak_in(&mut (0..longest).map(|j| (best + j) % n));
        Window {
            freeze: phase(best as f64 + 0.5),
            melt: phase((best + longest) as f64 + 0.5),
            peak,
            at: phase(k as f64 + 1.0),
        }
    }
}

/// The ices that come and go through the year (F246): at each place of the cycle's layout
/// its land's frost's window and its sea's ice's, and the phase of the snapshot's season.
#[derive(Clone, Debug)]
pub struct Cryosphere {
    land: Vec<Window>,
    sea: Vec<Window>,
    layout: (usize, usize),
    /// The phase of the year of the snapshot's season.
    pub snapshot: f64,
}

impl Cryosphere {
    /// Whether an ice, its land's frost or, `sea`, its sea's ice, lies anywhere in the year.
    pub fn holds(&self, sea: bool) -> bool {
        let windows = if sea { &self.sea } else { &self.land };
        windows.iter().any(|w| w.peak > 0.0)
    }

    /// The window of an ice in a direction (F246), its land's frost's or, `sea`, its sea's
    /// ice's: its peak read between the cycle's places, as its mass is, and its phases those
    /// of the places it lies at, by their weights, on the circle of the year; all year where
    /// the nearest of them holds it all year.
    pub fn at(&self, d: [f64; 3], sea: bool) -> Window {
        let (lon, lat) = lon_lat_from_direction(d);
        let windows = if sea { &self.sea } else { &self.land };
        let corners = reader(self.layout.0, self.layout.1, lat, lon);
        let peak: f64 = corners.iter().map(|&(p, w)| w * windows[p].peak).sum();
        let lying: Vec<(f64, Window)> = corners
            .iter()
            .filter(|&&(p, w)| w > 0.0 && windows[p].peak > 0.0)
            .map(|&(p, w)| (w, windows[p]))
            .collect();
        let Some(&(_, nearest)) = lying.iter().max_by(|a, b| a.0.total_cmp(&b.0)) else {
            return Window::default();
        };
        let circular = |among: &[(f64, Window)], pick: fn(&Window) -> f64| {
            let (s, c) = among.iter().fold((0.0, 0.0), |(s, c), (w, window)| {
                let angle = TAU * pick(window);
                (s + w * angle.sin(), c + w * angle.cos())
            });
            (s.atan2(c) / TAU).rem_euclid(1.0)
        };
        let (freeze, melt) = if nearest.perennial() {
            (0.0, 1.0)
        } else {
            let seasonal: Vec<(f64, Window)> = lying
                .iter()
                .copied()
                .filter(|(_, w)| !w.perennial())
                .collect();
            (
                circular(&seasonal, |w| w.freeze),
                circular(&seasonal, |w| w.melt),
            )
        };
        Window {
            freeze,
            melt,
            peak,
            at: circular(&lying, |w| w.at),
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
    /// Each place's frost and its sea's ice at the end of each step, kg/m², place after
    /// place (F246).
    masses: [Vec<f32>; 2],
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
    /// Each place's frost and its sea's ice at the end of each step, kg/m², place after
    /// place (F246).
    masses: [Vec<f32>; 2],
}

impl Year {
    /// A year without frost or sea.
    fn bare(mean: f64, max: f64, min: f64) -> Self {
        Self {
            mean,
            max,
            min,
            frost_share: 0.0,
            frost_now: 0.0,
            sea: [mean, min, max, 0.0, 0.0],
        }
    }
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
        let evaporation = match volatile("H2O") {
            Some(w) if graph.str("volatiles.working_liquid") == "water" && ventilation > 0.0 => {
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
        }
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

    /// The heat the air exchanges with a sea of water per kelvin (F245), W/m²/K: the
    /// ground's, and its evaporation's, the latent heat of the vapour the sea gives the air
    /// for each kelvin it stands above it, the air held near saturation, so that its vapour's
    /// pressure climbs with its temperature as the sea's does (Clausius and Clapeyron).
    fn sea_coupling(&self, t: f64) -> f64 {
        let latent = match volatile("H2O") {
            Some(w) if self.evaporation > 0.0 => {
                let l = if t >= w.t_triple { w.l_vap } else { w.l_sub };
                l * self.evaporation * vapour_pressure(w, t) * l / (w.r_specific * t * t)
            }
            _ => 0.0,
        };
        self.coupling(t) + latent
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

    /// The year's mean climate on a grid over caps of ice (F241): `caps` the share of each
    /// cell they cover and their albedo there, which each cell reflects at by its share.
    pub fn solve_over(&self, grid: CubeSphere, caps: Option<&(Field<f32>, Field<f32>)>) -> Climate {
        let sun = self.sky.annual_sunlight(grid);
        let n = grid.cell_count();
        let area: Vec<f64> = (0..n).map(|k| grid.cell_area(grid.cell(k))).collect();
        let links = links(grid);
        let flux: Vec<f64> = sun.flux.as_slice().iter().map(|&q| q as f64).collect();
        let height: Vec<f64> = sun.height.as_slice().iter().map(|&m| m as f64).collect();
        // The albedo's low-sun term about the sunlight-weighted mean height, so that without
        // ice the planet reflects its Bond albedo.
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
    /// grid, its sea (F245).
    pub fn cycle(
        &self,
        climate: &Climate,
        ground: &Ground,
        sea: Option<(&Ocean, &Field<f32>)>,
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
        let (n_lat, n_lon) = match self.sky.resonance {
            None => (CYCLE_LATITUDES, 1),
            Some(_) => CYCLE_RESONANT_GRID,
        };
        let paths: Vec<Vec<PathStep>> = match self.sky.resonance {
            Some(_) => tilts.iter().map(|&eps| self.resonant_steps(eps)).collect(),
            None => Vec::new(),
        };
        let ocean = sea.map(|s| s.0);
        let each: Vec<TiltYear> = if self.exchange > 0.0 {
            let layout = self.layout(climate, sea.map(|s| s.1));
            tilts
                .iter()
                .enumerate()
                .map(|(j, &eps)| match self.sky.resonance {
                    None => self.seasonal_years(ground, &layout, eps, ocean),
                    Some((_, orbits)) => {
                        self.seasonal_orbits(ground, &layout, (orbits, &paths[j]), ocean)
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
                    masses: [Vec::new(), Vec::new()],
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
        // The ices' windows: each place's masses through its cycle, the tilts' mean step by
        // step, the phase of the year counted from the northern spring equinox.
        let places = n_lat * n_lon;
        let steps = each
            .first()
            .map_or(0, |e| e.masses[0].len() / places.max(1));
        let cryosphere = (steps > 0).then(|| {
            let orbits = self.sky.resonance.map_or(1, |(_, q)| q) as f64;
            let e = self.sky.eccentricity;
            let equinox = snapshot_anomaly(0.0, self.sky.perihelion_season, e) / TAU / orbits;
            let phase = |t: f64| (t / steps as f64 - equinox).rem_euclid(1.0);
            let windows = |j: usize, cover: f64| -> Vec<Window> {
                (0..places)
                    .map(|p| {
                        let span = p * steps..(p + 1) * steps;
                        let masses: Vec<f32> = (0..steps)
                            .map(|k| {
                                each.iter()
                                    .map(|t| t.masses[j][span.start + k])
                                    .sum::<f32>()
                                    / each.len() as f32
                            })
                            .collect();
                        Window::of(&masses, cover, phase)
                    })
                    .collect()
            };
            Cryosphere {
                land: windows(0, FROST_COVER_MASS * std::f64::consts::LN_2),
                sea: windows(1, SEA_ICE_COVER_MASS * std::f64::consts::LN_2),
                layout: (n_lat, n_lon),
                snapshot: (self.snapshot / TAU / orbits - equinox).rem_euclid(1.0),
            }
        });
        let seas = ocean.map(|ocean| Seas {
            mean: field(|y| y.sea[0]),
            min: field(|y| y.sea[1]),
            max: field(|y| y.sea[2]),
            ice_now: field(|y| y.sea[3]),
            ice_share: field(|y| y.sea[4]),
            ocean: *ocean,
            mean_rows: rows.iter().map(|y| y.sea[0]).collect(),
            ice_rows: rows.iter().map(|y| y.sea[3]).collect(),
            cover_rows: rows
                .iter()
                .map(|y| 1.0 - (-y.sea[3].max(0.0) / SEA_ICE_COVER_MASS).exp())
                .collect(),
            n_lon,
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
            tilts,
            partials: each.iter().map(|e| e.partials.clone()).collect(),
            seas,
            cryosphere,
            air: each.into_iter().map(|e| e.air).collect(),
            layout: (n_lat, n_lon),
        }
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
    fn layout(&self, climate: &Climate, sea: Option<&Field<f32>>) -> Layout {
        let grid = climate.temperature.grid();
        let (n_lat, n_lon) = match self.sky.resonance {
            None => (CYCLE_LATITUDES, 1),
            Some(_) => CYCLE_RESONANT_GRID,
        };
        let (dphi, dlon) = (PI / (n_lat - 1) as f64, TAU / n_lon as f64);
        let lat_of = |i: usize| -PI / 2.0 + dphi * i as f64;
        let edges = |i: usize| {
            let lat = lat_of(i);
            (
                (lat - dphi / 2.0).max(-PI / 2.0),
                (lat + dphi / 2.0).min(PI / 2.0),
            )
        };
        let share_at = |lat: f64, lon: f64| {
            sea.map_or(0.0, |f| {
                f.get(f.grid().locate(direction_from_lon_lat(lon, lat)).0) as f64
            })
        };
        let places = (0..n_lat * n_lon)
            .map(|k| {
                let (i, j) = (k / n_lon, k % n_lon);
                let (lat, lon) = (lat_of(i), -PI + dlon * j as f64);
                let (lo, hi) = edges(i);
                let (c, _, _) = grid.locate(direction_from_lon_lat(lon, lat));
                let sea = match n_lon {
                    1 => {
                        (0..SEA_SHARE_LONGITUDES)
                            .map(|m| share_at(lat, TAU * m as f64 / SEA_SHARE_LONGITUDES as f64))
                            .sum::<f64>()
                            / SEA_SHARE_LONGITUDES as f64
                    }
                    _ => share_at(lat, lon),
                };
                Place {
                    lat,
                    lon,
                    weight: (hi.sin() - lo.sin()) / n_lon as f64,
                    steady: climate.temperature.get(c) as f64,
                    albedo: climate.ground_albedo.get(c) as f64,
                    snow: climate.snow_albedo.get(c) as f64,
                    sea: sea.clamp(0.0, 1.0),
                }
            })
            .collect();
        let links = (0..n_lat * n_lon)
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
            .collect();
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
        Layout {
            places,
            links,
            n_lon,
            modes: [
                by_mode(f64::cos),
                by_mode(f64::sin),
                by_longitude(f64::cos),
                by_longitude(f64::sin),
            ],
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
        let lockstep = self.lockstep(ground, &column, layout, timing, light, ocean, every);
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
                for record in &r.records {
                    let f = frost.map(|f| Frost {
                        point: record.point,
                        ..*f
                    });
                    let (flux, decl) = sun[record.k];
                    let (sd, cd) = decl.sin_cos();
                    let shed = self.shed(record.air);
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
                                    let sheds = (shed, self.frost_shed(record.air, f));
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
                            }
                        }
                    }
                    mean += sum / CYCLE_DAY_STEPS as f64;
                }
                let mean = mean / r.records.len().max(1) as f64;
                Year {
                    mean,
                    max: hi,
                    min: lo,
                    frost_share: r.frosted as f64 / steps as f64,
                    frost_now: r.now,
                    sea: r.sea_year(steps, (mean, hi, lo)),
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
        let lockstep = self.lockstep(ground, &column, layout, timing, light, ocean, 0);
        let years = lockstep
            .runs
            .iter()
            .map(|r| {
                let mean = r.sums[0] / steps as f64;
                Year {
                    mean,
                    max: r.hi,
                    min: r.lo,
                    frost_share: r.frosted as f64 / steps as f64,
                    frost_now: r.now,
                    sea: r.sea_year(steps, (mean, r.hi, r.lo)),
                }
            })
            .collect();
        TiltYear {
            years,
            range: lockstep.range,
            growth: lockstep.growth.map(|g| g / orbits.max(1) as f64),
            partials: lockstep.partials,
            air: lockstep.air,
            masses: lockstep.masses,
        }
    }

    /// Places run in step through their cycles under one atmosphere (F239, F245).
    ///
    /// The air over each place holds its column's heat (F236) and exchanges heat with the
    /// ground under it, its land's and its sea's by their shares of the place (`coupling`,
    /// `sea_coupling`), and with the air beside it by the transport, which the cycle answers
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
    /// The cycle is `steps` steps of `dt` with the snapshot's step, each place lit at a step
    /// as `light` gives, the ground recorded every `every` steps (never at 0).
    #[allow(clippy::too_many_arguments)]
    fn lockstep(
        &self,
        ground: &Ground,
        column: &Column,
        layout: &Layout,
        (steps, dt, snapshot): (usize, f64, usize),
        light: impl Fn(usize, &Place) -> f64 + Sync,
        ocean: Option<&Ocean>,
        every: usize,
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
            .map(|p| Run::new((column.len(), steps), p, ocean))
            .collect();
        let mut air: Vec<f64> = places.iter().map(|p| p.steady).collect();
        let mut series = vec![0.0; steps * n];
        let (mut diag, mut rhs) = (vec![0.0; n], vec![0.0; n]);
        let hold = ground.air / dt;
        // The transport through the air, less the share a sea's water carries across the
        // sea's share of each boundary.
        let water = ocean.map_or(0.0, |_| OCEAN_TRANSPORT_SHARE);
        let air_links: Vec<Vec<(usize, f64)>> = layout
            .links
            .iter()
            .enumerate()
            .map(|(p, links)| {
                links
                    .iter()
                    .map(|&(q, g)| {
                        let wet = places[p].sea.min(places[q].sea);
                        (q, self.transport * g * (1.0 - water * wet))
                    })
                    .collect()
            })
            .collect();
        let mut inflow = vec![0.0; n];
        let (mut held, mut load) = (frost.map_or(0.0, |f| f.column), 0.0);
        let mut range = (f64::MAX, f64::MIN);
        let mut growth = 0.0;
        let mut partials = Vec::with_capacity(steps);
        let mut last = vec![f64::NAN; 2 * n];
        for _ in 0..CYCLE_MAX_YEARS {
            range = (f64::MAX, f64::MIN);
            partials.clear();
            runs.iter_mut().for_each(Run::begin);
            for k in 0..steps {
                // The air, from the ground under it as the last step left it.
                for (q, (r, p)) in runs.iter().zip(places).enumerate() {
                    let land = self.coupling(r.t[0]) * (1.0 - p.sea);
                    let (wet, surface) = match &r.sea {
                        Some(s) => (self.sea_coupling(s.surface) * p.sea, s.surface),
                        None => (0.0, 0.0),
                    };
                    diag[q] = hold + land + wet;
                    rhs[q] = hold * air[q] + land * r.t[0] + wet * surface;
                }
                transport_solve(layout, &air_links, (&diag, &rhs), &mut air, true);
                series[k * n..(k + 1) * n].copy_from_slice(&air);
                if let Some(ocean) = ocean.filter(|_| water > 0.0) {
                    let step = (self.transport * water, dt);
                    sea_transport(layout, &mut runs, step, ocean, &mut inflow);
                }
                let highest = frost.map(|f| f.point_at(partial(held)));
                runs.par_iter_mut()
                    .zip(places.par_iter())
                    .zip(air.par_iter().zip(inflow.par_iter()))
                    .for_each(|((r, p), (&a, &brought))| {
                        let lit = light(k, p);
                        r.albedo = self.land_albedo((p.albedo, p.snow), r.t[0]);
                        let bare = (1.0 - r.albedo) * lit;
                        let shed = self.shed(a);
                        let heat = self.heat_flux;
                        match (frost, highest) {
                            (Some(f), Some(highest)) => {
                                let gains = (bare, (1.0 - f.albedo) * lit);
                                let rows = (&mut r.work, &mut r.trial);
                                column.trial(&r.t, rows, dt, gains, heat, shed, r.mass, highest);
                            }
                            _ => column.step(&mut r.t, &mut r.work, dt, bare, heat, shed),
                        }
                        if let (Some(sea), Some(ocean)) = (r.sea.as_mut(), ocean) {
                            let exchange = self.sea_coupling(sea.surface);
                            let water = (ocean, brought);
                            self.sea_step(sea, dt, (lit, p.albedo), (a, exchange), water);
                        }
                    });
                let point = match (frost, highest) {
                    (Some(f), Some(highest)) => {
                        let gas = (held - load, whole);
                        let point =
                            self.frost_point(column, (&runs, places, &air), gas, dt, f, highest);
                        let ((olr, slope), exchange) =
                            (self.olr_and_slope(point), self.coupling(point));
                        let (e, de) = f.emissivity(point);
                        runs.par_iter_mut().zip(air.par_iter()).for_each(|(r, &a)| {
                            let shed = (
                                e * olr + exchange * (point - a),
                                e * slope + de * olr + exchange,
                            );
                            column.finish(
                                &mut r.t,
                                &r.trial,
                                dt,
                                point,
                                shed,
                                f.latent,
                                &mut r.mass,
                            );
                        });
                        load = over(&runs, |r| r.mass);
                        let now = partial(held - load);
                        range = (range.0.min(now), range.1.max(now));
                        partials.push(now);
                        point
                    }
                    _ => f64::NAN,
                };
                runs.par_iter_mut()
                    .zip(air.par_iter())
                    .for_each(|(r, &a)| r.record(k, snapshot, every, point, a));
            }
            // The frost that lay all cycle and grew through it: a cap's growth.
            growth = over(&runs, Run::growth);
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
            if let Some(f) = frost {
                let mean_load = over(&runs, |r| r.frost_sum) / steps as f64;
                let kept = ((held - mean_load) / column_gas).clamp(0.0, 1.0);
                for r in runs.iter_mut() {
                    r.mass -= kept * r.growth();
                }
                load = over(&runs, |r| r.mass);
                held = (f.column + (mean_load - growth) * kept).max(load);
            }
            for r in runs.iter_mut() {
                let sums = std::mem::take(&mut r.sums);
                settle(&mut r.t, &sums, steps);
                r.sums = sums;
            }
        }
        // The surface pressure: the other gases' and the frost's gas's.
        let range = frost.map(|f| {
            let others = (self.pressure - f.partial).max(0.0);
            (others + range.0, others + range.1)
        });
        let masses = [
            runs.iter()
                .flat_map(|r| r.masses[0].iter().copied())
                .collect(),
            runs.iter()
                .flat_map(|r| r.masses[1].iter().copied())
                .collect(),
        ];
        Lockstep {
            runs,
            range,
            growth: frost.map(|_| growth),
            partials,
            air: series,
            masses,
        }
    }

    /// The frost point a step of the lockstep settles at (F239): the point at which the frost
    /// the places take through the step leaves the gas, `gas` kg/m² over the body before it,
    /// at the partial pressure whose point it is, each place's ground shedding its heat to
    /// space and to the air over it (F245). The frost taken rises with the point and the
    /// point of what is left falls, so the two cross once, between none and `highest`, the
    /// point of all the gas; found by false position, halving the side that stays, the low
    /// side, which leaves the air the gas the frost does not take.
    fn frost_point(
        &self,
        column: &Column,
        (runs, places, air): (&[Run], &[Place], &[f64]),
        (gas, whole): (f64, f64),
        dt: f64,
        frost: &Frost,
        highest: f64,
    ) -> f64 {
        let taking: Vec<(&Run, f64, f64)> = runs
            .iter()
            .zip(places)
            .zip(air)
            .filter(|((r, _), _)| r.mass > 0.0 || !r.trial.clear)
            .map(|((r, p), &a)| (r, p.weight / whole, a))
            .collect();
        let gap = |point: f64| {
            let (olr, exchange) = (self.olr(point), self.coupling(point));
            let olr = olr * frost.emissivity(point).0;
            let taken: f64 = taking
                .iter()
                .map(|&(r, w, a)| {
                    let radiated = olr + exchange * (point - a);
                    column
                        .deposit(&r.trial, dt, point, radiated, frost.latent, r.mass)
                        .1
                        * w
                })
                .sum();
            point - frost.point_at(frost.partial * (gas - taken).max(0.0) / frost.column.max(1e-30))
        };
        let (mut lo, mut hi) = (0.0, highest);
        let (mut g_lo, mut g_hi) = (gap(lo), gap(hi));
        if g_hi <= 0.0 {
            return hi;
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
        // A cell's air is the cycle's over its place, as much colder the higher the cell
        // stands as its ground is (F237): near its ground the air follows the ground, a thin
        // atmosphere's barely cooling with height; its height the one its pressure's share
        // gives.
        let (n_lat, n_lon) = cycle.layout;
        let cell_at = |lat: f64, lon: f64, thin: f64| {
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
        };
        let frost = cycle.frost.as_ref();
        // A fast spinner's traps are its latitudes' at a few heights, the frost's point the
        // only thing a height changes, read at each cell between them; a resonant body's are
        // each cell's own.
        let levels: Option<Vec<f64>> = match self.sky.resonance {
            None => {
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
            Some(_) => None,
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
            // Each place's value through the tilts' mean, and its warmest at any.
            let runs: Vec<(f64, f64)> = places
                .par_iter()
                .map(|cell| {
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
                })
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
                    }
                }
            }
            mean += sum / CYCLE_DAY_STEPS as f64;
        }
        Year::bare(mean / seasons.len().max(1) as f64, hi, lo)
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
        for _ in 0..CYCLE_MAX_YEARS {
            let mut sums = vec![0.0; t.len()];
            (hi, lo) = (f64::MIN, f64::MAX);
            for s in path {
                let cz = (sp * s.sin_decl + cp * s.cos_decl * (s.hour + lon).cos()).max(0.0);
                let gain = (1.0 - albedo) * scale * s.flux * cz + brought;
                column.step(&mut t, &mut work, dt, gain, self.heat_flux, olr);
                for (sum, v) in sums.iter_mut().zip(&t) {
                    *sum += v;
                }
                hi = hi.max(t[0]);
                lo = lo.min(t[0]);
            }
            mean = sums[0] / path.len() as f64;
            if (mean - last).abs() < CYCLE_SETTLED_K {
                break;
            }
            last = mean;
            settle(&mut t, &sums, path.len());
        }
        Year::bare(mean, hi, lo)
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
    /// thickness and its peak's phase, a field apiece.
    pub fn seasons(
        &self,
        cryosphere: &Cryosphere,
        mean: &Field<f32>,
        (height, liquid): (&Field<f32>, &Field<f32>),
        densities: (f64, f64),
    ) -> [Field<f32>; 4] {
        let src = mean.grid();
        let warmth = PaddedField::new(&src, mean.as_slice());
        let grid = height.grid();
        let windows: Vec<Window> = (0..grid.cell_count())
            .into_par_iter()
            .map(|k| {
                let c = grid.cell(k);
                let d = grid.center(c);
                if liquid.get(c) > 0.0 {
                    let w = cryosphere.at(d, true);
                    return Window {
                        peak: w.peak / densities.1,
                        ..w
                    };
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
                if peak <= 0.0 || read.peak <= 0.0 {
                    return Window::default();
                }
                Window {
                    peak: peak / densities.0,
                    ..read
                }
            })
            .collect();
        let field = |pick: fn(&Window) -> f64| {
            Field::from_fn(grid, |c, _| pick(&windows[grid.index(c)]) as f32)
        };
        [
            field(|w| w.freeze),
            field(|w| w.melt),
            field(|w| w.peak),
            field(|w| w.at),
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
fn read_table(values: &[f64], levels: &[f64], lat: f64, thin: f64, log: bool) -> f64 {
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
}

impl Work {
    fn new(n: usize) -> Self {
        Self {
            a: vec![0.0; n],
            b: vec![0.0; n],
            c: vec![0.0; n],
            d: vec![0.0; n],
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

/// A place of the cycle's layout (F239, F245): its latitude and longitude, rad, its share of
/// the sphere's area, its year's mean temperature, its ground's albedo without ice and snow,
/// and the share of it the sea covers.
#[derive(Clone, Copy, Debug)]
struct Place {
    lat: f64,
    lon: f64,
    weight: f64,
    steady: f64,
    albedo: f64,
    snow: f64,
    sea: f64,
}

/// The places the cycle runs at under an atmosphere (F245): each place's links to its
/// neighbours, the transport's conductance to each per unit of the place's area at unit
/// transport, W/m²/K; and the layout's longitudes, a fast spinner's one.
struct Layout {
    places: Vec<Place>,
    links: Vec<Vec<(usize, f64)>>,
    n_lon: usize,
    /// A resonant body's longitudes' Fourier modes: the cos and sin of each mode's phase at
    /// each longitude, mode after mode, and longitude after longitude.
    modes: [Vec<f64>; 4],
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
    /// The frost and the sea's ice at the end of each step of the cycle, kg/m² (F246).
    masses: [Vec<f32>; 2],
}

impl Run {
    fn new((n, steps): (usize, usize), place: &Place, ocean: Option<&Ocean>) -> Self {
        let sea = ocean.filter(|_| place.sea > 0.0).map(|o| {
            let t = place.steady.max(o.freezing);
            Sea {
                t,
                ice: 0.0,
                surface: t,
            }
        });
        Self {
            t: vec![place.steady; n],
            mass: 0.0,
            work: Work::new(n),
            trial: Trial::new(n),
            sea,
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
            masses: [vec![0.0; steps], vec![0.0; steps]],
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
        self.frost_sum += self.mass;
        self.frost_least = self.frost_least.min(self.mass);
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
        self.trial(t, (w, tr), dt, gains, bottom, bare, *mass, f.point);
        self.finish(t, tr, dt, f.point, frosted(f.point), f.latent, mass);
    }

    /// A step's trial under the frost (F239). While frost lies on the ground, or the ground
    /// would cool below the frost point, the ground stands at the point and the step's energy
    /// sublimes or condenses frost instead, at its latent heat (Leighton and Murray 1966);
    /// frost reflects at its own albedo, so the gains are the bare ground's and the frosted
    /// ground's. The trial takes bare ground's step left to itself and, unless that leaves it
    /// above `highest`, the column below a surface held at a point, which stands linear in the
    /// point, so that the point may be settled after.
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
        tr.gain = if mass > 0.0 { frosted } else { bare };
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
    /// left to sublime warms the ground past the point.
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
    ) {
        if *mass <= 0.0 && (tr.clear || tr.free[0] >= point) {
            t.copy_from_slice(&tr.free);
            return;
        }
        let (found, taken) = self.deposit(tr, dt, point, radiated, latent, *mass);
        for (v, (b, l)) in t.iter_mut().zip(tr.base.iter().zip(&tr.lift)) {
            *v = b + point * l;
        }
        *mass = (*mass + taken).max(0.0);
        if taken > found {
            let g = self.conductance.first().copied().unwrap_or(0.0);
            let respond = self.capacity[0] + dt * (slope + g);
            t[0] = point + (taken - found) * latent / respond;
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
fn reader(n_lat: usize, n_lon: usize, lat: f64, lon: f64) -> [(usize, f64); 4] {
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
        (eb.cycle(&climate, &Ground::from_graph(graph), sea), climate)
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
        for lat in [-60.0, 60.0] {
            assert!(at(&cycle.frost_share, lat, 0.0) > 0.4, "{lat}");
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
        let cycle = eb.cycle(&climate, &Ground::from_graph(&graph), None);
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
        let cycle = eb.cycle(&climate, &Ground::from_graph(world), Some((&ocean, &share)));
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
    /// (Sergeev et al. 2022).
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
        let (n_lat, _) = cycle.layout;
        let swing = |lat: f64| {
            let i = ((lat + 90.0) / 180.0 * (n_lat - 1) as f64).round() as usize;
            let air = &cycle.air[0];
            let (lo, hi) = (0..air.len() / n_lat)
                .map(|k| air[k * n_lat + i])
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
    /// round the cycle's end, its phases from the northern spring equinox; all year as 0 and
    /// 1; none where it never lies.
    #[test]
    fn an_ice_s_window_is_its_longest_run() {
        let phase = |t: f64| (t / 10.0 - 0.25).rem_euclid(1.0);
        // Lying from step 8 round to step 2, its most at step 0; a shorter run at 5.
        let masses = [9.0f32, 5.0, 2.0, 0.0, 0.0, 3.0, 0.0, 0.0, 1.0, 4.0];
        let w = Window::of(&masses, 0.5, phase);
        let near = |a: f64, b: f64| (a - b).abs() < 1e-9;
        assert!(near(w.freeze, 0.6) && near(w.melt, 0.1), "{w:?}");
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
    /// lying all year.
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

    /// Prints the Earth-like's sea ice through its year near its north pole, kg/m², every
    /// fourth step from perihelion, over its sea at every latitude (EXOGEN_SEA_SHARE).
    #[test]
    #[ignore = "diagnostic: the sea ice's year"]
    fn print_sea_series() {
        let earth = crate::presets::preset("Earth").unwrap().to_graph().0;
        let eb = EnergyBalance::from_graph(&earth);
        let climate = eb.solve(CubeSphere::new(CLIMATE_FACE_SIZE));
        let (ocean, share) = sea_of(&earth).unwrap();
        let layout = eb.layout(&climate, Some(&share));
        let ground = Ground::from_graph(&earth);
        let year = eb.seasonal_years(&ground, &layout, eb.sky.obliquity, Some(&ocean));
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
        let (n_lat, _) = cycle.layout;
        for lat in [-60.0f64, 0.0, 45.0, 60.0, 80.0] {
            let i = ((lat + 90.0) / 180.0 * (n_lat - 1) as f64).round() as usize;
            let air = &cycle.air[0];
            let steps = air.len() / n_lat;
            let row: Vec<String> = (0..steps)
                .step_by(steps / 12)
                .map(|k| format!("{:.0}", air[k * n_lat + i]))
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
        let layout = eb.layout(&climate, None);
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
}
