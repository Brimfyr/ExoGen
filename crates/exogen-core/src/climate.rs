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
//! liquid, ice and snow brighten the ground the year leaves below freezing (North 1975).
//!
//! The balance is solved by Newton's method from the zero-dimensional mean temperature, the
//! warm branch, each step's linear system by conjugate gradients; the local slope is floored
//! at a share of the radiation's own so that the ice's feedback, which can make it negative,
//! leaves the system positive definite.

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
    mean_heat_capacity, mean_molecular_mass, regolith_heat, saturation_temperature,
    species_columns, tau_dry, vapour_pressure, volatile, water_column, water_optical_depth,
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
    /// The heat capacity of the air over a square metre, J m⁻² K⁻¹, warming and cooling
    /// with the ground under it.
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

    /// The frost point at a partial pressure of the gas, K; none at none.
    pub fn point_at(&self, partial: f64) -> f64 {
        if partial <= 0.0 {
            0.0
        } else {
            saturation_temperature(self.volatile, partial)
        }
    }
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
}

impl Traps {
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
/// share of the datum's pressure its height keeps, its year's balance (temperature, albedo,
/// transport) and the heat its air exchanges with the ground per kelvin, W/m²/K.
#[derive(Clone, Copy, Debug)]
struct TrapCell {
    lat: f64,
    lon: f64,
    thin: f64,
    steady: f64,
    albedo: f64,
    brought: f64,
    exchange: f64,
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
}

/// An ice's column and the bare ground's beside it in the traps (F240), their frost, and
/// their sums through a cycle.
#[derive(Clone)]
struct Pair {
    air: Vec<f64>,
    air_mass: f64,
    ice: Vec<f64>,
    ice_mass: f64,
    sums: Vec<f64>,
    air_sums: Vec<f64>,
    work: Work,
    trial: Trial,
}

impl Pair {
    fn new(column: &Column, steady: f64) -> Self {
        let n = column.len();
        Self {
            air: vec![steady; n],
            air_mass: 0.0,
            ice: vec![steady; n],
            ice_mass: 0.0,
            sums: vec![0.0; n],
            air_sums: vec![0.0; n],
            work: Work::new(n),
            trial: Trial::new(n),
        }
    }

    /// A new cycle.
    fn begin(&mut self) {
        self.sums.fill(0.0);
        self.air_sums.fill(0.0);
    }

    /// The deep ground of both columns settled at their surfaces' means over a cycle of
    /// `steps`.
    fn settle(&mut self, steps: usize) {
        settle(&mut self.ice, &self.sums, steps);
        settle(&mut self.air, &self.air_sums, steps);
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

/// The year's surface temperatures through the days and seasons (F236), K, and the frost
/// (F239).
#[derive(Clone, Debug)]
pub struct Cycle {
    /// The mean over the year.
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
    /// the model does not yet keep (5B's caps), where it lays frost.
    pub cap_growth: Option<f64>,
    /// The tilts the year was run at, and at each the frost's gas's partial pressure at the
    /// datum through the year's steps, Pa, where it lays frost (F240).
    pub tilts: Vec<f64>,
    pub partials: Vec<Vec<f64>>,
}

/// A place's year: its mean, warmest and coldest temperature, the share of the year frost
/// lies on it and its frost at the snapshot's season.
#[derive(Clone, Copy, Debug)]
struct Year {
    mean: f64,
    max: f64,
    min: f64,
    frost_share: f64,
    frost_now: f64,
}

/// The frost's year at a tilt (F239): the places' years, the surface pressure's range, Pa,
/// the caps' growth in a year over the body, kg/m², and the frost's gas's partial pressure at
/// each step of the year.
struct FrostYear {
    years: Vec<Year>,
    range: (f64, f64),
    growth: f64,
    partials: Vec<f64>,
}

/// The frost's lockstep through its last cycle (F239): the places' runs, the surface
/// pressure's range, Pa, the frost the caps grow by in the cycle over the body, kg/m², and
/// the frost's gas's partial pressure after each step.
struct Lockstep {
    runs: Vec<Run>,
    range: (f64, f64),
    growth: f64,
    partials: Vec<f64>,
}

impl Year {
    /// A year without frost.
    fn bare(mean: f64, max: f64, min: f64) -> Self {
        Self {
            mean,
            max,
            min,
            frost_share: 0.0,
            frost_now: 0.0,
        }
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
        Self {
            sky: Sky::from_graph(graph),
            albedo: graph.f64("atmosphere.bond_albedo"),
            low_sun: LOW_SUN_ALBEDO
                * (graph.f64("atmosphere.column_mass") / LOW_SUN_COLUMN).clamp(0.0, 1.0),
            ice: graph.bool("climate.ice_albedo_feedback")
                && graph.str("volatiles.working_liquid") == "water",
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
        let base: Vec<f64> = (0..n)
            .map(|k| (self.albedo + self.low_sun * (mean_height - height[k])).clamp(0.0, 1.0))
            .collect();
        let albedo = |k: usize, t: f64| -> (f64, f64) {
            let (share, slope) = self.ice_share(t);
            if ICE_ALBEDO <= base[k] {
                (base[k], 0.0)
            } else {
                (
                    base[k] + (ICE_ALBEDO - base[k]) * share,
                    (ICE_ALBEDO - base[k]) * slope,
                )
            }
        };
        let d = self.transport;
        let mut t = vec![self.start.max(1.0); n];
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
            transport: Field::from_fn(grid, |c, _| brought[grid.index(c)] as f32),
            iterations,
        }
    }

    /// The year's surface temperatures through the days and seasons on the climate's grid
    /// (F236): at each place a column of ground under the star's path, its sunlight absorbed
    /// at the albedo the year's balance found and the air's heat brought by its transport. A
    /// body that spins many times a year is run through the year on the day's mean sunlight
    /// until its seasons repeat, then through a day at each season; a body in a spin-orbit
    /// resonance through the orbits it repeats in, until they repeat.
    pub fn cycle(&self, climate: &Climate, ground: &Ground) -> Cycle {
        let grid = climate.temperature.grid();
        let at = |lat: f64, lon: f64| -> (f64, f64, f64) {
            let (c, _, _) = grid.locate(direction_from_lon_lat(lon, lat));
            (
                climate.temperature.get(c) as f64,
                climate.albedo.get(c) as f64,
                climate.transport.get(c) as f64,
            )
        };
        let all = self.sky.obliquities();
        let tilts: Vec<f64> = if all.len() > CYCLE_SWING_SAMPLES {
            (0..CYCLE_SWING_SAMPLES)
                .map(|j| all[j * all.len() / CYCLE_SWING_SAMPLES])
                .collect()
        } else {
            all
        };
        let combine = |each: Vec<Year>| -> Year {
            let n = each.len().max(1) as f64;
            each.iter().fold(
                Year {
                    mean: 0.0,
                    max: f64::MIN,
                    min: f64::MAX,
                    frost_share: 0.0,
                    frost_now: 0.0,
                },
                |acc, e| Year {
                    mean: acc.mean + e.mean / n,
                    max: acc.max.max(e.max),
                    min: acc.min.min(e.min),
                    frost_share: acc.frost_share + e.frost_share / n,
                    frost_now: acc.frost_now + e.frost_now / n,
                },
            )
        };
        let n_lat = match self.sky.resonance {
            None => CYCLE_LATITUDES,
            Some(_) => CYCLE_RESONANT_GRID.0,
        };
        let lat_of = |i: usize| -PI / 2.0 + PI * i as f64 / (n_lat - 1) as f64;
        let mut pressure_range = None;
        let mut cap_growth = None;
        let mut partials = Vec::new();
        // A frost's runs, one a tilt: their places' years combined, the pressure's range the
        // widest and the caps' growth the mean.
        let mut gather = |runs: Vec<FrostYear>| -> Vec<Year> {
            pressure_range = runs
                .iter()
                .map(|r| r.range)
                .reduce(|a, b| (a.0.min(b.0), a.1.max(b.1)));
            cap_growth =
                Some(runs.iter().map(|r| r.growth).sum::<f64>() / runs.len().max(1) as f64);
            partials = runs.iter().map(|r| r.partials.clone()).collect();
            (0..runs.first().map_or(0, |r| r.years.len()))
                .map(|i| combine(runs.iter().map(|r| r.years[i]).collect()))
                .collect()
        };
        // A place of the frost's lockstep, its latitude's band of the sphere's area shared
        // among its longitudes.
        let frost_place = |i: usize, lon: f64, share: f64| {
            let (lat, half) = (lat_of(i), PI / (n_lat - 1) as f64 / 2.0);
            let band = (lat + half).min(PI / 2.0).sin() - (lat - half).max(-PI / 2.0).sin();
            let (steady, albedo, brought) = at(lat, lon);
            Place {
                lat,
                lon,
                weight: band * share,
                steady,
                albedo,
                brought,
            }
        };
        let (rows, n_lon): (Vec<Year>, Option<usize>) = match (self.sky.resonance, self.frost) {
            // A frost's caps draw on one atmosphere, so its places run in step.
            (None, Some(frost)) => {
                let places: Vec<Place> = (0..n_lat).map(|i| frost_place(i, 0.0, 1.0)).collect();
                let runs = tilts
                    .iter()
                    .map(|&eps| self.frosty_years(ground, &places, eps, &frost))
                    .collect();
                (gather(runs), None)
            }
            (None, None) => (
                (0..n_lat)
                    .into_par_iter()
                    .map(|i| {
                        let (steady, albedo, brought) = at(lat_of(i), 0.0);
                        combine(
                            tilts
                                .iter()
                                .map(|&eps| {
                                    self.zonal_cycle(
                                        ground,
                                        lat_of(i),
                                        eps,
                                        steady,
                                        albedo,
                                        brought,
                                    )
                                })
                                .collect(),
                        )
                    })
                    .collect(),
                None,
            ),
            (Some((_, orbits)), frost) => {
                let n_lon = CYCLE_RESONANT_GRID.1;
                let lon_of = |j: usize| -PI + TAU * j as f64 / n_lon as f64;
                let paths: Vec<Vec<PathStep>> = tilts
                    .iter()
                    .map(|&eps| self.sky.resonant_path(eps, orbits))
                    .collect();
                let rows = match frost {
                    Some(frost) => {
                        let share = 1.0 / n_lon as f64;
                        let places: Vec<Place> = (0..n_lat * n_lon)
                            .map(|k| frost_place(k / n_lon, lon_of(k % n_lon), share))
                            .collect();
                        let runs = paths
                            .iter()
                            .map(|path| self.frosty_orbits(ground, orbits, path, &places, &frost))
                            .collect();
                        gather(runs)
                    }
                    None => (0..n_lat * n_lon)
                        .into_par_iter()
                        .map(|k| {
                            let (i, j) = (k / n_lon, k % n_lon);
                            let lon = lon_of(j);
                            let (steady, albedo, brought) = at(lat_of(i), lon);
                            combine(
                                paths
                                    .iter()
                                    .map(|path| {
                                        let place = (lat_of(i), lon);
                                        let start = (steady, albedo, brought);
                                        self.resonant_cycle(ground, orbits, path, place, start)
                                    })
                                    .collect(),
                            )
                        })
                        .collect(),
                };
                (rows, Some(n_lon))
            }
        };
        let field = |pick: fn(&Year) -> f64| {
            let values: Vec<f64> = rows.iter().map(pick).collect();
            Field::from_fn(grid, |_, d| {
                let (lon, lat) = lon_lat_from_direction(d);
                let value = match n_lon {
                    None => read_row(&values, lat),
                    Some(n_lon) => read_grid(&values, n_lon, lat, lon),
                };
                value as f32
            })
        };
        Cycle {
            mean: field(|r| r.mean),
            max: field(|r| r.max),
            min: field(|r| r.min),
            frost_share: field(|r| r.frost_share),
            frost_now: field(|r| r.frost_now),
            frost: self.frost,
            pressure_range,
            cap_growth,
            tilts,
            partials,
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

    /// The years of a fast spinner whose atmosphere lays frost (F239): its latitudes run in
    /// step through the year on the day's mean sunlight, then each through a day at each
    /// season from the ground, its frost and the frost point as the season left them.
    fn frosty_years(
        &self,
        ground: &Ground,
        places: &[Place],
        eps: f64,
        frost: &Frost,
    ) -> FrostYear {
        let olr = |t: f64| self.olr_and_slope(t);
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
        let lockstep = self.frosty_cycles(&column, places, timing, light, frost, every);
        let dt = day / CYCLE_DAY_STEPS as f64;
        let years = lockstep
            .runs
            .par_iter()
            .zip(places.par_iter())
            .map(|(r, p)| {
                let mut work = Work::new(column.len());
                let mut trial = Trial::new(column.len());
                let (sp, cp) = p.lat.sin_cos();
                let (mut mean, mut hi, mut lo) = (0.0, f64::MIN, f64::MAX);
                for (m, state, frost_mass, point) in &r.records {
                    let f = Frost {
                        point: *point,
                        ..*frost
                    };
                    let (flux, decl) = sun[*m];
                    let (sd, cd) = decl.sin_cos();
                    let mut t = state.clone();
                    let mut mass = *frost_mass;
                    let mut sum = 0.0;
                    for d in 0..CYCLE_DAYS {
                        for h in 0..CYCLE_DAY_STEPS {
                            let hour = TAU * (h as f64 + 0.5) / CYCLE_DAY_STEPS as f64 - PI;
                            let light = flux * (sp * sd + cp * cd * hour.cos()).max(0.0);
                            let gains = (
                                (1.0 - p.albedo) * light + p.brought,
                                (1.0 - f.albedo) * light + p.brought.max(0.0),
                            );
                            let (rows, heat) = ((&mut work, &mut trial), self.heat_flux);
                            column.frosty_step(&mut t, rows, dt, gains, heat, olr, &f, &mut mass);
                            if d + 1 == CYCLE_DAYS {
                                sum += t[0];
                                hi = hi.max(t[0]);
                                lo = lo.min(t[0]);
                            }
                        }
                    }
                    mean += sum / CYCLE_DAY_STEPS as f64;
                }
                Year {
                    mean: mean / r.records.len().max(1) as f64,
                    max: hi,
                    min: lo,
                    frost_share: r.frosted as f64 / steps as f64,
                    frost_now: r.now,
                }
            })
            .collect();
        FrostYear {
            years,
            range: lockstep.range,
            growth: lockstep.growth,
            partials: lockstep.partials,
        }
    }

    /// The cycles of a resonant body whose atmosphere lays frost (F239): its places run in
    /// step through the orbits it repeats in.
    fn frosty_orbits(
        &self,
        ground: &Ground,
        orbits: u32,
        path: &[PathStep],
        places: &[Place],
        frost: &Frost,
    ) -> FrostYear {
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
        let lockstep = self.frosty_cycles(&column, places, timing, light, frost, 0);
        let years = lockstep
            .runs
            .iter()
            .map(|r| Year {
                mean: r.sums[0] / steps as f64,
                max: r.hi,
                min: r.lo,
                frost_share: r.frosted as f64 / steps as f64,
                frost_now: r.now,
            })
            .collect();
        FrostYear {
            years,
            range: lockstep.range,
            growth: lockstep.growth / orbits.max(1) as f64,
            partials: lockstep.partials,
        }
    }

    /// Places run in step through their cycles while their frost draws on one atmosphere
    /// (F239): the frost the cold ground takes is the atmosphere's own, so that its pressure,
    /// and the frost point with it, falls as the caps grow and rises as they sublime (Leighton
    /// and Murray 1966). At each step every place tries its step, and the point is the one at
    /// which the frost they take together leaves the gas at the pressure whose point it is: a
    /// thin atmosphere's pressure follows its frost's balance through the step. Frosted
    /// ground takes the heat the air brings it but gives it none: the cycle's transport is the
    /// year's steady one (F236), and ground its frost holds below its steady temperature
    /// would draw heat from its surroundings rather than shed it. The gas's
    /// whole mass is set so that its pressure over the cycle averages the body's where the
    /// frost allows it. Frost that lies all through the cycle and grows is a cap's growth,
    /// which the model does not yet keep where the gas is given back (5B's caps): it is
    /// counted, and taken off in the share given. The cycle is `steps` steps of `dt` with the
    /// snapshot's step, each place lit at a step as `light` gives, the ground recorded every
    /// `every` steps (never at 0).
    fn frosty_cycles(
        &self,
        column: &Column,
        places: &[Place],
        (steps, dt, snapshot): (usize, f64, usize),
        light: impl Fn(usize, &Place) -> f64 + Sync,
        frost: &Frost,
        every: usize,
    ) -> Lockstep {
        let olr = |t: f64| self.olr_and_slope(t);
        let whole = places.iter().map(|p| p.weight).sum::<f64>().max(1e-30);
        let over = |runs: &[Run], pick: fn(&Run) -> f64| {
            runs.iter()
                .zip(places)
                .map(|(r, p)| pick(r) * p.weight)
                .sum::<f64>()
                / whole
        };
        let partial = |air: f64| frost.partial * air.max(0.0) / frost.column.max(1e-30);
        let mut runs: Vec<Run> = places
            .iter()
            .map(|p| Run::new(column.len(), p.steady))
            .collect();
        let (mut held, mut load) = (frost.column, 0.0);
        let mut range = (f64::MAX, f64::MIN);
        let mut growth = 0.0;
        let mut partials = Vec::with_capacity(steps);
        let mut last = vec![f64::NAN; places.len()];
        for _ in 0..CYCLE_MAX_YEARS {
            range = (f64::MAX, f64::MIN);
            partials.clear();
            runs.iter_mut().for_each(Run::begin);
            for k in 0..steps {
                let highest = frost.point_at(partial(held));
                runs.par_iter_mut()
                    .zip(places.par_iter())
                    .for_each(|(r, p)| {
                        let lit = light(k, p);
                        let gains = (
                            (1.0 - p.albedo) * lit + p.brought,
                            (1.0 - frost.albedo) * lit + p.brought.max(0.0),
                        );
                        let (rows, heat) = ((&mut r.work, &mut r.trial), self.heat_flux);
                        column.trial(&r.t, rows, dt, gains, heat, olr, r.mass, highest);
                    });
                let air = (held - load, whole);
                let point = self.frost_point(column, &runs, places, air, dt, frost, highest);
                let shed = olr(point);
                runs.par_iter_mut().for_each(|r| {
                    column.finish(
                        &mut r.t,
                        &r.trial,
                        dt,
                        point,
                        shed,
                        frost.latent,
                        &mut r.mass,
                    );
                    r.record(k, snapshot, every, point);
                });
                load = over(&runs, |r| r.mass);
                let now = partial(held - load);
                range = (range.0.min(now), range.1.max(now));
                partials.push(now);
            }
            // The frost that lay all cycle and grew through it: a cap's growth.
            growth = over(&runs, Run::growth);
            let means: Vec<f64> = runs.iter().map(|r| r.sums[0] / steps as f64).collect();
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
            let mean_load = over(&runs, |r| r.frost_sum) / steps as f64;
            let kept = ((held - mean_load) / frost.column.max(1e-30)).clamp(0.0, 1.0);
            for r in runs.iter_mut() {
                r.mass -= kept * r.growth();
                let sums = std::mem::take(&mut r.sums);
                settle(&mut r.t, &sums, steps);
                r.sums = sums;
            }
            load = over(&runs, |r| r.mass);
            held = (frost.column + (mean_load - growth) * kept).max(load);
        }
        // The surface pressure: the other gases' and the frost's gas's.
        let others = (self.pressure - frost.partial).max(0.0);
        Lockstep {
            runs,
            range: (others + range.0, others + range.1),
            growth,
            partials,
        }
    }

    /// The frost point a step of the lockstep settles at (F239): the point at which the frost
    /// the places take through the step leaves the gas, `air` kg/m² over the body before it,
    /// at the partial pressure whose point it is. The frost taken rises with the point and the
    /// point of what is left falls, so the two cross once, between none and `highest`, the
    /// point of all the gas; found by false position, halving the side that stays. Frosted
    /// ground gives the air none of its heat (`frosty_cycles`), so that at no point does it
    /// take frost the air does not hold.
    #[allow(clippy::too_many_arguments)]
    fn frost_point(
        &self,
        column: &Column,
        runs: &[Run],
        places: &[Place],
        (air, whole): (f64, f64),
        dt: f64,
        frost: &Frost,
        highest: f64,
    ) -> f64 {
        let taking: Vec<(&Run, f64)> = runs
            .iter()
            .zip(places)
            .filter(|(r, _)| r.mass > 0.0 || !r.trial.clear)
            .map(|(r, p)| (r, p.weight / whole))
            .collect();
        let gap = |point: f64| {
            let radiated = self.olr(point);
            let taken: f64 = taking
                .iter()
                .map(|(r, w)| {
                    column
                        .deposit(&r.trial, dt, point, radiated, frost.latent, r.mass)
                        .1
                        * w
                })
                .sum();
            point - frost.point_at(frost.partial * (air - taken).max(0.0) / frost.column.max(1e-30))
        };
        let (mut lo, mut hi) = (0.0, highest);
        let (mut g_lo, mut g_hi) = (gap(lo), gap(hi));
        if g_hi <= 0.0 {
            return hi;
        }
        // The low side always leaves the air the gas the frost does not take.
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
    /// year was run at, beside the bare ground about it, whose air exchanges heat with the
    /// ice as F237's does with a slope (`pair_step`). The frost's own gas's ice lies all year from an unlimited store, held at
    /// its point at the cell's partial pressure, and gains or loses what the surface's
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
        let cell_at = |lat: f64, lon: f64, thin: f64| {
            let (c, _, _) = grid.locate(direction_from_lon_lat(lon, lat));
            TrapCell {
                lat,
                lon,
                thin,
                steady: climate.temperature.get(c) as f64,
                albedo: climate.albedo.get(c) as f64,
                brought: climate.transport.get(c) as f64,
                exchange: downscale.exchange,
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
        };
        for ice in ices {
            let own = frost.filter(|f| f.species == ice.species);
            let values: Vec<f64> = places
                .par_iter()
                .map(|cell| {
                    let n = cycle.tilts.len().max(1) as f64;
                    cycle
                        .tilts
                        .iter()
                        .enumerate()
                        .map(|(j, &eps)| {
                            let partials = cycle.partials.get(j).map_or(&[][..], |p| &p[..]);
                            match own {
                                Some(f) => self.ice_growth(ground, cell, eps, partials, f),
                                None => self.ice_vapour(ground, cell, eps, partials, frost, ice),
                            }
                        })
                        .sum::<f64>()
                        / n
                })
                .collect();
            // The vapour at the datum's air: the saturation pressure over the ice over the
            // share of that air's pressure its height keeps.
            let field = Field::from_fn(grid, |c, d| {
                let k = grid.index(c);
                let value = match &levels {
                    Some(levels) => {
                        let (_, lat) = lon_lat_from_direction(d);
                        read_table(&values, levels, lat, thins[k], own.is_none())
                    }
                    None => values[k],
                };
                match own {
                    Some(_) => value as f32,
                    None => (value / thins[k]) as f32,
                }
            });
            match own {
                Some(_) => traps.growth = Some((*ice, field)),
                None => traps.vapour.push((*ice, field)),
            }
        }
        traps
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
                let path = self.sky.resonant_path(eps, orbits);
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
    /// pressure, it gains or loses what the surface's balance leaves (F239), the deep ground
    /// settled under it.
    fn ice_growth(
        &self,
        ground: &Ground,
        cell: &TrapCell,
        eps: f64,
        partials: &[f64],
        frost: &Frost,
    ) -> f64 {
        let (column, _) = self.trap_column(ground, eps);
        let (light, _, dt, years) = self.trap_sun(cell, eps);
        let mut pair = Pair::new(&column, cell.steady);
        let mut gained = 0.0;
        for _ in 0..TRAP_YEARS {
            pair.ice_mass = TRAP_STORE;
            pair.begin();
            for (k, &lit) in light.iter().enumerate() {
                let f = cell.laid(Some(frost), partials, k);
                self.pair_step(&column, cell, &mut pair, dt, lit, f.as_ref(), frost.albedo);
            }
            gained = (pair.ice_mass - TRAP_STORE) / years;
            pair.settle(light.len());
        }
        gained
    }

    /// The year-mean saturation pressure over a cell's ice of another gas at a tilt (F240),
    /// Pa: the frost's gas laying frost on it where the cell's partial pressure's point
    /// allows (F239), run
    /// through its year until it repeats and then, a fast spinner's, through a day at each
    /// season.
    fn ice_vapour(
        &self,
        ground: &Ground,
        cell: &TrapCell,
        eps: f64,
        partials: &[f64],
        frost: Option<&Frost>,
        ice: &Ice,
    ) -> f64 {
        let (column, day) = self.trap_column(ground, eps);
        let (light, sun, dt, _) = self.trap_sun(cell, eps);
        let albedo = ice.albedo();
        let mut pair = Pair::new(&column, cell.steady);
        let mut last = f64::NAN;
        let mut seasons: Vec<(usize, Pair)> = Vec::new();
        let mut vapour = 0.0;
        let every = (light.len() / CYCLE_SEASONS).max(1);
        for _ in 0..CYCLE_MAX_YEARS {
            seasons.clear();
            vapour = 0.0;
            pair.begin();
            for (k, &lit) in light.iter().enumerate() {
                let f = cell.laid(frost, partials, k);
                self.pair_step(&column, cell, &mut pair, dt, lit, f.as_ref(), albedo);
                vapour += vapour_pressure(ice.volatile, pair.ice[0]);
                if !sun.is_empty() && k.is_multiple_of(every) {
                    seasons.push((k, pair.clone()));
                }
            }
            let surface = pair.sums[0] / light.len() as f64;
            if (surface - last).abs() < CYCLE_SETTLED_K {
                break;
            }
            last = surface;
            pair.settle(light.len());
        }
        if sun.is_empty() {
            return vapour / light.len() as f64;
        }
        // A fast spinner's days at each season, from the ground and its frost as the season
        // left them: the saturation pressure through the last of them.
        let (sp, cp) = cell.lat.sin_cos();
        let dt = day / CYCLE_DAY_STEPS as f64;
        let mut sum = 0.0;
        for (k, state) in &seasons {
            let (flux, decl) = sun[*k];
            let (sd, cd) = decl.sin_cos();
            let f = cell.laid(frost, partials, *k);
            let mut pair = state.clone();
            for d in 0..CYCLE_DAYS {
                for h in 0..CYCLE_DAY_STEPS {
                    let hour = TAU * (h as f64 + 0.5) / CYCLE_DAY_STEPS as f64 - PI;
                    let lit = flux * (sp * sd + cp * cd * hour.cos()).max(0.0);
                    self.pair_step(&column, cell, &mut pair, dt, lit, f.as_ref(), albedo);
                    if d + 1 == CYCLE_DAYS {
                        sum += vapour_pressure(ice.volatile, pair.ice[0]);
                    }
                }
            }
        }
        sum / (seasons.len().max(1) * CYCLE_DAY_STEPS) as f64
    }

    /// A cell's ice and the bare ground about it through a step (F240): the bare ground as
    /// the cycle takes it, the frost laying its frost on it where `frost` allows (F239); then
    /// the ice at `albedo`, the frost on it in turn, the air over it at the bare ground's
    /// temperature exchanging heat with it as F237's air does with a slope, so that a thick
    /// atmosphere holds the ice near its own temperature and a thin one leaves it its own
    /// sunlight.
    #[allow(clippy::too_many_arguments)]
    fn pair_step(
        &self,
        column: &Column,
        cell: &TrapCell,
        pair: &mut Pair,
        dt: f64,
        lit: f64,
        frost: Option<&Frost>,
        albedo: f64,
    ) {
        let olr = |t: f64| self.olr_and_slope(t);
        let heat = self.heat_flux;
        let bare = (1.0 - cell.albedo) * lit + cell.brought;
        let frosted = |f: &Frost| (1.0 - f.albedo) * lit + cell.brought.max(0.0);
        let rows = (&mut pair.work, &mut pair.trial);
        match frost {
            Some(f) => {
                let gains = (bare, frosted(f));
                column.frosty_step(
                    &mut pair.air,
                    rows,
                    dt,
                    gains,
                    heat,
                    olr,
                    f,
                    &mut pair.air_mass,
                );
            }
            None => column.step(&mut pair.air, rows.0, dt, bare, heat, olr),
        }
        let (air, exchange) = (pair.air[0], cell.exchange);
        let shed = move |t: f64| {
            let (r, slope) = self.olr_and_slope(t);
            (r + exchange * (t - air), slope + exchange)
        };
        let own = (1.0 - albedo) * lit + cell.brought;
        let rows = (&mut pair.work, &mut pair.trial);
        match frost {
            Some(f) => {
                let gains = (own, frosted(f));
                column.frosty_step(
                    &mut pair.ice,
                    rows,
                    dt,
                    gains,
                    heat,
                    shed,
                    f,
                    &mut pair.ice_mass,
                );
            }
            None => column.step(&mut pair.ice, rows.0, dt, own, heat, shed),
        }
        for (s, v) in pair.sums.iter_mut().zip(&pair.ice) {
            *s += v;
        }
        for (s, v) in pair.air_sums.iter_mut().zip(&pair.air) {
            *s += v;
        }
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
        let mu = mean_molecular_mass(composition);
        let t = graph.f64("atmosphere.surface_temperature_mean").max(1.0);
        let density = pressure * mu / (GAS_CONSTANT * t);
        Self {
            gases: [
                tau("CO2", TAU_CO2),
                tau("CH4", TAU_CH4),
                tau("N2", TAU_N2_CIA),
            ],
            water,
            scale_height: graph.f64("atmosphere.scale_height"),
            lapse: graph.f64("climate.lapse_rate"),
            exchange: density
                * mean_heat_capacity(composition)
                * AIR_EXCHANGE_COEFFICIENT
                * SURFACE_WIND_SPEED,
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

/// A column of ground: each layer's heat capacity per square metre, the surface's carrying
/// the air's, and the conductance from each layer to the one below.
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

/// A place of the frost's lockstep (F239): its latitude and longitude, rad, its share of the
/// sphere's area, and its year's balance: temperature, albedo and transport.
#[derive(Clone, Copy, Debug)]
struct Place {
    lat: f64,
    lon: f64,
    weight: f64,
    steady: f64,
    albedo: f64,
    brought: f64,
}

/// A place's ground and frost through the frost's lockstep (F239), and its last cycle.
struct Run {
    t: Vec<f64>,
    mass: f64,
    work: Work,
    trial: Trial,
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
    /// The step, the ground, its frost and the point at each step recorded.
    records: Vec<(usize, Vec<f64>, f64, f64)>,
}

impl Run {
    fn new(n: usize, steady: f64) -> Self {
        Self {
            t: vec![steady; n],
            mass: 0.0,
            work: Work::new(n),
            trial: Trial::new(n),
            sums: vec![0.0; n],
            hi: f64::MIN,
            lo: f64::MAX,
            frosted: 0,
            now: 0.0,
            frost_start: 0.0,
            frost_sum: 0.0,
            frost_least: f64::MAX,
            records: Vec::new(),
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
        self.records.clear();
    }

    /// The record of the cycle's `k`th step, just taken at `point`.
    fn record(&mut self, k: usize, snapshot: usize, every: usize, point: f64) {
        for (s, v) in self.sums.iter_mut().zip(&self.t) {
            *s += v;
        }
        (self.hi, self.lo) = (self.hi.max(self.t[0]), self.lo.min(self.t[0]));
        if self.mass > 0.0 {
            self.frosted += 1;
        }
        self.frost_sum += self.mass;
        self.frost_least = self.frost_least.min(self.mass);
        if k == snapshot {
            self.now = self.mass;
        }
        if every > 0 && k.is_multiple_of(every) {
            self.records.push((k, self.t.clone(), self.mass, point));
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
        let mut capacity: Vec<f64> = dz.iter().map(|d| ground.heat * d).collect();
        capacity[0] += ground.air;
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

    /// One step with the frost at its point (F239): the step's trial, ended at the point.
    #[allow(clippy::too_many_arguments)]
    fn frosty_step(
        &self,
        t: &mut [f64],
        (w, tr): (&mut Work, &mut Trial),
        dt: f64,
        gains: (f64, f64),
        bottom: f64,
        olr: impl Fn(f64) -> (f64, f64) + Copy,
        f: &Frost,
        mass: &mut f64,
    ) {
        self.trial(t, (w, tr), dt, gains, bottom, olr, *mass, f.point);
        self.finish(t, tr, dt, f.point, olr(f.point), f.latent, mass);
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

    fn cycle_of(name: &str) -> (Cycle, Climate) {
        let graph = crate::presets::preset(name).unwrap().to_graph().0;
        let eb = EnergyBalance::from_graph(&graph);
        let climate = eb.solve(CubeSphere::new(CLIMATE_FACE_SIZE));
        (eb.cycle(&climate, &Ground::from_graph(&graph)), climate)
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
            let (cycle, climate) = cycle_of(name);
            let took = start.elapsed().as_secs_f64();
            let (mean, equator, poles) = bands(&cycle.mean);
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

    /// F236: a thick atmosphere holds the day and the seasons: Venus's and Titan's means
    /// through them stand at their steady balances. The Earth-like's seasons swing as land's
    /// do until its ocean holds them (5B), its summer poles too warm, its mean 2 K under its
    /// steady balance and still within the gate's 5 K of Earth's.
    #[test]
    fn a_thick_atmosphere_keeps_its_steady_mean() {
        for (name, within) in [("Venus", 0.5), ("Titan", 0.5), ("Earth", 3.0)] {
            let (cycle, climate) = cycle_of(name);
            let (mean, ..) = bands(&cycle.mean);
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

    /// F239: a thin atmosphere's frost draws on its own gas: Pluto's frost holds no more
    /// than the 1.6 kg/m² its atmosphere weighs, and the pressure follows the frost's balance
    /// through the year. That balance holds it far under the 1 Pa measured (Stern et al.
    /// 2015): the model's nitrogen frost radiates as a black body where the ice's emissivity
    /// is some 0.6 (Stansberry et al. 1996), and Pluto's reservoir of nitrogen ice waits on
    /// 5B's caps.
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
        let (lo, hi) = cycle
            .pressure_range
            .expect("its pressure follows the frost");
        assert!(lo < hi && hi < f.partial, "{lo} to {hi}");
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
        let cycle = eb.cycle(&climate, &Ground::from_graph(&graph));
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
            let climate = eb.solve(CubeSphere::new(CLIMATE_FACE_SIZE));
            let ground = Ground::from_graph(&graph);
            let cycle = eb.cycle(&climate, &ground);
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
        let climate = eb.solve(CubeSphere::new(CLIMATE_FACE_SIZE));
        let ground = Ground::from_graph(&graph);
        let cycle = eb.cycle(&climate, &ground);
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

    fn downscale_of(name: &str) -> Downscale {
        Downscale::from_graph(&crate::presets::preset(name).unwrap().to_graph().0)
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
