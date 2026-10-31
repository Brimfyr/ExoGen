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
use crate::grid::{CubeSphere, direction_from_lon_lat, lon_lat_from_direction};
use crate::insolation::{PathStep, Sky, diurnal_share, read_grid, read_row};
use crate::params::ParamGraph;
use crate::params::derive::physics::{
    mean_heat_capacity, regolith_heat, species_columns, tau_dry, water_optical_depth,
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

/// The year's surface temperatures through the days and seasons (F236), K.
#[derive(Clone, Debug)]
pub struct Cycle {
    /// The mean over the year.
    pub mean: Field<f32>,
    /// The warmest the surface stands.
    pub max: Field<f32>,
    /// The coldest.
    pub min: Field<f32>,
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
        let combine = |each: Vec<(f64, f64, f64)>| -> (f64, f64, f64) {
            let n = each.len().max(1) as f64;
            each.iter().fold((0.0, f64::MIN, f64::MAX), |acc, e| {
                (acc.0 + e.0 / n, acc.1.max(e.1), acc.2.min(e.2))
            })
        };
        let n_lat = match self.sky.resonance {
            None => CYCLE_LATITUDES,
            Some(_) => CYCLE_RESONANT_GRID.0,
        };
        let lat_of = |i: usize| -PI / 2.0 + PI * i as f64 / (n_lat - 1) as f64;
        let (rows, n_lon): (Vec<(f64, f64, f64)>, Option<usize>) = match self.sky.resonance {
            None => (
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
            Some((_, orbits)) => {
                let n_lon = CYCLE_RESONANT_GRID.1;
                let paths: Vec<Vec<PathStep>> = tilts
                    .iter()
                    .map(|&eps| self.sky.resonant_path(eps, orbits))
                    .collect();
                (
                    (0..n_lat * n_lon)
                        .into_par_iter()
                        .map(|k| {
                            let (i, j) = (k / n_lon, k % n_lon);
                            let lon = -PI + TAU * j as f64 / n_lon as f64;
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
                    Some(n_lon),
                )
            }
        };
        let field = |pick: fn(&(f64, f64, f64)) -> f64| {
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
            mean: field(|r| r.0),
            max: field(|r| r.1),
            min: field(|r| r.2),
        }
    }

    /// A fast spinner's year at a latitude and tilt: the mean, the warmest and the coldest.
    fn zonal_cycle(
        &self,
        ground: &Ground,
        lat: f64,
        eps: f64,
        steady: f64,
        albedo: f64,
        brought: f64,
    ) -> (f64, f64, f64) {
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
        (mean / seasons.len().max(1) as f64, hi, lo)
    }

    /// A resonant body's cycle at a place (latitude, longitude), started from its year's
    /// balance (temperature, albedo, transport): the mean, the warmest and the coldest.
    fn resonant_cycle(
        &self,
        ground: &Ground,
        orbits: u32,
        path: &[PathStep],
        (lat, lon): (f64, f64),
        (steady, albedo, brought): (f64, f64, f64),
    ) -> (f64, f64, f64) {
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
        (mean, hi, lo)
    }
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
        olr: impl Fn(f64) -> (f64, f64),
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
