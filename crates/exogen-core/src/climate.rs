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

use rayon::prelude::*;

use crate::constants::*;
use crate::fields::Field;
use crate::grid::{CubeSphere, lon_lat_from_direction};
use crate::insolation::Sky;
use crate::params::ParamGraph;
use crate::params::derive::physics::{species_columns, tau_dry, water_optical_depth};

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
}

/// The year's mean climate on the climate grid.
#[derive(Clone, Debug)]
pub struct Climate {
    /// Surface temperature, K.
    pub temperature: Field<f32>,
    /// The sunlight absorbed, W/m².
    pub absorbed: Field<f32>,
    /// Newton's steps to converge.
    pub iterations: usize,
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
        }
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
        let absorbed: Vec<f64> = (0..n)
            .map(|k| flux[k] * (1.0 - albedo(k, t[k]).0))
            .collect();
        Climate {
            temperature: Field::from_fn(grid, |c, _| t[grid.index(c)] as f32),
            absorbed: Field::from_fn(grid, |c, _| absorbed[grid.index(c)] as f32),
            iterations,
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
