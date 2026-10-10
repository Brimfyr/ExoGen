// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! Phase 1 placeholder terrain process.
//!
//! It exists so that the loop, keyframes, scrubber, resampling, hydrosphere and export can
//! be built and tested end to end before the physical processes arrive. It is still built
//! from the pieces those processes will use: a low-degree primordial crust with an
//! isostatic elevation (PLAN.md 4.1), a seeded noise detail term scaled by the relief limit,
//! a slow time-varying perturbation so keyframes differ, patchy resurfacing that stamps
//! the surface-age layer, and diffusion. Nothing here is a model of a real process.
//!
//! Time variation is a blend between two noise fields, so the pattern changes in place. Two
//! rules learned here apply to every later process: a pattern must never be moved by
//! translating the 3D noise sample point, because a translation is not a rigid motion on the
//! sphere and compresses the pattern on one hemisphere while stretching it on the other; and
//! a pattern whose effect accumulates over steps must not drift at all, or it smears into
//! streaks. A true drift of a non-accumulating field is a rotation about an axis (`rotate`).

use rayon::prelude::*;

use crate::constants::{DENSITY_CRUST_ICE, DENSITY_CRUST_ROCK};
use crate::fields::{Field, LayerSet};
use crate::grid::CubeSphere;
use crate::noise::GradientNoise;
use crate::params::ParamGraph;
use crate::rng::Rng;
use crate::spharm::Harmonics;

/// Bumped when the behaviour changes, so caches are not reused across versions.
pub const VERSION: u64 = 2;

/// Geological unit labels written by this process.
pub const UNIT_PRIMORDIAL: u16 = 1;
pub const UNIT_RESURFACED: u16 = 2;

/// Two-level crust for mobile-lid bodies: thin oceanic crust and thick continental crust
/// covering a fraction of the area, which gives the bimodal hypsometry that lets an ocean
/// fill the basins and leave continents standing. Plate tectonics replaces it in Phase 6.
struct Continents {
    oceanic_thickness: f64,
    continental_thickness: f64,
    area_fraction: f64,
}

pub struct PlaceholderTerrain {
    grid: CubeSphere,
    noise: GradientNoise,
    noise_b: GradientNoise,
    harmonics: Harmonics,
    dichotomy_axis: [f64; 3],
    dichotomy_amplitude: f64,
    crust_mean: f64,
    crust_sigma: f64,
    continents: Option<Continents>,
    relief: f64,
    isostatic_factor: f64,
    age: f64,
}

impl PlaceholderTerrain {
    pub fn new(graph: &ParamGraph, grid: CubeSphere) -> Self {
        let seed = graph.f64("seed.terrain") as u64;
        let rng = Rng::new(seed);
        let icy = matches!(
            graph.str("body.composition_class"),
            "rock_ice" | "ice_dominated"
        );
        let (rho_mantle, rho_crust) = if icy {
            (1000.0, DENSITY_CRUST_ICE)
        } else {
            (3300.0, DENSITY_CRUST_ROCK)
        };
        let relief = 0.25 * graph.f64("interior.max_relief") * graph.f64("multiplier.relief");
        let continents = if graph.str("interior.tectonic_regime") == "mobile_lid" {
            Some(Continents {
                oceanic_thickness: 8_000.0,
                continental_thickness: 38_000.0,
                area_fraction: 0.30,
            })
        } else {
            None
        };
        Self {
            grid,
            noise: GradientNoise::new(seed ^ 0x9E37_79B9),
            noise_b: GradientNoise::new(seed ^ 0x7F4A_7C15),
            harmonics: Harmonics::random(&rng, 2, 8, 2.0),
            dichotomy_axis: rng.unit_direction(1, 0),
            dichotomy_amplitude: graph.f64("tectonics.dichotomy_amplitude"),
            crust_mean: graph.f64("tectonics.initial_crust_thickness"),
            crust_sigma: 0.12 * graph.f64("tectonics.initial_crust_thickness"),
            continents,
            relief,
            isostatic_factor: (rho_mantle - rho_crust) / rho_mantle,
            age: graph.f64("body.age"),
        }
    }

    /// Initial layered state at formation.
    pub fn init(&self) -> LayerSet {
        let grid = self.grid;
        let mut layers = LayerSet::new(grid);

        // Normalise the harmonic field to unit rms over the grid.
        let raw: Vec<f64> = (0..grid.cell_count())
            .into_par_iter()
            .map(|k| self.harmonics.eval(grid.center(grid.cell(k))))
            .collect();
        let rms = (raw.iter().map(|v| v * v).sum::<f64>() / raw.len() as f64)
            .sqrt()
            .max(1e-12);

        let axis = self.dichotomy_axis;
        let crust: Vec<f32> = match &self.continents {
            Some(cont) => {
                // Threshold at the quantile that leaves the requested continental area.
                let mut sorted = raw.clone();
                sorted.sort_by(|a, b| a.total_cmp(b));
                let thr = sorted[((1.0 - cont.area_fraction) * (sorted.len() - 1) as f64) as usize];
                let width = 0.35 * rms;
                (0..grid.cell_count())
                    .into_par_iter()
                    .map(|k| {
                        let s = ((raw[k] - thr) / width + 0.5).clamp(0.0, 1.0);
                        let smooth = s * s * (3.0 - 2.0 * s);
                        (cont.oceanic_thickness
                            + (cont.continental_thickness - cont.oceanic_thickness) * smooth)
                            as f32
                    })
                    .collect()
            }
            None => (0..grid.cell_count())
                .into_par_iter()
                .map(|k| {
                    let d = grid.center(grid.cell(k));
                    let dot = d[0] * axis[0] + d[1] * axis[1] + d[2] * axis[2];
                    let dichotomy = 0.5 * self.dichotomy_amplitude * (3.0 * dot).tanh();
                    (self.crust_mean + dichotomy + self.crust_sigma * raw[k] / rms) as f32
                })
                .collect(),
        };
        let mean_crust = crust.iter().map(|v| *v as f64).sum::<f64>() / crust.len() as f64;

        let bedrock: Vec<f32> = (0..grid.cell_count())
            .into_par_iter()
            .map(|k| {
                let d = grid.center(grid.cell(k));
                let isostatic = (crust[k] as f64 - mean_crust) * self.isostatic_factor;
                let detail = self.relief * self.noise.fbm(d, 3.0, 6, 2.0, 0.5);
                (isostatic + detail) as f32
            })
            .collect();

        layers
            .crust_thickness
            .as_mut_slice()
            .copy_from_slice(&crust);
        layers.bedrock.as_mut_slice().copy_from_slice(&bedrock);
        layers.unit_id = Field::filled(grid, UNIT_PRIMORDIAL);
        layers
    }

    /// Advance the state to time `t` (Gyr since formation) over an interval `dt`.
    pub fn step(&self, layers: &mut LayerSet, t: f64, dt: f64) {
        let grid = self.grid;
        let frac = dt / self.age;
        let progress = t / self.age;
        let relief = self.relief;
        // Blend weights that morph the pattern from one noise field to the other over the age.
        // The pattern changes in place: it is accumulated into the layers every step, so any
        // drift would smear into streaks.
        let phase = 0.5 * std::f64::consts::PI * progress;
        let (wa, wb) = (phase.cos().powi(2), phase.sin().powi(2));
        let morph = |d: [f64; 3], freq: f64, octaves: u32| {
            wa * self.noise.fbm(d, freq, octaves, 2.0, 0.5)
                + wb * self.noise_b.fbm(d, freq, octaves, 2.0, 0.5)
        };

        // Slow large-scale perturbation that morphs with time.
        let uplift: Vec<f32> = (0..grid.cell_count())
            .into_par_iter()
            .map(|k| {
                let d = grid.center(grid.cell(k));
                (frac * 0.6 * relief * morph(d, 2.5, 4)) as f32
            })
            .collect();
        layers
            .bedrock
            .as_mut_slice()
            .par_iter_mut()
            .zip(uplift.par_iter())
            .for_each(|(h, u)| *h += u);

        // Patchy resurfacing: where a morphing higher-frequency field is high, the surface is
        // renewed: the unit changes and the surface age is stamped with the current time.
        let threshold = 0.42f64;
        let marks: Vec<bool> = (0..grid.cell_count())
            .into_par_iter()
            .map(|k| {
                let d = grid.center(grid.cell(k));
                morph(d, 6.0, 3) > threshold
            })
            .collect();
        let age_slice = layers.surface_age.as_mut_slice();
        let unit_slice = layers.unit_id.as_mut_slice();
        age_slice
            .par_iter_mut()
            .zip(unit_slice.par_iter_mut())
            .zip(marks.par_iter())
            .for_each(|((a, u), &m)| {
                if m {
                    *a = t as f32;
                    *u = UNIT_RESURFACED;
                }
            });

        // Diffusion standing in for erosion.
        relax(&mut layers.bedrock, 0.05);
    }
}

/// Rotate a vector about a unit axis by an angle in radians (Rodrigues).
pub fn rotate(d: [f64; 3], axis: [f64; 3], angle: f64) -> [f64; 3] {
    let (s, c) = angle.sin_cos();
    let dot = axis[0] * d[0] + axis[1] * d[1] + axis[2] * d[2];
    let cross = [
        axis[1] * d[2] - axis[2] * d[1],
        axis[2] * d[0] - axis[0] * d[2],
        axis[0] * d[1] - axis[1] * d[0],
    ];
    [
        d[0] * c + cross[0] * s + axis[0] * dot * (1.0 - c),
        d[1] * c + cross[1] * s + axis[1] * dot * (1.0 - c),
        d[2] * c + cross[2] * s + axis[2] * dot * (1.0 - c),
    ]
}

/// One explicit relaxation step: h += k (mean of the four edge neighbours - h).
pub fn relax(field: &mut Field<f32>, k: f32) {
    let grid = field.grid();
    let src = field.as_slice().to_vec();
    field
        .as_mut_slice()
        .par_iter_mut()
        .enumerate()
        .for_each(|(idx, h)| {
            let c = grid.cell(idx);
            let mut sum = 0.0f32;
            for (di, dj) in [(1, 0), (-1, 0), (0, 1), (0, -1)] {
                sum += src[grid.index(grid.neighbor(c, di, dj))];
            }
            *h += k * (0.25 * sum - *h);
        });
}
