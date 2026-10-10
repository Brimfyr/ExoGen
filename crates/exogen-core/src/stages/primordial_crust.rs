// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! Primordial crust (PLAN.md 4.1): the initial layered state at formation.
//!
//! Crustal thickness is a low-degree spherical-harmonic field with an optional hemispheric
//! dichotomy term, or a two-level crust for mobile-lid bodies until plate tectonics
//! exists. Elevation follows from Airy isostasy, with a seeded noise term scaled by the
//! relief limit standing in for the finer structure that later processes produce.

use rayon::prelude::*;

use crate::constants::{DENSITY_CRUST_ICE, DENSITY_CRUST_ROCK};
use crate::fields::{Field, LayerSet};
use crate::grid::CubeSphere;
use crate::noise::GradientNoise;
use crate::params::ParamGraph;
use crate::rng::Rng;
use crate::spharm::Harmonics;
use crate::units::{UNIT_CONTINENTAL, UNIT_OCEANIC, UNIT_PRIMORDIAL};

/// Bumped when the behaviour changes, so caches are not reused across versions.
pub const VERSION: u64 = 4;

/// Two-level crust for mobile-lid bodies: thin oceanic crust and thick continental crust
/// covering a fraction of the area, which gives the bimodal hypsometry that lets an ocean
/// fill the basins and leave continents standing. Plate tectonics replaces it in Phase 6.
struct Continents {
    oceanic_thickness: f64,
    continental_thickness: f64,
    area_fraction: f64,
}

pub struct PrimordialCrust {
    grid: CubeSphere,
    noise: GradientNoise,
    harmonics: Harmonics,
    dichotomy_axis: [f64; 3],
    dichotomy_amplitude: f64,
    crust_mean: f64,
    crust_sigma: f64,
    continents: Option<Continents>,
    relief: f64,
    isostatic_factor: f64,
}

impl PrimordialCrust {
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
            harmonics: Harmonics::random(&rng, 2, 8, 2.0),
            dichotomy_axis: rng.unit_direction(1, 0),
            dichotomy_amplitude: graph.f64("tectonics.dichotomy_amplitude"),
            crust_mean: graph.f64("tectonics.initial_crust_thickness"),
            crust_sigma: 0.12 * graph.f64("tectonics.initial_crust_thickness"),
            continents,
            relief,
            isostatic_factor: (rho_mantle - rho_crust) / rho_mantle,
        }
    }

    /// Initial layered state at formation.
    pub fn init(&self) -> LayerSet {
        let grid = self.grid;
        let mut layers = LayerSet::new(grid);

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

        let units: Vec<u16> = match &self.continents {
            Some(cont) => {
                let mid = 0.5 * (cont.oceanic_thickness + cont.continental_thickness) as f32;
                crust
                    .iter()
                    .map(|&c| {
                        if c > mid {
                            UNIT_CONTINENTAL
                        } else {
                            UNIT_OCEANIC
                        }
                    })
                    .collect()
            }
            None => vec![UNIT_PRIMORDIAL; grid.cell_count()],
        };

        layers
            .crust_thickness
            .as_mut_slice()
            .copy_from_slice(&crust);
        layers.bedrock.as_mut_slice().copy_from_slice(&bedrock);
        layers.unit_id = Field::new(grid);
        layers.unit_id.as_mut_slice().copy_from_slice(&units);
        layers
    }
}
