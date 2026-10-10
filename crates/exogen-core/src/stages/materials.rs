// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! Phase 1 placeholder material model: an albedo per cell from the composition class, the
//! height, and the liquid layer. The spectral end-member model (PLAN.md 4.10) replaces it
//! in Phase 2. Colours are linear reflectance.

use rayon::prelude::*;

use crate::fields::{Field, LayerSet};
use crate::params::ParamGraph;

pub const VERSION: u64 = 1;

fn mix(a: [f32; 3], b: [f32; 3], t: f32) -> [f32; 3] {
    [
        a[0] + (b[0] - a[0]) * t,
        a[1] + (b[1] - a[1]) * t,
        a[2] + (b[2] - a[2]) * t,
    ]
}

fn scale(a: [f32; 3], s: f32) -> [f32; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

/// Linear-reflectance albedo of the surface.
pub fn albedo(graph: &ParamGraph, layers: &LayerSet, height: &Field<f32>) -> Field<[f32; 3]> {
    let base: [f32; 3] = match graph.str("body.composition_class") {
        "iron_rich" => [0.11, 0.11, 0.11],
        "rock" => [0.15, 0.13, 0.11],
        "rock_ice" => [0.50, 0.52, 0.55],
        _ => [0.62, 0.64, 0.68],
    };
    let (deep, shallow): ([f32; 3], [f32; 3]) = match graph.str("volatiles.working_liquid") {
        "water" | "ammonia_water" => ([0.012, 0.035, 0.085], [0.06, 0.22, 0.26]),
        "methane_ethane" => ([0.02, 0.02, 0.015], [0.07, 0.07, 0.05]),
        _ => ([0.03, 0.03, 0.03], [0.08, 0.08, 0.08]),
    };
    let (lo, hi) = height.min_max();
    let span = (hi - lo).max(1.0);
    let grid = layers.grid;
    let h = height.as_slice();
    let liquid = layers.liquid.as_slice();
    let age = layers.surface_age.as_slice();
    let data: Vec<[f32; 3]> = (0..grid.cell_count())
        .into_par_iter()
        .map(|k| {
            if liquid[k] > 0.0 {
                let t = (-liquid[k] / 60.0).exp();
                mix(deep, shallow, t)
            } else {
                let t = (h[k] - lo) / span;
                let brightness = 0.75 + 0.5 * t;
                // Younger surfaces are a little darker, standing in for weathering state.
                let youth = 1.0 - 0.15 * (age[k] > 0.0) as u8 as f32;
                scale(base, brightness * youth)
            }
        })
        .collect();
    let mut f = Field::new(grid);
    f.as_mut_slice().copy_from_slice(&data);
    f
}
