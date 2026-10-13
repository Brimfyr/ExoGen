// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! Sub-grid crater roughness (docs/parameters.md F54). The craters smaller than the
//! finish resolves are a population too, and a surface's age says how many it has
//! taken. Their relief goes under everything as a self-affine field: anchored at the
//! smallest resolved diameter to the depth of a crater that size and to how near that
//! size stands to saturation on each cell's surface, and falling toward finer scales
//! with the Hurst exponent of cratered terrain. An old surface is rough at every scale,
//! a renewed one smooth, and the contact between them is a gradient of roughness rather
//! than a line.

use rayon::prelude::*;

use crate::constants::*;
use crate::crater::dimensions;
use crate::crater::production::{geometric_saturation, survivors};
use crate::fields::LayerSet;
use crate::noise::GradientNoise;
use crate::stages::impacts::CraterInputs;

/// Add the roughness of the unresolved craters to the bedrock, by each cell's surface
/// age: `d_min` is the smallest diameter the finish stamps, `cell` the grid's cell, m.
/// `keep` is what a cell's old fine relief has kept where a renewal touched it without
/// taking it (F60), across a renewed terrain's edge and halo: the unresolved craters
/// were taken there in the same measure as the resolved ones.
pub fn apply(
    layers: &mut LayerSet,
    body_radius: f64,
    cell: f64,
    d_min: f64,
    inputs: &CraterInputs,
    present: f64,
    keep: Option<&[f32]>,
) {
    let grid = layers.grid;
    let depth = dimensions(d_min, &inputs.scaling).depth;
    if depth <= 0.0 || present <= 0.0 {
        return;
    }
    let d_km = d_min / 1e3;
    let n_eq = EMPIRICAL_SATURATION_FRACTION * geometric_saturation(d_km);
    let per_gyr = inputs.pf.cumulative_1gyr(d_km) * inputs.flux_scale;
    // The share of saturation a surface of an age has reached at the anchor scale.
    let share_at = |age: f64| -> f64 {
        if n_eq <= 0.0 {
            return 0.0;
        }
        let n_prod = per_gyr * inputs.flux.equivalent_gyr(age.min(present), present);
        survivors(n_prod, n_eq) / n_eq
    };
    let bins = 256usize;
    let table: Vec<f64> = (0..=bins)
        .map(|i| share_at(present * i as f64 / bins as f64).sqrt())
        .collect();
    // The field: octaves from the anchor scale down to two cells, each 2^-H as strong.
    let top_frequency = body_radius / d_min;
    let octaves = (d_min / (2.0 * cell)).log2().ceil().max(1.0) as u32;
    let gain = 2f64.powf(-ROUGHNESS_HURST);
    let noise = GradientNoise::new(inputs.seed ^ 0x5B61_0A9E);
    let raw: Vec<f32> = (0..grid.cell_count())
        .into_par_iter()
        .map(|k| noise.fbm(grid.center(grid.cell(k)), top_frequency, octaves, 2.0, gain) as f32)
        .collect();
    let rms = (raw.iter().map(|v| (*v as f64).powi(2)).sum::<f64>() / raw.len().max(1) as f64)
        .sqrt()
        .max(1e-9);
    let amplitude = ROUGHNESS_RMS_SHARE * depth / rms;
    let ages = layers.surface_age.as_slice();
    // The strength from each cell's age, then smoothed over a cell or so: the edge of a
    // young deposit on old ground is a short ramp of roughness, not a staircase of cells.
    let shares: Vec<f32> = (0..grid.cell_count())
        .into_par_iter()
        .map(|k| {
            let x = (ages[k] as f64).clamp(0.0, present) / present * bins as f64;
            let i = (x.floor() as usize).min(bins - 1);
            let f = x - i as f64;
            (table[i] + (table[i + 1] - table[i]) * f) as f32 * keep.map_or(1.0, |v| v[k])
        })
        .collect();
    let shares = crate::stages::episodic::smooth(&grid, &shares, ROUGHNESS_BLEND_PASSES);
    layers
        .bedrock
        .as_mut_slice()
        .par_iter_mut()
        .enumerate()
        .for_each(|(k, h)| {
            *h += amplitude as f32 * shares[k] * raw[k];
        });
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::grid::CubeSphere;

    /// An old surface takes the roughness of its unresolved craters; a renewed one
    /// stays smooth.
    #[test]
    fn old_surfaces_are_rough_and_young_ones_smooth() {
        let (g, _) = crate::presets::preset("Moon").unwrap().to_graph();
        let inputs = CraterInputs::from_graph(&g);
        let grid = CubeSphere::new(16);
        let body_radius = g.f64("body.radius");
        let cell = grid.delta() * body_radius;
        let present = g.f64("body.age");
        let mut layers = LayerSet::new(grid);
        let n = grid.cell_count();
        for (k, a) in layers.surface_age.as_mut_slice().iter_mut().enumerate() {
            *a = if k < n / 2 {
                0.0
            } else {
                present as f32 - 0.001
            };
        }
        apply(
            &mut layers,
            body_radius,
            cell,
            4.0 * cell,
            &inputs,
            present,
            None,
        );
        let h = layers.bedrock.as_slice();
        // The strength blends over a cell or so at a contact, so the halves are judged
        // away from the face edges where old and young meet.
        let interior = |k: usize| {
            let c = grid.cell(k);
            (4..12u32).contains(&c.i) && (4..12u32).contains(&c.j)
        };
        let rms = |old_half: bool| {
            let v: Vec<f64> = (0..n)
                .filter(|&k| (k < n / 2) == old_half && interior(k))
                .map(|k| (h[k] as f64).powi(2))
                .collect();
            (v.iter().sum::<f64>() / v.len().max(1) as f64).sqrt()
        };
        let (old, young) = (rms(true), rms(false));
        let depth = dimensions(4.0 * cell, &inputs.scaling).depth;
        assert!(
            old > 0.1 * depth && old < depth,
            "old rms {old} of depth {depth}"
        );
        assert!(young < 0.05 * old, "young rms {young} against old {old}");
    }
}
