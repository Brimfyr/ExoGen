// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! Wrinkle ridges on volcanic plains (PLAN.md 4.4).
//!
//! A pile of flood lava subsides as it cools and compacts, most where it is thickest, and
//! the basin it fills sags under its load. The plains shorten over the resulting slopes:
//! around the margins of a filled basin and over buried rims, where the fill thins
//! quickly, long sinuous ridges form roughly along the contours of the lava thickness,
//! with a broad arch and a narrow crenulated crest. They form within a few hundred
//! million years of the flows and are cut by whatever comes later.

use std::collections::HashMap;

use crate::constants::*;
use crate::crater::{CellWrite, StampOutput, Stamper};
use crate::fields::LayerSet;
use crate::geom::{cross, dot, normalize};
use crate::grid::CubeSphere;
use crate::rng::Rng;
use crate::stages::volcanism::is_lava_unit;

pub const VERSION: u64 = 1;

/// A wrinkle ridge, kept for the sidecar and for stamping at the finishing resolution.
#[derive(Clone, Debug)]
pub struct RidgeRecord {
    /// The trace, evenly spaced points along the crest.
    pub points: Vec<[f64; 3]>,
    /// Spacing of the trace points, metres.
    pub spacing: f64,
    pub length: f64,
    /// Height of the arch and of the crest riding on it, metres.
    pub height: f64,
    pub crest: f64,
    pub time_gyr: f64,
    pub id: u64,
}

impl RidgeRecord {
    pub fn midpoint(&self) -> [f64; 3] {
        self.points[self.points.len() / 2]
    }

    /// Angular radius of the cap holding the ridge and its arch.
    pub fn cap_radius(&self, body_radius: f64) -> f64 {
        (0.5 * self.length + 3.0 * RIDGE_ARCH_WIDTH) / body_radius
    }
}

/// Height at signed across-strike distance `d` from the crest line: a broad arch and a
/// narrow crest, the crest broadened by `blur` for degradation and grid smoothing.
pub fn ridge_profile(d: f64, height: f64, crest: f64, blur: f64) -> f64 {
    let arch = height * (-(d / RIDGE_ARCH_WIDTH).powi(2)).exp();
    let w = RIDGE_CREST_WIDTH.hypot(blur);
    let top = crest * RIDGE_CREST_WIDTH / w * (-(d / w).powi(2)).exp();
    arch + top
}

/// Gradient of a coarse field at a cell, as a tangent vector per metre.
fn gradient(grid: &CubeSphere, body_radius: f64, field: &[f32], idx: usize) -> [f64; 3] {
    let cell = grid.cell(idx);
    let p = grid.center(cell);
    let mut g = [0.0; 3];
    for (di, dj) in [(1, 0), (-1, 0), (0, 1), (0, -1)] {
        let n = grid.index(grid.neighbor(cell, di, dj));
        let q = grid.center(grid.cell(n));
        let v = [q[0] - p[0], q[1] - p[1], q[2] - p[2]];
        let dist2 = dot(v, v) * body_radius * body_radius;
        if dist2 <= 0.0 {
            continue;
        }
        let dv = (field[n] - field[idx]) as f64;
        // Half the difference along each axis pair, in the tangent plane.
        g[0] += 0.5 * dv * v[0] * body_radius / dist2;
        g[1] += 0.5 * dv * v[1] * body_radius / dist2;
        g[2] += 0.5 * dv * v[2] * body_radius / dist2;
    }
    let radial = dot(g, p);
    [
        g[0] - radial * p[0],
        g[1] - radial * p[1],
        g[2] - radial * p[2],
    ]
}

/// Lay out the ridges of the plains at the end of the loop from the cumulative lava
/// thickness and the time of the last flood at each cell. Seeds favour steep thickness
/// gradients; each ridge follows the thickness contours with some inertia and wobble,
/// and stops where the plains end.
pub fn plan_ridges(
    grid: &CubeSphere,
    body_radius: f64,
    units: &[u16],
    thickness: &[f32],
    last_time: &[f32],
    age_gyr: f64,
    seed: u64,
) -> Vec<RidgeRecord> {
    let rng = Rng::new(seed ^ 0x51D6_E5A1);
    let n = grid.cell_count();
    let cell_m = grid.delta() * body_radius;
    // Seed weights: lava cells, by gradient of the fill over a floor, so quiet plains
    // interiors still get some ridges.
    let mut weights = vec![0.0f64; n];
    let mut lava_area = 0.0;
    let mut grads = vec![0.0f64; n];
    // Plains proper: lava that was laid down here, not crater floors or ejecta on it.
    let plains = |k: usize| units[k] == crate::units::UNIT_LAVA && last_time[k].is_finite();
    for (k, gk) in grads.iter_mut().enumerate() {
        if !plains(k) {
            continue;
        }
        let a = grid.cell_area(grid.cell(k)) * body_radius * body_radius;
        lava_area += a;
        let g = gradient(grid, body_radius, thickness, k);
        *gk = dot(g, g).sqrt();
    }
    if lava_area <= 0.0 {
        return Vec::new();
    }
    let mean_grad = {
        let (mut s, mut c) = (0.0, 0usize);
        for (k, g) in grads.iter().enumerate() {
            if plains(k) {
                s += g;
                c += 1;
            }
        }
        if c > 0 { s / c as f64 } else { 0.0 }
    };
    let mut cumulative = Vec::with_capacity(n);
    let mut acc = 0.0;
    for k in 0..n {
        if plains(k) {
            let rel = if mean_grad > 0.0 {
                grads[k] / mean_grad
            } else {
                1.0
            };
            weights[k] = RIDGE_SEED_FLOOR + rel;
        }
        acc += weights[k];
        cumulative.push(acc);
    }
    let total_weight = acc;
    let target_length = RIDGE_LENGTH_PER_AREA * lava_area;
    let spacing = 0.5 * cell_m;
    let (l_min, l_max) = RIDGE_LENGTH;
    let mut out = Vec::new();
    let mut laid = 0.0;
    let mut j = 0u64;
    while laid < target_length && j < 200_000 {
        let draw = 8 * j;
        let length = l_min * (l_max / l_min).powf(rng.uniform(0, draw));
        let u = rng.uniform(0, draw + 1) * total_weight;
        let start = cumulative.partition_point(|&c| c < u).min(n - 1);
        let height = RIDGE_HEIGHT * (0.6 + 0.8 * rng.uniform(0, draw + 2));
        let crest = RIDGE_CREST_HEIGHT * (0.6 + 0.8 * rng.uniform(0, draw + 3));
        // Initial heading: along the contour, either way.
        let p0 = grid.center(grid.cell(start));
        let g0 = gradient(grid, body_radius, thickness, start);
        let mut heading = if dot(g0, g0) > 0.0 {
            normalize(cross(p0, g0))
        } else {
            let (east, north) = crate::geom::local_basis(p0);
            let phi = 2.0 * std::f64::consts::PI * rng.uniform(0, draw + 4);
            normalize([
                east[0] * phi.cos() + north[0] * phi.sin(),
                east[1] * phi.cos() + north[1] * phi.sin(),
                east[2] * phi.cos() + north[2] * phi.sin(),
            ])
        };
        if rng.uniform(0, draw + 5) < 0.5 {
            heading = [-heading[0], -heading[1], -heading[2]];
        }
        let steps = (length / spacing).ceil().max(2.0) as usize;
        let mut points = vec![p0];
        let mut p = p0;
        let wobble_phase = 2.0 * std::f64::consts::PI * rng.uniform(0, draw + 6);
        for i in 1..=steps {
            // Follow the contour with inertia, plus a slow wobble.
            let idx = grid.index(grid.locate(p).0);
            if !is_lava_unit(units[idx]) {
                break;
            }
            let g = gradient(grid, body_radius, thickness, idx);
            let mut along = cross(p, g);
            if dot(along, heading) < 0.0 {
                along = [-along[0], -along[1], -along[2]];
            }
            let gn = dot(g, g).sqrt();
            let follow = if mean_grad > 0.0 {
                (gn / mean_grad).min(1.0) * RIDGE_CONTOUR_FOLLOWING
            } else {
                0.0
            };
            let along = if gn > 0.0 { normalize(along) } else { heading };
            let mut h = [
                (1.0 - follow) * heading[0] + follow * along[0],
                (1.0 - follow) * heading[1] + follow * along[1],
                (1.0 - follow) * heading[2] + follow * along[2],
            ];
            // Wobble: a slow sinusoidal turn about the normal.
            let turn = RIDGE_WOBBLE_DEG.to_radians()
                * (2.0 * std::f64::consts::PI * i as f64 * spacing / RIDGE_WOBBLE_LENGTH
                    + wobble_phase)
                    .sin()
                * spacing
                / RIDGE_WOBBLE_LENGTH
                * 2.0
                * std::f64::consts::PI;
            let side = cross(p, h);
            h = [
                h[0] + turn * side[0],
                h[1] + turn * side[1],
                h[2] + turn * side[2],
            ];
            let radial = dot(h, p);
            heading = normalize([
                h[0] - radial * p[0],
                h[1] - radial * p[1],
                h[2] - radial * p[2],
            ]);
            let step = spacing / body_radius;
            p = normalize([
                p[0] + step * heading[0],
                p[1] + step * heading[1],
                p[2] + step * heading[2],
            ]);
            points.push(p);
        }
        j += 1;
        if points.len() < 4 {
            continue;
        }
        let actual = (points.len() - 1) as f64 * spacing;
        // The ridge stands some time after the last flow it crosses.
        let last = points
            .iter()
            .map(|p| last_time[grid.index(grid.locate(*p).0)])
            .filter(|t| t.is_finite())
            .fold(f32::NEG_INFINITY, f32::max);
        if !last.is_finite() {
            continue;
        }
        let time_gyr = (last as f64 + RIDGE_DELAY_GYR).min(age_gyr);
        laid += actual;
        out.push(RidgeRecord {
            points,
            spacing,
            length: actual,
            height,
            crest,
            time_gyr,
            id: j,
        });
    }
    out.sort_by(|a, b| a.time_gyr.total_cmp(&b.time_gyr).then(a.id.cmp(&b.id)));
    out
}

/// The height one ridge adds, with the crest broadened by `blur` for degradation.
pub fn ridge_compute(
    layers: &LayerSet,
    grid: &CubeSphere,
    body_radius: f64,
    rec: &RidgeRecord,
    blur: f64,
    stamper: &mut Stamper,
) -> StampOutput {
    let bedrock = layers.bedrock.as_slice();
    let cell = grid.delta() * body_radius;
    let width = 3.0 * RIDGE_ARCH_WIDTH + 2.0 * cell;
    let half_band = 0.5 * rec.spacing + cell;
    let mut best: HashMap<usize, (f64, f64)> = HashMap::new();
    let n = rec.points.len();
    let noise = crate::noise::GradientNoise::new(rec.id ^ 0x51D6_0000);
    for k in 0..n {
        let p = rec.points[k];
        let prev = rec.points[k.saturating_sub(1)];
        let next = rec.points[(k + 1).min(n - 1)];
        let chord = [next[0] - prev[0], next[1] - prev[1], next[2] - prev[2]];
        let radial = dot(chord, p);
        let t = normalize([
            chord[0] - radial * p[0],
            chord[1] - radial * p[1],
            chord[2] - radial * p[2],
        ]);
        if dot(t, t) < 0.5 {
            continue;
        }
        let normal = normalize(cross(p, t));
        let along_trace = k as f64 * rec.spacing;
        let taper = (std::f64::consts::PI * along_trace / rec.length)
            .sin()
            .max(0.0)
            .sqrt();
        if taper <= 0.0 {
            continue;
        }
        // The crest wanders within the arch.
        let wander = RIDGE_CREST_WANDER
            * RIDGE_ARCH_WIDTH
            * noise.fbm(p, body_radius / RIDGE_CREST_WANDER_LENGTH, 2, 2.0, 0.5);
        stamper.gather(grid, body_radius, p, width);
        for &(idx, _, _) in &stamper.footprint {
            let d = grid.center(grid.cell(idx));
            let along = dot(d, t) * body_radius;
            if along.abs() > half_band {
                continue;
            }
            let across = dot(d, normal) * body_radius;
            let u = taper
                * (rec.height * (-(across / RIDGE_ARCH_WIDTH).powi(2)).exp()
                    + ridge_profile(across - wander, 0.0, rec.crest, blur));
            let entry = best.entry(idx).or_insert((f64::INFINITY, 0.0));
            if along.abs() < entry.0 {
                *entry = (along.abs(), u);
            }
        }
    }
    let mut cells: Vec<CellWrite> = best
        .into_iter()
        .map(|(idx, (_, u))| CellWrite {
            idx: idx as u32,
            bedrock: bedrock[idx] + u as f32,
            surface: false,
            age: 0.0,
            unit: 0,
            deep: 0.0,
            melt: 0.0,
            cover: 0.0,
        })
        .collect();
    cells.sort_by_key(|c| c.idx);
    StampOutput {
        cells,
        composition: None,
    }
}

/// Total ridge length per unit of plains area for a laid-out set, 1/m.
pub fn ridge_density(ridges: &[RidgeRecord], lava_area: f64) -> f64 {
    if lava_area <= 0.0 {
        return 0.0;
    }
    ridges.iter().map(|r| r.length).sum::<f64>() / lava_area
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::angle_between;
    use crate::units::{UNIT_LAVA, UNIT_PRIMORDIAL};

    #[test]
    fn the_profile_is_an_arch_with_a_crest() {
        let h = ridge_profile(0.0, 200.0, 60.0, 0.0);
        assert!((h - 260.0).abs() < 1e-6);
        let far = ridge_profile(3.0 * RIDGE_ARCH_WIDTH, 200.0, 60.0, 0.0);
        assert!(far < 1.0, "{far}");
        // Blur spreads the crest without changing its volume much.
        let sharp: f64 = (-50..=50)
            .map(|i| ridge_profile(i as f64 * 100.0, 0.0, 60.0, 0.0) * 100.0)
            .sum();
        let soft: f64 = (-50..=50)
            .map(|i| ridge_profile(i as f64 * 100.0, 0.0, 60.0, 2_000.0) * 100.0)
            .sum();
        assert!((sharp - soft).abs() < 0.05 * sharp, "{sharp} vs {soft}");
        assert!(ridge_profile(0.0, 0.0, 60.0, 2_000.0) < 60.0);
    }

    /// A round patch of plains thickest in the middle: the ridges lie on it, run roughly
    /// along its contours, and together reach the planned length per area.
    #[test]
    fn ridges_follow_the_contours_of_a_filled_basin() {
        let grid = CubeSphere::new(64);
        let body_radius = 1_700e3;
        let n = grid.cell_count();
        let centre = normalize([1.0, 0.3, 0.2]);
        let mut units = vec![UNIT_PRIMORDIAL; n];
        let mut thickness = vec![0.0f32; n];
        let last_time = vec![3.0f32; n];
        let radius = 0.25;
        for k in 0..n {
            let a = angle_between(centre, grid.center(grid.cell(k)));
            if a < radius {
                units[k] = UNIT_LAVA;
                thickness[k] = (2_000.0 * (1.0 - (a / radius).powi(2))) as f32;
            }
        }
        let ridges = plan_ridges(&grid, body_radius, &units, &thickness, &last_time, 4.5, 7);
        assert!(ridges.len() > 10, "{} ridges", ridges.len());
        let lava_area: f64 = (0..n)
            .filter(|&k| units[k] == UNIT_LAVA)
            .map(|k| grid.cell_area(grid.cell(k)) * body_radius * body_radius)
            .sum();
        let density = ridge_density(&ridges, lava_area);
        assert!(
            density > 0.7 * RIDGE_LENGTH_PER_AREA && density < 1.3 * RIDGE_LENGTH_PER_AREA,
            "{density} vs {RIDGE_LENGTH_PER_AREA}"
        );
        // Along a contour the angular distance from the centre stays nearly constant.
        let mut concentric = 0;
        for r in &ridges {
            let d0 = angle_between(centre, r.points[0]);
            let d1 = angle_between(centre, *r.points.last().unwrap());
            let along = r.length / body_radius;
            if (d1 - d0).abs() < 0.5 * along {
                concentric += 1;
            }
            for p in &r.points {
                assert!(angle_between(centre, *p) < radius + 2.0 * grid.delta());
            }
            assert!((r.time_gyr - (3.0 + RIDGE_DELAY_GYR)).abs() < 1e-6);
        }
        assert!(
            concentric * 10 > ridges.len() * 7,
            "{concentric} of {}",
            ridges.len()
        );
    }
}
