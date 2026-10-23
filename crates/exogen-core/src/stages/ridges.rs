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
//!
//! The plains stand in the lid's stress field (F156): the shrinking of a whole lid
//! squeezes them, their own subsidence brings them past their strength, and the field
//! adds its compression in the lows and at the feet of the rises and takes it away on the
//! highs, so the ridges gather where the plains are squeezed most and are absent where the
//! field holds them in tension (Venus's ridged plains, 93 per cent of them below its mean
//! radius, Bilotti and Suppe 1999), and strike across the field's compression where it
//! decides their strike, around Tharsis as around Venus's highlands, following the lava's
//! contours where it does not.

use std::collections::HashMap;

use rayon::prelude::*;

use crate::constants::*;
use crate::crater::{CellWrite, StampOutput, Stamper};
use crate::fields::LayerSet;
use crate::geom::{cross, dot, normalize};
use crate::grid::CubeSphere;
use crate::params::ParamGraph;
use crate::rng::Rng;
use crate::stages::stress::{FieldHistory, Squeeze, shrinking};
use crate::stages::tectonics::Systems;
use crate::stages::volcanism::is_lava_unit;
use crate::thermal::ThermalHistory;

pub const VERSION: u64 = 4;

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
        (0.5 * self.length + 3.0 * RIDGE_ARCH_WIDTH * self.size()) / body_radius
    }

    /// The ridge's size against a ridge of RIDGE_HEIGHT, which its arch and crest widen
    /// with (F159): a higher ridge stands over a thrust rooted deeper.
    pub fn size(&self) -> f64 {
        (self.height / RIDGE_HEIGHT).max(RIDGE_MIN_SIZE)
    }
}

/// A ridge's height for its length (F159): RIDGE_HEIGHT at RIDGE_REFERENCE_LENGTH, growing
/// with the length as a thrust's displacement does (the scarps' exponent, F155), scattered
/// by `u` from 0.6 to 1.4 of that.
pub fn ridge_height(length: f64, u: f64) -> f64 {
    RIDGE_HEIGHT
        * (length / RIDGE_REFERENCE_LENGTH).powf(SCARP_HEIGHT_EXPONENT)
        * (0.6 + 0.8 * u.clamp(0.0, 1.0))
}

/// Height at signed across-strike distance `d` from the crest line: a broad arch and a
/// narrow crest, both widened by the ridge's `size` (F159), the crest broadened by `blur`
/// for degradation and grid smoothing.
pub fn ridge_profile(d: f64, height: f64, crest: f64, size: f64, blur: f64) -> f64 {
    let arch_width = RIDGE_ARCH_WIDTH * size;
    let crest_width = RIDGE_CREST_WIDTH * size;
    let arch = height * (-(d / arch_width).powi(2)).exp();
    let w = crest_width.hypot(blur);
    let top = crest * crest_width / w * (-(d / w).powi(2)).exp();
    arch + top
}

/// How far the lid squeezes the plains at each cell of the loop's grid, read at the time
/// the cell's ridges would stand by the field of the loop's last epoch then (F156, F157),
/// and the strike the field gives there. The squeeze is the overstress, the compression past
/// the lid's strength: what a whole lid's shrinking puts on it, the plains standing at their
/// strength at least, brought there by drives the field does not hold (a stated choice, with
/// nothing fitted: Venus's climate warming its plains after their renewal, Solomon et al.
/// 1999, which the model's atmosphere does not yet carry, and the plains' own cooling), and
/// the field's compression, the floods' own loads among it, where both horizontal stresses
/// press harder than the overburden; given as a share that saturates at the lid's strength,
/// past which the ridges stand as densely as they can and the shortening goes into their
/// slip. So the field decides where the plains fold and where they do not.
pub fn plains_squeeze(
    graph: &ParamGraph,
    thermal: &ThermalHistory,
    history: &FieldHistory,
    grid: CubeSphere,
    last_time: &[f32],
) -> (Vec<f64>, Vec<Squeeze>) {
    let whole = graph.str("interior.tectonic_regime") == "stagnant_lid";
    let shrink = |t: f64| {
        if whole {
            shrinking(thermal, history.body_radius, t)
        } else {
            0.0
        }
    };
    squeeze_plains(history, grid, last_time, graph.f64("body.age"), &shrink)
}

/// The plains' squeeze at each cell given the body's age, Gyr, and what its shrinking puts
/// on the lid at a time, Pa.
pub fn squeeze_plains(
    history: &FieldHistory,
    grid: CubeSphere,
    last_time: &[f32],
    age: f64,
    shrink: &(dyn Fn(f64) -> f64 + Sync),
) -> (Vec<f64>, Vec<Squeeze>) {
    let n = grid.cell_count();
    // Each plains cell at its time, grouped by the epochs read by then.
    let mut cells: Vec<(usize, usize, f64)> = (0..n)
        .filter(|&k| last_time[k].is_finite())
        .map(|k| {
            let t = (last_time[k] as f64 + RIDGE_DELAY_GYR).min(age);
            (history.read_by(t), k, t)
        })
        .collect();
    cells.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
    let mut squeezed = vec![0.0; n];
    let mut squeeze = vec![Squeeze::default(); n];
    let mut from = 0;
    while from < cells.len() {
        let epochs = cells[from].0;
        let to = cells[from..]
            .iter()
            .position(|c| c.0 != epochs)
            .map_or(cells.len(), |d| from + d);
        let field = history.field(cells[from].2);
        let strength = field.yield_compression;
        let read: Vec<(usize, f64, Squeeze)> = cells[from..to]
            .par_iter()
            .map(|&(_, k, t)| {
                let q = field.squeeze(grid.center(grid.cell(k)), t);
                // The plains stand at their strength at least, brought there by drives the
                // field does not hold (a stated choice): Venus's climate warming its plains
                // after their renewal (Solomon et al. 1999), the plains' own cooling.
                let squeeze = shrink(t).max(strength);
                let drive = squeeze - strength;
                // A ridge's thrust where both horizontal stresses press harder than the
                // overburden (F157).
                let over = if q.least <= squeeze {
                    (drive + q.stress).max(0.0)
                } else {
                    0.0
                };
                (k, over / (over + strength.max(1.0)), q)
            })
            .collect();
        for (k, share, q) in read {
            squeezed[k] = share;
            squeeze[k] = q;
        }
        from = to;
    }
    (squeezed, squeeze)
}

/// Lay out the ridges of the plains at the end of the loop from the time of the last flood
/// at each cell and the lid's squeeze on it (F156, F157). Seeds stand as the plains are
/// squeezed; half run beside a ridge drawn before near where they start, the thrusts' ridges
/// among them (`existing`), en-echelon along a span of it, as the maria's ridges stand in
/// systems (F160); the rest strike across the field's compression as far as the field
/// decides it and run their own way with some inertia where it does not, meandering about
/// their course; each stops where the plains end or are no longer squeezed, or, running its
/// own way, where it meets a ridge drawn before at least as long.
#[allow(clippy::too_many_arguments)]
pub fn plan_ridges(
    grid: &CubeSphere,
    body_radius: f64,
    units: &[u16],
    last_time: &[f32],
    squeezed: &[f64],
    squeeze: &[Squeeze],
    existing: &[RidgeRecord],
    age_gyr: f64,
    seed: u64,
) -> Vec<RidgeRecord> {
    let rng = Rng::new(seed ^ 0x51D6_E5A1);
    let n = grid.cell_count();
    let cell_m = grid.delta() * body_radius;
    // Plains proper: lava that was laid down here, not crater floors or ejecta on it, that
    // the lid squeezes; seeds as they are squeezed.
    let plains = |k: usize| {
        units[k] == crate::units::UNIT_LAVA && last_time[k].is_finite() && squeezed[k] > 0.0
    };
    let mut lava_area = 0.0;
    let mut cumulative = Vec::with_capacity(n);
    let mut acc = 0.0;
    for (k, &share) in squeezed.iter().enumerate() {
        if plains(k) {
            lava_area += grid.cell_area(grid.cell(k)) * body_radius * body_radius;
            acc += share;
        }
        cumulative.push(acc);
    }
    if lava_area <= 0.0 || acc <= 0.0 {
        return Vec::new();
    }
    let total_weight = acc;
    let target_length = RIDGE_LENGTH_PER_AREA * lava_area;
    let spacing = 0.5 * cell_m;
    let (l_min, l_max) = RIDGE_LENGTH;
    let reach = |height: f64| 3.0 * RIDGE_ARCH_WIDTH * (height / RIDGE_HEIGHT).max(RIDGE_MIN_SIZE);
    let mut systems = Systems::new(SYSTEM_INDEX_WIDTH, body_radius);
    for r in existing {
        systems.add(r.points.clone(), r.spacing, r.length, reach(r.height));
    }
    let on_plains = |p: [f64; 3]| {
        let idx = grid.index(grid.locate(p).0);
        is_lava_unit(units[idx]) && squeezed[idx] > 0.0
    };
    let mut out = Vec::new();
    let mut laid = 0.0;
    let mut j = 0u64;
    while laid < target_length && j < 200_000 {
        let draw = 16 * j;
        // Lengths from the thrusts' power law (F155, F159), most of them short.
        let length = crate::stages::tectonics::sample_length(rng.uniform(0, draw), l_min, l_max);
        let u = rng.uniform(0, draw + 1) * total_weight;
        let start = cumulative.partition_point(|&c| c < u).min(n - 1);
        // Its first heading: across the field's compression as far as it decides the
        // strike, its own way otherwise, either way along it.
        let p0 = grid.center(grid.cell(start));
        let (east, north) = crate::geom::local_basis(p0);
        let phi = 2.0 * std::f64::consts::PI * rng.uniform(0, draw + 4);
        let free = [
            east[0] * phi.cos() + north[0] * phi.sin(),
            east[1] * phi.cos() + north[1] * phi.sin(),
            east[2] * phi.cos() + north[2] * phi.sin(),
        ];
        let q = squeeze[start];
        let strike = if dot(q.strike, free) < 0.0 {
            [-q.strike[0], -q.strike[1], -q.strike[2]]
        } else {
            q.strike
        };
        let mut heading = normalize([
            q.decided * strike[0] + (1.0 - q.decided) * free[0],
            q.decided * strike[1] + (1.0 - q.decided) * free[1],
            q.decided * strike[2] + (1.0 - q.decided) * free[2],
        ]);
        if rng.uniform(0, draw + 5) < 0.5 {
            heading = [-heading[0], -heading[1], -heading[2]];
        }
        // Half beside a ridge drawn before near where it starts (F160), as far as the
        // squeezed plains reach along it.
        let beside = systems.beside(
            p0,
            length,
            rng.uniform(0, draw + 8),
            rng.uniform(0, draw + 9),
            rng.uniform(0, draw + 10),
            rng.uniform(0, draw + 11),
        );
        let (points, trace_spacing) = if let Some((trace, host_spacing)) = beside {
            let kept: Vec<[f64; 3]> = trace.into_iter().take_while(|&q| on_plains(q)).collect();
            (kept, host_spacing)
        } else {
            let steps = (length / spacing).ceil().max(2.0) as usize;
            let mut points = vec![p0];
            let mut p = p0;
            let wobble_phase = 2.0 * std::f64::consts::PI * rng.uniform(0, draw + 6);
            let second_phase = 2.0 * std::f64::consts::PI * rng.uniform(0, draw + 7);
            // Running its own way, it stops where it meets a ridge drawn before at least as long.
            let half = SYSTEM_ARREST_SHARE * reach(ridge_height(length, 0.5));
            for i in 1..=steps {
                let idx = grid.index(grid.locate(p).0);
                if !is_lava_unit(units[idx])
                    || squeezed[idx] <= 0.0
                    || (i > 1 && systems.meets(p, half, length))
                {
                    break;
                }
                // Toward the field's strike, as far as it decides it, with inertia.
                let q = squeeze[idx];
                let strike = if dot(q.strike, heading) < 0.0 {
                    [-q.strike[0], -q.strike[1], -q.strike[2]]
                } else {
                    q.strike
                };
                let ff = q.decided * THRUST_FIELD_FOLLOWING;
                let h = [
                    (1.0 - ff) * heading[0] + ff * strike[0],
                    (1.0 - ff) * heading[1] + ff * strike[1],
                    (1.0 - ff) * heading[2] + ff * strike[2],
                ];
                let radial = dot(h, p);
                heading = normalize([
                    h[0] - radial * p[0],
                    h[1] - radial * p[1],
                    h[2] - radial * p[2],
                ]);
                // Its meander about that course: the step's direction wobbles slowly, in two
                // harmonics, which the pull of the field does not take back.
                let s = i as f64 * spacing / RIDGE_WOBBLE_LENGTH;
                let (ratio, share) = RIDGE_WOBBLE_HARMONIC;
                let wobble = RIDGE_WOBBLE_DEG.to_radians()
                    * ((2.0 * std::f64::consts::PI * s + wobble_phase).sin()
                        + share * (2.0 * std::f64::consts::PI * ratio * s + second_phase).sin());
                let side = normalize(cross(p, heading));
                let (c, w) = (wobble.cos(), wobble.sin());
                let step = spacing / body_radius;
                p = normalize([
                    p[0] + step * (c * heading[0] + w * side[0]),
                    p[1] + step * (c * heading[1] + w * side[1]),
                    p[2] + step * (c * heading[2] + w * side[2]),
                ]);
                points.push(p);
            }
            (points, spacing)
        };
        j += 1;
        if points.len() < 4 {
            continue;
        }
        let actual = (points.len() - 1) as f64 * trace_spacing;
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
        // Its height for the length it ran, its crest with it (F159).
        let height = ridge_height(actual, rng.uniform(0, draw + 2));
        let crest =
            RIDGE_CREST_HEIGHT * height / RIDGE_HEIGHT * (0.6 + 0.8 * rng.uniform(0, draw + 3));
        laid += actual;
        systems.add(points.clone(), trace_spacing, actual, reach(height));
        out.push(RidgeRecord {
            points,
            spacing: trace_spacing,
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
    let size = rec.size();
    let arch_width = RIDGE_ARCH_WIDTH * size;
    let width = 3.0 * arch_width + 2.0 * cell;
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
            * arch_width
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
                * (rec.height * (-(across / arch_width).powi(2)).exp()
                    + ridge_profile(across - wander, 0.0, rec.crest, size, blur));
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
        stain: None,
        steep: None,
        ejecta_age: None,
        compositions: None,
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
        let h = ridge_profile(0.0, 200.0, 60.0, 1.0, 0.0);
        assert!((h - 260.0).abs() < 1e-6);
        let far = ridge_profile(3.0 * RIDGE_ARCH_WIDTH, 200.0, 60.0, 1.0, 0.0);
        assert!(far < 1.0, "{far}");
        // Blur spreads the crest without changing its volume much.
        let sharp: f64 = (-50..=50)
            .map(|i| ridge_profile(i as f64 * 100.0, 0.0, 60.0, 1.0, 0.0) * 100.0)
            .sum();
        let soft: f64 = (-50..=50)
            .map(|i| ridge_profile(i as f64 * 100.0, 0.0, 60.0, 1.0, 2_000.0) * 100.0)
            .sum();
        assert!((sharp - soft).abs() < 0.05 * sharp, "{sharp} vs {soft}");
        assert!(ridge_profile(0.0, 0.0, 60.0, 1.0, 2_000.0) < 60.0);
    }

    /// The lid's squeeze places the plains' ridges and strikes them (F156): on plains
    /// squeezed in their northern half alone with the field's strike east, the ridges stand
    /// in the north and run east and west but for their wobble.
    #[test]
    fn the_lid_s_squeeze_places_the_plains_ridges_and_strikes_them() {
        let grid = CubeSphere::new(64);
        let body_radius = 1_700e3;
        let n = grid.cell_count();
        let centre = normalize([1.0, 0.0, 0.0]);
        let mut units = vec![UNIT_PRIMORDIAL; n];
        let last_time = vec![3.0f32; n];
        let mut squeezed = vec![0.0; n];
        let mut squeeze = vec![Squeeze::default(); n];
        for k in 0..n {
            let p = grid.center(grid.cell(k));
            if angle_between(centre, p) < 0.5 {
                units[k] = UNIT_LAVA;
            }
            if p[2] > 0.0 {
                squeezed[k] = 0.5;
            }
            squeeze[k] = Squeeze {
                stress: 1e6,
                least: 0.0,
                strike: crate::geom::local_basis(p).0,
                decided: 1.0,
            };
        }
        let ridges = plan_ridges(
            &grid,
            body_radius,
            &units,
            &last_time,
            &squeezed,
            &squeeze,
            &[],
            4.5,
            11,
        );
        assert!(ridges.len() > 10, "{} ridges", ridges.len());
        let (mut steps, mut east) = (0usize, 0usize);
        for r in &ridges {
            for w in r.points.windows(2) {
                assert!(w[0][2] > -grid.delta(), "a ridge in the south: {}", w[0][2]);
                let d = normalize([w[1][0] - w[0][0], w[1][1] - w[0][1], w[1][2] - w[0][2]]);
                steps += 1;
                if dot(d, crate::geom::local_basis(w[0]).0).abs() > 0.87 {
                    east += 1;
                }
            }
        }
        assert!(
            east * 10 > steps * 6,
            "{east} of {steps} steps east and west"
        );
    }

    /// A basin filled with lava bends the lid, a shell, under its own load (F157): its top
    /// is squeezed through the fill, most outward, so that a thrust strikes round the basin,
    /// and stretched outward beyond the margin; and the ridges stand within the fill at the
    /// planned length per area of the squeezed plains, after the last flow.
    #[test]
    fn a_filled_basin_s_load_squeezes_its_plains_and_stretches_its_margin() {
        use crate::stages::stress::{FieldEpoch, Relief};
        let grid = CubeSphere::new(64);
        let (body_radius, gravity, elastic) = (1_737e3, 1.62, 60e3);
        let n = grid.cell_count();
        let centre = normalize([1.0, 0.3, 0.2]);
        let radius = 0.25;
        let mut units = vec![UNIT_PRIMORDIAL; n];
        let mut thickness = vec![0.0f32; n];
        let mut last_time = vec![f32::NAN; n];
        for k in 0..n {
            let a = angle_between(centre, grid.center(grid.cell(k)));
            if a < radius {
                units[k] = UNIT_LAVA;
                thickness[k] = 2_000.0;
                last_time[k] = 3.0;
            }
        }
        let flat = Relief {
            grid,
            height: vec![0.0; n],
        };
        let history = FieldHistory {
            rises: Vec::new(),
            limbs: Default::default(),
            body_radius,
            gravity,
            epochs: vec![
                FieldEpoch {
                    flow: None,
                    time: 0.0,
                    elastic,
                    relief: flat.clone(),
                    load: None,
                },
                FieldEpoch {
                    flow: None,
                    time: 3.0,
                    elastic,
                    relief: flat,
                    load: crate::stages::load::lay(
                        grid,
                        &thickness,
                        body_radius,
                        gravity,
                        elastic,
                        3_000.0,
                        None,
                    )
                    .map(|l| l.stress),
                },
            ],
        };
        let (squeezed, squeeze) = squeeze_plains(&history, grid, &last_time, 4.5, &|_| 0.0);
        // The fill squeezes the lid's top through the basin, past its strength, and
        // stretches it outward beyond the margin (as the Hankel transform's solution of a
        // disc load does: -47 MPa at the middle, -54 at two thirds out, +44 radially at
        // 1.2 radii).
        let field = history.field(3.2);
        let (east, _) = crate::geom::local_basis(centre);
        let along = |a: f64| {
            normalize([
                centre[0] * a.cos() + east[0] * a.sin(),
                centre[1] * a.cos() + east[1] * a.sin(),
                centre[2] * a.cos() + east[2] * a.sin(),
            ])
        };
        for a in [0.0, 0.5 * radius, 0.8 * radius] {
            assert!(squeezed[grid.index(grid.locate(along(a)).0)] > 0.5, "{a}");
        }
        let out_at = |p: [f64; 3]| {
            normalize([
                p[0] * dot(centre, p) - centre[0],
                p[1] * dot(centre, p) - centre[1],
                p[2] * dot(centre, p) - centre[2],
            ])
        };
        let beyond = along(1.2 * radius);
        let (s1, _, strike, _, _) = field.at(beyond, 3.2);
        assert!(s1 > 5e6, "{s1}");
        // Stretched outward there: the normal faults' strike runs round the basin.
        assert!(
            dot(strike, out_at(beyond)).abs() < 0.5,
            "{}",
            dot(strike, out_at(beyond))
        );
        // Within the fill the shell's squeeze is greatest outward, and a thrust strikes
        // round the basin, as the maria's ridges run in rings.
        for share in [0.6, 0.85] {
            let p = along(share * radius);
            let q = field.squeeze(p, 3.2);
            assert!(
                q.decided > 0.3 && dot(q.strike, out_at(p)).abs() < 0.5,
                "{share}: {q:?}"
            );
        }
        let ridges = plan_ridges(
            &grid,
            body_radius,
            &units,
            &last_time,
            &squeezed,
            &squeeze,
            &[],
            4.5,
            7,
        );
        assert!(ridges.len() > 10, "{} ridges", ridges.len());
        let squeezed_area: f64 = (0..n)
            .filter(|&k| units[k] == UNIT_LAVA && squeezed[k] > 0.0)
            .map(|k| grid.cell_area(grid.cell(k)) * body_radius * body_radius)
            .sum();
        let density = ridge_density(&ridges, squeezed_area);
        assert!(
            density > 0.7 * RIDGE_LENGTH_PER_AREA && density < 1.3 * RIDGE_LENGTH_PER_AREA,
            "{density} vs {RIDGE_LENGTH_PER_AREA}"
        );
        for r in &ridges {
            for p in &r.points {
                assert!(angle_between(centre, *p) < radius + 2.0 * grid.delta());
            }
            assert!((r.time_gyr - (3.0 + RIDGE_DELAY_GYR)).abs() < 1e-6);
        }
    }
}
