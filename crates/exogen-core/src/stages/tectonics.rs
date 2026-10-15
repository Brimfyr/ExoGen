// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! Global contraction expressed as lobate thrust scarps (PLAN.md 4.3, stagnant lid).
//!
//! The thermal history gives the cumulative radial contraction. The areal strain 2 ΔR / R
//! times the surface area is the shortening area the scarps take up: a scarp of length L
//! and height h shortens the crust by h / tan(dip) across its length. As the budget grows
//! through the loop, scarps are sampled from a power-law length distribution and stamped
//! at the step they form, so they cut what is older and are cut by what is younger. Each
//! has a steep face on the hanging-wall side and a gentle back slope, tapering to nothing
//! at its ends, along a gently curving trace.

use std::collections::HashMap;

use crate::constants::*;
use crate::crater::{CellWrite, StampOutput, Stamper};
use crate::fields::LayerSet;
use crate::geom::{azimuth_from, cross, dot, local_basis, normalize, offset};
use crate::grid::CubeSphere;
use crate::params::ParamGraph;
use crate::rng::Rng;
use crate::thermal::ThermalHistory;

pub const VERSION: u64 = 1;

/// A scarp that was emplaced, kept for the sidecar and for re-stamping its face at the
/// finishing resolution.
#[derive(Clone, Debug)]
pub struct ScarpRecord {
    /// The trace, evenly spaced points along the fault.
    pub points: Vec<[f64; 3]>,
    /// Spacing of the trace points, metres.
    pub spacing: f64,
    pub length: f64,
    pub height: f64,
    pub time_gyr: f64,
    pub id: u64,
}

impl ScarpRecord {
    pub fn midpoint(&self) -> [f64; 3] {
        self.points[self.points.len() / 2]
    }

    /// Angular radius of the cap holding the scarp and its back slope.
    pub fn cap_radius(&self, body_radius: f64) -> f64 {
        (0.5 * self.length + influence_width(self.height)) / body_radius
    }

    /// Strike at the midpoint, the trace heading there, as `geom::offset` measures it:
    /// counterclockwise from east.
    pub fn midpoint_heading(&self) -> f64 {
        let m = self.points.len() / 2;
        let p = self.points[m];
        let (east, north) = local_basis(p);
        azimuth_from(
            p,
            east,
            north,
            self.points[(m + 1).min(self.points.len() - 1)],
        )
    }
}

/// Height of a scarp from its length (Enterprise Rupes: 1,000 km long, 3 km high).
pub fn scarp_height(length: f64) -> f64 {
    SCARP_HEIGHT_COEFF_KM * (length / 1e3).powf(SCARP_HEIGHT_EXPONENT) * 1e3
}

/// Horizontal shortening across a thrust of this height.
pub fn scarp_shortening(height: f64) -> f64 {
    height / SCARP_DIP_DEG.to_radians().tan()
}

/// How far from the trace a scarp changes the surface.
fn influence_width(height: f64) -> f64 {
    3.0 * SCARP_BACK_SLOPE_FACTOR * height
}

/// Uplift at signed across-strike distance `d`, positive on the hanging wall: a step up
/// the face, then a back slope decaying away from the fault. The step is a tanh whose
/// central slope is the face slope, so its half-width is h / (2 tan(face slope)). `blur`
/// broadens the face, for degradation and for the coarse grid's smoothing.
pub fn scarp_profile(d: f64, height: f64, blur: f64) -> f64 {
    let face = (0.5 * height / SCARP_FACE_SLOPE_DEG.to_radians().tan()).hypot(blur);
    let step = 0.5 * (1.0 + (d / face).tanh());
    let back = SCARP_BACK_SLOPE_FACTOR * height;
    let decay = if d > 0.0 { (-d / back).exp() } else { 1.0 };
    height * step * decay
}

/// A length from the truncated power law N(>L) ∝ L^-a between the limits.
fn sample_length(u: f64, l_min: f64, l_max: f64) -> f64 {
    let a = SCARP_LENGTH_EXPONENT;
    let lo = l_min.powf(-a);
    let hi = l_max.powf(-a);
    (lo - u.clamp(0.0, 1.0) * (lo - hi)).powf(-1.0 / a)
}

/// A curving trace of `length` from `start` along `heading`: an overall turn of `turn`
/// with lobes from a few harmonics of heading wobble at the given `phases`, so the trace
/// bends and scallops the way lobate scarps do rather than running straight.
fn make_trace(
    start: [f64; 3],
    heading: f64,
    turn: f64,
    phases: [f64; SCARP_LOBE_HARMONICS],
    length: f64,
    spacing: f64,
    body_radius: f64,
) -> Vec<[f64; 3]> {
    let n = (length / spacing).ceil().max(2.0) as usize;
    let mut points = Vec::with_capacity(n + 1);
    let mut p = start;
    for i in 0..=n {
        points.push(p);
        let (east, north) = local_basis(p);
        let s = (i as f64 + 0.5) / n as f64;
        let mut head = heading + turn * s;
        for (k, phase) in phases.iter().enumerate() {
            let k = (k + 1) as f64;
            head += SCARP_LOBE_DEG.to_radians() / k
                * (2.0 * std::f64::consts::PI * k * s + phase).sin();
        }
        p = offset(p, east, north, spacing / body_radius, head);
    }
    points
}

/// The height change one scarp makes: the profile with `blur_fine`, minus the profile
/// with `blur_coarse` when restoring detail over a coarse stamp.
pub fn scarp_compute(
    layers: &LayerSet,
    grid: &CubeSphere,
    body_radius: f64,
    rec: &ScarpRecord,
    blur_fine: f64,
    blur_coarse: Option<f64>,
    stamper: &mut Stamper,
) -> StampOutput {
    let bedrock = layers.bedrock.as_slice();
    let cell = grid.delta() * body_radius;
    let width = influence_width(rec.height) + 2.0 * cell;
    // Bands a cell wider than the spacing, so the fan on the outside of a bend leaves no
    // gaps; the nearest point still wins where bands overlap.
    let half_band = 0.5 * rec.spacing + cell;
    // Each cell takes the value from the nearest trace point along strike.
    let mut best: HashMap<usize, (f64, f64)> = HashMap::new();
    let n = rec.points.len();
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
        stamper.gather(grid, body_radius, p, width);
        for &(idx, _, _) in &stamper.footprint {
            let d = grid.center(grid.cell(idx));
            let along = dot(d, t) * body_radius;
            if along.abs() > half_band {
                continue;
            }
            let across = dot(d, normal) * body_radius;
            let mut u = scarp_profile(across, rec.height, blur_fine);
            if let Some(c) = blur_coarse {
                u -= scarp_profile(across, rec.height, c);
            }
            let u = taper * u;
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
    }
}

/// Contraction scarps forming through the coarse loop.
pub struct ContractionProcess {
    grid: CubeSphere,
    body_radius: f64,
    area: f64,
    thermal: ThermalHistory,
    age: f64,
    rng: Rng,
    stamper: Stamper,
    /// Shortening area still to be taken up, carried between steps.
    carry: f64,
    l_min: f64,
    l_max: f64,
    stamp_min: f64,
    pub records: Vec<ScarpRecord>,
}

impl ContractionProcess {
    /// None for regimes that do not express contraction as scarps.
    pub fn new(graph: &ParamGraph, grid: CubeSphere, thermal: &ThermalHistory) -> Option<Self> {
        if graph.str("interior.tectonic_regime") != "stagnant_lid" {
            return None;
        }
        let body_radius = graph.f64("body.radius");
        let cell = 2.0 * std::f64::consts::PI * body_radius / (4.0 * grid.face_size() as f64);
        // The longest fault grows with the strain taken up: deeper faulting under higher
        // stress and the linking of segments into longer systems. Mercury's half per cent
        // gives thousand-kilometre scarps; the Moon's hundredth of a per cent, kilometres.
        let strain = graph.f64("tectonics.contraction_strain");
        let l_max = (SCARP_REFERENCE_MAX_LENGTH
            * (strain / SCARP_REFERENCE_STRAIN).powf(SCARP_LENGTH_STRAIN_EXPONENT))
        .min(SCARP_MAX_LENGTH_FRACTION * body_radius);
        if l_max <= SCARP_MIN_LENGTH {
            return None;
        }
        Some(Self {
            grid,
            body_radius,
            area: 4.0 * std::f64::consts::PI * body_radius * body_radius,
            thermal: thermal.clone(),
            age: graph.f64("body.age"),
            rng: Rng::new(graph.f64("seed.tectonics") as u64 ^ 0x5CA2_9F11),
            stamper: Stamper::new(&grid),
            carry: 0.0,
            l_min: SCARP_MIN_LENGTH,
            l_max,
            stamp_min: SCARP_LOOP_MIN_CELLS * cell,
            records: Vec::new(),
        })
    }

    /// Shortening area to take up between two times.
    fn budget(&self, t0: f64, t1: f64) -> f64 {
        let strain = |t: f64| 2.0 * self.thermal.at(t).contraction / self.body_radius;
        ((strain(t1) - strain(t0)) * self.area).max(0.0)
    }

    /// Scarps formed between two times, in order of formation, not yet emplaced. The
    /// budget is spent on the whole population; only those the loop grid can carry are
    /// returned.
    pub fn sample(&mut self, t0: f64, t1: f64) -> Vec<ScarpRecord> {
        let mut out = Vec::new();
        for (bin, b0, b1) in crate::stages::time_bins(self.age, t0, t1) {
            self.sample_bin(bin, b0, b1, &mut out);
        }
        out.sort_by(|a, b| a.time_gyr.total_cmp(&b.time_gyr).then(a.id.cmp(&b.id)));
        out
    }

    /// The scarps of one time bin, drawn with the bin as the random index.
    fn sample_bin(&mut self, bin: u64, t0: f64, t1: f64, out: &mut Vec<ScarpRecord>) {
        let step_index = bin;
        let mut budget = self.budget(t0, t1) + self.carry;
        let mut j = 0u64;
        while budget > 0.0 && j < 200_000 {
            let draw = 16 + j * 12;
            let length = sample_length(self.rng.uniform(step_index, draw), self.l_min, self.l_max);
            let height = scarp_height(length);
            budget -= length * scarp_shortening(height);
            if length >= self.stamp_min {
                let spacing = 0.5 * self.stamp_min / SCARP_LOOP_MIN_CELLS;
                let start = self.rng.unit_direction(step_index, draw + 1);
                let heading = 2.0 * std::f64::consts::PI * self.rng.uniform(step_index, draw + 3);
                let turn = SCARP_TURN_DEG.to_radians()
                    * (2.0 * self.rng.uniform(step_index, draw + 4) - 1.0);
                let time_gyr = t0 + (t1 - t0) * self.rng.uniform(step_index, draw + 5);
                let mut phases = [0.0; SCARP_LOBE_HARMONICS];
                for (k, phase) in phases.iter_mut().enumerate() {
                    *phase = 2.0
                        * std::f64::consts::PI
                        * self.rng.uniform(step_index, draw + 6 + k as u64);
                }
                let rec = ScarpRecord {
                    points: make_trace(
                        start,
                        heading,
                        turn,
                        phases,
                        length,
                        spacing,
                        self.body_radius,
                    ),
                    spacing,
                    length,
                    height,
                    time_gyr,
                    id: (step_index << 32) | j,
                };
                out.push(rec);
            }
            j += 1;
        }
        self.carry = budget.min(0.0);
    }

    /// Cut one scarp into the surface as it is now, and record it.
    pub fn apply(&mut self, layers: &mut LayerSet, rec: ScarpRecord) {
        let out = scarp_compute(
            layers,
            &self.grid,
            self.body_radius,
            &rec,
            0.0,
            None,
            &mut self.stamper,
        );
        out.apply(layers);
        self.records.push(rec);
    }

    /// Sample and emplace a step's scarps.
    pub fn step(&mut self, layers: &mut LayerSet, t0: f64, t1: f64, _step_index: u64) {
        for rec in self.sample(t0, t1) {
            self.apply(layers, rec);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_profile_steps_up_the_face_and_decays_down_the_back_slope() {
        let h = 1_000.0;
        assert!(scarp_profile(-50_000.0, h, 0.0).abs() < 1.0);
        let face = scarp_profile(2_000.0, h, 0.0);
        assert!(face > 0.5 * h, "{face}");
        let far = scarp_profile(100_000.0, h, 0.0);
        assert!(far < 0.1 * h, "{far}");
        // The central slope is the face slope.
        let face = 0.5 * h / SCARP_FACE_SLOPE_DEG.to_radians().tan();
        let slope = (scarp_profile(1.0, h, 0.0) - scarp_profile(-1.0, h, 0.0)) / 2.0;
        assert!(
            (slope - SCARP_FACE_SLOPE_DEG.to_radians().tan()).abs() < 0.02,
            "{slope}"
        );
        assert!(face > 1_000.0 && face < 3_000.0, "{face}");
        // Blur broadens the face without changing its height far from it.
        let blurred = scarp_profile(0.0, h, 10_000.0);
        assert!(blurred > 0.4 * h && blurred < 0.6 * h);
        assert!((scarp_shortening(1_000.0) - 2_144.5).abs() < 5.0);
        assert!((scarp_height(1.0e6) - 3_000.0).abs() < 50.0);
    }

    /// A scarp stamped on a flat grid: nothing on the footwall, most of the height just
    /// behind the face, little far behind it, and little near the ends of the trace.
    #[test]
    fn a_stamped_scarp_steps_up_behind_its_face_and_tapers_at_its_ends() {
        let grid = CubeSphere::new(256);
        let body_radius = 300e3;
        let mut layers = LayerSet::new(grid);
        let length = 120e3;
        let height = scarp_height(length);
        let spacing = 0.5 * grid.delta() * body_radius;
        let start = normalize([1.0, 0.2, 0.1]);
        let points = make_trace(
            start,
            0.3,
            0.2,
            [0.0; SCARP_LOBE_HARMONICS],
            length,
            spacing,
            body_radius,
        );
        let rec = ScarpRecord {
            points,
            spacing,
            length,
            height,
            time_gyr: 1.0,
            id: 1,
        };
        let mut stamper = Stamper::new(&grid);
        let out = scarp_compute(&layers, &grid, body_radius, &rec, 0.0, None, &mut stamper);
        assert!(out.cells.len() > 100);
        out.apply(&mut layers);
        let face = 0.5 * height / SCARP_FACE_SLOPE_DEG.to_radians().tan();
        let back = SCARP_BACK_SLOPE_FACTOR * height;
        let at = |k: usize, across: f64| -> f64 {
            let p = rec.points[k];
            let prev = rec.points[k - 1];
            let next = rec.points[k + 1];
            let chord = [next[0] - prev[0], next[1] - prev[1], next[2] - prev[2]];
            let radial = dot(chord, p);
            let t = normalize([
                chord[0] - radial * p[0],
                chord[1] - radial * p[1],
                chord[2] - radial * p[2],
            ]);
            let n = normalize(cross(p, t));
            let s = across / body_radius;
            let q = normalize([p[0] + s * n[0], p[1] + s * n[1], p[2] + s * n[2]]);
            layers.bedrock.as_slice()[grid.index(grid.locate(q).0)] as f64
        };
        let m = rec.points.len() / 2;
        let footwall = at(m, -4.0 * face);
        let crest = at(m, 4.0 * face);
        let far = at(m, 3.0 * back);
        assert!(footwall.abs() < 0.05 * height, "{footwall} of {height}");
        assert!(
            crest > 0.5 * height && crest < 0.75 * height,
            "{crest} of {height}"
        );
        assert!(far < 0.1 * height, "{far} of {height}");
        // Near the ends the taper has brought it most of the way down.
        let end = at(2, 4.0 * face);
        assert!(end < 0.5 * crest, "{end} vs {crest}");
    }

    #[test]
    fn lengths_follow_the_power_law_within_the_limits() {
        let (lo, hi) = (10_000.0, 1_000_000.0);
        let n = 20_000;
        let mut above_100 = 0;
        for i in 0..n {
            let l = sample_length((i as f64 + 0.5) / n as f64, lo, hi);
            assert!(l >= lo - 1e-6 && l <= hi + 1e-6);
            if l >= 100_000.0 {
                above_100 += 1;
            }
        }
        // N(>100 km) / N(>10 km) for a -1.5 slope, truncated at 1,000 km.
        let a = SCARP_LENGTH_EXPONENT;
        let expected = (100_000f64.powf(-a) - hi.powf(-a)) / (lo.powf(-a) - hi.powf(-a));
        let got = above_100 as f64 / n as f64;
        assert!((got - expected).abs() < 0.01, "{got} vs {expected}");
    }
}
