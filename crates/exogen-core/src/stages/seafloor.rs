// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! The ground under a liquid (PLAN 5C, the sea's floor): where and how deep the body's liquid
//! has stood over the finish's ground through its history (F411), and what that liquid kept
//! from its floor (F412).

use rayon::prelude::*;

use crate::crater::production::Oblique;
use crate::crater::{Bombardment, CraterSample, Screen};
use crate::fields::Field;
use crate::grid::CubeSphere;
use crate::grid::projection::{PaddedField, sample_bilinear};
use crate::params::ParamGraph;
use crate::stages::LoopOutput;

/// The liquid's past over the finish's ground (F411): at each of the loop's keyframes the liquid
/// its history kept standing then, filled over the finish's own ground as it stood then, the
/// finish's ground before its craters moved by as much as the loop's has moved since.
pub struct LiquidPast {
    coarse: CubeSphere,
    /// Each keyframe's time, Gyr since formation; the liquid's level then over the finish's
    /// datum, m, NaN where none stood; and how far the loop's ground stood then above where it
    /// stands at the present, m.
    frames: Vec<(f64, f64, PaddedField)>,
}

impl LiquidPast {
    /// The liquid's past from the loop's keyframes up to the present's, `keyframe_index`, over
    /// the finish's `ground` as it stands at the present, from the keyframe before `from`, Gyr
    /// since formation, on; none where no liquid ever stands.
    pub fn new(
        graph: &ParamGraph,
        loop_output: &LoopOutput,
        keyframe_index: usize,
        ground: &Field<f32>,
        from: f64,
    ) -> Option<Self> {
        let history = graph.curve("volatiles.liquid_history")?;
        if !history.points.iter().any(|p| p.1 > 0.0) {
            return None;
        }
        let radius = graph.f64("body.radius");
        let keyframes = &loop_output.keyframes[..=keyframe_index];
        let coarse = loop_output.grid;
        let now = keyframes[keyframe_index].layers.height();
        let fine = ground.grid();
        let area = 4.0 * std::f64::consts::PI * radius * radius;
        let frames = keyframes
            .iter()
            .enumerate()
            .filter(|(i, _)| *i == keyframe_index || keyframes[i + 1].time_gyr >= from)
            .map(|(_, k)| {
                let shift: Vec<f32> = k
                    .layers
                    .height()
                    .as_slice()
                    .iter()
                    .zip(now.as_slice())
                    .map(|(a, b)| a - b)
                    .collect();
                let shift = PaddedField::new(&coarse, &shift);
                let gel = history.eval(k.time_gyr);
                let level = if gel > 0.0 {
                    let mut then = ground.clone();
                    then.as_mut_slice()
                        .par_iter_mut()
                        .enumerate()
                        .for_each(|(c, h)| {
                            *h += sample_bilinear(&coarse, &shift, fine.center(fine.cell(c)));
                        });
                    then.fill_level(radius, gel * area)
                } else {
                    f64::NAN
                };
                (k.time_gyr, level, shift)
            })
            .collect();
        Some(Self { coarse, frames })
    }

    /// The keyframes' times, Gyr since formation.
    pub fn times(&self) -> impl Iterator<Item = f64> + '_ {
        self.frames.iter().map(|f| f.0)
    }

    /// The depth of the liquid at keyframe `i` over a place of the finish's ground in direction
    /// `d`, standing `ground` m at the present, m.
    fn depth_at(&self, i: usize, d: [f64; 3], ground: f64) -> f64 {
        let (_, level, shift) = &self.frames[i];
        if !level.is_finite() {
            return 0.0;
        }
        (level - ground - sample_bilinear(&self.coarse, shift, d) as f64).max(0.0)
    }

    /// The depth of the liquid at `t`, Gyr since formation, over a place of the finish's ground
    /// in direction `d`, standing `ground` m at the present, m: read between the keyframes.
    pub fn depth(&self, d: [f64; 3], ground: f64, t: f64) -> f64 {
        let n = self.frames.len();
        let k = self.frames.partition_point(|f| f.0 <= t);
        if k == 0 {
            return self.depth_at(0, d, ground);
        }
        if k == n {
            return self.depth_at(n - 1, d, ground);
        }
        let (t0, t1) = (self.frames[k - 1].0, self.frames[k].0);
        let w = ((t - t0) / (t1 - t0).max(1e-12)).clamp(0.0, 1.0);
        (1.0 - w) * self.depth_at(k - 1, d, ground) + w * self.depth_at(k, d, ground)
    }
}

/// The crater the floor takes under the liquid where and when it struck (F412): the share of
/// its dry diameter its impactor keeps through the liquid's depth then; none where the liquid
/// keeps it from the floor or leaves it under `d_min`, the smallest the finish stamps.
pub fn screened(
    c: CraterSample,
    past: &LiquidPast,
    screen: &Screen,
    ground: &Field<f32>,
    d_min: f64,
    oblique: &Oblique,
) -> Option<CraterSample> {
    let grid = ground.grid();
    let h = ground.as_slice()[grid.index(grid.locate(c.center).0)] as f64;
    let depth = past.depth(c.center, h, c.time_gyr);
    if depth <= 0.0 {
        return Some(c);
    }
    let diameter = c.diameter * screen.kept(c.diameter, depth);
    (diameter >= d_min).then(|| CraterSample {
        diameter,
        ellipticity: oblique.ellipticity(diameter, c.angle),
        ..c
    })
}

/// The craters at or above the finish's smallest that the liquid kept from a place's floor
/// since a time (F412): per km², what the bombardment made there between the keyframes after
/// that time, each span's by the share of them its liquid's depth kept off, between its ends'.
/// The share kept of the craters at or above `d_min` under a depth is the share of the
/// production at or above the dry diameter that leaves `d_min` on the floor, the kept share of
/// the diameter to the power of the production's slope there.
pub struct KeptOff<'a> {
    past: &'a LiquidPast,
    screen: &'a Screen,
    bombardment: &'a Bombardment,
    d_min: f64,
    slope: f64,
    /// The keyframes' times and the present, Gyr since formation, and the craters per km² at
    /// or above `d_min` made from each to the present.
    times: Vec<f64>,
    made: Vec<f64>,
}

impl<'a> KeptOff<'a> {
    pub fn new(
        past: &'a LiquidPast,
        screen: &'a Screen,
        bombardment: &'a Bombardment,
        d_min: f64,
        present: f64,
    ) -> Self {
        let d_km = d_min / 1e3;
        let (n0, n1) = (
            bombardment.cumulative_1gyr(d_km),
            bombardment.cumulative_1gyr(1.1 * d_km),
        );
        let slope = if n0 > 0.0 && n1 > 0.0 && n0 > n1 {
            (n0 / n1).ln() / 1.1f64.ln()
        } else {
            2.0
        };
        let times: Vec<f64> = past
            .times()
            .map(|t| t.min(present))
            .chain(std::iter::once(present))
            .collect();
        let made = times
            .iter()
            .map(|&t| bombardment.density_above(d_km, t, present))
            .collect();
        Self {
            past,
            screen,
            bombardment,
            d_min,
            slope,
            times,
            made,
        }
    }

    /// Craters per km² at or above `d_min` the liquid kept from the floor in direction `d`,
    /// standing `ground` m at the present, since `since`, Gyr since formation.
    pub fn since(&self, d: [f64; 3], ground: f64, since: f64) -> f64 {
        let off = |depth: f64| 1.0 - self.screen.kept(self.d_min, depth).powf(self.slope);
        let n = self.times.len();
        // The share kept off at each keyframe; the present's is its last keyframe's.
        let offs: Vec<f64> = (0..n)
            .map(|i| off(self.past.depth(d, ground, self.times[i])))
            .collect();
        let mut kept_off = 0.0;
        for i in 0..n - 1 {
            let (a, b) = (self.times[i], self.times[i + 1]);
            if b <= since || b <= a || (offs[i] <= 0.0 && offs[i + 1] <= 0.0) {
                continue;
            }
            kept_off += if since <= a {
                0.5 * (offs[i] + offs[i + 1]) * (self.made[i] - self.made[i + 1])
            } else {
                let at = off(self.past.depth(d, ground, since));
                0.5 * (at + offs[i + 1])
                    * self.bombardment.density_above(self.d_min / 1e3, since, b)
            };
        }
        kept_off.max(0.0)
    }
}
