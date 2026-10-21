//! What the slow processes have added to the ground by a time, for the finishing pass.
//!
//! The loop lifts the rises' domes, grows the shields and carves the canyon on the coarse
//! ground as they happen, a step at a time and before that step's floods, and bends the
//! ground under its loads at the end of a step (F162); a corona's
//! profile is the finish's (F146), so it adds nothing here. A flood refilled at the finish must see the ground as it found
//! it, so the growth since then is taken off the base it is refilled against and put
//! back on top of what it writes. Each process has a growth law the loop followed, and
//! this evaluates the same law at any point and time, through an index of the coarse
//! cells each feature touches so a point only sees the features over it.

use crate::crater::Stamper;
use crate::geom::angle_between;
use crate::grid::CubeSphere;
use crate::grid::projection::{PaddedField, sample_bilinear};
use crate::stages::plume::{
    Rise, Shield, Trough, TroughFrames, dome, rise_progress, shield_height, shield_progress,
    trough_depth_at, trough_progress,
};

#[derive(Clone, Copy)]
enum Deformer {
    Rise(usize),
    Shield(usize),
    /// A canyon and the range of its frames that reach the cell.
    Canyon(usize, u32, u32),
}

/// The slow processes' growth, indexed by the coarse cells each feature touches.
pub struct Deformation<'a> {
    coarse: CubeSphere,
    body_radius: f64,
    /// Length of a loop step, Gyr.
    step_gyr: f64,
    present: f64,
    rises: &'a [Rise],
    shields: &'a [Shield],
    canyons: &'a [Trough],
    frames: Vec<TroughFrames>,
    /// The features over each coarse cell, as ranges into `list`.
    offsets: Vec<u32>,
    list: Vec<Deformer>,
    /// How far the ground sank under each load, on the coarse grid, and when (F162).
    flexure: Vec<(f64, PaddedField)>,
}

impl<'a> Deformation<'a> {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        coarse: CubeSphere,
        body_radius: f64,
        step_gyr: f64,
        present: f64,
        rises: &'a [Rise],
        shields: &'a [Shield],
        canyons: &'a [Trough],
        flexure: &[(f64, Vec<f32>)],
    ) -> Self {
        let cell = coarse.delta() * body_radius;
        let n = coarse.cell_count();
        let mut stamper = Stamper::new(&coarse);
        let mut pairs: Vec<(u32, Deformer)> = Vec::new();
        for (i, r) in rises.iter().enumerate() {
            stamper.gather(&coarse, body_radius, r.center, r.radius + cell);
            pairs.extend(
                stamper
                    .footprint
                    .iter()
                    .map(|&(idx, _, _)| (idx as u32, Deformer::Rise(i))),
            );
        }
        for (j, s) in shields.iter().enumerate() {
            stamper.gather(&coarse, body_radius, s.center, s.extent() + cell);
            pairs.extend(
                stamper
                    .footprint
                    .iter()
                    .map(|&(idx, _, _)| (idx as u32, Deformer::Shield(j))),
            );
        }
        let mut frames = Vec::with_capacity(canyons.len());
        // Per coarse cell, the first and last frame of the canyon in hand that reach it.
        let mut range: Vec<(u32, u32)> = vec![(u32::MAX, 0); n];
        let mut touched: Vec<usize> = Vec::new();
        for (k, c) in canyons.iter().enumerate() {
            let f = TroughFrames::new(c);
            let reach = f.search(c, 1.0, cell, cell) + cell;
            for (fi, p) in f.points().enumerate() {
                stamper.gather(&coarse, body_radius, p, reach);
                for &(idx, _, _) in &stamper.footprint {
                    let r = &mut range[idx];
                    if r.0 == u32::MAX {
                        touched.push(idx);
                    }
                    r.0 = r.0.min(fi as u32);
                    r.1 = r.1.max(fi as u32 + 1);
                }
            }
            for &idx in &touched {
                let (lo, hi) = range[idx];
                pairs.push((idx as u32, Deformer::Canyon(k, lo, hi)));
                range[idx] = (u32::MAX, 0);
            }
            touched.clear();
            frames.push(f);
        }
        pairs.sort_by_key(|p| p.0);
        let mut offsets = vec![0u32; n + 1];
        for &(c, _) in &pairs {
            offsets[c as usize + 1] += 1;
        }
        for i in 0..n {
            offsets[i + 1] += offsets[i];
        }
        let list = pairs.into_iter().map(|p| p.1).collect();
        Self {
            coarse,
            body_radius,
            step_gyr,
            present,
            rises,
            shields,
            canyons,
            frames,
            offsets,
            list,
            flexure: flexure
                .iter()
                .map(|(t, sunk)| (*t, PaddedField::new(&coarse, sunk)))
                .collect(),
        }
    }

    /// Whether anything deforms at all.
    pub fn is_empty(&self) -> bool {
        self.list.is_empty() && self.flexure.is_empty()
    }

    /// End of the loop step a time falls in: the step's growth was on the ground before
    /// that step's floods.
    fn step_end(&self, t: f64) -> f64 {
        if self.step_gyr <= 0.0 {
            return t;
        }
        (self.step_gyr * (t / self.step_gyr - 1e-9).ceil()).clamp(0.0, self.present)
    }

    /// Height added at a point by the continuous processes up to `t_slow` and by the
    /// discrete events up to `t_events`, m.
    fn at(&self, d: [f64; 3], t_slow: f64, _t_events: f64) -> f64 {
        // The ground's subsidence under the loads laid by then: each at the end of its
        // step, after that step's floods, so a time within a step sees the loads before it,
        // and the present all of its own.
        let mut sunk = 0.0;
        for (time, field) in &self.flexure {
            let before = if t_slow >= self.present - 1e-9 {
                *time <= self.present + 1e-9
            } else {
                *time < t_slow - 1e-9
            };
            if before {
                sunk += sample_bilinear(&self.coarse, field, d) as f64;
            }
        }
        let c = self.coarse.index(self.coarse.locate(d).0);
        let (lo, hi) = (self.offsets[c] as usize, self.offsets[c + 1] as usize);
        if lo == hi {
            return -sunk;
        }
        let r_body = self.body_radius;
        let coarse_cell = self.coarse.delta() * r_body;
        let mut h = 0.0;
        for &deformer in &self.list[lo..hi] {
            match deformer {
                Deformer::Rise(i) => {
                    let rise = &self.rises[i];
                    let x = angle_between(rise.center, d) * r_body / rise.radius;
                    if x < 1.0 {
                        h += rise.height * dome(x) * rise_progress(rise, t_slow);
                    }
                }
                Deformer::Shield(j) => {
                    let s = &self.shields[j];
                    let p = shield_progress(t_slow, s);
                    if p > 0.0 {
                        h += shield_height(s, d, r_body, false, 0.0) * p;
                    }
                }
                Deformer::Canyon(k, lo, hi) => {
                    let c = &self.canyons[k];
                    let p = trough_progress(c, t_slow);
                    if p > 0.0 {
                        h -= trough_depth_at(
                            d,
                            r_body,
                            c,
                            &self.frames[k],
                            (lo as usize, hi as usize),
                            p,
                            coarse_cell,
                            coarse_cell,
                            0.0,
                        );
                    }
                }
            }
        }
        h - sunk
    }

    /// Height added at a point since a flood at `t` found the ground: the continuous
    /// growth after the flood's step and the events after it, m.
    pub fn since_flood(&self, d: [f64; 3], t: f64) -> f64 {
        let t1 = self.step_end(t);
        self.at(d, self.present, self.present) - self.at(d, t1, t1)
    }

    /// The same for the resurfacing event, which the loop applies before the coronae of
    /// its own step.
    pub fn since_event(&self, d: [f64; 3], t: f64) -> f64 {
        let t1 = self.step_end(t);
        self.at(d, self.present, self.present) - self.at(d, t1, t)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn growth_after_a_flood_is_what_the_loop_added_later() {
        let coarse = CubeSphere::new(32);
        let body_radius = 3.39e6;
        let center = [1.0, 0.0, 0.0];
        let shields = vec![Shield {
            center,
            radius: 300e3,
            height: 20e3,
            start_gyr: 1.0,
            calderas: Vec::new(),
            rift: 0.0,
            elongation: 1.0,
            growth_gyr: 1.0,
            id: 7,
        }];
        let present = 4.5;
        let steps = 90;
        let deform = Deformation::new(
            coarse,
            body_radius,
            present / steps as f64,
            present,
            &[],
            &shields,
            &[],
            &[],
        );
        // At the summit, a flood at 2 Gyr saw the growth to the end of its step; what
        // came after is the rest of the exponential approach to the full height.
        let t = 2.0;
        let t1 = (present / steps as f64) * (t * steps as f64 / present).ceil();
        let tau = shields[0].growth_gyr;
        let expect = 20e3 * ((-(t1 - 1.0) / tau).exp() - (-(present - 1.0) / tau).exp());
        let got = deform.since_flood(center, t);
        assert!((got - expect).abs() < 1.0, "{got} vs {expect}");
        // Nothing grew since the present, and nothing grows far from the shield.
        assert!(deform.since_flood(center, present).abs() < 1e-6);
        assert_eq!(deform.since_flood([0.0, 1.0, 0.0], 0.5), 0.0);
        // A flood before the shield started sees all of it as later growth.
        let all = 20e3 * (1.0 - (-(present - 1.0) / tau).exp());
        assert!((deform.since_flood(center, 0.5) - all).abs() < 1.0);
    }
}
