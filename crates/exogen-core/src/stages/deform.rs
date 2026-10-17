//! What the slow processes have added to the ground by a time, for the finishing pass.
//!
//! The loop lifts the rises' domes, grows the shields, carves the canyon and raises the
//! coronae's cores on the coarse ground as they happen, a step at a time and before
//! that step's floods. A flood refilled at the finish must see the ground as it found
//! it, so the growth since then is taken off the base it is refilled against and put
//! back on top of what it writes. Each process has a growth law the loop followed, and
//! this evaluates the same law at any point and time, through an index of the coarse
//! cells each feature touches so a point only sees the features over it.

use crate::constants::*;
use crate::crater::Stamper;
use crate::geom::angle_between;
use crate::grid::CubeSphere;
use crate::stages::episodic::{Corona, corona_core};
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
    Corona(usize),
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
    coronae: &'a [Corona],
    /// The features over each coarse cell, as ranges into `list`.
    offsets: Vec<u32>,
    list: Vec<Deformer>,
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
        coronae: &'a [Corona],
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
        for (k, c) in coronae.iter().enumerate() {
            stamper.gather(
                &coarse,
                body_radius,
                c.center,
                CORONA_RIM_POSITION * c.radius + cell,
            );
            pairs.extend(
                stamper
                    .footprint
                    .iter()
                    .map(|&(idx, _, _)| (idx as u32, Deformer::Corona(k))),
            );
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
            coronae,
            offsets,
            list,
        }
    }

    /// Whether anything deforms at all.
    pub fn is_empty(&self) -> bool {
        self.list.is_empty()
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
    fn at(&self, d: [f64; 3], t_slow: f64, t_events: f64) -> f64 {
        let c = self.coarse.index(self.coarse.locate(d).0);
        let (lo, hi) = (self.offsets[c] as usize, self.offsets[c + 1] as usize);
        if lo == hi {
            return 0.0;
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
                    let p = shield_progress(t_slow, s.start_gyr);
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
                Deformer::Corona(k) => {
                    let c = &self.coronae[k];
                    if c.time_gyr <= t_events {
                        let x = angle_between(c.center, d) * r_body / c.radius;
                        h += corona_core(x, c.centre_relief);
                    }
                }
            }
        }
        h
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
    use crate::constants::SHIELD_GROWTH_GYR;

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
        let expect = 20e3
            * ((-(t1 - 1.0) / SHIELD_GROWTH_GYR).exp()
                - (-(present - 1.0) / SHIELD_GROWTH_GYR).exp());
        let got = deform.since_flood(center, t);
        assert!((got - expect).abs() < 1.0, "{got} vs {expect}");
        // Nothing grew since the present, and nothing grows far from the shield.
        assert!(deform.since_flood(center, present).abs() < 1e-6);
        assert_eq!(deform.since_flood([0.0, 1.0, 0.0], 0.5), 0.0);
        // A flood before the shield started sees all of it as later growth.
        let all = 20e3 * (1.0 - (-(present - 1.0) / SHIELD_GROWTH_GYR).exp());
        assert!((deform.since_flood(center, 0.5) - all).abs() < 1.0);
    }
}
