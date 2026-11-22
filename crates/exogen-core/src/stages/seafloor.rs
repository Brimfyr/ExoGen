// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! The ground under a liquid (PLAN 5C, the sea's floor): where and how deep the body's liquid
//! has stood over the finish's ground through its history (F411), what that liquid kept from
//! its floor (F412), and the wedges its rivers' load builds on its floor (F414).

use std::cmp::Reverse;
use std::collections::BinaryHeap;

use rayon::prelude::*;

use crate::constants::{FORESET_SLOPE, FRONT_TAIL, WEDGE_ROUNDS};
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

/// The most a walk by the cells about each cell outruns the straight line it follows, along a
/// line at 22.5° to the grid's: cos 22.5° - sin 22.5° + √2 sin 22.5°.
const STEP_DETOUR: f64 = 1.082_392_2;

/// The cells about a cell the load spreads to.
const AROUND: [(i32, i32); 8] = [
    (1, 0),
    (-1, 0),
    (0, 1),
    (0, -1),
    (1, 1),
    (1, -1),
    (-1, 1),
    (-1, -1),
];

/// A distance and a cell, ordered by distance, then cell.
#[derive(Clone, Copy, PartialEq)]
struct Near(f64, u32);

impl Eq for Near {}

impl PartialOrd for Near {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Near {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.0.total_cmp(&other.0).then(self.1.cmp(&other.1))
    }
}

/// The depth a wedge's front lays over a floor `depth` m under the liquid's level, `x` m past
/// its topset's edge (F414), none past the front's reach: a front that advances while its load
/// spreads down it as by diffusion thins as e^(-x/L) (in its own frame, D f'' = -U f'), here
/// over the floor's depth, its length L that of FORESET_SLOPE at the topset's edge, its tail
/// cut where its gradient has fallen to e^-FRONT_TAIL of that and lowered to meet the floor
/// there, so that it falls at the measured slope at its break and flattens into the floor at
/// its foot, the continental rise's gentle gradients between.
fn front(x: f64, depth: f64) -> Option<f64> {
    if depth <= 0.0 {
        return None;
    }
    let tail = (-FRONT_TAIL).exp();
    let length = depth / (FORESET_SLOPE * (1.0 - tail));
    let u = x / length;
    (u < FRONT_TAIL).then(|| depth * ((-u).exp() - tail) / (1.0 - tail))
}

/// The reach of a wedge whose mouth stands `depth` m under the liquid's level, laying `load` m³
/// on flat floor at its mouth's ground (F414), m: a mound's, a cone of the front's slope, while
/// that stays under the level, negative, its top that far under it along the front; past that,
/// its topset's, the cone's top flattened at the level.
fn reach_of(load: f64, depth: f64) -> f64 {
    let pi = std::f64::consts::PI;
    let height = (3.0 * FORESET_SLOPE * FORESET_SLOPE * load / pi).cbrt();
    if height <= depth || depth <= 0.0 {
        (height - depth) / FORESET_SLOPE
    } else {
        let w = depth / FORESET_SLOPE;
        let c = w * w / 3.0 - load / (pi * depth);
        0.5 * (-w + (w * w - 4.0 * c).max(0.0).sqrt())
    }
}

/// A sea's search for the share its mouths' wedges are built at (F414): the share, the share
/// that laid less than its mouths brought and the one that laid at least that, each with what it
/// laid less what they brought, weighed as the false position keeps them (the Illinois rule);
/// which end was kept last; and whether it is found.
#[derive(Clone, Copy)]
struct Share {
    share: f64,
    low: (f64, f64),
    high: Option<(f64, f64)>,
    kept: i8,
    found: bool,
}

/// The rivers' load laid on the floor of the sea each enters (F414): each mouth's wedge built
/// out from it through the sea, its top at the liquid's level and its front falling at
/// FORESET_SLOPE to the floor; a load too small to reach the level a mound on the floor at its
/// mouth, one whose mouth stands over deeper floor beside it falling to the deep. A sea's
/// mouths' wedges stand together, the ground at each place taking the highest that reaches it:
/// each wedge's surface `level - FORESET_SLOPE * max(0, d - R)`, `d` its shortest way from its
/// mouth through the sea, on from each place whose ground it reaches and no further past one
/// whose ground stands over it, and `R` its reach, the topset's, or under the level a mound's
/// top's (`reach_of`), at a share of its load the same for all the sea's mouths: the share at
/// which the sea takes all they bring, from one, doubled or halved until it brackets it, then
/// between by the false position, the sea's wedges laid by what they bring's share of what they
/// laid. A sea that has no room left for its load, its wedges reaching as far as any can, is
/// filled to its level and the rest, where `hold`, held at its mouths by their shares for the
/// rivers to carry on, else spread over all it reached. `mouths` the load each cell of the sea
/// `sea` takes, m³, its floor at `ground`, m, under a liquid at `level`, m, on a grid of a body
/// of `radius` m. Returns the depth laid at each cell, m, each full sea's mouths with the rest
/// each held or spread, m³, and the cells of the seas it filled.
pub fn wedges(
    grid: CubeSphere,
    ground: &[f32],
    sea: &[bool],
    mouths: &[f64],
    (level, radius): (f64, f64),
    hold: bool,
) -> (Vec<f64>, Vec<(usize, f64)>, Vec<u32>) {
    let n = grid.cell_count();
    let clock = std::time::Instant::now();
    let sources: Vec<usize> = (0..n).filter(|&k| sea[k] && mouths[k] > 0.0).collect();
    let mut laid = vec![0.0f64; n];
    if sources.is_empty() {
        return (laid, Vec::new(), Vec::new());
    }
    let g = |k: usize| ground[k] as f64;
    let tangents = grid.center_tangents();
    let area = |k: usize| grid.cell_area(grid.cell(k)) * radius * radius;
    // Each mouth's sea: the sea's cells joined by their neighbours, from each mouth.
    let mut basin = vec![u32::MAX; n];
    let mut seas = 0u32;
    for &m in &sources {
        if basin[m] != u32::MAX {
            continue;
        }
        let mut stack = vec![m];
        basin[m] = seas;
        while let Some(k) = stack.pop() {
            let c = grid.cell(k);
            for (di, dj) in AROUND {
                let j = grid.index(grid.neighbor(c, di, dj));
                if sea[j] && basin[j] == u32::MAX {
                    basin[j] = seas;
                    stack.push(j);
                }
            }
        }
        seas += 1;
    }
    let seas = seas as usize;
    let of: Vec<usize> = sources.iter().map(|&m| basin[m] as usize).collect();
    let mut brought = vec![0.0f64; seas];
    for (i, &m) in sources.iter().enumerate() {
        brought[of[i]] += mouths[m];
    }
    // A share past which a sea's wedges reach as far as any can: the topset of its mouths'
    // widest across the body's half circumference.
    let widest = std::f64::consts::PI * radius;
    let mut shares = vec![
        Share {
            share: 1.0,
            low: (0.0, f64::NAN),
            high: None,
            kept: 0,
            found: false,
        };
        seas
    ];
    let mut key = vec![0.0f64; n];
    let mut owner = vec![0u32; n];
    // The way walked from each place's mouth and the surface the load comes to it from, each
    // mouth's centre, and its reach this round.
    let mut walked = vec![0.0f64; n];
    let mut came = vec![0.0f64; n];
    let centres: Vec<[f64; 3]> = sources
        .iter()
        .map(|&m| CubeSphere::center_from(grid.cell(m), &tangents))
        .collect();
    let mut reaches = vec![0.0f64; sources.len()];
    let (mut seen, mut done) = (vec![0u32; n], vec![0u32; n]);
    let mut heap: BinaryHeap<Reverse<Near>> = BinaryHeap::new();
    let mut reached: Vec<(u32, u32, f64)> = Vec::new();
    let mut volumes = vec![0.0f64; seas];
    let mut full = vec![false; seas];
    let (mut rounds, mut taken) = (0u32, 0u64);
    for round in 1..=WEDGE_ROUNDS {
        rounds = round;
        heap.clear();
        reached.clear();
        volumes.iter_mut().for_each(|v| *v = 0.0);
        // Each sea's widest reach of its mouths at its share, to tell when it can reach no
        // further.
        let mut most = vec![f64::NEG_INFINITY; seas];
        for (i, &m) in sources.iter().enumerate() {
            let s = of[i];
            let reach = reach_of(shares[s].share * mouths[m], level - g(m));
            most[s] = most[s].max(reach);
            reaches[i] = reach;
            let start = -reach;
            if seen[m] != round || start < key[m] {
                (seen[m], key[m], owner[m], walked[m]) = (round, start, i as u32, 0.0);
                came[m] = level;
                heap.push(Reverse(Near(start, m as u32)));
            }
        }
        while let Some(Reverse(Near(d, k))) = heap.pop() {
            let k = k as usize;
            if done[k] == round || d > key[k] {
                continue;
            }
            done[k] = round;
            taken += 1;
            // Its ground over the surface the load comes from holds it back; within its
            // front's reach, the front's depth over its floor, no higher than that surface.
            if g(k) > came[k] + 1e-6 {
                continue;
            }
            let Some(lies) = front(d.max(0.0), level - g(k)) else {
                continue;
            };
            let depth = lies.min(came[k] - g(k)).max(0.0);
            let s = g(k) + depth;
            let i = owner[k];
            reached.push((k as u32, i, depth));
            if depth > 0.0 {
                volumes[of[i as usize]] += depth * area(k);
            }
            let c = grid.cell(k);
            let here = CubeSphere::center_from(c, &tangents);
            for (di, dj) in AROUND {
                let to = grid.neighbor(c, di, dj);
                let m = grid.index(to);
                if !sea[m] || done[m] == round {
                    continue;
                }
                // The way walked from its mouth, step by step between the centres (their chord
                // their arc within a part in a million at a cell's reach), and the straight
                // way back to it: the walked way over the steps' most along a straight line
                // where that is longer, round a headland, else the straight way.
                let there = CubeSphere::center_from(to, &tangents);
                let step = ((here[0] - there[0]).powi(2)
                    + (here[1] - there[1]).powi(2)
                    + (here[2] - there[2]).powi(2))
                .sqrt();
                let way = walked[k] + radius * step;
                let from = centres[i as usize];
                let chord = ((from[0] - there[0]).powi(2)
                    + (from[1] - there[1]).powi(2)
                    + (from[2] - there[2]).powi(2))
                .sqrt();
                let straight = 2.0 * radius * (0.5 * chord).min(1.0).asin();
                let next = straight.max(way / STEP_DETOUR) - reaches[i as usize];
                if seen[m] != round || next < key[m] {
                    (seen[m], key[m], owner[m], walked[m]) = (round, next, i, way);
                    came[m] = s;
                    heap.push(Reverse(Near(next, m as u32)));
                }
            }
        }
        // Each sea's share from what its wedges laid: found where they laid at least what its
        // mouths brought, within a part in a billion of it or a share within a part in a million
        // of one that laid less; full where its widest wedge reaches as far as any can.
        let mut settled = true;
        for (s, sh) in shares.iter_mut().enumerate() {
            if sh.found {
                continue;
            }
            let (load, volume) = (brought[s], volumes[s]);
            let off = volume - load;
            if off >= 0.0 {
                if sh.kept == 1 {
                    sh.low.1 *= 0.5;
                }
                (sh.high, sh.kept) = (Some((sh.share, off)), 1);
            } else {
                if sh.kept == -1
                    && let Some(h) = sh.high.as_mut()
                {
                    h.1 *= 0.5;
                }
                (sh.low, sh.kept) = ((sh.share, off), -1);
            }
            let tight = |h: f64, lo: f64| h - lo <= 1e-6 * h;
            match sh.high {
                Some((h, _)) if off >= 0.0 && (off <= 1e-9 * load || tight(h, sh.low.0)) => {
                    sh.found = true;
                }
                None if most[s] >= widest => {
                    (sh.found, full[s]) = (true, true);
                }
                _ => {}
            }
            if sh.found {
                continue;
            }
            settled = false;
            sh.share = match sh.high {
                None => sh.share * 4.0,
                Some((h, _)) if tight(h, sh.low.0) => h,
                Some((h, above)) if sh.low.1.is_finite() => {
                    let (lo, below) = sh.low;
                    let guess = (lo * above - h * below) / (above - below);
                    if guess > lo && guess < h {
                        guess
                    } else {
                        0.5 * (lo + h)
                    }
                }
                // No share yet that laid less: a quarter of the least that laid more.
                Some((h, _)) => 0.25 * h,
            };
        }
        if settled {
            break;
        }
    }
    // Each sea's wedges laid by what its mouths brought's share of what they laid, the rest of a
    // full one's held at its mouths by their shares, or spread over all it reached.
    let scale: Vec<f64> = (0..seas)
        .map(|s| {
            if volumes[s] > brought[s] {
                brought[s] / volumes[s]
            } else {
                1.0
            }
        })
        .collect();
    let rest: Vec<f64> = (0..seas)
        .map(|s| (brought[s] - volumes[s]).max(0.0))
        .collect();
    let mut spread = vec![0.0f64; seas];
    for &(k, i, _) in &reached {
        spread[of[i as usize]] += area(k as usize);
    }
    let mut filled_cells = Vec::new();
    for &(k, i, depth) in &reached {
        let s = of[i as usize];
        let mut d = depth * scale[s];
        if rest[s] > 0.0 && !hold && spread[s] > 0.0 {
            d += rest[s] / spread[s];
        }
        laid[k as usize] = d;
        if full[s] {
            filled_cells.push(k);
        }
    }
    let mut rests = Vec::new();
    for (i, &m) in sources.iter().enumerate() {
        let s = of[i];
        if rest[s] > 0.0 {
            let part = rest[s] * mouths[m] / brought[s];
            if !hold && spread[s] <= 0.0 {
                laid[m] += part / area(m);
            }
            rests.push((m, part));
        }
    }
    if std::env::var_os("EXOGEN_DEBUG").is_some() {
        eprintln!(
            "wedges: {} mouths in {seas} seas, {rounds} rounds, {taken} cells taken, {} full; {:.2} s",
            sources.len(),
            full.iter().filter(|&&f| f).count(),
            clock.elapsed().as_secs_f64()
        );
    }
    if std::env::var_os("EXOGEN_DEBUG").is_some() {
        let given: f64 = brought.iter().sum();
        let held: f64 = rests.iter().map(|r| r.1).sum();
        let placed: f64 = (0..n)
            .filter(|&k| laid[k] > 0.0)
            .map(|k| laid[k] * area(k))
            .sum();
        let kept = placed + if hold { held } else { 0.0 };
        if (kept - given).abs() > 1e-6 * given.max(1.0) {
            eprintln!("wedges: given {given:.6e} m³, laid {placed:.6e}, held {held:.4e}");
        }
    }
    (laid, rests, filled_cells)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::dot;

    const RADIUS: f64 = 2e5;

    /// A grid whose cells' ground and sea are `f` of their centres.
    fn ground(f: impl Fn([f64; 3]) -> (f32, bool)) -> (CubeSphere, Vec<f32>, Vec<bool>) {
        let grid = CubeSphere::new(64);
        let (g, s) = (0..grid.cell_count())
            .map(|k| f(grid.center(grid.cell(k))))
            .unzip();
        (grid, g, s)
    }

    fn area(grid: CubeSphere, k: usize) -> f64 {
        grid.cell_area(grid.cell(k)) * RADIUS * RADIUS
    }

    fn apart(a: [f64; 3], b: [f64; 3]) -> f64 {
        RADIUS * dot(a, b).clamp(-1.0, 1.0).acos()
    }

    /// F414: a river's load builds a wedge into its sea from its mouth: all of its volume,
    /// none on the land; its top at the liquid's level at the mouth (within a metre, the sea's
    /// wedges laid by its load's share of what they laid) and no higher anywhere; its front
    /// falling from that top toward the floor, nowhere steeper than FORESET_SLOPE and nearly as
    /// steep as that at its topset's edge.
    #[test]
    fn a_river_s_load_builds_a_wedge_into_its_sea() {
        let (grid, floor, sea) = ground(|d| {
            if d[0] < 0.0 {
                (100.0, false)
            } else {
                (-2_000.0, true)
            }
        });
        let n = grid.cell_count();
        let mouth = grid.index(grid.locate(crate::geom::normalize([0.012, 1.0, 0.0])).0);
        assert!(sea[mouth]);
        let load = 3e13;
        let mut mouths = vec![0.0; n];
        mouths[mouth] = load;
        let (laid, rests, _) = wedges(grid, &floor, &sea, &mouths, (0.0, RADIUS), false);
        let over: f64 = rests.iter().map(|r| r.1).sum();
        let volume: f64 = (0..n).map(|k| laid[k] * area(grid, k)).sum();
        let surface = |k: usize| floor[k] as f64 + laid[k];
        let (mut steepest, mut topset, mut front) = (0.0f64, 0, 0);
        for k in (0..n).filter(|&k| laid[k] > 0.0) {
            if surface(k) > -1.0 {
                topset += 1;
            } else {
                front += 1;
            }
            let c = grid.cell(k);
            for (di, dj) in [(1, 0), (0, 1), (-1, 0), (0, -1)] {
                let m = grid.index(grid.neighbor(c, di, dj));
                if laid[m] <= 0.0 {
                    continue;
                }
                let run = apart(grid.center(c), grid.center(grid.cell(m)));
                steepest = steepest.max((surface(k) - surface(m)).abs() / (FORESET_SLOPE * run));
            }
        }
        println!(
            "laid {volume:.4e} m³ of {load:.4e}, {over:.3e} over the level; its topset {topset} cells, its front {front}; its steepest {steepest:.3} of the front's slope"
        );
        assert!((volume / load - 1.0).abs() < 1e-6, "{volume} {load}");
        assert!(over == 0.0);
        assert!((0..n).all(|k| sea[k] || laid[k] == 0.0));
        assert!((0..n).all(|k| surface(k) < 1e-6 || !sea[k]));
        assert!(surface(mouth).abs() < 1.0, "{}", surface(mouth));
        assert!(topset > 1 && front > topset, "{topset} {front}");
        assert!(steepest > 0.9 && steepest < 1.0 + 1e-6, "{steepest}");
    }

    /// F414: a wedge's front flattens into the floor at its foot: where it ends it stands over
    /// the floor by less than a twentieth of the floor's depth, where a straight front at
    /// FORESET_SLOPE ended a cell's fall of it, some 210 m, above the floor at its toe.
    #[test]
    fn a_front_flattens_into_the_floor_at_its_foot() {
        let probe = CubeSphere::new(64);
        let centre = probe.center(probe.locate([0.0, 1.0, 0.0]).0);
        let (grid, floor, sea) = ground(|_| (-2_000.0, true));
        let n = grid.cell_count();
        let mouth = grid.index(grid.locate(centre).0);
        let mut mouths = vec![0.0; n];
        mouths[mouth] = 3.45e13;
        let (laid, _, _) = wedges(grid, &floor, &sea, &mouths, (0.0, RADIUS), true);
        let mut step = 0.0f64;
        for k in (0..n).filter(|&k| laid[k] > 0.0) {
            let c = grid.cell(k);
            let toe = AROUND
                .iter()
                .any(|&(di, dj)| laid[grid.index(grid.neighbor(c, di, dj))] <= 0.0);
            if toe {
                step = step.max(laid[k]);
            }
        }
        println!("its foot stands {step:.1} m over the floor at most, the floor 2000 m deep");
        assert!(step < 0.05 * 2_000.0, "{step}");
    }

    /// F414: a wedge in an open sea over a flat floor is round: its toe as far from its mouth
    /// every way round within a cell, where a walk by the grid's steps alone fell
    /// short by up to 8% along lines at 22.5° to the grid's, its fan an octagon.
    #[test]
    fn a_wedge_on_a_flat_floor_is_round() {
        let probe = CubeSphere::new(64);
        let centre = probe.center(probe.locate([0.0, 1.0, 0.0]).0);
        let (grid, floor, sea) = ground(|_| (-2_000.0, true));
        let n = grid.cell_count();
        let mouth = grid.index(grid.locate(centre).0);
        let mut mouths = vec![0.0; n];
        mouths[mouth] = 3.45e13;
        let (laid, _, _) = wedges(grid, &floor, &sea, &mouths, (0.0, RADIUS), true);
        let (east, north) = crate::geom::local_basis(centre);
        // In narrow sectors, so that none holds both a line of the grid's and one at 22.5° to it.
        let mut far = [0.0f64; 48];
        for k in (0..n).filter(|&k| laid[k] > 1.0) {
            let d = grid.center(grid.cell(k));
            let off = [d[0] - centre[0], d[1] - centre[1], d[2] - centre[2]];
            let angle = dot(off, north).atan2(dot(off, east));
            let bin = ((angle / std::f64::consts::TAU + 1.0) * 48.0).floor() as usize % 48;
            far[bin] = far[bin].max(apart(d, centre));
        }
        let (least, most) = far
            .iter()
            .fold((f64::INFINITY, 0.0f64), |(a, b), &f| (a.min(f), b.max(f)));
        let cell = grid.delta() * RADIUS;
        println!(
            "its toe from {:.1} to {:.1} km off its mouth every way round, a cell {:.1} km",
            least / 1e3,
            most / 1e3,
            cell / 1e3
        );
        assert!(least > 10.0 * cell, "{least}");
        assert!(most - least < cell, "{least} {most}");
    }

    /// F414: a wedge's front spreads on over the ground it reaches and stops at a ridge its
    /// surface does not top: a mound in an open sea within a submerged ring its front reaches,
    /// its crest over the mound's surface there, lays nothing on the ring nor beyond it, where
    /// without the ring it lays beyond.
    #[test]
    fn a_wedge_stops_at_a_ridge_its_front_does_not_top() {
        let centre = [0.0, 1.0, 0.0];
        let ring = |d: [f64; 3]| (15e3..22.5e3).contains(&apart(d, centre));
        let (grid, floor, sea) = ground(|d| (if ring(d) { -1_200.0 } else { -2_000.0 }, true));
        let (_, open, _) = ground(|_| (-2_000.0, true));
        let n = grid.cell_count();
        let mouth = grid.index(grid.locate(centre).0);
        let load = 3e11;
        let mut mouths = vec![0.0; n];
        mouths[mouth] = load;
        let beyond = |laid: &[f64]| -> f64 {
            (0..n)
                .filter(|&k| apart(grid.center(grid.cell(k)), centre) >= 22.5e3)
                .map(|k| laid[k] * area(grid, k))
                .sum()
        };
        let (walled, _, _) = wedges(grid, &floor, &sea, &mouths, (0.0, RADIUS), false);
        let (free, _, _) = wedges(grid, &open, &sea, &mouths, (0.0, RADIUS), false);
        let on_ring: f64 = (0..n)
            .filter(|&k| ring(grid.center(grid.cell(k))))
            .map(|k| walled[k])
            .sum();
        let volume: f64 = (0..n).map(|k| walled[k] * area(grid, k)).sum();
        println!(
            "within the ring {volume:.4e} m³, on it {on_ring:.3} m, beyond it {:.3e} m³; without it {:.3e} beyond",
            beyond(&walled),
            beyond(&free)
        );
        assert!((volume / load - 1.0).abs() < 1e-6);
        assert!(beyond(&free) > 1e-3 * load);
        assert!(on_ring == 0.0 && beyond(&walled) == 0.0);
    }

    /// F414: where the sea a load reaches has no room for it, the sea is filled to its level
    /// and the rest lies over what it filled, none of it on the land about it.
    #[test]
    fn a_closed_sea_fills_and_holds_the_rest_over_it() {
        let centre = [0.0, 1.0, 0.0];
        let (grid, floor, sea) = ground(|d| {
            if apart(d, centre) < 12e3 {
                (-300.0, true)
            } else {
                (100.0, false)
            }
        });
        let n = grid.cell_count();
        let room: f64 = (0..n)
            .filter(|&k| sea[k])
            .map(|k| -(floor[k] as f64) * area(grid, k))
            .sum();
        let mouth = grid.index(grid.locate(centre).0);
        let mut mouths = vec![0.0; n];
        mouths[mouth] = 2.0 * room;
        let (laid, rests, _) = wedges(grid, &floor, &sea, &mouths, (0.0, RADIUS), false);
        let over: f64 = rests.iter().map(|r| r.1).sum();
        let volume: f64 = (0..n).map(|k| laid[k] * area(grid, k)).sum();
        println!("the sea's room {room:.4e} m³; laid {volume:.4e}, {over:.4e} over the level");
        assert!((volume / (2.0 * room) - 1.0).abs() < 1e-9);
        assert!((over / room - 1.0).abs() < 1e-9, "{over} {room}");
        assert!((0..n).all(|k| !sea[k] || floor[k] as f64 + laid[k] > 0.0));
        assert!((0..n).all(|k| sea[k] || laid[k] == 0.0));
        // Held, the rest waits at the mouth for the rivers, the sea filled to its level.
        let (held, rests, _) = wedges(grid, &floor, &sea, &mouths, (0.0, RADIUS), true);
        let volume: f64 = (0..n).map(|k| held[k] * area(grid, k)).sum();
        assert!(rests.len() == 1 && rests[0].0 == mouth, "{rests:?}");
        assert!((rests[0].1 / room - 1.0).abs() < 1e-9 && (volume / room - 1.0).abs() < 1e-9);
        assert!((0..n).all(|k| !sea[k] || (floor[k] as f64 + held[k]).abs() < 1e-6));
    }

    /// F414: a load whose mouth stands over deeper floor beside it, less than a wedge at the
    /// mouth's own ground would lay, falls to the deep beside it: all of it laid, none made,
    /// nothing over the mouth's ground.
    #[test]
    fn a_load_on_a_high_mouth_falls_to_the_deep_beside_it() {
        let probe = CubeSphere::new(64);
        let centre = probe.center(probe.locate([0.0, 1.0, 0.0]).0);
        let (grid, floor, sea) = ground(|d| {
            if apart(d, centre) < 1e3 {
                (-50.0, true)
            } else {
                (-3_000.0, true)
            }
        });
        let n = grid.cell_count();
        let mouth = grid.index(grid.locate(centre).0);
        assert!(floor[mouth] == -50.0);
        assert!((0..n).filter(|&k| floor[k] == -50.0).count() == 1);
        let mut mouths = vec![0.0; n];
        mouths[mouth] = 1e9;
        let (laid, rests, _) = wedges(grid, &floor, &sea, &mouths, (0.0, RADIUS), true);
        let volume: f64 = (0..n).map(|k| laid[k] * area(grid, k)).sum();
        println!(
            "laid {volume:.6e} m³ of 1e9, at the mouth {:.3} m, {} held",
            laid[mouth],
            rests.len()
        );
        assert!((volume / 1e9 - 1.0).abs() < 1e-9 && rests.is_empty());
        assert!(laid[mouth] < 1e-2, "{}", laid[mouth]);
    }

    /// F414: over a rough floor of many seas, small and large, every load its mouths take is
    /// laid, none made and none lost.
    #[test]
    fn many_mouths_lay_all_they_take() {
        let noise = crate::noise::GradientNoise::new(11);
        let (grid, floor, sea) = ground(|d| {
            let h = 3_000.0 * noise.fbm(d, 8.0, 4, 2.0, 0.5) - 300.0;
            (h as f32, h < 0.0)
        });
        let n = grid.cell_count();
        let mut mouths = vec![0.0; n];
        let mut taken = 0.0;
        for k in (0..n).filter(|&k| sea[k] && k % 7 == 0) {
            let c = grid.cell(k);
            let beside = AROUND
                .iter()
                .any(|&(di, dj)| !sea[grid.index(grid.neighbor(c, di, dj))]);
            if beside {
                mouths[k] =
                    1e9 * (1.0 + (k % 13) as f64) * (1.0 + (k % 5 == 0) as u8 as f64 * 300.0);
                taken += mouths[k];
            }
        }
        let (laid, rests, _) = wedges(grid, &floor, &sea, &mouths, (0.0, RADIUS), false);
        let over: f64 = rests.iter().map(|r| r.1).sum();
        let volume: f64 = (0..n).map(|k| laid[k] * area(grid, k)).sum();
        let negative = laid.iter().fold(0.0f64, |a, &b| a.min(b));
        println!(
            "taken {taken:.6e} m³, laid {volume:.6e}, {over:.4e} over the level; the least laid {negative}"
        );
        assert!(negative >= 0.0);
        assert!((volume / taken - 1.0).abs() < 1e-9, "{volume} {taken}");
    }
}
