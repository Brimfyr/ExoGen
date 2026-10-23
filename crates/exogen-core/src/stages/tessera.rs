// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! The tesserae's fabric from the strain their lid took (docs/parameters.md F187, F190).
//!
//! Where the mantle's flow dragged the young lid past its strength, the lid broke and went
//! with the flow (crust::drag), shortening where the flow converged and stretching where
//! it spread, each way by the flow's own curvature; and later, as the plateaus collapsed,
//! the lower crust's flow stretched it again (F188). A strong layer shortened buckles into
//! folds across the shortening, some five times as long as the layer is thick, and one
//! stretched necks into graben across the stretching, three and a half times; the lid has
//! two such layers, the brittle layer its faults cut and the elastic lid beneath, so each
//! strain makes folds or graben at two scales (the short and long folds and the ribbons of
//! Ghent and Hansen). A fold grows as the layer's length takes up the shortening, until its
//! limbs dip as steeply as the layer's friction lets them slip, or its relief weighs on the
//! layer as much as the layer's strength holds; a graben deepens with the stretching until
//! its relief weighs as much as the layer pulled apart holds. Each state of the flow made
//! its own sets, the later laid over the earlier.
//!
//! Each fold's crest and each graben is a line along the strain's strike, the lines of a
//! set as far apart as their wavelength (evenly spaced streamlines, Jobard and Lefer 1997):
//! a line runs on as long as the strain keeps its sense and no other line of its set
//! crowds it, so folds end where the shortening does or where they converge on their
//! neighbours, plunging out at their ends. A fold is a narrow crest over broad troughs, a
//! graben a flat floor between steep walls with flat blocks between; crossing sets bound
//! cells, basins and domes where the lid shortened both ways and blocks where it stretched
//! both ways, the tiles the tesserae are named for.

use std::collections::HashMap;

use rayon::prelude::*;

use crate::constants::*;
use crate::geom::{cross, dot, local_basis, normalize};
use crate::grid::CubeSphere;
use crate::grid::projection::{PaddedField, sample_bilinear};
use crate::rng::Rng;
use crate::stages::episodic::LidStrain;
use crate::stages::stress::brittle_yields;

/// A piece of a line: from `a` to `b`, its relief at each end, m, and its line.
#[derive(Clone, Copy, Debug)]
struct Segment {
    a: [f64; 3],
    b: [f64; 3],
    ra: f32,
    rb: f32,
    line: u32,
}

/// What a line is: a fold's crest or a graben, and its set's spacing, m.
#[derive(Clone, Copy, Debug)]
struct Kind {
    graben: bool,
    spacing: f64,
}

/// One family of lines: one principal strain of one state in one layer, its segments
/// by the cell of a grid about as fine as its reach.
struct Family {
    layer: usize,
    kinds: Vec<Kind>,
    segments: Vec<Segment>,
    index: HashMap<u32, Vec<u32>>,
    grid: CubeSphere,
}

/// The fabric of the tesserae: its lines, family by family.
pub struct TesseraFabric {
    body_radius: f64,
    families: Vec<Family>,
    /// The relief a layer's strength holds, m.
    hold: f64,
    /// A graben's walls, m wide at the least: the grid's cell.
    wall: f64,
    /// The highest relief any one line reaches, m, for the fabric's bound.
    pub bound: f64,
}

/// A line's reach about it, in its set's spacings: the trough between two fold crests
/// lies half a spacing out, and the fabric fades to its mean beyond.
const REACH: f64 = 0.75;

/// The strain of one state, read anywhere: its three components in the local east and
/// north basis.
struct Strain<'a> {
    coarse: &'a CubeSphere,
    parts: [PaddedField; 3],
}

impl Strain<'_> {
    /// The principal strain `which` (0 the greatest, 1 the least) at a direction and the
    /// strike across it, a tangent unit vector.
    fn principal(&self, d: [f64; 3], which: usize) -> (f64, [f64; 3]) {
        let [a, b, c] = [0, 1, 2].map(|i| sample_bilinear(self.coarse, &self.parts[i], d) as f64);
        let mean = 0.5 * (a + c);
        let radius = (0.25 * (a - c) * (a - c) + b * b).sqrt();
        let theta = 0.5 * (2.0 * b).atan2(a - c);
        let (east, north) = local_basis(d);
        let axis = |angle: f64| {
            let (s, cs) = angle.sin_cos();
            normalize([
                east[0] * cs + north[0] * s,
                east[1] * cs + north[1] * s,
                east[2] * cs + north[2] * s,
            ])
        };
        // The greatest runs along `theta` and its strike across it; the least the reverse.
        if which == 0 {
            (mean + radius, axis(theta + std::f64::consts::FRAC_PI_2))
        } else {
            (mean - radius, axis(theta))
        }
    }
}

/// A step of `h` radians from `p` along the tangent `dir`.
fn advance(p: [f64; 3], dir: [f64; 3], h: f64) -> [f64; 3] {
    normalize([p[0] + dir[0] * h, p[1] + dir[1] * h, p[2] + dir[2] * h])
}

/// The distance from `c` to the segment, m, and how far along it the nearest point is.
fn to_segment(c: [f64; 3], s: &Segment, radius: f64) -> (f64, f64) {
    let ab = [s.b[0] - s.a[0], s.b[1] - s.a[1], s.b[2] - s.a[2]];
    let ac = [c[0] - s.a[0], c[1] - s.a[1], c[2] - s.a[2]];
    let len2 = dot(ab, ab).max(1e-30);
    let t = (dot(ac, ab) / len2).clamp(0.0, 1.0);
    let q = [ac[0] - t * ab[0], ac[1] - t * ab[1], ac[2] - t * ab[2]];
    (dot(q, q).sqrt() * radius, t)
}

impl TesseraFabric {
    /// The fabric from the lid's strain sets on the coarse grid, drawn on a grid of cells
    /// `cell` m wide: the lines the grid cannot resolve are left out. `near` says which
    /// coarse cells the fabric is wanted about.
    pub fn new(
        coarse: &CubeSphere,
        sets: &[LidStrain],
        body_radius: f64,
        gravity: f64,
        cell: f64,
        seed: u64,
        near: &(dyn Fn(usize) -> bool + Sync),
    ) -> Self {
        let (yield_tension, yield_compression) = brittle_yields(gravity);
        let weight = DENSITY_CRUST_ROCK * gravity;
        let friction = ROCK_FRICTION;
        // Anderson's normal faults dip at 45 degrees plus half the friction angle.
        let dip = std::f64::consts::FRAC_PI_4 + 0.5 * friction.atan();
        let fold_cap = yield_compression / (weight * TESSERA_BUCKLING_RATIO);
        let graben_cap = yield_tension / (weight * GROOVE_NECKING_WAVELENGTH);
        // The grid draws a line's set where its spacing spans enough cells to read as
        // ridges, not a hatch of cells, fading in over a range.
        let drawn = |length: f64| {
            let (from, to) = TESSERA_MIN_CELLS;
            let x = ((length / cell - from) / (to - from)).clamp(0.0, 1.0);
            x * x * (3.0 - 2.0 * x)
        };
        let mut jobs = Vec::new();
        for (s, set) in sets.iter().enumerate() {
            for (l, &share) in [GROOVE_RIFT_LID_SHARE, 1.0].iter().enumerate() {
                let layer = share * set.elastic;
                let fold = TESSERA_BUCKLING_RATIO * layer;
                let graben = GROOVE_NECKING_WAVELENGTH * layer;
                if drawn(fold) <= 0.0 && drawn(graben) <= 0.0 {
                    continue;
                }
                for which in 0..2 {
                    jobs.push((s, l, which, fold, graben));
                }
            }
        }
        let families: Vec<Family> = jobs
            .par_iter()
            .map(|&(s, l, which, fold_length, graben_length)| {
                let set = &sets[s];
                let parts = [0, 1, 2].map(|i| {
                    let v: Vec<f32> = set.strain.iter().map(|t| t[i]).collect();
                    PaddedField::new(coarse, &v)
                });
                let strain = Strain { coarse, parts };
                let rng = Rng::new(seed ^ 0x7E55_E2A0 ^ ((s as u64) << 20) ^ ((l as u64) << 16));
                let id = which as u64;
                // Where a line may run: about the tesserae, where the state strained the lid.
                let inside = |d: [f64; 3]| {
                    let k = coarse.index(coarse.locate(d).0);
                    near(k) && set.strain[k] != [0.0; 3]
                };
                let relief_of = |value: f64, graben: bool| -> (f64, f64) {
                    if graben {
                        let h = (0.5 * value * graben_length * dip.tan()).min(graben_cap);
                        (h * drawn(graben_length), graben_length)
                    } else {
                        let h = ((fold_length / std::f64::consts::PI) * (-value).sqrt())
                            .min(friction * fold_length / std::f64::consts::TAU)
                            .min(fold_cap);
                        (h * drawn(fold_length), fold_length)
                    }
                };
                // The lines' points by the cell of a grid as fine as the closer spacing,
                // for how near a line runs to the others.
                let near_spacing = fold_length.min(graben_length);
                let face = ((std::f64::consts::FRAC_PI_2 * body_radius / near_spacing).ceil()
                    as u32)
                    .clamp(4, 4096);
                let crowd = CubeSphere::new(face);
                let mut points: HashMap<u32, Vec<[f64; 3]>> = HashMap::new();
                let crowded = |points: &HashMap<u32, Vec<[f64; 3]>>, p: [f64; 3], within: f64| {
                    let (c, _, _) = crowd.locate(p);
                    let reach = (within / near_spacing).ceil() as i32 + 1;
                    for di in -reach..=reach {
                        for dj in -reach..=reach {
                            let k = crowd.index(crowd.neighbor(c, di, dj)) as u32;
                            if let Some(list) = points.get(&k)
                                && list.iter().any(|&q| {
                                    crate::geom::angle_between(p, q) * body_radius < within
                                })
                            {
                                return true;
                            }
                        }
                    }
                    false
                };
                let mut lines: Vec<(Kind, Vec<[f64; 3]>, Vec<f32>)> = Vec::new();
                let mut seeds: Vec<[f64; 3]> = Vec::new();
                let mut next_seed = 0usize;
                let mut coarse_seed = 0usize;
                let mut draws = 0u64;
                loop {
                    // The next seed: one beside a line drawn, or a coarse cell's centre.
                    let p0 = if next_seed < seeds.len() {
                        next_seed += 1;
                        seeds[next_seed - 1]
                    } else {
                        let mut found = None;
                        while coarse_seed < coarse.cell_count() {
                            let k = coarse_seed;
                            coarse_seed += 1;
                            if near(k) && set.strain[k] != [0.0; 3] {
                                found = Some(coarse.center(coarse.cell(k)));
                                break;
                            }
                        }
                        match found {
                            Some(p) => p,
                            None => break,
                        }
                    };
                    if !inside(p0) {
                        continue;
                    }
                    let (v0, strike0) = strain.principal(p0, which);
                    if v0 == 0.0 {
                        continue;
                    }
                    let graben = v0 > 0.0;
                    let (h0, spacing) = relief_of(v0, graben);
                    if h0 <= 0.0 {
                        continue;
                    }
                    draws += 1;
                    let spread = 2f64
                        .powf(TESSERA_SPACING_SCATTER * (2.0 * rng.uniform(id, 3 * draws) - 1.0));
                    let gap = spacing * spread;
                    if crowded(&points, p0, 0.9 * gap) {
                        continue;
                    }
                    let scatter =
                        1.0 + TESSERA_RELIEF_SCATTER * (2.0 * rng.uniform(id, 3 * draws + 1) - 1.0);
                    // Traced both ways from the seed along the strike, a sixth of a spacing
                    // at a step.
                    let step = spacing / 6.0;
                    let h = step / body_radius;
                    // Its length drawn within the range folds of a set run, so that they
                    // stand en echelon and plunge out at their own ends.
                    let (short, long) = TESSERA_LINE_LENGTH;
                    let length = short + (long - short) * rng.uniform(id, 3 * draws + 2);
                    let longest = (length * spacing / step) as usize;
                    let mut halves: [Vec<([f64; 3], f32)>; 2] = [Vec::new(), Vec::new()];
                    for (half, sign) in [(0usize, 1.0f64), (1, -1.0)] {
                        let aligned = |d: [f64; 3], prev: [f64; 3]| {
                            if dot(d, prev) < 0.0 {
                                [-d[0], -d[1], -d[2]]
                            } else {
                                d
                            }
                        };
                        let mut p = p0;
                        let mut dir = [strike0[0] * sign, strike0[1] * sign, strike0[2] * sign];
                        for _ in 0..longest / 2 {
                            let (_, s1) = strain.principal(p, which);
                            let d1 = aligned(s1, dir);
                            // The strike turning sharply: the strain alike every way, no
                            // strike to follow.
                            if dot(d1, dir) < std::f64::consts::FRAC_1_SQRT_2 {
                                break;
                            }
                            let mid = advance(p, d1, 0.5 * h);
                            let (_, s2) = strain.principal(mid, which);
                            let d2 = aligned(s2, d1);
                            let q = advance(p, d2, h);
                            if !inside(q) {
                                break;
                            }
                            let (v, _) = strain.principal(q, which);
                            if v == 0.0 || (v > 0.0) != graben {
                                break;
                            }
                            if crowded(&points, q, 0.5 * gap) {
                                break;
                            }
                            let (hq, _) = relief_of(v, graben);
                            halves[half].push((q, (hq * scatter) as f32));
                            dir = normalize(cross(cross(q, d2), q));
                            p = q;
                        }
                    }
                    let [forward, backward] = halves;
                    let mut pts: Vec<[f64; 3]> = backward.iter().rev().map(|x| x.0).collect();
                    let mut rel: Vec<f32> = backward.iter().rev().map(|x| x.1).collect();
                    pts.push(p0);
                    rel.push((h0 * scatter) as f32);
                    pts.extend(forward.iter().map(|x| x.0));
                    rel.extend(forward.iter().map(|x| x.1));
                    let n = pts.len();
                    // A line shorter than its spacing is no fold.
                    if n < 6 {
                        continue;
                    }
                    // Its ends plunge out over a spacing.
                    for (i, r) in rel.iter_mut().enumerate() {
                        let x = (i.min(n - 1 - i) as f64 / 6.0).min(1.0);
                        *r *= (x * x * (3.0 - 2.0 * x)) as f32;
                    }
                    for &q in &pts {
                        let k = crowd.index(crowd.locate(q).0) as u32;
                        points.entry(k).or_default().push(q);
                    }
                    // Seeds a spacing to either side, every spacing along.
                    for i in (0..n).step_by(6) {
                        let q = pts[i];
                        let other = if i + 1 < n { pts[i + 1] } else { pts[i - 1] };
                        let t = normalize([other[0] - q[0], other[1] - q[1], other[2] - q[2]]);
                        let side = normalize(cross(q, t));
                        for sgn in [1.0, -1.0] {
                            let dir = [side[0] * sgn, side[1] * sgn, side[2] * sgn];
                            seeds.push(advance(q, dir, gap / body_radius));
                        }
                    }
                    lines.push((Kind { graben, spacing }, pts, rel));
                }
                // The segments, by the cell of a grid as coarse as the farthest a line
                // reaches, so that a point finds every line that reaches it among the cells
                // about its own.
                let widest = fold_length.max(graben_length) * REACH;
                let face = ((std::f64::consts::FRAC_PI_2 * body_radius / widest).floor() as u32)
                    .clamp(4, 4096);
                let grid = CubeSphere::new(face);
                let mut kinds = Vec::with_capacity(lines.len());
                let mut segments = Vec::new();
                let mut index: HashMap<u32, Vec<u32>> = HashMap::new();
                for (li, (kind, pts, rel)) in lines.into_iter().enumerate() {
                    kinds.push(kind);
                    for i in 0..pts.len().saturating_sub(1) {
                        let sid = segments.len() as u32;
                        segments.push(Segment {
                            a: pts[i],
                            b: pts[i + 1],
                            ra: rel[i],
                            rb: rel[i + 1],
                            line: li as u32,
                        });
                        let k = grid.index(grid.locate(pts[i]).0) as u32;
                        index.entry(k).or_default().push(sid);
                    }
                }
                Family {
                    layer: l,
                    kinds,
                    segments,
                    index,
                    grid,
                }
            })
            .filter(|f| !f.segments.is_empty())
            .collect();
        let bound = families
            .iter()
            .flat_map(|f| f.segments.iter().map(|s| s.ra.max(s.rb) as f64))
            .fold(0.0f64, f64::max);
        Self {
            body_radius,
            families,
            hold: fold_cap,
            wall: cell,
            bound,
        }
    }

    pub fn is_empty(&self) -> bool {
        self.families.is_empty()
    }

    /// The fabric's relief at a direction, m: each family's nearest line's profile, the
    /// layers' reliefs each no higher than their strength holds.
    pub fn at(&self, d: [f64; 3]) -> f64 {
        let mut layer = [0.0f64; 2];
        let r = self.body_radius;
        let fade = |x: f64| {
            let x = x.clamp(0.0, 1.0);
            1.0 - x * x * (3.0 - 2.0 * x)
        };
        for f in &self.families {
            let (c, _, _) = f.grid.locate(d);
            let mut best = (f64::INFINITY, 0.0f64, None::<Kind>);
            for di in -1..=1 {
                for dj in -1..=1 {
                    let k = f.grid.index(f.grid.neighbor(c, di, dj)) as u32;
                    let Some(list) = f.index.get(&k) else {
                        continue;
                    };
                    for &sid in list {
                        let s = &f.segments[sid as usize];
                        let (dist, t) = to_segment(d, s, r);
                        if dist < best.0 {
                            let relief = s.ra as f64 + (s.rb as f64 - s.ra as f64) * t;
                            best = (dist, relief, Some(f.kinds[s.line as usize]));
                        }
                    }
                }
            }
            let (dist, relief, Some(kind)) = best else {
                continue;
            };
            let spacing = kind.spacing;
            if dist >= REACH * spacing {
                continue;
            }
            let value = if kind.graben {
                // A flat floor across part of the spacing between walls a cell wide.
                let half = 0.5 * TESSERA_GRABEN_FLOOR * spacing;
                -relief * fade((dist - half) / self.wall + 0.5)
            } else {
                // A narrow crest over broad troughs, fading to the mean beyond them.
                let x = (dist / spacing).min(0.5);
                let profile = 4.0 / std::f64::consts::PI - 2.0 * (std::f64::consts::PI * x).sin();
                let beyond = if dist > 0.5 * spacing {
                    fade((dist - 0.5 * spacing) / ((REACH - 0.5) * spacing))
                } else {
                    1.0
                };
                relief * profile * beyond
            };
            layer[f.layer] += value;
        }
        // Each layer's relief, every set it took laid together, weighs on it no more than
        // its strength holds: past that the relief its sets would add gives way.
        layer
            .iter()
            .map(|&v| self.hold * (v / self.hold).tanh())
            .sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::offset;

    /// F187, F190: a lid shortened one way folds across that way, its crests lines running
    /// along the other, as far apart as five times its layer, no higher than its strength
    /// holds; one stretched one way has graben across it, flat blocks between.
    #[test]
    fn a_strained_lid_folds_across_its_shortening_and_necks_across_its_stretching() {
        let coarse = CubeSphere::new(64);
        let r = 6.05e6;
        let n = coarse.cell_count();
        let elastic = 8e3;
        let fold = |strain: [f32; 3]| LidStrain {
            time: 0.0,
            elastic,
            strain: vec![strain; n],
            pressed: vec![0.0; n],
            pulled: vec![0.0; n],
        };
        let at = coarse.center(coarse.cell(coarse.index(crate::grid::Cell {
            face: 0,
            i: 32,
            j: 32,
        })));
        let (east, north) = local_basis(at);
        let cell = 2e3;
        let near = |k: usize| crate::geom::angle_between(coarse.center(coarse.cell(k)), at) < 0.05;
        // Shortened east-west: the relief varies east-west and hardly north-south.
        let fabric =
            TesseraFabric::new(&coarse, &[fold([-0.2, 0.0, 0.0])], r, 8.87, cell, 7, &near);
        // The relief's total variation along a line from the centre at an azimuth.
        let walk = |fabric: &TesseraFabric, phi: f64| -> f64 {
            let mut sum = 0.0;
            let mut last = fabric.at(at);
            for i in 1..400 {
                let h = fabric.at(offset(at, east, north, i as f64 * 250.0 / r, phi));
                sum += (h - last).abs();
                last = h;
            }
            sum
        };
        let across = walk(&fabric, 0.0);
        let along = walk(&fabric, std::f64::consts::FRAC_PI_2);
        assert!(across > 2.0 * along, "across {across}, along {along}");
        let (_, yc) = brittle_yields(8.87);
        let cap = yc / (DENSITY_CRUST_ROCK * 8.87 * TESSERA_BUCKLING_RATIO);
        assert!(
            fabric.bound <= cap * (1.0 + TESSERA_RELIEF_SCATTER) + 1.0,
            "{}",
            fabric.bound
        );
        // Three sets shortened every way, laid together: each layer's relief no more than
        // its strength holds (F189).
        let sets = [
            fold([-0.3, 0.05, -0.2]),
            fold([-0.25, -0.1, -0.3]),
            fold([-0.2, 0.0, -0.2]),
        ];
        let fabric = TesseraFabric::new(&coarse, &sets, r, 8.87, cell, 7, &near);
        let highest = (0..400)
            .map(|i| {
                fabric
                    .at(offset(at, east, north, i as f64 * 250.0 / r, 0.7))
                    .abs()
            })
            .fold(0.0f64, f64::max);
        assert!(
            highest <= 2.0 * cap + 1.0,
            "{highest} against {cap} a layer"
        );
        // Stretched: graben, flat blocks standing between them at the mean.
        let fabric = TesseraFabric::new(&coarse, &[fold([0.3, 0.0, 0.0])], r, 8.87, cell, 7, &near);
        let mut hs: Vec<f64> = (0..400)
            .map(|i| fabric.at(offset(at, east, north, i as f64 * 250.0 / r, 0.0)))
            .collect();
        hs.sort_by(f64::total_cmp);
        assert!(
            hs[4] < -50.0 && hs[395].abs() < 1.0,
            "{} {}",
            hs[4],
            hs[395]
        );
    }
}
