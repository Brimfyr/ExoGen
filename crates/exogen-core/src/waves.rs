// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! The waves the wind raises on the liquid (F428).
//!
//! At each step of the cycle's year each place's wind blows over its liquid: its mean, and the
//! weather's gusts about it at the body's mean wind (M76), taken as a normal spread whose mean
//! square is that wind's, sampled at the Gauss-Hermite nodes. A wind of speed U10 at 10 m
//! raises waves over the open liquid upwind of each place, its fetch X, by the fetch-limited
//! growth of the Coastal Engineering Manual (EM 1110-2-1100, Part II, Chapter 2, eq. II-2-36:
//! g H_m0 / u*² = 4.13e-2 (g X / u*²)^(1/2), g T_p / u* = 0.651 (g X / u*²)^(1/3)) up to its
//! fully developed sea (eq. II-2-37: 2.115e2 and 2.398e2), the friction velocity u*² = C_D
//! U10², C_D = 0.001 (1.1 + 0.035 U10), and over no more than the fetch the wind's time
//! blowing from one way gives (eq. II-2-38): the weather's gusts blow from one way for the time
//! the time their eddies, half the Rhines wavelength of their own wind across as the wind's
//! belts are (F267), take to turn their wind a quarter of the compass as they pass, a day or so
//! in the Earth's middle latitudes, the longer the slower the body spins; a storm's sea does
//! not grow to the fully developed sea of its strongest wind. The fetch is the straight way
//! over the liquid, as the manual recommends, walked over a grid a quarter of the finish's
//! across, its cells the liquid where it covers more than half of them; the waves grow over
//! the open water among its ice, the fetch by the share of the place's sea open of its ice,
//! that of its ice's mass (SEA_ICE_COVER_MASS), as the moisture reads it, so that a sea its ice
//! covers through a season raises next to none in it. The law is the manual's, fitted on the
//! Earth's air over its water, for every body
//! (a stated choice: the share of the air's momentum a liquid takes up grows with the air's
//! density over the liquid's, which differs on Titan).
//!
//! Each place keeps: the year's mean significant height; the height exceeded 12 hours a year
//! and its period (Hallermeier 1981's H_e and T_e, as Udo, Ranasinghe and Takeda 2020 give
//! them), the same share of the time on any body; and for the waves' carrying along the coast
//! (Ashton and Murray 2006; Nienhuis, Ashton and Giosan 2015), each direction's share of their
//! energy, H^(12/5) T^(1/5), by the direction they come from, and their net height and period,
//! (mean H^(12/5))^(5/12) and (mean T^(1/5))^5.

use crate::constants::*;
use crate::geom::dot;
use crate::grid::{CubeSphere, lon_lat_from_direction};
use rayon::prelude::*;
use std::f64::consts::{PI, TAU};

/// The Gauss-Hermite nodes and weights of a normal spread of unit variance (probabilists'):
/// the physicists' nine-point rule's nodes times sqrt(2), its weights over sqrt(pi).
const NODES: [(f64, f64); 9] = [
    (-4.512_745_863_399_783, 2.234_584_400_774_658e-5),
    (-3.205_429_002_856_47, 2.789_141_321_231_769e-3),
    (-2.076_847_978_677_83, 4.991_640_676_521_788e-2),
    (-1.023_255_663_789_133, 0.244_097_502_894_939),
    (0.0, 0.406_349_206_349_206),
    (1.023_255_663_789_133, 0.244_097_502_894_939),
    (2.076_847_978_677_83, 4.991_640_676_521_788e-2),
    (3.205_429_002_856_47, 2.789_141_321_231_769e-3),
    (4.512_745_863_399_783, 2.234_584_400_774_658e-5),
];

/// The waves through the year on a grid a quarter of the finish's across (F428).
#[derive(Clone, Debug)]
pub struct Waves {
    pub grid: CubeSphere,
    /// Each cell the liquid, more than half of it.
    pub sea: Vec<bool>,
    /// At each cell: the year's mean significant height, m; the height exceeded 12 hours a
    /// year and its period, m and s; the net height and period of the waves' carrying, m and
    /// s; and the share of their energy coming from each of WAVE_DIRECTIONS directions,
    /// clockwise from north.
    pub mean: Vec<f32>,
    pub exceeded: Vec<f32>,
    pub exceeded_period: Vec<f32>,
    pub net_height: Vec<f32>,
    pub net_period: Vec<f32>,
    pub energy: Vec<[f32; WAVE_DIRECTIONS]>,
    /// The fetch at each cell from each direction, m.
    pub fetch: Vec<[f32; WAVE_DIRECTIONS]>,
}

/// The waves' year in brief (F428): the mean significant height over the liquid, m; and at the
/// coasts, the liquid's cells beside the land, the median and the most of the height exceeded
/// 12 hours a year, m.
#[derive(Clone, Copy, Debug, Default)]
pub struct Summary {
    pub mean: f64,
    pub coast_exceeded_median: f64,
    pub coast_exceeded_most: f64,
}

/// The significant height and the peak period of the waves a wind of `u10` m/s at 10 m raises
/// over a fetch of `fetch` m under gravity `g` (F428): the manual's fetch-limited growth to
/// its fully developed sea.
pub fn grown(u10: f64, fetch: f64, g: f64) -> (f64, f64) {
    if u10 <= 0.0 || fetch <= 0.0 {
        return (0.0, 0.0);
    }
    let drag = WAVE_DRAG.0 * (WAVE_DRAG.1 + WAVE_DRAG.2 * u10);
    let ustar2 = drag * u10 * u10;
    let ustar = ustar2.sqrt();
    let x = g * fetch / ustar2;
    let height = (WAVE_GROWTH.0 * x.sqrt()).min(WAVE_DEVELOPED.0) * ustar2 / g;
    let period = (WAVE_GROWTH.1 * x.cbrt()).min(WAVE_DEVELOPED.1) * ustar / g;
    (height, period)
}

/// The fetch a wind of `u10` m/s blowing for `duration` s raises its waves over, m (F428): the
/// manual's duration-limited growth by its equivalent fetch, g X / u*² = WAVE_DURATION (g t /
/// u*)^(3/2) (eq. II-2-38).
pub fn duration_fetch(u10: f64, duration: f64, g: f64) -> f64 {
    if !duration.is_finite() {
        return f64::INFINITY;
    }
    let drag = WAVE_DRAG.0 * (WAVE_DRAG.1 + WAVE_DRAG.2 * u10);
    let ustar2 = drag * u10 * u10;
    if ustar2 <= 0.0 {
        return f64::INFINITY;
    }
    WAVE_DURATION * (g * duration / ustar2.sqrt()).powf(1.5) * ustar2 / g
}

/// How long a gust of the weather blows from one way at latitude `lat` on a body spinning at
/// `omega` rad/s of radius `radius` m, its gusts' speed `gust` m/s, s (F428): its eddy, half
/// the Rhines wavelength of its own wind across, pi (2 U / beta)^(1/2), as the wind's belts are
/// (F267), passes in the time it takes to cross it, its wind turning through half the compass
/// as it does, so that it blows within a quarter turn of one way for that time over
/// WAVE_EDDY_TURNS; unbounded where the spin's beta is none.
pub fn gust_duration(lat: f64, omega: f64, radius: f64, gust: f64) -> f64 {
    let speed = gust.max(WIND_LEAST_SPEED);
    let beta = 2.0 * omega * lat.cos().abs() / radius.max(1.0);
    if beta <= 0.0 {
        return f64::INFINITY;
    }
    std::f64::consts::PI * (2.0 * speed / beta).sqrt() / speed / WAVE_EDDY_TURNS
}

/// The fetch at which a wind of `u10` m/s raises its fully developed sea, m (F428).
pub fn developed_fetch(u10: f64, g: f64) -> f64 {
    let drag = WAVE_DRAG.0 * (WAVE_DRAG.1 + WAVE_DRAG.2 * u10);
    let ustar2 = drag * u10 * u10;
    (WAVE_DEVELOPED.0 / WAVE_GROWTH.0).powi(2) * ustar2 / g
}

/// The direction from which a wind blowing `east` and `north` comes, rad clockwise from north.
fn from_azimuth(east: f64, north: f64) -> f64 {
    (-east).atan2(-north).rem_euclid(TAU)
}

/// The east and north unit vectors at a direction `p`.
fn local_frame(p: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    let mut e = [-p[1], p[0], 0.0];
    let len = (e[0] * e[0] + e[1] * e[1]).sqrt();
    if len < 1e-12 {
        e = [1.0, 0.0, 0.0];
    } else {
        e = [e[0] / len, e[1] / len, 0.0];
    }
    let n = [
        p[1] * e[2] - p[2] * e[1],
        p[2] * e[0] - p[0] * e[2],
        p[0] * e[1] - p[1] * e[0],
    ];
    (e, n)
}

/// The fetch at each cell of `grid` where `sea`, from each direction, m (F428): the straight
/// way over the liquid toward it, along the great circle, in steps of a cell's width, to the
/// first cell not the liquid's, at most its `reach` m; none on the land.
pub fn fetches(
    grid: CubeSphere,
    sea: &[bool],
    radius: f64,
    reach: &[f64],
) -> Vec<[f32; WAVE_DIRECTIONS]> {
    let step = grid.delta();
    (0..grid.cell_count())
        .into_par_iter()
        .map(|k| {
            let mut out = [0.0f32; WAVE_DIRECTIONS];
            if !sea[k] {
                return out;
            }
            let reach = reach[k];
            let most = ((reach / radius) / step).ceil().max(1.0) as usize;
            let p = grid.center(grid.cell(k));
            let (e, n) = local_frame(p);
            for (b, x) in out.iter_mut().enumerate() {
                let az = TAU * b as f64 / WAVE_DIRECTIONS as f64;
                let t = [
                    az.cos() * n[0] + az.sin() * e[0],
                    az.cos() * n[1] + az.sin() * e[1],
                    az.cos() * n[2] + az.sin() * e[2],
                ];
                let mut walked = most;
                for s in 1..=most {
                    let a = s as f64 * step;
                    let q = [
                        p[0] * a.cos() + t[0] * a.sin(),
                        p[1] * a.cos() + t[1] * a.sin(),
                        p[2] * a.cos() + t[2] * a.sin(),
                    ];
                    if !sea[grid.index(grid.locate(q).0)] {
                        walked = s;
                        break;
                    }
                }
                *x = (((walked as f64 - 0.5) * step * radius).min(reach)) as f32;
            }
            out
        })
        .collect()
}

/// The fetch from the direction `azimuth`, rad clockwise from north, read between the
/// directions' fetches.
fn fetch_from(fetch: &[f32; WAVE_DIRECTIONS], azimuth: f64) -> f64 {
    let x = azimuth.rem_euclid(TAU) / TAU * WAVE_DIRECTIONS as f64;
    let i = (x.floor() as usize) % WAVE_DIRECTIONS;
    let f = x - x.floor();
    fetch[i] as f64 * (1.0 - f) + fetch[(i + 1) % WAVE_DIRECTIONS] as f64 * f
}

/// The year's wind over a place: its steps' mean winds, east and north, m/s, the weather's
/// gusts' speed, m/s, how long a gust blows from one way, s, and the share of its liquid open
/// of ice, at each step.
pub struct Year {
    pub east: Vec<f64>,
    pub north: Vec<f64>,
    pub gust: Vec<f64>,
    pub duration: Vec<f64>,
    pub open: Vec<f64>,
}

/// The height a spread of `samples`, each a height, its period and its share of the time,
/// exceeds for the `share` of the time, and its period: the samples gathered into
/// WAVE_HEIGHT_BINS bins of height, each at its samples' mean height and period, read between
/// the bins, each at the middle of its share of the time from the highest down.
fn exceeded_by(samples: &[(f64, f64, f64)], share: f64) -> (f64, f64) {
    let top = samples.iter().fold(0.0f64, |a, s| a.max(s.0));
    if top <= 0.0 {
        return (0.0, 0.0);
    }
    let mut bins = vec![(0.0f64, 0.0f64, 0.0f64); WAVE_HEIGHT_BINS];
    for &(h, t, w) in samples {
        let i = ((h / top * WAVE_HEIGHT_BINS as f64) as usize).min(WAVE_HEIGHT_BINS - 1);
        bins[i].0 += w * h;
        bins[i].1 += w * t;
        bins[i].2 += w;
    }
    let mut above = 0.0;
    let mut last: Option<(f64, f64, f64)> = None;
    for &(wh, wt, w) in bins.iter().rev() {
        if w <= 0.0 {
            continue;
        }
        let (h, t) = (wh / w, wt / w);
        let middle = above + 0.5 * w;
        if middle >= share {
            return match last {
                None => (h, t),
                Some((lh, lt, lm)) => {
                    let f = ((share - lm) / (middle - lm).max(1e-300)).clamp(0.0, 1.0);
                    (lh + (h - lh) * f, lt + (t - lt) * f)
                }
            };
        }
        last = Some((h, t, middle));
        above += w;
    }
    last.map_or((0.0, 0.0), |(h, t, _)| (h, t))
}

/// The waves at a place through its `year` over its `fetch` under gravity `g` (F428): the
/// year's mean significant height, the height exceeded the share WAVE_EXCEEDED_HOURS of a year
/// of hours and its period, the net height and period of the waves' carrying, and their
/// energy's share by the direction they come from.
pub fn at_place(
    year: &Year,
    fetch: &[f32; WAVE_DIRECTIONS],
    g: f64,
) -> (f64, f64, f64, f64, f64, [f32; WAVE_DIRECTIONS]) {
    let steps = year.east.len().max(1);
    let mut samples: Vec<(f64, f64, f64)> = Vec::with_capacity(steps * NODES.len().pow(2) + 1);
    let (mut mean, mut h24, mut t02) = (0.0, 0.0, 0.0);
    let mut energy = [0.0f64; WAVE_DIRECTIONS];
    for k in 0..steps {
        let open = year.open[k].clamp(0.0, 1.0);
        for &(gx, wx) in &NODES {
            for &(gy, wy) in &NODES {
                let w = wx * wy / steps as f64;
                let (u, v) = (
                    year.east[k] + year.gust[k] * gx / std::f64::consts::SQRT_2,
                    year.north[k] + year.gust[k] * gy / std::f64::consts::SQRT_2,
                );
                let speed = u.hypot(v);
                let az = from_azimuth(u, v);
                let reach =
                    (open * fetch_from(fetch, az)).min(duration_fetch(speed, year.duration[k], g));
                let (h, t) = grown(speed, reach, g);
                samples.push((h, t, w));
                mean += w * h;
                let (hp, tp) = if h > 0.0 {
                    (h.powf(2.4), t.powf(0.2))
                } else {
                    (0.0, 0.0)
                };
                let e = hp * tp;
                h24 += w * hp;
                t02 += w * tp;
                let b = ((az / TAU * WAVE_DIRECTIONS as f64).round() as usize) % WAVE_DIRECTIONS;
                energy[b] += w * e;
            }
        }
    }
    let (exceeded, period) = exceeded_by(&samples, WAVE_EXCEEDED_HOURS / (365.25 * 24.0));
    let total: f64 = energy.iter().sum();
    let mut shares = [0.0f32; WAVE_DIRECTIONS];
    if total > 0.0 {
        for (s, e) in shares.iter_mut().zip(energy) {
            *s = (e / total) as f32;
        }
    }
    (
        mean,
        exceeded,
        period,
        h24.powf(1.0 / 2.4),
        t02.powi(5),
        shares,
    )
}

/// The waves on the liquid of a finish's `fine` grid where `sea`, a body of `radius` m under
/// gravity `g` spinning at `omega` rad/s, through the cycle's year of `steps` steps at its
/// places of layout `(n_lat,
/// n_lon)`: their mean winds, east and north, step after step, place after place within each,
/// the weather's gusts at each step, and the mass of each place's sea's ice at each step,
/// place after place, step after step, kg/m² (F428). None without wind.
pub fn raise(
    fine: CubeSphere,
    sea: &[bool],
    (radius, g, omega): (f64, f64, f64),
    (n_lat, n_lon, steps): (usize, usize, usize),
    (east, north, gust): (&[f32], &[f32], &[f64]),
    ice: &[f32],
) -> Option<Waves> {
    let places = n_lat * n_lon;
    if steps == 0 || east.len() != steps * places || gust.len() != steps || !sea.iter().any(|&s| s)
    {
        return None;
    }
    let grid = CubeSphere::new((fine.face_size() / WAVE_FACE_SHARE).max(WAVE_LEAST_FACE));
    // The liquid's share of each cell of the coarser grid.
    let mut share = vec![(0.0f64, 0.0f64); grid.cell_count()];
    for (k, &wet) in sea.iter().enumerate() {
        let c = grid.index(grid.locate(fine.center(fine.cell(k))).0);
        share[c].1 += 1.0;
        if wet {
            share[c].0 += 1.0;
        }
    }
    let coarse_sea: Vec<bool> = share.iter().map(|&(s, n)| n > 0.0 && s / n > 0.5).collect();
    // The fetch past which the strongest of the year's winds raises nothing more.
    let strongest = (0..steps)
        .flat_map(|k| (0..places).map(move |p| (k, p)))
        .map(|(k, p)| {
            (east[k * places + p] as f64).hypot(north[k * places + p] as f64)
                + gust[k] * NODES[NODES.len() - 1].0
        })
        .fold(0.0f64, f64::max);
    // No wave reads past the fetch its strongest wind's sea develops over, nor past what that
    // wind blowing for its longest time from one way at the place gives.
    let developed = developed_fetch(strongest, g).min(PI * radius);
    let calmest = gust.iter().fold(f64::INFINITY, |a, &u| a.min(u));
    let reach: Vec<f64> = (0..grid.cell_count())
        .map(|c| {
            let (_, lat) = lon_lat_from_direction(grid.center(grid.cell(c)));
            let longest = gust_duration(lat, omega, radius, calmest);
            developed.min(duration_fetch(strongest, longest, g))
        })
        .collect();
    let fetch = fetches(grid, &coarse_sea, radius, &reach);
    let has_ice = ice.len() == steps * places;
    let read = |values: &[f32], k: usize, at: &[(usize, f64); 4]| -> f64 {
        at.iter()
            .map(|&(p, w)| w * values[k * places + p] as f64)
            .sum()
    };
    let rows: Vec<_> = (0..grid.cell_count())
        .into_par_iter()
        .map(|c| {
            if !coarse_sea[c] {
                return (0.0, 0.0, 0.0, 0.0, 0.0, [0.0f32; WAVE_DIRECTIONS]);
            }
            let (lon, lat) = lon_lat_from_direction(grid.center(grid.cell(c)));
            let at = crate::climate::reader(n_lat, n_lon, lat, lon);
            let mut year = Year {
                east: Vec::with_capacity(steps),
                north: Vec::with_capacity(steps),
                gust: gust.to_vec(),
                duration: gust
                    .iter()
                    .map(|&u| gust_duration(lat, omega, radius, u))
                    .collect(),
                open: Vec::with_capacity(steps),
            };
            for k in 0..steps {
                year.east.push(read(east, k, &at));
                year.north.push(read(north, k, &at));
                let mass = if has_ice {
                    at.iter()
                        .map(|&(p, w)| w * ice[p * steps + k] as f64)
                        .sum::<f64>()
                } else {
                    0.0
                };
                year.open.push((-mass / SEA_ICE_COVER_MASS).exp());
            }
            at_place(&year, &fetch[c], g)
        })
        .collect();
    let mut w = Waves {
        grid,
        sea: coarse_sea,
        mean: Vec::with_capacity(rows.len()),
        exceeded: Vec::with_capacity(rows.len()),
        exceeded_period: Vec::with_capacity(rows.len()),
        net_height: Vec::with_capacity(rows.len()),
        net_period: Vec::with_capacity(rows.len()),
        energy: Vec::with_capacity(rows.len()),
        fetch,
    };
    for (m, e, ep, nh, np, en) in rows {
        w.mean.push(m as f32);
        w.exceeded.push(e as f32);
        w.exceeded_period.push(ep as f32);
        w.net_height.push(nh as f32);
        w.net_period.push(np as f32);
        w.energy.push(en);
    }
    Some(w)
}

impl Waves {
    /// The year in brief (F428).
    pub fn summary(&self) -> Summary {
        let g = self.grid;
        let (mut area, mut mean) = (0.0, 0.0);
        let mut coast = Vec::new();
        for k in (0..g.cell_count()).filter(|&k| self.sea[k]) {
            let a = g.cell_area(g.cell(k));
            area += a;
            mean += a * self.mean[k] as f64;
            let c = g.cell(k);
            let beside_land = (-1..=1)
                .flat_map(|di| (-1..=1).map(move |dj| (di, dj)))
                .any(|(di, dj)| !self.sea[g.index(g.neighbor(c, di, dj))]);
            if beside_land {
                coast.push(self.exceeded[k] as f64);
            }
        }
        coast.sort_by(f64::total_cmp);
        Summary {
            mean: if area > 0.0 { mean / area } else { 0.0 },
            coast_exceeded_median: coast.get(coast.len() / 2).copied().unwrap_or(0.0),
            coast_exceeded_most: coast.last().copied().unwrap_or(0.0),
        }
    }

    /// A quantity of the coarser grid's, `values` at each of its cells, read on the finish's
    /// `fine` grid where `sea` (F428): between the four of its cells about each fine cell's
    /// centre, bilinear over those of them the liquid's, else the liquid's nearest; none where
    /// not `sea`.
    pub fn on(&self, fine: CubeSphere, sea: &[bool], values: &[f32]) -> Vec<f32> {
        let g = self.grid;
        (0..fine.cell_count())
            .into_par_iter()
            .map(|k| {
                if !sea[k] {
                    return 0.0;
                }
                let d = fine.center(fine.cell(k));
                let (c, x, y) = g.locate(d);
                let (fx, fy) = (x - c.i as f64 - 0.5, y - c.j as f64 - 0.5);
                let (si, sj) = (if fx < 0.0 { -1 } else { 1 }, if fy < 0.0 { -1 } else { 1 });
                let (ax, ay) = (fx.abs(), fy.abs());
                let (mut sum, mut weight) = (0.0f64, 0.0f64);
                for (di, dj, w) in [
                    (0, 0, (1.0 - ax) * (1.0 - ay)),
                    (si, 0, ax * (1.0 - ay)),
                    (0, sj, (1.0 - ax) * ay),
                    (si, sj, ax * ay),
                ] {
                    let m = g.index(g.neighbor(c, di, dj));
                    if self.sea[m] && w > 0.0 {
                        sum += w * values[m] as f64;
                        weight += w;
                    }
                }
                if weight > 0.0 {
                    (sum / weight) as f32
                } else {
                    self.cell_at(d).map_or(0.0, |m| values[m])
                }
            })
            .collect()
    }

    /// A cell of the coarser grid at the direction `d`, the liquid's nearest where its own is
    /// the land's: its own, else the first of the liquid's cells about it, ring by ring, up to
    /// WAVE_COAST_RINGS rings.
    pub fn cell_at(&self, d: [f64; 3]) -> Option<usize> {
        let g = self.grid;
        let c = g.locate(d).0;
        let k = g.index(c);
        if self.sea[k] {
            return Some(k);
        }
        let mut best: Option<(f64, usize)> = None;
        for r in 1..=WAVE_COAST_RINGS as i32 {
            for di in -r..=r {
                for dj in -r..=r {
                    if di.abs() != r && dj.abs() != r {
                        continue;
                    }
                    let m = g.index(g.neighbor(c, di, dj));
                    if self.sea[m] {
                        let near = dot(g.center(g.cell(m)), d);
                        if best.is_none_or(|b| near > b.0) {
                            best = Some((near, m));
                        }
                    }
                }
            }
            if best.is_some() {
                break;
            }
        }
        best.map(|b| b.1)
    }

    /// The waves' carrying along a shore at cell `c` of the coarser grid, its seaward normal
    /// facing `facing`, rad clockwise from north, under gravity `g`: m³/s of the ground, toward
    /// the left looking out to sea (F429). Ashton and Murray's (2006) deep-water form of the
    /// CERC formula, K1 H^(12/5) T^(1/5) cos^(6/5)(x) sin(x), x the angle the waves come from
    /// clockwise of the normal, over the waves' energy by direction at their net height and
    /// period (Nienhuis, Ashton and Giosan 2015, eqs. 1 to 4); none from the land's side. K1
    /// WAVE_CARRY at the Earth's gravity, and at another as Komar's (1998) g^(3/2) times the
    /// deep-water form's g^(1/10).
    pub fn carried(&self, c: usize, facing: f64, g: f64) -> f64 {
        let k1 = WAVE_CARRY.0 * (g / WAVE_CARRY.1).powf(1.6);
        let base =
            k1 * (self.net_height[c] as f64).powf(2.4) * (self.net_period[c] as f64).powf(0.2);
        let mut q = 0.0;
        for (b, &e) in self.energy[c].iter().enumerate() {
            if e <= 0.0 {
                continue;
            }
            let from = TAU * b as f64 / WAVE_DIRECTIONS as f64;
            let x = (from - facing + PI).rem_euclid(TAU) - PI;
            if x.abs() >= std::f64::consts::FRAC_PI_2 {
                continue;
            }
            q += e as f64 * x.cos().powf(1.2) * x.sin();
        }
        base * q
    }

    /// Hallermeier's (1981) depth of closure at cell `c` of the coarser grid under gravity `g`,
    /// m: 2.28 H_e - 68.5 H_e² / (g T_e²), H_e the height exceeded 12 hours a year and T_e its
    /// period (Udo, Ranasinghe and Takeda 2020, eq. 2), none where it would be less.
    pub fn closure(&self, c: usize, g: f64) -> f64 {
        let (h, t) = (self.exceeded[c] as f64, self.exceeded_period[c] as f64);
        if h <= 0.0 || t <= 0.0 {
            return 0.0;
        }
        (WAVE_CLOSURE.0 * h - WAVE_CLOSURE.1 * h * h / (g * t * t)).max(0.0)
    }

    /// How the load a mouth at cell `m` of the finish's `fine` grid brings the sea `sea` in a
    /// step is laid along its coast (F429): the cells of the coast it goes to and their shares,
    /// none where its river builds its own lobe. Its coast faces the sea's side of the cells
    /// about it; the waves' carrying away from it along each flank as the flank turns from it,
    /// at most `left` and `right`; its fluvial dominance ratio R, its load `rate` m³/s over the
    /// two together (Nienhuis, Ashton and Giosan 2015). At R of one or more its river's lobe;
    /// under it each flank turned until the waves carry away from it R of the most they can,
    /// the river's own share of what the flank carries, at least WAVE_FLANK_LEAST, and the
    /// delta a cuspate one, its flanks from its apex at those angles, standing as deep as the
    /// sea about the mouth, its floor at `ground` under the liquid's `level`, the depth of
    /// closure at least (Ashton and Murray's shoreface, down to where it meets the shelf), its
    /// area holding the span's load, `span`
    /// m³: its apex L out, sqrt(2 V / (D (cot a_l + cot a_r))), each flank L cot a along the
    /// coast. Along the coast's cells, the sea's beside the land, the way along them from the
    /// mouth, the load laid from each by the delta's reach there, 1 - s / W of its flank's.
    #[allow(clippy::too_many_arguments)]
    pub fn delta(
        &self,
        fine: CubeSphere,
        sea: &[bool],
        coast: &[bool],
        (ground, level): (&[f32], f64),
        m: usize,
        (rate, span): (f64, f64),
        (radius, g): (f64, f64),
    ) -> Option<Vec<(usize, f64)>> {
        let p = fine.center(fine.cell(m));
        let c = self.cell_at(p)?;
        // The coast's seaward normal: toward the sea's side of the cells about the mouth.
        let (e, n) = local_frame(p);
        let r = WAVE_FACING_CELLS as i32;
        let cm = fine.cell(m);
        let (mut se, mut sn) = (0.0, 0.0);
        let (mut deep, mut wet) = (0.0, 0.0);
        for di in -r..=r {
            for dj in -r..=r {
                let k = fine.index(fine.neighbor(cm, di, dj));
                let q = fine.center(fine.cell(k));
                let v = [q[0] - p[0], q[1] - p[1], q[2] - p[2]];
                let s = if sea[k] { 1.0 } else { -1.0 };
                if sea[k] {
                    deep += (level - ground[k] as f64).max(0.0);
                    wet += 1.0;
                }
                se += s * dot(v, e);
                sn += s * dot(v, n);
            }
        }
        if se.hypot(sn) <= 0.0 {
            return None;
        }
        let facing = se.atan2(sn).rem_euclid(TAU);
        // The waves' carrying away from the mouth along each flank as it turns from the coast.
        let turns: Vec<f64> = (0..WAVE_FLANK_SAMPLES)
            .map(|i| std::f64::consts::FRAC_PI_2 * i as f64 / (WAVE_FLANK_SAMPLES - 1) as f64)
            .collect();
        let left: Vec<f64> = turns
            .iter()
            .map(|&a| self.carried(c, facing - a, g))
            .collect();
        let right: Vec<f64> = turns
            .iter()
            .map(|&a| -self.carried(c, facing + a, g))
            .collect();
        let most = |v: &[f64]| v.iter().fold(0.0f64, |a, &b| a.max(b));
        let (most_l, most_r) = (most(&left), most(&right));
        if most_l + most_r <= 0.0 {
            return None;
        }
        let ratio = rate / (most_l + most_r);
        if ratio >= 1.0 {
            return None;
        }
        // Each flank turned until the waves carry away R of the most they can along it.
        let angle = |v: &[f64], top: f64| -> Option<f64> {
            if top <= 0.0 {
                return None;
            }
            let want = ratio * top;
            for i in 0..v.len() {
                if v[i] >= want {
                    let a = if i == 0 {
                        turns[0]
                    } else {
                        let f = ((want - v[i - 1]) / (v[i] - v[i - 1]).max(1e-300)).clamp(0.0, 1.0);
                        turns[i - 1] + f * (turns[i] - turns[i - 1])
                    };
                    return Some(a.max(WAVE_FLANK_LEAST));
                }
            }
            None
        };
        let (a_l, a_r) = (angle(&left, most_l), angle(&right, most_r));
        let cot = |a: Option<f64>| a.map_or(0.0, |a| 1.0 / a.tan());
        let (cot_l, cot_r) = (cot(a_l), cot(a_r));
        let depth = self
            .closure(c, g)
            .max(if wet > 0.0 { deep / wet } else { 0.0 });
        if depth <= 0.0 || cot_l + cot_r <= 0.0 {
            return None;
        }
        let apex = (2.0 * span / (depth * (cot_l + cot_r))).sqrt();
        let (w_l, w_r) = (apex * cot_l, apex * cot_r);
        let reach = w_l.max(w_r);
        // The way along the coast from the mouth, each cell's side by the coast's tangent.
        let along = facing + std::f64::consts::FRAC_PI_2;
        let to_right = [
            along.cos() * n[0] + along.sin() * e[0],
            along.cos() * n[1] + along.sin() * e[1],
            along.cos() * n[2] + along.sin() * e[2],
        ];
        let mut dist: std::collections::HashMap<usize, f32> = std::collections::HashMap::new();
        let mut heap = std::collections::BinaryHeap::new();
        dist.insert(m, 0.0);
        heap.push(std::cmp::Reverse((Lowest(0.0, m as u32), m)));
        let mut shares = Vec::new();
        while let Some(std::cmp::Reverse((Lowest(s, _), k))) = heap.pop() {
            if dist.get(&k).is_some_and(|&d| d < s) {
                continue;
            }
            let s = s as f64;
            let q = fine.center(fine.cell(k));
            let v = [q[0] - p[0], q[1] - p[1], q[2] - p[2]];
            let w_side = if dot(v, to_right) >= 0.0 { w_r } else { w_l };
            let w = if k == m {
                1.0
            } else if w_side > 0.0 {
                (1.0 - s / w_side).max(0.0)
            } else {
                0.0
            };
            if w > 0.0 {
                shares.push((k, w));
            }
            let ck = fine.cell(k);
            for (di, dj) in AROUND {
                let j = fine.index(fine.neighbor(ck, di, dj));
                if !coast[j] {
                    continue;
                }
                let step = dot(q, fine.center(fine.cell(j))).clamp(-1.0, 1.0).acos() * radius;
                let next = s + step;
                let near = next as f32;
                if next <= reach && dist.get(&j).is_none_or(|&d| near < d) {
                    dist.insert(j, near);
                    heap.push(std::cmp::Reverse((Lowest(next as f32, j as u32), j)));
                }
            }
        }
        let total: f64 = shares.iter().map(|s| s.1).sum();
        if shares.len() < 2 || total <= 0.0 {
            return None;
        }
        Some(shares.into_iter().map(|(k, w)| (k, w / total)).collect())
    }

    /// The load each cell of the sea `sea` on the finish's `fine` grid, its floor at `ground`
    /// under the liquid's `level`, takes in a step of
    /// `step` s, `mouths` m³, laid along the coasts by the waves (F429), the span's load each
    /// mouth's the span's `steps` steps of it: each mouth's as its delta lays it (`delta`).
    /// Returns the mouths, those whose loads the waves laid along their coasts and the volume
    /// they moved, m³.
    pub fn spread(
        &self,
        fine: CubeSphere,
        sea: &[bool],
        (ground, level): (&[f32], f64),
        mouths: &mut [f64],
        (step, steps): (f64, usize),
        (radius, g): (f64, f64),
    ) -> (usize, usize, f64) {
        let n = fine.cell_count();
        let coast: Vec<bool> = (0..n)
            .into_par_iter()
            .map(|k| {
                sea[k] && {
                    let c = fine.cell(k);
                    AROUND
                        .iter()
                        .any(|&(di, dj)| !sea[fine.index(fine.neighbor(c, di, dj))])
                }
            })
            .collect();
        let sources: Vec<usize> = (0..n).filter(|&k| sea[k] && mouths[k] > 0.0).collect();
        let plans: Vec<Option<Vec<(usize, f64)>>> = sources
            .par_iter()
            .map(|&m| {
                self.delta(
                    fine,
                    sea,
                    &coast,
                    (ground, level),
                    m,
                    (mouths[m] / step, mouths[m] * steps as f64),
                    (radius, g),
                )
            })
            .collect();
        let mut laid = mouths.to_vec();
        let (mut waved, mut moved) = (0, 0.0);
        for (&m, plan) in sources.iter().zip(plans) {
            let Some(plan) = plan else { continue };
            let load = mouths[m];
            laid[m] -= load;
            for (k, w) in plan {
                laid[k] += load * w;
                if k != m {
                    moved += load * w;
                }
            }
            waved += 1;
        }
        mouths.copy_from_slice(&laid);
        (sources.len(), waved, moved)
    }
}

/// The cells about a cell the waves carry the load along (F429).
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

/// A way's length and a cell, the shortest first, then the lower index (F429).
#[derive(Clone, Copy, PartialEq)]
struct Lowest(f32, u32);

impl Eq for Lowest {}

impl PartialOrd for Lowest {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Lowest {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.0.total_cmp(&other.0).then(self.1.cmp(&other.1))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// F428: the manual's growth rises as the fetch's square root, its period as its cube
    /// root, and meets the fully developed sea at the fetch `developed_fetch` gives, no higher
    /// past it; a 10 m/s wind's fully developed sea some 3.1 m and 9.3 s.
    #[test]
    fn waves_grow_with_their_fetch_to_the_fully_developed_sea() {
        let g = 9.81;
        let (h1, t1) = grown(10.0, 1.0e4, g);
        let (h4, t4) = grown(10.0, 4.0e4, g);
        assert!((h4 / h1 - 2.0).abs() < 1e-9, "{h1} {h4}");
        assert!((t4 / t1 - 4f64.cbrt()).abs() < 1e-9, "{t1} {t4}");
        let full = developed_fetch(10.0, g);
        let (hf, _) = grown(10.0, full, g);
        let (hp, tp) = grown(10.0, 10.0 * full, g);
        println!(
            "10 m/s: developed past {:.0} km at {hf:.3} m, {tp:.2} s",
            full / 1e3
        );
        assert!((hf - hp).abs() < 1e-9 && (hp - 3.126).abs() < 0.01 && (tp - 9.31).abs() < 0.01);
        assert_eq!(grown(0.0, 1.0e5, g), (0.0, 0.0));
    }

    /// F428: a storm's sea grows only as long as it blows from one way: the weather's gusts in
    /// the Earth's middle latitudes, at some 7 m/s, turn within a day or two, over which a
    /// 23 m/s gale over an open ocean raises a sea short of its fully developed one; a body
    /// without spin's gusts blow from one way for ever.
    #[test]
    fn a_storm_s_sea_grows_as_long_as_it_blows_from_one_way() {
        let (g, day) = (9.81, 86_400.0);
        let earth = gust_duration(45f64.to_radians(), 7.292e-5, 6.371e6, 7.0);
        let (storm, _) = grown(23.0, duration_fetch(23.0, earth, g), g);
        let (full, _) = grown(23.0, 1.0e9, g);
        println!(
            "the Earth's gusts blow from one way {:.2} days; a 23 m/s gale raises {storm:.1} m in that time, {full:.1} m fully developed",
            earth / day
        );
        assert!((0.5 * day..2.0 * day).contains(&earth), "{earth}");
        assert!(storm < 0.6 * full && storm > 5.0, "{storm} {full}");
        assert!(gust_duration(0.5, 0.0, 6.371e6, 7.0).is_infinite());
        assert!(duration_fetch(10.0, f64::INFINITY, g).is_infinite());
    }

    /// A straight coast along the meridian through [0, 1, 0], the sea east of it, on a fine
    /// grid, its waves all from one direction of WAVE_DIRECTIONS at 2 m and 8 s net, 4 m and
    /// 10 s exceeded 12 hours a year; and a mouth on the coast at [0, 1, 0].
    fn straight_coast(from: usize) -> (CubeSphere, Vec<bool>, Waves, usize) {
        let fine = CubeSphere::new(128);
        let sea: Vec<bool> = (0..fine.cell_count())
            .map(|k| fine.center(fine.cell(k))[0] < 0.0)
            .collect();
        let grid = CubeSphere::new(32);
        let n = grid.cell_count();
        let coarse: Vec<bool> = (0..n).map(|k| grid.center(grid.cell(k))[0] < 0.0).collect();
        let mut energy = [0.0f32; WAVE_DIRECTIONS];
        energy[from] = 1.0;
        let waves = Waves {
            grid,
            sea: coarse,
            mean: vec![1.5; n],
            exceeded: vec![4.0; n],
            exceeded_period: vec![10.0; n],
            net_height: vec![2.0; n],
            net_period: vec![8.0; n],
            energy: vec![energy; n],
            fetch: vec![[1.0e6; WAVE_DIRECTIONS]; n],
        };
        let mouth = (0..fine.cell_count())
            .filter(|&k| sea[k])
            .max_by(|&a, &b| {
                let (pa, pb) = (fine.center(fine.cell(a)), fine.center(fine.cell(b)));
                pa[1].total_cmp(&pb[1]).then(b.cmp(&a))
            })
            .unwrap();
        (fine, sea, waves, mouth)
    }

    /// F429: waves square onto a straight coast spread a small river's load along it, as much
    /// to each side; a river bringing more than they can carry away keeps its lobe; waves from
    /// the north-east carry the load to the south, down the coast.
    #[test]
    fn the_waves_lay_a_river_s_load_along_its_coast() {
        let (g, radius, year) = (9.81, 6.371e6, 3.156e7);
        let share_south = |from: usize, rate: f64| {
            let (fine, sea, waves, mouth) = straight_coast(from);
            let mut mouths = vec![0.0f64; fine.cell_count()];
            let ground: Vec<f32> = sea
                .iter()
                .map(|&s| if s { -100.0 } else { 100.0 })
                .collect();
            let step = 2.0e5 * year;
            mouths[mouth] = rate * step;
            let (_, waved, moved) = waves.spread(
                fine,
                &sea,
                (&ground, 0.0),
                &mut mouths,
                (step, 10),
                (radius, g),
            );
            let z = fine.center(fine.cell(mouth))[2];
            let (mut south, mut north, mut cells) = (0.0, 0.0, 0);
            for (k, &v) in mouths.iter().enumerate().filter(|(_, v)| **v > 0.0) {
                cells += 1;
                let dz = fine.center(fine.cell(k))[2] - z;
                if dz < -1e-9 {
                    south += v;
                } else if dz > 1e-9 {
                    north += v;
                }
            }
            println!(
                "waves from bin {from}, {rate} m³/s: {waved} laid along the coast, {moved:.3e} m³ moved, over {cells} cells, {:.2} of it south of the mouth, {:.2} north",
                south / (rate * 2.0e5 * year),
                north / (rate * 2.0e5 * year)
            );
            (
                waved,
                south / (rate * 2.0e5 * year),
                north / (rate * 2.0e5 * year),
                cells,
            )
        };
        let east = WAVE_DIRECTIONS / 4;
        let (waved, south, north, cells) = share_south(east, 0.05);
        assert!(waved == 1 && cells > 10, "{waved} {cells}");
        assert!(
            (south - north).abs() < 0.1 * (south + north),
            "{south} {north}"
        );
        let (waved, ..) = share_south(east, 1.0);
        assert_eq!(waved, 0);
        let (waved, south, north, _) = share_south(east / 2, 0.05);
        assert!(waved == 1 && south > 2.0 * north, "{south} {north}");
    }

    /// F428: across a round sea the fetch from its centre is its radius in every direction,
    /// within a cell, and none on its land.
    #[test]
    fn the_fetch_across_a_round_sea_is_its_radius() {
        let grid = CubeSphere::new(64);
        let centre = [0.0f64, 0.6, 0.8];
        let r = 0.3;
        let sea: Vec<bool> = (0..grid.cell_count())
            .map(|k| {
                dot(grid.center(grid.cell(k)), centre)
                    .clamp(-1.0, 1.0)
                    .acos()
                    < r
            })
            .collect();
        let radius = 1.0e6;
        let f = fetches(grid, &sea, radius, &vec![1.0e7; grid.cell_count()]);
        let k = grid.index(grid.locate(centre).0);
        let cell = grid.delta() * radius;
        let worst = f[k]
            .iter()
            .map(|&x| (x as f64 - r * radius).abs())
            .fold(0.0, f64::max);
        println!(
            "the fetch from the centre {:?} m, the radius {:.0} m, a cell {cell:.0} m",
            f[k],
            r * radius
        );
        assert!(worst < 1.5 * cell, "{worst}");
        assert!(
            (0..grid.cell_count())
                .filter(|&k| !sea[k])
                .all(|k| f[k] == [0.0; WAVE_DIRECTIONS])
        );
    }

    /// F428: a steady wind from the west over a long fetch raises its fully developed sea all
    /// the year, all its energy from the west; ice over half the year halves the mean and
    /// leaves the height exceeded 12 hours a year as it was; ice open only by its leads, a
    /// hundredth of it, the waves growing over that much of the fetch, next to none; the
    /// weather's gusts raise that height over the mean.
    #[test]
    fn a_place_s_waves_follow_its_wind_its_ice_and_its_weather() {
        let g = 9.81;
        let steps = 12;
        let year = |open: &dyn Fn(usize) -> f64, gust: f64| Year {
            east: vec![10.0; steps],
            north: vec![0.0; steps],
            gust: vec![gust; steps],
            duration: vec![f64::INFINITY; steps],
            open: (0..steps).map(open).collect(),
        };
        let fetch = [1.0e7f32; WAVE_DIRECTIONS];
        let (mean, exceeded, _, _, _, energy) = at_place(&year(&|_| 1.0, 0.0), &fetch, g);
        let west = 3 * WAVE_DIRECTIONS / 4;
        println!(
            "steady: mean {mean:.3} m, exceeded {exceeded:.3} m, the west's share {:.3}",
            energy[west]
        );
        assert!((mean - 3.126).abs() < 0.01 && (exceeded - mean).abs() < 1e-6);
        assert!((energy[west] - 1.0).abs() < 1e-6);
        let seasons = |k: usize| if k < steps / 2 { 1.0 } else { 0.0 };
        let (half, half_exceeded, ..) = at_place(&year(&seasons, 0.0), &fetch, g);
        assert!((half / mean - 0.5).abs() < 1e-6 && (half_exceeded - exceeded).abs() < 1e-3);
        let short = [1.0e5f32; WAVE_DIRECTIONS];
        let (leads, ..) = at_place(&year(&|_| 0.01, 0.0), &short, g);
        println!(
            "seasonal ice: mean {half:.3} m, exceeded {half_exceeded:.3} m; leads: mean {leads:.3} m"
        );
        assert!(leads < 0.15 * mean, "{leads}");
        let (gusty, gusty_exceeded, ..) = at_place(&year(&|_| 1.0, 6.0), &fetch, g);
        println!("gusts of 6 m/s: mean {gusty:.3} m, exceeded {gusty_exceeded:.3} m");
        assert!(gusty_exceeded > 1.5 * gusty, "{gusty} {gusty_exceeded}");
    }
}
