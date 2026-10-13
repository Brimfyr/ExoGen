// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! Partial resurfacing of an ice shell (PLAN.md 4.3 and Phase 4): grooved terrain and a
//! focused polar terrain.
//!
//! A shell that has stretched renews part of its surface as lanes of new ice: swaths of
//! parallel grooves, bright and clean, cutting the old dark cratered terrain into
//! polygons, as two thirds of Ganymede were renewed. A lane is a swath about a centre
//! line of random length and width, traced along the tidal field where the shell cracks
//! (so it curves as the field turns, as Enceladus's sulci do) and along a great circle
//! where nothing steers it (Ganymede's), closing in a lens at its ends; lanes are drawn
//! until they cover the resurfaced share of the surface. Inside them the old relief is
//! smoothed away at their time and sinks by the extension's subsidence, spread over the
//! lithosphere's flexural parameter so that no lane ends on a step; the surface age
//! resets, and the finish cuts the grooves. A shell over an ocean whose tidal heating focuses under a
//! pole is thin there, floats lower there, and renews that polar terrain continuously, so
//! it carries no craters, most of the cracks, and folds between them, as Enceladus's
//! south polar terrain does.

use std::collections::HashMap;

use rayon::prelude::*;

use crate::constants::*;
use crate::crater::{CellWrite, StampOutput, Stamper, with_stamper};
use crate::fields::{Field, LayerSet};
use crate::geom::{angle_between, cross, dot, normalize};
use crate::grid::CubeSphere;
use crate::noise::GradientNoise;
use crate::params::ParamGraph;
use crate::rng::Rng;
use crate::stages::ice::{FocusEpisode, ShellMap};
use crate::stages::lineae::{Cracks, Linea};
use crate::units::*;

pub const VERSION: u64 = 12;

/// A lane of renewed terrain: a swath about a centre line traced across the shell,
/// along the tidal field where the shell cracks, so that it curves as the field turns
/// (Enceladus's sulci), and along a great circle where nothing steers it (Ganymede's);
/// its width wobbles along it and closes in a lens at its ends.
#[derive(Clone, Debug)]
pub struct Lane {
    /// The centre line: points `spacing` apart, from one end to the other.
    pub points: Vec<[f64; 3]>,
    pub spacing: f64,
    /// Half the centre line's length, m.
    pub half_length: f64,
    pub half_width: f64,
    /// Whether the lane is grooved or smooth bright terrain.
    pub grooved: bool,
    /// How clean the lane's ice is, 0 to 1: lanes differ in the dirt the new ice took up.
    pub clean: f64,
    /// The width's wander along the lane, a self-affine profile sampled from one end to
    /// the other, as a factor on the half width.
    pub width_profile: Vec<f32>,
    /// When the lane formed, Gyr since formation.
    pub time_gyr: f64,
    pub id: u64,
}

/// A cell near a lane: its signed distance across the centre line, m, its signed
/// distance along it from the midpoint, m, and the centre line's tangent there.
#[derive(Clone, Copy, Debug)]
pub struct LaneSample {
    pub idx: usize,
    pub across: f64,
    pub along: f64,
    pub tangent: [f64; 3],
}

impl Lane {
    /// The points of a centre line along the great circle with pole `axis` through
    /// `center`, `half_length` each way at `spacing`.
    pub fn great_circle(
        axis: [f64; 3],
        center: [f64; 3],
        half_length: f64,
        spacing: f64,
        body_radius: f64,
    ) -> Vec<[f64; 3]> {
        let t = normalize(cross(axis, center));
        let n = (half_length / spacing.max(1.0)).ceil().max(1.0) as i64;
        (-n..=n)
            .map(|i| {
                let (s, c) = (i as f64 * spacing / body_radius).sin_cos();
                normalize([
                    center[0] * c + t[0] * s,
                    center[1] * c + t[1] * s,
                    center[2] * c + t[2] * s,
                ])
            })
            .collect()
    }

    /// A lane about a centre line, its edges wandering in a self-affine profile from
    /// three widths down to half a width.
    pub fn from_points(
        points: Vec<[f64; 3]>,
        spacing: f64,
        half_width: f64,
        grooved: bool,
        clean: f64,
        time_gyr: f64,
        id: u64,
    ) -> Self {
        let half_length = 0.5 * points.len().saturating_sub(1) as f64 * spacing;
        let noise = GradientNoise::new(id ^ 0x6700_D1F7);
        let samples = 128usize;
        let width = 2.0 * half_width;
        let gain = 2f64.powf(-ROUGHNESS_HURST);
        let raw: Vec<f64> = (0..=samples)
            .map(|i| {
                let along = -half_length + 2.0 * half_length * i as f64 / samples as f64;
                noise.fbm([along / (3.0 * width), 0.3, 0.7], 1.0, 3, 2.0, gain)
            })
            .collect();
        let rms = (raw.iter().map(|v| v * v).sum::<f64>() / raw.len() as f64)
            .sqrt()
            .max(1e-9);
        let width_profile = raw
            .iter()
            .map(|v| (1.0 + GROOVE_WIDTH_WOBBLE * v / rms).max(0.3) as f32)
            .collect();
        Self {
            points,
            spacing,
            half_length,
            half_width,
            grooved,
            clean,
            width_profile,
            time_gyr,
            id,
        }
    }

    /// The lane's midpoint.
    pub fn center(&self) -> [f64; 3] {
        self.points[self.points.len() / 2]
    }

    /// The centre line's tangent at one of its points.
    fn tangent_at(&self, i: usize) -> [f64; 3] {
        let n = self.points.len();
        let a = self.points[i.saturating_sub(1)];
        let b = self.points[(i + 1).min(n - 1)];
        let p = self.points[i];
        let chord = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let radial = dot(chord, p);
        normalize([
            chord[0] - radial * p[0],
            chord[1] - radial * p[1],
            chord[2] - radial * p[2],
        ])
    }

    /// Signed distance along the centre line from the midpoint to one of its points, m.
    fn along_of(&self, i: usize) -> f64 {
        (i as f64 - 0.5 * (self.points.len() - 1) as f64) * self.spacing
    }

    /// The nearest point of the centre line to a direction: the signed distance across
    /// the line, m, the signed distance along it from the midpoint, m, and the tangent.
    pub fn locate(&self, d: [f64; 3], body_radius: f64) -> (f64, f64, [f64; 3]) {
        let (mut best, mut best_dot) = (0usize, f64::NEG_INFINITY);
        for (i, p) in self.points.iter().enumerate() {
            let x = dot(*p, d);
            if x > best_dot {
                best_dot = x;
                best = i;
            }
        }
        let p = self.points[best];
        let t = self.tangent_at(best);
        let normal = normalize(cross(p, t));
        let across = dot(d, normal).clamp(-1.0, 1.0).asin() * body_radius;
        let along = self.along_of(best) + dot(d, t).clamp(-1.0, 1.0).asin() * body_radius;
        (across, along, t)
    }

    /// The lane's half width at a signed distance along it, m: the mean times the
    /// wander profile there, closing over a width at each end as the bounding faults die
    /// out.
    pub fn half_width_at(&self, along: f64) -> f64 {
        let s = (0.5 + 0.5 * along / self.half_length.max(1.0)).clamp(0.0, 1.0);
        let n = self.width_profile.len();
        let x = s * (n - 1) as f64;
        let i = (x.floor() as usize).min(n.saturating_sub(2));
        let f = (x - i as f64).clamp(0.0, 1.0);
        let wander =
            self.width_profile[i] as f64 * (1.0 - f) + self.width_profile[i + 1] as f64 * f;
        let to_end = (self.half_length - along.abs()).max(0.0);
        let taper = (to_end / (GROOVE_END_TAPER_WIDTHS * 2.0 * self.half_width).max(1.0))
            .clamp(0.0, 1.0)
            .powf(GROOVE_LANE_TAPER_POWER);
        self.half_width * wander * taper
    }

    /// How much of a cell the lane covers at an across and an along, 0 to 1: a one-cell
    /// ramp across the lane's edges and ends, for edges drawn without stair steps.
    pub fn coverage_at(&self, across: f64, along: f64, cell: f64) -> f64 {
        let ramp = GROOVE_EDGE_CELLS * cell.max(1.0);
        let side = (self.half_width_at(along) - across.abs()) / ramp + 0.5;
        let end = (self.half_length - along.abs()) / ramp + 0.5;
        side.clamp(0.0, 1.0) * end.clamp(0.0, 1.0)
    }

    pub fn inside(&self, d: [f64; 3], body_radius: f64) -> bool {
        let (across, along, _) = self.locate(d, body_radius);
        along.abs() < self.half_length && across.abs() < self.half_width_at(along)
    }

    /// Angular radius of a cap round the midpoint that holds the lane.
    pub fn cap_radius(&self, body_radius: f64) -> f64 {
        let c = self.center();
        let reach = self
            .points
            .iter()
            .map(|p| angle_between(c, *p))
            .fold(0.0, f64::max);
        reach + (self.half_width * (1.0 + GROOVE_WIDTH_WOBBLE) + 2.0 * self.spacing) / body_radius
    }

    /// The cells within `margin` of the lane with their across and along: a disc round
    /// each point of the centre line, each cell kept with the point it lies nearest
    /// along, in cell order.
    pub fn sample(
        &self,
        grid: &CubeSphere,
        body_radius: f64,
        margin: f64,
        stamper: &mut Stamper,
    ) -> Vec<LaneSample> {
        let half = 0.5 * self.spacing;
        let reach = (self.half_width * (1.0 + GROOVE_WIDTH_WOBBLE) + margin).hypot(half + margin);
        let mut best: HashMap<usize, (f64, LaneSample)> = HashMap::new();
        for i in 0..self.points.len() {
            let p = self.points[i];
            let t = self.tangent_at(i);
            let normal = normalize(cross(p, t));
            let along0 = self.along_of(i);
            stamper.gather(grid, body_radius, p, reach);
            for &(idx, _, _) in &stamper.footprint {
                let d = grid.center(grid.cell(idx));
                let local = dot(d, t).clamp(-1.0, 1.0).asin() * body_radius;
                if local.abs() > half + margin {
                    continue;
                }
                let across = dot(d, normal).clamp(-1.0, 1.0).asin() * body_radius;
                let sample = LaneSample {
                    idx,
                    across,
                    along: along0 + local,
                    tangent: t,
                };
                let e = best.entry(idx).or_insert((f64::INFINITY, sample));
                if local.abs() < e.0 {
                    *e = (local.abs(), sample);
                }
            }
        }
        let mut out: Vec<LaneSample> = best.into_values().map(|(_, s)| s).collect();
        out.sort_by_key(|s| s.idx);
        out
    }
}

/// The fissures' points bucketed on the sphere, for the folds to follow the nearest one.
struct FabricIndex {
    scale: f64,
    buckets: HashMap<(i32, i32, i32), Vec<FabricPoint>>,
}

/// A fissure point and the tangent there.
type FabricPoint = ([f64; 3], [f64; 3]);

impl FabricIndex {
    fn new(cracks: &[Linea], range_rad: f64) -> Self {
        let scale = 1.0 / range_rad.max(1e-6);
        let mut buckets: HashMap<(i32, i32, i32), Vec<FabricPoint>> = HashMap::new();
        for l in cracks {
            let n = l.points.len();
            for i in 0..n {
                let a = l.points[i.saturating_sub(1)];
                let b = l.points[(i + 1).min(n - 1)];
                let p = l.points[i];
                let chord = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
                let radial = dot(chord, p);
                let t = normalize([
                    chord[0] - radial * p[0],
                    chord[1] - radial * p[1],
                    chord[2] - radial * p[2],
                ]);
                let key = (
                    (p[0] * scale).floor() as i32,
                    (p[1] * scale).floor() as i32,
                    (p[2] * scale).floor() as i32,
                );
                buckets.entry(key).or_default().push((p, t));
            }
        }
        Self { scale, buckets }
    }

    /// The distance across the fissures within a bucket's reach of a direction, rad: a
    /// soft minimum, over `softness` rad, of the perpendicular distances to each point's
    /// local line, each point counted by a window on its straight distance of `window`
    /// rad; and the true straight minimum. Measured across the line rather than to the
    /// point, the ridges run parallel past a fissure's end instead of wrapping round it.
    fn soft_distance(&self, d: [f64; 3], softness: f64, window: f64) -> Option<(f64, f64)> {
        let (i, j, k) = (
            (d[0] * self.scale).floor() as i32,
            (d[1] * self.scale).floor() as i32,
            (d[2] * self.scale).floor() as i32,
        );
        let softness = softness.max(1e-9);
        let window = window.max(1e-9);
        let (mut sum, mut min) = (0.0f64, f64::INFINITY);
        for di in -1..=1 {
            for dj in -1..=1 {
                for dk in -1..=1 {
                    let Some(v) = self.buckets.get(&(i + di, j + dj, k + dk)) else {
                        continue;
                    };
                    for &(p, t) in v {
                        let dist = angle_between(p, d);
                        let n = normalize(cross(p, t));
                        let perp = dot([d[0] - p[0], d[1] - p[1], d[2] - p[2]], n).abs();
                        let w = (-(dist / window).powi(2)).exp();
                        sum += w * (-perp / softness).exp();
                        min = min.min(dist);
                    }
                }
            }
        }
        (sum > 1e-300).then(|| (-softness * sum.ln(), min))
    }
}

/// How much of a structure of width `w`, m, survives gardening diffusion of `kappa_t`,
/// m², as a factor on its height: it spreads and lowers as the craters do.
fn degraded(w: f64, kappa_t: f64) -> f64 {
    (w * w / (w * w + 4.0 * kappa_t.max(0.0))).sqrt()
}

/// The flexural parameter of an elastic ice lithosphere `te` thick over an ocean, m: the
/// length over which a load's deflection spreads (Turcotte and Schubert, (4 D / rho g)^1/4
/// with D = E te^3 / 12 (1 - nu^2)).
pub fn flexural_parameter(te: f64, gravity: f64) -> f64 {
    let rigidity =
        ICE_YOUNG_MODULUS * te.max(1.0).powi(3) / (12.0 * (1.0 - ICE_POISSON * ICE_POISSON));
    (4.0 * rigidity / (DENSITY_ICE * gravity.max(1e-3))).powf(0.25)
}

/// The strike of the shell's extension at a point (F44): the direction lanes run in
/// where no crack field steers them, the regional pattern the stress alone does not
/// resolve, as an angle from east that turns smoothly over the body (a field of a few
/// cells about it; not the contours of a scalar, which would close in loops about its
/// extrema) and drifts through the episode by `phase` (0 at the start, 1 at the end) so
/// that generations cross-cut.
pub fn extension_heading(noise: &GradientNoise, d: [f64; 3], phase: f64) -> [f64; 3] {
    let drift = GROOVE_HEADING_DRIFT * phase;
    let v = noise.fbm(
        [d[0] + drift, d[1] - 0.5 * drift, d[2]],
        GROOVE_HEADING_FREQUENCY,
        2,
        2.0,
        0.5,
    ) / 0.6;
    let theta = std::f64::consts::PI * v.clamp(-1.0, 1.0);
    let (east, north) = crate::geom::local_basis(d);
    normalize([
        east[0] * theta.cos() + north[0] * theta.sin(),
        east[1] * theta.cos() + north[1] * theta.sin(),
        east[2] * theta.cos() + north[2] * theta.sin(),
    ])
}

/// Rodrigues rotation of a vector about a unit axis.
fn rotate(v: [f64; 3], k: [f64; 3], angle: f64) -> [f64; 3] {
    let (s, c) = angle.sin_cos();
    let kv = cross(k, v);
    let kd = dot(k, v) * (1.0 - c);
    [
        v[0] * c + kv[0] * s + k[0] * kd,
        v[1] * c + kv[1] * s + k[1] * kd,
        v[2] * c + kv[2] * s + k[2] * kd,
    ]
}

/// A coarse field from its values.
pub fn field_from(grid: CubeSphere, values: &[f32]) -> Field<f32> {
    let mut f = Field::new(grid);
    f.as_mut_slice().copy_from_slice(values);
    f
}

/// The renewed polar terrain on a grid: the thinned shell's core and the ground within
/// the spread width of the fissures that opened in the thinned zone, wherever the shell
/// is thinned at all; so the outline follows the stripes and lobes where they end.
pub fn polar_mask(
    grid: &CubeSphere,
    body_radius: f64,
    shell: &ShellMap,
    cracks: &[Linea],
) -> Vec<bool> {
    let n = grid.cell_count();
    if !shell.has_active_terrain() {
        return vec![false; n];
    }
    let mut mask: Vec<bool> = (0..n)
        .map(|k| shell.profile(grid.center(grid.cell(k))) > POLAR_CORE_PROFILE)
        .collect();
    with_stamper(grid, |stamper| {
        for l in cracks.iter().filter(|l| shell.thin(l.midpoint())) {
            let every = ((0.25 * POLAR_SPREAD_WIDTH / l.spacing.max(1.0)) as usize).max(1);
            for p in l.points.iter().step_by(every) {
                // The terrain follows the fissure out of the core as a tongue that
                // narrows to nothing where the taper falls to the outer profile, and
                // ends where the discs along the fissure would no longer meet.
                let reach = POLAR_SPREAD_WIDTH
                    * ((shell.profile(*p) - POLAR_OUTER_PROFILE)
                        / (POLAR_CORE_PROFILE - POLAR_OUTER_PROFILE))
                        .clamp(0.0, 1.0);
                if reach < 0.5 * every as f64 * l.spacing {
                    continue;
                }
                stamper.gather(grid, body_radius, *p, reach);
                for &(k, _, _) in &stamper.footprint {
                    if shell.profile(grid.center(grid.cell(k))) > POLAR_OUTER_PROFILE {
                        mask[k] = true;
                    }
                }
            }
        }
    });
    mask
}

/// Signed distance to a mask's edge in cells, positive inside, by breadth-first search
/// from the edge out to a limit; beyond it the limit itself.
pub fn edge_distance(grid: &CubeSphere, mask: &[bool], limit: usize) -> Vec<f32> {
    let n = mask.len();
    let neighbours = |k: usize| {
        let c = grid.cell(k);
        [(1, 0), (-1, 0), (0, 1), (0, -1)].map(|(di, dj)| grid.index(grid.neighbor(c, di, dj)))
    };
    let mut dist = vec![f32::NAN; n];
    let mut queue = std::collections::VecDeque::new();
    for k in 0..n {
        if neighbours(k).iter().any(|&nb| mask[nb] != mask[k]) {
            dist[k] = if mask[k] { 0.5 } else { -0.5 };
            queue.push_back(k);
        }
    }
    while let Some(k) = queue.pop_front() {
        let next = dist[k].abs() + 1.0;
        if next > limit as f32 {
            continue;
        }
        for nb in neighbours(k) {
            if dist[nb].is_nan() {
                dist[nb] = if mask[nb] { next } else { -next };
                queue.push_back(nb);
            }
        }
    }
    for (k, d) in dist.iter_mut().enumerate() {
        if d.is_nan() {
            *d = if mask[k] {
                limit as f32
            } else {
                -(limit as f32)
            };
        }
    }
    // The search counts steps along the grid, so its contours are diamonds; two passes of
    // smoothing round them.
    crate::stages::episodic::smooth(grid, &dist, 2)
}

/// How far from the polar terrain's edge, in cells, the finish keeps a distance: the
/// margin ridge's reach, or the margin ramp's, whichever is further.
pub fn polar_edge_limit(cell: f64) -> usize {
    ((4.0 * POLAR_MARGIN_WIDTH).max(POLAR_MARGIN_RAMP) / cell.max(1.0)).ceil() as usize + 2
}

/// How much of a renewed terrain a cell takes at a signed distance from its edge, in
/// cells: a smooth ramp over the margin width centred on the edge, 0.5 at the edge
/// itself, so the floor and the folds fade into the old surface and the old craters sink
/// part way across it instead of ending on a cut.
pub fn margin_coverage(edge_cells: f32, cell: f64) -> f64 {
    let x = (edge_cells as f64 * cell / POLAR_MARGIN_RAMP + 0.5).clamp(0.0, 1.0);
    x * x * (3.0 - 2.0 * x)
}

/// A mask eroded by one cell: a cell stays only with all four neighbours. The coarse loop
/// renews these; the finish draws the true edge around them at its own resolution, so
/// the coarse cells never reach past it.
fn eroded(grid: &CubeSphere, mask: &[bool]) -> Vec<bool> {
    (0..mask.len())
        .into_par_iter()
        .map(|k| {
            mask[k] && {
                let c = grid.cell(k);
                [(1, 0), (-1, 0), (0, 1), (0, -1)]
                    .iter()
                    .all(|&(di, dj)| mask[grid.index(grid.neighbor(c, di, dj))])
            }
        })
        .collect()
}

/// The shell's partial resurfacing through the loop.
pub struct IceShellProcess {
    grid: CubeSphere,
    body_radius: f64,
    pub shell: ShellMap,
    mean_shell: f64,
    /// The lanes of grooved terrain, each with its own time within the episode.
    pub lanes: Vec<Lane>,
    /// The share of the surface the lanes cover.
    pub lane_coverage: f64,
    /// The equilibrium surface of the thinned shell, which relaxation is not to take for
    /// a load.
    equilibrium: Option<Vec<f32>>,
    /// The polar terrain's renewed surface as the last relaxation step left it, and the
    /// lanes' floor at their time: the old relief smoothed and dropped.
    pub polar_floor: Option<Vec<f32>>,
    pub lane_floor: Option<Vec<f32>>,
    active: Vec<bool>,
    depressed: bool,
    /// The thermal history, for the elastic lithosphere at each lane's time.
    thermal: crate::thermal::ThermalHistory,
    gravity: f64,
    /// When the present episode of focusing began, Gyr since formation: 0 where the
    /// thinning is steady, the episode's length before the present where it renews.
    present_start: f64,
    /// The past episodes of focused heating (F53), their shell maps centred where the
    /// body left their terrain, the terrain each renewed, and its floor when it ended.
    pub episodes: Vec<FocusEpisode>,
    episode_maps: Vec<ShellMap>,
    episode_masks: Vec<Vec<bool>>,
    pub episode_floors: Vec<Option<Vec<f32>>>,
    episode_sunk: Vec<bool>,
    episode_rebounded: Vec<bool>,
    /// The terrain renewed at the current step: the present pole's while its episode
    /// lasts, and each past episode's through its own.
    active_now: Vec<bool>,
    /// The share of the surface the lanes, the present polar terrain and the past
    /// episodes' terrains renewed between them, by area: what is not old cratered ground,
    /// whatever craters it has taken since.
    pub renewed_fraction: f64,
}

impl IceShellProcess {
    pub fn new(graph: &ParamGraph, grid: CubeSphere) -> Option<Self> {
        if !crate::stages::ice::is_ice_surface(graph) {
            return None;
        }
        let body_radius = graph.f64("body.radius");
        let shell = ShellMap::new(graph);
        let n = grid.cell_count();
        let fraction = graph
            .f64("tectonics.ice_resurfaced_fraction")
            .clamp(0.0, 0.95);
        let mut lanes: Vec<Lane> = Vec::new();
        let mut lane_coverage = 0.0;
        // Where the shell cracks, the tidal field steers the lanes as it steers the
        // cracks; where it does not, nothing does, and a lane follows a great circle.
        let cracks = Cracks::new(graph);
        let thermal = crate::thermal::ThermalHistory::from_graph(graph);
        let areas: Vec<f64> = (0..n).map(|k| grid.cell_area(grid.cell(k))).collect();
        let total: f64 = areas.iter().sum();
        let mut covered = vec![false; n];
        if fraction > 0.0 {
            // The episode: its centre and its length.
            let age = graph.f64("body.age");
            let te = (age - graph.f64("tectonics.ice_resurfacing_age")).max(0.0);
            let duration = graph.f64("tectonics.ice_resurfacing_duration").max(0.0);
            // Lanes until they cover the share, their lengths and widths log-uniform, the
            // coverage counted by area on the coarse grid.
            let rng = Rng::new(graph.f64("seed.tectonics") as u64 ^ 0x6700_0E5A);
            let belts = GradientNoise::new(graph.f64("seed.tectonics") as u64 ^ 0x6700_BE17);
            let strikes = GradientNoise::new(graph.f64("seed.tectonics") as u64 ^ 0x6700_4EAD);
            let cell = grid.delta() * body_radius;
            let mut area = 0.0;
            // Widths from the shell that breaks: the elastic thickness at the episode.
            let elastic = thermal.at(te).elastic.max(1e3);
            let (w_lo, w_hi) = (
                GROOVE_LANE_WIDTH_ELASTIC.0 * elastic,
                GROOVE_LANE_WIDTH_ELASTIC.1 * elastic,
            );
            let (o_lo, o_hi) = GROOVE_NEIGHBOUR_OFFSET;
            let tau = 2.0 * std::f64::consts::PI;
            with_stamper(&grid, |stamper| {
                for k in 0..GROOVE_MAX_LANES {
                    if area >= fraction * total {
                        break;
                    }
                    let half_width = (0.5 * w_lo * (w_hi / w_lo).powf(rng.uniform(k, 100)))
                        .min(0.5 * GROOVE_LANE_MAX_WIDTH_SHARE * body_radius);
                    let spacing = half_width.min(GROOVE_LANE_POINT_SPACING_RAD * body_radius);
                    // A set grows beside one of the lanes laid lately, while the
                    // extension is there, not beside any lane ever laid, which would
                    // crowd the crowded regions further.
                    let parent = (!lanes.is_empty()
                        && rng.uniform(k, 105) < GROOVE_NEIGHBOUR_SHARE)
                        .then(|| {
                            let recent = lanes.len().min(GROOVE_SET_RECENT);
                            let back = (rng.uniform(k, 106) * recent as f64) as usize;
                            lanes.len() - 1 - back.min(recent - 1)
                        });
                    // A lane of a set belongs to its parent's generation and is no longer
                    // than it; a lane on its own takes its time from the episode.
                    let (time_gyr, half_length) = match parent {
                        Some(i) => (
                            (lanes[i].time_gyr
                                + (rng.uniform(k, 111) - 0.5) * GROOVE_SET_TIME_SPREAD * duration)
                                .clamp(0.0, age),
                            lanes[i].half_length
                                * (GROOVE_SET_LENGTH_SHARE
                                    + (1.0 - GROOVE_SET_LENGTH_SHARE) * rng.uniform(k, 101)),
                        ),
                        None => (
                            (te + (rng.uniform(k, 111) - 0.5) * duration).clamp(0.0, age),
                            (half_width
                                * GROOVE_LANE_ASPECT.0
                                * (GROOVE_LANE_ASPECT.1 / GROOVE_LANE_ASPECT.0)
                                    .powf(rng.uniform(k, 101)))
                            .min(0.5 * GROOVE_LANE_MAX_LENGTH_RAD * body_radius),
                        ),
                    };
                    // The start and the heading: beside an earlier lane, one to three
                    // widths to its side somewhere along it, its heading turned a little;
                    // or anywhere in the belts of extension, by rejection, along the
                    // strike of the field of extension at the lane's time.
                    let (start, heading) = match parent {
                        Some(i) => {
                            let p = &lanes[i];
                            let j = ((rng.uniform(k, 110) * p.points.len() as f64) as usize)
                                .min(p.points.len() - 1);
                            let (pj, tj) = (p.points[j], p.tangent_at(j));
                            let nj = normalize(cross(pj, tj));
                            let side = if rng.uniform(k, 108) < 0.5 { -1.0 } else { 1.0 };
                            let offset = side * 2.0 * p.half_width / body_radius
                                * (o_lo + (o_hi - o_lo) * rng.uniform(k, 109));
                            let start = normalize([
                                pj[0] + nj[0] * offset,
                                pj[1] + nj[1] * offset,
                                pj[2] + nj[2] * offset,
                            ]);
                            let turn = GROOVE_NEIGHBOUR_TURN_DEG.to_radians()
                                * (2.0 * rng.uniform(k, 107) - 1.0);
                            (start, rotate(tj, start, turn))
                        }
                        None => {
                            let axis = rng.unit_direction(k, 1);
                            let u = normalize(cross(
                                axis,
                                if axis[0].abs() < 0.9 {
                                    [1.0, 0.0, 0.0]
                                } else {
                                    [0.0, 1.0, 0.0]
                                },
                            ));
                            let v = cross(axis, u);
                            let mut center = u;
                            for attempt in 0..16u64 {
                                let ang = tau * rng.uniform(k, 3 + 2 * attempt);
                                let (sa, ca) = ang.sin_cos();
                                center = normalize([
                                    u[0] * ca + v[0] * sa,
                                    u[1] * ca + v[1] * sa,
                                    u[2] * ca + v[2] * sa,
                                ]);
                                let w = 0.5
                                    + 0.5
                                        * belts.fbm(
                                            center,
                                            GROOVE_BELT_FREQUENCY / tau,
                                            2,
                                            2.0,
                                            0.5,
                                        );
                                if rng.uniform(k, 4 + 2 * attempt) < w.powf(GROOVE_BELT_POWER) {
                                    break;
                                }
                            }
                            let phase = if duration > 0.0 {
                                ((time_gyr - (te - 0.5 * duration)) / duration).clamp(0.0, 1.0)
                            } else {
                                0.5
                            };
                            (center, extension_heading(&strikes, center, phase))
                        }
                    };
                    // The centre line: along the field at the lane's time where the
                    // shell cracks and the tension at the start opens it, else along the
                    // great circle of the heading.
                    let traced = cracks
                        .as_ref()
                        .and_then(|c| {
                            let axis = crate::stages::lineae::mean_axis_at(
                                time_gyr,
                                age,
                                c.nsr_period_gyr,
                            );
                            let (tension, phase) = c.shell.peak_tension(axis, start);
                            (tension >= ICE_CRACK_STRESS).then(|| {
                                c.trace(
                                    axis,
                                    start,
                                    phase,
                                    spacing,
                                    2.0 * half_length,
                                    k ^ 0x6700_7ACE,
                                )
                            })
                        })
                        .filter(|pts| pts.len() >= 3);
                    let points = traced.unwrap_or_else(|| {
                        let axis = normalize(cross(start, heading));
                        Lane::great_circle(axis, start, half_length, spacing, body_radius)
                    });
                    let lane = Lane::from_points(
                        points,
                        spacing,
                        half_width,
                        rng.uniform(k, 102) >= GROOVE_SMOOTH_LANE_SHARE,
                        GROOVE_CLEAN_RANGE.0
                            + (GROOVE_CLEAN_RANGE.1 - GROOVE_CLEAN_RANGE.0) * rng.uniform(k, 112),
                        time_gyr,
                        k,
                    );
                    for sm in lane.sample(&grid, body_radius, 2.0 * cell, stamper) {
                        if !covered[sm.idx] && lane.coverage_at(sm.across, sm.along, cell) >= 0.5 {
                            covered[sm.idx] = true;
                            area += areas[sm.idx];
                        }
                    }
                    lanes.push(lane);
                }
            });
            lane_coverage = area / total.max(1e-300);
            if std::env::var_os("EXOGEN_DEBUG").is_some() {
                eprintln!(
                    "grooved terrain: {} lanes cover {lane_coverage:.3} of the surface for a share of {fraction:.3}, over {duration:.2} Gyr about {te:.2} Gyr, elastic {elastic:.0} m",
                    lanes.len()
                );
            }
        }
        let mean_shell = graph.f64("interior.ice_shell_thickness");
        // The renewed terrain follows the fissures of the thinned shell, traced here at the
        // finishing grid's step so that the loop and the finish draw the same outline;
        // and the past episodes' terrains follow theirs, carried to where the body left
        // them (F53).
        let age = graph.f64("body.age");
        let fissures = if shell.has_active_terrain() {
            let seed = graph.f64("seed.cratering") as u64;
            let fine = CubeSphere::new(graph.f64("detail.finishing_resolution") as u32);
            cracks
                .as_ref()
                .map(|c| c.generate(age, &fine, seed))
                .unwrap_or_default()
        } else {
            Vec::new()
        };
        let active = if shell.has_active_terrain() {
            eroded(&grid, &polar_mask(&grid, body_radius, &shell, &fissures))
        } else {
            vec![false; n]
        };
        let episodes = crate::stages::ice::focusing_episodes(graph, &shell);
        let episode_maps: Vec<ShellMap> =
            episodes.iter().map(|e| shell.with_pole(e.center)).collect();
        let episode_masks: Vec<Vec<bool>> = episode_maps
            .iter()
            .map(|m| eroded(&grid, &polar_mask(&grid, body_radius, m, &fissures)))
            .collect();
        let present_start = if shell.has_active_terrain() {
            (age - FOCUSING_EPISODE_DURATION_GYR).max(0.0)
        } else {
            0.0
        };
        let count = episodes.len();
        let renewed_fraction = (0..n)
            .filter(|&k| covered[k] || active[k] || episode_masks.iter().any(|m| m[k]))
            .map(|k| areas[k])
            .sum::<f64>()
            / total.max(1e-300);
        Some(Self {
            grid,
            body_radius,
            shell,
            mean_shell,
            lanes,
            lane_coverage,
            equilibrium: None,
            polar_floor: None,
            lane_floor: None,
            active,
            depressed: false,
            thermal,
            gravity: graph.f64("body.gravity"),
            present_start,
            episodes,
            episode_maps,
            episode_masks,
            episode_floors: vec![None; count],
            episode_sunk: vec![false; count],
            episode_rebounded: vec![false; count],
            active_now: vec![false; n],
            renewed_fraction,
        })
    }

    /// What the relaxation step takes from the shell: the thinned shell's equilibrium
    /// surface and the active polar terrain, renewed at every step down to the renewal
    /// length. The loop relaxes with the mean shell (see `relaxation::Local`); the local
    /// thickness acts through the replayed craters' survival, the crack density and the
    /// renewal.
    pub fn relaxation(&mut self, t: f64) -> crate::stages::relaxation::Local<'_> {
        // The terrain renewed at this step: the present pole's through its episode, each
        // past episode's through its own.
        let present = self.shell.has_active_terrain() && t >= self.present_start;
        let mut any = false;
        for k in 0..self.active_now.len() {
            let mut on = present && self.active[k];
            for (e, mask) in self.episodes.iter().zip(&self.episode_masks) {
                if t >= e.start && t <= e.end + 1e-9 && mask[k] {
                    on = true;
                }
            }
            self.active_now[k] = on;
            any |= on;
        }
        crate::stages::relaxation::Local {
            equilibrium: self.equilibrium.as_deref(),
            // The thinned zone keeps nothing up to its own scale while it is hot.
            renew: any.then(|| {
                (
                    self.active_now.as_slice(),
                    self.smooth_length()
                        .max(POLAR_RENEWAL_LENGTH_SHARE * self.body_radius),
                )
            }),
        }
    }

    /// The renewed polar surface the relaxation step produced.
    pub fn set_floor(&mut self, floor: Vec<f32>) {
        self.polar_floor = Some(floor);
    }

    /// Whether a point lies in a lane.
    pub fn in_lanes(&self, d: [f64; 3]) -> bool {
        self.lanes.iter().any(|l| l.inside(d, self.body_radius))
    }

    fn smooth_length(&self) -> f64 {
        GROOVE_SMOOTH_LENGTH.min(0.05 * self.body_radius)
    }

    /// Before the step's events: the polar depression once, and the lanes at their time.
    pub fn before(&mut self, layers: &mut LayerSet, t0: f64, t1: f64) {
        let grid = self.grid;
        let n = grid.cell_count();
        let cell = grid.delta() * self.body_radius;
        // The thinned shell floats lower while the focusing lasts (F43, F53): the present
        // zone's depression from its episode's start (from the first step where the
        // thinning is steady), each past episode's through its own and back up when it
        // ends. The shell floats to its equilibrium at once, and the relaxation keeps
        // that equilibrium rather than taking the depression for a load.
        let mut changed = false;
        let mut shift = |layers: &mut LayerSet, map: &ShellMap, sign: f32| {
            let mean = self.mean_shell;
            let bedrock = layers.bedrock.as_mut_slice();
            for (k, h) in bedrock.iter_mut().enumerate() {
                let d = grid.center(grid.cell(k));
                *h += sign * map.depression(d, mean) as f32;
            }
            changed = true;
        };
        if !self.depressed && self.shell.pole.is_some() && t1 >= self.present_start {
            self.depressed = true;
            shift(layers, &self.shell, -1.0);
        }
        for k in 0..self.episodes.len() {
            if !self.episode_sunk[k] && t1 >= self.episodes[k].start {
                self.episode_sunk[k] = true;
                shift(layers, &self.episode_maps[k], -1.0);
            }
            if self.episode_sunk[k] && !self.episode_rebounded[k] && t0 >= self.episodes[k].end {
                self.episode_rebounded[k] = true;
                shift(layers, &self.episode_maps[k], 1.0);
            }
        }
        if changed {
            let mean = self.mean_shell;
            let mut eq = vec![0.0f32; n];
            let mut any = false;
            if self.depressed {
                for (k, e) in eq.iter_mut().enumerate() {
                    *e -= self.shell.depression(grid.center(grid.cell(k)), mean) as f32;
                }
                any = true;
            }
            for (i, map) in self.episode_maps.iter().enumerate() {
                if self.episode_sunk[i] && !self.episode_rebounded[i] {
                    for (k, e) in eq.iter_mut().enumerate() {
                        *e -= map.depression(grid.center(grid.cell(k)), mean) as f32;
                    }
                    any = true;
                }
            }
            self.equilibrium = any.then_some(eq);
        }
        // The lanes whose time falls in this step, youngest last: each takes the relief as
        // it is now, smoothed and dropped, and the cells a younger lane also covers end up
        // with the younger lane's time.
        let mut due: Vec<usize> = (0..self.lanes.len())
            .filter(|&i| {
                let t = self.lanes[i].time_gyr;
                t > t0 && t <= t1
            })
            .collect();
        if !due.is_empty() {
            due.sort_by(|&a, &b| self.lanes[a].time_gyr.total_cmp(&self.lanes[b].time_gyr));
            let smoothed: Vec<f32> =
                crate::stages::relaxation::low_pass(&layers.bedrock, cell, self.smooth_length());
            let body_radius = self.body_radius;
            // The cells each due lane holds, the youngest lane's time winning.
            let mut inside: Vec<Option<(f32, f32)>> = vec![None; n];
            with_stamper(&grid, |stamper| {
                for &i in &due {
                    let lane = &self.lanes[i];
                    for sm in lane.sample(&grid, body_radius, 2.0 * cell, stamper) {
                        if lane.coverage_at(sm.across, sm.along, cell) >= 0.5 {
                            inside[sm.idx] = Some((lane.time_gyr as f32, lane.clean as f32));
                        }
                    }
                }
            });
            let mask: Vec<bool> = inside.iter().map(|t| t.is_some()).collect();
            let core = eroded(&grid, &mask);
            // The lanes' drop, the extension's subsidence, spreads over the lithosphere's
            // flexural parameter at the time, so a lane narrower than it sinks little and
            // none ends on a step.
            let te = self.thermal.at(0.5 * (t0 + t1)).elastic;
            let alpha = flexural_parameter(te, self.gravity);
            let mask_field: Vec<f32> = core.iter().map(|&m| if m { 1.0 } else { 0.0 }).collect();
            let ramp =
                crate::stages::relaxation::low_pass(&field_from(grid, &mask_field), cell, alpha);
            let floor = self.lane_floor.get_or_insert_with(|| smoothed.clone());
            let bedrock = layers.bedrock.as_mut_slice();
            let ages = layers.surface_age.as_mut_slice();
            let units = layers.unit_id.as_mut_slice();
            let deep = layers.deep.as_mut_slice();
            let melt = layers.melt.as_mut_slice();
            let cover = layers.cover.as_mut_slice();
            let clean = layers.clean.as_mut_slice();
            for k in 0..n {
                let drop = (GROOVE_LANE_DROP * ramp[k] as f64) as f32;
                if core[k]
                    && let Some((t, cl)) = inside[k]
                {
                    bedrock[k] = smoothed[k] - drop;
                    floor[k] = bedrock[k];
                    ages[k] = t;
                    units[k] = UNIT_GROOVED;
                    deep[k] = 0.0;
                    melt[k] = 0.0;
                    cover[k] = 0.0;
                    clean[k] = cl;
                } else if drop > 0.0 {
                    bedrock[k] -= drop;
                }
            }
        }
    }

    /// After the step's events: the polar terrain renews itself, erasing what the step
    /// brought it.
    pub fn after(&mut self, layers: &mut LayerSet, t1: f64) {
        let Some(floor) = &self.polar_floor else {
            return;
        };
        let bedrock = layers.bedrock.as_mut_slice();
        let ages = layers.surface_age.as_mut_slice();
        let units = layers.unit_id.as_mut_slice();
        let deep = layers.deep.as_mut_slice();
        let melt = layers.melt.as_mut_slice();
        let cover = layers.cover.as_mut_slice();
        let clean = layers.clean.as_mut_slice();
        let now = (t1 - POLAR_RENEWAL_GYR).max(0.0) as f32;
        for (k, &active) in self.active_now.iter().enumerate() {
            if active {
                bedrock[k] = floor[k];
                ages[k] = now;
                units[k] = UNIT_POLAR_TERRAIN;
                deep[k] = 0.0;
                melt[k] = 0.0;
                cover[k] = 0.0;
                clean[k] = 1.0;
            }
        }
        // A past episode's floor as it ended, for the finish to lay back.
        for k in 0..self.episodes.len() {
            if self.episode_floors[k].is_none()
                && self.episode_sunk[k]
                && t1 >= self.episodes[k].end
            {
                self.episode_floors[k] = Some(floor.clone());
            }
        }
    }
}

/// The lanes and the polar terrain on the resolved grid, before the finish stamps
/// anything. The coarse loop drew them cell by cell; here their edges are redrawn at the
/// finishing resolution, so the crater population that follows reads the ages they leave
/// and the colour follows the true edge. Cells something younger has taken keep it.
#[allow(clippy::too_many_arguments)]
pub fn rasterize(
    layers: &mut LayerSet,
    body_radius: f64,
    cell: f64,
    coarse_cell: f64,
    lanes: &[Lane],
    lane_floor: Option<&[f32]>,
) {
    let grid = layers.grid;
    let n = grid.cell_count();
    if !lanes.is_empty() {
        // Each lane's coverage of the cells near it, out to a coarse cell beyond its
        // edge, where it covers anything or the coarse layers left a trace of it, with
        // the lane's time.
        let old = layers.clean.as_slice();
        let margin = coarse_cell.max(2.0 * cell);
        let hits: Vec<Vec<(u32, f32)>> = lanes
            .par_iter()
            .map(|lane| {
                with_stamper(&grid, |stamper| {
                    lane.sample(&grid, body_radius, margin, stamper)
                        .into_iter()
                        .filter_map(|sm| {
                            let c = lane.coverage_at(sm.across, sm.along, cell);
                            (c > 0.0 || old[sm.idx] > 0.0).then_some((sm.idx as u32, c as f32))
                        })
                        .collect()
                })
            })
            .collect();
        let mut coverage = vec![0.0f32; n];
        let mut touched = vec![false; n];
        let mut times = vec![0.0f32; n];
        let mut cleans = vec![1.0f32; n];
        for (lane, h) in lanes.iter().zip(&hits) {
            for &(k, c) in h {
                let k = k as usize;
                touched[k] = true;
                if c > coverage[k] {
                    coverage[k] = c;
                    cleans[k] = lane.clean as f32;
                }
                if c >= 0.5 {
                    times[k] = times[k].max(lane.time_gyr as f32);
                }
            }
        }
        redraw(
            layers,
            &coverage,
            &touched,
            UNIT_GROOVED,
            &times,
            &cleans,
            false,
        );
        if let Some(floor) = lane_floor {
            lay_floor(layers, &coverage, floor);
        }
    }
}

/// A renewed polar terrain on the resolved grid, the present one or a past episode's:
/// its outline redrawn at this resolution from its shell map and the fissures, its unit
/// and time set, and its floor laid back. Returns the signed distance to its edge, cells.
#[allow(clippy::too_many_arguments)]
pub fn rasterize_polar(
    layers: &mut LayerSet,
    body_radius: f64,
    cell: f64,
    shell: &ShellMap,
    cracks: &[Linea],
    time: f64,
    floor: Option<&[f32]>,
    raise_only: bool,
) -> Vec<f32> {
    let grid = layers.grid;
    let n = grid.cell_count();
    let mask = polar_mask(&grid, body_radius, shell, cracks);
    let limit = polar_edge_limit(cell);
    let dist = edge_distance(&grid, &mask, limit);
    // The margin ramp across the edge.
    let coverage: Vec<f32> = dist
        .iter()
        .map(|&d| margin_coverage(d, cell) as f32)
        .collect();
    let touched: Vec<bool> = dist.iter().map(|&d| d > -(limit as f32)).collect();
    let times = vec![time as f32; n];
    redraw(
        layers,
        &coverage,
        &touched,
        UNIT_POLAR_TERRAIN,
        &times,
        &vec![1.0f32; n],
        raise_only,
    );
    if let Some(floor) = floor {
        lay_floor(layers, &coverage, floor);
    }
    if std::env::var_os("EXOGEN_DEBUG").is_some() {
        let inside = coverage.iter().filter(|&&c| c >= 0.5).count();
        let partial = coverage.iter().filter(|&&c| c > 0.0 && c < 0.5).count();
        eprintln!(
            "polar terrain on the fine grid at {time:.3} Gyr: {inside} cells inside, {partial} on the edge, {} touched",
            touched.iter().filter(|&&t| t).count()
        );
    }
    dist
}

/// The renewed terrain's floor laid under it before anything is stamped, by the
/// coverage, so that a crater formed since sits on the renewed surface rather than on
/// the old relief the coarse loop left outside its cells.
fn lay_floor(layers: &mut LayerSet, coverage: &[f32], floor: &[f32]) {
    for ((h, &c), &f) in layers
        .bedrock
        .as_mut_slice()
        .iter_mut()
        .zip(coverage)
        .zip(floor)
    {
        if c > 0.0 {
            *h += (f - *h) * c;
        }
    }
}

/// One renewed terrain's surface fields at the resolved edge: the clean substrate by the
/// coverage, the unit and age where the terrain holds more than half the cell and nothing
/// younger has taken it.
fn redraw(
    layers: &mut LayerSet,
    coverage: &[f32],
    touched: &[bool],
    unit: u16,
    times: &[f32],
    cleans: &[f32],
    raise_only: bool,
) {
    let ages = layers.surface_age.as_mut_slice();
    let units = layers.unit_id.as_mut_slice();
    let deep = layers.deep.as_mut_slice();
    let melt = layers.melt.as_mut_slice();
    let cover = layers.cover.as_mut_slice();
    let clean = layers.clean.as_mut_slice();
    for (k, &within) in touched.iter().enumerate() {
        if !within {
            continue;
        }
        let c = coverage[k];
        let time = times[k] as f64;
        let cl = c * cleans[k];
        clean[k] = if raise_only { clean[k].max(cl) } else { cl };
        if c >= 0.5 {
            if ages[k] as f64 <= time + 1e-3 {
                ages[k] = time as f32;
                units[k] = unit;
                deep[k] = 0.0;
                melt[k] = 0.0;
                cover[k] = 0.0;
            }
        } else if units[k] == unit {
            // The coarse cell reached past the edge.
            units[k] = UNIT_PRIMORDIAL;
            clean[k] = c * cleans[k];
        }
    }
}

fn write(idx: usize, bedrock: f64) -> CellWrite {
    CellWrite {
        idx: idx as u32,
        bedrock: bedrock as f32,
        surface: false,
        age: 0.0,
        unit: 0,
        deep: 0.0,
        melt: 0.0,
        cover: 0.0,
    }
}

/// A lane on the resolved grid at its time: the lane's surface back to its floor, the
/// smoothed coarse relief, with parallel grooves along the lane a few kilometres apart,
/// fading at the lane's edges and turning slowly along their length in sets. Cells
/// something younger has taken are left alone.
#[allow(clippy::too_many_arguments)]
pub fn groove_compute(
    layers: &LayerSet,
    grid: &CubeSphere,
    body_radius: f64,
    floor: &[f32],
    lane: &Lane,
    time_gyr: f64,
    cell: f64,
    stamper: &mut Stamper,
) -> StampOutput {
    let bedrock = layers.bedrock.as_slice();
    let ages = layers.surface_age.as_slice();
    let units = layers.unit_id.as_slice();
    let noise = GradientNoise::new(lane.id ^ 0x6700_0000);
    let rough_noise = GradientNoise::new(lane.id ^ 0x6700_0C4A);
    let rough_gain = 2f64.powf(-ROUGHNESS_HURST);
    let set_frequency = body_radius / GROOVE_SET_LENGTH;
    let spacing = GROOVE_SPACING.max(GROOVE_MIN_CELLS * cell);
    let samples = lane.sample(grid, body_radius, 2.0 * cell, stamper);
    let mut cells: Vec<CellWrite> = samples
        .par_iter()
        .filter_map(|sm| {
            let k = sm.idx;
            let c = lane.coverage_at(sm.across, sm.along, cell);
            if c <= 0.0 {
                return None;
            }
            // A cell the lane holds that something younger has taken keeps that.
            let held = units[k] != UNIT_GROOVED || ages[k] as f64 > time_gyr + 1e-3;
            if c >= 0.5 && held {
                return None;
            }
            let d = grid.center(grid.cell(k));
            // The tectonic roughness of a swath deformed at every scale, under the
            // grooves.
            let rough = GROOVE_ROUGHNESS
                * rough_noise.fbm(d, body_radius / GROOVE_ROUGH_LENGTH, 4, 2.0, rough_gain)
                / 0.4;
            let target = if lane.grooved {
                // Grooves measured in from the nearer edge, so that they follow the
                // edge's wander and the outermost sits just inside it as the bounding
                // trough; their spacing and strength vary with the set field, they
                // wander so that they are sinuous, and they come in segments that die
                // out and relay, so that no groove runs the lane's length unbroken.
                let half_width = lane.half_width_at(sm.along).max(1.0);
                let inward = (half_width - sm.across.abs()).max(0.0);
                let sets = noise.fbm(d, set_frequency, 2, 2.0, 0.5);
                let local_spacing = spacing * (1.0 + GROOVE_SPACING_VARIATION * sets);
                let wander = GROOVE_WANDER
                    * noise.fbm(
                        [d[0] + 3.1, d[1] - 1.7, d[2] + 4.3],
                        body_radius / GROOVE_WANDER_LENGTH,
                        3,
                        2.0,
                        0.5,
                    );
                let segment = {
                    let x = (GROOVE_SEGMENT_LEVEL
                        + noise.fbm(
                            [d[0] - 7.7, d[1] + 2.9, d[2] - 1.1],
                            body_radius / GROOVE_SEGMENT_LENGTH,
                            2,
                            2.0,
                            0.5,
                        ))
                    .clamp(0.0, 1.0);
                    x * x * (3.0 - 2.0 * x)
                };
                let strength = GROOVE_SLACK + (1.0 - GROOVE_SLACK) * (sets.abs() / 0.3).min(1.0);
                let phase = 2.0 * std::f64::consts::PI * inward / local_spacing + wander;
                let grooves = GROOVE_AMPLITUDE * strength * segment * phase.sin();
                let trough = -GROOVE_MARGIN_TROUGH * (-(inward / (0.5 * spacing)).powi(2)).exp();
                floor[k] as f64 + grooves + trough + rough
            } else {
                floor[k] as f64 + 0.5 * rough
            };
            Some(write(
                k,
                bedrock[k] as f64 + (target - bedrock[k] as f64) * c,
            ))
        })
        .collect();
    cells.sort_by_key(|c| c.idx);
    StampOutput {
        cells,
        composition: None,
    }
}

/// The polar terrain on the resolved grid at its time: its surface back to its floor,
/// the smoothed coarse relief, with folds between the cracks parallel to the nearest
/// one, and the margin ridge where the spreading terrain meets the old; a past terrain's
/// structures spread and lowered by the gardening since (`kappa_t`, m²).
#[allow(clippy::too_many_arguments)]
pub fn polar_compute(
    layers: &LayerSet,
    grid: &CubeSphere,
    body_radius: f64,
    floor: &[f32],
    shell: &ShellMap,
    edge: &[f32],
    cracks: &[Linea],
    _fold_direction: [f64; 3],
    time_gyr: f64,
    kappa_t: f64,
    cell: f64,
    stamper: &mut Stamper,
) -> StampOutput {
    let Some(pole) = shell.pole else {
        return StampOutput::default();
    };
    let bedrock = layers.bedrock.as_slice();
    let ages = layers.surface_age.as_slice();
    let units = layers.unit_id.as_slice();
    let limit = polar_edge_limit(cell) as f32;
    let reach = shell.radius_rad;
    let spacing = POLAR_FOLD_SPACING.max(GROOVE_MIN_CELLS * cell);
    let fabric = FabricIndex::new(cracks, POLAR_FABRIC_RANGE / body_radius);
    let ridge_keep = degraded(POLAR_MARGIN_WIDTH, kappa_t);
    let belt_keep = degraded(0.5 * POLAR_BELT_SPACING.max(3.0 * cell), kappa_t);
    let fold_keep = degraded(0.5 * spacing, kappa_t);
    let noise = GradientNoise::new(0x5027_F01D);
    let belt_noise = GradientNoise::new(0x5027_BE17);
    let rough_noise = GradientNoise::new(0x5027_0C4A);
    let rough_gain = 2f64.powf(-ROUGHNESS_HURST);
    let set_frequency = body_radius / GROOVE_SET_LENGTH.min(0.2 * body_radius);
    stamper.gather(grid, body_radius, pole, reach * body_radius);
    let mut cells: Vec<CellWrite> = stamper
        .footprint
        .par_iter()
        .filter_map(|&(k, _, _)| {
            let d = grid.center(grid.cell(k));
            let dk = edge[k];
            if dk <= -limit {
                return None;
            }
            // The margin ridge straddles the edge, higher and lower along it as the
            // belt's segments are; inside, the folds fade in from it.
            let seg = belt_noise.fbm(d, body_radius / POLAR_BELT_SEGMENT, 2, 2.0, 0.5);
            let x = dk as f64 * cell / POLAR_MARGIN_WIDTH;
            let ridge = POLAR_MARGIN_HEIGHT * ridge_keep * (0.75 + 0.25 * seg) * (-x * x).exp();
            let c = margin_coverage(dk, cell);
            let held = units[k] != UNIT_POLAR_TERRAIN || ages[k] as f64 > time_gyr + 1e-3;
            if c <= 0.0 || (c >= 0.5 && held) {
                if ridge < 1.0 {
                    return None;
                }
                return Some(write(k, bedrock[k] as f64 + ridge));
            }
            // Inside: folds parallel to the nearest fissure, measured from the point on
            // it nearest this cell so that they curve with it, the field's direction
            // where none is near; towards the edge, the belt where the spreading is
            // taken up, folds parallel to the edge instead.
            let inward = dk as f64 * cell;
            let taper = (inward / POLAR_FOLD_TAPER).clamp(0.0, 1.0);
            let sets = noise.fbm(d, set_frequency, 2, 2.0, 0.5);
            // The folds parallel the nearest fissure: their phase is the distance to it,
            // continuous across the terrain and curving with it, and they fade with that
            // distance so that each fissure carries its own belt of ridges; beyond the
            // fabric's range the plain carries its roughness alone.
            let (across, near) = match fabric.soft_distance(
                d,
                POLAR_FABRIC_SOFTNESS / body_radius,
                0.5 * POLAR_FABRIC_RANGE / body_radius,
            ) {
                Some((soft, min)) => {
                    let dist = min * body_radius;
                    let x =
                        ((POLAR_FABRIC_RANGE - dist) / (0.5 * POLAR_FABRIC_RANGE)).clamp(0.0, 1.0);
                    (soft * body_radius, x * x * (3.0 - 2.0 * x))
                }
                None => (0.0, 0.0),
            };
            // The sets curve: the distance is warped by a smooth field.
            let warp = POLAR_FOLD_WARP
                * noise.fbm(
                    [d[0] - 5.5, d[1] + 2.2, d[2] - 8.8],
                    body_radius / POLAR_FOLD_WARP_LENGTH,
                    2,
                    2.0,
                    0.5,
                );
            let across = across + warp;
            // Fields of ridges: the spacing and the strength vary with the set field,
            // the ridges wander so that they are sinuous, and they come in segments that
            // die out and relay, so the fabric is not one grating.
            let local_spacing = spacing * (1.0 + POLAR_FOLD_SPACING_VARIATION * sets);
            let wander = POLAR_FOLD_WANDER
                * noise.fbm(d, body_radius / POLAR_FOLD_WANDER_LENGTH, 3, 2.0, 0.5);
            let segment = {
                let x = (POLAR_FOLD_SEGMENT_LEVEL
                    + noise.fbm(
                        [d[0] + 7.1, d[1] - 3.3, d[2] + 1.9],
                        body_radius / POLAR_FOLD_SEGMENT,
                        2,
                        2.0,
                        0.5,
                    ))
                .clamp(0.0, 1.0);
                x * x * (3.0 - 2.0 * x)
            };
            let phase = 2.0 * std::f64::consts::PI * across / local_spacing + wander;
            let strength = (sets.abs() / POLAR_FOLD_SET_FULL).min(1.0);
            let folds = POLAR_FOLD_AMPLITUDE
                * fold_keep
                * (POLAR_FOLD_SLACK + (1.0 - POLAR_FOLD_SLACK) * strength)
                * segment
                * taper
                * near
                * phase.sin();
            // The tectonic roughness of a plain deformed at every scale, under the folds.
            let rough = POLAR_TECTONIC_RMS
                * fold_keep
                * taper
                * rough_noise.fbm(d, body_radius / POLAR_TECTONIC_LENGTH, 4, 2.0, rough_gain)
                / 0.4;
            // The belt decays inward from the edge; its folds come in segments that
            // shift against their neighbours, en echelon, never as concentric rings.
            let belt_w = (-inward / POLAR_BELT_WIDTH).exp();
            let belt = POLAR_BELT_HEIGHT
                * belt_keep
                * belt_w
                * (0.6 + 0.4 * seg)
                * (2.0 * std::f64::consts::PI * inward / POLAR_BELT_SPACING.max(3.0 * cell)
                    + POLAR_BELT_OFFSET * seg)
                    .sin();
            let target = floor[k] as f64 + folds * (1.0 - belt_w) + rough + belt;
            Some(write(
                k,
                bedrock[k] as f64 + (target - bedrock[k] as f64) * c + ridge,
            ))
        })
        .collect();
    cells.sort_by_key(|c| c.idx);
    if std::env::var_os("EXOGEN_DEBUG").is_some() {
        let moved = cells
            .iter()
            .filter(|c| (c.bedrock - bedrock[c.idx as usize]).abs() > 1.0)
            .count();
        eprintln!(
            "polar terrain stamp: {} cells written, {moved} moved by over a metre",
            cells.len()
        );
    }
    StampOutput {
        cells,
        composition: None,
    }
}

/// The shares of the surface, by area, that are grooved terrain and polar terrain.
pub fn terrain_shares(layers: &LayerSet) -> (f64, f64) {
    let grid = layers.grid;
    let units = layers.unit_id.as_slice();
    let (mut grooved, mut polar, mut total) = (0.0, 0.0, 0.0);
    for (k, &u) in units.iter().enumerate() {
        let a = grid.cell_area(grid.cell(k));
        total += a;
        if u == UNIT_GROOVED {
            grooved += a;
        } else if u == UNIT_POLAR_TERRAIN {
            polar += a;
        }
    }
    (grooved / total.max(1e-300), polar / total.max(1e-300))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_lane_is_a_swath_along_part_of_a_great_circle() {
        let r = 2.6e6;
        let points = Lane::great_circle([0.0, 0.0, 1.0], [1.0, 0.0, 0.0], 1000e3, 50e3, r);
        let lane = Lane::from_points(points, 50e3, 100e3, true, 1.0, 0.0, 1);
        assert!((lane.half_length - 1000e3).abs() < 1.0);
        assert!(lane.inside([1.0, 0.0, 0.0], r));
        // Off the circle by atan(0.1): 259 km, outside the width however it wobbles.
        let off = normalize([1.0, 0.0, 0.1]);
        let (across, along, _) = lane.locate(off, r);
        assert!((across.abs() - 0.1f64.atan() * r).abs() < 1.0);
        assert!(along.abs() < 1.0);
        assert!(!lane.inside(off, r));
        // The width wanders by a quarter or so, never to nothing.
        let mid = lane.half_width_at(0.0);
        assert!((50e3..=160e3).contains(&mid), "{mid}");
        assert!(lane.half_width_at(150e3) <= 160e3 && lane.half_width_at(150e3) >= 40e3);
        // The lens: the width closes at the ends.
        assert!(lane.half_width_at(990e3) < 0.5 * lane.half_width_at(0.0));
        // Along the circle by 0.3 rad: 780 km, inside the length; by 0.5 rad, outside.
        let near = [0.3f64.cos(), 0.3f64.sin(), 0.0];
        let (_, along, tangent) = lane.locate(near, r);
        assert!((along - 0.3 * r).abs() < 1.0);
        assert!(dot(tangent, [-(0.3f64.sin()), 0.3f64.cos(), 0.0]) > 0.999);
        assert!(lane.inside(near, r));
        let far = [0.5f64.cos(), -(0.5f64.sin()), 0.0];
        assert!(!lane.inside(far, r));
    }
}
