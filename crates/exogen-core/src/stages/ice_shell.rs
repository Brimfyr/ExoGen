// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! Partial resurfacing of an ice shell (PLAN.md 4.3 and Phase 4): grooved terrain and a
//! focused polar terrain.
//!
//! A shell that has stretched renews part of its surface as lanes of new ice: swaths of
//! parallel grooves, bright and clean, cutting the old dark cratered terrain into
//! polygons, as two thirds of Ganymede were renewed. The lanes are great-circle swaths of
//! random length and width drawn until they cover the resurfaced share of the surface;
//! inside them the old relief is smoothed away at their time, the surface age resets, and
//! the finish cuts the grooves. A shell over an ocean whose tidal heating focuses under a
//! pole is thin there, floats lower there, and renews that polar terrain continuously, so
//! it carries no craters, most of the cracks, and folds between them, as Enceladus's
//! south polar terrain does.

use rayon::prelude::*;

use crate::constants::*;
use crate::crater::{CellWrite, StampOutput, Stamper, with_stamper};
use crate::fields::{Field, LayerSet};
use crate::geom::{angle_between, cross, dot, normalize};
use crate::grid::CubeSphere;
use crate::noise::GradientNoise;
use crate::params::ParamGraph;
use crate::rng::Rng;
use crate::stages::ice::ShellMap;
use crate::units::*;

pub const VERSION: u64 = 2;

/// A lane of grooved terrain: a swath along part of a great circle.
#[derive(Clone, Debug)]
pub struct Lane {
    /// Pole of the great circle the lane follows.
    pub axis: [f64; 3],
    /// The lane's midpoint, on the circle.
    pub center: [f64; 3],
    pub half_length: f64,
    pub half_width: f64,
    /// Whether the lane is grooved or smooth bright terrain.
    pub grooved: bool,
    /// Phases of the two waves that wobble the width along the lane.
    pub wobble_phase: [f64; 2],
    pub id: u64,
}

impl Lane {
    /// Distance from the lane's centre line, m.
    pub fn across(&self, d: [f64; 3], body_radius: f64) -> f64 {
        dot(self.axis, d).clamp(-1.0, 1.0).asin().abs() * body_radius
    }

    /// Distance along the centre line from the midpoint, m, signed.
    pub fn signed_along(&self, d: [f64; 3], body_radius: f64) -> f64 {
        let t = cross(self.axis, self.center);
        dot(d, t).atan2(dot(d, self.center)) * body_radius
    }

    /// Distance along the centre line from the midpoint, m.
    pub fn along(&self, d: [f64; 3], body_radius: f64) -> f64 {
        self.signed_along(d, body_radius).abs()
    }

    /// The lane's half width at a signed distance along it, m: the mean wobbled by two
    /// waves three and seven widths long.
    pub fn half_width_at(&self, along: f64) -> f64 {
        let w = 2.0 * self.half_width;
        let turn = 2.0 * std::f64::consts::PI * along;
        let wobble = 0.5
            * ((turn / (3.0 * w) + self.wobble_phase[0]).sin()
                + (turn / (7.0 * w) + self.wobble_phase[1]).sin());
        self.half_width * (1.0 + GROOVE_WIDTH_WOBBLE * wobble)
    }

    pub fn inside(&self, d: [f64; 3], body_radius: f64) -> bool {
        let along = self.signed_along(d, body_radius);
        along.abs() < self.half_length && self.across(d, body_radius) < self.half_width_at(along)
    }

    /// How much of a cell at a direction the lane covers, 0 to 1: a one-cell ramp across
    /// the lane's edges and ends, for edges drawn without stair steps.
    pub fn coverage(&self, d: [f64; 3], body_radius: f64, cell: f64) -> f64 {
        let along = self.signed_along(d, body_radius);
        let cell = cell.max(1.0);
        let side = (self.half_width_at(along) - self.across(d, body_radius)) / cell + 0.5;
        let end = (self.half_length - along.abs()) / cell + 0.5;
        side.clamp(0.0, 1.0) * end.clamp(0.0, 1.0)
    }

    /// Angular radius of a cap round the midpoint that holds the lane.
    pub fn cap_radius(&self, body_radius: f64) -> f64 {
        (self.half_length + self.half_width * (1.0 + GROOVE_WIDTH_WOBBLE)) / body_radius
    }
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

/// Passes of the four-neighbour smoothing that blur over a length on a grid.
fn smoothing_passes(length: f64, cell: f64) -> usize {
    ((4.0 * (length / cell.max(1.0)).powi(2)).ceil() as usize).clamp(1, 400)
}

/// The shell's partial resurfacing through the loop.
pub struct IceShellProcess {
    grid: CubeSphere,
    body_radius: f64,
    pub shell: ShellMap,
    mean_shell: f64,
    /// When the lanes of grooved terrain formed, if any.
    pub groove_time: Option<f64>,
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
    grooves_done: bool,
}

impl IceShellProcess {
    pub fn new(graph: &ParamGraph, grid: CubeSphere) -> Option<Self> {
        if !crate::stages::ice::is_ice_surface(graph) {
            return None;
        }
        let body_radius = graph.f64("body.radius");
        let shell = ShellMap::new(graph);
        let n = grid.cell_count();
        let centers: Vec<[f64; 3]> = (0..n).map(|k| grid.center(grid.cell(k))).collect();
        let fraction = graph
            .f64("tectonics.ice_resurfaced_fraction")
            .clamp(0.0, 0.95);
        let mut lanes: Vec<Lane> = Vec::new();
        let mut groove_time = None;
        let mut lane_coverage = 0.0;
        if fraction > 0.0 {
            // Lanes until they cover the share: parts of great circles at random, their
            // lengths and widths log-uniform, the coverage counted by area on the coarse
            // grid.
            let rng = Rng::new(graph.f64("seed.tectonics") as u64 ^ 0x6700_0E5A);
            let belts = GradientNoise::new(graph.f64("seed.tectonics") as u64 ^ 0x6700_BE17);
            let areas: Vec<f64> = (0..n).map(|k| grid.cell_area(grid.cell(k))).collect();
            let total: f64 = areas.iter().sum();
            let mut covered = vec![false; n];
            let mut area = 0.0;
            let (w_lo, w_hi) = GROOVE_LANE_WIDTH;
            let (l_lo, l_hi) = GROOVE_LANE_LENGTH_RAD;
            let (o_lo, o_hi) = GROOVE_NEIGHBOUR_OFFSET;
            let tau = 2.0 * std::f64::consts::PI;
            for k in 0..GROOVE_MAX_LANES {
                if area >= fraction * total {
                    break;
                }
                let parent = (!lanes.is_empty() && rng.uniform(k, 105) < GROOVE_NEIGHBOUR_SHARE)
                    .then(|| {
                        let i = (rng.uniform(k, 106) * lanes.len() as f64) as usize;
                        lanes[i.min(lanes.len() - 1)].clone()
                    });
                let (axis, center) = match parent {
                    // Beside an earlier lane: its direction turned a little, one to three
                    // widths to its side, somewhere along its length.
                    Some(p) => {
                        let turn = GROOVE_NEIGHBOUR_TURN_DEG.to_radians()
                            * (2.0 * rng.uniform(k, 107) - 1.0);
                        let tangent = rotate(cross(p.axis, p.center), p.center, turn);
                        let side = if rng.uniform(k, 108) < 0.5 { -1.0 } else { 1.0 };
                        let offset = side * 2.0 * p.half_width / body_radius
                            * (o_lo + (o_hi - o_lo) * rng.uniform(k, 109));
                        let shift =
                            (2.0 * rng.uniform(k, 110) - 1.0) * 0.5 * p.half_length / body_radius;
                        let c = normalize([
                            p.center[0] + offset * p.axis[0],
                            p.center[1] + offset * p.axis[1],
                            p.center[2] + offset * p.axis[2],
                        ]);
                        let t = normalize([
                            tangent[0] - c[0] * dot(tangent, c),
                            tangent[1] - c[1] * dot(tangent, c),
                            tangent[2] - c[2] * dot(tangent, c),
                        ]);
                        let (ss, cs) = shift.sin_cos();
                        let center = normalize([
                            c[0] * cs + t[0] * ss,
                            c[1] * cs + t[1] * ss,
                            c[2] * cs + t[2] * ss,
                        ]);
                        (normalize(cross(center, t)), center)
                    }
                    // Or anywhere, in the belts of extension, by rejection.
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
                            let w =
                                0.5 + 0.5 * belts.fbm(center, GROOVE_BELT_FREQUENCY, 2, 2.0, 0.5);
                            if rng.uniform(k, 4 + 2 * attempt) < w * w {
                                break;
                            }
                        }
                        (axis, center)
                    }
                };
                let lane = Lane {
                    axis,
                    center,
                    half_length: 0.5 * body_radius * l_lo * (l_hi / l_lo).powf(rng.uniform(k, 101)),
                    half_width: (0.5 * w_lo * (w_hi / w_lo).powf(rng.uniform(k, 100)))
                        .min(0.3 * body_radius),
                    grooved: rng.uniform(k, 102) >= GROOVE_SMOOTH_LANE_SHARE,
                    wobble_phase: [tau * rng.uniform(k, 103), tau * rng.uniform(k, 104)],
                    id: k,
                };
                for (i, c) in centers.iter().enumerate() {
                    if !covered[i] && lane.inside(*c, body_radius) {
                        covered[i] = true;
                        area += areas[i];
                    }
                }
                lanes.push(lane);
            }
            lane_coverage = area / total.max(1e-300);
            if std::env::var_os("EXOGEN_DEBUG").is_some() {
                eprintln!(
                    "grooved terrain: {} lanes cover {lane_coverage:.3} of the surface for a share of {fraction:.3}",
                    lanes.len()
                );
            }
            groove_time =
                Some((graph.f64("body.age") - graph.f64("tectonics.ice_resurfacing_age")).max(0.0));
        }
        let mean_shell = graph.f64("interior.ice_shell_thickness");
        let equilibrium = shell.pole.is_some().then(|| {
            centers
                .iter()
                .map(|&c| -shell.depression(c, mean_shell) as f32)
                .collect()
        });
        let active: Vec<bool> = centers.iter().map(|&c| shell.active(c)).collect();
        let active = eroded(&grid, &active);
        Some(Self {
            grid,
            body_radius,
            shell,
            mean_shell,
            groove_time,
            lanes,
            lane_coverage,
            equilibrium,
            polar_floor: None,
            lane_floor: None,
            active,
            depressed: false,
            grooves_done: false,
        })
    }

    /// What the relaxation step takes from the shell: the thinned shell's equilibrium
    /// surface and the active polar terrain, renewed at every step down to the renewal
    /// length. The loop relaxes with the mean shell (see `relaxation::Local`); the local
    /// thickness acts through the replayed craters' survival, the crack density and the
    /// renewal.
    pub fn relaxation(&self) -> crate::stages::relaxation::Local<'_> {
        crate::stages::relaxation::Local {
            equilibrium: self.equilibrium.as_deref(),
            renew: self
                .shell
                .has_active_terrain()
                .then(|| (self.active.as_slice(), self.smooth_length())),
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
        if !self.depressed && self.shell.pole.is_some() {
            // The thinner shell floats lower, once: the depression is its equilibrium,
            // which the relaxation keeps.
            self.depressed = true;
            let bedrock = layers.bedrock.as_mut_slice();
            for (k, h) in bedrock.iter_mut().enumerate() {
                let d = grid.center(grid.cell(k));
                *h -= self.shell.depression(d, self.mean_shell) as f32;
            }
        }
        if let Some(te) = self.groove_time
            && !self.grooves_done
            && te > t0
            && te <= t1
        {
            self.grooves_done = true;
            let passes = smoothing_passes(self.smooth_length(), cell);
            let smoothed: Vec<f32> =
                crate::stages::episodic::smooth(&grid, layers.bedrock.as_slice(), passes)
                    .into_iter()
                    .map(|h| h - GROOVE_LANE_DROP as f32)
                    .collect();
            let inside: Vec<bool> = (0..n)
                .into_par_iter()
                .map(|k| self.in_lanes(grid.center(grid.cell(k))))
                .collect();
            let inside = eroded(&grid, &inside);
            let bedrock = layers.bedrock.as_mut_slice();
            let ages = layers.surface_age.as_mut_slice();
            let units = layers.unit_id.as_mut_slice();
            let deep = layers.deep.as_mut_slice();
            let melt = layers.melt.as_mut_slice();
            let cover = layers.cover.as_mut_slice();
            let clean = layers.clean.as_mut_slice();
            for k in 0..n {
                if inside[k] {
                    bedrock[k] = smoothed[k];
                    ages[k] = te as f32;
                    units[k] = UNIT_GROOVED;
                    deep[k] = 0.0;
                    melt[k] = 0.0;
                    cover[k] = 0.0;
                    clean[k] = 1.0;
                }
            }
            self.lane_floor = Some(smoothed);
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
        for (k, &active) in self.active.iter().enumerate() {
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
    lanes: &[Lane],
    groove_time: Option<f64>,
    shell: &ShellMap,
    polar: bool,
    polar_time: f64,
) {
    let grid = layers.grid;
    let n = grid.cell_count();
    if let Some(te) = groove_time
        && !lanes.is_empty()
    {
        // Each lane's coverage of the cells in its cap, where it covers anything or the
        // coarse layers left a trace of it.
        let old = layers.clean.as_slice();
        let hits: Vec<Vec<(u32, f32)>> = lanes
            .par_iter()
            .map(|lane| {
                with_stamper(&grid, |stamper| {
                    stamper.gather(
                        &grid,
                        body_radius,
                        lane.center,
                        lane.cap_radius(body_radius) * body_radius,
                    );
                    stamper
                        .footprint
                        .iter()
                        .filter_map(|&(k, _, _)| {
                            let c = lane.coverage(grid.center(grid.cell(k)), body_radius, cell);
                            (c > 0.0 || old[k] > 0.0).then_some((k as u32, c as f32))
                        })
                        .collect()
                })
            })
            .collect();
        let mut coverage = vec![0.0f32; n];
        let mut touched = vec![false; n];
        for h in &hits {
            for &(k, c) in h {
                let k = k as usize;
                touched[k] = true;
                coverage[k] = coverage[k].max(c);
            }
        }
        redraw(layers, &coverage, &touched, UNIT_GROOVED, te, false);
    }
    if polar && let Some(pole) = shell.pole {
        let mut coverage = vec![0.0f32; n];
        let mut touched = vec![false; n];
        with_stamper(&grid, |stamper| {
            stamper.gather(&grid, body_radius, pole, shell.radius_rad * body_radius);
            for &(k, _, _) in &stamper.footprint {
                touched[k] = true;
                coverage[k] =
                    shell.active_coverage(grid.center(grid.cell(k)), body_radius, cell) as f32;
            }
        });
        redraw(
            layers,
            &coverage,
            &touched,
            UNIT_POLAR_TERRAIN,
            polar_time,
            true,
        );
        if std::env::var_os("EXOGEN_DEBUG").is_some() {
            let inside = coverage.iter().filter(|&&c| c >= 0.5).count();
            let partial = coverage.iter().filter(|&&c| c > 0.0 && c < 0.5).count();
            eprintln!(
                "polar terrain on the fine grid: {inside} cells inside, {partial} on the edge, {} touched",
                touched.iter().filter(|&&t| t).count()
            );
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
    time: f64,
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
        clean[k] = if raise_only { clean[k].max(c) } else { c };
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
            clean[k] = c;
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
    let set_frequency = 2.0 * std::f64::consts::PI * body_radius / GROOVE_SET_LENGTH;
    let spacing = GROOVE_SPACING.max(GROOVE_MIN_CELLS * cell);
    let along_axis = normalize(cross(lane.axis, lane.center));
    stamper.gather(
        grid,
        body_radius,
        lane.center,
        lane.cap_radius(body_radius) * body_radius,
    );
    let mut cells: Vec<CellWrite> = stamper
        .footprint
        .par_iter()
        .filter_map(|&(k, _, _)| {
            let d = grid.center(grid.cell(k));
            let c = lane.coverage(d, body_radius, cell);
            if c <= 0.0 {
                return None;
            }
            // A cell the lane holds that something younger has taken keeps that.
            let held = units[k] != UNIT_GROOVED || ages[k] as f64 > time_gyr + 1e-3;
            if c >= 0.5 && held {
                return None;
            }
            let target = if lane.grooved {
                let across = lane.across(d, body_radius);
                let half_width = lane.half_width_at(lane.signed_along(d, body_radius));
                let edge = ((half_width - across) / (0.15 * half_width)).clamp(0.0, 1.0);
                let along = dot(d, along_axis) * body_radius;
                let sets = noise.fbm(d, set_frequency, 2, 2.0, 0.5);
                let phase = 2.0 * std::f64::consts::PI * across / spacing
                    + GROOVE_SET_TURN * sets * along / spacing;
                let amplitude = GROOVE_AMPLITUDE * (0.6 + 0.4 * sets.abs()) * edge;
                floor[k] as f64 + amplitude * phase.sin()
            } else {
                floor[k] as f64
            };
            Some(write(
                k,
                bedrock[k] as f64 + (target - bedrock[k] as f64) * c,
            ))
        })
        .collect();
    cells.sort_by_key(|c| c.idx);
    StampOutput { cells }
}

/// The polar terrain on the resolved grid at the end of the sequence: its surface back
/// to its floor, the smoothed coarse relief, with folds between the cracks along their
/// direction, and the margin ridge where the spreading terrain meets the old.
#[allow(clippy::too_many_arguments)]
pub fn polar_compute(
    layers: &LayerSet,
    grid: &CubeSphere,
    body_radius: f64,
    floor: &[f32],
    shell: &ShellMap,
    fold_direction: [f64; 3],
    time_gyr: f64,
    cell: f64,
    stamper: &mut Stamper,
) -> StampOutput {
    let Some(pole) = shell.pole else {
        return StampOutput::default();
    };
    let bedrock = layers.bedrock.as_slice();
    let ages = layers.surface_age.as_slice();
    let units = layers.unit_id.as_slice();
    let reach = shell.radius_rad.min(
        shell.active_radius_rad() * (1.0 + POLAR_MARGIN_WOBBLE)
            + 4.0 * POLAR_MARGIN_WIDTH / body_radius,
    );
    let normal = normalize(cross(pole, fold_direction));
    let spacing = POLAR_FOLD_SPACING.max(GROOVE_MIN_CELLS * cell);
    let noise = GradientNoise::new(0x5027_F01D);
    let set_frequency =
        2.0 * std::f64::consts::PI * body_radius / GROOVE_SET_LENGTH.min(0.2 * body_radius);
    stamper.gather(grid, body_radius, pole, reach * body_radius);
    let mut cells: Vec<CellWrite> = stamper
        .footprint
        .par_iter()
        .filter_map(|&(k, _, _)| {
            let d = grid.center(grid.cell(k));
            let angle = angle_between(pole, d);
            if angle >= reach {
                return None;
            }
            let edge = shell.active_radius_at(d);
            let x = (angle - edge) * body_radius / POLAR_MARGIN_WIDTH;
            let ridge = POLAR_MARGIN_HEIGHT * (-x * x).exp();
            let c = shell.active_coverage(d, body_radius, cell);
            let held = units[k] != UNIT_POLAR_TERRAIN || ages[k] as f64 > time_gyr + 1e-3;
            if c <= 0.0 || (c >= 0.5 && held) {
                if ridge < 1.0 {
                    return None;
                }
                return Some(write(k, bedrock[k] as f64 + ridge));
            }
            let taper = ((edge - angle) / (0.15 * edge)).clamp(0.0, 1.0);
            let sets = noise.fbm(d, set_frequency, 2, 2.0, 0.5);
            let across = dot(d, normal) * body_radius;
            let along = dot(d, fold_direction) * body_radius;
            let phase = 2.0 * std::f64::consts::PI * across / spacing
                + GROOVE_SET_TURN * sets * along / spacing;
            let folds = POLAR_FOLD_AMPLITUDE * (0.6 + 0.4 * sets.abs()) * taper * phase.sin();
            let target = floor[k] as f64 + folds;
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
    StampOutput { cells }
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
        let lane = Lane {
            axis: [0.0, 0.0, 1.0],
            center: [1.0, 0.0, 0.0],
            half_length: 1000e3,
            half_width: 100e3,
            grooved: true,
            wobble_phase: [0.0, 0.0],
            id: 1,
        };
        assert!(lane.inside([1.0, 0.0, 0.0], r));
        // Off the circle by atan(0.1): 259 km, outside the width however it wobbles.
        let off = normalize([1.0, 0.0, 0.1]);
        assert!((lane.across(off, r) - 0.1f64.atan() * r).abs() < 1.0);
        assert!(!lane.inside(off, r));
        assert!((lane.half_width_at(0.0) - 100e3).abs() < 1.0);
        assert!(lane.half_width_at(150e3) <= 125e3 && lane.half_width_at(150e3) >= 75e3);
        // Along the circle by 0.3 rad: 780 km, inside the length; by 0.5 rad, outside.
        let near = [0.3f64.cos(), 0.3f64.sin(), 0.0];
        assert!((lane.along(near, r) - 0.3 * r).abs() < 1.0);
        assert!(lane.inside(near, r));
        let far = [0.5f64.cos(), -(0.5f64.sin()), 0.0];
        assert!(!lane.inside(far, r));
    }

    #[test]
    fn smoothing_passes_follow_the_length_over_the_cell() {
        assert_eq!(smoothing_passes(100e3, 50e3), 16);
        assert_eq!(smoothing_passes(10e3, 50e3), 1);
    }
}
