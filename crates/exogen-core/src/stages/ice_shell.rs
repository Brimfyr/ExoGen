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
//! where nothing steers it (Ganymede's), bounded by its two faults; lanes are drawn
//! until they cover the resurfaced share of the surface. Inside them the old relief is
//! smoothed away at their time and sinks by the extension's subsidence, spread over the
//! lithosphere's flexural parameter so that no lane ends on a step; the surface age
//! resets, and the finish cuts the grooves. A shell over an ocean whose tidal heating focuses under a
//! pole is thin there, floats lower there, and renews that polar terrain continuously, so
//! it carries no craters, most of the cracks, and folds between them, as Enceladus's
//! south polar terrain does.

use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashMap};

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
use crate::stages::lineae::{Cracks, Linea, LineaKind, crack_direction, mean_axis_at};
use crate::units::*;

pub const VERSION: u64 = 21;

/// A lane of renewed terrain (F44, F61): a swath about a centre line traced across the
/// shell, along the tidal field where the shell cracks, so that it curves as the field
/// turns (Enceladus's sulci), and along the strike of the extension where nothing else
/// steers it (Ganymede's). Its two edges are its bounding faults: each a trace of fault
/// segments that wanders as a fault trace does and steps where one segment relays to
/// the next. An end stops at a cut along an older lane's strike where the lane met one
/// that took up its extension; an end that has run its length dies out in the old
/// ground, its bounding faults converging and its grooves fading (F142).
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
    /// The bounding faults' traces: each edge's distance from the centre line as a factor
    /// on the half width, sampled from the first point to the last; the edge on the side
    /// of negative across, then the other.
    pub edges: [Vec<f32>; 2],
    /// Where the end at the first point, and the end at the last, stops; an end without
    /// a stop closes as a wedge over a few widths.
    pub stops: [Option<LaneStop>; 2],
    /// The spacing of the lane's grooves, m.
    pub groove_spacing: f64,
    /// The strike of each later pulse's grooves against the lane's own, rad (F63).
    pub turns: Vec<f64>,
    /// When the lane formed, Gyr since formation.
    pub time_gyr: f64,
    pub id: u64,
}

/// Where a lane's end stops inside an older lane (F61): the great circle it is cut off
/// along, as that circle's pole on the lane's own side, and how far from the end the cut
/// reaches, m.
#[derive(Clone, Copy, Debug)]
pub struct LaneStop {
    pub pole: [f64; 3],
    pub reach: f64,
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

/// What a lane makes of a cell (F61): its cover, where all the old relief gives way to
/// its floor; the relief it lays on that floor, m; the relief it adds to the ground
/// beyond its faults, m, the flank's uplift and the fractures of its damage zone; the
/// cleanliness of its ice where it covers; and of the ground its fractures break.
#[derive(Clone, Copy, Debug, Default)]
pub struct LaneFabric {
    pub cover: f64,
    pub relief: f64,
    pub outer: f64,
    pub tone: f64,
    pub fracture: f64,
}

/// A lane's fields of unresolved variation: where its groove sets change, how its
/// grooves wander and die out, and the roughness under them.
pub struct LaneNoise {
    grooves: GradientNoise,
    rough: GradientNoise,
    /// The draws that belong to one groove or one pit: its depth, its place, its size.
    draws: Rng,
}

/// One pulse's set of grooves at a cell: the phase across the set in turns, with the
/// groove the cell lies in and how far across it; the spacing there and the distance
/// along the set, m; how far down the groove's trough the cell lies, 0 to 1; how present
/// the groove's segment is there, 0 to 1; the groove's depth against the middling one;
/// and its relief there in units of the middling amplitude.
pub(crate) struct GrooveSet {
    pub(crate) turns: f64,
    pub(crate) stripe: f64,
    pub(crate) frac: f64,
    pub(crate) local: f64,
    pub(crate) along: f64,
    pub(crate) trough: f64,
    pub(crate) gate: f64,
    pub(crate) depth: f64,
    pub(crate) shape: f64,
    /// The level the set's troughs cut down from, in units of the middling amplitude (F132).
    pub(crate) plateau: f64,
}

/// A fault's depth against the middling one of its population, from a uniform draw
/// (F63): a power-law tail, so that half are shallower than the middling one by up to
/// three tenths and a few are much deeper, bounded softly at the deepest a lid holds.
fn fault_depth(u: f64) -> f64 {
    let raw = (u.max(1e-6) / 0.5).powf(-1.0 / GROOVE_DEPTH_TAIL);
    if raw <= 1.0 {
        raw
    } else {
        let room = GROOVE_DEPTH_MAX - 1.0;
        1.0 + room * (1.0 - (-(raw - 1.0) / room).exp())
    }
}

/// A groove's index among its set's, as a key for its own draws.
fn stripe_key(stripe: f64) -> u64 {
    (stripe as i64 + (1 << 30)) as u64
}

/// One pulse's set of grooves at a cell (F61, F63): a wave across the faults' zone,
/// struck `turn` rad from its own strike, whose spacing follows the set field and which
/// sways and wanders so that the grooves fan and bend; each groove lives in segments
/// along itself; and each has a depth of its own, drawn from the heavy tail faults have,
/// so that most are shallow and a few are master troughs several times as deep, the extra
/// depth going down and not up. `across` and `along` are m across the zone and along it;
/// where the faults run round a closed line `period` m long (a corona's belt, F133), the
/// slots along them are a whole number round it and close up without a seam. A `skew`
/// makes each trough a tilt block's: its wall towards rising `across` short and steep, the
/// scarp, the other long, the back slope (none for a lane's grooves).
#[allow(clippy::too_many_arguments)]
pub(crate) fn fault_set(
    across: f64,
    along: f64,
    d: [f64; 3],
    body_radius: f64,
    spacing: f64,
    sets: f64,
    level: f64,
    pulse: usize,
    turn: f64,
    grooves: &GradientNoise,
    draws: &Rng,
    period: Option<f64>,
    skew: f64,
) -> GrooveSet {
    let tau = 2.0 * std::f64::consts::PI;
    let g = grooves;
    let shift = 17.3 * pulse as f64;
    let (sin, cos) = turn.sin_cos();
    let (across, along) = (across * cos + along * sin, along * cos - across * sin);
    let local = spacing * (1.0 + GROOVE_SPACING_VARIATION * sets);
    let sway = GROOVE_SWAY
        * g.fbm(
            [d[0] + 9.2 + shift, d[1] + 0.4, d[2] - 6.1],
            body_radius / (GROOVE_SWAY_SPACINGS * spacing),
            2,
            2.0,
            0.5,
        );
    let wander = GROOVE_WANDER
        * g.fbm(
            [d[0] + 3.1, d[1] - 1.7 + shift, d[2] + 4.3],
            body_radius / (GROOVE_WANDER_SPACINGS * spacing),
            3,
            2.0,
            0.5,
        );
    let turns = across / local + (sway + wander) / tau;
    // The grooves' segments (F132): each slot along a groove holds a segment with
    // GROOVE_SEGMENT_PRESENCE, fewer the farther beyond the bounding fault, its line
    // stepped across the groove's and its throw falling to its tips. Each segment
    // cuts its own trough, a hollow one spacing across its line, into the plateau, so
    // that two beside each other in a relay both show, and the ground between them
    // and in their gaps stands at the plateau, with no seam where one gives way.
    let presence = GROOVE_SEGMENT_PRESENCE * {
        let x = (level / GROOVE_SEGMENT_LEVEL).clamp(0.0, 1.0);
        x * x * (3.0 - 2.0 * x)
    };
    let (slot, slots, along) = match period {
        Some(p) if p > 0.0 => {
            let n = (p / (GROOVE_SEGMENT_SPACINGS * spacing)).round().max(1.0);
            (p / n, Some(n as i64), along.rem_euclid(p))
        }
        _ => (GROOVE_SEGMENT_SPACINGS * spacing, None, along),
    };
    let (s0, j0) = (turns.floor() as i64, (along / slot).floor() as i64);
    let mut troughs = 0.0;
    // The segment that holds the cell most: its stripe, where across its trough the
    // cell lies, its taper, its depth and its hollow there.
    let mut best = (0.0, s0, 0.0, 0.0, 1.0, 0.0);
    for stripe in s0 - 1..=s0 + 1 {
        for j in j0 - 1..=j0 + 1 {
            // Round a closed line a slot's draws are its own whichever way it is reached.
            let jk = slots.map_or(j, |n| j.rem_euclid(n));
            let key = stripe_key(stripe as f64) ^ ((jk as u64) << 34) ^ ((pulse as u64) << 60);
            let draw = |k: u64| draws.uniform(key, 60 + k);
            if draw(0) >= presence {
                continue;
            }
            let centre = (j as f64 + draw(1)) * slot;
            let (short, long) = GROOVE_SEGMENT_LENGTH;
            let half = 0.5 * slot * (short + (long - short) * draw(2));
            let t = (along - centre) / half;
            if t.abs() >= 1.0 {
                continue;
            }
            let x = turns - (stripe as f64 + 0.75 + GROOVE_STEP * (2.0 * draw(3) - 1.0));
            if x.abs() >= 0.5 {
                continue;
            }
            let hollow = 0.5 + 0.5 * (tau * x + skew * (tau * x).sin()).cos();
            let taper = 1.0 - t * t;
            let depth = fault_depth(draw(4).max(1e-3));
            troughs += (1.0 + depth) * hollow * taper;
            if hollow * taper > best.0 {
                best = (hollow * taper, stripe, x, taper, depth, hollow);
            }
        }
    }
    let plateau = presence * (1.0 + GROOVE_DEPTH_MEAN) / 3.0;
    let (_, stripe, x, taper, depth, hollow) = best;
    GrooveSet {
        turns,
        stripe: stripe as f64,
        frac: (0.75 + x).rem_euclid(1.0),
        local,
        along,
        trough: hollow,
        gate: taper,
        depth,
        shape: plateau - troughs,
        plateau,
    }
}

impl LaneNoise {
    pub fn new(lane: &Lane) -> Self {
        Self {
            grooves: GradientNoise::new(lane.id ^ 0x6700_0000),
            rough: GradientNoise::new(lane.id ^ 0x6700_0C4A),
            draws: Rng::new(lane.id ^ 0x6700_D0A5),
        }
    }
}

/// How many samples a bounding fault's trace holds along a lane.
const LANE_EDGE_SAMPLES: usize = 256;

/// A bounding fault's trace along a lane, as a factor on the half width: the self-affine
/// wander of a fault trace, and the steps where one fault segment relays to the next,
/// each a jog of a groove spacing or so across a relay ramp several spacings long, so
/// that the trace runs straight and bends rather than notches, the offsets reverting to
/// the mean trace so that the lane keeps its width.
fn edge_trace(
    id: u64,
    side: usize,
    half_length: f64,
    half_width: f64,
    groove_spacing: f64,
) -> Vec<f32> {
    let salt = (side as u64 + 1) << 44;
    let noise = GradientNoise::new(id ^ 0x6700_D1F7 ^ salt);
    let rng = Rng::new(id ^ 0x6700_E0E5 ^ salt);
    let length = (2.0 * half_length).max(1.0);
    let half_width = half_width.max(1.0);
    let gain = 2f64.powf(-GROOVE_EDGE_HURST);
    let base = (GROOVE_EDGE_WANDER_WIDTHS * 2.0 * half_width)
        .min(length)
        .max(1.0);
    let raw: Vec<f64> = (0..=LANE_EDGE_SAMPLES)
        .map(|i| {
            let s = length * i as f64 / LANE_EDGE_SAMPLES as f64;
            noise.fbm([s / base, 0.3, 0.7], 1.0, 5, 2.0, gain)
        })
        .collect();
    let rms = (raw.iter().map(|v| v * v).sum::<f64>() / raw.len() as f64)
        .sqrt()
        .max(1e-9);
    let segment = (GROOVE_EDGE_SEGMENT_SPACINGS * groove_spacing).max(1.0);
    let count = (length / segment).ceil() as usize + 2;
    let (j_lo, j_hi) = GROOVE_EDGE_JOG;
    let mut offsets = Vec::with_capacity(count);
    let mut offset = 0.0f64;
    for i in 0..count as u64 {
        let sign = if rng.uniform(i, 0) < 0.5 { -1.0 } else { 1.0 };
        let jog = sign * groove_spacing * (j_lo + (j_hi - j_lo) * rng.uniform(i, 1));
        offset = GROOVE_EDGE_REVERT * offset + jog;
        offsets.push(offset);
    }
    let riser = (GROOVE_EDGE_RELAY_SPACINGS * groove_spacing / segment).clamp(0.05, 0.6);
    let limit = GROOVE_EDGE_STEP_LIMIT * half_width;
    raw.iter()
        .enumerate()
        .map(|(i, v)| {
            let s = length * i as f64 / LANE_EDGE_SAMPLES as f64;
            let x = s / segment;
            let j = (x.floor() as usize).min(count - 2);
            let t = ((x - j as f64 - (1.0 - riser)) / riser).clamp(0.0, 1.0);
            let t = t * t * (3.0 - 2.0 * t);
            let step = (offsets[j] * (1.0 - t) + offsets[j + 1] * t).clamp(-limit, limit);
            (1.0 + GROOVE_EDGE_WANDER * v / rms + step / half_width)
                .clamp(1.0 - GROOVE_EDGE_REACH, 1.0 + GROOVE_EDGE_REACH) as f32
        })
        .collect()
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

    /// A lane about a centre line with grooves of the spacing a lane of its width
    /// carries on a lid of the middling width in lid thicknesses, closing at both ends.
    pub fn from_points(
        points: Vec<[f64; 3]>,
        spacing: f64,
        half_width: f64,
        grooved: bool,
        clean: f64,
        time_gyr: f64,
        id: u64,
    ) -> Self {
        let (lo, hi) = GROOVE_LANE_WIDTH_ELASTIC;
        let lid = 2.0 * half_width / (lo * hi).sqrt();
        Self::structured(
            points,
            spacing,
            half_width,
            grooved,
            clean,
            time_gyr,
            id,
            GROOVE_NECKING_WAVELENGTH * GROOVE_RIFT_LID_SHARE * lid,
            [None, None],
        )
    }

    /// A lane about a centre line, with its grooves' spacing and where each end stops.
    #[allow(clippy::too_many_arguments)]
    pub fn structured(
        points: Vec<[f64; 3]>,
        spacing: f64,
        half_width: f64,
        grooved: bool,
        clean: f64,
        time_gyr: f64,
        id: u64,
        groove_spacing: f64,
        stops: [Option<LaneStop>; 2],
    ) -> Self {
        let half_length = 0.5 * points.len().saturating_sub(1) as f64 * spacing;
        let edges = [
            edge_trace(id, 0, half_length, half_width, groove_spacing),
            edge_trace(id, 1, half_length, half_width, groove_spacing),
        ];
        Self {
            points,
            spacing,
            half_length,
            half_width,
            grooved,
            clean,
            edges,
            stops,
            groove_spacing,
            turns: Vec::new(),
            time_gyr,
            id,
        }
    }

    /// The lane's midpoint.
    pub fn center(&self) -> [f64; 3] {
        self.points[self.points.len() / 2]
    }

    /// The centre line's tangent at one of its points.
    pub fn tangent_at(&self, i: usize) -> [f64; 3] {
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

    /// How far the lane's width has closed at a signed distance along it, 0 to 1: a
    /// wedge over a few widths at an end where the bounding faults die out, nothing at
    /// an end that stops at a cut.
    fn taper_at(&self, along: f64) -> f64 {
        if self.stops[(along >= 0.0) as usize].is_some() {
            return 1.0;
        }
        let to_end = (self.half_length - along.abs()).max(0.0);
        let reach = (GROOVE_END_TAPER_WIDTHS * 2.0 * self.half_width)
            .min(GROOVE_END_TAPER_SHARE * self.half_length)
            .max(1.0);
        (to_end / reach)
            .clamp(0.0, 1.0)
            .powf(GROOVE_LANE_TAPER_POWER)
    }

    /// One edge's distance from the centre line at a signed distance along the lane, m:
    /// the bounding fault's trace on that side (0 the side of negative across, 1 the
    /// other), closing toward an end that tapers.
    pub fn edge_at(&self, along: f64, side: usize) -> f64 {
        let s = (0.5 + 0.5 * along / self.half_length.max(1.0)).clamp(0.0, 1.0);
        let profile = &self.edges[side];
        let n = profile.len();
        let x = s * (n - 1) as f64;
        let i = (x.floor() as usize).min(n.saturating_sub(2));
        let f = (x - i as f64).clamp(0.0, 1.0);
        let trace = profile[i] as f64 * (1.0 - f) + profile[i + 1] as f64 * f;
        self.half_width * trace * self.taper_at(along)
    }

    /// The lane's half width at a signed distance along it, m: the mean of its two edges.
    pub fn half_width_at(&self, along: f64) -> f64 {
        0.5 * (self.edge_at(along, 0) + self.edge_at(along, 1))
    }

    /// The lid the lane broke, m: the layer whose necking spaced its grooves, over the
    /// share of the lid that necks in a rift.
    pub fn lid(&self) -> f64 {
        self.groove_spacing / (GROOVE_NECKING_WAVELENGTH * GROOVE_RIFT_LID_SHARE)
    }

    /// The width of the lane's border zone on a grid of a cell, m. A bounding fault cuts
    /// the whole lid, and the lane meets the old ground across about a lid's thickness:
    /// a physical width, the same fault on any grid; the grid's own ramp of a few cells
    /// is only its floor, where the grid cannot draw it narrower.
    pub fn border(&self, cell: f64) -> f64 {
        (GROOVE_BORDER_LIDS * self.lid()).max(GROOVE_EDGE_CELLS * cell.max(1.0))
    }

    /// How much of a cell the lane covers at an across and an along, 0 to 1: a ramp
    /// across its border zone at its bounding faults and its ends.
    pub fn coverage_at(&self, across: f64, along: f64, cell: f64) -> f64 {
        let ramp = self.border(cell);
        let edge = self.edge_at(along, (across >= 0.0) as usize);
        let side = (edge - across.abs()) / ramp + 0.5;
        let end = (self.half_length - along.abs()) / ramp + 0.5;
        side.clamp(0.0, 1.0) * end.clamp(0.0, 1.0)
    }

    /// Where a cell stands against the cut its end of the lane stops at: the signed
    /// distance to the cut, m, positive on the lane's side; nothing where the end has no
    /// cut or the cell is out of its reach.
    fn cut_at(&self, d: [f64; 3], along: f64, body_radius: f64) -> Option<f64> {
        let end = (along >= 0.0) as usize;
        let stop = self.stops[end]?;
        if self.half_length - along.abs() > stop.reach {
            return None;
        }
        Some(dot(d, stop.pole).clamp(-1.0, 1.0).asin() * body_radius)
    }

    /// The lane's cover of a cell, 0 to 1: the ramp across its border zone at its
    /// bounding faults, cut off where an end stops.
    pub fn cover_at(&self, sm: &LaneSample, d: [f64; 3], body_radius: f64, cell: f64) -> f64 {
        let stop = match self.cut_at(d, sm.along, body_radius) {
            Some(dist) => (dist / self.border(cell) + 0.5).clamp(0.0, 1.0),
            None => 1.0,
        };
        stop * self.coverage_at(sm.across, sm.along, cell)
    }

    pub fn inside(&self, d: [f64; 3], body_radius: f64) -> bool {
        let (across, along, _) = self.locate(d, body_radius);
        along.abs() < self.half_length
            && across.abs() < self.edge_at(along, (across >= 0.0) as usize)
    }

    /// The spacing of the lane's grooves on a grid of a cell, m: no finer than the grid
    /// draws.
    pub fn grooves_on(&self, cell: f64) -> f64 {
        self.groove_spacing.max(GROOVE_MIN_CELLS * cell)
    }

    /// The length the lane levels the relief it renews over, m (F141): its own width,
    /// which its stretching reworks, or the lanes' smoothing length if that is more.
    pub fn removal_length(&self, body_radius: f64) -> f64 {
        lane_smooth_length(body_radius).max(2.0 * self.half_width)
    }

    /// How far past its edges the lane reaches on a grid of a cell, m: the damage zone
    /// its fractures break or the outer half of its border zone, whichever is wider.
    pub fn reach(&self, cell: f64) -> f64 {
        (GROOVE_DAMAGE_SPACINGS * self.grooves_on(cell)).max(0.5 * self.border(cell)) + 2.0 * cell
    }

    /// Angular radius of a cap round the midpoint that holds the lane and its damage
    /// zone on a grid of a cell.
    pub fn cap_radius(&self, body_radius: f64, cell: f64) -> f64 {
        let c = self.center();
        let reach = self
            .points
            .iter()
            .map(|p| angle_between(c, *p))
            .fold(0.0, f64::max);
        reach
            + (self.half_width * (1.0 + GROOVE_EDGE_REACH) + self.reach(cell) + 2.0 * self.spacing)
                / body_radius
    }

    /// What the lane makes of a cell (F61). Inside its bounding faults, parallel grooves
    /// from the necking of its lid: a wave across the lane whose spacing and strength
    /// change from set to set, which sways and wanders at several scales so that the
    /// grooves fan and bend, and of which each groove lives in segments along itself, so
    /// that it dies out and its neighbours relay it; on them the fine set of tilted
    /// fault blocks, and under them the roughness of a swath deformed at every scale. At each edge the bounding fault: the lane's side down,
    /// the flank outside up, both decaying from the fault, with a throw that comes and
    /// goes along the trace as the fault's segments do and dies with the faults at a
    /// tapering end. Beyond the edge the damage zone: the same grooves as open fractures
    /// cut into the old ground, each in fewer and shorter segments the farther out, so
    /// that the lane's strain thins out through discrete breaks rather than fading; and
    /// where `loose` m of loose cover lie on that ground, the fractures drain it and
    /// show as chains of pits (F64). The tone follows the structure: sets differ,
    /// troughs gather the dark lag, and a fracture's walls show clean ice in the old
    /// ground.
    pub fn fabric_at(
        &self,
        sm: &LaneSample,
        d: [f64; 3],
        body_radius: f64,
        cell: f64,
        loose: f64,
        noise: &LaneNoise,
    ) -> LaneFabric {
        let tau = 2.0 * std::f64::consts::PI;
        let spacing = self.grooves_on(cell);
        let damage = GROOVE_DAMAGE_SPACINGS * spacing;
        let side = (sm.across >= 0.0) as usize;
        let edge = self.edge_at(sm.along, side);
        let past = sm.along.abs() - self.half_length;
        // The distance in from the bounding fault, negative beyond it and beyond an end.
        let inward = if past > 0.0 {
            -(sm.across.abs() - edge).max(0.0).hypot(past)
        } else {
            edge - sm.across.abs()
        };
        // The cut the end stops at, inside an older lane, leaves everything beyond it to
        // that lane.
        let (mut stop, mut beyond) = (1.0, 1.0);
        if let Some(dist) = self.cut_at(d, sm.along, body_radius) {
            stop = (dist / self.border(cell) + 0.5).clamp(0.0, 1.0);
            beyond = stop;
            if stop <= 0.0 {
                return LaneFabric::default();
            }
        }
        if inward <= -damage {
            return LaneFabric::default();
        }
        let cover = stop * self.coverage_at(sm.across, sm.along, cell);
        // The relief is that of the lane's own grooves, whatever spacing the grid
        // draws them at.
        let amplitude = GROOVE_RELIEF_SHARE * self.groove_spacing;
        let g = &noise.grooves;
        let set_frequency = body_radius / (GROOVE_SET_SPACINGS * spacing);
        let sets = g.fbm(d, set_frequency, 2, 2.0, 0.5);
        let tones = g.fbm(
            [d[0] - 2.3, d[1] + 5.9, d[2] + 0.7],
            set_frequency,
            2,
            2.0,
            0.5,
        );
        // Each groove's segments along itself; fewer and shorter the farther beyond the
        // fault.
        let out = (-inward / damage).clamp(0.0, 1.0);
        let level = GROOVE_SEGMENT_LEVEL - (GROOVE_SEGMENT_LEVEL + GROOVE_DAMAGE_LEVEL_DROP) * out;
        let strength = GROOVE_SLACK + (1.0 - GROOVE_SLACK) * (sets.abs() / 0.3).min(1.0);
        // The first pulse's set runs the whole lane. Each later pulse reworked part of
        // it (F63): in the patches it took, each of its grooves stands where its own
        // segment does, and across its own width cuts whatever stood there before, so
        // that the sets overlap groove by groove and never interfere as two waves would.
        let base = self.set_at(sm, d, body_radius, spacing, sets, level, 0, 0.0, noise);
        // The set's own segments and plateau (F132).
        let mut shape = base.shape;
        let mut lag = base.gate * base.trough * base.trough * base.depth.min(2.0);
        // The set that stands at the cell, and how fully.
        let mut top = (base.turns, 1.0);
        for (j, &turn) in self.turns.iter().enumerate() {
            let pulse = j + 1;
            let shift = 23.7 * pulse as f64;
            let patch = {
                let field = g.fbm(
                    [d[0] + shift, d[1] - 4.4, d[2] + 1.6],
                    body_radius / (GROOVE_PULSE_SPACINGS * spacing),
                    2,
                    2.0,
                    0.5,
                );
                let x = ((field - GROOVE_PULSE_LEVEL) / GROOVE_PULSE_EDGE + 0.5).clamp(0.0, 1.0);
                x * x * (3.0 - 2.0 * x)
            };
            if patch <= 0.0 {
                continue;
            }
            let set = self.set_at(sm, d, body_radius, spacing, sets, level, pulse, turn, noise);
            // The groove's own footprint: its segment along it, and its width across.
            // Across its own trough, where its segment stands (F132).
            let own = patch * set.gate * set.trough;
            shape += (set.shape - shape) * own;
            lag += (set.trough * set.trough * set.depth.min(2.0) - lag) * own;
            top = if own > 0.5 {
                (set.turns, 2.0 * own - 1.0)
            } else {
                (top.0, top.1 * (1.0 - 2.0 * own))
            };
        }
        // The grooves die out over the last stretch in from the bounding faults, which they
        // do not cross, leaving the plateau they cut (F132).
        let margin = {
            let x = (inward / (GROOVE_MARGIN_SPACINGS * spacing)).clamp(0.0, 1.0);
            x * x * (3.0 - 2.0 * x)
        };
        shape = base.plateau + margin * (shape - base.plateau);
        let rough_octaves = ((GROOVE_ROUGH_SPACINGS * spacing) / (2.0 * cell.max(1.0)))
            .log2()
            .ceil()
            .clamp(1.0, 10.0) as u32;
        let rough = GROOVE_ROUGH_SHARE
            * amplitude
            * noise.rough.fbm(
                d,
                body_radius / (GROOVE_ROUGH_SPACINGS * spacing),
                rough_octaves,
                2.0,
                2f64.powf(-ROUGHNESS_HURST),
            )
            / 0.4;
        // The fine set: the tilted fault blocks that ride on the necks (Pappalardo et
        // al. 1998; Patel et al. 1999), where the grid resolves them. They follow the
        // grooves of the set that stands at the cell with a wander of their own, face
        // the lane's axis from either side, a long back slope and a short scarp, and
        // each lives in its own segments.
        let fine = if self.groove_spacing / GROOVE_FINE_RATIO >= GROOVE_MIN_CELLS * cell {
            let turns = GROOVE_FINE_RATIO * top.0
                + GROOVE_WANDER
                    * g.fbm(
                        [d[0] - 5.3, d[1] + 8.1, d[2] - 0.9],
                        body_radius * GROOVE_FINE_RATIO / (GROOVE_WANDER_SPACINGS * spacing),
                        2,
                        2.0,
                        0.5,
                    )
                    / tau;
            let stripe = turns.floor();
            let phase = tau * (turns - stripe);
            let facing = (sm.across / spacing).clamp(-1.0, 1.0);
            let run = g.fbm(
                [
                    stripe * 7.2331 + 0.71,
                    sm.along * GROOVE_FINE_RATIO / (GROOVE_SEGMENT_SPACINGS * spacing),
                    0.47,
                ],
                1.0,
                2,
                2.0,
                0.5,
            );
            let x = (level + run).clamp(0.0, 1.0);
            GROOVE_FINE_RELIEF
                * amplitude
                * margin
                * top.1
                * x
                * x
                * (3.0 - 2.0 * x)
                * (phase + facing * GROOVE_FINE_SKEW * phase.sin()).sin()
        } else {
            0.0
        };
        // The bounding fault's throw here: its segments come and go along the trace, and
        // it dies with the faults at an end that tapers.
        let throw = {
            let run = noise.rough.fbm(
                [
                    side as f64 * 7.7 + 1.3,
                    sm.along / (GROOVE_EDGE_SEGMENT_SPACINGS * spacing),
                    0.59,
                ],
                1.0,
                2,
                2.0,
                0.5,
            );
            let x = (GROOVE_BOUNDARY_LEVEL + run).clamp(0.0, 1.0);
            let tip = if past > 0.0 {
                0.0
            } else {
                self.taper_at(sm.along)
            };
            GROOVE_BOUNDARY_THROW * amplitude * x * x * (3.0 - 2.0 * x) * tip
        };
        let down = -0.5 * throw * (-inward.max(0.0) / spacing).exp();
        let up = 0.5 * throw * (inward.min(0.0) / (GROOVE_FLANK_SPACINGS * spacing)).exp();
        // A smooth lane's grooves are subdued, and it breaks the ground beside it less.
        let (depth, breaking) = if self.grooved {
            (1.0, 1.0)
        } else {
            (GROOVE_SMOOTH_RELIEF, 0.5)
        };
        let grooves = depth * (strength * (amplitude * shape + fine) + rough);
        // Beyond the fault the first set's grooves go on as open fractures. Under loose
        // cover a fracture drains what lies on it and shows as a chain of pits (F64):
        // the cover fills the fracture by its depth over the fracture's, and where it
        // fills, the pits stand in its place, each as wide as the cover's slope of
        // repose makes it and as deep as the cover, where the grid can draw them.
        let fracture_depth = GROOVE_FRACTURE_SHARE * amplitude;
        let open = breaking * base.gate * base.trough.powi(3);
        let filled = (loose / fracture_depth.max(1.0)).clamp(0.0, 1.0);
        let pit_width = PIT_WIDTH_COVERS * loose.min(fracture_depth);
        let pits = if filled > 0.0 && pit_width >= PIT_MIN_CELLS * cell {
            let axis = (base.frac - 0.75) * base.local;
            breaking * self.pits_at(axis, base.along, base.stripe, pit_width, base.gate, noise)
        } else {
            0.0
        };
        let crack = open * (1.0 - filled) + pits * filled;
        let dug =
            fracture_depth * open * (1.0 - filled) + loose.min(fracture_depth) * pits * filled;
        let relief = grooves + down;
        let outer = beyond * (1.0 - out) * (up - dug);
        let set_tone = 1.0 - GROOVE_TONE_SETS * (0.5 + 0.5 * (tones / 0.5).clamp(-1.0, 1.0));
        let lag = depth * GROOVE_TONE_TROUGH * strength * lag.clamp(0.0, 1.0);
        LaneFabric {
            cover,
            relief,
            outer,
            tone: self.clean * set_tone * (1.0 - lag),
            fracture: beyond * GROOVE_DAMAGE_CLEAN * self.clean * (1.0 - out) * crack,
        }
    }

    /// One pulse's set of grooves at a cell (F61, F63): see `fault_set`.
    #[allow(clippy::too_many_arguments)]
    fn set_at(
        &self,
        sm: &LaneSample,
        d: [f64; 3],
        body_radius: f64,
        spacing: f64,
        sets: f64,
        level: f64,
        pulse: usize,
        turn: f64,
        noise: &LaneNoise,
    ) -> GrooveSet {
        fault_set(
            sm.across,
            sm.along,
            d,
            body_radius,
            spacing,
            sets,
            level,
            pulse,
            turn,
            &noise.grooves,
            &noise.draws,
            None,
            0.0,
        )
    }

    /// A chain of pits along a fracture (F64), 0 to 1 at a cell `axis` m off the
    /// fracture's line and `along` m along it: pits `width` m across, one every pit
    /// and a half or so, each a little off its place and of its own size, and standing
    /// only where the fracture's own segment does, so that the chain breaks as the
    /// fracture does.
    fn pits_at(
        &self,
        axis: f64,
        along: f64,
        stripe: f64,
        width: f64,
        gate: f64,
        noise: &LaneNoise,
    ) -> f64 {
        if axis.abs() > width {
            return 0.0;
        }
        let step = PIT_SPACING * width;
        let key = stripe_key(stripe);
        let place = along / step + noise.draws.uniform(key, 50);
        let nearest = place.floor() as i64;
        let mut best = 0.0f64;
        for m in nearest - 1..=nearest + 1 {
            let pit = key ^ ((m as u64) << 32);
            if noise.draws.uniform(pit, 51) > gate {
                continue;
            }
            let centre = (m as f64 + 0.5 + PIT_JITTER * (noise.draws.uniform(pit, 52) - 0.5)
                - noise.draws.uniform(key, 50))
                * step;
            let radius = 0.5 * width * (0.7 + 0.6 * noise.draws.uniform(pit, 53));
            let rho = axis.hypot(along - centre) / radius;
            if rho < 1.0 {
                best = best.max((1.0 - rho * rho).powf(1.5));
            }
        }
        best
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
        let reach = (self.half_width * (1.0 + GROOVE_EDGE_REACH) + margin).hypot(half + margin);
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

/// The fissures' points bucketed on the sphere, for the folds to follow the nearest one,
/// and the same points thinned along each fissure for the fissures' direction about a
/// place, each standing for the points about it.
struct FabricIndex {
    scale: f64,
    buckets: HashMap<(i32, i32, i32), Vec<FabricPoint>>,
    thinned: HashMap<(i32, i32, i32), Vec<(FabricPoint, f64)>>,
}

/// A fissure point and the tangent there.
type FabricPoint = ([f64; 3], [f64; 3]);

impl FabricIndex {
    /// The points of the direction's average stand a tenth of its smoothing length apart.
    const THINNING: f64 = 0.1;

    fn key(&self, d: [f64; 3]) -> (i32, i32, i32) {
        (
            (d[0] * self.scale).floor() as i32,
            (d[1] * self.scale).floor() as i32,
            (d[2] * self.scale).floor() as i32,
        )
    }

    fn new<'a>(cracks: impl IntoIterator<Item = &'a Linea>, range_rad: f64) -> Self {
        let mut index = Self {
            scale: 1.0 / range_rad.max(1e-6),
            buckets: HashMap::new(),
            thinned: HashMap::new(),
        };
        let gap = Self::THINNING * POLAR_FOLD_ORIENTATION_SHARE * range_rad;
        for l in cracks {
            let n = l.points.len();
            let (mut run, mut count) = (f64::INFINITY, 0.0);
            let mut last: Option<FabricPoint> = None;
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
                index.buckets.entry(index.key(p)).or_default().push((p, t));
                if i > 0 {
                    run += angle_between(l.points[i - 1], p);
                }
                count += 1.0;
                if run >= gap {
                    if let Some(q) = last {
                        index
                            .thinned
                            .entry(index.key(q.0))
                            .or_default()
                            .push((q, count - 1.0));
                        count = 1.0;
                    }
                    last = Some((p, t));
                    run = 0.0;
                }
            }
            if let Some(q) = last {
                index
                    .thinned
                    .entry(index.key(q.0))
                    .or_default()
                    .push((q, count));
            }
        }
        index
    }

    /// The fissures' direction about a direction: their tangents within a bucket's reach,
    /// carried to the point's tangent plane and averaged as doubled angles with a
    /// Gaussian weight of `sigma` rad, so that the field turns smoothly between fissures.
    fn orientation(&self, d: [f64; 3], sigma: f64) -> Option<[f64; 3]> {
        let (i, j, k) = self.key(d);
        let (east, north) = crate::geom::local_basis(d);
        let inv = 1.0 / (sigma.max(1e-9) * sigma.max(1e-9));
        let (mut c2, mut s2) = (0.0f64, 0.0f64);
        for di in -1..=1 {
            for dj in -1..=1 {
                for dk in -1..=1 {
                    let Some(v) = self.thinned.get(&(i + di, j + dj, k + dk)) else {
                        continue;
                    };
                    for &((p, t), count) in v {
                        // The chord stands for the arc: the two differ by a part in a few
                        // hundred at the reach of the weight.
                        let chord2 = 2.0 - 2.0 * dot(p, d);
                        if chord2 * inv > 18.42 {
                            continue;
                        }
                        let w = count * (-0.5 * chord2 * inv).exp();
                        let (a, b) = (dot(t, east), dot(t, north));
                        let r2 = a * a + b * b;
                        if r2 <= 0.0 {
                            continue;
                        }
                        c2 += w * (a * a - b * b) / r2;
                        s2 += w * 2.0 * a * b / r2;
                    }
                }
            }
        }
        if c2.hypot(s2) < 1e-9 {
            return None;
        }
        let theta = 0.5 * s2.atan2(c2);
        Some(normalize([
            east[0] * theta.cos() + north[0] * theta.sin(),
            east[1] * theta.cos() + north[1] * theta.sin(),
            east[2] * theta.cos() + north[2] * theta.sin(),
        ]))
    }

    /// The nearest fissure point within a bucket's reach of a direction, and the tangent
    /// there.
    fn nearest(&self, d: [f64; 3]) -> Option<([f64; 3], [f64; 3])> {
        let (i, j, k) = self.key(d);
        let mut best: Option<(f64, [f64; 3], [f64; 3])> = None;
        for di in -1..=1 {
            for dj in -1..=1 {
                for dk in -1..=1 {
                    let Some(v) = self.buckets.get(&(i + di, j + dj, k + dk)) else {
                        continue;
                    };
                    for &(p, t) in v {
                        // Nearest by the cosine, which falls as the angle grows.
                        let far = -dot(p, d);
                        if best.is_none_or(|b| far < b.0) {
                            best = Some((far, p, t));
                        }
                    }
                }
            }
        }
        best.map(|b| (b.1, b.2))
    }
}

/// One oriented segment of a renewed plain's fabric (F43): a wave windowed about its
/// centre, laid along the fissure nearest the centre and in phase with that fissure's
/// line.
struct FoldKernel {
    center: [f64; 3],
    tangent: [f64; 3],
    normal: [f64; 3],
    wavenumber: f64,
    phase: f64,
    strength: f64,
}

/// The fabric of a renewed plain as a sum of oriented segments: one kernel per cell of a
/// lattice at the segment scale over the terrain, each a wave windowed over half the
/// segment length along the ridges and a spacing or so across, laid along the fissure
/// nearest its centre and in phase with that fissure's line, so that the ridges parallel
/// each fissure, run on past its ends and relay from one fissure's set to the next
/// without a seam or a loop; the sum is normalized so that the amplitude holds where
/// kernels overlap and fades where they thin.
struct FoldKernels {
    scale: f64,
    body_radius: f64,
    sigma_along: f64,
    sigma_across: f64,
    kernels: HashMap<(i32, i32, i32), FoldKernel>,
}

impl FoldKernels {
    fn key(d: [f64; 3], scale: f64) -> (i32, i32, i32) {
        (
            (d[0] * scale).floor() as i32,
            (d[1] * scale).floor() as i32,
            (d[2] * scale).floor() as i32,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn new(
        fabric: &FabricIndex,
        shell: &ShellMap,
        footprint: &[(usize, f64, f64)],
        grid: &CubeSphere,
        body_radius: f64,
        spacing: f64,
        sets: &GradientNoise,
        set_frequency: f64,
        seed: u64,
    ) -> Self {
        let lattice = POLAR_FOLD_KERNEL_LATTICE * POLAR_FOLD_SEGMENT;
        // A fold's wavelength goes as the lid it buckles, and the lid thins toward the
        // terrain's centre: the spacing is the terrain's at the edge of its core, finer
        // inside and coarser in the tongues and the halo, never below what the grid
        // draws.
        let edge_factor = (1.0 - shell.thinning * POLAR_CORE_PROFILE).max(0.05);
        let finest = GROOVE_MIN_CELLS * grid.delta() * body_radius;
        let scale = body_radius / lattice;
        let reach = POLAR_FOLD_KERNEL_REACH;
        let rng = Rng::new(seed ^ 0x5027_4E12);
        let mut cells: std::collections::HashSet<(i32, i32, i32)> =
            std::collections::HashSet::new();
        for &(k, _, _) in footprint {
            cells.insert(Self::key(grid.center(grid.cell(k)), scale));
        }
        let mut keys: std::collections::HashSet<(i32, i32, i32)> = std::collections::HashSet::new();
        for (i, j, l) in cells {
            for di in -reach..=reach {
                for dj in -reach..=reach {
                    for dl in -reach..=reach {
                        keys.insert((i + di, j + dj, l + dl));
                    }
                }
            }
        }
        let kernels = keys
            .into_iter()
            .filter_map(|(i, j, l)| {
                let id = ((i as i64 & 0xF_FFFF) as u64)
                    | (((j as i64 & 0xF_FFFF) as u64) << 20)
                    | (((l as i64 & 0xF_FFFF) as u64) << 40);
                let center = normalize([
                    (i as f64 + rng.uniform(id, 0)) / scale,
                    (j as f64 + rng.uniform(id, 1)) / scale,
                    (l as f64 + rng.uniform(id, 2)) / scale,
                ]);
                let (p, t) = fabric.nearest(center)?;
                let dist = angle_between(p, center) * body_radius;
                let near = {
                    let x =
                        ((POLAR_FABRIC_RANGE - dist) / (0.5 * POLAR_FABRIC_RANGE)).clamp(0.0, 1.0);
                    x * x * (3.0 - 2.0 * x)
                };
                if near <= 0.0 {
                    return None;
                }
                // The segment runs along the fissures' direction about its centre, the
                // nearest fissure's tangent where no field forms.
                let field = fabric
                    .orientation(
                        center,
                        POLAR_FOLD_ORIENTATION_SHARE * POLAR_FABRIC_RANGE / body_radius,
                    )
                    .unwrap_or(t);
                let radial = dot(field, center);
                let tangent = normalize([
                    field[0] - radial * center[0],
                    field[1] - radial * center[1],
                    field[2] - radial * center[2],
                ]);
                // An earlier pulse (F63): where its patches stand, a kernel may keep the
                // strike of the field as it stood then, and lies in phase with its own
                // neighbours, not with the fissure's line.
                let patch = {
                    let x = (POLAR_PULSE_LEVEL
                        + sets.fbm(
                            [center[0] - 3.9, center[1] + 6.2, center[2] + 2.5],
                            scale / POLAR_PULSE_PATCH,
                            2,
                            2.0,
                            0.5,
                        ))
                    .clamp(0.0, 1.0);
                    x * x * (3.0 - 2.0 * x)
                };
                let earlier = rng.uniform(id, 5) < patch;
                let tangent = if earlier {
                    let turn = POLAR_PULSE_TURN_DEG.to_radians()
                        * (sets.fbm(
                            [center[0] + 8.8, center[1] - 0.6, center[2] - 4.1],
                            scale / POLAR_PULSE_PATCH,
                            2,
                            2.0,
                            0.5,
                        ) / 0.5)
                            .clamp(-1.0, 1.0);
                    rotate(tangent, center, turn)
                } else {
                    tangent
                };
                let normal = normalize(cross(center, tangent));
                let n_p = if earlier {
                    normal
                } else {
                    normalize(cross(p, t))
                };
                let s = sets.fbm(center, set_frequency, 2, 2.0, 0.5);
                let local_spacing = (spacing
                    * (shell.factor(center) / edge_factor)
                    * (1.0 + POLAR_FOLD_SPACING_VARIATION * s)
                    * 2f64.powf(POLAR_FOLD_SPACING_JITTER * (2.0 * rng.uniform(id, 4) - 1.0)))
                .max(finest);
                let wavenumber = 2.0 * std::f64::consts::PI / local_spacing;
                // In phase with the fissure's line: the centre's distance across it.
                let phase = wavenumber
                    * dot([center[0] - p[0], center[1] - p[1], center[2] - p[2]], n_p)
                    * body_radius;
                // A share of the segments are absent, so the fabric dies out and relays.
                let present = {
                    let x = (POLAR_FOLD_SEGMENT_LEVEL
                        + sets.fbm(
                            [center[0] + 7.1, center[1] - 3.3, center[2] + 1.9],
                            body_radius / POLAR_FOLD_SEGMENT,
                            2,
                            2.0,
                            0.5,
                        ))
                    .clamp(0.0, 1.0);
                    x * x * (3.0 - 2.0 * x)
                };
                // Each packet of folds has a depth of its own, from the heavy tail fault
                // throws have (F63), narrowed as the mean of a packet's few folds is, and
                // capped where the tail turns deep: the deep folds are the train's (F65).
                let strength = (POLAR_FOLD_SLACK
                    + (1.0 - POLAR_FOLD_SLACK) * (s.abs() / POLAR_FOLD_SET_FULL).min(1.0))
                    * near
                    * present
                    * fault_depth(rng.uniform(id, 3).max(POLAR_FOLD_TAIL_CAP))
                        .powf(POLAR_FOLD_PACKET_SPREAD)
                    * if earlier { POLAR_PULSE_RELIEF } else { 1.0 };
                Some((
                    (i, j, l),
                    FoldKernel {
                        center,
                        tangent,
                        normal,
                        wavenumber,
                        phase,
                        strength,
                    },
                ))
            })
            .collect();
        Self {
            scale,
            body_radius,
            sigma_along: 0.5 * POLAR_FOLD_SEGMENT,
            sigma_across: POLAR_FOLD_KERNEL_WIDTH * spacing,
            kernels,
        }
    }

    /// The fabric at a direction, in units of the fold amplitude.
    fn fabric(&self, d: [f64; 3]) -> f64 {
        let (i, j, l) = Self::key(d, self.scale);
        let reach = POLAR_FOLD_KERNEL_REACH;
        let (mut num, mut den) = (0.0f64, 0.0f64);
        for di in -reach..=reach {
            for dj in -reach..=reach {
                for dl in -reach..=reach {
                    let Some(k) = self.kernels.get(&(i + di, j + dj, l + dl)) else {
                        continue;
                    };
                    let off = [d[0] - k.center[0], d[1] - k.center[1], d[2] - k.center[2]];
                    let across = dot(off, k.normal) * self.body_radius;
                    let along = dot(off, k.tangent) * self.body_radius;
                    let w = (-0.5
                        * ((across / self.sigma_across).powi(2)
                            + (along / self.sigma_along).powi(2)))
                    .exp();
                    if w < 1e-4 {
                        continue;
                    }
                    num += w * k.strength * (k.wavenumber * across + k.phase).sin();
                    den += w * w;
                }
            }
        }
        if den <= 0.0 {
            return 0.0;
        }
        (num / den.sqrt().max(POLAR_FOLD_KERNEL_FLOOR)).clamp(-GROOVE_DEPTH_MAX, GROOVE_DEPTH_MAX)
    }
}

/// A hinge of the first-order train of deep folds (F65): its lattice cell, for a fixed
/// order; where it stands; and its depth against the middling fold (F63).
#[derive(Clone, Copy, Debug)]
struct DeepNucleus {
    key: (i32, i32, i32),
    center: [f64; 3],
    depth: f64,
}

/// What a growing hinge reads where it stands (F65): the folds' direction, where a
/// fabric stands; the fabric's strength, 0 to 1; the deep folds' wavelength, m, which
/// goes as the lid that folds; and the distance to the nearest of the terrain's own
/// vents, m.
#[derive(Clone, Copy, Debug)]
struct DeepSite {
    direction: Option<[f64; 3]>,
    strength: f64,
    wavelength: f64,
    fissure: f64,
}

/// How an arm of a deep fold ended: still growing; free, out of time or of fabric or at a
/// vent, where it plunges; or linked into a neighbour.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum ArmEnd {
    Growing,
    Free,
    Linked,
}

/// One arm of a deep fold: its hinge as it grew one way from the nucleus, the time each
/// point was reached (in wavelengths per unit of depth), the wavelength there, and, once
/// grown, the fold's height there in units of the fold amplitude.
#[derive(Clone, Debug)]
struct DeepArm {
    line: usize,
    depth: f64,
    points: Vec<[f64; 3]>,
    times: Vec<f64>,
    wavelengths: Vec<f64>,
    heights: Vec<f64>,
    heading: [f64; 3],
    /// The folds' direction where the arm's tip stands.
    direction: Option<[f64; 3]>,
    end: ArmEnd,
}

/// A piece of a deep fold's hinge as the relief reads it: the fold it belongs to, its
/// ends, the fold's height at each in units of the fold amplitude, and its wavelength at
/// each, m.
#[derive(Clone, Copy, Debug)]
struct DeepSegment {
    line: usize,
    a: [f64; 3],
    b: [f64; 3],
    height: [f64; 2],
    wavelength: [f64; 2],
}

/// The hinges' points bucketed on the sphere: each entry an arm and a point's index.
type HingePoints = HashMap<(i32, i32, i32), Vec<(u32, u32)>>;

/// The deep folds of a renewed plain (F65): a first-order train whose hinges nucleate one
/// to a wavelength, each grown both ways along the folds' direction from where it
/// nucleated, at a rate that goes with its depth, until it links with a neighbour whose
/// hinge passes within about half a wavelength (bending into it where the two run alike,
/// a Y where the neighbour runs on; abutting it where they meet across), runs out of
/// fabric or of time, or meets one of its terrain's vents, and there its end plunges.
/// Each hinge comes up when its fold has amplified enough, the deeper sooner, so that the
/// deep hinges grow long and the shallow fill the gaps and link into them; once up they
/// grow together, in the order of the time each step is reached, and who links with whom
/// is decided by who got there first. Where the train converges its hinges merge, and
/// where they run apart they end en echelon. A hinge stands highest where it nucleated
/// and lower where it reached last, so that a linkage leaves a saddle (Grasemann and
/// Schmalholz 2012).
struct DeepFolds {
    scale: f64,
    body_radius: f64,
    /// How far the folds' surface is softened where two hinges are nearly as near, m.
    blend: f64,
    segments: Vec<DeepSegment>,
    buckets: HashMap<(i32, i32, i32), Vec<u32>>,
    /// The arms as they grew.
    arms: Vec<DeepArm>,
}

/// A vector carried into the tangent plane at `p` and normalized.
fn tangential(p: [f64; 3], v: [f64; 3]) -> [f64; 3] {
    let r = dot(v, p);
    normalize([v[0] - r * p[0], v[1] - r * p[1], v[2] - r * p[2]])
}

/// The point `angle` rad from `p` along the great circle that leaves it along `heading`.
fn step_along(p: [f64; 3], heading: [f64; 3], angle: f64) -> [f64; 3] {
    let (s, c) = angle.sin_cos();
    normalize([
        p[0] * c + heading[0] * s,
        p[1] * c + heading[1] * s,
        p[2] * c + heading[2] * s,
    ])
}

/// A fold's profile across its hinge at `s` widths (the Ricker profile): a crest, and a
/// trough each side that takes up what the crest stands above.
fn ricker(s: f64) -> f64 {
    let s2 = s * s;
    (1.0 - s2) * (-0.5 * s2).exp()
}

impl DeepFolds {
    /// A profile is read out to this many widths, where it is under a thousandth.
    const REACH_WIDTHS: f64 = 4.5;
    /// The folds' surface is softened over this share of a fold's width where two hinges
    /// are nearly as near, so that the trough between them has no crease.
    const BLEND_WIDTHS: f64 = 0.25;
    /// A tip whose heading parts from its neighbour's hinge by more than this sine (about
    /// three degrees) is running away from it.
    const PARTING: f64 = 0.05;
    /// The steepest flank of the Ricker profile, per width.
    const STEEPEST: f64 = 1.3801;
    /// A join is not taken where the neighbour ends more than this share of the way
    /// short of where the tip would meet it.
    const SHORTFALL: f64 = 0.2;

    fn key(d: [f64; 3], scale: f64) -> (i32, i32, i32) {
        (
            (d[0] * scale).floor() as i32,
            (d[1] * scale).floor() as i32,
            (d[2] * scale).floor() as i32,
        )
    }

    /// The width of a deep fold's Ricker profile, m: a narrow crest, its troughs close by
    /// and the train's broad troughs between.
    fn width(wavelength: f64) -> f64 {
        POLAR_DEEP_CREST * wavelength
    }

    fn insert(points: &mut HingePoints, d: [f64; 3], arm: usize, k: usize, scale: f64) {
        points
            .entry(Self::key(d, scale))
            .or_default()
            .push((arm as u32, k as u32));
    }

    /// The nearest point within `reach` rad of `d` on the hinge of a line other than
    /// `line`: its arm and index.
    fn nearest_point(
        points: &HingePoints,
        arms: &[DeepArm],
        d: [f64; 3],
        line: usize,
        reach: f64,
        scale: f64,
    ) -> Option<(usize, usize)> {
        let (i, j, l) = Self::key(d, scale);
        let mut best: Option<(f64, usize, usize)> = None;
        for di in -1..=1 {
            for dj in -1..=1 {
                for dl in -1..=1 {
                    let Some(v) = points.get(&(i + di, j + dj, l + dl)) else {
                        continue;
                    };
                    for &(a, k) in v {
                        let (a, k) = (a as usize, k as usize);
                        if arms[a].line == line {
                            continue;
                        }
                        let dist = angle_between(arms[a].points[k], d);
                        if dist < reach && best.is_none_or(|b| dist < b.0) {
                            best = Some((dist, a, k));
                        }
                    }
                }
            }
        }
        best.map(|b| (b.1, b.2))
    }

    /// Where the arc from `p` to `q` first meets the hinge of a line other than the two in
    /// `ignore`, the pieces of hinge sought among the points within `reach` rad of the
    /// arc's middle (the longest piece of hinge is shorter than that).
    fn crossing(
        points: &HingePoints,
        arms: &[DeepArm],
        p: [f64; 3],
        q: [f64; 3],
        ignore: [usize; 2],
        reach: f64,
        scale: f64,
    ) -> Option<[f64; 3]> {
        // The gnomonic plane about `p`, where arcs of great circles are straight.
        let (e1, e2) = crate::geom::local_basis(p);
        let flat = |v: [f64; 3]| {
            let s = dot(v, p);
            (s > 0.0).then(|| (dot(v, e1) / s, dot(v, e2) / s))
        };
        let end = flat(q)?;
        let mid = normalize([p[0] + q[0], p[1] + q[1], p[2] + q[2]]);
        let near = 0.5 * angle_between(p, q) + reach;
        let (i, j, l) = Self::key(mid, scale);
        let mut first: Option<f64> = None;
        for di in -1..=1 {
            for dj in -1..=1 {
                for dl in -1..=1 {
                    let Some(v) = points.get(&(i + di, j + dj, l + dl)) else {
                        continue;
                    };
                    for &(a, k) in v {
                        let (a, k) = (a as usize, k as usize);
                        let hinge = &arms[a].points;
                        if ignore.contains(&arms[a].line) || angle_between(hinge[k], mid) > near {
                            continue;
                        }
                        for (u, w) in [(k.wrapping_sub(1), k), (k, k + 1)] {
                            if u >= hinge.len() || w >= hinge.len() {
                                continue;
                            }
                            let (Some(c), Some(d)) = (flat(hinge[u]), flat(hinge[w])) else {
                                continue;
                            };
                            let r = end;
                            let s = (d.0 - c.0, d.1 - c.1);
                            let denom = r.0 * s.1 - r.1 * s.0;
                            if denom.abs() < 1e-18 {
                                continue;
                            }
                            let t = (c.0 * s.1 - c.1 * s.0) / denom;
                            let v = (c.0 * r.1 - c.1 * r.0) / denom;
                            if t > 1e-9 && t <= 1.0 && (0.0..=1.0).contains(&v) {
                                first = Some(first.map_or(t, |f: f64| f.min(t)));
                            }
                        }
                    }
                }
            }
        }
        first.map(|t| {
            normalize([
                p[0] + t * end.0 * e1[0] + t * end.1 * e2[0],
                p[1] + t * end.0 * e1[1] + t * end.1 * e2[1],
                p[2] + t * end.0 * e1[2] + t * end.1 * e2[2],
            ])
        })
    }

    /// Adds points to an arm, reached at time `t`.
    fn extend(
        arms: &mut [DeepArm],
        points: &mut HingePoints,
        id: usize,
        path: &[[f64; 3]],
        t: f64,
        wavelength: f64,
        scale: f64,
    ) {
        let arm = &mut arms[id];
        for &pt in path {
            arm.points.push(pt);
            arm.times.push(t);
            arm.wavelengths.push(wavelength);
            Self::insert(points, pt, id, arm.points.len() - 1, scale);
        }
    }

    /// An arm's direction at its point `k`, the way it grew.
    fn tangent(arm: &DeepArm, k: usize) -> [f64; 3] {
        let n = arm.points.len();
        if n < 2 {
            return arm.heading;
        }
        let (a, b) = (
            arm.points[k.saturating_sub(1)],
            arm.points[(k + 1).min(n - 1)],
        );
        tangential(arm.points[k], [b[0] - a[0], b[1] - a[1], b[2] - a[2]])
    }

    /// The point `lead` rad along a fold's hinge from point `k` of arm `a`, going the way
    /// the arm grew (`forward`) or back toward its nucleus and on along the arm that grew
    /// the other way, and the hinge's direction there in the sense of travel; at the end
    /// of a hinge, that end, the arm whose end it is and how far short of the lead it
    /// stopped, rad.
    fn walk(
        arms: &[DeepArm],
        a: usize,
        k: usize,
        forward: bool,
        lead: f64,
    ) -> ([f64; 3], [f64; 3], Option<usize>, f64) {
        let (mut arm, mut k, mut forward, mut left) = (a, k, forward, lead);
        loop {
            let points = &arms[arm].points;
            let next = if forward {
                (k + 1 < points.len()).then_some(k + 1)
            } else if k > 0 {
                Some(k - 1)
            } else {
                // Back through the nucleus onto the arm that grew the other way.
                let other = arm ^ 1;
                if other < arms.len()
                    && arms[other].line == arms[arm].line
                    && arms[other].points.len() > 1
                {
                    arm = other;
                    forward = true;
                    continue;
                }
                None
            };
            let Some(n) = next else {
                break;
            };
            let step = angle_between(points[k], points[n]);
            if step >= left {
                let f = left / step.max(1e-12);
                let (u, v) = (points[k], points[n]);
                let at = normalize([
                    u[0] + f * (v[0] - u[0]),
                    u[1] + f * (v[1] - u[1]),
                    u[2] + f * (v[2] - u[2]),
                ]);
                return (
                    at,
                    tangential(at, [v[0] - u[0], v[1] - u[1], v[2] - u[2]]),
                    None,
                    0.0,
                );
            }
            left -= step;
            k = n;
        }
        let points = &arms[arm].points;
        let at = points[k];
        let way = if forward {
            if k > 0 {
                let u = points[k - 1];
                [at[0] - u[0], at[1] - u[1], at[2] - u[2]]
            } else {
                arms[arm].heading
            }
        } else if k + 1 < points.len() {
            let u = points[k + 1];
            [at[0] - u[0], at[1] - u[1], at[2] - u[2]]
        } else {
            let h = arms[arm].heading;
            [-h[0], -h[1], -h[2]]
        };
        (at, tangential(at, way), forward.then_some(arm), left)
    }

    /// A tip that came within linking distance of a neighbour's hinge joins it. Where the
    /// two run alike it bends into the neighbour where its own heading meets it, or, if it
    /// runs alongside, no more steeply than the least merging angle, along a cubic that
    /// leaves the tip along its heading and arrives along the neighbour; beyond the
    /// junction the two are one fold: a Y where the neighbour runs on; where the neighbour
    /// ends there, its end, which then neither plunges nor, if it was coming the other way,
    /// grows on. A neighbour's tip coming the other way past the junction stops too, the
    /// two one fold. Where they meet across the tip runs straight on into the fold ahead.
    /// The tip passes a neighbour that it is running away from, that ends short of where
    /// it would join it or lies behind it, and one it could join only across a third fold:
    /// returns whether the tip linked.
    #[allow(clippy::too_many_arguments)]
    fn link(
        arms: &mut [DeepArm],
        points: &mut HingePoints,
        id: usize,
        (host, k): (usize, usize),
        q: [f64; 3],
        heading: [f64; 3],
        t: f64,
        wavelength: f64,
        body_radius: f64,
        scale: f64,
    ) -> bool {
        let hp = arms[host].points[k];
        let th = Self::tangent(&arms[host], k);
        let forward = dot(th, heading) >= 0.0;
        let along = if forward {
            th
        } else {
            [-th[0], -th[1], -th[2]]
        };
        let angle = dot(along, heading).clamp(-1.0, 1.0).acos();
        let ahead = |p: [f64; 3]| dot([p[0] - q[0], p[1] - q[1], p[2] - q[2]], heading) > 0.0;
        let (line, host_line) = (arms[id].line, arms[host].line);
        let reach = 0.5 / scale;
        if angle > POLAR_DEEP_ABUT_DEG.to_radians() {
            // Across: the tip runs straight on until it meets a fold.
            let far = step_along(q, heading, 2.0 * POLAR_DEEP_LINK * wavelength / body_radius);
            let Some(x) = Self::crossing(points, arms, q, far, [line, line], reach, scale) else {
                return false;
            };
            Self::extend(arms, points, id, &[q, x], t, wavelength, scale);
            arms[id].end = ArmEnd::Linked;
            return true;
        }
        let normal = normalize(cross(hp, along));
        let side = dot([q[0] - hp[0], q[1] - hp[1], q[2] - hp[2]], normal);
        // The sine of the angle at which the tip closes on the neighbour; hinges that run
        // apart do not link.
        let closing = -dot(heading, normal) * side.signum();
        if closing < -Self::PARTING {
            return false;
        }
        let approach = closing
            .clamp(0.0, 1.0)
            .asin()
            .max(POLAR_DEEP_MERGE_DEG.to_radians());
        let lead = side.abs() / approach.tan();
        let (joint, tj, end, short) = Self::walk(arms, host, k, forward, lead);
        if short > Self::SHORTFALL * lead || !ahead(joint) {
            return false;
        }
        let mut path = vec![q];
        let chord = angle_between(q, joint);
        let step = POLAR_DEEP_STEP * wavelength / body_radius;
        let n = ((chord / step).ceil() as usize).max(1);
        for s in 1..=n {
            let u = s as f64 / n as f64;
            let (u2, u3) = (u * u, u * u * u);
            let (h00, h10, h01, h11) = (
                2.0 * u3 - 3.0 * u2 + 1.0,
                u3 - 2.0 * u2 + u,
                3.0 * u2 - 2.0 * u3,
                u3 - u2,
            );
            path.push(normalize([
                h00 * q[0] + h10 * chord * heading[0] + h01 * joint[0] + h11 * chord * tj[0],
                h00 * q[1] + h10 * chord * heading[1] + h01 * joint[1] + h11 * chord * tj[1],
                h00 * q[2] + h10 * chord * heading[2] + h01 * joint[2] + h11 * chord * tj[2],
            ]));
        }
        // A join that would cross a third fold on its way is not taken.
        for i in 1..path.len() {
            let ignore = [line, host_line];
            if Self::crossing(points, arms, path[i - 1], path[i], ignore, reach, scale).is_some() {
                return false;
            }
        }
        Self::extend(arms, points, id, &path, t, wavelength, scale);
        arms[id].end = ArmEnd::Linked;
        let head_on = |a: &DeepArm| a.end == ArmEnd::Growing && dot(a.heading, heading) < 0.0;
        match end {
            Some(e) => {
                if arms[e].end == ArmEnd::Free || head_on(&arms[e]) {
                    arms[e].end = ArmEnd::Linked;
                }
            }
            None => {
                let tip = *arms[host].points.last().expect("an arm holds its nucleus");
                if head_on(&arms[host])
                    && angle_between(tip, joint) < 2.0 * POLAR_DEEP_LINK * wavelength / body_radius
                {
                    arms[host].end = ArmEnd::Linked;
                }
            }
        }
        true
    }

    /// Grows the deep folds from their nuclei through the field the hinges read, the
    /// lid's heterogeneity deflecting them where `wander` is given.
    fn grow(
        nuclei: &[DeepNucleus],
        field: &dyn Fn([f64; 3]) -> DeepSite,
        wander: Option<&GradientNoise>,
        body_radius: f64,
    ) -> Self {
        let sites: Vec<DeepSite> = nuclei.iter().map(|n| field(n.center)).collect();
        let mut lengths: Vec<f64> = sites
            .iter()
            .map(|s| s.wavelength)
            .filter(|&w| w > 0.0)
            .collect();
        if lengths.is_empty() {
            return Self::from_arms(Vec::new(), body_radius);
        }
        lengths.sort_by(f64::total_cmp);
        let longest = lengths[lengths.len() - 1];
        let typical = lengths[lengths.len() / 2];
        // The neighbour test reaches no further than a bucket, twice the longest distance
        // it looks at a nucleus.
        let bucket = 2.0 * POLAR_DEEP_LINK.max(POLAR_DEEP_SPACE) * longest / body_radius;
        let scale = 1.0 / bucket;
        let wander_frequency = body_radius / (POLAR_DEEP_WANDER_LENGTH * typical);
        let max_turn = POLAR_DEEP_TURN_DEG.to_radians();
        let mut points: HingePoints = HashMap::new();
        let mut arms: Vec<DeepArm> = Vec::new();
        // Events in time: a nucleus coming up (0) or an arm's next step (1).
        let mut heap: BinaryHeap<Reverse<(u64, u8, usize)>> = BinaryHeap::new();
        for (i, n) in nuclei.iter().enumerate() {
            let onset = POLAR_DEEP_ONSET * (GROOVE_DEPTH_MAX / n.depth.max(1e-6)).ln().max(0.0);
            if onset < POLAR_DEEP_REACH {
                heap.push(Reverse((onset.to_bits(), 0, i)));
            }
        }
        let mut lines = 0;
        while let Some(Reverse((bits, kind, id))) = heap.pop() {
            if kind == 0 {
                let (n, site, onset) = (nuclei[id], sites[id], f64::from_bits(bits));
                let Some(o) = site.direction else {
                    continue;
                };
                if site.strength < POLAR_DEEP_STOP
                    || site.fissure < POLAR_DEEP_FISSURE_NUCLEUS * site.wavelength
                {
                    continue;
                }
                // A hinge that comes up beside a fold already grown is part of it.
                let space = (POLAR_DEEP_SPACE * site.wavelength / body_radius).min(bucket);
                if Self::nearest_point(&points, &arms, n.center, usize::MAX, space, scale).is_some()
                {
                    continue;
                }
                for sign in [1.0, -1.0] {
                    let arm = arms.len();
                    Self::insert(&mut points, n.center, arm, 0, scale);
                    arms.push(DeepArm {
                        line: lines,
                        depth: n.depth,
                        points: vec![n.center],
                        times: vec![onset],
                        wavelengths: vec![site.wavelength],
                        heights: Vec::new(),
                        heading: [o[0] * sign, o[1] * sign, o[2] * sign],
                        direction: site.direction,
                        end: ArmEnd::Growing,
                    });
                    heap.push(Reverse((bits, 1, arm)));
                }
                lines += 1;
                continue;
            }
            if arms[id].end != ArmEnd::Growing {
                continue;
            }
            let arm = &arms[id];
            let p = *arm.points.last().expect("an arm holds its nucleus");
            let lambda = *arm.wavelengths.last().expect("a wavelength per point");
            let (heading, line) = (arm.heading, arm.line);
            // A hinge grows while the shortening lasts, the deeper the faster.
            let t = f64::from_bits(bits) + POLAR_DEEP_STEP / arm.depth;
            if t > POLAR_DEEP_REACH {
                arms[id].end = ArmEnd::Free;
                continue;
            }
            let Some(o) = arms[id].direction else {
                arms[id].end = ArmEnd::Free;
                continue;
            };
            let o = if dot(o, heading) < 0.0 {
                [-o[0], -o[1], -o[2]]
            } else {
                o
            };
            // The lid's unresolved heterogeneity deflects the hinge as it grows.
            let target = match wander {
                Some(noise) => rotate(
                    o,
                    p,
                    POLAR_DEEP_WANDER_DEG.to_radians()
                        * (noise.fbm(p, wander_frequency, 2, 2.0, 0.5) / 0.5).clamp(-1.0, 1.0),
                ),
                None => o,
            };
            let turn = dot(cross(heading, target), p)
                .atan2(dot(heading, target))
                .clamp(-max_turn, max_turn);
            let heading = tangential(p, rotate(heading, p, turn));
            let q = step_along(p, heading, POLAR_DEEP_STEP * lambda / body_radius);
            let site = field(q);
            if site.direction.is_none()
                || site.strength < POLAR_DEEP_STOP
                || site.fissure < POLAR_DEEP_FISSURE * site.wavelength
            {
                arms[id].end = ArmEnd::Free;
                continue;
            }
            if let Some(x) = Self::crossing(&points, &arms, p, q, [line, line], 0.5 * bucket, scale)
            {
                Self::extend(&mut arms, &mut points, id, &[x], t, site.wavelength, scale);
                arms[id].end = ArmEnd::Linked;
                continue;
            }
            let link = (POLAR_DEEP_LINK * site.wavelength / body_radius).min(bucket);
            if let Some(host) = Self::nearest_point(&points, &arms, q, line, link, scale)
                && Self::link(
                    &mut arms,
                    &mut points,
                    id,
                    host,
                    q,
                    heading,
                    t,
                    site.wavelength,
                    body_radius,
                    scale,
                )
            {
                continue;
            }
            let arm = &mut arms[id];
            arm.points.push(q);
            arm.times.push(t);
            arm.wavelengths.push(site.wavelength);
            arm.heading = tangential(q, heading);
            arm.direction = site.direction;
            let k = arm.points.len() - 1;
            Self::insert(&mut points, q, id, k, scale);
            heap.push(Reverse((t.to_bits(), 1, id)));
        }
        Self::from_arms(arms, body_radius)
    }

    /// The grown arms as segments to read: each point's height from the time the hinge
    /// reached it, a free end plunging.
    fn from_arms(mut arms: Vec<DeepArm>, body_radius: f64) -> Self {
        for arm in &mut arms {
            arm.heights = arm
                .times
                .iter()
                .map(|&t| arm.depth * (1.0 - POLAR_DEEP_SADDLE * (t / POLAR_DEEP_REACH).min(1.0)))
                .collect();
        }
        // A free end plunges: the end of an arm that did not link, and the nucleus of a
        // fold whose other arm never grew.
        let lone: Vec<bool> = (0..arms.len())
            .map(|a| {
                arms.get(a ^ 1)
                    .is_none_or(|o| o.line != arms[a].line || o.points.len() < 2)
            })
            .collect();
        for (a, arm) in arms.iter_mut().enumerate() {
            if arm.end != ArmEnd::Linked {
                Self::plunge(arm, true, body_radius);
            }
            if lone[a] {
                Self::plunge(arm, false, body_radius);
            }
        }
        let mut segments = Vec::new();
        let mut reaches = Vec::new();
        let mut lengths = Vec::new();
        for arm in &arms {
            lengths.push(arm.wavelengths[0]);
            for i in 1..arm.points.len() {
                let segment = DeepSegment {
                    line: arm.line,
                    a: arm.points[i - 1],
                    b: arm.points[i],
                    height: [arm.heights[i - 1], arm.heights[i]],
                    wavelength: [arm.wavelengths[i - 1], arm.wavelengths[i]],
                };
                let widest = Self::width(segment.wavelength[0].max(segment.wavelength[1]));
                reaches.push(
                    0.5 * angle_between(segment.a, segment.b)
                        + Self::REACH_WIDTHS * widest / body_radius,
                );
                segments.push(segment);
            }
        }
        lengths.sort_by(f64::total_cmp);
        let blend = lengths
            .get(lengths.len() / 2)
            .map_or(1.0, |&l| Self::BLEND_WIDTHS * Self::width(l));
        // Each segment goes into every bucket its reach touches, so that a direction reads
        // only its own.
        let mut sorted = reaches.clone();
        sorted.sort_by(f64::total_cmp);
        let scale = sorted.get(sorted.len() / 2).map_or(1.0, |&r| 1.0 / r);
        let mut buckets: HashMap<(i32, i32, i32), Vec<u32>> = HashMap::new();
        for (i, s) in segments.iter().enumerate() {
            let mid = normalize([s.a[0] + s.b[0], s.a[1] + s.b[1], s.a[2] + s.b[2]]);
            let r = reaches[i];
            let lo = Self::key([mid[0] - r, mid[1] - r, mid[2] - r], scale);
            let hi = Self::key([mid[0] + r, mid[1] + r, mid[2] + r], scale);
            for x in lo.0..=hi.0 {
                for y in lo.1..=hi.1 {
                    for z in lo.2..=hi.2 {
                        buckets.entry((x, y, z)).or_default().push(i as u32);
                    }
                }
            }
        }
        Self {
            scale,
            body_radius,
            blend,
            segments,
            buckets,
            arms,
        }
    }

    /// Lowers an arm's heights over the taper from its end (`from_end`) or its nucleus, to
    /// nothing there.
    fn plunge(arm: &mut DeepArm, from_end: bool, body_radius: f64) {
        let n = arm.points.len();
        let mut run = 0.0;
        for step in 0..n {
            let i = if from_end { n - 1 - step } else { step };
            if step > 0 {
                let before = if from_end { i + 1 } else { i - 1 };
                run += angle_between(arm.points[i], arm.points[before]) * body_radius;
            }
            let x = run / (POLAR_DEEP_TAPER * arm.wavelengths[i]);
            if x >= 1.0 {
                break;
            }
            arm.heights[i] *= x * x * (3.0 - 2.0 * x);
        }
    }

    /// The deep folds at a direction: their relief in units of the fold amplitude, and how
    /// fully a deep fold holds the ground, 0 to 1, for the fabric's folds to give way to
    /// it. The folds are one surface: a fold's profile about the distance to the nearest
    /// hinge, softened where two are nearly as near, with the height and wavelength of
    /// the hinges in the same measure, so that two hinges merge where they meet, the one
    /// that runs on unbroken, and a trough between two folds is shared.
    fn at(&self, d: [f64; 3]) -> (f64, f64) {
        // The nearest point of each fold about (up to a few, the nearest kept): its
        // distance, m, and the fold's height and wavelength there.
        const FOLDS: usize = 8;
        let mut nearest = [(usize::MAX, 0.0f64, 0.0f64, 0.0f64); FOLDS];
        let mut count = 0;
        if let Some(v) = self.buckets.get(&Self::key(d, self.scale)) {
            for &k in v {
                let s = &self.segments[k as usize];
                let ab = [s.b[0] - s.a[0], s.b[1] - s.a[1], s.b[2] - s.a[2]];
                let ad = [d[0] - s.a[0], d[1] - s.a[1], d[2] - s.a[2]];
                let len2 = dot(ab, ab);
                let t = if len2 > 0.0 {
                    (dot(ad, ab) / len2).clamp(0.0, 1.0)
                } else {
                    0.0
                };
                let off = [ad[0] - t * ab[0], ad[1] - t * ab[1], ad[2] - t * ab[2]];
                let x = dot(off, off).sqrt() * self.body_radius;
                let wavelength = s.wavelength[0] + t * (s.wavelength[1] - s.wavelength[0]);
                if x > Self::REACH_WIDTHS * Self::width(wavelength) {
                    continue;
                }
                let height = s.height[0] + t * (s.height[1] - s.height[0]);
                let entry = (s.line, x, height, wavelength);
                match nearest[..count].iter_mut().find(|e| e.0 == s.line) {
                    Some(e) => {
                        if x < e.1 {
                            *e = entry;
                        }
                    }
                    None if count < FOLDS => {
                        nearest[count] = entry;
                        count += 1;
                    }
                    None => {
                        let far = (0..FOLDS)
                            .max_by(|&a, &b| nearest[a].1.total_cmp(&nearest[b].1))
                            .expect("a full list");
                        if x < nearest[far].1 {
                            nearest[far] = entry;
                        }
                    }
                }
            }
        }
        if count == 0 {
            return (0.0, 0.0);
        }
        let nearest = &nearest[..count];
        let x0 = nearest.iter().map(|e| e.1).fold(f64::INFINITY, f64::min);
        let (mut sum, mut height, mut wavelength) = (0.0, 0.0, 0.0);
        for &(_, x, h, w) in nearest {
            let p = (-(x - x0) / self.blend).exp();
            sum += p;
            height += p * h;
            wavelength += p * w;
        }
        let (height, wavelength) = (height / sum, wavelength / sum);
        let x = (x0 - self.blend * sum.ln()).max(0.0);
        (
            height * ricker(x / Self::width(wavelength)),
            height.min(1.0) * (-(0.5 * x / Self::width(wavelength)).powi(2)).exp(),
        )
    }
}

/// The fold fabric's strength at a signed distance `dist` from a renewed terrain's edge,
/// m: in full inside past the folds' taper, the halo's share at the edge, and fading
/// across the halo beyond it.
fn fabric_strength(dist: f64, halo: f64) -> f64 {
    if dist >= 0.0 {
        let x = (dist / POLAR_FOLD_TAPER).clamp(0.0, 1.0);
        POLAR_HALO_PEAK + (1.0 - POLAR_HALO_PEAK) * x * x * (3.0 - 2.0 * x)
    } else {
        let x = (1.0 + dist / halo.max(1.0)).clamp(0.0, 1.0);
        POLAR_HALO_PEAK * x * x
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

/// Signed distance to a mask's edge, m, positive inside: the great-circle distance to
/// the nearest cell on the same side that touches the edge, plus the half cell between
/// that cell and the edge, found by carrying each edge cell outward, nearest first, out
/// to `limit` m; beyond it the limit itself. The distance is the true one, so its
/// contours round the corners of the cells rather than step along them.
pub fn edge_distance_m(grid: &CubeSphere, body_radius: f64, mask: &[bool], limit: f64) -> Vec<f32> {
    use std::cmp::Reverse;
    use std::collections::BinaryHeap;
    let n = mask.len();
    let half = (0.5 * grid.delta() * body_radius) as f32;
    let limit = (limit as f32).max(half);
    let mut dist = vec![f32::INFINITY; n];
    let mut heap: BinaryHeap<(Reverse<u32>, u32, u32)> = BinaryHeap::new();
    for k in 0..n {
        let c = grid.cell(k);
        let on_edge = [(1, 0), (-1, 0), (0, 1), (0, -1)]
            .iter()
            .any(|&(di, dj)| mask[grid.index(grid.neighbor(c, di, dj))] != mask[k]);
        if on_edge {
            dist[k] = half;
            heap.push((Reverse(half.to_bits()), k as u32, k as u32));
        }
    }
    while let Some((Reverse(bits), k, source)) = heap.pop() {
        let (k, source) = (k as usize, source as usize);
        if f32::from_bits(bits) > dist[k] {
            continue;
        }
        let origin = grid.center(grid.cell(source));
        let c = grid.cell(k);
        for di in -1..=1 {
            for dj in -1..=1 {
                if di == 0 && dj == 0 {
                    continue;
                }
                let nb = grid.index(grid.neighbor(c, di, dj));
                if mask[nb] != mask[source] {
                    continue;
                }
                let d =
                    (angle_between(origin, grid.center(grid.cell(nb))) * body_radius) as f32 + half;
                if d < dist[nb] && d <= limit {
                    dist[nb] = d;
                    heap.push((Reverse(d.to_bits()), nb as u32, source as u32));
                }
            }
        }
    }
    dist.iter()
        .zip(mask)
        .map(|(&d, &inside)| {
            let d = d.min(limit);
            if inside { d } else { -d }
        })
        .collect()
}

/// How far from a polar terrain's edge its distance is kept, m: the margin ridge's
/// reach, the margin ramp, the folds' taper or the halo, whichever is furthest.
pub fn polar_reach(halo: f64) -> f64 {
    (4.0 * POLAR_MARGIN_WIDTH)
        .max(POLAR_MARGIN_RAMP)
        .max(POLAR_FOLD_TAPER)
        .max(halo)
}

/// How much of a cell a renewed terrain holds at a signed distance from its edge, m: a
/// smooth ramp over the margin width centred on the edge, a half at the edge itself.
pub fn margin_coverage(dist: f64) -> f64 {
    let x = (dist / POLAR_MARGIN_RAMP + 0.5).clamp(0.0, 1.0);
    x * x * (3.0 - 2.0 * x)
}

/// A mask eroded by one cell: a cell stays only with all four neighbours. The coarse loop
/// gives a lane's surface to these alone; the finish draws the true edge around them at
/// its own resolution, so the coarse cells never reach past it.
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
    /// The coronae over the upwellings of an ice mantle's convection (F133), each with its
    /// own time within the heating episode.
    pub coronae: Vec<crate::stages::corona::Corona>,
    /// The equilibrium surface of the thinned shell, which relaxation is not to take for
    /// a load.
    equilibrium: Option<Vec<f32>>,
    /// The surface as the last renewing relaxation step left it, which the terrain
    /// renewed at that step takes for its floor; and how much of each cell's relief the
    /// lanes have renewed at their times, 0 to 1, which the finish completes out to
    /// their resolved edges.
    polar_floor: Option<Vec<f32>>,
    pub lane_core: Option<Vec<f32>>,
    /// How much of its relief below the renewal length each cell gives up to the present
    /// polar terrain when it renews (F60): all of it inside and a ramp across the edge.
    active: Vec<f32>,
    depressed: bool,
    /// The thermal history, for the elastic lithosphere at each lane's time.
    thermal: crate::thermal::ThermalHistory,
    gravity: f64,
    /// When the present episode of focusing began, Gyr since formation: 0 where the
    /// thinning is steady, the episode's length before the present where it renews.
    pub present_start: f64,
    /// The past episodes of focused heating (F53), and the weights of the terrain each
    /// renewed, about where the body left it.
    pub episodes: Vec<FocusEpisode>,
    episode_weights: Vec<Vec<f32>>,
    /// The weights of the terrain renewed at the current step, the present pole's while
    /// its episode lasts and a past episode's at the step its end falls in, with the time
    /// the surface it renews takes, and whether any renews at this step.
    active_now: Vec<f32>,
    renewal_time: Vec<f32>,
    renewing: bool,
    /// The factor the thinning of the terrain renewed at the current step puts on the
    /// shell under each cell, and how long it lasted within the step, Gyr (F66).
    thinned_now: Vec<f32>,
    thinned_for: f64,
    /// The share of the surface the lanes, the present polar terrain and the past
    /// episodes' terrains renewed between them, by area: what is not old cratered ground,
    /// whatever craters it has taken since.
    pub renewed_fraction: f64,
    /// The relief the body formed with, on the loop's grid (F141), and its bands shorter
    /// than each length a renewal has asked for: the texture a renewed cell keeps under
    /// the new terrain.
    primordial: Option<Vec<f32>>,
    fine_cache: Vec<(f64, Vec<f32>)>,
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
            let (lid_lo, lid_hi) = GROOVE_LANE_WIDTH_ELASTIC;
            // The oldest lane to hold each cell, as its index and one: where a lane meets
            // an older one.
            let mut oldest = vec![0u32; n];
            let (o_lo, o_hi) = GROOVE_NEIGHBOUR_OFFSET;
            // The lid a lane breaks is the lid where and when it opens (F61): the elastic
            // lid of the time, thinner by the shell's own thinning where a focusing
            // episode, a past one or the present, has the shell thinned at that time.
            // With it comes the thinning's profile there, which says whether the ground
            // is being renewed.
            let episodes = crate::stages::ice::focusing_episodes(graph, &shell);
            let episode_maps: Vec<ShellMap> =
                episodes.iter().map(|e| shell.with_pole(e.center)).collect();
            let present_start = if shell.has_active_terrain() {
                (age - FOCUSING_EPISODE_DURATION_GYR).max(0.0)
            } else {
                f64::INFINITY
            };
            let lid_at = |d: [f64; 3], t: f64| -> (f64, f64) {
                let mut thinnest = (1.0, 0.0);
                for (e, m) in episodes.iter().zip(&episode_maps) {
                    if t >= e.start && t <= e.end {
                        let f = m.factor(d);
                        if f < thinnest.0 {
                            thinnest = (f, m.profile(d));
                        }
                    }
                }
                if t >= present_start {
                    let f = shell.factor(d);
                    if f < thinnest.0 {
                        thinnest = (f, shell.profile(d));
                    }
                }
                (thermal.at(t).elastic.max(1e3) * thinnest.0, thinnest.1)
            };
            // A lane must be able to be the least of its lengths within the greatest:
            // the widest a lane may be on this body.
            let widest = GROOVE_LANE_MAX_LENGTH_RAD * body_radius / GROOVE_LANE_ASPECT.0;
            // The times the lanes may take, and within them the episodes' own. Where
            // the lid of the time is too thick everywhere for the narrowest lane, only
            // an episode's thinned zone can hold one, and the times are drawn from the
            // episodes, which is the same as drawing from all and keeping those.
            let (t_lo, t_hi) = (
                (te - 0.5 * duration).max(0.0),
                (te + 0.5 * duration).min(age),
            );
            let mut windows: Vec<(f64, f64)> = episodes
                .iter()
                .map(|e| (e.start.max(t_lo), e.end.min(t_hi)))
                .filter(|w| w.1 > w.0)
                .collect();
            if present_start.is_finite() && t_hi > present_start {
                windows.push((present_start.max(t_lo), t_hi));
            }
            let thinnest_lid = (0..=16)
                .map(|i| {
                    thermal
                        .at(t_lo + (t_hi - t_lo) * i as f64 / 16.0)
                        .elastic
                        .max(1e3)
                })
                .fold(f64::INFINITY, f64::min);
            let only_thinned = lid_lo * thinnest_lid > widest && !windows.is_empty();
            let mut misses = 0usize;
            let tau = 2.0 * std::f64::consts::PI;
            with_stamper(&grid, |stamper| {
                for k in 0..GROOVE_MAX_LANES {
                    if area >= fraction * total {
                        break;
                    }
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
                    // The lane's time, where it opens and its width in lid thicknesses
                    // are drawn together, and drawn again until the lid there and then
                    // lets a lane that wide be a lane: a draw on a lid too thick is a
                    // rift that never formed. The first draw of each lane uses the
                    // streams it always did, so a body whose lid lets every lane form is
                    // not changed by the rule.
                    let mut made = None;
                    for attempt in 0..GROOVE_LANE_ATTEMPTS {
                        let st = |base: u64| base + 4096 * attempt;
                        let lids = lid_lo * (lid_hi / lid_lo).powf(rng.uniform(k, st(100)));
                        // A lane of a set belongs to its parent's generation; a lane on
                        // its own takes its time from the episode.
                        let time_gyr = match parent {
                            Some(i) => (lanes[i].time_gyr
                                + (rng.uniform(k, st(111)) - 0.5)
                                    * GROOVE_SET_TIME_SPREAD
                                    * duration)
                                .clamp(0.0, age),
                            None if only_thinned => {
                                let w = windows[((rng.uniform(k, st(113)) * windows.len() as f64)
                                    as usize)
                                    .min(windows.len() - 1)];
                                w.0 + (w.1 - w.0) * rng.uniform(k, st(111))
                            }
                            None => {
                                (te + (rng.uniform(k, st(111)) - 0.5) * duration).clamp(0.0, age)
                            }
                        };
                        // Where it opens: beside an earlier lane, somewhere along it; or
                        // anywhere in the belts of extension, by rejection.
                        let (anchor, beside) = match parent {
                            Some(i) => {
                                let p = &lanes[i];
                                let j = (((0.1 + 0.8 * rng.uniform(k, st(110)))
                                    * p.points.len() as f64)
                                    as usize)
                                    .min(p.points.len() - 1);
                                (p.points[j], Some(p.tangent_at(j)))
                            }
                            None => {
                                let axis = rng.unit_direction(k, st(1));
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
                                for inner in 0..16u64 {
                                    let ang = tau * rng.uniform(k, st(3 + 2 * inner));
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
                                    if rng.uniform(k, st(4 + 2 * inner)) < w.powf(GROOVE_BELT_POWER)
                                    {
                                        break;
                                    }
                                }
                                (center, None)
                            }
                        };
                        // A lane wider than a lane can be is not formed; nor is one that
                        // opens in ground being renewed about it, which keeps no trace.
                        let (lid, profile) = lid_at(anchor, time_gyr);
                        let half_width = 0.5 * lids * lid;
                        if 2.0 * half_width > widest || profile >= POLAR_CORE_PROFILE {
                            continue;
                        }
                        made = Some((time_gyr, anchor, beside, lid, half_width, attempt));
                        break;
                    }
                    let Some((time_gyr, anchor, beside, lid, half_width, attempt)) = made else {
                        misses += 1;
                        if misses > GROOVE_LANE_MISSES {
                            break;
                        }
                        continue;
                    };
                    misses = 0;
                    let st = |base: u64| base + 4096 * attempt;
                    // Its grooves are the necking of the lid it broke (F61).
                    let groove_spacing = GROOVE_NECKING_WAVELENGTH * GROOVE_RIFT_LID_SHARE * lid;
                    let spacing = half_width.min(GROOVE_LANE_POINT_SPACING_RAD * body_radius);
                    // A lane of a set is no longer than its parent.
                    let half_length = match parent {
                        Some(i) => {
                            lanes[i].half_length
                                * (GROOVE_SET_LENGTH_SHARE
                                    + (1.0 - GROOVE_SET_LENGTH_SHARE) * rng.uniform(k, st(101)))
                        }
                        None => (half_width
                            * GROOVE_LANE_ASPECT.0
                            * (GROOVE_LANE_ASPECT.1 / GROOVE_LANE_ASPECT.0)
                                .powf(rng.uniform(k, st(101))))
                        .min(0.5 * GROOVE_LANE_MAX_LENGTH_RAD * body_radius),
                    };
                    // The start and the heading: beside the parent, to one side of it with
                    // the two overlapping, its heading turned a little; or where it
                    // opened, along the strike of the field of extension at its time.
                    let (start, heading) = match (parent, beside) {
                        (Some(i), Some(tj)) => {
                            let p = &lanes[i];
                            let nj = normalize(cross(anchor, tj));
                            let side = if rng.uniform(k, st(108)) < 0.5 {
                                -1.0
                            } else {
                                1.0
                            };
                            let offset = side * (p.half_width + half_width) / body_radius
                                * (o_lo + (o_hi - o_lo) * rng.uniform(k, st(109)));
                            let start = normalize([
                                anchor[0] + nj[0] * offset,
                                anchor[1] + nj[1] * offset,
                                anchor[2] + nj[2] * offset,
                            ]);
                            let turn = GROOVE_NEIGHBOUR_TURN_DEG.to_radians()
                                * (2.0 * rng.uniform(k, st(107)) - 1.0);
                            (start, rotate(tj, start, turn))
                        }
                        _ => {
                            let phase = if duration > 0.0 {
                                ((time_gyr - (te - 0.5 * duration)) / duration).clamp(0.0, 1.0)
                            } else {
                                0.5
                            };
                            (anchor, extension_heading(&strikes, anchor, phase))
                        }
                    };
                    // The centre line is traced on past the lane's length, as far as an
                    // end may look for a lane to stop against, no farther than the
                    // longest a lane may be, with room for an oblique cut.
                    let long = (GROOVE_JUNCTION_REACH * half_length)
                        .min(0.5 * GROOVE_LANE_MAX_LENGTH_RAD * body_radius)
                        .max(half_length)
                        + 2.0 * half_width
                        + 3.0 * cell;
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
                                c.trace(axis, start, phase, spacing, 2.0 * long, k ^ 0x6700_7ACE)
                            })
                        })
                        .filter(|pts| pts.len() >= 3);
                    let tidal = traced.is_some();
                    let points = traced.unwrap_or_else(|| {
                        let axis = normalize(cross(start, heading));
                        Lane::great_circle(axis, start, long, spacing, body_radius)
                    });
                    // The ends (F61, F142). Each way from where the lane began it runs
                    // its length. Where it enters an older lane that crosses it, that
                    // lane more often than not takes up its extension: the end stops at
                    // a cut along the older lane's strike, a little inside it. A lane of
                    // a set ends where its parent does, on the same cut. An end that has
                    // run its length looks farther for one; finding none, it dies out in
                    // the old ground, its faults converging and its grooves fading.
                    let nucleus = points
                        .iter()
                        .enumerate()
                        .max_by(|a, b| dot(*a.1, start).total_cmp(&dot(*b.1, start)))
                        .map_or(0, |(i, _)| i);
                    let run = (half_length / spacing).ceil().max(1.0) as usize;
                    let crossing_cos = GROOVE_JUNCTION_ANGLE_DEG.to_radians().cos();
                    let inherited: Vec<LaneStop> = parent
                        .map(|i| lanes[i].stops.iter().flatten().copied().collect())
                        .unwrap_or_default();
                    let mut ends = [nucleus, nucleus];
                    let mut stops: [Option<LaneStop>; 2] = [None, None];
                    for (e, sign) in [-1i64, 1].into_iter().enumerate() {
                        let avail = if sign < 0 {
                            nucleus
                        } else {
                            points.len() - 1 - nucleus
                        };
                        if avail == 0 {
                            continue;
                        }
                        let at =
                            |step: usize| points[(nucleus as i64 + sign * step as i64) as usize];
                        let ahead = |step: usize| {
                            let (p, back) = (at(step), at(step - 1));
                            normalize([p[0] - back[0], p[1] - back[1], p[2] - back[2]])
                        };
                        let mut passed: Vec<u32> = Vec::new();
                        let own = oldest[grid.index(grid.locate(points[nucleus]).0)];
                        if own != 0 {
                            passed.push(own);
                        }
                        // The step the cut crosses the centre line at and the cut's pole.
                        let mut found: Option<(usize, [f64; 3])> = None;
                        for step in 1..=avail {
                            let p = at(step);
                            if let Some(s) = inherited.iter().find(|s| {
                                dot(points[nucleus], s.pole) > 0.0 && dot(p, s.pole) <= 0.0
                            }) {
                                found = Some((step, s.pole));
                                break;
                            }
                            let j = oldest[grid.index(grid.locate(p).0)];
                            if j == 0 || passed.contains(&j) {
                                continue;
                            }
                            let other = &lanes[(j - 1) as usize];
                            let (_, _, strike) = other.locate(p, body_radius);
                            let cos = dot(ahead(step), strike).abs();
                            let far = step as f64 * spacing
                                >= GROOVE_JUNCTION_MIN_WIDTHS * 2.0 * half_width;
                            if other.time_gyr < time_gyr
                                && cos < crossing_cos
                                && far
                                && (step > run
                                    || rng.uniform(k ^ ((j as u64) << 24), 113)
                                        < GROOVE_JUNCTION_STOP)
                            {
                                let sin = (1.0 - cos * cos).sqrt().max(0.3);
                                let inset = (cell + GROOVE_JUNCTION_INSET * other.half_width) / sin;
                                let cut = (step + (inset / spacing).round() as usize).min(avail);
                                found = Some((cut, normalize(cross(at(cut), strike))));
                                break;
                            }
                            passed.push(j);
                        }
                        let Some((cut, pole)) = found else {
                            // No older lane in reach: the end dies out in the old ground
                            // where its run ends, with no cut.
                            ends[e] = (nucleus as i64 + sign * run.min(avail) as i64) as usize;
                            continue;
                        };
                        // The centre line runs on past the cut until the lane's far
                        // corner has passed it.
                        let q = at(cut);
                        let cos = dot(ahead(cut), normalize(cross(pole, q))).abs().min(0.95);
                        let sin = (1.0 - cos * cos).sqrt();
                        let beyond = half_width * (1.0 + GROOVE_EDGE_REACH) * cos / sin
                            + GROOVE_EDGE_STEP_LIMIT * half_width
                            + 3.0 * cell;
                        let pole = if dot(points[nucleus], pole) < 0.0 {
                            [-pole[0], -pole[1], -pole[2]]
                        } else {
                            pole
                        };
                        stops[e] = Some(LaneStop {
                            pole,
                            reach: 2.0 * beyond + 2.0 * half_width,
                        });
                        let steps = (cut + (beyond / spacing).ceil() as usize).min(avail);
                        ends[e] = (nucleus as i64 + sign * steps as i64) as usize;
                    }
                    let points: Vec<[f64; 3]> = points[ends[0]..=ends[1]].to_vec();
                    if points.len() < 3 {
                        continue;
                    }
                    let lane = Lane::structured(
                        points,
                        spacing,
                        half_width,
                        rng.uniform(k, 102) >= GROOVE_SMOOTH_LANE_SHARE,
                        GROOVE_CLEAN_RANGE.0
                            + (GROOVE_CLEAN_RANGE.1 - GROOVE_CLEAN_RANGE.0) * rng.uniform(k, 112),
                        time_gyr,
                        k,
                        groove_spacing,
                        stops,
                    );
                    // The strike of each later pulse (F63): the field the lane follows,
                    // as it stands later in the lane's working life, against its strike
                    // when the lane opened. The life is a set's spread of the episode,
                    // or what is left of the focusing episode that thinned the lid.
                    let mut lane = lane;
                    let life = windows
                        .iter()
                        .find(|w| time_gyr >= w.0 && time_gyr <= w.1)
                        .map_or(f64::INFINITY, |w| w.1 - time_gyr)
                        .min(GROOVE_SET_TIME_SPREAD * duration);
                    let strike_at = |t: f64| -> [f64; 3] {
                        if tidal && let Some(c) = cracks.as_ref() {
                            let axis = mean_axis_at(t, age, c.nsr_period_gyr);
                            let (_, phase) = c.shell.peak_tension(axis, start);
                            return crack_direction(
                                start,
                                &c.shell.total(axis, start, phase),
                                None,
                            );
                        }
                        let phase = if duration > 0.0 {
                            ((t - (te - 0.5 * duration)) / duration).clamp(0.0, 1.0)
                        } else {
                            0.5
                        };
                        extension_heading(&strikes, start, phase)
                    };
                    let opened = strike_at(time_gyr);
                    lane.turns = (1..GROOVE_PULSES)
                        .map(|j| {
                            let later = strike_at(
                                time_gyr + life * j as f64 / (GROOVE_PULSES - 1).max(1) as f64,
                            );
                            // A strike has no sense: the turn is the lesser way round.
                            let turn = dot(cross(opened, later), start).atan2(dot(opened, later));
                            let half = std::f64::consts::FRAC_PI_2;
                            if turn > half {
                                turn - 2.0 * half
                            } else if turn < -half {
                                turn + 2.0 * half
                            } else {
                                turn
                            }
                        })
                        .collect();
                    let index = lanes.len() as u32 + 1;
                    for sm in lane.sample(&grid, body_radius, 2.0 * cell, stamper) {
                        let d = grid.center(grid.cell(sm.idx));
                        if lane.cover_at(&sm, d, body_radius, cell) >= 0.5 {
                            if !covered[sm.idx] {
                                covered[sm.idx] = true;
                                area += areas[sm.idx];
                            }
                            let held = oldest[sm.idx];
                            if held == 0 || lanes[(held - 1) as usize].time_gyr > lane.time_gyr {
                                oldest[sm.idx] = index;
                            }
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
        // Each renewed terrain's weights on this grid (F60): all of the relief below the
        // renewal length inside it and a ramp across its edge; the ground beside it relaxes
        // under the thinned shell instead (F66).
        let halo = POLAR_HALO_SHELLS * mean_shell;
        let weights_of = |map: &ShellMap| -> Vec<f32> {
            let mask = polar_mask(&grid, body_radius, map, &fissures);
            edge_distance_m(&grid, body_radius, &mask, polar_reach(halo))
                .iter()
                .map(|&d| margin_coverage(d as f64) as f32)
                .collect()
        };
        let active = if shell.has_active_terrain() {
            weights_of(&shell)
        } else {
            vec![0.0; n]
        };
        let episodes = crate::stages::ice::focusing_episodes(graph, &shell);
        let episode_weights: Vec<Vec<f32>> = episodes
            .iter()
            .map(|e| weights_of(&shell.with_pole(e.center)))
            .collect();
        let present_start = if shell.has_active_terrain() {
            (age - FOCUSING_EPISODE_DURATION_GYR).max(0.0)
        } else {
            0.0
        };
        // The coronae over an ice mantle's upwellings (F133), and the cells they hold.
        let coronae = crate::stages::corona::ice_coronae(graph, &thermal);
        let corona_cell = grid.delta() * body_radius;
        with_stamper(&grid, |stamper| {
            for c in &coronae {
                for (idx, cover, _) in c.covers(&grid, body_radius, corona_cell, stamper) {
                    if cover >= 0.5 {
                        covered[idx] = true;
                    }
                }
            }
        });
        let renewed_fraction = (0..n)
            .filter(|&k| {
                covered[k] || active[k] >= 0.5 || episode_weights.iter().any(|w| w[k] >= 0.5)
            })
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
            coronae,
            equilibrium: None,
            polar_floor: None,
            lane_core: None,
            active,
            depressed: false,
            thermal,
            gravity: graph.f64("body.gravity"),
            present_start,
            episodes,
            episode_weights,
            active_now: vec![0.0; n],
            renewal_time: vec![0.0; n],
            renewing: false,
            thinned_now: vec![1.0; n],
            thinned_for: 0.0,
            renewed_fraction,
            primordial: None,
            fine_cache: Vec::new(),
        })
    }

    /// The relief the body formed with (F141), whose texture the renewals keep.
    pub fn set_primordial(&mut self, bedrock: &[f32]) {
        self.primordial = Some(bedrock.to_vec());
        self.fine_cache.clear();
    }

    /// The base relief's bands shorter than a length, m, on the loop's grid: the relief
    /// less its low-pass over that length, kept per level of the pyramid the low-pass
    /// reads.
    fn fine_index(&mut self, length: f64) -> Option<usize> {
        let p = self.primordial.as_ref()?;
        let cell = self.grid.delta() * self.body_radius;
        let mut level_length = 4.0 * cell;
        while level_length < length && level_length < self.body_radius * 8.0 {
            level_length *= 2.0;
        }
        if let Some(i) = self.fine_cache.iter().position(|(l, _)| *l == level_length) {
            return Some(i);
        }
        let low =
            crate::stages::relaxation::low_pass(&field_from(self.grid, p), cell, level_length);
        let fine: Vec<f32> = p.iter().zip(&low).map(|(a, b)| a - b).collect();
        self.fine_cache.push((level_length, fine));
        Some(self.fine_cache.len() - 1)
    }

    fn fine_relief(&mut self, length: f64) -> Option<Vec<f32>> {
        let i = self.fine_index(length)?;
        Some(self.fine_cache[i].1.clone())
    }

    /// What the relaxation step takes from the shell: the thinned shell's equilibrium
    /// surface, and the terrain that renews at this step with its weights (F60), the
    /// present pole's while its episode lasts and a past episode's at the step its end
    /// falls in, each once, whatever the step's length. The loop relaxes with the mean
    /// shell (see `relaxation::Local`); the local thickness acts through the replayed
    /// craters' survival, the crack density and the renewal.
    pub fn relaxation(&mut self, t0: f64, t1: f64) -> crate::stages::relaxation::Local<'_> {
        let grid = self.grid;
        self.active_now.iter_mut().for_each(|w| *w = 0.0);
        self.thinned_now.iter_mut().for_each(|f| *f = 1.0);
        self.thinned_for = 0.0;
        self.renewing = false;
        for (e, weights) in self.episodes.iter().zip(&self.episode_weights) {
            if e.end > t0 && e.end <= t1 {
                self.renewing = true;
                for (k, &w) in weights.iter().enumerate() {
                    if w > self.active_now[k] {
                        self.active_now[k] = w;
                        self.renewal_time[k] = e.end as f32;
                    }
                }
                // The episode's thinned shell under the ground beside its terrain, for as
                // long as the episode lasted (F66).
                let map = self.shell.with_pole(e.center);
                for (k, f) in self.thinned_now.iter_mut().enumerate() {
                    *f = f.min(map.factor(grid.center(grid.cell(k))) as f32);
                }
                self.thinned_for = self.thinned_for.max(e.end - e.start);
            }
        }
        if self.shell.has_active_terrain() && t1 >= self.present_start {
            self.renewing = true;
            let now = (t1 - POLAR_RENEWAL_GYR).max(0.0) as f32;
            for (k, &w) in self.active.iter().enumerate() {
                if w > 0.0 && w >= self.active_now[k] {
                    self.active_now[k] = w;
                    self.renewal_time[k] = now;
                }
            }
            for (k, f) in self.thinned_now.iter_mut().enumerate() {
                *f = f.min(self.shell.factor(grid.center(grid.cell(k))) as f32);
            }
            self.thinned_for = self.thinned_for.max(t1 - t0.max(self.present_start));
        }
        // A renewed terrain keeps nothing up to its own scale while it is hot, the zone's
        // radius, and takes the base relief's texture below the scale of its folds (F141).
        let length = self
            .smooth_length()
            .max(self.shell.radius_rad * self.body_radius);
        let cell = self.grid.delta() * self.body_radius;
        let fine = if self.renewing {
            self.fine_index(FABRIC_TEXTURE_SHARE * POLAR_FOLD_SPACING.max(GROOVE_MIN_CELLS * cell))
        } else {
            None
        };
        crate::stages::relaxation::Local {
            equilibrium: self.equilibrium.as_deref(),
            renew: self
                .renewing
                .then_some((self.active_now.as_slice(), length)),
            halo: self
                .renewing
                .then_some((self.thinned_now.as_slice(), self.thinned_for)),
            fine: fine.map(|i| self.fine_cache[i].1.as_slice()),
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
        // The thinned shell floats lower while the focusing lasts (F43): the present
        // zone's depression from its episode's start (from the first step where the
        // thinning is steady). The shell floats to its equilibrium at once, the floor the
        // relaxation step left goes down with it, and the relaxation keeps that
        // equilibrium rather than taking the depression for a load. A past episode's
        // depression came and went within its own step and leaves nothing.
        if !self.depressed && self.shell.pole.is_some() && t1 >= self.present_start {
            self.depressed = true;
            let mean = self.mean_shell;
            let eq: Vec<f32> = (0..n)
                .map(|k| -(self.shell.depression(grid.center(grid.cell(k)), mean) as f32))
                .collect();
            for (h, e) in layers.bedrock.as_mut_slice().iter_mut().zip(&eq) {
                *h += e;
            }
            if let Some(floor) = &mut self.polar_floor {
                for (h, e) in floor.iter_mut().zip(&eq) {
                    *h += e;
                }
            }
            self.equilibrium = Some(eq);
        }
        // The lanes whose time falls in this step, youngest last: each takes the relief as
        // it is now, smoothed and dropped, by its cover of each cell across its border
        // zone (F60); a cell a younger lane also takes ends up with the younger lane's
        // time.
        let mut due: Vec<usize> = (0..self.lanes.len())
            .filter(|&i| {
                let t = self.lanes[i].time_gyr;
                t > t0 && t <= t1
            })
            .collect();
        if !due.is_empty() {
            due.sort_by(|&a, &b| self.lanes[a].time_gyr.total_cmp(&self.lanes[b].time_gyr));
            // Each lane's floor (F141): the relief smoothed over what the lane levels, its
            // own width or the lanes' smoothing length, with the base relief's texture
            // below the scale of its grooves under that.
            let body_radius = self.body_radius;
            let longest = due
                .iter()
                .map(|&i| self.lanes[i].removal_length(body_radius))
                .fold(0.0, f64::max);
            let levels = crate::stages::relaxation::low_passes(&layers.bedrock, cell, longest);
            let fines: Vec<Option<Vec<f32>>> = due
                .iter()
                .map(|&i| {
                    let scale = FABRIC_TEXTURE_SHARE * self.lanes[i].grooves_on(cell);
                    self.fine_relief(scale)
                })
                .collect();
            let mut weight = vec![0.0f32; n];
            let mut target = vec![0.0f32; n];
            let mut held: Vec<Option<(f32, f32)>> = vec![None; n];
            with_stamper(&grid, |stamper| {
                for (n_due, &i) in due.iter().enumerate() {
                    let lane = &self.lanes[i];
                    let length = lane.removal_length(body_radius);
                    let fine = fines[n_due].as_deref();
                    for sm in lane.sample(&grid, body_radius, lane.reach(cell), stamper) {
                        let d = grid.center(grid.cell(sm.idx));
                        let c = lane.cover_at(&sm, d, body_radius, cell) as f32;
                        if c > weight[sm.idx] {
                            weight[sm.idx] = c;
                            target[sm.idx] =
                                crate::stages::relaxation::blend_levels(&levels, length, sm.idx)
                                    + fine.map_or(0.0, |f| f[sm.idx]);
                        }
                        if c >= 0.5 {
                            held[sm.idx] = Some((lane.time_gyr as f32, lane.clean as f32));
                        }
                    }
                }
            });
            // The surface is the lanes' on the cells wholly inside one, a cell in from
            // the edge, so that the loop's cells never reach past the edge the finish
            // draws; the relief is renewed by the weight, out to the edge and across it.
            let inside: Vec<bool> = weight.iter().map(|&w| w >= 0.5).collect();
            let core = eroded(&grid, &inside);
            // The lanes' drop, the extension's subsidence, spreads over the lithosphere's
            // flexural parameter at the time, so a lane narrower than it sinks little and
            // none ends on a step.
            let te = self.thermal.at(0.5 * (t0 + t1)).elastic;
            let alpha = flexural_parameter(te, self.gravity);
            let ramp = crate::stages::relaxation::low_pass(&field_from(grid, &weight), cell, alpha);
            let renewed = self.lane_core.get_or_insert_with(|| vec![0.0; n]);
            let bedrock = layers.bedrock.as_mut_slice();
            let ages = layers.surface_age.as_mut_slice();
            let units = layers.unit_id.as_mut_slice();
            let deep = layers.deep.as_mut_slice();
            let melt = layers.melt.as_mut_slice();
            let cover = layers.cover.as_mut_slice();
            let clean = layers.clean.as_mut_slice();
            for k in 0..n {
                let w = weight[k];
                if w > 0.0 {
                    bedrock[k] += (target[k] - bedrock[k]) * w;
                    renewed[k] = renewed[k].max(w);
                }
                bedrock[k] -= (GROOVE_LANE_DROP * ramp[k] as f64) as f32;
                if core[k]
                    && let Some((t, cl)) = held[k]
                {
                    ages[k] = t;
                    units[k] = UNIT_GROOVED;
                    deep[k] = 0.0;
                    melt[k] = 0.0;
                    cover[k] = 0.0;
                    clean[k] = cl;
                }
            }
        }
        // The coronae whose time falls in this step (F133), youngest last: each renews the
        // relief inside its outline, smoothed over its own length, by its cover of each
        // cell, as a lane does, and holds the cells wholly inside it.
        let mut due: Vec<usize> = (0..self.coronae.len())
            .filter(|&i| {
                let t = self.coronae[i].time_gyr;
                t > t0 && t <= t1
            })
            .collect();
        if !due.is_empty() {
            due.sort_by(|&a, &b| {
                self.coronae[a]
                    .time_gyr
                    .total_cmp(&self.coronae[b].time_gyr)
            });
            let body_radius = self.body_radius;
            let mut weight = vec![0.0f32; n];
            let mut target = vec![0.0f32; n];
            let mut held: Vec<Option<(f32, f32)>> = vec![None; n];
            // Each corona's floor: the relief smoothed over its own length, with the base
            // relief's texture below the scale of its faults under that (F141).
            let fines: Vec<Option<Vec<f32>>> = due
                .iter()
                .map(|&i| {
                    let scale = FABRIC_TEXTURE_SHARE * self.coronae[i].faults_on(cell);
                    self.fine_relief(scale)
                })
                .collect();
            with_stamper(&grid, |stamper| {
                for (n_due, &i) in due.iter().enumerate() {
                    let c = &self.coronae[i];
                    let smoothed: Vec<f32> = crate::stages::relaxation::low_pass(
                        &layers.bedrock,
                        cell,
                        c.smooth_length(),
                    );
                    let fine = fines[n_due].as_deref();
                    for (idx, w, core_share) in c.covers(&grid, body_radius, cell, stamper) {
                        let w = w as f32;
                        if w > weight[idx] {
                            weight[idx] = w;
                            target[idx] = smoothed[idx] + fine.map_or(0.0, |f| f[idx]);
                        }
                        if w >= 0.5 {
                            // The core's ice is clean, renewed from below; the belt is the
                            // old crust broken into blocks (F143).
                            held[idx] = Some((c.time_gyr as f32, (c.clean * core_share) as f32));
                        }
                    }
                }
            });
            let inside: Vec<bool> = weight.iter().map(|&w| w >= 0.5).collect();
            let core = eroded(&grid, &inside);
            let renewed = self.lane_core.get_or_insert_with(|| vec![0.0; n]);
            let bedrock = layers.bedrock.as_mut_slice();
            let ages = layers.surface_age.as_mut_slice();
            let units = layers.unit_id.as_mut_slice();
            let deep = layers.deep.as_mut_slice();
            let melt = layers.melt.as_mut_slice();
            let cover = layers.cover.as_mut_slice();
            let clean = layers.clean.as_mut_slice();
            for k in 0..n {
                let w = weight[k];
                if w > 0.0 {
                    bedrock[k] += (target[k] - bedrock[k]) * w;
                    renewed[k] = renewed[k].max(w);
                }
                if core[k]
                    && let Some((t, cl)) = held[k]
                {
                    ages[k] = t;
                    units[k] = UNIT_CORONA;
                    deep[k] = 0.0;
                    melt[k] = 0.0;
                    cover[k] = 0.0;
                    clean[k] = cl;
                }
            }
        }
    }

    /// After the step's events: the terrain that renewed at this step takes the floor the
    /// relaxation left it, erasing what the step brought it, with the time of its
    /// renewal.
    pub fn after(&mut self, layers: &mut LayerSet) {
        if !self.renewing {
            return;
        }
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
        for (k, &w) in self.active_now.iter().enumerate() {
            if w >= 0.5 {
                bedrock[k] = floor[k];
                ages[k] = self.renewal_time[k];
                units[k] = UNIT_POLAR_TERRAIN;
                deep[k] = 0.0;
                melt[k] = 0.0;
                cover[k] = 0.0;
                clean[k] = 1.0;
            }
        }
    }
}

/// What the renewed terrains leave on the resolved grid (F60, F61), composed in the
/// order they formed: for each cell, the latest time a lane or a polar terrain took more
/// than half of it, with that terrain's unit; the clean ice on it, each terrain's laid
/// over what was there by its cover, and brightened where a lane's fractures break the
/// old ground; the lanes' cover of the cell with the time of the lane that holds it,
/// for laying the lanes' floor; and what the old fine relief keeps where renewals have
/// touched the cell without taking it, across their edges and their halos.
pub struct Renewal {
    pub time: Vec<f32>,
    unit: Vec<u16>,
    clean: Vec<f32>,
    lane_cover: Vec<f32>,
    lane_time: Vec<f32>,
    /// The length the terrain that holds the cell levels the relief over, m, and the
    /// scale of its fabric, m, below which the base relief keeps its texture (F141).
    lane_length: Vec<f32>,
    lane_fabric: Vec<f32>,
    pub keep: Vec<f32>,
}

/// A lane's hold of a cell: its cover, the tone of its ice there, and the clean ice its
/// fractures show beyond its edge.
#[derive(Clone, Copy, Debug)]
pub struct LaneHit {
    pub(crate) k: u32,
    pub(crate) cover: f32,
    pub(crate) tone: f32,
    pub(crate) fracture: f32,
}

impl Renewal {
    pub fn new(cells: usize) -> Self {
        Self {
            time: vec![f32::NEG_INFINITY; cells],
            unit: vec![UNIT_NONE; cells],
            clean: vec![0.0; cells],
            lane_cover: vec![0.0; cells],
            lane_time: vec![f32::NEG_INFINITY; cells],
            lane_length: vec![0.0; cells],
            lane_fabric: vec![0.0; cells],
            keep: vec![1.0; cells],
        }
    }

    /// A terrain's hold of a cell, taken in the order the terrains formed: its ice laid
    /// over what was there by its cover; the cell its own, with its time and unit, where
    /// it takes more than half, and the old fine relief less by its weight where it does
    /// not.
    fn take(&mut self, k: usize, cover: f32, weight: f32, tone: f32, time: f32, unit: u16) {
        self.clean[k] += (tone - self.clean[k]) * cover;
        if cover >= 0.5 {
            self.time[k] = time;
            self.unit[k] = unit;
            self.keep[k] = 1.0;
        } else {
            self.keep[k] *= 1.0 - weight.clamp(0.0, 1.0);
        }
    }

    /// A lane's hold of the cells about it.
    pub fn take_lane(&mut self, lane: &Lane, hits: &[LaneHit], body_radius: f64, cell: f64) {
        self.take_hits(
            lane.time_gyr,
            UNIT_GROOVED,
            hits,
            lane.removal_length(body_radius),
            lane.grooves_on(cell),
        );
    }

    /// A terrain's hold of the cells about it, a lane's or a corona's (F133), each a
    /// floor of its own with its cover: levelled over `length` m, with the base relief's
    /// texture below `fabric` m, the scale of its grooves or faults (F141).
    pub fn take_hits(
        &mut self,
        time_gyr: f64,
        unit: u16,
        hits: &[LaneHit],
        length: f64,
        fabric: f64,
    ) {
        let time = time_gyr as f32;
        for h in hits {
            let k = h.k as usize;
            if h.cover > 0.0 {
                self.take(k, h.cover, h.cover, h.tone, time, unit);
                // The lane that holds the cell for its floor: the latest of those that
                // take more than half of it, else the one that takes most.
                if h.cover >= 0.5 || (self.lane_cover[k] < 0.5 && h.cover > self.lane_cover[k]) {
                    self.lane_time[k] = time;
                    self.lane_length[k] = length as f32;
                    self.lane_fabric[k] = fabric as f32;
                }
                self.lane_cover[k] = self.lane_cover[k].max(h.cover);
            }
            if h.fracture > self.clean[k] {
                self.clean[k] = h.fracture;
            }
        }
    }

    /// A polar terrain's hold of the cells within its reach, from the signed distance to
    /// its edge, m: its clean ice and its weight on the old relief, across the margin's
    /// ramp.
    pub fn take_polar(&mut self, edge: &[f32], time: f64) {
        for (k, &d) in edge.iter().enumerate() {
            let cover = margin_coverage(d as f64) as f32;
            if cover <= 0.0 {
                continue;
            }
            self.take(k, cover, cover, 1.0, time as f32, UNIT_POLAR_TERRAIN);
        }
    }

    /// Write the terrains' surface fields: the clean substrate everywhere; the unit and
    /// the time where a terrain holds the cell and nothing younger has taken it since. A
    /// cell the coarse loop gave a terrain past its resolved edge goes back to the old
    /// ground.
    pub fn apply(&self, layers: &mut LayerSet) {
        let ages = layers.surface_age.as_mut_slice();
        let units = layers.unit_id.as_mut_slice();
        let deep = layers.deep.as_mut_slice();
        let melt = layers.melt.as_mut_slice();
        let cover = layers.cover.as_mut_slice();
        let clean = layers.clean.as_mut_slice();
        for k in 0..self.time.len() {
            clean[k] = self.clean[k];
            if self.time[k] > f32::NEG_INFINITY {
                if ages[k] <= self.time[k] + 1e-3 {
                    ages[k] = self.time[k];
                    units[k] = self.unit[k];
                    deep[k] = 0.0;
                    melt[k] = 0.0;
                    cover[k] = 0.0;
                }
            } else if matches!(units[k], UNIT_GROOVED | UNIT_POLAR_TERRAIN | UNIT_CORONA) {
                units[k] = UNIT_PRIMORDIAL;
            }
        }
    }
}

/// The lanes on the resolved grid, before the finish stamps anything: each lane's hold
/// of the cells about it, for the renewal to take in the order the lanes formed. The
/// coarse loop drew the lanes cell by cell; here their edges, their tone and their
/// damage zones are drawn at the finishing resolution, so that the crater population
/// that follows reads the ages they leave and the colour follows the true edge.
pub fn lane_hits(
    grid: &CubeSphere,
    body_radius: f64,
    cell: f64,
    lanes: &[Lane],
    loose: Option<&[f32]>,
) -> Vec<Vec<LaneHit>> {
    lanes
        .par_iter()
        .map(|lane| {
            let noise = LaneNoise::new(lane);
            with_stamper(grid, |stamper| {
                lane.sample(grid, body_radius, lane.reach(cell), stamper)
                    .into_iter()
                    .filter_map(|sm| {
                        let d = grid.center(grid.cell(sm.idx));
                        let snow = loose.map_or(0.0, |l| l[sm.idx] as f64);
                        let f = lane.fabric_at(&sm, d, body_radius, cell, snow, &noise);
                        (f.cover > 0.0 || f.fracture > 0.004).then_some(LaneHit {
                            k: sm.idx as u32,
                            cover: f.cover as f32,
                            tone: f.tone as f32,
                            fracture: f.fracture as f32,
                        })
                    })
                    .collect()
            })
        })
        .collect()
}

/// The length a lane smooths the relief it renews over, m.
pub fn lane_smooth_length(body_radius: f64) -> f64 {
    GROOVE_SMOOTH_LENGTH.min(0.05 * body_radius)
}

/// The lanes' floor on the resolved grid (F60), before anything is stamped. The loop
/// renewed each lane's relief by its cover on its own cells (`renewed`, how much of each
/// cell's relief it took, resampled and smoothed over the loop's cell); here the same
/// renewal is completed out to the lane's resolved edge, and in the lanes too narrow for
/// the loop's cells: where a lane is the latest thing to have renewed a cell, what the
/// loop left of the relief shorter than the lanes' smoothing length gives way to the
/// mean of the ground about it, until the cell has given up its relief by the lane's
/// resolved cover. The mean holds the regional relief on both sides of the edge, so
/// the floor stands on what surrounds it and no step is laid; nothing is kept from
/// the loop but the relief itself. A crater younger than the lane keeps its ground,
/// floor and thick ejecta alike, and is left out of the mean; so does a cell
/// something younger renewed. Each cell's floor is levelled over the length its terrain
/// levels, a lane's width or a corona's, and under it the base relief keeps its texture
/// (F141): its octaves below the scale of the terrain's fabric where the floor is laid
/// here, and of those only the ones below the loop's cells, which the loop's grid could
/// not hold, where the loop laid it (`texture`: the relief the body formed with, and the
/// loop's cell, m).
pub fn lay_lane_floor(
    layers: &mut LayerSet,
    body_radius: f64,
    renewed: Option<&[f32]>,
    renewal: &Renewal,
    texture: Option<(&crate::stages::primordial_crust::PrimordialCrust, f64)>,
) {
    let grid = layers.grid;
    let n = grid.cell_count();
    let cell = grid.delta() * body_radius;
    let units = layers.unit_id.as_slice();
    let ages = layers.surface_age.as_slice();
    let mut weight = vec![0.0f32; n];
    let mut taken = vec![0.0f32; n];
    let mut counted = vec![1.0f32; n];
    let mut any = false;
    for k in 0..n {
        let c = renewal.lane_cover[k];
        if c <= 0.0 {
            continue;
        }
        let time = renewal.lane_time[k];
        if is_crater_unit(units[k]) && ages[k] > time + 1e-3 {
            counted[k] = 0.0;
            continue;
        }
        if renewal.time[k] > time + 1e-3 {
            continue;
        }
        // What is left to take, of what the loop left, for the cell to have given up
        // its relief by the lane's cover; less than that where the loop left little.
        let q = renewed.map_or(0.0, |r| r[k]).clamp(0.0, 1.0);
        weight[k] = ((c - q) / (1.0 - q).max(GROOVE_FLOOR_RESIDUE as f32)).clamp(0.0, 1.0);
        taken[k] = q * (1.0 - weight[k]);
        any |= weight[k] > 0.0 || taken[k] > 0.0;
    }
    if !any {
        return;
    }
    let bedrock = layers.bedrock.as_mut_slice();
    let relief: Vec<f32> = (0..n).map(|k| bedrock[k] * counted[k]).collect();
    let floor_length =
        |k: usize| (renewal.lane_length[k] as f64).max(lane_smooth_length(body_radius));
    let longest = (0..n)
        .filter(|&k| weight[k] > 0.0)
        .map(floor_length)
        .fold(lane_smooth_length(body_radius), f64::max);
    let sums = crate::stages::relaxation::low_passes(&field_from(grid, &relief), cell, longest);
    let shares = crate::stages::relaxation::low_passes(&field_from(grid, &counted), cell, longest);
    // The textures, by the first octave below each terrain's fabric: all of them, and
    // those below the loop's cells.
    let mut textures: HashMap<u32, (Vec<f32>, Vec<f32>)> = HashMap::new();
    for k in 0..n {
        if weight[k] <= 0.0 && taken[k] <= 0.0 {
            continue;
        }
        let (all, sub) = match texture {
            Some((crust, loop_cell)) => {
                let first = crust.first_octave(
                    FABRIC_TEXTURE_SHARE * renewal.lane_fabric[k] as f64,
                    body_radius,
                );
                let (all, sub) = textures.entry(first).or_insert_with(|| {
                    let below_loop = first.max(crust.first_octave(2.0 * loop_cell, body_radius));
                    (
                        crust.fine_relief_from(&grid, first),
                        crust.fine_relief_from(&grid, below_loop),
                    )
                });
                (all[k], sub[k])
            }
            None => (0.0, 0.0),
        };
        if weight[k] > 0.0 {
            let length = floor_length(k);
            let sum = crate::stages::relaxation::blend_levels(&sums, length, k);
            let share = crate::stages::relaxation::blend_levels(&shares, length, k);
            if share > 0.05 {
                bedrock[k] += (sum / share + all - bedrock[k]) * weight[k];
            }
        }
        if taken[k] > 0.0 {
            bedrock[k] += sub * taken[k];
        }
    }
}

/// A renewed polar terrain on the resolved grid, the present one or a past episode's:
/// its outline drawn at this resolution from its shell map and the fissures, as the
/// signed distance to its edge, m, out to the halo's reach.
pub fn polar_edge(
    grid: &CubeSphere,
    body_radius: f64,
    shell: &ShellMap,
    cracks: &[Linea],
    halo: f64,
) -> Vec<f32> {
    let mask = polar_mask(grid, body_radius, shell, cracks);
    if std::env::var_os("EXOGEN_DEBUG").is_some() {
        let inside = mask.iter().filter(|&&m| m).count();
        eprintln!("polar terrain on the fine grid: {inside} cells inside");
    }
    edge_distance_m(grid, body_radius, &mask, polar_reach(halo))
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

/// A lane's fabric on the resolved grid at its time (F61). Inside its bounding faults
/// the lane renews the ground: whatever fine relief stood there, an older lane's
/// grooves or a small crater's bowl, gives way by the lane's cover to the floor (`base`,
/// the relief before anything was stamped) with the lane's own grooves on it. Beyond the
/// faults the flank's uplift and the damage zone's fractures are added to the ground as
/// it stands. A cell something younger renewed is left alone; a younger crater cuts and
/// buries the fabric when its own turn comes.
#[allow(clippy::too_many_arguments)]
pub fn groove_compute(
    layers: &LayerSet,
    grid: &CubeSphere,
    body_radius: f64,
    lane: &Lane,
    time_gyr: f64,
    cell: f64,
    renewed: &[f32],
    base: &[f32],
    loose: Option<&[f32]>,
    stamper: &mut Stamper,
) -> StampOutput {
    let bedrock = layers.bedrock.as_slice();
    let noise = LaneNoise::new(lane);
    let samples = lane.sample(grid, body_radius, lane.reach(cell), stamper);
    let mut cells: Vec<CellWrite> = samples
        .par_iter()
        .filter_map(|sm| {
            let k = sm.idx;
            // A cell something younger renewed keeps that.
            if renewed[k] as f64 > time_gyr + 1e-3 {
                return None;
            }
            let d = grid.center(grid.cell(k));
            let snow = loose.map_or(0.0, |l| l[k] as f64);
            let f = lane.fabric_at(sm, d, body_radius, cell, snow, &noise);
            if f.cover <= 0.0 && f.outer.abs() < 0.5 {
                return None;
            }
            let h = bedrock[k] as f64;
            Some(write(
                k,
                h + (base[k] as f64 + f.relief - h) * f.cover + f.outer * (1.0 - f.cover),
            ))
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

/// A polar terrain's fabric on the resolved grid at its time (F43, F60), over the relief
/// the loop left it: folds along its fissures and the roughness of a plain deformed at
/// every scale, in full inside and fading across the halo beyond its edge; the belt
/// where the spreading is taken up, about the edge; and the margin ridge where the
/// spreading terrain meets the old. A past terrain's structures are spread and lowered
/// by the gardening since (`kappa_t`, m²). The terrain renews the ground by its weight:
/// whatever fine relief stood there gives way to the floor the loop left (`base`, the
/// relief before anything was stamped), wholly inside and less across its edge and its
/// halo, so that an older crater beside it is softened and one inside it is gone. A cell
/// something younger renewed is left alone; a younger crater cuts and buries the fabric
/// when its own turn comes. Under the floor the base relief keeps the octaves below the
/// loop's cells (`fine`, F141), which the loop's floor could not hold.
#[allow(clippy::too_many_arguments)]
pub fn polar_compute(
    layers: &LayerSet,
    grid: &CubeSphere,
    body_radius: f64,
    shell: &ShellMap,
    edge: &[f32],
    cracks: &[Linea],
    time_gyr: f64,
    kappa_t: f64,
    cell: f64,
    halo: f64,
    renewed: &[f32],
    base: &[f32],
    fine: Option<&[f32]>,
    stamper: &mut Stamper,
) -> StampOutput {
    let Some(pole) = shell.pole else {
        return StampOutput::default();
    };
    let bedrock = layers.bedrock.as_slice();
    let limit = polar_reach(halo);
    let reach = shell.radius_rad;
    let spacing = POLAR_FOLD_SPACING.max(GROOVE_MIN_CELLS * cell);
    let belt_spacing = POLAR_BELT_SPACING.max(3.0 * cell);
    let fabric = FabricIndex::new(cracks, POLAR_FABRIC_RANGE / body_radius);
    // The terrain's own vents, which opened in its episode, bound its deep folds (F65).
    let vents = FabricIndex::new(
        cracks.iter().filter(|l| {
            l.kind == LineaKind::Vent
                && (l.time_gyr - time_gyr).abs() < FOCUSING_EPISODE_DURATION_GYR + 1e-3
        }),
        POLAR_FABRIC_RANGE / body_radius,
    );
    let clock = std::time::Instant::now();
    let ridge_keep = degraded(POLAR_MARGIN_WIDTH, kappa_t);
    let belt_keep = degraded(0.5 * belt_spacing, kappa_t);
    let fold_keep = degraded(0.5 * spacing, kappa_t);
    let deep_keep = degraded(2.0 * POLAR_DEEP_CREST * POLAR_DEEP_ORDER * spacing, kappa_t);
    // A deep fold of the middling depth is as steep as the fabric's folds (F65).
    let deep_unit = POLAR_DEEP_RELIEF
        * POLAR_FOLD_AMPLITUDE
        * std::f64::consts::TAU
        * POLAR_DEEP_CREST
        * POLAR_DEEP_ORDER
        / DeepFolds::STEEPEST;
    let noise = GradientNoise::new(0x5027_F01D);
    let belt_noise = GradientNoise::new(0x5027_BE17);
    let rough_noise = GradientNoise::new(0x5027_0C4A);
    let rough_gain = 2f64.powf(-ROUGHNESS_HURST);
    let set_frequency = body_radius / GROOVE_SET_LENGTH.min(0.2 * body_radius);
    stamper.gather(grid, body_radius, pole, reach * body_radius);
    let gathered = clock.elapsed().as_secs_f64();
    let kernels = FoldKernels::new(
        &fabric,
        shell,
        &stamper.footprint,
        grid,
        body_radius,
        spacing,
        &noise,
        set_frequency,
        time_gyr.to_bits(),
    );
    let kerneled = clock.elapsed().as_secs_f64();
    // The deep folds (F65) grow from the fabric's deepest hinges through the same field:
    // the fissures' direction, the fabric's strength, a wavelength that goes as the lid.
    let edge_factor = (1.0 - shell.thinning * POLAR_CORE_PROFILE).max(0.05);
    let finest = GROOVE_MIN_CELLS * grid.delta() * body_radius;
    let sigma = POLAR_FOLD_ORIENTATION_SHARE * POLAR_FABRIC_RANGE / body_radius;
    let site = |p: [f64; 3]| {
        let dist = edge[grid.index(grid.locate(p).0)] as f64;
        let nearest = fabric.nearest(p);
        let from = |n: Option<([f64; 3], [f64; 3])>| {
            n.map_or(f64::INFINITY, |(f, _)| angle_between(f, p) * body_radius)
        };
        let near = {
            let x =
                ((POLAR_FABRIC_RANGE - from(nearest)) / (0.5 * POLAR_FABRIC_RANGE)).clamp(0.0, 1.0);
            x * x * (3.0 - 2.0 * x)
        };
        DeepSite {
            direction: fabric
                .orientation(p, sigma)
                .or(nearest.map(|(_, t)| t))
                .map(|t| tangential(p, t)),
            strength: fabric_strength(dist, halo) * near,
            wavelength: POLAR_DEEP_ORDER * (spacing * shell.factor(p) / edge_factor).max(finest),
            fissure: from(vents.nearest(p)),
        }
    };
    let nuclei = {
        let lattice = body_radius / (POLAR_DEEP_ORDER * spacing);
        let rng = Rng::new(time_gyr.to_bits() ^ 0x5027_DEE1);
        let mut keys = std::collections::HashSet::new();
        for &(k, _, _) in &stamper.footprint {
            if f64::from(edge[k]) > -halo {
                keys.insert(FoldKernels::key(grid.center(grid.cell(k)), lattice));
            }
        }
        let mut nuclei: Vec<DeepNucleus> = keys
            .into_iter()
            .map(|(i, j, l)| {
                let id = ((i as i64 & 0xF_FFFF) as u64)
                    | (((j as i64 & 0xF_FFFF) as u64) << 20)
                    | (((l as i64 & 0xF_FFFF) as u64) << 40);
                DeepNucleus {
                    key: (i, j, l),
                    center: normalize([
                        (i as f64 + rng.uniform(id, 0)) / lattice,
                        (j as f64 + rng.uniform(id, 1)) / lattice,
                        (l as f64 + rng.uniform(id, 2)) / lattice,
                    ]),
                    depth: fault_depth(rng.uniform(id, 3)),
                }
            })
            .collect();
        // In a fixed order, whatever order the lattice came in.
        nuclei.sort_by_key(|n| n.key);
        nuclei
    };
    let wander = GradientNoise::new(0x5027_DEE9 ^ time_gyr.to_bits());
    let laid = clock.elapsed().as_secs_f64();
    let deep = DeepFolds::grow(&nuclei, &site, Some(&wander), body_radius);
    let grown = clock.elapsed().as_secs_f64();
    if std::env::var_os("EXOGEN_DEBUG").is_some() {
        let ends = |e: ArmEnd| deep.arms.iter().filter(|a| a.end == e).count();
        eprintln!(
            "deep folds at {time_gyr:.4} Gyr: {} nuclei, {} arms, {} linked, {} free; gather {gathered:.1} s, kernels {:.1} s, nuclei {:.1} s, growth {:.1} s",
            nuclei.len(),
            deep.arms.len(),
            ends(ArmEnd::Linked),
            ends(ArmEnd::Free),
            kerneled - gathered,
            laid - kerneled,
            grown - laid
        );
    }
    if let Some(dir) = std::env::var_os("EXOGEN_DEEP_DUMP") {
        let arms: Vec<serde_json::Value> = deep
            .arms
            .iter()
            .map(|a| {
                let pts: Vec<[f64; 2]> = a
                    .points
                    .iter()
                    .map(|p| {
                        let (lon, lat) = crate::grid::lon_lat_from_direction(*p);
                        [lon.to_degrees(), lat.to_degrees()]
                    })
                    .collect();
                serde_json::json!({
                    "line": a.line,
                    "end": format!("{:?}", a.end),
                    "depth": a.depth,
                    "heights": a.heights,
                    "pts": pts,
                })
            })
            .collect();
        let path = std::path::Path::new(&dir).join(format!("deep_{time_gyr:.5}.json"));
        let _ = std::fs::write(path, serde_json::to_vec(&arms).unwrap_or_default());
    }
    let mut cells: Vec<CellWrite> = stamper
        .footprint
        .par_iter()
        .filter_map(|&(k, _, _)| {
            let dist = edge[k] as f64;
            if dist <= 1.0 - limit {
                return None;
            }
            // A cell something younger renewed keeps that.
            if renewed[k] as f64 > time_gyr + 1e-3 {
                return None;
            }
            let d = grid.center(grid.cell(k));
            // The margin ridge straddles the edge where the belt's arcs stand: where the
            // spreading meets the old ground head on it folds and piles, and between
            // the arcs the fabric's fractures run out into the old ground instead.
            let seg = belt_noise.fbm(d, body_radius / POLAR_BELT_SEGMENT, 2, 2.0, 0.5);
            let arc = {
                let t = (0.5 + POLAR_BELT_ARC_CONTRAST * seg).clamp(0.0, 1.0);
                t * t * (3.0 - 2.0 * t)
            };
            let x = dist / POLAR_MARGIN_WIDTH;
            let ridge = POLAR_MARGIN_HEIGHT * ridge_keep * arc * (-x * x).exp();
            let strength = fabric_strength(dist, halo);
            // The fabric of folds, a sum of oriented segments each laid along the
            // fissures' direction about its centre and in phase with the nearest one's
            // line, giving way in part to the deep folds where they stand (F65), and
            // under it the tectonic roughness of a plain deformed at every scale; none
            // past the halo.
            let (folds, rough) = if strength > 0.0 {
                let (deep_relief, hold) = deep.at(d);
                (
                    strength
                        * (POLAR_FOLD_AMPLITUDE
                            * fold_keep
                            * kernels.fabric(d)
                            * (1.0 - POLAR_DEEP_REPLACE * hold)
                            + deep_unit * deep_keep * deep_relief),
                    POLAR_TECTONIC_RMS
                        * fold_keep
                        * strength
                        * rough_noise.fbm(
                            d,
                            body_radius / POLAR_TECTONIC_LENGTH,
                            4,
                            2.0,
                            rough_gain,
                        )
                        / 0.4,
                )
            } else {
                (0.0, 0.0)
            };
            // The belt where the spreading is taken up, about the edge and decaying from
            // it: its folds parallel the edge in segments that shift against their
            // neighbours, en echelon, never as concentric rings; the fissures' folds give
            // way to it there.
            let inside = margin_coverage(dist);
            let belt_w = (-dist.abs() / POLAR_BELT_WIDTH).exp();
            let belt = POLAR_BELT_HEIGHT
                * belt_keep
                * belt_w
                * inside
                * arc
                * (2.0 * std::f64::consts::PI * dist / belt_spacing + POLAR_BELT_OFFSET * seg)
                    .sin();
            let relief =
                folds * (1.0 - belt_w * arc * (2.0 * inside).min(1.0)) + rough + belt + ridge;
            let weight = margin_coverage(dist);
            let h = bedrock[k] as f64;
            let floor = base[k] as f64 + fine.map_or(0.0, |f| f[k] as f64);
            (weight > 0.0 || relief.abs() >= 0.5)
                .then(|| write(k, h + (floor - h) * weight + relief))
        })
        .collect();
    cells.sort_by_key(|c| c.idx);
    if std::env::var_os("EXOGEN_DEBUG").is_some() {
        eprintln!(
            "polar terrain stamp: {} cells written, cells {:.1} s",
            cells.len(),
            clock.elapsed().as_secs_f64() - grown
        );
    }
    StampOutput {
        cells,
        composition: None,
        stain: None,
        steep: None,
        ejecta_age: None,
        compositions: None,
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

    /// A straight test lane along the equator, 2,000 km long and 200 km wide, with
    /// grooves 7 km apart, and a cell `across` and `along` it.
    fn equator_lane(stops: [Option<LaneStop>; 2]) -> (Lane, f64) {
        let r = 2.6e6;
        let points = Lane::great_circle([0.0, 0.0, 1.0], [1.0, 0.0, 0.0], 1000e3, 50e3, r);
        (
            Lane::structured(points, 50e3, 100e3, true, 1.0, 0.0, 1, 7e3, stops),
            r,
        )
    }

    fn on_equator_lane(r: f64, across: f64, along: f64) -> ([f64; 3], LaneSample) {
        let (b, c) = (along / r, across / r);
        (
            normalize([b.cos() * c.cos(), b.sin() * c.cos(), c.sin()]),
            LaneSample {
                idx: 0,
                across,
                along,
                tangent: [-(b.sin()), b.cos(), 0.0],
            },
        )
    }

    /// A lane meets the old ground across a border zone as wide as the lid it broke (F61):
    /// the same width on any grid that can draw it, and the grid's own ramp where it
    /// cannot.
    #[test]
    fn a_lane_meets_the_old_ground_across_its_lid() {
        let (lane, _) = equator_lane([None, None]);
        let lid = lane.lid();
        assert!((lid - 7e3 / (GROOVE_NECKING_WAVELENGTH * GROOVE_RIFT_LID_SHARE)).abs() < 1.0);
        assert!((lane.border(1e3) - GROOVE_BORDER_LIDS * lid).abs() < 1.0);
        assert!((lane.border(2e3) - lane.border(1e3)).abs() < 1.0);
        assert!((lane.border(20e3) - GROOVE_EDGE_CELLS * 20e3).abs() < 1.0);
        let edge = lane.edge_at(0.0, 1);
        let border = lane.border(1e3);
        assert!((lane.coverage_at(edge, 0.0, 1e3) - 0.5).abs() < 1e-6);
        assert!(lane.coverage_at(edge - 0.6 * border, 0.0, 1e3) > 0.999);
        assert!(lane.coverage_at(edge + 0.6 * border, 0.0, 1e3) < 1e-6);
        let quarter = lane.coverage_at(edge + 0.25 * border, 0.0, 1e3);
        assert!((quarter - 0.25).abs() < 1e-6, "{quarter}");
    }

    /// An end that stops is cut off along its stop's great circle and does not taper,
    /// and the lane leaves the ground beyond a cut inside an older lane alone; an end
    /// without a stop closes as a wedge.
    #[test]
    fn a_stopped_end_is_cut_along_its_stop() {
        let r: f64 = 2.6e6;
        // A cut square across the lane 800 km along it, its pole on the lane's side.
        let a = 800e3 / r;
        let stop = LaneStop {
            pole: [a.sin(), -(a.cos()), 0.0],
            reach: 600e3,
        };
        let (lane, r) = (equator_lane([None, Some(stop)]).0, r);
        let noise = LaneNoise::new(&lane);
        let cell = 4e3;
        let (d, sm) = on_equator_lane(r, 0.0, 700e3);
        assert!(lane.cover_at(&sm, d, r, cell) > 0.999);
        let (d, sm) = on_equator_lane(r, 0.0, 900e3);
        assert!(lane.cover_at(&sm, d, r, cell) < 1e-6);
        let f = lane.fabric_at(&sm, d, r, cell, 0.0, &noise);
        assert!(f.cover == 0.0 && f.outer == 0.0 && f.fracture == 0.0);
        // Far from the stop the cut does not reach.
        let (d, sm) = on_equator_lane(r, 0.0, -700e3);
        assert!(lane.cover_at(&sm, d, r, cell) > 0.999);
        assert!(lane.edge_at(790e3, 1) > 0.39 * lane.half_width);
        assert!(lane.half_width_at(-990e3) < 0.2 * lane.half_width);
    }

    /// A lane's grooves keep their relief whatever spacing the grid draws them at: a
    /// coarse grid draws them farther apart, not taller.
    #[test]
    fn grooves_keep_their_relief_on_any_grid() {
        let (lane, r) = equator_lane([None, None]);
        let noise = LaneNoise::new(&lane);
        assert!(lane.grooves_on(1e3) == 7e3 && lane.grooves_on(5e3) == 15e3);
        let rms = |cell: f64| {
            let (mut sum, mut n) = (0.0, 0usize);
            for i in 0..400 {
                for j in 0..40 {
                    let (d, sm) =
                        on_equator_lane(r, -40e3 + 2e3 * j as f64, -400e3 + 2e3 * i as f64);
                    let f = lane.fabric_at(&sm, d, r, cell, 0.0, &noise);
                    sum += f.relief * f.relief;
                    n += 1;
                }
            }
            (sum / n as f64).sqrt()
        };
        let (fine, coarse) = (rms(1e3), rms(5e3));
        let amplitude = GROOVE_RELIEF_SHARE * 7e3;
        assert!(fine > 0.2 * amplitude && fine < amplitude, "{fine}");
        assert!(
            coarse / fine > 0.6 && coarse / fine < 1.7,
            "coarse {coarse} against fine {fine}"
        );
    }

    /// F132: a groove is a chain of fault segments. Along one line it comes and goes, its
    /// segments stepped across within GROOVE_STEP of a spacing, so that where one ends the
    /// next takes up beside it; the troughs meet the plateau and each other with no seam;
    /// and the grooves die out toward the bounding faults.
    #[test]
    fn a_groove_comes_and_goes() {
        let (lane, r) = equator_lane([None, None]);
        let noise = LaneNoise::new(&lane);
        let spacing = lane.grooves_on(1e3);
        let n = 800;
        let (mut held, mut jump) = (0usize, 0.0f64);
        let mut steps: Vec<f64> = Vec::new();
        let mut last: Option<f64> = None;
        for i in 0..n {
            let (d, sm) = on_equator_lane(r, 0.0, -400e3 + 1e3 * i as f64);
            let set = lane.set_at(
                &sm,
                d,
                r,
                spacing,
                0.0,
                GROOVE_SEGMENT_LEVEL,
                0,
                0.0,
                &noise,
            );
            if set.gate > 0.1 {
                held += 1;
            }
            if set.gate > 0.5 {
                let x = if set.frac >= 0.25 {
                    set.frac - 0.75
                } else {
                    set.frac + 0.25
                };
                steps.push(set.turns - x - set.stripe - 0.75);
            }
            if let Some(prev) = last {
                jump = jump.max((set.shape - prev).abs());
            }
            last = Some(set.shape);
        }
        let share = held as f64 / n as f64;
        assert!(
            share > 0.3 && share < 0.97,
            "a segment holds {share:.2} of the line"
        );
        assert!(
            jump < 0.5,
            "a jump of {jump:.2} between neighbours a seventh of a spacing apart"
        );
        let (lo, hi) = steps
            .iter()
            .fold((f64::MAX, f64::MIN), |(a, b), &s| (a.min(s), b.max(s)));
        assert!(
            lo >= -GROOVE_STEP - 1e-9 && hi <= GROOVE_STEP + 1e-9 && hi - lo > 0.2,
            "steps from {lo:.2} to {hi:.2}"
        );
        // The grooves' relief half a spacing in from the bounding fault, against the middle.
        let rms = |across: f64| {
            let (mut sum, mut m) = (0.0, 0usize);
            for i in 0..400 {
                let (d, sm) = on_equator_lane(r, across, -400e3 + 2e3 * i as f64);
                let f = lane.fabric_at(&sm, d, r, 1e3, 0.0, &noise);
                if f.cover > 0.99 {
                    sum += f.relief * f.relief;
                    m += 1;
                }
            }
            (sum / m.max(1) as f64).sqrt()
        };
        let edge = lane.edge_at(0.0, 1);
        let (near, middle) = (rms(edge - 0.5 * spacing), rms(0.0));
        assert!(
            near < 0.7 * middle,
            "{near:.0} m by the fault against {middle:.0} m"
        );
    }

    /// Terrains are composed in the order they formed (F60): a younger one's ice lies
    /// over an older one's by its cover, the latest to take more than half a cell holds
    /// it, and one that only touches a cell takes its share of the old fine relief.
    #[test]
    fn a_younger_terrain_lies_over_an_older_one() {
        let mut r = Renewal::new(2);
        r.take(0, 1.0, 1.0, 0.6, 1.0, UNIT_GROOVED);
        r.take(0, 1.0, 1.0, 1.0, 2.0, UNIT_POLAR_TERRAIN);
        assert_eq!(r.time[0], 2.0);
        assert_eq!(r.unit[0], UNIT_POLAR_TERRAIN);
        assert!((r.clean[0] - 1.0).abs() < 1e-6);
        // A duller lane's edge over clean ice dulls it by its cover and holds nothing.
        r.take(0, 0.25, 0.25, 0.6, 3.0, UNIT_GROOVED);
        assert!((r.clean[0] - 0.9).abs() < 1e-6);
        assert_eq!(r.time[0], 2.0);
        assert!((r.keep[0] - 0.75).abs() < 1e-6);
        // A halo takes its weight of the old relief without laying ice.
        r.take(1, 0.0, 0.4, 1.0, 1.0, UNIT_POLAR_TERRAIN);
        assert_eq!(r.clean[1], 0.0);
        assert_eq!(r.time[1], f32::NEG_INFINITY);
        assert!((r.keep[1] - 0.6).abs() < 1e-6);
    }

    /// The finish completes a lane's renewal where the loop's cells did not reach (F60):
    /// the short relief gives way inside the lane, stays outside it, and stays under a
    /// crater younger than the lane.
    #[test]
    fn the_lanes_renewal_is_completed_on_the_resolved_grid() {
        let grid = CubeSphere::new(32);
        let n = grid.cell_count();
        let mut layers = LayerSet::new(grid);
        // A checkerboard of relief, a cell wide.
        for (k, h) in layers.bedrock.as_mut_slice().iter_mut().enumerate() {
            let c = grid.cell(k);
            *h = if (c.i + c.j).is_multiple_of(2) {
                100.0
            } else {
                -100.0
            };
        }
        let before = layers.bedrock.as_slice().to_vec();
        // A lane holds one face since 2 Gyr; one of its cells carries a crater of 3 Gyr.
        let cell_at = |face: u8, i: u32, j: u32| {
            (0..n)
                .find(|&k| {
                    let c = grid.cell(k);
                    c.face == face && c.i == i && c.j == j
                })
                .unwrap()
        };
        let mut renewal = Renewal::new(n);
        for k in 0..n {
            if grid.cell(k).face == 0 {
                renewal.lane_cover[k] = 1.0;
                renewal.lane_time[k] = 2.0;
                renewal.time[k] = 2.0;
            }
        }
        let young = cell_at(0, 16, 16);
        layers.unit_id.as_mut_slice()[young] = UNIT_CRATER_FLOOR;
        layers.surface_age.as_mut_slice()[young] = 3.0;
        lay_lane_floor(&mut layers, 252e3, None, &renewal, None);
        let h = layers.bedrock.as_slice();
        let inside = cell_at(0, 8, 8);
        let outside = cell_at(3, 16, 16);
        assert!(h[inside].abs() < 30.0, "{}", h[inside]);
        assert_eq!(h[outside], before[outside]);
        assert_eq!(h[young], before[young]);
        // Where the loop already took the relief, the finish takes no more.
        let mut again = LayerSet::new(grid);
        again.bedrock.as_mut_slice().copy_from_slice(&before);
        let taken = vec![1.0f32; n];
        lay_lane_floor(&mut again, 252e3, Some(&taken), &renewal, None);
        assert_eq!(again.bedrock.as_slice()[inside], before[inside]);
    }

    /// A lane's ends (F61, F142): a good share stop at a cut against an older lane that
    /// crosses them, the rest die out in the old ground, tapering to a point.
    #[test]
    fn lane_ends_stop_at_older_lanes_or_die_out() {
        let (g, _) = crate::presets::preset("Ganymede").unwrap().to_graph();
        let p = IceShellProcess::new(&g, CubeSphere::new(64)).unwrap();
        let ends = 2 * p.lanes.len();
        let cuts = p
            .lanes
            .iter()
            .flat_map(|l| l.stops.iter().flatten())
            .count();
        let share = cuts as f64 / ends.max(1) as f64;
        eprintln!(
            "ganymede: {} lanes, {ends} ends: {share:.2} stop against an older lane, the rest taper",
            p.lanes.len()
        );
        assert!(share > 0.25 && share < 0.9, "{share}");
        // A tapering end closes to nothing; a cut end keeps its width to the cut.
        let lane = p
            .lanes
            .iter()
            .find(|l| l.stops[1].is_none())
            .expect("a lane with a free end");
        assert!(lane.half_width_at(lane.half_length) < 0.05 * lane.half_width);
        assert!(lane.half_width_at(0.0) > 0.5 * lane.half_width);
    }

    /// A lane takes the lid where and when it opens (F61). On Enceladus the lid of the
    /// time is too thick for any lane, so every lane opens in a zone a focusing episode
    /// has thinned, during that episode, on a lid thinner than the mean; and none is
    /// wider than a lane can be.
    #[test]
    fn a_lane_opens_where_the_lid_lets_it() {
        // A less heated Enceladus, its eccentricity 0.8 of the real one's, has a lid too
        // thick for any lane on the lid of the time: its lanes open only where and when its
        // episodes thinned the shell, on the thinned lid.
        let (g, p) = enceladus_heated(0.8, 64);
        assert!(!p.lanes.is_empty());
        let widest = GROOVE_LANE_MAX_LENGTH_RAD * g.f64("body.radius") / GROOVE_LANE_ASPECT.0;
        let mean_lid = g.f64("interior.elastic_thickness");
        assert!(
            GROOVE_LANE_WIDTH_ELASTIC.0 * mean_lid > widest,
            "{mean_lid}"
        );
        let in_episode = |p: &IceShellProcess, t: f64| {
            p.episodes.iter().any(|e| t >= e.start && t <= e.end) || t >= p.present_start
        };
        for l in &p.lanes {
            assert!(2.0 * l.half_width <= widest * 1.0001, "{}", l.half_width);
            let t = l.time_gyr;
            assert!(
                in_episode(&p, t),
                "a lane at {t} Gyr, outside every episode"
            );
            assert!(l.lid() < 0.8 * mean_lid, "lid {}", l.lid());
            // The lane is as wide as its own lid makes it.
            let lids = 2.0 * l.half_width / l.lid();
            assert!(
                (GROOVE_LANE_WIDTH_ELASTIC.0 * 0.999..=GROOVE_LANE_WIDTH_ELASTIC.1 * 1.001)
                    .contains(&lids),
                "{lids}"
            );
        }
        // The real one's lid, over the shell its porous core keeps thin (F68), lets lanes
        // open on the lid of the time too, each as wide as its own lid makes it.
        let (g1, p1) = enceladus_heated(1.0, 64);
        let mean1 = g1.f64("interior.elastic_thickness");
        assert!(GROOVE_LANE_WIDTH_ELASTIC.0 * mean1 < widest, "{mean1}");
        assert!(p1.lanes.iter().any(|l| !in_episode(&p1, l.time_gyr)));
        for l in &p1.lanes {
            let lids = 2.0 * l.half_width / l.lid();
            assert!(
                (GROOVE_LANE_WIDTH_ELASTIC.0 * 0.999..=GROOVE_LANE_WIDTH_ELASTIC.1 * 1.001)
                    .contains(&lids),
                "{lids}"
            );
        }
        let span = |p: &IceShellProcess, of: &dyn Fn(&Lane) -> f64| {
            let v: Vec<f64> = p.lanes.iter().map(of).collect();
            (
                v.iter().cloned().fold(f64::INFINITY, f64::min),
                v.iter().cloned().fold(0.0, f64::max),
            )
        };
        for (name, p) in [("less heated", &p), ("enceladus", &p1)] {
            let (lid, width) = (
                span(p, &|l| l.lid() / 1e3),
                span(p, &|l| 2.0 * l.half_width / 1e3),
            );
            eprintln!(
                "{name}: {} lanes on lids of {:.1} to {:.1} km, {:.0} to {:.0} km wide",
                p.lanes.len(),
                lid.0,
                lid.1,
                width.0,
                width.1
            );
        }
    }

    /// Enceladus with its eccentricity scaled, and so its tides by the square: a body less
    /// or more heated than the real one, on the path between them.
    fn enceladus_heated(eccentricity: f64, face: u32) -> (ParamGraph, IceShellProcess) {
        let (mut g, _) = crate::presets::preset("Enceladus").unwrap().to_graph();
        let e = g.f64("orbit.eccentricity") * eccentricity;
        g.set("orbit.eccentricity", e).unwrap();
        let p = IceShellProcess::new(&g, CubeSphere::new(face)).unwrap();
        (g, p)
    }

    /// Between a body whose lid lets no lane form outside its thinned zones and one whose
    /// lid lets every lane form lie bodies that get some of each (principle 15): as a less
    /// heated Enceladus (its eccentricity 0.8 of the real one's) is made larger its lid
    /// takes up less of it, and lanes open on the lid of the time, outside any focusing
    /// episode, a few and then many; a fully heated one half again as large makes most of
    /// its lanes on the lid of the time and some in its episodes.
    #[test]
    fn lanes_come_in_by_degrees_as_the_lid_thins_against_the_body() {
        let mut outside_counts = Vec::new();
        let mut inside_counts = Vec::new();
        for (scale, heat) in [(1.0, 0.8), (1.2, 0.8), (1.5, 0.8), (1.5, 1.0)] {
            let (mut g, _) = crate::presets::preset("Enceladus").unwrap().to_graph();
            let e = g.f64("orbit.eccentricity") * heat;
            g.set("orbit.eccentricity", e).unwrap();
            let r = g.f64("body.radius") * scale;
            g.set("body.radius", r).unwrap();
            g.set("detail.quality", "preview").unwrap();
            let p = IceShellProcess::new(&g, CubeSphere::new(48)).unwrap();
            let in_episode = |t: f64| {
                p.episodes.iter().any(|e| t >= e.start && t <= e.end) || t >= p.present_start
            };
            let outside = p.lanes.iter().filter(|l| !in_episode(l.time_gyr)).count();
            let widest = GROOVE_LANE_MAX_LENGTH_RAD * r / GROOVE_LANE_ASPECT.0;
            for l in &p.lanes {
                assert!(2.0 * l.half_width <= widest * 1.0001);
            }
            eprintln!(
                "x{scale}: radius {:.0} km, shell {:.3} of it, lid {:.1} km, widest lane {:.0} km, share to lanes {:.3}, {} lanes, {} outside episodes",
                r / 1e3,
                g.f64("interior.ice_shell_thickness") / r,
                g.f64("interior.elastic_thickness") / 1e3,
                widest / 1e3,
                g.f64("tectonics.ice_resurfaced_fraction"),
                p.lanes.len(),
                outside,
            );
            // A body with no share for lanes has none.
            assert_eq!(
                p.lanes.is_empty(),
                g.f64("tectonics.ice_resurfaced_fraction") <= 0.0
            );
            outside_counts.push(outside);
            inside_counts.push(p.lanes.len() - outside);
        }
        assert_eq!(outside_counts[0], 0);
        assert!(inside_counts[0] > 0);
        // A few, then more than twice as many (12 and 34 with the ends tapering, F142).
        assert!(outside_counts[1] > 0 && outside_counts[2] > 2 * outside_counts[1]);
        // The fully heated one is a world in between: lanes of both kinds.
        assert!(
            inside_counts[3] > 0 && outside_counts[3] > inside_counts[3],
            "{inside_counts:?}"
        );
    }

    /// A groove's depth has a heavy tail (F63): the middling one at one, half of them a
    /// little shallower, an eighth much deeper, none past the bound.
    #[test]
    fn groove_depths_have_a_heavy_tail() {
        assert!((fault_depth(0.5) - 1.0).abs() < 1e-9);
        assert!(fault_depth(1.0) > 0.7 && fault_depth(1.0) < 1.0);
        // One groove in eight is more than 1.7 times as deep as the middling one.
        assert!(fault_depth(0.125) > 1.7 && fault_depth(0.125) < 2.0);
        assert!(
            fault_depth(1e-9) <= GROOVE_DEPTH_MAX && fault_depth(1e-9) > 0.95 * GROOVE_DEPTH_MAX
        );
        let mut last = f64::INFINITY;
        for i in 1..=100 {
            let v = fault_depth(i as f64 / 100.0);
            assert!(v <= last);
            last = v;
        }
    }

    /// A later pulse reworks part of a lane (F63): its grooves cut the first set where
    /// they stand and leave it elsewhere.
    #[test]
    fn a_later_pulse_cuts_the_first_set_in_patches() {
        let (plain, r) = equator_lane([None, None]);
        let mut worked = plain.clone();
        worked.turns = vec![0.3, -0.25];
        let noise = LaneNoise::new(&plain);
        let (mut same, mut cut, mut n) = (0usize, 0usize, 0usize);
        for i in 0..300 {
            for j in 0..30 {
                let (d, sm) = on_equator_lane(r, -30e3 + 2e3 * j as f64, -300e3 + 2e3 * i as f64);
                let a = plain.fabric_at(&sm, d, r, 1e3, 0.0, &noise).relief;
                let b = worked.fabric_at(&sm, d, r, 1e3, 0.0, &noise).relief;
                n += 1;
                if (a - b).abs() < 1.0 {
                    same += 1;
                } else if (a - b).abs() > 20.0 {
                    cut += 1;
                }
            }
        }
        let (same, cut) = (same as f64 / n as f64, cut as f64 / n as f64);
        assert!(same > 0.15, "untouched share {same}");
        assert!(cut > 0.1, "reworked share {cut}");
    }

    /// Under loose cover a fracture beyond a lane's fault shows as a chain of pits (F64):
    /// pits along its line where its segment stands, none where the cover is too thin
    /// for the grid to draw them, and the open fracture itself where there is no cover.
    #[test]
    fn a_covered_fracture_shows_as_pits() {
        let (lane, r) = equator_lane([None, None]);
        let noise = LaneNoise::new(&lane);
        let width = 1000.0;
        let along_line = |gate: f64| -> Vec<f64> {
            (0..4000)
                .map(|i| lane.pits_at(0.0, 50.0 * i as f64, 3.0, width, gate, &noise))
                .collect()
        };
        let share = |v: &[f64]| v.iter().filter(|&&x| x > 0.0).count() as f64 / v.len() as f64;
        // A pit every pit and a half: most of the line is pit, with ground between.
        let whole = along_line(1.0);
        let pitted = share(&whole);
        assert!(pitted > 0.45 && pitted < 0.8, "{pitted}");
        assert!(whole.iter().cloned().fold(0.0, f64::max) > 0.9);
        assert!(whole.contains(&0.0));
        // The chain breaks as the fracture does, and nothing stands off its line.
        let broken = share(&along_line(0.3));
        assert!(
            broken > 0.05 && broken < 0.6 * pitted,
            "{broken} against {pitted}"
        );
        assert_eq!(
            lane.pits_at(1.1 * width, 12_345.0, 3.0, width, 1.0, &noise),
            0.0
        );
        // In the damage zone, thin cover fills its share of the open fractures and makes
        // no pits a grid of 250 m can draw: the walls show less by that share.
        let edge = lane.edge_at(0.0, 1);
        // Over a band of lines, since the fractures beyond the fault are segments that
        // come and go and step across (F132).
        let shown = |loose: f64| -> f64 {
            (0..2000 * 11)
                .map(|i| {
                    let across = edge + 5.9e3 + 100.0 * (i % 11) as f64;
                    let (d, sm) = on_equator_lane(r, across, -200e3 + 200.0 * (i / 11) as f64);
                    lane.fabric_at(&sm, d, r, 250.0, loose, &noise).fracture
                })
                .sum()
        };
        let (open, thin) = (shown(0.0), shown(20.0));
        assert!(open > 0.0);
        let fracture_depth = GROOVE_FRACTURE_SHARE * GROOVE_RELIEF_SHARE * 7e3;
        assert!(
            (thin / open - (1.0 - 20.0 / fracture_depth)).abs() < 1e-6,
            "{}",
            thin / open
        );
        // Thick cover fills them, and pits the grid can draw stand in their place.
        let thick = shown(400.0);
        assert!(thick > 0.0 && thick < open, "{thick} against {open}");
    }

    /// A later pulse's grooves are struck along the field as it stands later in the
    /// lane's working life (F63): on Ganymede the strike field of extension drifts, and
    /// the later sets are turned by degrees to tens of degrees; on Enceladus the tidal
    /// field barely turns within an episode, and they run nearly parallel to the first.
    #[test]
    fn later_pulses_take_the_field_as_it_turns() {
        let turns_of = |name: &str| -> Vec<f64> {
            let (g, _) = crate::presets::preset(name).unwrap().to_graph();
            let p = IceShellProcess::new(&g, CubeSphere::new(64)).unwrap();
            assert!(p.lanes.iter().all(|l| l.turns.len() == GROOVE_PULSES - 1));
            let mut v: Vec<f64> = p
                .lanes
                .iter()
                .flat_map(|l| l.turns.iter().map(|t| t.abs().to_degrees()))
                .collect();
            v.sort_by(f64::total_cmp);
            v
        };
        let ganymede = turns_of("Ganymede");
        // A less heated Enceladus's lanes all open in its episodes, within which its tidal
        // field barely turns.
        let enceladus: Vec<f64> = {
            let (_, p) = enceladus_heated(0.8, 64);
            assert!(p.lanes.iter().all(|l| l.turns.len() == GROOVE_PULSES - 1));
            let mut v: Vec<f64> = p
                .lanes
                .iter()
                .flat_map(|l| l.turns.iter().map(|t| t.abs().to_degrees()))
                .collect();
            v.sort_by(f64::total_cmp);
            v
        };
        let at = |v: &[f64], q: f64| v[((v.len() - 1) as f64 * q) as usize];
        eprintln!(
            "later pulses' turns, degrees: ganymede median {:.1}, nine tenths under {:.1}; enceladus median {:.2}, largest {:.2}",
            at(&ganymede, 0.5),
            at(&ganymede, 0.9),
            at(&enceladus, 0.5),
            at(&enceladus, 1.0)
        );
        assert!((1.0..45.0).contains(&at(&ganymede, 0.5)));
        assert!(at(&enceladus, 1.0) < 5.0);
    }

    /// A point on a flat patch about [1, 0, 0], km east and north of it, and back.
    fn patch_point(east_km: f64, north_km: f64, r: f64) -> [f64; 3] {
        normalize([1.0, east_km * 1e3 / r, north_km * 1e3 / r])
    }

    fn patch_coords(p: [f64; 3], r: f64) -> (f64, f64) {
        (p[1] / p[0] * r / 1e3, p[2] / p[0] * r / 1e3)
    }

    /// Two deep folds 30 km apart on one strike, the second offset across it by a share
    /// of the wavelength, grown in a uniform field with no wander.
    fn two_hinges(offset: f64) -> (DeepFolds, f64) {
        let r = 2.52e5;
        let wavelength = 7.5e3;
        let field = move |p: [f64; 3]| {
            let (x, y) = patch_coords(p, r);
            DeepSite {
                direction: Some(tangential(p, [0.0, 1.0, 0.0])),
                strength: if x.abs() < 80.0 && y.abs() < 40.0 {
                    1.0
                } else {
                    0.0
                },
                wavelength,
                fissure: f64::INFINITY,
            }
        };
        let nuclei = [
            DeepNucleus {
                key: (0, 0, 0),
                center: patch_point(-15.0, 0.0, r),
                depth: 2.0,
            },
            DeepNucleus {
                key: (1, 0, 0),
                center: patch_point(15.0, offset * wavelength / 1e3, r),
                depth: 2.0,
            },
        ];
        (DeepFolds::grow(&nuclei, &field, None, r), r)
    }

    #[test]
    fn a_renewed_terrain_takes_nothing_beyond_its_ramp() {
        // F66: no halo is drawn from the outline; the ground beside a terrain relaxes by
        // its own law under the thinned shell instead.
        assert_eq!(margin_coverage(-0.5 * POLAR_MARGIN_RAMP), 0.0);
        assert!((margin_coverage(0.0) - 0.5).abs() < 1e-12);
        assert_eq!(margin_coverage(0.5 * POLAR_MARGIN_RAMP), 1.0);
        // A cell beyond the ramp keeps its old relief whole; one on the ramp gives up part.
        let mut renewal = Renewal::new(2);
        renewal.take_polar(&[-20e3, -0.25 * POLAR_MARGIN_RAMP as f32], 1.0);
        assert_eq!(renewal.keep[0], 1.0);
        assert!(renewal.keep[1] < 1.0);
    }

    #[test]
    fn hinges_link_within_half_a_wavelength() {
        // A third of a wavelength apart, the two link into one fold, lower where they
        // joined than where either nucleated (Grasemann and Schmalholz 2012).
        let (near, _) = two_hinges(0.3);
        let linked: Vec<&DeepArm> = near
            .arms
            .iter()
            .filter(|a| a.end == ArmEnd::Linked)
            .collect();
        assert!(!linked.is_empty());
        for a in &linked {
            let joint = *a.heights.last().unwrap();
            assert!(joint < a.depth && joint > 0.5 * a.depth, "{joint}");
        }
        // Four fifths of a wavelength apart, they pass each other en echelon and end
        // free, plunging to nothing.
        let (far, r) = two_hinges(0.8);
        assert!(far.arms.iter().all(|a| a.end == ArmEnd::Free));
        let (tip, _) = patch_coords(*far.arms[0].points.last().unwrap(), r);
        assert!(
            tip > 15.0,
            "the first fold's tip, at {tip} km, passes the second's nucleus"
        );
        for a in &far.arms {
            assert_eq!(*a.heights.last().unwrap(), 0.0);
        }
    }

    #[test]
    fn a_deep_fold_is_a_crest_between_two_troughs() {
        // Across the first fold 15 km west of its nucleus, where the second has no arm.
        let (far, r) = two_hinges(0.8);
        let at = |north_km: f64| far.at(patch_point(-30.0, north_km, r));
        let (crest, hold) = at(0.0);
        assert!(crest > 1.0, "{crest}");
        assert!((hold - 1.0).abs() < 1e-6, "{hold}");
        // The troughs either side take up what the crest stands above.
        let out = 3f64.sqrt() * DeepFolds::width(7.5e3) / 1e3;
        for side in [-out, out] {
            let ratio = at(side).0 / crest;
            assert!((ratio + 2.0 * (-1.5f64).exp()).abs() < 0.02, "{ratio}");
        }
        assert!(at(-30.0).0.abs() < 1e-3);
    }

    #[test]
    fn deep_folds_join_into_a_network() {
        let r = 2.52e5;
        let wavelength = 7.5e3;
        // A field whose lines converge toward the middle as they run east.
        let field = move |p: [f64; 3]| {
            let (x, y) = patch_coords(p, r);
            DeepSite {
                direction: Some(tangential(p, [0.0, 1.0, -0.3 * y / 45.0])),
                strength: if x.abs() < 70.0 && y.abs() < 45.0 {
                    1.0
                } else {
                    0.0
                },
                wavelength,
                fissure: f64::INFINITY,
            }
        };
        let rng = Rng::new(11);
        let mut nuclei = Vec::new();
        for i in -7i32..=7 {
            for j in -5i32..=5 {
                let id = ((i + 64) * 256 + (j + 64)) as u64;
                if rng.uniform(id, 3) > 0.4 {
                    continue;
                }
                let x = i as f64 * 9.0 + 6.0 * (rng.uniform(id, 0) - 0.5);
                let y = j as f64 * 9.0 + 6.0 * (rng.uniform(id, 1) - 0.5);
                nuclei.push(DeepNucleus {
                    key: (i, j, 0),
                    center: patch_point(x, y, r),
                    depth: fault_depth(rng.uniform(id, 2)),
                });
            }
        }
        let noise = GradientNoise::new(5);
        let folds = DeepFolds::grow(&nuclei, &field, Some(&noise), r);
        // Of the arms that ended well inside the field, most linked into a neighbour.
        let inside = |p: [f64; 3]| {
            let (x, y) = patch_coords(p, r);
            x.abs() < 62.5 && y.abs() < 37.5
        };
        let ended: Vec<&DeepArm> = folds
            .arms
            .iter()
            .filter(|a| inside(*a.points.last().unwrap()))
            .collect();
        let linked = ended.iter().filter(|a| a.end == ArmEnd::Linked).count();
        eprintln!(
            "{} nuclei, {} arms, {} ended inside, {linked} linked",
            nuclei.len(),
            folds.arms.len(),
            ended.len()
        );
        assert!(linked * 2 >= ended.len(), "{linked} of {}", ended.len());
        // No two folds cross: where they meet, they join.
        type Piece = (usize, (f64, f64), (f64, f64));
        let pieces: Vec<Piece> = folds
            .arms
            .iter()
            .flat_map(|a| {
                a.points
                    .windows(2)
                    .map(move |w| (a.line, patch_coords(w[0], r), patch_coords(w[1], r)))
            })
            .collect();
        let orient = |p: (f64, f64), q: (f64, f64), s: (f64, f64)| {
            (q.0 - p.0) * (s.1 - p.1) - (q.1 - p.1) * (s.0 - p.0)
        };
        let mut crossings = 0;
        for (i, a) in pieces.iter().enumerate() {
            for b in &pieces[i + 1..] {
                if a.0 == b.0 {
                    continue;
                }
                let d = [
                    orient(b.1, b.2, a.1),
                    orient(b.1, b.2, a.2),
                    orient(a.1, a.2, b.1),
                    orient(a.1, a.2, b.2),
                ];
                if d[0] * d[1] < 0.0 && d[2] * d[3] < 0.0 && d.iter().all(|v| v.abs() > 1e-9) {
                    crossings += 1;
                }
            }
        }
        assert_eq!(crossings, 0);
    }
}
