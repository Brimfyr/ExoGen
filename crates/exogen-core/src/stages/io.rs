// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! The Io regime (PLAN.md 4.4): a rocky body heated past about half a watt per square
//! metre.
//!
//! Resurfacing outruns impacts by orders of magnitude, so the surface is always young
//! and carries no craters. What it carries is in a steady state, and the loop lays it
//! down afresh at every step on the body's own base terrain, so every epoch has an Io of
//! its own. The volcanic centres are paterae, volcano-tectonic depressions with irregular
//! fault-bounded outlines, a share of them active with floors of fresh silicate; the
//! active ones and a few quiet ones feed compound flow fields, lobes of different ages;
//! and every centre throws plumes from time to time, the gas-rich Pele-class ones far
//! and red with sulfur, the Prometheus-class ones nearer and white with sulfur dioxide
//! frost. Mountains are tilted crustal blocks raised by the subsidence of the crust
//! under the burial, standing a few hundred million years.
//!
//! The colour is not painted: the surface is what the plumes have left on it (F84).
//! Every cell carries a top layer that the fallout buries as it falls and that ages
//! between falls: red sulfur reverts to yellow, frost sublimes as fast as the
//! temperature of its latitude allows, and where the atmosphere is too thin to shield
//! the ground the plasma alters what lies on top to red-brown. A lava flow or a patera
//! floor is black when it is laid and is coated by whatever falls on it after. The
//! history leaves the dark silicate share in `deep`, the red sulfur in `melt`, the
//! altered share in `stain` and the visible frost in `cover`, which the materials stage
//! lights. Between the stamps and the history, a lava cell carries its coverage in
//! `cover` and the years since it was laid in `melt`.

use crate::constants::*;
use crate::crater::{CellWrite, StampOutput, Stamper};
use crate::fields::{Field, LayerSet};
use crate::geom::{angle_between, cross, dot, local_basis, normalize};
use crate::grid::CubeSphere;
use crate::grid::projection::{PaddedField, sample_bilinear};
use crate::noise::GradientNoise;
use crate::params::ParamGraph;
use crate::rng::Rng;
use crate::units::*;
use rayon::prelude::*;

pub const VERSION: u64 = 24;

/// A patera: a volcano-tectonic depression with a flat floor, its outline the collapses
/// that made it and the faults that cut it (F88).
#[derive(Clone, Debug)]
pub struct Patera {
    pub center: [f64; 3],
    /// Its first collapse's radius, the depth of its floor below the rim, m.
    pub radius: f64,
    pub depth: f64,
    /// A lava lake on the floor at present.
    pub active: bool,
    /// Years since the floor last held a lake: it is coated by what has fallen on it
    /// since (F84).
    pub floor_age_yr: f64,
    /// Gas-rich, Pele-class, whose plumes are hot and red; the rest are Prometheus-class.
    pub pele: bool,
    /// The azimuth of its long axis, rad from east towards north at its centre, along
    /// the crust's fabric there, and its collapses' long axes over their short (F90).
    pub axis: f64,
    pub ellipticity: f64,
    /// Its collapses, each an ellipse about a point offset from its centre (east and
    /// north, m) with the radius of the circle of its area, m; and the faults that cut
    /// straight sides into it, each the azimuth of its normal and its distance from the
    /// centre, m.
    pub collapses: Vec<[f64; 3]>,
    pub faults: Vec<[f64; 2]>,
    /// The mountain it sank beside, whose foot is its margin there (F90).
    pub abuts: Option<Mountain>,
    pub birth_gyr: f64,
    pub death_gyr: f64,
    pub id: u64,
}

/// A mountain (F87): a crustal block thrust up along a fault across the crust's
/// compression.
#[derive(Clone, Debug)]
pub struct Mountain {
    pub center: [f64; 3],
    /// The fault's strike, rad from east towards north at the centre; the block's length
    /// along it and its width across it, and its height, m.
    pub strike: f64,
    pub length: f64,
    pub width: f64,
    pub height: f64,
    /// The side of the fault the block stands on, +1 or -1, and the share of its height
    /// its back slope loses across it.
    pub vergence: f64,
    pub tilt: f64,
    /// Its form: how level its top is (a mesa's or a plateau's against a massif's or a
    /// ridge's), and how ridged along its strike (a lineated mountain's).
    pub level: f64,
    pub ridging: f64,
    pub birth_gyr: f64,
    pub death_gyr: f64,
    pub id: u64,
}

/// One lobe of a flow field: a tube about a centre line traced downslope from the
/// patera's rim (F59), laid down some years ago.
#[derive(Clone, Debug)]
pub struct Lobe {
    /// The centre line from the rim outward, unit directions, `spacing` m apart.
    pub points: Vec<[f64; 3]>,
    pub spacing: f64,
    /// The mean half width, m, and its factor at each point.
    pub half_width: f64,
    pub widths: Vec<f32>,
    pub age_yr: f64,
}

/// An eruption of a patera's lava (F117): a flood like any other body's, from a vent on
/// the patera's rim where the ground falls away, this volume, this many years before the
/// time the surface is drawn at.
#[derive(Clone, Copy, Debug)]
pub struct IoEruption {
    pub vent: [f64; 3],
    pub volume: f64,
    pub age_yr: f64,
    pub id: u64,
    /// How fast its lava poured out while it erupted, m³/s (F117, F201).
    pub effusion: f64,
    /// The years its field grew over, its last lava `age_yr` years ago (F201).
    pub span_yr: f64,
}

/// The orders of the harmonics of a plume deposit's range with azimuth.
const RANGE_ORDERS: (usize, usize) = (2, 6);

/// One eruption of a plume over the history the surface remembers (F84): a vent
/// throwing its fallout to the ballistic range of its eruption speed, from `start_yr`
/// to `end_yr` years before the present, `end_yr` 0 while it erupts.
#[derive(Clone, Debug)]
pub struct PlumeEpisode {
    pub center: [f64; 3],
    /// The ballistic range v² / g, m.
    pub radius: f64,
    pub pele: bool,
    pub start_yr: f64,
    pub end_yr: f64,
    /// The fallout at the ring's crest, optical depths a year.
    pub crest: f64,
    /// The deposit's elongation across an elongated vent: its aspect, and the azimuth of
    /// its long axis in the vent's local frame.
    pub aspect: f64,
    pub axis: f64,
    /// The silicate share of the inner fall, and the frost share of the ring.
    pub dark: f64,
    pub frost: f64,
    /// The variation of the range with azimuth, as harmonics: each an amplitude times
    /// the cosine and sine of its phase.
    range: [[f64; 2]; RANGE_ORDERS.1 - RANGE_ORDERS.0 + 1],
    /// Where its mottling lies in the noise.
    mottle: [f64; 3],
    /// The azimuths its pyroclasts jet along from the cracks, rad (F208); none for an
    /// even fall.
    jets: Vec<f64>,
    /// One of the small eruptions that lay the background fallout (F89).
    pub small: bool,
    pub id: u64,
}

impl Patera {
    /// Its floor's area, m²: its collapses' areas, overlapping where they coalesced.
    pub fn floor_area(&self) -> f64 {
        std::f64::consts::PI * self.collapses.iter().map(|c| c[2] * c[2]).sum::<f64>()
    }
    /// The farthest its outline reaches from its centre, m.
    pub fn extent(&self) -> f64 {
        let long = self.ellipticity.sqrt();
        self.collapses
            .iter()
            .map(|c| c[0].hypot(c[1]) + c[2] * long)
            .fold(self.radius * long, f64::max)
            * (1.0 + 4.0 * IO_PATERA_WANDER)
    }
    pub fn alive(&self, t: f64) -> bool {
        self.birth_gyr - IO_ALIVE_TOLERANCE_GYR <= t && t < self.death_gyr - IO_ALIVE_TOLERANCE_GYR
    }
}
impl Mountain {
    /// The farthest the block and its slopes reach from its centre, m.
    pub fn extent(&self) -> f64 {
        let scarp = mountain_scarp(self.height, f64::MIN_POSITIVE);
        (MOUNTAIN_FAULT_SHARE * self.length).hypot(MOUNTAIN_BACK_REACH * self.width)
            + scarp * (1.0 + 3.0 * IO_MOUNTAIN_APRON_REACH)
            + mountain_heave(self.height)
    }
    /// How far a drawing of the block on a grid of `cell` reaches, its scarp's debris
    /// with it, m.
    pub fn reach(&self, cell: f64) -> f64 {
        self.extent() + 3.0 * IO_MOUNTAIN_APRON_REACH * mountain_scarp(self.height, cell)
    }
    /// A point's place in the block's frame: along its fault from its middle, and across
    /// it towards the block, m.
    pub fn frame(&self, r: f64, phi: f64) -> (f64, f64) {
        let (x, y) = (r * phi.cos(), r * phi.sin());
        let (c, s) = (self.strike.cos(), self.strike.sin());
        (x * c + y * s, self.vergence * (y * c - x * s))
    }
    /// The block's height in a direction on the sphere, m, as a grid of `cell` draws it.
    pub fn height_at(&self, d: [f64; 3], body_radius: f64, cell: f64) -> f64 {
        let (east, north) = local_basis(self.center);
        let r = angle_between(self.center, d) * body_radius;
        let phi = crate::geom::azimuth_from(self.center, east, north, d);
        let (u, v) = self.frame(r, phi);
        mountain_height(self, u, v, cell)
    }
    pub fn alive(&self, t: f64) -> bool {
        self.birth_gyr - IO_ALIVE_TOLERANCE_GYR <= t && t < self.death_gyr - IO_ALIVE_TOLERANCE_GYR
    }
}
/// An eruption as the history reads it at every cell: its reach, its oval and range
/// harmonics, and its fallout's crest, profile and make-up, from `start` to `end` years
/// from the present (negative, within the history).
struct Plume {
    cos_cap: f64,
    reach: f64,
    /// The ring's width beyond its crest, as a share of the range.
    outer: f64,
    jets: Vec<f64>,
    radius: f64,
    axis: (f64, f64),
    aspect: f64,
    range: [[f64; 2]; RANGE_ORDERS.1 - RANGE_ORDERS.0 + 1],
    crest: f64,
    inner: f64,
    frost: f64,
    red: f64,
    dark: f64,
    mottle_frequency: f64,
    mottle: [f64; 3],
    start: f64,
    end: f64,
}

impl Plume {
    fn new(e: &PlumeEpisode, body_radius: f64, window: f64) -> Self {
        let (red, inner, outer) = if e.pele {
            (1.0, IO_RING_INNER_PELE, IO_RING_OUTER_WIDTH_PELE)
        } else {
            (
                IO_PROMETHEUS_RED,
                IO_RING_INNER_PROMETHEUS,
                IO_RING_OUTER_WIDTH,
            )
        };
        let (start, end) = (-e.start_yr.min(window), -e.end_yr);
        // The deposit reaches as far as its ring's tail lays a fall that could show over
        // the eruption, as a share of its range, at most one and a half ranges; a small
        // eruption's, a disc with no ring, to its range (F89).
        let total = e.crest * (end - start).max(0.0) / NEGLIGIBLE_FALL;
        let reach = if e.small {
            1.0
        } else {
            (1.0 + outer * total.max(1.0).ln().sqrt()).min(1.5)
        };
        let wobble: f64 = e.range.iter().map(|h| h[0].hypot(h[1])).sum();
        let cap = reach * e.radius * e.aspect.sqrt() * (1.0 + wobble) / body_radius;
        Self {
            cos_cap: cap.min(std::f64::consts::PI).cos(),
            reach,
            outer,
            jets: e.jets.clone(),
            radius: e.radius,
            axis: (e.axis.cos(), e.axis.sin()),
            aspect: e.aspect,
            range: e.range,
            crest: e.crest,
            inner,
            frost: e.frost,
            red,
            dark: e.dark,
            // The fallout's mass varies on a fifth of the range, unresolved.
            mottle_frequency: body_radius / (0.15 * e.radius),
            mottle: e.mottle,
            start,
            end,
        }
    }

    /// The share of the pyroclasts at an azimuth (F208): IO_PELE_JET_SHARE along the
    /// cracks' jets, the rest even, the mean of all azimuths one; even where there are no
    /// jets.
    fn jet_at(&self, azimuth: f64) -> f64 {
        if self.jets.is_empty() {
            return 1.0;
        }
        let tau = std::f64::consts::TAU;
        let w = IO_PELE_JET_WIDTH;
        let sum: f64 = self
            .jets
            .iter()
            .map(|&j| {
                let d = (azimuth - j + 0.5 * tau).rem_euclid(tau) - 0.5 * tau;
                (-(d / w).powi(2)).exp()
            })
            .sum();
        let jetted = sum * tau / (self.jets.len() as f64 * w * std::f64::consts::PI.sqrt());
        1.0 - IO_PELE_JET_SHARE + IO_PELE_JET_SHARE * jetted
    }

    /// The deposit's range at an azimuth, as a share of the ballistic range, from the
    /// cosines and sines of the azimuth's multiples (index k for order k): the vent's
    /// oval and the scatter of the eruption speed with azimuth.
    fn range_at(&self, cos_k: &[f64], sin_k: &[f64]) -> f64 {
        let (ca, sa) = self.axis;
        let c = cos_k[1] * ca + sin_k[1] * sa;
        let s = sin_k[1] * ca - cos_k[1] * sa;
        let oval = 1.0 / (c * c / self.aspect + self.aspect * s * s).sqrt();
        let mut wobble = 1.0;
        for (j, h) in self.range.iter().enumerate() {
            let k = RANGE_ORDERS.0 + j;
            wobble += h[0] * cos_k[k] - h[1] * sin_k[k];
        }
        oval * wobble.max(0.3)
    }
}

/// Where a ray at an azimuth from a patera's centre leaves the last of its collapses
/// (F88, F90): ellipses along its axis, overlapping, so that their union is scalloped
/// where they meet, as nested and coalesced calderas are.
fn outline_radius(phi: f64, axis: f64, ellipticity: f64, collapses: &[[f64; 3]]) -> f64 {
    let (ca, sa) = (axis.cos(), axis.sin());
    // The ray along the axis and across it.
    let (du, dv) = ((phi - axis).cos(), (phi - axis).sin());
    let long = ellipticity.sqrt();
    let mut edge: f64 = 0.0;
    for c in collapses {
        let (cu, cv) = (c[0] * ca + c[1] * sa, c[1] * ca - c[0] * sa);
        let (a, b) = (c[2] * long, c[2] / long);
        let qa = (du / a).powi(2) + (dv / b).powi(2);
        let qb = du * cu / (a * a) + dv * cv / (b * b);
        let qc = (cu / a).powi(2) + (cv / b).powi(2) - 1.0;
        let disc = qb * qb - qa * qc;
        if disc >= 0.0 {
            edge = edge.max((qb + disc.sqrt()) / qa);
        }
    }
    edge
}

/// How far a patera's margin wanders off its collapses' arcs at an azimuth, as a factor
/// of their radius (F118): at every scale from a quarter of its outline down.
fn margin_wander(p: &Patera, phi: f64) -> f64 {
    let wander = 0.7 * crate::crater::stamp::crenulation(p.id ^ 0xA7E2, phi)
        + 0.3 * crate::crater::stamp::lobes(p.id ^ 0xA7E3, phi, 40.0);
    1.0 + IO_PATERA_WANDER * wander
}

/// Radius of a patera's outline at an azimuth (F88): where a ray from its centre leaves
/// the last of its collapses, unless a fault that cut a straight side into it comes
/// first; the faults' sides end at sharp corners.
pub fn patera_radius(phi: f64, p: &Patera) -> f64 {
    let (dx, dy) = (phi.cos(), phi.sin());
    let mut edge = outline_radius(phi, p.axis, p.ellipticity, &p.collapses) * margin_wander(p, phi);
    for f in &p.faults {
        let n = dx * f[0].cos() + dy * f[0].sin();
        if n > 1e-9 {
            edge = edge.min(f[1] / n);
        }
    }
    edge
}

/// The years since a patera floor's lava was laid at a point (F117): a lake's crust in
/// patches from its last overturn to IO_FLOW_ACTIVE_YR old, a quiet floor's patches up to
/// twice older or younger than its age; the plumes' fallout coats each as old as it is.
fn floor_age(p: &Patera, d: [f64; 3], patches: &GradientNoise, body_radius: f64) -> f64 {
    let scale = body_radius / (IO_FLOOR_PATCH_SHARE * p.radius).max(1.0);
    let n = (0.5 * (1.0 + patches.fbm(d, scale, 3, 2.0, 0.5))).clamp(0.0, 1.0);
    if p.active {
        let (lo, hi) = IO_FLOW_ACTIVE_YR;
        lo * (hi / lo).powf(n)
    } else {
        p.floor_age_yr * 2f64.powf(2.0 * n - 1.0)
    }
}

/// Depth of a flat-floored pit at a distance from its centre, with walls of a width.
pub fn pit_profile(r: f64, radius: f64, depth: f64, wall: f64) -> f64 {
    depth * 0.5 * (1.0 - ((r - radius) / wall).tanh())
}

/// Azimuths at which a patera beside a mountain finds the mountain's foot (F90).
const FOOT_AZIMUTHS: usize = 360;

/// Where a patera that sank beside a mountain meets it (F90): along each of
/// FOOT_AZIMUTHS rays from its centre, the distance to where the block, as a grid of
/// `cell` draws it, first stands above its share of its height, found in steps of `step`
/// and refined between them; infinite where the ray leaves the outline first. Paterae
/// that abut mountains have arcuate or scalloped margins, and some wrap round them
/// (Radebaugh et al. 2001: Gish Bar, Monan): the margin is the block's foot, bowed and
/// cut back as the foot is.
fn foot_edges(p: &Patera, m: &Mountain, body_radius: f64, cell: f64, step: f64) -> Vec<f64> {
    let (east, north) = local_basis(p.center);
    let top = MOUNTAIN_UNIT_SHARE * m.height;
    let rises = |t: f64, phi: f64| {
        let d = crate::geom::offset(p.center, east, north, t / body_radius, phi);
        m.height_at(d, body_radius, cell) > top
    };
    (0..FOOT_AZIMUTHS)
        .map(|k| {
            let phi = std::f64::consts::TAU * k as f64 / FOOT_AZIMUTHS as f64;
            let reach = patera_radius(phi, p);
            let mut t = 0.0;
            while t <= reach {
                if rises(t, phi) {
                    let (mut lo, mut hi) = ((t - step).max(0.0), t);
                    for _ in 0..4 {
                        let mid = 0.5 * (lo + hi);
                        if rises(mid, phi) {
                            hi = mid;
                        } else {
                            lo = mid;
                        }
                    }
                    return hi;
                }
                t += step;
            }
            f64::INFINITY
        })
        .collect()
}

/// The margin a mountain gives a patera at an azimuth, between the rays of its table.
fn foot_edge(edges: &[f64], phi: f64) -> f64 {
    let n = edges.len();
    let x = phi.rem_euclid(std::f64::consts::TAU) / std::f64::consts::TAU * n as f64;
    let i = (x.floor() as usize) % n;
    let j = (i + 1) % n;
    let f = x - x.floor();
    match (edges[i].is_finite(), edges[j].is_finite()) {
        (true, true) => edges[i] + (edges[j] - edges[i]) * f,
        (true, false) if f < 0.5 => edges[i],
        (false, true) if f >= 0.5 => edges[j],
        _ => f64::INFINITY,
    }
}

/// A patera as a grid draws it (F88, F90): its wall's half width, and where the
/// mountain it sank beside, drawn as a grid of `mountain_cell` draws it, gives it its
/// margin.
struct PateraDrawing<'a> {
    patera: &'a Patera,
    wall: f64,
    mountain_cell: f64,
    foot: Option<Vec<f64>>,
}

impl<'a> PateraDrawing<'a> {
    fn new(patera: &'a Patera, body_radius: f64, cell: f64, mountain_cell: f64) -> Self {
        Self {
            patera,
            wall: patera_wall(patera.depth, cell),
            mountain_cell,
            foot: patera
                .abuts
                .as_ref()
                .map(|m| foot_edges(patera, m, body_radius, mountain_cell, 0.5 * cell)),
        }
    }

    /// The depth the patera takes from a cell at a distance and azimuth from its centre,
    /// in direction `d`, m, and whether the cell is its floor. The mountain's debris
    /// below its scarp sank with the floor, or lies under its lava, so the floor is level.
    fn cut(&self, r: f64, phi: f64, d: [f64; 3], body_radius: f64) -> (f64, bool) {
        let p = self.patera;
        let edge = self.edge(phi);
        let pit = pit_profile(r, edge, p.depth, self.wall);
        let debris = match &p.abuts {
            Some(m) if pit > 0.01 * p.depth => {
                let h = m.height_at(d, body_radius, self.mountain_cell);
                h.min(MOUNTAIN_UNIT_SHARE * m.height) * pit / p.depth
            }
            _ => 0.0,
        };
        (pit + debris, r < edge)
    }

    fn needs_direction(&self) -> bool {
        self.foot.is_some()
    }

    /// The floor's edge at azimuth `phi`, m from the centre: its outline, against the
    /// mountain's foot where it sank beside one.
    fn edge(&self, phi: f64) -> f64 {
        let edge = patera_radius(phi, self.patera);
        match &self.foot {
            Some(foot) => edge.min(foot_edge(foot, phi)),
            None => edge,
        }
    }
}

/// The half width of a patera's wall for its depth on a grid of `cell`, m (F88): steep
/// as Chaac Patera's, 2.7 km high at 70 degrees (Radebaugh et al. 2001), at least a cell.
pub fn patera_wall(depth: f64, cell: f64) -> f64 {
    (0.5 * depth / IO_PATERA_WALL_SLOPE_DEG.to_radians().tan()).max(cell)
}

/// The share of a block's height above which its ground is the block's own, not its
/// debris or the plains' embayments.
const MOUNTAIN_UNIT_SHARE: f64 = 0.15;

/// The fault's half length as a share of the block's exposed length: the throw falls as
/// a cosine squared to the tips, and the block is exposed where it stands above the
/// plains' burial.
const MOUNTAIN_FAULT_SHARE: f64 = 0.7;

/// How far across a block, as a share of its width from its scarp, it and its relief
/// reach.
const MOUNTAIN_BACK_REACH: f64 = 1.6;

/// The width of a block's scarp, worn back by landslides to its slope, for a height, m,
/// at least two cells, the narrowest the grid draws without stepping.
fn mountain_scarp(height: f64, cell: f64) -> f64 {
    (height / IO_MOUNTAIN_SCARP_SLOPE_DEG.to_radians().tan()).max(2.0 * cell)
}

/// How far a block's thrust has carried it out over the plains where its throw is
/// greatest, m: the heave of a throw of its height on its fault.
fn mountain_heave(height: f64) -> f64 {
    height / IO_MOUNTAIN_FAULT_DIP_DEG.to_radians().tan()
}

fn smoothstep(x: f64) -> f64 {
    let x = x.clamp(0.0, 1.0);
    x * x * (3.0 - 2.0 * x)
}

/// The spread of one octave of the gradient noise over a plane through it, and how far
/// apart its zero crossings lie along a line, in its lattice's cells (measured).
const NOISE_SPREAD: f64 = 0.275;
const NOISE_CROSSING_SPACING: f64 = 1.26;

/// How many times longer along its strike a lineated block's ridges run than they lie
/// apart, before they bend, branch or end (a stated choice).
const RIDGE_STRETCH: f64 = 8.0;

/// The noise a block's unresolved relief is drawn from, offset for each block.
fn relief_noise() -> &'static GradientNoise {
    static NOISE: std::sync::OnceLock<GradientNoise> = std::sync::OnceLock::new();
    NOISE.get_or_init(|| GradientNoise::new(0x10_B10C_0087))
}

/// A point of a plane through the noise, `a` and `b` along it and `depth` across it, the
/// plane slanted to the noise's lattice so that the lattice's axes leave no grain.
fn slanted(a: f64, b: f64, depth: f64) -> [f64; 3] {
    [
        0.866 * a - 0.5 * b,
        0.383 * a + 0.663 * b - 0.643 * depth,
        0.321 * a + 0.557 * b + 0.766 * depth,
    ]
}

/// The height of a mountain at a point of its frame (F87), `u` along its fault from its
/// middle and `v` across it towards the block, m, on a grid of `cell`. The block is the
/// crust thrust up along its fault: its uplift is the fault's throw, falling to its tips
/// as a fault's slip does, times its fall across the block to its down-dip edge, a
/// plateau's that holds level, a tilted block's that falls steadily, or a massif's that
/// rounds; the thrust has carried it out over the plains by the throw's heave, so that
/// it bows out towards them. Its surface carries its own relief, unresolved, the more on
/// its dissected margins, and the plains bury it to a share of its height, so its back
/// and its tips are where its surface dips under them. Its front is the fault's scarp,
/// cut back in arcs by landslides and worn to its slope, their debris spread below. A
/// lineated block is ridged along its strike where the grid resolves it.
pub fn mountain_height(m: &Mountain, u: f64, v: f64, cell: f64) -> f64 {
    // The fault runs past the block's exposed length, its throw dying out under the
    // plains at its tips.
    let a = MOUNTAIN_FAULT_SHARE * m.length;
    let s = u / a;
    if s.abs() >= 1.0 {
        return 0.0;
    }
    let throw = (0.5 * std::f64::consts::PI * s).cos().powi(2);
    // The block as it stood before the thrust carried it out by its heave, the farthest
    // where its throw is greatest: its trace bows out towards the plains, as a thrust's
    // does where its displacement is greatest (the bow-and-arrow rule).
    let v = v + mountain_heave(m.height) * throw;
    let scarp = mountain_scarp(m.height, cell);
    let rng = Rng::new(m.id ^ 0x5CA1_10B5);
    // The scarp's line: arcs where landslides have cut back into the block, meeting in
    // rounded cusps, shallower towards the tips with the throw.
    let count = ((m.length / IO_MOUNTAIN_SCALLOP_LENGTH).round() as u64).max(1);
    let soft = 0.3 * scarp;
    let mut cut = 0.0;
    for k in 0..count {
        let centre = -a + 2.0 * a * (k as f64 + 0.25 + 0.5 * rng.uniform(k, 1)) / count as f64;
        let half = a / count as f64 * (0.7 + 0.8 * rng.uniform(k, 2));
        let depth = IO_MOUNTAIN_SCALLOP_DEPTH * m.width * rng.uniform(k, 3);
        let x = (u - centre) / half;
        cut += (depth * (1.0 - x * x).max(0.0) / soft).exp() - 1.0;
    }
    let front = soft * cut.ln_1p() * throw;
    // The surface across the block, from its fault's line; the landslides took away what
    // stood in front of the scarp.
    let noise = relief_noise();
    let o = (m.id % 4096) as f64 * 3.7;
    // The block's own relief, of unit spread: octaves from half its width down to what
    // the grid resolves, each IO_MOUNTAIN_RELIEF_GAIN of the one before.
    let mut octaves = 2;
    while octaves < 5 && 0.5 * m.width / 2f64.powi(octaves) >= 4.0 * cell {
        octaves += 1;
    }
    let spread = NOISE_SPREAD
        * (0..octaves)
            .map(|k| IO_MOUNTAIN_RELIEF_GAIN.powi(2 * k))
            .sum::<f64>()
            .sqrt();
    let relief_at = |u: f64, v: f64| {
        let (mut sum, mut amplitude, mut f) = (0.0, 1.0, 2.0 / m.width);
        for k in 0..octaves {
            let shift = o + 17.0 * k as f64;
            sum += amplitude
                * noise.sample(slanted(
                    u * f + shift,
                    v * f + 0.5 * shift,
                    0.37 + 1.3 * k as f64,
                ));
            amplitude *= IO_MOUNTAIN_RELIEF_GAIN;
            f *= 2.0;
        }
        sum / spread
    };
    // A plateau's back, an erosional scarp that wanders along the strike.
    let f = 1.0 / m.width;
    let wander = (noise.sample(slanted(u * f + o, 0.0, 7.7))
        + 0.5 * noise.sample(slanted(2.0 * u * f, o, 9.0)))
        / (NOISE_SPREAD * 1.25f64.sqrt());
    let back = 1.0 + 0.13 * wander;
    // A lineated block's ridges along its strike, where the grid resolves them: the folds
    // and slices of the thrust sheet, crests that wander, branch and end, the noise's zero
    // crossings drawn out along the strike, in patches over the block; its surface's own
    // relief, buried with it.
    let ridged = m.ridging > 0.0 && IO_MOUNTAIN_RIDGE_SPACING >= 3.0 * cell;
    let ridge = |v: f64| {
        if !ridged {
            return 0.0;
        }
        let across = NOISE_CROSSING_SPACING / IO_MOUNTAIN_RIDGE_SPACING;
        let n = noise.sample(slanted(u * across / RIDGE_STRETCH + o, v * across, 5.3));
        let crest = (1.0 - (n / NOISE_SPREAD).powi(2)).max(0.0).powi(2);
        let patches =
            smoothstep(0.5 + 0.5 * noise.sample(slanted(u * f, v * f + o, 9.1)) / NOISE_SPREAD);
        m.ridging * IO_MOUNTAIN_RIDGE_RELIEF * patches * crest
    };
    let surface = |v: f64| {
        let x = v / m.width;
        let plateau = smoothstep((back - x) / 0.35);
        let tilted = (1.0 - x).max(0.0);
        let round = (1.0 - x * x).max(0.0);
        let tilt = (m.tilt / IO_MOUNTAIN_TILT).clamp(0.0, 1.0);
        let fall = m.level * plateau + (1.0 - m.level) * (tilt * tilted + (1.0 - tilt) * round);
        let lift = throw * fall;
        // The relief is the raised crust's, dissected the more where it stands low, and
        // none where the crust was not raised.
        let amplitude = (IO_MOUNTAIN_RELIEF * (1.0 - lift) + IO_MOUNTAIN_TOP_RELIEF * lift)
            * smoothstep(lift / IO_MOUNTAIN_BURIAL);
        let rough = amplitude * relief_at(u, v) + throw * ridge(v);
        ((lift + rough - IO_MOUNTAIN_BURIAL) / (1.0 - IO_MOUNTAIN_BURIAL)).max(0.0)
    };
    if v < front {
        // The scarp, down to the top of its debris at its foot, and the debris below,
        // thinning out to nothing within its reach.
        let crest = m.height * surface(front);
        let toe = IO_MOUNTAIN_APRON * crest;
        let foot = front - scarp;
        return if v > foot {
            toe + (crest - toe) * smoothstep((v - foot) / scarp)
        } else {
            let reach = IO_MOUNTAIN_APRON_REACH * scarp;
            let d = (foot - v) / reach;
            toe * (-d).exp() * smoothstep((3.0 - d) / 1.0)
        };
    }
    let x = v / m.width;
    if x >= MOUNTAIN_BACK_REACH {
        return 0.0;
    }
    // Everything the block holds fades out before the end of its reach.
    let fade = smoothstep((MOUNTAIN_BACK_REACH - x) / 0.4);
    m.height * surface(v) * fade
}

/// The crust's fabric at a point (F87): the strike a thrust takes there, about 30 or 60
/// degrees from the tide's most compressive direction at its most compressive phase
/// (Seeger and de Kleer 2025: Io's mountains form most often so offset from it), which
/// of the four the draw picks; without a tide, any.
fn fabric_strike(
    tide: Option<&crate::stages::tidal::TidalShell>,
    p: [f64; 3],
    pick: f64,
    free: f64,
) -> f64 {
    let Some(shell) = tide else {
        return std::f64::consts::TAU * free;
    };
    let axis = [1.0, 0.0, 0.0];
    let mut best = (f64::INFINITY, 0.0);
    for k in 0..TIDAL_PHASE_SAMPLES {
        let phase = std::f64::consts::TAU * k as f64 / TIDAL_PHASE_SAMPLES as f64;
        let (_, s2, psi) = shell.diurnal(axis, p, phase).principal();
        if s2 < best.0 {
            best = (s2, psi + std::f64::consts::FRAC_PI_2);
        }
    }
    let offsets = [
        std::f64::consts::FRAC_PI_6,
        -std::f64::consts::FRAC_PI_6,
        std::f64::consts::FRAC_PI_3,
        -std::f64::consts::FRAC_PI_3,
    ];
    best.1 + offsets[((pick * 4.0) as usize).min(3)]
}

/// A fall too thin to show over its whole eruption, optical depths.
const NEGLIGIBLE_FALL: f64 = 0.02;

/// Bins of the years a flow's lava was laid, as shares of the history, for the plumes at
/// its fronts (F202): finer towards the present, where the frost they lay still shows.
const FRONT_AGE_BINS: [f64; 8] = [0.0, 0.002, 0.006, 0.02, 0.06, 0.2, 0.6, 1.0];

/// The length over which the plumes at the flows' fronts share their eruption speed, m
/// (F202): a plume keeps its form as its source wanders with its front (Kieffer et al.
/// 2000), and a field's pieces erupt alike, its neighbours' differently.
const FRONT_SPEED_LENGTH: f64 = 600e3;

/// How far back the history reaches, years (F84): a few burial times of the background
/// fallout, past which nothing a plume laid can still show.
pub fn history_years() -> f64 {
    IO_HISTORY_BURIALS / IO_BACKGROUND_FALLOUT
}

/// Fallout of a plume at a fraction `x` of its range (F84), as the crest's share: the
/// ring where the gas condenses at the canopy and falls back, ending beyond the range of
/// the fastest over `outer` of the range (F208) and rising gently within; and, as a share
/// of `inner`, the pyroclasts' thinner fall within it.
pub fn deposit_profile(x: f64, outer: f64) -> (f64, f64) {
    let width = if x < 1.0 { IO_RING_WIDTH } else { outer };
    let ring = (-((x - 1.0) / width).powi(2)).exp();
    let fill = ((1.0 - x) / 0.5).clamp(0.0, 1.0);
    (ring, fill)
}

/// A small eruption's fall at a fraction `x` of its range (F89), as its centre's share:
/// thinning smoothly to nothing at the range.
pub fn small_profile(x: f64) -> f64 {
    (1.0 - x * x).max(0.0).powi(2)
}

/// The area under a deposit's profile, as the integral of the profile times x over x: a
/// deposit of range R spreads its mass over 2π R² times this.
pub fn profile_moment(inner: f64, outer: f64) -> f64 {
    let n = 4000;
    let dx = 1.5 / n as f64;
    (0..n)
        .map(|i| {
            let x = (i as f64 + 0.5) * dx;
            let (ring, fill) = deposit_profile(x, outer);
            (ring + inner * fill) * x * dx
        })
        .sum()
}

/// Where a Pele-class eruption's gas rises from its lava lake (F208): `count` cracks
/// along the margin, drawn by length, each facing out across the margin there; the
/// share of the gas facing each azimuth, each crack's spread over the angle its gas
/// covers, as the harmonics of orders RANGE_ORDERS of its departure from even; and the
/// cracks' azimuths, along which the pyroclasts jet.
fn margin_facing(p: &Patera, count: usize, rng: &Rng, draw: u64) -> MarginFacing {
    let n = 360;
    let tau = std::f64::consts::TAU;
    let points: Vec<(f64, f64)> = (0..n)
        .map(|i| {
            let phi = tau * i as f64 / n as f64;
            let r = patera_radius(phi, p);
            (r * phi.cos(), r * phi.sin())
        })
        .collect();
    let pieces: Vec<(f64, f64)> = (0..n)
        .map(|i| {
            let (x0, y0) = points[i];
            let (x1, y1) = points[(i + 1) % n];
            let (dx, dy) = (x1 - x0, y1 - y0);
            (dx.hypot(dy), (-dx).atan2(dy))
        })
        .collect();
    let total: f64 = pieces
        .iter()
        .map(|p| p.0)
        .sum::<f64>()
        .max(f64::MIN_POSITIVE);
    let jets: Vec<f64> = (0..count.max(1) as u64)
        .map(|c| {
            let mut at = rng.uniform(draw, c) * total;
            let mut psi = pieces[0].1;
            for &(ds, normal) in &pieces {
                psi = normal;
                if at < ds {
                    break;
                }
                at -= ds;
            }
            psi
        })
        .collect();
    let mut range = [[0.0; 2]; RANGE_ORDERS.1 - RANGE_ORDERS.0 + 1];
    for (j, h) in range.iter_mut().enumerate() {
        let k = (RANGE_ORDERS.0 + j) as f64;
        let (a, b) = jets.iter().fold((0.0, 0.0), |(a, b), &psi| {
            (a + (k * psi).cos(), b + (k * psi).sin())
        });
        let spread = (-0.5 * (k * IO_PELE_GAS_SPREAD).powi(2)).exp();
        let scale = 2.0 / jets.len() as f64 * spread * IO_PELE_MARGIN_SHAPE;
        // The range harmonics' convention: h[0] cos kφ − h[1] sin kφ.
        *h = [scale * a, -scale * b];
    }
    MarginFacing { range, jets }
}

struct MarginFacing {
    range: [[f64; 2]; RANGE_ORDERS.1 - RANGE_ORDERS.0 + 1],
    jets: Vec<f64>,
}

/// The plumes of the Pele-class paterae over the history the surface remembers (F84;
/// the Prometheus class rise at the flows' fronts, `front_falls`, F202). Each patera's
/// eruptions are drawn from the present back: one erupting now if the patera is active
/// and its class's share of active centres throws a plume, part way through; then quiet
/// spells and eruptions alternating at the duty its class keeps, the spells memoryless
/// and the eruptions from a few months to decades long. Each eruption has its own speed
/// and so its range and its own mass; its gas rises from the cracks along its lava
/// lake's margin, so that its ring reaches farther where more of the margin faces and
/// its pyroclasts jet from the cracks (F208).
pub fn plume_history(paterae: &[Patera], gravity: f64, seed: u64) -> Vec<PlumeEpisode> {
    let rng = Rng::new(seed ^ 0x0E91_5DE5_0000_0084);
    let (d_lo, d_hi) = IO_PLUME_DURATION_YR;
    let mean_duration = (d_hi - d_lo) / (d_hi / d_lo).ln();
    let moment = profile_moment(IO_RING_INNER_PELE, IO_RING_OUTER_WIDTH_PELE);
    let tau = 2.0 * std::f64::consts::PI;
    let window = history_years();
    let mut out = Vec::new();
    for p in paterae.iter().filter(|p| p.pele) {
        let idx = p.id;
        let share = IO_PLUME_FRACTION_PELE;
        let duty = IO_PATERA_ACTIVE_FRACTION * share;
        let quiet = mean_duration * (1.0 - duty) / duty;
        let duration = |draw: u64| d_lo * (d_hi / d_lo).powf(rng.uniform(idx, draw));
        let mut spans: Vec<(f64, f64)> = Vec::new();
        let mut t = 0.0;
        if p.active && rng.uniform(idx, 0) < share {
            let elapsed = duration(1) * rng.uniform(idx, 2);
            spans.push((elapsed, 0.0));
            t = elapsed;
        }
        let mut draw = 3u64;
        while spans.len() < IO_PLUME_EPISODES_MAX {
            t -= quiet * (1.0 - rng.uniform(idx, draw)).ln();
            if t >= window {
                break;
            }
            let d = duration(draw + 1);
            spans.push((t + d, t));
            t += d;
            draw += 2;
        }
        for (j, &(start, end)) in spans.iter().enumerate() {
            let base = 1_000 + 64 * j as u64;
            let speed =
                IO_PLUME_SPEED_PELE * (IO_PLUME_SPEED_SCATTER * rng.normal(idx, base)).exp();
            let radius = speed * speed / gravity.max(0.01);
            let mass =
                IO_PLUME_MASS_PELE * (IO_PLUME_MASS_SCATTER * rng.normal(idx, base + 2)).exp();
            let crest = mass * YEAR / (tau * radius * radius * moment * IO_DEPOSIT_OPTICAL_MASS);
            // Each eruption cracks the lake's margin afresh.
            let (j_lo, j_hi) = IO_PELE_JETS;
            let count = (j_lo * (j_hi / j_lo).powf(rng.uniform(idx, base + 7))).round() as usize;
            let facing = margin_facing(p, count, &rng, idx ^ (base << 20));
            out.push(PlumeEpisode {
                center: p.center,
                radius,
                pele: true,
                start_yr: start,
                end_yr: end,
                crest,
                aspect: 1.0,
                axis: 0.0,
                dark: IO_DARK_PELE,
                frost: IO_PELE_FROST,
                range: facing.range,
                mottle: [
                    1e3 * rng.uniform(idx, base + 40),
                    1e3 * rng.uniform(idx, base + 41),
                    1e3 * rng.uniform(idx, base + 42),
                ],
                jets: facing.jets,
                small: false,
                id: idx ^ (base << 32),
            });
        }
    }
    out
}

/// The weight of a latitude in the regime's heat flow, against the mean over the
/// sphere: the paterae's, and the small eruptions'.
fn heat_weight(sin_lat: f64) -> f64 {
    let w = IO_PATERA_POLAR_WEIGHT + (1.0 - IO_PATERA_POLAR_WEIGHT) * (1.0 - sin_lat * sin_lat);
    w / (IO_PATERA_POLAR_WEIGHT + (1.0 - IO_PATERA_POLAR_WEIGHT) * 2.0 / 3.0)
}

/// Io's small eruptions over the history the surface remembers (F89): the vents too
/// small to be paterae and the stealth plumes, which carry much of Io's heat unseen
/// (Veeder et al. 2012: the known volcanoes 60 percent; de Pater et al. 2020: SO
/// emission spread widely, mostly away from known volcanoes). Their ranges follow a power
/// law, as many per factor of range as cover equal areas, each rising too low for a
/// canopy, so that its frost and sulfur, some with silicate, fall in a smooth bump about
/// the vent, at a place the heat flow favours and a time within the history, at a rate
/// of the order of the surface changes Galileo saw; those under two cells of `cell`
/// across are left out, their fall a small share of the background's.
pub fn small_eruptions(seed: u64, cell: f64) -> Vec<PlumeEpisode> {
    let (r_min, r_max) = IO_SMALL_ERUPTION_RANGE;
    let b = IO_SMALL_ERUPTION_SLOPE;
    let r_lo = (2.0 * cell).clamp(r_min, r_max);
    let window = history_years();
    let count = (IO_SMALL_ERUPTIONS_PER_YR * window * (r_min / r_lo).powf(b)).round() as u64;
    let mut out = Vec::new();
    if r_lo >= r_max {
        return out;
    }
    let tau = 2.0 * std::f64::consts::PI;
    let rng = Rng::new(seed ^ 0x5_3A11_E2B7_0089);
    let (d_lo, d_hi) = IO_SMALL_ERUPTION_DURATION_YR;
    let (depth_med, depth_spread) = IO_SMALL_ERUPTION_DEPTH;
    let (frost_lo, frost_hi) = IO_PROMETHEUS_FROST;
    let (lo_b, hi_b) = (r_lo.powf(-b), r_max.powf(-b));
    for idx in 0..count {
        // Where the heat flow favours: draw until the latitude weight accepts.
        let mut center = rng.unit_direction(idx, 0);
        for attempt in 0..12u64 {
            let w = heat_weight(center[2]) / heat_weight(0.0);
            if rng.uniform(idx ^ crate::rng::ACCEPT_STREAM, attempt) < w {
                break;
            }
            center = rng.unit_direction(idx ^ crate::rng::ACCEPT_STREAM, 100 + 2 * attempt);
        }
        // Draws 0 and 1 place it.
        let radius = (lo_b - rng.uniform(idx, 30) * (lo_b - hi_b)).powf(-1.0 / b);
        let duration = d_lo * (d_hi / d_lo).powf(rng.uniform(idx, 2));
        let depth = depth_med * (depth_spread * rng.normal(idx, 3)).exp();
        let frost = frost_lo + (frost_hi - frost_lo) * rng.uniform(idx, 5);
        // Within the history, some still erupting.
        let t = (window + duration) * rng.uniform(idx, 6) - duration;
        let (start, end) = ((t + duration).min(window), t.max(0.0));
        let crest = depth / duration;
        if start <= end {
            continue;
        }
        out.push(PlumeEpisode {
            center,
            radius,
            pele: false,
            start_yr: start,
            end_yr: end,
            crest,
            aspect: 1.0 + IO_VENT_ELONGATION * rng.uniform(idx, 7),
            axis: tau * rng.uniform(idx, 8),
            dark: if rng.uniform(idx, 9) < IO_DARK_CONTAMINATED {
                IO_DARK_PROMETHEUS
            } else {
                0.0
            },
            frost,
            range: [[0.0; 2]; RANGE_ORDERS.1 - RANGE_ORDERS.0 + 1],
            mottle: [
                1e3 * rng.uniform(idx, 20),
                1e3 * rng.uniform(idx, 21),
                1e3 * rng.uniform(idx, 22),
            ],
            jets: Vec::new(),
            small: true,
            id: 0x5_0000_0000_0000 + idx,
        });
    }
    out
}

/// Io's background fallout at every cell (F89), optical depths a year: the far fallout
/// of the plumes, spread about each volcanic centre alive at the epoch over
/// IO_FAR_FALLOUT_REACH of its class's range and weighted by its plumes' output, the
/// share of the time it throws one times their mass, so that the plains are buried
/// fastest among the active centres and slowest far from them; scaled to the background's
/// mean over the surface, its mass clumped by a factor whose logarithm spreads by
/// IO_FAR_FALLOUT_SPREAD from IO_FAR_FALLOUT_MOTTLE down, unresolved.
pub fn background_fall(
    grid: CubeSphere,
    body_radius: f64,
    gravity: f64,
    paterae: &[Patera],
    seed: u64,
) -> Vec<f32> {
    let uniform = vec![IO_BACKGROUND_FALLOUT as f32; grid.cell_count()];
    if paterae.is_empty() {
        return uniform;
    }
    // The far fallout varies over hundreds of kilometres: drawn on a coarse grid.
    let far = far_fallout(
        CubeSphere::new(grid.face_size().min(64)),
        body_radius,
        gravity,
        paterae,
    );
    let mean = far.area_mean();
    if mean <= 0.0 {
        return uniform;
    }
    let far = far.resample(grid);
    // Its clumpiness, unresolved: a log-normal factor of the noise at unit spread, from
    // the mottling's length down by three octaves.
    let noise = GradientNoise::new(seed ^ 0xFA2_FA11_0089);
    let frequency = body_radius / IO_FAR_FALLOUT_MOTTLE;
    let clumps = Field::from_fn(grid, |_, d| noise.fbm(d, frequency, 4, 2.0, 0.6) as f32);
    let values = clumps.as_slice();
    let n = values.len() as f64;
    let mean = values.iter().map(|&v| v as f64).sum::<f64>() / n;
    let spread = (values
        .iter()
        .map(|&v| (v as f64 - mean).powi(2))
        .sum::<f64>()
        / n)
        .sqrt()
        .max(1e-9);
    let mottled = Field::from_fn(grid, |c, _| {
        let k = grid.index(c);
        let unit = (values[k] as f64 - mean) / spread;
        far.as_slice()[k] * (IO_FAR_FALLOUT_SPREAD * unit).exp() as f32
    });
    let scale = IO_BACKGROUND_FALLOUT / mottled.area_mean();
    mottled
        .as_slice()
        .iter()
        .map(|&v| (v as f64 * scale) as f32)
        .collect()
}

/// The fall of the plumes at the flows' fronts over one bin of years (F202), optical
/// depths a year at each cell, from `start` to `end` years from the present (negative):
/// the ring's condensate, the pyroclasts within it, and their silicate.
struct FrontFall {
    start: f64,
    end: f64,
    ring: Vec<f32>,
    fill: Vec<f32>,
    dark: Vec<f32>,
    /// The plumes' gas settling across their regions (F205).
    gas: Vec<f32>,
}

/// The plumes at the young flows' fronts (F202). A Prometheus-class plume rises where a
/// flow runs over the plains' SO2, not at its vent (Kieffer et al. 2000: Prometheus's
/// plume wandered 75 to 95 km west with its flow field over twenty years while its vent
/// stayed, its form unchanged; Milazzo et al. 2001), so its source is the lava laid in
/// each bin of years, IO_FRONT_SO2 of SO2 for each cubic metre, thrown at the class's
/// speed with its scatter and falling in the class's ring about each piece of the
/// source: the deposits wrap the flows and trail their growth, as Io's white diffuse
/// deposits surround its flows in irregular units (Williams et al. 2011). The rings are
/// broad and drawn on a coarse grid.
fn front_falls(
    grid: CubeSphere,
    body_radius: f64,
    gravity: f64,
    units: &[u16],
    coverage: &[f32],
    laid_yr: &[f32],
    window: f64,
) -> Vec<FrontFall> {
    let coarse = CubeSphere::new(grid.face_size().min(128));
    let range = IO_PLUME_SPEED_PROMETHEUS * IO_PLUME_SPEED_PROMETHEUS / gravity.max(0.01);
    let moment = profile_moment(IO_RING_INNER_PROMETHEUS, IO_RING_OUTER_WIDTH);
    // Each piece of a source erupts at the speed its place draws, with the class's
    // scatter, a unit spread of the noise; the range goes as the speed squared.
    let speeds = GradientNoise::new(0xF207_5EED);
    let frequency = body_radius / FRONT_SPEED_LENGTH;
    let radius_at = |d: [f64; 3]| {
        let z = (2.5 * speeds.fbm(d, frequency, 2, 2.0, 0.5)).clamp(-2.0, 2.0);
        range * (2.0 * IO_PLUME_SPEED_SCATTER * z).exp()
    };
    let farthest = 1.4;
    let bins: Vec<(f64, f64)> = FRONT_AGE_BINS
        .windows(2)
        .map(|b| (b[0] * window, b[1] * window))
        .collect();
    // Each bin's lava, m², in each coarse cell, and the share of each cell it covers.
    let mut lava = vec![vec![0.0f64; coarse.cell_count()]; bins.len()];
    let mut laid: Vec<Field<f32>> = bins.iter().map(|_| Field::<f32>::new(grid)).collect();
    for k in 0..units.len() {
        if units[k] != UNIT_IO_FLOW {
            continue;
        }
        let age = laid_yr[k] as f64;
        let Some(b) = bins.iter().position(|&(lo, hi)| age >= lo && age < hi) else {
            continue;
        };
        let c = grid.cell(k);
        let share = coverage[k].clamp(0.0, 1.0);
        let area = grid.cell_area(c) * body_radius * body_radius;
        lava[b][coarse.index(coarse.locate(grid.center(c)).0)] += area * share as f64;
        laid[b].as_mut_slice()[k] = share;
    }
    let tau = std::f64::consts::TAU;
    let dark_share = IO_DARK_CONTAMINATED * IO_DARK_PROMETHEUS;
    bins.par_iter()
        .zip(lava.par_iter())
        .zip(laid.par_iter())
        .map(|((&(lo, hi), lava), laid)| {
            let n = coarse.cell_count();
            let (mut ring, mut fill, mut dark) =
                (vec![0.0f64; n], vec![0.0f64; n], vec![0.0f64; n]);
            let mut stamper = Stamper::new(&coarse);
            let years = (hi - lo).max(1e-9);
            for (s, &area) in lava.iter().enumerate() {
                if area <= 0.0 {
                    continue;
                }
                // Its SO2 a year over the bin that condenses in the canopy, in optical
                // depths over a square metre.
                let rate = area * FLOOD_SHEET_THICKNESS * IO_FRONT_SO2 * IO_RING_GAS_SHARE
                    / (years * IO_DEPOSIT_OPTICAL_MASS);
                let center = coarse.center(coarse.cell(s));
                let radius = radius_at(center);
                let per = rate / (tau * radius * radius * moment);
                stamper.gather(&coarse, body_radius, center, farthest * radius);
                for &(idx, r, _) in &stamper.footprint {
                    let x = r / radius;
                    if x >= farthest {
                        continue;
                    }
                    let (crest, within) = deposit_profile(x, IO_RING_OUTER_WIDTH);
                    let within = IO_RING_INNER_PROMETHEUS * within;
                    ring[idx] += per * crest;
                    fill[idx] += per * within;
                    dark[idx] += per * within * dark_share * (-(x / IO_DARK_REACH).powi(2)).exp();
                }
            }
            let fine = |v: Vec<f64>| {
                let mut f = Field::<f32>::new(coarse);
                for (o, x) in f.as_mut_slice().iter_mut().zip(v) {
                    *o = x as f32;
                }
                f.resample(grid).as_slice().to_vec()
            };
            // The rest of the gas settles across the region, its mass kept by the
            // smoothing.
            let settles = FLOOD_SHEET_THICKNESS * IO_FRONT_SO2 * (1.0 - IO_RING_GAS_SHARE)
                / (years * IO_DEPOSIT_OPTICAL_MASS);
            let gas = smoothed_over(laid, IO_PLUME_GAS_LENGTH, body_radius)
                .into_iter()
                .map(|v| (v as f64 * settles) as f32)
                .collect();
            FrontFall {
                start: -hi,
                end: -lo,
                ring: fine(ring),
                fill: fine(fill),
                dark: fine(dark),
                gas,
            }
        })
        .collect()
}

/// The far fallout of the plumes on a grid, in proportion (F89): about each volcanic
/// centre over IO_FAR_FALLOUT_REACH of its class's range, weighted by its plumes' output,
/// the share of the time it throws one times their mass.
fn far_fallout(grid: CubeSphere, body_radius: f64, gravity: f64, paterae: &[Patera]) -> Field<f32> {
    let g = gravity.max(0.01);
    let sources: Vec<([f64; 3], f64, f64)> = paterae
        .iter()
        .map(|p| {
            let (share, mass, speed) = if p.pele {
                (
                    IO_PLUME_FRACTION_PELE,
                    IO_PLUME_MASS_PELE,
                    IO_PLUME_SPEED_PELE,
                )
            } else {
                (
                    IO_PLUME_FRACTION_PROMETHEUS,
                    IO_PLUME_MASS_PROMETHEUS,
                    IO_PLUME_SPEED_PROMETHEUS,
                )
            };
            let reach = IO_FAR_FALLOUT_REACH * speed * speed / g;
            (p.center, IO_PATERA_ACTIVE_FRACTION * share * mass, reach)
        })
        .collect();
    Field::from_fn(grid, |_, d| {
        sources
            .iter()
            .map(|&(c, w, reach)| {
                let r = angle_between(c, d) * body_radius;
                w * (-r / reach).exp() / (reach * reach)
            })
            .sum::<f64>() as f32
    })
}

/// How the top layer of a cell ages, per year (F84): red sulfur's reversion to yellow,
/// the frost's sublimation in optical depths, and the plasma's alteration where the
/// atmosphere does not shield the ground.
#[derive(Clone, Copy, Debug)]
struct Kinetics {
    revert: f64,
    sublime: f64,
    alter: f64,
}

/// The fallout reaching a cell, optical depths a year: the sulfur and silicate that mix
/// into the top layer, with its red and dark shares as rates, and the frost.
#[derive(Clone, Copy, Debug, Default)]
struct Fall {
    solid: f64,
    red: f64,
    dark: f64,
    frost: f64,
}

impl Fall {
    fn add(&mut self, o: &Fall, share: f64) {
        self.solid += share * o.solid;
        self.red += share * o.red;
        self.dark += share * o.dark;
        self.frost += share * o.frost;
    }
}

/// A share after `dt` years of dx/dt = source - rate x, in closed form, exact as the
/// rate goes to nothing.
fn relax(x: f64, source: f64, rate: f64, dt: f64) -> f64 {
    let y = rate * dt;
    let grown = if y > 1e-8 {
        -(-y).exp_m1() / y
    } else {
        1.0 - 0.5 * y
    };
    x * (-y).exp() + source * dt * grown
}

/// The cover of sulfur and radiation products over a frost column beyond which the
/// frost no longer shows, optical depths.
const LAG_CAP: f64 = 30.0;

/// The slopes shorter than this, m, the frost reads as a spread of its temperature by
/// day, and so of its sublimation, over these factors of the mean rate (stated choices:
/// a few kelvin either way, the plains' slopes of a few degrees over their relief of a
/// few tens of kilometres), which blends a frost line across the width its ground gives
/// it, in steps fine enough not to show where a broad fall crosses them (F205); longer
/// slopes, a mountain's flanks, it reads one by one.
const FROST_SLOPE_LENGTH: f64 = 30e3;
const FROST_SPREAD: [f64; 16] = [
    0.333, 0.386, 0.447, 0.517, 0.599, 0.693, 0.803, 0.929, 1.076, 1.246, 1.442, 1.670, 1.933,
    2.238, 2.591, 3.000,
];

/// A field smoothed over a length, m: averaged onto coarser grids until their cells are
/// half the length, blurred there, and brought back grid by grid with a blur at each, so
/// that no coarse cell shows.
fn smoothed_over(field: &Field<f32>, length: f64, body_radius: f64) -> Vec<f32> {
    let mut grids = vec![field.grid()];
    let mut f = field.clone();
    while f.grid().delta() * body_radius < 0.5 * length {
        match crate::stages::relaxation::halve(&f) {
            Some(c) => {
                grids.push(c.grid());
                f = c;
            }
            None => break,
        }
    }
    let mut f = crate::stages::relaxation::blurred(&f, 2);
    for g in grids.iter().rev().skip(1) {
        f = crate::stages::relaxation::blurred(&f.resample(*g), 2);
    }
    f.as_slice().to_vec()
}

/// The top layer of a cell (F84): the share of the ground that is still bare lava; the
/// dark silicate, red and altered shares of the sulfur over it; and the frost column
/// over all, optical depths as it fell, at each of the cell's spread of sublimation
/// rates, under the cover of sulfur and radiation products laid on it since it fell;
/// and how much of the frost fell within about IO_FROST_COARSENING_YR, still fine, the
/// rest coarse (F203).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
struct Top {
    lava: f64,
    dark: f64,
    red: f64,
    altered: f64,
    frost: [f64; 16],
    lag: f64,
    fine: f64,
}

impl Top {
    /// The layer in balance with a steady fall, the frost column what `span` years of
    /// its excess over the sublimation laid.
    fn steady(fall: Fall, k: Kinetics, span: f64) -> Self {
        let s = fall.solid;
        let frost = FROST_SPREAD.map(|f| ((fall.frost - f * k.sublime) * span).max(0.0));
        Self {
            lava: 0.0,
            dark: if s > 0.0 { fall.dark / s } else { 0.0 },
            red: if s + k.revert > 0.0 {
                fall.red / (s + k.revert)
            } else {
                0.0
            },
            altered: if s + k.alter > 0.0 {
                k.alter / (s + k.alter)
            } else {
                0.0
            },
            frost,
            lag: if fall.frost > 0.0 {
                ((s + k.alter) / fall.frost).min(LAG_CAP)
            } else {
                LAG_CAP
            },
            fine: fall.frost * IO_FROST_COARSENING_YR,
        }
    }

    /// The layer after `dt` years of a steady fall, in closed form: each share relaxes
    /// towards its balance at the rate of its burial and its ageing together.
    fn age(&mut self, dt: f64, fall: Fall, k: Kinetics) {
        if dt <= 0.0 {
            return;
        }
        let s = fall.solid;
        // Rough lava hides only under several optical depths of fallout.
        self.lava *= (-s / IO_LAVA_COAT_DEPTHS * dt).exp();
        self.dark = relax(self.dark, fall.dark, s, dt);
        self.red = relax(self.red, fall.red, s + k.revert, dt);
        self.altered = relax(self.altered, k.alter, s + k.alter, dt);
        // The frost column grows by its fall and sublimes at its latitude's rate; the
        // sulfur falling with it and the radiation products cover it, and fresh frost
        // covers them in turn.
        for (column, spread) in self.frost.iter_mut().zip(FROST_SPREAD) {
            let grow = fall.frost - spread * k.sublime;
            if *column > 0.0 || grow > 0.0 {
                *column = (*column + grow * dt).max(0.0);
            }
        }
        self.lag = relax(self.lag, s + k.alter, fall.frost, dt).min(LAG_CAP);
        // Fresh frost falls fine and coarsens (F203).
        self.fine = relax(self.fine, fall.frost, 1.0 / IO_FROST_COARSENING_YR, dt);
    }

    /// Fresh lava laid over a share of the cell.
    fn lava_laid(&mut self, share: f64) {
        let keep = 1.0 - share;
        self.lava = self.lava * keep + share;
        self.dark *= keep;
        self.red *= keep;
        self.altered *= keep;
        self.frost = self.frost.map(|f| f * keep);
        self.lag *= keep;
        self.fine *= keep;
    }

    fn blend(&self, other: &Top, share: f64) -> Top {
        let mix = |a: f64, b: f64| a + (b - a) * share;
        Top {
            lava: mix(self.lava, other.lava),
            dark: mix(self.dark, other.dark),
            red: mix(self.red, other.red),
            altered: mix(self.altered, other.altered),
            frost: std::array::from_fn(|i| mix(self.frost[i], other.frost[i])),
            lag: mix(self.lag, other.lag),
            fine: mix(self.fine, other.fine),
        }
    }

    /// The share of the ground the frost hides, over the cell's spread of sublimation:
    /// its fine frost in full, its coarse as many times less as its grains are larger,
    /// the column subliming its coarse frost first (F203).
    fn visible_frost(&self) -> f64 {
        let columns: f64 = self
            .frost
            .iter()
            .map(|&f| {
                let fine = self.fine.min(f);
                1.0 - (-(fine + (f - fine) / IO_FROST_COARSE_GRAIN)).exp()
            })
            .sum();
        columns / self.frost.len() as f64 * (-self.lag).exp()
    }
}

/// A patera's plumes as one group for the search: its vent, the local frame there, the
/// farthest any of its deposits reaches, and its episodes.
struct Vent {
    center: [f64; 3],
    east: [f64; 3],
    north: [f64; 3],
    cos_cap: f64,
    episodes: std::ops::Range<usize>,
}

/// Latitude and longitude bins of the vents' search.
const SEARCH_BINS: (usize, usize) = (90, 180);

fn search_bin(d: [f64; 3]) -> usize {
    let lat = d[2].clamp(-1.0, 1.0).asin();
    let lon = d[1].atan2(d[0]);
    let row = (((lat + std::f64::consts::FRAC_PI_2) / std::f64::consts::PI * SEARCH_BINS.0 as f64)
        as usize)
        .min(SEARCH_BINS.0 - 1);
    let col = (((lon + std::f64::consts::PI) / std::f64::consts::TAU * SEARCH_BINS.1 as f64)
        as usize)
        .min(SEARCH_BINS.1 - 1);
    row * SEARCH_BINS.1 + col
}

/// Io's surface from what its plumes have left on it (F84). Every cell's top layer
/// starts at its balance with the history's mean fall, or as lava when a flow or a
/// floor was laid before the history began, and runs through the eruptions that reached
/// it in their order: each fall buries what lay on top by its optical depth, the red
/// sulfur reverts, the frost sublimes at the rate the vapour pressure of its own
/// temperature gives (the frost's albedo, the latitude and the cell's slope towards the
/// pole), and the plasma alters what lies on top where the column of the sublimation
/// atmosphere, which the same vapour pressure sets, is too thin to stop it. A flow or a
/// floor laid during the history starts black over its coverage at its time. A flow
/// spreading over the plains drives off the sulfur dioxide in the ground it covers,
/// which condenses again round it as a halo of frost (Williams et al. 2011: white
/// diffuse deposits, often round lava flow margins, from SO2 driven off by warm flows and
/// condensed again; Prometheus's plume rises where its flows advance over SO2-rich
/// plains, Kieffer et al. 2000): each flow's release spread over the halo's reach,
/// landing round the young lava rather than on it. The small eruptions fall with the
/// plumes, and all of it under the plumes' far fallout about the volcanic centres,
/// `background` at each cell (F89).
#[allow(clippy::too_many_arguments)]
pub fn deposit_history(
    layers: &mut LayerSet,
    body_radius: f64,
    gravity: f64,
    insolation: f64,
    episodes: &[PlumeEpisode],
    small: &[PlumeEpisode],
    background: &[f32],
    slope: Option<&Field<f32>>,
) {
    let grid = layers.grid;
    let window = history_years();
    let all: Vec<&PlumeEpisode> = episodes.iter().chain(small.iter()).collect();
    let plumes: Vec<Plume> = all
        .iter()
        .map(|e| Plume::new(e, body_radius, window))
        .collect();
    let smalls: Vec<bool> = all.iter().map(|e| e.small).collect();
    // The vents, each patera's episodes being consecutive.
    let mut vents: Vec<Vent> = Vec::new();
    let mut i = 0;
    while i < all.len() {
        let center = all[i].center;
        let mut j = i;
        let mut cos_cap: f64 = 1.0;
        while j < all.len() && all[j].center == center {
            cos_cap = cos_cap.min(plumes[j].cos_cap);
            j += 1;
        }
        let (east, north) = local_basis(center);
        vents.push(Vent {
            center,
            east,
            north,
            cos_cap,
            episodes: i..j,
        });
        i = j;
    }
    // Each bin lists the vents whose deposits may reach into it.
    let mut bins: Vec<Vec<u32>> = vec![Vec::new(); SEARCH_BINS.0 * SEARCH_BINS.1];
    let row_height = std::f64::consts::PI / SEARCH_BINS.0 as f64;
    let col_width = std::f64::consts::TAU / SEARCH_BINS.1 as f64;
    for (v, vent) in vents.iter().enumerate() {
        let cap = vent.cos_cap.clamp(-1.0, 1.0).acos();
        let lat = vent.center[2].clamp(-1.0, 1.0).asin();
        let lon = vent.center[1].atan2(vent.center[0]);
        let lo = ((lat - cap + std::f64::consts::FRAC_PI_2) / row_height).floor() - 1.0;
        let hi = ((lat + cap + std::f64::consts::FRAC_PI_2) / row_height).floor() + 1.0;
        let all = lat.abs() + cap >= std::f64::consts::FRAC_PI_2 - row_height;
        let half = if all {
            std::f64::consts::PI
        } else {
            (cap.sin() / lat.cos()).min(1.0).asin() + col_width
        };
        let c0 = ((lon - half + std::f64::consts::PI) / col_width).floor() as i64;
        let c1 = ((lon + half + std::f64::consts::PI) / col_width).floor() as i64;
        for row in (lo.max(0.0) as usize)..=(hi.min(SEARCH_BINS.0 as f64 - 1.0) as usize) {
            if all || c1 - c0 + 1 >= SEARCH_BINS.1 as i64 {
                for col in 0..SEARCH_BINS.1 {
                    bins[row * SEARCH_BINS.1 + col].push(v as u32);
                }
            } else {
                for c in c0..=c1 {
                    let col = c.rem_euclid(SEARCH_BINS.1 as i64) as usize;
                    bins[row * SEARCH_BINS.1 + col].push(v as u32);
                }
            }
        }
    }
    let noise = GradientNoise::new(0x10_DE90_5175);
    let frost_k = (insolation * (1.0 - IO_FROST_ALBEDO) / SIGMA_SB)
        .max(1.0)
        .powf(0.25);
    let sublime_equator = IO_FROST_SUBLIMATION / IO_DEPOSIT_OPTICAL_MASS;
    let units = layers.unit_id.as_slice();
    let tilts = slope.map(|s| smoothed_over(s, FROST_SLOPE_LENGTH, body_radius));
    let coverage = layers.cover.as_slice();
    let laid_yr = layers.melt.as_slice();
    // The plumes at the young flows' fronts (F202).
    let front_clock = std::time::Instant::now();
    let fronts = front_falls(grid, body_radius, gravity, units, coverage, laid_yr, window);
    let front_frost = 0.5 * (IO_PROMETHEUS_FROST.0 + IO_PROMETHEUS_FROST.1);
    let front_mottle = body_radius
        / (0.15 * IO_PLUME_SPEED_PROMETHEUS * IO_PLUME_SPEED_PROMETHEUS / gravity.max(0.01));
    let debug = std::env::var_os("EXOGEN_DEBUG").is_some();
    if debug {
        eprintln!(
            "io deposit history: the fronts' plumes and the slopes {:.1} s",
            front_clock.elapsed().as_secs_f64()
        );
    }
    let clock = std::time::Instant::now();
    // Probes of cells' histories (EXOGEN_IO_TRACE=lon,lat;lon,lat in degrees).
    let traced: Vec<usize> = std::env::var("EXOGEN_IO_TRACE")
        .map(|v| {
            v.split(';')
                .filter_map(|pair| {
                    let mut it = pair.split(',').map(|x| x.trim().parse::<f64>());
                    match (it.next(), it.next()) {
                        (Some(Ok(lon)), Some(Ok(lat))) => Some(
                            grid.index(
                                grid.locate(crate::grid::direction_from_lon_lat(
                                    lon.to_radians(),
                                    lat.to_radians(),
                                ))
                                .0,
                            ),
                        ),
                        _ => None,
                    }
                })
                .collect()
        })
        .unwrap_or_default();
    let tops: Vec<Top> = (0..grid.cell_count())
        .into_par_iter()
        .map_init(
            || {
                (
                    Vec::<(f64, f64, Fall, bool)>::new(),
                    Vec::<(f64, u8, usize)>::new(),
                    Vec::<usize>::new(),
                )
            },
            |(falls, marks, active), k| {
                let d = grid.center(grid.cell(k));
                // The far fallout of the plumes about the volcanic centres (F89).
                let background = Fall {
                    solid: background[k] as f64,
                    ..Fall::default()
                };
                let expected = background;
                // The frost's temperature by day against the equator's sets its vapour
                // pressure: the frost sublimes at its own slope's, and the atmosphere the
                // sublimation supports shields the ground by the latitude's, being
                // spread over far more than a slope.
                let pressure = |tilt: f64| {
                    let cos_lat = crate::stages::ice::effective_cos_lat(d, tilt).max(0.02);
                    (SO2_SUBLIMATION_K * (1.0 / frost_k - 1.0 / (frost_k * cos_lat.powf(0.25))))
                        .exp()
                };
                let kinetics = Kinetics {
                    revert: 1.0 / IO_RED_FADE_YR,
                    sublime: sublime_equator
                        * pressure(tilts.as_ref().map_or(0.0, |t| t[k] as f64)),
                    alter: (-IO_SHIELD_COLUMNS * pressure(0.0)).exp() / IO_ALTERATION_YR,
                };
                // The eruptions whose fallout reached the cell, and when.
                falls.clear();
                let mut cos_k = [0.0f64; RANGE_ORDERS.1 + 1];
                let mut sin_k = [0.0f64; RANGE_ORDERS.1 + 1];
                for &v in &bins[search_bin(d)] {
                    let vent = &vents[v as usize];
                    let c = dot(vent.center, d);
                    if c < vent.cos_cap {
                        continue;
                    }
                    let r = angle_between(vent.center, d) * body_radius;
                    let rel = [
                        d[0] - vent.center[0],
                        d[1] - vent.center[1],
                        d[2] - vent.center[2],
                    ];
                    let (x, y) = (dot(rel, vent.east), dot(rel, vent.north));
                    let len = x.hypot(y);
                    let (c1, s1) = if len > 1e-15 {
                        (x / len, y / len)
                    } else {
                        (1.0, 0.0)
                    };
                    cos_k[0] = 1.0;
                    sin_k[0] = 0.0;
                    for n in 1..=RANGE_ORDERS.1 {
                        cos_k[n] = cos_k[n - 1] * c1 - sin_k[n - 1] * s1;
                        sin_k[n] = sin_k[n - 1] * c1 + cos_k[n - 1] * s1;
                    }
                    for (e, &is_small) in plumes[vent.episodes.clone()]
                        .iter()
                        .zip(&smalls[vent.episodes.clone()])
                    {
                        if c < e.cos_cap || e.end <= e.start {
                            continue;
                        }
                        let x = r / (e.radius * e.range_at(&cos_k, &sin_k));
                        if x >= e.reach {
                            continue;
                        }
                        let (ring, fill) = deposit_profile(x, e.outer);
                        // A small eruption rises too low for a canopy: its fall is a smooth
                        // bump about the vent, frost and sulfur together (F89).
                        let (ring, fill) = if is_small {
                            (0.0, small_profile(x))
                        } else {
                            (ring, e.inner * fill)
                        };
                        if e.crest * (ring + fill) * (e.end - e.start) < NEGLIGIBLE_FALL {
                            continue;
                        }
                        let (f, m) = (e.mottle_frequency, &e.mottle);
                        let mottle =
                            noise.sample([d[0] * f + m[0], d[1] * f + m[1], d[2] * f + m[2]]);
                        let mass = e.crest * (1.0 - IO_DEPOSIT_MOTTLING * (0.5 + 0.5 * mottle));
                        // The ring's gas condensate, frost and the plume's red or yellow
                        // sulfur; the pyroclasts within, sulfur and silicate.
                        let ring = mass * ring;
                        let fill = mass * fill;
                        let dark =
                            e.dark * (-(x / IO_DARK_REACH).powi(2)).exp() * e.jet_at(s1.atan2(c1));
                        let fall = if is_small {
                            let solid = fill * (1.0 - e.frost);
                            Fall {
                                solid,
                                red: solid * e.red,
                                dark: solid * dark,
                                frost: fill * e.frost,
                            }
                        } else {
                            let condensed = ring * (1.0 - e.frost);
                            Fall {
                                solid: condensed + fill,
                                red: condensed * e.red,
                                dark: fill * dark,
                                frost: ring * e.frost,
                            }
                        };
                        falls.push((e.start, e.end, fall, is_small));
                    }
                }
                // The plumes at the young flows' fronts about the cell (F202), their mass
                // mottled as an eruption's is.
                let mottle = 1.0
                    - IO_DEPOSIT_MOTTLING
                        * (0.5
                            + 0.5
                                * noise.sample([
                                    d[0] * front_mottle + 37.0,
                                    d[1] * front_mottle + 11.0,
                                    d[2] * front_mottle + 5.0,
                                ]));
                for front in &fronts {
                    let ring = front.ring[k] as f64 * mottle;
                    let fill = front.fill[k] as f64 * mottle;
                    let gas = front.gas[k] as f64;
                    if (ring + fill + gas) * (front.end - front.start) < NEGLIGIBLE_FALL {
                        continue;
                    }
                    let condensed = ring * (1.0 - front_frost);
                    falls.push((
                        front.start,
                        front.end,
                        Fall {
                            solid: condensed + fill,
                            red: condensed * IO_PROMETHEUS_RED,
                            dark: front.dark[k] as f64 * mottle,
                            frost: ring * front_frost + gas,
                        },
                        false,
                    ));
                }
                // Where the cell stood when the history began: what fell before has aged,
                // its red reverted and its frost gone or covered, and it was buried at the
                // history's mean rate, the small eruptions' left to the background.
                let mut mean = expected;
                for (start, end, fall, _) in falls.iter().filter(|f| !f.3) {
                    let w = (end - start) / window;
                    mean.add(
                        &Fall {
                            solid: fall.solid * w,
                            dark: fall.dark * w,
                            ..Fall::default()
                        },
                        1.0,
                    );
                }
                let mut top = Top::steady(mean, kinetics, window);
                // A flow or a floor: laid before the history began, it has been coated
                // since by the mean fall; laid during it, it is laid in its turn.
                let lava = matches!(units[k], UNIT_IO_FLOW | UNIT_PATERA_FLOOR);
                let share = if lava {
                    (coverage[k] as f64).clamp(0.0, 1.0)
                } else {
                    0.0
                };
                let laid = -(laid_yr[k] as f64).max(0.0);
                if share > 0.0 && laid <= -window {
                    let mut fresh = Top {
                        lava: 1.0,
                        ..Top::default()
                    };
                    fresh.age(-window - laid, mean, kinetics);
                    top = top.blend(&fresh, share);
                }
                // The eruptions' starts and ends, and the lava, in order.
                marks.clear();
                for (n, (start, end, _, _)) in falls.iter().enumerate() {
                    marks.push((*start, 0, n));
                    marks.push((*end, 1, n));
                }
                if share > 0.0 && laid > -window {
                    marks.push((laid, 2, 0));
                }
                marks.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
                let mut now = -window;
                // The fall at each moment is the sum of the eruptions then reaching the
                // cell, summed afresh at each start and end.
                let mut fall = background;
                active.clear();
                for &(t, kind, n) in marks.iter() {
                    top.age(t - now, fall, kinetics);
                    now = now.max(t);
                    match kind {
                        0 => active.push(n),
                        1 => active.retain(|&a| a != n),
                        _ => top.lava_laid(share),
                    }
                    if kind < 2 {
                        fall = background;
                        for &a in active.iter() {
                            fall.add(&falls[a].2, 1.0);
                        }
                    }
                }
                top.age(-now, fall, kinetics);
                if traced.contains(&k) {
                    eprintln!(
                        "io trace: cell {k} unit {} coverage {:.3} laid {:.1} yr, {kinetics:?}",
                        units[k], coverage[k], laid_yr[k]
                    );
                    for (start, end, fall, small) in falls.iter() {
                        eprintln!(
                            "io trace:   fall {start:.1} to {end:.1} yr (small {small}): {fall:?}"
                        );
                    }
                    for (b, front) in fronts.iter().enumerate() {
                        eprintln!(
                            "io trace:   front bin {b}: ring {:.3e} fill {:.3e} a year",
                            front.ring[k], front.fill[k]
                        );
                    }
                    eprintln!(
                        "io trace:   top {top:?}, frost shows {:.3}",
                        top.visible_frost()
                    );
                }
                top
            },
        )
        .collect();
    if debug {
        eprintln!(
            "io deposit history: the cells {:.1} s",
            clock.elapsed().as_secs_f64()
        );
        // The lava still showing, by area: what Williams et al. (2011) would map as flows.
        let mut shares = [0.0f64; 4];
        let mut total = 0.0;
        for (k, top) in tops.iter().enumerate() {
            let a = grid.cell_area(grid.cell(k));
            total += a;
            for (s, lim) in shares.iter_mut().zip([0.05, 0.2, 0.5, 0.8]) {
                if top.lava > lim {
                    *s += a;
                }
            }
        }
        eprintln!(
            "io deposit history: lava showing over 0.05/0.2/0.5/0.8 of the ground: {:.3} {:.3} {:.3} {:.3}",
            shares[0] / total,
            shares[1] / total,
            shares[2] / total,
            shares[3] / total
        );
    }
    let deep = layers.deep.as_mut_slice();
    let melt = layers.melt.as_mut_slice();
    let cover = layers.cover.as_mut_slice();
    let stain = layers.stain.as_mut_slice();
    for (k, top) in tops.iter().enumerate() {
        // Frost does not hold on the lava still showing (F206): dark ground takes in some
        // three times the sunlight the frost and the plains do and stands tens of kelvin
        // warmer by day, where the frost's vapour pressure is thousands of times its own,
        // so the frost that falls on it by night is gone by noon; only where the fallout
        // has mantled the lava bright does the frost lie.
        let bare = top.lava;
        deep[k] = (top.lava + (1.0 - top.lava) * top.dark) as f32;
        melt[k] = top.red as f32;
        stain[k] = top.altered as f32;
        cover[k] = (top.visible_frost() * (1.0 - bare)) as f32;
    }
    if std::env::var_os("EXOGEN_DEBUG").is_some() {
        // The frost showing by band of latitude, and over the whole (Geissler et al.
        // 1999: white and grey SO2 over 27 percent of Io, mostly at low latitudes).
        let edges = [0.0, 10.0, 20.0, 30.0, 45.0, 60.0, 90.0];
        let mut sums = [[0.0f64; 2]; 6];
        let mut all = [0.0f64; 2];
        for (k, &c) in cover.iter().enumerate() {
            let cell = grid.cell(k);
            let a = grid.cell_area(cell);
            let lat = grid.center(cell)[2]
                .clamp(-1.0, 1.0)
                .asin()
                .to_degrees()
                .abs();
            let b = edges
                .windows(2)
                .position(|e| lat >= e[0] && lat < e[1])
                .unwrap_or(5);
            sums[b][0] += c as f64 * a;
            sums[b][1] += a;
            all[0] += c as f64 * a;
            all[1] += a;
        }
        let bands: Vec<String> = sums
            .iter()
            .zip(edges.windows(2))
            .map(|(s, e)| format!("{}-{}: {:.3}", e[0], e[1], s[0] / s[1]))
            .collect();
        eprintln!(
            "io frost showing {:.3}; by latitude {}",
            all[0] / all[1],
            bands.join(", ")
        );
    }
}

/// The Io regime's features through the loop: a steady state, laid down afresh at every
/// step on the base terrain, with the full timeline kept so any epoch can be finished.
pub struct IoProcess {
    grid: CubeSphere,
    body_radius: f64,
    age: f64,
    dt: f64,
    rng: Rng,
    n_paterae: usize,
    n_mountains: usize,
    /// How deep the plains are buried since their surface formed, m (F86).
    burial_depth: f64,
    /// The tide that orients the crust's faults (F87).
    tide: Option<crate::stages::tidal::TidalShell>,
    /// The terrain before the regime's features, which every step starts from, and the
    /// finish too.
    pub base: Option<LayerSet>,
    pub paterae: Vec<Patera>,
    pub mountains: Vec<Mountain>,
    stamper: Stamper,
}

impl IoProcess {
    /// None unless the body is a rocky world heated past the regime's threshold.
    pub fn new(graph: &ParamGraph, grid: CubeSphere) -> Option<Self> {
        if !is_io_regime(graph) {
            return None;
        }
        Some(Self {
            grid,
            body_radius: graph.f64("body.radius"),
            age: graph.f64("body.age"),
            dt: 0.0,
            rng: Rng::new(graph.f64("seed.volcanism") as u64 ^ 0x10_5EED_0001),
            n_paterae: graph.f64("volcanism.paterae_count").round().max(0.0) as usize,
            n_mountains: graph.f64("tectonics.mountain_count").round().max(0.0) as usize,
            burial_depth: graph.f64("volcanism.resurfacing_rate") * IO_SURFACE_AGE_GYR * 1e3,
            tide: crate::stages::tidal::TidalShell::new(graph),
            base: None,
            paterae: Vec::new(),
            mountains: Vec::new(),
            stamper: Stamper::new(&grid),
        })
    }

    fn log_uniform(&self, index: u64, draw: u64, range: (f64, f64)) -> f64 {
        range.0 * (range.1 / range.0).powf(self.rng.uniform(index, draw))
    }

    /// The regime's surface afresh, as every step lays it: the plains' age and unit,
    /// nothing on them but the lava of its floods.
    pub fn clear_surface(layers: &mut LayerSet, t: f64) {
        let now = (t - IO_SURFACE_AGE_GYR) as f32;
        let ages = layers.surface_age.as_mut_slice();
        let units = layers.unit_id.as_mut_slice();
        let deep = layers.deep.as_mut_slice();
        let melt = layers.melt.as_mut_slice();
        let cover = layers.cover.as_mut_slice();
        let stain = layers.stain.as_mut_slice();
        for k in 0..ages.len() {
            // The floods' lava the step laid stays (F117): the finish's flood passes take
            // it up as they do on any body.
            if units[k] == UNIT_LAVA {
                continue;
            }
            ages[k] = now;
            units[k] = UNIT_PRIMORDIAL;
            deep[k] = 0.0;
            melt[k] = 0.0;
            cover[k] = 0.0;
            stain[k] = 0.0;
        }
    }

    /// The mountains through the body's life (F87): each of the standing count is a chain
    /// of blocks, one raised as the last is buried, drawn from the seed alone; each a
    /// block thrust up along a fault striking with the crust's fabric, its length, width
    /// and height log-normal (Schenk et al. 2001: Io's mountains 157 km long on average,
    /// the longest 570 km, 6.3 km high on average, the highest 17.5 km), its form a
    /// lineated block's, a layered plateau's or mesa's, or a massif's.
    fn raise_mountains(&mut self) {
        let log_normal = |median: f64, sigma: f64, k: u64, draw: u64, max: f64| {
            (median * (sigma * self.rng.normal(k, draw)).exp()).min(max)
        };
        for slot in 0..self.n_mountains as u64 {
            let mut generation = 0u64;
            let first = 5_000 + slot * 1_000;
            let mut birth =
                -self.log_uniform(first, 7, IO_MOUNTAIN_LIFETIME_GYR) * self.rng.uniform(first, 8);
            while birth < self.age && generation < 900 {
                let k = first + generation;
                let center = self.rng.unit_direction(k, 9);
                let (l_med, l_sigma) = IO_MOUNTAIN_LENGTH;
                let length = log_normal(l_med, l_sigma, k, 20, IO_MOUNTAIN_LENGTH_MAX);
                let (a_med, a_sigma) = IO_MOUNTAIN_ASPECT;
                let width = length / log_normal(a_med, a_sigma, k, 22, 8.0).max(1.0);
                let (h_med, h_sigma) = IO_MOUNTAIN_HEIGHT;
                let height = log_normal(h_med, h_sigma, k, 24, IO_MOUNTAIN_HEIGHT_MAX);
                let life = self.log_uniform(k, 7, IO_MOUNTAIN_LIFETIME_GYR);
                // Its form, by the shares of Io's mountain units.
                let form = self.rng.uniform(k, 26);
                let (level, ridging, tilt) = if form < IO_MOUNTAIN_FORMS.0 {
                    (0.3, 1.0, IO_MOUNTAIN_TILT)
                } else if form < IO_MOUNTAIN_FORMS.0 + IO_MOUNTAIN_FORMS.1 {
                    (0.9, 0.0, 0.3 * IO_MOUNTAIN_TILT)
                } else {
                    (0.0, 0.0, 0.5 * IO_MOUNTAIN_TILT)
                };
                self.mountains.push(Mountain {
                    center,
                    strike: fabric_strike(
                        self.tide.as_ref(),
                        center,
                        self.rng.uniform(k, 27),
                        self.rng.uniform(k, 28),
                    ),
                    length,
                    width,
                    height,
                    vergence: if self.rng.uniform(k, 29) < 0.5 {
                        -1.0
                    } else {
                        1.0
                    },
                    tilt: tilt * self.rng.uniform(k, 5),
                    level,
                    ridging,
                    birth_gyr: birth,
                    death_gyr: birth + life,
                    id: k,
                });
                birth += life;
                generation += 1;
            }
        }
    }

    /// The features alive through a step: paterae and flows are renewed at every kept
    /// step, since the surface is younger than the step.
    fn populate(&mut self, t1: f64, lay: bool) {
        // The time bin the step ends in is the random index, so the same surface comes
        // out at every quality.
        let step = (t1 / self.age * crate::stages::LOOP_TIME_BINS as f64).round() as u64;
        let death = t1 + self.dt.max(1e-6);
        // Paterae for the steps whose surface is kept; their flows are the finish's
        // (F85).
        if !lay {
            return;
        }
        for k in 0..self.n_paterae as u64 {
            let idx = 1_000_000 + step * 100_003 + k;
            // Low latitudes favoured: draw until the latitude weight accepts.
            let mut center = self.rng.unit_direction(idx, 1);
            for attempt in 0..12u64 {
                let lat = center[2].asin();
                let w = IO_PATERA_POLAR_WEIGHT + (1.0 - IO_PATERA_POLAR_WEIGHT) * lat.cos().powi(2);
                if self.rng.uniform(idx ^ crate::rng::ACCEPT_STREAM, attempt) < w {
                    break;
                }
                center = self
                    .rng
                    .unit_direction(idx ^ crate::rng::ACCEPT_STREAM, 100 + 2 * attempt);
            }
            // Its size, log-normal (Radebaugh et al. 2001: 417 paterae, mean 41.0 km
            // across, median 34.5, from 2.5 to 203 km).
            let (d_med, d_sigma) = IO_PATERA_DIAMETER;
            let (d_lo, d_hi) = IO_PATERA_DIAMETER_RANGE;
            let diameter = (d_med * (d_sigma * self.rng.normal(idx, 60)).exp()).clamp(d_lo, d_hi);
            let radius = 0.5 * diameter;
            let pele = self.rng.uniform(idx, 43) < IO_PELE_CLASS_FRACTION;
            let active = self.rng.uniform(idx, 42) < IO_PATERA_ACTIVE_FRACTION;
            let floor_age_yr = if active {
                0.0
            } else {
                self.log_uniform(idx, 44, IO_FLOOR_AGE_YR)
            };
            // A share sink at the foot of a mountain's scarp, along its fault, where magma
            // rises along the fault (Radebaugh et al. 2001: 13 percent beside mountains);
            // the rest lie along the crust's fabric where they are.
            let near = if self.rng.uniform(idx, 62) < IO_PATERA_MOUNTAIN_SHARE {
                let standing: Vec<&Mountain> =
                    self.mountains.iter().filter(|m| m.alive(t1)).collect();
                let total: f64 = standing.iter().map(|m| m.length).sum();
                let mut pick = self.rng.uniform(idx, 63) * total;
                let mut chosen = None;
                for m in &standing {
                    if pick < m.length {
                        chosen = Some(*m);
                        break;
                    }
                    pick -= m.length;
                }
                chosen
            } else {
                None
            };
            // Its elongation along its axis (Radebaugh et al. 2001: ellipticities from 1
            // to 5.73, 70 percent above 1.25), half its collapses' own and half their
            // spread along it.
            let (e_spread, e_max) = IO_PATERA_ELLIPTICITY;
            let elongation = (e_spread * self.rng.normal(idx, 68)).abs().exp().min(e_max);
            let ellipticity = elongation.sqrt();
            let (center, axis, abuts) = match near {
                Some(m) => {
                    let u = m.length * 0.3 * (2.0 * self.rng.uniform(idx, 64) - 1.0);
                    let scarp = mountain_scarp(m.height, 0.0);
                    let throw = (0.5 * std::f64::consts::PI * u
                        / (MOUNTAIN_FAULT_SHARE * m.length))
                        .cos()
                        .powi(2);
                    // Its centre half its width out from the scarp's foot, so that it
                    // sinks against the block, whose foot is its margin there.
                    let v = -(scarp
                        + mountain_heave(m.height) * throw
                        + 0.5 * radius / elongation.sqrt());
                    let (c, s) = (m.strike.cos(), m.strike.sin());
                    let (x, y) = (u * c - m.vergence * v * s, u * s + m.vergence * v * c);
                    let (e, n) = local_basis(m.center);
                    let azimuth = y.atan2(x);
                    let at =
                        crate::geom::offset(m.center, e, n, x.hypot(y) / self.body_radius, azimuth);
                    // Its long axis along the block's foot.
                    let normal_world = [
                        -m.vergence * s * e[0] + m.vergence * c * n[0],
                        -m.vergence * s * e[1] + m.vergence * c * n[1],
                        -m.vergence * s * e[2] + m.vergence * c * n[2],
                    ];
                    let (pe, pn) = local_basis(at);
                    let normal = dot(normal_world, pn).atan2(dot(normal_world, pe));
                    (at, normal + std::f64::consts::FRAC_PI_2, Some(m.clone()))
                }
                None => (
                    center,
                    fabric_strike(
                        self.tide.as_ref(),
                        center,
                        self.rng.uniform(idx, 65),
                        self.rng.uniform(idx, 66),
                    ),
                    None,
                ),
            };
            // Its collapses: a chain along its axis, as calderas coalesce along a fault
            // (Radebaugh et al. 2001: Monan's worm shape), each an ellipse of the square
            // root of its elongation as long as its radius, their centres spread as far
            // as the elongation has it reach; and up to two more nested about its centre.
            let (ca, sa) = (axis.cos(), axis.sin());
            let at = |along: f64, aside: f64, r: f64| {
                [along * ca - aside * sa, along * sa + aside * ca, r]
            };
            let link = radius / ellipticity.sqrt();
            let reach = radius * (ellipticity - 1.0);
            let links = 1 + (2.0 * reach / (1.2 * radius)).ceil() as u64;
            let mut collapses = Vec::new();
            for c in 0..links {
                let along = if links > 1 {
                    -reach + 2.0 * reach * c as f64 / (links - 1) as f64
                } else {
                    0.0
                };
                // Each link of its own size, wandering across the axis by a third of
                // its width.
                let jitter = 0.15 * radius * (2.0 * self.rng.uniform(idx, 104 + 2 * c) - 1.0);
                let size = link * (0.85 + 0.3 * self.rng.uniform(idx, 112 + c));
                let aside = 0.35 * size / ellipticity.sqrt()
                    * (2.0 * self.rng.uniform(idx, 105 + 2 * c) - 1.0);
                collapses.push(at(along + jitter, aside, size));
            }
            let count = 1 + (3.0 * self.rng.uniform(idx, 67).powi(2)) as u64;
            for c in 1..count {
                let r = link * (0.4 + 0.5 * self.rng.uniform(idx, 70 + 4 * c));
                let along = (if self.rng.uniform(idx, 71 + 4 * c) < 0.5 {
                    -1.0
                } else {
                    1.0
                }) * r
                    * (0.2 + 0.6 * self.rng.uniform(idx, 72 + 4 * c));
                let aside =
                    0.25 * r / ellipticity * (2.0 * self.rng.uniform(idx, 73 + 4 * c) - 1.0);
                collapses.push(at(along, aside, r));
            }
            // It collapses piecemeal (F117): more collapses of comparable size, set off
            // from its centre in any direction, break out of the first, so that its margin
            // is a run of arcs of different curvature meeting at corners.
            let (n_lo, n_hi) = IO_PATERA_PIECES;
            let pieces = n_lo + self.rng.bits(idx, 199) % (n_hi - n_lo + 1);
            let turn = std::f64::consts::TAU * self.rng.uniform(idx, 198);
            for m in 0..pieces {
                // Spread round the patera, each in its own share of the turn.
                let phi = turn
                    + std::f64::consts::TAU
                        * (m as f64 + 0.5 + 0.8 * (self.rng.uniform(idx, 200 + 3 * m) - 0.5))
                        / pieces as f64;
                let (s_lo, s_hi) = IO_PATERA_PIECE_SIZE;
                let size = link * (s_lo + (s_hi - s_lo) * self.rng.uniform(idx, 201 + 3 * m));
                let (o_lo, o_hi) = IO_PATERA_PIECE_OFFSET;
                let d = link * (o_lo + (o_hi - o_lo) * self.rng.uniform(idx, 202 + 3 * m));
                collapses.push([d * phi.cos(), d * phi.sin(), size]);
            }
            // Straight sides where faults of the fabric cut it (Radebaugh et al. 2001: 42
            // percent with straight or irregular margins ending at sharp corners).
            let mut faults = Vec::new();
            if abuts.is_none() && self.rng.uniform(idx, 90) < IO_PATERA_FAULTED_SHARE {
                let sides = 1 + self.rng.bits(idx, 91) % 3;
                for f in 0..sides {
                    // Each a fault along one of the fabric's strikes, on either side.
                    let strike = fabric_strike(
                        self.tide.as_ref(),
                        center,
                        self.rng.uniform(idx, 92 + 4 * f),
                        self.rng.uniform(idx, 93 + 4 * f),
                    );
                    let side = std::f64::consts::PI * (self.rng.bits(idx, 94 + 4 * f) % 2) as f64;
                    let normal = strike + std::f64::consts::FRAC_PI_2 + side;
                    // Cutting across the outline as far out as it reaches that way.
                    let out = outline_radius(normal, axis, ellipticity, &collapses);
                    faults.push([
                        normal,
                        out * (0.55 + 0.35 * self.rng.uniform(idx, 95 + 4 * f)),
                    ]);
                }
            }
            self.paterae.push(Patera {
                center,
                radius,
                depth: IO_PATERA_DEPTH * (0.5 + self.rng.uniform(idx, 41)),
                active,
                floor_age_yr,
                pele,
                axis,
                ellipticity,
                collapses,
                faults,
                abuts,
                birth_gyr: t1,
                death_gyr: death,
                id: idx,
            });
        }
    }

    /// The surface at the end of a step: the base terrain and the mountains standing, on
    /// the coarse grid; the paterae are `cut_paterae`'s, after the step's floods (F117). Only the steps whose surface the
    /// loop keeps as a keyframe are laid; the rest keep the timeline going. The flows
    /// and the plumes' fallout are the finish's (F84, F85): they bear no relief the loop
    /// needs, and the finish lays them at its own resolution.
    pub fn step(&mut self, layers: &mut LayerSet, t0: f64, t1: f64, lay: bool) {
        match &self.base {
            None => {
                self.dt = t1 - t0;
                // The regime buries its plains as fast as it renews them (F86).
                let cell = self.grid.delta() * self.body_radius;
                let buried = bury(&layers.bedrock, cell, self.body_radius, self.burial_depth);
                layers.bedrock.as_mut_slice().copy_from_slice(&buried);
                self.base = Some(layers.clone());
            }
            Some(base) => {
                if lay {
                    *layers = base.clone();
                }
            }
        }
        if self.mountains.is_empty() {
            self.raise_mountains();
        }
        self.populate(t1, lay);
        if !lay {
            return;
        }
        let grid = self.grid;
        let body_radius = self.body_radius;
        let cell = grid.delta() * body_radius;
        Self::clear_surface(layers, t1);
        for m in self.mountains.iter().filter(|m| m.alive(t1)) {
            self.stamper
                .gather(&grid, body_radius, m.center, m.reach(2.0 * cell));
            let bedrock = layers.bedrock.as_mut_slice();
            let units = layers.unit_id.as_mut_slice();
            let cover = layers.cover.as_mut_slice();
            for &(idx, r, phi) in &self.stamper.footprint {
                // The loop's block is smooth, its scarp at least four of its cells wide;
                // the finish draws the block afresh at its own resolution.
                let (u, v) = m.frame(r, phi);
                let h = mountain_height(m, u, v, 2.0 * cell);
                if h > 0.0 {
                    bedrock[idx] += h as f32;
                    if h > MOUNTAIN_UNIT_SHARE * m.height {
                        units[idx] = UNIT_MOUNTAIN;
                        cover[idx] = 0.0;
                    }
                }
            }
        }
    }

    /// The paterae alive at the end of a kept step, cut into its surface after the floods
    /// of the step (F117): a patera is born at the step's end, after the lava of the
    /// centuries before it, and its collapse goes through that lava, its floor lava of
    /// its own.
    pub fn cut_paterae(&mut self, layers: &mut LayerSet, t1: f64) {
        let grid = self.grid;
        let body_radius = self.body_radius;
        let cell = grid.delta() * body_radius;
        for p in self.paterae.iter().filter(|p| p.alive(t1)) {
            // The loop's mountains are smooth (F87).
            let drawing = PateraDrawing::new(p, body_radius, cell, 2.0 * cell);
            self.stamper.gather(
                &grid,
                body_radius,
                p.center,
                p.extent() + 3.0 * drawing.wall,
            );
            let bedrock = layers.bedrock.as_mut_slice();
            let units = layers.unit_id.as_mut_slice();
            let deep = layers.deep.as_mut_slice();
            let melt = layers.melt.as_mut_slice();
            let cover = layers.cover.as_mut_slice();
            for &(idx, r, phi) in &self.stamper.footprint {
                let d = if drawing.needs_direction() {
                    grid.center(grid.cell(idx))
                } else {
                    p.center
                };
                let (cut, floor) = drawing.cut(r, phi, d, body_radius);
                bedrock[idx] -= cut as f32;
                if floor {
                    units[idx] = UNIT_PATERA_FLOOR;
                    deep[idx] = 1.0;
                    melt[idx] = p.floor_age_yr as f32;
                    cover[idx] = 1.0;
                }
            }
        }
    }
}

/// How a flow's path is traced and its lobe's width drawn (F59): the shortest spacing of
/// its points, how strongly the fall line steers it, over what slope and what length of
/// path, how far it wanders over what length, how far above the lowest ground it has
/// reached it may climb before it ponds, and how much its width wobbles.
#[derive(Clone, Copy, Debug)]
pub struct PathRules {
    pub min_spacing: f64,
    pub slope_follow: f64,
    pub slope_scale: f64,
    pub steer_length: f64,
    pub wander_rad: f64,
    pub wander_length: f64,
    pub fill_depth: f64,
    pub width_wobble: f64,
}

/// A flow's path (F59): from a point along a heading, downslope over the ground with
/// the inertia of a flow and a wander, deflecting along a contour but never turning
/// back, stopping where `blocked` says something stands, where it would have to climb out of the
/// ground it has ponded in, or at its length; a lobe about it this half wide on level
/// ground, widening by `toe` from its source to its end, narrowing where the ground is
/// steep and spreading where it flattens (F101). Returns the lobe and how it ended: 0
/// at its length, 1 at a mountain, 2 ponded.
#[allow(clippy::too_many_arguments)]
pub fn trace_flow(
    grid: &CubeSphere,
    body_radius: f64,
    rules: &PathRules,
    idx: u64,
    draw: u64,
    start: [f64; 3],
    heading: [f64; 3],
    length: f64,
    half_width: f64,
    toe: f64,
    ground: &PaddedField,
    blocked: &dyn Fn([f64; 3]) -> bool,
    age_yr: f64,
) -> (Lobe, usize) {
    let r = body_radius;
    let mut p = start;
    let radial = dot(heading, p);
    let mut h = normalize([
        heading[0] - radial * p[0],
        heading[1] - radial * p[1],
        heading[2] - radial * p[2],
    ]);
    let spacing = (length / 64.0).max(rules.min_spacing);
    let steps = (length / spacing).ceil() as usize;
    let noise = GradientNoise::new(crate::rng::hash_words(&[idx, draw, 0x0F10_7ACE]));
    let delta = (0.5 * spacing / r).max(1e-7);
    let mut points = vec![p];
    let mut slopes = Vec::with_capacity(steps + 1);
    let mut lowest = sample_bilinear(grid, ground, p) as f64;
    let mut end = 0usize;
    for i in 0..steps {
        // The ground's fall line at the point.
        let (e, n) = local_basis(p);
        let at = |u: f64, v: f64| {
            normalize([
                p[0] + delta * (u * e[0] + v * n[0]),
                p[1] + delta * (u * e[1] + v * n[1]),
                p[2] + delta * (u * e[2] + v * n[2]),
            ])
        };
        let hx = (sample_bilinear(grid, ground, at(1.0, 0.0))
            - sample_bilinear(grid, ground, at(-1.0, 0.0))) as f64
            / (2.0 * delta * r);
        let hy = (sample_bilinear(grid, ground, at(0.0, 1.0))
            - sample_bilinear(grid, ground, at(0.0, -1.0))) as f64
            / (2.0 * delta * r);
        let grad = hx.hypot(hy);
        slopes.push(grad);
        let down = if grad > 1e-12 {
            normalize([
                -(hx * e[0] + hy * n[0]),
                -(hx * e[1] + hy * n[1]),
                -(hx * e[2] + hy * n[2]),
            ])
        } else {
            h
        };
        // The fall line turns the front at a rate per distance, over the steering length
        // at full slope; on ground steeper than FLOW_STEEP_SLOPE faster as the slope
        // grows, as gravity's push across the front does (F101): on a shield's flank a
        // breakout swings downhill within a few kilometres.
        let follow = rules.slope_follow
            * (grad / rules.slope_scale).tanh()
            * (grad / FLOW_STEEP_SLOPE).clamp(1.0, FLOW_STEER_MAX)
            * (spacing / rules.steer_length).min(1.0);
        // The fall line steers the flow sideways and speeds it on, but never turns it
        // back: a front deflects along a contour and ponds rather than reverse.
        let along = dot(down, h).max(0.0);
        let lateral = {
            let a = dot(down, h);
            [down[0] - a * h[0], down[1] - a * h[1], down[2] - a * h[2]]
        };
        // The wander: a turn per step that sums to the wander over its length.
        let turn = rules.wander_rad
            * (spacing / rules.wander_length)
            * noise.fbm(
                [i as f64 * spacing / rules.wander_length, 0.3, 0.7],
                1.0,
                2,
                2.0,
                0.5,
            )
            / 0.6;
        let side = normalize(cross(p, h));
        let ahead = 1.0 + follow * along;
        let mut d = [
            ahead * (h[0] * turn.cos() + side[0] * turn.sin()) + follow * lateral[0],
            ahead * (h[1] * turn.cos() + side[1] * turn.sin()) + follow * lateral[1],
            ahead * (h[2] * turn.cos() + side[2] * turn.sin()) + follow * lateral[2],
        ];
        let radial = dot(d, p);
        d = normalize([
            d[0] - radial * p[0],
            d[1] - radial * p[1],
            d[2] - radial * p[2],
        ]);
        let theta = spacing / r;
        let q = normalize([
            p[0] * theta.cos() + d[0] * theta.sin(),
            p[1] * theta.cos() + d[1] * theta.sin(),
            p[2] * theta.cos() + d[2] * theta.sin(),
        ]);
        // Something in the way stops the flow.
        let blocked = blocked(q);
        if blocked {
            end = 1;
            break;
        }
        // Lava ponds: it crosses a bump no higher than its own pile above the lowest
        // ground it has reached, and stops where it would have to climb more.
        let height = sample_bilinear(grid, ground, q) as f64;
        if height > lowest + rules.fill_depth {
            end = 2;
            break;
        }
        lowest = lowest.min(height);
        let radial = dot(d, q);
        h = normalize([
            d[0] - radial * q[0],
            d[1] - radial * q[1],
            d[2] - radial * q[2],
        ]);
        p = q;
        points.push(p);
    }
    // The width along the lobe: a channel at the vent widening to the toe, wandering, and
    // narrowing on steep ground, where the flow runs fast between its levees, as one over
    // the slope past FLOW_NARROWING_SLOPE (Hulme 1974). The lobe's half width is its mean.
    let count = points.len();
    while slopes.len() < count {
        slopes.push(slopes.last().copied().unwrap_or(0.0));
    }
    let narrowing: Vec<f64> = slopes
        .iter()
        .map(|&g| 1.0 / (1.0 + (g / FLOW_NARROWING_SLOPE).powi(2)).sqrt())
        .collect();
    let mean_narrowing = narrowing.iter().sum::<f64>() / count as f64;
    let widths = (0..count)
        .map(|i| {
            let s = i as f64 / (count.max(2) - 1) as f64;
            let wander = 1.0
                + rules.width_wobble
                    * noise.fbm(
                        [i as f64 * spacing / (3.0 * half_width), 1.3, 2.1],
                        1.0,
                        2,
                        2.0,
                        0.5,
                    )
                    / 0.6;
            (wander.max(0.4) * (1.0 - toe + toe * s) * narrowing[i] / mean_narrowing) as f32
        })
        .collect();
    (
        Lobe {
            points,
            spacing,
            half_width: half_width * mean_narrowing,
            widths,
            age_yr,
        },
        end,
    )
}

/// The flow fields of the paterae alive at a time, over the years their lava can show
/// (F117, F201, audits M8 and M20): the lava the heat budget erupts, `rate` m³ a year,
/// shared among the paterae by their floors' areas (a stated choice: Veeder et al. 2012
/// read the paterae's heat from their dark floors, Loki's the largest), each erupting a
/// share IO_PATERA_ACTIVE_FRACTION of the time at the rate that keeps its share. Each
/// patera builds its fields one after another back from the present, a field taking as
/// many years of its share as its volume needs and growing lobe by lobe through them
/// (`volcanism::lobe_cover_laid`): its last lava an active patera's erupting now, a
/// quiet one's a quiet spell ago, memoryless at the duty its eruptions keep, each lasting
/// IO_FLOW_SPAN_YR; until its fields are older than the fallout about the patera takes
/// to hide rough lava IO_FLOW_SHOW_DEPTHS times over, or than the plains, at most
/// IO_FLOW_FIELDS_MAX. Each field is a flood from the patera's rim where the ground falls
/// away, run down over the ground, filling and spilling over the lows it comes to as a
/// flood does on any body (F107, F109) and drawn as its flows at the finish (F94, F110),
/// as large as lays a flood's sheet over the share of its reach squared that a compound
/// field covers (F85), the reach log-uniform (IO_FLOW_REACH); its vent stands past the
/// pit's wall as the loop draws it, so that the lava runs away from the patera rather
/// than back into it. Oldest first.
pub fn eruptions(
    graph: &ParamGraph,
    paterae: &[Patera],
    bedrock: &Field<f32>,
    t: f64,
    rate: f64,
) -> Vec<IoEruption> {
    let grid = bedrock.grid();
    let body_radius = graph.f64("body.radius");
    let rng = Rng::new(graph.f64("seed.volcanism") as u64 ^ 0x10_5EED_0001);
    let log_uniform = |index: u64, draw: u64, range: (f64, f64)| {
        range.0 * (range.1 / range.0).powf(rng.uniform(index, draw))
    };
    let cell = grid.delta() * body_radius;
    let smoothed = crate::stages::relaxation::low_pass(bedrock, cell, IO_FLOW_GROUND_SMOOTH);
    let ground = PaddedField::new(&grid, &smoothed);
    let alive: Vec<&Patera> = paterae.iter().filter(|p| p.alive(t)).collect();
    let floors: f64 = alive.iter().map(|p| p.floor_area()).sum();
    if alive.is_empty() || floors <= 0.0 || rate <= 0.0 {
        return Vec::new();
    }
    // The fallout about each patera, which hides its old flows.
    let far = far_fallout(
        CubeSphere::new(grid.face_size().min(64)),
        body_radius,
        graph.f64("body.gravity"),
        &alive.iter().map(|p| (*p).clone()).collect::<Vec<_>>(),
    );
    let far_mean = far.area_mean();
    let (s_lo, s_hi) = IO_FLOW_SPAN_YR;
    let mean_span = (s_hi - s_lo) / (s_hi / s_lo).ln();
    let duty = IO_PATERA_ACTIVE_FRACTION;
    let quiet = mean_span * (1.0 - duty) / duty;
    let mut out = Vec::new();
    for p in alive {
        let idx = p.id;
        // Its share of the budget, m³ a year, and the rate it pours out while it erupts.
        let output = rate * p.floor_area() / floors;
        let effusion = output / (duty * YEAR);
        let fall = if far_mean > 0.0 {
            let c = far.grid().locate(p.center).0;
            IO_BACKGROUND_FALLOUT * far.get(c) as f64 / far_mean
        } else {
            IO_BACKGROUND_FALLOUT
        };
        let window = (IO_FLOW_SHOW_DEPTHS * IO_LAVA_COAT_DEPTHS / fall.max(1e-12))
            .min(IO_SURFACE_AGE_GYR * 1e9);
        // Its last lava: an active patera's erupting now, a quiet one's a quiet spell ago.
        let mut age_yr = if p.active {
            log_uniform(idx, 46, IO_FLOW_ACTIVE_YR)
        } else {
            -quiet * (1.0 - rng.uniform(idx, 46)).ln()
        };
        for g in 0..IO_FLOW_FIELDS_MAX as u64 {
            if age_yr >= window {
                break;
            }
            let reach = log_uniform(idx, 49 + 90 * g, IO_FLOW_REACH).max(2.0 * p.extent());
            let field = idx ^ ((g + 1) << 48);
            // Its vent on the rim where the ground falls away, past the pit's wall.
            let (east, north) = local_basis(p.center);
            let delta = 0.5 * p.radius / body_radius;
            let at = |u: f64, v: f64| {
                normalize([
                    p.center[0] + delta * (u * east[0] + v * north[0]),
                    p.center[1] + delta * (u * east[1] + v * north[1]),
                    p.center[2] + delta * (u * east[2] + v * north[2]),
                ])
            };
            let gx = (sample_bilinear(&grid, &ground, at(1.0, 0.0))
                - sample_bilinear(&grid, &ground, at(-1.0, 0.0))) as f64;
            let gy = (sample_bilinear(&grid, &ground, at(0.0, 1.0))
                - sample_bilinear(&grid, &ground, at(0.0, -1.0))) as f64;
            let downhill = if gx.hypot(gy) > 1e-9 {
                (-gy).atan2(-gx)
            } else {
                std::f64::consts::TAU * rng.uniform(field, 3)
            };
            let azimuth = downhill + IO_FLOW_SOURCE_SPREAD * rng.normal(field, 4);
            let out_to = p.extent() + patera_wall(p.depth, 2.0 * cell) + cell;
            let volume = IO_FLOW_FIELD_AREA * reach * reach * FLOOD_SHEET_THICKNESS;
            // It grows over as many years of the patera's share as its volume needs.
            let span_yr = volume / output;
            out.push(IoEruption {
                vent: crate::geom::offset(p.center, east, north, out_to / body_radius, azimuth),
                volume,
                age_yr,
                id: field,
                effusion,
                span_yr,
            });
            age_yr += span_yr;
        }
    }
    out.sort_by(|a, b| b.age_yr.total_cmp(&a.age_yr).then(a.id.cmp(&b.id)));
    out
}

/// A flood's cells as Io's lava (F117): the fresh silicate flow, or a patera's floor where
/// it ponded on one (Williams et al. 2011 map a patera's floor as one unit, its lavas
/// with it); its coverage varying with the flow units within its lobes (F85); dark
/// silicate as far as it covers; and its age in years, from which the plumes' fallout
/// coats it (F84): its field's last lava `age_yr` years ago, the field growing over
/// `span_yr` years before, each cell's lobe laid when the flood's pass says (F201). Over
/// lava already laid the stronger coverage wins, and the younger age where the new lava
/// covers the cell well.
pub fn as_io_flow(
    out: &mut StampOutput,
    layers: &LayerSet,
    grid: &CubeSphere,
    body_radius: f64,
    id: u64,
    age_yr: f64,
    span_yr: f64,
) {
    let units = layers.unit_id.as_slice();
    let covered = layers.cover.as_slice();
    let laid = layers.melt.as_slice();
    let noise = GradientNoise::new(id ^ 0xF10_0000);
    for c in out.cells.iter_mut() {
        if !c.surface || c.unit != UNIT_LAVA {
            continue;
        }
        let idx = c.idx as usize;
        let d = grid.center(grid.cell(idx));
        let texture = 1.0
            - IO_FLOW_TEXTURE
                * 0.5
                * (1.0 + noise.fbm(d, body_radius / IO_FLOW_TEXTURE_LENGTH, 2, 2.0, 0.5));
        let cover = (c.cover as f64 * texture).clamp(0.0, 1.0) as f32;
        let age_yr = age_yr + span_yr * (1.0 - (c.melt as f64).clamp(0.0, 1.0));
        let floor = units[idx] == UNIT_PATERA_FLOOR;
        let (cover, age) = if matches!(units[idx], UNIT_IO_FLOW | UNIT_PATERA_FLOOR) {
            (
                cover.max(covered[idx]),
                if cover >= 0.5 {
                    (age_yr as f32).min(laid[idx])
                } else {
                    laid[idx]
                },
            )
        } else {
            (cover, age_yr as f32)
        };
        c.unit = if floor {
            UNIT_PATERA_FLOOR
        } else {
            UNIT_IO_FLOW
        };
        c.cover = cover;
        c.deep = if floor { 1.0 } else { cover };
        c.melt = age;
    }
    out.composition = None;
}

/// A point of a path and its heading there.
pub(crate) fn path_point(points: &[[f64; 3]], k: usize) -> ([f64; 3], [f64; 3]) {
    let k = k.min(points.len() - 1);
    let (a, b) = if k + 1 < points.len() {
        (points[k], points[k + 1])
    } else {
        (points[k.saturating_sub(1)], points[k])
    };
    let p = points[k];
    let d = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let radial = dot(d, p);
    (
        p,
        normalize([
            d[0] - radial * p[0],
            d[1] - radial * p[1],
            d[2] - radial * p[2],
        ]),
    )
}

/// A heading at a point turned by an angle about the vertical there.
pub(crate) fn rotate_about(p: [f64; 3], heading: [f64; 3], angle: f64) -> [f64; 3] {
    let side = normalize(cross(p, heading));
    let (s, c) = angle.sin_cos();
    normalize([
        heading[0] * c + side[0] * s,
        heading[1] * c + side[1] * s,
        heading[2] * c + side[2] * s,
    ])
}

/// The plains buried (F86): the relief of the base at every scale of its Gaussian
/// pyramid, each band shrunk by the burial depth, m, so that only what stands higher
/// than the burial shows; the broadest level whole.
pub fn bury(bedrock: &Field<f32>, cell: f64, body_radius: f64, depth: f64) -> Vec<f32> {
    let levels = crate::stages::relaxation::low_passes(
        bedrock,
        cell,
        0.5 * std::f64::consts::PI * body_radius,
    );
    let mut finer: Vec<f32> = bedrock.as_slice().to_vec();
    let mut out = vec![0.0f32; finer.len()];
    let depth = depth as f32;
    for (_, low) in &levels {
        for k in 0..out.len() {
            let band = finer[k] - low[k];
            out[k] += band.signum() * (band.abs() - depth).max(0.0);
        }
        finer.copy_from_slice(low);
    }
    for k in 0..out.len() {
        out[k] += finer[k];
    }
    out
}

/// Whether the body is in the Io regime: rock, heated past the threshold.
pub fn is_io_regime(graph: &ParamGraph) -> bool {
    graph.str("body.composition_class") == "rock"
        && graph.f64("interior.tidal_flux") > IO_REGIME_TIDAL_FLUX
}

/// Finishing pass of a patera: the pit at this resolution, cut into the buried plains the
/// finish starts from, against the foot of the mountain it sank beside, and the floor's
/// unit, its lava's coverage and the years since it was laid.
pub fn patera_compute(
    layers: &LayerSet,
    grid: &CubeSphere,
    body_radius: f64,
    p: &Patera,
    stamper: &mut Stamper,
) -> StampOutput {
    let bedrock = layers.bedrock.as_slice();
    let cell = grid.delta() * body_radius;
    let drawing = PateraDrawing::new(p, body_radius, cell, cell);
    stamper.gather(grid, body_radius, p.center, p.extent() + 3.0 * drawing.wall);
    let time = layers.surface_age.as_slice();
    let units = layers.unit_id.as_slice();
    let covered = layers.cover.as_slice();
    let patches = GradientNoise::new(p.id ^ 0xF100_A7E5);
    // The floor's share of each cell, from where its edge crosses the cell (F131).
    let (east, north) = crate::geom::local_basis(p.center);
    let inside = |e: [f64; 3]| {
        drawing.edge(crate::geom::azimuth_from(p.center, east, north, e))
            - crate::geom::angle_between(p.center, e) * body_radius
    };
    let cells = stamper
        .footprint
        .iter()
        .map(|&(idx, r, phi)| {
            let d = if drawing.needs_direction() {
                grid.center(grid.cell(idx))
            } else {
                p.center
            };
            let (cut, _) = drawing.cut(r, phi, d, body_radius);
            let share =
                crate::fields::level_share(&inside, grid.center(grid.cell(idx)), 0.0, grid.delta())
                    as f32;
            let mut write = CellWrite {
                idx: idx as u32,
                bedrock: (bedrock[idx] as f64 - cut) as f32,
                surface: false,
                age: 0.0,
                unit: 0,
                deep: 0.0,
                melt: 0.0,
                cover: 0.0,
            };
            // A floor cell carries its lava's share of the cell, over a flow's lava on the
            // rest, as a flow's margin does, so that the deposit history coats each by its
            // share.
            if share > 0.0 {
                let under = if matches!(units[idx], UNIT_IO_FLOW | UNIT_PATERA_FLOOR) {
                    covered[idx]
                } else {
                    0.0
                };
                write.surface = true;
                write.age = time[idx];
                write.unit = UNIT_PATERA_FLOOR;
                write.deep = share;
                write.melt =
                    floor_age(p, grid.center(grid.cell(idx)), &patches, body_radius) as f32;
                write.cover = share + (1.0 - share) * under;
            }
            write
        })
        .collect();
    StampOutput {
        cells,
        composition: None,
        stain: None,
        steep: None,
        ejecta_age: None,
        base: None,
        compositions: None,
    }
}

/// Finishing pass of a mountain (F87): the block at this resolution, its scarp and ridges
/// with it, raised on the buried plains the finish starts from.
pub fn mountain_compute(
    layers: &LayerSet,
    grid: &CubeSphere,
    body_radius: f64,
    m: &Mountain,
    stamper: &mut Stamper,
) -> StampOutput {
    let bedrock = layers.bedrock.as_slice();
    let ages = layers.surface_age.as_slice();
    let melt = layers.melt.as_slice();
    let cover = layers.cover.as_slice();
    let cell = grid.delta() * body_radius;
    stamper.gather(grid, body_radius, m.center, m.reach(cell));
    let cells = stamper
        .footprint
        .iter()
        .map(|&(idx, r, phi)| {
            let (u, v) = m.frame(r, phi);
            let fine = mountain_height(m, u, v, cell);
            let mut write = CellWrite {
                idx: idx as u32,
                bedrock: (bedrock[idx] as f64 + fine) as f32,
                surface: false,
                age: 0.0,
                unit: 0,
                deep: 0.0,
                melt: 0.0,
                cover: 0.0,
            };
            if fine > MOUNTAIN_UNIT_SHARE * m.height {
                write.surface = true;
                write.age = ages[idx];
                write.unit = UNIT_MOUNTAIN;
                write.melt = melt[idx];
                write.cover = cover[idx];
            }
            write
        })
        .collect();
    StampOutput {
        cells,
        composition: None,
        stain: None,
        steep: None,
        ejecta_age: None,
        base: None,
        compositions: None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn profiles_have_the_right_shape() {
        let d = pit_profile(0.0, 20e3, 1e3, 1e3);
        assert!((d - 1e3).abs() < 1.0);
        assert!(pit_profile(30e3, 20e3, 1e3, 1e3) < 1.0);
        // A thrust block (F87): highest near its scarp, lower towards its back, falling to
        // nothing past its tips and well beyond its scarp's foot, its debris below.
        let m = Mountain {
            center: [1.0, 0.0, 0.0],
            strike: 0.0,
            length: 150e3,
            width: 60e3,
            height: 6e3,
            vergence: 1.0,
            tilt: 0.4,
            level: 0.9,
            ridging: 0.0,
            birth_gyr: 0.0,
            death_gyr: 1.0,
            id: 1,
        };
        let cell = 1e3;
        let near = mountain_height(&m, 0.0, 15e3, cell);
        let far = mountain_height(&m, 0.0, 55e3, cell);
        assert!(near > far && near > 3e3, "{near} vs {far}");
        assert!(mountain_height(&m, 0.0, -40e3, cell) < 0.1 * near);
        assert!(mountain_height(&m, 140e3, 20e3, cell) == 0.0);
        assert!(mountain_height(&m, 76e3, 20e3, cell) < 0.3 * near);
        assert!(mountain_height(&m, 70e3, 20e3, cell) < mountain_height(&m, 0.0, 20e3, cell));
        let scarp = mountain_scarp(6e3, cell);
        let apron = mountain_height(&m, 0.0, -mountain_heave(6e3) - scarp - 1e3, cell);
        assert!(apron > 0.0 && apron < 0.5 * near, "{apron}");
        // The frame: along the strike and across towards the block.
        let (u, v) = m.frame(10e3, std::f64::consts::FRAC_PI_2);
        assert!(u.abs() < 1e-6 && (v - 10e3).abs() < 1e-6, "{u} {v}");
        // The ring peaks at the range and ends more steeply outside than it rises
        // inside; the pyroclasts fall within.
        let total = |x: f64| {
            let (ring, fill) = deposit_profile(x, IO_RING_OUTER_WIDTH);
            ring + 0.3 * fill
        };
        assert!(total(1.0) > total(0.5) && total(0.5) > total(1.3));
        let w = IO_RING_OUTER_WIDTH;
        assert!(deposit_profile(1.15, w).0 < deposit_profile(0.85, w).0);
        assert!(deposit_profile(1.1, w).1 == 0.0 && deposit_profile(0.2, w).1 == 1.0);
        let moment = profile_moment(0.3, w);
        assert!(moment > 0.2 && moment < 0.6, "{moment}");
    }

    /// Writes three blocks' heights, a lineated, a plateau and a massif, to the file
    /// EXOGEN_BLOCK_OUT names, for a look at their forms: each 400 by 260 cells of 1 km,
    /// little-endian f32 rows after the counts.
    #[test]
    #[ignore]
    fn write_blocks() {
        let Some(path) = std::env::var_os("EXOGEN_BLOCK_OUT") else {
            return;
        };
        let forms = [
            (0.5, 0.3, 1.0, 5_123),
            (0.0, 0.9, 0.0, 7_001),
            (0.0, 0.0, 0.0, 9_377),
        ];
        let (nx, ny, cell) = (400usize, 260usize, 1e3);
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&(forms.len() as u32).to_le_bytes());
        bytes.extend_from_slice(&(nx as u32).to_le_bytes());
        bytes.extend_from_slice(&(ny as u32).to_le_bytes());
        for &(tilt, level, ridging, id) in &forms {
            let m = Mountain {
                center: [1.0, 0.0, 0.0],
                strike: 0.0,
                length: 260e3,
                width: 90e3,
                height: 5e3,
                vergence: 1.0,
                tilt,
                level,
                ridging,
                birth_gyr: 0.0,
                death_gyr: 1.0,
                id,
            };
            for row in 0..ny {
                let v = (0.65 * ny as f64 - row as f64) * cell;
                for col in 0..nx {
                    let u = (col as f64 - nx as f64 / 2.0) * cell;
                    let h = mountain_height(&m, u, v, 2.8e3) as f32;
                    bytes.extend_from_slice(&h.to_le_bytes());
                }
            }
        }
        std::fs::write(path, bytes).unwrap();
    }

    /// A patera's outline (F88): the union of its nested collapses, cut straight by its
    /// faults.
    #[test]
    fn a_patera_outline_is_its_collapses_and_faults() {
        let p = Patera {
            center: [1.0, 0.0, 0.0],
            radius: 30e3,
            depth: 1e3,
            active: false,
            floor_age_yr: 100.0,
            pele: false,
            axis: 0.0,
            ellipticity: 1.0,
            collapses: vec![[0.0, 0.0, 30e3], [15e3, 0.0, 20e3]],
            faults: vec![[std::f64::consts::FRAC_PI_2, 20e3]],
            abuts: None,
            birth_gyr: 0.0,
            death_gyr: 1.0,
            id: 3,
        };
        // Along the axis the second collapse reaches farther; the fault cuts the north
        // straight. The margin wanders off the collapses' arcs (F118), not off the fault.
        let arc = |phi: f64, p: &Patera| patera_radius(phi, p) / margin_wander(p, phi);
        assert!((arc(0.0, &p) - 35e3).abs() < 1.0, "{}", arc(0.0, &p));
        assert!((arc(std::f64::consts::PI, &p) - 30e3).abs() < 1.0);
        assert!((patera_radius(std::f64::consts::FRAC_PI_2, &p) - 20e3).abs() < 1.0);
        assert!(arc(-std::f64::consts::FRAC_PI_2, &p) > 29e3);
        for k in 0..360 {
            let phi = (k as f64).to_radians();
            assert!(patera_radius(phi, &p) <= p.extent());
            assert!((margin_wander(&p, phi) - 1.0).abs() < 4.0 * IO_PATERA_WANDER);
        }
        // An elongated collapse keeps its area: twice as long along its axis as the
        // circle of its area, half as wide across it (F90).
        let oval = Patera {
            ellipticity: 4.0,
            collapses: vec![[0.0, 0.0, 30e3]],
            faults: Vec::new(),
            axis: std::f64::consts::FRAC_PI_4,
            ..p.clone()
        };
        let along = arc(std::f64::consts::FRAC_PI_4, &oval);
        let across = arc(3.0 * std::f64::consts::FRAC_PI_4, &oval);
        assert!(
            (along - 60e3).abs() < 1.0 && (across - 15e3).abs() < 1.0,
            "{along} {across}"
        );
        assert!((oval.extent() - 60e3 * (1.0 + 4.0 * IO_PATERA_WANDER)).abs() < 1.0);
    }

    /// Io's small eruptions (F89): as many as their rate gives above the grid, their sizes
    /// the same east and west; and the background is the plumes' far fallout, heaviest
    /// among the volcanic centres, its mean the background's.
    #[test]
    fn small_eruptions_and_the_background() {
        let cell = 20e3;
        let small = small_eruptions(7, cell);
        let (r_min, r_max) = IO_SMALL_ERUPTION_RANGE;
        let want = IO_SMALL_ERUPTIONS_PER_YR * history_years() * (r_min / (2.0 * cell)).powi(2);
        assert!(
            (small.len() as f64 / want - 1.0).abs() < 0.05,
            "{} vs {want}",
            small.len()
        );
        let side = |east: bool| {
            let sizes: Vec<f64> = small
                .iter()
                .filter(|e| (e.center[1] >= 0.0) == east)
                .map(|e| e.radius)
                .collect();
            sizes.iter().sum::<f64>() / sizes.len() as f64
        };
        let (east, west) = (side(true), side(false));
        assert!((east / west - 1.0).abs() < 0.08, "{east} vs {west}");
        assert!(
            small
                .iter()
                .all(|e| e.radius >= 2.0 * cell && e.radius <= r_max)
        );
        // The background about one centre.
        let grid = CubeSphere::new(32);
        let patera = Patera {
            center: [1.0, 0.0, 0.0],
            radius: 20e3,
            depth: 1e3,
            active: true,
            floor_age_yr: 0.0,
            pele: false,
            axis: 0.0,
            ellipticity: 1.0,
            collapses: vec![[0.0, 0.0, 20e3]],
            faults: Vec::new(),
            abuts: None,
            birth_gyr: 0.0,
            death_gyr: 1.0,
            id: 1,
        };
        let fall = background_fall(grid, 1_821.6e3, 1.796, &[patera], 3);
        let field = Field::from_fn(grid, |c, _| fall[grid.index(c)]);
        assert!((field.area_mean() / IO_BACKGROUND_FALLOUT - 1.0).abs() < 0.01);
        let at = |d: [f64; 3]| fall[grid.index(grid.locate(d).0)];
        assert!(at([1.0, 0.0, 0.0]) > 3.0 * at([-1.0, 0.0, 0.0]));
    }

    /// A patera that sank beside a mountain (F90): its margin towards the block is the
    /// block's foot, and away from it its own outline.
    #[test]
    fn a_patera_beside_a_mountain_ends_at_its_foot() {
        let r = 1_821e3;
        let m = Mountain {
            center: [1.0, 0.0, 0.0],
            strike: 0.0,
            length: 150e3,
            width: 60e3,
            height: 5e3,
            vergence: 1.0,
            tilt: 0.0,
            level: 0.9,
            ridging: 0.0,
            birth_gyr: 0.0,
            death_gyr: 1.0,
            id: 11,
        };
        // South of the block's front, which faces south.
        let (east, north) = local_basis(m.center);
        let at = crate::geom::offset(
            m.center,
            east,
            north,
            25e3 / r,
            -std::f64::consts::FRAC_PI_2,
        );
        let p = Patera {
            center: at,
            radius: 20e3,
            depth: 1e3,
            active: false,
            floor_age_yr: 100.0,
            pele: false,
            axis: 0.0,
            ellipticity: 1.0,
            collapses: vec![[0.0, 0.0, 20e3]],
            faults: Vec::new(),
            abuts: Some(m.clone()),
            birth_gyr: 0.0,
            death_gyr: 1.0,
            id: 12,
        };
        let edges = foot_edges(&p, &m, r, 1e3, 500.0);
        let towards = foot_edge(&edges, std::f64::consts::FRAC_PI_2);
        let away = foot_edge(&edges, -std::f64::consts::FRAC_PI_2);
        assert!(towards < 20e3 && towards > 2e3, "{towards}");
        assert!(away.is_infinite(), "{away}");
        // The block stands above its share there, and not short of it.
        let d = |t: f64| {
            crate::geom::offset(
                at,
                local_basis(at).0,
                local_basis(at).1,
                t / r,
                std::f64::consts::FRAC_PI_2,
            )
        };
        assert!(m.height_at(d(towards + 300.0), r, 1e3) > MOUNTAIN_UNIT_SHARE * m.height);
        assert!(m.height_at(d(towards - 600.0), r, 1e3) <= MOUNTAIN_UNIT_SHARE * m.height);
    }

    /// The top layer (F84): a red deposit stays red while it falls and reverts within a
    /// few years after; a thick frost outlasts a thin one; the frost survives longer in
    /// the cold; lava is hidden by its coat; and the closed form agrees with small steps.
    #[test]
    fn the_top_layer_ages_and_is_buried() {
        let warm = Kinetics {
            revert: 1.0 / IO_RED_FADE_YR,
            sublime: 60.0,
            alter: 0.0,
        };
        let cold = Kinetics {
            sublime: 0.6,
            ..warm
        };
        let red_fall = Fall {
            solid: 4.0,
            red: 4.0,
            ..Fall::default()
        };
        let mut top = Top::default();
        top.age(10.0, red_fall, warm);
        assert!(top.red > 0.85, "red while it falls: {}", top.red);
        top.age(10.0, Fall::default(), warm);
        assert!(top.red < 0.01, "reverted after: {}", top.red);
        // Frost from a plume: a long eruption's column outlasts a short one's.
        let frost_fall = Fall {
            solid: 20.0,
            frost: 140.0,
            ..Fall::default()
        };
        let lasting = |years: f64, k: Kinetics| {
            let mut t = Top::default();
            t.age(years, frost_fall, k);
            let mut after = 0.0f64;
            while t.visible_frost() > 0.5 && after < 1e4 {
                t.age(1.0, Fall::default(), k);
                after += 1.0;
            }
            after
        };
        let short = lasting(1.0, warm);
        let long = lasting(10.0, warm);
        assert!(long > 3.0 * short && long > 5.0, "{short} vs {long}");
        assert!(lasting(1.0, cold) > 10.0 * short.max(1.0));
        // Lava hides under its coat, slower than a smooth plain: once its hollows fill,
        // some IO_LAVA_COAT_DEPTHS of fall (F201).
        let mut lava = Top {
            lava: 1.0,
            ..Top::default()
        };
        lava.age(1.0, red_fall, warm);
        assert!(lava.lava > 0.8, "{}", lava.lava);
        lava.age(0.5 * IO_LAVA_COAT_DEPTHS / red_fall.solid, red_fall, warm);
        assert!(lava.lava > 0.5, "{}", lava.lava);
        lava.age(5.0 * IO_LAVA_COAT_DEPTHS / red_fall.solid, red_fall, warm);
        assert!(lava.lava < 0.01, "{}", lava.lava);
        // A fall of nothing but rounding leaves the cover over the frost as it was.
        let mut covered = Top {
            lag: 14.0,
            ..Top::default()
        };
        covered.age(
            100.0,
            Fall {
                frost: 1e-18,
                ..Fall::default()
            },
            warm,
        );
        assert!((covered.lag - 14.0).abs() < 1e-9, "{}", covered.lag);
        // Closed form against small steps.
        let fall = Fall {
            solid: 2.0,
            red: 0.5,
            dark: 0.3,
            frost: 30.0,
        };
        let k = Kinetics {
            revert: 0.5,
            sublime: 20.0,
            alter: 0.05,
        };
        let mut one = Top {
            lava: 0.7,
            dark: 0.2,
            red: 0.1,
            altered: 0.4,
            frost: [3.0; 16],
            lag: 2.0,
            fine: 1.0,
        };
        let mut many = one;
        one.age(5.0, fall, k);
        for _ in 0..5000 {
            many.age(1e-3, fall, k);
        }
        let close = |a: f64, b: f64| (a - b).abs() < 1e-6 * (1.0 + a.abs());
        assert!(close(one.red, many.red) && close(one.altered, many.altered));
        assert!(close(one.frost[1], many.frost[1]) && close(one.lag, many.lag));
        assert!(close(one.lava, many.lava) && close(one.dark, many.dark));
        assert!(close(one.fine, many.fine));
        // A frost coarsens as it ages and shows less for its mass (F203): a thick column
        // laid in the cold shows whole while it is fresh, and decades later, unsublimed,
        // lets the ground through.
        let frozen = Kinetics {
            sublime: 0.0,
            ..warm
        };
        let mut aged = Top::default();
        aged.age(
            1.0,
            Fall {
                frost: 5.0,
                ..Fall::default()
            },
            frozen,
        );
        let fresh = aged.visible_frost();
        aged.age(30.0, Fall::default(), frozen);
        assert!(fresh > 0.95, "{fresh}");
        assert!(aged.visible_frost() < 0.4, "{}", aged.visible_frost());
    }

    /// Io's eruptions (F117, F201): oldest first, each from past its patera's pit on the
    /// side the ground falls away, as large as a flood's sheet over the field's share of
    /// its reach squared; each patera's fields back to back from its last lava, an active
    /// patera's within its first years, each growing over the years its share of the
    /// budget takes to erupt it, at the rate that keeps the share at the paterae's duty.
    #[test]
    fn io_s_eruptions_leave_their_paterae_downhill() {
        let (g, _) = crate::presets::preset("Io").unwrap().to_graph();
        let grid = CubeSphere::new(128);
        let r = g.f64("body.radius");
        // Ground falling away along the body's y axis.
        let bedrock = Field::from_fn(grid, |_, d| (-2000.0 * d[1]) as f32);
        let center = normalize([1.0, 0.0, 0.0]);
        let patera = |active: bool, id: u64| Patera {
            center,
            radius: 15e3,
            depth: 1e3,
            active,
            floor_age_yr: 0.0,
            pele: false,
            axis: 0.0,
            ellipticity: 1.0,
            collapses: vec![[0.0, 0.0, 15e3]],
            faults: Vec::new(),
            abuts: None,
            birth_gyr: 0.0,
            death_gyr: 10.0,
            id,
        };
        let paterae: Vec<Patera> = (0..40).map(|k| patera(k % 2 == 0, 7 + k)).collect();
        let rate = 4.8e10;
        let all = eruptions(&g, &paterae, &bedrock, 4.5, rate);
        assert!(all.len() > 60, "{} eruptions", all.len());
        // Each patera's share, the same floors here, and its fields back to back.
        let output = rate / paterae.len() as f64;
        let mut by_patera: std::collections::HashMap<u64, Vec<&IoEruption>> =
            std::collections::HashMap::new();
        for e in &all {
            assert!((e.span_yr - e.volume / output).abs() < 1e-6 * e.span_yr);
            assert!((e.effusion - output / (IO_PATERA_ACTIVE_FRACTION * YEAR)).abs() < 1e-6);
            by_patera.entry(e.id & 0xFFFF_FFFF).or_default().push(e);
        }
        for (id, mut fields) in by_patera {
            fields.sort_by(|a, b| a.age_yr.total_cmp(&b.age_yr));
            if id % 2 == 1 {
                assert!(fields[0].age_yr <= IO_FLOW_ACTIVE_YR.1);
            }
            for w in fields.windows(2) {
                let end = w[0].age_yr + w[0].span_yr;
                assert!(
                    (w[1].age_yr - end).abs() < 1e-6 * end,
                    "{id}: {end} then {}",
                    w[1].age_yr
                );
            }
            assert!(fields.len() <= IO_FLOW_FIELDS_MAX);
        }
        for w in all.windows(2) {
            assert!(w[0].age_yr >= w[1].age_yr);
        }
        let cell = grid.delta() * r;
        let (lo, hi) = IO_FLOW_REACH;
        let per_area = IO_FLOW_FIELD_AREA * FLOOD_SHEET_THICKNESS;
        let mut downhill = 0;
        for e in &all {
            let dist = angle_between(center, e.vent) * r;
            assert!(dist > 15e3 + cell, "a vent {dist} m out");
            if e.vent[1] > center[1] {
                downhill += 1;
            }
            assert!(e.volume >= per_area * lo * lo * 0.99 && e.volume <= per_area * hi * hi * 1.01);
            assert!(e.age_yr >= 0.0 && e.age_yr <= IO_SURFACE_AGE_GYR * 1e9);
        }
        assert!(
            downhill as f64 > 0.8 * all.len() as f64,
            "{downhill} of {} downhill",
            all.len()
        );
    }

    /// The plains buried (F86): relief lower than the burial depth goes at every scale,
    /// while the body's broadest undulation stands.
    #[test]
    fn burial_fills_the_small_relief_and_keeps_the_broad() {
        let grid = CubeSphere::new(64);
        let r = 1.8216e6;
        let noise = GradientNoise::new(5);
        let relief = Field::from_fn(grid, |_, d| {
            (3000.0 * d[2] + 150.0 * noise.sample([d[0] * 30.0, d[1] * 30.0, d[2] * 30.0])) as f32
        });
        let buried = bury(&relief, grid.delta() * r, r, 1000.0);
        let roughness = |h: &[f32]| {
            let mut sum = 0.0;
            for k in 0..h.len() {
                let n = grid.index(grid.neighbor(grid.cell(k), 1, 0));
                sum += (h[k] - h[n]).abs() as f64;
            }
            sum / h.len() as f64
        };
        let before = roughness(relief.as_slice());
        let after = roughness(&buried);
        assert!(after < 0.3 * before, "{after} vs {before}");
        let (mut north, mut south, mut n_count, mut s_count) = (0.0, 0.0, 0.0, 0.0);
        for (k, h) in buried.iter().enumerate() {
            let z = grid.center(grid.cell(k))[2];
            if z > 0.7 {
                north += *h as f64;
                n_count += 1.0;
            } else if z < -0.7 {
                south += *h as f64;
                s_count += 1.0;
            }
        }
        let broad = north / n_count - south / s_count;
        assert!(broad > 2000.0, "{broad}");
    }

    /// The plumes at the flows' fronts (F202): lava laid some years ago throws its SO2 in
    /// a ring about it, within the bin of years it was laid in, the mass kept, the fall
    /// heaviest near the class's range from the lava rather than on it.
    #[test]
    fn a_young_flow_throws_its_so2_in_a_ring() {
        let grid = CubeSphere::new(128);
        let r = 1_821.6e3;
        let g = 1.796;
        let n = grid.cell_count();
        let source = normalize([1.0, 0.2, 0.1]);
        let (mut units, mut coverage, mut laid) = (vec![0u16; n], vec![0f32; n], vec![0f32; n]);
        let mut area = 0.0;
        for k in 0..n {
            let d = grid.center(grid.cell(k));
            if angle_between(d, source) * r < 20e3 {
                units[k] = UNIT_IO_FLOW;
                coverage[k] = 1.0;
                laid[k] = 5.0;
                area += grid.cell_area(grid.cell(k)) * r * r;
            }
        }
        let window = history_years();
        let fronts = front_falls(grid, r, g, &units, &coverage, &laid, window);
        let bin = fronts
            .iter()
            .position(|f| -f.start > 5.0 && -f.end <= 5.0)
            .unwrap();
        let mut fallen = 0.0;
        for (b, f) in fronts.iter().enumerate() {
            let years = f.end - f.start;
            let sum: f64 = (0..n)
                .map(|k| {
                    (f.ring[k] + f.fill[k] + f.gas[k]) as f64 * grid.cell_area(grid.cell(k)) * r * r
                })
                .sum();
            if b != bin {
                assert!(sum == 0.0, "bin {b} took {sum}");
            }
            fallen += sum * years * IO_DEPOSIT_OPTICAL_MASS;
        }
        let thrown = area * FLOOD_SHEET_THICKNESS * IO_FRONT_SO2;
        // The gas that does not condense in the canopy settles across the region (F205).
        let settled: f64 = fronts[bin]
            .gas
            .iter()
            .enumerate()
            .map(|(k, &g)| g as f64 * grid.cell_area(grid.cell(k)) * r * r)
            .sum::<f64>()
            * (fronts[bin].end - fronts[bin].start)
            * IO_DEPOSIT_OPTICAL_MASS;
        assert!(
            (settled / thrown - (1.0 - IO_RING_GAS_SHARE)).abs() < 0.05,
            "{settled:e} of {thrown:e} kg settled"
        );
        assert!(
            (fallen / thrown - 1.0).abs() < 0.1,
            "{fallen:e} of {thrown:e} kg"
        );
        // The ring's crest lies near the range, the lava itself under the thin fill.
        let range = IO_PLUME_SPEED_PROMETHEUS * IO_PLUME_SPEED_PROMETHEUS / g;
        let f = &fronts[bin];
        let (mut at_source, mut crest, mut crest_r) = (0.0f32, 0.0f32, 0.0);
        for k in 0..n {
            let d = angle_between(grid.center(grid.cell(k)), source) * r;
            let fall = f.ring[k] + f.fill[k];
            if d < 10e3 {
                at_source = at_source.max(fall);
            }
            if fall > crest {
                crest = fall;
                crest_r = d;
            }
        }
        assert!(
            crest > 5.0 * at_source,
            "crest {crest} at the source {at_source}"
        );
        assert!(
            crest_r > 0.5 * range && crest_r < 1.8 * range,
            "crest at {crest_r} m, range {range} m"
        );
    }

    /// A Pele-class plume rises from the cracks along its lava lake's margin (F208): an
    /// elongated patera's cracks face across its long axis, so its ring reaches farther
    /// across than along; and its pyroclasts jet along the cracks, their mean over every
    /// azimuth kept.
    #[test]
    fn a_pele_class_ring_follows_its_lake_s_margin() {
        let patera = Patera {
            center: [1.0, 0.0, 0.0],
            radius: 20e3,
            depth: 1e3,
            active: true,
            floor_age_yr: 1.0,
            pele: true,
            axis: 0.0,
            ellipticity: 4.0,
            collapses: vec![[0.0, 0.0, 20e3]],
            faults: Vec::new(),
            abuts: None,
            birth_gyr: 0.0,
            death_gyr: 10.0,
            id: 11,
        };
        let rng = Rng::new(5);
        let facing = margin_facing(&patera, 400, &rng, 3);
        let range = |phi: f64| {
            1.0 + facing
                .range
                .iter()
                .enumerate()
                .map(|(j, h)| {
                    let k = (RANGE_ORDERS.0 + j) as f64;
                    h[0] * (k * phi).cos() - h[1] * (k * phi).sin()
                })
                .sum::<f64>()
        };
        let (along, across) = (range(0.0), range(0.5 * std::f64::consts::PI));
        assert!(across > along + 0.1, "across {across:.3}, along {along:.3}");
        let episode = PlumeEpisode {
            center: patera.center,
            radius: 500e3,
            pele: true,
            start_yr: 10.0,
            end_yr: 0.0,
            crest: 1.0,
            aspect: 1.0,
            axis: 0.0,
            dark: IO_DARK_PELE,
            frost: IO_PELE_FROST,
            range: facing.range,
            mottle: [0.0; 3],
            jets: facing.jets[..12].to_vec(),
            small: false,
            id: 1,
        };
        let plume = Plume::new(&episode, 1.82e6, history_years());
        let n = 3600;
        let shares: Vec<f64> = (0..n)
            .map(|i| plume.jet_at(std::f64::consts::TAU * i as f64 / n as f64))
            .collect();
        let mean = shares.iter().sum::<f64>() / n as f64;
        let peak = shares.iter().cloned().fold(0.0, f64::max);
        assert!((mean - 1.0).abs() < 0.02, "mean {mean:.3}");
        assert!(peak > 1.3, "peak {peak:.3}");
    }

    /// The history of a patera's plumes (F84): an active Pele-class patera with a plume
    /// erupting now, the rest at their class's duty, every eruption within the history;
    /// the Prometheus class throws none from its centre, its plumes rising at the flows'
    /// fronts (F202).
    #[test]
    fn plumes_erupt_at_their_class_duty() {
        let mut paterae = Vec::new();
        for k in 0..400u64 {
            paterae.push(Patera {
                center: [1.0, 0.0, 0.0],
                radius: 20e3,
                depth: 1e3,
                active: k % 5 == 0,
                floor_age_yr: 10.0,
                pele: k % 20 == 0,
                axis: 0.0,
                ellipticity: 1.0,
                collapses: vec![[0.0, 0.0, 20e3]],
                faults: Vec::new(),
                abuts: None,
                birth_gyr: 0.0,
                death_gyr: 1.0,
                id: 7_000 + k,
            });
        }
        let episodes = plume_history(&paterae, 1.8, 3);
        let window = history_years();
        let mut erupting = [0usize; 2];
        let mut busy = [0.0f64; 2];
        for e in &episodes {
            assert!(e.start_yr > e.end_yr && e.end_yr >= 0.0);
            assert!(e.end_yr < window);
            let class = e.pele as usize;
            if e.end_yr == 0.0 {
                erupting[class] += 1;
            }
            busy[class] += e.start_yr.min(window) - e.end_yr;
        }
        let duty = |class: usize, count: f64| busy[class] / (count * window);
        let pele = duty(1, 20.0);
        let want = IO_PATERA_ACTIVE_FRACTION;
        assert!(
            (pele / (want * IO_PLUME_FRACTION_PELE) - 1.0).abs() < 0.5,
            "Pele-class duty {pele}"
        );
        assert_eq!(busy[0], 0.0, "Prometheus-class plumes at their centres");
        assert!(erupting[1] > 3, "{erupting:?}");
        // Pele-class rings reach hundreds of kilometres (Geissler et al. 2004).
        let mean = episodes.iter().map(|e| e.radius).sum::<f64>() / episodes.len() as f64;
        assert!(mean > 400e3, "{mean}");
    }
}
