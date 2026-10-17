// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! Plume loading on a stagnant lid (PLAN.md 4.3): a volcanic rise with its graben,
//! canyon and shields.
//!
//! A mantle plume under a stagnant lid builds a broad dome, Tharsis-like, over the first
//! couple of billion years. The flexure around the load puts the flanks in radial
//! extension, cut by sets of radial graben, and the periphery in concentric compression,
//! which the ridge layout picks up through the load field. On a large rise one graben
//! grows into a canyon system by collapse of its walls: a chain of overlapping troughs
//! with scalloped walls and a floor of landslide debris. Shield volcanoes then grow on
//! the rise, their size set by the magma supply and the relief the lithosphere can hold,
//! with nested summit calderas; they keep growing to the present, so their flanks are
//! young. The loop carries the dome, the canyon and the shields on the coarse grid; the
//! finishing pass restores the graben, the canyon walls and the caldera detail.

use std::collections::HashMap;

use crate::constants::*;
use crate::crater::{CellWrite, StampOutput, Stamper};
use crate::fields::LayerSet;
use crate::geom::{angle_between, azimuth_from, cross, dot, local_basis, normalize, offset};
use crate::grid::CubeSphere;
use crate::grid::projection::PaddedField;
use crate::noise::GradientNoise;
use crate::params::ParamGraph;
use crate::rng::Rng;
use crate::stages::io::{Lobe, PathRules};
use crate::units::{UNIT_LAVA, UNIT_TESSERA};

pub const VERSION: u64 = 5;

/// A volcanic rise: a broad dome growing between two times.
#[derive(Clone, Debug)]
pub struct Rise {
    pub center: [f64; 3],
    pub radius: f64,
    pub height: f64,
    pub start_gyr: f64,
    pub end_gyr: f64,
}

/// A shield volcano on a rise, growing from its start to the present.
#[derive(Clone, Debug)]
pub struct Shield {
    pub center: [f64; 3],
    pub radius: f64,
    pub height: f64,
    pub start_gyr: f64,
    /// Its summit's collapses, oldest first (F102).
    pub calderas: Vec<Caldera>,
    /// The azimuth of its rift zone, rad from east towards north at its centre, along
    /// which its flank vents open (F101).
    pub rift: f64,
    /// Its edifice's axis along the rift over its axis across (F101).
    pub elongation: f64,
    pub id: u64,
}

/// One collapse of a shield's summit (F102): a pit sunk over the magma chamber as an
/// eruption drained it, about as wide as the chamber.
#[derive(Clone, Copy, Debug)]
pub struct Caldera {
    /// Its centre's offset from the shield's, east and north, m.
    pub east: f64,
    pub north: f64,
    /// The radius of the circle of its area, m; its long axis over its short, and the
    /// long axis's azimuth, rad from east towards north.
    pub radius: f64,
    pub ellipticity: f64,
    pub axis: f64,
    /// Its floor's depth below the summit, m.
    pub depth: f64,
}

/// One eruption that shows on a shield (F101): its flow field's lobes, when it erupted,
/// its lava's composition and whether it came from the summit.
pub struct ShieldFlow {
    pub lobes: Vec<Lobe>,
    pub time_gyr: f64,
    pub composition: f32,
    pub summit: bool,
}

/// A trough along a trace: graben on the flanks and the canyon.
#[derive(Clone, Debug)]
pub struct Trough {
    pub points: Vec<[f64; 3]>,
    pub spacing: f64,
    pub length: f64,
    pub width: f64,
    pub depth: f64,
    /// The share of the depth the flanks stand raised by, the faults' footwalls rebounding
    /// (F69's ring graben); none for the rifts and canyons.
    pub flank: f64,
    /// Formation, or for the canyon the start and end of its growth.
    pub start_gyr: f64,
    pub end_gyr: f64,
    pub id: u64,
}

impl Trough {
    pub fn midpoint(&self) -> [f64; 3] {
        self.points[self.points.len() / 2]
    }

    pub fn cap_radius(&self, body_radius: f64) -> f64 {
        (0.5 * self.length + 1.5 * self.width) / body_radius
    }
}

impl Shield {
    pub fn cap_radius(&self, body_radius: f64) -> f64 {
        1.05 * SHIELD_FLOW_REACH * self.extent() / body_radius
    }

    /// A point's distance from the summit, m, and its azimuth there, rad from east
    /// towards north.
    pub fn polar(&self, d: [f64; 3], body_radius: f64) -> (f64, f64) {
        let (east, north) = local_basis(self.center);
        (
            angle_between(self.center, d) * body_radius,
            azimuth_from(self.center, east, north, d),
        )
    }

    /// The edifice's outline toward an azimuth (F101): the radius of its core, longer
    /// along its rift zone, whose fissures feed it more lava there, by its elongation with
    /// the area kept; and the lobing of its foot, where its flows ran out farther in some
    /// directions than others: harmonics of orders SHIELD_OUTLINE_ORDERS, their amplitudes
    /// the order to the -0.75, SHIELD_OUTLINE_LOBING rms.
    fn outline(&self, phi: f64) -> (f64, f64) {
        let a = phi - self.rift;
        let e = self.elongation;
        let core = self.radius / (a.cos().powi(2) / e + e * a.sin().powi(2)).sqrt();
        let rng = Rng::new(self.id ^ 0x0071_1E5E_0000);
        let (lo, hi) = SHIELD_OUTLINE_ORDERS;
        let (mut sum, mut power) = (0.0, 0.0);
        for k in lo..=hi {
            let amplitude = (k as f64).powf(-0.75);
            let phase = 2.0 * std::f64::consts::PI * rng.uniform(k as u64, 0);
            sum += amplitude * (k as f64 * phi + phase).cos();
            power += 0.5 * amplitude * amplitude;
        }
        let lobing = (SHIELD_OUTLINE_LOBING * sum / power.sqrt())
            .clamp(-3.0 * SHIELD_OUTLINE_LOBING, 3.0 * SHIELD_OUTLINE_LOBING);
        (core, lobing)
    }

    /// A point's share of the way from the summit to the edifice's foot (F101): the
    /// core's ellipse on the upper flanks, the lobing coming in as the square of the
    /// distance, so that the summit stays round and the foot is lobed.
    fn reach_share(&self, r: f64, phi: f64) -> f64 {
        let (core, lobing) = self.outline(phi);
        let u = r / core;
        u / (1.0 + lobing * u * u)
    }

    /// The edifice's radius toward an azimuth, to its foot, m.
    pub fn radius_toward(&self, phi: f64) -> f64 {
        let (core, lobing) = self.outline(phi);
        if lobing.abs() < 1e-9 {
            core
        } else {
            core * (1.0 - (1.0 - 4.0 * lobing).sqrt()) / (2.0 * lobing)
        }
    }

    /// The farthest the edifice reaches from its summit, m.
    pub fn extent(&self) -> f64 {
        self.radius * self.elongation.sqrt() * (1.0 + 5.0 * SHIELD_OUTLINE_LOBING)
    }
}

/// Dome shape, 1 at the centre to 0 at the radius.
pub fn dome(x: f64) -> f64 {
    if x >= 1.0 { 0.0 } else { (1.0 - x * x).powi(2) }
}

/// Shield flank shape, 1 at the centre to 0 at the radius: gentle near the summit,
/// steepest low on the flank.
pub fn shield_shape(x: f64) -> f64 {
    if x >= 1.0 {
        0.0
    } else {
        (1.0 - x.powf(SHIELD_FLANK_EXPONENT)).powf(SHIELD_FLANK_POWER)
    }
}

/// Smooth step from 0 at `a` to 1 at `b`.
fn smooth(t: f64, a: f64, b: f64) -> f64 {
    let x = ((t - a) / (b - a)).clamp(0.0, 1.0);
    x * x * (3.0 - 2.0 * x)
}

/// Growth of a rise's dome by a time, 0 to 1.
pub fn rise_progress(rise: &Rise, t: f64) -> f64 {
    smooth(t, rise.start_gyr, rise.end_gyr)
}

/// Growth of a trough by a time, 0 to 1.
pub fn trough_progress(trough: &Trough, t: f64) -> f64 {
    smooth(t, trough.start_gyr, trough.end_gyr)
}

/// Shield growth: most of it early, a tail to the present.
pub fn shield_progress(t: f64, start: f64) -> f64 {
    if t <= start {
        0.0
    } else {
        1.0 - (-(t - start) / SHIELD_GROWTH_GYR).exp()
    }
}

/// Caldera depth at a point, metres below the summit surface (F102): each collapse a pit
/// with a flat floor at its own depth inside a scarp, its outline scalloped by the slumps
/// of its wall; a younger pit, deeper, cuts the older ones it overlaps.
fn caldera_depth(shield: &Shield, d: [f64; 3], body_radius: f64, blur: f64) -> f64 {
    let (east, north) = local_basis(shield.center);
    let v = [
        d[0] - shield.center[0],
        d[1] - shield.center[1],
        d[2] - shield.center[2],
    ];
    let (x, y) = (dot(v, east) * body_radius, dot(v, north) * body_radius);
    let mut depth: f64 = 0.0;
    for (k, c) in shield.calderas.iter().enumerate() {
        let (dx, dy) = (x - c.east, y - c.north);
        let (ca, sa) = (c.axis.cos(), c.axis.sin());
        let (u, w) = (dx * ca + dy * sa, dy * ca - dx * sa);
        let long = c.ellipticity.sqrt();
        let r = ((u / long).powi(2) + (w * long).powi(2)).sqrt();
        let scallop = crate::crater::stamp::crenulation(
            shield.id.wrapping_mul(0x9E37_79B9) ^ k as u64,
            w.atan2(u),
        );
        let edge = c.radius * (1.0 + CALDERA_SCALLOP * scallop);
        let wall = (CALDERA_WALL_WIDTH * c.radius).hypot(blur);
        depth = depth.max(c.depth * 0.5 * (1.0 - ((r - edge) / wall).tanh()));
    }
    depth
}

/// A shield's summit collapses (F102): each time an eruption drains the chamber the roof
/// sinks into it, a pit about as wide as the chamber; the vent wanders between collapses
/// and the chamber shrinks as the supply wanes, so each pit is offset from the last,
/// smaller, and sunk deeper into the floors the later lava ponded in. Olympus Mons's
/// summit holds six, the youngest the smallest and the deepest (Mouginis-Mark et al.
/// 1991; Mouginis-Mark 2021).
pub fn summit_collapses(rng: &Rng, draw: u64, radius: f64, height: f64) -> Vec<Caldera> {
    let tau = 2.0 * std::f64::consts::PI;
    let n = (1 + (CALDERA_COUNT_MAX as f64 * rng.uniform(draw, 0).powi(2)).floor() as usize)
        .min(CALDERA_COUNT_MAX);
    let deepest = (CALDERA_DEPTH_FRACTION * height).min(CALDERA_MAX_DEPTH);
    let (c_lo, c_hi) = CALDERA_RADIUS;
    let mut size = radius * c_lo * (c_hi / c_lo).powf(rng.uniform(draw, 1));
    let mut depth = deepest * (0.35 + 0.25 * rng.uniform(draw, 2));
    let mut at = (0.0, 0.0);
    let mut out = Vec::with_capacity(n);
    for c in 0..n as u64 {
        let k = 3 + 8 * c;
        if c > 0 {
            let phi = tau * rng.uniform(draw, k);
            let step = CALDERA_WANDER * size * rng.uniform(draw, k + 1);
            at = (at.0 + step * phi.cos(), at.1 + step * phi.sin());
            let (s0, s1) = CALDERA_SHRINK;
            size *= s0 + (s1 - s0) * rng.uniform(draw, k + 2);
            depth += (deepest - depth) * (0.25 + 0.5 * rng.uniform(draw, k + 3));
        }
        out.push(Caldera {
            east: at.0,
            north: at.1,
            radius: size,
            ellipticity: 1.0 + CALDERA_ELONGATION * rng.uniform(draw, k + 4),
            axis: tau * rng.uniform(draw, k + 5),
            depth,
        });
    }
    out
}

/// A shield's lava composition at a time (F101): its lavas come from one reservoir under
/// its summit, which mixes the batches that feed it, so they share the source there, the
/// drift by the time and the shield's own batch, an eruption's batch varying about it by
/// SHIELD_BATCH_SPREAD of a flood's.
pub fn shield_composition(shield: &Shield, time_gyr: f64, eruption: Option<u64>) -> f32 {
    let rng = Rng::new(shield.id ^ 0x5B1E_BA7C_0000);
    let batch =
        rng.normal(0, 0) + eruption.map_or(0.0, |j| SHIELD_BATCH_SPREAD * rng.normal(1 + j, 0));
    crate::stages::volcanism::flood_composition(time_gyr, shield.center, batch)
}

/// The time by which a shield had grown to a share of its edifice.
fn progress_time(shield: &Shield, progress: f64) -> f64 {
    shield.start_gyr - SHIELD_GROWTH_GYR * (1.0 - progress).ln()
}

/// The edifice's mean height over its area as a share of its summit's.
fn shield_mean_shape() -> f64 {
    let n = 256;
    (0..n)
        .map(|i| {
            let x = (i as f64 + 0.5) / n as f64;
            2.0 * x * shield_shape(x) / n as f64
        })
        .sum()
}

/// The growth of a shield, shares of its edifice, over which the eruptions that can
/// show on it came (F101): since the edifice last grew on average by SHIELD_FLOW_COVER
/// sheets of flows, which lie over the rest; none before it starts.
fn shield_span(shield: &Shield, present: f64) -> Option<(f64, f64)> {
    let now = shield_progress(present, shield.start_gyr);
    if now <= 0.0 {
        return None;
    }
    let laid = SHIELD_FLOW_COVER * FLOOD_SHEET_THICKNESS / (shield.height * shield_mean_shape());
    Some(((now - laid).max(0.0), now))
}

/// The times of the eruptions that can show on a shield (F101).
pub fn shield_window(shield: &Shield, present: f64) -> Option<(f64, f64)> {
    shield_span(shield, present)
        .map(|(from, _)| (progress_time(shield, from).min(present), present))
}

/// The eruptions that show on a shield (F101). Its lava comes out at circumferential
/// fissures round its summit caldera and on its flanks, and at fissures low on its rift
/// zone (Chadwick and Howard 1991), whose eruptions are the largest and run onto the
/// plains as aprons. Each eruption's lava runs down the edifice as a
/// flood's does (F94), a compound flow field: a trunk from its vent away from the summit,
/// widening as a breakout does, its branches and breakouts, until its lobes cover a field
/// as Io's do, its run-out long and half as wide, FLOOD_LOBE_COVER times over. A flow
/// narrows on the steep flanks and spreads where the ground flattens, on the summit and at
/// the foot (`trace_flow`). The eruptions follow one another until their lobes cover the
/// shield SHIELD_FLOW_COVER times over, the later lying on the earlier.
pub fn shield_flows(
    shield: &Shield,
    grid: &CubeSphere,
    body_radius: f64,
    ground: &PaddedField,
    present: f64,
) -> Vec<ShieldFlow> {
    let Some((from, now)) = shield_span(shield, present) else {
        return Vec::new();
    };
    let tau = 2.0 * std::f64::consts::PI;
    let cell = grid.delta() * body_radius;
    let rng = Rng::new(shield.id ^ 0x5F10_0000_0000);
    let center = shield.center;
    let (east, north) = local_basis(center);
    let rules = |hw: f64| PathRules {
        min_spacing: 0.5 * cell,
        slope_follow: IO_FLOW_SLOPE_FOLLOW,
        slope_scale: IO_FLOW_SLOPE_SCALE,
        steer_length: FLOOD_STEER_LENGTH,
        wander_rad: FLOOD_WANDER_RAD,
        wander_length: FLOOD_WANDER_WIDTHS * hw,
        fill_depth: FLOOD_SHEET_THICKNESS,
        width_wobble: IO_FLOW_WIDTH_WOBBLE,
    };
    let outside = |q: [f64; 3], hw: f64| {
        let (r, phi) = shield.polar(q, body_radius);
        r > SHIELD_FLOW_REACH * shield.radius_toward(phi) - hw
    };
    let target = SHIELD_FLOW_COVER * std::f64::consts::PI * shield.radius.powi(2);
    let chamber = shield
        .calderas
        .first()
        .map_or(CALDERA_RADIUS.0 * shield.radius, |c| c.radius);
    let (m_lo, m_hi) = SHIELD_RIM_VENTS;
    let mut laid: Vec<(Vec<Lobe>, bool)> = Vec::new();
    let mut area = 0.0;
    let mut i = 0u64;
    while area < target && laid.len() < SHIELD_ERUPTIONS_MAX {
        let k = 8 * i;
        let summit = rng.uniform(k, 0) < SHIELD_SUMMIT_VENT_SHARE;
        let rift = !summit && rng.uniform(k, 1) < SHIELD_RIFT_VENT_SHARE;
        let (dist, phi) = if summit {
            // A circumferential fissure round the caldera.
            (
                chamber * (m_lo + (m_hi - m_lo) * rng.uniform(k, 3)),
                tau * rng.uniform(k, 2),
            )
        } else if rift {
            // A fissure on the rift zone either side of the summit, low on the flank,
            // where the dikes that feed the largest eruptions break out.
            let phi = shield.rift
                + if rng.uniform(k, 2) < 0.5 {
                    0.0
                } else {
                    std::f64::consts::PI
                }
                + SHIELD_RIFT_SPREAD * (2.0 * rng.uniform(k, 5) - 1.0);
            let (lo, hi) = SHIELD_RIFT_VENTS;
            (
                shield.radius_toward(phi) * (lo + (hi - lo) * rng.uniform(k, 3)),
                phi,
            )
        } else {
            // A fissure anywhere on the upper flanks.
            let phi = tau * rng.uniform(k, 2);
            let (lo, hi) = SHIELD_FLANK_VENTS;
            (
                shield.radius_toward(phi) * (lo + (hi - lo) * rng.uniform(k, 3)),
                phi,
            )
        };
        let vent = offset(center, east, north, dist / body_radius, phi);
        // The lava sets off away from the summit.
        let heading = [
            vent[0] - center[0],
            vent[1] - center[1],
            vent[2] - center[2],
        ];
        // The rift zone's eruptions are the largest and run the farthest, onto the plains
        // as aprons; the fissures round the summit and on the flanks build the edifice.
        let (l_lo, l_hi) = if rift {
            SHIELD_RIFT_FLOW_LENGTH
        } else {
            SHIELD_FISSURE_FLOW_LENGTH
        };
        let run = shield.radius * l_lo * (l_hi / l_lo).powf(rng.uniform(k, 6));
        let field = FLOOD_LOBE_COVER * IO_FLOW_FIELD_AREA * run * run;
        let id = crate::rng::hash_words(&[shield.id, i, 0x5F1E_1D00]);
        let lobes = crate::stages::volcanism::flow_field(
            grid,
            body_radius,
            ground,
            &Rng::new(id),
            id,
            vent,
            heading,
            run,
            field,
            IO_FLOW_TOE,
            &rules,
            &outside,
        );
        area += lobes
            .iter()
            .map(crate::stages::volcanism::lobe_area)
            .sum::<f64>();
        laid.push((lobes, summit));
        i += 1;
    }
    let n = laid.len().max(1) as f64;
    laid.into_iter()
        .enumerate()
        .map(|(j, (lobes, summit))| {
            // Eruptions come as the supply does, evenly through the growth.
            let time =
                progress_time(shield, from + (now - from) * (j as f64 + 0.5) / n).min(present);
            ShieldFlow {
                lobes,
                time_gyr: time,
                composition: shield_composition(shield, time, Some(j as u64)),
                summit,
            }
        })
        .collect()
}

/// Height of a shield's edifice at a point, fully grown, with or without the calderas
/// cut into the summit. The loop builds the edifice without them, since they are too
/// small for its grid; the finishing pass cuts them.
pub fn shield_height(
    shield: &Shield,
    d: [f64; 3],
    body_radius: f64,
    with_calderas: bool,
    blur: f64,
) -> f64 {
    let (r, phi) = shield.polar(d, body_radius);
    let h = shield.height * shield_shape(shield.reach_share(r, phi));
    if with_calderas {
        (h - caldera_depth(shield, d, body_radius, blur)).max(0.0)
    } else {
        h
    }
}

/// Trough cross-profile: depth of the floor at a signed across distance, flat-floored
/// with walls of the given width, broadened by `blur`.
pub fn trough_profile(across: f64, width: f64, depth: f64, wall: f64, blur: f64) -> f64 {
    let w = wall.hypot(blur);
    let edge = 0.5 * width;
    let inside = 0.5 * (1.0 - ((across.abs() - edge) / w).tanh());
    depth * inside
}

/// The raised flanks of a trough whose faults' footwalls rebound: 1 at the walls,
/// falling away over a share of the trough's width outside them, nothing inside.
pub fn flank_profile(across: f64, width: f64, wall: f64, blur: f64) -> f64 {
    let w = wall.hypot(blur);
    let edge = 0.5 * width;
    let outside = 0.5 * (1.0 + ((across.abs() - edge) / w).tanh());
    let beyond = (across.abs() - edge).max(0.0);
    outside * (-beyond / (TROUGH_FLANK_DECAY * width).max(1.0)).exp()
}

/// The dome, canyon and shields through the coarse loop.
pub struct PlumeProcess {
    grid: CubeSphere,
    body_radius: f64,
    pub rises: Vec<Rise>,
    pub shields: Vec<Shield>,
    pub canyons: Vec<Trough>,
    pub graben: Vec<Trough>,
    /// Cumulative dome uplift per coarse cell, m.
    pub load: Vec<f32>,
    rise_height: f64,
    /// The body's age, Gyr.
    present: f64,
    /// Shield height so far per coarse cell, m.
    shield_built: Vec<f32>,
    /// Canyon depth carved so far per coarse cell, m.
    carved: Vec<f32>,
    /// The share of that the shields' lava has filled since, m.
    filled: Vec<f32>,
    stamper: Stamper,
}

impl PlumeProcess {
    /// None when the body raises no plume.
    pub fn new(graph: &ParamGraph, grid: CubeSphere) -> Option<Self> {
        let count = graph.f64("tectonics.plume_rise_count").round() as usize;
        if count == 0
            || !matches!(
                graph.str("interior.tectonic_regime"),
                "stagnant_lid" | "episodic_lid"
            )
        {
            return None;
        }
        let body_radius = graph.f64("body.radius");
        let episodic = graph.str("interior.tectonic_regime") == "episodic_lid";
        // On an episodic lid the rifts and graben come after the event, which would
        // otherwise fill them; on a stagnant lid they follow the dome.
        let age = graph.f64("body.age");
        let event_time = if episodic {
            (age - graph.f64("tectonics.resurfacing_event_age")).max(0.0)
        } else {
            0.0
        };
        let rng = Rng::new(graph.f64("seed.tectonics") as u64 ^ 0x9D1F_E5A7);
        let rise_radius = 0.5 * graph.f64("tectonics.plume_rise_diameter");
        let rise_height = graph.f64("tectonics.plume_rise_height");
        let mut rises = Vec::new();
        for i in 0..count as u64 {
            let mut center = rng.unit_direction(i, 1);
            // Rises keep apart.
            for attempt in 0..16u64 {
                if rises.iter().all(|r: &Rise| {
                    angle_between(r.center, center) > 2.5 * rise_radius / body_radius
                }) {
                    break;
                }
                center = rng.unit_direction(i, 10 + 2 * attempt);
            }
            rises.push(Rise {
                center,
                radius: rise_radius * (0.8 + 0.4 * rng.uniform(i, 4)),
                height: rise_height * (0.8 + 0.4 * rng.uniform(i, 5)),
                start_gyr: RISE_START_GYR,
                end_gyr: RISE_END_GYR,
            });
        }
        // Shields: a few large ones on the first rise, smaller ones on each.
        let shield_count = graph.f64("volcanism.shield_count").round() as usize;
        let relief_limit = graph.f64("interior.max_relief");
        let mut shields = Vec::new();
        let mut sid = 0u64;
        for (ri, rise) in rises.iter().enumerate() {
            let n_here = if ri == 0 {
                shield_count
            } else {
                (shield_count / 2).max(1)
            };
            let (east, north) = local_basis(rise.center);
            for k in 0..n_here as u64 {
                let draw = 100 + sid * 16;
                // The largest first, then down a power law.
                let rank = k as f64 + 1.0;
                let diameter = SHIELD_MAX_DIAMETER_FRACTION * 2.0 * rise.radius / rank.powf(0.5);
                let radius = 0.5 * diameter;
                let height = (radius * SHIELD_SLOPE_DEG.to_radians().tan()).min(relief_limit);
                let mut center = rise.center;
                for attempt in 0..24u64 {
                    let dist = rise.radius
                        * (SHIELD_PLACEMENT.0
                            + (SHIELD_PLACEMENT.1 - SHIELD_PLACEMENT.0)
                                * rng.uniform(draw, 1 + 3 * attempt));
                    let phi = 2.0 * std::f64::consts::PI * rng.uniform(draw, 2 + 3 * attempt);
                    center = offset(rise.center, east, north, dist / body_radius, phi);
                    let clear = shields.iter().all(|s: &Shield| {
                        angle_between(s.center, center) * body_radius > 1.1 * (s.radius + radius)
                    });
                    if clear {
                        break;
                    }
                }
                let calderas = summit_collapses(&rng, draw ^ 0xCA1D_E7A5_0000, radius, height);
                // Its rift zone runs along the rise's radial graben, where the dikes rise
                // that feed its flank vents.
                let (se, sn) = local_basis(center);
                let rift = azimuth_from(center, se, sn, rise.center);
                shields.push(Shield {
                    center,
                    radius,
                    height,
                    start_gyr: SHIELD_START_GYR.0
                        + (SHIELD_START_GYR.1 - SHIELD_START_GYR.0) * rng.uniform(draw, 14),
                    calderas,
                    rift,
                    elongation: 1.0 + (SHIELD_ELONGATION_MAX - 1.0) * rng.uniform(draw, 15),
                    id: sid,
                });
                sid += 1;
            }
        }
        // Rifting: the uplift stretches the lithosphere, and the extension is taken up
        // by a population of normal-fault troughs, from chasmata kilometres deep to
        // graben at the grid's limit, drawn from a power law until the budget is spent.
        // Faults cluster beside faults, so the big ones gather en-echelon neighbours and
        // the small ones form swarms.
        let cell = grid.delta() * body_radius;
        let extension = graph.f64("tectonics.rift_extension");
        let mut canyons = Vec::new();
        let mut graben = Vec::new();
        for (ri, rise) in rises.iter().enumerate() {
            let (east, north) = local_basis(rise.center);
            let budget = extension
                * (rise.height / rise_height.max(1.0))
                * 2.0
                * std::f64::consts::PI
                * RIFT_TRANSECT_FRACTION
                * rise.radius;
            let l_min = RIFT_MIN_LENGTH;
            let l_max = (RIFT_MAX_LENGTH_FRACTION * 2.0 * rise.radius).max(2.0 * l_min);
            let (t_lo, t_hi) = if episodic {
                (event_time, age)
            } else {
                (rise.start_gyr, RIFT_END_GYR.min(age))
            };
            let depth_cap = (CANYON_DEPTH_FRACTION * rise.height).max(RIFT_MIN_DEPTH);
            let spacing = 0.5 * cell;
            let mut faults: Vec<Trough> = Vec::new();
            let mut spent = 0.0;
            let mut k = 0u64;
            while spent < budget && k < RIFT_MAX_FAULTS {
                let draw = 5_000 + ri as u64 * 1_000_000 + k * 16;
                k += 1;
                // Length from the truncated power law by inverse transform.
                let u = rng.uniform(draw, 1).clamp(1e-9, 1.0 - 1e-9);
                let ratio = (l_min / l_max).powf(RIFT_LENGTH_SLOPE);
                let mut length = l_min * (1.0 - u * (1.0 - ratio)).powf(-1.0 / RIFT_LENGTH_SLOPE);
                let scatter = (RIFT_DEPTH_SCATTER * rng.normal(draw, 2)).exp();
                let depth =
                    (RIFT_DEPTH_PER_LENGTH * length * scatter).clamp(RIFT_MIN_DEPTH, depth_cap);
                let width = (RIFT_WIDTH_PER_DEPTH * depth).max(RIFT_MIN_WIDTH);
                // Beside an existing fault, chosen by length, or on its own.
                let host = if !faults.is_empty() && rng.uniform(draw, 3) < RIFT_CLUSTER_PROBABILITY
                {
                    let total: f64 = faults.iter().map(|f| f.length).sum();
                    let mut pick = rng.uniform(draw, 4) * total;
                    let mut chosen = faults.len() - 1;
                    for (i, f) in faults.iter().enumerate() {
                        if pick < f.length {
                            chosen = i;
                            break;
                        }
                        pick -= f.length;
                    }
                    Some(chosen)
                } else {
                    None
                };
                let points = match host {
                    Some(h) => {
                        let hf = &faults[h];
                        let n_points = hf.points.len();
                        let span = ((length / hf.spacing).round() as usize).clamp(2, n_points);
                        let i0 = ((n_points - span) as f64 * rng.uniform(draw, 5)) as usize;
                        let side = if rng.uniform(draw, 6) < 0.5 {
                            -1.0
                        } else {
                            1.0
                        };
                        let (g_lo, g_hi) = RIFT_CLUSTER_OFFSET;
                        let gap = g_lo + (g_hi - g_lo) * rng.uniform(draw, 7).powi(2);
                        length = span as f64 * hf.spacing;
                        offset_trace(
                            &hf.points[i0..i0 + span],
                            side * gap * hf.width.max(width),
                            body_radius,
                        )
                    }
                    None => {
                        let phi = 2.0 * std::f64::consts::PI * rng.uniform(draw, 5);
                        // The extension is greatest on the rise: the longest faults
                        // start on its flank, the short ones anywhere out past its foot.
                        let (r_lo, r_hi) = RIFT_START_FRACTION;
                        let r_hi = r_hi - (r_hi - r_lo) * 0.75 * (length / l_max).min(1.0);
                        let r0 = rise.radius * (r_lo + (r_hi - r_lo) * rng.uniform(draw, 6));
                        let start = offset(rise.center, east, north, r0 / body_radius, phi);
                        radial_trace(start, rise.center, length, spacing, body_radius, &rng, draw)
                    }
                };
                let time = t_lo + (t_hi - t_lo) * rng.uniform(draw, 8);
                let carried = width >= RIFT_LOOP_WIDTH_CELLS * cell;
                let (start_gyr, end_gyr) = if carried {
                    (
                        (time - 0.5 * RIFT_CARVING_GYR).max(t_lo),
                        (time + 0.5 * RIFT_CARVING_GYR).min(age),
                    )
                } else {
                    (time, time)
                };
                spent += RIFT_HEAVE_PER_DEPTH * depth * length;
                faults.push(Trough {
                    points,
                    spacing,
                    length,
                    width,
                    depth,
                    flank: 0.0,
                    start_gyr,
                    end_gyr,
                    id: 1_000 + ri as u64 * 1_000_000 + k,
                });
            }
            for f in faults {
                if f.width >= RIFT_LOOP_WIDTH_CELLS * cell {
                    canyons.push(f);
                } else {
                    graben.push(f);
                }
            }
        }
        Some(Self {
            grid,
            body_radius,
            rises,
            shields,
            canyons,
            graben,
            load: vec![0.0; grid.cell_count()],
            rise_height: rise_height.max(1.0),
            present: graph.f64("body.age"),
            shield_built: vec![0.0; grid.cell_count()],
            carved: vec![0.0; grid.cell_count()],
            filled: vec![0.0; grid.cell_count()],
            stamper: Stamper::new(&grid),
        })
    }

    /// Multiplier on the vent weights: eruptions favour the rise.
    pub fn vent_boost(&self) -> Vec<f64> {
        self.load
            .iter()
            .map(|&l| 1.0 + PLUME_VENT_BOOST * (l as f64 / self.rise_height).clamp(0.0, 1.0))
            .collect()
    }

    /// Dome uplift, shield growth and canyon carving between two times.
    pub fn step(&mut self, layers: &mut LayerSet, t0: f64, t1: f64) {
        let grid = self.grid;
        let body_radius = self.body_radius;
        // Dome: lift everything, which carries craters up with it.
        for rise in &self.rises {
            let dp =
                smooth(t1, rise.start_gyr, rise.end_gyr) - smooth(t0, rise.start_gyr, rise.end_gyr);
            if dp <= 0.0 {
                continue;
            }
            self.stamper
                .gather(&grid, body_radius, rise.center, rise.radius);
            let bedrock = layers.bedrock.as_mut_slice();
            for &(idx, r, _) in &self.stamper.footprint {
                let dh = rise.height * dp * dome(r / rise.radius);
                bedrock[idx] += dh as f32;
                self.load[idx] += dh as f32;
            }
        }
        // Shields: flows pile up, renewing the surface where they add more than a veneer.
        let t_mid = 0.5 * (t0 + t1);
        for shield in &self.shields {
            let dp = shield_progress(t1, shield.start_gyr) - shield_progress(t0, shield.start_gyr);
            if dp <= 0.0 {
                continue;
            }
            self.stamper
                .gather(&grid, body_radius, shield.center, shield.extent());
            let bedrock = layers.bedrock.as_mut_slice();
            let ages = layers.surface_age.as_mut_slice();
            let units = layers.unit_id.as_mut_slice();
            let deep = layers.deep.as_mut_slice();
            let melt = layers.melt.as_mut_slice();
            let cover = layers.cover.as_mut_slice();
            // The eruptions that show from the window's start are drawn in the finishing
            // pass, so the surface they lie on is dated no later than that (F101); the
            // finishing pass lays its lava's composition at its own resolution.
            let renewed = shield_window(shield, self.present).map_or(t_mid, |w| t_mid.min(w.0));
            for &(idx, _, _) in &self.stamper.footprint {
                let d = grid.center(grid.cell(idx));
                let full = shield_height(shield, d, body_radius, false, 0.0);
                let dh = full * dp;
                if dh <= 0.0 {
                    continue;
                }
                // The lava ponds in a canyon it comes to, filling it as fast as the
                // shield thickens round it (F103).
                let fill = (self.carved[idx] - self.filled[idx]).clamp(0.0, dh as f32);
                self.filled[idx] += fill;
                bedrock[idx] += dh as f32 + fill;
                self.shield_built[idx] += dh as f32;
                if dh > SHIELD_RENEWAL_THICKNESS {
                    ages[idx] = renewed as f32;
                    units[idx] = UNIT_LAVA;
                    deep[idx] = 0.0;
                    melt[idx] = 0.0;
                    cover[idx] = 1.0;
                }
            }
        }
        if std::env::var_os("EXOGEN_DEBUG").is_some()
            && let (Some(rise), Some(shield)) = (self.rises.first(), self.shields.first())
        {
            {
                let rc = grid.index(grid.locate(rise.center).0);
                let sc = grid.index(grid.locate(shield.center).0);
                let b = layers.bedrock.as_slice();
                eprintln!(
                    "plume t={t1:.2} rise load {:.0} bedrock {:.0} | shield built {:.0} bedrock {:.0}",
                    self.load[rc], b[rc], self.shield_built[sc], b[sc]
                );
            }
        }
        // Canyon: the trough deepens and widens through its growth.
        for canyon in &self.canyons {
            let p1 = smooth(t1, canyon.start_gyr, canyon.end_gyr);
            let p0 = smooth(t0, canyon.start_gyr, canyon.end_gyr);
            if p1 <= p0 {
                continue;
            }
            let cell = grid.delta() * body_radius;
            let target =
                trough_depths(&grid, body_radius, canyon, p1, cell, 0.0, &mut self.stamper);
            let bedrock = layers.bedrock.as_mut_slice();
            let ages = layers.surface_age.as_mut_slice();
            for (idx, depth) in target {
                let more = depth - self.carved[idx] as f64;
                if more <= 0.0 {
                    continue;
                }
                bedrock[idx] -= more as f32;
                self.carved[idx] += more as f32;
                if more > CANYON_RENEWAL_THICKNESS {
                    ages[idx] = t_mid as f32;
                }
            }
        }
    }
}

/// A trace shifted sideways by a distance: each point moved along the local normal.
fn offset_trace(points: &[[f64; 3]], lateral: f64, body_radius: f64) -> Vec<[f64; 3]> {
    let n = points.len();
    (0..n)
        .map(|k| {
            let p = points[k];
            let prev = points[k.saturating_sub(1)];
            let next = points[(k + 1).min(n - 1)];
            let chord = [next[0] - prev[0], next[1] - prev[1], next[2] - prev[2]];
            let radial = dot(chord, p);
            let t = normalize([
                chord[0] - radial * p[0],
                chord[1] - radial * p[1],
                chord[2] - radial * p[2],
            ]);
            let normal = normalize(cross(p, t));
            let s = lateral / body_radius;
            normalize([
                p[0] + s * normal[0],
                p[1] + s * normal[1],
                p[2] + s * normal[2],
            ])
        })
        .collect()
}

/// A trace running outward from `start`, away from `center`, with a slow wobble.
pub(crate) fn radial_trace(
    start: [f64; 3],
    center: [f64; 3],
    length: f64,
    spacing: f64,
    body_radius: f64,
    rng: &Rng,
    draw: u64,
) -> Vec<[f64; 3]> {
    let n = (length / spacing).ceil().max(2.0) as usize;
    let mut points = Vec::with_capacity(n + 1);
    let mut p = start;
    let noise = GradientNoise::new(draw ^ 0x7A11_CE00);
    let turn_scale = TROUGH_WOBBLE_DEG.to_radians() * (0.5 + rng.uniform(draw, 8));
    for i in 0..=n {
        points.push(p);
        // Outward from the centre in the tangent plane at p.
        let radial = dot(center, p);
        let mut out = normalize([
            p[0] * radial - center[0],
            p[1] * radial - center[1],
            p[2] * radial - center[2],
        ]);
        if dot(out, out) < 0.5 {
            let (east, _) = local_basis(p);
            out = east;
        }
        let s = i as f64 * spacing;
        let turn = turn_scale
            * (noise.fbm([s / TROUGH_WOBBLE_LENGTH, 0.5, 0.25], 1.0, 2, 2.0, 0.5) / 0.6)
                .clamp(-1.0, 1.0);
        let side = normalize(cross(p, out));
        let heading = normalize([
            out[0] * turn.cos() + side[0] * turn.sin(),
            out[1] * turn.cos() + side[1] * turn.sin(),
            out[2] * turn.cos() + side[2] * turn.sin(),
        ]);
        let step = spacing / body_radius;
        p = normalize([
            p[0] + step * heading[0],
            p[1] + step * heading[1],
            p[2] + step * heading[2],
        ]);
    }
    points
}

/// A trough's trace resolved into frames: each point with its tangent and normal, the
/// taper and the width pulse there. Shared by the stamp over a grid and the depth at a
/// point.
pub struct TroughFrames {
    frames: Vec<Frame>,
    width: f64,
    flank: f64,
    noise: GradientNoise,
}

struct Frame {
    p: [f64; 3],
    t: [f64; 3],
    normal: [f64; 3],
    /// Depth and width factors towards the ends of the trace.
    taper: f64,
    narrow: f64,
    pulse: f64,
}

/// A trough's dimensions at a growth fraction, m.
struct TroughSize {
    width: f64,
    depth: f64,
    wall: f64,
    search: f64,
    half_band: f64,
}

impl TroughFrames {
    pub fn new(trough: &Trough) -> Self {
        let n = trough.points.len();
        let noise = GradientNoise::new(trough.id ^ 0xCA17_0000);
        let mut frames = Vec::with_capacity(n);
        for k in 0..n {
            let p = trough.points[k];
            let prev = trough.points[k.saturating_sub(1)];
            let next = trough.points[(k + 1).min(n - 1)];
            let chord = [next[0] - prev[0], next[1] - prev[1], next[2] - prev[2]];
            let radial = dot(chord, p);
            let t = normalize([
                chord[0] - radial * p[0],
                chord[1] - radial * p[1],
                chord[2] - radial * p[2],
            ]);
            if dot(t, t) < 0.5 {
                continue;
            }
            let normal = normalize(cross(p, t));
            let along_trace = k as f64 * trough.spacing;
            let taper = (std::f64::consts::PI * along_trace / trough.length)
                .sin()
                .max(0.0);
            if taper <= 0.0 {
                continue;
            }
            let narrow = taper.sqrt();
            let taper = taper.powf(TROUGH_TAPER_POWER);
            // Chasmata: the width swells and narrows along the trace, and each wall is
            // scalloped.
            let pulse = (1.0
                + TROUGH_WIDTH_PULSE
                    * (noise.fbm(
                        [along_trace / TROUGH_PULSE_LENGTH, 0.37, 0.61],
                        1.0,
                        2,
                        2.0,
                        0.5,
                    ) / 0.6))
                .clamp(0.5, 1.6);
            frames.push(Frame {
                p,
                t,
                normal,
                taper,
                narrow,
                pulse,
            });
        }
        Self {
            frames,
            width: trough.width,
            flank: trough.flank,
            noise,
        }
    }

    /// The trace points that carry a frame, in frame order.
    pub fn points(&self) -> impl Iterator<Item = [f64; 3]> + '_ {
        self.frames.iter().map(|f| f.p)
    }

    pub fn len(&self) -> usize {
        self.frames.len()
    }

    pub fn is_empty(&self) -> bool {
        self.frames.is_empty()
    }

    fn size(trough: &Trough, progress: f64, wall_floor: f64, cell: f64) -> TroughSize {
        let width = trough.width * (0.3 + 0.7 * progress);
        let depth = trough.depth * progress;
        // Walls no narrower than the grid that carries them.
        let wall = (TROUGH_WALL_FRACTION * width).max(wall_floor);
        TroughSize {
            width,
            depth,
            wall,
            search: 0.5 * width
                + 3.0 * wall
                + 2.0 * cell
                + if trough.flank != 0.0 { width } else { 0.0 },
            half_band: 0.5 * trough.spacing + cell,
        }
    }

    /// How far from a trace point the trough reaches at a growth fraction, m.
    pub fn search(&self, trough: &Trough, progress: f64, wall_floor: f64, cell: f64) -> f64 {
        Self::size(trough, progress, wall_floor, cell).search
    }

    /// The depth a frame gives a point, with the point's distance along strike, or
    /// nothing when the point is outside the frame's band.
    fn depth_from(
        &self,
        d: [f64; 3],
        body_radius: f64,
        size: &TroughSize,
        blur: f64,
        f: &Frame,
    ) -> Option<(f64, f64)> {
        let along = dot(d, f.t) * body_radius;
        if along.abs() > size.half_band {
            return None;
        }
        let across = dot(d, f.normal) * body_radius;
        let scallop_frequency = body_radius / TROUGH_SCALLOP_LENGTH.max(self.width);
        let scallop = 1.0
            + TROUGH_SCALLOP
                * self.noise.fbm(
                    [d[0] + across.signum() * 3.0, d[1], d[2]],
                    scallop_frequency,
                    2,
                    2.0,
                    0.5,
                );
        let width = size.width * f.pulse * scallop * f.narrow;
        let mut u = trough_profile(across, width, size.depth, size.wall, blur);
        if self.flank != 0.0 {
            u -= self.flank * size.depth * flank_profile(across, width, size.wall, blur);
        }
        let u = f.taper * u;
        Some((along.abs(), u))
    }
}

/// Depth of a trough at each cell it touches, at a growth fraction, with scalloped
/// walls and a tapered floor: nearest trace point along strike.
pub fn trough_depths(
    grid: &CubeSphere,
    body_radius: f64,
    trough: &Trough,
    progress: f64,
    wall_floor: f64,
    blur: f64,
    stamper: &mut Stamper,
) -> Vec<(usize, f64)> {
    let cell = grid.delta() * body_radius;
    let frames = TroughFrames::new(trough);
    let size = TroughFrames::size(trough, progress, wall_floor, cell);
    let mut best: HashMap<usize, (f64, f64)> = HashMap::new();
    for f in &frames.frames {
        stamper.gather(grid, body_radius, f.p, size.search);
        for &(idx, _, _) in &stamper.footprint {
            let d = grid.center(grid.cell(idx));
            if let Some((along, u)) = frames.depth_from(d, body_radius, &size, blur, f) {
                let entry = best.entry(idx).or_insert((f64::INFINITY, 0.0));
                if along < entry.0 {
                    *entry = (along, u);
                }
            }
        }
    }
    let mut out: Vec<(usize, f64)> = best.into_iter().map(|(i, (_, u))| (i, u)).collect();
    out.sort_by_key(|e| e.0);
    out
}

/// Depth of a trough at a point, m, as a grid of cell `cell` with walls no narrower
/// than `wall_floor` carries it: the same profile the stamp lays down, at any point.
#[allow(clippy::too_many_arguments)]
pub fn trough_depth_at(
    d: [f64; 3],
    body_radius: f64,
    trough: &Trough,
    frames: &TroughFrames,
    frame_range: (usize, usize),
    progress: f64,
    wall_floor: f64,
    cell: f64,
    blur: f64,
) -> f64 {
    let size = TroughFrames::size(trough, progress, wall_floor, cell);
    let cos_reach = ((size.search + 0.5 * cell) / body_radius).cos();
    let mut best = (f64::INFINITY, 0.0);
    let (lo, hi) = frame_range;
    for f in &frames.frames[lo.min(frames.frames.len())..hi.min(frames.frames.len())] {
        if dot(f.p, d) < cos_reach {
            continue;
        }
        if let Some((along, u)) = frames.depth_from(d, body_radius, &size, blur, f)
            && along < best.0
        {
            best = (along, u);
        }
    }
    best.1
}

/// The shields' lava laid on a point since a time, m (F103): each shield's growth there
/// since then, which fills a trough formed before it, as it ponds in the trough rather
/// than draping it.
pub fn shield_burial(
    shields: &[Shield],
    d: [f64; 3],
    body_radius: f64,
    since: f64,
    present: f64,
) -> f64 {
    shields
        .iter()
        .map(|s| {
            let grown = shield_progress(present, s.start_gyr) - shield_progress(since, s.start_gyr);
            if grown <= 0.0 || angle_between(s.center, d) * body_radius > s.extent() {
                0.0
            } else {
                shield_height(s, d, body_radius, false, 0.0) * grown
            }
        })
        .sum()
}

/// Finishing pass of a trough: the fine profile minus the version the loop carried,
/// with the loop's wall width and the coarse smear, or the whole trough when the loop
/// did not carry it. Inside the floor the fine relief of older features is removed, as a
/// crater's cavity removes it. A shield that grew over the trough since has filled it by
/// the lava it laid there (F103), the loop's coarse trough as well, and the trough shows
/// only where it is deeper than that, its raised flanks only where higher.
#[allow(clippy::too_many_arguments)]
pub fn trough_compute(
    layers: &LayerSet,
    grid: &CubeSphere,
    body_radius: f64,
    base: &[f32],
    trough: &Trough,
    time_gyr: f32,
    blur_fine: f64,
    coarse: Option<(f64, f64)>,
    shields: &[Shield],
    present: f64,
    stamper: &mut Stamper,
) -> StampOutput {
    let bedrock = layers.bedrock.as_slice();
    let ages = layers.surface_age.as_slice();
    let cell = grid.delta() * body_radius;
    let fine = trough_depths(grid, body_radius, trough, 1.0, cell, blur_fine, stamper);
    let coarse: HashMap<usize, f64> = match coarse {
        Some((coarse_cell, blur)) => {
            trough_depths(grid, body_radius, trough, 1.0, coarse_cell, blur, stamper)
                .into_iter()
                .collect()
        }
        None => HashMap::new(),
    };
    let cells = fine
        .into_iter()
        .map(|(idx, depth)| {
            let buried = shield_burial(
                shields,
                grid.center(grid.cell(idx)),
                body_radius,
                time_gyr as f64,
                present,
            );
            let smeared = (coarse.get(&idx).copied().unwrap_or(0.0) - buried).max(0.0);
            let depth = depth.signum() * (depth.abs() - buried).max(0.0);
            // How deep inside the trough this cell is, 0 on the rim to 1 on the floor.
            let inside = (depth / trough.depth.max(1.0)).clamp(0.0, 1.0);
            // What the finish has added since the base is kept outside the floor and
            // removed inside it, as a crater's cavity removes what it cuts. (Taking the
            // current ground as the reference as well as keeping the difference counted
            // that difference twice, and overlapping graben doubled each other.)
            let kept = (bedrock[idx] - base[idx]) as f64 * (1.0 - inside);
            let mut write = CellWrite {
                idx: idx as u32,
                bedrock: (base[idx] as f64 + kept - (depth - smeared)) as f32,
                surface: false,
                age: 0.0,
                unit: 0,
                deep: 0.0,
                melt: 0.0,
                cover: 0.0,
            };
            if inside > 0.5 && ages[idx] <= time_gyr + 1e-6 {
                write.surface = true;
                write.age = time_gyr;
                write.unit = layers.unit_id.as_slice()[idx];
                write.deep = layers.deep.as_slice()[idx];
                write.melt = layers.melt.as_slice()[idx];
                write.cover = layers.cover.as_slice()[idx];
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
        compositions: None,
    }
}

/// Finishing pass of a shield (F101, F102): the calderas cut into the summit, and the
/// eruptions that show laid over the edifice, whose surface the loop dated to their
/// start, each flow field a sheet with fronts where its lobes end; the ground takes the
/// topmost eruption's lava, its age and composition, and the caldera floors are the lava
/// the last summit eruption ponded in them. Elsewhere on the edifice the loop's lava keeps
/// its age and takes the shield's composition when the eruptions began, out to the
/// edifice's own outline, over the tessera a resurfacing may have made of it. Only the
/// eruptions after `after` are laid: those before a resurfacing lie under it.
#[allow(clippy::too_many_arguments)]
pub fn shield_compute(
    layers: &LayerSet,
    grid: &CubeSphere,
    body_radius: f64,
    shield: &Shield,
    ground: Option<&PaddedField>,
    present: f64,
    after: f64,
    stamper: &mut Stamper,
) -> StampOutput {
    let flows: Vec<ShieldFlow> = ground.map_or_else(Vec::new, |g| {
        shield_flows(shield, grid, body_radius, g, present)
            .into_iter()
            .filter(|f| f.time_gyr > after)
            .collect()
    });
    let covers: Vec<HashMap<usize, f32>> = flows
        .iter()
        .enumerate()
        .map(|(j, f)| {
            crate::stages::volcanism::lobe_cover(
                &f.lobes,
                grid,
                body_radius,
                shield.id ^ ((j as u64) << 24),
                stamper,
            )
        })
        .collect();
    let ponded = flows
        .iter()
        .rev()
        .find(|f| f.summit)
        .or(flows.last())
        .map(|f| (f.time_gyr as f32, f.composition));
    let bedrock = layers.bedrock.as_slice();
    let ages = layers.surface_age.as_slice();
    let units = layers.unit_id.as_slice();
    // The edifice's lava under the eruptions that show, as it was when they began.
    let older = shield_window(shield, present).map(|w| shield_composition(shield, w.0, None));
    stamper.gather(
        grid,
        body_radius,
        shield.center,
        shield.cap_radius(body_radius) * body_radius,
    );
    let mut cells = Vec::with_capacity(stamper.footprint.len());
    let mut compositions = Vec::with_capacity(stamper.footprint.len());
    for &(idx, _, _) in &stamper.footprint {
        let d = grid.center(grid.cell(idx));
        let fine = shield_height(shield, d, body_radius, true, 0.0);
        let coarse = shield_height(shield, d, body_radius, false, 0.0);
        let mut stack = 0.0;
        let mut top: Option<usize> = None;
        for (j, c) in covers.iter().enumerate() {
            if let Some(&v) = c.get(&idx) {
                stack += FLOOD_SHEET_THICKNESS * v as f64;
                if v >= 0.5 {
                    top = Some(j);
                }
            }
        }
        let mut write = CellWrite {
            idx: idx as u32,
            bedrock: (bedrock[idx] as f64 + (fine - coarse) + stack) as f32,
            surface: false,
            age: 0.0,
            unit: 0,
            deep: 0.0,
            melt: 0.0,
            cover: 0.0,
        };
        let mut composition = 0.0;
        let lava = match top {
            Some(j) => Some((flows[j].time_gyr as f32, flows[j].composition)),
            None if coarse - fine > 1.0 => ponded,
            None if coarse > 0.0 && matches!(units[idx], UNIT_LAVA | UNIT_TESSERA) => {
                older.map(|c| (ages[idx], c))
            }
            None => None,
        };
        if let Some((age, comp)) = lava {
            write.surface = true;
            write.age = age;
            write.unit = UNIT_LAVA;
            write.cover = 1.0;
            composition = comp;
        }
        cells.push(write);
        compositions.push(composition);
    }
    StampOutput {
        cells,
        composition: None,
        compositions: Some(compositions),
        stain: None,
        steep: None,
        ejecta_age: None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shapes_fall_from_the_centre_to_the_edge() {
        assert!((dome(0.0) - 1.0).abs() < 1e-12 && dome(1.0) == 0.0);
        assert!(dome(0.5) > 0.5 && dome(0.5) < 0.6);
        assert!((shield_shape(0.0) - 1.0).abs() < 1e-12 && shield_shape(1.0) == 0.0);
        let s = shield_shape(0.5);
        assert!(s > 0.6 && s < 0.95, "{s}");
        // A trough is full depth on the floor and nothing beyond the walls.
        let d = trough_profile(0.0, 200e3, 7e3, 20e3, 0.0);
        assert!((d - 7e3).abs() < 1.0);
        assert!(trough_profile(200e3, 200e3, 7e3, 20e3, 0.0) < 10.0);
        let mid = trough_profile(100e3, 200e3, 7e3, 20e3, 0.0);
        assert!((mid - 3.5e3).abs() < 1.0);
    }

    #[test]
    fn a_shield_has_a_summit_caldera() {
        let shield = Shield {
            center: [1.0, 0.0, 0.0],
            radius: 300e3,
            height: 20e3,
            start_gyr: 1.0,
            calderas: vec![Caldera {
                east: 0.0,
                north: 0.0,
                radius: 30e3,
                ellipticity: 1.0,
                axis: 0.0,
                depth: 2.5e3,
            }],
            rift: 0.0,
            elongation: 1.0,
            id: 1,
        };
        let r = 3.39e6;
        let summit = shield_height(&shield, [1.0, 0.0, 0.0], r, true, 0.0);
        assert!((summit - 17.5e3).abs() < 1.0, "{summit}");
        let (east, north) = local_basis(shield.center);
        let rim = offset(shield.center, east, north, 45e3 / r, 0.3);
        let h_rim = shield_height(&shield, rim, r, true, 0.0);
        assert!(h_rim > 19.0e3 && h_rim < 20.0e3, "{h_rim}");
        let flank = offset(shield.center, east, north, 150e3 / r, 1.0);
        let h_flank = shield_height(&shield, flank, r, true, 0.0);
        assert!(h_flank < h_rim && h_flank > 10e3, "{h_flank}");
        assert!(
            shield_height(
                &shield,
                offset(shield.center, east, north, shield.extent() / r + 1e-4, 0.0),
                r,
                true,
                0.0
            ) == 0.0
        );
    }

    #[test]
    fn a_summit_sinks_in_collapses_each_smaller_and_deeper() {
        // F102: over many shields, a summit sinks in one to CALDERA_COUNT_MAX collapses,
        // each no larger and no shallower than the one before, wandering from the first.
        let rng = Rng::new(9);
        let mut most = 0;
        for s in 0..200u64 {
            let calderas = summit_collapses(&rng, s * 64, 300e3, 20e3);
            most = most.max(calderas.len());
            assert!((1..=CALDERA_COUNT_MAX).contains(&calderas.len()));
            for pair in calderas.windows(2) {
                assert!(pair[1].radius <= pair[0].radius + 1e-6);
                assert!(pair[1].depth >= pair[0].depth - 1e-6);
            }
            let first = calderas[0].radius;
            for c in &calderas {
                assert!(c.east.hypot(c.north) < 3.0 * first);
                assert!(c.depth <= CALDERA_MAX_DEPTH + 1e-6);
            }
        }
        assert!(most >= 4, "the most collapses on a summit was {most}");
    }

    #[test]
    fn a_shield_s_flows_run_down_its_flanks() {
        // F101: the eruptions that show on a shield start round its summit or on its
        // flanks, their trunks run down it, and together their lobes cover it more than
        // once.
        let grid = CubeSphere::new(256);
        let body_radius = 3.39e6;
        let center = normalize([0.3, 1.0, 0.2]);
        let shield = Shield {
            center,
            radius: 300e3,
            height: 20e3,
            start_gyr: 1.0,
            calderas: summit_collapses(&Rng::new(3), 0, 300e3, 20e3),
            rift: 0.7,
            elongation: 1.3,
            id: 5,
        };
        let height: Vec<f32> = (0..grid.cell_count())
            .map(|k| {
                shield_height(&shield, grid.center(grid.cell(k)), body_radius, false, 0.0) as f32
            })
            .collect();
        let ground = PaddedField::new(&grid, &height);
        let flows = shield_flows(&shield, &grid, body_radius, &ground, 4.5);
        assert!(flows.len() > 5, "{} flows", flows.len());
        let mut downhill = 0;
        for f in &flows {
            let p = &f.lobes[0].points;
            let h = |q: [f64; 3]| shield_height(&shield, q, body_radius, false, 0.0);
            if h(p[p.len() - 1]) < h(p[0]) {
                downhill += 1;
            }
            for q in f.lobes.iter().flat_map(|l| &l.points) {
                let (r, phi) = shield.polar(*q, body_radius);
                assert!(r <= SHIELD_FLOW_REACH * shield.radius_toward(phi) + 1e3);
            }
        }
        assert!(
            downhill * 10 >= flows.len() * 9,
            "{downhill} of {} run downhill",
            flows.len()
        );
        let area: f64 = flows
            .iter()
            .flat_map(|f| &f.lobes)
            .map(crate::stages::volcanism::lobe_area)
            .sum();
        assert!(area > std::f64::consts::PI * shield.radius.powi(2));
        assert!(flows.windows(2).all(|w| w[1].time_gyr > w[0].time_gyr));
    }

    #[test]
    fn a_shield_buries_the_troughs_it_grows_over() {
        // F103: a shield fills a trough formed before it by the lava it laid there since:
        // near its summit kilometres, at its foot a few hundred metres, and over a trough
        // formed near the present hardly any; nothing beyond its foot.
        let shield = Shield {
            center: [1.0, 0.0, 0.0],
            radius: 300e3,
            height: 20e3,
            start_gyr: 1.0,
            calderas: Vec::new(),
            rift: 0.0,
            elongation: 1.0,
            id: 3,
        };
        let shields = std::slice::from_ref(&shield);
        let r = 3.39e6;
        let (east, north) = local_basis(shield.center);
        let at = |dist: f64| offset(shield.center, east, north, dist / r, 0.3);
        let summit = shield_burial(shields, at(10e3), r, 0.9, 4.5);
        let foot = shield_burial(shields, at(0.99 * shield.radius_toward(0.3)), r, 0.9, 4.5);
        let late = shield_burial(shields, at(10e3), r, 4.4, 4.5);
        assert!(summit > 18e3 && summit < 20e3, "{summit}");
        assert!(foot > 0.0 && foot < 1.5e3, "{foot}");
        assert!(late < 0.01 * summit, "{late}");
        assert_eq!(
            shield_burial(shields, at(1.1 * shield.extent()), r, 0.9, 4.5),
            0.0
        );
    }

    #[test]
    fn a_shield_stretches_along_its_rift_and_lobes_at_its_foot() {
        // F101: the edifice is longer along its rift, round near its summit and lobed at
        // its foot; the eruptions that show came after its last few sheets were laid.
        let shield = Shield {
            center: [1.0, 0.0, 0.0],
            radius: 300e3,
            height: 20e3,
            start_gyr: 1.0,
            calderas: Vec::new(),
            rift: 0.0,
            elongation: 1.3,
            id: 2,
        };
        let pi = std::f64::consts::PI;
        let along = shield.radius_toward(0.0) + shield.radius_toward(pi);
        let across = shield.radius_toward(0.5 * pi) + shield.radius_toward(1.5 * pi);
        assert!(along > 1.1 * across, "{along} {across}");
        // Near the summit the share of the way to the foot follows the ellipse alone; at
        // the foot it varies with the lobing.
        let core = |phi: f64| shield.outline(phi).0;
        let (mut near, mut foot) = (0.0f64, 0.0f64);
        for k in 0..64 {
            let phi = 2.0 * pi * k as f64 / 64.0;
            near = near.max((shield.reach_share(0.1 * core(phi), phi) - 0.1).abs());
            foot = foot.max((shield.radius_toward(phi) / core(phi) - 1.0).abs());
            assert!(shield.radius_toward(phi) <= shield.extent());
        }
        assert!(near < 0.002 && foot > 0.03, "{near} {foot}");
        let (window, now) = shield_window(&shield, 4.5).unwrap();
        assert!(window > shield.start_gyr && window < now && now == 4.5);
    }
}
