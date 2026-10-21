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
use crate::crater::stamp::burial_factor;
use crate::crater::{CellWrite, StampOutput, Stamper};
use crate::fields::LayerSet;
use crate::geom::{
    angle_between, azimuth_from, cross, dot, local_basis, normalize, offset, rotate,
};
use crate::grid::CubeSphere;
use crate::grid::projection::PaddedField;
use crate::noise::GradientNoise;
use crate::params::{ParamGraph, Reads};
use crate::rng::Rng;
use crate::stages::io::{Lobe, PathRules};
use crate::thermal::ThermalHistory;
use crate::units::{UNIT_LAVA, UNIT_TESSERA};

pub const VERSION: u64 = 21;

/// A volcanic rise: a broad dome growing between two times.
#[derive(Clone, Debug)]
pub struct Rise {
    pub center: [f64; 3],
    pub radius: f64,
    pub height: f64,
    pub start_gyr: f64,
    pub end_gyr: f64,
    /// Its latitude when it began to rise, before the body turned under it (F104), rad.
    pub paleo_latitude: f64,
    /// The axis it grows along, a unit vector in the tangent plane at its centre (F114).
    pub axis: [f64; 3],
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
    /// The time over which it grows by a factor e of what it has still to grow, Gyr: its
    /// magma comes as fast as its mantle melts as it starts (F120).
    pub growth_gyr: f64,
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

/// What made a trough (F111): a graben, its faults meeting at the base of a faulted layer
/// near the surface; a rift, cutting deeper, as wide as its throw allows; a ring graben
/// round a basin (F69); a fracture round a corona.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TroughKind {
    Graben,
    Rift,
    Ring,
    Fracture,
}

/// One segment of a fault (F111): where it starts and ends along the trace, m, its
/// greatest throw, m, and for a rift's half-graben the side of its master fault, +1 or -1
/// along the trace's normal; 0 for a graben's symmetric floor.
#[derive(Clone, Copy, Debug)]
pub struct Segment {
    pub start: f64,
    pub end: f64,
    pub throw: f64,
    pub polarity: f64,
}

/// The side a trough's floor tilts down toward at a distance along it (F111), -1 to 1: a
/// rift segment's master fault's, turning over a transfer zone about each link; 0 for a
/// symmetric floor.
pub fn trough_polarity(trough: &Trough, s: f64) -> f64 {
    if trough.segments.is_empty() {
        return 0.0;
    }
    let n = trough.segments.len();
    let j = trough
        .segments
        .iter()
        .position(|g| s <= g.end)
        .unwrap_or(n - 1);
    let here = trough.segments[j].polarity;
    let zone = 0.5 * RIFT_TRANSFER_WIDTHS * trough.width;
    // Toward a link with a segment of the other polarity, the floor turns linearly to
    // level at the link itself.
    let blend = |other: f64, distance: f64| {
        if zone <= 0.0 || distance >= zone {
            here
        } else {
            let mid = 0.5 * (here + other);
            mid + (here - mid) * distance / zone
        }
    };
    let g = &trough.segments[j];
    if j > 0 && s - g.start < zone {
        blend(trough.segments[j - 1].polarity, s - g.start)
    } else if j + 1 < n && g.end - s < zone {
        blend(trough.segments[j + 1].polarity, g.end - s)
    } else {
        here
    }
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
    pub kind: TroughKind,
    /// A fault's segments in order along the trace (F111); none for a trough drawn as one
    /// profile tapering to its ends, a ring's arc or a fracture.
    pub segments: Vec<Segment>,
}

/// A trough's depth at a distance along its trace at its full growth, m. A fault's is its
/// segments' throw (F111), dipping by SEGMENT_LINK_DEFICIT where segments link and running
/// on unbroken from one segment's throw to the next's. A graben's faults are confined to
/// its layer, so its throw stays up along its length and falls only toward its tips, as a
/// flat-topped profile (F112); a rift's falls toward its ends as its whole profile does. A
/// trough without segments tapers to its ends as a sine to TROUGH_TAPER_POWER.
pub fn trough_throw(trough: &Trough, s: f64) -> f64 {
    let whole = (std::f64::consts::PI * s / trough.length.max(1.0))
        .sin()
        .max(0.0)
        .powf(TROUGH_TAPER_POWER);
    if trough.segments.is_empty() {
        return trough.depth * whole;
    }
    let n = trough.segments.len();
    let j = trough
        .segments
        .iter()
        .position(|g| s <= g.end)
        .unwrap_or(n - 1);
    let g = &trough.segments[j];
    let x = (2.0 * (s - g.start) / (g.end - g.start).max(1.0) - 1.0).clamp(-1.0, 1.0);
    let ellipse = (1.0 - x * x).max(0.0).sqrt();
    // The segment's throw turns toward the mean of its own and its neighbour's at a link.
    let throw = if x < 0.0 && j > 0 {
        g.throw + 0.5 * (trough.segments[j - 1].throw - g.throw) * -x
    } else if x > 0.0 && j + 1 < n {
        g.throw + 0.5 * (trough.segments[j + 1].throw - g.throw) * x
    } else {
        g.throw
    };
    let outer = (j == 0 && x < 0.0) || (j + 1 == n && x > 0.0);
    let link = if outer {
        1.0
    } else {
        1.0 - SEGMENT_LINK_DEFICIT * (1.0 - ellipse)
    };
    if trough.kind == TroughKind::Rift {
        return throw * whole * link;
    }
    // Confined to its layer, a graben keeps its throw along its length, rising from each
    // tip as fast as twice its throw per length allows.
    let tip = (trough.depth / (2.0 * SEGMENT_THROW_PER_LENGTH))
        .min(0.5 * trough.length)
        .max(1.0);
    let d = s.min(trough.length - s).max(0.0);
    let rise = if d >= tip {
        1.0
    } else {
        let u = d / tip;
        (1.0 - (1.0 - u) * (1.0 - u)).sqrt()
    };
    throw * rise * link
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

/// The fluid Love number of a body with a moment of inertia factor (the Darwin-Radau
/// relation): 3/2 for a uniform body, about 0.94 for the Earth.
pub fn fluid_love_number(moment_of_inertia: f64) -> f64 {
    let s = 2.5 * (1.0 - 1.5 * moment_of_inertia);
    (4.0 - s * s) / (1.0 + s * s)
}

/// The share of a load at degree two that an elastic lithosphere `te` thick holds up
/// rather than isostasy, by its membrane strength (Turcotte et al. 1981; its bending
/// strength is negligible at degree two): a few percent for Mars's.
pub fn degree_two_support(te: f64, gravity: f64, body_radius: f64) -> f64 {
    let n = 6.0;
    let xi = ROCK_YOUNG_MODULUS * te.max(0.0)
        / (DENSITY_ROCK * gravity.max(1e-3) * body_radius * body_radius)
        * (n - 2.0)
        / (n - (1.0 - ROCK_POISSON));
    xi / (1.0 + xi)
}

/// The rotational bulge an elastic lithosphere holds frozen as the body turns under it,
/// kg m²: the excess of the polar moment of inertia over the equatorial, (k_f / 3)
/// ω² R⁵ / G for a fluid body, times the share the lithosphere supports (Willemann 1984;
/// Matsuyama et al. 2006).
pub fn frozen_bulge(
    body_radius: f64,
    gravity: f64,
    angular_velocity: f64,
    moment_of_inertia: f64,
    te: f64,
) -> f64 {
    let fluid =
        fluid_love_number(moment_of_inertia) / 3.0 * angular_velocity.powi(2) * body_radius.powi(5)
            / G;
    fluid * degree_two_support(te, gravity, body_radius)
}

/// True polar wander (F104): the turn that brings the body's axis of greatest moment of
/// inertia to its spin axis, from the rises' mass and the frozen bulge about the axis
/// they began under. A plume holds its rise up rather than isostasy, so its mass turns
/// the body against only the bulge the elastic lithosphere holds frozen, and a large
/// rise is carried to the equator (Matsuyama et al. 2006: a dynamic planet's stability
/// is the remnant bulge's; the Earth's two great upwellings lie under its equator; Bouley
/// et al. 2016 on Tharsis's turn of Mars); several rise to a great circle. Returns the
/// axis and angle of the rotation that carries the body-fixed directions of before into
/// the frame of the present spin axis, the shortest turn that does.
pub fn polar_wander(rises: &[Rise], body_radius: f64, bulge: f64) -> ([f64; 3], f64) {
    let pole = [0.0, 0.0, 1.0];
    // The traceless inertia tensor: the bulge about the old pole, less each rise's mass
    // m R² (u uᵀ - I/3) summed over its dome.
    let mut tensor = [[0.0f64; 3]; 3];
    for (i, row) in tensor.iter_mut().enumerate() {
        for (j, t) in row.iter_mut().enumerate() {
            *t = bulge * (pole[i] * pole[j] - if i == j { 1.0 / 3.0 } else { 0.0 });
        }
    }
    let (rings, spokes) = (48, 48);
    for rise in rises {
        let (east, north) = local_basis(rise.center);
        let dr = rise.radius / rings as f64;
        let dphi = 2.0 * std::f64::consts::PI / spokes as f64;
        for a in 0..rings {
            let r = (a as f64 + 0.5) * dr;
            let mass = DENSITY_CRUST_ROCK * rise.height * dome(r / rise.radius) * r * dr * dphi;
            for b in 0..spokes {
                let u = offset(rise.center, east, north, r / body_radius, b as f64 * dphi);
                for (i, row) in tensor.iter_mut().enumerate() {
                    for (j, t) in row.iter_mut().enumerate() {
                        let iso = if i == j { 1.0 / 3.0 } else { 0.0 };
                        *t -= mass * body_radius * body_radius * (u[i] * u[j] - iso);
                    }
                }
            }
        }
    }
    let (values, vectors) = crate::geom::symmetric_eigen(tensor);
    // The new spin axis: the greatest axis, or where two tie, the one of their plane
    // nearest the old pole.
    let spin = if values[0] - values[1] <= 1e-9 * values[0].abs().max(f64::MIN_POSITIVE) {
        let (a, b) = (dot(pole, vectors[0]), dot(pole, vectors[1]));
        let v = [
            a * vectors[0][0] + b * vectors[1][0],
            a * vectors[0][1] + b * vectors[1][1],
            a * vectors[0][2] + b * vectors[1][2],
        ];
        if dot(v, v) > 1e-24 {
            normalize(v)
        } else {
            vectors[0]
        }
    } else {
        vectors[0]
    };
    let spin = if dot(spin, pole) < 0.0 {
        [-spin[0], -spin[1], -spin[2]]
    } else {
        spin
    };
    let axis = cross(spin, pole);
    let size = dot(axis, axis).sqrt();
    if size < 1e-12 {
        return (pole, 0.0);
    }
    let axis = [axis[0] / size, axis[1] / size, axis[2] / size];
    (axis, dot(spin, pole).clamp(-1.0, 1.0).acos())
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
pub fn shield_progress(t: f64, shield: &Shield) -> f64 {
    if t <= shield.start_gyr {
        0.0
    } else {
        1.0 - (-(t - shield.start_gyr) / shield.growth_gyr.max(1e-6)).exp()
    }
}

/// The share of a rise's dikes its shields' loads gather by a time (F119), and so the share
/// of its eruptions its shields take (F108). A rising dike bends toward a volcano's load,
/// from several of its widths away where the load's pressure outweighs the regional
/// differential stress, a feedback between the edifice and the ground it draws dikes from
/// (Muller, Ito and Martel 2001). While a rise's dome grows it stretches its lithosphere,
/// its faults and dike swarms taking up the stretching (F105, F113, F116), and its dikes
/// run along that stress to erupt as fissure floods; as the growth slows the stress falls
/// and the loads gather the dikes. The stress follows the growth's rate, as the faults do,
/// and the shields gather the share by which it has fallen from its peak (a stated choice,
/// for the comparison of pressures whose magnitudes the loop does not carry): none while
/// the growth quickens, all once it has stopped (Carr and Head 2010: Mars's volcanism
/// gathered at the large shields of Tharsis and Elysium in the Amazonian, after the
/// widespread plains volcanism of the Hesperian).
pub fn rise_focus(rise: &Rise, t: f64) -> f64 {
    let span = rise.end_gyr - rise.start_gyr;
    if span <= 0.0 {
        return if t >= rise.end_gyr { 1.0 } else { 0.0 };
    }
    let x = ((t - rise.start_gyr) / span).clamp(0.0, 1.0);
    if x <= 0.5 {
        0.0
    } else {
        // One less the growth's rate against its peak, 4 x (1 - x).
        (2.0 * x - 1.0).powi(2)
    }
}

/// When a rise's shield starts (F119): once its rise's shields gather as much of its dikes
/// as the larger shields before it hold, a share of their volume, so that the largest
/// starts as the dome's growth turns and the rest follow as the stretching wanes.
pub fn shield_start(rise: &Rise, before: f64) -> f64 {
    let x = 0.5 * (1.0 + before.clamp(0.0, 1.0).sqrt());
    rise.start_gyr + x * (rise.end_gyr - rise.start_gyr)
}

/// Where a dike swarm starts (F113, F116): its source's centre and radius, how far its
/// dikes reach, its random key, and the volcano that feeds it, if one does.
struct SwarmSource<'a> {
    center: [f64; 3],
    radius: f64,
    reach: f64,
    key: u64,
    volcano: Option<&'a Shield>,
}

/// The rise a point belongs to: the nearest in units of each rise's radius.
fn rise_of(rises: &[Rise], d: [f64; 3]) -> usize {
    rises
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| {
            (angle_between(a.center, d) / a.radius)
                .total_cmp(&(angle_between(b.center, d) / b.radius))
        })
        .map_or(0, |(i, _)| i)
}

/// Where a rise's eruptions go in a step (F108): its dome's height over the tallest
/// rise's at each cell it covers, and the share of its eruptions its shields take.
pub struct RiseVent {
    pub dome: Vec<(usize, f32)>,
    pub capture: f64,
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
    shield.start_gyr - shield.growth_gyr * (1.0 - progress).ln()
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
    let now = shield_progress(present, shield);
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
            0.0,
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

/// The rises a body's plumes raise (F104, F120): where they stand once the body has turned
/// under them, their size and when they grow. None when the body raises no plume.
pub fn place_rises<P: Reads>(graph: &P, thermal: &ThermalHistory) -> Option<Vec<Rise>> {
    let count = graph.num("tectonics.plume_rise_count").round() as usize;
    if count == 0
        || !matches!(
            graph.text("interior.tectonic_regime"),
            "stagnant_lid" | "episodic_lid"
        )
    {
        return None;
    }
    let body_radius = graph.num("body.radius");
    let age = graph.num("body.age");
    // A rise grows as its plume's mantle melts beneath its lid, since the lid last
    // formed (F120): a stagnant lid when the body formed, an episodic lid at its
    // overturn, which foundered the lid and the domes on it.
    let renewal = if graph.text("interior.tectonic_regime") == "episodic_lid" {
        (age - graph.num("tectonics.resurfacing_event_age")).max(0.0)
    } else {
        0.0
    };
    let clock = thermal.melt_clock(4.0 * std::f64::consts::PI * body_radius * body_radius);
    let since = clock.at(renewal);
    let rise_start = clock.time_of(since + RISE_START_MELT);
    let rise_end = clock.time_of(since + RISE_END_MELT);
    let rng = Rng::new(graph.num("seed.tectonics") as u64 ^ 0x9D1F_E5A7);
    let rise_radius = 0.5 * graph.num("tectonics.plume_rise_diameter");
    let rise_height = graph.num("tectonics.plume_rise_height");
    let mut rises = Vec::new();
    for i in 0..count as u64 {
        let mut center = rng.unit_direction(i, 1);
        // Rises keep apart.
        for attempt in 0..16u64 {
            if rises
                .iter()
                .all(|r: &Rise| angle_between(r.center, center) > 2.5 * rise_radius / body_radius)
            {
                break;
            }
            center = rng.unit_direction(i, 10 + 2 * attempt);
        }
        rises.push(Rise {
            center,
            radius: rise_radius * (0.8 + 0.4 * rng.uniform(i, 4)),
            height: rise_height * (0.8 + 0.4 * rng.uniform(i, 5)),
            start_gyr: rise_start,
            end_gyr: rise_end,
            paleo_latitude: center[2].clamp(-1.0, 1.0).asin(),
            axis: [0.0, 0.0, 0.0],
        });
    }
    // The body turns under its rises as they grow, carrying them toward its equator
    // (F104); everything on them is placed where they stand now.
    {
        let mid = (0.5 * (rise_start + rise_end)).min(age);
        let lithosphere_now = graph.num("interior.lithosphere_thickness").max(1.0);
        let te = graph.num("interior.elastic_thickness")
            * graph
                .series("interior.lithosphere_history")
                .map_or(lithosphere_now, |c| c.eval(mid))
            / lithosphere_now;
        let bulge = frozen_bulge(
            body_radius,
            graph.num("body.gravity"),
            graph.num("rotation.angular_velocity"),
            graph.num("body.moment_of_inertia_factor"),
            te,
        );
        let (axis, angle) = polar_wander(&rises, body_radius, bulge);
        for rise in &mut rises {
            rise.center = normalize(rotate(rise.center, axis, angle));
        }
    }
    Some(rises)
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
    pub fn new(graph: &ParamGraph, grid: CubeSphere, thermal: &ThermalHistory) -> Option<Self> {
        let mut rises = place_rises(graph, thermal)?;
        let body_radius = graph.f64("body.radius");
        let age = graph.f64("body.age");
        let clock = thermal.melt_clock(4.0 * std::f64::consts::PI * body_radius * body_radius);
        let rng = Rng::new(graph.f64("seed.tectonics") as u64 ^ 0x9D1F_E5A7);
        let rise_height = graph.f64("tectonics.plume_rise_height");
        // Each rise grows along an axis (F114), its azimuth drawn as its first dikes break
        // the symmetry of its crest.
        let axes: Vec<(f64, Axial)> = rises
            .iter_mut()
            .enumerate()
            .map(|(ri, rise)| {
                let (east, north) = local_basis(rise.center);
                let phi = std::f64::consts::TAU * rng.uniform(0x0A21_5000 + ri as u64, 0);
                let axis = normalize([
                    east[0] * phi.cos() + north[0] * phi.sin(),
                    east[1] * phi.cos() + north[1] * phi.sin(),
                    east[2] * phi.cos() + north[2] * phi.sin(),
                ]);
                rise.axis = axis;
                (
                    phi,
                    Axial {
                        axis,
                        stress: RISE_AXIAL_STRESS,
                    },
                )
            })
            .collect();
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
            // The largest first, then down a power law.
            let size = |k: usize| {
                let rank = k as f64 + 1.0;
                let diameter = SHIELD_MAX_DIAMETER_FRACTION * 2.0 * rise.radius / rank.powf(0.5);
                let radius = 0.5 * diameter;
                (
                    radius,
                    (radius * SHIELD_SLOPE_DEG.to_radians().tan()).min(relief_limit),
                )
            };
            // Each starts once the rise's shields gather what the larger ones hold (F119).
            let volume = |k: usize| {
                let (r, h) = size(k);
                r * r * h
            };
            let total: f64 = (0..n_here).map(volume).sum();
            let mut before = 0.0;
            for k in 0..n_here as u64 {
                let draw = 100 + sid * 16;
                let (radius, height) = size(k as usize);
                let mut center = rise.center;
                for attempt in 0..24u64 {
                    let dist = rise.radius
                        * (SHIELD_PLACEMENT.0
                            + (SHIELD_PLACEMENT.1 - SHIELD_PLACEMENT.0)
                                * rng.uniform(draw, 1 + 3 * attempt));
                    // On the rise's axis, on either side of its centre (F114).
                    let side = if rng.uniform(draw, 2 + 3 * attempt) < 0.5 {
                        0.0
                    } else {
                        std::f64::consts::PI
                    };
                    let spread = SHIELD_AXIS_SPREAD_DEG.to_radians()
                        * rng.normal(draw ^ 0xA715_0000_0000, 2 * attempt);
                    let phi = axes[ri].0 + side + spread;
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
                let start_gyr = shield_start(rise, before / total.max(f64::MIN_POSITIVE));
                shields.push(Shield {
                    center,
                    radius,
                    height,
                    start_gyr,
                    calderas,
                    rift,
                    elongation: 1.0 + (SHIELD_ELONGATION_MAX - 1.0) * rng.uniform(draw, 15),
                    growth_gyr: SHIELD_GROWTH_MELT / clock.rate(start_gyr.min(age)).max(1e-9),
                    id: sid,
                });
                before += volume(k as usize);
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
        // Where faults start: as often as the ground is in tension, over the area round a
        // rise (F105), by the total stress (F114).
        let starts = StartTable::new(RISE_AXIAL_STRESS);
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
            // Faults form as the rise's stress grows (F105, F116): as its dome rises, the
            // dome's bending and load being that stress, at the pace of its growth and none
            // once it has stopped, on every lid, the dome growing since its lid last formed
            // (F120); a rise still growing has stretched by its growth so far, and one whose
            // plume has yet to rise not at all.
            if rise.start_gyr >= age {
                continue;
            }
            let budget = budget * rise_progress(rise, age);
            let (t_lo, t_hi) = (rise.start_gyr, rise.end_gyr.min(age));
            let time_cdf = Tabulated::new(t_lo, t_hi, 256, |t| {
                // The rate of the dome's smoothstep growth.
                let x = ((t - rise.start_gyr) / (rise.end_gyr - rise.start_gyr)).clamp(0.0, 1.0);
                6.0 * x * (1.0 - x)
            });
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
                // The faulted layer near the surface: a graben is as wide as the layer
                // sets it, and a fault whose throw over its whole length would exceed the
                // layer cuts deeper, a rift as wide as its throw (F111).
                let layer_width = (GRABEN_WIDTH
                    * (GRABEN_WIDTH_SPREAD * rng.normal(draw ^ 0x6A4E_0000_0000, 0)).exp())
                .clamp(GRABEN_WIDTH_RANGE.0, GRABEN_WIDTH_RANGE.1);
                let layer = (FAULT_LAYER_DEPTH
                    * (FAULT_LAYER_SPREAD * rng.normal(draw ^ 0x6A4E_0000_0000, 2)).exp())
                .clamp(FAULT_LAYER_RANGE.0, FAULT_LAYER_RANGE.1);
                let zone_depth = |length: f64| {
                    (RIFT_DEPTH_PER_LENGTH * length * scatter).clamp(RIFT_MIN_DEPTH, depth_cap)
                };
                let class_width = |length: f64| {
                    let d = zone_depth(length);
                    if d >= layer {
                        (RIFT_WIDTH_PER_DEPTH * d).max(RIFT_MIN_WIDTH)
                    } else {
                        layer_width
                    }
                };
                let width = class_width(length);
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
                let (points, breaks) = match host {
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
                        // Graben in a set stand apart by a few thicknesses of their faulted
                        // layer, wall to wall (F112); chasmata by a width or a few.
                        let distance = if hf.kind == TroughKind::Graben {
                            let spread =
                                GRABEN_SPACING_SPREAD * rng.normal(draw ^ 0x5BAC_0000_0000, 0);
                            let gap = layer * (GRABEN_SPACING_LAYERS * spread.exp()).max(0.25);
                            0.5 * (hf.width + width) + gap
                        } else {
                            let (g_lo, g_hi) = RIFT_CLUSTER_OFFSET;
                            let gap = g_lo + (g_hi - g_lo) * rng.uniform(draw, 7).powi(2);
                            gap * hf.width.max(width)
                        };
                        length = span as f64 * hf.spacing;
                        // The neighbour follows its host's segments along the span it runs
                        // beside.
                        let (s0, s1) = (i0 as f64 * hf.spacing, (i0 + span) as f64 * hf.spacing);
                        let breaks: Vec<f64> = hf
                            .segments
                            .iter()
                            .map(|g| g.end)
                            .filter(|&e| e > s0 && e < s1)
                            .map(|e| e - s0)
                            .collect();
                        (
                            offset_trace(&hf.points[i0..i0 + span], side * distance, body_radius),
                            breaks,
                        )
                    }
                    None => {
                        // A fault starts where the ground is in tension, as often as the
                        // tension is strong, and runs as the stress field leads it (F105),
                        // in segments (F111).
                        let (x, a) = starts.sample(
                            rng.uniform(draw, 20),
                            rng.uniform(draw, 22),
                            rng.uniform(draw, 23),
                        );
                        // The axis lies at `a` from the direction outward at the start.
                        let phi = axes[ri].0 - a;
                        let start =
                            offset(rise.center, east, north, x * rise.radius / body_radius, phi);
                        let (trace, breaks) = stress_trace(
                            rise,
                            start,
                            length,
                            spacing,
                            body_radius,
                            &rng,
                            draw,
                            width,
                            None,
                            RIFT_MIN_TENSION,
                            Some(axes[ri].1),
                        );
                        length = (trace.len().max(2) - 1) as f64 * spacing;
                        let breaks = breaks.into_iter().filter(|&b| b < length).collect();
                        (trace, breaks)
                    }
                };
                // A fault that died out early slipped as its own length allows; it is a
                // rift if that throw exceeds the faulted layer, else a graben whose
                // segments each throw as their own length allows (F111).
                let zone = zone_depth(length);
                let rift = zone >= layer;
                let width = if rift {
                    (RIFT_WIDTH_PER_DEPTH * zone).max(RIFT_MIN_WIDTH)
                } else {
                    layer_width
                };
                let mut bounds = vec![0.0];
                bounds.extend(breaks);
                bounds.push(length);
                // A rift's segments are half-graben whose master fault's side flips now
                // and then from one to the next.
                let mut polarity = if rng.uniform(draw ^ 0x5E67_0000_0000, 1) < 0.5 {
                    -1.0
                } else {
                    1.0
                };
                let segments: Vec<Segment> = bounds
                    .windows(2)
                    .enumerate()
                    .map(|(j, w)| {
                        // A graben confined to its layer throws as far as its segments'
                        // length allows, linked or not, along its whole length (F112), within
                        // its layer and its width, each segment a little more or less.
                        let throw = if rift {
                            zone
                        } else {
                            let cap = layer.min(GRABEN_DEPTH_PER_WIDTH * width);
                            let mean_segment = length / (bounds.len() - 1).max(1) as f64;
                            let plateau = (SEGMENT_THROW_PER_LENGTH * mean_segment * scatter)
                                .clamp(RIFT_MIN_DEPTH, cap);
                            let own = (SEGMENT_THROW_SCATTER
                                * rng.normal(draw ^ 0x5E67_0000_0000, 2 * j as u64))
                            .exp();
                            (plateau * own).clamp(RIFT_MIN_DEPTH, cap)
                        };
                        if j > 0
                            && rng.uniform(draw ^ 0x5E67_0000_0000, 2 * j as u64 + 1)
                                < RIFT_POLARITY_FLIP
                        {
                            polarity = -polarity;
                        }
                        Segment {
                            start: w[0],
                            end: w[1],
                            throw,
                            polarity: if rift { polarity } else { 0.0 },
                        }
                    })
                    .collect();
                let depth = segments.iter().map(|g| g.throw).fold(0.0, f64::max);
                let time = time_cdf.sample(rng.uniform(draw, 8));
                let carried = width >= RIFT_LOOP_WIDTH_CELLS * cell;
                let (start_gyr, end_gyr) = if carried {
                    // Carved over the melt its mantle makes about its time (F120).
                    let m = clock.at(time);
                    (
                        clock.time_of(m - 0.5 * RIFT_CARVING_MELT).max(t_lo),
                        clock.time_of(m + 0.5 * RIFT_CARVING_MELT).min(age),
                    )
                } else {
                    (time, time)
                };
                // The extension a fault takes up: its heave along its length, each segment's
                // mean throw over it.
                let mean = 1.0 - SEGMENT_LINK_DEFICIT * (1.0 - std::f64::consts::FRAC_PI_4);
                spent += RIFT_HEAVE_PER_DEPTH
                    * mean
                    * segments
                        .iter()
                        .map(|g| g.throw * (g.end - g.start))
                        .sum::<f64>();
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
                    kind: if rift {
                        TroughKind::Rift
                    } else {
                        TroughKind::Graben
                    },
                    segments,
                });
            }
            // Radiating dike swarms (F113): the plume's head and each volcano on the rise
            // send dikes out laterally through the crust, each on its own heading within
            // its swarm's fan and on along the stress trajectories as far as its magma
            // carries it, magma-driven whatever the tension; a graben opens above each, as
            // wide as twice the depth to the dike's top, subsiding as a graben over its
            // faulted layer does.
            // A swarm comes with the magma that feeds it (F116): the plume head's as the
            // dome rises, a volcano's as it grows.
            let mut sources: Vec<SwarmSource> = Vec::new();
            for m in 0..DIKE_SWARMS_PER_RISE {
                let key = 0x0D1C_0000 + ri as u64 * 0x1000 + m;
                let (q_lo, q_hi) = DIKE_RISE_REACH;
                let reach = rise.radius * q_lo * (q_hi / q_lo).powf(rng.uniform(key, 0));
                sources.push(SwarmSource {
                    center: rise.center,
                    radius: rise.radius,
                    reach,
                    key,
                    volcano: None,
                });
            }
            for s in shields.iter().filter(|s| rise_of(&rises, s.center) == ri) {
                let key = 0x0D1D_0000 + s.id;
                let (q_lo, q_hi) = DIKE_SHIELD_REACH;
                let reach = s.radius * q_lo * (q_hi / q_lo).powf(rng.uniform(key, 0));
                sources.push(SwarmSource {
                    center: s.center,
                    radius: s.radius,
                    reach,
                    key,
                    volcano: Some(s),
                });
            }
            for &SwarmSource {
                center: source,
                radius: source_radius,
                reach,
                key,
                volcano,
            } in &sources
            {
                let time = match volcano {
                    Some(v) if v.start_gyr < age => {
                        let grown = shield_progress(age, v);
                        progress_time(v, rng.uniform(key, 1) * grown).min(age)
                    }
                    _ => time_cdf.sample(rng.uniform(key, 1)),
                };
                // The swarm fans along the rise's axis, to one side or the other of its
                // source: its magma opens its way most easily across the axis (F114).
                let (se, sn) = local_basis(source);
                let along = axis_at(rise.center, axes[ri].1.axis, source);
                let side = if rng.uniform(key, 2) < 0.5 {
                    0.0
                } else {
                    std::f64::consts::PI
                };
                let axis = dot(along, sn).atan2(dot(along, se)) + side;
                let (w_lo, w_hi) = DIKE_SWARM_SECTOR_DEG;
                let sector = (w_lo * (w_hi / w_lo).powf(rng.uniform(key, 3))).to_radians();
                let (n_lo, n_hi) = DIKES_PER_SWARM;
                let count = (n_lo * (n_hi / n_lo).powf(rng.uniform(key, 4))).round() as u64;
                for d in 0..count {
                    let draw = key * 0x100 + d;
                    let phi = axis + sector * (rng.uniform(draw, 0) - 0.5);
                    let (l_lo, l_hi) = DIKE_LENGTH_SHARE;
                    let length = reach * l_lo * (l_hi / l_lo).powf(rng.uniform(draw, 1));
                    let (o_lo, o_hi) = DIKE_START_SHARE;
                    let from = source_radius * (o_lo + (o_hi - o_lo) * rng.uniform(draw, 2));
                    let start = offset(source, se, sn, from / body_radius, phi);
                    let along = dot(source, start);
                    let out = normalize([
                        start[0] * along - source[0],
                        start[1] * along - source[1],
                        start[2] * along - source[2],
                    ]);
                    let width = (GRABEN_WIDTH
                        * (GRABEN_WIDTH_SPREAD * rng.normal(draw ^ 0x6A4E_0000_0000, 0)).exp())
                    .clamp(GRABEN_WIDTH_RANGE.0, GRABEN_WIDTH_RANGE.1);
                    let layer = (FAULT_LAYER_DEPTH
                        * (FAULT_LAYER_SPREAD * rng.normal(draw ^ 0x6A4E_0000_0000, 2)).exp())
                    .clamp(FAULT_LAYER_RANGE.0, FAULT_LAYER_RANGE.1);
                    let (trace, breaks) = stress_trace(
                        rise,
                        start,
                        length,
                        spacing,
                        body_radius,
                        &rng,
                        draw,
                        width,
                        Some(out),
                        f64::NEG_INFINITY,
                        Some(axes[ri].1),
                    );
                    let length = (trace.len().max(2) - 1) as f64 * spacing;
                    let mut bounds = vec![0.0];
                    bounds.extend(breaks.into_iter().filter(|&b| b < length));
                    bounds.push(length);
                    let cap = layer.min(GRABEN_DEPTH_PER_WIDTH * width);
                    let mean_segment = length / (bounds.len() - 1).max(1) as f64;
                    let scatter = (RIFT_DEPTH_SCATTER * rng.normal(draw, 4)).exp();
                    let plateau = (SEGMENT_THROW_PER_LENGTH * mean_segment * scatter)
                        .clamp(RIFT_MIN_DEPTH, cap);
                    let segments: Vec<Segment> = bounds
                        .windows(2)
                        .enumerate()
                        .map(|(j, w)| {
                            let own = (SEGMENT_THROW_SCATTER
                                * rng.normal(draw ^ 0x5E67_0000_0000, 2 * j as u64))
                            .exp();
                            Segment {
                                start: w[0],
                                end: w[1],
                                throw: (plateau * own).clamp(RIFT_MIN_DEPTH, cap),
                                polarity: 0.0,
                            }
                        })
                        .collect();
                    faults.push(Trough {
                        points: trace,
                        spacing,
                        length,
                        width,
                        depth: segments.iter().map(|g| g.throw).fold(0.0, f64::max),
                        flank: 0.0,
                        start_gyr: time,
                        end_gyr: time,
                        id: 0x0D1C_0000_0000_0000 | draw,
                        kind: TroughKind::Graben,
                        segments,
                    });
                }
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

    /// Where each rise's eruptions go for the vents a step draws between two times
    /// (F108): its dome as it stands at the step's end, and the share of its eruptions
    /// its shields take at the step's middle.
    pub fn rise_vents(&mut self, t0: f64, t1: f64) -> Vec<RiseVent> {
        let grid = self.grid;
        let body_radius = self.body_radius;
        let mut out = Vec::with_capacity(self.rises.len());
        for (ri, rise) in self.rises.iter().enumerate() {
            let grown = rise_progress(rise, t1);
            self.stamper
                .gather(&grid, body_radius, rise.center, rise.radius);
            let dome_cells: Vec<(usize, f32)> = self
                .stamper
                .footprint
                .iter()
                .map(|&(idx, r, _)| {
                    let d = rise.height * grown * dome(r / rise.radius) / self.rise_height;
                    (idx, d as f32)
                })
                .filter(|&(_, d)| d > 0.0)
                .collect();
            let mine: Vec<&Shield> = self
                .shields
                .iter()
                .filter(|s| rise_of(&self.rises, s.center) == ri)
                .collect();
            out.push(RiseVent {
                dome: dome_cells,
                capture: if mine.is_empty() {
                    0.0
                } else {
                    rise_focus(rise, 0.5 * (t0 + t1))
                },
            });
        }
        out
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
        // Shields: flows pile up, renewing the surface where they add more than a veneer,
        // and burying what they flow over (F115).
        let t_mid = 0.5 * (t0 + t1);
        for shield in &self.shields {
            let dp = shield_progress(t1, shield) - shield_progress(t0, shield);
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
            // The ground the edifice grows on, in rings about its centre: the surface less
            // the lava the shields have piled there and less the canyons' open depth, which
            // their lava fills as it comes (F103). The lava buries that ground's relief,
            // measured against its ring's mean, as a deposit buries what it falls on: the
            // lows fill first and the highs drown as it thickens (F115). A crater the
            // shield grows over fills and vanishes under it, rather than riding up on it.
            let extent = shield.extent();
            let ring_of = |r: f64| {
                ((r / extent * SHIELD_BURIAL_RINGS as f64) as usize).min(SHIELD_BURIAL_RINGS - 1)
            };
            let mut ring_sum = [0.0f64; SHIELD_BURIAL_RINGS];
            let mut ring_count = [0.0f64; SHIELD_BURIAL_RINGS];
            for &(idx, r, _) in &self.stamper.footprint {
                let under =
                    bedrock[idx] - self.shield_built[idx] + (self.carved[idx] - self.filled[idx]);
                ring_sum[ring_of(r)] += under as f64;
                ring_count[ring_of(r)] += 1.0;
            }
            for &(idx, r, _) in &self.stamper.footprint {
                let d = grid.center(grid.cell(idx));
                let full = shield_height(shield, d, body_radius, false, 0.0);
                let dh = full * dp;
                if dh <= 0.0 {
                    continue;
                }
                let under =
                    bedrock[idx] - self.shield_built[idx] + (self.carved[idx] - self.filled[idx]);
                let relief = under as f64 - ring_sum[ring_of(r)] / ring_count[ring_of(r)].max(1.0);
                let drowned = relief * (1.0 - burial_factor(relief, dh));
                // The lava ponds in a canyon it comes to, filling it as fast as the
                // shield thickens round it (F103).
                let fill = (self.carved[idx] - self.filled[idx]).clamp(0.0, dh as f32);
                self.filled[idx] += fill;
                bedrock[idx] += (dh - drowned) as f32 + fill;
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

/// A distribution over an interval tabulated from its density, sampled by inverse transform.
struct Tabulated {
    lo: f64,
    hi: f64,
    cumulative: Vec<f64>,
}

impl Tabulated {
    fn new(lo: f64, hi: f64, bins: usize, density: impl Fn(f64) -> f64) -> Self {
        let width = (hi - lo) / bins as f64;
        let mut cumulative = Vec::with_capacity(bins + 1);
        let mut sum = 0.0;
        cumulative.push(0.0);
        for k in 0..bins {
            sum += density(lo + (k as f64 + 0.5) * width).max(0.0) * width;
            cumulative.push(sum);
        }
        Self { lo, hi, cumulative }
    }

    /// The point below which a share `u` of the distribution lies; uniform over the
    /// interval if the density is nowhere positive.
    fn sample(&self, u: f64) -> f64 {
        let total = *self.cumulative.last().unwrap_or(&0.0);
        let bins = self.cumulative.len() - 1;
        if total <= 0.0 || bins == 0 {
            return self.lo + (self.hi - self.lo) * u;
        }
        let target = u.clamp(0.0, 1.0) * total;
        let k = self
            .cumulative
            .partition_point(|&c| c <= target)
            .clamp(1, bins);
        let (c0, c1) = (self.cumulative[k - 1], self.cumulative[k]);
        let f = if c1 > c0 {
            (target - c0) / (c1 - c0)
        } else {
            0.5
        };
        self.lo + (self.hi - self.lo) * ((k - 1) as f64 + f) / bins as f64
    }
}

/// The surface stress about a rise at a distance from its centre in rise radii (F105):
/// radial and hoop, tension positive, in units of the tension at its crest. The
/// lithosphere bends over the plume as the dome does: tension in every direction on the
/// crest, in hoop on the flanks, where radial graben open, and compression at the foot,
/// radial the most, where the ridges stand concentric; the bending stresses go as the
/// dome's curvatures, -(w'' + ν w'/r) radially and -(w'/r + ν w'') in hoop. And the
/// rise's load stretches the ground beyond its foot in hoop and squeezes it radially,
/// as a load on a spherical shell does (Banerdt et al. 1992), by RIFT_FAR_HOOP of the
/// crest's tension at the foot, falling as the square of the distance.
pub fn rise_stress(x: f64) -> (f64, f64) {
    let nu = ROCK_POISSON;
    let crest = 4.0 * (1.0 + nu);
    let (radial, hoop) = if x < 1.0 {
        // dome(x) = (1 - x²)²: w' = -4x (1 - x²), w'' = -4 (1 - 3x²).
        let w1x = -4.0 * (1.0 - x * x);
        let w2 = -4.0 * (1.0 - 3.0 * x * x);
        (-(w2 + nu * w1x) / crest, -(w1x + nu * w2) / crest)
    } else {
        (0.0, 0.0)
    };
    let far = RIFT_FAR_HOOP * smooth(x, 0.7, 1.0) / x.max(1.0).powi(2);
    (radial - far, hoop + far)
}

/// The stress a rise's growth along its axis adds (F114): the axis, a unit vector in the
/// tangent plane at the rise's centre, and the tension across it and the compression along
/// it, a share of the crest's tension, fading beyond the foot as the load's own stress
/// does, as the inverse square of the distance.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Axial {
    pub axis: [f64; 3],
    pub stress: f64,
}

/// A rise's axis carried from its centre along the great circle to a point (F114): a unit
/// vector in the tangent plane there, keeping its shares along and across that circle.
pub(crate) fn axis_at(center: [f64; 3], axis: [f64; 3], p: [f64; 3]) -> [f64; 3] {
    let along = dot(center, p);
    // At the centre, toward the point; at the point, away from the centre.
    let toward = [
        p[0] - along * center[0],
        p[1] - along * center[1],
        p[2] - along * center[2],
    ];
    let out = [
        p[0] * along - center[0],
        p[1] * along - center[1],
        p[2] * along - center[2],
    ];
    if dot(toward, toward) < 1e-20 || dot(out, out) < 1e-20 {
        let r = dot(axis, p);
        return normalize([axis[0] - r * p[0], axis[1] - r * p[1], axis[2] - r * p[2]]);
    }
    let (toward, out) = (normalize(toward), normalize(out));
    let round = normalize(cross(p, out));
    let (a, b) = (dot(axis, toward), dot(axis, round));
    normalize([
        out[0] * a + round[0] * b,
        out[1] * a + round[1] * b,
        out[2] * a + round[2] * b,
    ])
}

/// The total surface stress about a rise at a point (F105, F114), tension positive in
/// units of the crest's tension, in the frame of the direction outward from the rise's
/// centre and the hoop direction round it there: radial, hoop and their shear, and the
/// frame's two unit vectors. The axis is carried from the centre along the great circle to
/// the point, keeping its shares along and across that circle.
pub(crate) fn total_stress(
    rise: &Rise,
    axial: Option<Axial>,
    p: [f64; 3],
    body_radius: f64,
) -> (f64, f64, f64, [f64; 3], [f64; 3]) {
    let c = rise.center;
    let x = angle_between(c, p) * body_radius / rise.radius;
    let (radial, hoop) = rise_stress(x);
    let along = dot(c, p);
    let out = [
        p[0] * along - c[0],
        p[1] * along - c[1],
        p[2] * along - c[2],
    ];
    let out = if dot(out, out) < 1e-20 {
        local_basis(p).0
    } else {
        normalize(out)
    };
    let round = normalize(cross(p, out));
    let Some(a) = axial else {
        return (radial, hoop, 0.0, out, round);
    };
    let carried = axis_at(c, a.axis, p);
    let (ca, sa) = (dot(carried, out), dot(carried, round));
    let (cos2, sin2) = (ca * ca - sa * sa, 2.0 * ca * sa);
    let s = a.stress / x.max(1.0).powi(2);
    (radial - s * cos2, hoop + s * cos2, -s * sin2, out, round)
}

/// The principal stresses of a surface stress given in a frame, the most tensile first,
/// and the strike across the most tensile direction in the frame's components: the
/// direction a fault or a dike runs.
fn principal(radial: f64, hoop: f64, shear: f64) -> (f64, f64, f64, f64) {
    let mean = 0.5 * (radial + hoop);
    let half = (0.25 * (radial - hoop).powi(2) + shear * shear).sqrt();
    // The most tensile lies at this angle from the radial direction, toward the hoop.
    let theta = 0.5 * (2.0 * shear).atan2(radial - hoop);
    (mean + half, mean - half, -theta.sin(), theta.cos())
}

/// Where faults start about a rise (F105, F114): as often as the ground is in tension,
/// over the distance from its centre in rise radii and the angle of its axis from the
/// direction outward there.
struct StartTable {
    cumulative: Vec<f64>,
    nx: usize,
    na: usize,
}

impl StartTable {
    fn new(axial: f64) -> Self {
        let (nx, na) = (256, 72);
        let dx = RIFT_REACH / nx as f64;
        let da = std::f64::consts::TAU / na as f64;
        let mut cumulative = Vec::with_capacity(nx * na + 1);
        cumulative.push(0.0);
        let mut sum = 0.0;
        for ix in 0..nx {
            let x = (ix as f64 + 0.5) * dx;
            let (radial, hoop) = rise_stress(x);
            let s = axial / x.max(1.0).powi(2);
            for ia in 0..na {
                let a = 2.0 * (ia as f64 + 0.5) * da;
                let (t1, _, _, _) =
                    principal(radial - s * a.cos(), hoop + s * a.cos(), -s * a.sin());
                sum += (t1 - RIFT_MIN_TENSION).max(0.0) * x * dx * da;
                cumulative.push(sum);
            }
        }
        Self { cumulative, nx, na }
    }

    /// A start's distance in rise radii and the angle of the axis from the direction
    /// outward there, rad, for three uniform draws.
    fn sample(&self, u: f64, v: f64, w: f64) -> (f64, f64) {
        let total = *self.cumulative.last().unwrap_or(&0.0);
        let cells = self.cumulative.len() - 1;
        let k = if total <= 0.0 {
            ((u * cells as f64) as usize).min(cells - 1)
        } else {
            let target = u.clamp(0.0, 1.0) * total;
            self.cumulative
                .partition_point(|&c| c <= target)
                .saturating_sub(1)
                .min(cells - 1)
        };
        let (ix, ia) = (k / self.na, k % self.na);
        (
            (ix as f64 + v) * RIFT_REACH / self.nx as f64,
            (ia as f64 + w) * std::f64::consts::TAU / self.na as f64,
        )
    }
}

/// A fault's trace through a rise's stress field (F105), from a point, in straight segments
/// (F111): it strikes across the most tensile direction of the total stress (F114), along
/// the radius where the hoop is the more tensile and round the rise where the radial is,
/// and keeps its own heading where the tension is about the same every way, as on the
/// crest. Each segment runs straight, SEGMENT_LENGTH_WIDTHS of the trough's `width` long;
/// where it ends the next turns toward the stress field by SEGMENT_STEER of the way as far
/// as the stress is anisotropic and by a few degrees more or less, and now and then steps
/// sideways to one side, en-echelon, over a short oblique link. The fault dies out where
/// the tension falls below RIFT_MIN_TENSION of the crest's. Returns the trace and the
/// distances along it where each segment after the first begins.
#[allow(clippy::too_many_arguments)]
pub(crate) fn stress_trace(
    rise: &Rise,
    start: [f64; 3],
    length: f64,
    spacing: f64,
    body_radius: f64,
    rng: &Rng,
    draw: u64,
    width: f64,
    first: Option<[f64; 3]>,
    min_tension: f64,
    axial: Option<Axial>,
) -> (Vec<[f64; 3]>, Vec<f64>) {
    let n = (length / spacing).ceil().max(2.0) as usize;
    let key = draw ^ 0x5E61_0000_0000;
    let (l_lo, l_hi) = SEGMENT_LENGTH_WIDTHS;
    let segment = |j: u64| width * l_lo * (l_hi / l_lo).powf(rng.uniform(key, 4 * j));
    // The total stress at a point, its principal values and the strike across the most
    // tensile direction there, and the frame's outward direction.
    let field = |p: [f64; 3]| {
        let (radial, hoop, shear, out, round) = total_stress(rise, axial, p, body_radius);
        let (s1, s2, so, sk) = principal(radial, hoop, shear);
        let strike = [
            out[0] * so + round[0] * sk,
            out[1] * so + round[1] * sk,
            out[2] * so + round[2] * sk,
        ];
        (s1, s2, strike, out, round)
    };
    let tangent = |v: [f64; 3], p: [f64; 3]| {
        let radial = dot(v, p);
        normalize([
            v[0] - radial * p[0],
            v[1] - radial * p[1],
            v[2] - radial * p[2],
        ])
    };
    // The first heading: across the most tensile direction, outward where that runs along
    // the radius, or the fault's own where the stress is about the same every way.
    let (s1, s2, strike, out, round) = field(start);
    let own = 2.0 * std::f64::consts::PI * rng.uniform(draw, 9);
    let flip = |v: [f64; 3]| [-v[0], -v[1], -v[2]];
    // A dike leaves its source on a heading of its own (F113).
    let mut heading = if let Some(h) = first {
        tangent(h, start)
    } else if s1 - s2 <= RIFT_LABYRINTH * s1.abs() {
        normalize([
            out[0] * own.cos() + round[0] * own.sin(),
            out[1] * own.cos() + round[1] * own.sin(),
            out[2] * own.cos() + round[2] * own.sin(),
        ])
    } else if dot(strike, out).abs() >= dot(strike, round).abs() {
        if dot(strike, out) >= 0.0 {
            strike
        } else {
            flip(strike)
        }
    } else if rng.uniform(draw, 10) < 0.5 {
        strike
    } else {
        flip(strike)
    };
    // A fault's en-echelon steps all go to one side.
    let side = if rng.uniform(draw, 11) < 0.5 {
        -1.0
    } else {
        1.0
    };
    let (s_lo, s_hi) = SEGMENT_STEP_WIDTHS;
    let link_angle = SEGMENT_LINK_ANGLE_DEG.to_radians();
    let mut points = Vec::with_capacity(n + 1);
    let mut breaks = Vec::new();
    let mut j = 0u64;
    let mut next_break = segment(0);
    // The steps left on a sideways link between segments, and its heading.
    let mut link = 0usize;
    let mut link_heading = heading;
    let mut p = start;
    for i in 0..=n {
        points.push(p);
        let (s1, s2, strike, _, _) = field(p);
        if i >= 2 && s1 < min_tension {
            break;
        }
        let s = i as f64 * spacing;
        if link == 0 && s >= next_break && i < n {
            j += 1;
            breaks.push(s);
            // The next segment turns toward the stress field, and a little either way.
            let along = if dot(strike, heading) < 0.0 {
                flip(strike)
            } else {
                strike
            };
            let isotropic = s1 - s2 <= RIFT_LABYRINTH * s1.abs();
            let anisotropy = ((s1 - s2) / (s1.abs() + s2.abs()).max(1e-9)).clamp(0.0, 1.0);
            let k = if isotropic {
                0.0
            } else {
                SEGMENT_STEER * anisotropy
            };
            let steered = tangent(
                [
                    heading[0] + k * (along[0] - heading[0]),
                    heading[1] + k * (along[1] - heading[1]),
                    heading[2] + k * (along[2] - heading[2]),
                ],
                p,
            );
            let turn = SEGMENT_STRIKE_SPREAD_DEG.to_radians() * rng.normal(key, 4 * j + 1);
            heading = crate::stages::io::rotate_about(p, steered, turn);
            // Now and then it steps sideways to the next segment, en-echelon.
            if rng.uniform(key, 4 * j + 3) < SEGMENT_STEP_PROBABILITY {
                let step = width * (s_lo + (s_hi - s_lo) * rng.uniform(key, 4 * j + 2));
                link = (step / (spacing * link_angle.sin())).ceil() as usize;
                link_heading = crate::stages::io::rotate_about(p, heading, side * link_angle);
            }
            next_break = s + link as f64 * spacing + segment(j);
        }
        let dir = if link > 0 {
            link -= 1;
            link_heading
        } else {
            heading
        };
        let step = spacing / body_radius;
        p = normalize([
            p[0] + step * dir[0],
            p[1] + step * dir[1],
            p[2] + step * dir[2],
        ]);
        // The headings stay in the tangent plane as the trace moves on: straight segments
        // along great circles.
        heading = tangent(heading, p);
        link_heading = tangent(link_heading, p);
    }
    (points, breaks)
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
    /// The trough's full depth, and the landslides off its walls in order along it
    /// (F106), with the farthest any reaches along the wall from its centre, m.
    depth: f64,
    landslides: Vec<Landslide>,
    slide_reach: f64,
    /// A debris lobe's volume over its peak thickness, its alcove's half width and its
    /// run-out.
    lobe_shape: f64,
    /// A graben's walls are rough fault scarps (F111, F112), and a fault's trace wanders.
    graben: bool,
    segmented: bool,
}

/// A debris lobe's volume over its peak thickness, its alcove's half width and its
/// run-out (F106): the integral of its thickness across, (1 - x²)^1.5 to 3π/8, times its
/// shape along, (1 - u²) times its widening and its rounded toe.
fn lobe_shape() -> f64 {
    let n = 256;
    let along: f64 = (0..n)
        .map(|k| {
            let u = (k as f64 + 0.5) / n as f64;
            (1.0 - u * u) * (1.0 + LANDSLIDE_SPREAD * u) * (1.0 - u.powi(4)).sqrt() / n as f64
        })
        .sum();
    3.0 * std::f64::consts::PI / 8.0 * along
}

/// One landslide off a canyon's wall (F106): an alcove cut back into the plateau and the
/// debris it shed, a lobe running out across the floor. Its sizes are in metres at the
/// trough's full depth, and grow with the trough.
#[derive(Clone, Copy, Debug)]
struct Landslide {
    /// Its centre along the trace, m, and the wall it fell from: the side of the normal.
    along: f64,
    side: f64,
    /// The alcove's half width along the wall and how far it cut back, m; the wall's
    /// height it fell, m; and its run-out over that height.
    half: f64,
    back: f64,
    fall: f64,
    reach: f64,
}

struct Frame {
    p: [f64; 3],
    t: [f64; 3],
    normal: [f64; 3],
    /// Its distance along the trace, m.
    along: f64,
    /// Depth and width factors towards the ends of the trace.
    taper: f64,
    narrow: f64,
    pulse: f64,
    /// The side its floor tilts down toward, -1 to 1 (F111).
    polarity: f64,
}

/// A trough's dimensions at a growth fraction, m.
struct TroughSize {
    width: f64,
    depth: f64,
    wall: f64,
    search: f64,
    half_band: f64,
    /// The cell of the grid that draws it, m.
    cell: f64,
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
            let sine = (std::f64::consts::PI * along_trace / trough.length)
                .sin()
                .max(0.0);
            if sine <= 0.0 {
                continue;
            }
            // A fault's depth follows its segments' throw (F111); a trough drawn as one
            // profile, a ring's arc or ridge whose depth may be negative, tapers as a sine.
            let taper = if trough.segments.is_empty() {
                sine.powf(TROUGH_TAPER_POWER)
            } else {
                trough_throw(trough, along_trace) / trough.depth.max(1.0)
            };
            if taper <= 0.0 {
                continue;
            }
            // A fault's trough keeps its width, narrowing only to its tips, where its
            // faults meet (F111); one drawn as one profile narrows all the way.
            let narrow = if trough.segments.is_empty() {
                sine.sqrt()
            } else {
                (along_trace.min(trough.length - along_trace)
                    / (TROUGH_TIP_WIDTHS * trough.width).max(1.0))
                .clamp(0.0, 1.0)
                .sqrt()
            };
            // Chasmata: a rift's width swells and narrows along the trace, and each wall is
            // scalloped; a graben keeps its width between straight walls.
            let pulse = if trough.kind == TroughKind::Graben {
                1.0
            } else {
                (1.0 + TROUGH_WIDTH_PULSE
                    * (noise.fbm(
                        [along_trace / TROUGH_PULSE_LENGTH, 0.37, 0.61],
                        1.0,
                        2,
                        2.0,
                        0.5,
                    ) / 0.6))
                    .clamp(0.5, 1.6)
            };
            frames.push(Frame {
                p,
                t,
                normal,
                along: along_trace,
                taper,
                narrow,
                pulse,
                polarity: trough_polarity(trough, along_trace),
            });
        }
        // The walls retreat in landslides where they stand high enough (F106): along each
        // wall, alcoves sized by the wall's height there, one after another.
        let mut landslides = Vec::new();
        if trough.depth >= LANDSLIDE_MIN_WALL && trough.length > 0.0 {
            let rng = Rng::new(trough.id ^ 0x1A2D_5D1D_0000);
            for (w, side) in [1.0f64, -1.0].into_iter().enumerate() {
                let mut s = 0.0;
                let mut k = 0u64;
                loop {
                    let draw = (w as u64) << 40 | k;
                    k += 1;
                    let fall = trough_throw(trough, s);
                    let (h_lo, h_hi) = LANDSLIDE_HALF_WIDTH;
                    let half = fall.max(LANDSLIDE_MIN_WALL)
                        * h_lo
                        * (h_hi / h_lo).powf(rng.uniform(draw, 0));
                    let (g_lo, g_hi) = LANDSLIDE_GAP;
                    s += half * (1.0 + 2.0 * (g_lo + (g_hi - g_lo) * rng.uniform(draw, 1)));
                    if s >= trough.length {
                        break;
                    }
                    let fall = trough_throw(trough, s);
                    if fall >= LANDSLIDE_MIN_WALL {
                        let (b_lo, b_hi) = LANDSLIDE_BACK;
                        let (r_lo, r_hi) = LANDSLIDE_REACH;
                        landslides.push(Landslide {
                            along: s,
                            side,
                            half,
                            back: half * (b_lo + (b_hi - b_lo) * rng.uniform(draw, 2)),
                            fall,
                            reach: r_lo + (r_hi - r_lo) * rng.uniform(draw, 3),
                        });
                    }
                    s += half;
                }
            }
            landslides.sort_by(|a, b| a.along.total_cmp(&b.along));
        }
        let slide_reach = landslides
            .iter()
            .map(|l| l.half * (1.0 + LANDSLIDE_SPREAD))
            .fold(0.0, f64::max);
        Self {
            frames,
            width: trough.width,
            flank: trough.flank,
            noise,
            depth: trough.depth,
            landslides,
            slide_reach,
            lobe_shape: lobe_shape(),
            graben: trough.kind == TroughKind::Graben,
            segmented: !trough.segments.is_empty(),
        }
    }

    /// How far a fault's trace wanders sideways off its segments at a distance along it, m
    /// (F112): self-affine, FAULT_MEANDER of its width over wavelengths from
    /// FAULT_MEANDER_WIDTHS widths down; none for a trough without segments.
    fn meander(&self, s: f64) -> f64 {
        if !self.segmented {
            return 0.0;
        }
        let wavelength = (FAULT_MEANDER_WIDTHS * self.width).max(1.0);
        FAULT_MEANDER
            * self.width
            * self.noise.fbm(
                [s / wavelength, 0.71, 0.13],
                1.0,
                3,
                2.0,
                2f64.powf(-FAULT_HURST),
            )
    }

    /// What the landslides near a point do to it (F106), at the trough's growth: how far
    /// the alcoves cut the wall back on its side, m, and how thick their debris lies on
    /// it, m. `along` is its distance along the trace and `across` off it, the trough
    /// `width` wide with walls `wall` wide there. A grid of cell `cell` draws only the
    /// landslides wider than its cell, leaving the finer to a finer grid: the loop's
    /// sampled alcoves stood in blocky steps where the finish drew the trough anew.
    fn landslides_at(
        &self,
        along: f64,
        across: f64,
        width: f64,
        wall: f64,
        scale: f64,
        cell: f64,
    ) -> (f64, f64) {
        if self.landslides.is_empty() || scale <= 0.0 {
            return (0.0, 0.0);
        }
        let first = self
            .landslides
            .partition_point(|l| l.along < along - self.slide_reach);
        let side = if across >= 0.0 { 1.0 } else { -1.0 };
        let foot = (0.5 * width - wall).max(0.0);
        let (mut cut, mut debris) = (0.0f64, 0.0f64);
        for l in &self.landslides[first..] {
            if l.along > along + self.slide_reach {
                break;
            }
            let ds = along - l.along;
            let (half, back, fall) = (l.half * scale, l.back * scale, l.fall * scale);
            if half < cell {
                continue;
            }
            if l.side == side && ds.abs() < half {
                cut = cut.max(back * (1.0 - (ds / half).powi(2)));
            }
            // The debris: the alcove's volume, a half cone, in a lobe from the foot of its
            // wall across the floor, as far as the floor goes; it widens toward its rounded
            // toe and thins to its edges, and drapes the wall's lower slope above the foot.
            // The volume spreads over the lobe and the drape together, so where the walls
            // close in on the floor, at a trough's tips or where the grid widens the walls,
            // the slope takes it rather than a lobe with no floor to run out on (a lobe a
            // metre long had made the debris tens of kilometres thick and emptied the
            // trough at each landslide).
            let run = (l.reach * fall).min(2.0 * foot);
            let into = foot - l.side * across;
            if into < -wall || into > run {
                continue;
            }
            let u = if run > 0.0 {
                (into / run).max(0.0)
            } else {
                0.0
            };
            let spread = half * (1.0 + LANDSLIDE_SPREAD * u) * (1.0 - u.powi(4)).max(0.0).sqrt();
            if ds.abs() >= spread {
                continue;
            }
            let volume = std::f64::consts::FRAC_PI_2 * half * back * fall / 3.0;
            let draped = 3.0 * std::f64::consts::PI / 8.0 * 0.5 * wall;
            let peak = volume / (half * (self.lobe_shape * run + draped));
            let drape = if into < 0.0 {
                smooth(into, -wall, 0.0)
            } else {
                1.0
            };
            debris +=
                peak * drape * (1.0 - u * u) * (1.0 - (ds / spread).powi(2)).max(0.0).powf(1.5);
        }
        (cut, debris)
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
                + if trough.flank != 0.0 { width } else { 0.0 }
                + if trough.segments.is_empty() {
                    0.0
                } else {
                    (1.5 * FAULT_MEANDER + GRABEN_WALL_ROUGHNESS) * trough.width
                }
                + if trough.depth >= LANDSLIDE_MIN_WALL {
                    LANDSLIDE_BACK.1 * LANDSLIDE_HALF_WIDTH.1 * depth
                } else {
                    0.0
                },
            half_band: 0.5 * trough.spacing + cell,
            cell,
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
        // A fault's trace wanders off its segments' great circles, rough at every scale
        // (F112).
        let across = dot(d, f.normal) * body_radius - self.meander(f.along + along);
        let scallop = if self.graben {
            // A graben's walls are rough fault traces, each standing in or out on its own,
            // so the graben pinches and swells (F112).
            1.0 + GRABEN_WALL_ROUGHNESS
                * self.noise.fbm(
                    [d[0] + across.signum() * 3.0, d[1], d[2]],
                    body_radius / (GRABEN_WALL_WIDTHS * self.width).max(1.0),
                    3,
                    2.0,
                    2f64.powf(-FAULT_HURST),
                )
        } else {
            1.0 + TROUGH_SCALLOP
                * self.noise.fbm(
                    [d[0] + across.signum() * 3.0, d[1], d[2]],
                    body_radius / TROUGH_SCALLOP_LENGTH.max(self.width),
                    2,
                    2.0,
                    0.5,
                )
        };
        let width = size.width * f.pulse * scallop * f.narrow;
        // The landslides cut the wall back on their side and lay their debris on the floor,
        // by the trough's own walls, not the grid's: a grid that widens the walls draws the
        // same landslides, only blurred.
        let (cut, debris) = self.landslides_at(
            f.along + along,
            across,
            width,
            TROUGH_WALL_FRACTION * size.width,
            size.depth / self.depth.max(1.0),
            size.cell,
        );
        let mut u = trough_profile(across, width + 2.0 * cut, size.depth, size.wall, blur);
        // A half-graben's floor tilts down toward its master fault, the far wall standing
        // lower (F111).
        if f.polarity != 0.0 && u > 0.0 {
            let toward = (f.polarity * across / (0.5 * width).max(1.0)).clamp(-1.0, 1.0);
            let tilt = f.polarity.abs() * RIFT_TILT;
            u *= 1.0 - tilt * 0.5 * (1.0 - toward);
        }
        if self.flank != 0.0 {
            u -= self.flank * size.depth * flank_profile(across, width, size.wall, blur);
        }
        let u = f.taper * u;
        let u = if u > 0.0 { (u - debris).max(0.0) } else { u };
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
pub fn shield_burial<'a>(
    shields: impl IntoIterator<Item = &'a Shield>,
    d: [f64; 3],
    body_radius: f64,
    since: f64,
    present: f64,
) -> f64 {
    shields
        .into_iter()
        .map(|s| {
            let grown = shield_progress(present, s) - shield_progress(since, s);
            if grown <= 0.0 || angle_between(s.center, d) * body_radius > s.extent() {
                0.0
            } else {
                shield_height(s, d, body_radius, false, 0.0) * grown
            }
        })
        .sum()
}

/// The share of a feature's relief that stands after lava this thick has buried it, the
/// feature this deep or high (F115): the lava fills it from below as it thickens, taking
/// its relief across its whole extent, and it is gone once the lava is as thick.
pub fn lava_keep(lava: f64, relief: f64) -> f64 {
    (1.0 - lava / relief.max(1.0)).clamp(0.0, 1.0)
}

/// Bury a crater's stamp under the lava the shields laid where it fell since it formed
/// (F115), as a trough is buried (F103): each cell's change kept by the share of the
/// crater's `relief` the lava there has not reached; where it has reached most of it, the
/// crater's surface no longer shows. `reach` is the stamp's radius about `center`.
#[allow(clippy::too_many_arguments)]
pub fn bury_under_shields(
    out: &mut StampOutput,
    layers: &LayerSet,
    grid: &CubeSphere,
    shields: &[Shield],
    body_radius: f64,
    center: [f64; 3],
    reach: f64,
    relief: f64,
    since: f64,
    present: f64,
) {
    let near: Vec<&Shield> = shields
        .iter()
        .filter(|s| {
            s.start_gyr < present
                && angle_between(s.center, center) * body_radius < s.extent() + reach
        })
        .collect();
    if near.is_empty() {
        return;
    }
    let bedrock = layers.bedrock.as_slice();
    let steeps = layers.steep.as_slice();
    let ejecta_ages = layers.ejecta_age.as_slice();
    for (j, cell) in out.cells.iter_mut().enumerate() {
        let idx = cell.idx as usize;
        let d = grid.center(grid.cell(idx));
        let lava = shield_burial(near.iter().copied(), d, body_radius, since, present);
        if lava <= 0.0 {
            continue;
        }
        let old = bedrock[idx] as f64;
        let change = cell.bedrock as f64 - old;
        let keep = lava_keep(lava, relief);
        cell.bedrock = (old + change * keep) as f32;
        if keep < 0.5 {
            cell.surface = false;
            if let Some(s) = &mut out.steep {
                s[j] = steeps[idx];
            }
            if let Some(e) = &mut out.ejecta_age {
                e[j] = ejecta_ages[idx];
            }
        }
    }
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
            // The lava fills the trough as it fills a crater (F115): what stands is the
            // share of its depth the lava has not reached, across its whole width.
            let keep = lava_keep(buried, trough.depth);
            let smeared = coarse.get(&idx).copied().unwrap_or(0.0) * keep;
            let depth = depth * keep;
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
    let deep = layers.deep.as_slice();
    let melt = layers.melt.as_slice();
    let cover = layers.cover.as_slice();
    // The summit's pond where the fine summit lies below the coarse, and the edifice's
    // own outline, each by the share of a cell its margin covers (F131).
    let pond = |e: [f64; 3]| {
        shield_height(shield, e, body_radius, false, 0.0)
            - shield_height(shield, e, body_radius, true, 0.0)
    };
    let edifice = |e: [f64; 3]| shield_height(shield, e, body_radius, false, 0.0);
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
        // The topmost eruption that holds the cell, else the topmost that touches it, and
        // the share of the cell the eruptions cover together (F131).
        let (mut whole, mut touched): (Option<usize>, Option<usize>) = (None, None);
        let mut bare = 1.0f64;
        for (j, c) in covers.iter().enumerate() {
            if let Some(&v) = c.get(&idx) {
                stack += FLOOD_SHEET_THICKNESS * v as f64;
                if v > 0.0 {
                    touched = Some(j);
                    bare *= 1.0 - (v as f64).clamp(0.0, 1.0);
                }
                if v >= 0.5 {
                    whole = Some(j);
                }
            }
        }
        let top = whole.or(touched);
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
        let delta = grid.delta();
        let (lava, share) = match top {
            Some(j) => (
                Some((flows[j].time_gyr as f32, flows[j].composition)),
                1.0 - bare,
            ),
            None => {
                let pond_share = if coarse - fine > 0.0 {
                    crate::fields::level_share(&pond, d, 1.0, delta)
                } else {
                    0.0
                };
                if pond_share > 0.0 {
                    (ponded, pond_share)
                } else if coarse > 0.0 && matches!(units[idx], UNIT_LAVA | UNIT_TESSERA) {
                    (
                        older.map(|c| (ages[idx], c)),
                        crate::fields::level_share(&edifice, d, 1.0, delta),
                    )
                } else {
                    (None, 0.0)
                }
            }
        };
        if let Some((age, comp)) = lava
            && share > 0.0
        {
            let share = share as f32;
            let held = share >= 0.5;
            write.surface = true;
            write.age = if held { age } else { ages[idx] };
            write.unit = if held { UNIT_LAVA } else { units[idx] };
            write.deep = (1.0 - share) * deep[idx];
            write.melt = (1.0 - share) * melt[idx];
            write.cover = share + (1.0 - share) * cover[idx];
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
            growth_gyr: 1.0,
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
            growth_gyr: 1.0,
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
    #[ignore = "diagnostic: lists Mars's troughs passing a point (EXOGEN_PT=lon,lat,km,face)"]
    fn troughs_near_a_point() {
        let spec = std::env::var("EXOGEN_PT").unwrap_or_else(|_| "-140.62,-3.57,120,512".into());
        let v: Vec<f64> = spec.split(',').map(|x| x.parse().unwrap()).collect();
        let (lon, lat) = (v[0].to_radians(), v[1].to_radians());
        let p = [lat.cos() * lon.cos(), lat.cos() * lon.sin(), lat.sin()];
        let (g, _) = crate::presets::preset("Mars").unwrap().to_graph();
        let r = g.f64("body.radius");
        let age = g.f64("body.age");
        let pp = PlumeProcess::new(
            &g,
            CubeSphere::new(v[3] as u32),
            &ThermalHistory::from_graph(&g),
        )
        .unwrap();
        for (kind, list) in [("canyon", &pp.canyons), ("graben", &pp.graben)] {
            for t in list.iter() {
                let near = t
                    .points
                    .iter()
                    .map(|q| angle_between(*q, p) * r / 1e3)
                    .fold(f64::INFINITY, f64::min);
                if near < v[2] {
                    println!(
                        "{kind} {:?} id {:x} time {:.3} ({:.2} Ga) start {:.3} width {:.1} km depth {:.0} m length {:.0} km nearest {:.0} km",
                        t.kind,
                        t.id,
                        t.end_gyr,
                        age - t.end_gyr,
                        t.start_gyr,
                        t.width / 1e3,
                        t.depth,
                        t.length / 1e3,
                        near
                    );
                }
            }
        }
    }

    #[test]
    #[ignore = "diagnostic: writes a canyon's relief with its landslides to a file"]
    fn write_canyon_patch() {
        let r = 3.39e6;
        let (east, north) = local_basis([1.0, 0.0, 0.0]);
        let points: Vec<[f64; 3]> = (0..=200)
            .map(|k| {
                offset(
                    [1.0, 0.0, 0.0],
                    east,
                    north,
                    (k as f64 * 10e3 - 1000e3) / r,
                    0.0,
                )
            })
            .collect();
        let canyon = Trough {
            points,
            spacing: 10e3,
            length: 2000e3,
            width: 150e3,
            depth: 6e3,
            flank: 0.0,
            start_gyr: 1.0,
            end_gyr: 2.0,
            id: 11,
            kind: TroughKind::Rift,
            segments: Vec::new(),
        };
        let frames = TroughFrames::new(&canyon);
        let n = 600usize;
        let half = 300e3;
        let mut out = String::from("P2\n600 600\n65535\n");
        let cell = 1e3;
        for y in 0..n {
            for x in 0..n {
                let e = (x as f64 / (n - 1) as f64 * 2.0 - 1.0) * half;
                let no = (1.0 - y as f64 / (n - 1) as f64 * 2.0) * half;
                let d = normalize([1.0, e / r, no / r]);
                let depth = trough_depth_at(
                    d,
                    r,
                    &canyon,
                    &frames,
                    (0, frames.len()),
                    1.0,
                    0.0,
                    cell,
                    0.0,
                );
                let v = ((1.0 - depth / 6.5e3).clamp(0.0, 1.0) * 65535.0) as u32;
                out.push_str(&format!("{v} "));
            }
            out.push('\n');
        }
        std::fs::write(std::env::temp_dir().join("exogen_canyon_patch.pgm"), out).unwrap();
    }

    #[test]
    fn a_canyon_s_walls_fall_in_landslides() {
        // F106: a deep canyon's walls carry alcoves along both sides, a shallow trough's
        // none, and a landslide's debris on the floor holds the volume its alcove lost.
        let r = 3.39e6;
        let (east, north) = local_basis([1.0, 0.0, 0.0]);
        let points: Vec<[f64; 3]> = (0..=200)
            .map(|k| offset([1.0, 0.0, 0.0], east, north, k as f64 * 10e3 / r, 0.3))
            .collect();
        let canyon = Trough {
            points: points.clone(),
            spacing: 10e3,
            length: 2000e3,
            width: 150e3,
            depth: 6e3,
            flank: 0.0,
            start_gyr: 1.0,
            end_gyr: 2.0,
            id: 11,
            kind: TroughKind::Rift,
            segments: Vec::new(),
        };
        let frames = TroughFrames::new(&canyon);
        assert!(frames.landslides.len() > 20, "{}", frames.landslides.len());
        assert!(frames.landslides.iter().any(|l| l.side > 0.0));
        assert!(frames.landslides.iter().any(|l| l.side < 0.0));
        let shallow = Trough {
            depth: 500.0,
            id: 12,
            ..canyon.clone()
        };
        assert!(TroughFrames::new(&shallow).landslides.is_empty());
        // One landslide's debris, summed over the floor about it, against its alcove.
        let l = frames
            .landslides
            .iter()
            .find(|l| l.along > 500e3 && l.side > 0.0)
            .unwrap();
        let (width, wall) = (150e3, TROUGH_WALL_FRACTION * 150e3);
        let only = TroughFrames {
            frames: Vec::new(),
            width,
            flank: 0.0,
            noise: GradientNoise::new(1),
            depth: canyon.depth,
            landslides: vec![*l],
            slide_reach: l.half * (1.0 + LANDSLIDE_SPREAD),
            lobe_shape: lobe_shape(),
            graben: false,
            segmented: false,
        };
        let step = 250.0;
        let mut volume = 0.0;
        let mut a = l.along - 2.0 * l.half;
        while a < l.along + 2.0 * l.half {
            let mut x = -0.5 * width;
            while x < 0.5 * width {
                volume += only.landslides_at(a, x, width, wall, 1.0, 0.0).1 * step * step;
                x += step;
            }
            a += step;
        }
        let alcove = std::f64::consts::FRAC_PI_2 * l.half * l.back * l.fall / 3.0;
        assert!((volume / alcove - 1.0).abs() < 0.1, "{volume} {alcove}");
    }

    #[test]
    fn a_canyon_keeps_its_depth_along_it_at_the_loop_s_cells() {
        // F106: the landslides follow the canyon's own walls whatever the grid draws them
        // at, and their debris spreads over its lobe and the wall's foot together, so with
        // walls a loop cell wide, which leave the floor no room, the canyon keeps its depth
        // at each landslide rather than filling to its rim there.
        let r = 3.39e6;
        let (east, north) = local_basis([1.0, 0.0, 0.0]);
        let points: Vec<[f64; 3]> = (0..=100)
            .map(|k| offset([1.0, 0.0, 0.0], east, north, k as f64 * 10e3 / r, 0.3))
            .collect();
        let canyon = Trough {
            points: points.clone(),
            spacing: 10e3,
            length: 1000e3,
            width: 40e3,
            depth: 2e3,
            flank: 0.0,
            start_gyr: 1.0,
            end_gyr: 2.0,
            id: 13,
            kind: TroughKind::Rift,
            segments: Vec::new(),
        };
        let frames = TroughFrames::new(&canyon);
        assert!(!frames.landslides.is_empty());
        // Walls a loop cell wide, which leave the floor no room, and the landslides drawn,
        // as a grid fine enough for them would.
        let (wall_floor, cell) = (21e3, 1e3);
        let depths: Vec<f64> = points[25..75]
            .iter()
            .map(|&p| {
                trough_depth_at(
                    p,
                    r,
                    &canyon,
                    &frames,
                    (0, frames.len()),
                    1.0,
                    wall_floor,
                    cell,
                    0.0,
                )
            })
            .collect();
        let deepest = depths.iter().cloned().fold(0.0, f64::max);
        let shallowest = depths.iter().cloned().fold(f64::INFINITY, f64::min);
        assert!(
            deepest > 0.0 && shallowest > 0.6 * deepest,
            "{shallowest} of {deepest}"
        );
    }

    #[test]
    fn a_rise_s_faults_follow_its_stress() {
        // F105: the crest is in tension every way alike, the flank in hoop tension with
        // radial compression, the foot in compression, and beyond it the load's hoop
        // tension falls away; a fault on the flank runs along the radius, and the starts
        // are drawn as the tension is.
        let (r0, h0) = rise_stress(0.0);
        assert!((r0 - 1.0).abs() < 1e-9 && (h0 - 1.0).abs() < 1e-9);
        let (r1, h1) = rise_stress(0.75);
        assert!(h1 > 0.0 && r1 < 0.0 && h1 > r1, "{r1} {h1}");
        let (r2, h2) = rise_stress(0.95);
        assert!(r2 < 0.0 && h2 < 0.05 && r2 < h2, "{r2} {h2}");
        let (r3, h3) = rise_stress(1.5);
        assert!(h3 > 0.0 && r3 < 0.0 && h3 < 0.1, "{r3} {h3}");
        let half = Tabulated::new(0.0, 2.0, 64, |x| if x > 1.0 { 1.0 } else { 0.0 });
        assert!((half.sample(0.0) - 1.0).abs() < 0.04 && (half.sample(0.5) - 1.5).abs() < 0.04);
        let rise = Rise {
            center: [1.0, 0.0, 0.0],
            radius: 2.85e6,
            height: 8e3,
            start_gyr: 0.3,
            end_gyr: 1.5,
            paleo_latitude: 0.0,
            axis: [0.0, 1.0, 0.0],
        };
        let r = 3.39e6;
        let (east, north) = local_basis(rise.center);
        let start = offset(rise.center, east, north, 0.55 * rise.radius / r, 1.0);
        let (trace, _) = stress_trace(
            &rise,
            start,
            800e3,
            10e3,
            r,
            &Rng::new(4),
            7,
            3e3,
            None,
            RIFT_MIN_TENSION,
            None,
        );
        assert!(trace.len() > 20);
        let (out, _) = (
            normalize([
                trace[1][0] - trace[0][0],
                trace[1][1] - trace[0][1],
                trace[1][2] - trace[0][2],
            ]),
            0,
        );
        // The outward direction in the tangent plane at the start.
        let c = rise.center;
        let radial = normalize([
            start[0] * dot(c, start) - c[0],
            start[1] * dot(c, start) - c[1],
            start[2] * dot(c, start) - c[2],
        ]);
        assert!(dot(out, radial).abs() > 0.9, "{}", dot(out, radial));
        let end = trace[trace.len() - 1];
        let x_end = angle_between(rise.center, end) * r / rise.radius;
        let x_start = angle_between(rise.center, start) * r / rise.radius;
        assert!((x_end - x_start).abs() > 0.15, "{x_start} {x_end}");
    }

    #[test]
    fn a_fault_grows_in_straight_segments_that_link() {
        // F111: a fault's trace runs straight within each segment and turns or steps where
        // segments meet; a graben keeps its width to its tips, and its throw peaks within
        // each segment, dips where they link and falls to nothing at its ends.
        let rise = Rise {
            center: [1.0, 0.0, 0.0],
            radius: 2.85e6,
            height: 8e3,
            start_gyr: 0.3,
            end_gyr: 1.5,
            paleo_latitude: 0.0,
            axis: [0.0, 1.0, 0.0],
        };
        let r = 3.39e6;
        let (east, north) = local_basis(rise.center);
        let start = offset(rise.center, east, north, 0.55 * rise.radius / r, 1.0);
        let (width, spacing) = (3e3, 1e3);
        let (trace, breaks) = stress_trace(
            &rise,
            start,
            600e3,
            spacing,
            r,
            &Rng::new(5),
            11,
            width,
            None,
            RIFT_MIN_TENSION,
            None,
        );
        assert!(breaks.len() >= 3, "{} segments", breaks.len() + 1);
        // Within a segment the trace keeps to a great circle: successive steps turn by far
        // less than a degree, except at a break and on the link after it.
        let link = (SEGMENT_STEP_WIDTHS.1 * width
            / (spacing * SEGMENT_LINK_ANGLE_DEG.to_radians().sin()))
        .ceil();
        let chord = |k: usize| {
            let (a, b) = (trace[k], trace[k + 1]);
            let m = normalize([a[0] + b[0], a[1] + b[1], a[2] + b[2]]);
            let v = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
            let radial = dot(v, m);
            normalize([
                v[0] - radial * m[0],
                v[1] - radial * m[1],
                v[2] - radial * m[2],
            ])
        };
        let mut bent = 0;
        for k in 1..trace.len() - 1 {
            let s = k as f64 * spacing;
            let near = breaks
                .iter()
                .any(|&b| s >= b - spacing && s <= b + (link + 1.0) * spacing);
            let turn = dot(chord(k - 1), chord(k))
                .clamp(-1.0, 1.0)
                .acos()
                .to_degrees();
            if !near && turn > 0.2 {
                bent += 1;
            }
        }
        assert_eq!(bent, 0, "a segment bends");
        let length = (trace.len() - 1) as f64 * spacing;
        let mut bounds = vec![0.0];
        bounds.extend(breaks.iter().copied());
        bounds.push(length);
        let segments: Vec<Segment> = bounds
            .windows(2)
            .enumerate()
            .map(|(j, w)| Segment {
                start: w[0],
                end: w[1],
                throw: 200.0 + 10.0 * (j % 5) as f64,
                polarity: 0.0,
            })
            .collect();
        let trough = Trough {
            points: trace.clone(),
            spacing,
            length,
            width,
            depth: segments.iter().map(|g| g.throw).fold(0.0, f64::max),
            flank: 0.0,
            start_gyr: 1.0,
            end_gyr: 1.0,
            id: 5,
            kind: TroughKind::Graben,
            segments: segments.clone(),
        };
        // Away from its tips, a graben keeps each segment's throw at its middle and dips by
        // SEGMENT_LINK_DEFICIT at a link, between the two segments' throws; at its tips it
        // falls to nothing (F112).
        let tip = trough.depth / (2.0 * SEGMENT_THROW_PER_LENGTH);
        let inside = |s: f64| s > tip && s < length - tip;
        for g in &segments {
            let s = 0.5 * (g.start + g.end);
            if inside(s) {
                let mid = trough_throw(&trough, s);
                assert!((mid - g.throw).abs() < 1e-6, "{mid} vs {}", g.throw);
            }
        }
        for (k, &b) in breaks.iter().enumerate() {
            if inside(b) {
                let dip = trough_throw(&trough, b);
                let expect = (1.0 - SEGMENT_LINK_DEFICIT)
                    * 0.5
                    * (segments[k].throw + segments[k + 1].throw);
                assert!((dip - expect).abs() < 1.0, "{dip} vs {expect}");
            }
        }
        assert!(trough_throw(&trough, 0.0) < 1e-6 && trough_throw(&trough, length) < 1e-6);
        // The profile runs on unbroken across a link.
        for &b in &breaks {
            let (a, c) = (
                trough_throw(&trough, b - 1.0),
                trough_throw(&trough, b + 1.0),
            );
            assert!((a - c).abs() < 5.0, "{a} {c}");
        }
        let frames = TroughFrames::new(&trough);
        for f in &frames.frames {
            if f.along > 2.0 * width && f.along < length - 2.0 * width {
                assert!(
                    (f.narrow - 1.0).abs() < 1e-9 && f.pulse == 1.0,
                    "{} {}",
                    f.narrow,
                    f.pulse
                );
            }
        }
        // A rift's half-graben: its floor stands deeper toward its master fault, by
        // RIFT_TILT at the far wall, and turns over at a link where the polarity flips.
        let mut rift = trough.clone();
        rift.kind = TroughKind::Rift;
        for (j, g) in rift.segments.iter_mut().enumerate() {
            g.polarity = if j % 2 == 0 { 1.0 } else { -1.0 };
        }
        let mid = 0.5 * (rift.segments[0].start + rift.segments[0].end);
        assert_eq!(trough_polarity(&rift, mid), 1.0);
        assert!(trough_polarity(&rift, rift.segments[0].end).abs() < 1e-9);
        let frames = TroughFrames::new(&rift);
        let f = frames
            .frames
            .iter()
            .min_by(|a, b| (a.along - mid).abs().total_cmp(&(b.along - mid).abs()))
            .unwrap();
        let size = TroughFrames::size(&rift, 1.0, 0.0, 100.0);
        let side = |s: f64| {
            let q = normalize([
                f.p[0] + s * 0.3 * width / r * f.normal[0],
                f.p[1] + s * 0.3 * width / r * f.normal[1],
                f.p[2] + s * 0.3 * width / r * f.normal[2],
            ]);
            frames.depth_from(q, r, &size, 0.0, f).unwrap().1
        };
        let (master, far) = (side(1.0), side(-1.0));
        assert!(master > far && far > 0.0, "{master} {far}");
    }

    #[test]
    fn graben_are_as_wide_as_their_layer_and_long_faults_open_rifts() {
        // F111: on Mars's rise a graben's width is its faulted layer's, whatever its
        // length, and its throw stays within the layer; the faults whose throw over their
        // length would exceed the layer are rifts, as wide as their throw, the long ones.
        let (g, _) = crate::presets::preset("Mars").unwrap().to_graph();
        let p =
            PlumeProcess::new(&g, CubeSphere::new(64), &ThermalHistory::from_graph(&g)).unwrap();
        let all: Vec<&Trough> = p.canyons.iter().chain(&p.graben).collect();
        let graben: Vec<&&Trough> = all
            .iter()
            .filter(|t| t.kind == TroughKind::Graben)
            .collect();
        let rifts: Vec<&&Trough> = all.iter().filter(|t| t.kind == TroughKind::Rift).collect();
        assert!(
            graben.len() > 1000 && !rifts.is_empty(),
            "{} graben, {} rifts",
            graben.len(),
            rifts.len()
        );
        for t in &graben {
            assert!(
                t.width >= GRABEN_WIDTH_RANGE.0 && t.width <= GRABEN_WIDTH_RANGE.1,
                "{}",
                t.width
            );
            assert!(t.depth <= FAULT_LAYER_RANGE.1 + 1e-6, "{}", t.depth);
            assert!(t.depth <= (GRABEN_DEPTH_PER_WIDTH * t.width).max(RIFT_MIN_DEPTH) + 1e-6);
            assert!(!t.segments.is_empty());
        }
        for t in &rifts {
            assert!(
                t.width >= RIFT_WIDTH_PER_DEPTH * FAULT_LAYER_RANGE.0 - 1.0,
                "{}",
                t.width
            );
        }
        let median = |v: &mut Vec<f64>| {
            v.sort_by(f64::total_cmp);
            v[v.len() / 2]
        };
        let rift_length = median(&mut rifts.iter().map(|t| t.length).collect());
        let graben_length = median(&mut graben.iter().map(|t| t.length).collect());
        assert!(
            rift_length > 3.0 * graben_length,
            "{rift_length} vs {graben_length}"
        );
    }

    /// Diagnostic: a preset's rises and shields with their times (EXOGEN_PRESET, Venus by
    /// default); run with `--ignored --nocapture`.
    #[test]
    #[ignore]
    fn print_rises_and_shields() {
        let name = std::env::var("EXOGEN_PRESET").unwrap_or_else(|_| "Venus".into());
        let (g, _) = crate::presets::preset(&name).unwrap().to_graph();
        let thermal = ThermalHistory::from_graph(&g);
        let p = PlumeProcess::new(&g, CubeSphere::new(32), &thermal).unwrap();
        for (i, r) in p.rises.iter().enumerate() {
            println!(
                "rise {i}: {:.3} to {:.3} Gyr, radius {:.0} km",
                r.start_gyr,
                r.end_gyr,
                r.radius / 1e3
            );
        }
        for s in &p.shields {
            println!(
                "shield {} on rise {}: start {:.3}, grows over {:.3} Gyr, radius {:.0} km",
                s.id,
                rise_of(&p.rises, s.center),
                s.start_gyr,
                s.growth_gyr,
                s.radius / 1e3
            );
        }
    }

    #[test]
    fn a_rise_grows_after_its_lid_last_formed_as_fast_as_its_mantle_melts() {
        // F120, F121: Mars's rise keeps Tharsis's window; Venus's grows after its overturn,
        // the largest shield on each starting before the present and growing faster, its
        // mantle melting faster.
        let run = |name: &str| {
            let (g, _) = crate::presets::preset(name).unwrap().to_graph();
            let thermal = ThermalHistory::from_graph(&g);
            let p = PlumeProcess::new(&g, CubeSphere::new(32), &thermal).unwrap();
            (g, thermal, p)
        };
        let (_, _, mars) = run("Mars");
        let rise = &mars.rises[0];
        assert!(
            (rise.start_gyr - 0.3).abs() < 0.02 && (rise.end_gyr - 1.5).abs() < 0.02,
            "Mars's rise {} to {}",
            rise.start_gyr,
            rise.end_gyr
        );
        let (g, thermal, venus) = run("Venus");
        let age = g.f64("body.age");
        let overturn = age - g.f64("tectonics.resurfacing_event_age");
        assert!(!venus.rises.is_empty());
        for (ri, r) in venus.rises.iter().enumerate() {
            assert!(
                r.start_gyr > overturn && r.start_gyr < age,
                "overturn {overturn}, rise {} to {}",
                r.start_gyr,
                r.end_gyr
            );
            let largest = venus
                .shields
                .iter()
                .filter(|s| rise_of(&venus.rises, s.center) == ri)
                .max_by(|a, b| a.radius.total_cmp(&b.radius))
                .unwrap();
            assert!(
                largest.start_gyr > overturn && largest.start_gyr < age,
                "{}",
                largest.start_gyr
            );
        }
        for s in &venus.shields {
            assert!(s.start_gyr > overturn, "{}", s.start_gyr);
            assert!(
                s.growth_gyr < mars.shields[0].growth_gyr,
                "{} vs {}",
                s.growth_gyr,
                mars.shields[0].growth_gyr
            );
        }
        // The melt clock runs both ways, and on past the present.
        let r = g.f64("body.radius");
        let clock = thermal.melt_clock(4.0 * std::f64::consts::PI * r * r);
        for t in [0.5, 2.0, 4.2, 5.0] {
            assert!((clock.time_of(clock.at(t)) - t).abs() < 1e-9, "{t}");
        }
    }

    #[test]
    fn a_rise_grows_along_its_axis() {
        // F114: tension across the axis and compression along it. On the crest the tension
        // is still nearly the same every way; beyond the foot the ground along the axis is
        // in tension farther out than the load alone puts it, and the ground across it is
        // not; on Mars the faults gather along the axis and the shields line up on it.
        let rise = Rise {
            center: [1.0, 0.0, 0.0],
            radius: 2.85e6,
            height: 8e3,
            start_gyr: 0.3,
            end_gyr: 1.5,
            paleo_latitude: 0.0,
            axis: [0.0, 1.0, 0.0],
        };
        let r = 3.39e6;
        let (east, north) = local_basis(rise.center);
        let axial = Axial {
            axis: east,
            stress: RISE_AXIAL_STRESS,
        };
        let at = |x: f64, phi: f64| offset(rise.center, east, north, x * rise.radius / r, phi);
        let tension = |p: [f64; 3]| {
            let (radial, hoop, shear, _, _) = total_stress(&rise, Some(axial), p, r);
            principal(radial, hoop, shear)
        };
        let (c1, c2, _, _) = tension(at(0.01, 0.3));
        assert!((c1 - 1.0 - RISE_AXIAL_STRESS).abs() < 0.01 && c1 - c2 <= RIFT_LABYRINTH * c1);
        let (along, _, _, _) = tension(at(1.4, 0.0));
        let (across, _, _, _) = tension(at(1.4, std::f64::consts::FRAC_PI_2));
        assert!(
            along > 1.5 * RIFT_MIN_TENSION && across < RIFT_MIN_TENSION,
            "{along} {across}"
        );
        // Along the axis the faults reach past where the load's tension alone would end
        // them, and die out short of twice the rise's radius.
        let alone = |x: f64| {
            let (radial, hoop) = rise_stress(x);
            radial.max(hoop)
        };
        let (beyond, _, _, _) = tension(at(1.6, 0.0));
        let (farther, _, _, _) = tension(at(2.0, 0.0));
        assert!(
            alone(1.6) < RIFT_MIN_TENSION
                && beyond > RIFT_MIN_TENSION
                && farther < RIFT_MIN_TENSION,
            "{} {beyond} {farther}",
            alone(1.6)
        );
        // On Mars the faults beyond the foot gather along the axis, and the shields lie on it.
        let (g, _) = crate::presets::preset("Mars").unwrap().to_graph();
        let p =
            PlumeProcess::new(&g, CubeSphere::new(64), &ThermalHistory::from_graph(&g)).unwrap();
        let rise = &p.rises[0];
        let (re, rn) = local_basis(rise.center);
        let phi_axis = dot(rise.axis, rn).atan2(dot(rise.axis, re));
        let off_axis = |d: [f64; 3]| {
            let phi = azimuth_from(rise.center, re, rn, d) - phi_axis;
            let wrapped = phi.rem_euclid(std::f64::consts::PI);
            wrapped.min(std::f64::consts::PI - wrapped)
        };
        for s in &p.shields {
            assert!(off_axis(s.center) < 3.5 * SHIELD_AXIS_SPREAD_DEG.to_radians());
        }
        let (mut near, mut far) = (0, 0);
        for t in p.graben.iter().chain(&p.canyons) {
            let m = t.midpoint();
            if angle_between(rise.center, m) * g.f64("body.radius") < 1.2 * rise.radius {
                continue;
            }
            let o = off_axis(m);
            if o < 30f64.to_radians() {
                near += 1;
            } else if o > 60f64.to_radians() {
                far += 1;
            }
        }
        assert!(near > 3 * far, "{near} along the axis, {far} across it");
    }

    #[test]
    fn a_rise_s_dike_swarms_radiate_long_graben() {
        // F113: the plume and the volcanoes on Mars's rise send dikes out in swarms, each
        // leaving its source and running on away from it, a graben above each; some run
        // thousands of kilometres, as Tharsis's radial graben systems do.
        let (g, _) = crate::presets::preset("Mars").unwrap().to_graph();
        let p =
            PlumeProcess::new(&g, CubeSphere::new(64), &ThermalHistory::from_graph(&g)).unwrap();
        let r = g.f64("body.radius");
        let dikes: Vec<&Trough> = p
            .graben
            .iter()
            .chain(&p.canyons)
            .filter(|t| t.id & 0xFFFF_0000_0000_0000 == 0x0D1C_0000_0000_0000)
            .collect();
        assert!(dikes.len() > 50, "{} dike graben", dikes.len());
        let longest = dikes.iter().map(|t| t.length).fold(0.0, f64::max);
        assert!(longest > 2000e3, "{longest}");
        let mut sources: Vec<[f64; 3]> = p.rises.iter().map(|x| x.center).collect();
        sources.extend(p.shields.iter().map(|s| s.center));
        for t in &dikes {
            assert_eq!(t.kind, TroughKind::Graben);
            assert!(t.width <= GRABEN_WIDTH_RANGE.1);
            let k = t.points.len().min(20) - 1;
            let gain = sources
                .iter()
                .map(|&c| angle_between(c, t.points[k]) - angle_between(c, t.points[0]))
                .fold(f64::NEG_INFINITY, f64::max);
            assert!(
                gain * r > 0.5 * k as f64 * t.spacing,
                "a dike runs back toward its source"
            );
        }
    }

    #[test]
    fn a_rise_turns_the_body_toward_its_equator() {
        // F104: the fluid Love number is 3/2 for a uniform body and about 0.93 for the
        // Earth's moment of inertia; Mars's lithosphere holds a few percent of a degree-two
        // load; a rise near the pole is carried near the equator against a small frozen
        // bulge and hardly at all against a large one; with no bulge one rise goes to the
        // equator, and the eigen solver finds a known matrix's axes.
        assert!((fluid_love_number(0.4) - 1.5).abs() < 1e-12);
        assert!((fluid_love_number(0.3307) - 0.93).abs() < 0.01);
        let support = degree_two_support(33e3, 3.73, 3.39e6);
        assert!(support > 0.005 && support < 0.05, "{support}");
        let (east, north) = local_basis([0.0, 0.0, 1.0]);
        let rise = Rise {
            center: offset([0.0, 0.0, 1.0], east, north, 7f64.to_radians(), 0.4),
            radius: 2.85e6,
            height: 8e3,
            start_gyr: 0.3,
            end_gyr: 1.5,
            paleo_latitude: 83f64.to_radians(),
            axis: [0.0, 1.0, 0.0],
        };
        let r = 3.39e6;
        let latitude = |bulge: f64| {
            let (axis, angle) = polar_wander(std::slice::from_ref(&rise), r, bulge);
            normalize(rotate(rise.center, axis, angle))[2]
                .asin()
                .to_degrees()
        };
        let mars = frozen_bulge(r, 3.73, 7.09e-5, 0.33, 33e3);
        assert!(latitude(mars).abs() < 15.0, "{}", latitude(mars));
        assert!(latitude(1e3 * mars) > 75.0, "{}", latitude(1e3 * mars));
        assert!(latitude(0.0).abs() < 1e-6, "{}", latitude(0.0));
        let (c, s) = (0.6f64, 0.8f64);
        let m = [
            [2.0 * c * c + s * s, (2.0 - 1.0) * c * s, 0.0],
            [(2.0 - 1.0) * c * s, 2.0 * s * s + c * c, 0.0],
            [0.0, 0.0, -3.0],
        ];
        let (values, vectors) = crate::geom::symmetric_eigen(m);
        assert!((values[0] - 2.0).abs() < 1e-12 && (values[1] - 1.0).abs() < 1e-12);
        assert!((dot(vectors[0], [c, s, 0.0]).abs() - 1.0).abs() < 1e-12);
        assert!((values[2] + 3.0).abs() < 1e-12);
    }

    #[test]
    fn a_rise_s_shields_gather_its_dikes_as_its_stretching_wanes() {
        // F119: none while the dome's growth quickens, then the share by which its rate has
        // fallen from its peak, all once it has stopped; each shield starts when that share
        // reaches what the larger ones hold, the largest as the growth turns.
        let rise = Rise {
            center: [1.0, 0.0, 0.0],
            radius: 2.5e6,
            height: 7e3,
            start_gyr: 0.3,
            end_gyr: 1.5,
            paleo_latitude: 0.0,
            axis: [0.0, 0.0, 1.0],
        };
        let at = |t: f64| rise_focus(&rise, t);
        assert_eq!(at(0.2), 0.0);
        assert!(at(0.9) < 1e-12, "{}", at(0.9));
        assert!((at(1.2) - 0.25).abs() < 1e-12, "{}", at(1.2));
        assert_eq!(at(1.6), 1.0);
        assert!((shield_start(&rise, 0.0) - 0.9).abs() < 1e-12);
        for share in [0.1, 0.57, 0.9] {
            let t = shield_start(&rise, share);
            assert!((at(t) - share).abs() < 1e-12, "{share}: {t} {}", at(t));
        }
        assert!((shield_start(&rise, 1.0) - 1.5).abs() < 1e-12);
        // A rise raised at once gathers all after it.
        let sudden = Rise {
            start_gyr: 4.0,
            end_gyr: 4.0,
            ..rise
        };
        assert_eq!(rise_focus(&sudden, 3.9), 0.0);
        assert_eq!(rise_focus(&sudden, 4.1), 1.0);
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
            growth_gyr: 1.0,
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
    fn a_shield_buries_the_craters_it_grows_over() {
        // F115: a crater on the ground a shield grows over fills with its lava rather than
        // riding up on it. Runs of Mars's plume, one on flat ground and others with bowls
        // 1 km deep stamped before anything grew, one at mid-flank on its largest shield and
        // one beyond that shield's foot: at the present, the difference between the runs is
        // what is left of each bowl. The flank's bowl is tried at eight places round the
        // shield, each in its own run, since the grabens and canyons crossing the flank
        // decide what is left at any one place.
        let (g, _) = crate::presets::preset("Mars").unwrap().to_graph();
        let grid = CubeSphere::new(128);
        let r = g.f64("body.radius");
        let age = g.f64("body.age");
        let thermal = ThermalHistory::from_graph(&g);
        let mut plain = PlumeProcess::new(&g, grid, &thermal).unwrap();
        let shield = plain
            .shields
            .iter()
            .max_by(|a, b| a.height.total_cmp(&b.height))
            .unwrap()
            .clone();
        let (east, north) = local_basis(shield.center);
        // Beyond its foot, where no other shield reaches either.
        let beyond = (0..16)
            .map(|k| {
                let phi = k as f64 * std::f64::consts::PI / 8.0;
                offset(
                    shield.center,
                    east,
                    north,
                    (shield.extent() + 200e3) / r,
                    phi,
                )
            })
            .find(|q| {
                plain
                    .shields
                    .iter()
                    .all(|s| angle_between(s.center, *q) * r > s.extent() + 100e3)
            })
            .unwrap();
        let bowl = |c: [f64; 3], q: [f64; 3]| {
            let x = angle_between(c, q) * r / 60e3;
            if x < 1.0 { -1e3 * (1.0 - x * x) } else { 0.0 }
        };
        let steps = 90;
        let run = |process: &mut PlumeProcess, layers: &mut LayerSet| {
            for k in 0..steps {
                let (t0, t1) = (
                    age * k as f64 / steps as f64,
                    age * (k + 1) as f64 / steps as f64,
                );
                process.step(layers, t0, t1);
            }
        };
        let mut flat = LayerSet::new(grid);
        run(&mut plain, &mut flat);
        let depth = |layers: &LayerSet, c: [f64; 3]| {
            let i = grid.index(grid.locate(c).0);
            (flat.bedrock.as_slice()[i] - layers.bedrock.as_slice()[i]) as f64
        };
        let mut kept = Vec::new();
        for k in 0..8 {
            let flank = offset(
                shield.center,
                east,
                north,
                0.5 * shield.radius / r,
                0.7 + k as f64 * std::f64::consts::FRAC_PI_4,
            );
            let mut pitted = LayerSet::new(grid);
            for (i, b) in pitted.bedrock.as_mut_slice().iter_mut().enumerate() {
                let q = grid.center(grid.cell(i));
                *b = (bowl(flank, q) + bowl(beyond, q)) as f32;
            }
            // The bowls' depths before anything grew: flat ground less the bowl.
            let i_on = grid.index(grid.locate(flank).0);
            let i_off = grid.index(grid.locate(beyond).0);
            let (on, off) = (
                -pitted.bedrock.as_slice()[i_on] as f64,
                -pitted.bedrock.as_slice()[i_off] as f64,
            );
            let mut cratered = PlumeProcess::new(&g, grid, &thermal).unwrap();
            run(&mut cratered, &mut pitted);
            let (on_now, off_now) = (depth(&pitted, flank), depth(&pitted, beyond));
            assert!(on > 500.0, "on the flank {on} m");
            assert!(
                (off_now - off).abs() < 0.02 * off,
                "beyond the foot {off} m, now {off_now} m"
            );
            kept.push(on_now / on);
        }
        // Riding up on the shield a bowl would keep its depth; buried, it keeps a tenth or
        // so where no graben crosses it.
        kept.sort_by(|a, b| a.total_cmp(b));
        let median = 0.5 * (kept[3] + kept[4]);
        assert!(median < 0.25, "kept {kept:?}");
    }

    #[test]
    fn faults_and_swarms_come_with_what_drives_them() {
        // F116: on Mars the rise's faults form while its dome rises, and the plume head's
        // dikes with them; each volcano's dikes come while it grows.
        let (g, _) = crate::presets::preset("Mars").unwrap().to_graph();
        let p =
            PlumeProcess::new(&g, CubeSphere::new(64), &ThermalHistory::from_graph(&g)).unwrap();
        let rise = &p.rises[0];
        let age = g.f64("body.age");
        let (mut faults, mut plume_dikes, mut volcano_dikes) = (0, 0, 0);
        for t in p.graben.iter().chain(&p.canyons) {
            let time = 0.5 * (t.start_gyr + t.end_gyr);
            let rising = time >= rise.start_gyr - 1e-9 && time <= rise.end_gyr + 1e-9;
            if t.id & 0xFFFF_0000_0000_0000 == 0x0D1C_0000_0000_0000 {
                let key = (t.id & 0x0000_FFFF_FFFF_FFFF) >> 8;
                if key & 0xFFFF_0000 == 0x0D1D_0000 {
                    let s = p
                        .shields
                        .iter()
                        .find(|s| s.id == key - 0x0D1D_0000)
                        .unwrap();
                    assert!(time >= s.start_gyr - 1e-9 && time <= age + 1e-9, "{time}");
                    volcano_dikes += 1;
                } else {
                    assert!(rising, "a plume dike at {time}");
                    plume_dikes += 1;
                }
            } else {
                assert!(rising, "a fault at {time}");
                faults += 1;
            }
        }
        assert!(
            faults > 100 && plume_dikes > 20 && volcano_dikes > 5,
            "{faults} faults, {plume_dikes} plume dikes, {volcano_dikes} volcano dikes"
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
            growth_gyr: 1.0,
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
