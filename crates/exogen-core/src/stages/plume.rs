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
use crate::geom::{
    angle_between, azimuth_from, cross, dot, local_basis, normalize, offset, rotate,
};
use crate::grid::CubeSphere;
use crate::grid::projection::PaddedField;
use crate::noise::GradientNoise;
use crate::params::ParamGraph;
use crate::rng::Rng;
use crate::stages::io::{Lobe, PathRules};
use crate::units::{UNIT_LAVA, UNIT_TESSERA};

pub const VERSION: u64 = 8;

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
                paleo_latitude: center[2].clamp(-1.0, 1.0).asin(),
            });
        }
        // The body turns under its rises as they grow, carrying them toward its equator
        // (F104); everything on them is placed where they stand now.
        {
            let mid = 0.5 * (RISE_START_GYR + RISE_END_GYR);
            let lithosphere_now = graph.f64("interior.lithosphere_thickness").max(1.0);
            let te = graph.f64("interior.elastic_thickness")
                * graph
                    .curve("interior.lithosphere_history")
                    .map_or(lithosphere_now, |c| c.eval(mid))
                / lithosphere_now;
            let bulge = frozen_bulge(
                body_radius,
                graph.f64("body.gravity"),
                graph.f64("rotation.angular_velocity"),
                graph.f64("body.moment_of_inertia_factor"),
                te,
            );
            let (axis, angle) = polar_wander(&rises, body_radius, bulge);
            for rise in &mut rises {
                rise.center = normalize(rotate(rise.center, axis, angle));
            }
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
        // Where faults start: as often as the ground is in tension, over the area round a
        // rise (F105).
        let start_cdf = Tabulated::new(0.0, RIFT_REACH, 512, |x| {
            let (sr, st) = rise_stress(x);
            (sr.max(st) - RIFT_MIN_TENSION).max(0.0) * x
        });
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
            // Faults form as the rise's load grows (F105): on a stagnant lid from the rise's
            // start at the pace of the body's volcanism, which builds the load, to the
            // present; on an episodic lid after the resurfacing event.
            let (t_lo, t_hi) = if episodic {
                (event_time, age)
            } else {
                (rise.start_gyr.min(age), age)
            };
            let time_cdf = Tabulated::new(t_lo, t_hi, 256, |t| {
                if episodic {
                    1.0
                } else {
                    graph
                        .curve("volcanism.history")
                        .map_or(1.0, |c| c.eval(t).max(0.0))
                }
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
                        // A fault starts where the ground is in tension, as often as the
                        // tension is strong, and runs as the stress field leads it (F105).
                        let x = start_cdf.sample(rng.uniform(draw, 20));
                        let phi = 2.0 * std::f64::consts::PI * rng.uniform(draw, 21);
                        let start =
                            offset(rise.center, east, north, x * rise.radius / body_radius, phi);
                        let trace =
                            stress_trace(rise, start, length, spacing, body_radius, &rng, draw);
                        length = (trace.len().max(2) - 1) as f64 * spacing;
                        trace
                    }
                };
                // A fault that died out early slipped as its own length allows.
                let depth =
                    (RIFT_DEPTH_PER_LENGTH * length * scatter).clamp(RIFT_MIN_DEPTH, depth_cap);
                let width = (RIFT_WIDTH_PER_DEPTH * depth).max(RIFT_MIN_WIDTH);
                let time = time_cdf.sample(rng.uniform(draw, 8));
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

/// A fault's trace through a rise's stress field (F105), from a point: it strikes across
/// the most tensile direction, along the radius where the hoop is the more tensile and
/// round the rise where the radial is, turning toward that by RIFT_STEER a step as far as
/// the stress is anisotropic and keeping its own heading where the tension is about the
/// same every way, as on the crest; it wanders a little, and dies out where the tension
/// falls below RIFT_MIN_TENSION of the crest's.
pub(crate) fn stress_trace(
    rise: &Rise,
    start: [f64; 3],
    length: f64,
    spacing: f64,
    body_radius: f64,
    rng: &Rng,
    draw: u64,
) -> Vec<[f64; 3]> {
    let n = (length / spacing).ceil().max(2.0) as usize;
    let noise = GradientNoise::new(draw ^ 0x57E5_5000);
    let turn_scale = TROUGH_WOBBLE_DEG.to_radians() * (0.5 + rng.uniform(draw, 8));
    let directions = |p: [f64; 3]| {
        // The radial direction outward from the centre and the hoop direction round it,
        // in the tangent plane at p.
        let c = rise.center;
        let radial_dot = dot(c, p);
        let out = [
            p[0] * radial_dot - c[0],
            p[1] * radial_dot - c[1],
            p[2] * radial_dot - c[2],
        ];
        let out = if dot(out, out) < 1e-20 {
            local_basis(p).0
        } else {
            normalize(out)
        };
        (out, normalize(cross(p, out)))
    };
    let distance = |p: [f64; 3]| angle_between(rise.center, p) * body_radius / rise.radius;
    // The first heading: across the most tensile direction, or the fault's own where the
    // stress is about the same every way.
    let (sr, st) = rise_stress(distance(start));
    let (out, round) = directions(start);
    let own = 2.0 * std::f64::consts::PI * rng.uniform(draw, 9);
    let mut heading = if (st - sr).abs() <= RIFT_LABYRINTH * sr.max(st).abs() {
        normalize([
            out[0] * own.cos() + round[0] * own.sin(),
            out[1] * own.cos() + round[1] * own.sin(),
            out[2] * own.cos() + round[2] * own.sin(),
        ])
    } else if st >= sr {
        out
    } else if rng.uniform(draw, 10) < 0.5 {
        round
    } else {
        [-round[0], -round[1], -round[2]]
    };
    let mut points = Vec::with_capacity(n + 1);
    let mut p = start;
    for i in 0..=n {
        points.push(p);
        let (sr, st) = rise_stress(distance(p));
        if i >= 2 && sr.max(st) < RIFT_MIN_TENSION {
            break;
        }
        let (out, round) = directions(p);
        let along = if st >= sr { out } else { round };
        let along = if dot(along, heading) < 0.0 {
            [-along[0], -along[1], -along[2]]
        } else {
            along
        };
        let anisotropy = ((st - sr).abs() / (st.abs() + sr.abs()).max(1e-9)).clamp(0.0, 1.0);
        let k = RIFT_STEER * anisotropy;
        let s = i as f64 * spacing;
        let turn = turn_scale
            * (noise.fbm([s / TROUGH_WOBBLE_LENGTH, 0.5, 0.25], 1.0, 2, 2.0, 0.5) / 0.6)
                .clamp(-1.0, 1.0)
            / n.max(1) as f64
            * 4.0;
        let side = normalize(cross(p, heading));
        let steered = [
            heading[0] + k * (along[0] - heading[0]),
            heading[1] + k * (along[1] - heading[1]),
            heading[2] + k * (along[2] - heading[2]),
        ];
        let steered = [
            steered[0] + turn * side[0],
            steered[1] + turn * side[1],
            steered[2] + turn * side[2],
        ];
        let radial = dot(steered, p);
        heading = normalize([
            steered[0] - radial * p[0],
            steered[1] - radial * p[1],
            steered[2] - radial * p[2],
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
                along: along_trace,
                taper,
                narrow,
                pulse,
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
                    let taper = (std::f64::consts::PI * s / trough.length).sin().max(0.0);
                    let fall = trough.depth * taper.powf(TROUGH_TAPER_POWER);
                    let (h_lo, h_hi) = LANDSLIDE_HALF_WIDTH;
                    let half = fall.max(LANDSLIDE_MIN_WALL)
                        * h_lo
                        * (h_hi / h_lo).powf(rng.uniform(draw, 0));
                    let (g_lo, g_hi) = LANDSLIDE_GAP;
                    s += half * (1.0 + 2.0 * (g_lo + (g_hi - g_lo) * rng.uniform(draw, 1)));
                    if s >= trough.length {
                        break;
                    }
                    let taper = (std::f64::consts::PI * s / trough.length).sin().max(0.0);
                    let fall = trough.depth * taper.powf(TROUGH_TAPER_POWER);
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
        }
    }

    /// What the landslides near a point do to it (F106), at the trough's growth: how far
    /// the alcoves cut the wall back on its side, m, and how thick their debris lies on
    /// it, m. `along` is its distance along the trace and `across` off it, the trough
    /// `width` wide with walls `wall` wide there.
    fn landslides_at(
        &self,
        along: f64,
        across: f64,
        width: f64,
        wall: f64,
        scale: f64,
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
            if l.side == side && ds.abs() < half {
                cut = cut.max(back * (1.0 - (ds / half).powi(2)));
            }
            // The debris: the alcove's volume, a half cone, in a lobe from the foot of its
            // wall across the floor, as far as the floor goes; it widens toward its rounded
            // toe and thins to its edges, and drapes the wall's lower slope above the foot.
            let run = (l.reach * fall).min(2.0 * foot).max(1.0);
            let into = foot - l.side * across;
            if into < -wall || into > run {
                continue;
            }
            let u = (into / run).max(0.0);
            let spread = half * (1.0 + LANDSLIDE_SPREAD * u) * (1.0 - u.powi(4)).max(0.0).sqrt();
            if ds.abs() >= spread {
                continue;
            }
            let volume = std::f64::consts::FRAC_PI_2 * half * back * fall / 3.0;
            let peak = volume / (self.lobe_shape * half * run);
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
                + if trough.depth >= LANDSLIDE_MIN_WALL {
                    LANDSLIDE_BACK.1 * LANDSLIDE_HALF_WIDTH.1 * depth
                } else {
                    0.0
                },
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
        // The landslides cut the wall back on their side and lay their debris on the floor.
        let (cut, debris) = self.landslides_at(
            f.along + along,
            across,
            width,
            size.wall,
            size.depth / self.depth.max(1.0),
        );
        let mut u = trough_profile(across, width + 2.0 * cut, size.depth, size.wall, blur);
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
        };
        let step = 250.0;
        let mut volume = 0.0;
        let mut a = l.along - 2.0 * l.half;
        while a < l.along + 2.0 * l.half {
            let mut x = -0.5 * width;
            while x < 0.5 * width {
                volume += only.landslides_at(a, x, width, wall, 1.0).1 * step * step;
                x += step;
            }
            a += step;
        }
        let alcove = std::f64::consts::FRAC_PI_2 * l.half * l.back * l.fall / 3.0;
        assert!((volume / alcove - 1.0).abs() < 0.1, "{volume} {alcove}");
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
        };
        let r = 3.39e6;
        let (east, north) = local_basis(rise.center);
        let start = offset(rise.center, east, north, 0.55 * rise.radius / r, 1.0);
        let trace = stress_trace(&rise, start, 800e3, 10e3, r, &Rng::new(4), 7);
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
