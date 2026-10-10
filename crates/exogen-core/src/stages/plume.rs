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
use crate::geom::{angle_between, cross, dot, local_basis, normalize, offset};
use crate::grid::CubeSphere;
use crate::noise::GradientNoise;
use crate::params::ParamGraph;
use crate::rng::Rng;
use crate::thermal::ThermalHistory;
use crate::units::UNIT_LAVA;

pub const VERSION: u64 = 1;

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
    /// Summit calderas: offset from the centre as a fraction of the radius along an
    /// azimuth, radius and depth, metres.
    pub calderas: Vec<Caldera>,
    pub id: u64,
}

#[derive(Clone, Copy, Debug)]
pub struct Caldera {
    pub offset: f64,
    pub azimuth: f64,
    pub radius: f64,
    pub depth: f64,
}

/// A trough along a trace: graben on the flanks and the canyon.
#[derive(Clone, Debug)]
pub struct Trough {
    pub points: Vec<[f64; 3]>,
    pub spacing: f64,
    pub length: f64,
    pub width: f64,
    pub depth: f64,
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
        1.3 * self.radius / body_radius
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

/// Shield growth: most of it early, a tail to the present.
fn shield_progress(t: f64, start: f64) -> f64 {
    if t <= start {
        0.0
    } else {
        1.0 - (-(t - start) / SHIELD_GROWTH_GYR).exp()
    }
}

/// Caldera depth at a point, metres below the summit surface: nested flat-floored pits.
fn caldera_depth(shield: &Shield, d: [f64; 3], body_radius: f64, blur: f64) -> f64 {
    let (east, north) = local_basis(shield.center);
    let mut depth: f64 = 0.0;
    for c in &shield.calderas {
        let centre = offset(
            shield.center,
            east,
            north,
            c.offset * shield.radius / body_radius,
            c.azimuth,
        );
        let r = angle_between(centre, d) * body_radius;
        let wall = (CALDERA_WALL_WIDTH * c.radius).hypot(blur);
        let w = 0.5 * (1.0 - ((r - c.radius) / wall).tanh());
        depth = depth.max(c.depth * w);
    }
    depth
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
    let r = angle_between(shield.center, d) * body_radius;
    let h = shield.height * shield_shape(r / shield.radius);
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
    /// Shield height so far per coarse cell, m.
    shield_built: Vec<f32>,
    /// Canyon depth carved so far per coarse cell, m.
    carved: Vec<f32>,
    stamper: Stamper,
}

impl PlumeProcess {
    /// None when the body raises no plume.
    pub fn new(graph: &ParamGraph, grid: CubeSphere, thermal: &ThermalHistory) -> Option<Self> {
        let count = graph.f64("tectonics.plume_count").round() as usize;
        if count == 0 || graph.str("interior.tectonic_regime") != "stagnant_lid" {
            return None;
        }
        let body_radius = graph.f64("body.radius");
        let age = graph.f64("body.age");
        let rng = Rng::new(graph.f64("seed.tectonics") as u64 ^ 0x9D1F_E5A7);
        let rise_radius = 0.5 * graph.f64("tectonics.rise_diameter");
        let rise_height = graph.f64("tectonics.rise_height");
        let mut rises = Vec::new();
        for i in 0..count as u64 {
            let mut center = rng.unit_direction(i, 1);
            // Rises keep apart.
            for _ in 0..16 {
                if rises.iter().all(|r: &Rise| {
                    angle_between(r.center, center) > 2.5 * rise_radius / body_radius
                }) {
                    break;
                }
                center = rng.unit_direction(i, 2 + rng.bits(i, 3) % 64);
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
        let relief_limit = thermal.at(age).relief_limit;
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
                let n_cal = 1 + (rng.bits(draw, 4) % 3) as usize;
                let mut calderas = Vec::with_capacity(n_cal);
                for c in 0..n_cal as u64 {
                    calderas.push(Caldera {
                        offset: if c == 0 {
                            0.0
                        } else {
                            CALDERA_RADIUS_FRACTION * (0.6 + 0.8 * rng.uniform(draw, 5 + 3 * c))
                        },
                        azimuth: 2.0 * std::f64::consts::PI * rng.uniform(draw, 6 + 3 * c),
                        radius: CALDERA_RADIUS_FRACTION
                            * radius
                            * (if c == 0 {
                                1.0
                            } else {
                                0.4 + 0.4 * rng.uniform(draw, 7 + 3 * c)
                            }),
                        depth: (CALDERA_DEPTH_FRACTION * height).min(CALDERA_MAX_DEPTH),
                    });
                }
                shields.push(Shield {
                    center,
                    radius,
                    height,
                    start_gyr: SHIELD_START_GYR.0
                        + (SHIELD_START_GYR.1 - SHIELD_START_GYR.0) * rng.uniform(draw, 14),
                    calderas,
                    id: sid,
                });
                sid += 1;
            }
        }
        // Canyons and graben on the first rise's flank.
        let cell = grid.delta() * body_radius;
        let canyon_count = graph.f64("tectonics.canyon_count").round() as usize;
        let mut canyons = Vec::new();
        let mut graben = Vec::new();
        if let Some(rise) = rises.first() {
            let (east, north) = local_basis(rise.center);
            for k in 0..canyon_count as u64 {
                let draw = 5_000 + k * 8;
                let phi = 2.0 * std::f64::consts::PI * rng.uniform(draw, 1);
                let length = CANYON_LENGTH_FRACTION * 2.0 * rise.radius;
                let start = offset(
                    rise.center,
                    east,
                    north,
                    CANYON_START_FRACTION * rise.radius / body_radius,
                    phi,
                );
                let points = radial_trace(
                    start,
                    rise.center,
                    length,
                    0.5 * cell,
                    body_radius,
                    &rng,
                    draw,
                );
                canyons.push(Trough {
                    points,
                    spacing: 0.5 * cell,
                    length,
                    width: CANYON_WIDTH,
                    depth: CANYON_DEPTH.min(0.9 * rise_height + CANYON_DEPTH * 0.5),
                    start_gyr: CANYON_START_GYR,
                    end_gyr: CANYON_END_GYR,
                    id: k,
                });
            }
            let n_graben = (GRABEN_PER_RISE * (rise.radius / 2.5e6).clamp(0.5, 2.0)) as u64;
            for k in 0..n_graben {
                let draw = 9_000 + k * 8;
                let phi = 2.0 * std::f64::consts::PI * rng.uniform(draw, 1);
                let r0 = rise.radius
                    * (GRABEN_START_FRACTION.0
                        + (GRABEN_START_FRACTION.1 - GRABEN_START_FRACTION.0)
                            * rng.uniform(draw, 2));
                let (l_lo, l_hi) = GRABEN_LENGTH;
                let length = l_lo * (l_hi / l_lo).powf(rng.uniform(draw, 3));
                let start = offset(rise.center, east, north, r0 / body_radius, phi);
                let spacing = 0.5 * cell;
                let points =
                    radial_trace(start, rise.center, length, spacing, body_radius, &rng, draw);
                let (w_lo, w_hi) = GRABEN_WIDTH;
                let width = w_lo * (w_hi / w_lo).powf(rng.uniform(draw, 4));
                let time = rise.start_gyr + (rise.end_gyr - rise.start_gyr) * rng.uniform(draw, 5);
                graben.push(Trough {
                    points,
                    spacing,
                    length,
                    width,
                    depth: GRABEN_DEPTH * (0.5 + rng.uniform(draw, 6)),
                    start_gyr: time,
                    end_gyr: time,
                    id: 1_000 + k,
                });
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
            shield_built: vec![0.0; grid.cell_count()],
            carved: vec![0.0; grid.cell_count()],
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
                .gather(&grid, body_radius, shield.center, shield.radius);
            let bedrock = layers.bedrock.as_mut_slice();
            let ages = layers.surface_age.as_mut_slice();
            let units = layers.unit_id.as_mut_slice();
            let deep = layers.deep.as_mut_slice();
            let melt = layers.melt.as_mut_slice();
            for &(idx, _, _) in &self.stamper.footprint {
                let d = grid.center(grid.cell(idx));
                let full = shield_height(shield, d, body_radius, false, 0.0);
                let dh = full * dp;
                if dh <= 0.0 {
                    continue;
                }
                bedrock[idx] += dh as f32;
                self.shield_built[idx] += dh as f32;
                if dh > SHIELD_RENEWAL_THICKNESS {
                    ages[idx] = t_mid as f32;
                    units[idx] = UNIT_LAVA;
                    deep[idx] = 0.0;
                    melt[idx] = 0.0;
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

/// A trace running outward from `start`, away from `center`, with a slow wobble.
fn radial_trace(
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
    let phase = 2.0 * std::f64::consts::PI * rng.uniform(draw, 7);
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
        let turn =
            turn_scale * (2.0 * std::f64::consts::PI * s / TROUGH_WOBBLE_LENGTH + phase).sin();
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
    let noise = GradientNoise::new(trough.id ^ 0xCA17_0000);
    let scallop_frequency = body_radius / TROUGH_SCALLOP_LENGTH.max(trough.width);
    let width_now = trough.width * (0.3 + 0.7 * progress);
    let depth_now = trough.depth * progress;
    // Walls no narrower than the grid that carries them.
    let wall = (TROUGH_WALL_FRACTION * width_now).max(wall_floor);
    let search = 0.5 * width_now + 3.0 * wall + 2.0 * cell;
    let half_band = 0.5 * trough.spacing + cell;
    let mut best: HashMap<usize, (f64, f64)> = HashMap::new();
    let n = trough.points.len();
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
        let taper = taper.powf(TROUGH_TAPER_POWER);
        // Chasmata: the width pulses along the trace, and each wall is scalloped.
        let pulse = 1.0
            + TROUGH_WIDTH_PULSE
                * (2.0 * std::f64::consts::PI * along_trace / TROUGH_PULSE_LENGTH).sin();
        stamper.gather(grid, body_radius, p, search);
        for &(idx, _, _) in &stamper.footprint {
            let d = grid.center(grid.cell(idx));
            let along = dot(d, t) * body_radius;
            if along.abs() > half_band {
                continue;
            }
            let across = dot(d, normal) * body_radius;
            let scallop = 1.0
                + TROUGH_SCALLOP
                    * noise.fbm(
                        [d[0] + across.signum() * 3.0, d[1], d[2]],
                        scallop_frequency,
                        2,
                        2.0,
                        0.5,
                    );
            let width = width_now * pulse * scallop;
            let u = taper * trough_profile(across, width, depth_now, wall, blur);
            let entry = best.entry(idx).or_insert((f64::INFINITY, 0.0));
            if along.abs() < entry.0 {
                *entry = (along.abs(), u);
            }
        }
    }
    let mut out: Vec<(usize, f64)> = best.into_iter().map(|(i, (_, u))| (i, u)).collect();
    out.sort_by_key(|e| e.0);
    out
}

/// Finishing pass of a trough: the fine profile minus the version the loop carried,
/// with the loop's wall width and the coarse smear, or the whole trough when the loop
/// did not carry it. Inside the floor the fine relief of older features is removed, as a
/// crater's cavity removes it.
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
    stamper: &mut Stamper,
) -> StampOutput {
    let bedrock = layers.bedrock.as_slice();
    let ages = layers.surface_age.as_slice();
    let cell = grid.delta() * body_radius;
    let fine = trough_depths(grid, body_radius, trough, 1.0, cell, blur_fine, stamper);
    let blur_coarse = coarse.map(|c| c.1);
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
            let smeared = coarse.get(&idx).copied().unwrap_or(0.0);
            // How deep inside the trough this cell is, 0 on the rim to 1 on the floor.
            let inside = (depth / trough.depth.max(1.0)).clamp(0.0, 1.0);
            let kept = if blur_coarse.is_some() {
                // The loop carried the trough: keep older fine relief outside the floor.
                (bedrock[idx] - base[idx]) as f64 * (1.0 - inside)
            } else {
                (bedrock[idx] - base[idx]) as f64
            };
            let reference = if blur_coarse.is_some() {
                base[idx] as f64
            } else {
                bedrock[idx] as f64
            };
            let mut write = CellWrite {
                idx: idx as u32,
                bedrock: (reference + kept - (depth - smeared)) as f32,
                surface: false,
                age: 0.0,
                unit: 0,
                deep: 0.0,
                melt: 0.0,
            };
            if inside > 0.5 && ages[idx] <= time_gyr + 1e-6 {
                write.surface = true;
                write.age = time_gyr;
                write.unit = layers.unit_id.as_slice()[idx];
                write.deep = layers.deep.as_slice()[idx];
                write.melt = layers.melt.as_slice()[idx];
            }
            write
        })
        .collect();
    StampOutput { cells }
}

/// Finishing pass of a shield: the calderas cut into the summit and the flank texture,
/// neither of which the coarse grid carried.
pub fn shield_compute(
    layers: &LayerSet,
    grid: &CubeSphere,
    body_radius: f64,
    shield: &Shield,
    stamper: &mut Stamper,
) -> StampOutput {
    let bedrock = layers.bedrock.as_slice();
    stamper.gather(grid, body_radius, shield.center, shield.radius * 1.05);
    let noise = GradientNoise::new(shield.id ^ 0x5B1E_1D00);
    let flow_frequency = body_radius / SHIELD_FLOW_TEXTURE_LENGTH;
    let cells = stamper
        .footprint
        .iter()
        .map(|&(idx, r, phi)| {
            let d = grid.center(grid.cell(idx));
            let fine = shield_height(shield, d, body_radius, true, 0.0);
            let coarse = shield_height(shield, d, body_radius, false, 0.0);
            let x = r / shield.radius;
            // Radial flow ridges down the flank, fading at the summit and the base.
            let flank = (4.0 * x * (1.0 - x)).clamp(0.0, 1.0);
            let streak = noise.fbm(
                [
                    d[0] * flow_frequency + 11.0 * phi.cos(),
                    d[1] * flow_frequency + 11.0 * phi.sin(),
                    d[2] * flow_frequency,
                ],
                1.0,
                3,
                2.0,
                0.5,
            );
            let texture = SHIELD_FLOW_TEXTURE * shield.height * flank * streak * 0.01;
            CellWrite {
                idx: idx as u32,
                bedrock: (bedrock[idx] as f64 + (fine - coarse) + texture) as f32,
                surface: false,
                age: 0.0,
                unit: 0,
                deep: 0.0,
                melt: 0.0,
            }
        })
        .collect();
    StampOutput { cells }
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
                offset: 0.0,
                azimuth: 0.0,
                radius: 30e3,
                depth: 2.5e3,
            }],
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
                offset(shield.center, east, north, 310e3 / r, 0.0),
                r,
                true,
                0.0
            ) == 0.0
        );
    }
}
