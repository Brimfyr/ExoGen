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
use crate::geom::{angle_between, azimuth_from, cross, dot, local_basis, normalize};
use crate::grid::CubeSphere;
use crate::grid::projection::{PaddedField, sample_bilinear};
use crate::noise::GradientNoise;
use crate::params::ParamGraph;
use crate::rng::Rng;
use crate::units::*;
use rayon::prelude::*;
use std::collections::HashMap;

pub const VERSION: u64 = 10;

/// A patera: a volcano-tectonic depression with a flat floor and a fault-bounded,
/// irregular outline.
#[derive(Clone, Debug)]
pub struct Patera {
    pub center: [f64; 3],
    /// Mean radius and the depth of the floor below the rim, m.
    pub radius: f64,
    pub depth: f64,
    /// A lava lake on the floor at present.
    pub active: bool,
    /// Years since the floor last held a lake: it is coated by what has fallen on it
    /// since (F84).
    pub floor_age_yr: f64,
    /// Gas-rich, Pele-class, whose plumes are hot and red; the rest are Prometheus-class.
    pub pele: bool,
    /// The outline: a polygon with this many sides at this rotation, under noise.
    pub sides: u32,
    pub rotation: f64,
    pub birth_gyr: f64,
    pub death_gyr: f64,
    pub id: u64,
}

/// A mountain: a tilted crustal block.
#[derive(Clone, Debug)]
pub struct Mountain {
    pub center: [f64; 3],
    pub radius: f64,
    pub height: f64,
    /// Direction the block tilts down towards, radians in the local frame, and how much
    /// of the height the tilt moves across the block.
    pub tilt_azimuth: f64,
    pub tilt: f64,
    /// Elongation along the tilt azimuth, and the seed of the outline noise.
    pub aspect: f64,
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

/// A compound flow field spreading from a patera.
#[derive(Clone, Debug)]
pub struct FlowField {
    pub center: [f64; 3],
    pub reach: f64,
    pub lobes: Vec<Lobe>,
    pub birth_gyr: f64,
    pub death_gyr: f64,
    pub id: u64,
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
    pub id: u64,
}

impl Patera {
    pub fn cap_radius(&self, body_radius: f64) -> f64 {
        1.5 * self.radius / body_radius
    }
    pub fn alive(&self, t: f64) -> bool {
        self.birth_gyr - IO_ALIVE_TOLERANCE_GYR <= t && t < self.death_gyr - IO_ALIVE_TOLERANCE_GYR
    }
}
impl Mountain {
    pub fn cap_radius(&self, body_radius: f64) -> f64 {
        1.2 * self.radius * self.aspect.sqrt() * (1.0 + IO_MOUNTAIN_OUTLINE) / body_radius
    }
    pub fn alive(&self, t: f64) -> bool {
        self.birth_gyr - IO_ALIVE_TOLERANCE_GYR <= t && t < self.death_gyr - IO_ALIVE_TOLERANCE_GYR
    }
}
impl FlowField {
    pub fn cap_radius(&self, body_radius: f64) -> f64 {
        let mut far: f64 = 0.0;
        for l in &self.lobes {
            let hw = l.half_width * (1.0 + IO_FLOW_LOBING + IO_FLOW_WIDTH_WOBBLE) * 1.5;
            for p in &l.points {
                far = far.max(angle_between(self.center, *p) + hw / body_radius);
            }
        }
        far.max(2.0 * IO_FLOW_MIN_HALF_WIDTH / body_radius)
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
        let (red, inner) = if e.pele {
            (1.0, IO_RING_INNER_PELE)
        } else {
            (IO_PROMETHEUS_RED, IO_RING_INNER_PROMETHEUS)
        };
        let (start, end) = (-e.start_yr.min(window), -e.end_yr);
        // The deposit reaches as far as its ring's tail lays a fall that could show over
        // the eruption, as a share of its range, at most one and a half ranges.
        let total = e.crest * (end - start).max(0.0) / NEGLIGIBLE_FALL;
        let reach = (1.0 + IO_RING_OUTER_WIDTH * total.max(1.0).ln().sqrt()).min(1.5);
        let wobble: f64 = e.range.iter().map(|h| h[0].hypot(h[1])).sum();
        let cap = reach * e.radius * e.aspect.sqrt() * (1.0 + wobble) / body_radius;
        Self {
            cos_cap: cap.min(std::f64::consts::PI).cos(),
            reach,
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

/// Radius of a patera's outline at an azimuth: a polygon's edge, as the faults that
/// bound a collapse leave it, under noise.
pub fn patera_radius(phi: f64, p: &Patera, noise: &GradientNoise) -> f64 {
    let sector = 2.0 * std::f64::consts::PI / p.sides as f64;
    let a = (phi - p.rotation).rem_euclid(sector) - 0.5 * sector;
    let polygon = (0.5 * sector).cos() / a.cos();
    let wobble = 1.0
        + IO_PATERA_OUTLINE * noise.fbm([phi.cos() * 1.3, phi.sin() * 1.3, 0.37], 2.0, 2, 2.0, 0.5);
    p.radius * (IO_PATERA_POLYGON * polygon + (1.0 - IO_PATERA_POLYGON)) * wobble
}

/// Depth of a flat-floored pit at a distance from its centre, with walls of a width.
pub fn pit_profile(r: f64, radius: f64, depth: f64, wall: f64) -> f64 {
    depth * 0.5 * (1.0 - ((r - radius) / wall).tanh())
}

/// The block's radius at an azimuth: an ellipse along the tilt azimuth with a noisy
/// outline.
pub fn block_radius(phi: f64, m: &Mountain, noise: &GradientNoise) -> f64 {
    let t = phi - m.tilt_azimuth;
    let ellipse = 1.0 / (t.cos().powi(2) + (m.aspect * t.sin()).powi(2)).sqrt() * m.aspect.sqrt();
    let wobble = 1.0
        + IO_MOUNTAIN_OUTLINE
            * noise.fbm([phi.cos() * 1.5, phi.sin() * 1.5, 0.3], 2.0, 2, 2.0, 0.5);
    m.radius * ellipse * wobble
}

/// Height of a tilted block at a point: a plateau sloping down towards the tilt azimuth,
/// dropping to nothing at its noisy edge.
pub fn block_profile(r: f64, phi: f64, m: &Mountain, noise: &GradientNoise) -> f64 {
    let edge = block_radius(phi, m, noise);
    if r >= edge {
        return 0.0;
    }
    let x = r / edge;
    let tilt = 1.0 - m.tilt * 0.5 * (1.0 + x * (phi - m.tilt_azimuth).cos());
    // The upthrown side ends in a fault scarp; the downthrown slopes under the plains.
    let up = 0.5 * (1.0 - (phi - m.tilt_azimuth).cos());
    let power = IO_MOUNTAIN_APRON_POWER + (IO_MOUNTAIN_SCARP_POWER - IO_MOUNTAIN_APRON_POWER) * up;
    let rim = (1.0 - x.powf(power)).max(0.0);
    m.height * tilt * rim
}

/// A fall too thin to show over its whole eruption, optical depths.
const NEGLIGIBLE_FALL: f64 = 0.02;

/// Bins of a young flow's age as shares of the history, for its frost halo: each bin's
/// halo takes the release-weighted mean age of the lava in it.
const HALO_AGE_BINS: [f64; 4] = [0.0, 0.02, 0.1, 1.0];

/// How far back the history reaches, years (F84): a few burial times of the background
/// fallout, past which nothing a plume laid can still show.
pub fn history_years() -> f64 {
    IO_HISTORY_BURIALS / IO_BACKGROUND_FALLOUT
}

/// Fallout of a plume at a fraction `x` of its range (F84), as the crest's share: the
/// ring where the gas condenses at the canopy and falls back, ending steeply beyond the
/// range of the fastest and rising gently within; and, as a share of `inner`, the
/// pyroclasts' thinner fall within it.
pub fn deposit_profile(x: f64) -> (f64, f64) {
    let width = if x < 1.0 {
        IO_RING_WIDTH
    } else {
        IO_RING_OUTER_WIDTH
    };
    let ring = (-((x - 1.0) / width).powi(2)).exp();
    let fill = ((1.0 - x) / 0.5).clamp(0.0, 1.0);
    (ring, fill)
}

/// The area under a deposit's profile, as the integral of the profile times x over x: a
/// deposit of range R spreads its mass over 2π R² times this.
pub fn profile_moment(inner: f64) -> f64 {
    let n = 4000;
    let dx = 1.5 / n as f64;
    (0..n)
        .map(|i| {
            let x = (i as f64 + 0.5) * dx;
            let (ring, fill) = deposit_profile(x);
            (ring + inner * fill) * x * dx
        })
        .sum()
}

/// The plumes of every patera over the history the surface remembers (F84). Each
/// patera's eruptions are drawn from the present back: one erupting now if the patera is
/// active and its class's share of active centres throws a plume, part way through; then
/// quiet spells and eruptions alternating at the duty its class keeps, the spells
/// memoryless and the eruptions from a few months to decades long. Each eruption has its
/// own speed and so its range, its own mass and its vent's elongation.
pub fn plume_history(paterae: &[Patera], gravity: f64, seed: u64) -> Vec<PlumeEpisode> {
    let rng = Rng::new(seed ^ 0x0E91_5DE5_0000_0084);
    let (d_lo, d_hi) = IO_PLUME_DURATION_YR;
    let mean_duration = (d_hi - d_lo) / (d_hi / d_lo).ln();
    let moments = [
        profile_moment(IO_RING_INNER_PROMETHEUS),
        profile_moment(IO_RING_INNER_PELE),
    ];
    let tau = 2.0 * std::f64::consts::PI;
    let window = history_years();
    let mut out = Vec::new();
    for p in paterae {
        let idx = p.id;
        let share = if p.pele {
            IO_PLUME_FRACTION_PELE
        } else {
            IO_PLUME_FRACTION_PROMETHEUS
        };
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
            let (speed, mass, elongation) = if p.pele {
                (
                    IO_PLUME_SPEED_PELE,
                    IO_PLUME_MASS_PELE,
                    IO_VENT_ELONGATION.0,
                )
            } else {
                (
                    IO_PLUME_SPEED_PROMETHEUS,
                    IO_PLUME_MASS_PROMETHEUS,
                    IO_VENT_ELONGATION.1,
                )
            };
            let speed = speed * (IO_PLUME_SPEED_SCATTER * rng.normal(idx, base)).exp();
            let radius = speed * speed / gravity.max(0.01);
            let mass = mass * (IO_PLUME_MASS_SCATTER * rng.normal(idx, base + 2)).exp();
            let crest = mass * YEAR
                / (tau * radius * radius * moments[p.pele as usize] * IO_DEPOSIT_OPTICAL_MASS);
            let dark = if p.pele {
                IO_DARK_PELE
            } else if rng.uniform(idx, base + 4) < IO_DARK_CONTAMINATED {
                IO_DARK_PROMETHEUS
            } else {
                0.0
            };
            let harmonic = |amplitude: f64, draw: u64| {
                let phase = tau * rng.uniform(idx, draw);
                [amplitude * phase.cos(), amplitude * phase.sin()]
            };
            let n_range = RANGE_ORDERS.1 - RANGE_ORDERS.0 + 1;
            let a_range = IO_RING_EDGE_WOBBLE * (2.0 / n_range as f64).sqrt();
            let mut range = [[0.0; 2]; RANGE_ORDERS.1 - RANGE_ORDERS.0 + 1];
            for (k, h) in range.iter_mut().enumerate() {
                *h = harmonic(a_range, base + 10 + k as u64);
            }
            out.push(PlumeEpisode {
                center: p.center,
                radius,
                pele: p.pele,
                start_yr: start,
                end_yr: end,
                crest,
                aspect: 1.0 + elongation * rng.uniform(idx, base + 5),
                // The deposit reaches farthest across the vent, which runs along the
                // patera's long side.
                axis: p.rotation + 0.5 * std::f64::consts::PI,
                dark,
                frost: if p.pele {
                    IO_PELE_FROST
                } else {
                    let (lo, hi) = IO_PROMETHEUS_FROST;
                    lo + (hi - lo) * rng.uniform(idx, base + 6)
                },
                range,
                mottle: [
                    1e3 * rng.uniform(idx, base + 40),
                    1e3 * rng.uniform(idx, base + 41),
                    1e3 * rng.uniform(idx, base + 42),
                ],
                id: idx ^ (base << 32),
            });
        }
    }
    out
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
/// it; longer slopes, a mountain's flanks, it reads one by one.
const FROST_SLOPE_LENGTH: f64 = 30e3;
const FROST_SPREAD: [f64; 6] = [0.333, 0.517, 0.803, 1.246, 1.934, 3.0];

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
/// over all, optical depths, at each of the cell's spread of sublimation rates, under
/// the cover of sulfur and radiation products laid on it since it fell.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
struct Top {
    lava: f64,
    dark: f64,
    red: f64,
    altered: f64,
    frost: [f64; 6],
    lag: f64,
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
        }
    }

    /// The share of the ground the frost hides, over the cell's spread of sublimation.
    fn visible_frost(&self) -> f64 {
        let columns: f64 = self.frost.iter().map(|f| 1.0 - (-f).exp()).sum();
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
/// landing round the young lava rather than on it.
pub fn deposit_history(
    layers: &mut LayerSet,
    body_radius: f64,
    insolation: f64,
    episodes: &[PlumeEpisode],
    slope: Option<&Field<f32>>,
) {
    let grid = layers.grid;
    let window = history_years();
    let plumes: Vec<Plume> = episodes
        .iter()
        .map(|e| Plume::new(e, body_radius, window))
        .collect();
    // The vents, each patera's episodes being consecutive.
    let mut vents: Vec<Vent> = Vec::new();
    let mut i = 0;
    while i < episodes.len() {
        let center = episodes[i].center;
        let mut j = i;
        let mut cos_cap: f64 = 1.0;
        while j < episodes.len() && episodes[j].center == center {
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
    let background = Fall {
        solid: IO_BACKGROUND_FALLOUT,
        ..Fall::default()
    };
    let units = layers.unit_id.as_slice();
    let tilts = slope.map(|s| smoothed_over(s, FROST_SLOPE_LENGTH, body_radius));
    let coverage = layers.cover.as_slice();
    let laid_yr = layers.melt.as_slice();
    // The young flows' halos: their release in bins of age, and the release-weighted
    // age, each spread over the halo's reach; the lava itself in each bin takes no halo
    // of its own.
    let halo_clock = std::time::Instant::now();
    let halos: Vec<(Vec<f32>, Vec<f32>, Vec<f32>)> = HALO_AGE_BINS
        .windows(2)
        .map(|bin| {
            let mut lava = Field::<f32>::new(grid);
            let mut aged = Field::<f32>::new(grid);
            {
                let (l, a) = (lava.as_mut_slice(), aged.as_mut_slice());
                for k in 0..l.len() {
                    let age = laid_yr[k] as f64 / window;
                    if units[k] == UNIT_IO_FLOW && age >= bin[0] && age < bin[1] {
                        l[k] = coverage[k].clamp(0.0, 1.0);
                        a[k] = l[k] * laid_yr[k];
                    }
                }
            }
            (
                smoothed_over(&lava, IO_HALO_RANGE, body_radius),
                smoothed_over(&aged, IO_HALO_RANGE, body_radius),
                lava.as_slice().to_vec(),
            )
        })
        .collect();
    let halo_rate = IO_HALO_RELEASE / (IO_HALO_DURATION_YR * IO_DEPOSIT_OPTICAL_MASS);
    let debug = std::env::var_os("EXOGEN_DEBUG").is_some();
    if debug {
        eprintln!(
            "io deposit history: the halos and the slopes {:.1} s",
            halo_clock.elapsed().as_secs_f64()
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
                    Vec::<(f64, f64, Fall)>::new(),
                    Vec::<(f64, u8, usize)>::new(),
                    Vec::<usize>::new(),
                )
            },
            |(falls, marks, active), k| {
                let d = grid.center(grid.cell(k));
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
                    for e in &plumes[vent.episodes.clone()] {
                        if c < e.cos_cap || e.end <= e.start {
                            continue;
                        }
                        let x = r / (e.radius * e.range_at(&cos_k, &sin_k));
                        if x >= e.reach {
                            continue;
                        }
                        let (ring, fill) = deposit_profile(x);
                        let fill = e.inner * fill;
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
                        let condensed = ring * (1.0 - e.frost);
                        let dark = e.dark * (-(x / IO_DARK_REACH).powi(2)).exp();
                        falls.push((
                            e.start,
                            e.end,
                            Fall {
                                solid: condensed + fill,
                                red: condensed * e.red,
                                dark: fill * dark,
                                frost: ring * e.frost,
                            },
                        ));
                    }
                }
                // The halos of the young flows about the cell.
                for (spread, spread_age, lava) in &halos {
                    let s = spread[k] as f64;
                    let rate = halo_rate * s * (1.0 - lava[k] as f64);
                    if rate * IO_HALO_DURATION_YR < NEGLIGIBLE_FALL {
                        continue;
                    }
                    let start = -(spread_age[k] as f64 / s).clamp(0.0, window);
                    let end = (start + IO_HALO_DURATION_YR).min(0.0);
                    if end > start {
                        falls.push((
                            start,
                            end,
                            Fall {
                                frost: rate,
                                ..Fall::default()
                            },
                        ));
                    }
                }
                // Where the cell stood when the history began: what fell before has aged,
                // its red reverted and its frost gone or covered, and it was buried at the
                // history's mean rate.
                let mut mean = background;
                for (start, end, fall) in falls.iter() {
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
                for (n, (start, end, _)) in falls.iter().enumerate() {
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
                    for (start, end, fall) in falls.iter() {
                        eprintln!("io trace:   fall {start:.1} to {end:.1} yr: {fall:?}");
                    }
                    for (b, (spread, spread_age, lava)) in halos.iter().enumerate() {
                        eprintln!(
                            "io trace:   halo bin {b}: spread {:.3e} age sum {:.3e} own {:.3}",
                            spread[k], spread_age[k], lava[k]
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
    }
    let deep = layers.deep.as_mut_slice();
    let melt = layers.melt.as_mut_slice();
    let cover = layers.cover.as_mut_slice();
    let stain = layers.stain.as_mut_slice();
    for (k, top) in tops.iter().enumerate() {
        deep[k] = (top.lava + (1.0 - top.lava) * top.dark) as f32;
        melt[k] = top.red as f32;
        stain[k] = top.altered as f32;
        cover[k] = top.visible_frost() as f32;
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
    /// The terrain before the regime's features, which every step starts from.
    base: Option<LayerSet>,
    pub paterae: Vec<Patera>,
    pub mountains: Vec<Mountain>,
    pub flows: Vec<FlowField>,
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
            base: None,
            paterae: Vec::new(),
            mountains: Vec::new(),
            flows: Vec::new(),
            stamper: Stamper::new(&grid),
        })
    }

    fn log_uniform(&self, index: u64, draw: u64, range: (f64, f64)) -> f64 {
        range.0 * (range.1 / range.0).powf(self.rng.uniform(index, draw))
    }

    /// A lobe's centre line (F59): from the patera's rim at an azimuth, downslope over
    /// the ground with the inertia of a flow and a wander, deflecting along a contour
    /// but never turning back, stopping where a mountain stands, where it would have to
    /// climb out of the ground it has ponded in, or at its length. Returns the lobe and
    /// how it ended: 0 at its length, 1 at a mountain, 2 ponded.
    #[allow(clippy::too_many_arguments)]
    fn trace_lobe(
        &self,
        idx: u64,
        draw: u64,
        center: [f64; 3],
        rim: f64,
        reach: f64,
        ground: &PaddedField,
        standing: &[(Mountain, GradientNoise)],
        age_yr: f64,
    ) -> (Lobe, usize) {
        let r = self.body_radius;
        let grid = &self.grid;
        let (east, north) = local_basis(center);
        let a0 = 2.0 * std::f64::consts::PI * self.rng.uniform(idx, draw);
        let mut p = crate::geom::offset(center, east, north, rim / r, a0);
        let outward = |p: [f64; 3]| {
            let c = dot(center, p);
            normalize([
                c * p[0] - center[0],
                c * p[1] - center[1],
                c * p[2] - center[2],
            ])
        };
        let mut h = outward(p);
        let length = reach * (0.3 + 0.7 * self.rng.uniform(idx, draw + 2));
        let spacing = (length / 64.0).max(IO_FLOW_POINT_SPACING);
        let steps = (length / spacing).ceil() as usize;
        let (w_lo, w_hi) = IO_FLOW_WIDTH_SHARE;
        let half_width = (reach * (w_lo + (w_hi - w_lo) * self.rng.uniform(idx, draw + 1)))
            .max(IO_FLOW_MIN_HALF_WIDTH);
        let noise = GradientNoise::new(idx ^ (draw << 40) ^ 0x0F10_7ACE);
        let delta = (0.5 * spacing / r).max(1e-7);
        let mut points = vec![p];
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
            let down = if grad > 1e-12 {
                normalize([
                    -(hx * e[0] + hy * n[0]),
                    -(hx * e[1] + hy * n[1]),
                    -(hx * e[2] + hy * n[2]),
                ])
            } else {
                h
            };
            let follow = IO_FLOW_SLOPE_FOLLOW * (grad / IO_FLOW_SLOPE_SCALE).tanh();
            // The fall line steers the flow sideways and speeds it on, but never turns it
            // back: a front deflects along a contour and ponds rather than reverse.
            let along = dot(down, h).max(0.0);
            let lateral = {
                let a = dot(down, h);
                [down[0] - a * h[0], down[1] - a * h[1], down[2] - a * h[2]]
            };
            // The wander: a turn per step that sums to the wander over its length.
            let turn = IO_FLOW_WANDER_RAD
                * (spacing / IO_FLOW_WANDER_LENGTH)
                * noise.fbm(
                    [i as f64 * spacing / IO_FLOW_WANDER_LENGTH, 0.3, 0.7],
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
            // A mountain in the way stops the flow.
            let blocked = standing.iter().any(|(m, mn)| {
                let dist = angle_between(m.center, q) * r;
                dist < m.radius * m.aspect.sqrt() * (1.0 + IO_MOUNTAIN_OUTLINE) && {
                    let (me, mnn) = local_basis(m.center);
                    dist < block_radius(azimuth_from(m.center, me, mnn, q), m, mn)
                }
            });
            if blocked {
                end = 1;
                break;
            }
            // Lava ponds: it crosses a bump no higher than its own pile above the lowest
            // ground it has reached, and stops where it would have to climb more.
            let height = sample_bilinear(grid, ground, q) as f64;
            if height > lowest + IO_FLOW_FILL_DEPTH {
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
        // The width along the lobe: a channel at the vent widening to the toe, wandering.
        let count = points.len();
        let widths = (0..count)
            .map(|i| {
                let s = i as f64 / (count.max(2) - 1) as f64;
                let wander = 1.0
                    + IO_FLOW_WIDTH_WOBBLE
                        * noise.fbm(
                            [i as f64 * spacing / (3.0 * half_width), 1.3, 2.1],
                            1.0,
                            2,
                            2.0,
                            0.5,
                        )
                        / 0.6;
                (wander.max(0.4) * (1.0 - IO_FLOW_TOE + IO_FLOW_TOE * s)) as f32
            })
            .collect();
        (
            Lobe {
                points,
                spacing,
                half_width,
                widths,
                age_yr,
            },
            end,
        )
    }

    /// The regime's surface afresh, as every step lays it: the plains' age and unit,
    /// nothing on them.
    pub fn clear_surface(layers: &mut LayerSet, t: f64) {
        let now = (t - IO_SURFACE_AGE_GYR) as f32;
        let ages = layers.surface_age.as_mut_slice();
        let units = layers.unit_id.as_mut_slice();
        let deep = layers.deep.as_mut_slice();
        let melt = layers.melt.as_mut_slice();
        let cover = layers.cover.as_mut_slice();
        let stain = layers.stain.as_mut_slice();
        for k in 0..ages.len() {
            ages[k] = now;
            units[k] = UNIT_PRIMORDIAL;
            deep[k] = 0.0;
            melt[k] = 0.0;
            cover[k] = 0.0;
            stain[k] = 0.0;
        }
    }

    /// The mountains through the body's life: each of the standing count is a chain of
    /// blocks, one raised as the last is buried, drawn from the seed alone.
    fn raise_mountains(&mut self) {
        for slot in 0..self.n_mountains as u64 {
            let mut generation = 0u64;
            let first = 5_000 + slot * 1_000;
            let mut birth =
                -self.log_uniform(first, 7, IO_MOUNTAIN_LIFETIME_GYR) * self.rng.uniform(first, 8);
            while birth < self.age && generation < 900 {
                let k = first + generation;
                let (lo, hi) = IO_MOUNTAIN_DIAMETER;
                let diameter = lo * (hi / lo).powf(self.rng.uniform(k, 1));
                let (h_lo, h_hi) = IO_MOUNTAIN_HEIGHT;
                let life = self.log_uniform(k, 7, IO_MOUNTAIN_LIFETIME_GYR);
                self.mountains.push(Mountain {
                    center: self.rng.unit_direction(k, 9),
                    radius: 0.5 * diameter,
                    height: h_lo * (h_hi / h_lo).powf(self.rng.uniform(k, 3)),
                    tilt_azimuth: 2.0 * std::f64::consts::PI * self.rng.uniform(k, 4),
                    tilt: IO_MOUNTAIN_TILT * self.rng.uniform(k, 5),
                    aspect: 1.0 + (IO_MOUNTAIN_ASPECT - 1.0) * self.rng.uniform(k, 6),
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
    fn populate(&mut self, t1: f64, lay: bool, bedrock: &Field<f32>) {
        // The time bin the step ends in is the random index, so the same surface comes
        // out at every quality.
        let step = (t1 / self.age * crate::stages::LOOP_TIME_BINS as f64).round() as u64;
        let death = t1 + self.dt.max(1e-6);
        // Paterae, with their flows, for the steps whose surface is kept.
        if !lay {
            return;
        }
        // The ground the flows run over, the regional slope of the base, and the
        // mountains that stop them.
        let cell = self.grid.delta() * self.body_radius;
        let smoothed = crate::stages::relaxation::low_pass(bedrock, cell, IO_FLOW_GROUND_SMOOTH);
        let ground = PaddedField::new(&self.grid, &smoothed);
        let mut ends = [0usize; 3];
        let mut lobe_lengths = (0usize, 0.0f64);
        let standing: Vec<(Mountain, GradientNoise)> = self
            .mountains
            .iter()
            .filter(|m| m.alive(t1))
            .map(|m| (m.clone(), GradientNoise::new(m.id ^ 0x0B10_C000)))
            .collect();
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
            let (lo, hi) = IO_PATERA_DIAMETER;
            let diameter = lo * (hi / lo).powf(self.rng.uniform(idx, 40).powf(IO_PATERA_SIZE_SKEW));
            let radius = 0.5 * diameter;
            let pele = self.rng.uniform(idx, 43) < IO_PELE_CLASS_FRACTION;
            let active = self.rng.uniform(idx, 42) < IO_PATERA_ACTIVE_FRACTION;
            let floor_age_yr = if active {
                0.0
            } else {
                self.log_uniform(idx, 44, IO_FLOOR_AGE_YR)
            };
            let (s_lo, s_hi) = IO_PATERA_SIDES;
            self.paterae.push(Patera {
                center,
                radius,
                depth: IO_PATERA_DEPTH * (0.5 + self.rng.uniform(idx, 41)),
                active,
                floor_age_yr,
                pele,
                sides: s_lo + (self.rng.bits(idx, 45) % (s_hi - s_lo + 1) as u64) as u32,
                rotation: 2.0 * std::f64::consts::PI * self.rng.uniform(idx, 46),
                birth_gyr: t1,
                death_gyr: death,
                id: idx,
            });
            // Flow fields: the active paterae's and a few quiet ones' now, and every
            // patera's older generations, lobes of all ages coated by what has fallen on
            // them since (F56, F84).
            for g in 0..IO_FLOW_GENERATIONS {
                let current = g == 0;
                if current && !(active || self.rng.uniform(idx, 47) < IO_QUIET_FLOW_FRACTION) {
                    continue;
                }
                let base = 60 + 40 * g as u64;
                let (n_lo, n_hi) = IO_FLOW_LOBES;
                let n =
                    n_lo + (self.rng.bits(idx, 48 + g as u64) % (n_hi - n_lo + 1) as u64) as u32;
                let (a_lo, a_hi) = IO_FLOW_AGE_YR;
                let age_range = if current {
                    (a_lo, a_hi)
                } else {
                    (IO_FLOW_EARLIER_YR, a_hi)
                };
                let reach = self
                    .log_uniform(idx, 49 + 90 * g as u64, IO_FLOW_REACH)
                    .max(2.0 * radius);
                let rim = radius * (1.0 + IO_PATERA_OUTLINE);
                let lobes = (0..n)
                    .map(|j| {
                        let d = base + 4 * j as u64;
                        let age_yr = self.log_uniform(idx, d + 3, age_range);
                        let (lobe, end) =
                            self.trace_lobe(idx, d, center, rim, reach, &ground, &standing, age_yr);
                        ends[end] += 1;
                        lobe_lengths.0 += 1;
                        lobe_lengths.1 += (lobe.points.len() - 1) as f64 * lobe.spacing;
                        lobe
                    })
                    .collect();
                self.flows.push(FlowField {
                    center,
                    reach,
                    lobes,
                    birth_gyr: t1,
                    death_gyr: death,
                    id: idx + 10_000_000 + 30_000_000 * g as u64,
                });
            }
        }
        if std::env::var_os("EXOGEN_DEBUG").is_some() && lobe_lengths.0 > 0 {
            eprintln!(
                "io flows at {t1:.3} Gyr: {} lobes, mean length {:.0} km; ended at their length {}, at a mountain {}, ponded {}",
                lobe_lengths.0,
                lobe_lengths.1 / lobe_lengths.0 as f64 / 1e3,
                ends[0],
                ends[1],
                ends[2]
            );
        }
    }

    /// The surface at the end of a step: the base terrain, the mountains standing, the
    /// paterae and the flows, all on the coarse grid. Only the steps whose surface the
    /// loop keeps as a keyframe are laid; the rest keep the timeline going. The plumes'
    /// fallout is the finish's (F84), since the loop's surfaces are not lit.
    pub fn step(&mut self, layers: &mut LayerSet, t0: f64, t1: f64, lay: bool) {
        match &self.base {
            None => {
                self.dt = t1 - t0;
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
        self.populate(t1, lay, &layers.bedrock);
        if !lay {
            return;
        }
        let grid = self.grid;
        let body_radius = self.body_radius;
        let cell = grid.delta() * body_radius;
        Self::clear_surface(layers, t1);
        for m in self.mountains.iter().filter(|m| m.alive(t1)) {
            let noise = GradientNoise::new(m.id ^ 0x0B10_C000);
            self.stamper.gather(
                &grid,
                body_radius,
                m.center,
                m.radius * m.aspect.sqrt() * (1.0 + IO_MOUNTAIN_OUTLINE),
            );
            let bedrock = layers.bedrock.as_mut_slice();
            let units = layers.unit_id.as_mut_slice();
            let cover = layers.cover.as_mut_slice();
            for &(idx, r, phi) in &self.stamper.footprint {
                let h = block_profile(r, phi, m, &noise);
                if h > 0.0 {
                    bedrock[idx] += h as f32;
                    units[idx] = UNIT_MOUNTAIN;
                    cover[idx] = 0.0;
                }
            }
        }
        for p in self.paterae.iter().filter(|p| p.alive(t1)) {
            let noise = GradientNoise::new(p.id ^ 0x9A7E_2A00);
            let wall = (IO_PATERA_WALL_FRACTION * p.radius).max(cell);
            self.stamper.gather(
                &grid,
                body_radius,
                p.center,
                p.radius * (1.0 + IO_PATERA_OUTLINE) + 3.0 * wall,
            );
            let bedrock = layers.bedrock.as_mut_slice();
            let units = layers.unit_id.as_mut_slice();
            let deep = layers.deep.as_mut_slice();
            let melt = layers.melt.as_mut_slice();
            let cover = layers.cover.as_mut_slice();
            for &(idx, r, phi) in &self.stamper.footprint {
                let edge = patera_radius(phi, p, &noise);
                bedrock[idx] -= pit_profile(r, edge, p.depth, wall) as f32;
                if r < edge {
                    units[idx] = UNIT_PATERA_FLOOR;
                    deep[idx] = 1.0;
                    melt[idx] = p.floor_age_yr as f32;
                    cover[idx] = 1.0;
                }
            }
        }
        let flows: Vec<FlowField> = self.flows.iter().filter(|f| f.alive(t1)).cloned().collect();
        for f in &flows {
            let out = flow_compute(layers, &grid, body_radius, f, &mut self.stamper);
            out.apply(layers);
        }
    }
}

/// Whether the body is in the Io regime: rock, heated past the threshold.
pub fn is_io_regime(graph: &ParamGraph) -> bool {
    graph.str("body.composition_class") == "rock"
        && graph.f64("interior.tidal_flux") > IO_REGIME_TIDAL_FLUX
}

/// Finishing pass of a patera: the pit's fine profile minus the version the loop
/// carried, and the floor's unit, its lava's coverage and the years since it was laid.
pub fn patera_compute(
    layers: &LayerSet,
    grid: &CubeSphere,
    coarse_cell: f64,
    body_radius: f64,
    p: &Patera,
    stamper: &mut Stamper,
) -> StampOutput {
    let bedrock = layers.bedrock.as_slice();
    let cell = grid.delta() * body_radius;
    let noise = GradientNoise::new(p.id ^ 0x9A7E_2A00);
    let wall_fine = (IO_PATERA_WALL_FRACTION * p.radius).max(cell);
    let wall_coarse = (IO_PATERA_WALL_FRACTION * p.radius).max(coarse_cell);
    let smear = (0.5 * coarse_cell).hypot(wall_coarse);
    stamper.gather(
        grid,
        body_radius,
        p.center,
        p.radius * (1.0 + IO_PATERA_OUTLINE) + 3.0 * smear,
    );
    let time = layers.surface_age.as_slice();
    let cells = stamper
        .footprint
        .iter()
        .map(|&(idx, r, phi)| {
            let edge = patera_radius(phi, p, &noise);
            let fine = pit_profile(r, edge, p.depth, wall_fine);
            let coarse = pit_profile(r, edge, p.depth, smear);
            let mut write = CellWrite {
                idx: idx as u32,
                bedrock: (bedrock[idx] as f64 - (fine - coarse)) as f32,
                surface: false,
                age: 0.0,
                unit: 0,
                deep: 0.0,
                melt: 0.0,
                cover: 0.0,
            };
            if r < edge {
                write.surface = true;
                write.age = time[idx];
                write.unit = UNIT_PATERA_FLOOR;
                write.deep = 1.0;
                write.melt = p.floor_age_yr as f32;
                write.cover = 1.0;
            }
            write
        })
        .collect();
    StampOutput {
        cells,
        composition: None,
        stain: None,
    }
}

/// Finishing pass of a mountain: ridged texture over the block.
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
    stamper.gather(
        grid,
        body_radius,
        m.center,
        m.radius * m.aspect.sqrt() * (1.0 + IO_MOUNTAIN_OUTLINE),
    );
    let outline = GradientNoise::new(m.id ^ 0x0B10_C000);
    let noise = GradientNoise::new(m.id ^ 0x3B10_C400);
    let f = body_radius / IO_MOUNTAIN_TEXTURE_LENGTH;
    let cells = stamper
        .footprint
        .iter()
        .map(|&(idx, r, phi)| {
            let d = grid.center(grid.cell(idx));
            let x = r / block_radius(phi, m, &outline);
            let block = block_profile(r, phi, m, &outline);
            let ridges = noise.ridged([d[0] * f, d[1] * f, d[2] * f], 1.0, 3, 2.0, 0.5) - 0.5;
            let texture = IO_MOUNTAIN_TEXTURE * m.height * ridges * (1.0 - x.powi(4)).max(0.0);
            let mut write = CellWrite {
                idx: idx as u32,
                bedrock: (bedrock[idx] as f64 + texture) as f32,
                surface: false,
                age: 0.0,
                unit: 0,
                deep: 0.0,
                melt: 0.0,
                cover: 0.0,
            };
            if block > 0.0 {
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
    }
}

/// A flow field over the ground: each lobe a tube about its centre line, its margin
/// lobate at the width's scale and soft over a share of the width; the strongest
/// coverage wins where lobes cross, with the years since the youngest lobe that covers
/// the cell well was laid, and a little relief goes with it. What falls on the lava
/// after coats it (F84).
pub fn flow_compute(
    layers: &LayerSet,
    grid: &CubeSphere,
    body_radius: f64,
    f: &FlowField,
    stamper: &mut Stamper,
) -> StampOutput {
    let bedrock = layers.bedrock.as_slice();
    let units = layers.unit_id.as_slice();
    let ages = layers.surface_age.as_slice();
    let laid = layers.melt.as_slice();
    let covered = layers.cover.as_slice();
    let cell = grid.delta() * body_radius;
    let noise = GradientNoise::new(f.id ^ 0xF10_0000);
    // Per cell: the strongest coverage, the age of the lobe giving it, and the age of
    // the youngest lobe covering at least half the cell.
    let mut best: HashMap<usize, (f64, f64, f64)> = HashMap::new();
    for lobe in &f.lobes {
        for (i, p) in lobe.points.iter().enumerate() {
            let hw = lobe.half_width * lobe.widths[i] as f64;
            let edge = (IO_FLOW_EDGE * hw).max(1.5 * cell);
            stamper.gather(grid, body_radius, *p, hw * (1.0 + IO_FLOW_LOBING) + edge);
            for &(idx, r, _) in &stamper.footprint {
                let d = grid.center(grid.cell(idx));
                let lobing =
                    1.0 + IO_FLOW_LOBING * noise.fbm(d, body_radius / (1.5 * hw), 2, 2.0, 0.5);
                let c = ((hw * lobing - r) / edge).clamp(0.0, 1.0);
                if c <= 0.0 {
                    continue;
                }
                let e = best
                    .entry(idx)
                    .or_insert((0.0, f64::INFINITY, f64::INFINITY));
                if c > e.0 {
                    e.0 = c;
                    e.1 = lobe.age_yr;
                }
                if c >= 0.5 {
                    e.2 = e.2.min(lobe.age_yr);
                }
            }
        }
    }
    let mut cells: Vec<CellWrite> = best
        .into_iter()
        .filter_map(|(idx, (cover, strongest, youngest))| {
            if units[idx] == UNIT_MOUNTAIN {
                return None;
            }
            // The flow units within a lobe overlap and differ: the coverage varies.
            let d = grid.center(grid.cell(idx));
            let texture = 1.0
                - IO_FLOW_TEXTURE
                    * 0.5
                    * (1.0 + noise.fbm(d, body_radius / IO_FLOW_TEXTURE_LENGTH, 2, 2.0, 0.5));
            let cover = cover * texture;
            if cover <= 0.0 {
                return None;
            }
            let age_yr = if youngest.is_finite() {
                youngest
            } else {
                strongest
            };
            // Over lava already laid, the younger wins the age and the stronger the
            // coverage.
            let (cover, age_yr) = if matches!(units[idx], UNIT_IO_FLOW | UNIT_PATERA_FLOOR) {
                (
                    (cover as f32).max(covered[idx]),
                    if cover >= 0.5 {
                        (age_yr as f32).min(laid[idx])
                    } else {
                        laid[idx]
                    },
                )
            } else {
                (cover as f32, age_yr as f32)
            };
            Some(CellWrite {
                idx: idx as u32,
                bedrock: bedrock[idx] + (IO_FLOW_THICKNESS as f32) * cover,
                surface: true,
                age: ages[idx],
                unit: UNIT_IO_FLOW,
                deep: cover,
                melt: age_yr,
                cover,
            })
        })
        .collect();
    cells.sort_by_key(|c| c.idx);
    StampOutput {
        cells,
        composition: None,
        stain: None,
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
        let m = Mountain {
            center: [1.0, 0.0, 0.0],
            radius: 60e3,
            height: 10e3,
            tilt_azimuth: 0.0,
            tilt: 0.4,
            aspect: 1.0,
            birth_gyr: 0.0,
            death_gyr: 1.0,
            id: 1,
        };
        let noise = GradientNoise::new(1);
        let top = block_profile(0.0, 0.0, &m, &noise);
        assert!((top - 8e3).abs() < 1.0, "{top}");
        // Tilted down towards azimuth 0: lower there than opposite.
        let near = block_profile(0.5 * m.radius, 0.0, &m, &noise);
        let far = block_profile(0.5 * m.radius, std::f64::consts::PI, &m, &noise);
        assert!(near < far, "{near} vs {far}");
        assert!(block_profile(2.0 * m.radius, 0.3, &m, &noise) == 0.0);
        // The ring peaks at the range and ends more steeply outside than it rises
        // inside; the pyroclasts fall within.
        let total = |x: f64| {
            let (ring, fill) = deposit_profile(x);
            ring + 0.3 * fill
        };
        assert!(total(1.0) > total(0.5) && total(0.5) > total(1.3));
        assert!(deposit_profile(1.15).0 < deposit_profile(0.85).0);
        assert!(deposit_profile(1.1).1 == 0.0 && deposit_profile(0.2).1 == 1.0);
        let moment = profile_moment(0.3);
        assert!(moment > 0.2 && moment < 0.6, "{moment}");
    }

    #[test]
    fn a_patera_outline_is_a_noisy_polygon() {
        let p = Patera {
            center: [1.0, 0.0, 0.0],
            radius: 30e3,
            depth: 1e3,
            active: false,
            floor_age_yr: 100.0,
            pele: false,
            sides: 6,
            rotation: 0.0,
            birth_gyr: 0.0,
            death_gyr: 1.0,
            id: 3,
        };
        let noise = GradientNoise::new(3);
        let mut lo = f64::INFINITY;
        let mut hi = 0.0f64;
        for k in 0..360 {
            let r = patera_radius((k as f64).to_radians(), &p, &noise);
            lo = lo.min(r);
            hi = hi.max(r);
        }
        assert!(lo > 0.6 * p.radius && hi < 1.4 * p.radius, "{lo} {hi}");
        assert!(hi - lo > 0.05 * p.radius);
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
        // Lava hides under its coat, slower than a smooth plain.
        let mut lava = Top {
            lava: 1.0,
            ..Top::default()
        };
        lava.age(1.0, red_fall, warm);
        assert!(lava.lava > 0.8, "{}", lava.lava);
        lava.age(50.0, red_fall, warm);
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
            frost: [3.0; 6],
            lag: 2.0,
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
    }

    /// The history of a patera's plumes (F84): an active Pele-class patera with a plume
    /// erupting now, the rest at their class's duty, every eruption within the history.
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
                sides: 6,
                rotation: 0.0,
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
        let prometheus = duty(0, 380.0);
        let want = IO_PATERA_ACTIVE_FRACTION;
        assert!(
            (pele / (want * IO_PLUME_FRACTION_PELE) - 1.0).abs() < 0.5,
            "Pele-class duty {pele}"
        );
        assert!(
            (prometheus / (want * IO_PLUME_FRACTION_PROMETHEUS) - 1.0).abs() < 0.3,
            "Prometheus-class duty {prometheus}"
        );
        assert!(erupting[0] + erupting[1] > 3, "{erupting:?}");
        // Pele-class rings are far wider than Prometheus-class ones.
        let mean = |pele: bool| {
            let r: Vec<f64> = episodes
                .iter()
                .filter(|e| e.pele == pele)
                .map(|e| e.radius)
                .collect();
            r.iter().sum::<f64>() / r.len() as f64
        };
        assert!(mean(true) > 2.5 * mean(false));
    }
}
