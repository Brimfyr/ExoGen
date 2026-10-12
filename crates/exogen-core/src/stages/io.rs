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
//! active ones and a few quiet ones feed compound flow fields, lobes of different ages
//! from black fresh lava to older flows going brown under the sulfur rain; the gas-rich
//! Pele-class centres and some of the Prometheus-class ones throw plumes whose fallout
//! settles in rings at the ballistic range of the eruption speed, red sulfur from the
//! hot gas-rich plumes and white sulfur dioxide frost from the rest, and the fallout of
//! a quiet centre's last plumes fades as it is buried. Mountains are tilted crustal
//! blocks raised by the subsidence of the crust under the burial, standing a few hundred
//! million years.
//!
//! The colour is not painted. Each cell carries what the eruptions left on it, the fresh
//! silicate exposure in `deep`, the red sulfur fallout in `melt` and the frost fallout in
//! `cover` (on a flow or a floor, `cover` is the lava's own coverage), and the materials
//! stage lights those by the temperature the latitude's insolation gives, which decides
//! where frost survives the day and where radiation keeps sulfur red.

use crate::constants::*;
use crate::crater::{CellWrite, StampOutput, Stamper};
use crate::fields::LayerSet;
use crate::geom::{angle_between, azimuth_from, local_basis};
use crate::grid::CubeSphere;
use crate::noise::GradientNoise;
use crate::params::ParamGraph;
use crate::rng::Rng;
use crate::units::*;

pub const VERSION: u64 = 4;

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
    /// Years since the floor last held a lake: it goes from black silicate to a sulfur
    /// coat.
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

/// One lobe of a flow field: a sector of the field's reach, laid down some years ago.
#[derive(Clone, Debug)]
pub struct Lobe {
    pub azimuth: f64,
    pub spread: f64,
    pub reach_fraction: f64,
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

/// A plume's fallout round a patera: a ring at the ballistic range, red sulfur from a
/// Pele-class plume and frost from a Prometheus-class one, at a density that fades as
/// the fallout is buried after the plume dies.
#[derive(Clone, Debug)]
pub struct PlumeRing {
    pub center: [f64; 3],
    pub radius: f64,
    pub pele: bool,
    /// Densities of the red sulfur and the frost in the fallout, 0 to 1, each faded by
    /// its own time since the plume died.
    pub red: f64,
    pub white: f64,
    pub birth_gyr: f64,
    pub death_gyr: f64,
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
        1.1 * self.reach / body_radius
    }
    pub fn alive(&self, t: f64) -> bool {
        self.birth_gyr - IO_ALIVE_TOLERANCE_GYR <= t && t < self.death_gyr - IO_ALIVE_TOLERANCE_GYR
    }
}
impl PlumeRing {
    pub fn cap_radius(&self, body_radius: f64) -> f64 {
        (1.0 + 3.0 * IO_RING_WIDTH) * self.radius / body_radius
    }
    pub fn alive(&self, t: f64) -> bool {
        self.birth_gyr - IO_ALIVE_TOLERANCE_GYR <= t && t < self.death_gyr - IO_ALIVE_TOLERANCE_GYR
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
    let rim = (1.0 - x.powi(8)).max(0.0);
    m.height * tilt * rim
}

/// A flow field's coverage and freshness at a point: the strongest lobe reaching it,
/// and the freshest lava among the lobes that reach it.
fn flow_at(d: [f64; 3], f: &FlowField, body_radius: f64, noise: &GradientNoise) -> (f64, f64) {
    let r = angle_between(f.center, d) * body_radius;
    if r > f.reach {
        return (0.0, 0.0);
    }
    let (east, north) = local_basis(f.center);
    let phi = azimuth_from(f.center, east, north, d);
    // The vent's own pool at the centre; the lobes reach out of a connected body.
    let mut cover: f64 = ((IO_FLOW_VENT * f.reach - r) / (IO_FLOW_EDGE * f.reach)).clamp(0.0, 1.0);
    let mut fresh: f64 = 0.0;
    for (j, lobe) in f.lobes.iter().enumerate() {
        let mut da = (phi - lobe.azimuth).rem_euclid(2.0 * std::f64::consts::PI);
        if da > std::f64::consts::PI {
            da -= 2.0 * std::f64::consts::PI;
        }
        let sector = (-(da / lobe.spread).powi(2)).exp();
        let lobing = 1.0
            + IO_FLOW_LOBING
                * noise.fbm(
                    [
                        phi.cos() * 2.0 + j as f64 * 3.1,
                        phi.sin() * 2.0,
                        r / f.reach * 2.0,
                    ],
                    1.0,
                    2,
                    2.0,
                    0.5,
                );
        let extent =
            f.reach * (IO_FLOW_BODY + (1.0 - IO_FLOW_BODY) * lobe.reach_fraction * sector) * lobing;
        // An old lobe is coated over by the plains' fallout and fades.
        let visible = (-lobe.age_yr / IO_FLOW_COAT_YR).exp();
        let c = ((extent - r) / (IO_FLOW_EDGE * f.reach)).clamp(0.0, 1.0) * visible;
        if c > 0.0 {
            cover = cover.max(c);
            fresh = fresh.max(c * (-lobe.age_yr / IO_FLOW_DARK_YR).exp());
        }
    }
    // The flow units within a lobe overlap and differ: the surface varies along it.
    let texture = 1.0
        - IO_FLOW_TEXTURE
            * 0.5
            * (1.0 + noise.fbm(d, body_radius / IO_FLOW_TEXTURE_LENGTH, 2, 2.0, 0.5));
    (cover * texture, fresh)
}

/// Fallout density of a plume at a fraction of its ring radius: the ring where the
/// ballistic shells land, over a thinner fall inside it.
pub fn ring_profile(x: f64, inner: f64) -> f64 {
    let ring = (-((x - 1.0) / IO_RING_WIDTH).powi(2)).exp();
    let fill = inner * ((1.0 - x) / 0.5).clamp(0.0, 1.0);
    (ring + fill).min(1.0)
}

/// The Io regime's features through the loop: a steady state, laid down afresh at every
/// step on the base terrain, with the full timeline kept so any epoch can be finished.
pub struct IoProcess {
    grid: CubeSphere,
    body_radius: f64,
    gravity: f64,
    /// The day's peak temperature at the subsolar point, K, from the insolation.
    subsolar_k: f64,
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
    pub rings: Vec<PlumeRing>,
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
            gravity: graph.f64("body.gravity"),
            subsolar_k: {
                let absorbed =
                    graph.f64("orbit.insolation") * (1.0 - graph.f64("atmosphere.bond_albedo"));
                (absorbed.max(0.0) / SIGMA_SB).powf(0.25)
            },
            age: graph.f64("body.age"),
            dt: 0.0,
            rng: Rng::new(graph.f64("seed.volcanism") as u64 ^ 0x10_5EED_0001),
            n_paterae: graph.f64("volcanism.paterae_count").round().max(0.0) as usize,
            n_mountains: graph.f64("tectonics.mountain_count").round().max(0.0) as usize,
            base: None,
            paterae: Vec::new(),
            mountains: Vec::new(),
            flows: Vec::new(),
            rings: Vec::new(),
            stamper: Stamper::new(&grid),
        })
    }

    fn log_uniform(&self, index: u64, draw: u64, range: (f64, f64)) -> f64 {
        range.0 * (range.1 / range.0).powf(self.rng.uniform(index, draw))
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
                    center: self.rng.unit_direction(k, 2),
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

    /// The features alive through a step: paterae, flows and fallout are renewed at
    /// every kept step, since the surface is younger than the step.
    fn populate(&mut self, t1: f64, lay: bool) {
        // The time bin the step ends in is the random index, so the same surface comes
        // out at every quality.
        let step = (t1 / self.age * crate::stages::LOOP_TIME_BINS as f64).round() as u64;
        let death = t1 + self.dt.max(1e-6);
        // Paterae, with their flows and fallout, for the steps whose surface is kept.
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
                if self.rng.uniform(idx, 2 + attempt) < w {
                    break;
                }
                center = self.rng.unit_direction(idx, 20 + attempt);
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
            // patera's older generations, lobes of all ages fading into the plains as
            // the fallout coats them (F56).
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
                    (IO_FLOW_DARK_YR, a_hi)
                };
                let lobes = (0..n)
                    .map(|j| {
                        let d = base + 4 * j as u64;
                        let (sp_lo, sp_hi) = IO_FLOW_SPREAD;
                        Lobe {
                            azimuth: 2.0 * std::f64::consts::PI * self.rng.uniform(idx, d),
                            spread: sp_lo + (sp_hi - sp_lo) * self.rng.uniform(idx, d + 1),
                            reach_fraction: 0.3 + 0.7 * self.rng.uniform(idx, d + 2),
                            age_yr: self.log_uniform(idx, d + 3, age_range),
                        }
                    })
                    .collect();
                self.flows.push(FlowField {
                    center,
                    reach: self
                        .log_uniform(idx, 49 + 90 * g as u64, IO_FLOW_REACH)
                        .max(2.0 * radius),
                    lobes,
                    birth_gyr: t1,
                    death_gyr: death,
                    id: idx + 10_000_000 + 30_000_000 * g as u64,
                });
            }
            // Fallout: a plume now, or the fading fallout of the last ones, each
            // generation older than the last (F56).
            let plume_fraction = if pele {
                IO_PLUME_FRACTION_PELE
            } else {
                IO_PLUME_FRACTION_PROMETHEUS
            };
            for g in 0..IO_PLUME_GENERATIONS {
                let plume_age_yr = if g == 0 && active && self.rng.uniform(idx, 50) < plume_fraction
                {
                    0.0
                } else {
                    self.log_uniform(idx, 51 + 10 * g as u64, IO_PLUME_AGE_YR)
                };
                // Red sulfur fades in a couple of years; frost lasts until it is buried,
                // longer where the cold keeps it.
                let cos_lat = (1.0 - center[2] * center[2]).max(0.0).sqrt();
                let t_peak = self.subsolar_k * cos_lat.max(0.02).powf(0.25);
                let persist = 1.0 / (1.0 + ((t_peak - IO_FROST_HALF_K) / IO_FROST_WIDTH_K).exp());
                let frost_fade = IO_FALLOUT_FADE_YR * (1.0 + IO_FROST_COLD_FACTOR * persist);
                let (red_share, white_share) = if pele {
                    (1.0, IO_PELE_FROST)
                } else {
                    (IO_PROMETHEUS_RED, 1.0)
                };
                let red = red_share * (-plume_age_yr / IO_RED_FADE_YR).exp();
                let white = white_share * (-plume_age_yr / frost_fade).exp();
                if red.max(white) > 0.05 {
                    let speed = if pele {
                        IO_PLUME_SPEED_PELE
                    } else {
                        IO_PLUME_SPEED_PROMETHEUS
                    } * (IO_PLUME_SPEED_SCATTER * self.rng.normal(idx, 52)).exp();
                    self.rings.push(PlumeRing {
                        center,
                        radius: speed * speed / self.gravity.max(0.01),
                        pele,
                        red,
                        white,
                        birth_gyr: t1,
                        death_gyr: death,
                        id: idx + 20_000_000 + 30_000_000 * g as u64,
                    });
                }
            }
        }
    }

    /// The surface at the end of a step: the base terrain, the mountains standing, the
    /// paterae, the fallout and the flows, all on the coarse grid. Only the steps whose
    /// surface the loop keeps as a keyframe are laid; the rest keep the timeline going.
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
        self.populate(t1, lay);
        if !lay {
            return;
        }
        let grid = self.grid;
        let body_radius = self.body_radius;
        let cell = grid.delta() * body_radius;
        let now = (t1 - IO_SURFACE_AGE_GYR) as f32;
        {
            let ages = layers.surface_age.as_mut_slice();
            let units = layers.unit_id.as_mut_slice();
            let deep = layers.deep.as_mut_slice();
            let melt = layers.melt.as_mut_slice();
            let cover = layers.cover.as_mut_slice();
            for k in 0..ages.len() {
                ages[k] = now;
                units[k] = UNIT_PRIMORDIAL;
                deep[k] = 0.0;
                melt[k] = 0.0;
                cover[k] = 0.0;
            }
        }
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
            let cover = layers.cover.as_mut_slice();
            let fresh = floor_freshness(p) as f32;
            for &(idx, r, phi) in &self.stamper.footprint {
                let edge = patera_radius(phi, p, &noise);
                bedrock[idx] -= pit_profile(r, edge, p.depth, wall) as f32;
                if r < edge {
                    units[idx] = if p.active {
                        UNIT_IO_FLOW
                    } else {
                        UNIT_PATERA_FLOOR
                    };
                    deep[idx] = fresh;
                    cover[idx] = 1.0;
                }
            }
        }
        let rings: Vec<PlumeRing> = self.rings.iter().filter(|r| r.alive(t1)).cloned().collect();
        for r in &rings {
            let out = ring_compute(layers, &grid, body_radius, r, &mut self.stamper);
            out.apply(layers);
        }
        let flows: Vec<FlowField> = self.flows.iter().filter(|f| f.alive(t1)).cloned().collect();
        for f in &flows {
            let out = flow_compute(layers, &grid, body_radius, f, &mut self.stamper);
            out.apply(layers);
        }
    }
}

/// How much of a patera's floor is still fresh silicate.
pub fn floor_freshness(p: &Patera) -> f64 {
    if p.active {
        1.0
    } else {
        (-p.floor_age_yr / IO_FLOOR_COAT_YR).exp()
    }
}

/// Whether the body is in the Io regime: rock, heated past the threshold.
pub fn is_io_regime(graph: &ParamGraph) -> bool {
    graph.str("body.composition_class") == "rock"
        && graph.f64("interior.tidal_flux") > IO_REGIME_TIDAL_FLUX
}

/// Finishing pass of a patera: the pit's fine profile minus the version the loop
/// carried, and the floor's unit and freshness.
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
    let melt = layers.melt.as_slice();
    let fresh = floor_freshness(p) as f32;
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
                write.unit = if p.active {
                    UNIT_IO_FLOW
                } else {
                    UNIT_PATERA_FLOOR
                };
                write.deep = fresh;
                write.melt = melt[idx] * (1.0 - fresh);
                write.cover = 1.0;
            }
            write
        })
        .collect();
    StampOutput {
        cells,
        composition: None,
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
    }
}

/// A flow field over the ground: its lobes' coverage and freshness, a little relief;
/// over an earlier flow the stronger coverage and the fresher lava win, and the lava
/// buries the fallout under it.
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
    let deep = layers.deep.as_slice();
    let melt = layers.melt.as_slice();
    let covered = layers.cover.as_slice();
    stamper.gather(grid, body_radius, f.center, f.reach);
    let noise = GradientNoise::new(f.id ^ 0xF10_0000);
    let cells = stamper
        .footprint
        .iter()
        .filter_map(|&(idx, _, _)| {
            if units[idx] == UNIT_MOUNTAIN {
                return None;
            }
            let d = grid.center(grid.cell(idx));
            let (cover, fresh) = flow_at(d, f, body_radius, &noise);
            if cover <= 0.0 {
                return None;
            }
            let (cover, fresh) = if matches!(units[idx], UNIT_IO_FLOW | UNIT_PATERA_FLOOR) {
                (
                    (cover as f32).max(covered[idx]),
                    (fresh as f32).max(deep[idx]),
                )
            } else {
                (cover as f32, fresh as f32)
            };
            Some(CellWrite {
                idx: idx as u32,
                bedrock: bedrock[idx] + (IO_FLOW_THICKNESS as f32) * cover,
                surface: true,
                age: ages[idx],
                unit: UNIT_IO_FLOW,
                deep: fresh,
                melt: melt[idx] * (1.0 - cover),
                cover,
            })
        })
        .collect();
    StampOutput {
        cells,
        composition: None,
    }
}

/// A plume's fallout over the ground: red sulfur or frost by the plume's class, at the
/// ring's density and the fallout's, mottled, laid over the plains and the mountains
/// and over cooled lava, where it joins what earlier plumes left.
pub fn ring_compute(
    layers: &LayerSet,
    grid: &CubeSphere,
    body_radius: f64,
    ring: &PlumeRing,
    stamper: &mut Stamper,
) -> StampOutput {
    let bedrock = layers.bedrock.as_slice();
    let units = layers.unit_id.as_slice();
    let ages = layers.surface_age.as_slice();
    let deep = layers.deep.as_slice();
    let melt = layers.melt.as_slice();
    let cover = layers.cover.as_slice();
    stamper.gather(
        grid,
        body_radius,
        ring.center,
        (1.0 + 3.0 * IO_RING_WIDTH) * ring.radius,
    );
    let noise = GradientNoise::new(ring.id ^ 0x21_0600);
    let mottle_frequency = body_radius / (0.15 * ring.radius);
    let inner = if ring.pele {
        IO_RING_INNER_PELE
    } else {
        IO_RING_INNER_PROMETHEUS
    };
    let cells = stamper
        .footprint
        .iter()
        .filter_map(|&(idx, r, phi)| {
            // Hot lava takes no fallout.
            if deep[idx] > 0.5 {
                return None;
            }
            let d = grid.center(grid.cell(idx));
            let edge =
                IO_RING_EDGE_WOBBLE * noise.fbm([phi.cos(), phi.sin(), 0.7], 3.0, 2, 2.0, 0.5);
            let x = r / ring.radius / (1.0 + edge);
            let mottle = 1.0
                - IO_DEPOSIT_MOTTLING * (0.5 + 0.5 * noise.fbm(d, mottle_frequency, 2, 2.0, 0.5));
            let fall = ring_profile(x, inner) * mottle;
            if fall < 0.02 {
                return None;
            }
            let red = (ring.red * fall) as f32;
            let white = (ring.white * fall) as f32;
            // On a flow or a floor the coverage is the lava's own; frost goes on plains
            // and mountains.
            let lava = matches!(units[idx], UNIT_IO_FLOW | UNIT_PATERA_FLOOR);
            let new_cover = if lava {
                cover[idx]
            } else {
                cover[idx].max(white)
            };
            let new_melt = melt[idx].max(red);
            if new_cover <= cover[idx] && new_melt <= melt[idx] {
                return None;
            }
            Some(CellWrite {
                idx: idx as u32,
                bedrock: bedrock[idx],
                surface: true,
                age: ages[idx],
                unit: units[idx],
                deep: deep[idx],
                melt: new_melt,
                cover: new_cover,
            })
        })
        .collect();
    StampOutput {
        cells,
        composition: None,
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
        // The ring peaks at the ring radius and the inner fall is thinner.
        assert!(ring_profile(1.0, 0.3) > ring_profile(0.5, 0.3));
        assert!(ring_profile(0.5, 0.3) > ring_profile(1.6, 0.3));
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
        assert!((floor_freshness(&p) - (-100.0f64 / IO_FLOOR_COAT_YR).exp()).abs() < 1e-9);
    }
}
