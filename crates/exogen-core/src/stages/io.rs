// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! The Io regime (PLAN.md 4.4): a rocky body heated past about a watt per square metre.
//!
//! Resurfacing outruns impacts by orders of magnitude, so the surface is always young
//! and carries no craters. What it carries instead: hundreds of paterae, flat-floored
//! volcanic depressions tens to two hundred kilometres across, concentrated at low
//! latitudes where the tidal heat comes out; dark silicate flow fields spreading from
//! the active ones; the ring deposits of plumes, red sulfur near the vent and white
//! sulfur dioxide frost further out; and isolated mountains, tilted crustal blocks up
//! to seventeen kilometres high raised by the subsidence of the crust under the burial.
//! Everything here is laid down in the last step of the loop, since the oldest surface
//! is a fraction of a million years old, and refined at the finish.

use crate::constants::*;
use crate::crater::{CellWrite, StampOutput, Stamper};
use crate::fields::LayerSet;
use crate::geom::{angle_between, azimuth_from, local_basis};
use crate::grid::CubeSphere;
use crate::noise::GradientNoise;
use crate::params::ParamGraph;
use crate::rng::Rng;
use crate::units::*;

pub const VERSION: u64 = 2;

/// A patera: a volcanic depression with a flat floor, active when its floor is a lake of
/// fresh silicate lava.
#[derive(Clone, Debug)]
pub struct Patera {
    pub center: [f64; 3],
    pub radius: f64,
    pub depth: f64,
    pub active: bool,
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
    pub id: u64,
}

/// A flow field spreading from a patera.
#[derive(Clone, Debug)]
pub struct FlowField {
    pub center: [f64; 3],
    pub reach: f64,
    pub id: u64,
}

/// A plume's ring deposits round a patera.
#[derive(Clone, Debug)]
pub struct PlumeRing {
    pub center: [f64; 3],
    pub radius: f64,
    /// The one Pele-class plume, whose ring is red; the rest are mostly frost.
    pub giant: bool,
    pub id: u64,
}

impl Patera {
    pub fn cap_radius(&self, body_radius: f64) -> f64 {
        1.3 * self.radius / body_radius
    }
}
impl Mountain {
    pub fn cap_radius(&self, body_radius: f64) -> f64 {
        1.2 * self.radius / body_radius
    }
}
impl FlowField {
    pub fn cap_radius(&self, body_radius: f64) -> f64 {
        1.1 * self.reach / body_radius
    }
}
impl PlumeRing {
    pub fn cap_radius(&self, body_radius: f64) -> f64 {
        1.5 * self.radius / body_radius
    }
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
    let outline =
        1.0 + IO_MOUNTAIN_OUTLINE * noise.fbm([phi.cos(), phi.sin(), 0.5], 2.5, 2, 2.0, 0.5);
    m.radius * ellipse * outline
}

/// Height of a tilted block at a distance and azimuth from its centre.
pub fn block_profile(r: f64, phi: f64, m: &Mountain, noise: &GradientNoise) -> f64 {
    let x = r / block_radius(phi, m, noise);
    if x >= 1.0 {
        return 0.0;
    }
    let plateau = (1.0 - x.powi(4)).powf(0.5);
    let tilt = 1.0 - m.tilt * 0.5 * (1.0 + x * (phi - m.tilt_azimuth).cos());
    m.height * plateau * tilt
}

/// Coverage of a lobate flow field at a point, 1 inside, fading to 0 across a soft
/// margin at a noisy radius round the source.
fn flow_cover(d: [f64; 3], f: &FlowField, body_radius: f64, noise: &GradientNoise) -> f64 {
    let r = angle_between(f.center, d) * body_radius;
    if r > f.reach {
        return 0.0;
    }
    let (east, north) = local_basis(f.center);
    let phi = azimuth_from(f.center, east, north, d);
    let lobes = noise.fbm(
        [phi.cos() * 1.7, phi.sin() * 1.7, r / f.reach * 2.0],
        1.0,
        2,
        2.0,
        0.5,
    );
    let edge = 0.55 + IO_FLOW_LOBING * lobes;
    ((edge - r / f.reach) / IO_DEPOSIT_EDGE).clamp(0.0, 1.0)
}

/// The Io regime's features, laid down in the last loop step.
pub struct IoProcess {
    grid: CubeSphere,
    body_radius: f64,
    age: f64,
    done: bool,
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
        let body_radius = graph.f64("body.radius");
        let age = graph.f64("body.age");
        let rng = Rng::new(graph.f64("seed.volcanism") as u64 ^ 0x10_5EED_0001);
        let n_paterae = graph.f64("volcanism.paterae_count").round() as u64;
        let mut paterae = Vec::new();
        for k in 0..n_paterae {
            // Low latitudes favoured: draw until the latitude weight accepts.
            let mut center = rng.unit_direction(k, 1);
            for attempt in 0..12u64 {
                let lat = center[2].asin();
                let w = IO_PATERA_POLAR_WEIGHT + (1.0 - IO_PATERA_POLAR_WEIGHT) * lat.cos().powi(2);
                if rng.uniform(k, 2 + attempt) < w {
                    break;
                }
                center = rng.unit_direction(k, 20 + attempt);
            }
            let (lo, hi) = IO_PATERA_DIAMETER;
            let diameter = lo * (hi / lo).powf(rng.uniform(k, 40).powf(IO_PATERA_SIZE_SKEW));
            paterae.push(Patera {
                center,
                radius: 0.5 * diameter,
                depth: IO_PATERA_DEPTH * (0.5 + rng.uniform(k, 41)),
                active: rng.uniform(k, 42) < IO_PATERA_ACTIVE_FRACTION,
                id: k,
            });
        }
        let n_mountains = graph.f64("tectonics.mountain_count").round() as u64;
        let mut mountains = Vec::new();
        for k in 0..n_mountains {
            let (lo, hi) = IO_MOUNTAIN_DIAMETER;
            let diameter = lo * (hi / lo).powf(rng.uniform(5_000 + k, 1));
            let (h_lo, h_hi) = IO_MOUNTAIN_HEIGHT;
            mountains.push(Mountain {
                center: rng.unit_direction(5_000 + k, 2),
                radius: 0.5 * diameter,
                height: h_lo * (h_hi / h_lo).powf(rng.uniform(5_000 + k, 3)),
                tilt_azimuth: 2.0 * std::f64::consts::PI * rng.uniform(5_000 + k, 4),
                tilt: IO_MOUNTAIN_TILT * rng.uniform(5_000 + k, 5),
                aspect: 1.0 + (IO_MOUNTAIN_ASPECT - 1.0) * rng.uniform(5_000 + k, 6),
                id: 5_000 + k,
            });
        }
        // Flow fields from the active paterae and a few quiet ones; plume rings round a
        // share of the active ones.
        let mut flows = Vec::new();
        let mut rings: Vec<PlumeRing> = Vec::new();
        for p in &paterae {
            let draw = 9_000 + p.id;
            let flowing = p.active || rng.uniform(draw, 1) < IO_QUIET_FLOW_FRACTION;
            if flowing {
                let (lo, hi) = IO_FLOW_REACH;
                flows.push(FlowField {
                    center: p.center,
                    reach: (lo * (hi / lo).powf(rng.uniform(draw, 2))).max(2.0 * p.radius),
                    id: draw,
                });
            }
            if p.active && rng.uniform(draw, 3) < IO_PLUME_FRACTION {
                // The first plume drawn is the giant one.
                let giant = rings.is_empty();
                let (lo, hi) = if giant {
                    IO_GIANT_RING_RADIUS
                } else {
                    IO_PLUME_RING_RADIUS
                };
                rings.push(PlumeRing {
                    center: p.center,
                    radius: lo * (hi / lo).powf(rng.uniform(draw, 4)),
                    giant,
                    id: draw + 20_000,
                });
            }
        }
        Some(Self {
            grid,
            body_radius,
            age,
            done: false,
            paterae,
            mountains,
            flows,
            rings,
            stamper: Stamper::new(&grid),
        })
    }

    /// In the last step: the paterae and mountains on the coarse grid, and a surface
    /// everywhere as young as the resurfacing makes it.
    pub fn step(&mut self, layers: &mut LayerSet, _t0: f64, t1: f64) {
        if self.done || t1 < self.age - 1e-9 {
            return;
        }
        self.done = true;
        let grid = self.grid;
        let body_radius = self.body_radius;
        let cell = grid.delta() * body_radius;
        let now = (self.age - IO_SURFACE_AGE_GYR) as f32;
        {
            let ages = layers.surface_age.as_mut_slice();
            let units = layers.unit_id.as_mut_slice();
            for (a, u) in ages.iter_mut().zip(units.iter_mut()) {
                *a = now;
                *u = UNIT_PRIMORDIAL;
            }
        }
        for m in &self.mountains {
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
                    cover[idx] = 1.0;
                }
            }
        }
        for p in &self.paterae {
            let wall = (IO_PATERA_WALL_FRACTION * p.radius).max(cell);
            self.stamper
                .gather(&grid, body_radius, p.center, p.radius + 3.0 * wall);
            let bedrock = layers.bedrock.as_mut_slice();
            let units = layers.unit_id.as_mut_slice();
            let cover = layers.cover.as_mut_slice();
            for &(idx, r, _) in &self.stamper.footprint {
                let depth = pit_profile(r, p.radius, p.depth, wall);
                bedrock[idx] -= depth as f32;
                if r < p.radius {
                    units[idx] = if p.active {
                        UNIT_IO_FLOW
                    } else {
                        UNIT_PATERA_FLOOR
                    };
                    cover[idx] = 1.0;
                }
            }
        }
    }
}

/// Whether the body is in the Io regime: rock, heated past the threshold.
pub fn is_io_regime(graph: &ParamGraph) -> bool {
    graph.str("body.composition_class") == "rock"
        && graph.f64("interior.tidal_flux") > IO_REGIME_TIDAL_FLUX
}

/// Finishing pass of a patera: the walls sharpened from the loop's cell-wide ones, the
/// floor's unit.
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
    let wall_fine = (IO_PATERA_WALL_FRACTION * p.radius).max(cell);
    let wall_coarse = (IO_PATERA_WALL_FRACTION * p.radius).max(coarse_cell);
    let smear = (0.5 * coarse_cell).hypot(wall_coarse);
    stamper.gather(grid, body_radius, p.center, p.radius + 3.0 * smear);
    let time = layers.surface_age.as_slice();
    let cells = stamper
        .footprint
        .iter()
        .map(|&(idx, r, _)| {
            let fine = pit_profile(r, p.radius, p.depth, wall_fine);
            let coarse = pit_profile(r, p.radius, p.depth, smear);
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
            if r < p.radius + 0.5 * wall_fine {
                write.surface = true;
                write.age = time[idx];
                write.unit = if r < p.radius {
                    if p.active {
                        UNIT_IO_FLOW
                    } else {
                        UNIT_PATERA_FLOOR
                    }
                } else {
                    UNIT_PRIMORDIAL
                };
                write.cover = if r < p.radius { 1.0 } else { 0.0 };
            }
            write
        })
        .collect();
    StampOutput { cells }
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
                write.cover = 1.0;
            }
            write
        })
        .collect();
    StampOutput { cells }
}

/// Finishing pass of a flow field: dark fresh silicate over the plains, a little relief.
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
    stamper.gather(grid, body_radius, f.center, f.reach);
    let noise = GradientNoise::new(f.id ^ 0xF10_0000);
    let covered = layers.cover.as_slice();
    let cells = stamper
        .footprint
        .iter()
        .filter_map(|&(idx, _, _)| {
            let d = grid.center(grid.cell(idx));
            if units[idx] == UNIT_MOUNTAIN {
                return None;
            }
            let cover = flow_cover(d, f, body_radius, &noise) as f32;
            if cover <= 0.0 {
                return None;
            }
            // A flow over an earlier flow keeps the stronger coverage.
            let cover = if units[idx] == UNIT_IO_FLOW {
                cover.max(covered[idx])
            } else {
                cover
            };
            Some(CellWrite {
                idx: idx as u32,
                bedrock: bedrock[idx] + (IO_FLOW_THICKNESS as f32) * cover,
                surface: true,
                age: ages[idx],
                unit: UNIT_IO_FLOW,
                deep: 0.0,
                melt: 0.0,
                cover,
            })
        })
        .collect();
    StampOutput { cells }
}

/// Finishing pass of a plume's rings: red deposits inside the ring radius, white frost
/// beyond it, over whatever plains they fall on.
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
    let existing = layers.cover.as_slice();
    stamper.gather(
        grid,
        body_radius,
        ring.center,
        (IO_RING_WHITE_OUTER + 2.0 * IO_DEPOSIT_EDGE) * ring.radius,
    );
    let noise = GradientNoise::new(ring.id ^ 0x21_0600);
    let mottle_frequency = body_radius / (0.15 * ring.radius);
    let cells = stamper
        .footprint
        .iter()
        .filter_map(|&(idx, r, phi)| {
            // Over the plains, or over another ring's deposit where this one is denser.
            let over_deposit = matches!(units[idx], UNIT_PLUME_RED | UNIT_PLUME_WHITE);
            if units[idx] != UNIT_PRIMORDIAL && !over_deposit {
                return None;
            }
            let d = grid.center(grid.cell(idx));
            let x = r / ring.radius;
            let edge =
                IO_RING_EDGE_WOBBLE * noise.fbm([phi.cos(), phi.sin(), 0.7], 3.0, 2, 2.0, 0.5);
            let soft = |a: f64| (a / IO_DEPOSIT_EDGE).clamp(0.0, 1.0);
            // Red from the inner edge to the ring, white frost fading out beyond it, both
            // mottled.
            let red = soft(x - IO_RING_RED_INNER - edge) * soft(1.0 + edge - x);
            let white = soft(x - 1.0 - edge + IO_DEPOSIT_EDGE)
                * ((IO_RING_WHITE_OUTER + edge - x) / (IO_RING_WHITE_OUTER - 1.0)).clamp(0.0, 1.0);
            let mottle = 1.0
                - IO_DEPOSIT_MOTTLING * (0.5 + 0.5 * noise.fbm(d, mottle_frequency, 2, 2.0, 0.5));
            let (red, white) = if ring.giant {
                (red, white * IO_GIANT_RING_WHITE)
            } else {
                (red * IO_SMALL_RING_RED, white)
            };
            let (unit, cover) = if red >= white {
                (UNIT_PLUME_RED, red * mottle)
            } else {
                (UNIT_PLUME_WHITE, white * mottle)
            };
            if cover < 0.02 || (over_deposit && cover <= existing[idx] as f64) {
                return None;
            }
            Some(CellWrite {
                idx: idx as u32,
                bedrock: bedrock[idx],
                surface: true,
                age: ages[idx],
                unit,
                deep: 0.0,
                melt: 0.0,
                cover: cover as f32,
            })
        })
        .collect();
    StampOutput { cells }
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
    }
}
