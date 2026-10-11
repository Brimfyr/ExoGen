// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! What an ice shell does with its heat and its ocean's salts (PLAN.md 4.3 and Phase
//! 4): chaos where a thin shell melts through, cryovolcanic domes and salt-bright crater
//! floors on a shell over rock, and a polar plume's fallout of fresh ice.
//!
//! Chaos is a population of patches, their number from the share of the surface the
//! shell's heat disrupts and their sizes from a distribution, drawn in a heat field; the
//! matrix inside a patch drops and breaks into hummocks around rafts of the old surface,
//! and the ocean's salts stain it. On a shell over rock the ice magma that does not
//! flood extrudes as domes of brine-rich ice, which relax like any other load, and a
//! crater deep enough to reach the brines gets salts on its floor that darken with age.
//! A plume lays fresh ice down over the body, fastest near its pole: it buries lags and
//! salts before they can show and keeps the surface in frost.

use rayon::prelude::*;

use crate::constants::*;
use crate::crater::{CellWrite, CraterRecord, StampOutput, Stamper, dimensions};
use crate::fields::LayerSet;
use crate::geom::angle_between;
use crate::grid::CubeSphere;
use crate::noise::GradientNoise;
use crate::params::ParamGraph;
use crate::rng::Rng;
use crate::units::*;

pub const VERSION: u64 = 1;

/// A polar plume's fallout of fresh ice.
pub struct Fallout {
    pole: [f64; 3],
    /// The rate at the pole, m/Gyr.
    rate: f64,
    angle: f64,
}

impl Fallout {
    pub fn new(graph: &ParamGraph) -> Option<Self> {
        let rate = graph.f64("volcanism.plume_fallout_rate");
        if rate <= 0.0 {
            return None;
        }
        let pole = crate::stages::ice::ShellMap::new(graph).pole?;
        Some(Self {
            pole,
            rate,
            angle: PLUME_FALLOUT_ANGLE_DEG.to_radians(),
        })
    }

    /// The fallout rate at a direction, m/Gyr.
    pub fn rate_at(&self, d: [f64; 3]) -> f64 {
        let x = angle_between(self.pole, d) / self.angle;
        self.rate * (PLUME_FALLOUT_FAR + (1.0 - PLUME_FALLOUT_FAR) * (-x).exp())
    }

    /// How long the fallout takes to bury a lag or a salt at a direction, Gyr.
    pub fn burial_gyr(&self, d: [f64; 3]) -> f64 {
        PLUME_BURIAL_DEPTH / self.rate_at(d).max(1e-30)
    }

    /// The frost cover the fallout keeps at a direction, 0 to 1.
    pub fn frost_at(&self, d: [f64; 3]) -> f64 {
        (self.rate_at(d) / PLUME_FROST_RATE).min(1.0)
    }

    /// How much a fissure at a direction is one of the plume's warm vents, 0 to 1.
    pub fn vent_weight(&self, d: [f64; 3]) -> f64 {
        (-angle_between(self.pole, d) / self.angle).exp()
    }
}

/// A patch of chaos: the shell melted through or disrupted from below.
#[derive(Clone, Debug)]
pub struct Chaos {
    pub center: [f64; 3],
    pub radius: f64,
    pub time_gyr: f64,
    pub id: u64,
}

impl Chaos {
    pub fn cap_radius(&self, body_radius: f64) -> f64 {
        (1.0 + CHAOS_EDGE_WOBBLE + 0.1) * self.radius / body_radius
    }
}

/// The chaos present at the epoch: patches drawn per time bin within the retention age,
/// their number from the disrupted share of the surface over the mean patch area, their
/// places in the shell's heat field.
pub fn chaos_population(graph: &ParamGraph, present: f64, seed: u64) -> Vec<Chaos> {
    let fraction = graph.f64("tectonics.chaos_fraction");
    if fraction <= 0.0 {
        return Vec::new();
    }
    let body_radius = graph.f64("body.radius");
    let age = graph.f64("body.age");
    let retention = graph.f64("cratering.retention_age").max(1e-6);
    let (lo, hi) = CHAOS_DIAMETER;
    let mean_area = std::f64::consts::FRAC_PI_4 * (hi * hi - lo * lo) / (2.0 * (hi / lo).ln());
    let expected = fraction * 4.0 * std::f64::consts::PI * body_radius * body_radius / mean_area;
    let rng = Rng::new(seed ^ 0xC4A0_5000);
    let heat = GradientNoise::new(graph.f64("seed.tectonics") as u64 ^ 0xC4A0_4EA7);
    let window_start = (present - retention).max(0.0);
    let mut out = Vec::new();
    for (bin, b0, b1) in crate::stages::time_bins(age, 0.0, present) {
        if b1 <= window_start {
            continue;
        }
        let b0 = b0.max(window_start);
        let share = ((b1 - b0) / retention).min(1.0);
        let (count, _) = crate::crater::production::poisson(&rng, bin, 0, expected * share);
        for j in 0..count as u64 {
            let draw = 16 + j * 40;
            // Where: in the heat field, by rejection.
            let mut center = None;
            for attempt in 0..16u64 {
                let p = rng.unit_direction(bin, draw + 1 + 2 * attempt);
                let field = 0.5
                    + 0.5
                        * heat.fbm(
                            p,
                            CHAOS_HEAT_FREQUENCY / (2.0 * std::f64::consts::PI),
                            2,
                            2.0,
                            0.5,
                        );
                let w = 1.0 - CHAOS_HEAT_CONTRAST * (1.0 - field.clamp(0.0, 1.0));
                if rng.uniform(bin, draw + 2 + 2 * attempt) < w {
                    center = Some(p);
                    break;
                }
            }
            let Some(center) = center else { continue };
            let diameter = lo * (hi / lo).powf(rng.uniform(bin, draw + 34));
            out.push(Chaos {
                center,
                radius: 0.5 * diameter,
                time_gyr: b0 + (b1 - b0) * rng.uniform(bin, draw + 35),
                id: (bin << 32) | j,
            });
        }
    }
    out
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

/// A chaos patch on the resolved grid: inside a wobbling outline the matrix drops and
/// breaks into hummocks, stained with the ocean's salts, around rafts of the old surface
/// that keep what was on them.
pub fn chaos_compute(
    layers: &LayerSet,
    grid: &CubeSphere,
    body_radius: f64,
    chaos: &Chaos,
    stamper: &mut Stamper,
) -> StampOutput {
    let bedrock = layers.bedrock.as_slice();
    let noise = GradientNoise::new(chaos.id ^ 0xC4A0_0000);
    // Noise frequencies are lattice cells per radian: the radius over the feature length.
    let outline_frequency = body_radius / (2.0 * chaos.radius);
    let raft_frequency = body_radius / CHAOS_RAFT_SCALE;
    let hummock_frequency = body_radius / CHAOS_HUMMOCK_SPACING;
    // The share of a patch standing as rafts sets the threshold on the raft field.
    let raft_threshold = 0.55 * (1.0 - 2.0 * CHAOS_RAFT_SHARE);
    let time = chaos.time_gyr as f32;
    stamper.gather(
        grid,
        body_radius,
        chaos.center,
        chaos.cap_radius(body_radius) * body_radius,
    );
    let mut cells: Vec<CellWrite> = stamper
        .footprint
        .par_iter()
        .filter_map(|&(k, r, _)| {
            let d = grid.center(grid.cell(k));
            let edge = chaos.radius
                * (1.0 + CHAOS_EDGE_WOBBLE * noise.fbm(d, outline_frequency, 2, 2.0, 0.5));
            if r >= edge {
                return None;
            }
            let raft = noise.fbm(d, raft_frequency, 1, 2.0, 0.5) > raft_threshold;
            if raft {
                return None;
            }
            let taper = ((edge - r) / (0.1 * edge)).clamp(0.0, 1.0);
            let hummocks = CHAOS_HUMMOCK * noise.fbm(d, hummock_frequency, 2, 2.0, 0.5);
            let mut w = write(k, bedrock[k] as f64 + taper * (hummocks - CHAOS_DROP));
            w.surface = true;
            w.age = time;
            w.unit = UNIT_CHAOS;
            w.cover = 1.0;
            Some(w)
        })
        .collect();
    cells.sort_by_key(|c| c.idx);
    StampOutput { cells }
}

/// A cryovolcanic dome: brine-rich ice extruded onto a shell over rock.
#[derive(Clone, Debug)]
pub struct Cryodome {
    pub center: [f64; 3],
    pub radius: f64,
    /// Height as built, m; it relaxes with the shell.
    pub height: f64,
    pub time_gyr: f64,
    pub id: u64,
}

impl Cryodome {
    pub fn cap_radius(&self, body_radius: f64) -> f64 {
        1.3 * self.radius / body_radius
    }
}

/// The domes built through the body's life: per time bin, the dome share of the ice
/// magma budget spent in domes of volumes drawn log-uniform as a share of the body's
/// cube, placed at random.
pub fn cryodome_population(graph: &ParamGraph, present: f64, seed: u64) -> Vec<Cryodome> {
    if graph.str("interior.tectonic_regime") != "ice_shell_rock"
        || !crate::stages::ice::is_ice_surface(graph)
    {
        return Vec::new();
    }
    // The magma rate is kept in m³/s; the budget here is in km³ per Gyr.
    let rate = graph.f64("volcanism.rate") * GYR / 1e9 * CRYODOME_SHARE;
    if rate <= 0.0 {
        return Vec::new();
    }
    let body_radius = graph.f64("body.radius");
    let age = graph.f64("body.age");
    let r_km3 = (body_radius / 1e3).powi(3);
    let (lo, hi) = (
        CRYODOME_VOLUME_SHARE.0 * r_km3,
        CRYODOME_VOLUME_SHARE.1 * r_km3,
    );
    let rng = Rng::new(seed ^ 0xD0CE_5000);
    let mut carry = 0.0;
    // The next dome's volume is drawn ahead and waits until the budget covers it, so
    // that a slow budget builds domes of every size rather than only the smallest.
    let mut pending: Option<f64> = None;
    let mut out = Vec::new();
    for (bin, b0, b1) in crate::stages::time_bins(age, 0.0, present) {
        carry += rate * (b1 - b0);
        let mut j = 0u64;
        loop {
            let draw = 16 + j * 8;
            let volume =
                *pending.get_or_insert_with(|| lo * (hi / lo).powf(rng.uniform(bin, draw)));
            if carry < volume || j >= 1000 {
                break;
            }
            pending = None;
            carry -= volume;
            // V = (2π/3) aspect r³ for a bell of height 2 aspect r.
            let radius =
                (3.0 * volume * 1e9 / (2.0 * std::f64::consts::PI * CRYODOME_ASPECT)).cbrt();
            out.push(Cryodome {
                center: rng.unit_direction(bin, draw + 1),
                radius,
                height: 2.0 * CRYODOME_ASPECT * radius,
                time_gyr: b0 + (b1 - b0) * rng.uniform(bin, draw + 2),
                id: (bin << 32) | j,
            });
            j += 1;
        }
    }
    out
}

/// A dome on the resolved grid at its height now: a bell with streaked flanks, its
/// surface the dome's own, salted.
pub fn cryodome_compute(
    layers: &LayerSet,
    grid: &CubeSphere,
    body_radius: f64,
    dome: &Cryodome,
    height_now: f64,
    stamper: &mut Stamper,
) -> StampOutput {
    let bedrock = layers.bedrock.as_slice();
    let ages = layers.surface_age.as_slice();
    let noise = GradientNoise::new(dome.id ^ 0xD0CE_0000);
    let streak_frequency = body_radius / (0.2 * dome.radius);
    let time = dome.time_gyr as f32;
    stamper.gather(grid, body_radius, dome.center, dome.radius);
    let mut cells: Vec<CellWrite> = stamper
        .footprint
        .par_iter()
        .map(|&(k, r, _)| {
            let d = grid.center(grid.cell(k));
            let x = (r / dome.radius).min(1.0);
            let bell = (1.0 - x * x).powf(1.5);
            let streaks = CRYODOME_FLANK_ROUGHNESS
                * bell
                * (1.0 - bell)
                * noise.fbm(d, streak_frequency, 2, 2.0, 0.5);
            let mut w = write(k, bedrock[k] as f64 + height_now * (bell + streaks));
            if bell > 0.05 && ages[k] <= time + 1e-6 {
                w.surface = true;
                w.age = time;
                w.unit = UNIT_CRYODOME;
                w.cover = CRYODOME_SALT as f32;
            }
            w
        })
        .collect();
    cells.sort_by_key(|c| c.idx);
    StampOutput { cells }
}

/// After the stamps: the domes' substrate is partly clean ice, and the floors of
/// craters deep enough to reach the brines carry salts in patches.
pub fn post_pass(
    layers: &mut LayerSet,
    body_radius: f64,
    domes: &[Cryodome],
    craters: &[CraterRecord],
    scaling: &crate::crater::Scaling,
    present: f64,
) {
    let grid = layers.grid;
    if !domes.is_empty() {
        let units = layers.unit_id.as_slice().to_vec();
        let clean = layers.clean.as_mut_slice();
        crate::crater::with_stamper(&grid, |stamper| {
            for dome in domes {
                stamper.gather(&grid, body_radius, dome.center, dome.radius);
                for &(k, _, _) in &stamper.footprint {
                    if units[k] == UNIT_CRYODOME {
                        clean[k] = clean[k].max(CRYODOME_CLEAN as f32);
                    }
                }
            }
        });
    }
    let faculae: Vec<&CraterRecord> = craters
        .iter()
        .filter(|c| {
            c.visible
                && present - c.time_gyr < 4.0 * FACULA_FADE_GYR
                && dimensions(c.diameter, scaling).depth >= FACULA_BRINE_DEPTH
        })
        .collect();
    if faculae.is_empty() {
        return;
    }
    let units = layers.unit_id.as_slice().to_vec();
    let ages = layers.surface_age.as_slice().to_vec();
    let cover = layers.cover.as_mut_slice();
    let noise = GradientNoise::new(0xFAC0_1A00);
    let patch_frequency = body_radius / FACULA_PATCH_SCALE;
    crate::crater::with_stamper(&grid, |stamper| {
        for c in faculae {
            let dim = dimensions(c.diameter, scaling);
            stamper.gather(
                &grid,
                body_radius,
                c.center,
                FACULA_FLOOR_SHARE * dim.radius,
            );
            for &(k, r, _) in &stamper.footprint {
                let floor = units[k] == UNIT_CRATER_FLOOR || units[k] == UNIT_DEEP_FLOOR;
                if !floor || (ages[k] as f64 - c.time_gyr).abs() > 1e-3 {
                    continue;
                }
                let d = grid.center(grid.cell(k));
                let taper = 1.0 - r / (FACULA_FLOOR_SHARE * dim.radius);
                let patch = noise.fbm(d, patch_frequency, 2, 2.0, 0.5) + 0.6 * taper;
                if patch > 0.45 {
                    cover[k] = cover[k].max(1.0);
                }
            }
        }
    });
}

/// The share of the surface that is chaos matrix.
pub fn chaos_share(layers: &LayerSet) -> f64 {
    let grid = layers.grid;
    let units = layers.unit_id.as_slice();
    let (mut chaos, mut total) = (0.0, 0.0);
    for (k, &u) in units.iter().enumerate() {
        let a = grid.cell_area(grid.cell(k));
        total += a;
        if u == UNIT_CHAOS {
            chaos += a;
        }
    }
    chaos / total.max(1e-300)
}

/// The leading and trailing apexes of a locked moon: the sub-parent point lies on +x and
/// the orbit runs prograde about +z, so the moon moves towards -y and its trailing
/// hemisphere faces +y.
pub const TRAILING_APEX: [f64; 3] = [0.0, 1.0, 0.0];

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fallout_falls_off_from_the_pole_and_buries_quickly() {
        let f = Fallout {
            pole: [0.0, 0.0, -1.0],
            rate: 1e5,
            angle: 25f64.to_radians(),
        };
        assert!((f.rate_at([0.0, 0.0, -1.0]) - 1e5).abs() < 1.0);
        let equator = f.rate_at([1.0, 0.0, 0.0]);
        assert!(equator > 1e3 && equator < 1e4, "{equator}");
        assert!(f.burial_gyr([1.0, 0.0, 0.0]) < 1e-5);
        assert!(f.frost_at([1.0, 0.0, 0.0]) >= 1.0);
    }
}
