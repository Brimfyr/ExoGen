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

/// The share of the body a plume's fallout keeps under frost (F47, F52): the frost cover
/// at each distance from the pole, integrated over the sphere.
pub fn mean_frost_share(rate_pole: f64) -> f64 {
    if rate_pole <= 0.0 {
        return 0.0;
    }
    let angle = PLUME_FALLOUT_ANGLE_DEG.to_radians();
    let n = 360;
    let mut sum = 0.0;
    for k in 0..n {
        let theta = (k as f64 + 0.5) / n as f64 * std::f64::consts::PI;
        let rate =
            rate_pole * (PLUME_FALLOUT_FAR + (1.0 - PLUME_FALLOUT_FAR) * (-theta / angle).exp());
        sum += (1.0 - (-rate / PLUME_FROST_RATE).exp()) * theta.sin();
    }
    0.5 * sum * std::f64::consts::PI / n as f64
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
        // A soft saturation: the fallout wins against the lag and the gardening by
        // degrees, so the far side keeps some of its own contrast under a thinner frost.
        1.0 - (-self.rate_at(d) / PLUME_FROST_RATE).exp()
    }

    /// How much a fissure at a direction is one of the plume's warm vents, 0 to 1.
    pub fn vent_weight(&self, d: [f64; 3]) -> f64 {
        (-angle_between(self.pole, d) / self.angle).exp()
    }
}

/// The plumes' fallout mantles the ground it lies on (F74). Snow laid over billions of
/// years creeps down the slopes of what it covers, filling the hollows and thinning on the
/// crests, so relief smaller than its thickness is buried and what stands larger is
/// rounded by no more than it: Enceladus's old craters are muted (Kirchoff and Schenk
/// 2009). Its loose top, what drains into fractures, is a few hundred metres (Martin et al.
/// 2023); the whole deposit, sintered below that, mantles. The ground is drawn toward its
/// own low-pass over MANTLE_SPREAD times the deposit's thickness `cover`, m, by no more
/// than the thickness either way; the deposit's mean, which only raises the ground, is
/// left out. The mantle buries only what is older than it: the finish lays it before the
/// unresolved craters, which count only those dug since the cover above them grew
/// thinner than they are deep (LooseCover::buried_before).
pub fn mantle(layers: &mut LayerSet, body_radius: f64, cover: &[f32]) {
    let grid = layers.grid;
    let cell = grid.delta() * body_radius;
    let thickest = cover.iter().cloned().fold(0.0f32, f32::max) as f64;
    if thickest <= 0.0 {
        return;
    }
    if std::env::var_os("EXOGEN_DEBUG").is_some() {
        let mut sorted: Vec<f32> = cover.iter().copied().filter(|h| *h > 0.0).collect();
        sorted.sort_by(f32::total_cmp);
        let at = |q: f64| sorted[((sorted.len() - 1) as f64 * q) as usize];
        eprintln!(
            "mantle: fallout {:.0} m at the 10th percentile, {:.0} at the median, {:.0} at the 90th, {:.0} at most",
            at(0.1),
            at(0.5),
            at(0.9),
            thickest
        );
    }
    let lows =
        crate::stages::relaxation::low_passes(&layers.bedrock, cell, MANTLE_SPREAD * thickest);
    if lows.is_empty() {
        return;
    }
    layers
        .bedrock
        .as_mut_slice()
        .par_iter_mut()
        .enumerate()
        .for_each(|(k, z)| {
            let h = cover[k] as f64;
            if h <= 0.0 {
                return;
            }
            let here = *z as f64;
            let want = MANTLE_SPREAD * h;
            // The low-pass at that length, between the pyramid's levels about it.
            let low = match lows.iter().position(|(len, _)| *len >= want) {
                Some(0) => here + (want / lows[0].0) * (lows[0].1[k] as f64 - here),
                Some(i) => {
                    let (l0, l1) = (lows[i - 1].0, lows[i].0);
                    let t = (want / l0).ln() / (l1 / l0).ln();
                    lows[i - 1].1[k] as f64 * (1.0 - t) + lows[i].1[k] as f64 * t
                }
                None => lows[lows.len() - 1].1[k] as f64,
            };
            *z = (here + (low - here).clamp(-h, h)) as f32;
        });
}

/// The loose cover the plumes' fallout has laid on the ground (F64): the snow of the
/// present plume and of each past episode's, each falling off from its own pole as the
/// fallout does (F47), summed over the time a surface has stood.
pub struct LooseCover {
    /// Each plume's pole in the present frame, and when it ran, Gyr.
    plumes: Vec<([f64; 3], f64, f64)>,
    /// The rate at a plume's pole, m/Gyr.
    rate: f64,
    angle: f64,
    /// The depth a crack stays open to, m: what it can drain.
    crack_depth: f64,
}

impl LooseCover {
    /// The plumes that have run by `present`: one at the centre of each episode of
    /// focusing (F53) while it lasted, each venting as the present one does.
    pub fn new(
        graph: &ParamGraph,
        episodes: &[crate::stages::ice::FocusEpisode],
        present: f64,
    ) -> Option<Self> {
        let rate = graph.f64("volcanism.plume_fallout_rate");
        if rate <= 0.0 {
            return None;
        }
        let shell = crate::stages::ice::ShellMap::new(graph);
        let pole = shell.pole?;
        let mut plumes: Vec<([f64; 3], f64, f64)> = episodes
            .iter()
            .filter(|e| e.start < present)
            .map(|e| (e.center, e.start, e.end.min(present)))
            .collect();
        let start = (graph.f64("body.age") - FOCUSING_EPISODE_DURATION_GYR).max(0.0);
        if shell.has_active_terrain() && present > start {
            plumes.push((pole, start, present));
        }
        // A crack stays open down to where the weight of the ice matches the tension
        // that opened it.
        let crack_depth =
            ICE_CRACK_STRESS / (DENSITY_CRUST_ICE * graph.f64("body.gravity").max(1e-6));
        Some(Self::of_plumes(plumes, rate, crack_depth))
    }

    /// The cover of given plumes: each one's pole and the times it ran between, Gyr, all
    /// at one rate at the pole, m/Gyr; and the depth a crack can drain, m.
    pub fn of_plumes(plumes: Vec<([f64; 3], f64, f64)>, rate: f64, crack_depth: f64) -> Self {
        Self {
            plumes,
            rate,
            angle: PLUME_FALLOUT_ANGLE_DEG.to_radians(),
            crack_depth,
        }
    }

    /// How much of a cover a crack drains, m: no more than the crack is deep.
    pub fn drained(&self, cover: f64) -> f64 {
        cover.min(self.crack_depth)
    }

    /// The time before which what was dug at a direction lies buried by `until`, Gyr:
    /// after it no more than `depth` of cover fell. `since` when even all the cover
    /// since then is thinner (F74).
    pub fn buried_before(&self, d: [f64; 3], since: f64, until: f64, depth: f64) -> f64 {
        if self.thickness_between(d, since, until) <= depth {
            return since;
        }
        let (mut lo, mut hi) = (since, until);
        for _ in 0..40 {
            let mid = 0.5 * (lo + hi);
            if self.thickness_between(d, mid, until) > depth {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    }

    /// The cover laid at a direction between two times, m: what lies on a surface that
    /// formed at `since` by the time `until`.
    pub fn thickness_between(&self, d: [f64; 3], since: f64, until: f64) -> f64 {
        self.plumes
            .iter()
            .map(|&(pole, t0, t1)| {
                let span = t1.min(until) - t0.max(since);
                if span <= 0.0 {
                    return 0.0;
                }
                let x = angle_between(pole, d) / self.angle;
                self.rate * (PLUME_FALLOUT_FAR + (1.0 - PLUME_FALLOUT_FAR) * (-x).exp()) * span
            })
            .sum()
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
                if rng.uniform(bin ^ crate::rng::ACCEPT_STREAM, draw + attempt) < w {
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
    // About CHAOS_EDGE_WAVES waves round the patch's circumference.
    let outline_frequency =
        body_radius * CHAOS_EDGE_WAVES / (2.0 * std::f64::consts::PI * chaos.radius);
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
    StampOutput {
        cells,
        composition: None,
        stain: None,
    }
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
                time_gyr: b0 + (b1 - b0) * rng.uniform(bin, draw + 3),
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
    StampOutput {
        cells,
        composition: None,
        stain: None,
    }
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

    /// Europa's chaos is placed over every longitude (the acceptance draw has a stream
    /// of its own).
    #[test]
    fn europa_s_chaos_reaches_every_longitude() {
        let (g, _) = crate::presets::preset("Europa").unwrap().to_graph();
        let chaos = chaos_population(&g, g.f64("body.age"), 11);
        let mut sectors = [0usize; 6];
        for c in &chaos {
            let lon = c.center[1].atan2(c.center[0]).to_degrees() + 180.0;
            sectors[((lon / 60.0) as usize).min(5)] += 1;
        }
        let mean = chaos.len() as f64 / 6.0;
        assert!(chaos.len() > 300);
        assert!(
            sectors.iter().all(|&n| n as f64 > 0.6 * mean),
            "{sectors:?} of {}",
            chaos.len()
        );
    }

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
        // The frost saturates softly: nearly full at the equator, part of the way on
        // the far side, where the floor of the fallout is the frost threshold.
        assert!(f.frost_at([1.0, 0.0, 0.0]) > 0.9);
        let far = f.frost_at([0.0, 0.0, 1.0]);
        assert!(far > 0.5 && far < 0.75, "{far}");
    }

    /// The loose cover is the fallout of every plume that has run (F64): the present
    /// episode's at the present pole at least, the far share of every episode's far
    /// from the poles, and nothing on a surface that has only just formed.
    #[test]
    fn loose_cover_is_the_sum_of_the_plumes_fallout() {
        let (g, _) = crate::presets::preset("Enceladus").unwrap().to_graph();
        let shell = crate::stages::ice::ShellMap::new(&g);
        let episodes = crate::stages::ice::focusing_episodes(&g, &shell);
        let present = g.f64("body.age");
        let cover = LooseCover::new(&g, &episodes, present).unwrap();
        let pole = shell.pole.unwrap();
        let rate = g.f64("volcanism.plume_fallout_rate");
        let one = rate * FOCUSING_EPISODE_DURATION_GYR;
        let at_pole = cover.thickness_between(pole, 0.0, present);
        assert!(at_pole >= 0.999 * one, "{at_pole} against {one}");
        assert_eq!(cover.thickness_between(pole, present, present), 0.0);
        // The least any old surface carries: the far share of every plume.
        let floor = PLUME_FALLOUT_FAR * one * (episodes.len() + 1) as f64;
        let mut least = f64::INFINITY;
        let mut most: f64 = 0.0;
        for k in 0..2000u64 {
            let d = Rng::new(7).unit_direction(k, 0);
            let h = cover.thickness_between(d, 0.0, present);
            least = least.min(h);
            most = most.max(h);
        }
        eprintln!(
            "enceladus: {} past episodes; loose cover {least:.0} to {most:.0} m on old ground, {at_pole:.0} m at the pole",
            episodes.len()
        );
        assert!(least >= 0.999 * floor, "{least} against {floor}");
        // The same order as the pit chains of its cratered plains have it (Martin et al.
        // 2017: 250 m on average, to 700 m), and two to three times over: nine plumes'
        // fallout, where the real body's history is not known.
        assert!(least > 100.0 && most < 3000.0, "{least} to {most}");
        // A crack there drains a few hundred metres of it.
        assert!(
            (300.0..500.0).contains(&cover.drained(most)),
            "{}",
            cover.drained(most)
        );
    }
}
