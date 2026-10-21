// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! Flood volcanism on stagnant lids (PLAN.md 4.4).
//!
//! The thermal history gives the magma production through time. A fraction of it erupts,
//! less as the lithosphere thickens and once global contraction has put the crust into
//! compression, from vents that favour thin crust and low ground. Each eruption floods the
//! terrain from its vent by priority fill: the lava runs downhill and fills one flat lake
//! in the low it reaches, taking in ground in order of height and crossing passes into
//! deeper ground as the level reaches them, so basin floors flood first, a basin fed many
//! times fills as one lake, and older crater rims stand out of the fill where they rise
//! above it. A flood resets the
//! surface age and unit of what it covers. The loop carries the coarse fill; the finishing
//! pass sharpens each flood's margins against the resolved terrain.

use rayon::prelude::*;
use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashMap, HashSet};

use crate::constants::*;
use crate::crater::{CellWrite, StampOutput, Stamper};
use crate::fields::{Field, LayerSet};
use crate::geom::{angle_between, dot, normalize};
use crate::grid::CubeSphere;
use crate::grid::projection::{PaddedField, SparseField, sample_bilinear_sparse};
use crate::noise::GradientNoise;
use crate::params::ParamGraph;
use crate::rng::Rng;
use crate::stages::deform::Deformation;
use crate::stages::io::{Lobe, PathRules, path_point, rotate_about, trace_flow};
use crate::stages::plume::RiseVent;
use crate::thermal::ThermalHistory;
use crate::units::*;

pub const VERSION: u64 = 21;

/// The drift of a body's lavas from their material's mean at a time, -1 to 1 (F55): a
/// smooth function of the eruptive history with the composition period, so units close in
/// time are alike and units an age apart differ, written by each eruption into the
/// surface it lays and kept by the ejecta of that surface.
pub fn composition_drift(time_gyr: f64) -> f32 {
    let noise = GradientNoise::new(0x1A7A_71A7);
    let t = time_gyr / LAVA_COMPOSITION_PERIOD_GYR;
    (noise.fbm([t, 0.3, 0.7], 1.0, 1, 2.0, 0.5) / 0.6).clamp(-1.0, 1.0) as f32
}

/// The composition of the mantle source under a vent, -1 to 1 (F94): a smooth field of
/// LAVA_SOURCE_ORDER cycles about the body, for the unevenness of the mantle that its
/// magma ocean's cumulates left.
pub fn source_composition(vent: [f64; 3]) -> f64 {
    let noise = GradientNoise::new(0x5A0C_E5ED);
    (noise.fbm(vent, LAVA_SOURCE_ORDER, 2, 2.0, 0.5) / 0.6).clamp(-1.0, 1.0)
}

/// The composition of a flood's lava, -1 to 1 (F55, F94): its eruptive time's drift, its
/// source's and its own batch's, `batch` a unit normal draw, a share of each.
pub fn flood_composition(time_gyr: f64, vent: [f64; 3], batch: f64) -> f32 {
    (LAVA_DRIFT_SHARE * composition_drift(time_gyr) as f64
        + LAVA_SOURCE_SHARE * source_composition(vent)
        + LAVA_BATCH_SHARE * batch)
        .clamp(-1.0, 1.0) as f32
}

/// The composition of the lava an overturn lays everywhere at once (F120): its time's
/// drift, one lava over the whole plains; the mantle's source field read at every place had
/// painted its noise across them (each basin's own lava, from the source under it, is for
/// the plains' units to come, audit M10).
pub fn overturn_composition(time_gyr: f64) -> f32 {
    (LAVA_DRIFT_SHARE * composition_drift(time_gyr) as f64).clamp(-1.0, 1.0) as f32
}

/// The share of a rocky body's magma that its crust lets through (F125): magma rising from
/// `depth`, the stiff lid's base, gains buoyancy through the mantle part of its path and
/// loses it through a crust `thickness` thick and lighter than itself, and erupts by what
/// the loss leaves of the gain; a crust denser than the magma takes nothing.
pub fn crust_buoyancy(magma: f64, crust: f64, thickness: f64, depth: f64) -> f64 {
    let loss = (magma - crust).max(0.0) * thickness;
    let gain = (DENSITY_ROCK - magma).max(1.0) * (depth - thickness).max(1.0);
    (1.0 - loss / gain).clamp(0.0, 1.0)
}

/// A coarse cell a flood covered: the thickness added and what the surface was before.
#[derive(Clone, Copy, Debug)]
pub struct FloodCell {
    pub idx: u32,
    pub thickness: f32,
    pub prev_unit: u16,
    pub prev_age: f32,
    /// The composition of the lava it found, if it found lava (F55).
    pub prev_composition: f32,
}

/// One pond of a flood: the level it filled to and the coarse cells it covered; or a
/// sheet the lava left on the slope it ran down (M15), its thickness over the ground,
/// whose level is its highest surface.
#[derive(Clone, Debug)]
pub struct Pond {
    pub level: f32,
    pub cells: Vec<FloodCell>,
    pub sheet: bool,
}

/// What one run of an eruption laid (M15): a sheet down a slope or a lake, the cells it
/// covered with the surface each found and the lava laid there, and its level: a lake's,
/// or a sheet's highest surface.
struct Run {
    level: f64,
    sheet: bool,
    cells: Vec<(usize, f32, f32)>,
}

/// An eruption drawn for a step, not yet emplaced.
#[derive(Clone, Copy, Debug)]
pub struct FloodEvent {
    pub vent: usize,
    pub volume: f64,
    pub time_gyr: f64,
    pub id: u64,
    /// How fast it pours out, m³/s (F117).
    pub effusion: f64,
}

/// What a flood did to one coarse cell over all its ponds: the thickness added, the last
/// lake's level, or none where only sheets covered it (M15), and the surface before.
#[derive(Clone, Copy, Debug)]
pub struct Covered {
    pub thickness: f32,
    pub level: f32,
    pub sheet: bool,
    pub prev_unit: u16,
    pub prev_age: f32,
    pub prev_composition: f32,
}

/// A flood that was emplaced, kept for the sidecar and for sharpening its margins at the
/// finishing resolution.
#[derive(Clone, Debug)]
pub struct FloodRecord {
    pub vent: [f64; 3],
    /// Erupted volume, m³.
    pub volume: f64,
    pub time_gyr: f64,
    /// Area covered, m².
    pub area: f64,
    pub ponds: Vec<Pond>,
    /// The lava's composition, -1 to 1 (F55, F94).
    pub composition: f32,
    pub id: u64,
    /// How fast it poured out, m³/s, which sets the form of its flows (F117).
    pub effusion: f64,
}

/// How compound a lava's flows are for how fast it pours out (F117): 1 at or below the
/// slow rate, 0 at or above the fast, between them by the logarithm.
pub fn compound_share(effusion: f64) -> f64 {
    let (slow, fast) = LAVA_COMPOUND_RATE;
    ((fast.ln() - effusion.max(1e-9).ln()) / (fast.ln() - slow.ln())).clamp(0.0, 1.0)
}

impl FloodRecord {
    /// Coarse cells covered, with their total thickness, final level and the surface
    /// the flood found. A sheet has no level of its own, lying its thickness over the
    /// ground (M15): a cell only sheets covered has none.
    pub fn cells(&self) -> HashMap<u32, Covered> {
        let mut out: HashMap<u32, Covered> = HashMap::new();
        for pond in &self.ponds {
            for c in &pond.cells {
                let e = out.entry(c.idx).or_insert(Covered {
                    thickness: 0.0,
                    level: f32::NAN,
                    sheet: true,
                    prev_unit: c.prev_unit,
                    prev_age: c.prev_age,
                    prev_composition: c.prev_composition,
                });
                e.thickness += c.thickness;
                if !pond.sheet {
                    e.level = pond.level;
                    e.sheet = false;
                }
            }
        }
        out
    }

    /// A cap holding every covered cell: centre direction and angular radius, without
    /// margin.
    pub fn cap(&self, grid: &CubeSphere) -> ([f64; 3], f64) {
        let cells = self.cells();
        if cells.is_empty() {
            return (self.vent, 0.0);
        }
        // Summed in the cells' order, so that the centre does not depend on the map's.
        let mut keys: Vec<u32> = cells.keys().copied().collect();
        keys.sort_unstable();
        let mut sum = [0.0; 3];
        for &idx in &keys {
            let d = grid.center(grid.cell(idx as usize));
            sum[0] += d[0];
            sum[1] += d[1];
            sum[2] += d[2];
        }
        let center = if sum[0] * sum[0] + sum[1] * sum[1] + sum[2] * sum[2] < 1e-12 {
            self.vent
        } else {
            normalize(sum)
        };
        let radius = keys
            .iter()
            .map(|&idx| crate::geom::angle_between(center, grid.center(grid.cell(idx as usize))))
            .fold(0.0, f64::max);
        (center, radius)
    }
}

/// Priority of a frontier cell: height first, then distance from the seed, so a lake on
/// level ground grows outward from where the lava arrives rather than in index order.
#[derive(PartialEq)]
struct Key(f32, f32, usize);
impl Eq for Key {}
impl PartialOrd for Key {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for Key {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.0
            .total_cmp(&other.0)
            .then(self.1.total_cmp(&other.1))
            .then(self.2.cmp(&other.2))
    }
}

/// Flood volcanism through the coarse loop.
/// A magmatic province: a region of the mantle more fertile than the rest.
#[derive(Clone, Debug)]
pub struct Province {
    pub center: [f64; 3],
    /// Gaussian angular radius, rad.
    pub sigma: f64,
}

pub struct FloodProcess {
    grid: CubeSphere,
    body_radius: f64,
    thermal: ThermalHistory,
    rng: Rng,
    /// Cell areas, m².
    areas: Vec<f64>,
    /// Magmatic provinces and their weight per coarse cell, which the vents favour.
    pub provinces: Vec<Province>,
    province_weight: Vec<f64>,
    /// Each rise's own share of that weight per coarse cell (F108).
    rise_province: Vec<Vec<f32>>,
    /// Cumulative lava thickness and the time of the last flood, per coarse cell, for
    /// laying out the ridges at the end.
    pub thickness: Vec<f32>,
    pub last_time: Vec<f32>,
    /// Visit marks for the pond search, by generation.
    visited: Vec<u32>,
    generation: u32,
    heap: BinaryHeap<Reverse<Key>>,
    /// Extrusive volume carried between steps, m³.
    carry: f64,
    /// On an episodic lid: when the surface was renewed, before which no flows are kept.
    event_time: Option<f64>,
    /// The share of the ice magma budget that floods rather than builds domes.
    flood_share: f64,
    /// A rocky crust's buoyancy for its magma (F125): the magma's and the crust's densities,
    /// kg/m³, and the crust's thickness, m. None for cryolava.
    buoyancy: Option<(f64, f64, f64)>,
    age: f64,
    pub records: Vec<FloodRecord>,
}

const NEIGHBOURS: [(i32, i32); 4] = [(1, 0), (-1, 0), (0, 1), (0, -1)];

impl FloodProcess {
    /// None for regimes whose volcanism is handled elsewhere.
    /// The plume rises are provinces in their own right: a plume head is the most
    /// fertile mantle there is, and its rise draws the eruptions before any other region.
    pub fn new(
        graph: &ParamGraph,
        grid: CubeSphere,
        thermal: &ThermalHistory,
        rises: &[crate::stages::plume::Rise],
    ) -> Option<Self> {
        if crate::stages::io::is_io_regime(graph) {
            return None;
        }
        // A rocky crust holds its magma down by how much lighter than the magma it is
        // (F125): the magma by the mantle's iron, the crust by what it is.
        let feo = graph.f64("body.mantle_feo");
        let crust_density = match graph.str("materials.bedrock") {
            "anorthosite" => CRUST_DENSITY_ANORTHOSITE,
            "felsic" => CRUST_DENSITY_FELSIC,
            _ => CRUST_DENSITY_MAFIC_BASE + CRUST_DENSITY_MAFIC_PER_FEO * feo,
        };
        let rocky = Some((
            MAGMA_DENSITY_BASE + MAGMA_DENSITY_PER_FEO * feo,
            crust_density,
            graph.f64("tectonics.initial_crust_thickness"),
        ));
        let buoyancy = match graph.str("interior.tectonic_regime") {
            "stagnant_lid" | "episodic_lid" => rocky,
            _ => None,
        };
        let event_time = match graph.str("interior.tectonic_regime") {
            "stagnant_lid" => None,
            "episodic_lid" => Some(
                (graph.f64("body.age") - graph.f64("tectonics.resurfacing_event_age")).max(0.0),
            ),
            // Cryolava on an ice shell over rock, from the ice magma budget: slush that
            // is buoyant through the shell. Water under a shell over an ocean is denser
            // than the ice and does not flood the surface; what that shell does with
            // its heat, chaos, fissures and fallout, is its own process.
            "ice_shell_rock" => None,
            _ => return None,
        };
        let body_radius = graph.f64("body.radius");
        let areas = (0..grid.cell_count())
            .map(|k| grid.cell_area(grid.cell(k)) * body_radius * body_radius)
            .collect();
        let rng = Rng::new(graph.f64("seed.volcanism") as u64 ^ 0x7F10_0D5A);
        let n_provinces = graph.f64("volcanism.province_count").round().max(0.0) as usize;
        let mut provinces: Vec<Province> = rises
            .iter()
            .map(|r| Province {
                center: r.center,
                sigma: r.radius / body_radius,
            })
            .collect();
        for k in provinces.len()..n_provinces {
            provinces.push(Province {
                center: rng.unit_direction(0xB10 + k as u64, 1),
                sigma: MAGMA_PROVINCE_SIGMA_RAD,
            });
        }
        let province_weight = (0..grid.cell_count())
            .map(|k| {
                if provinces.is_empty() {
                    return 1.0;
                }
                let d = grid.center(grid.cell(k));
                provinces.iter().fold(MAGMA_PROVINCE_FLOOR, |w, p| {
                    w + (-(angle_between(p.center, d) / p.sigma).powi(2)).exp()
                })
            })
            .collect();
        let rise_province = provinces
            .iter()
            .take(rises.len())
            .map(|p| {
                (0..grid.cell_count())
                    .map(|k| {
                        let d = grid.center(grid.cell(k));
                        (-(angle_between(p.center, d) / p.sigma).powi(2)).exp() as f32
                    })
                    .collect()
            })
            .collect();
        Some(Self {
            grid,
            body_radius,
            thermal: thermal.clone(),
            rng,
            areas,
            provinces,
            province_weight,
            rise_province,
            thickness: vec![0.0; grid.cell_count()],
            last_time: vec![f32::NEG_INFINITY; grid.cell_count()],
            visited: vec![0; grid.cell_count()],
            generation: 0,
            heap: BinaryHeap::new(),
            carry: 0.0,
            event_time,
            flood_share: if graph.str("interior.tectonic_regime") == "ice_shell_rock" {
                1.0 - CRYODOME_SHARE
            } else {
                1.0
            },
            buoyancy,
            age: graph.f64("body.age"),
            records: Vec::new(),
        })
    }

    /// Io's lava as floods (F117): no provinces nor rises, its eruptions the paterae's
    /// (`io::eruptions`), emplaced as any body's floods are.
    pub fn for_io(graph: &ParamGraph, grid: CubeSphere, thermal: &ThermalHistory) -> Self {
        let body_radius = graph.f64("body.radius");
        let n = grid.cell_count();
        Self {
            grid,
            body_radius,
            thermal: thermal.clone(),
            rng: Rng::new(graph.f64("seed.volcanism") as u64 ^ 0x7F10_0D5A),
            areas: (0..n)
                .map(|k| grid.cell_area(grid.cell(k)) * body_radius * body_radius)
                .collect(),
            provinces: Vec::new(),
            province_weight: vec![1.0; n],
            rise_province: Vec::new(),
            thickness: vec![0.0; n],
            last_time: vec![f32::NEG_INFINITY; n],
            visited: vec![0; n],
            generation: 0,
            heap: BinaryHeap::new(),
            carry: 0.0,
            event_time: None,
            flood_share: 1.0,
            buoyancy: None,
            age: graph.f64("body.age"),
            records: Vec::new(),
        }
    }

    /// Take the resurfacing event's fill into the thickness and last-flow records.
    /// The density of the lava a rocky body's floods lay, kg/m³ (F125); none for cryolava.
    pub fn magma_density(&self) -> Option<f64> {
        self.buoyancy.map(|b| b.0)
    }

    pub fn note_resurfacing(&mut self, fill: &[f32], time_gyr: f64) {
        for (k, &f) in fill.iter().enumerate() {
            if f > 0.0 {
                self.thickness[k] += f;
                self.last_time[k] = time_gyr as f32;
            }
        }
    }

    /// Fraction of the magma that reaches the surface at a time: the extrusive share,
    /// cut back as the lithosphere thickens and as contraction closes the pathways. On a
    /// rocky body the share is what the crust's buoyancy leaves (F125): magma rising from the
    /// stiff lid's base gains buoyancy through the mantle part of its path and loses it
    /// through a crust lighter than itself, and erupts by what the loss leaves of the gain,
    /// so that a light anorthositic crust holds most of a dense basalt down while the lid is
    /// thin.
    pub fn eruptive_fraction(&self, t: f64) -> f64 {
        let s = self.thermal.at(t);
        let strain = 2.0 * s.contraction / self.body_radius;
        let ascent = 1.0
            / (1.0
                + ((s.lithosphere - VOLCANISM_ASCENT_LITHOSPHERE) / VOLCANISM_ASCENT_WIDTH).exp());
        let share = match self.buoyancy {
            Some((magma, crust, thickness)) => {
                VOLCANISM_EXTRUSIVE_FRACTION_ROCK
                    * crust_buoyancy(magma, crust, thickness, s.melt_depth)
            }
            None => VOLCANISM_EXTRUSIVE_FRACTION,
        };
        share * ascent * (-(strain / VOLCANISM_SHUTOFF_STRAIN).powi(2)).exp()
    }

    /// Extrusive volume between two times, m³, by midpoint rule over sub-intervals. On an
    /// episodic lid nothing before the event is kept, and a small share erupts after it.
    pub fn budget(&self, t0: f64, t1: f64) -> f64 {
        let (t0, share) = match self.event_time {
            Some(event) if t1 <= event => return 0.0,
            Some(event) => (t0.max(event), EPISODIC_POST_EVENT_FRACTION),
            None => (t0, 1.0),
        };
        let n = 8;
        let dt = (t1 - t0) / n as f64;
        let mut v = 0.0;
        for i in 0..n {
            let t = t0 + (i as f64 + 0.5) * dt;
            v += self.thermal.at(t).magma_rate * self.eruptive_fraction(t) * dt * GYR;
        }
        v * share * self.flood_share
    }

    /// Vent weights, cumulative, and the share of the floods' weight a rise's shields leave
    /// them: thin crust and low ground, and any rise's boost. A rise's shields take the
    /// share of its eruptions that those started hold (F108): what the rise adds to the
    /// weights, by its fertile mantle and its dome, goes by that share to its shields,
    /// whose edifices and flows are their own, so that a rise's late lava builds its
    /// shields rather than flooding its dome and canyons.
    fn vent_weights(
        &self,
        layers: &LayerSet,
        boost: Option<&[f64]>,
        rises: &[RiseVent],
    ) -> (Vec<f64>, f64) {
        let crust = layers.crust_thickness.as_slice();
        let bedrock = layers.bedrock.as_slice();
        let n = bedrock.len();
        let mut heights: Vec<f32> = bedrock.to_vec();
        heights.sort_by(f32::total_cmp);
        let median = heights[n / 2] as f64;
        let crust_min = crust.iter().cloned().fold(f32::INFINITY, f32::min) as f64;
        let base: Vec<f64> = (0..n)
            .map(|k| {
                (-(crust[k] as f64 - crust_min) / VOLCANISM_VENT_CRUST_SCALE).exp()
                    * (-((bedrock[k] as f64 - median) / VOLCANISM_VENT_HEIGHT_SCALE)
                        .clamp(-4.0, 4.0))
                    .exp()
                    * self.areas[k]
            })
            .collect();
        let lift = |k: usize| boost.map_or(1.0, |b| b[k]);
        let full: Vec<f64> = (0..n)
            .map(|k| base[k] * self.province_weight[k] * lift(k))
            .collect();
        let mut w = full.clone();
        for (r, rise) in rises.iter().enumerate() {
            let Some(own) = self.rise_province.get(r) else {
                continue;
            };
            if rise.capture <= 0.0 {
                continue;
            }
            let mut dome = vec![0.0f64; n];
            for &(idx, d) in &rise.dome {
                dome[idx] = d as f64;
            }
            // What the rise adds at each cell: the weight less what it would be without
            // the rise's province and dome.
            for k in 0..n {
                let without = (self.province_weight[k] - own[k] as f64)
                    * (lift(k) - PLUME_VENT_BOOST * dome[k]).max(1.0);
                w[k] -= rise.capture * (full[k] - base[k] * without).max(0.0);
            }
        }
        let whole: f64 = full.iter().sum();
        let mut acc = 0.0;
        let cumulative: Vec<f64> = w
            .into_iter()
            .map(|v| {
                acc += v;
                acc
            })
            .collect();
        let kept = if whole > 0.0 { acc / whole } else { 1.0 };
        (cumulative, kept)
    }

    /// Eruptions between two times, in order of formation, not yet emplaced. Vents are
    /// drawn against the surface as it is at the start of the step.
    pub fn sample(
        &mut self,
        layers: &LayerSet,
        boost: Option<&[f64]>,
        rises: &[RiseVent],
        t0: f64,
        t1: f64,
    ) -> Vec<FloodEvent> {
        let bins = crate::stages::time_bins(self.age, t0, t1);
        // What the shields take is theirs, not the floods' (F108).
        let kept = if rises.iter().any(|r| r.capture > 0.0) {
            Some(self.vent_weights(layers, boost, rises))
        } else {
            None
        };
        let share = kept.as_ref().map_or(1.0, |k| k.1);
        if std::env::var_os("EXOGEN_DEBUG").is_some() {
            let full: f64 = bins.iter().map(|&(_, b0, b1)| self.budget(b0, b1)).sum();
            eprintln!(
                "floods {t0:.3} to {t1:.3}: budget {:.4e} km3, the shields take {:.4}",
                full / 1e9,
                1.0 - share
            );
        }
        let step_budget: f64 = bins
            .iter()
            .map(|&(_, b0, b1)| self.budget(b0, b1))
            .sum::<f64>()
            * share;
        if step_budget + self.carry < VOLCANISM_EVENT_VOLUME.0 {
            self.carry += step_budget;
            return Vec::new();
        }
        let (cumulative, share) = kept.unwrap_or_else(|| self.vent_weights(layers, boost, rises));
        let total = *cumulative.last().unwrap_or(&0.0);
        if total <= 0.0 {
            return Vec::new();
        }
        let mut out = Vec::new();
        for (bin, b0, b1) in bins {
            self.sample_bin(&cumulative, total, share, bin, b0, b1, &mut out);
        }
        out.sort_by(|a, b| a.time_gyr.total_cmp(&b.time_gyr).then(a.id.cmp(&b.id)));
        out
    }

    /// The eruptions of one time bin, drawn with the bin as the random index against
    /// the vent weights of the step.
    #[allow(clippy::too_many_arguments)]
    fn sample_bin(
        &mut self,
        cumulative: &[f64],
        total: f64,
        share: f64,
        bin: u64,
        t0: f64,
        t1: f64,
        out: &mut Vec<FloodEvent>,
    ) {
        let step_index = bin;
        let mut budget = self.budget(t0, t1) * share + self.carry;
        let (lo, hi) = VOLCANISM_EVENT_VOLUME;
        let mut j = 0u64;
        while budget >= lo && j < 100_000 {
            let draw = 16 + j * 4;
            let volume = lo * (hi / lo).powf(self.rng.uniform(step_index, draw));
            let volume = volume.min(budget);
            let u = self.rng.uniform(step_index, draw + 1) * total;
            let vent = cumulative
                .partition_point(|&c| c < u)
                .min(cumulative.len() - 1);
            let time_gyr = t0 + (t1 - t0) * self.rng.uniform(step_index, draw + 2);
            let (e_lo, e_hi) = FLOOD_EFFUSION_RATE;
            out.push(FloodEvent {
                vent,
                volume,
                time_gyr,
                id: (step_index << 32) | j,
                effusion: e_lo * (e_hi / e_lo).powf(self.rng.uniform(step_index, draw + 3)),
            });
            budget -= volume;
            j += 1;
        }
        self.carry = budget.max(0.0);
    }

    /// Erupt one event onto the surface as it is now, and record it.
    pub fn apply(&mut self, layers: &mut LayerSet, ev: FloodEvent) {
        let mut rec = self.erupt(layers, ev.vent, ev.volume, ev.time_gyr, ev.id);
        rec.effusion = ev.effusion;
        if !rec.ponds.is_empty() {
            self.records.push(rec);
        }
    }

    /// Sample and emplace a step's eruptions.
    pub fn step(&mut self, layers: &mut LayerSet, t0: f64, t1: f64, _step_index: u64) {
        for ev in self.sample(layers, None, &[], t0, t1) {
            self.apply(layers, ev);
        }
    }

    /// Steepest descent over a surface from a cell to a local minimum, or to the last cell
    /// within `within` of the vent: the cells in order, from the first to the low.
    fn descend(
        &self,
        height: impl Fn(usize) -> f32,
        within: impl Fn(usize) -> bool,
        mut idx: usize,
    ) -> Vec<usize> {
        let mut path = vec![idx];
        for _ in 0..self.grid.cell_count() {
            let cell = self.grid.cell(idx);
            let mut best = (height(idx), idx);
            for (di, dj) in NEIGHBOURS {
                let n = self.grid.index(self.grid.neighbor(cell, di, dj));
                let h = height(n);
                if h < best.0 {
                    best = (h, n);
                }
            }
            if best.1 == idx || !within(best.1) {
                break;
            }
            idx = best.1;
            path.push(idx);
        }
        path
    }

    /// Flood the terrain from a vent with a volume. The lava runs down by the steepest way
    /// to the nearest low, leaving a sheet a flow's thickness deep on each cell it crosses
    /// as far as its volume goes (audit M15: running down without it, the lava gathered in
    /// the lows and left the slopes bare), and fills a flat lake with the rest: cells join
    /// in order of height from the low, until the volume cannot lift the lake to the next
    /// cell. A lake that fills to its lowest pass spills the rest over it: the next cell
    /// would lie more than a flow's thickness below the lake, on ground that falls away
    /// beyond the pass, so the lava runs on down from there to the next low, leaving its
    /// sheet on the way, and fills a lake of its own, and so on (audit W7: one lake taking
    /// in the ground beyond its pass stood level over it, a cliff where it stopped). A lake
    /// whose spill reaches back to an earlier one rises with it as one. A basin fed many
    /// times fills as one lake.
    fn erupt(
        &mut self,
        layers: &mut LayerSet,
        vent: usize,
        volume: f64,
        time_gyr: f64,
        id: u64,
    ) -> FloodRecord {
        let grid = self.grid;
        let vent_dir = grid.center(grid.cell(vent));
        let drape = FLOOD_SHEET_THICKNESS as f32;
        // Flows run only so far from their vent.
        let body_radius = self.body_radius;
        let within = |n: usize| {
            crate::geom::angle_between(vent_dir, grid.center(grid.cell(n))) * body_radius
                <= FLOOD_MAX_REACH
        };
        // The runs, each a sheet down a slope or a lake, in order; the lava surface they
        // leave, which a later run of the eruption finds.
        let mut runs: Vec<Run> = Vec::new();
        let mut surface: HashMap<usize, f32> = HashMap::new();
        let mut remaining = volume;
        let mut start = vent;
        for _ in 0..FLOOD_MAX_SPILLS {
            let path = {
                let bedrock = layers.bedrock.as_slice();
                let height = |i: usize| surface.get(&i).copied().unwrap_or(bedrock[i]);
                self.descend(height, within, start)
            };
            let seed = *path.last().unwrap_or(&start);
            // The lava runs down to the low a flow's thickness deep, leaving its sheet on
            // each cell it crosses and, as wide as its volume spreads it, on the cells round
            // its path it need not climb to, nearest first, as far as its volume goes.
            {
                let bedrock = layers.bedrock.as_slice();
                let half_width = 0.5
                    * (remaining / (FLOOD_SHEET_THICKNESS * FLOOD_SHEET_ASPECT)).sqrt()
                    / body_radius;
                let mut sheet: Vec<(usize, f32, f32)> = Vec::new();
                'path: for &c in &path[..path.len() - 1] {
                    let centre = grid.center(grid.cell(c));
                    let h_c = surface.get(&c).copied().unwrap_or(bedrock[c]);
                    let top = h_c + FLOOD_SHEET_THICKNESS as f32;
                    // The cells it spreads to from here, through those it can reach.
                    let mut swath = vec![c];
                    let mut seen: HashSet<usize> = HashSet::from([c]);
                    let mut at = 0;
                    while at < swath.len() {
                        let cell = grid.cell(swath[at]);
                        at += 1;
                        for (di, dj) in NEIGHBOURS {
                            let n = grid.index(grid.neighbor(cell, di, dj));
                            if seen.insert(n)
                                && angle_between(centre, grid.center(grid.cell(n))) <= half_width
                                && within(n)
                                && surface.get(&n).copied().unwrap_or(bedrock[n]) <= top
                            {
                                swath.push(n);
                            }
                        }
                    }
                    for n in swath {
                        if surface.contains_key(&n) {
                            continue;
                        }
                        let a = self.areas[n];
                        let t = (remaining / a).min(FLOOD_SHEET_THICKNESS);
                        if t <= FLOOD_MIN_THICKNESS {
                            break 'path;
                        }
                        let h = bedrock[n];
                        sheet.push((n, h, t as f32));
                        surface.insert(n, h + t as f32);
                        remaining -= t * a;
                    }
                }
                if !sheet.is_empty() {
                    let level = sheet
                        .iter()
                        .map(|&(_, h, t)| h + t)
                        .fold(f32::NEG_INFINITY, f32::max);
                    runs.push(Run {
                        level: level as f64,
                        sheet: true,
                        cells: sheet,
                    });
                }
            }
            if remaining <= FLOOD_MIN_THICKNESS * self.areas[seed] {
                break;
            }
            let bedrock = layers.bedrock.as_slice();
            let height = |i: usize| surface.get(&i).copied().unwrap_or(bedrock[i]);
            self.generation = self.generation.wrapping_add(1);
            if self.generation == 0 {
                self.visited.iter_mut().for_each(|v| *v = 0);
                self.generation = 1;
            }
            let generation = self.generation;
            let seed_dir = grid.center(grid.cell(seed));
            self.heap.clear();
            self.visited[seed] = generation;
            self.heap.push(Reverse(Key(height(seed), 0.0, seed)));
            // The cells taken, each with the surface it had and the lake's highest
            // ground when it joined: the pass it was reached over.
            let mut popped: Vec<(usize, f32, f32)> = Vec::new();
            let (mut area_sum, mut vol_sum) = (0.0f64, 0.0f64);
            let mut top = f32::NEG_INFINITY;
            // The lowest ground next to the lake past the flows' reach: the lake rises no
            // more than a flow's thickness above it.
            let mut outlet = f32::INFINITY;
            let mut next: Option<f32> = None;
            let mut spill: Option<usize> = None;
            let mut full = false;
            while let Some(Reverse(Key(h, _, idx))) = self.heap.pop() {
                // The lowest cell on the frontier. Below the lake's highest ground by more
                // than a flow's thickness, it lies past a pass: the lake is full to it and
                // the rest spills on from here.
                if h < top - drape {
                    spill = Some(idx);
                    break;
                }
                // A cell joins if the volume can lift everything taken so far, and the
                // cell, to a flow's thickness above it: a flow advances only where it
                // stands that thick, so the lake ends in a flow's front. Nor does it rise
                // more than that above the ground past the flows' reach.
                let a = self.areas[idx];
                let lift = (h + drape) as f64;
                if lift * (area_sum + a) - (vol_sum + a * h as f64) > remaining
                    || lift > (outlet + drape) as f64
                {
                    next = Some(h + drape);
                    break;
                }
                area_sum += a;
                vol_sum += a * h as f64;
                top = top.max(h);
                popped.push((idx, h, top));
                if popped.len() >= FLOOD_MAX_CELLS {
                    full = true;
                    break;
                }
                let cell = grid.cell(idx);
                for (di, dj) in NEIGHBOURS {
                    let n = grid.index(grid.neighbor(cell, di, dj));
                    if self.visited[n] != generation {
                        self.visited[n] = generation;
                        if !within(n) {
                            outlet = outlet.min(height(n));
                            continue;
                        }
                        let dist = crate::geom::angle_between(seed_dir, grid.center(grid.cell(n)));
                        self.heap.push(Reverse(Key(height(n), dist as f32, n)));
                    }
                }
            }
            if popped.is_empty() {
                // Too little to pond, or at the flows' reach: what it can, a sheet at most,
                // stops in the low.
                let t = (remaining / self.areas[seed]).min(FLOOD_SHEET_THICKNESS);
                if t > FLOOD_MIN_THICKNESS {
                    let h = height(seed);
                    runs.push(Run {
                        level: (h as f64) + t,
                        sheet: true,
                        cells: vec![(seed, h, t as f32)],
                    });
                }
                break;
            }
            // The level that spends the volume over the cells taken, every one of them
            // below it: full to its pass where it spills; where it stopped, at the volume
            // or at the most cells a lake takes, no more than a flow's thickness above the
            // ground at its edge; and no more than that above the ground past the flows'
            // reach.
            let mut level = (remaining + vol_sum) / area_sum;
            if spill.is_some() {
                level = level.min(top as f64);
            } else if let Some(h) = next {
                level = level.min(h.max(top) as f64);
            }
            if full && let Some(Reverse(Key(h, _, _))) = self.heap.peek() {
                level = level.min((h + drape).max(top) as f64);
            }
            level = level.min((outlet + drape) as f64);
            let lake: Vec<(usize, f32, f32)> = popped
                .iter()
                .filter(|&&(_, h, pass)| pass as f64 <= level && (h as f64) < level)
                .map(|&(idx, h, _)| (idx, h, (level - h as f64) as f32))
                .collect();
            let used: f64 = lake
                .iter()
                .map(|&(idx, _, t)| t as f64 * self.areas[idx])
                .sum();
            remaining -= used;
            for &(idx, _, _) in &lake {
                surface.insert(idx, level as f32);
            }
            if !lake.is_empty() {
                runs.push(Run {
                    level,
                    sheet: false,
                    cells: lake,
                });
            }
            match spill {
                Some(cell) if remaining > 0.0 => start = cell,
                _ => break,
            }
        }
        // Lay the runs' lava, each over what the ones before it left.
        let mut ponds = Vec::new();
        let mut placed = 0.0;
        {
            let bedrock = layers.bedrock.as_mut_slice();
            let ages = layers.surface_age.as_mut_slice();
            let units = layers.unit_id.as_mut_slice();
            let melt = layers.melt.as_mut_slice();
            let cover = layers.cover.as_mut_slice();
            let composition = layers.composition.as_mut_slice();
            let lava_age = layers.lava_age.as_mut_slice();
            let drift = flood_composition(time_gyr, vent_dir, self.rng.normal(id, 0xBA7C));
            let mut laid: HashMap<usize, (f64, u16, f32, f32)> = HashMap::new();
            for run in &runs {
                let mut recorded = Vec::new();
                for &(idx, _, t) in &run.cells {
                    let t = t as f64;
                    if t <= FLOOD_MIN_THICKNESS {
                        continue;
                    }
                    // The pile subsides by a share of its thickness: what is recorded and
                    // added is the height left after that.
                    let added = (1.0 - LAVA_SUBSIDENCE) * t;
                    placed += t * self.areas[idx];
                    let e =
                        laid.entry(idx)
                            .or_insert((0.0, units[idx], ages[idx], composition[idx]));
                    e.0 += t;
                    recorded.push(FloodCell {
                        idx: idx as u32,
                        thickness: added as f32,
                        prev_unit: e.1,
                        prev_age: e.2,
                        prev_composition: e.3,
                    });
                    self.thickness[idx] += t as f32;
                    self.last_time[idx] = time_gyr as f32;
                    bedrock[idx] += added as f32;
                    ages[idx] = time_gyr as f32;
                    units[idx] = UNIT_LAVA;
                    melt[idx] = 0.0;
                    cover[idx] = (e.0 / LAVA_COVER_DEPTH).min(1.0) as f32;
                    composition[idx] = drift;
                    lava_age[idx] = time_gyr as f32;
                }
                if !recorded.is_empty() {
                    ponds.push(Pond {
                        level: run.level as f32,
                        cells: recorded,
                        sheet: run.sheet,
                    });
                }
            }
        }
        let mut rec = FloodRecord {
            vent: vent_dir,
            volume: placed,
            time_gyr,
            area: 0.0,
            ponds,
            composition: flood_composition(time_gyr, vent_dir, self.rng.normal(id, 0xBA7C)),
            id,
            effusion: FLOOD_EFFUSION_RATE.0,
        };
        // Summed in the cells' order, so that the area does not depend on the map's.
        let mut keys: Vec<u32> = rec.cells().keys().copied().collect();
        keys.sort_unstable();
        rec.area = keys.iter().map(|&i| self.areas[i as usize]).sum();
        rec
    }
}

/// A coarse cell the fine fill may reach: the level, and the previous surface when the
/// cell is inside the coarse footprint.
type Reach = (f32, Option<(u16, f32, f32)>);

/// A fine cell a flood's pass visits: its index, the level the lava could reach there and,
/// inside the coarse footprint, the surface the flood found.
type Visit = (usize, f32, Option<(u16, f32, f32)>);

/// The finishing pass's lava bookkeeping. The base carries every flood's coarse fill,
/// smeared over a coarse cell at its margin. Each flood's pass takes its own smear out
/// where it goes, so the base minus what remains is the ground as the flood found it,
/// and later floods' smears ride on top until their turn.
pub struct LavaSmear {
    /// Every flood's smear, summed, on the fine grid.
    pub total: Vec<f32>,
    /// The smear of floods whose pass has run.
    pub done: Vec<f32>,
    /// Cells a flood's pass has written.
    pub touched: Vec<bool>,
    /// The time of the latest loop crater whose cavity covers each cell at the finishing
    /// resolution, or none: inside the rim of a crater younger than a flood the base is
    /// the crater stamped on the lava, and is right as it is; beyond it the crater's
    /// restore buries the pass's relief by its ejecta's thickness (M15).
    pub crater_time: Vec<f32>,
    /// Counts over all passes: cells written as lava, left above the level, left
    /// outside the outline, skipped under a younger crater.
    pub stats: [usize; 4],
}

impl LavaSmear {
    pub fn new(
        coarse: CubeSphere,
        grid: CubeSphere,
        floods: &[FloodRecord],
        extra: Option<&[f32]>,
    ) -> Self {
        let mut thickness = Field::<f32>::new(coarse);
        {
            let t = thickness.as_mut_slice();
            for f in floods {
                for (&idx, c) in &f.cells() {
                    t[idx as usize] += c.thickness;
                }
            }
            if let Some(e) = extra {
                for (v, x) in t.iter_mut().zip(e) {
                    *v += x;
                }
            }
        }
        let n = grid.cell_count();
        Self {
            total: thickness.resample(grid).as_slice().to_vec(),
            done: vec![0.0; n],
            touched: vec![false; n],
            crater_time: vec![f32::NEG_INFINITY; n],
            stats: [0; 4],
        }
    }

    /// Record a flood's pass: the smear it took out, the cells it wrote, and its counts.
    pub fn note(&mut self, taken: &[(u32, f32)], written: &StampOutput, stats: [usize; 4]) {
        for &(idx, s) in taken {
            self.done[idx as usize] += s;
        }
        for c in &written.cells {
            self.touched[c.idx as usize] = true;
        }
        for (a, b) in self.stats.iter_mut().zip(stats) {
            *a += b;
        }
    }
}

/// Sharpen one flood's margins at the finishing resolution. The base carries the coarse
/// fill smeared over a coarse cell at its edge; this removes the smeared thickness and
/// refills to the pond level wherever the resolved terrain lies below it, inside a
/// margin that follows the smeared outline of the coarse footprint with lobes of noise
/// along it, so a flow front over flat ground is lobate rather than blocky. What the
/// lava does not reach gets the surface it had before the flood; what later craters
/// resurfaced is left alone. Returns the writes and the smear taken out per cell, for
/// the bookkeeping.
#[allow(clippy::too_many_arguments)]
pub fn flood_compute(
    layers: &LayerSet,
    grid: &CubeSphere,
    coarse: &CubeSphere,
    rec: &FloodRecord,
    smear: &LavaSmear,
    deform: &Deformation,
    flows: Option<FloodFlows>,
    base: &[f32],
    stamper: &mut Stamper,
) -> (StampOutput, Vec<(u32, f32)>, [usize; 4]) {
    let ratio = grid.face_size() / coarse.face_size();
    if ratio == 0 || !grid.face_size().is_multiple_of(coarse.face_size()) {
        return (StampOutput::default(), Vec::new(), [0; 4]);
    }
    let mut stats = [0usize; 4];
    let covered = rec.cells();
    // The coarse thickness, sampled exactly as the base's resampling did, and the
    // footprint indicator whose smeared outline bounds the fine fill.
    let padded = SparseField::new(
        *coarse,
        covered.iter().map(|(&idx, c)| (idx as usize, c.thickness)),
    );
    // The footprint's indicator, each coarse cell half its own and a sixteenth each of
    // its eight neighbours', so that its contour rounds the staircase of the coarse
    // cells while a single flooded cell keeps its area; sampled bilinearly between, it
    // is the flood's outline at the fine resolution, which the lobes then wander.
    let padded_ind = {
        let mut cells: HashSet<usize> = HashSet::new();
        for &idx in covered.keys() {
            let c = coarse.cell(idx as usize);
            for di in -1..=1 {
                for dj in -1..=1 {
                    cells.insert(coarse.index(coarse.neighbor(c, di, dj)));
                }
            }
        }
        SparseField::new(
            *coarse,
            cells.into_iter().map(|i| {
                let c = coarse.cell(i);
                let mut v = 0.0f32;
                for di in -1..=1 {
                    for dj in -1..=1 {
                        let n = coarse.index(coarse.neighbor(c, di, dj)) as u32;
                        if covered.contains_key(&n) {
                            v += if di == 0 && dj == 0 { 0.5 } else { 1.0 / 16.0 };
                        }
                    }
                }
                (i, v)
            }),
        )
    };
    let lobes = GradientNoise::new(rec.id ^ 0x0F10_0D00);
    let lobe_frequency = 1.0 / (FLOOD_MARGIN_LOBE_CELLS * coarse.delta());
    // A flood that shows at the end is drawn as its flows (F94): its lobes, traced over
    // the ground from where its lava came in, lay its sheet, and where they do not reach
    // the ground keeps the surface it had, unless it lies so far below the level that the
    // lava ponds there. A buried flood keeps the coarse outline its successors hide.
    let lobed: Option<HashMap<usize, f32>> = flows.map(|f| {
        let lobes = flood_lobes(
            rec,
            coarse,
            grid,
            f.body_radius,
            f.ground,
            f.seed,
            compound_share(rec.effusion),
        );
        lobe_cover(&lobes, grid, f.body_radius, rec.id, stamper)
    });
    let soft = flows.is_some_and(|f| f.soft);
    let hollow_steps = flows.map_or(0, |f| {
        (f.hollow_reach / (grid.delta() * f.body_radius)).ceil() as usize
    });
    // Where the loop filled the flood deeper than its sheets, the lava pond fills the low
    // whatever its lobes do: judged by the flood's own coarse thickness, sampled smoothly,
    // since the ground pieced together from the coarse smears is not good to a sheet.
    let ponded_thickness =
        ((1.0 - LAVA_SUBSIDENCE) * FLOOD_PONDED_SHEETS * FLOOD_SHEET_THICKNESS) as f32;
    // The sheets' own thickness (M15): a sheet lies its thickness over the ground it
    // crosses, which the coarse thickness sampled smoothly follows along a slope.
    let sheets = SparseField::new(
        *coarse,
        covered
            .iter()
            .filter(|(_, c)| c.sheet)
            .map(|(&idx, c)| (idx as usize, c.thickness)),
    );
    // Coarse cells to visit: the footprint and its ring, with the level the lava could
    // reach there and, inside the footprint, the surface the flood found. A cell only
    // sheets covered, or a ring cell beside none but them, has no level: its lava lies its
    // thickness over the ground.
    let mut levels: HashMap<usize, Reach> = HashMap::new();
    for (&idx, c) in &covered {
        levels.insert(
            idx as usize,
            (c.level, Some((c.prev_unit, c.prev_age, c.prev_composition))),
        );
        let cell = coarse.cell(idx as usize);
        for di in -1..=1 {
            for dj in -1..=1 {
                let n = coarse.index(coarse.neighbor(cell, di, dj));
                if !covered.contains_key(&(n as u32)) {
                    let e = levels.entry(n).or_insert((c.level, None));
                    e.0 = e.0.max(c.level);
                }
            }
        }
    }
    let bedrock = layers.bedrock.as_slice();
    let ages = layers.surface_age.as_slice();
    let units = layers.unit_id.as_slice();
    let deep = layers.deep.as_slice();
    let melt = layers.melt.as_slice();
    let cover = layers.cover.as_slice();
    let found = layers.composition.as_slice();
    let time = rec.time_gyr as f32;
    // A hollow the lava ponds in is finer relief, a crater the loop never saw: the fine
    // ground lies below the floor the lava's level stands over, `floor`, while the ground the
    // loop saw there, the base less what came after, does not (audit W7, M15). Where the
    // loop's own ground runs below the floor too, the cell lies on an open slope beyond the
    // low, down the coarse resampling's ramp toward lower ground or down a slope its flows
    // crossed, which the lava runs down rather than ponds in; the gaps earlier floods' fine
    // outlines left there stand below the base as well, but are no hollow either.
    let loop_ground = |idx: usize, remaining: f32, later: f32| base[idx] - remaining - later;
    let mut cells = Vec::with_capacity(levels.len() * (ratio * ratio) as usize);
    let mut taken = Vec::with_capacity(levels.len() * (ratio * ratio) as usize);
    let mut coarse_cells: Vec<(&usize, &Reach)> = levels.iter().collect();
    coarse_cells.sort_by_key(|(k, _)| **k);
    // The fine cells to visit: those of the coarse footprint and its ring, and those the
    // flood's lobes reach beyond them (F94), at the flood's highest level.
    let mut visits: Vec<Visit> = Vec::with_capacity(levels.len() * (ratio * ratio) as usize);
    for (&cidx, &(level, previous)) in coarse_cells {
        let c = coarse.cell(cidx);
        for a in 0..ratio {
            for b in 0..ratio {
                let fine = crate::grid::Cell {
                    face: c.face,
                    i: c.i * ratio + a,
                    j: c.j * ratio + b,
                };
                visits.push((grid.index(fine), level, previous));
            }
        }
    }
    // The cells of the coarse footprint and its ring come first; those its lobes reach
    // beyond follow.
    let in_ring = visits.len();
    if let Some(lobed) = &lobed {
        let seen: HashSet<usize> = visits.iter().map(|v| v.0).collect();
        let top = rec
            .ponds
            .iter()
            .map(|p| p.level)
            .fold(f32::NEG_INFINITY, f32::max);
        let mut beyond: Vec<usize> = lobed
            .iter()
            .filter(|(i, c)| **c > 0.0 && !seen.contains(*i))
            .map(|(i, _)| *i)
            .collect();
        beyond.sort_unstable();
        visits.extend(beyond.into_iter().map(|i| (i, top, None)));
    }
    // A cell to follow through the flood passes, for debugging: EXOGEN_DEBUG_CELL=lon,lat.
    let watched: Option<usize> = std::env::var("EXOGEN_DEBUG_CELL").ok().and_then(|v| {
        let mut it = v.split(',').map(|x| x.trim().parse::<f64>().ok());
        let (lon, lat) = (it.next()??, it.next()??);
        let dir = crate::grid::direction_from_lon_lat(lon.to_radians(), lat.to_radians());
        Some(grid.index(grid.locate(dir).0))
    });
    // What each visited cell offers the flood: the ground as the flood found it, the level
    // its lava could reach there, and whether its flows reached it.
    struct Work {
        idx: usize,
        level: f32,
        previous: Option<(u16, f32, f32)>,
        own: f32,
        ground: f32,
        carry: f32,
        younger: bool,
        near: bool,
        hollow: bool,
        reached: bool,
        ponded: bool,
        lobe: f32,
        // The level the loop's pond stood at in the cell's coarse cell, if any, and whether
        // the cell lies in one of a soft flood's lakes.
        pond: f32,
        lake: bool,
        // The share of the cell the flood's margin covers, and whether the cell lies on
        // its margin outside what it holds whole (F131).
        share: f32,
        margin: bool,
    }
    let mut work: Vec<Work> = visits
        .into_par_iter()
        .enumerate()
        .map(|(k, (idx, level, previous))| {
            let d = grid.center(grid.cell(idx));
            let own = sample_bilinear_sparse(coarse, &padded, d);
            let mut w = Work {
                idx,
                level,
                previous,
                own,
                ground: bedrock[idx],
                carry: 0.0,
                younger: smear.crater_time[idx] > time + 1e-6,
                near: k < in_ring,
                hollow: false,
                reached: false,
                ponded: false,
                lobe: lobed
                    .as_ref()
                    .map_or(0.0, |l| l.get(&idx).copied().unwrap_or(0.0)),
                pond: level,
                lake: false,
                share: 1.0,
                margin: false,
            };
            // A buried flood's sheet lies its thickness, sampled smoothly, wherever the loop
            // laid it: its margin is hidden under the flows after it, and its outline cut a
            // step a sheet deep along the coarse cells, a staircase on the plains (M15).
            w.reached = match &lobed {
                // A soft flood's lobes reach where they go, and its lakes are the fine
                // ground's closed lows (below), not its coarse thickness's contour, which runs
                // in steps along the coarse cells.
                Some(_) if soft => w.lobe > 0.0,
                // A lobe holds the cells it covers more than half of, and lays its share of
                // the rest (F131).
                Some(_) => {
                    let ponds = own > ponded_thickness;
                    if !ponds {
                        w.share = w.lobe.clamp(0.0, 1.0);
                    }
                    ponds || w.lobe >= 0.5
                }
                None if level.is_nan() => true,
                None => {
                    // The outline by the share of the cell beyond its margin (F131).
                    let outline = |e: [f64; 3]| {
                        let mut o = sample_bilinear_sparse(coarse, &padded_ind, e) as f64;
                        if o > 0.0 && o < 1.0 {
                            o += FLOOD_MARGIN_WOBBLE * lobes.fbm(e, lobe_frequency, 3, 2.0, 0.5);
                        }
                        o
                    };
                    w.share = crate::fields::level_share(
                        &outline,
                        d,
                        FLOOD_MARGIN_THRESHOLD,
                        grid.delta(),
                    ) as f32;
                    w.share >= 0.5
                }
            };
            w.margin = !w.reached && !soft && w.share > 0.0 && w.share < 0.5;
            // The ground as the flood found it: the base less the smears not yet
            // taken out and less what the slow processes, the domes, shields,
            // canyon and coronae, have added since; both ride on top of whatever
            // is written.
            let remaining = smear.total[idx] - smear.done[idx];
            let later = deform.since_flood(d, rec.time_gyr) as f32;
            let ground = bedrock[idx] - remaining - later;
            let sheet = level.is_nan();
            let level = if sheet {
                ground + sample_bilinear_sparse(coarse, &sheets, d) / (1.0 - LAVA_SUBSIDENCE as f32)
            } else {
                level
            };
            w.level = level;
            if w.younger {
                return w;
            }
            // A hollow finer than the loop's cells, a crater the loop never saw, lies
            // deeper below the level than the ground the loop flooded there (F95). Both
            // grounds carry the same ramps between the coarse cells, so only the finer
            // relief tells, where the loop's own ground stands at the floor or above.
            let floor = level - own / (1.0 - LAVA_SUBSIDENCE as f32) - FLOOD_SHEET_THICKNESS as f32;
            w.hollow = !soft
                && w.near
                && !sheet
                && ground < floor
                && loop_ground(idx, remaining, later) >= floor;
            // The ground is the ground the flood found, on every body: what stands above
            // the level ends the flood there, a rim, an older floor or a slope. Levelled to
            // a sheet below the level where the units, resampled from the loop's cells,
            // still said this lava (as the resampling's ramps were once read), slopes had
            // stood in terraces a coarse cell square, on a shield's flank and in old
            // craters on Mars as on Io's plains (F117, F120).
            w.ground = ground;
            w.carry = remaining - own + later;
            w
        })
        .collect();
    // A soft flood's lakes (F117): its lava ponds in the fine ground's closed lows that its
    // flows run into within its cells and their ring, each filled to where it spills
    // (priority flood, Barnes et al. 2014) and no higher than the loop's pond stood, the
    // lava draining away over the ring's edge wherever the ground there lies lower. The
    // loop's pond stood in a coarse low that the fine ground, smoothed over the loop's
    // cell, holds in part or not at all: filled to the pond's level wherever the ground lay
    // below it, a lake had stood perched on a slope and ended in steps along the coarse
    // cells. Elsewhere its lobes lay their sheets over the ground they cross. The finer
    // hollows of F95 are these lakes on a soft flood: the gaps between the lobes of the
    // floods before it, which the loop saw filled with their coarse lava, had passed for
    // craters it never saw and filled to the level of the lava beside them, across the
    // earlier flood's coarse footprint.
    if soft {
        let index: HashMap<usize, usize> = work
            .iter()
            .enumerate()
            .filter(|(_, w)| w.near)
            .map(|(k, w)| (w.idx, k))
            .collect();
        let around = |idx: usize| {
            let cell = grid.cell(idx);
            (-1..=1)
                .flat_map(move |di| (-1..=1).map(move |dj| (di, dj)))
                .filter(|&d| d != (0, 0))
                .map(move |(di, dj)| grid.index(grid.neighbor(cell, di, dj)))
        };
        let mut spill = vec![f32::INFINITY; work.len()];
        let mut heap = BinaryHeap::new();
        for (k, w) in work.iter().enumerate() {
            if w.near && around(w.idx).any(|n| !index.contains_key(&n)) {
                spill[k] = w.ground;
                heap.push(Reverse(Key(w.ground, 0.0, k)));
            }
        }
        while let Some(Reverse(Key(at, _, k))) = heap.pop() {
            if at > spill[k] {
                continue;
            }
            for n in around(work[k].idx) {
                if let Some(&j) = index.get(&n) {
                    let up = work[j].ground.max(at);
                    if up < spill[j] {
                        spill[j] = up;
                        heap.push(Reverse(Key(up, 0.0, j)));
                    }
                }
            }
        }
        // The closed lows, and those of them the flows run into, each a lake at one level:
        // where it spills, and no higher than the loop's pond stood in any of its cells (a
        // lake at each cell's own pond level had stood in steps along the coarse cells).
        let low = |w: &Work, top: f32| w.near && !w.younger && w.ground < top;
        let mut seen = vec![false; work.len()];
        for start in 0..work.len() {
            if seen[start] || !low(&work[start], spill[start]) {
                continue;
            }
            seen[start] = true;
            let mut lake = vec![start];
            let mut entered = false;
            let mut at = 0;
            while at < lake.len() {
                let k = lake[at];
                entered |= work[k].lobe > 0.0;
                for n in around(work[k].idx) {
                    if let Some(&j) = index.get(&n) {
                        entered |= work[j].lobe > 0.0;
                        if !seen[j] && low(&work[j], spill[j]) {
                            seen[j] = true;
                            lake.push(j);
                        }
                    }
                }
                at += 1;
            }
            if entered {
                let level = lake
                    .iter()
                    .map(|&k| work[k].pond)
                    .fold(f32::NEG_INFINITY, f32::max)
                    .min(spill[start]);
                for k in lake {
                    let w = &mut work[k];
                    if w.ground < level {
                        w.level = level;
                        w.reached = true;
                        w.lake = true;
                    }
                }
            }
        }
    }
    taken.extend(work.iter().map(|w| (w.idx as u32, w.own)));
    // The lava ponds in the finer hollows it comes to, filling each to its level from
    // where it comes in over the hollow's edge (F95): from the lava lying in the low, over
    // the cells its flows covered that stand no higher than they drape, into the hollows
    // beside them. A crater whose rim the lava breached or overtopped fills; one whose rim
    // stands above it keeps its bowl, whatever a lobe traced over the smoothed ground drew
    // across it, since the flow went round it. A hollow the visited cells cut across goes
    // on past them, its ground there above the level only where the finer relief does not
    // reach below it, since the coarse ground beyond the ring stood above the level too;
    // the lava follows it as far as the loop's largest unseen crater reaches.
    if lobed.is_some() && !soft {
        let drape = FLOOD_SHEET_THICKNESS as f32;
        // Lava lying in the low or its ring: the pond and the flows next to it, and the
        // sheets running down to it, not a lobe's sheet far beyond. A sheet fills the
        // finer hollows it crosses to its surface; down an open slope the loop's own
        // ground falls below that, and the slope stays as it is (the strips its flows
        // missed time after time had stood as trenches between them).
        let lava =
            |w: &Work| w.near && !w.younger && !w.hollow && w.reached && w.ground < w.level + drape;
        let mut index: HashMap<usize, usize> =
            work.iter().enumerate().map(|(k, w)| (w.idx, k)).collect();
        let neighbours = |idx: usize| {
            let cell = grid.cell(idx);
            (-1..=1)
                .flat_map(move |di| (-1..=1).map(move |dj| (di, dj)))
                .map(move |(di, dj)| grid.index(grid.neighbor(cell, di, dj)))
        };
        // The body of the lava: its connected stretches at least a coarse cell in size. A
        // sliver a lobe drew inside a rim the flow could not cross is no part of it.
        let mut seen = vec![false; work.len()];
        let mut queue: Vec<usize> = Vec::new();
        let body = (ratio * ratio) as usize;
        for start in 0..work.len() {
            if seen[start] || !lava(&work[start]) {
                continue;
            }
            seen[start] = true;
            let mut stretch = vec![start];
            let mut at = 0;
            while at < stretch.len() {
                for n in neighbours(work[stretch[at]].idx) {
                    if let Some(&j) = index.get(&n)
                        && !seen[j]
                        && lava(&work[j])
                    {
                        seen[j] = true;
                        stretch.push(j);
                    }
                }
                at += 1;
            }
            if stretch.len() >= body && stretch.iter().any(|&k| work[k].previous.is_some()) {
                queue.extend(stretch);
            }
        }
        let mut seen = vec![false; work.len()];
        for &k in &queue {
            seen[k] = true;
        }
        let mut steps = vec![0usize; work.len()];
        for w in work.iter_mut() {
            if w.hollow {
                w.reached = false;
            }
        }
        while let Some(k) = queue.pop() {
            let level = work[k].level;
            for n in neighbours(work[k].idx) {
                {
                    if let Some(&j) = index.get(&n) {
                        if seen[j] {
                            continue;
                        }
                        if work[j].hollow && !work[j].younger {
                            work[j].ponded = true;
                            work[j].reached = true;
                        } else if !lava(&work[j]) {
                            // A cell beside the pond lying below its floor ponds too,
                            // whatever the level of its own coarse cell, a flood's
                            // sheet or the ring round it: the fill follows the fine
                            // relief, not the coarse cells (the edge of a coarse cell
                            // had stopped it in a step a cell long). The loop's own
                            // ground there stands at the floor, so it is finer relief.
                            let wj = &work[j];
                            let floor = level
                                - wj.own / (1.0 - LAVA_SUBSIDENCE as f32)
                                - FLOOD_SHEET_THICKNESS as f32;
                            if !(work[k].ponded || lava(&work[k]))
                                || wj.younger
                                || wj.ground >= floor
                                || base[wj.idx] - (wj.carry + wj.own) < floor
                            {
                                continue;
                            }
                            work[j].level = level;
                            work[j].hollow = true;
                            work[j].ponded = true;
                            work[j].reached = true;
                        }
                        seen[j] = true;
                        queue.push(j);
                        continue;
                    }
                    if !(work[k].ponded || lava(&work[k]))
                        || steps[k] >= hollow_steps
                        || smear.crater_time[n] > time + 1e-6
                    {
                        continue;
                    }
                    let d = grid.center(grid.cell(n));
                    let own = sample_bilinear_sparse(coarse, &padded, d);
                    let remaining = smear.total[n] - smear.done[n];
                    let later = deform.since_flood(d, rec.time_gyr) as f32;
                    let ground = bedrock[n] - remaining - later;
                    let floor =
                        level - own / (1.0 - LAVA_SUBSIDENCE as f32) - FLOOD_SHEET_THICKNESS as f32;
                    if ground >= floor || loop_ground(n, remaining, later) < floor {
                        continue;
                    }
                    index.insert(n, work.len());
                    steps.push(steps[k] + 1);
                    seen.push(true);
                    queue.push(work.len());
                    taken.push((n as u32, own));
                    work.push(Work {
                        idx: n,
                        level,
                        previous: None,
                        own,
                        ground,
                        carry: remaining - own + later,
                        younger: false,
                        near: false,
                        hollow: true,
                        reached: true,
                        ponded: true,
                        lobe: 0.0,
                        pond: level,
                        lake: false,
                        share: 1.0,
                        margin: false,
                    });
                }
            }
        }
    }
    // The composition each written cell takes: the flood's where its lava lies, the lava it
    // found where it gives the surface back, whose coarse cells would otherwise carry the
    // flood's tint over its whole footprint.
    let mut compositions: Vec<f32> = Vec::with_capacity(work.len());
    for w in work {
        let Work {
            idx,
            level,
            previous,
            own,
            ground,
            carry,
            younger,
            hollow,
            reached,
            ponded,
            lobe,
            lake,
            share: margin_share,
            margin,
            ..
        } = w;
        if younger {
            stats[3] += 1;
            // Inside the rim of a crater younger than the flood the base is the crater
            // stamped on the lava, right as it is in relief: the crater's restore removes
            // what the pass laid there, and beyond the rim it buries the pass's relief by
            // its ejecta's thickness, so that the flood's sheets and gaps fade under the
            // ejecta rather than stop at a ring (M15). What the crater dug into is the
            // lava where the flood's flows reached (F94), each flood over the last in
            // their order as beside the crater, and where none did, the surface the floods
            // found.
            let flowed = match &lobed {
                Some(l) => l.get(&idx).copied().unwrap_or(0.0) >= 0.5,
                None => previous.is_some(),
            };
            if flowed {
                cells.push(CellWrite {
                    idx: idx as u32,
                    bedrock: bedrock[idx],
                    surface: true,
                    age: time,
                    unit: UNIT_LAVA,
                    deep: deep[idx],
                    melt: 0.0,
                    cover: 1.0,
                });
                compositions.push(rec.composition);
            } else if let Some((prev_unit, prev_age, prev_composition)) = previous
                && !smear.touched[idx]
                && units[idx] == UNIT_LAVA
            {
                cells.push(CellWrite {
                    idx: idx as u32,
                    bedrock: bedrock[idx],
                    surface: true,
                    age: prev_age,
                    unit: prev_unit,
                    deep: deep[idx],
                    melt: melt[idx],
                    cover: cover[idx],
                });
                compositions.push(prev_composition);
            }
            continue;
        }
        // Outside the low the loop filled, a lobe of a flood drawn as its flows lays a
        // sheet on the ground it crosses rather than filling to the low's level, unless
        // the lava ponds in a hollow there; outside its lakes, a soft flood's lobes lay
        // their sheets wherever they go.
        let level = if soft && !lake && !ponded {
            ground + FLOOD_SHEET_THICKNESS as f32
        } else if lobed.is_some() && previous.is_none() && !ponded {
            level.min(ground + FLOOD_SHEET_THICKNESS as f32)
        } else {
            level
        };
        if watched == Some(idx) {
            eprintln!(
                "watch flood {} t {:.4}: level {:.1} own {:.1} ground {:.1} hollow {} ponded {} previous {:?} touched {} unit {} lobe {:?} reached {:?} lake {}",
                rec.id,
                rec.time_gyr,
                level,
                own,
                ground,
                hollow,
                ponded,
                previous,
                smear.touched[idx],
                units[idx],
                lobed.as_ref().map(|l| l.get(&idx).copied().unwrap_or(0.0)),
                reached,
                lake,
            );
        }
        // A flow drapes what stands less than its own thickness above its level.
        let drape = FLOOD_SHEET_THICKNESS as f32;
        if ground < level + drape && (reached || margin) {
            stats[0] += 1;
            // A soft flood's lobe margin (F117) covers its cell by the lobe's share and
            // lays that share of its sheet; a lake or a pond it fills whole, a lobe running
            // into a lake keeping its sheet until the lake is deeper. A hard flood's lobes
            // and outline lay the share of the cell their margin covers (F131).
            let share = if soft && !ponded && !lake {
                lobe.clamp(0.0, 1.0)
            } else if ponded || lake || hollow {
                1.0
            } else {
                margin_share
            };
            let fill = if lake {
                (level - ground).max(FLOOD_SHEET_THICKNESS as f32 * lobe.clamp(0.0, 1.0))
            } else {
                (level - ground).max(0.0) * share
            };
            // A high the flow drapes carries a veneer, the sheet's thickness less
            // what the high stands above the level: thin lava over it, but lava,
            // not the ground it hides (the ghost craters of the maria show in
            // relief, not in albedo).
            let sheet = fill.max((drape - (ground - level).max(0.0)) * share);
            let surface = ground + (1.0 - LAVA_SUBSIDENCE as f32) * fill;
            compositions.push(rec.composition);
            // Its time and unit where it holds the cell, the ground's where it lays only
            // its margin's share.
            let held = soft || share >= 0.5;
            cells.push(CellWrite {
                idx: idx as u32,
                bedrock: surface + carry,
                surface: true,
                age: if held { time } else { ages[idx] },
                unit: if held { UNIT_LAVA } else { units[idx] },
                deep: deep[idx],
                melt: if held { 0.0 } else { melt[idx] * (1.0 - share) },
                cover: if soft {
                    share
                } else {
                    (sheet / LAVA_COVER_DEPTH as f32).clamp(0.0, 1.0)
                },
            });
        } else {
            if ground >= level + drape {
                stats[1] += 1;
            } else {
                stats[2] += 1;
            }
            // Not reached: take the smear off. Inside the coarse footprint the
            // base says lava where no pass has written yet; give it back the
            // surface the flood found, unless a crater formed since then covers the
            // cell, whose coarse cell the loop's record of that surface missed.
            let (unit, age, composition) = match previous {
                Some((u, a, c))
                    if !smear.touched[idx]
                        && units[idx] == UNIT_LAVA
                        && a >= smear.crater_time[idx] =>
                {
                    (u, a, c)
                }
                _ => (units[idx], ages[idx], found[idx]),
            };
            compositions.push(composition);
            cells.push(CellWrite {
                idx: idx as u32,
                bedrock: ground + carry,
                surface: unit != units[idx] || age != ages[idx],
                age,
                unit,
                deep: deep[idx],
                melt: melt[idx],
                cover: cover[idx],
            });
        }
    }
    debug_assert_eq!(cells.len(), compositions.len());
    (
        StampOutput {
            cells,
            composition: Some(rec.composition),
            stain: None,
            steep: None,
            ejecta_age: None,
            compositions: Some(compositions),
        },
        taken,
        stats,
    )
}

/// What a flood's pass needs to draw it as its flows (F94): the ground its lobes are
/// traced over, the body's radius and the seed of its lobes.
#[derive(Clone, Copy)]
pub struct FloodFlows<'a> {
    pub ground: &'a PaddedField,
    pub body_radius: f64,
    pub seed: u64,
    /// How far past the cells a flood's pass visits its lava may follow a hollow the
    /// loop never saw, m: the radius of the largest crater the loop left to the finish.
    pub hollow_reach: f64,
    /// Whether a lobe's margin covers its cells by its share (F117): Io's flows, whose
    /// colours the coverage sets, fade out at their edges as they did (F85); elsewhere
    /// a lava unit covers a cell whole where its lobes reach half of it.
    pub soft: bool,
}

/// Which floods can show at the end (F94): those among the last FLOOD_VISIBLE_DEPTH floods
/// to cover at least one coarse cell of their footprint. A successor's lobes leave gaps,
/// so a flood can show through the last few; the rest lie hidden.
pub fn visible_floods(floods: &[FloodRecord]) -> Vec<bool> {
    let mut order: Vec<usize> = (0..floods.len()).collect();
    order.sort_by(|&a, &b| {
        floods[a]
            .time_gyr
            .total_cmp(&floods[b].time_gyr)
            .then(a.cmp(&b))
    });
    let mut last: HashMap<u32, Vec<usize>> = HashMap::new();
    for &k in &order {
        // Each cell once, however many of the flood's lakes covered it.
        for idx in floods[k].cells().into_keys() {
            let stack = last.entry(idx).or_default();
            stack.push(k);
            if stack.len() > FLOOD_VISIBLE_DEPTH {
                stack.remove(0);
            }
        }
    }
    let mut visible = vec![false; floods.len()];
    for stack in last.values() {
        for &k in stack {
            visible[k] = true;
        }
    }
    visible
}

/// A lobe's area, m²: its length times its mean width.
pub(crate) fn lobe_area(lobe: &Lobe) -> f64 {
    (lobe.points.len().max(1) - 1) as f64 * lobe.spacing * 2.0 * lobe.half_width
}

/// The flows of one flood at the finishing resolution (F94): for each run of it, a sheet
/// down a slope or a lake (F109), a compound flow field from the run's cell nearest the
/// vent, where its lava came in, its trunk heading across the run toward its middle,
/// until its lobes cover the run FLOOD_LOBE_COVER times over. The lobes stay within two
/// coarse cells of the flood's own cells, inside the cap its pass may write: bound by the
/// cap alone, a cascade's flows, starting at its vent on the cap's edge, stopped before
/// they began.
pub fn flood_lobes(
    rec: &FloodRecord,
    coarse: &CubeSphere,
    grid: &CubeSphere,
    body_radius: f64,
    ground: &PaddedField,
    seed: u64,
    compound: f64,
) -> Vec<Lobe> {
    let covered = rec.cells();
    if covered.is_empty() {
        return Vec::new();
    }
    let r = body_radius;
    let cell = grid.delta() * r;
    let coarse_cell = coarse.delta() * r;
    // The flood's coarse cells and those within two of them, each with how far it lies
    // inside that area's edge, in coarse cells: a lobe goes where its margin stays inside.
    let mut area: HashSet<usize> = covered.keys().map(|&i| i as usize).collect();
    for _ in 0..2 {
        let ring: Vec<usize> = area
            .iter()
            .flat_map(|&i| {
                let c = coarse.cell(i);
                NEIGHBOURS
                    .iter()
                    .map(move |&(di, dj)| coarse.index(coarse.neighbor(c, di, dj)))
            })
            .filter(|n| !area.contains(n))
            .collect();
        area.extend(ring);
    }
    let mut inside: HashMap<usize, u32> = HashMap::with_capacity(area.len());
    let mut front: Vec<usize> = area
        .iter()
        .copied()
        .filter(|&i| {
            let c = coarse.cell(i);
            NEIGHBOURS
                .iter()
                .any(|&(di, dj)| !area.contains(&coarse.index(coarse.neighbor(c, di, dj))))
        })
        .collect();
    for &i in &front {
        inside.insert(i, 0);
    }
    let mut depth = 0;
    while !front.is_empty() {
        depth += 1;
        let mut next = Vec::new();
        for &i in &front {
            let c = coarse.cell(i);
            for (di, dj) in NEIGHBOURS {
                let n = coarse.index(coarse.neighbor(c, di, dj));
                if area.contains(&n) && !inside.contains_key(&n) {
                    inside.insert(n, depth);
                    next.push(n);
                }
            }
        }
        front = next;
    }
    let outside = |q: [f64; 3], hw: f64| {
        let margin = hw * (1.0 + IO_FLOW_LOBING) + (IO_FLOW_EDGE * hw).max(1.5 * cell);
        match inside.get(&coarse.index(coarse.locate(q).0)) {
            Some(&d) => (d as f64 + 0.5) * coarse_cell < margin,
            None => true,
        }
    };
    let mix = |a: f64, b: f64| a + (b - a) * compound;
    let rules = |half_width: f64| PathRules {
        min_spacing: 0.5 * cell,
        slope_follow: IO_FLOW_SLOPE_FOLLOW,
        slope_scale: IO_FLOW_SLOPE_SCALE,
        steer_length: mix(FLOOD_STEER_LENGTH, IO_FLOW_STEER_LENGTH),
        wander_rad: mix(FLOOD_WANDER_RAD, IO_FLOW_WANDER_RAD),
        wander_length: mix(FLOOD_WANDER_WIDTHS * half_width, IO_FLOW_WANDER_LENGTH),
        fill_depth: FLOOD_SHEET_THICKNESS,
        width_wobble: IO_FLOW_WIDTH_WOBBLE,
    };
    let mut lobes = Vec::new();
    for (k, pond) in rec.ponds.iter().enumerate() {
        let cells: Vec<[f64; 3]> = pond
            .cells
            .iter()
            .map(|c| coarse.center(coarse.cell(c.idx as usize)))
            .collect();
        let Some(&entry) = cells
            .iter()
            .min_by(|a, b| angle_between(rec.vent, **a).total_cmp(&angle_between(rec.vent, **b)))
        else {
            continue;
        };
        let middle = {
            let sum = cells
                .iter()
                .fold([0.0; 3], |s, d| [s[0] + d[0], s[1] + d[1], s[2] + d[2]]);
            if sum[0] * sum[0] + sum[1] * sum[1] + sum[2] * sum[2] < 1e-18 {
                entry
            } else {
                normalize(sum)
            }
        };
        let run_area: f64 = pond
            .cells
            .iter()
            .map(|c| coarse.cell_area(coarse.cell(c.idx as usize)) * r * r)
            .sum();
        let reach = cells
            .iter()
            .map(|&d| angle_between(entry, d))
            .fold(0.0, f64::max)
            * r
            + coarse_cell;
        let id = rec.id ^ ((k as u64) << 52);
        let rng = Rng::new(seed ^ id ^ 0xF100_D10B);
        // The trunk heads across the run toward its middle.
        let toward = {
            let v = [
                middle[0] - entry[0],
                middle[1] - entry[1],
                middle[2] - entry[2],
            ];
            let radial = dot(v, entry);
            let t = [
                v[0] - radial * entry[0],
                v[1] - radial * entry[1],
                v[2] - radial * entry[2],
            ];
            if t[0] * t[0] + t[1] * t[1] + t[2] * t[2] > 1e-18 {
                normalize(t)
            } else {
                // Entered at the middle: any heading.
                let a = if entry[2].abs() < 0.9 {
                    [0.0, 0.0, 1.0]
                } else {
                    [1.0, 0.0, 0.0]
                };
                rotate_about(
                    entry,
                    normalize(crate::geom::cross(entry, a)),
                    std::f64::consts::TAU * rng.uniform(2, 0),
                )
            }
        };
        lobes.extend(flow_field(
            grid,
            r,
            ground,
            &rng,
            id,
            entry,
            toward,
            reach,
            mix(FLOOD_LOBE_COVER, IO_FLOW_OVERLAP) * run_area,
            0.0,
            compound,
            &rules,
            &outside,
        ));
    }
    lobes
}

/// A compound flow field (F94): from its source a trunk along a heading, a share of the
/// field's reach long, up to FLOOD_BRANCHES.1 branches off it, and breakouts on the lobes
/// already laid, at their fronts and along them, until the lobes cover `target` m². Each
/// lobe is a path over the ground (F59), a share of the reach long and a share of its own
/// length wide, wandering at its width's scale and ponding where it would climb more than
/// a sheet; the trunk widens by `trunk_toe` from its source, a breakout by IO_FLOW_TOE.
/// A path stops where `outside` says its margin, at its half width, would leave the
/// field's bounds.
#[allow(clippy::too_many_arguments)]
pub(crate) fn flow_field(
    grid: &CubeSphere,
    body_radius: f64,
    ground: &PaddedField,
    rng: &Rng,
    id: u64,
    source: [f64; 3],
    heading: [f64; 3],
    reach: f64,
    target: f64,
    trunk_toe: f64,
    compound: f64,
    rules: &dyn Fn(f64) -> PathRules,
    outside: &dyn Fn([f64; 3], f64) -> bool,
) -> Vec<Lobe> {
    let r = body_radius;
    let cell = grid.delta() * r;
    // A compound field's breakouts are short and broad, a sheet flow's lobes long
    // (F117): its shares by how compound the eruption is.
    let blend =
        |a: (f64, f64), b: (f64, f64)| (a.0 + (b.0 - a.0) * compound, a.1 + (b.1 - a.1) * compound);
    // The breakouts' half widths as shares of their lengths.
    let (w_lo, w_hi) = blend(
        FLOOD_LOBE_WIDTH,
        (
            0.5 * IO_FLOW_BREAKOUT_ASPECT.0,
            0.5 * IO_FLOW_BREAKOUT_ASPECT.1,
        ),
    );
    let width = |index: u64, length: f64| {
        (length * (w_lo + (w_hi - w_lo) * rng.uniform(index, 0))).max(1.5 * cell)
    };
    // The trunk and its branches: a sheet flow's as its lobes are, a compound field's the
    // narrow insulated paths its breakouts spread from (F85).
    let (f_lo, f_hi) = FLOOD_LOBE_WIDTH;
    let path_width = |index: u64, length: f64| {
        let sheet = length * (f_lo + (f_hi - f_lo) * rng.uniform(index, 0));
        (sheet + (0.5 * IO_FLOW_TRUNK_WIDTH - sheet) * compound).max(1.5 * cell)
    };
    let (t_lo, t_hi) = FLOOD_TRUNK_SHARE;
    let trunk_length = reach * (t_lo + (t_hi - t_lo) * rng.uniform(3, 0));
    let hw = path_width(4, trunk_length);
    let mut lobes = vec![
        trace_flow(
            grid,
            r,
            &rules(hw),
            id,
            5,
            source,
            heading,
            trunk_length,
            hw,
            trunk_toe,
            ground,
            &|q| outside(q, hw),
            0.0,
        )
        .0,
    ];
    let mut area = lobe_area(&lobes[0]);
    // Branches from the trunk.
    let (b_lo, b_hi) = FLOOD_BRANCHES;
    let branches = b_lo + (rng.bits(6, 0) % (b_hi - b_lo + 1) as u64) as u32;
    for b in 0..branches as u64 {
        let trunk = &lobes[0];
        if trunk.points.len() < 4 {
            break;
        }
        let s = 0.25 + 0.55 * rng.uniform(10 + b, 0);
        let k = ((trunk.points.len() - 1) as f64 * s) as usize;
        let (p, tangent) = path_point(&trunk.points, k);
        let side = if rng.uniform(10 + b, 1) < 0.5 {
            -1.0
        } else {
            1.0
        };
        let turn = side * (20f64 + 30.0 * rng.uniform(10 + b, 2)).to_radians();
        let remaining = (trunk.points.len() - 1 - k) as f64 * trunk.spacing;
        let length = remaining * (0.4 + 0.5 * rng.uniform(10 + b, 3));
        if length < 4.0 * cell {
            continue;
        }
        let hw = path_width(20 + b, length);
        let heading = rotate_about(p, tangent, turn);
        let lobe = trace_flow(
            grid,
            r,
            &rules(hw),
            id,
            30 + b,
            p,
            heading,
            length,
            hw,
            0.0,
            ground,
            &|q| outside(q, hw),
            0.0,
        )
        .0;
        area += lobe_area(&lobe);
        lobes.push(lobe);
    }
    // Breakouts on the lobes already laid, until they cover the target.
    let (l_lo, l_hi) = blend(FLOOD_LOBE_LENGTH, IO_FLOW_BREAKOUT_LENGTH);
    let mut j = 0u64;
    while area < target && lobes.len() < FLOOD_LOBES_MAX && j < 4 * FLOOD_LOBES_MAX as u64 {
        let base = 1_000 + 16 * j;
        j += 1;
        let lengths: Vec<f64> = lobes
            .iter()
            .map(|l| (l.points.len().max(1) - 1) as f64 * l.spacing)
            .collect();
        let total: f64 = lengths.iter().sum();
        if total <= 0.0 {
            break;
        }
        let mut pick = rng.uniform(base, 0) * total;
        let mut which = 0;
        for (n, len) in lengths.iter().enumerate() {
            which = n;
            if pick < *len {
                break;
            }
            pick -= len;
        }
        let path = &lobes[which];
        if path.points.len() < 2 {
            continue;
        }
        let s = if rng.uniform(base, 1) < IO_FLOW_FRONT_SHARE {
            1.0
        } else {
            rng.uniform(base, 2).sqrt()
        };
        let k = ((path.points.len() - 1) as f64 * s).round() as usize;
        let (p, tangent) = path_point(&path.points, k);
        let heading = rotate_about(p, tangent, IO_FLOW_BREAKOUT_SPREAD * rng.normal(base, 3));
        let length = reach * l_lo * (l_hi / l_lo).powf(rng.uniform(base, 4));
        let hw = width(base + 5, length);
        let lobe = trace_flow(
            grid,
            r,
            &rules(hw),
            id,
            base + 6,
            p,
            heading,
            length,
            hw,
            IO_FLOW_TOE,
            ground,
            &|q| outside(q, hw),
            0.0,
        )
        .0;
        area += lobe_area(&lobe);
        lobes.push(lobe);
    }
    lobes
}

/// The share of each cell a flood's lobes cover (F94), the strongest where they cross: a
/// tube about each lobe's path, as Io's flows are drawn (F59), its margin lobate at its
/// width's scale and soft over a share of its width or a cell and a half.
pub(crate) fn lobe_cover(
    lobes: &[Lobe],
    grid: &CubeSphere,
    body_radius: f64,
    id: u64,
    stamper: &mut Stamper,
) -> HashMap<usize, f32> {
    let cell = grid.delta() * body_radius;
    let noise = GradientNoise::new(id ^ 0xF10D_0000);
    let mut cover: HashMap<usize, f32> = HashMap::new();
    for lobe in lobes {
        if lobe.points.len() < 2 {
            continue;
        }
        let widest = lobe.widths.iter().cloned().fold(0.0f32, f32::max) as f64;
        let reach = lobe.half_width * widest * (1.0 + IO_FLOW_LOBING)
            + (IO_FLOW_EDGE * lobe.half_width * widest).max(1.5 * cell);
        stamper.gather_corridor(grid, body_radius, &lobe.points, reach);
        let last = (lobe.widths.len() - 1) as f64;
        for &(idx, r, at) in &stamper.footprint {
            let at = at.clamp(0.0, last);
            let i = at.floor() as usize;
            let f = at - i as f64;
            let w = lobe.widths[i] as f64 * (1.0 - f)
                + lobe.widths[(i + 1).min(lobe.widths.len() - 1)] as f64 * f;
            let hw = lobe.half_width * w;
            let edge = (IO_FLOW_EDGE * hw).max(1.5 * cell);
            let d = grid.center(grid.cell(idx));
            let lobing = 1.0 + IO_FLOW_LOBING * noise.fbm(d, body_radius / (1.5 * hw), 2, 2.0, 0.5);
            let c = ((hw * lobing - r) / edge).clamp(0.0, 1.0) as f32;
            if c > 0.0 {
                let e = cover.entry(idx).or_insert(0.0);
                if c > *e {
                    *e = c;
                }
            }
        }
    }
    cover
}

/// Area fraction of the surface carrying a lava unit.
pub fn lava_fraction(layers: &LayerSet) -> f64 {
    let grid = layers.grid;
    let units = layers.unit_id.as_slice();
    let mut lava = 0.0;
    let mut total = 0.0;
    for (k, &unit) in units.iter().enumerate() {
        let a = grid.cell_area(grid.cell(k));
        total += a;
        if is_lava_unit(unit) {
            lava += a;
        }
    }
    lava / total
}

/// Whether a unit is lava or something formed on it.
pub fn is_lava_unit(unit: u16) -> bool {
    matches!(unit, UNIT_LAVA | UNIT_LAVA_CRATER_FLOOR | UNIT_LAVA_EJECTA)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::presets::preset;

    #[test]
    fn a_light_crust_holds_a_dense_magma_down() {
        // F125: the Moon's iron-rich magma, denser than its anorthositic crust, loses through
        // the crust much of the buoyancy it gains below it while the lid is thin, less as the
        // lid thickens; Mercury's iron-poor magma, lighter than its crust, rises through it.
        let moon = |depth: f64| crust_buoyancy(2_905.0, 2_750.0, 40e3, depth);
        assert!(moon(60e3) < 0.5 && moon(200e3) > moon(60e3) && moon(200e3) < 1.0);
        assert_eq!(crust_buoyancy(2_705.0, 2_871.0, 40e3, 80e3), 1.0);
    }

    fn process(grid: CubeSphere) -> FloodProcess {
        let (g, _) = preset("Moon").unwrap().to_graph();
        let thermal = ThermalHistory::from_graph(&g);
        let mut p = FloodProcess::new(&g, grid, &thermal, &[]).unwrap();
        p.body_radius = 1_000e3;
        p.areas = (0..grid.cell_count())
            .map(|k| grid.cell_area(grid.cell(k)) * 1e12)
            .collect();
        p
    }

    /// Diagnostic: the eruptive history of a preset, run with `--ignored --nocapture`.
    #[test]
    #[ignore]
    fn print_eruptive_history() {
        for name in ["Moon", "Mercury", "Mars", "Venus", "Io"] {
            let (g, _) = preset(name).unwrap().to_graph();
            let thermal = ThermalHistory::from_graph(&g);
            let grid = CubeSphere::new(8);
            let Some(p) = FloodProcess::new(&g, grid, &thermal, &[]) else {
                println!("{name}: its volcanism is handled elsewhere");
                continue;
            };
            println!(
                "{name}: t  lithosphere_km  elastic_km  heat_mW  magma_m3_s  strain  eruptive  budget_km3_per_0.1Gyr"
            );
            let mut t = 0.0;
            while t < g.f64("body.age") {
                let s = thermal.at(t);
                println!(
                    "  {t:.1}  {:.0}  {:.0}  {:.0}  {:.2}  {:.2e}  {:.3}  {:.3e}",
                    s.lithosphere / 1e3,
                    s.elastic / 1e3,
                    s.heat_flux * 1e3,
                    s.magma_rate,
                    2.0 * s.contraction / p.body_radius,
                    p.eruptive_fraction(t),
                    p.budget(t, t + 0.1) / 1e9
                );
                t += 0.5;
            }
        }
    }

    /// From a vent a few cells up a slope the lava runs down to a pit, leaving a flow's sheet
    /// on each cell it crosses (M15), and the pit takes the rest; the volume is kept.
    #[test]
    fn lava_runs_downhill_leaving_its_sheet_and_fills_the_pit_with_the_rest() {
        let grid = CubeSphere::new(32);
        let mut p = process(grid);
        let mut layers = LayerSet::new(grid);
        let cell = |i, j| grid.index(crate::grid::Cell { face: 0, i, j });
        // A pit of one cell, 1 km deep, and a vent four cells away up the slope to it.
        let (pit, vent) = (cell(16, 16), cell(20, 16));
        for k in 0..grid.cell_count() {
            let d =
                crate::geom::angle_between(grid.center(grid.cell(k)), grid.center(grid.cell(pit)));
            layers.bedrock.as_mut_slice()[k] = (d * 1e4) as f32;
        }
        layers.bedrock.as_mut_slice()[pit] = -1_000.0;
        let volume = 0.5 * 1_000.0 * p.areas[pit];
        let rec = p.erupt(&mut layers, vent, volume, 1.0, 1);
        let covered = rec.cells();
        // The way down: a sheet a flow thick on each cell from the vent to the pit, and on
        // the cells beside it the flow need not climb to; the pit holds the rest as a lake.
        for i in 17..=20 {
            let cov = covered
                .get(&(cell(i, 16) as u32))
                .expect("a cell on the way has no lava");
            assert!(cov.sheet && cov.level.is_nan(), "{cov:?}");
        }
        let sheet = (1.0 - LAVA_SUBSIDENCE) * FLOOD_SHEET_THICKNESS;
        let mut on_the_way = 0.0;
        for (&idx, cov) in &covered {
            if idx as usize != pit {
                assert!(cov.sheet && cov.thickness as f64 <= sheet + 1e-3, "{cov:?}");
                on_the_way +=
                    cov.thickness as f64 / (1.0 - LAVA_SUBSIDENCE) * p.areas[idx as usize];
            }
        }
        let pond = covered.get(&(pit as u32)).expect("pit not flooded");
        assert!(!pond.sheet);
        let rest = (volume - on_the_way) / p.areas[pit];
        let expected = -1_000.0 + (1.0 - LAVA_SUBSIDENCE) * rest;
        let b = layers.bedrock.as_slice()[pit];
        assert!((b as f64 - expected).abs() < 5.0, "{b} vs {expected}");
        assert_eq!(layers.unit_id.as_slice()[pit], UNIT_LAVA);
        assert_eq!(layers.surface_age.as_slice()[pit], 1.0);
        // What stands is the lava less its subsidence; the record keeps the volume.
        let filled: f64 = covered
            .iter()
            .map(|(&idx, c)| c.thickness as f64 * p.areas[idx as usize])
            .sum();
        let standing = (1.0 - LAVA_SUBSIDENCE) * volume;
        assert!(
            (filled - standing).abs() < 0.02 * volume,
            "{filled} vs {standing}"
        );
        assert!(
            (rec.volume - volume).abs() < 0.02 * volume,
            "{} vs {volume}",
            rec.volume
        );
    }

    /// Two pits joined by a pass: the lake fed at the shallower one crosses the pass once
    /// the level reaches it and the deeper pit takes the rest, so the first pit ends at
    /// or below the pass and the deeper one holds lava.
    #[test]
    fn a_full_pit_spills_into_the_next_one() {
        let grid = CubeSphere::new(32);
        let mut p = process(grid);
        let mut layers = LayerSet::new(grid);
        let cell = |i, j| grid.index(crate::grid::Cell { face: 0, i, j });
        let (a, b) = (cell(16, 16), cell(16, 20));
        {
            let h = layers.bedrock.as_mut_slice();
            for (k, v) in h.iter_mut().enumerate() {
                let da = crate::geom::angle_between(
                    grid.center(grid.cell(k)),
                    grid.center(grid.cell(a)),
                );
                let db = crate::geom::angle_between(
                    grid.center(grid.cell(k)),
                    grid.center(grid.cell(b)),
                );
                *v = (da.min(db) * 1e4) as f32;
            }
            h[a] = -500.0;
            h[b] = -2_000.0;
        }
        // The pass towards b is two cells out; the pond around a holds everything
        // below it within that distance.
        let (pass, capacity) = {
            let h = layers.bedrock.as_slice();
            let pass = h[cell(16, 18)];
            let mut capacity = 0.0;
            for i in 14..=18 {
                for j in 14..=18 {
                    let k = cell(i, j);
                    if h[k] < pass {
                        capacity += (pass - h[k]) as f64 * p.areas[k];
                    }
                }
            }
            (pass, capacity)
        };
        let volume = capacity + 1_000.0 * p.areas[a];
        let rec = p.erupt(&mut layers, a, volume, 2.0, 2);
        let covered = rec.cells();
        assert!(covered.contains_key(&(a as u32)) && covered.contains_key(&(b as u32)));
        let h = layers.bedrock.as_slice();
        assert!(
            h[a] > -500.0 && h[a] <= pass + 1.0,
            "a at {} vs pass {pass}",
            h[a]
        );
        assert!(h[b] > -2_000.0 && h[b] < pass, "b at {}", h[b]);
        let filled: f64 = covered.values().map(|c| c.thickness as f64).sum::<f64>();
        assert!(filled > 0.0);
        assert!(
            (rec.volume - volume).abs() < 0.02 * volume,
            "{} vs {volume}",
            rec.volume
        );
    }

    /// A lake stands no more than a flow's thickness above the open ground beyond it
    /// (audit W7): over rough ground, one that fills to its pass spills the rest on into
    /// the next low rather than standing level over the ground beyond the pass, so each
    /// flood ends in flow fronts, not cliffs.
    #[test]
    fn a_flood_over_rough_ground_ends_in_flow_fronts_not_cliffs() {
        let grid = CubeSphere::new(32);
        let mut p = process(grid);
        let mut layers = LayerSet::new(grid);
        let noise = GradientNoise::new(7);
        for (k, h) in layers.bedrock.as_mut_slice().iter_mut().enumerate() {
            *h = (3_000.0 * noise.fbm(grid.center(grid.cell(k)), 4.0, 4, 2.0, 0.5)) as f32;
        }
        let rng = Rng::new(11);
        let mut spilled = 0;
        for e in 0..60u64 {
            let vent = grid.index(grid.locate(rng.unit_direction(e, 0)).0);
            let before = layers.bedrock.as_slice().to_vec();
            let volume = 1e11 * 100f64.powf(rng.uniform(e, 1));
            let rec = p.erupt(&mut layers, vent, volume, 1.0, e);
            if rec.ponds.len() > 1 {
                spilled += 1;
            }
            let covered = rec.cells();
            for (&idx, c) in &covered {
                let t = c.thickness as f64 / (1.0 - LAVA_SUBSIDENCE);
                let cell = grid.cell(idx as usize);
                for (di, dj) in NEIGHBOURS {
                    let n = grid.index(grid.neighbor(cell, di, dj));
                    if covered.contains_key(&(n as u32)) {
                        continue;
                    }
                    let front = (c.level as f64 - before[n] as f64).min(t);
                    assert!(
                        front <= FLOOD_SHEET_THICKNESS + 1.0,
                        "flood {e}: lava {t:.0} m thick at level {} stands {front:.0} m over the ground beyond",
                        c.level
                    );
                }
            }
        }
        assert!(spilled > 0, "no lake spilled");
    }

    /// A rise's shields take the share of its eruptions that those started hold (F108):
    /// what the rise adds to the vent weights goes by that share to its shields, out of the
    /// floods' weights, and the weights away from the rise stay as they were.
    #[test]
    fn a_rise_s_shields_take_its_eruptions_as_they_start() {
        let grid = CubeSphere::new(16);
        let n = grid.cell_count();
        let mut p = process(grid);
        let layers = LayerSet::new(grid);
        let center = normalize([1.0, 0.2, 0.1]);
        let sigma = 0.5;
        let angle = |k: usize| angle_between(center, grid.center(grid.cell(k)));
        let own: Vec<f32> = (0..n)
            .map(|k| (-(angle(k) / sigma).powi(2)).exp() as f32)
            .collect();
        let dome_at = |k: usize| crate::stages::plume::dome(angle(k) / sigma);
        let boost: Vec<f64> = (0..n)
            .map(|k| 1.0 + PLUME_VENT_BOOST * dome_at(k))
            .collect();
        let dome_cells: Vec<(usize, f32)> = (0..n)
            .filter(|&k| dome_at(k) > 0.0)
            .map(|k| (k, dome_at(k) as f32))
            .collect();
        let per_cell = |c: Vec<f64>| {
            let mut prev = 0.0;
            c.into_iter()
                .map(|v| {
                    let w = v - prev;
                    prev = v;
                    w
                })
                .collect::<Vec<f64>>()
        };
        // The weights with no rise at all: the background a rise adds to.
        p.province_weight = vec![MAGMA_PROVINCE_FLOOR; n];
        let background = per_cell(p.vent_weights(&layers, None, &[]).0);
        p.province_weight = own
            .iter()
            .map(|&g| MAGMA_PROVINCE_FLOOR + g as f64)
            .collect();
        p.rise_province = vec![own];
        let weights = |capture: f64| {
            let rise = RiseVent {
                dome: dome_cells.clone(),
                capture,
            };
            let (cumulative, kept) = p.vent_weights(&layers, Some(&boost), &[rise]);
            (per_cell(cumulative), kept)
        };
        let ((before, whole), (after, kept)) = (weights(0.0), weights(0.6));
        let total = |w: &[f64]| w.iter().sum::<f64>();
        assert!((whole - 1.0).abs() < 1e-12);
        let left = total(&before) - 0.6 * (total(&before) - total(&background));
        assert!((total(&after) - left).abs() < 1e-9 * left);
        assert!((kept - left / total(&before)).abs() < 1e-9);
        let far = (0..n)
            .max_by(|&a, &b| angle(a).total_cmp(&angle(b)))
            .unwrap();
        assert!((after[far] - before[far]).abs() < 1e-9 * before[far]);
        let flank = (0..n)
            .min_by(|&a, &b| (angle(a) - 0.25).abs().total_cmp(&(angle(b) - 0.25).abs()))
            .unwrap();
        let stays = before[flank] - 0.6 * (before[flank] - background[flank]);
        assert!((after[flank] - stays).abs() < 1e-9 * before[flank]);
        assert!(
            after[flank] < 0.45 * before[flank],
            "{} vs {}",
            after[flank],
            before[flank]
        );
    }

    /// A flood over flat ground: a disc of coarse cells about a vent.
    fn disc_flood(
        coarse: &CubeSphere,
        center: [f64; 3],
        radius: f64,
        time: f64,
        id: u64,
    ) -> FloodRecord {
        let cells: Vec<FloodCell> = (0..coarse.cell_count())
            .filter(|&k| angle_between(center, coarse.center(coarse.cell(k))) < radius)
            .map(|k| FloodCell {
                idx: k as u32,
                thickness: 50.0,
                prev_unit: UNIT_PRIMORDIAL,
                prev_age: 0.0,
                prev_composition: 0.0,
            })
            .collect();
        let area = cells
            .iter()
            .map(|c| coarse.cell_area(coarse.cell(c.idx as usize)) * 1e12)
            .sum();
        FloodRecord {
            vent: center,
            volume: 50.0 * area,
            time_gyr: time,
            area,
            ponds: vec![Pond {
                level: 50.0,
                cells,
                sheet: false,
            }],
            composition: 0.0,
            id,
            effusion: FLOOD_EFFUSION_RATE.0,
        }
    }

    #[test]
    fn a_flood_fills_a_crater_it_comes_into_and_not_one_whose_rim_holds_it_out() {
        // F95: a crater finer than the loop's cells in a flood's low. Lava lying round it
        // pours over a rim the flood breaches and fills the bowl to the level; a rim
        // standing above the lava keeps the bowl dry, whatever lobe crosses it.
        let coarse = CubeSphere::new(64);
        let grid = CubeSphere::new(512);
        let body_radius = 1_000e3;
        let center = normalize([0.3, 1.0, 0.2]);
        let flat = vec![0.0f32; grid.cell_count()];
        let ground = PaddedField::new(&grid, &flat);
        let crater = 15e3;
        let rise_at_floor = |rim: f32| -> f32 {
            let mut flood = disc_flood(&coarse, center, 0.12, 1.0, 7);
            flood.area = flood
                .cells()
                .keys()
                .map(|&i| coarse.cell_area(coarse.cell(i as usize)) * body_radius * body_radius)
                .sum();
            let footprint: HashSet<u32> = flood.cells().keys().copied().collect();
            let smear = LavaSmear::new(coarse, grid, std::slice::from_ref(&flood), None);
            let mut layers = LayerSet::new(grid);
            let mut relief = vec![0.0f32; grid.cell_count()];
            for (k, relief) in relief.iter_mut().enumerate() {
                let d = grid.center(grid.cell(k));
                let x = angle_between(center, d) * body_radius / crater;
                *relief = if x < 1.0 {
                    -2_000.0 * (1.0 - (x * x) as f32)
                } else if x < 1.6 {
                    rim * (1.0 - ((x - 1.0) / 0.6) as f32)
                } else {
                    0.0
                };
                layers.bedrock.as_mut_slice()[k] = smear.total[k] + *relief;
                // The base says lava over the loop's footprint, and the crater the finish
                // stamped on it says crater.
                let in_footprint = footprint.contains(&(coarse.index(coarse.locate(d).0) as u32));
                layers.unit_id.as_mut_slice()[k] = if x < 1.6 {
                    UNIT_EJECTA
                } else if in_footprint {
                    UNIT_LAVA
                } else {
                    UNIT_PRIMORDIAL
                };
            }
            let deform = Deformation::new(coarse, body_radius, 0.01, 4.5, &[], &[], &[]);
            let mut stamper = Stamper::new(&grid);
            let flows = FloodFlows {
                ground: &ground,
                body_radius,
                seed: 3,
                hollow_reach: 40e3,
                soft: false,
            };
            let (out, _, _) = flood_compute(
                &layers,
                &grid,
                &coarse,
                &flood,
                &smear,
                &deform,
                Some(flows),
                &smear.total,
                &mut stamper,
            );
            out.apply(&mut layers);
            let mut rises: Vec<f32> = (0..grid.cell_count())
                .filter(|&k| {
                    angle_between(center, grid.center(grid.cell(k))) * body_radius < 0.6 * crater
                })
                .map(|k| layers.bedrock.as_slice()[k] - relief[k])
                .collect();
            rises.sort_by(|a, b| a.total_cmp(b));
            rises[0]
        };
        let breached = rise_at_floor(0.0);
        assert!(
            breached > 1_000.0,
            "a breached crater's floor rose {breached} m"
        );
        let held = rise_at_floor(500.0);
        assert!(
            held.abs() < 1.0,
            "a crater whose rim holds the lava out rose {held} m"
        );
    }

    #[test]
    fn a_soft_flood_ponds_only_in_the_closed_lows_its_lobes_run_into() {
        // F117: Io's flood ponded at 50 m in the loop. On an open slope of the fine ground
        // that level is no lake's: its lobes lay their sheets over the ground, with no lake
        // perched below the level and no terrace where the units still say lava. A closed
        // bowl its lobes run into fills to where it spills.
        let coarse = CubeSphere::new(64);
        let grid = CubeSphere::new(512);
        let body_radius = 1_000e3;
        let center = normalize([0.3, 1.0, 0.2]);
        let across = normalize([center[1], -center[0], 0.0]);
        let thickness = |relief: &dyn Fn([f64; 3]) -> f32| -> Vec<(f32, f32)> {
            let flood = disc_flood(&coarse, center, 0.12, 1.0, 7);
            let footprint: HashSet<u32> = flood.cells().keys().copied().collect();
            let smear = LavaSmear::new(coarse, grid, std::slice::from_ref(&flood), None);
            let mut layers = LayerSet::new(grid);
            let ground: Vec<f32> = (0..grid.cell_count())
                .map(|k| relief(grid.center(grid.cell(k))))
                .collect();
            for (k, g) in ground.iter().enumerate() {
                let d = grid.center(grid.cell(k));
                layers.bedrock.as_mut_slice()[k] = smear.total[k] + g;
                layers.unit_id.as_mut_slice()[k] =
                    if footprint.contains(&(coarse.index(coarse.locate(d).0) as u32)) {
                        UNIT_LAVA
                    } else {
                        UNIT_PRIMORDIAL
                    };
            }
            let base = layers.bedrock.as_slice().to_vec();
            let padded = PaddedField::new(&grid, &ground);
            let deform = Deformation::new(coarse, body_radius, 0.01, 4.5, &[], &[], &[]);
            let mut stamper = Stamper::new(&grid);
            let flows = FloodFlows {
                ground: &padded,
                body_radius,
                seed: 3,
                hollow_reach: 40e3,
                soft: true,
            };
            let (out, _, _) = flood_compute(
                &layers,
                &grid,
                &coarse,
                &flood,
                &smear,
                &deform,
                Some(flows),
                &base,
                &mut stamper,
            );
            out.apply(&mut layers);
            ground
                .iter()
                .zip(layers.bedrock.as_slice())
                .map(|(&g, &b)| (g, b - g))
                .collect()
        };
        let sheet = ((1.0 - LAVA_SUBSIDENCE) * FLOOD_SHEET_THICKNESS) as f32;
        // A slope of 1 in 200 across the flood: 600 m down from its middle to its edge.
        let slope = thickness(&|d| (0.005 * body_radius * dot(d, across)) as f32);
        let thickest = slope.iter().map(|c| c.1).fold(0.0f32, f32::max);
        assert!(
            thickest > 0.5 * sheet && thickest < sheet + 0.5,
            "on a slope the lava is up to {thickest} m thick, its sheet {sheet} m"
        );
        let cut = slope.iter().map(|c| c.1).fold(0.0f32, f32::min);
        assert!(cut > -0.5, "the slope was cut {cut} m down");
        // A bowl 100 m deep and 15 km across at the vent, on flat ground.
        let bowl = thickness(&|d| {
            let x = angle_between(center, d) * body_radius / 15e3;
            (-100.0 * (-x * x).exp()) as f32
        });
        let (deepest, filled) = bowl
            .iter()
            .copied()
            .fold((0.0f32, 0.0f32), |a, c| if c.0 < a.0 { c } else { a });
        assert!(
            (filled - (1.0 - LAVA_SUBSIDENCE as f32) * -deepest).abs() < 1.0,
            "the bowl's floor at {deepest} m took {filled} m of lava"
        );
        let plain = bowl
            .iter()
            .filter(|c| c.0 > -1.0)
            .map(|c| c.1)
            .fold(0.0f32, f32::max);
        assert!(plain < sheet + 0.5, "the plain round it took {plain} m");
    }

    #[test]
    fn a_flood_shows_through_its_last_few_successors() {
        // F94: six floods over the same ground; only the last FLOOD_VISIBLE_DEPTH can show.
        let coarse = CubeSphere::new(16);
        let center = normalize([1.0, 0.2, 0.1]);
        let floods: Vec<FloodRecord> = (0..6)
            .map(|k| disc_flood(&coarse, center, 0.2, k as f64 * 0.1, k as u64 + 1))
            .collect();
        let visible = visible_floods(&floods);
        let shown = visible.iter().filter(|&&v| v).count();
        assert_eq!(shown, FLOOD_VISIBLE_DEPTH, "{visible:?}");
        assert!(visible[5] && !visible[0], "{visible:?}");
    }

    #[test]
    fn a_cascade_s_flows_cover_its_runs() {
        // F110: a flood running as a sheet from its vent down to a lake draws flows over
        // both runs, each from its cell nearest the vent, and they stay within two coarse
        // cells of the flood's own cells. Bound by the flood's cap and started at the vent,
        // on the cap's edge, they had stopped before they began.
        let coarse = CubeSphere::new(32);
        let grid = CubeSphere::new(128);
        let body_radius = 1_000e3;
        let flat = vec![0.0f32; grid.cell_count()];
        let ground = PaddedField::new(&grid, &flat);
        let vent = normalize([0.3, 1.0, 0.2]);
        let (east, north) = crate::geom::local_basis(vent);
        let step = coarse.delta();
        let mut path: Vec<u32> = Vec::new();
        for k in 0..14 {
            let p = crate::geom::offset(vent, east, north, k as f64 * step, 0.4);
            let c = coarse.locate(p).0;
            for (di, dj) in [(0, 0), (1, 0), (-1, 0), (0, 1), (0, -1)] {
                let n = coarse.index(coarse.neighbor(c, di, dj)) as u32;
                if !path.contains(&n) {
                    path.push(n);
                }
            }
        }
        let end = crate::geom::offset(vent, east, north, 16.0 * step, 0.4);
        let lake = disc_flood(&coarse, end, 3.0 * step, 1.0, 9).ponds.remove(0);
        let sheet = Pond {
            level: 0.0,
            cells: path
                .iter()
                .map(|&idx| FloodCell {
                    idx,
                    thickness: 47.5,
                    prev_unit: UNIT_PRIMORDIAL,
                    prev_age: 0.0,
                    prev_composition: 0.0,
                })
                .collect(),
            sheet: true,
        };
        let mut flood = FloodRecord {
            vent,
            volume: 0.0,
            time_gyr: 1.0,
            area: 0.0,
            ponds: vec![sheet, lake],
            composition: 0.0,
            id: 9,
            effusion: FLOOD_EFFUSION_RATE.0,
        };
        flood.area = flood
            .cells()
            .keys()
            .map(|&i| coarse.cell_area(coarse.cell(i as usize)) * body_radius * body_radius)
            .sum();
        let lobes = flood_lobes(&flood, &coarse, &grid, body_radius, &ground, 3, 0.0);
        let mut stamper = Stamper::new(&grid);
        let cover = lobe_cover(&lobes, &grid, body_radius, flood.id, &mut stamper);
        let ratio = grid.face_size() / coarse.face_size();
        let share = |pond: &Pond| {
            let (mut n, mut hit) = (0usize, 0usize);
            for c in &pond.cells {
                let cell = coarse.cell(c.idx as usize);
                for a in 0..ratio {
                    for b in 0..ratio {
                        let fine = crate::grid::Cell {
                            face: cell.face,
                            i: cell.i * ratio + a,
                            j: cell.j * ratio + b,
                        };
                        n += 1;
                        if cover.get(&grid.index(fine)).copied().unwrap_or(0.0) >= 0.5 {
                            hit += 1;
                        }
                    }
                }
            }
            hit as f64 / n as f64
        };
        let (on_sheet, on_lake) = (share(&flood.ponds[0]), share(&flood.ponds[1]));
        assert!(
            on_sheet > 0.2,
            "the flows cover {on_sheet:.2} of the sheet's run"
        );
        assert!(on_lake > 0.2, "the flows cover {on_lake:.2} of the lake");
        let centres: Vec<[f64; 3]> = flood
            .cells()
            .keys()
            .map(|&i| coarse.center(coarse.cell(i as usize)))
            .collect();
        for (&idx, &c) in &cover {
            if c < 0.5 {
                continue;
            }
            let d = grid.center(grid.cell(idx));
            let near = centres
                .iter()
                .map(|&q| angle_between(q, d))
                .fold(f64::INFINITY, f64::min);
            assert!(
                near <= 3.0 * step,
                "a flow strays {:.1} coarse cells",
                near / step
            );
        }
    }

    #[test]
    fn a_flood_spreads_as_lobes_that_cover_its_low_and_stay_in_its_pass() {
        // F94: on flat ground the lobes fan out from the vent until they cover the low
        // FLOOD_LOBE_COVER times over, inside the cap the pass may write.
        let coarse = CubeSphere::new(32);
        let grid = CubeSphere::new(128);
        let body_radius = 1_000e3;
        let flat = vec![0.0f32; grid.cell_count()];
        let ground = PaddedField::new(&grid, &flat);
        let center = normalize([0.3, 1.0, 0.2]);
        let mut flood = disc_flood(&coarse, center, 0.12, 1.0, 7);
        flood.area = flood
            .cells()
            .keys()
            .map(|&i| coarse.cell_area(coarse.cell(i as usize)) * body_radius * body_radius)
            .sum();
        let lobes = flood_lobes(&flood, &coarse, &grid, body_radius, &ground, 3, 0.0);
        assert!(lobes.len() > 3, "{} lobes", lobes.len());
        let area: f64 = lobes.iter().map(lobe_area).sum();
        assert!(
            area >= FLOOD_LOBE_COVER * flood.area || lobes.len() == FLOOD_LOBES_MAX,
            "lobes cover {:.3e} m² of a {:.3e} m² low",
            area,
            flood.area
        );
        let (middle, cap) = flood.cap(&coarse);
        let limit = cap + 2.0 * coarse.delta();
        for l in &lobes {
            for p in &l.points {
                assert!(
                    angle_between(middle, *p) <= limit + 1e-9,
                    "a lobe leaves the pass"
                );
            }
        }
    }
}
