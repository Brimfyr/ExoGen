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
use crate::thermal::ThermalHistory;
use crate::units::*;

pub const VERSION: u64 = 14;

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

/// A coarse cell a flood covered: the thickness added and what the surface was before.
#[derive(Clone, Copy, Debug)]
pub struct FloodCell {
    pub idx: u32,
    pub thickness: f32,
    pub prev_unit: u16,
    pub prev_age: f32,
}

/// One pond of a flood: the level it filled to and the coarse cells it covered.
#[derive(Clone, Debug)]
pub struct Pond {
    pub level: f32,
    pub cells: Vec<FloodCell>,
}

/// An eruption drawn for a step, not yet emplaced.
#[derive(Clone, Copy, Debug)]
pub struct FloodEvent {
    pub vent: usize,
    pub volume: f64,
    pub time_gyr: f64,
    pub id: u64,
}

/// What a flood did to one coarse cell over all its ponds.
#[derive(Clone, Copy, Debug)]
pub struct Covered {
    pub thickness: f32,
    pub level: f32,
    pub prev_unit: u16,
    pub prev_age: f32,
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
}

impl FloodRecord {
    /// Coarse cells covered, with their total thickness, final level and the surface
    /// the flood found.
    pub fn cells(&self) -> HashMap<u32, Covered> {
        let mut out: HashMap<u32, Covered> = HashMap::new();
        for pond in &self.ponds {
            for c in &pond.cells {
                let e = out.entry(c.idx).or_insert(Covered {
                    thickness: 0.0,
                    level: pond.level,
                    prev_unit: c.prev_unit,
                    prev_age: c.prev_age,
                });
                e.thickness += c.thickness;
                e.level = pond.level;
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
        Some(Self {
            grid,
            body_radius,
            thermal: thermal.clone(),
            rng,
            areas,
            provinces,
            province_weight,
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
            age: graph.f64("body.age"),
            records: Vec::new(),
        })
    }

    /// Take the resurfacing event's fill into the thickness and last-flow records.
    pub fn note_resurfacing(&mut self, fill: &[f32], time_gyr: f64) {
        for (k, &f) in fill.iter().enumerate() {
            if f > 0.0 {
                self.thickness[k] += f;
                self.last_time[k] = time_gyr as f32;
            }
        }
    }

    /// Fraction of the magma that reaches the surface at a time: the extrusive share,
    /// cut back as the lithosphere thickens and as contraction closes the pathways.
    pub fn eruptive_fraction(&self, t: f64) -> f64 {
        let s = self.thermal.at(t);
        let strain = 2.0 * s.contraction / self.body_radius;
        let ascent = 1.0
            / (1.0
                + ((s.lithosphere - VOLCANISM_ASCENT_LITHOSPHERE) / VOLCANISM_ASCENT_WIDTH).exp());
        VOLCANISM_EXTRUSIVE_FRACTION * ascent * (-(strain / VOLCANISM_SHUTOFF_STRAIN).powi(2)).exp()
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

    /// Vent weights: thin crust and low ground, and any rise's boost.
    fn vent_weights(&self, layers: &LayerSet, boost: Option<&[f64]>) -> Vec<f64> {
        let crust = layers.crust_thickness.as_slice();
        let bedrock = layers.bedrock.as_slice();
        let n = bedrock.len();
        let mut heights: Vec<f32> = bedrock.to_vec();
        heights.sort_by(f32::total_cmp);
        let median = heights[n / 2] as f64;
        let crust_min = crust.iter().cloned().fold(f32::INFINITY, f32::min) as f64;
        let mut acc = 0.0;
        (0..n)
            .map(|k| {
                let w = (-(crust[k] as f64 - crust_min) / VOLCANISM_VENT_CRUST_SCALE).exp()
                    * (-((bedrock[k] as f64 - median) / VOLCANISM_VENT_HEIGHT_SCALE)
                        .clamp(-4.0, 4.0))
                    .exp()
                    * self.areas[k]
                    * self.province_weight[k]
                    * boost.map_or(1.0, |b| b[k]);
                acc += w;
                acc
            })
            .collect()
    }

    /// Eruptions between two times, in order of formation, not yet emplaced. Vents are
    /// drawn against the surface as it is at the start of the step.
    pub fn sample(
        &mut self,
        layers: &LayerSet,
        boost: Option<&[f64]>,
        t0: f64,
        t1: f64,
    ) -> Vec<FloodEvent> {
        let bins = crate::stages::time_bins(self.age, t0, t1);
        let step_budget: f64 = bins.iter().map(|&(_, b0, b1)| self.budget(b0, b1)).sum();
        if step_budget + self.carry < VOLCANISM_EVENT_VOLUME.0 {
            self.carry += step_budget;
            return Vec::new();
        }
        let cumulative = self.vent_weights(layers, boost);
        let total = *cumulative.last().unwrap_or(&0.0);
        if total <= 0.0 {
            return Vec::new();
        }
        let mut out = Vec::new();
        for (bin, b0, b1) in bins {
            self.sample_bin(&cumulative, total, bin, b0, b1, &mut out);
        }
        out.sort_by(|a, b| a.time_gyr.total_cmp(&b.time_gyr).then(a.id.cmp(&b.id)));
        out
    }

    /// The eruptions of one time bin, drawn with the bin as the random index against
    /// the vent weights of the step.
    fn sample_bin(
        &mut self,
        cumulative: &[f64],
        total: f64,
        bin: u64,
        t0: f64,
        t1: f64,
        out: &mut Vec<FloodEvent>,
    ) {
        let step_index = bin;
        let mut budget = self.budget(t0, t1) + self.carry;
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
            out.push(FloodEvent {
                vent,
                volume,
                time_gyr,
                id: (step_index << 32) | j,
            });
            budget -= volume;
            j += 1;
        }
        self.carry = budget.max(0.0);
    }

    /// Erupt one event onto the surface as it is now, and record it.
    pub fn apply(&mut self, layers: &mut LayerSet, ev: FloodEvent) {
        let rec = self.erupt(layers, ev.vent, ev.volume, ev.time_gyr, ev.id);
        if !rec.ponds.is_empty() {
            self.records.push(rec);
        }
    }

    /// Sample and emplace a step's eruptions.
    pub fn step(&mut self, layers: &mut LayerSet, t0: f64, t1: f64, _step_index: u64) {
        for ev in self.sample(layers, None, t0, t1) {
            self.apply(layers, ev);
        }
    }

    /// Steepest descent from a cell to a local minimum.
    fn descend(&self, bedrock: &[f32], mut idx: usize) -> usize {
        for _ in 0..self.grid.cell_count() {
            let cell = self.grid.cell(idx);
            let mut best = (bedrock[idx], idx);
            for (di, dj) in NEIGHBOURS {
                let n = self.grid.index(self.grid.neighbor(cell, di, dj));
                if bedrock[n] < best.0 {
                    best = (bedrock[n], n);
                }
            }
            if best.1 == idx {
                return idx;
            }
            idx = best.1;
        }
        idx
    }

    /// Flood the terrain from a vent with a volume. The lava runs down to the nearest
    /// low and fills one flat lake there: cells join in order of height from the seed,
    /// over passes into deeper ground when the level reaches them, until the volume
    /// cannot lift the lake to the next cell. A basin fed many times fills as one lake.
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
        let bedrock = layers.bedrock.as_slice();
        let seed = self.descend(bedrock, vent);
        self.generation = self.generation.wrapping_add(1);
        if self.generation == 0 {
            self.visited.iter_mut().for_each(|v| *v = 0);
            self.generation = 1;
        }
        let generation = self.generation;
        let seed_dir = grid.center(grid.cell(seed));
        self.heap.clear();
        self.visited[seed] = generation;
        self.heap.push(Reverse(Key(bedrock[seed], 0.0, seed)));
        let mut popped: Vec<(usize, f32)> = Vec::new();
        let (mut area_sum, mut vol_sum) = (0.0f64, 0.0f64);
        while let Some(Reverse(Key(h, _, idx))) = self.heap.pop() {
            // The lowest cell on the frontier: if the volume cannot lift everything
            // taken so far to its height, or would spread thinner than a flow does,
            // the lake stops here.
            let a = self.areas[idx];
            if h as f64 * (area_sum + a) - (vol_sum + a * h as f64) > volume
                || (area_sum + a) * FLOOD_SHEET_THICKNESS > volume
            {
                break;
            }
            area_sum += a;
            vol_sum += a * h as f64;
            popped.push((idx, h));
            if popped.len() >= FLOOD_MAX_CELLS {
                break;
            }
            let cell = grid.cell(idx);
            for (di, dj) in NEIGHBOURS {
                let n = grid.index(grid.neighbor(cell, di, dj));
                if self.visited[n] != generation {
                    self.visited[n] = generation;
                    let centre = grid.center(grid.cell(n));
                    // Flows run only so far from their vent.
                    if crate::geom::angle_between(vent_dir, centre) * self.body_radius
                        > FLOOD_MAX_REACH
                    {
                        continue;
                    }
                    let dist = crate::geom::angle_between(seed_dir, centre);
                    self.heap.push(Reverse(Key(bedrock[n], dist as f32, n)));
                }
            }
        }
        // Ground beyond a pass joins below the pass, so order by height before finding
        // the level: the largest k whose first k cells the volume lifts to the k-th
        // height, then the level that spends the rest.
        popped.sort_by(|a, b| a.1.total_cmp(&b.1).then(a.0.cmp(&b.0)));
        let (mut area_sum, mut vol_sum) = (0.0f64, 0.0f64);
        let mut level = popped.first().map_or(0.0, |p| p.1 as f64);
        for k in 0..popped.len() {
            let h = popped[k].1 as f64;
            let a = self.areas[popped[k].0];
            if h * (area_sum + a) - (vol_sum + a * h) > volume {
                break;
            }
            area_sum += a;
            vol_sum += a * h;
            level = (volume + vol_sum) / area_sum;
            if k + 1 < popped.len() {
                level = level.min(popped[k + 1].1 as f64);
            }
        }
        let mut recorded = Vec::new();
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
            for &(idx, h) in &popped {
                let t = level - h as f64;
                if t <= FLOOD_MIN_THICKNESS {
                    continue;
                }
                // The pile subsides by a share of its thickness: what is recorded and
                // added is the height left after that.
                let added = (1.0 - LAVA_SUBSIDENCE) * t;
                placed += t * self.areas[idx];
                recorded.push(FloodCell {
                    idx: idx as u32,
                    thickness: added as f32,
                    prev_unit: units[idx],
                    prev_age: ages[idx],
                });
                self.thickness[idx] += t as f32;
                self.last_time[idx] = time_gyr as f32;
                bedrock[idx] = (h as f64 + added) as f32;
                ages[idx] = time_gyr as f32;
                units[idx] = UNIT_LAVA;
                melt[idx] = 0.0;
                cover[idx] = (t / LAVA_COVER_DEPTH).min(1.0) as f32;
                composition[idx] = drift;
                lava_age[idx] = time_gyr as f32;
            }
        }
        let ponds = if recorded.is_empty() {
            Vec::new()
        } else {
            vec![Pond {
                level: level as f32,
                cells: recorded,
            }]
        };
        let mut rec = FloodRecord {
            vent: vent_dir,
            volume: placed,
            time_gyr,
            area: 0.0,
            ponds,
            composition: flood_composition(time_gyr, vent_dir, self.rng.normal(id, 0xBA7C)),
            id,
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
type Reach = (f32, Option<(u16, f32)>);

/// A fine cell a flood's pass visits: its index, the level the lava could reach there and,
/// inside the coarse footprint, the surface the flood found.
type Visit = (usize, f32, Option<(u16, f32)>);

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
    /// The time of the latest loop crater whose continuous ejecta covers each cell at
    /// the finishing resolution, or none: under a crater younger than a flood the base is
    /// the crater stamped on the lava, and is right as it is.
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
        let lobes = flood_lobes(rec, coarse, grid, f.body_radius, f.ground, f.seed);
        lobe_cover(&lobes, grid, f.body_radius, rec.id, stamper)
    });
    let hollow_steps = flows.map_or(0, |f| {
        (f.hollow_reach / (grid.delta() * f.body_radius)).ceil() as usize
    });
    // Where the loop filled the flood deeper than its sheets, the lava pond fills the low
    // whatever its lobes do: judged by the flood's own coarse thickness, sampled smoothly,
    // since the ground pieced together from the coarse smears is not good to a sheet.
    let ponded_thickness =
        ((1.0 - LAVA_SUBSIDENCE) * FLOOD_PONDED_SHEETS * FLOOD_SHEET_THICKNESS) as f32;
    // Coarse cells to visit: the footprint and its ring, with the level the lava could
    // reach there and, inside the footprint, the surface the flood found.
    let mut levels: HashMap<usize, Reach> = HashMap::new();
    for (&idx, c) in &covered {
        levels.insert(idx as usize, (c.level, Some((c.prev_unit, c.prev_age))));
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
    let time = rec.time_gyr as f32;
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
        previous: Option<(u16, f32)>,
        own: f32,
        ground: f32,
        carry: f32,
        younger: bool,
        near: bool,
        hollow: bool,
        reached: bool,
        ponded: bool,
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
            };
            if w.younger {
                return w;
            }
            // The ground as the flood found it: the base less the smears not yet
            // taken out and less what the slow processes, the domes, shields,
            // canyon and coronae, have added since; both ride on top of whatever
            // is written.
            let remaining = smear.total[idx] - smear.done[idx];
            let later = deform.since_flood(d, rec.time_gyr) as f32;
            let ground = bedrock[idx] - remaining - later;
            // A hollow finer than the loop's cells, a crater the loop never saw, lies
            // deeper below the level than the ground the loop flooded there (F95). Both
            // grounds carry the same ramps between the coarse cells, so only the finer
            // relief tells; in the ring round the low, whose coarse ground stood above the
            // level, a ramp only raises the ground, so a hollow there is the finer relief
            // too.
            w.hollow = w.near
                && ground
                    < level - own / (1.0 - LAVA_SUBSIDENCE as f32) - FLOOD_SHEET_THICKNESS as f32;
            // Inside the flood's own footprint, where the base still says this lava
            // and no finer pass has written, the loop's surface was the flood's own
            // sheet at its level: what stands above the level there is the
            // resampling's ramp toward a higher neighbour, not ground the flood
            // found, and the sheet lies level over it, at least its own thickness
            // deep. Elsewhere what stands above the level is a rim or an older
            // floor, and the flood ends against it. The lobed outline bounds the
            // fill either way; what it leaves inside the coarse footprint keeps the
            // loop's surface, the older lava's own.
            let own_floor = previous.is_some() && units[idx] == UNIT_LAVA && !smear.touched[idx];
            w.ground = if own_floor {
                ground.min(level - FLOOD_SHEET_THICKNESS as f32)
            } else {
                ground
            };
            w.carry = remaining - own + later;
            w.reached = match &lobed {
                Some(lobed) => {
                    own > ponded_thickness || lobed.get(&idx).copied().unwrap_or(0.0) >= 0.5
                }
                None => {
                    let mut outline = sample_bilinear_sparse(coarse, &padded_ind, d) as f64;
                    if outline > 0.0 && outline < 1.0 {
                        outline += FLOOD_MARGIN_WOBBLE * lobes.fbm(d, lobe_frequency, 3, 2.0, 0.5);
                    }
                    outline >= FLOOD_MARGIN_THRESHOLD
                }
            };
            w
        })
        .collect();
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
    if lobed.is_some() {
        let drape = FLOOD_SHEET_THICKNESS as f32;
        // Lava lying in the low or its ring: the pond and the flows next to it, not a
        // lobe's sheet far beyond.
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
                            continue;
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
                    if ground
                        >= level
                            - own / (1.0 - LAVA_SUBSIDENCE as f32)
                            - FLOOD_SHEET_THICKNESS as f32
                    {
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
                    });
                }
            }
        }
    }
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
            ..
        } = w;
        if younger {
            stats[3] += 1;
            // Under a crater younger than the flood the base is the crater stamped on the
            // lava, right as it is in relief; what the crater dug into is the lava where
            // the flood's flows reached (F94), each flood over the last in their order as
            // beside the crater, and where none did, the surface the floods found.
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
            } else if let Some((prev_unit, prev_age)) = previous
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
            }
            continue;
        }
        // Outside the low the loop filled, a lobe of a flood drawn as its flows lays a
        // sheet on the ground it crosses rather than filling to the low's level, unless
        // the lava ponds in a hollow there.
        let level = if lobed.is_some() && previous.is_none() && !ponded {
            level.min(ground + FLOOD_SHEET_THICKNESS as f32)
        } else {
            level
        };
        if watched == Some(idx) {
            eprintln!(
                "watch flood {} t {:.4}: level {:.1} own {:.1} ground {:.1} hollow {} ponded {} previous {:?} touched {} unit {} lobe {:?} reached {:?}",
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
            );
        }
        // A flow drapes what stands less than its own thickness above its level.
        let drape = FLOOD_SHEET_THICKNESS as f32;
        if ground < level + drape && reached {
            stats[0] += 1;
            let fill = (level - ground).max(0.0);
            // A high the flow drapes carries a veneer, the sheet's thickness less
            // what the high stands above the level: thin lava over it, but lava,
            // not the ground it hides (the ghost craters of the maria show in
            // relief, not in albedo).
            let sheet = fill.max(drape - (ground - level).max(0.0));
            let surface = ground + (1.0 - LAVA_SUBSIDENCE as f32) * fill;
            cells.push(CellWrite {
                idx: idx as u32,
                bedrock: surface + carry,
                surface: true,
                age: time,
                unit: UNIT_LAVA,
                deep: deep[idx],
                melt: 0.0,
                cover: (sheet / LAVA_COVER_DEPTH as f32).clamp(0.0, 1.0),
            });
        } else {
            if ground >= level + drape {
                stats[1] += 1;
            } else {
                stats[2] += 1;
            }
            // Not reached: take the smear off. Inside the coarse footprint the
            // base says lava where no pass has written yet; give it back the
            // surface the flood found.
            let (unit, age) = match previous {
                Some((u, a)) if !smear.touched[idx] && units[idx] == UNIT_LAVA => (u, a),
                _ => (units[idx], ages[idx]),
            };
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
    (
        StampOutput {
            cells,
            composition: Some(rec.composition),
            stain: None,
            steep: None,
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
        for pond in &floods[k].ponds {
            for c in &pond.cells {
                let stack = last.entry(c.idx).or_default();
                stack.push(k);
                if stack.len() > FLOOD_VISIBLE_DEPTH {
                    stack.remove(0);
                }
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
fn lobe_area(lobe: &Lobe) -> f64 {
    (lobe.points.len().max(1) - 1) as f64 * lobe.spacing * 2.0 * lobe.half_width
}

/// The flows of one flood at the finishing resolution (F94): from the coarse cell of its
/// low nearest its vent, where its lava came in, a trunk across the low toward its middle,
/// up to FLOOD_BRANCHES.1 branches, and breakouts on the lobes already laid, at their
/// fronts and along them, until the lobes cover the low FLOOD_LOBE_COVER times over. Each
/// lobe is a path over the ground (F59), a share of the low's reach long and a share of
/// its own length wide, wandering at its width's scale and ponding where it would climb
/// more than a sheet.
pub fn flood_lobes(
    rec: &FloodRecord,
    coarse: &CubeSphere,
    grid: &CubeSphere,
    body_radius: f64,
    ground: &PaddedField,
    seed: u64,
) -> Vec<Lobe> {
    let covered = rec.cells();
    let Some(entry_idx) = covered.keys().copied().min_by(|&a, &b| {
        let da = angle_between(rec.vent, coarse.center(coarse.cell(a as usize)));
        let db = angle_between(rec.vent, coarse.center(coarse.cell(b as usize)));
        da.total_cmp(&db).then(a.cmp(&b))
    }) else {
        return Vec::new();
    };
    let r = body_radius;
    let cell = grid.delta() * r;
    let entry = coarse.center(coarse.cell(entry_idx as usize));
    let (middle, cap) = rec.cap(coarse);
    // The lobes stay inside the flood's pass, whose cap reaches 2.5 coarse cells past its
    // footprint: a path stops where its lobe's margin would cross 2 of them.
    let limit = cap + 2.0 * coarse.delta();
    let outside = |q: [f64; 3], hw: f64| {
        angle_between(middle, q)
            > limit - (hw * (1.0 + IO_FLOW_LOBING) + (IO_FLOW_EDGE * hw).max(1.5 * cell)) / r
    };
    let reach = covered
        .keys()
        .map(|&i| angle_between(entry, coarse.center(coarse.cell(i as usize))))
        .fold(0.0, f64::max)
        * r
        + coarse.delta() * r;
    let rng = Rng::new(seed ^ rec.id ^ 0xF100_D10B);
    let rules = |half_width: f64| PathRules {
        min_spacing: 0.5 * cell,
        slope_follow: IO_FLOW_SLOPE_FOLLOW,
        slope_scale: IO_FLOW_SLOPE_SCALE,
        steer_length: FLOOD_STEER_LENGTH,
        wander_rad: FLOOD_WANDER_RAD,
        wander_length: FLOOD_WANDER_WIDTHS * half_width,
        fill_depth: FLOOD_SHEET_THICKNESS,
        width_wobble: IO_FLOW_WIDTH_WOBBLE,
    };
    let (w_lo, w_hi) = FLOOD_LOBE_WIDTH;
    let width = |index: u64, length: f64| {
        (length * (w_lo + (w_hi - w_lo) * rng.uniform(index, 0))).max(1.5 * cell)
    };
    // The trunk heads across the low toward its middle.
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
    let (t_lo, t_hi) = FLOOD_TRUNK_SHARE;
    let trunk_length = reach * (t_lo + (t_hi - t_lo) * rng.uniform(3, 0));
    let hw = width(4, trunk_length);
    let mut lobes = vec![
        trace_flow(
            grid,
            r,
            &rules(hw),
            rec.id,
            5,
            entry,
            toward,
            trunk_length,
            hw,
            0.0,
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
        let hw = width(20 + b, length);
        let heading = rotate_about(p, tangent, turn);
        let lobe = trace_flow(
            grid,
            r,
            &rules(hw),
            rec.id,
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
    // Breakouts on the lobes already laid, until they cover the low.
    let target = FLOOD_LOBE_COVER * rec.area;
    let (l_lo, l_hi) = FLOOD_LOBE_LENGTH;
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
            rec.id,
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
fn lobe_cover(
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
        for name in ["Moon", "Mercury", "Mars"] {
            let (g, _) = preset(name).unwrap().to_graph();
            let thermal = ThermalHistory::from_graph(&g);
            let grid = CubeSphere::new(8);
            let p = FloodProcess::new(&g, grid, &thermal, &[]).unwrap();
            println!(
                "{name}: t  lithosphere_km  magma_m3_s  strain  eruptive  budget_km3_per_0.1Gyr"
            );
            let mut t = 0.0;
            while t < g.f64("body.age") {
                let s = thermal.at(t);
                println!(
                    "  {t:.1}  {:.0}  {:.2}  {:.2e}  {:.3}  {:.3e}",
                    s.lithosphere / 1e3,
                    s.magma_rate,
                    2.0 * s.contraction / p.body_radius,
                    p.eruptive_fraction(t),
                    p.budget(t, t + 0.1) / 1e9
                );
                t += 0.5;
            }
        }
    }

    #[test]
    fn lava_runs_downhill_fills_the_pit_and_keeps_its_volume() {
        let grid = CubeSphere::new(32);
        let mut p = process(grid);
        let mut layers = LayerSet::new(grid);
        // A pit of one cell, 1 km deep, and a vent a few cells away on flat ground.
        let pit = grid.index(crate::grid::Cell {
            face: 0,
            i: 16,
            j: 16,
        });
        let vent = grid.index(crate::grid::Cell {
            face: 0,
            i: 20,
            j: 16,
        });
        layers.bedrock.as_mut_slice()[pit] = -1_000.0;
        // Slope towards the pit so the lava finds it.
        for k in 0..grid.cell_count() {
            let d =
                crate::geom::angle_between(grid.center(grid.cell(k)), grid.center(grid.cell(pit)));
            if k != pit {
                layers.bedrock.as_mut_slice()[k] = (d * 1e4) as f32;
            }
        }
        let area_pit = p.areas[pit];
        let volume = 0.5 * 1_000.0 * area_pit;
        let rec = p.erupt(&mut layers, vent, volume, 1.0, 1);
        let covered = rec.cells();
        assert!(covered.contains_key(&(pit as u32)), "pit not flooded");
        let filled: f64 = covered
            .iter()
            .map(|(&idx, c)| c.thickness as f64 * p.areas[idx as usize])
            .sum();
        // What stands is the fill less its subsidence; the record keeps the volume.
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
        // The pit took it all: half full less subsidence, flat, lava unit, new age.
        let b = layers.bedrock.as_slice()[pit];
        let expected = -1_000.0 + (1.0 - LAVA_SUBSIDENCE) * 500.0;
        assert!((b as f64 - expected).abs() < 5.0, "{b} vs {expected}");
        assert_eq!(layers.unit_id.as_slice()[pit], UNIT_LAVA);
        assert_eq!(layers.surface_age.as_slice()[pit], 1.0);
        assert_eq!(covered.len(), 1);
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
            ponds: vec![Pond { level: 50.0, cells }],
            composition: 0.0,
            id,
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
            let deform = Deformation::new(coarse, body_radius, 0.01, 4.5, &[], &[], &[], &[]);
            let mut stamper = Stamper::new(&grid);
            let flows = FloodFlows {
                ground: &ground,
                body_radius,
                seed: 3,
                hollow_reach: 40e3,
            };
            let (out, _, _) = flood_compute(
                &layers,
                &grid,
                &coarse,
                &flood,
                &smear,
                &deform,
                Some(flows),
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
        let lobes = flood_lobes(&flood, &coarse, &grid, body_radius, &ground, 3);
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
