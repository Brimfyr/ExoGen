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

use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashMap};

use crate::constants::*;
use crate::crater::{CellWrite, StampOutput};
use crate::fields::{Field, LayerSet};
use crate::geom::normalize;
use crate::grid::CubeSphere;
use crate::grid::projection::{SparseField, sample_bilinear_sparse};
use crate::noise::GradientNoise;
use crate::params::ParamGraph;
use crate::rng::Rng;
use crate::thermal::ThermalHistory;
use crate::units::*;

pub const VERSION: u64 = 2;

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
        let mut sum = [0.0; 3];
        for &idx in cells.keys() {
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
        let radius = cells
            .keys()
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
pub struct FloodProcess {
    grid: CubeSphere,
    body_radius: f64,
    thermal: ThermalHistory,
    rng: Rng,
    /// Cell areas, m².
    areas: Vec<f64>,
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
    pub records: Vec<FloodRecord>,
}

const NEIGHBOURS: [(i32, i32); 4] = [(1, 0), (-1, 0), (0, 1), (0, -1)];

impl FloodProcess {
    /// None for regimes whose volcanism is handled elsewhere.
    pub fn new(graph: &ParamGraph, grid: CubeSphere, thermal: &ThermalHistory) -> Option<Self> {
        if graph.str("interior.tectonic_regime") != "stagnant_lid" {
            return None;
        }
        let body_radius = graph.f64("body.radius");
        let areas = (0..grid.cell_count())
            .map(|k| grid.cell_area(grid.cell(k)) * body_radius * body_radius)
            .collect();
        Some(Self {
            grid,
            body_radius,
            thermal: thermal.clone(),
            rng: Rng::new(graph.f64("seed.volcanism") as u64 ^ 0x7F10_0D5A),
            areas,
            thickness: vec![0.0; grid.cell_count()],
            last_time: vec![f32::NEG_INFINITY; grid.cell_count()],
            visited: vec![0; grid.cell_count()],
            generation: 0,
            heap: BinaryHeap::new(),
            carry: 0.0,
            records: Vec::new(),
        })
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

    /// Extrusive volume between two times, m³, by midpoint rule over sub-intervals.
    pub fn budget(&self, t0: f64, t1: f64) -> f64 {
        let n = 8;
        let dt = (t1 - t0) / n as f64;
        let mut v = 0.0;
        for i in 0..n {
            let t = t0 + (i as f64 + 0.5) * dt;
            v += self.thermal.at(t).magma_rate * self.eruptive_fraction(t) * dt * GYR;
        }
        v
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
        step_index: u64,
    ) -> Vec<FloodEvent> {
        let mut budget = self.budget(t0, t1) + self.carry;
        if budget < VOLCANISM_EVENT_VOLUME.0 {
            self.carry = budget;
            return Vec::new();
        }
        let cumulative = self.vent_weights(layers, boost);
        let total = *cumulative.last().unwrap_or(&0.0);
        if total <= 0.0 {
            return Vec::new();
        }
        let (lo, hi) = VOLCANISM_EVENT_VOLUME;
        let mut j = 0u64;
        let mut out = Vec::new();
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
        out.sort_by(|a, b| a.time_gyr.total_cmp(&b.time_gyr).then(a.id.cmp(&b.id)));
        out
    }

    /// Erupt one event onto the surface as it is now, and record it.
    pub fn apply(&mut self, layers: &mut LayerSet, ev: FloodEvent) {
        let rec = self.erupt(layers, ev.vent, ev.volume, ev.time_gyr, ev.id);
        if !rec.ponds.is_empty() {
            self.records.push(rec);
        }
    }

    /// Sample and emplace a step's eruptions.
    pub fn step(&mut self, layers: &mut LayerSet, t0: f64, t1: f64, step_index: u64) {
        for ev in self.sample(layers, None, t0, t1, step_index) {
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
            let deep = layers.deep.as_mut_slice();
            let melt = layers.melt.as_mut_slice();
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
                deep[idx] = 0.0;
                melt[idx] = 0.0;
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
            id,
        };
        rec.area = rec.cells().keys().map(|&i| self.areas[i as usize]).sum();
        rec
    }
}

/// A coarse cell the fine fill may reach: the level, and the previous surface when the
/// cell is inside the coarse footprint.
type Reach = (f32, Option<(u16, f32)>);

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
    pub fn new(coarse: CubeSphere, grid: CubeSphere, floods: &[FloodRecord]) -> Self {
        let mut thickness = Field::<f32>::new(coarse);
        {
            let t = thickness.as_mut_slice();
            for f in floods {
                for (&idx, c) in &f.cells() {
                    t[idx as usize] += c.thickness;
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
pub fn flood_compute(
    layers: &LayerSet,
    grid: &CubeSphere,
    coarse: &CubeSphere,
    rec: &FloodRecord,
    smear: &LavaSmear,
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
    let padded_ind = SparseField::new(*coarse, covered.keys().map(|&idx| (idx as usize, 1.0)));
    let lobes = GradientNoise::new(rec.id ^ 0x0F10_0D00);
    let lobe_frequency = 1.0 / (FLOOD_MARGIN_LOBE_CELLS * coarse.delta());
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
    let time = rec.time_gyr as f32;
    let mut cells = Vec::with_capacity(levels.len() * (ratio * ratio) as usize);
    let mut taken = Vec::with_capacity(levels.len() * (ratio * ratio) as usize);
    let mut coarse_cells: Vec<(&usize, &Reach)> = levels.iter().collect();
    coarse_cells.sort_by_key(|(k, _)| **k);
    for (&cidx, &(level, previous)) in coarse_cells {
        let c = coarse.cell(cidx);
        for a in 0..ratio {
            for b in 0..ratio {
                let fine = crate::grid::Cell {
                    face: c.face,
                    i: c.i * ratio + a,
                    j: c.j * ratio + b,
                };
                let idx = grid.index(fine);
                let d = grid.center(fine);
                let own = sample_bilinear_sparse(coarse, &padded, d);
                taken.push((idx as u32, own));
                // A younger crater sits here: it was stamped on top of this lava in
                // the loop, so the base is right as it is.
                if smear.crater_time[idx] > time + 1e-6 {
                    stats[3] += 1;
                    continue;
                }
                // The ground as the flood found it, and the smears of floods still to
                // come, which ride on top of whatever is written.
                let remaining = smear.total[idx] - smear.done[idx];
                let ground = bedrock[idx] - remaining;
                let carry = remaining - own;
                let mut outline = sample_bilinear_sparse(coarse, &padded_ind, d) as f64;
                if (outline - FLOOD_MARGIN_THRESHOLD).abs() < FLOOD_MARGIN_WOBBLE {
                    outline += FLOOD_MARGIN_WOBBLE * lobes.fbm(d, lobe_frequency, 2, 2.0, 0.5);
                }
                // A flow drapes what stands less than its own thickness above its level.
                let drape = FLOOD_SHEET_THICKNESS as f32;
                if ground < level + drape && outline >= FLOOD_MARGIN_THRESHOLD {
                    stats[0] += 1;
                    let fill = (level - ground).max(0.0);
                    let surface = ground + (1.0 - LAVA_SUBSIDENCE as f32) * fill;
                    cells.push(CellWrite {
                        idx: idx as u32,
                        bedrock: surface + carry,
                        surface: true,
                        age: time,
                        unit: UNIT_LAVA,
                        deep: 0.0,
                        melt: 0.0,
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
                    });
                }
            }
        }
    }
    (StampOutput { cells }, taken, stats)
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
        let mut p = FloodProcess::new(&g, grid, &thermal).unwrap();
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
            let p = FloodProcess::new(&g, grid, &thermal).unwrap();
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
}
