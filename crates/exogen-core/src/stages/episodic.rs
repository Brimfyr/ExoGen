// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! The episodic lid (PLAN.md 4.3, Venus type).
//!
//! The lid's regions founder once in a while, each renewing its share of the surface (F296).
//! Each overturn floods its regions' plains: every depression fills most of the way to its
//! spill point, which takes the old craters with it and leaves smooth lowlands under the
//! regional slopes. The regions stand on the old crust's smoothed height as their strength
//! ranks them, the weakest on the lowest ground and the strongest on the plateaus of crust
//! the mantle's flow thickened over its downwellings (F184), compensated on their roots: the
//! strongest that outlast the overturns that renew the plains about them, the tesserae (F300),
//! deformed into a fabric of intersecting ridges.
//! Afterwards the surface ages again under a sparse, un-degraded crater population, with
//! flow fields, coronae over the mantle's diapirs (F145, the corona process) and
//! steep-sided pancake domes on the plains.

use rayon::prelude::*;
use std::cmp::Reverse;
use std::collections::BinaryHeap;

use crate::constants::*;
use crate::crater::{CellWrite, StampOutput, Stamper};
use crate::fields::LayerSet;
use crate::geom::{angle_between, local_basis, offset};
use crate::grid::CubeSphere;
use crate::grid::projection::{PaddedField, sample_bilinear};
use crate::noise::GradientNoise;
use crate::params::ParamGraph;
use crate::rng::Rng;
use crate::stages::deform::Deformation;
use crate::stages::plume::{Trough, TroughKind, radial_trace};
use crate::stages::tessera::TesseraFabric;
use crate::stages::volcanism::{LavaSmear, overturn_composition};
use crate::units::{UNIT_LAVA, UNIT_TESSERA};

pub const VERSION: u64 = 17;

/// The overturns' resurfacing, kept for the finishing pass (F296).
#[derive(Clone, Debug)]
pub struct Resurfacing {
    /// When the oldest standing region of the lid formed, Gyr since formation.
    pub time_gyr: f64,
    /// Thickness the overturns added per coarse cell, m.
    pub fill: Vec<f32>,
    /// The surface each cell's plains were last filled towards, per coarse cell: the terrain
    /// with its depressions filled.
    pub filled: Vec<f32>,
    /// Smoothed height, on which the lid's regions stand as their strength ranks them, the
    /// tesserae its high end.
    pub score: Vec<f32>,
    pub threshold: f32,
    /// The score at each region's upper edge, the weakest's first; the first of the
    /// tesserae, the regions' count where there are none (F300); and when each region last
    /// formed, Gyr since formation.
    pub edges: Vec<f32>,
    pub tesserae_from: usize,
    pub formed: Vec<f32>,
    /// The body's surface gravity, m/s².
    pub gravity: f64,
    /// The strain the lid the event keeps took from the mantle's flow, a set for each state
    /// of the flow (F187).
    pub strain: Vec<LidStrain>,
}

impl Resurfacing {
    /// The lid's region at a score.
    pub fn region_of(&self, score: f32) -> usize {
        self.edges
            .partition_point(|&e| e < score)
            .min(self.edges.len().saturating_sub(1))
    }

    /// When the ground at a direction last formed, Gyr since formation: its region's, by its
    /// score read on the coarse grid.
    pub fn formed_at(&self, coarse: &CubeSphere, score: &PaddedField, d: [f64; 3]) -> f32 {
        self.formed[self.region_of(sample_bilinear(coarse, score, d))]
    }

    /// When the plains at a direction last formed, Gyr since formation: within the
    /// tesserae, those at their edge, which its lava ponded from.
    pub fn plains_at(&self, coarse: &CubeSphere, score: &PaddedField, d: [f64; 3]) -> f32 {
        let region = self.region_of(sample_bilinear(coarse, score, d));
        self.formed[region.min(self.tesserae_from.saturating_sub(1))]
    }
}

/// The score at each of `regions` equal shares of the surface's upper edge, the lowest's
/// first and the highest's unbounded (F296).
fn region_edges(grid: &CubeSphere, score: &[f32], regions: usize) -> Vec<f32> {
    let mut cells: Vec<(f32, f64)> = score
        .iter()
        .enumerate()
        .map(|(k, &s)| (s, grid.cell_area(grid.cell(k))))
        .collect();
    cells.sort_by(|a, b| a.0.total_cmp(&b.0));
    let total: f64 = cells.iter().map(|c| c.1).sum();
    let mut edges = Vec::with_capacity(regions);
    let mut acc = 0.0;
    let mut i = 0;
    for r in 0..regions - 1 {
        let upto = (r + 1) as f64 / regions as f64 * total;
        while i < cells.len() - 1 && acc + cells[i].1 < upto {
            acc += cells[i].1;
            i += 1;
        }
        edges.push(cells[i].0);
    }
    edges.push(f32::INFINITY);
    edges
}

/// The strain one state of the mantle's flow put on the lid an overturn keeps (F187).
#[derive(Clone, Debug)]
pub struct LidStrain {
    /// The time it was taken at, Gyr, and the lid's elastic thickness then, m.
    pub time: f64,
    pub elastic: f64,
    /// Per coarse cell, the horizontal strain in the cell's east and north basis,
    /// [ee, en, nn], extension positive.
    pub strain: Vec<[f32; 3]>,
    /// Per coarse cell, the most the drag pressed and pulled the lid past its strength, Pa.
    pub pressed: Vec<f32>,
    pub pulled: Vec<f32>,
}

/// A pancake dome: flat-topped and steep-sided.
#[derive(Clone, Debug)]
pub struct Dome {
    pub center: [f64; 3],
    pub radius: f64,
    pub height: f64,
    pub time_gyr: f64,
    pub id: u64,
}

impl Dome {
    pub fn cap_radius(&self, body_radius: f64) -> f64 {
        1.3 * self.radius / body_radius
    }
}

#[derive(PartialEq)]
struct Key(f32, usize);
impl Eq for Key {}
impl PartialOrd for Key {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for Key {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.0.total_cmp(&other.0).then(self.1.cmp(&other.1))
    }
}

/// The terrain with every depression filled to its spill point (priority flood from the
/// lowest cell, Barnes et al. 2014).
pub fn fill_depressions(grid: &CubeSphere, heights: &[f32]) -> Vec<f32> {
    fill_depressions_to(grid, heights, |_| false)
}

/// The terrain with every depression filled to its spill point, a liquid on it spilling
/// onto the cells `outlet` names wherever they stand lower than its surface, as well as
/// towards the lowest cell: each outlet keeps its own height (F301).
pub fn fill_depressions_to(
    grid: &CubeSphere,
    heights: &[f32],
    outlet: impl Fn(usize) -> bool,
) -> Vec<f32> {
    let n = grid.cell_count();
    let mut filled = vec![f32::INFINITY; n];
    let mut visited = vec![false; n];
    let start = (0..n)
        .min_by(|&a, &b| heights[a].total_cmp(&heights[b]))
        .unwrap_or(0);
    let mut heap = BinaryHeap::new();
    for k in (0..n).filter(|&k| k == start || outlet(k)) {
        heap.push(Reverse(Key(heights[k], k)));
        visited[k] = true;
    }
    while let Some(Reverse(Key(level, idx))) = heap.pop() {
        filled[idx] = level;
        let cell = grid.cell(idx);
        for (di, dj) in [(1, 0), (-1, 0), (0, 1), (0, -1)] {
            let m = grid.index(grid.neighbor(cell, di, dj));
            if !visited[m] {
                visited[m] = true;
                heap.push(Reverse(Key(heights[m].max(level), m)));
            }
        }
    }
    filled
}

/// Neighbour-averaged smoothing, `passes` times.
pub fn smooth(grid: &CubeSphere, values: &[f32], passes: usize) -> Vec<f32> {
    let mut a = values.to_vec();
    let mut b = vec![0.0f32; a.len()];
    for _ in 0..passes {
        b.par_iter_mut().enumerate().for_each(|(idx, out)| {
            let cell = grid.cell(idx);
            let mut s = a[idx] * 4.0;
            for (di, dj) in [(1, 0), (-1, 0), (0, 1), (0, -1)] {
                s += a[grid.index(grid.neighbor(cell, di, dj))];
            }
            *out = s / 8.0;
        });
        std::mem::swap(&mut a, &mut b);
    }
    a
}

/// Pancake dome profile at a fraction of the radius: flat top, steep sides.
pub fn dome_profile(x: f64) -> f64 {
    if x >= 1.0 {
        0.0
    } else {
        (1.0 - x.powi(6)).powf(0.6)
    }
}

/// The episodic lid through the loop: its overturns, then the coronae and domes.
pub struct EpisodicProcess {
    grid: CubeSphere,
    body_radius: f64,
    /// When the oldest standing region of the lid formed, Gyr since formation; the
    /// overturns since, each with the regions it renewed, and the next to come (F296).
    window: f64,
    overturns: Vec<(f64, Vec<u16>)>,
    next: usize,
    /// The lid's regions, the first of the strongest, the tesserae, and when each formed
    /// last, Gyr since formation (F296).
    regions: usize,
    tesserae_from: usize,
    formed: Vec<f64>,
    /// The passes of the one-cell kernel that smooth the old crust's height over the
    /// flexural parameter of the elastic lid the overturn broke (audit W3).
    tessera_smoothing: usize,
    gravity: f64,
    pub resurfacing: Option<Resurfacing>,
    /// The coronae over the mantle's diapirs since the lid formed (F145).
    pub coronae: Vec<crate::stages::corona::Corona>,
    pub domes: Vec<Dome>,
    /// Radial fracture systems round some of the coronae, stamped as graben.
    pub novae: Vec<Trough>,
    stamper: Stamper,
}

impl EpisodicProcess {
    pub fn new(
        graph: &ParamGraph,
        grid: CubeSphere,
        rises: &[crate::stages::plume::Rise],
        limbs: &crate::stages::stress::Limbs,
    ) -> Option<Self> {
        if graph.str("interior.tectonic_regime") != "episodic_lid" {
            return None;
        }
        let body_radius = graph.f64("body.radius");
        let age = graph.f64("body.age");
        // The overturns since the oldest standing region of the lid formed, which laid every
        // region's surface (F296).
        let mantle = crate::mantle::evolve(&crate::mantle::MantleInputs::from_graph(graph));
        let window = mantle.oldest_standing()?;
        let overturns: Vec<(f64, Vec<u16>)> = mantle
            .overturns
            .iter()
            .zip(&mantle.renewed)
            .filter(|&(&o, _)| o >= window)
            .map(|(&o, r)| (o, r.clone()))
            .collect();
        let regions = mantle.regions;
        // The tesserae: its strongest regions that outlasted every overturn since they
        // formed, older than the plains about them (F300); none where the last overturn
        // renewed the whole lid.
        let tesserae_from = regions - mantle.tesserae();
        let formed = mantle.formed_by(age);
        let rng = Rng::new(graph.f64("seed.tectonics") as u64 ^ 0xE915_0D1C);
        // The coronae over the mantle's diapirs (F145).
        let coronae = crate::stages::corona::rock_coronae(graph, &mantle, rises, limbs);
        let n_domes = graph.f64("volcanism.dome_count").round() as u64;
        let mut domes = Vec::new();
        for k in 0..n_domes {
            let (lo, hi) = DOME_DIAMETER;
            let diameter = lo * (hi / lo).powf(rng.uniform(1_000 + k, 1));
            domes.push(Dome {
                center: rng.unit_direction(1_000 + k, 5),
                radius: 0.5 * diameter,
                height: DOME_HEIGHT * (0.6 + 0.8 * rng.uniform(1_000 + k, 3)),
                time_gyr: (window + (age - window) * rng.uniform(1_000 + k, 4)).min(age),
                id: 1_000 + k,
            });
        }
        // Novae: the young coronae, whose heads, still fed and rising, stretch their domes'
        // crests past what breaks the lid, send fractures radiating out beyond their rims
        // (F185), the first stage of a corona (Stofan et al. 2001); an older corona's dome
        // has subsided and its floor's lava flooded them.
        let cell = grid.delta() * body_radius;
        let mut novae = Vec::new();
        for c in &coronae {
            let rising = c.head.as_ref().is_some_and(|h| {
                h.lift_tension > CORONA_FAULT_STRENGTH && age <= c.time_gyr + h.feeding
            });
            if !rising {
                continue;
            }
            let (east, north) = local_basis(c.center);
            let (n_lo, n_hi) = NOVA_FRACTURES;
            let n = n_lo + ((n_hi - n_lo) as f64 * rng.uniform(c.id, 8)) as u64;
            for j in 0..n {
                let draw = 50_000 + c.id * 16 + j;
                let phi = 2.0 * std::f64::consts::PI * (j as f64 + rng.uniform(draw, 1)) / n as f64;
                let r0 = c.outline_radius() * (0.3 + 0.6 * rng.uniform(draw, 2));
                let (l_lo, l_hi) = NOVA_LENGTH;
                let length = c.outline_radius() * (l_lo + (l_hi - l_lo) * rng.uniform(draw, 3));
                let start = offset(c.center, east, north, r0 / body_radius, phi);
                let spacing = 0.5 * cell;
                let points =
                    radial_trace(start, c.center, length, spacing, body_radius, &rng, draw);
                let (w_lo, w_hi) = NOVA_WIDTH;
                novae.push(Trough {
                    points,
                    spacing,
                    length,
                    width: w_lo * (w_hi / w_lo).powf(rng.uniform(draw, 4)),
                    depth: NOVA_DEPTH * (0.6 + 0.8 * rng.uniform(draw, 5)),
                    flank: 0.0,
                    tip_rate: 2.0 * SEGMENT_THROW_PER_LENGTH,
                    tilt: RIFT_TILT,
                    start_gyr: c.time_gyr,
                    end_gyr: c.time_gyr,
                    id: 2_000_000 + draw,
                    kind: TroughKind::Fracture,
                    segments: Vec::new(),
                });
            }
        }
        // The relief the lid held up on its strength, shorter than its flexural parameter,
        // is not the plateaus of crust that stand on their roots (audit W3): the elastic
        // lid the overturn broke, down to ELASTIC_ISOTHERM_ROCK, and the passes whose
        // spread, half the square root of their number in cells, reaches its parameter.
        let before = (window - 1e-3).max(0.0);
        let ts = mantle.surface_temperature;
        let elastic = mantle.lid_at(before)
            * ((ELASTIC_ISOTHERM_ROCK - ts) / (mantle.temperature_at(before) - ts).max(1.0))
                .clamp(0.0, 1.0);
        let parameter = crate::stages::load::flexural_parameter(elastic, graph.f64("body.gravity"));
        let tessera_smoothing = (4.0 * (parameter / cell).powi(2)).round().max(1.0) as usize;
        Some(Self {
            grid,
            body_radius,
            window,
            overturns,
            next: 0,
            regions,
            tesserae_from,
            formed,
            tessera_smoothing,
            gravity: graph.f64("body.gravity"),
            resurfacing: None,
            coronae,
            domes,
            novae,
            stamper: Stamper::new(&grid),
        })
    }

    /// When the oldest standing region of the lid formed, Gyr since formation (F296).
    pub fn window(&self) -> f64 {
        self.window
    }

    /// The time over which the mantle's flow strained the lid the tesserae keep (F187,
    /// F300): from when they formed, all at once, to when the oldest plains standing about
    /// them were laid; none without tesserae.
    pub fn strain_window(&self) -> Option<(f64, f64)> {
        let (plains, tesserae) = self.formed.split_at(self.tesserae_from);
        let from = tesserae.iter().copied().reduce(f64::max)?;
        let to = plains.iter().copied().reduce(f64::min)?;
        (to > from + 1e-9).then_some((from, to))
    }

    /// The coronae and domes on the ground as it stands, each on ground its region's last
    /// overturn laid before it (F296), and what the overturns did.
    pub fn standing(
        self,
    ) -> (
        Option<Resurfacing>,
        Vec<crate::stages::corona::Corona>,
        Vec<Dome>,
        Vec<Trough>,
    ) {
        let Self {
            grid,
            resurfacing,
            mut coronae,
            mut domes,
            novae,
            ..
        } = self;
        if let Some(r) = &resurfacing {
            let score = PaddedField::new(&grid, &r.score);
            coronae.retain(|c| c.time_gyr as f32 >= r.formed_at(&grid, &score, c.center));
            domes.retain(|d| d.time_gyr as f32 >= r.formed_at(&grid, &score, d.center));
        }
        (resurfacing, coronae, domes, novae)
    }

    /// Apply what falls in the step: the overturns, and the coarse part of the coronae.
    /// Returns each overturn's fill and time as it happens, for the lava ledger. The lid's
    /// regions are ranked, at the first overturn since its oldest standing region formed, on
    /// the height less what the plume rises' domes and the mantle's flow hold up (`held`),
    /// which is young dynamic topography, not the crust: the plateaus of crust the flow
    /// thickened stand highest (F184), the strongest regions on them, and the tesserae among
    /// them keep their ground as it is (F296).
    pub fn step(
        &mut self,
        layers: &mut LayerSet,
        t0: f64,
        t1: f64,
        held: Option<&[f32]>,
    ) -> Vec<(Vec<f32>, f64)> {
        let grid = self.grid;
        let body_radius = self.body_radius;
        let mut fills = Vec::new();
        while self.next < self.overturns.len() && self.overturns[self.next].0 <= t1 {
            let (time, renewing) = self.overturns[self.next].clone();
            self.next += 1;
            let heights = layers.bedrock.as_slice().to_vec();
            let n = grid.cell_count();
            if self.resurfacing.is_none() {
                let old_crust: Vec<f32> = match held {
                    Some(l) => heights.iter().zip(l).map(|(h, l)| h - l).collect(),
                    None => heights.clone(),
                };
                let score = smooth(&grid, &old_crust, self.tessera_smoothing);
                let edges = region_edges(&grid, &score, self.regions);
                self.resurfacing = Some(Resurfacing {
                    time_gyr: time,
                    fill: vec![0.0; n],
                    filled: heights.clone(),
                    threshold: edges[self.tesserae_from - 1],
                    score,
                    edges,
                    tesserae_from: self.tesserae_from,
                    formed: vec![0.0; self.regions],
                    gravity: self.gravity,
                    strain: Vec::new(),
                });
            }
            let rec = self.resurfacing.as_mut().expect("the lid's regions");
            let mut renewed = vec![false; self.regions];
            for &r in &renewing {
                renewed[r as usize] = true;
                rec.formed[r as usize] = time as f32;
            }
            // The lava fills the depressions on the plains it lays, spilling onto the ground
            // beyond them wherever that stands lower than its surface (F301): it never stands
            // above the old lid about it, its edge onlapping that lid's ground.
            let lays_none = |k: usize| {
                let region = rec.region_of(rec.score[k]);
                !renewed[region] || region >= self.tesserae_from
            };
            let filled = fill_depressions_to(&grid, &heights, lays_none);
            let mut fill = vec![0.0f32; n];
            {
                let bedrock = layers.bedrock.as_mut_slice();
                let ages = layers.surface_age.as_mut_slice();
                let units = layers.unit_id.as_mut_slice();
                let deep = layers.deep.as_mut_slice();
                let melt = layers.melt.as_mut_slice();
                let cover = layers.cover.as_mut_slice();
                let composition = layers.composition.as_mut_slice();
                let lava_age = layers.lava_age.as_mut_slice();
                let lava = overturn_composition(time);
                for k in 0..n {
                    let region = rec.region_of(rec.score[k]);
                    if !renewed[region] {
                        continue;
                    }
                    if region >= self.tesserae_from {
                        units[k] = UNIT_TESSERA;
                        ages[k] = time as f32;
                        cover[k] = 0.0;
                        continue;
                    }
                    let t = (EPISODIC_FILL_FRACTION as f32 * (filled[k] - heights[k])).max(0.0)
                        + EPISODIC_VENEER as f32;
                    fill[k] = t;
                    rec.fill[k] += t;
                    rec.filled[k] = filled[k];
                    bedrock[k] += t;
                    ages[k] = time as f32;
                    units[k] = UNIT_LAVA;
                    deep[k] = 0.0;
                    melt[k] = 0.0;
                    cover[k] = 1.0;
                    // Its lava rose everywhere it renewed, each place's from the mantle under
                    // it.
                    composition[k] = lava;
                    lava_age[k] = time as f32;
                }
            }
            fills.push((fill, time));
        }
        // Coronae (F145): the lava their cores flood renews the coarse ground at their
        // time; their profile and fabric are the finish's, on top of what the loop left.
        let cell = grid.delta() * body_radius;
        for c in &self.coronae {
            if c.time_gyr <= t0 || c.time_gyr > t1 {
                continue;
            }
            let covers = c.covers(&grid, body_radius, cell, &mut self.stamper);
            let ages = layers.surface_age.as_mut_slice();
            let units = layers.unit_id.as_mut_slice();
            let cover = layers.cover.as_mut_slice();
            for (idx, w, core) in covers {
                if w * core >= 0.5 && units[idx] != UNIT_TESSERA {
                    ages[idx] = c.time_gyr as f32;
                    units[idx] = UNIT_LAVA;
                    cover[idx] = 1.0;
                }
            }
        }
        fills
    }
}

/// Finishing pass of the event at the resolved resolution (F186): the old ground with the
/// tessera's fabric from the strain its lid took (F187), over the tesserae and sinking
/// beneath the plains beyond them, and the
/// event's lava ponded against it at the tesserae's edge: the plains onlap the tesserae,
/// their lava running up the fabric's valleys and leaving its ridges standing as inliers,
/// each cell tessera or plains by the share of it the lava left above its surface. Returns
/// the writes and the smear taken out per cell.
#[allow(clippy::too_many_arguments)]
pub fn resurfacing_compute(
    layers: &LayerSet,
    grid: &CubeSphere,
    coarse: &CubeSphere,
    body_radius: f64,
    rec: &Resurfacing,
    smear: &LavaSmear,
    seed: u64,
    deform: &Deformation,
) -> (StampOutput, Vec<(u32, f32)>) {
    let padded_fill = PaddedField::new(coarse, &rec.fill);
    let padded_filled = PaddedField::new(coarse, &rec.filled);
    let padded_score = PaddedField::new(coarse, &rec.score);
    // Each plains' fill level read from the coarse cells its own overturn laid alone (F302):
    // a young plain's spill point never blended with the older plains' about it, nor theirs
    // with its, which drew a lip inside the young plains' edge and a moat outside it. Each
    // overturn's level and its cells' mask, read together, their ratio the level.
    let plains_time =
        |s: f32| rec.formed[rec.region_of(s).min(rec.tesserae_from.saturating_sub(1))];
    let times: Vec<f32> = rec.score.iter().map(|&s| plains_time(s)).collect();
    let mut overturns = times.clone();
    overturns.sort_by(f32::total_cmp);
    overturns.dedup();
    let levels: Vec<(f32, PaddedField, PaddedField)> = overturns
        .iter()
        .map(|&t| {
            let mask: Vec<f32> = times
                .iter()
                .map(|&p| if p == t { 1.0 } else { 0.0 })
                .collect();
            let held: Vec<f32> = rec.filled.iter().zip(&mask).map(|(f, m)| f * m).collect();
            (
                t,
                PaddedField::new(coarse, &held),
                PaddedField::new(coarse, &mask),
            )
        })
        .collect();
    let target_at = |e: [f64; 3], time: f32| -> f32 {
        if let Some((_, held, mask)) = levels.iter().find(|l| l.0 == time) {
            let w = sample_bilinear(coarse, mask, e);
            if w > 1e-3 {
                return sample_bilinear(coarse, held, e) / w;
            }
        }
        sample_bilinear(coarse, &padded_filled, e)
    };
    // How far the old ground's smoothed height stands below the tesserae's edge, m:
    // beyond them the old ground sinking beneath the plains, which thicken away from the
    // edge in an onlap wedge, and within them the lava ponded at the edge's level below the
    // ground, filling the fabric's valleys deeper than that.
    let below = |e: [f64; 3]| rec.threshold - sample_bilinear(coarse, &padded_score, e);
    // The fabric about the tesserae, as deep below their edge as all its waves at their
    // highest could stand.
    let near_band = |fabric_bound: f64| 4.0 * fabric_bound * rec.strain.len().max(1) as f64;
    let (_, yield_compression) = crate::stages::stress::brittle_yields(rec.gravity);
    let highest = (1.0 + TESSERA_RELIEF_SCATTER) * yield_compression
        / (DENSITY_CRUST_ROCK * rec.gravity * TESSERA_BUCKLING_RATIO);
    let band = near_band(highest) as f32;
    let fabric = TesseraFabric::new(
        coarse,
        &rec.strain,
        body_radius,
        rec.gravity,
        grid.delta() * body_radius,
        seed,
        &|k| rec.threshold - rec.score[k] < band,
    );
    // The fabric over the tesserae, then beyond them as far below their edge as its highest
    // ridges stand above its median.
    let draw_within = |drawn: &mut Vec<f32>, from: f32, to: f32| {
        if fabric.is_empty() {
            return;
        }
        drawn.par_iter_mut().enumerate().for_each(|(idx, v)| {
            let d = grid.center(grid.cell(idx));
            let b = below(d);
            if b >= from && b < to {
                *v = fabric.at(d) as f32;
            }
        });
    };
    let mut drawn = vec![0.0f32; grid.cell_count()];
    draw_within(&mut drawn, f32::NEG_INFINITY, 0.0);
    // Centred so that half of the tesserae's fabric stands above the plains at the edge;
    // beyond its highest ridges' reach all of it is buried.
    let median = {
        let mut inside: Vec<f32> = (0..grid.cell_count())
            .step_by(7)
            .filter(|&idx| below(grid.center(grid.cell(idx))) < 0.0)
            .map(|idx| drawn[idx])
            .collect();
        if inside.is_empty() {
            0.0
        } else {
            let mid = inside.len() / 2;
            *inside.select_nth_unstable_by(mid, |a, b| a.total_cmp(b)).1
        }
    };
    let reach = (drawn.iter().fold(0.0f32, |m, &v| m.max(v)) - median).min(band);
    draw_within(&mut drawn, 0.0, reach);
    let bedrock = layers.bedrock.as_slice();
    let ages = layers.surface_age.as_slice();
    let deep = layers.deep.as_slice();
    let melt = layers.melt.as_slice();
    let n = grid.cell_count();
    let mut cells = Vec::with_capacity(n);
    // The overturns' lava rose everywhere they renewed, each place's from the mantle under
    // it: its composition at this resolution, not the coarse cells' tint of the floods after
    // it.
    let mut compositions = Vec::with_capacity(n);
    let mut taken = Vec::with_capacity(n);
    for idx in 0..n {
        let d = grid.center(grid.cell(idx));
        let own = sample_bilinear(coarse, &padded_fill, d);
        taken.push((idx as u32, own));
        // When its region's plains, and its tesserae, last formed (F296).
        let time = rec.plains_at(coarse, &padded_score, d);
        let formed = rec.formed_at(coarse, &padded_score, d);
        if smear.crater_time[idx] > time + 1e-6 {
            continue;
        }
        // The ground as its last overturn found it: the base less the smears not yet taken
        // out and less what the rises, shields and coronae have added since; both ride on
        // top of whatever is written.
        let remaining = smear.total[idx] - smear.done[idx];
        let later = deform.since_event(d, time as f64) as f32;
        let ground = bedrock[idx] - remaining - later;
        let carry = remaining - own + later;
        // The fabric less the burial at a cell; beyond its highest ridges' reach all of it
        // is buried.
        let exposed = |k: usize, e: [f64; 3]| -> (f32, f32) {
            let b = below(e);
            if b >= reach {
                (-b, b)
            } else {
                (drawn[k] - median - b, b)
            }
        };
        let plains_depth = |e: [f64; 3]| -> f32 {
            let target = target_at(e, time);
            (EPISODIC_FILL_FRACTION as f32 * (target - ground)).max(0.0) + EPISODIC_VENEER as f32
        };
        // A lid whose last overturn renewed every region keeps no tesserae (F300): plains
        // everywhere.
        if rec.tesserae_from >= rec.edges.len() {
            cells.push(CellWrite {
                idx: idx as u32,
                bedrock: ground + plains_depth(d) + carry,
                surface: true,
                age: time,
                unit: UNIT_LAVA,
                deep: 0.0,
                melt: 0.0,
                cover: 1.0,
            });
            compositions.push(overturn_composition(time as f64));
            continue;
        }
        // The share of the cell the lava left standing above its surface, by where the
        // fabric crosses it across the cell, read from the cells beside it (F131): a contact
        // drawn by share, not by whole cells.
        let (x, b) = exposed(idx, d);
        let v = x - plains_depth(d);
        let held = {
            let cell = grid.cell(idx);
            let side = |di: i32, dj: i32| {
                let k = grid.index(grid.neighbor(cell, di, dj));
                let e = grid.center(grid.cell(k));
                exposed(k, e).0 - plains_depth(e)
            };
            let gx = 0.5 * (side(1, 0) - side(-1, 0));
            let gy = 0.5 * (side(0, 1) - side(0, -1));
            let across = (gx * gx + gy * gy).sqrt();
            if across < 1e-6 {
                if v > 0.0 { 1.0 } else { 0.0 }
            } else {
                (0.5 + v / across).clamp(0.0, 1.0)
            }
        };
        // The old ground with its fabric, and the lava's surface: the plains' where it
        // laid them, ponded at the edge's level within the tesserae.
        let old = ground + x + b.min(0.0);
        let pond = ground + plains_depth(d) + b.min(0.0);
        let surface = old.max(pond) + carry;
        let whole = held >= 0.5;
        cells.push(CellWrite {
            idx: idx as u32,
            bedrock: surface,
            surface: true,
            age: if whole { ages[idx].min(formed) } else { time },
            unit: if whole { UNIT_TESSERA } else { UNIT_LAVA },
            deep: held * deep[idx],
            melt: held * melt[idx],
            cover: 1.0 - held,
        });
        compositions.push(if held < 1.0 {
            overturn_composition(time as f64)
        } else {
            0.0
        });
    }
    (
        StampOutput {
            deposit: None,
            cells,
            composition: None,
            stain: None,
            steep: None,
            ejecta_age: None,
            base: None,
            compositions: Some(compositions),
        },
        taken,
    )
}

/// Finishing pass of a pancake dome: the whole edifice, with a fractured top.
pub fn dome_compute(
    layers: &LayerSet,
    grid: &CubeSphere,
    body_radius: f64,
    dome: &Dome,
    stamper: &mut Stamper,
) -> StampOutput {
    let bedrock = layers.bedrock.as_slice();
    let ages = layers.surface_age.as_slice();
    let units = layers.unit_id.as_slice();
    let deep = layers.deep.as_slice();
    let melt = layers.melt.as_slice();
    let cover = layers.cover.as_slice();
    stamper.gather(grid, body_radius, dome.center, dome.radius);
    let noise = GradientNoise::new(dome.id ^ 0xD03E_0000);
    // The dome's lava where its edifice stands, by the share of each cell (F131).
    let edifice =
        |e: [f64; 3]| dome_profile(angle_between(dome.center, e) * body_radius / dome.radius);
    let f = body_radius / (0.25 * dome.radius);
    let time = dome.time_gyr as f32;
    let cells = stamper
        .footprint
        .iter()
        .map(|&(idx, r, _)| {
            let d = grid.center(grid.cell(idx));
            let x = r / dome.radius;
            let top = dome_profile(x);
            let cracks = DOME_TOP_ROUGHNESS
                * top
                * noise.fbm([d[0] * f, d[1] * f, d[2] * f], 1.0, 2, 2.0, 0.5);
            let h = dome.height * (top + cracks);
            let mut write = CellWrite {
                idx: idx as u32,
                bedrock: (bedrock[idx] as f64 + h) as f32,
                surface: false,
                age: 0.0,
                unit: 0,
                deep: 0.0,
                melt: 0.0,
                cover: 0.0,
            };
            if ages[idx] <= time + 1e-6 {
                let share = crate::fields::level_share(&edifice, d, 0.05, grid.delta()) as f32;
                if share > 0.0 {
                    lava_by_share(
                        &mut write, share, time, ages[idx], units[idx], deep[idx], melt[idx],
                        cover[idx],
                    );
                }
            }
            write
        })
        .collect();
    StampOutput {
        deposit: None,
        cells,
        composition: None,
        stain: None,
        steep: None,
        ejecta_age: None,
        base: None,
        compositions: None,
    }
}

/// Lava laid over `share` of a cell (F131): its cover over what was there by that share,
/// its time and unit where it takes more than half, the ground's own otherwise.
#[allow(clippy::too_many_arguments)]
pub fn lava_by_share(
    write: &mut CellWrite,
    share: f32,
    time: f32,
    age: f32,
    unit: u16,
    deep: f32,
    melt: f32,
    cover: f32,
) {
    let whole = share >= 0.5;
    write.surface = true;
    write.age = if whole { time } else { age };
    write.unit = if whole { UNIT_LAVA } else { unit };
    write.deep = (1.0 - share) * deep;
    write.melt = (1.0 - share) * melt;
    write.cover = share + (1.0 - share) * cover;
}

/// Fraction of the surface that is tessera.
pub fn tessera_fraction(layers: &LayerSet) -> f64 {
    let grid = layers.grid;
    let units = layers.unit_id.as_slice();
    let (mut t, mut total) = (0.0, 0.0);
    for (k, &u) in units.iter().enumerate() {
        let a = grid.cell_area(grid.cell(k));
        total += a;
        if u == UNIT_TESSERA {
            t += a;
        }
    }
    t / total
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn depression_filling_levels_pits_and_keeps_slopes() {
        let grid = CubeSphere::new(16);
        let n = grid.cell_count();
        // A tilted plane with a pit.
        let mut h: Vec<f32> = (0..n)
            .map(|k| (grid.center(grid.cell(k))[2] * 1000.0) as f32)
            .collect();
        let pit = grid.index(crate::grid::Cell {
            face: 0,
            i: 8,
            j: 8,
        });
        let before = h[pit];
        h[pit] = before - 500.0;
        let filled = fill_depressions(&grid, &h);
        // The pit fills to its lowest neighbour, the spill point.
        let cell = grid.cell(pit);
        let spill = [(1, 0), (-1, 0), (0, 1), (0, -1)]
            .iter()
            .map(|&(di, dj)| h[grid.index(grid.neighbor(cell, di, dj))])
            .fold(f32::INFINITY, f32::min);
        assert!(
            (filled[pit] - spill).abs() < 1e-3,
            "{} vs {spill}",
            filled[pit]
        );
        assert!(filled[pit] > before - 200.0);
        for k in 0..n {
            assert!(filled[k] >= h[k] - 1e-3);
            if k != pit {
                assert!((filled[k] - h[k]).abs() < 1e-3);
            }
        }
        // A liquid laid beside ground it does not cover spills onto it where that stands
        // lower than its spill point (F301): with one of the pit's neighbours lowered and
        // named an outlet, the pit fills only to it, and the outlet keeps its height.
        let low = grid.index(grid.neighbor(cell, 1, 0));
        h[low] = before - 450.0;
        let filled = fill_depressions_to(&grid, &h, |k| k == low);
        assert!(
            (filled[pit] - h[low]).abs() < 1e-3,
            "{} vs {}",
            filled[pit],
            h[low]
        );
        assert!((filled[low] - h[low]).abs() < 1e-3);
        // Without the outlet the two fill together to the lowest of their rim, two cells
        // from the pit at most, under 400 m down the plane.
        let closed = fill_depressions(&grid, &h);
        assert!(closed[pit] > h[low] + 1.0, "{} vs {}", closed[pit], h[low]);
        assert!((closed[pit] - closed[low]).abs() < 1e-3);
    }

    #[test]
    fn profiles_have_the_right_shape() {
        assert!((dome_profile(0.0) - 1.0).abs() < 1e-9);
        assert!(dome_profile(0.6) > 0.9);
        assert!(dome_profile(0.95) < 0.5);
        assert!(dome_profile(1.0) == 0.0);
    }
}
