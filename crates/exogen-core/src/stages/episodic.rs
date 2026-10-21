// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! The episodic lid (PLAN.md 4.3, Venus type).
//!
//! The lid founders once in a while and the whole surface is renewed. The event floods
//! the plains: every depression fills most of the way to its spill point, which takes
//! the old craters with it and leaves smooth lowlands under the regional slopes, while
//! the highest-standing terrain, a small fraction, survives as tessera, deformed into a
//! fabric of intersecting ridges. Afterwards the surface ages again under a sparse,
//! un-degraded crater population, with flow fields, coronae over the mantle's diapirs
//! (F145, the corona process) and steep-sided pancake domes on the plains.

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
use crate::stages::volcanism::{LavaSmear, overturn_composition};
use crate::units::{UNIT_LAVA, UNIT_TESSERA};

pub const VERSION: u64 = 7;

/// The resurfacing event, kept for the finishing pass.
#[derive(Clone, Debug)]
pub struct Resurfacing {
    pub time_gyr: f64,
    /// Thickness the event added per coarse cell, m.
    pub fill: Vec<f32>,
    /// The surface the plains were filled towards, per coarse cell: the terrain with its
    /// depressions filled.
    pub filled: Vec<f32>,
    /// Smoothed height, whose high end is the tessera.
    pub score: Vec<f32>,
    pub threshold: f32,
    /// Surface age the tessera keep, Gyr since formation.
    pub tessera_age: f32,
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
    let n = grid.cell_count();
    let mut filled = vec![f32::INFINITY; n];
    let mut visited = vec![false; n];
    let start = (0..n)
        .min_by(|&a, &b| heights[a].total_cmp(&heights[b]))
        .unwrap_or(0);
    let mut heap = BinaryHeap::new();
    heap.push(Reverse(Key(heights[start], start)));
    visited[start] = true;
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

/// The event's score threshold: the value a fraction of the area lies above.
fn upper_quantile(grid: &CubeSphere, score: &[f32], fraction: f64) -> f32 {
    let mut cells: Vec<(f32, f64)> = score
        .iter()
        .enumerate()
        .map(|(k, &s)| (s, grid.cell_area(grid.cell(k))))
        .collect();
    cells.sort_by(|a, b| b.0.total_cmp(&a.0));
    let total: f64 = cells.iter().map(|c| c.1).sum();
    let mut acc = 0.0;
    for (s, a) in cells {
        acc += a;
        if acc >= fraction * total {
            return s;
        }
    }
    f32::NEG_INFINITY
}

/// Tessera fabric: two sets of ridges crossing, as ridged noise, m.
pub fn tessera_texture(d: [f64; 3], noise: &GradientNoise, body_radius: f64) -> f64 {
    let f = body_radius / TESSERA_FABRIC_LENGTH;
    let a = noise.ridged([d[0] * f, d[1] * f * 0.35, d[2] * f], 1.0, 3, 2.0, 0.5);
    let b = noise.ridged(
        [d[1] * f * 1.3 + 7.0, d[2] * f * 0.4 + 3.0, d[0] * f * 1.1],
        1.0,
        3,
        2.0,
        0.5,
    );
    let g = body_radius / TESSERA_FOLD_LENGTH;
    let fold = noise.ridged(
        [d[2] * g + 11.0, d[0] * g * 0.5, d[1] * g * 1.2],
        1.0,
        2,
        2.0,
        0.5,
    );
    TESSERA_RELIEF * (a.max(b) - 0.5) + TESSERA_FOLD_RELIEF * (fold - 0.5)
}

/// The plateau uplift of the tessera at a smoothed-height score, m: full above the
/// threshold, ramping in over the edge width below it.
pub fn plateau_lift(score: f32, threshold: f32) -> f32 {
    let x = ((score - threshold) / (2.0 * TESSERA_EDGE_WIDTH as f32) + 0.5).clamp(0.0, 1.0);
    TESSERA_PLATEAU_UPLIFT as f32 * x * x * (3.0 - 2.0 * x)
}

/// Pancake dome profile at a fraction of the radius: flat top, steep sides.
pub fn dome_profile(x: f64) -> f64 {
    if x >= 1.0 {
        0.0
    } else {
        (1.0 - x.powi(6)).powf(0.6)
    }
}

/// The episodic lid through the loop: the event, then the coronae and domes.
pub struct EpisodicProcess {
    grid: CubeSphere,
    body_radius: f64,
    event_time: f64,
    /// The surface age the tessera keep, Gyr since formation (F121).
    tessera_age: f64,
    tessera_fraction: f64,
    done: bool,
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
    ) -> Option<Self> {
        if graph.str("interior.tectonic_regime") != "episodic_lid" {
            return None;
        }
        let body_radius = graph.f64("body.radius");
        let age = graph.f64("body.age");
        let event_time = (age - graph.f64("tectonics.resurfacing_event_age")).max(0.0);
        // The tessera keep the surface the overturn before the event laid (F121), or one
        // TESSERA_AGE_LEAD_GYR older than the event where the mantle gives none before it.
        let mantle = crate::mantle::evolve(&crate::mantle::MantleInputs::from_graph(graph));
        let tessera_age = mantle
            .overturns
            .iter()
            .rev()
            .find(|&&t| t < event_time - 1e-6)
            .copied()
            .unwrap_or((event_time - TESSERA_AGE_LEAD_GYR).max(0.0));
        let rng = Rng::new(graph.f64("seed.tectonics") as u64 ^ 0xE915_0D1C);
        // The coronae over the mantle's diapirs (F145).
        let coronae = crate::stages::corona::rock_coronae(graph, &mantle, rises);
        let n_domes = graph.f64("volcanism.dome_count").round() as u64;
        let mut domes = Vec::new();
        for k in 0..n_domes {
            let (lo, hi) = DOME_DIAMETER;
            let diameter = lo * (hi / lo).powf(rng.uniform(1_000 + k, 1));
            domes.push(Dome {
                center: rng.unit_direction(1_000 + k, 5),
                radius: 0.5 * diameter,
                height: DOME_HEIGHT * (0.6 + 0.8 * rng.uniform(1_000 + k, 3)),
                time_gyr: (event_time + (age - event_time) * rng.uniform(1_000 + k, 4)).min(age),
                id: 1_000 + k,
            });
        }
        // Novae: a share of the coronae send fractures radiating out beyond their rims.
        let cell = grid.delta() * body_radius;
        let mut novae = Vec::new();
        for c in &coronae {
            if rng.uniform(c.id, 7) >= NOVA_FRACTION {
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
                    start_gyr: c.time_gyr,
                    end_gyr: c.time_gyr,
                    id: 2_000_000 + draw,
                    kind: TroughKind::Fracture,
                    segments: Vec::new(),
                });
            }
        }
        Some(Self {
            grid,
            body_radius,
            event_time,
            tessera_age,
            tessera_fraction: graph.f64("tectonics.tessera_fraction"),
            done: false,
            resurfacing: None,
            coronae,
            domes,
            novae,
            stamper: Stamper::new(&grid),
        })
    }

    pub fn event_time(&self) -> f64 {
        self.event_time
    }

    /// Apply what falls in the step: the event, and the coarse part of the coronae.
    /// Returns the event's fill when it happens, for the lava ledger. The tessera are
    /// chosen on the height without the plume rises' domes, which are young dynamic
    /// topography, not the old thick crust.
    pub fn step(
        &mut self,
        layers: &mut LayerSet,
        t0: f64,
        t1: f64,
        dome: Option<&[f32]>,
    ) -> Option<Vec<f32>> {
        let grid = self.grid;
        let body_radius = self.body_radius;
        let mut event_fill = None;
        if !self.done && self.event_time > t0 && self.event_time <= t1 {
            self.done = true;
            let heights = layers.bedrock.as_slice().to_vec();
            let old_crust: Vec<f32> = match dome {
                Some(l) => heights.iter().zip(l).map(|(h, l)| h - l).collect(),
                None => heights.clone(),
            };
            let score = smooth(&grid, &old_crust, TESSERA_SMOOTHING_PASSES);
            let threshold = upper_quantile(&grid, &score, self.tessera_fraction);
            let filled = fill_depressions(&grid, &heights);
            let n = grid.cell_count();
            let mut fill = vec![0.0f32; n];
            let tessera_age = self.tessera_age as f32;
            {
                let bedrock = layers.bedrock.as_mut_slice();
                let ages = layers.surface_age.as_mut_slice();
                let units = layers.unit_id.as_mut_slice();
                let deep = layers.deep.as_mut_slice();
                let melt = layers.melt.as_mut_slice();
                let cover = layers.cover.as_mut_slice();
                let composition = layers.composition.as_mut_slice();
                let lava_age = layers.lava_age.as_mut_slice();
                let lava = overturn_composition(self.event_time);
                for k in 0..n {
                    let lift = plateau_lift(score[k], threshold);
                    if score[k] > threshold {
                        fill[k] = lift;
                        bedrock[k] += lift;
                        units[k] = UNIT_TESSERA;
                        ages[k] = ages[k].min(tessera_age);
                        cover[k] = 0.0;
                        continue;
                    }
                    let t = (EPISODIC_FILL_FRACTION as f32 * (filled[k] - heights[k])).max(0.0)
                        + EPISODIC_VENEER as f32
                        + lift;
                    fill[k] = t;
                    bedrock[k] += t;
                    ages[k] = self.event_time as f32;
                    units[k] = UNIT_LAVA;
                    deep[k] = 0.0;
                    melt[k] = 0.0;
                    cover[k] = 1.0;
                    // Its lava rose everywhere, each place's from the mantle under it.
                    composition[k] = lava;
                    lava_age[k] = self.event_time as f32;
                }
            }
            event_fill = Some(fill.clone());
            self.resurfacing = Some(Resurfacing {
                time_gyr: self.event_time,
                fill,
                filled,
                score,
                threshold,
                tessera_age,
            });
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
        event_fill
    }
}

/// Finishing pass of the event: the tessera boundary at the resolved resolution, the
/// plains refilled against the resolved ground, the tessera fabric. Returns the writes
/// and the smear taken out per cell.
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
    let noise = GradientNoise::new(seed ^ 0x7E55_E2A0);
    let bedrock = layers.bedrock.as_slice();
    let ages = layers.surface_age.as_slice();
    let deep = layers.deep.as_slice();
    let melt = layers.melt.as_slice();
    let time = rec.time_gyr as f32;
    let n = grid.cell_count();
    let mut cells = Vec::with_capacity(n);
    // The event's lava rose everywhere, each place's from the mantle under it: its
    // composition at this resolution, not the coarse cells' tint of the floods after it.
    let mut compositions = Vec::with_capacity(n);
    let lava = overturn_composition(rec.time_gyr);
    let mut taken = Vec::with_capacity(n);
    for idx in 0..n {
        let d = grid.center(grid.cell(idx));
        let own = sample_bilinear(coarse, &padded_fill, d);
        taken.push((idx as u32, own));
        if smear.crater_time[idx] > time + 1e-6 {
            continue;
        }
        // The ground as the event found it: the base less the smears not yet taken out
        // and less what the rises, shields and coronae have added since; both ride on
        // top of whatever is written.
        let remaining = smear.total[idx] - smear.done[idx];
        let later = deform.since_event(d, rec.time_gyr) as f32;
        let ground = bedrock[idx] - remaining - later;
        let carry = remaining - own + later;
        let score = sample_bilinear(coarse, &padded_score, d);
        let lift = plateau_lift(score, rec.threshold);
        // The tessera's share of the cell, from where its score's margin crosses it, and
        // the plains' the rest (F131): a contact drawn by share, not a staircase of cells.
        let held = crate::fields::level_share(
            &|e: [f64; 3]| sample_bilinear(coarse, &padded_score, e) as f64,
            d,
            rec.threshold as f64,
            grid.delta(),
        ) as f32;
        let tessera = if held > 0.0 {
            ground + lift + tessera_texture(d, &noise, body_radius) as f32 + carry
        } else {
            0.0
        };
        let plains = if held < 1.0 {
            let target = sample_bilinear(coarse, &padded_filled, d);
            ground
                + (EPISODIC_FILL_FRACTION as f32 * (target - ground)).max(0.0)
                + EPISODIC_VENEER as f32
                + lift
                + carry
        } else {
            0.0
        };
        let whole = held >= 0.5;
        cells.push(CellWrite {
            idx: idx as u32,
            bedrock: held * tessera + (1.0 - held) * plains,
            surface: true,
            age: if whole {
                ages[idx].min(rec.tessera_age)
            } else {
                time
            },
            unit: if whole { UNIT_TESSERA } else { UNIT_LAVA },
            deep: held * deep[idx],
            melt: held * melt[idx],
            cover: 1.0 - held,
        });
        compositions.push(if held < 1.0 { lava } else { 0.0 });
    }
    (
        StampOutput {
            cells,
            composition: None,
            stain: None,
            steep: None,
            ejecta_age: None,
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
        cells,
        composition: None,
        stain: None,
        steep: None,
        ejecta_age: None,
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
    }

    #[test]
    fn profiles_have_the_right_shape() {
        assert!((dome_profile(0.0) - 1.0).abs() < 1e-9);
        assert!(dome_profile(0.6) > 0.9);
        assert!(dome_profile(0.95) < 0.5);
        assert!(dome_profile(1.0) == 0.0);
        let noise = GradientNoise::new(3);
        let v: Vec<f64> = (0..200)
            .map(|k| {
                let a = k as f64 * 0.01;
                tessera_texture([a.cos(), a.sin(), 0.2], &noise, 6.05e6)
            })
            .collect();
        let max = v.iter().cloned().fold(f64::MIN, f64::max);
        let min = v.iter().cloned().fold(f64::MAX, f64::min);
        assert!(max - min > 0.3 * TESSERA_RELIEF, "{min}..{max}");
    }
}
