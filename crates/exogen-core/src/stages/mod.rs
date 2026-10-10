// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! The stage runner (PLAN.md 3.4, ADR 0002, ADR 0003).
//!
//! A run has two parts. The coarse forward loop steps the body through its age at a low
//! resolution and keeps keyframes. The finishing pass takes the keyframe at the output
//! epoch, resamples it to the finishing resolution and produces the fields that export
//! reads. Both parts are cached by a hash of the parameters and seeds they depend on, so
//! changing a late parameter does not rerun the loop.
//!
//! Phase 1 ships one placeholder process in the loop and a placeholder material model. The
//! structure is what later phases fill in.

pub mod materials;
pub mod placeholder_terrain;

use std::collections::HashMap;
use std::sync::Arc;
use std::time::{Duration, Instant};

use crate::fields::{Field, LayerSet};
use crate::grid::CubeSphere;
use crate::params::ParamGraph;
use crate::params::derive::physics::volatile;
use crate::rng::hash_words;

/// Loop steps by quality: preview, standard, high.
pub const LOOP_STEPS: [usize; 3] = [30, 90, 90];
/// Keyframes kept in memory, at most. The last step is always kept.
pub const MAX_KEYFRAMES: usize = 24;

#[derive(Clone, Debug)]
pub struct Progress {
    pub stage: &'static str,
    pub fraction: f32,
    pub message: String,
}

/// Receives progress and answers whether the run should stop.
pub trait ProgressSink: Sync {
    fn report(&self, _p: Progress) {}
    fn cancelled(&self) -> bool {
        false
    }
}

pub struct NoProgress;
impl ProgressSink for NoProgress {}

#[derive(Debug, thiserror::Error)]
pub enum RunError {
    #[error("run cancelled")]
    Cancelled,
}

/// The coarse state at one loop step.
#[derive(Clone)]
pub struct Keyframe {
    pub step: usize,
    pub time_gyr: f64,
    pub layers: LayerSet,
}

pub struct LoopOutput {
    pub grid: CubeSphere,
    pub keyframes: Vec<Keyframe>,
    pub steps: usize,
    pub age_gyr: f64,
    pub hash: u64,
}

impl LoopOutput {
    /// Index of the stored keyframe nearest to a time since formation.
    pub fn keyframe_at(&self, time_gyr: f64) -> usize {
        let mut best = 0;
        let mut best_d = f64::INFINITY;
        for (i, k) in self.keyframes.iter().enumerate() {
            let d = (k.time_gyr - time_gyr).abs();
            if d < best_d {
                best_d = d;
                best = i;
            }
        }
        best
    }
}

pub struct FinishOutput {
    pub grid: CubeSphere,
    pub layers: LayerSet,
    pub height: Field<f32>,
    pub albedo: Field<[f32; 3]>,
    /// Liquid surface level relative to the mean radius, metres; NaN without a liquid.
    pub sea_level: f64,
    pub ocean_fraction: f64,
    pub epoch_time_gyr: f64,
    pub keyframe_index: usize,
    pub hash: u64,
}

pub struct PipelineResult {
    pub loop_output: Arc<LoopOutput>,
    pub finish: Arc<FinishOutput>,
    pub timings: Vec<(&'static str, Duration)>,
    pub loop_cached: bool,
    pub finish_cached: bool,
}

/// Runs pipelines and caches their stage outputs in memory.
#[derive(Default)]
pub struct Runner {
    loop_cache: HashMap<u64, Arc<LoopOutput>>,
    finish_cache: HashMap<u64, Arc<FinishOutput>>,
}

const CACHE_ENTRIES: usize = 4;

impl Runner {
    pub fn new() -> Self {
        Self::default()
    }

    /// Everything the coarse loop reads.
    pub fn loop_key(graph: &ParamGraph) -> u64 {
        hash_words(&[
            graph.hash_of_group("star"),
            graph.hash_of_group("orbit"),
            graph.hash_of_group("body"),
            graph.hash_of_group("rotation"),
            graph.hash_of_group("interior"),
            graph.hash_of_group("tectonics"),
            graph.hash_of_group("atmosphere"),
            graph.hash_of_group("volatiles"),
            graph.hash_of(&[
                "seed.terrain",
                "detail.loop_resolution",
                "detail.quality",
                "multiplier.relief",
            ]),
            placeholder_terrain::VERSION,
        ])
    }

    /// Everything the finishing pass reads beyond the loop output.
    pub fn finish_key(graph: &ParamGraph, loop_hash: u64) -> u64 {
        hash_words(&[
            loop_hash,
            graph.hash_of(&[
                "detail.finishing_resolution",
                "output.epoch",
                "volatiles.inventory",
                "volatiles.working_liquid",
                "body.radius",
                "body.composition_class",
            ]),
            materials::VERSION,
        ])
    }

    pub fn run(
        &mut self,
        graph: &ParamGraph,
        sink: &dyn ProgressSink,
    ) -> Result<PipelineResult, RunError> {
        let mut timings = Vec::new();

        let lk = Self::loop_key(graph);
        let (loop_output, loop_cached) = match self.loop_cache.get(&lk) {
            Some(l) => (l.clone(), true),
            None => {
                let t = Instant::now();
                let l = Arc::new(run_coarse_loop(graph, sink)?);
                timings.push(("coarse loop", t.elapsed()));
                trim(&mut self.loop_cache);
                self.loop_cache.insert(lk, l.clone());
                (l, false)
            }
        };

        let fk = Self::finish_key(graph, loop_output.hash);
        let (finish, finish_cached) = match self.finish_cache.get(&fk) {
            Some(f) => (f.clone(), true),
            None => {
                let t = Instant::now();
                let face = graph.f64("detail.finishing_resolution") as u32;
                let epoch_time = (loop_output.age_gyr - graph.f64("output.epoch"))
                    .clamp(0.0, loop_output.age_gyr);
                let index = loop_output.keyframe_at(epoch_time);
                let f = Arc::new(finish(graph, &loop_output, face, index, sink)?);
                timings.push(("finishing", t.elapsed()));
                trim(&mut self.finish_cache);
                self.finish_cache.insert(fk, f.clone());
                (f, false)
            }
        };

        Ok(PipelineResult {
            loop_output,
            finish,
            timings,
            loop_cached,
            finish_cached,
        })
    }
}

fn trim<T>(cache: &mut HashMap<u64, T>) {
    if cache.len() >= CACHE_ENTRIES {
        cache.clear();
    }
}

fn quality_index(graph: &ParamGraph) -> usize {
    match graph.str("detail.quality") {
        "preview" => 0,
        "high" => 2,
        _ => 1,
    }
}

/// The coarse forward loop with the Phase 1 placeholder process.
pub fn run_coarse_loop(
    graph: &ParamGraph,
    sink: &dyn ProgressSink,
) -> Result<LoopOutput, RunError> {
    let face = graph.f64("detail.loop_resolution") as u32;
    let grid = CubeSphere::new(face);
    let age = graph.f64("body.age");
    let steps = LOOP_STEPS[quality_index(graph)];
    let dt = age / steps as f64;
    let keep_every = steps.div_ceil(MAX_KEYFRAMES).max(1);

    let process = placeholder_terrain::PlaceholderTerrain::new(graph, grid);
    let mut layers = process.init();
    let mut keyframes = Vec::with_capacity(MAX_KEYFRAMES + 1);
    keyframes.push(Keyframe {
        step: 0,
        time_gyr: 0.0,
        layers: layers.clone(),
    });

    for s in 1..=steps {
        if sink.cancelled() {
            return Err(RunError::Cancelled);
        }
        let t = dt * s as f64;
        process.step(&mut layers, t, dt);
        if s % keep_every == 0 || s == steps {
            keyframes.push(Keyframe {
                step: s,
                time_gyr: t,
                layers: layers.clone(),
            });
        }
        sink.report(Progress {
            stage: "coarse loop",
            fraction: s as f32 / steps as f32,
            message: format!("step {s} of {steps}, {t:.2} Gyr"),
        });
    }

    let hashes: Vec<u64> = keyframes.iter().map(|k| k.layers.hash()).collect();
    let hash = hash_words(&hashes);
    Ok(LoopOutput {
        grid,
        keyframes,
        steps,
        age_gyr: age,
        hash,
    })
}

/// The finishing pass for one keyframe at a given resolution.
pub fn finish(
    graph: &ParamGraph,
    loop_output: &LoopOutput,
    face: u32,
    keyframe_index: usize,
    sink: &dyn ProgressSink,
) -> Result<FinishOutput, RunError> {
    let kf = &loop_output.keyframes[keyframe_index];
    let grid = CubeSphere::new(face);
    sink.report(Progress {
        stage: "finishing",
        fraction: 0.1,
        message: format!("resampling to face size {face}"),
    });
    let mut layers = kf.layers.resample(grid);
    if sink.cancelled() {
        return Err(RunError::Cancelled);
    }

    sink.report(Progress {
        stage: "finishing",
        fraction: 0.5,
        message: "hydrosphere".to_string(),
    });
    let (sea_level, ocean_fraction) = hydrosphere_fill(graph, &mut layers);
    let height = layers.height();

    sink.report(Progress {
        stage: "finishing",
        fraction: 0.8,
        message: "materials".to_string(),
    });
    let albedo = materials::albedo(graph, &layers, &height);
    let hash = hash_words(&[layers.hash(), height.hash(), sea_level.to_bits()]);
    sink.report(Progress {
        stage: "finishing",
        fraction: 1.0,
        message: "done".to_string(),
    });
    Ok(FinishOutput {
        grid,
        layers,
        height,
        albedo,
        sea_level,
        ocean_fraction,
        epoch_time_gyr: kf.time_gyr,
        keyframe_index,
        hash,
    })
}

/// Fill the working liquid to its inventory volume (docs/parameters.md F26). Returns the
/// liquid level relative to the mean radius and the covered area fraction.
pub fn hydrosphere_fill(graph: &ParamGraph, layers: &mut LayerSet) -> (f64, f64) {
    let species = match graph.str("volatiles.working_liquid") {
        "water" | "ammonia_water" => "H2O",
        "methane_ethane" => "CH4",
        "nitrogen" => "N2",
        "sulfur_dioxide" => "SO2",
        "carbon_dioxide" => "CO2",
        _ => {
            layers.liquid = Field::new(layers.grid);
            return (f64::NAN, 0.0);
        }
    };
    let mut gel = graph
        .get("volatiles.inventory")
        .as_map()
        .and_then(|m| m.get(species).copied())
        .unwrap_or(0.0);
    if species == "CH4" {
        gel += graph
            .get("volatiles.inventory")
            .as_map()
            .and_then(|m| m.get("C2H6").copied())
            .unwrap_or(0.0);
    }
    let _ = volatile(species);
    let radius = graph.f64("body.radius");
    if gel <= 0.0 {
        layers.liquid = Field::new(layers.grid);
        return (f64::NAN, 0.0);
    }
    let solid = layers.height();
    let volume = gel * 4.0 * std::f64::consts::PI * radius * radius;
    let level = solid.fill_level(radius, volume);
    let grid = layers.grid;
    let mut covered = 0.0;
    let mut total = 0.0;
    let mut liquid = Field::new(grid);
    for c in grid.cells() {
        let depth = (level - solid.get(c) as f64).max(0.0);
        liquid.set(c, depth as f32);
        let a = grid.cell_area(c);
        total += a;
        if depth > 0.0 {
            covered += a;
        }
    }
    layers.liquid = liquid;
    (level, covered / total)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::presets::preset;

    fn earth_preview() -> ParamGraph {
        let (mut g, _) = preset("Earth").unwrap().to_graph();
        g.set("detail.quality", "preview").unwrap();
        g
    }

    #[test]
    fn the_same_seed_gives_identical_output_on_different_thread_counts() {
        let g = earth_preview();
        let run_with = |threads: usize| {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap();
            pool.install(|| {
                let mut r = Runner::new();
                let out = r.run(&g, &NoProgress).unwrap();
                (out.loop_output.hash, out.finish.hash)
            })
        };
        assert_eq!(run_with(1), run_with(4));
    }

    #[test]
    fn keyframes_change_over_time_and_the_cache_is_used() {
        let g = earth_preview();
        let mut r = Runner::new();
        let out = r.run(&g, &NoProgress).unwrap();
        assert!(!out.loop_cached && !out.finish_cached);
        let first = &out.loop_output.keyframes[0];
        let last = out.loop_output.keyframes.last().unwrap();
        assert!(last.time_gyr > first.time_gyr);
        assert_ne!(first.layers.hash(), last.layers.hash());
        let again = r.run(&g, &NoProgress).unwrap();
        assert!(again.loop_cached && again.finish_cached);

        let mut g2 = g.clone();
        g2.set("output.epoch", 2.0).unwrap();
        let other = r.run(&g2, &NoProgress).unwrap();
        assert!(
            other.loop_cached,
            "changing the epoch must not rerun the loop"
        );
        assert!(!other.finish_cached);
        assert!(other.finish.epoch_time_gyr < out.finish.epoch_time_gyr);
    }

    #[test]
    fn earth_gets_an_ocean_covering_most_of_the_surface() {
        let g = earth_preview();
        let mut r = Runner::new();
        let out = r.run(&g, &NoProgress).unwrap();
        assert!(out.finish.sea_level.is_finite());
        assert!(
            out.finish.ocean_fraction > 0.5 && out.finish.ocean_fraction < 0.9,
            "liquid fraction {}",
            out.finish.ocean_fraction
        );
    }
}
