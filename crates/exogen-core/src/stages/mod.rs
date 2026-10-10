// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! The stage runner (PLAN.md 3.4, ADR 0002, ADR 0003).
//!
//! A run has two parts. The coarse forward loop steps the body through its age at a low
//! resolution, emplacing the craters above the loop threshold at the steps they form, and
//! keeps keyframes. The finishing pass takes the keyframe at the output epoch, resamples it,
//! restores the resolved morphology of the loop craters, replays the sub-threshold crater
//! population in age order with each crater pre-degraded for its age, fills the working
//! liquid, computes optical maturity and rays, and assigns materials. Both parts are cached
//! by a hash of the parameters and seeds they depend on.

pub mod episodic;
pub mod impacts;
pub mod materials;
pub mod plume;
pub mod primordial_crust;
pub mod ridges;
pub mod tectonics;
pub mod volcanism;

use std::collections::HashMap;
use std::sync::Arc;
use std::time::{Duration, Instant};

use rayon::prelude::*;

use crate::constants::*;
use crate::crater::production::poisson;
use crate::crater::{
    CraterRecord, CraterTexture, EdgeShape, PopulationConfig, ProfileTable, RestoreParams,
    StampOutput, StampParams, Stamper, dimensions, ejecta_modulation, restore_compute,
    sample_population, stamp_compute, with_stamper,
};
use crate::fields::{Field, LayerSet};
use crate::geom::{dot, local_basis, offset};
use crate::grid::CubeSphere;
use crate::params::ParamGraph;
use crate::params::derive::physics::volatile;
use crate::rng::{Rng, hash_words};
use crate::thermal::ThermalHistory;
use crate::units::UNIT_SECONDARY_FIELD;

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
    /// Craters emplaced inside the loop, in order of formation.
    pub craters: Vec<CraterRecord>,
    pub loop_threshold: f64,
    /// The interior through the body's life, as the loop used it.
    pub thermal: ThermalHistory,
    /// Contraction scarps emplaced inside the loop, in order of formation.
    pub scarps: Vec<tectonics::ScarpRecord>,
    /// Lava floods emplaced inside the loop, in order of formation.
    pub floods: Vec<volcanism::FloodRecord>,
    /// Wrinkle ridges laid out on the plains at the end of the loop.
    pub ridges: Vec<ridges::RidgeRecord>,
    /// Plume rises, their shields, canyons and graben.
    pub rises: Vec<plume::Rise>,
    pub shields: Vec<plume::Shield>,
    pub canyons: Vec<plume::Trough>,
    pub graben: Vec<plume::Trough>,
    /// An episodic lid's resurfacing event, coronae and domes.
    pub resurfacing: Option<episodic::Resurfacing>,
    pub coronae: Vec<episodic::Corona>,
    pub domes: Vec<episodic::Dome>,
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
    /// Optical maturity, 0 fresh to 1 fully weathered.
    pub maturity: Field<f32>,
    /// Ray brightening, 0 to 1.
    pub rays: Field<f32>,
    /// Liquid surface level relative to the mean radius, metres; NaN without a liquid.
    pub sea_level: f64,
    pub ocean_fraction: f64,
    pub epoch_time_gyr: f64,
    pub keyframe_index: usize,
    /// Loop craters present at the epoch plus the replayed craters at or above eight
    /// finishing cells, for the sidecar and the statistics.
    pub craters: Vec<CraterRecord>,
    /// Total replayed craters stamped, including those below the recorded size.
    pub replayed_count: usize,
    /// Scarps present at the epoch.
    pub scarps: Vec<tectonics::ScarpRecord>,
    /// Floods present at the epoch.
    pub floods: Vec<volcanism::FloodRecord>,
    /// Wrinkle ridges present at the epoch.
    pub ridges: Vec<ridges::RidgeRecord>,
    /// Plume rises, and the shields, canyons and graben present at the epoch.
    pub rises: Vec<plume::Rise>,
    pub shields: Vec<plume::Shield>,
    pub canyons: Vec<plume::Trough>,
    pub graben: Vec<plume::Trough>,
    /// The resurfacing event if it lies before the epoch, and the coronae and domes
    /// present at it.
    pub resurfacing: Option<episodic::Resurfacing>,
    pub coronae: Vec<episodic::Corona>,
    pub domes: Vec<episodic::Dome>,
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
            graph.hash_of_group("cratering"),
            graph.hash_of_group("seeds"),
            graph.hash_of(&[
                "detail.loop_resolution",
                "detail.quality",
                "multiplier.relief",
                "multiplier.crater_density",
                "materials.deep_bedrock_depth",
            ]),
            primordial_crust::VERSION,
            impacts::VERSION,
            tectonics::VERSION,
            volcanism::VERSION,
            ridges::VERSION,
            plume::VERSION,
            episodic::VERSION,
        ])
    }

    /// Everything the finishing pass reads beyond the loop output.
    pub fn finish_key(graph: &ParamGraph, loop_hash: u64) -> u64 {
        hash_words(&[
            loop_hash,
            graph.hash_of_group("cratering"),
            graph.hash_of_group("materials"),
            graph.hash_of_group("surface"),
            graph.hash_of_group("seeds"),
            graph.hash_of(&[
                "detail.finishing_resolution",
                "output.epoch",
                "volatiles.inventory",
                "volatiles.working_liquid",
                "body.radius",
                "body.composition_class",
                "multiplier.crater_density",
            ]),
            materials::VERSION,
            FINISH_VERSION,
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

/// Version of the finishing pass itself, bumped when its behaviour changes.
pub const FINISH_VERSION: u64 = 18;

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

/// The coarse forward loop.
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

    let thermal = ThermalHistory::from_graph(graph);
    let crust = primordial_crust::PrimordialCrust::new(graph, grid);
    let mut layers = crust.init();
    let mut impacts = impacts::ImpactProcess::new(graph, grid);
    let mut contraction = tectonics::ContractionProcess::new(graph, grid, &thermal);
    let mut floods = volcanism::FloodProcess::new(graph, grid, &thermal);
    let mut plume = plume::PlumeProcess::new(graph, grid);
    let mut episodic = episodic::EpisodicProcess::new(graph, grid);
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
        let t0 = dt * (s - 1) as f64;
        let t1 = dt * s as f64;
        // The slow processes first: the dome lifts, the shields grow, the canyon deepens;
        // an episodic lid's event renews the plains.
        if let Some(p) = &mut plume {
            p.step(&mut layers, t0, t1);
        }
        if let Some(e) = &mut episodic
            && let Some(fill) = e.step(&mut layers, t0, t1)
            && let Some(f) = &mut floods
        {
            f.note_resurfacing(&fill, e.event_time());
        }
        let boost = plume.as_ref().map(|p| p.vent_boost());
        // Everything that happens in the step, in order of time: craters, scarps and
        // eruptions interleaved as they come, so what is older is cut by what is younger.
        let craters = impacts.sample(t0, t1, s as u64);
        let scarps = contraction
            .as_mut()
            .map_or_else(Vec::new, |c| c.sample(t0, t1, s as u64));
        let eruptions = floods.as_mut().map_or_else(Vec::new, |f| {
            f.sample(&layers, boost.as_deref(), t0, t1, s as u64)
        });
        let mut order: Vec<(f64, u8, usize)> =
            Vec::with_capacity(craters.len() + scarps.len() + eruptions.len());
        order.extend(craters.iter().enumerate().map(|(i, c)| (c.time_gyr, 0, i)));
        order.extend(scarps.iter().enumerate().map(|(i, c)| (c.time_gyr, 1, i)));
        order.extend(
            eruptions
                .iter()
                .enumerate()
                .map(|(i, e)| (e.time_gyr, 2, i)),
        );
        order.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)).then(a.2.cmp(&b.2)));
        let mut craters: Vec<Option<CraterRecord>> = craters.into_iter().map(Some).collect();
        let mut scarps: Vec<Option<tectonics::ScarpRecord>> =
            scarps.into_iter().map(Some).collect();
        for (_, kind, i) in order {
            match kind {
                0 => impacts.apply(&mut layers, craters[i].take().unwrap()),
                1 => {
                    if let Some(c) = &mut contraction {
                        c.apply(&mut layers, scarps[i].take().unwrap());
                    }
                }
                _ => {
                    if let Some(f) = &mut floods {
                        f.apply(&mut layers, eruptions[i]);
                    }
                }
            }
        }
        if s % keep_every == 0 || s == steps {
            keyframes.push(Keyframe {
                step: s,
                time_gyr: t1,
                layers: layers.clone(),
            });
        }
        sink.report(Progress {
            stage: "coarse loop",
            fraction: s as f32 / steps as f32,
            message: format!(
                "step {s} of {steps}, {t1:.2} Gyr, {} craters",
                impacts.records.len()
            ),
        });
    }

    let hashes: Vec<u64> = keyframes.iter().map(|k| k.layers.hash()).collect();
    let hash = hash_words(&hashes);
    let loop_threshold = impacts.loop_threshold();
    let craters = impacts.into_records(&layers);
    let scarps = contraction.map(|c| c.records).unwrap_or_default();
    let (resurfacing, coronae, domes) = match episodic {
        Some(e) => (e.resurfacing, e.coronae, e.domes),
        None => (None, Vec::new(), Vec::new()),
    };
    let (rises, shields, canyons, graben, load) = match plume {
        Some(p) => (p.rises, p.shields, p.canyons, p.graben, Some(p.load)),
        None => (Vec::new(), Vec::new(), Vec::new(), Vec::new(), None),
    };
    let (floods, ridges) = match floods {
        Some(f) => {
            // The ridges follow the lava thickness and, around a rise, its load.
            let mut field = f.thickness.clone();
            if let Some(load) = &load {
                for (v, l) in field.iter_mut().zip(load) {
                    *v += RIDGE_LOAD_WEIGHT as f32 * l;
                }
            }
            let ridges = ridges::plan_ridges(
                &grid,
                graph.f64("body.radius"),
                layers.unit_id.as_slice(),
                &field,
                &f.last_time,
                age,
                graph.f64("seed.volcanism") as u64,
            );
            (f.records, ridges)
        }
        None => (Vec::new(), Vec::new()),
    };
    Ok(LoopOutput {
        grid,
        keyframes,
        steps,
        age_gyr: age,
        craters,
        loop_threshold,
        thermal,
        scarps,
        floods,
        ridges,
        rises,
        shields,
        canyons,
        graben,
        resurfacing,
        coronae,
        domes,
        hash,
    })
}

/// What one stamp of the finishing pass produces: the writes, a record for the sidecar
/// when the crater is large enough, and the lava smear a flood's pass took out.
type WaveOutput = (
    StampOutput,
    Option<CraterRecord>,
    Option<(Vec<(u32, f32)>, [usize; 4])>,
);

/// The finishing pass for one keyframe at a given resolution.
pub fn finish(
    graph: &ParamGraph,
    loop_output: &LoopOutput,
    face: u32,
    keyframe_index: usize,
    sink: &dyn ProgressSink,
) -> Result<FinishOutput, RunError> {
    let kf = &loop_output.keyframes[keyframe_index];
    let present = kf.time_gyr;
    let grid = CubeSphere::new(face);
    let inputs = impacts::CraterInputs::from_graph(graph);
    let body_radius = inputs.body_radius;
    let cell = 2.0 * std::f64::consts::PI * body_radius / (4.0 * face as f64);

    sink.report(Progress {
        stage: "finishing",
        fraction: 0.05,
        message: format!("resampling to face size {face}"),
    });
    let mut layers = kf.layers.resample(grid);
    if sink.cancelled() {
        return Err(RunError::Cancelled);
    }

    // Loop craters that exist at this epoch. Their coarse form is in the resampled base;
    // their resolved detail is restored below, in age order among the replayed craters.
    let loop_craters: Vec<CraterRecord> = loop_output
        .craters
        .iter()
        .filter(|c| c.time_gyr <= present + 1e-9)
        .cloned()
        .collect();
    let scarps: Vec<tectonics::ScarpRecord> = loop_output
        .scarps
        .iter()
        .filter(|s| s.time_gyr <= present + 1e-9)
        .cloned()
        .collect();
    let floods: Vec<volcanism::FloodRecord> = loop_output
        .floods
        .iter()
        .filter(|f| f.time_gyr <= present + 1e-9)
        .cloned()
        .collect();
    let ridges: Vec<ridges::RidgeRecord> = loop_output
        .ridges
        .iter()
        .filter(|r| r.time_gyr <= present + 1e-9)
        .cloned()
        .collect();
    let shields: Vec<plume::Shield> = loop_output
        .shields
        .iter()
        .filter(|s| s.start_gyr <= present + 1e-9)
        .cloned()
        .collect();
    let canyons: Vec<plume::Trough> = loop_output
        .canyons
        .iter()
        .filter(|c| c.start_gyr <= present + 1e-9)
        .cloned()
        .collect();
    let graben: Vec<plume::Trough> = loop_output
        .graben
        .iter()
        .filter(|g| g.start_gyr <= present + 1e-9)
        .cloned()
        .collect();
    let resurfacing = loop_output
        .resurfacing
        .clone()
        .filter(|r| r.time_gyr <= present + 1e-9);
    let coronae: Vec<episodic::Corona> = loop_output
        .coronae
        .iter()
        .filter(|c| c.time_gyr <= present + 1e-9)
        .cloned()
        .collect();
    let domes: Vec<episodic::Dome> = loop_output
        .domes
        .iter()
        .filter(|d| d.time_gyr <= present + 1e-9)
        .cloned()
        .collect();
    let base: Vec<f32> = layers.bedrock.as_slice().to_vec();
    let base_deep: Vec<f32> = layers.deep.as_slice().to_vec();
    let mut smear = volcanism::LavaSmear::new(
        loop_output.grid,
        grid,
        &floods,
        resurfacing.as_ref().map(|r| r.fill.as_slice()),
    );
    if !floods.is_empty() || resurfacing.is_some() {
        with_stamper(&grid, |stamper| {
            for c in &loop_craters {
                let dim = dimensions(c.diameter, &inputs.scaling);
                let radius = 1.05 * CONTINUOUS_EJECTA_MULTIPLE * dim.radius * c.ellipticity.sqrt();
                stamper.gather(&grid, body_radius, c.center, radius);
                let t = c.time_gyr as f32;
                for &(idx, _, _) in &stamper.footprint {
                    smear.crater_time[idx] = smear.crater_time[idx].max(t);
                }
            }
        });
    }
    let coarse_cell =
        2.0 * std::f64::consts::PI * body_radius / (4.0 * loop_output.grid.face_size() as f64);

    // Statistical sample of the sub-threshold population, from the surface ages the loop left.
    let d_min = graph
        .f64("cratering.min_diameter")
        .max(FINISHING_MIN_CELLS * cell);
    let cfg = PopulationConfig {
        pf: inputs.pf,
        flux: inputs.flux,
        flux_scale: inputs.flux_scale,
        body_radius,
        d_min,
        d_max: loop_output.loop_threshold,
        bins_per_octave: BINS_PER_OCTAVE,
        saturation_fraction: inputs.saturation_fraction,
        retention_age: inputs.retention_age,
        oblique_fraction: inputs.oblique_fraction,
        seed: inputs.seed,
    };
    sink.report(Progress {
        stage: "finishing",
        fraction: 0.15,
        message: "sampling the crater population".to_string(),
    });
    let population = sample_population(&grid, &layers.surface_age, present, &cfg);
    let secondaries = sample_secondaries(&loop_craters, &inputs, d_min, present);
    let replay: Vec<(crate::crater::CraterSample, bool)> = population
        .into_iter()
        .map(|c| (c, false))
        .chain(secondaries.into_iter().map(|c| (c, true)))
        .collect();

    // One age-ordered sequence of loop craters and replayed craters. A loop crater goes
    // before anything else of the same time, so its own secondaries land on it; craters it
    // destroyed are stamped before it and removed by its restore.
    enum Event {
        Loop(usize),
        Scarp(usize),
        Flood(usize),
        Ridge(usize),
        Shield(usize),
        Canyon(usize),
        Graben(usize),
        Resurfacing,
        Corona(usize),
        Dome(usize),
        Replay(usize),
    }
    let mut events: Vec<(f64, u8, u64, Event)> = Vec::with_capacity(
        loop_craters.len()
            + scarps.len()
            + floods.len()
            + ridges.len()
            + shields.len()
            + canyons.len()
            + graben.len()
            + replay.len(),
    );
    for (k, sh) in shields.iter().enumerate() {
        events.push((sh.start_gyr, 0, sh.id ^ 0x5B1E_0000_0000, Event::Shield(k)));
    }
    for (k, c) in canyons.iter().enumerate() {
        events.push((
            c.end_gyr.min(present),
            0,
            c.id ^ 0xCA17_0000_0000,
            Event::Canyon(k),
        ));
    }
    for (k, g) in graben.iter().enumerate() {
        events.push((g.start_gyr, 0, g.id ^ 0x6A8B_0000_0000, Event::Graben(k)));
    }
    if let Some(r) = &resurfacing {
        events.push((r.time_gyr, 0, 0xE915_0000_0000, Event::Resurfacing));
    }
    for (k, c) in coronae.iter().enumerate() {
        events.push((c.time_gyr, 0, c.id ^ 0xC020_0000_0000, Event::Corona(k)));
    }
    for (k, d) in domes.iter().enumerate() {
        events.push((d.time_gyr, 0, d.id ^ 0xD03E_0000_0000, Event::Dome(k)));
    }
    for (k, r) in ridges.iter().enumerate() {
        events.push((r.time_gyr, 0, r.id ^ 0x51D6_0000_0000, Event::Ridge(k)));
    }
    for (k, f) in floods.iter().enumerate() {
        events.push((f.time_gyr, 0, f.id ^ 0x7F10_0000_0000, Event::Flood(k)));
    }
    for (k, s) in scarps.iter().enumerate() {
        events.push((s.time_gyr, 0, s.id ^ 0x5CA2_0000_0000, Event::Scarp(k)));
    }
    for (i, c) in loop_craters.iter().enumerate() {
        events.push((c.time_gyr, 0, c.id, Event::Loop(i)));
    }
    for (j, (c, _)) in replay.iter().enumerate() {
        events.push((c.time_gyr, 1, c.key, Event::Replay(j)));
    }
    events.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)).then(a.2.cmp(&b.2)));
    sink.report(Progress {
        stage: "finishing",
        fraction: 0.3,
        message: format!("stamping {} craters", events.len()),
    });
    if sink.cancelled() {
        return Err(RunError::Cancelled);
    }

    let record_min = 8.0 * cell;
    let mut records = loop_craters.clone();
    let total = events.len().max(1);
    let smooth_sigma = 0.5 * coarse_cell;

    // Footprint caps: centre and angular radius with a cell of margin, for wave building.
    let caps: Vec<([f64; 3], f64)> = events
        .par_iter()
        .map(|(_, _, _, event)| {
            let crater_cap = |center: [f64; 3], diameter: f64, ellipticity: f64| {
                let dim = dimensions(diameter, &inputs.scaling);
                (
                    center,
                    dim.ejecta_radius * ellipticity.sqrt() / body_radius + grid.delta(),
                )
            };
            match event {
                Event::Loop(li) => {
                    let c = &loop_craters[*li];
                    crater_cap(c.center, c.diameter, c.ellipticity)
                }
                Event::Replay(j) => {
                    let c = &replay[*j].0;
                    crater_cap(c.center, c.diameter, c.ellipticity)
                }
                Event::Scarp(k) => {
                    let s = &scarps[*k];
                    (s.midpoint(), s.cap_radius(body_radius) + 2.0 * grid.delta())
                }
                Event::Flood(k) => {
                    let (c, r) = floods[*k].cap(&loop_output.grid);
                    (c, r + 2.5 * loop_output.grid.delta() + grid.delta())
                }
                Event::Ridge(k) => {
                    let r = &ridges[*k];
                    (r.midpoint(), r.cap_radius(body_radius) + 2.0 * grid.delta())
                }
                Event::Shield(k) => {
                    let sh = &shields[*k];
                    (sh.center, sh.cap_radius(body_radius) + 2.0 * grid.delta())
                }
                Event::Canyon(k) => {
                    let c = &canyons[*k];
                    (c.midpoint(), c.cap_radius(body_radius) + 2.0 * grid.delta())
                }
                Event::Graben(k) => {
                    let g = &graben[*k];
                    (g.midpoint(), g.cap_radius(body_radius) + 2.0 * grid.delta())
                }
                Event::Resurfacing => ([0.0, 0.0, 1.0], std::f64::consts::PI),
                Event::Corona(k) => {
                    let c = &coronae[*k];
                    (c.center, c.cap_radius(body_radius) + 2.0 * grid.delta())
                }
                Event::Dome(k) => {
                    let d = &domes[*k];
                    (d.center, d.cap_radius(body_radius) + 2.0 * grid.delta())
                }
            }
        })
        .collect();
    let waves = build_waves(&caps);

    let mut done = 0usize;
    let mut next_report = 0usize;
    for wave in &waves {
        // Stamps in a wave have disjoint footprints: compute them in parallel against the
        // state before the wave, then apply them in order.
        let outputs: Vec<Option<WaveOutput>> = wave
            .par_iter()
            .map(|&ei| {
                with_stamper(&grid, |stamper| match &events[ei].3 {
                    Event::Loop(li) => {
                        let c = &loop_craters[*li];
                        let dim = dimensions(c.diameter, &inputs.scaling);
                        // The detail the coarse grid could not carry: the resolved profile
                        // minus its coarse-cell-smoothed version.
                        let fresh = ProfileTable::new(dim, 0.0);
                        let kappa_t = inputs.diffusivity_m2_per_myr
                            * inputs.flux.equivalent_gyr(c.time_gyr, present)
                            * 1e3;
                        let degraded = ProfileTable::new(dim, kappa_t);
                        let smooth =
                            ProfileTable::new(dim, kappa_t + 0.5 * smooth_sigma * smooth_sigma);
                        let mut detail = degraded;
                        for (v, s) in detail.values.iter_mut().zip(&smooth.values) {
                            *v -= s;
                        }
                        let deep =
                            impacts::reaches_deep_bedrock(c.diameter, inputs.deep_bedrock_depth);
                        let params = RestoreParams {
                            time_gyr: c.time_gyr as f32,
                            unit_floor: impacts::floor_unit(dim.class, deep, c.on_lava),
                            unit_ejecta: impacts::ejecta_unit(deep, c.on_lava),
                            azimuth: c.azimuth,
                            angle: c.angle,
                            ellipticity: c.ellipticity,
                            detail_weight: 1.0,
                            seed: c.id,
                            deep,
                            fine_blur: detail.sigma,
                            coarse_blur: smooth.sigma,
                        };
                        let out = restore_compute(
                            &layers,
                            &grid,
                            body_radius,
                            c.center,
                            &base,
                            &base_deep,
                            &fresh,
                            &detail,
                            &params,
                            stamper,
                        );
                        Some((out, None, None))
                    }
                    Event::Flood(k) => {
                        let (out, taken, stats) = volcanism::flood_compute(
                            &layers,
                            &grid,
                            &loop_output.grid,
                            &floods[*k],
                            &smear,
                        );
                        Some((out, None, Some((taken, stats))))
                    }
                    Event::Shield(k) => {
                        let out = plume::shield_compute(
                            &layers,
                            &grid,
                            body_radius,
                            &shields[*k],
                            stamper,
                        );
                        Some((out, None, None))
                    }
                    Event::Canyon(k) => {
                        let c = &canyons[*k];
                        let kappa_t = inputs.diffusivity_m2_per_myr
                            * inputs.flux.equivalent_gyr(c.end_gyr.min(present), present)
                            * 1e3;
                        let out = plume::trough_compute(
                            &layers,
                            &grid,
                            body_radius,
                            &base,
                            c,
                            c.end_gyr.min(present) as f32,
                            (4.0 * kappa_t).sqrt(),
                            Some((
                                coarse_cell,
                                (4.0 * kappa_t + 2.0 * smooth_sigma * smooth_sigma).sqrt(),
                            )),
                            stamper,
                        );
                        Some((out, None, None))
                    }
                    Event::Resurfacing => {
                        let r = resurfacing.as_ref().expect("event present");
                        let (out, taken) = episodic::resurfacing_compute(
                            &layers,
                            &grid,
                            &loop_output.grid,
                            body_radius,
                            r,
                            &smear,
                            inputs.seed,
                        );
                        Some((out, None, Some((taken, [0; 4]))))
                    }
                    Event::Corona(k) => {
                        let c = &coronae[*k];
                        let kappa_t = inputs.diffusivity_m2_per_myr
                            * inputs.flux.equivalent_gyr(c.time_gyr, present)
                            * 1e3;
                        let out = episodic::corona_compute(
                            &layers,
                            &grid,
                            body_radius,
                            c,
                            (4.0 * kappa_t).sqrt(),
                            stamper,
                        );
                        Some((out, None, None))
                    }
                    Event::Dome(k) => {
                        let out = episodic::dome_compute(
                            &layers,
                            &grid,
                            body_radius,
                            &domes[*k],
                            stamper,
                        );
                        Some((out, None, None))
                    }
                    Event::Graben(k) => {
                        let g = &graben[*k];
                        let kappa_t = inputs.diffusivity_m2_per_myr
                            * inputs.flux.equivalent_gyr(g.start_gyr, present)
                            * 1e3;
                        let out = plume::trough_compute(
                            &layers,
                            &grid,
                            body_radius,
                            &base,
                            g,
                            g.start_gyr as f32,
                            (4.0 * kappa_t).sqrt(),
                            None,
                            stamper,
                        );
                        Some((out, None, None))
                    }
                    Event::Ridge(k) => {
                        let r = &ridges[*k];
                        let kappa_t = inputs.diffusivity_m2_per_myr
                            * inputs.flux.equivalent_gyr(r.time_gyr, present)
                            * 1e3;
                        let out = ridges::ridge_compute(
                            &layers,
                            &grid,
                            body_radius,
                            r,
                            (4.0 * kappa_t).sqrt(),
                            stamper,
                        );
                        Some((out, None, None))
                    }
                    Event::Scarp(k) => {
                        // The face the coarse grid could not carry: the profile at its
                        // degraded width minus the coarse-smoothed one.
                        let s = &scarps[*k];
                        let kappa_t = inputs.diffusivity_m2_per_myr
                            * inputs.flux.equivalent_gyr(s.time_gyr, present)
                            * 1e3;
                        let fine = (4.0 * kappa_t).sqrt();
                        let coarse = (4.0 * kappa_t + 2.0 * smooth_sigma * smooth_sigma).sqrt();
                        let out = tectonics::scarp_compute(
                            &layers,
                            &grid,
                            body_radius,
                            s,
                            fine,
                            Some(coarse),
                            stamper,
                        );
                        Some((out, None, None))
                    }
                    Event::Replay(j) => {
                        let (c, secondary) = &replay[*j];
                        let dim = dimensions(c.diameter, &inputs.scaling);
                        // Degradation by the gardening that followed the crater: the
                        // diffusivity at the present flux times the equivalent flux since.
                        let kappa_t = inputs.diffusivity_m2_per_myr
                            * inputs.flux.equivalent_gyr(c.time_gyr, present)
                            * 1e3;
                        let table = ProfileTable::new(dim, kappa_t);
                        if table.is_erased() {
                            return None;
                        }
                        let deep = !*secondary
                            && impacts::reaches_deep_bedrock(c.diameter, inputs.deep_bedrock_depth);
                        let on_lava = impacts::target_is_lava(&layers, &grid, c.center);
                        let params = StampParams {
                            time_gyr: c.time_gyr as f32,
                            unit_floor: if *secondary && !on_lava {
                                UNIT_SECONDARY_FIELD
                            } else {
                                impacts::floor_unit(dim.class, deep, on_lava)
                            },
                            unit_ejecta: impacts::ejecta_unit(deep, on_lava),
                            azimuth: c.azimuth,
                            angle: c.angle,
                            ellipticity: c.ellipticity,
                            amplitude: if *secondary {
                                SECONDARY_DEPTH_FACTOR
                            } else {
                                1.0
                            },
                            replace_floor: true,
                            set_age: true,
                            bury: true,
                            seed: c.key,
                            deep,
                        };
                        let out = stamp_compute(
                            &layers,
                            &grid,
                            body_radius,
                            c.center,
                            &table,
                            &params,
                            stamper,
                        );
                        let record =
                            (c.diameter >= record_min && !*secondary).then_some(CraterRecord {
                                center: c.center,
                                diameter: c.diameter,
                                time_gyr: c.time_gyr,
                                azimuth: c.azimuth,
                                angle: c.angle,
                                ellipticity: c.ellipticity,
                                class: dim.class,
                                from_loop: false,
                                relief: (dim.depth + dim.rim_height) as f32,
                                visible: true,
                                on_lava,
                                id: c.key,
                            });
                        Some((out, record, None))
                    }
                })
            })
            .collect();
        for (out, record, taken) in outputs.into_iter().flatten() {
            if let Some((t, stats)) = &taken {
                smear.note(t, &out, *stats);
            }
            out.apply(&mut layers);
            if let Some(r) = record {
                records.push(r);
            }
        }
        done += wave.len();
        if done >= next_report {
            if sink.cancelled() {
                return Err(RunError::Cancelled);
            }
            sink.report(Progress {
                stage: "finishing",
                fraction: 0.3 + 0.4 * done as f32 / total as f32,
                message: format!("stamping craters {done} of {total}"),
            });
            next_report = done + 10_000;
        }
    }
    let replayed_count = replay.len();
    if std::env::var_os("EXOGEN_DEBUG").is_some() {
        eprintln!(
            "flood passes: lava {} above {} outside {} skipped {}",
            smear.stats[0], smear.stats[1], smear.stats[2], smear.stats[3]
        );
    }

    sink.report(Progress {
        stage: "finishing",
        fraction: 0.72,
        message: "hydrosphere".to_string(),
    });
    let (sea_level, ocean_fraction) = hydrosphere_fill(graph, &mut layers);
    let height = layers.height();

    sink.report(Progress {
        stage: "finishing",
        fraction: 0.8,
        message: "weathering and rays".to_string(),
    });
    let tau_sw = graph.f64("surface.space_weathering_timescale");
    let maturity = Field::from_fn(grid, |c, _| {
        let exposure = (present - layers.surface_age.get(c) as f64).max(0.0);
        (1.0 - (-exposure / tau_sw).exp()) as f32
    });
    let rays = ray_field(&grid, body_radius, &records, &inputs, present, cell);

    sink.report(Progress {
        stage: "finishing",
        fraction: 0.92,
        message: "materials".to_string(),
    });
    let albedo = materials::albedo(graph, &layers, &maturity, &rays);
    let hash = hash_words(&[
        layers.hash(),
        height.hash(),
        maturity.hash(),
        rays.hash(),
        sea_level.to_bits(),
    ]);
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
        maturity,
        rays,
        sea_level,
        ocean_fraction,
        epoch_time_gyr: present,
        keyframe_index,
        craters: records,
        replayed_count,
        scarps,
        floods,
        ridges,
        rises: loop_output.rises.clone(),
        shields,
        canyons,
        graben,
        resurfacing,
        coronae,
        domes,
        hash,
    })
}

/// Secondary craters of the loop craters: clusters at a few primary radii, formed at the
/// primary's time, sized as a small fraction of the primary.
fn sample_secondaries(
    primaries: &[CraterRecord],
    inputs: &impacts::CraterInputs,
    d_min: f64,
    present: f64,
) -> Vec<crate::crater::CraterSample> {
    let rng = Rng::new(inputs.seed ^ 0x5EC0_7DA2);
    let mut out = Vec::new();
    let (f_lo, f_hi) = SECONDARY_SIZE_FRACTION;
    let (r_lo, r_hi) = SECONDARY_DISTANCE_MULTIPLE;
    let slope = 3.0f64;
    for p in primaries {
        if p.time_gyr > present {
            continue;
        }
        let expected =
            SECONDARIES_PER_100KM * (inputs.secondary_fraction / 0.2) * (p.diameter / 1.0e5);
        let (count, _) = poisson(&rng, p.id, 0, expected);
        if count == 0 {
            continue;
        }
        let (east, north) = local_basis(p.center);
        let radius = 0.5 * p.diameter;
        let groups = (count / SECONDARY_GROUP_SIZE).max(1);
        for j in 0..count as u64 {
            let group = j % groups as u64;
            let rank = (j / groups as u64) as f64;
            let base = 8 + group * 4;
            let dist0 = radius * (r_lo + (r_hi - r_lo) * rng.uniform(p.id, base));
            let phi0 = 2.0 * std::f64::consts::PI * rng.uniform(p.id, base + 1);
            let chain = rng.uniform(p.id, base + 2) < SECONDARY_CHAIN_FRACTION;
            let draw = 1_000 + j * 4;
            // A chain runs outward along the radial from the primary; a cluster scatters
            // about its centre.
            let (dist, phi) = if chain {
                let spacing =
                    radius * SECONDARY_CHAIN_SPACING * (0.7 + 0.6 * rng.uniform(p.id, draw));
                let lateral = 0.03 * radius * (rng.uniform(p.id, draw + 1) - 0.5);
                let d = dist0 + rank * spacing;
                (d, phi0 + lateral / d.max(1.0))
            } else {
                let jitter_r = 0.3 * radius * (rng.uniform(p.id, draw) - 0.5);
                let jitter_phi =
                    0.5 * (rng.uniform(p.id, draw + 1) - 0.5) * radius / dist0.max(1.0);
                (dist0 + jitter_r, phi0 + jitter_phi)
            };
            let center = offset(p.center, east, north, dist / inputs.body_radius, phi);
            let u = rng.uniform(p.id, draw + 2);
            let q = (f_hi / f_lo).powf(-slope);
            let diameter = p.diameter * f_lo * (1.0 - u * (1.0 - q)).powf(-1.0 / slope);
            if diameter < d_min {
                continue;
            }
            out.push(crate::crater::CraterSample {
                center,
                diameter,
                time_gyr: p.time_gyr,
                azimuth: phi,
                angle: SECONDARY_IMPACT_ANGLE,
                ellipticity: 1.2,
                key: (p.id << 16) | j,
            });
        }
    }
    out
}

/// Largest number of stamps in one wave and the largest share of the sphere their caps
/// may cover, which bounds the memory of a wave's pending writes.
const WAVE_MAX_STAMPS: usize = 2048;
const WAVE_MAX_AREA: f64 = 0.25;
/// Rayed craters computed per parallel batch.
const RAY_BATCH: usize = 64;

/// Group consecutive stamps into waves of non-overlapping footprints. Stamps in one wave
/// touch disjoint cells, so they can be computed in parallel and applied in any order,
/// while every overlapping pair keeps its order because a stamp never joins a wave that
/// holds an earlier stamp it overlaps.
fn build_waves(caps: &[([f64; 3], f64)]) -> Vec<Vec<usize>> {
    let mut waves = Vec::new();
    let mut current: Vec<usize> = Vec::new();
    let mut area = 0.0;
    for (i, &(c, rho)) in caps.iter().enumerate() {
        let share = 0.5 * (1.0 - rho.min(std::f64::consts::PI).cos());
        let close = current.len() >= WAVE_MAX_STAMPS
            || area + share > WAVE_MAX_AREA
            || current.iter().any(|&j| {
                let (cj, rj) = caps[j];
                dot(c, cj) > (rho + rj).min(std::f64::consts::PI).cos()
            });
        if close && !current.is_empty() {
            waves.push(std::mem::take(&mut current));
            area = 0.0;
        }
        current.push(i);
        area += share;
    }
    if !current.is_empty() {
        waves.push(current);
    }
    waves
}

fn ray_field(
    grid: &CubeSphere,
    body_radius: f64,
    craters: &[CraterRecord],
    inputs: &impacts::CraterInputs,
    present: f64,
    cell: f64,
) -> Field<f32> {
    let mut rays = Field::new(*grid);
    let out = rays.as_mut_slice();
    // Each crater's ray cells are computed in parallel and merged by maximum, which is
    // order-independent, so the result does not depend on the thread count.
    for batch in craters.chunks(RAY_BATCH) {
        let parts: Vec<Vec<(u32, f32)>> = batch
            .par_iter()
            .map(|c| {
                with_stamper(grid, |stamper| {
                    ray_cells(grid, body_radius, c, inputs, present, cell, stamper)
                })
            })
            .collect();
        for part in parts {
            for (idx, v) in part {
                let o = &mut out[idx as usize];
                if v > *o {
                    *o = v;
                }
            }
        }
    }
    rays
}

/// The ray brightening one crater adds: full inside the rim, fading to half at the
/// continuous-ejecta edge and mottled there, then rays of constant linear width whose
/// strength follows the oblique ejecta asymmetry, all scaled by freshness.
fn ray_cells(
    grid: &CubeSphere,
    body_radius: f64,
    c: &CraterRecord,
    inputs: &impacts::CraterInputs,
    present: f64,
    cell: f64,
    stamper: &mut Stamper,
) -> Vec<(u32, f32)> {
    let lifetime = inputs.ray_lifetime.max(1e-6);
    let age = present - c.time_gyr;
    if age < 0.0 || age >= lifetime || c.diameter < RAY_MIN_CELLS * cell {
        return Vec::new();
    }
    let tau = 2.0 * std::f64::consts::PI;
    let rng = Rng::new(inputs.seed ^ 0x9A7F_3C11);
    let (n_lo, n_hi) = RAYS_PER_CRATER;
    let (l_lo, l_hi) = RAY_LENGTH_MULTIPLE;
    let (w_lo, w_hi) = RAY_WIDTH_FRACTION;
    let freshness = (1.0 - age / lifetime).powi(2);
    let radius = 0.5 * c.diameter;
    let edge_shape = EdgeShape::new(c.id);
    let texture = (c.diameter >= RAY_TEXTURE_MIN_CELLS * cell)
        .then(|| CraterTexture::new(c.id, &dimensions(c.diameter, &inputs.scaling)));
    let n = n_lo + ((n_hi - n_lo) as f64 * rng.uniform(c.id, 0)) as u32;
    // Rays keep a roughly constant linear width, so their angular width shrinks with
    // distance, and large craters' rays are shorter in radii than small craters' rays.
    let length_scale = (RAY_REFERENCE_DIAMETER / c.diameter).sqrt().min(1.0);
    let ray_params: Vec<(f64, f64, f64, f64)> = (0..n)
        .map(|i| {
            let d = 4 + 3 * i as u64;
            let phi = tau * rng.uniform(c.id, d);
            let width = radius * (w_lo + (w_hi - w_lo) * rng.uniform(c.id, d + 1));
            let length = radius * (l_lo + (l_hi - l_lo) * rng.uniform(c.id, d + 2)) * length_scale;
            let strength = ejecta_modulation(phi, c.azimuth, c.angle);
            (phi, width, length, strength)
        })
        .collect();
    let max_len = ray_params
        .iter()
        .map(|p| p.2)
        .fold(0.0, f64::max)
        .max(CONTINUOUS_EJECTA_MULTIPLE * radius);
    stamper.gather(grid, body_radius, c.center, max_len);
    let ejecta_edge = CONTINUOUS_EJECTA_MULTIPLE * radius;
    let mut cells = Vec::with_capacity(stamper.footprint.len());
    for &(idx, r, phi) in &stamper.footprint {
        let edge = ejecta_edge * edge_shape.factor(phi);
        let mut v: f64 = if r <= radius {
            1.0
        } else if r <= edge {
            (1.0 - 0.5 * (r - radius) / (edge - radius))
                * ejecta_modulation(phi, c.azimuth, c.angle).min(1.0)
                * texture.as_ref().map_or(1.0, |t| {
                    1.0 + EJECTA_BRIGHTNESS_TEXTURE * t.noise_floor(phi, r / radius)
                })
        } else {
            0.0
        };
        if r > edge {
            let falloff = (edge / r).powf(1.2);
            for &(p, w, l, s) in &ray_params {
                if r > l {
                    continue;
                }
                let dphi = (phi - p + std::f64::consts::PI).rem_euclid(tau) - std::f64::consts::PI;
                let taper = ((l - r) / (0.3 * l)).clamp(0.0, 1.0);
                let across = dphi * r / w;
                v = v.max((-(across * across)).exp() * taper * falloff * s);
            }
        }
        let value = (v.min(1.0) * freshness) as f32;
        if value > 0.0 {
            cells.push((idx as u32, value));
        }
    }
    cells
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
    let depths: Vec<f32> = solid
        .as_slice()
        .par_iter()
        .map(|h| (level - *h as f64).max(0.0) as f32)
        .collect();
    let mut covered = 0.0;
    let mut total = 0.0;
    for (k, c) in grid.cells().enumerate() {
        let a = grid.cell_area(c);
        total += a;
        if depths[k] > 0.0 {
            covered += a;
        }
    }
    let mut liquid = Field::new(grid);
    liquid.as_mut_slice().copy_from_slice(&depths);
    layers.liquid = liquid;
    (level, covered / total)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::presets::preset;

    fn preset_graph(name: &str, quality: &str) -> ParamGraph {
        let (mut g, _) = preset(name).unwrap().to_graph();
        g.set("detail.quality", quality).unwrap();
        g
    }

    #[test]
    fn the_same_seed_gives_identical_output_on_different_thread_counts() {
        let g = preset_graph("Moon", "preview");
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
        let g = preset_graph("Moon", "preview");
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
        let g = preset_graph("Earth", "preview");
        let mut r = Runner::new();
        let out = r.run(&g, &NoProgress).unwrap();
        assert!(out.finish.sea_level.is_finite());
        assert!(
            out.finish.ocean_fraction > 0.5 && out.finish.ocean_fraction < 0.9,
            "liquid fraction {}",
            out.finish.ocean_fraction
        );
    }

    #[test]
    fn the_moon_is_heavily_cratered_and_earth_is_not() {
        let mut r = Runner::new();
        let moon = r
            .run(&preset_graph("Moon", "preview"), &NoProgress)
            .unwrap();
        let earth = r
            .run(&preset_graph("Earth", "preview"), &NoProgress)
            .unwrap();
        assert!(
            moon.loop_output.craters.len() > 100,
            "loop craters {}",
            moon.loop_output.craters.len()
        );
        assert!(
            moon.finish.replayed_count > 5_000,
            "replayed {}",
            moon.finish.replayed_count
        );
        // Earth retains craters for only about 100 Myr, so it has far fewer per area.
        let moon_area = 4.0 * std::f64::consts::PI * 1.7374e6f64.powi(2);
        let earth_area = 4.0 * std::f64::consts::PI * 6.371e6f64.powi(2);
        let moon_density = moon.finish.replayed_count as f64 / moon_area;
        let earth_density = earth.finish.replayed_count as f64 / earth_area;
        assert!(
            earth_density < 0.05 * moon_density,
            "earth {earth_density:e} moon {moon_density:e}"
        );
        // Rays exist only on young craters, so only a small fraction of the Moon is rayed.
        let rayed = moon
            .finish
            .rays
            .as_slice()
            .iter()
            .filter(|v| **v > 0.05)
            .count() as f64
            / moon.finish.rays.len() as f64;
        assert!(rayed > 0.001 && rayed < 0.25, "rayed fraction {rayed}");
        // Mature highlands are dark; fresh ejecta is brighter.
        let (lo, hi) = {
            let a = moon.finish.albedo.as_slice();
            a.iter()
                .fold((1.0f32, 0.0f32), |(lo, hi), c| (lo.min(c[1]), hi.max(c[1])))
        };
        assert!(lo < 0.2 && hi > 0.25, "albedo range {lo} to {hi}");
    }

    /// Nothing older than a loop crater survives inside it: the finishing pass stamps in
    /// age order and the restore removes what the impact destroyed, so every cell inside
    /// the rim dates from the crater or from something younger.
    #[test]
    fn nothing_older_than_a_loop_crater_survives_inside_it() {
        let mut r = Runner::new();
        let moon = r
            .run(&preset_graph("Moon", "preview"), &NoProgress)
            .unwrap();
        let grid = moon.finish.grid;
        let ages = moon.finish.layers.surface_age.as_slice();
        let mut stamper = Stamper::new(&grid);
        let (mut checked, mut older) = (0usize, 0usize);
        for c in moon.finish.craters.iter().filter(|c| c.from_loop) {
            let inner = 0.45 * c.diameter / c.ellipticity.sqrt();
            stamper.gather(&grid, 1.7374e6, c.center, inner);
            for &(idx, _, _) in &stamper.footprint {
                checked += 1;
                if (ages[idx] as f64) < c.time_gyr - 1e-3 {
                    older += 1;
                }
            }
        }
        assert!(checked > 1_000, "checked {checked} cells");
        assert_eq!(
            older, 0,
            "{older} of {checked} cells older than their crater"
        );
    }

    /// Craters erase each other: a share of the loop craters lose their relief to later
    /// cavities and ejecta, while most basins, which nothing relaxes yet, stay visible.
    #[test]
    fn later_impacts_erase_some_loop_craters_and_most_basins_stay_visible() {
        let mut r = Runner::new();
        let moon = r
            .run(&preset_graph("Moon", "preview"), &NoProgress)
            .unwrap();
        let all = &moon.loop_output.craters;
        let erased = all.iter().filter(|c| !c.visible).count();
        assert!(all.len() > 50, "loop craters {}", all.len());
        let fraction = erased as f64 / all.len() as f64;
        assert!(
            (0.02..0.5).contains(&fraction),
            "erased {erased} of {} loop craters",
            all.len()
        );
        let basins: Vec<_> = all.iter().filter(|c| c.diameter >= 300_000.0).collect();
        let visible = basins.iter().filter(|c| c.visible).count();
        assert!(basins.len() >= 20, "basins {}", basins.len());
        assert!(
            visible * 2 >= basins.len(),
            "visible basins {visible} of {}",
            basins.len()
        );
        // A visible basin keeps kilometres of relief.
        assert!(
            basins
                .iter()
                .filter(|c| c.visible)
                .all(|c| c.relief > 1_000.0)
        );
    }

    /// Mercury takes up its contraction in scarps: many, up to hundreds of kilometres long
    /// and kilometres high, carrying most of the shortening budget; the Moon far less.
    #[test]
    fn mercury_takes_up_its_contraction_in_scarps() {
        let mut r = Runner::new();
        let mercury = r
            .run(&preset_graph("Mercury", "preview"), &NoProgress)
            .unwrap();
        let g = preset_graph("Mercury", "preview");
        let radius = g.f64("body.radius");
        let area = 4.0 * std::f64::consts::PI * radius * radius;
        let budget = 2.0 * g.f64("interior.contraction_total") / radius * area;
        let scarps = &mercury.loop_output.scarps;
        let taken: f64 = scarps
            .iter()
            .map(|s| s.length * tectonics::scarp_shortening(s.height))
            .sum();
        assert!(scarps.len() > 100, "{} scarps", scarps.len());
        let longest = scarps.iter().map(|s| s.length).fold(0.0, f64::max);
        let highest = scarps.iter().map(|s| s.height).fold(0.0, f64::max);
        assert!(
            longest > 300_000.0 && highest > 1_500.0,
            "{longest} m, {highest} m"
        );
        // Stamped scarps carry most of the budget; the rest is in scarps too short for
        // the loop grid.
        assert!(
            taken > 0.3 * budget && taken <= 1.05 * budget,
            "taken {taken} of {budget}"
        );
        let moon = r
            .run(&preset_graph("Moon", "preview"), &NoProgress)
            .unwrap();
        let moon_taken: f64 = moon
            .loop_output
            .scarps
            .iter()
            .map(|s| s.length * tectonics::scarp_shortening(s.height))
            .sum();
        let moon_radius = preset_graph("Moon", "preview").f64("body.radius");
        let per_area = |t: f64, r: f64| t / (4.0 * std::f64::consts::PI * r * r);
        assert!(per_area(moon_taken, moon_radius) < 0.3 * per_area(taken, radius));
    }

    /// The Moon's magma floods its thin-crust lowlands: a sixth of the surface or so,
    /// mostly on the thin-crust side, and the finishing pass keeps that extent.
    #[test]
    fn the_moon_floods_its_thin_crust_basins() {
        let mut r = Runner::new();
        let g = preset_graph("Moon", "preview");
        let out = r.run(&g, &NoProgress).unwrap();
        let layers = &out.loop_output.keyframes.last().unwrap().layers;
        let fraction = volcanism::lava_fraction(layers);
        assert!(
            fraction > 0.05 && fraction < 0.35,
            "lava fraction {fraction}"
        );
        assert!(out.loop_output.floods.len() > 20);
        let crust = layers.crust_thickness.as_slice();
        let units = layers.unit_id.as_slice();
        let mut sorted: Vec<f32> = crust.to_vec();
        sorted.sort_by(f32::total_cmp);
        let median = sorted[sorted.len() / 2];
        let (mut thin, mut thick) = (0.0, 0.0);
        for k in 0..units.len() {
            if volcanism::is_lava_unit(units[k]) {
                if crust[k] < median {
                    thin += 1.0;
                } else {
                    thick += 1.0;
                }
            }
        }
        assert!(thin > 1.5 * thick, "thin {thin} thick {thick}");
        assert!(median > 1e4 && median < 1e5, "crust median {median} m");
        let fine = volcanism::lava_fraction(&out.finish.layers);
        assert!((fine - fraction).abs() < 0.05, "{fine} vs {fraction}");
        // The maria carry wrinkle ridges, all of them on lava.
        let ridges = &out.loop_output.ridges;
        assert!(ridges.len() > 20, "{} ridges", ridges.len());
        let g = out.loop_output.grid;
        for r in ridges {
            let idx = g.index(g.locate(r.midpoint()).0);
            assert!(volcanism::is_lava_unit(units[idx]));
        }
    }

    /// Mars raises one plume: a rise kilometres high with shields on it, a canyon
    /// kilometres deep on its flank, graben, and lava concentrated around it.
    #[test]
    fn mars_raises_a_rise_with_shields_and_a_canyon() {
        let mut r = Runner::new();
        let g = preset_graph("Mars", "preview");
        let out = r.run(&g, &NoProgress).unwrap();
        let lo = &out.loop_output;
        assert_eq!(lo.rises.len(), 1);
        assert!(lo.shields.len() >= 4, "{} shields", lo.shields.len());
        assert_eq!(lo.canyons.len(), 1);
        assert!(lo.graben.len() > 20, "{} graben", lo.graben.len());
        let grid = out.finish.grid;
        let radius = g.f64("body.radius");
        let height = out.finish.height.as_slice();
        let at = |d: [f64; 3]| height[grid.index(grid.locate(d).0)] as f64;
        let rise = &lo.rises[0];
        let (east, north) = crate::geom::local_basis(rise.center);
        // The rise stands above the ground beyond it.
        let far: f64 = (0..8)
            .map(|k| {
                let phi = k as f64 * std::f64::consts::FRAC_PI_4;
                at(crate::geom::offset(
                    rise.center,
                    east,
                    north,
                    1.2 * rise.radius / radius,
                    phi,
                ))
            })
            .sum::<f64>()
            / 8.0;
        let top = at(rise.center);
        assert!(top - far > 0.3 * rise.height, "rise {top} vs far {far}");
        // The tallest shield reaches most of its height above the dome around it.
        let big = lo
            .shields
            .iter()
            .max_by(|a, b| a.height.total_cmp(&b.height))
            .unwrap();
        let (e2, n2) = crate::geom::local_basis(big.center);
        let foot: f64 = (0..8)
            .map(|k| {
                let phi = k as f64 * std::f64::consts::FRAC_PI_4;
                at(crate::geom::offset(
                    big.center,
                    e2,
                    n2,
                    1.1 * big.radius / radius,
                    phi,
                ))
            })
            .sum::<f64>()
            / 8.0;
        let summit = at(crate::geom::offset(
            big.center,
            e2,
            n2,
            0.15 * big.radius / radius,
            0.0,
        ));
        assert!(
            summit - foot > 0.5 * big.height,
            "shield {summit} vs foot {foot}"
        );
        // The canyon floor lies kilometres below its rims.
        let canyon = &lo.canyons[0];
        let mid = canyon.midpoint();
        let next = canyon.points[canyon.points.len() / 2 + 1];
        let t = crate::geom::normalize([next[0] - mid[0], next[1] - mid[1], next[2] - mid[2]]);
        let n = crate::geom::normalize(crate::geom::cross(mid, t));
        let side = |s: f64| {
            let d = canyon.width * s / radius;
            at(crate::geom::normalize([
                mid[0] + d * n[0],
                mid[1] + d * n[1],
                mid[2] + d * n[2],
            ]))
        };
        let rims = 0.5 * (side(1.0) + side(-1.0));
        assert!(
            rims - at(mid) > 3_000.0,
            "canyon floor {} rims {rims}",
            at(mid)
        );
        // The Moon raises none.
        let moon = r
            .run(&preset_graph("Moon", "preview"), &NoProgress)
            .unwrap();
        assert!(moon.loop_output.rises.is_empty() && moon.loop_output.shields.is_empty());
    }

    /// Venus renews its surface: most of it is plains the age of the event, a few per
    /// cent is older tessera, coronae and domes dot the plains, and the craters are few.
    #[test]
    fn venus_renews_its_surface_and_keeps_tessera() {
        let mut r = Runner::new();
        let g = preset_graph("Venus", "preview");
        let out = r.run(&g, &NoProgress).unwrap();
        let lo = &out.loop_output;
        assert!(!lo.rises.is_empty() && !lo.shields.is_empty());
        let event = lo.resurfacing.as_ref().expect("an event").time_gyr;
        assert!((event - g.f64("tectonics.resurfacing_time")).abs() < 1e-9);
        assert!(lo.coronae.len() > 50 && !lo.domes.is_empty());
        let tessera = episodic::tessera_fraction(&out.finish.layers);
        assert!(tessera > 0.04 && tessera < 0.14, "tessera {tessera}");
        // Plains carry the event's age; tessera is older.
        let layers = &out.finish.layers;
        let units = layers.unit_id.as_slice();
        let ages = layers.surface_age.as_slice();
        let (mut plains, mut plains_at_event) = (0usize, 0usize);
        let mut tessera_older = 0usize;
        for k in 0..units.len() {
            if units[k] == crate::units::UNIT_LAVA {
                plains += 1;
                if (ages[k] as f64 - event).abs() < 0.3 {
                    plains_at_event += 1;
                }
            } else if units[k] == crate::units::UNIT_TESSERA && (ages[k] as f64) < event - 0.1 {
                tessera_older += 1;
            }
        }
        assert!(
            plains_at_event * 10 > plains * 7,
            "{plains_at_event} of {plains}"
        );
        assert!(tessera_older > 0);
        // Few craters, none older than the event.
        assert!(
            out.finish.craters.len() < 2_000,
            "{} craters",
            out.finish.craters.len()
        );
        assert!(
            out.finish
                .craters
                .iter()
                .all(|c| c.time_gyr >= event - 1e-6)
        );
    }

    #[test]
    fn moon_height_spectrum_has_a_terrain_like_slope() {
        let mut r = Runner::new();
        let moon = r
            .run(&preset_graph("Moon", "preview"), &NoProgress)
            .unwrap();
        let width = 1024usize;
        let row = crate::grid::projection::to_equirectangular(
            &moon.finish.grid,
            moon.finish.height.as_slice(),
            width,
            width / 2,
        );
        let equator: Vec<f64> = row[(width / 4) * width..(width / 4 + 1) * width]
            .iter()
            .map(|v| *v as f64)
            .collect();
        let slope = spectral_slope(&equator, 4, 128);
        assert!(slope < -1.3 && slope > -3.8, "spectral slope {slope}");
    }

    /// Log-log slope of the periodogram of a periodic row between two wavenumbers.
    fn spectral_slope(row: &[f64], k_lo: usize, k_hi: usize) -> f64 {
        let n = row.len();
        let mean = row.iter().sum::<f64>() / n as f64;
        let mut xs = Vec::new();
        let mut ys = Vec::new();
        let mut k = k_lo;
        while k <= k_hi {
            let (mut re, mut im) = (0.0, 0.0);
            for (i, v) in row.iter().enumerate() {
                let ang = 2.0 * std::f64::consts::PI * k as f64 * i as f64 / n as f64;
                re += (v - mean) * ang.cos();
                im -= (v - mean) * ang.sin();
            }
            let power = (re * re + im * im).max(1e-30);
            xs.push((k as f64).ln());
            ys.push(power.ln());
            k = (k as f64 * 1.25).ceil() as usize;
        }
        let mx = xs.iter().sum::<f64>() / xs.len() as f64;
        let my = ys.iter().sum::<f64>() / ys.len() as f64;
        let num: f64 = xs.iter().zip(&ys).map(|(x, y)| (x - mx) * (y - my)).sum();
        let den: f64 = xs.iter().map(|x| (x - mx).powi(2)).sum();
        num / den
    }
}
