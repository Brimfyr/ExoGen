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

pub mod cryo;
pub mod deform;
pub mod episodic;
pub mod fabric;
pub mod ice;
pub mod ice_shell;
pub mod impacts;
pub mod io;
pub mod lineae;
pub mod materials;
pub mod plume;
pub mod primordial_crust;
pub mod relaxation;
pub mod ridges;
pub mod rings;
pub mod roughness;
pub mod tectonics;
pub mod tidal;
pub mod volcanism;

use std::collections::HashMap;
use std::sync::Arc;
use std::time::{Duration, Instant};

use rayon::prelude::*;

use crate::constants::*;
use crate::crater::production::poisson;
use crate::crater::{
    CraterRecord, CraterTexture, Curtain, EdgeShape, PopulationConfig, ProfileTable, Rays,
    RestoreParams, StampOutput, StampParams, Stamper, dimensions, ejecta_coverage,
    ejecta_modulation, restore_compute, sample_population, stamp_compute, with_stamper,
};
use crate::fields::{Field, LayerSet};
use crate::geom::{angle_between, azimuth_from, dot, local_basis, offset};
use crate::grid::CubeSphere;
use crate::params::ParamGraph;
use crate::params::derive::physics::volatile;
use crate::rng::{ACCEPT_STREAM, Rng, hash_words};
use crate::thermal::ThermalHistory;
use crate::units::UNIT_SECONDARY_FIELD;

/// Loop steps by quality: preview, standard, high.
pub const LOOP_STEPS: [usize; 3] = [30, 90, 90];
/// The loop's fixed time bins. Every sampler draws per bin rather than per step, with the
/// bin as its random index, so the population of craters, scarps, eruptions and the Io
/// surface is the same at every quality and the preview predicts the export. The count is
/// a multiple of every step count in `LOOP_STEPS`, so a bin never straddles a step.
pub const LOOP_TIME_BINS: usize = 900;

/// The bins a step covers, as (index, start, end) in order.
pub fn time_bins(age: f64, t0: f64, t1: f64) -> Vec<(u64, f64, f64)> {
    let w = age / LOOP_TIME_BINS as f64;
    let b0 = (t0 / w).round() as usize;
    let b1 = (t1 / w).round() as usize;
    (b0..b1)
        .map(|b| (b as u64, b as f64 * w, (b + 1) as f64 * w))
        .collect()
}
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
    /// Lava floods emplaced inside the loop, in order of formation, and the magmatic
    /// provinces their vents favoured.
    pub floods: Vec<volcanism::FloodRecord>,
    pub provinces: Vec<volcanism::Province>,
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
    /// The Io regime's paterae and mountains through the body's life, each with its
    /// birth and death; the finish takes those alive at the epoch.
    pub paterae: Vec<io::Patera>,
    pub mountains: Vec<io::Mountain>,
    /// The Io regime's buried plains, the bedrock every step lays its features on, on
    /// the coarse grid.
    pub io_base: Option<Field<f32>>,
    /// An ice shell's lanes of grooved terrain and when they formed.
    pub lanes: Vec<ice_shell::Lane>,
    /// The share of the surface the lanes cover once all have formed.
    pub lane_coverage: f64,
    /// How much of each cell's relief the lanes renewed at their times, 0 to 1, on the
    /// coarse grid: the finish completes the renewal out to the lanes' resolved edges.
    pub lane_core: Option<Field<f32>>,
    /// The past episodes of focused heating (F53), and when the present one began.
    pub episodes: Vec<ice::FocusEpisode>,
    pub polar_start: f64,
    /// The share of the surface the lanes and the polar terrains, present and past,
    /// renewed between them, by area.
    pub renewed_fraction: f64,
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
    /// Floods present at the epoch, and the magmatic provinces.
    pub floods: Vec<volcanism::FloodRecord>,
    pub provinces: Vec<volcanism::Province>,
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
    /// The Io regime's features, all of the present surface, and the plumes of its
    /// paterae over the history the surface remembers (F84).
    pub paterae: Vec<io::Patera>,
    pub mountains: Vec<io::Mountain>,
    pub io_flows: Vec<io::FlowField>,
    pub plume_episodes: Vec<io::PlumeEpisode>,
    /// The small eruptions laid one by one over that history (F89).
    pub small_eruption_count: usize,
    /// An ice shell's cracks present at the epoch, its chaos, and the domes a shell over
    /// rock has built, with the share of the surface that is chaos matrix.
    pub lineae: Vec<lineae::Linea>,
    pub chaos: Vec<cryo::ChaosGeneration>,
    pub cryodomes: Vec<cryo::Cryodome>,
    pub chaos_fraction: f64,
    /// The lanes of grooved terrain present at the epoch, and the shares of the surface
    /// that are grooved terrain and renewed polar terrain.
    pub lanes: Vec<ice_shell::Lane>,
    /// The share of the surface the lanes cover, and the shares that still show as
    /// grooved terrain and as polar terrain, craters since having taken the rest.
    pub lane_coverage: f64,
    pub grooved_fraction: f64,
    pub polar_fraction: f64,
    /// The past episodes of focused heating whose terrains the surface carries (F53).
    pub episodes: Vec<ice::FocusEpisode>,
    /// The share of the surface the lanes and the polar terrains, present and past,
    /// renewed between them, by area, whatever craters they have taken since.
    pub renewed_fraction: f64,
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
            io::VERSION,
            ice_shell::VERSION,
            cryo::VERSION,
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
pub const FINISH_VERSION: u64 = 80;

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
    let mut plume = plume::PlumeProcess::new(graph, grid);
    let mut floods = volcanism::FloodProcess::new(
        graph,
        grid,
        &thermal,
        plume.as_ref().map_or(&[], |p| p.rises.as_slice()),
    );
    let mut episodic = episodic::EpisodicProcess::new(graph, grid);
    let mut io = io::IoProcess::new(graph, grid);
    let ice_shell = ice_shell_history(graph, &thermal);
    let mut ice_shell_process = ice_shell::IceShellProcess::new(graph, grid);
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
        // An ice shell relaxes what it carries, the longest wavelengths first, with the
        // shell as the thermal history has it now.
        if let Some(shell_at) = &ice_shell {
            let local = ice_shell_process
                .as_mut()
                .map_or(relaxation::Local::default(), |p| p.relaxation(t0, t1));
            let floor = relaxation::relax(
                &mut layers.bedrock,
                graph.f64("body.radius"),
                dt,
                &shell_at(t1),
                &local,
            );
            if let (Some(p), Some(f)) = (&mut ice_shell_process, floor) {
                p.set_floor(f);
            }
        }
        // The shell's partial resurfacing: the lanes of grooved terrain at their time.
        if let Some(p) = &mut ice_shell_process {
            p.before(&mut layers, t0, t1);
        }
        // The slow processes first: the dome lifts, the shields grow, the canyon deepens;
        // an episodic lid's event renews the plains.
        if let Some(p) = &mut plume {
            p.step(&mut layers, t0, t1);
        }
        if let Some(e) = &mut episodic
            && let Some(fill) = e.step(
                &mut layers,
                t0,
                t1,
                plume.as_ref().map(|p| p.load.as_slice()),
            )
            && let Some(f) = &mut floods
        {
            f.note_resurfacing(&fill, e.event_time());
        }
        if let Some(i) = &mut io {
            i.step(&mut layers, t0, t1, s % keep_every == 0 || s == steps);
        }
        let boost = plume.as_ref().map(|p| p.vent_boost());
        // Everything that happens in the step, in order of time: craters, scarps and
        // eruptions interleaved as they come, so what is older is cut by what is younger.
        let craters = impacts.sample(t0, t1);
        let scarps = contraction
            .as_mut()
            .map_or_else(Vec::new, |c| c.sample(t0, t1));
        let eruptions = floods
            .as_mut()
            .map_or_else(Vec::new, |f| f.sample(&layers, boost.as_deref(), t0, t1));
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
        // The polar terrain renews itself after everything the step brought it.
        if let Some(p) = &mut ice_shell_process {
            p.after(&mut layers);
            if std::env::var_os("EXOGEN_DEBUG").is_some() && s % 10 == 0 {
                let b = layers.bedrock.as_slice();
                let (lo, hi) = b
                    .iter()
                    .fold((f32::MAX, f32::MIN), |(lo, hi), &v| (lo.min(v), hi.max(v)));
                eprintln!("step {s}: bedrock {lo:.0} to {hi:.0} m");
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
    let (paterae, mountains, io_base) = match io {
        Some(i) => (i.paterae, i.mountains, i.base.map(|b| b.bedrock)),
        None => (Vec::new(), Vec::new(), None),
    };
    let (resurfacing, coronae, domes, novae) = match episodic {
        Some(e) => (e.resurfacing, e.coronae, e.domes, e.novae),
        None => (None, Vec::new(), Vec::new(), Vec::new()),
    };
    let (rises, shields, canyons, mut graben, load) = match plume {
        Some(p) => (p.rises, p.shields, p.canyons, p.graben, Some(p.load)),
        None => (Vec::new(), Vec::new(), Vec::new(), Vec::new(), None),
    };
    // The novae's fractures are stamped as the graben are.
    graben.extend(novae);
    let (floods, ridges, provinces) = match floods {
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
            (f.records, ridges, f.provinces)
        }
        None => (Vec::new(), Vec::new(), Vec::new()),
    };
    let (lanes, lane_coverage, lane_core, episodes, polar_start, renewed_fraction) =
        match ice_shell_process {
            Some(p) => (
                p.lanes,
                p.lane_coverage,
                p.lane_core.map(|v| ice_shell::field_from(grid, &v)),
                p.episodes,
                p.present_start,
                p.renewed_fraction,
            ),
            None => (Vec::new(), 0.0, None, Vec::new(), 0.0, 0.0),
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
        provinces,
        ridges,
        rises,
        shields,
        canyons,
        graben,
        resurfacing,
        coronae,
        domes,
        paterae,
        mountains,
        io_base,
        lanes,
        lane_coverage,
        lane_core,
        episodes,
        polar_start,
        renewed_fraction,
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
    let finish_clock = Instant::now();

    sink.report(Progress {
        stage: "finishing",
        fraction: 0.05,
        message: format!("resampling to face size {face}"),
    });
    let mut layers = kf.layers.resample(grid);
    // The coarse bedrock resampled bilinearly carries the facets of the loop's cells on
    // every slope the loop relaxed or built; smoothed over the loop's cell it carries the
    // same relief without them, and the passes lay their resolved detail on it after.
    {
        let coarse_cell_m = kf.layers.grid.delta() * body_radius;
        let smoothed = relaxation::low_pass(&layers.bedrock, cell, coarse_cell_m);
        layers.bedrock.as_mut_slice().copy_from_slice(&smoothed);
    }
    if sink.cancelled() {
        return Err(RunError::Cancelled);
    }
    // The Io regime lays its surface afresh at every step, and the finish is a step:
    // the loop's coarse stamps of its features are cleared, so that each is drawn once
    // at this resolution and no coarse halo of its unit or coverage outlives it.
    if io::is_io_regime(graph) {
        io::IoProcess::clear_surface(&mut layers, present);
    }
    let shell_map = ice::ShellMap::new(graph);
    let shell_history = ice_shell_history(graph, &loop_output.thermal);
    let cracks = lineae::Cracks::new(graph).map(|c| c.with_freezing(&loop_output.thermal));
    let lineae: Vec<lineae::Linea> = cracks
        .as_ref()
        .map(|c| c.generate(present, &grid, inputs.seed))
        .unwrap_or_default();
    // Chaos where a thin shell melts through, and the domes a shell over rock extrudes,
    // each at the height its relaxation has left it.
    let chaos = cryo::chaos_population(graph, present, inputs.seed);
    // The ridged plains' fabric renews what the drawn lineae and the chaos matrix do not
    // (F81).
    let fabric: Vec<fabric::FabricSlice> = match &cracks {
        Some(_) => fabric::fabric_slices(
            graph,
            present,
            fabric::drawn_coverage(&lineae, body_radius),
            chaos.iter().map(|c| c.matrix_share()).sum(),
        ),
        None => Vec::new(),
    };
    let cryodomes = cryo::cryodome_population(graph, present, inputs.seed);
    let dome_heights: Vec<f64> = cryodomes
        .iter()
        .map(|d| {
            d.height
                * match &shell_history {
                    Some(shell_at) => relaxation::survival_since(
                        2.0 * d.radius,
                        d.time_gyr,
                        present,
                        &**shell_at,
                        1.0,
                    ),
                    None => 1.0,
                }
        })
        .collect();
    // The lanes of grooved terrain formed by the epoch, and the polar terrain's renewal,
    // which stands at the end of the sequence with its folds along the cracks.
    let lanes: Vec<ice_shell::Lane> = loop_output
        .lanes
        .iter()
        .filter(|l| l.time_gyr <= present + 1e-9)
        .cloned()
        .collect();
    let polar = shell_map.has_active_terrain() && present + 1e-9 >= loop_output.polar_start;
    let polar_time = (present - POLAR_RENEWAL_GYR).max(0.0);
    // The past episodes of focused heating (F53) that ended by this epoch: each terrain's
    // shell map centred where the body left it.
    let episodes: Vec<ice::FocusEpisode> = loop_output
        .episodes
        .iter()
        .filter(|e| e.end <= present + 1e-9)
        .cloned()
        .collect();
    let episode_maps: Vec<ice::ShellMap> = episodes
        .iter()
        .map(|e| shell_map.with_pole(e.center))
        .collect();
    // The renewed terrains on the resolved grid (F60, F61). The loop renewed their relief
    // on its own grid, with the ramp across each edge and the halo beyond it, and every
    // crater since stands on that; here their outlines are drawn at this resolution and
    // composed in the order they formed, so that the crater population that follows
    // reads the ages they leave, the colour follows the true edge, and a younger terrain
    // lies over an older one.
    let halo = POLAR_HALO_SHELLS * graph.f64("interior.ice_shell_thickness");
    // The loose cover the plumes' fallout has laid on each cell's surface since it
    // formed (F64), which an open fracture drains into pits: on the ground as the loop
    // left it, for the lanes' damage zones, and again once the terrains hold their
    // cells.
    let loose = cryo::LooseCover::new(graph, &loop_output.episodes, present);
    let cover_on = |layers: &LayerSet| -> Option<Vec<f32>> {
        loose.as_ref().map(|l| {
            let ages = layers.surface_age.as_slice();
            (0..grid.cell_count())
                .into_par_iter()
                .map(|k| {
                    l.thickness_between(grid.center(grid.cell(k)), ages[k] as f64, present) as f32
                })
                .collect()
        })
    };
    let loose_before = cover_on(&layers);
    let mut renewal = ice_shell::Renewal::new(grid.cell_count());
    let fossil_edges: Vec<Vec<f32>> = episode_maps
        .iter()
        .map(|map| ice_shell::polar_edge(&grid, body_radius, map, &lineae, halo))
        .collect();
    let polar_edge: Vec<f32> = if polar && shell_map.pole.is_some() {
        ice_shell::polar_edge(&grid, body_radius, &shell_map, &lineae, halo)
    } else {
        Vec::new()
    };
    let polar = !polar_edge.is_empty();
    {
        let lane_hits =
            ice_shell::lane_hits(&grid, body_radius, cell, &lanes, loose_before.as_deref());
        // Each terrain with its time: a lane (0), a past episode's terrain (1), the
        // present polar terrain (2).
        let mut order: Vec<(f64, u8, usize)> = (0..lanes.len())
            .map(|i| (lanes[i].time_gyr, 0u8, i))
            .collect();
        order.extend((0..episodes.len()).map(|k| (episodes[k].end, 1u8, k)));
        if polar {
            order.push((polar_time, 2u8, 0));
        }
        order.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)).then(a.2.cmp(&b.2)));
        for (time, kind, i) in order {
            match kind {
                0 => renewal.take_lane(&lanes[i], &lane_hits[i]),
                1 => renewal.take_polar(&fossil_edges[i], time),
                _ => renewal.take_polar(&polar_edge, time),
            }
        }
    }
    // The lanes' floor where the loop's cells did not reach their edges: what the loop
    // renewed of each cell, resampled and smoothed over its cell so that no facet of
    // the coarse grid shows in a lane, and the rest taken here.
    if !lanes.is_empty() {
        let core: Option<Vec<f32>> = loop_output.lane_core.as_ref().map(|c| {
            let fine = c.resample(grid);
            relaxation::low_pass(&fine, cell, loop_output.grid.delta() * body_radius)
        });
        ice_shell::lay_lane_floor(&mut layers, body_radius, core.as_deref(), &renewal);
    }
    renewal.apply(&mut layers);
    let renewed: &[f32] = &renewal.time;
    drop(loose_before);
    let loose_cover = cover_on(&layers);

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
    // The Io regime's features alive at the epoch.
    let paterae: Vec<io::Patera> = loop_output
        .paterae
        .iter()
        .filter(|p| p.alive(present))
        .cloned()
        .collect();
    let mountains: Vec<io::Mountain> = loop_output
        .mountains
        .iter()
        .filter(|m| m.alive(present))
        .cloned()
        .collect();
    // The flow fields of the paterae alive at the epoch, traced over the finish's ground
    // (F85).
    let io_flows = if io::is_io_regime(graph) {
        io::flow_fields(graph, &paterae, &mountains, &layers.bedrock)
    } else {
        Vec::new()
    };
    // The Io regime's mountains and paterae are drawn whole at this resolution on the
    // buried plains the loop laid its coarse versions on (F87, F88), smoothed over the
    // loop's cell as the rest of the bedrock is, so nothing of the coarse ones is left.
    if let Some(io_base) = &loop_output.io_base {
        let coarse_cell_m = loop_output.grid.delta() * body_radius;
        let smoothed = relaxation::low_pass(&io_base.resample(grid), cell, coarse_cell_m);
        layers.bedrock.as_mut_slice().copy_from_slice(&smoothed);
    }
    // The small eruptions between the paterae over the same history, laid one by one where
    // the grid resolves them (F89).
    let small_eruptions = if io::is_io_regime(graph) {
        io::small_eruptions(graph.f64("seed.volcanism") as u64, cell)
    } else {
        Vec::new()
    };
    // The plumes the paterae alive at the epoch have thrown over the history the surface
    // remembers (F84).
    let plume_episodes = if io::is_io_regime(graph) {
        io::plume_history(
            &paterae,
            graph.f64("body.gravity"),
            graph.f64("seed.volcanism") as u64,
        )
    } else {
        Vec::new()
    };
    // An ice shell's cracks within the surface's retention age, drawn per time bin.
    let base: Vec<f32> = layers.bedrock.as_slice().to_vec();
    let base_deep: Vec<f32> = layers.deep.as_slice().to_vec();
    let mut smear = volcanism::LavaSmear::new(
        loop_output.grid,
        grid,
        &floods,
        resurfacing.as_ref().map(|r| r.fill.as_slice()),
    );
    // What the slow processes added after each flood, so a refill sees the ground the
    // flood found and not the dome or shield that grew under it since.
    let deform = deform::Deformation::new(
        loop_output.grid,
        body_radius,
        loop_output.age_gyr / loop_output.steps.max(1) as f64,
        present,
        &loop_output.rises,
        &shields,
        &canyons,
        &coronae,
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
    if std::env::var_os("EXOGEN_DEBUG").is_some() {
        eprintln!(
            "finish phases: resampling and terrains {:.1} s",
            finish_clock.elapsed().as_secs_f64()
        );
    }
    let phase_clock = Instant::now();
    let population = sample_population(&grid, &layers.surface_age, present, &cfg);
    let secondaries = if std::env::var_os("EXOGEN_SKIP_SECONDARIES").is_some() {
        Vec::new()
    } else {
        sample_secondaries(&loop_craters, &inputs, d_min, present)
    };
    let replay: Vec<(crate::crater::CraterSample, bool)> = population
        .into_iter()
        .map(|c| (c, false))
        .chain(secondaries.into_iter().map(|c| (c, true)))
        .collect();
    // On an ice shell a crater relaxes over its own life by wavelength (F70): each
    // wavelength keeps the share the shell's viscosity lets it keep, and the share each
    // past episode's thinned shell let it keep while that lasted (F66).
    let episode_shells: Vec<relaxation::Shell> = match &shell_history {
        Some(shell_at) => episodes
            .iter()
            .map(|e| shell_at(0.5 * (e.start + e.end)))
            .collect(),
        None => Vec::new(),
    };
    let crater_kept = |wavelength: f64, time: f64, center: [f64; 3]| -> f64 {
        let Some(shell_at) = &shell_history else {
            return 1.0;
        };
        let mut kept = relaxation::survival_since(
            wavelength,
            time,
            present,
            &**shell_at,
            shell_map.factor(center),
        );
        for ((e, map), shell) in episodes.iter().zip(&episode_maps).zip(&episode_shells) {
            if e.end > time {
                kept *= relaxation::thinned_survival(
                    wavelength,
                    e.end - e.start.max(time),
                    shell,
                    map.factor(center),
                );
            }
        }
        kept
    };
    // A crater relaxes by wavelength when the longest band its table carries loses any
    // of its relief.
    let band_relaxation = std::env::var_os("EXOGEN_SKIP_BAND_RELAX").is_none();
    let relaxes = |radius: f64, time: f64, center: [f64; 3]| -> bool {
        band_relaxation
            && shell_history.is_some()
            && crater_kept(
                RELAX_BAND_WAVELENGTH * RELAX_LONGEST_BLUR * radius,
                time,
                center,
            ) < 0.999
    };
    // What a replayed crater keeps at its diameter: its whole relief's share where it
    // does not relax by band, and the debug summary's measure.
    let relaxed: Vec<f32> = replay
        .par_iter()
        .map(|(c, _)| crater_kept(c.diameter, c.time_gyr, c.center) as f32)
        .collect();
    // The ring faults of the basins whose collapse drew the substrate in and broke the
    // lithosphere (F69): every crater, of the loop or replayed, with the lid and the
    // substrate as they were at its time and place.
    let ring_systems: Vec<rings::RingSystem> = match &shell_history {
        Some(_) if std::env::var_os("EXOGEN_SKIP_RINGS").is_some() => Vec::new(),
        Some(shell_at) => {
            let mut basins: Vec<([f64; 3], f64, f64, u64)> = loop_craters
                .iter()
                .map(|c| (c.center, c.diameter, c.time_gyr, c.id))
                .collect();
            basins.extend(
                replay
                    .iter()
                    .filter(|(_, secondary)| !*secondary)
                    .map(|(c, _)| (c.center, c.diameter, c.time_gyr, c.key)),
            );
            basins
                .par_iter()
                .filter_map(|&(center, diameter, time, seed)| {
                    let factor = shell_map.factor(center);
                    let shell = shell_at(time);
                    let local = relaxation::Shell {
                        thickness: shell.thickness * factor,
                        elastic: shell.elastic * factor,
                        ..shell
                    };
                    let transient =
                        crate::crater::morphology::transient_diameter(diameter, &inputs.scaling);
                    let lithosphere = rings::collapse_lithosphere(&local);
                    if rings::melosh_gamma(transient, lithosphere) < RING_FAILURE {
                        return None;
                    }
                    let inflow = 1.0 - relaxation::survival(diameter, RING_INFLOW_GYR, &local);
                    let kept = |l: f64| crater_kept(l, time, center);
                    rings::ring_system(
                        &rings::Basin {
                            center,
                            diameter,
                            time_gyr: time,
                            seed,
                            lithosphere,
                            inflow,
                            kept: &kept,
                        },
                        &inputs.scaling,
                        body_radius,
                        cell,
                    )
                })
                .collect()
        }
        None => Vec::new(),
    };
    if std::env::var_os("EXOGEN_DEBUG").is_some() {
        for r in &ring_systems {
            eprintln!(
                "ring system at {:.2} Gyr: transient {:.0} km on a {:.1} km lithosphere, gamma {:.1}, graben from {:.0} to {:.0} km, {} rings in {} arcs",
                r.time_gyr,
                r.transient / 1e3,
                r.lithosphere / 1e3,
                r.gamma,
                r.inner / 1e3,
                r.reach / 1e3,
                r.rings,
                r.faults.len()
            );
        }
    }
    if std::env::var_os("EXOGEN_DEBUG").is_some() && shell_history.is_some() {
        let median = |lo: f64, hi: f64| -> (f32, usize) {
            let mut v: Vec<f32> = replay
                .iter()
                .zip(&relaxed)
                .filter(|((c, _), _)| c.diameter >= lo && c.diameter < hi)
                .map(|(_, r)| *r)
                .collect();
            v.sort_by(f32::total_cmp);
            (v.get(v.len() / 2).copied().unwrap_or(f32::NAN), v.len())
        };
        let (s, ns) = median(30e3, 70e3);
        let (l, nl) = median(150e3, 2_000e3);
        eprintln!(
            "replay survival: median {s:.3} for {ns} craters of 30 to 70 km, {l:.3} for {nl} over 150 km"
        );
    }

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
        Ring(usize),
        Resurfacing,
        Corona(usize),
        Dome(usize),
        Patera(usize),
        Mountain(usize),
        IoFlow(usize),
        Linea(usize),
        Groove(usize),
        PolarTerrain,
        FossilTerrain(usize),
        Chaos(usize),
        Fabric(usize),
        Cryodome(usize),
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
            + ring_systems.len()
            + lineae.len()
            + lanes.len()
            + 1
            + chaos.len()
            + fabric.len()
            + cryodomes.len()
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
    // A basin's ring faults break as it collapses, after its cavity and before anything
    // younger; its arcs are drawn together, in parallel.
    for (k, r) in ring_systems.iter().enumerate() {
        events.push((r.time_gyr, 2, r.id ^ 0x2196_0000_0000, Event::Ring(k)));
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
    for (k, p) in paterae.iter().enumerate() {
        events.push((present, 0, p.id, Event::Patera(k)));
    }
    for (k, m) in mountains.iter().enumerate() {
        events.push((present, 0, m.id, Event::Mountain(k)));
    }
    for (k, f) in io_flows.iter().enumerate() {
        events.push((present, 0, f.id, Event::IoFlow(k)));
    }
    for (k, l) in lineae.iter().enumerate() {
        events.push((l.time_gyr, 0, l.id ^ 0x11CE_0000_0000, Event::Linea(k)));
    }
    for (k, l) in lanes.iter().enumerate() {
        events.push((l.time_gyr, 0, l.id ^ 0x6700_0000_0000, Event::Groove(k)));
    }
    if polar {
        events.push((polar_time, 0, 0x5027_0000_0000, Event::PolarTerrain));
    }
    for (k, e) in episodes.iter().enumerate() {
        events.push((
            e.end,
            0,
            0x5027_F000_0000 + k as u64,
            Event::FossilTerrain(k),
        ));
    }
    for (k, c) in chaos.iter().enumerate() {
        events.push((c.time_gyr, 0, c.id ^ 0xC4A0_0000_0000, Event::Chaos(k)));
    }
    for (k, f) in fabric.iter().enumerate() {
        events.push((f.time_gyr(), 0, f.id, Event::Fabric(k)));
    }
    for (k, d) in cryodomes.iter().enumerate() {
        events.push((d.time_gyr, 0, d.id ^ 0xD0CE_0000_0000, Event::Cryodome(k)));
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
    if std::env::var_os("EXOGEN_DEBUG").is_some() {
        eprintln!(
            "finish phases: population, relaxation and rings {:.1} s",
            phase_clock.elapsed().as_secs_f64()
        );
    }
    let stamp_clock = Instant::now();
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
                Event::Ring(k) => {
                    let (c, r) = ring_systems[*k].cap(body_radius);
                    (c, r + 2.0 * grid.delta())
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
                Event::Patera(k) => {
                    let p = &paterae[*k];
                    let reach = p.extent() + 3.0 * io::patera_wall(p.depth, cell);
                    (p.center, reach / body_radius + 2.0 * grid.delta())
                }
                Event::Mountain(k) => {
                    let m = &mountains[*k];
                    (m.center, m.reach(cell) / body_radius + 2.0 * grid.delta())
                }
                Event::IoFlow(k) => {
                    let f = &io_flows[*k];
                    (f.center, f.cap_radius(body_radius) + 2.0 * grid.delta())
                }
                Event::Linea(k) => {
                    let l = &lineae[*k];
                    (
                        l.midpoint(),
                        l.cap_radius(body_radius, cell) + 2.0 * grid.delta(),
                    )
                }
                Event::Groove(k) => {
                    let l = &lanes[*k];
                    (
                        l.center(),
                        l.cap_radius(body_radius, cell) + 2.0 * grid.delta(),
                    )
                }
                Event::PolarTerrain => (
                    shell_map.pole.unwrap_or([0.0, 0.0, 1.0]),
                    shell_map.radius_rad + 2.0 * grid.delta(),
                ),
                Event::FossilTerrain(k) => (
                    episodes[*k].center,
                    shell_map.radius_rad + 2.0 * grid.delta(),
                ),
                // A generation of chaos reads the whole shell's warmth, and a slice of
                // the fabric crosses the whole body.
                Event::Chaos(_) | Event::Fabric(_) => {
                    ([0.0, 0.0, 1.0], std::f64::consts::PI + 2.0 * grid.delta())
                }
                Event::Cryodome(k) => {
                    let d = &cryodomes[*k];
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
                        if std::env::var_os("EXOGEN_SKIP_LOOP_RESTORE").is_some() {
                            return None;
                        }
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
                        // On an ice shell the detail relaxes by wavelength over the
                        // crater's life (F70), and its terraces, rings and peak keep the
                        // share of a short wavelength; the loop relaxed the rest.
                        let mut feature_weight = 1.0;
                        if relaxes(dim.radius, c.time_gyr, c.center) {
                            let kept = |l: f64| crater_kept(l, c.time_gyr, c.center);
                            let full = ProfileTable::relaxed(dim, kappa_t, &kept);
                            let coarse = ProfileTable::relaxed(
                                dim,
                                kappa_t + 0.5 * smooth_sigma * smooth_sigma,
                                &kept,
                            );
                            detail.values = full.values;
                            for (v, s) in detail.values.iter_mut().zip(&coarse.values) {
                                *v -= s;
                            }
                            feature_weight = kept(RELAX_FEATURE_WAVELENGTH * c.diameter);
                        } else {
                            for (v, s) in detail.values.iter_mut().zip(&smooth.values) {
                                *v -= s;
                            }
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
                            detail_weight: feature_weight,
                            seed: c.id,
                            deep,
                            fine_blur: detail.sigma,
                            coarse_blur: smooth.sigma,
                            maturity_timescale: inputs.maturity_timescale,
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
                            &deform,
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
                            &deform,
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
                    Event::Patera(k) => {
                        let out =
                            io::patera_compute(&layers, &grid, body_radius, &paterae[*k], stamper);
                        Some((out, None, None))
                    }
                    Event::Mountain(k) => {
                        let out = io::mountain_compute(
                            &layers,
                            &grid,
                            body_radius,
                            &mountains[*k],
                            stamper,
                        );
                        Some((out, None, None))
                    }
                    Event::IoFlow(k) => {
                        let out =
                            io::flow_compute(&layers, &grid, body_radius, &io_flows[*k], stamper);
                        Some((out, None, None))
                    }
                    Event::Groove(k) => {
                        let out = ice_shell::groove_compute(
                            &layers,
                            &grid,
                            body_radius,
                            &lanes[*k],
                            lanes[*k].time_gyr,
                            cell,
                            renewed,
                            &base,
                            loose_cover.as_deref(),
                            stamper,
                        );
                        Some((out, None, None))
                    }
                    Event::PolarTerrain => {
                        let out = ice_shell::polar_compute(
                            &layers,
                            &grid,
                            body_radius,
                            &shell_map,
                            &polar_edge,
                            &lineae,
                            polar_time,
                            0.0,
                            cell,
                            halo,
                            renewed,
                            &base,
                            stamper,
                        );
                        Some((out, None, None))
                    }
                    Event::FossilTerrain(k) => {
                        // Gardening since the episode spreads and lowers its structures,
                        // as it degrades the craters of the time.
                        let kappa_t = inputs.diffusivity_m2_per_myr
                            * inputs.flux.equivalent_gyr(episodes[*k].end, present)
                            * 1e3;
                        let out = ice_shell::polar_compute(
                            &layers,
                            &grid,
                            body_radius,
                            &episode_maps[*k],
                            &fossil_edges[*k],
                            &lineae,
                            episodes[*k].end,
                            kappa_t,
                            cell,
                            halo,
                            renewed,
                            &base,
                            stamper,
                        );
                        Some((out, None, None))
                    }
                    Event::Chaos(k) => {
                        let out = cryo::chaos_compute(&layers, &grid, cell, present, &chaos[*k]);
                        Some((out, None, None))
                    }
                    Event::Fabric(k) => {
                        let out = fabric::fabric_compute(
                            &layers,
                            &grid,
                            cell,
                            present,
                            &fabric[*k],
                            cracks.as_ref()?,
                            inputs.seed,
                        );
                        Some((out, None, None))
                    }
                    Event::Cryodome(k) => {
                        let out = cryo::cryodome_compute(
                            &layers,
                            &grid,
                            body_radius,
                            &cryodomes[*k],
                            dome_heights[*k],
                            stamper,
                        );
                        Some((out, None, None))
                    }
                    Event::Linea(k) => {
                        let l = &lineae[*k];
                        let kappa_t = inputs.diffusivity_m2_per_myr
                            * inputs.flux.equivalent_gyr(l.time_gyr, present)
                            * 1e3;
                        let out = lineae::linea_compute(
                            &layers,
                            &grid,
                            body_radius,
                            l,
                            (4.0 * kappa_t).sqrt(),
                            loose.as_ref(),
                            present,
                            stamper,
                        );
                        Some((out, None, None))
                    }
                    Event::Ring(k) => {
                        let r = &ring_systems[*k];
                        let kappa_t = inputs.diffusivity_m2_per_myr
                            * inputs.flux.equivalent_gyr(r.time_gyr, present)
                            * 1e3;
                        let out = rings::system_compute(
                            &layers,
                            &grid,
                            body_radius,
                            r,
                            (4.0 * kappa_t).sqrt(),
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
                        // A secondary carries its parent's time: where the ground was
                        // renewed from below since, it is gone with the parent.
                        if *secondary {
                            let k = grid.index(grid.locate(c.center).0);
                            if layers.clean.as_slice()[k] > 0.0
                                && layers.surface_age.as_slice()[k] as f64 > c.time_gyr + 1e-6
                            {
                                return None;
                            }
                        }
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
                        // On an ice shell the crater relaxes by wavelength (F70): its table
                        // carries the bowl's relaxation, and its terraces, rings and peak
                        // keep the share of a short wavelength.
                        let (table, kept, features) = if relaxes(dim.radius, c.time_gyr, c.center) {
                            let kept = |l: f64| crater_kept(l, c.time_gyr, c.center);
                            (
                                ProfileTable::relaxed(dim, kappa_t, &kept),
                                1.0,
                                kept(RELAX_FEATURE_WAVELENGTH * c.diameter),
                            )
                        } else {
                            let whole = relaxed[*j] as f64;
                            (table, whole, whole)
                        };
                        let depth_factor = if *secondary {
                            SECONDARY_DEPTH_FACTOR
                        } else {
                            1.0
                        };
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
                            amplitude: kept * depth_factor,
                            feature_amplitude: features * depth_factor,
                            replace_floor: true,
                            set_age: true,
                            bury: true,
                            seed: c.key,
                            deep,
                            maturity_timescale: inputs.maturity_timescale,
                            renewal_guard: false,
                            terraced: true,
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
    // On an ice shell over rock the domes' substrate is partly clean ice, and the deep
    // young craters' floors carry salts.
    if ice::is_ice_surface(graph) {
        cryo::post_pass(
            &mut layers,
            body_radius,
            &cryodomes,
            &records,
            &inputs.scaling,
            present,
        );
    }
    // The plumes' fallout mantles the ground it has lain on since each surface formed
    // (F74), before the unresolved craters are laid: those dug before the cover above them
    // grew thinner than they are deep lie buried, the younger ones pit its surface.
    let mantle = loose
        .as_ref()
        .filter(|_| std::env::var_os("EXOGEN_SKIP_MANTLE").is_none());
    let buried: Option<Vec<f32>> = match (mantle, cover_on(&layers)) {
        (Some(l), Some(cover)) => {
            cryo::mantle(&mut layers, body_radius, &cover);
            let depth = dimensions(d_min, &inputs.scaling).depth;
            let ages = layers.surface_age.as_slice();
            Some(
                (0..grid.cell_count())
                    .into_par_iter()
                    .map(|k| {
                        l.buried_before(grid.center(grid.cell(k)), ages[k] as f64, present, depth)
                            as f32
                    })
                    .collect(),
            )
        }
        _ => None,
    };
    // The craters the finish does not resolve, as the relief their population leaves on
    // each surface by its age (F54): rough where old, smooth where renewed.
    roughness::apply(
        &mut layers,
        body_radius,
        cell,
        d_min,
        &inputs,
        present,
        Some(&renewal.keep),
        buried.as_deref(),
    );
    let (sea_level, ocean_fraction) = hydrosphere_fill(graph, &mut layers);
    let height = layers.height();
    // On ice, each slope sits at its own mean temperature: frost and lag follow it.
    let slope = (ice::is_ice_surface(graph) || io::is_io_regime(graph))
        .then(|| ice::poleward_slope(&height, body_radius));
    // Io's surface is what its plumes have left on it, each slope at its own temperature
    // (F84).
    if io::is_io_regime(graph) {
        let clock = std::time::Instant::now();
        let background = io::background_fall(
            grid,
            body_radius,
            graph.f64("body.gravity"),
            &paterae,
            graph.f64("seed.volcanism") as u64,
        );
        io::deposit_history(
            &mut layers,
            body_radius,
            graph.f64("orbit.insolation"),
            &plume_episodes,
            &small_eruptions,
            &background,
            slope.as_ref(),
        );
        if std::env::var_os("EXOGEN_DEBUG").is_some() {
            eprintln!(
                "io deposit history: {} eruptions, {} erupting now, {} small ones, background {:.4} to {:.4} optical depths a year, {:.1} s",
                plume_episodes.len(),
                plume_episodes.iter().filter(|e| e.end_yr == 0.0).count(),
                small_eruptions.len(),
                background.iter().cloned().fold(f32::MAX, f32::min),
                background.iter().cloned().fold(f32::MIN, f32::max),
                clock.elapsed().as_secs_f64()
            );
        }
    }

    sink.report(Progress {
        stage: "finishing",
        fraction: 0.8,
        message: "weathering and rays".to_string(),
    });
    if std::env::var_os("EXOGEN_DEBUG").is_some() {
        eprintln!(
            "finish phases: the event sequence {:.1} s",
            stamp_clock.elapsed().as_secs_f64()
        );
    }
    // Maturity: space weathering on rock; on ice, the sublimation lag the latitude's
    // temperature and the surface's age have left.
    let tau_sw = graph.f64("surface.space_weathering_timescale");
    let (maturity, lag) = if ice::is_ice_surface(graph) {
        let (maturity, lag) = ice::lag_fields(graph, &layers, present, slope.as_ref());
        (maturity, Some(lag))
    } else {
        let maturity = Field::from_fn(grid, |c, _| {
            let exposure = (present - layers.surface_age.get(c) as f64).max(0.0);
            (1.0 - (-exposure / tau_sw).exp()) as f32
        });
        (maturity, None)
    };
    // The skin a ray's secondaries dig through (F91): the Moon's, or an icy surface's lag
    // where it is thicker, which leaves them less of the area they would excavate.
    let skin = lag.map(|lag| {
        Field::from_fn(grid, |c, _| {
            let depth = (lag.get(c) as f64).max(RAY_SKIN_DEPTH_M);
            (depth / RAY_SKIN_DEPTH_M).powf(-RAY_SKIN_SLOPE) as f32
        })
    });
    // Steep ground sheds its regolith and stays fresher than flat ground (F55): the
    // bright walls and scarps of every airless body, from the ground's own slopes.
    let maturity = {
        let steep = ice::slope_magnitude(&height, body_radius);
        let (start, repose) = (
            SLOPE_SHED_START_DEG.to_radians(),
            REPOSE_SLOPE_DEG.to_radians(),
        );
        let mut m = maturity;
        m.as_mut_slice()
            .par_iter_mut()
            .zip(steep.as_slice())
            .for_each(|(v, &s)| {
                let shed = ((s as f64 - start) / (repose - start)).clamp(0.0, 1.0);
                *v *= 1.0 - (SLOPE_FRESHNESS * shed) as f32;
            });
        m
    };
    let rays = ray_field(
        &grid,
        body_radius,
        &records,
        &inputs,
        present,
        cell,
        skin.as_ref(),
    );

    sink.report(Progress {
        stage: "finishing",
        fraction: 0.92,
        message: "materials".to_string(),
    });
    let albedo = materials::albedo(graph, &layers, &maturity, &rays, present, slope.as_ref());
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
    let (grooved_fraction, polar_fraction) = ice_shell::terrain_shares(&layers);
    // The share the lanes formed by this epoch cover, out of all of them.
    let lane_coverage =
        loop_output.lane_coverage * lanes.len() as f64 / loop_output.lanes.len().max(1) as f64;
    let chaos_fraction = cryo::chaos_share(&layers);
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
        provinces: loop_output.provinces.clone(),
        ridges,
        rises: loop_output.rises.clone(),
        shields,
        canyons,
        graben,
        resurfacing,
        coronae,
        domes,
        paterae,
        mountains,
        io_flows,
        plume_episodes,
        small_eruption_count: small_eruptions.len(),
        lineae,
        chaos,
        cryodomes,
        chaos_fraction,
        lanes,
        lane_coverage,
        grooved_fraction,
        polar_fraction,
        episodes,
        renewed_fraction: loop_output.renewed_fraction,
        hash,
    })
}

/// The ice shell through time for relaxation, when the surface is ice: its thickness and
/// basal temperature from the thermal history, the surface temperature from the climate.
fn ice_shell_history(
    graph: &ParamGraph,
    thermal: &ThermalHistory,
) -> Option<Box<dyn Fn(f64) -> relaxation::Shell + Sync>> {
    if !ice::is_ice_surface(graph) {
        return None;
    }
    let thermal = thermal.clone();
    let surface_temperature = graph.f64("atmosphere.surface_temperature_mean");
    let gravity = graph.f64("body.gravity");
    let rock = crate::stages::materials::non_ice_fraction(graph.str("materials.bedrock"));
    let viscosity_scale = 10f64.powf(ICE_ROCK_VISCOSITY_DECADES * rock);
    Some(Box::new(move |t: f64| {
        let state = thermal.at(t);
        relaxation::Shell {
            thickness: state.lithosphere,
            elastic: state.elastic,
            surface_temperature,
            basal_temperature: state.mantle_temperature,
            gravity,
            viscosity_scale,
        }
    }))
}

/// A launch azimuth drawn in proportion to an ejecta curtain's mass: candidates on one
/// stream, each accepted on another by its mass against the largest.
fn jet_azimuth(curtain: &Curtain, rng: &Rng, index: u64, slot: u64) -> f64 {
    let top = curtain.max_mass().max(1e-12);
    let mut phi = 0.0;
    for attempt in 0..64u64 {
        let draw = slot * 64 + attempt;
        phi = 2.0 * std::f64::consts::PI * rng.uniform(index ^ 0x7E75_0000_0000_0000, draw);
        if rng.uniform(index.wrapping_add(ACCEPT_STREAM), draw) * top <= curtain.mass(phi) {
            break;
        }
    }
    phi
}

/// Secondary craters of the loop craters: groups at a few primary radii, formed at the
/// primary's time, sized as a small fraction of the primary. Each group lands along a jet
/// of the primary's ejecta curtain, drawn by the jets' mass, on the jet's bent path
/// (F91), so the secondaries lie in the rays.
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
        let curtain = Curtain::new(p.id, p.azimuth, p.angle);
        let groups = (count / SECONDARY_GROUP_SIZE).max(1);
        let launches: Vec<f64> = (0..groups as u64)
            .map(|g| jet_azimuth(&curtain, &rng, p.id, g))
            .collect();
        for j in 0..count as u64 {
            let group = j % groups as u64;
            let rank = (j / groups as u64) as f64;
            let base = 8 + group * 4;
            let dist0 = radius * (r_lo + (r_hi - r_lo) * rng.uniform(p.id, base));
            let launch = launches[group as usize];
            let path = |d: f64| curtain.path_azimuth(launch, (d / radius).max(1.0).ln());
            let chain = rng.uniform(p.id, base + 2) < SECONDARY_CHAIN_FRACTION;
            let draw = 1_000 + j * 4;
            // A chain runs outward along the radial from the primary; a cluster scatters
            // about its centre.
            let (dist, phi) = if chain {
                let spacing =
                    radius * SECONDARY_CHAIN_SPACING * (0.7 + 0.6 * rng.uniform(p.id, draw));
                let lateral = 0.03 * radius * (rng.uniform(p.id, draw + 1) - 0.5);
                let d = dist0 + rank * spacing;
                (d, path(d) + lateral / d.max(1.0))
            } else {
                let jitter_r = 0.3 * radius * (rng.uniform(p.id, draw) - 0.5);
                let jitter_phi =
                    0.5 * (rng.uniform(p.id, draw + 1) - 0.5) * radius / dist0.max(1.0);
                let d = dist0 + jitter_r;
                (d, path(d) + jitter_phi)
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
/// A crater whose rays could reach more cells than this is drawn over the whole grid in
/// parallel instead of being gathered: its footprint would hold too much.
const RAY_GATHER_MAX_CELLS: f64 = 5.0e5;

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
    skin: Option<&Field<f32>>,
) -> Field<f32> {
    let skin = skin.map(|s| s.as_slice());
    let mut rays = Field::new(*grid);
    let out = rays.as_mut_slice();
    let lifetime = inputs.ray_lifetime.max(1e-6);
    let systems: Vec<RaySystem> = craters
        .iter()
        .filter(|c| {
            let age = present - c.time_gyr;
            age >= 0.0 && age < lifetime && c.diameter >= RAY_MIN_CELLS * cell
        })
        .collect::<Vec<_>>()
        .par_iter()
        .map(|c| RaySystem::new(c, inputs, present, cell))
        .collect();
    // A system whose rays could cover many cells is drawn over the whole grid in
    // parallel; the rest are gathered, a batch of systems at a time.
    let cap_cells = |reach: f64| {
        let theta = (reach / body_radius).min(std::f64::consts::PI);
        2.0 * std::f64::consts::PI * body_radius * body_radius * (1.0 - theta.cos()) / (cell * cell)
    };
    let (large, small): (Vec<&RaySystem>, Vec<&RaySystem>) = systems
        .iter()
        .partition(|s| cap_cells(s.reach) > RAY_GATHER_MAX_CELLS);
    let clock = std::time::Instant::now();
    // Each system's share joins those under it as independent shares do, 1 - (1 - a)(1 -
    // b), whatever their order.
    for batch in small.chunks(RAY_BATCH) {
        let parts: Vec<Vec<(u32, f32)>> = batch
            .par_iter()
            .map(|s| with_stamper(grid, |stamper| s.gathered(grid, body_radius, skin, stamper)))
            .collect();
        for part in parts {
            for (idx, v) in part {
                let o = &mut out[idx as usize];
                *o = 1.0 - (1.0 - *o) * (1.0 - v);
            }
        }
    }
    for s in &large {
        let (east, north) = local_basis(s.center);
        let theta_max = (s.reach / body_radius).min(std::f64::consts::PI) + grid.delta();
        out.par_iter_mut().enumerate().for_each(|(k, o)| {
            let d = grid.center(grid.cell(k));
            let theta = angle_between(s.center, d);
            if theta > theta_max {
                return;
            }
            let v = s.value(
                theta * body_radius,
                azimuth_from(s.center, east, north, d),
                skin.map_or(1.0, |s| s[k] as f64),
            );
            if v > 0.0 {
                *o = 1.0 - (1.0 - *o) * (1.0 - v);
            }
        });
    }
    if std::env::var_os("EXOGEN_DEBUG").is_some() {
        eprintln!(
            "rays: {} systems, {} over the whole grid, {:.1} s",
            systems.len(),
            large.len(),
            clock.elapsed().as_secs_f64()
        );
    }
    rays
}

/// One young crater's brightening (F91): full inside the rim; its continuous blanket,
/// falling with the deposit's cover to half at its lobed edge and on past it, mottled;
/// and beyond, its rays, the share of the ground its ejecta curtain's fragments
/// excavated from under the weathered skin. All of it fades with the crater's age.
struct RaySystem {
    center: [f64; 3],
    radius: f64,
    azimuth: f64,
    angle: f64,
    freshness: f64,
    edge: EdgeShape,
    ejecta_edge: f64,
    halo: f64,
    texture: Option<CraterTexture>,
    rays: Rays,
    reach: f64,
}

impl RaySystem {
    fn new(c: &CraterRecord, inputs: &impacts::CraterInputs, present: f64, cell: f64) -> Self {
        let lifetime = inputs.ray_lifetime.max(1e-6);
        let age = present - c.time_gyr;
        let radius = 0.5 * c.diameter;
        let dim = dimensions(c.diameter, &inputs.scaling);
        let rays = Rays::new(Curtain::new(c.id, c.azimuth, c.angle), radius);
        let halo = EJECTA_RADIUS_MULTIPLE * radius;
        let reach = rays.reach().max(halo);
        Self {
            center: c.center,
            radius,
            azimuth: c.azimuth,
            angle: c.angle,
            freshness: (1.0 - age / lifetime).powi(2),
            edge: EdgeShape::new(c.id, dim.ejecta_radius, cell),
            ejecta_edge: CONTINUOUS_EJECTA_MULTIPLE * radius,
            halo,
            texture: (c.diameter >= RAY_TEXTURE_MIN_CELLS * cell)
                .then(|| CraterTexture::new(c.id, &dim)),
            rays,
            reach,
        }
    }

    /// The brightening at distance `r`, m, and azimuth `phi` from the centre, where the
    /// skin leaves the share `skin` of the area the rays' secondaries would excavate.
    fn value(&self, r: f64, phi: f64, skin: f64) -> f32 {
        if r <= self.radius {
            return self.freshness as f32;
        }
        let edge = self.ejecta_edge * self.edge.factor_at(phi, r / self.radius);
        // The blanket is bright where its deposit covers the ground (F51), fading with the
        // cover past its lobed edge rather than ending there, mottled.
        let blanket = (ejecta_coverage(r, self.radius, edge, self.halo)
            * ejecta_modulation(phi, self.azimuth, self.angle).min(1.0)
            * self.texture.as_ref().map_or(1.0, |t| {
                1.0 + EJECTA_BRIGHTNESS_TEXTURE * t.noise_floor(phi, r / self.radius)
            }))
        .clamp(0.0, 1.0);
        let ray = self.rays.excavated(r, phi, skin);
        ((1.0 - (1.0 - blanket) * (1.0 - ray)) * self.freshness) as f32
    }

    /// The cells the system brightens, gathered within its reach.
    fn gathered(
        &self,
        grid: &CubeSphere,
        body_radius: f64,
        skin: Option<&[f32]>,
        stamper: &mut Stamper,
    ) -> Vec<(u32, f32)> {
        stamper.gather(grid, body_radius, self.center, self.reach);
        stamper
            .footprint
            .iter()
            .filter_map(|&(idx, r, phi)| {
                let v = self.value(r, phi, skin.map_or(1.0, |s| s[idx] as f64));
                (v > 0.0).then_some((idx as u32, v))
            })
            .collect()
    }
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
    #[test]
    fn time_bins_tile_every_quality_the_same_way() {
        use super::{LOOP_STEPS, LOOP_TIME_BINS, time_bins};
        for &steps in &LOOP_STEPS {
            assert_eq!(LOOP_TIME_BINS % steps, 0, "{steps} steps");
        }
        let age = 4.52;
        let mut seen = Vec::new();
        let dt = age / 30.0;
        for s in 1..=30 {
            seen.extend(time_bins(age, dt * (s - 1) as f64, dt * s as f64));
        }
        assert_eq!(seen.len(), LOOP_TIME_BINS);
        assert!(seen.iter().enumerate().all(|(i, b)| b.0 == i as u64));
        let dt = age / 90.0;
        let fine: Vec<_> = (1..=90)
            .flat_map(|s| time_bins(age, dt * (s - 1) as f64, dt * s as f64))
            .collect();
        assert_eq!(fine.len(), LOOP_TIME_BINS);
        assert!(
            fine.iter()
                .zip(&seen)
                .all(|(a, b)| a.0 == b.0 && (a.1 - b.1).abs() < 1e-12)
        );
    }

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
        // A rift population: at least one trough the loop carried, and many graben.
        assert!(!lo.canyons.is_empty(), "no loop-carried canyon");
        assert!(lo.graben.len() > 20, "{} graben", lo.graben.len());
        let grid = out.finish.grid;
        let radius = g.f64("body.radius");
        let height = out.finish.height.as_slice();
        let at = |d: [f64; 3]| height[grid.index(grid.locate(d).0)] as f64;
        let rise = &lo.rises[0];
        let (east, north) = crate::geom::local_basis(rise.center);
        // The rise stands above the ground beyond it: the mean over its crown against the
        // mean round a ring beyond its foot, many samples each so that craters, shields
        // and the canyon do not decide it.
        let ring_mean = |frac: f64, n: usize| -> f64 {
            (0..n)
                .map(|k| {
                    let phi = k as f64 * 2.0 * std::f64::consts::PI / n as f64;
                    at(crate::geom::offset(
                        rise.center,
                        east,
                        north,
                        frac * rise.radius / radius,
                        phi,
                    ))
                })
                .sum::<f64>()
                / n as f64
        };
        let crown = 0.5 * (ring_mean(0.1, 16) + ring_mean(0.2, 32));
        let far = ring_mean(1.2, 64);
        assert!(
            crown - far > 0.3 * rise.height,
            "crown {crown} vs far {far}, rise {}",
            rise.height
        );
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
        // Some trough the loop carried still stands kilometres deep below the plateau
        // either side, the rims being the highest ground within a few widths past the
        // en-echelon neighbours; a later flood may have filled another.
        let relief = |canyon: &plume::Trough| {
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
            let rim = |sign: f64| {
                [1.0, 1.5, 2.0, 2.5]
                    .iter()
                    .map(|&s| side(sign * s))
                    .fold(f64::NEG_INFINITY, f64::max)
            };
            0.5 * (rim(1.0) + rim(-1.0)) - at(mid)
        };
        let deepest = lo
            .canyons
            .iter()
            .map(relief)
            .fold(f64::NEG_INFINITY, f64::max);
        assert!(deepest > 3_000.0, "deepest standing trough {deepest} m");
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
        let expected = g.f64("body.age") - g.f64("tectonics.resurfacing_event_age");
        assert!((event - expected).abs() < 1e-9);
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

    /// Europa's cracks open perpendicular to the most tensile tidal stress of their moment.
    #[test]
    fn europa_cracks_follow_the_tidal_stress() {
        let mut r = Runner::new();
        let g = preset_graph("Europa", "preview");
        let out = r.run(&g, &NoProgress).unwrap();
        let cracks = lineae::Cracks::new(&g).expect("Europa cracks");
        let lineae = &out.finish.lineae;
        assert!(lineae.len() > 200, "{} lineae", lineae.len());
        let cycloids = lineae
            .iter()
            .filter(|l| l.kind == lineae::LineaKind::Cycloid)
            .count();
        assert!(cycloids > 0, "no cycloids");
        // Every ridge runs within twenty degrees of the direction perpendicular to the
        // most tensile stress at its midpoint, in the field of its own time.
        let present = out.finish.epoch_time_gyr;
        let mut within = 0;
        let mut total = 0;
        for l in lineae
            .iter()
            .filter(|l| l.kind != lineae::LineaKind::Cycloid)
        {
            let n = l.points.len();
            if n < 5 {
                continue;
            }
            let (a, b) = (l.points[n / 2 - 1], l.points[n / 2 + 1]);
            let mid = l.points[n / 2];
            let chord = crate::geom::normalize([b[0] - a[0], b[1] - a[1], b[2] - a[2]]);
            let axis = lineae::mean_axis_at(l.time_gyr, present, cracks.nsr_period_gyr);
            let stress = cracks.shell.total(axis, mid, l.phase);
            let d = lineae::crack_direction(mid, &stress, Some(chord));
            let cos = crate::geom::dot(chord, d).abs().min(1.0);
            if cos.acos().to_degrees() < 20.0 {
                within += 1;
            }
            total += 1;
        }
        assert!(total > 50, "{total} ridges measured");
        assert!(
            within as f64 > 0.8 * total as f64,
            "{within} of {total} ridges follow the field"
        );
        let share = lineae::crack_fraction(&out.finish.layers);
        assert!(share > 0.005, "crack share {share}");
    }

    /// Europa's thin shell breaks into chaos where its heat melts it through, and its
    /// plains are as young as its surface: the fabric of narrow ridges renews them (F81),
    /// so the surface's mean age is its retention age and little of it is older than
    /// three of them.
    #[test]
    fn europa_breaks_into_chaos_and_renews_its_plains() {
        let mut r = Runner::new();
        let g = preset_graph("Europa", "preview");
        let out = r.run(&g, &NoProgress).unwrap();
        assert!(
            out.finish.chaos.len() == CHAOS_GENERATIONS,
            "{} chaos generations",
            out.finish.chaos.len()
        );
        let share = out.finish.chaos_fraction;
        assert!(
            (0.08..=0.4).contains(&share),
            "chaos matrix share {share:.3}"
        );
        let layers = &out.finish.layers;
        let grid = layers.grid;
        let present = g.f64("body.age");
        let tau = g.f64("cratering.retention_age");
        let (mut total, mut mean, mut old) = (0.0, 0.0, 0.0);
        for (k, &age) in layers.surface_age.as_slice().iter().enumerate() {
            let a = grid.cell_area(grid.cell(k));
            let before = present - age as f64;
            total += a;
            mean += a * before;
            if before > 3.0 * tau {
                old += a;
            }
        }
        let (mean, old) = (mean / total, old / total);
        assert!(
            (0.5 * tau..=1.5 * tau).contains(&mean),
            "mean surface age {mean:.4} Gyr against a retention of {tau}"
        );
        assert!(
            old < 0.1,
            "{old:.3} of the surface older than three retention ages"
        );
    }

    /// Ceres builds domes of brine-rich ice and brightens the floors of its deep young
    /// craters with salts.
    #[test]
    fn ceres_builds_domes_and_bright_floors() {
        let mut r = Runner::new();
        let g = preset_graph("Ceres", "preview");
        let out = r.run(&g, &NoProgress).unwrap();
        assert!(
            out.finish.cryodomes.len() >= 5,
            "{} domes",
            out.finish.cryodomes.len()
        );
        let layers = &out.finish.layers;
        let units = layers.unit_id.as_slice();
        let cover = layers.cover.as_slice();
        let salted = units
            .iter()
            .zip(cover)
            .filter(|(u, c)| {
                (**u == crate::units::UNIT_CRATER_FLOOR || **u == crate::units::UNIT_DEEP_FLOOR)
                    && **c > 0.0
            })
            .count();
        assert!(salted > 0, "no salted crater floors");
    }

    /// Enceladus's plume keeps its surface in fresh frost and buries the salts in its
    /// cracks before they can show.
    #[test]
    fn enceladus_fallout_keeps_it_bright() {
        let mut r = Runner::new();
        let g = preset_graph("Enceladus", "preview");
        let out = r.run(&g, &NoProgress).unwrap();
        let m = out.finish.maturity.as_slice();
        let mean = m.iter().map(|&v| v as f64).sum::<f64>() / m.len() as f64;
        assert!(mean < 0.1, "mean maturity {mean:.3}");
        let units = out.finish.layers.unit_id.as_slice();
        let albedo = out.finish.albedo.as_slice();
        let (mut sum, mut n) = ([0.0f64; 3], 0usize);
        for (k, &u) in units.iter().enumerate() {
            if u == crate::units::UNIT_LINEA || u == crate::units::UNIT_BAND {
                for (ch, acc) in sum.iter_mut().enumerate() {
                    *acc += albedo[k][ch] as f64;
                }
                n += 1;
            }
        }
        assert!(n > 0, "no cracks");
        let (r, b) = (sum[0] / n as f64, sum[2] / n as f64);
        assert!(r - b < 0.03, "cracks are brown: r {r:.3} b {b:.3}");
    }

    /// Ganymede has two terrains: the old dark cratered one and the lanes of bright
    /// grooved terrain that cut it, which carry far fewer craters.
    #[test]
    fn ganymede_has_two_terrains() {
        let mut r = Runner::new();
        let g = preset_graph("Ganymede", "preview");
        let out = r.run(&g, &NoProgress).unwrap();
        let grooved = out.finish.lane_coverage;
        assert!(
            (0.45..=0.85).contains(&grooved),
            "{} lanes cover {grooved:.2} of the surface",
            out.finish.lanes.len()
        );
        let body_radius = g.f64("body.radius");
        let (mut inside, mut outside) = (0usize, 0usize);
        for c in out.finish.craters.iter().filter(|c| c.visible) {
            if out
                .finish
                .lanes
                .iter()
                .any(|l| l.inside(c.center, body_radius))
            {
                inside += 1;
            } else {
                outside += 1;
            }
        }
        assert!(outside >= 10, "{outside} craters on the dark terrain");
        let density_in = inside as f64 / grooved;
        let density_out = outside as f64 / (1.0 - grooved);
        eprintln!(
            "ganymede: {inside} recorded craters on the {grooved:.2} of the surface the lanes cover, {outside} on the rest; density ratio {:.2}",
            density_in / density_out
        );
        assert!(
            density_in < 0.5 * density_out,
            "{inside} craters on {grooved:.2} of the surface against {outside} on the rest"
        );
    }

    /// Enceladus carries the terrains of its past episodes of focusing (F53): several of
    /// them, which with the present one and the lanes take a good share of the surface
    /// from the cratered plains.
    #[test]
    fn enceladus_keeps_the_terrains_of_its_past_episodes() {
        let mut r = Runner::new();
        let g = preset_graph("Enceladus", "preview");
        let out = r.run(&g, &NoProgress).unwrap();
        let n = out.finish.episodes.len();
        assert!(n >= 3, "{n} past episodes");
        let renewed = out.finish.renewed_fraction;
        eprintln!(
            "enceladus: {n} past episodes; renewed {renewed:.3} of the surface, polar terrain units {:.3}, lanes {:.3}",
            out.finish.polar_fraction, out.finish.lane_coverage
        );
        // Kirchoff and Schenk 2009 map 43% of the real body as cratered terrain.
        assert!((0.35..=0.75).contains(&renewed), "renewed {renewed:.3}");
    }

    /// Ganymede's ice relaxes its large craters into palimpsests and keeps its small ones
    /// shallow but present, as the dark terrain's craters are.
    #[test]
    fn ganymede_relaxes_its_large_craters_first() {
        let mut r = Runner::new();
        let g = preset_graph("Ganymede", "preview");
        let out = r.run(&g, &NoProgress).unwrap();
        let grid = out.finish.grid;
        let radius = g.f64("body.radius");
        let height = out.finish.height.as_slice();
        let at = |d: [f64; 3]| height[grid.index(grid.locate(d).0)] as f64;
        let ring_mean = |center: [f64; 3], r_m: f64| -> f64 {
            let (east, north) = crate::geom::local_basis(center);
            (0..16)
                .map(|k| {
                    let phi = k as f64 * std::f64::consts::PI / 8.0;
                    at(crate::geom::offset(center, east, north, r_m / radius, phi))
                })
                .sum::<f64>()
                / 16.0
        };
        let depth_ratio = |lo: f64, hi: f64| -> f64 {
            let mut ratios: Vec<f64> = out
                .finish
                .craters
                .iter()
                .filter(|c| c.diameter >= lo && c.diameter < hi)
                .map(|c| {
                    let rim = ring_mean(c.center, 0.5 * c.diameter);
                    // The flat floor between a central pit's rim and the wall's foot
                    // (F72).
                    let floor = ring_mean(c.center, 0.22 * c.diameter);
                    (rim - floor) / c.diameter
                })
                .collect();
            assert!(ratios.len() >= 10, "{} craters in {lo}..{hi}", ratios.len());
            ratios.sort_by(f64::total_cmp);
            ratios[ratios.len() / 2]
        };
        let small = depth_ratio(30e3, 70e3);
        let large = depth_ratio(150e3, 2_000e3);
        eprintln!("ganymede relaxation: median depth ratio small {small:.4}, large {large:.4}");
        let (n_loop, n_all) = out
            .finish
            .craters
            .iter()
            .filter(|c| c.diameter >= 30e3 && c.diameter < 70e3)
            .fold((0usize, 0usize), |(l, a), c| {
                (l + c.from_loop as usize, a + 1)
            });
        eprintln!("ganymede relaxation: {n_loop} of {n_all} small craters came from the loop");
        // The dark terrain's 60 to 70 km craters are about a kilometre deep on the real
        // body (Schenk 2002), a depth ratio near 0.012: relaxed, not flattened.
        assert!(small > 0.008, "small craters keep their depth: {small}");
        assert!(large < 0.25 * small, "large {large} vs small {small}");
    }

    /// Io keeps no craters, carries hundreds of paterae and dozens of mountains, and its
    /// plains are sulfur yellow.
    #[test]
    fn io_is_craterless_with_paterae_and_mountains() {
        let mut r = Runner::new();
        let g = preset_graph("Io", "preview");
        let out = r.run(&g, &NoProgress).unwrap();
        let lo = &out.loop_output;
        assert!(lo.paterae.len() > 100, "{} paterae", lo.paterae.len());
        assert!(lo.mountains.len() > 20, "{} mountains", lo.mountains.len());
        assert!(!out.finish.io_flows.is_empty() && !out.finish.plume_episodes.is_empty());
        assert!(
            out.finish.craters.is_empty(),
            "{} craters",
            out.finish.craters.len()
        );
        let units = out.finish.layers.unit_id.as_slice();
        let flows = units
            .iter()
            .filter(|&&u| u == crate::units::UNIT_IO_FLOW)
            .count() as f64
            / units.len() as f64;
        assert!(flows > 0.005 && flows < 0.3, "flow share {flows}");
        let albedo = out.finish.albedo.as_slice();
        let mean: [f64; 3] = albedo.iter().fold([0.0; 3], |a, c| {
            [a[0] + c[0] as f64, a[1] + c[1] as f64, a[2] + c[2] as f64]
        });
        let n = albedo.len() as f64;
        let (rr, gg, bb) = (mean[0] / n, mean[1] / n, mean[2] / n);
        assert!(
            rr > gg && gg > bb && rr > 0.3,
            "mean colour {rr:.2} {gg:.2} {bb:.2}"
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
