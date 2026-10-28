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

pub mod corona;
pub mod crust;
pub mod cryo;
pub mod deform;
pub mod episodic;
pub mod fabric;
pub mod ice;
pub mod ice_shell;
pub mod impacts;
pub mod io;
pub mod lineae;
pub mod load;
pub mod materials;
pub mod mega;
pub mod plume;
pub mod primordial_crust;
pub mod rifting;
pub mod relaxation;
pub mod ridges;
pub mod rings;
pub mod roughness;
pub mod strain;
pub mod stress;
pub mod tectonics;
pub mod tessera;
pub mod tidal;
pub mod volcanism;

use std::collections::HashMap;
use std::sync::Arc;
use std::time::{Duration, Instant};

use rayon::prelude::*;

use crate::constants::*;
use crate::crater::curtain::{
    Cluster, Grain, largest_secondary, secondaries_above, secondary_size_at,
};
use crate::crater::production::poisson;
use crate::crater::{
    CraterRecord, CraterTexture, Curtain, EdgeShape, PopulationConfig, ProfileTable, Rays,
    RestoreParams, StampOutput, StampParams, Stamper, dimensions, ejecta_coverage,
    ejecta_modulation, restore_compute, sample_population, stamp_compute, with_stamper,
};
use crate::fields::{Field, LayerSet};
use crate::geom::{angle_between, azimuth_from, local_basis, offset};
use crate::grid::CubeSphere;
use crate::params::ParamGraph;
use crate::params::derive::physics::volatile;
use crate::rng::{Rng, hash_words, mix64};
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
    /// The mantle's convection beneath the lid, as it stands under the map (F173).
    pub support: Option<crate::convection::Support>,
    /// Plume rises, their shields, canyons and graben.
    pub rises: Vec<plume::Rise>,
    pub shields: Vec<plume::Shield>,
    pub canyons: Vec<plume::Trough>,
    pub graben: Vec<plume::Trough>,
    /// An episodic lid's resurfacing event, its coronae over the mantle's diapirs (F145)
    /// and its domes.
    pub resurfacing: Option<episodic::Resurfacing>,
    pub coronae: Vec<corona::Corona>,
    pub domes: Vec<episodic::Dome>,
    /// How far the ground sank under each load the loop laid, per coarse cell, and when
    /// (F162), which the finish's floods see after them.
    pub flexure: Vec<(f64, Vec<f32>)>,
    /// The intrusions laid since the last epoch at each, m per coarse cell, and when (F174).
    pub intrusions: Vec<(f64, Vec<f32>)>,
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
    /// The coronae over an ice mantle's upwellings and when they formed (F133).
    pub ice_coronae: Vec<corona::Corona>,
    /// The share of the surface the lanes, the coronae and the polar terrains, present
    /// and past, renewed between them, by area.
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
    pub coronae: Vec<corona::Corona>,
    pub domes: Vec<episodic::Dome>,
    /// The Io regime's features, all of the present surface, and the plumes of its
    /// paterae over the history the surface remembers (F84).
    pub paterae: Vec<io::Patera>,
    pub mountains: Vec<io::Mountain>,
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
    /// The coronae over an ice mantle's upwellings present at the epoch (F133).
    pub ice_coronae: Vec<corona::Corona>,
    /// The share of the surface the lanes, the coronae and the polar terrains, present
    /// and past, renewed between them, by area, whatever craters they have taken since.
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
            relaxation::VERSION,
            impacts::VERSION,
            tectonics::VERSION,
            volcanism::VERSION,
            ridges::VERSION,
            plume::VERSION,
            crate::convection::VERSION,
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
pub const FINISH_VERSION: u64 = 160;

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
    // The ground before the rises, low-passed over the lid's flexural wavelength, for the
    // stress its potential energy puts on the lid (F150).
    let relief = {
        let elastic = graph.f64("interior.elastic_thickness");
        let wavelength =
            stress::StressField::flexural_wavelength(elastic, graph.f64("body.gravity"));
        let cell = grid.delta() * graph.f64("body.radius");
        stress::Relief {
            grid,
            height: relaxation::low_pass(&layers.bedrock, cell, wavelength),
        }
    };
    let mut plume = plume::PlumeProcess::new(graph, grid, &thermal, Some(relief));
    // The mantle's convection beneath the lid (F173), turned with the body as its rises
    // turned it (F104), and the ground its flow holds up on each coarse cell so far.
    let support = crate::convection::Convection::from_graph(graph, &thermal).map(|convection| {
        crate::convection::Support {
            convection,
            turn: plume.as_ref().map_or_else(Default::default, |p| p.turn),
        }
    });
    let mut lifted = vec![0.0f32; grid.cell_count()];
    // The lid's field through the loop (F154, F157), which the thrusts read from its other
    // side (F155): the ground and the lid read at epochs, and the lava laid since each load
    // bending the lid of its time.
    let body_radius = graph.f64("body.radius");
    let gravity = graph.f64("body.gravity");
    let cell = grid.delta() * body_radius;
    let mut history = stress::FieldHistory::new(
        graph,
        plume.as_ref().map(|p| p.rises.clone()).unwrap_or_default(),
        plume.as_ref().map(|p| p.limbs.clone()).unwrap_or_default(),
    );
    let first_lid = thermal.at(0.0).elastic;
    // The stress the mantle's flow bends the lid with as it stands at a time (F182).
    let flow_at = |t: f64, elastic: f64| {
        let drag = thermal
            .mantle
            .as_ref()
            .map_or(0.0, |m| crate::convection::drag_stress(m, t));
        support
            .as_ref()
            .filter(|_| std::env::var_os("EXOGEN_SKIP_FLOW").is_none())
            .map(|s| {
                load::bent(
                    &|d| s.at(t, d),
                    &|d| drag * s.pattern(t, d),
                    body_radius,
                    gravity,
                    elastic,
                )
            })
    };
    history.push(stress::FieldEpoch {
        time: 0.0,
        elastic: first_lid,
        relief: stress::Relief::read(&layers.bedrock, None, cell, first_lid, gravity),
        load: None,
        flow: flow_at(0.0, first_lid),
    });
    let mut contraction = tectonics::ContractionProcess::new(graph, grid, &thermal, &history);
    // The lid's tension faults, opened through the loop as its field stretches it (F178),
    // and those no rise carves.
    let mut rifts = rifting::RiftProcess::new(graph, grid, &thermal, &history);
    // The lid's own faults, the wide ones carved through the loop (F147, F153).
    let mut strain = strain::StrainProcess::new(graph, grid, &thermal);
    let mut floods = volcanism::FloodProcess::new(
        graph,
        grid,
        &thermal,
        plume.as_ref().map_or(&[], |p| p.rises.as_slice()),
    );
    let no_limbs = stress::Limbs::default();
    let mut episodic = episodic::EpisodicProcess::new(
        graph,
        grid,
        plume.as_ref().map_or(&[], |p| p.rises.as_slice()),
        plume.as_ref().map_or(&no_limbs, |p| &p.limbs),
    );
    // The strain the lid the last overturn keeps took from the mantle's flow since the
    // overturn before (F187), which the tesserae's fabric records.
    let mut lid_strain = episodic.as_ref().and_then(|e| {
        let to = e.event_time();
        let from = thermal
            .mantle
            .as_ref()?
            .overturns
            .iter()
            .rev()
            .find(|&&o| o < to - 1e-9)
            .copied()
            .unwrap_or(0.0);
        Some(crust::StrainRecord::new(from, to))
    });
    let mut io = io::IoProcess::new(graph, grid);
    // Io's lava erupts as floods, from its paterae (F117).
    let mut io_floods = io
        .as_ref()
        .map(|_| volcanism::FloodProcess::for_io(graph, grid, &thermal));
    let ice_shell = ice_shell_history(graph, &thermal);
    let rock_crust = rock_crust_history(graph, &thermal);
    // The relief the ice shell holds, band by band, from step to step (F134).
    let mut relax_bands = relaxation::Bands::default();
    let mut ice_shell_process = ice_shell::IceShellProcess::new(graph, grid);
    // The renewals keep the base relief's texture under the new terrain (F141), as it
    // stands: the relief the body formed with, relaxed as the ground relaxes and levelled
    // by its craters' ejecta as the ground is (F196), its own bands held as the surface's.
    let mut base_bands = relaxation::Bands::default();
    if ice_shell_process.is_some() {
        layers.base = Some(layers.bedrock.clone());
    }
    // The shields grow from the lava their rises' eruptions bring them, scheduled against
    // the ground the loop begins on; the floods lose what they take (F164).
    let taken: Vec<f64> = match (&mut plume, &floods) {
        (Some(p), Some(f)) => p.schedule(f, &layers, steps, dt),
        _ => Vec::new(),
    };
    let mut keyframes = Vec::with_capacity(MAX_KEYFRAMES + 1);
    keyframes.push(Keyframe {
        step: 0,
        time_gyr: 0.0,
        layers: layers.clone(),
    });
    // The lava laid by the last load the field holds, per coarse cell, and when the field is
    // read next (F157); and how far the ground sank under each load, when (F162).
    let mut laid = vec![0.0f32; grid.cell_count()];
    let mut next_epoch = LOAD_EPOCH_GYR;
    let mut flexure: Vec<(f64, Vec<f32>)> = Vec::new();
    // The magma that stays below (F174): what has intruded since the last epoch, all of it,
    // and each epoch's, for the finish.
    let mut pending = vec![0.0f32; grid.cell_count()];
    let mut intruded = vec![0.0f32; grid.cell_count()];
    let mut intrusions: Vec<(f64, Vec<f32>)> = Vec::new();
    // The intrusions the last load held.
    let mut intruded_laid = vec![0.0f32; grid.cell_count()];

    for s in 1..=steps {
        if sink.cancelled() {
            return Err(RunError::Cancelled);
        }
        let t0 = dt * (s - 1) as f64;
        let t1 = dt * s as f64;
        // The mantle's flow holds the ground up over its upwellings and draws it down over
        // its downwellings, as its pattern stands at the step's end (F173).
        if let Some(sup) = &support {
            let now: Vec<f32> = (0..grid.cell_count())
                .into_par_iter()
                .map(|k| sup.at(t1, grid.center(grid.cell(k))) as f32)
                .collect();
            let bedrock = layers.bedrock.as_mut_slice();
            for ((b, l), v) in bedrock.iter_mut().zip(lifted.iter_mut()).zip(&now) {
                *b += *v - *l;
                *l = *v;
            }
        }
        // An ice shell relaxes what it carries, the longest wavelengths first, with the
        // shell as the thermal history has it now; the base relief kept under the
        // renewals relaxes with it, and the renewals take it as it then stands (F196).
        if let Some(shell_at) = &ice_shell {
            if let (Some(p), Some(base)) = (&mut ice_shell_process, layers.base.as_mut()) {
                relaxation::relax(
                    &mut base_bands,
                    base,
                    graph.f64("body.radius"),
                    dt,
                    &shell_at(t1),
                    &relaxation::Local::default(),
                );
                p.set_base(base.as_slice());
            }
            let local = ice_shell_process
                .as_mut()
                .map_or(relaxation::Local::default(), |p| p.relaxation(t0, t1));
            let floor = relaxation::relax(
                &mut relax_bands,
                &mut layers.bedrock,
                graph.f64("body.radius"),
                dt,
                &shell_at(t1),
                &local,
            );
            if let (Some(p), Some(f)) = (&mut ice_shell_process, floor) {
                p.set_floor(f);
            }
        } else if let Some(crust_at) = &rock_crust {
            // The mantle's flow carries the crust beneath the stiff lid toward its
            // downwellings (F184), before the lower crust's own flow spreads what piles up.
            if let (Some(sup), Some(mantle)) = (&support, &thermal.mantle) {
                let moved = crust::drag(
                    grid,
                    layers.crust_thickness.as_mut_slice(),
                    layers.bedrock.as_mut_slice(),
                    sup,
                    mantle,
                    t0,
                    t1,
                    lid_strain.as_mut(),
                );
                if moved.carried > 0 && std::env::var_os("EXOGEN_DEBUG").is_some() {
                    eprintln!(
                        "crust {t0:.3} to {t1:.3}: the flow carried crust in {} cell-substeps, by up to {:.0} m, the thickest {:.1} km, the mean {:.1} km, {:.2} km of it foundered",
                        moved.carried,
                        moved.most,
                        moved.thickest / 1e3,
                        moved.mean / 1e3,
                        moved.foundered / 1e3
                    );
                }
            }
            // A rocky lid's lower crust flows where it is warm, the relief relaxing the
            // longest wavelengths last (F170), each cell under its own crust, the lid as it
            // stood through the step (F184).
            let crust: Vec<f32> = layers.crust_thickness.as_slice().to_vec();
            let thickest = crust.iter().fold(0.0f32, |m, &c| m.max(c)) as f64;
            let shell = crust_at(t1, thickest);
            let substeps: Vec<(f64, f64)> = (0..RELAX_SUBSTEPS)
                .map(|i| {
                    let sub = (t1 - t0) / RELAX_SUBSTEPS as f64;
                    (t0 + (i as f64 + 0.5) * sub, sub)
                })
                .collect();
            // The rises the plumes hold up are not the crust's to relax: the relief relaxes
            // toward their domes.
            let body_radius = graph.f64("body.radius");
            let domes: Option<Vec<f32>> = plume.as_ref().filter(|p| !p.rises.is_empty()).map(|p| {
                (0..grid.cell_count())
                    .into_par_iter()
                    .map(|k| {
                        let d = grid.center(grid.cell(k));
                        p.rises
                            .iter()
                            .map(|rise| {
                                let x = angle_between(rise.center, d) * body_radius / rise.radius;
                                rise.height * plume::dome(x) * plume::rise_progress(rise, t1)
                            })
                            .sum::<f64>() as f32
                    })
                    .collect()
            });
            // Nor is the ground the mantle's flow holds up (F173).
            let equilibrium = match domes {
                Some(mut d) if support.is_some() => {
                    for (x, l) in d.iter_mut().zip(&lifted) {
                        *x += l;
                    }
                    Some(d)
                }
                None if support.is_some() => Some(lifted.clone()),
                d => d,
            };
            let before: Vec<f32> = layers.bedrock.as_slice().to_vec();
            // Within the lid the tesserae keep, the flow that relaxes the relief is kept for
            // the strain it puts on the lid (F188).
            let collapsing = lid_strain.as_ref().is_some_and(|r| r.covers(t0, t1));
            if collapsing {
                relax_bands.flow = Some(vec![0.0; grid.cell_count()]);
            }
            relaxation::relax(
                &mut relax_bands,
                &mut layers.bedrock,
                body_radius,
                dt,
                &shell,
                &relaxation::Local {
                    equilibrium: equilibrium.as_deref(),
                    rock: Some(relaxation::RockCrust {
                        crust: &crust,
                        shell_at: &**crust_at,
                        steps: &substeps,
                    }),
                    ..Default::default()
                },
            );
            // The crust keeps the relief compensated: what the flow took from under a high
            // it moved to under a low (F170, F184).
            let airy = ((DENSITY_ROCK - DENSITY_CRUST_ROCK) / DENSITY_ROCK) as f32;
            if let (Some(record), Some(flow)) = (lid_strain.as_mut(), relax_bands.flow.take()) {
                let relaxed: Vec<f32> = layers
                    .bedrock
                    .as_slice()
                    .iter()
                    .zip(&before)
                    .map(|(h, b)| h - b)
                    .collect();
                record.add_collapse(grid, body_radius, &flow, &relaxed, &crust, t1);
            }
            let thickness = layers.crust_thickness.as_mut_slice();
            for ((c, h), b) in thickness.iter_mut().zip(layers.bedrock.as_slice()).zip(&before) {
                *c += (h - b) / airy;
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
        if let Some(s) = &mut strain {
            s.step(&mut layers, grid, graph.f64("body.radius"), t0, t1);
        }
        let mut renewed = false;
        let held = held_up(
            plume.as_ref().map(|p| p.load.as_slice()),
            support.as_ref().map(|_| lifted.as_slice()),
        );
        if let Some(e) = &mut episodic
            && let Some(fill) = e.step(&mut layers, t0, t1, held.as_deref())
            && let Some(f) = &mut floods
        {
            f.note_resurfacing(&fill, e.event_time());
            renewed = true;
        }
        if let Some(i) = &mut io {
            let lay = s % keep_every == 0 || s == steps;
            i.step(&mut layers, t0, t1, lay);
            // Io's lava in the surface the step keeps: its paterae's eruptions of the last
            // few centuries of millennia, floods as on any body (F117).
            if lay && let Some(f) = &mut io_floods {
                let rate = f.extrusion_rate(t1);
                let eruptions = io::eruptions(graph, &i.paterae, &layers.bedrock, t1, rate);
                if s == steps && std::env::var_os("EXOGEN_DEBUG").is_some() {
                    eprintln!(
                        "io eruptions: {:.1} km³ a year, {} fields, {:.0} km³",
                        rate * 1e-9,
                        eruptions.len(),
                        eruptions.iter().map(|e| e.volume).sum::<f64>() * 1e-9
                    );
                }
                for (n, e) in eruptions.into_iter().enumerate() {
                    f.apply(
                        &mut layers,
                        volcanism::FloodEvent {
                            vent: grid.index(grid.locate(e.vent).0),
                            volume: e.volume,
                            time_gyr: t1 - e.age_yr * 1e-9,
                            id: ((s as u64) << 32) | n as u64,
                            effusion: e.effusion,
                            span_yr: e.span_yr,
                            composition: None,
                        },
                    );
                }
            }
            if lay {
                i.cut_paterae(&mut layers, t1);
            }
        }
        let boost = plume.as_ref().map(|p| p.vent_boost());
        let rise_vents = plume
            .as_mut()
            .map_or_else(Vec::new, |p| p.rise_vents(t0, t1));
        // Everything that happens in the step, in order of time: craters, scarps and
        // eruptions interleaved as they come, so what is older is cut by what is younger.
        let craters = impacts.sample(t0, t1);
        let scarps = contraction
            .as_mut()
            .map_or_else(Vec::new, |c| c.sample(t0, t1, &history));
        let mut eruptions = floods.as_mut().map_or_else(Vec::new, |f| {
            f.sample(
                &layers,
                boost.as_deref(),
                &rise_vents,
                t0,
                t1,
                taken.get(s - 1).copied(),
                held.as_deref(),
            )
        });
        // The lava the shields' flows sent past their feet, ponding beyond them (F163).
        if floods.is_some()
            && let Some(p) = &mut plume
        {
            let foot = p.foot_eruptions(t0, t1);
            if std::env::var_os("EXOGEN_DEBUG").is_some() && !foot.is_empty() {
                let volume: f64 = foot.iter().map(|e| e.volume).sum();
                eprintln!(
                    "shields' feet {t0:.3} to {t1:.3}: {} eruptions, {:.4e} km3",
                    foot.len(),
                    volume / 1e9
                );
            }
            eruptions.extend(foot);
        }
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
        // The melt that stays below thickens the crust beneath where it rose from: the
        // intrusive share of the step's melt, and the magma the full edifices could not
        // lift, beneath their rises (F174).
        if let Some(f) = &floods {
            f.intrude(f.step_intrusive(t0, t1), 0.5 * (t0 + t1), &mut pending);
            if let Some(p) = &plume {
                for (rise, volume) in p.surplus(t0, t1) {
                    f.intrude_rise(rise, volume, &mut pending);
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
        // The lid's field read anew (F157): every LOAD_EPOCH_GYR, at a renewal and at the
        // end, the ground as it stands and the lid as thick as it is then, and the rock lava
        // laid since the last load, once there is enough of it, bending that lid.
        if t1 >= next_epoch - 1e-9 || renewed || s == steps {
            while next_epoch <= t1 + 1e-9 {
                next_epoch += LOAD_EPOCH_GYR;
            }
            // The crust the intrusions thickened stands as much higher, its weight bending
            // the lid with the lava's (F174).
            if pending.iter().any(|&v| v > 0.0) {
                let bedrock = layers.bedrock.as_mut_slice();
                let crust = layers.crust_thickness.as_mut_slice();
                for k in 0..pending.len() {
                    bedrock[k] += pending[k];
                    crust[k] += pending[k];
                    intruded[k] += pending[k];
                }
                let n = pending.len();
                intrusions.push((t1, std::mem::replace(&mut pending, vec![0.0; n])));
            }
            let elastic = thermal.at(t1).elastic;
            let held = held_up(
                plume.as_ref().map(|p| p.load.as_slice()),
                support.as_ref().map(|_| lifted.as_slice()),
            );
            let relief =
                stress::Relief::read(&layers.bedrock, held.as_deref(), cell, elastic, gravity);
            let mut load = None;
            if let Some(f) = &floods
                && let Some(density) = f.magma_density()
            {
                // The lava the floods and the shields have laid (F161), less what the last
                // load held.
                let shields = plume.as_ref().map(|p| p.shield_built());
                let aureoles = plume.as_ref().map(|p| p.aureole_laid());
                let now: Vec<f32> = (0..grid.cell_count())
                    .map(|k| {
                        f.thickness[k]
                            + shields.map_or(0.0, |s| s[k])
                            + aureoles.map_or(0.0, |a| a[k])
                    })
                    .collect();
                let mut sunk_all: Option<Vec<f32>> = None;
                let mut stress_all: Option<load::LoadStress> = None;
                let mut datum = 0.0f64;
                let mean_of = |v: &[f32]| -> f64 {
                    let (mut sum, mut area) = (0.0, 0.0);
                    for (k, &x) in v.iter().enumerate() {
                        let a = grid.cell_area(grid.cell(k));
                        sum += a * x as f64;
                        area += a;
                    }
                    sum / area.max(1e-30)
                };
                let increment: Vec<f32> = now
                    .iter()
                    .zip(&laid)
                    .map(|(t, l)| (t - l).max(0.0))
                    .collect();
                let most = increment.iter().fold(0.0f32, |m, &v| m.max(v));
                if most as f64 >= LOAD_MIN_INCREMENT
                    && let Some(laid_load) = load::lay(
                        grid,
                        &increment,
                        body_radius,
                        gravity,
                        elastic,
                        density,
                        shields,
                    )
                {
                    // The ground bends under the load (F162): it sinks, but beneath a shield
                    // the edifice's lava fills the sag as far as it is thick, so that its
                    // surface stays as it built it and only the ground about it sinks into a
                    // moat, rising in a bulge beyond.
                    let w = laid_load.deflection_on(grid);
                    sunk_all = Some(load::sink(&w, shields));
                    stress_all = Some(laid_load.stress);
                    datum += mean_of(&increment);
                    laid.copy_from_slice(&now);
                }
                // The magma that stayed below weighs on the lid as well (F174), a load of its
                // own: the edifices standing on it sink with the ground, their lava filling
                // only the sags their own weight makes.
                let deep: Vec<f32> = intruded
                    .iter()
                    .zip(&intruded_laid)
                    .map(|(t, l)| (t - l).max(0.0))
                    .collect();
                let most = deep.iter().fold(0.0f32, |m, &v| m.max(v));
                if most as f64 >= LOAD_MIN_INCREMENT
                    && let Some(laid_load) =
                        load::lay(grid, &deep, body_radius, gravity, elastic, density, None)
                {
                    let w = laid_load.deflection_on(grid);
                    match &mut sunk_all {
                        Some(s) => s.iter_mut().zip(&w).for_each(|(a, b)| *a += b),
                        None => sunk_all = Some(w),
                    }
                    match &mut stress_all {
                        Some(s) => s.add(&laid_load.stress),
                        None => stress_all = Some(laid_load.stress),
                    }
                    datum += mean_of(&deep);
                    intruded_laid.copy_from_slice(&intruded);
                }
                // What the loads laid evenly over the whole body bends nothing but makes the
                // body that much larger, its mean radius, the heights' datum, with it (F177).
                if let Some(sunk) = &mut sunk_all {
                    for s in sunk.iter_mut() {
                        *s += datum as f32;
                    }
                }
                if let Some(sunk) = sunk_all {
                    let bedrock = layers.bedrock.as_mut_slice();
                    for (b, s) in bedrock.iter_mut().zip(&sunk) {
                        *b -= s;
                    }
                    flexure.push((t1, sunk));
                }
                load = stress_all;
            }
            history.push(stress::FieldEpoch {
                time: t1,
                elastic,
                relief,
                load,
                flow: flow_at(t1, elastic),
            });
        }
        // The faults the field opened or grew in the step (F178, F179): a rift wide enough
        // for the loop's grid is carved with the rise's canyons as it grows, the rest drawn
        // by the finish.
        if let Some(r) = &mut rifts {
            for f in r.step(&history, t1) {
                if let Some(p) = &mut plume {
                    match p.canyons.iter_mut().find(|c| c.id == f.id) {
                        Some(c) => *c = f,
                        None => p.canyons.push(f),
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
    let (scarps, thrust_ridges) = match contraction {
        Some(c) => (c.records, c.ridges),
        None => (Vec::new(), Vec::new()),
    };
    if std::env::var_os("EXOGEN_DEBUG").is_some() && !scarps.is_empty() {
        eprintln!("thrusts: {}", tectonics::strike_summary(&scarps));
    }
    let contraction_records = &scarps;
    let (paterae, mountains, io_base) = match io {
        Some(i) => (i.paterae, i.mountains, i.base.map(|b| b.bedrock)),
        None => (Vec::new(), Vec::new(), None),
    };
    let (mut resurfacing, coronae, domes, novae) = match episodic {
        Some(e) => (e.resurfacing, e.coronae, e.domes, e.novae),
        None => (None, Vec::new(), Vec::new(), Vec::new()),
    };
    // The tesserae keep the strain their lid took, each set with the lid of its time.
    if let (Some(r), Some(record)) = (resurfacing.as_mut(), lid_strain) {
        r.strain = record
            .sets
            .into_iter()
            .chain(record.collapse)
            .filter(|s| s.strain.iter().any(|v| v.iter().any(|&x| x != 0.0)))
            .map(|s| episodic::LidStrain {
                time: s.time,
                elastic: thermal.at(s.time).elastic,
                strain: s.strain,
                pressed: s.pressed,
                pulled: s.pulled,
            })
            .collect();
    }
    let rift_graben = rifts.map_or_else(Vec::new, |r| r.finish(plume.is_some()));
    let (rises, shields, mut canyons, mut graben) = match plume {
        Some(p) => (p.rises, p.shields, p.canyons, p.graben),
        None => (Vec::new(), Vec::new(), Vec::new(), Vec::new()),
    };
    // The novae's fractures are stamped as the graben are, and so are the faults the lid's
    // growth opens (F147): the wide ones carved through the loop as canyons (F153).
    graben.extend(novae);
    if let Some(s) = strain {
        canyons.extend(s.canyons);
        graben.extend(s.graben);
    }
    graben.extend(rift_graben);
    let (floods, mut ridges, provinces) = match floods {
        Some(f) => {
            // The plains' ridges stand where the lid squeezes them, read when each place's
            // would stand, and follow the lava's thickness where the field leaves them free
            // (F156).
            let (squeezed, squeeze) =
                ridges::plains_squeeze(graph, &thermal, &history, grid, &f.last_time);
            let ridges = ridges::plan_ridges(
                &grid,
                graph.f64("body.radius"),
                layers.unit_id.as_slice(),
                &f.last_time,
                &squeezed,
                &squeeze,
                &thrust_ridges,
                age,
                graph.f64("seed.volcanism") as u64,
            );
            // Where the ridges and thrusts stand against the mantle's flow (F185), for the
            // debug log: their length over its downwellings and upwellings against the area.
            if let (Some(sup), true) = (&support, std::env::var_os("EXOGEN_DEBUG").is_some()) {
                // Each against the pattern of its own time.
                let class = |d: [f64; 3], t: f64| {
                    let p = sup.pattern(t, d);
                    if p < -1.0 { 0 } else if p > 1.0 { 2 } else { 1 }
                };
                let mut area = [0.0f64; 3];
                for k in 0..grid.cell_count() {
                    let c = grid.cell(k);
                    area[class(grid.center(c), age)] += grid.cell_area(c);
                }
                let whole: f64 = area.iter().sum();
                let tally = |traces: &mut dyn Iterator<Item = (&Vec<[f64; 3]>, f64, f64)>| {
                    let mut l = [0.0f64; 3];
                    for (points, spacing, t) in traces {
                        for &q in points.iter() {
                            l[class(q, t)] += spacing;
                        }
                    }
                    let all: f64 = l.iter().sum::<f64>().max(1.0);
                    format!(
                        "{:.0} km: down {:.2}, between {:.2}, up {:.2}",
                        all / 1e3,
                        l[0] / all,
                        l[1] / all,
                        l[2] / all
                    )
                };
                eprintln!(
                    "flow classes by area: down {:.2}, between {:.2}, up {:.2}; plains ridges {}; thrusts {}",
                    area[0] / whole,
                    area[1] / whole,
                    area[2] / whole,
                    tally(&mut ridges.iter().map(|r| (&r.points, r.spacing, r.time_gyr))),
                    tally(
                        &mut contraction_records
                            .iter()
                            .map(|r| (&r.points, r.spacing, r.time_gyr))
                    )
                );
            }
            (f.records, ridges, f.provinces)
        }
        None => (
            io_floods.map_or_else(Vec::new, |f| f.records),
            Vec::new(),
            Vec::new(),
        ),
    };
    // The thrusts beneath the plains stand among the plains' own ridges (F155).
    ridges.extend(thrust_ridges);
    ridges.sort_by(|a, b| a.time_gyr.total_cmp(&b.time_gyr).then(a.id.cmp(&b.id)));
    let (lanes, lane_coverage, lane_core, episodes, polar_start, ice_coronae, renewed_fraction) =
        match ice_shell_process {
            Some(p) => (
                p.lanes,
                p.lane_coverage,
                p.lane_core.map(|v| ice_shell::field_from(grid, &v)),
                p.episodes,
                p.present_start,
                p.coronae,
                p.renewed_fraction,
            ),
            None => (Vec::new(), 0.0, None, Vec::new(), 0.0, Vec::new(), 0.0),
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
        support,
        rises,
        shields,
        canyons,
        graben,
        resurfacing,
        coronae,
        domes,
        flexure,
        intrusions,
        paterae,
        mountains,
        io_base,
        lanes,
        lane_coverage,
        lane_core,
        episodes,
        polar_start,
        ice_coronae,
        renewed_fraction,
        hash,
    })
}

/// The ground the plumes' domes and the mantle's flow hold up (F150, F173), which the lid
/// does not bear: their sum where either stands.
fn held_up(load: Option<&[f32]>, lifted: Option<&[f32]>) -> Option<Vec<f32>> {
    match (load, lifted) {
        (None, None) => None,
        (Some(a), None) => Some(a.to_vec()),
        (None, Some(b)) => Some(b.to_vec()),
        (Some(a), Some(b)) => Some(a.iter().zip(b).map(|(x, y)| x + y).collect()),
    }
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
    // The base relief as the loop left it (F196), for the floors laid here; the finish's
    // craters, laid after, have no renewal to keep it for.
    let kept_base = layers.base.take();
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
    // The coronae over an ice mantle's upwellings formed by the epoch (F133).
    let ice_coronae: Vec<corona::Corona> = loop_output
        .ice_coronae
        .iter()
        .chain(&loop_output.coronae)
        .filter(|c| c.time_gyr <= present + 1e-9)
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
    // Each terrain's margin belt (F194), under the lid of its time.
    let belt_of = |map: &ice::ShellMap, time: f64| {
        ice_shell::MarginBelt::new(
            map,
            loop_output.thermal.at(time).elastic,
            graph.f64("body.gravity"),
            body_radius,
        )
    };
    let fossil_belts: Vec<ice_shell::MarginBelt> = episode_maps
        .iter()
        .zip(&episodes)
        .map(|(map, e)| belt_of(map, e.end))
        .collect();
    let polar_belt = belt_of(&shell_map, polar_time);
    let fossil_edges: Vec<Vec<f32>> = episode_maps
        .iter()
        .zip(&fossil_belts)
        .map(|(map, belt)| ice_shell::polar_edge(&grid, body_radius, map, &lineae, halo, belt))
        .collect();
    let polar_edge: Vec<f32> = if polar && shell_map.pole.is_some() {
        ice_shell::polar_edge(&grid, body_radius, &shell_map, &lineae, halo, &polar_belt)
    } else {
        Vec::new()
    };
    let polar = !polar_edge.is_empty();
    {
        let lane_hits =
            ice_shell::lane_hits(&grid, body_radius, cell, &lanes, loose_before.as_deref());
        let corona_hits = corona::corona_hits(&grid, body_radius, cell, &ice_coronae);
        // Each terrain with its time: a lane (0), a past episode's terrain (1), the
        // present polar terrain (2), a corona (3).
        let mut order: Vec<(f64, u8, usize)> = (0..lanes.len())
            .map(|i| (lanes[i].time_gyr, 0u8, i))
            .collect();
        order.extend((0..episodes.len()).map(|k| (episodes[k].end, 1u8, k)));
        if polar {
            order.push((polar_time, 2u8, 0));
        }
        order.extend((0..ice_coronae.len()).map(|k| (ice_coronae[k].time_gyr, 3u8, k)));
        order.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)).then(a.2.cmp(&b.2)));
        for (time, kind, i) in order {
            match kind {
                0 => renewal.take_lane(&lanes[i], &lane_hits[i], body_radius, cell),
                1 => renewal.take_polar(&fossil_edges[i], time),
                2 => renewal.take_polar(&polar_edge, time),
                _ => renewal.take_hits(
                    time,
                    if ice_coronae[i].head.is_some() {
                        crate::units::UNIT_LAVA
                    } else {
                        crate::units::UNIT_CORONA
                    },
                    &corona_hits[i],
                    ice_coronae[i].smooth_length(),
                    ice_coronae[i].faults_on(cell),
                ),
            }
        }
    }
    // The lanes' and the coronae's floor where the loop's cells did not reach their
    // edges: what the loop renewed of each cell, resampled and smoothed over its cell so
    // that no facet of the coarse grid shows in a lane, and the rest taken here.
    // The base relief under the renewed terrains (F141): kept below the length each
    // levels over, and, where the loop laid the floor, its octaves below the loop's
    // cells, which that grid could not hold.
    let crust = primordial_crust::PrimordialCrust::new(graph, loop_output.grid);
    let loop_cell = loop_output.grid.delta() * body_radius;
    let fine_below_loop: Option<Vec<f32>> = (polar || !episode_maps.is_empty()).then(|| {
        let folds = POLAR_FOLD_SPACING.max(GROOVE_MIN_CELLS * cell);
        crust.below_loop(&grid, folds, body_radius)
    });
    if !lanes.is_empty() || !ice_coronae.is_empty() {
        let core: Option<Vec<f32>> = loop_output.lane_core.as_ref().map(|c| {
            let fine = c.resample(grid);
            relaxation::low_pass(&fine, cell, loop_cell)
        });
        ice_shell::lay_lane_floor(
            &mut layers,
            body_radius,
            core.as_deref(),
            &renewal,
            Some(&crust),
            kept_base.as_ref(),
        );
    }
    renewal.apply(&mut layers);
    let renewed: &[f32] = &renewal.time;
    drop(loose_before);
    let loose_cover = cover_on(&layers);

    // Loop craters that exist at this epoch. Their coarse form is in the resampled base;
    // their resolved detail is restored below, in age order among the replayed craters. A
    // mega-impact's form is the crust and the ground it holds up, all in the base, and it
    // is not restored (F168).
    let (megas, loop_craters): (Vec<CraterRecord>, Vec<CraterRecord>) = loop_output
        .craters
        .iter()
        .filter(|c| c.time_gyr <= present + 1e-9)
        .cloned()
        .partition(|c| mega::is_mega(c.diameter, body_radius));
    // The loop craters a later mega-impact's ejecta renewed are gone: their restore would
    // bring them back over it.
    let loop_craters: Vec<CraterRecord> = loop_craters
        .into_iter()
        .filter(|c| {
            !megas
                .iter()
                .any(|m| m.time_gyr > c.time_gyr && mega::renews(m, c.center, body_radius))
        })
        .collect();
    let scarps: Vec<tectonics::ScarpRecord> = loop_output
        .scarps
        .iter()
        .filter(|s| s.time_gyr <= present + 1e-9)
        .cloned()
        .collect();
    // On Io the floods are the epoch's alone: each kept step erupts its own over the
    // surface it lays afresh (F117).
    let io_regime = io::is_io_regime(graph);
    let floods: Vec<volcanism::FloodRecord> = loop_output
        .floods
        .iter()
        .filter(|f| f.time_gyr <= present + 1e-9 && (!io_regime || present - f.time_gyr < 1e-3))
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
    // The faults inside the wide troughs (F152): the floors of the canyons and graben the
    // grid resolves break into tilted blocks, drawn as graben just after their trough.
    let graben: Vec<plume::Trough> = {
        let mut all = graben;
        let wide: Vec<&plume::Trough> = canyons.iter().chain(all.iter()).collect();
        let cell = grid.delta() * body_radius;
        let interior = plume::interior_faults(
            &wide,
            cell,
            body_radius,
            graph.f64("seed.tectonics") as u64 ^ 0x1B10,
        );
        all.extend(interior);
        all
    };
    let resurfacing = loop_output
        .resurfacing
        .clone()
        .filter(|r| r.time_gyr <= present + 1e-9);
    // A rock lid's coronae over its mantle's diapirs (F145) take the corona process's
    // path with the icy ones.
    let rock_coronae: Vec<corona::Corona> = loop_output
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
    // The floods' coarse lava, each flood's smeared over a coarse cell at its margin, for its
    // pass to take up where it goes.
    let mut smear = volcanism::LavaSmear::new(
        loop_output.grid,
        grid,
        &floods,
        resurfacing.as_ref().map(|r| r.fill.as_slice()),
    );
    // The Io regime's mountains and paterae are drawn whole at this resolution on the
    // buried plains the loop laid its coarse versions on (F87, F88), smoothed over the
    // loop's cell as the rest of the bedrock is, so nothing of the coarse ones is left.
    if let Some(io_base) = &loop_output.io_base {
        let coarse_cell_m = loop_output.grid.delta() * body_radius;
        // The floods' coarse lava lies on the smoothed plains just as their passes take it
        // up (F117), so that the ground each finds is the plains: smoothed with them, the
        // lava had left its coarse cells' facets where a pass took it up, and on the flat
        // plains its lakes' edges ran in steps along them.
        let smoothed = relaxation::low_pass(&io_base.resample(grid), cell, coarse_cell_m);
        for ((b, s), l) in layers
            .bedrock
            .as_mut_slice()
            .iter_mut()
            .zip(&smoothed)
            .zip(&smear.total)
        {
            *b = s + l;
        }
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
    // The floods that show at the end are drawn as their flows (F94), traced over the base
    // smoothed over a front's width.
    // Io's floods are all young, years apart, and its successors' lobes leave gaps a
    // buried flood shows through: every one is drawn as its flows (F117).
    let flood_visible = if io_regime {
        vec![true; floods.len()]
    } else {
        volcanism::visible_floods(&floods)
    };
    let flood_ground = (!floods.is_empty() || !shields.is_empty()).then(|| {
        let smoothed = relaxation::low_pass(&layers.bedrock, cell, FLOOD_GROUND_SMOOTH);
        crate::grid::projection::PaddedField::new(&grid, &smoothed)
    });
    let flood_seed = graph.f64("seed.volcanism") as u64;
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
        &loop_output.flexure,
        loop_output.support.as_ref(),
        &loop_output.intrusions,
    );
    if !floods.is_empty() || resurfacing.is_some() {
        with_stamper(&grid, |stamper| {
            // Inside its rim: beyond, its restore buries the passes' relief by the
            // ejecta's thickness, so that the floods fade under its ejecta (M15).
            for c in &loop_craters {
                let dim = dimensions(c.diameter, &inputs.scaling);
                let radius = 1.05 * dim.radius * c.ellipticity.sqrt();
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
        bombardment: inputs.bombardment.clone(),
        body_radius,
        d_min,
        d_max: loop_output.loop_threshold,
        bins_per_octave: BINS_PER_OCTAVE,
        saturation_fraction: inputs.saturation_fraction,
        retention_age: inputs.retention_age,
        oblique: inputs.oblique,
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
    // The secondaries the loop's craters and the finish's own throw (F129), each with the
    // crater it came from, whose rays draw its fresh ground where it lies (F128).
    let secondaries = if std::env::var_os("EXOGEN_SKIP_SECONDARIES").is_some() {
        Vec::new()
    } else {
        let primaries: Vec<Primary> = loop_craters
            .iter()
            .map(|c| Primary {
                id: c.id,
                center: c.center,
                diameter: c.diameter,
                time_gyr: c.time_gyr,
                azimuth: c.azimuth,
                angle: c.angle,
            })
            .chain(
                population
                    .iter()
                    .filter(|c| throws_visible(c.diameter, d_min))
                    .map(|c| Primary {
                        id: c.key,
                        center: c.center,
                        diameter: c.diameter,
                        time_gyr: c.time_gyr,
                        azimuth: c.azimuth,
                        angle: c.angle,
                    }),
            )
            .collect();
        sample_secondaries(
            &primaries,
            &inputs,
            d_min,
            present,
            grid.delta() * body_radius,
        )
    };
    let mut thrown: HashMap<u64, Vec<crate::crater::CraterSample>> = HashMap::new();
    for (c, parent) in &secondaries {
        thrown.entry(*parent).or_default().push(c.clone());
    }
    let replay: Vec<(crate::crater::CraterSample, bool)> = population
        .into_iter()
        .map(|c| (c, false))
        .chain(secondaries.into_iter().map(|(c, _)| (c, true)))
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
                    // The share of a graben no later crater as wide has struck out: the
                    // equilibrium the craters themselves come to, at the graben's width
                    // (F97).
                    let survives = |w: f64| {
                        let d_km = w / 1e3;
                        let later = inputs.bombardment.density_above(d_km, time, present);
                        let equilibrium = inputs.saturation_fraction
                            * crate::crater::production::geometric_saturation(d_km);
                        (-later / equilibrium.max(1e-30)).exp()
                    };
                    rings::ring_system(
                        &rings::Basin {
                            center,
                            diameter,
                            time_gyr: time,
                            seed,
                            lithosphere,
                            inflow,
                            kept: &kept,
                            survives: &survives,
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
        Dome(usize),
        Patera(usize),
        Mountain(usize),
        Linea(usize),
        Groove(usize),
        IceCorona(usize),
        PolarTerrain,
        FossilTerrain(usize),
        Chaos(usize),
        Fabric(usize),
        Cryodome(usize),
        Mega(usize),
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
            + ice_coronae.len()
            + 1
            + chaos.len()
            + fabric.len()
            + cryodomes.len()
            + replay.len(),
    );
    // A shield is drawn when the flows that show on it begin (F101), over the craters they
    // bury, or just after a resurfacing they outlast, which buries the ones before it.
    let shield_times: Vec<f64> = shields
        .iter()
        .map(|sh| {
            let start = plume::shield_window(sh, present).map_or(sh.start_gyr, |w| w.0);
            match &resurfacing {
                Some(r) if r.time_gyr >= start && r.time_gyr < present => r.time_gyr + 1e-6,
                _ => start,
            }
        })
        .collect();
    for (k, sh) in shields.iter().enumerate() {
        events.push((
            shield_times[k],
            0,
            sh.id ^ 0x5B1E_0000_0000,
            Event::Shield(k),
        ));
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
    for (k, d) in domes.iter().enumerate() {
        events.push((d.time_gyr, 0, d.id ^ 0xD03E_0000_0000, Event::Dome(k)));
    }
    for (k, p) in paterae.iter().enumerate() {
        // After the floods (F117).
        events.push((present, 2, p.id, Event::Patera(k)));
    }
    for (k, m) in mountains.iter().enumerate() {
        events.push((present, 0, m.id, Event::Mountain(k)));
    }
    for (k, l) in lineae.iter().enumerate() {
        events.push((l.time_gyr, 0, l.id ^ 0x11CE_0000_0000, Event::Linea(k)));
    }
    for (k, l) in lanes.iter().enumerate() {
        events.push((l.time_gyr, 0, l.id ^ 0x6700_0000_0000, Event::Groove(k)));
    }
    for (k, c) in ice_coronae.iter().enumerate() {
        events.push((c.time_gyr, 0, c.id ^ 0xC020_1CE0_0000, Event::IceCorona(k)));
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
        // Io's lava runs over the paterae and mountains the step laid (F117), in the order
        // it erupted.
        if io_regime {
            events.push((present, 1, f.id, Event::Flood(k)));
        } else {
            events.push((f.time_gyr, 0, f.id ^ 0x7F10_0000_0000, Event::Flood(k)));
        }
    }
    for (k, s) in scarps.iter().enumerate() {
        events.push((s.time_gyr, 0, s.id ^ 0x5CA2_0000_0000, Event::Scarp(k)));
    }
    for (i, c) in loop_craters.iter().enumerate() {
        events.push((c.time_gyr, 0, c.id, Event::Loop(i)));
    }
    for (i, m) in megas.iter().enumerate() {
        events.push((m.time_gyr, 0, m.id, Event::Mega(i)));
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
                Event::IceCorona(k) => {
                    let c = &ice_coronae[*k];
                    (
                        c.center,
                        c.cap_radius(body_radius, cell) + 2.0 * grid.delta(),
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
                // A mega-impact's ejecta reaches round the body.
                Event::Mega(_) => ([0.0, 0.0, 1.0], std::f64::consts::PI + 2.0 * grid.delta()),
            }
        })
        .collect();
    let waves = build_waves(&caps);

    let mut done = 0usize;
    let mut next_report = 0usize;
    // The craters' records by their events' order, whatever wave each came in.
    let mut stamped: Vec<(usize, CraterRecord)> = Vec::new();
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
                            * inputs.bombardment.equivalent_gyr(c.time_gyr, present)
                            * 1e3;
                        let degraded = ProfileTable::new(dim, kappa_t);
                        let smooth =
                            ProfileTable::new(dim, kappa_t + 0.5 * smooth_sigma * smooth_sigma);
                        // The crater as it stands, whose walls shed where steep (F98).
                        let mut walls = degraded.clone();
                        let mut detail = degraded;
                        // On an ice shell the detail relaxes by wavelength over the
                        // crater's life (F70), and its terraces, rings and peak keep the
                        // share of a short wavelength; the loop relaxed the rest.
                        let mut feature_weight = 1.0;
                        if relaxes(dim.radius, c.time_gyr, c.center) {
                            let kept = |l: f64| crater_kept(l, c.time_gyr, c.center);
                            let full = ProfileTable::relaxed(dim, kappa_t, &kept);
                            walls.values = full.values.clone();
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
                            shed: inputs.shed,
                        };
                        let mut out = restore_compute(
                            &layers,
                            &grid,
                            body_radius,
                            c.center,
                            &base,
                            &base_deep,
                            &fresh,
                            &detail,
                            &walls,
                            &params,
                            stamper,
                        );
                        // A shield that grew over it since buried it (F115).
                        plume::bury_under_shields(
                            &mut out,
                            &layers,
                            &grid,
                            &shields,
                            body_radius,
                            c.center,
                            fresh.dim.ejecta_radius * c.ellipticity.sqrt(),
                            fresh.dim.depth + fresh.dim.rim_height,
                            c.time_gyr,
                            present,
                        );
                        Some((out, None, None))
                    }
                    Event::Flood(k) => {
                        let flows =
                            flood_ground
                                .as_ref()
                                .filter(|_| flood_visible[*k])
                                .map(|ground| volcanism::FloodFlows {
                                    ground,
                                    body_radius,
                                    seed: flood_seed,
                                    hollow_reach: 0.5 * graph.f64("cratering.loop_threshold"),
                                    soft: io_regime,
                                });
                        let (mut out, taken, stats) = volcanism::flood_compute(
                            &layers,
                            &grid,
                            &loop_output.grid,
                            &floods[*k],
                            &smear,
                            &deform,
                            flows,
                            &base,
                            stamper,
                        );
                        if io_regime {
                            io::as_io_flow(
                                &mut out,
                                &layers,
                                &grid,
                                body_radius,
                                floods[*k].id,
                                (present - floods[*k].time_gyr) * 1e9,
                                floods[*k].span_yr,
                            );
                        }
                        Some((out, None, Some((taken, stats))))
                    }
                    Event::Shield(k) => {
                        let out = plume::shield_compute(
                            &layers,
                            &grid,
                            body_radius,
                            &shields[*k],
                            flood_ground.as_ref(),
                            present,
                            shield_times[*k],
                            stamper,
                        );
                        Some((out, None, None))
                    }
                    Event::Canyon(k) => {
                        let c = &canyons[*k];
                        let kappa_t = inputs.diffusivity_m2_per_myr
                            * inputs
                                .bombardment
                                .equivalent_gyr(c.end_gyr.min(present), present)
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
                            &shields,
                            present,
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
                    Event::IceCorona(k) => {
                        let out = corona::corona_compute(
                            &layers,
                            &grid,
                            body_radius,
                            &ice_coronae[*k],
                            cell,
                            present,
                            renewed,
                            &base,
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
                            fine_below_loop.as_deref(),
                            &polar_belt,
                            stamper,
                        );
                        Some((out, None, None))
                    }
                    Event::FossilTerrain(k) => {
                        // Gardening since the episode spreads and lowers its structures,
                        // as it degrades the craters of the time.
                        let kappa_t = inputs.diffusivity_m2_per_myr
                            * inputs.bombardment.equivalent_gyr(episodes[*k].end, present)
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
                            fine_below_loop.as_deref(),
                            &fossil_belts[*k],
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
                    Event::Mega(i) => {
                        let floor = inputs.mega_floor_crust;
                        let out = mega::renew_compute(
                            &layers,
                            &grid,
                            body_radius,
                            &megas[*i],
                            &base,
                            floor / MEGA_FLOOR_CRUST_SHARE - floor,
                        );
                        Some((out, None, None))
                    }
                    Event::Linea(k) => {
                        let l = &lineae[*k];
                        let kappa_t = inputs.diffusivity_m2_per_myr
                            * inputs.bombardment.equivalent_gyr(l.time_gyr, present)
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
                            * inputs.bombardment.equivalent_gyr(r.time_gyr, present)
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
                            * inputs.bombardment.equivalent_gyr(g.start_gyr, present)
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
                            &shields,
                            present,
                            stamper,
                        );
                        Some((out, None, None))
                    }
                    Event::Ridge(k) => {
                        let r = &ridges[*k];
                        let kappa_t = inputs.diffusivity_m2_per_myr
                            * inputs.bombardment.equivalent_gyr(r.time_gyr, present)
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
                        // degraded width minus the coarse-smoothed one; a scarp the loop
                        // did not carry, whole (F155).
                        let s = &scarps[*k];
                        let kappa_t = inputs.diffusivity_m2_per_myr
                            * inputs.bombardment.equivalent_gyr(s.time_gyr, present)
                            * 1e3;
                        let fine = (4.0 * kappa_t).sqrt();
                        let coarse = (4.0 * kappa_t + 2.0 * smooth_sigma * smooth_sigma).sqrt();
                        let out = tectonics::scarp_compute(
                            &layers,
                            &grid,
                            body_radius,
                            s,
                            fine,
                            s.carried.then_some(coarse),
                            stamper,
                        );
                        Some((out, None, None))
                    }
                    Event::Replay(j) => {
                        let (c, secondary) = &replay[*j];
                        // A secondary carries its parent's time: where the ground was
                        // renewed since, by ice from below, lava or a younger crater, it is
                        // gone with the parent, as the population's craters are (F129).
                        if *secondary {
                            let k = grid.index(grid.locate(c.center).0);
                            if layers.surface_age.as_slice()[k] as f64 > c.time_gyr + 1e-6 {
                                return None;
                            }
                        }
                        let dim = dimensions(c.diameter, &inputs.scaling);
                        // Degradation by the gardening that followed the crater: the
                        // diffusivity at the present flux times the equivalent flux since.
                        let kappa_t = inputs.diffusivity_m2_per_myr
                            * inputs.bombardment.equivalent_gyr(c.time_gyr, present)
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
                            levelling: inputs.levelling,
                            seed: c.key,
                            deep,
                            maturity_timescale: inputs.maturity_timescale,
                            renewal_guard: false,
                            terraced: true,
                            shed: inputs.shed,
                        };
                        let mut out = stamp_compute(
                            &layers,
                            &grid,
                            body_radius,
                            c.center,
                            &table,
                            &params,
                            stamper,
                        );
                        // A shield that grew over it since buried it (F115).
                        plume::bury_under_shields(
                            &mut out,
                            &layers,
                            &grid,
                            &shields,
                            body_radius,
                            c.center,
                            dim.ejecta_radius * c.ellipticity.sqrt(),
                            dim.depth + dim.rim_height,
                            c.time_gyr,
                            present,
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
        for (&ei, output) in wave.iter().zip(outputs) {
            let Some((out, record, taken)) = output else {
                continue;
            };
            if let Some((t, stats)) = &taken {
                smear.note(t, &out, *stats);
            }
            out.apply(&mut layers);
            if let Some(r) = record {
                stamped.push((ei, r));
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
    stamped.sort_by_key(|&(ei, _)| ei);
    records.extend(stamped.into_iter().map(|(_, r)| r));
    if std::env::var_os("EXOGEN_DEBUG").is_some() {
        eprintln!("finish: {} events in {} waves", events.len(), waves.len());
    }
    // Io's floods leave no generic lava (F117): a margin cell that took a flood's coarse
    // lava in the resampling and that no pass reached keeps the plains.
    if io_regime {
        let now = (present - IO_SURFACE_AGE_GYR) as f32;
        let units = layers.unit_id.as_mut_slice();
        let ages = layers.surface_age.as_mut_slice();
        let deep = layers.deep.as_mut_slice();
        let melt = layers.melt.as_mut_slice();
        let cover = layers.cover.as_mut_slice();
        for k in 0..units.len() {
            if units[k] == crate::units::UNIT_LAVA {
                units[k] = crate::units::UNIT_PRIMORDIAL;
                ages[k] = now;
                deep[k] = 0.0;
                melt[k] = 0.0;
                cover[k] = 0.0;
            }
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
    let roughness_share = roughness::apply(
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
            graph.f64("body.gravity"),
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
    let (maturity, skin) = if ice::is_ice_surface(graph) {
        let (maturity, skin) = ice::lag_fields(
            graph,
            &layers,
            present,
            slope.as_ref(),
            Some(&roughness_share),
        );
        (maturity, Some(skin))
    } else {
        let maturity = Field::from_fn(grid, |c, _| {
            let exposure = (present - layers.exposure_age(grid.index(c)) as f64).max(0.0);
            (1.0 - (-exposure / tau_sw).exp()) as f32
        });
        (maturity, None)
    };
    // Steep ground sheds its regolith and stays fresher than flat ground (F55): the
    // bright walls and scarps of every airless body, from the ground's own slopes and,
    // finer than the grid, the craters' walls at the scale they are steep (F98). An icy
    // surface's lag sheds from gentler slopes, the ice subliming from under it, and a
    // slope that has shed it is as clean as the ice; a silicate regolith matures in place
    // and stays only a share fresher.
    let maturity = {
        let steep = ice::slope_magnitude(&height, body_radius);
        let walls = layers.steep.as_slice();
        let (start, repose) = inputs.shed;
        let freshness = if ice::is_ice_surface(graph) {
            SLOPE_FRESHNESS_ICE
        } else {
            SLOPE_FRESHNESS
        };
        let mut m = maturity;
        m.as_mut_slice()
            .par_iter_mut()
            .zip(steep.as_slice())
            .zip(walls)
            .for_each(|((v, &s), &w)| {
                let shed = ((s as f64 - start) / (repose - start))
                    .clamp(0.0, 1.0)
                    .max(w as f64);
                *v *= 1.0 - (freshness * shed) as f32;
            });
        m
    };
    let rays = ray_field(
        &grid,
        body_radius,
        &records,
        &inputs,
        present,
        skin.as_ref(),
        &thrown,
    );
    if std::env::var_os("EXOGEN_DEBUG").is_some() {
        // The share of the surface the rays brighten past a tenth and a half of the fresh
        // contrast.
        let (mut tenth, mut half, mut whole) = (0.0, 0.0, 0.0);
        for (k, &v) in rays.as_slice().iter().enumerate() {
            let a = grid.cell_area(grid.cell(k));
            whole += a;
            if v > 0.1 {
                tenth += a;
            }
            if v > 0.5 {
                half += a;
            }
        }
        eprintln!(
            "rays brighten {:.1}% of the surface past a tenth of the fresh contrast, {:.2}% past a half",
            100.0 * tenth / whole,
            100.0 * half / whole
        );
    }

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
        craters: records.into_iter().chain(megas).collect(),
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
        coronae: rock_coronae,
        domes,
        paterae,
        mountains,
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
        ice_coronae,
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
    // The ice that flows is the shell under the crust, which large craters excavate (F126):
    // under a shell that floats on an ocean the frozen ocean, otherwise the crust's own mix.
    let rock = crate::stages::materials::Library::from_graph(graph)
        .non_ice(graph.str("materials.deep_bedrock"));
    let viscosity_scale = 10f64.powf(ICE_ROCK_VISCOSITY_DECADES * rock);
    let strength = graph.f64("interior.crust_strength");
    Some(Box::new(move |t: f64| {
        let state = thermal.at(t);
        relaxation::Shell {
            thickness: state.lithosphere,
            // The lid that holds the relief up is the one that keeps its stress for as
            // long as the relief stands (F197).
            elastic: relaxation::load_lid(
                state.lithosphere,
                surface_temperature,
                state.mantle_temperature,
                viscosity_scale,
                t,
            ),
            surface_temperature,
            basal_temperature: state.mantle_temperature,
            gravity,
            viscosity_scale,
            rock: false,
            strength,
        }
    }))
}

/// A rocky lid's crust through time for relaxation (F170), as a function of the time and
/// the crust's thickness: its base's temperature on the lid's conductive profile, from the
/// surface to the mantle at the lithosphere's base, and the elastic lid that holds up the
/// short loads. A stagnant or an episodic lid's (a stated scope): an icy surface's shell
/// relaxes it; a mobile lid's crust is its plates' (Phase 6); Io's heat pipe is not the
/// conductive lid its thermal history describes. An episodic lid's highlands stand on the
/// crust the mantle's flow thickened over its downwellings (F184), which relaxes as the
/// rest does.
fn rock_crust_history(
    graph: &ParamGraph,
    thermal: &ThermalHistory,
) -> Option<Box<dyn Fn(f64, f64) -> relaxation::Shell + Sync>> {
    if ice::is_ice_surface(graph)
        || !matches!(
            graph.str("interior.tectonic_regime"),
            "stagnant_lid" | "episodic_lid"
        )
        || io::is_io_regime(graph)
    {
        return None;
    }
    let thermal = thermal.clone();
    let surface_temperature = graph.f64("atmosphere.surface_temperature_mean");
    let gravity = graph.f64("body.gravity");
    let strength = graph.f64("interior.crust_strength");
    Some(Box::new(move |t: f64, crust: f64| {
        let state = thermal.at(t);
        let share = (crust / state.lithosphere.max(1.0)).min(1.0);
        relaxation::Shell {
            thickness: crust,
            elastic: state.elastic,
            surface_temperature,
            basal_temperature: surface_temperature
                + (state.mantle_temperature - surface_temperature) * share,
            gravity,
            viscosity_scale: 1.0,
            rock: true,
            strength,
        }
    }))
}

/// A crater that throws secondaries: its id, centre, diameter, time and impact.
struct Primary {
    id: u64,
    center: [f64; 3],
    diameter: f64,
    time_gyr: f64,
    azimuth: f64,
    angle: f64,
}

/// Whether a crater `diameter` m across throws secondaries of `d_min` m and more in any
/// number to count: its largest secondary at least half that size (F129).
fn throws_visible(diameter: f64, d_min: f64) -> bool {
    largest_secondary(diameter, CONTINUOUS_EJECTA_MULTIPLE) >= 0.5 * d_min
}

/// The secondaries the craters `primaries` throw (F129), each with the crater it came from:
/// as many of `d_min` m and more as Orientale's census gives for the crater's largest
/// secondary, each sized from that count, landing in one of the clusters its rays'
/// fragments make (F128) beyond its continuous blanket, among those that hold secondaries
/// its size, chosen by the ground they excavate, and lying in that cluster's footprint, so
/// that the larger lie nearer, along the jets, and a cluster that takes several holds them
/// as a group. Each strikes low along its path, at the primary's time.
fn sample_secondaries(
    primaries: &[Primary],
    inputs: &impacts::CraterInputs,
    d_min: f64,
    present: f64,
    cell: f64,
) -> Vec<(crate::crater::CraterSample, u64)> {
    let rng = Rng::new(inputs.seed ^ 0x5EC0_7DA2);
    let scale = inputs.secondary_fraction / SECONDARY_FRACTION_REFERENCE;
    let body_radius = inputs.body_radius;
    let per_crater: Vec<Vec<(crate::crater::CraterSample, u64)>> = primaries
        .par_iter()
        .map(|p| {
            if p.time_gyr > present || !throws_visible(p.diameter, d_min) {
                return Vec::new();
            }
            let largest = largest_secondary(p.diameter, CONTINUOUS_EJECTA_MULTIPLE);
            let expected = scale * secondaries_above(d_min, largest);
            let (count, _) = poisson(&rng, p.id, 0, expected);
            let count = count.min(SECONDARY_MAX_PER_CRATER);
            if count == 0 {
                return Vec::new();
            }
            // Past the cap the smallest are dropped: the sizes come from the top of the count.
            let pool = expected.min(SECONDARY_MAX_PER_CRATER as f64) / scale.max(1e-12);
            let radius = 0.5 * p.diameter;
            let thrown = throw(p.id, p.diameter, p.azimuth, p.angle, inputs, cell);
            // The clusters beyond the continuous blanket, those that hold the largest
            // secondaries first, with the ground they excavate summed: the secondaries land
            // in the clusters the rays show, not the solid ray past the cap, which would
            // follow the curtain's strongest jet alone.
            let mut sites: Vec<(f64, &Cluster)> = thrown
                .clusters
                .iter()
                .filter(|c| c.r >= CONTINUOUS_EJECTA_MULTIPLE * radius)
                .map(|c| (c.width / RAY_CLUSTER_WIDTH, c))
                .collect();
            if sites.is_empty() {
                return Vec::new();
            }
            sites.sort_by(|a, b| b.0.total_cmp(&a.0));
            let cumulative: Vec<f64> = sites
                .iter()
                .scan(0.0, |sum, (_, c)| {
                    *sum += c.amount;
                    Some(*sum)
                })
                .collect();
            let (east, north) = local_basis(p.center);
            (0..count as u64)
                .map(|j| {
                    let draw = 16 + j * 8;
                    let u = 1.0 - rng.uniform(p.id, draw);
                    let diameter = secondary_size_at(u * pool, largest).max(d_min);
                    // The clusters that hold secondaries this size, or the nearest if none.
                    let k = sites.partition_point(|(own, _)| *own >= diameter).max(1);
                    let pick = rng.uniform(p.id, draw + 1) * cumulative[k - 1];
                    let c = sites[cumulative[..k].partition_point(|&a| a < pick).min(k - 1)].1;
                    // Where in its footprint: along the radial and across it, opening
                    // downrange, as the cluster lays its ground (F128).
                    let at = offset(p.center, east, north, c.r / body_radius, c.phi);
                    let (ce, cn) = local_basis(at);
                    let down = azimuth_from(at, ce, cn, p.center) + std::f64::consts::PI;
                    let (su, sv) = (c.length / 2.5, c.width / 2.5);
                    let along = su * rng.normal(p.id, draw + 2);
                    let across = sv
                        * (1.0 + RAY_CLUSTER_SPLAY * (along / su).max(0.0))
                        * rng.normal(p.id, draw + 4);
                    let center = offset(
                        at,
                        ce,
                        cn,
                        along.hypot(across) / body_radius,
                        down + across.atan2(along),
                    );
                    let (se, sn) = local_basis(center);
                    (
                        crate::crater::CraterSample {
                            center,
                            diameter,
                            time_gyr: p.time_gyr,
                            azimuth: azimuth_from(center, se, sn, p.center) + std::f64::consts::PI,
                            angle: SECONDARY_IMPACT_ANGLE,
                            ellipticity: 1.2,
                            key: mix64(p.id ^ mix64(j + 1)),
                        },
                        p.id,
                    )
                })
                .collect()
        })
        .collect();
    if std::env::var_os("EXOGEN_DEBUG").is_some() {
        let throwers = per_crater.iter().filter(|v| !v.is_empty()).count();
        let total: usize = per_crater.iter().map(Vec::len).sum();
        eprintln!("secondaries: {total} from {throwers} craters");
    }
    per_crater.into_iter().flatten().collect()
}

/// What a crater throws (F128, F129): its rays, the edge and halo of its ejecta, how far its
/// rays reach and the clusters their fragments' secondaries make, the same for its rays as
/// for the secondaries it stamps.
struct Throw {
    rays: Rays,
    edge: EdgeShape,
    halo: f64,
    reach: f64,
    clusters: Vec<Cluster>,
}

fn throw(
    id: u64,
    diameter: f64,
    azimuth: f64,
    angle: f64,
    inputs: &impacts::CraterInputs,
    cell: f64,
) -> Throw {
    let radius = 0.5 * diameter;
    let dim = dimensions(diameter, &inputs.scaling);
    let rays = Rays::new(Curtain::new(id, azimuth, angle), radius);
    // The blanket's deposit runs out along each jet until it is too thin to show, the
    // strongest jets furthest, rather than ending on the footprint's circle (F96): where a
    // lag hides the rays beyond, as on Callisto, a common edge would draw a circle.
    let edge = EdgeShape::new(id, dim.ejecta_radius, cell);
    let thinnest = -(1.0 - RAY_VISIBLE).log2();
    let halo = (EJECTA_RADIUS_MULTIPLE * radius)
        .max(CONTINUOUS_EJECTA_MULTIPLE * radius * (edge.strongest_jet() / thinnest).cbrt());
    // No ray reaches past where its fragments still dig below the skin (F212), nor past a
    // quarter of the way round (F127).
    let reach = rays
        .fresh_reach()
        .min(std::f64::consts::FRAC_PI_2 * inputs.body_radius)
        .max(halo);
    let least = RAY_CLUSTER_DRAW_CELLS * cell;
    let body_radius = inputs.body_radius;
    let clusters = rays.clusters(id, reach, least, &|r| ray_coherence(r / body_radius));
    Throw {
        rays,
        edge,
        halo,
        reach,
        clusters,
    }
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

/// The coarse grid of buckets the waves judge overlap on: a stamp's cap reaches the
/// buckets whose centres lie within its radius and a bucket and a half of its centre, so
/// that two caps that overlap always reach a bucket in common.
const WAVE_BUCKET_FACE: u32 = 128;

/// Group the stamps into waves of non-overlapping footprints. Stamps in one wave touch
/// disjoint cells, so they can be computed in parallel and applied in any order, and each
/// stamp joins the first wave with room after every wave that holds an earlier stamp it
/// overlaps, so that every overlapping pair keeps its order while a stamp never waits for
/// one it does not touch. Overlap is judged on a coarse grid of buckets, each holding the
/// latest wave of the stamps whose caps reach it: a stamp may wait for one beside it in a
/// bucket it shares, never less than its own overlaps ask.
fn build_waves(caps: &[([f64; 3], f64)]) -> Vec<Vec<usize>> {
    let buckets = CubeSphere::new(WAVE_BUCKET_FACE);
    let margin = 1.5 * buckets.delta();
    let mut latest: Vec<Option<usize>> = vec![None; buckets.cell_count()];
    let mut stamper = Stamper::new(&buckets);
    let mut waves: Vec<Vec<usize>> = Vec::new();
    let mut areas: Vec<f64> = Vec::new();
    for (i, &(c, rho)) in caps.iter().enumerate() {
        stamper.gather(&buckets, 1.0, c, rho + margin);
        let mut w = stamper
            .footprint
            .iter()
            .filter_map(|&(b, _, _)| latest[b])
            .max()
            .map_or(0, |last| last + 1);
        let share = 0.5 * (1.0 - rho.min(std::f64::consts::PI).cos());
        while w < waves.len()
            && (waves[w].len() >= WAVE_MAX_STAMPS
                || (!waves[w].is_empty() && areas[w] + share > WAVE_MAX_AREA))
        {
            w += 1;
        }
        if w == waves.len() {
            waves.push(Vec::new());
            areas.push(0.0);
        }
        waves[w].push(i);
        areas[w] += share;
        for &(b, _, _) in &stamper.footprint {
            latest[b] = Some(latest[b].map_or(w, |last| last.max(w)));
        }
    }
    waves
}

fn ray_field(
    grid: &CubeSphere,
    body_radius: f64,
    craters: &[CraterRecord],
    inputs: &impacts::CraterInputs,
    present: f64,
    skin: Option<&Field<f32>>,
    thrown: &HashMap<u64, Vec<crate::crater::CraterSample>>,
) -> Field<f32> {
    let cell = grid.delta() * body_radius;
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
        .map(|c| {
            let own = thrown.get(&c.id).map_or(&[][..], |v| v.as_slice());
            RaySystem::new(c, inputs, present, cell, own)
        })
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
        let mut clustered = vec![0.0f32; grid.cell_count()];
        s.lay_clusters_par(grid, |k, v| clustered[k] += v);
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
                clustered[k] as f64,
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

/// The share of a crater's ray fragments that still land on it at `theta`, radians, round
/// the body (F127): all near the crater, thinning from RAY_COHERENT_K of the circular
/// orbital speed to none at a quarter of the way round, beyond which a jet's fragments
/// scatter over the far hemisphere or escape; a ray thins, and breaks into its clusters,
/// as it ends (F128).
fn ray_coherence(theta: f64) -> f64 {
    if theta >= std::f64::consts::FRAC_PI_2 {
        return 0.0;
    }
    let t = (0.5 * theta).tan();
    let k = 2.0 * t / (1.0 + t);
    let x = ((1.0 - k) / (1.0 - RAY_COHERENT_K)).clamp(0.0, 1.0);
    x * x * (3.0 - 2.0 * x)
}

/// One young crater's brightening (F91): full inside the rim; its continuous blanket,
/// falling with the deposit's cover to half at its lobed edge and on past it, mottled;
/// and beyond, its rays, the share of the ground its ejecta curtain's fragments
/// excavated from under the weathered skin. All of it fades with the crater's age, and on
/// ice the lag that hides its rays hides its floor and blanket as it grows (F127).
struct RaySystem {
    body_radius: f64,
    /// The clusters its rays' secondaries make, and the secondaries the finish stamped
    /// (F128).
    clusters: Vec<crate::crater::curtain::Cluster>,
    /// Where the stamped secondaries' own clusters begin among them, which are laid whole;
    /// the crater's id, which seeds its clusters' grain; and the narrowest cluster drawn, m
    /// (F130).
    halos_from: usize,
    seed: u64,
    least: f64,
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
    fn new(
        c: &CraterRecord,
        inputs: &impacts::CraterInputs,
        present: f64,
        cell: f64,
        secondaries: &[crate::crater::CraterSample],
    ) -> Self {
        let lifetime = inputs.ray_lifetime.max(1e-6);
        let age = present - c.time_gyr;
        let radius = 0.5 * c.diameter;
        let dim = dimensions(c.diameter, &inputs.scaling);
        // Its rays are the clusters its fragments' secondaries make, and the secondaries the
        // finish stamped are clusters of their own, their fresh ground where they lie, a
        // crater's width round it and twice that downrange (F128, F129).
        let Throw {
            rays,
            edge,
            halo,
            reach,
            mut clusters,
        } = throw(c.id, c.diameter, c.azimuth, c.angle, inputs, cell);
        let least = RAY_CLUSTER_DRAW_CELLS * cell;
        let halos_from = clusters.len();
        let (east, north) = local_basis(c.center);
        for s in secondaries {
            let width = (2.0 * s.diameter).max(least);
            clusters.push(crate::crater::curtain::Cluster {
                r: angle_between(c.center, s.center) * inputs.body_radius,
                phi: azimuth_from(c.center, east, north, s.center),
                width,
                length: 2.0 * width,
                amount: 2.0 * width * width,
            });
        }
        Self {
            body_radius: inputs.body_radius,
            clusters,
            halos_from,
            seed: c.id,
            least,
            center: c.center,
            radius,
            azimuth: c.azimuth,
            angle: c.angle,
            freshness: (1.0 - age / lifetime).powi(2),
            edge,
            ejecta_edge: CONTINUOUS_EJECTA_MULTIPLE * radius,
            halo,
            texture: (c.diameter >= RAY_TEXTURE_MIN_CELLS * cell)
                .then(|| CraterTexture::new(c.id, &dim, cell)),
            rays,
            reach,
        }
    }

    /// The area the clusters excavate per unit area, cell by cell, handed to `add` (F128).
    fn lay_clusters(
        &self,
        grid: &CubeSphere,
        stamper: &mut Stamper,
        mut add: impl FnMut(usize, f32),
    ) {
        let mut out = Vec::new();
        for (i, c) in self.clusters.iter().enumerate() {
            out.clear();
            self.lay_cluster(grid, stamper, i, c, &mut out);
            for &(k, v) in &out {
                add(k, v);
            }
        }
    }

    /// As `lay_clusters`, a chunk of clusters at a time in parallel, for a system whose rays
    /// cover much of the body.
    fn lay_clusters_par(&self, grid: &CubeSphere, mut add: impl FnMut(usize, f32)) {
        const CHUNK: usize = 128;
        for (n, chunk) in self.clusters.chunks(CHUNK).enumerate() {
            let parts: Vec<Vec<(usize, f32)>> = chunk
                .par_iter()
                .enumerate()
                .map(|(j, c)| {
                    with_stamper(grid, |stamper| {
                        let mut out = Vec::new();
                        self.lay_cluster(grid, stamper, n * CHUNK + j, c, &mut out);
                        out
                    })
                })
                .collect();
            for part in parts {
                for (k, v) in part {
                    add(k, v);
                }
            }
        }
    }

    /// The ground the `i`th cluster excavates per unit area, cell by cell, into `out` (F128,
    /// F130): its grain, each part a Gaussian along the radial whose width opens downrange
    /// by RAY_CLUSTER_SPLAY, laying its whole share on the cells it covers; a stamped
    /// secondary's own cluster is laid whole.
    fn lay_cluster(
        &self,
        grid: &CubeSphere,
        stamper: &mut Stamper,
        i: usize,
        c: &Cluster,
        out: &mut Vec<(usize, f32)>,
    ) {
        let (east, north) = local_basis(self.center);
        let r2 = self.body_radius * self.body_radius;
        let at = offset(self.center, east, north, c.r / self.body_radius, c.phi);
        let (ce, cn) = local_basis(at);
        let down = azimuth_from(at, ce, cn, self.center) + std::f64::consts::PI;
        let grains = if i < self.halos_from {
            c.grain(self.least, mix64(self.seed ^ mix64(i as u64 + 1)))
        } else {
            vec![Grain {
                u: 0.0,
                v: 0.0,
                su: c.length / 2.5,
                sv: c.width / 2.5,
                amount: c.amount,
            }]
        };
        let reach = grains
            .iter()
            .map(|g| g.u.hypot(g.v) + 3.0 * g.su.max(g.sv * (1.0 + 3.0 * RAY_CLUSTER_SPLAY)))
            .fold(0.0, f64::max);
        stamper.gather(grid, self.body_radius, at, reach);
        // The cells in the cluster's frame, along the radial outward and across it, ordered
        // along it, so that each part reads only the cells in its reach.
        let mut cells: Vec<(f64, f64, usize, f64)> = stamper
            .footprint
            .iter()
            .map(|&(idx, r, phi)| {
                (
                    r * (phi - down).cos(),
                    r * (phi - down).sin(),
                    idx,
                    grid.cell_area(grid.cell(idx)),
                )
            })
            .collect();
        cells.sort_by(|a, b| a.0.total_cmp(&b.0));
        let mut weights: Vec<(usize, f64)> = Vec::new();
        for g in &grains {
            let lo = cells.partition_point(|x| x.0 < g.u - 3.0 * g.su);
            let hi = cells.partition_point(|x| x.0 <= g.u + 3.0 * g.su);
            weights.clear();
            let mut sum = 0.0;
            for &(u, v, idx, area) in &cells[lo..hi] {
                let du = u - g.u;
                let spread = g.sv * (1.0 + RAY_CLUSTER_SPLAY * (du / g.su).max(0.0));
                let dv = v - g.v;
                if dv.abs() > 4.3 * spread {
                    continue;
                }
                let k = (-0.5 * ((du / g.su).powi(2) + (dv / spread).powi(2))).exp();
                if k > 1e-4 {
                    weights.push((idx, k));
                    sum += k * area;
                }
            }
            if sum <= 0.0 {
                let point = offset(
                    at,
                    ce,
                    cn,
                    g.u.hypot(g.v) / self.body_radius,
                    down + g.v.atan2(g.u),
                );
                let k = grid.index(grid.locate(point).0);
                out.push((k, (g.amount / (grid.cell_area(grid.cell(k)) * r2)) as f32));
                continue;
            }
            for &(idx, k) in &weights {
                out.push((idx, (g.amount * k / (sum * r2)) as f32));
            }
        }
    }

    /// The brightening at distance `r`, m, and azimuth `phi` from the centre, where the
    /// skin leaves the share `skin` of the area the rays' secondaries would excavate and
    /// the clusters there excavate `clustered` per unit area.
    fn value(&self, r: f64, phi: f64, skin: f64, clustered: f64) -> f32 {
        if r <= self.radius {
            return (self.freshness * skin) as f32;
        }
        let edge = self.ejecta_edge * self.edge.factor_at(phi, r / self.radius);
        // The blanket is bright where its deposit covers the ground (F51), fading with the
        // cover past its lobed edge rather than ending there, mottled.
        let blanket = (ejecta_coverage(r, self.radius, edge, self.halo, self.edge.jet(phi))
            * ejecta_modulation(phi, self.azimuth, self.angle).min(1.0)
            * self.texture.as_ref().map_or(1.0, |t| {
                1.0 + EJECTA_BRIGHTNESS_TEXTURE * t.noise_floor(phi, r / self.radius)
            }))
        .clamp(0.0, 1.0)
            * skin;
        // The rays are their clusters; where the jets excavate more than clusters carry,
        // the rest is a solid ray (F128).
        let landing = ray_coherence(r / self.body_radius);
        let solid = (landing * self.rays.area(r, phi) - RAY_CLUSTER_CAP).max(0.0);
        // What they excavate is fresh by the share of it from beneath the skin.
        let ray = (1.0 - (-(solid + clustered) * skin).exp()) * self.rays.fresh_share(r, phi, skin);
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
        let mut clustered: HashMap<usize, f32> = HashMap::new();
        self.lay_clusters(grid, stamper, |k, v| {
            *clustered.entry(k).or_insert(0.0) += v
        });
        stamper.gather(grid, body_radius, self.center, self.reach);
        stamper
            .footprint
            .iter()
            .filter_map(|&(idx, r, phi)| {
                let v = self.value(
                    r,
                    phi,
                    skin.map_or(1.0, |s| s[idx] as f64),
                    clustered.get(&idx).copied().unwrap_or(0.0) as f64,
                );
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
    use crate::geom::dot;
    use crate::presets::preset;

    fn preset_graph(name: &str, quality: &str) -> ParamGraph {
        let (mut g, _) = preset(name).unwrap().to_graph();
        g.set("detail.quality", quality).unwrap();
        g
    }

    /// Prints each body's elastic lid through time and the loads its floods lay: the lava's
    /// thickness over the plains, when it was laid, and the flexural parameter then.
    #[test]
    #[ignore = "diagnostic: the lid through time and the floods' loads"]
    fn print_loads() {
        let names =
            std::env::var("EXOGEN_BODY").unwrap_or_else(|_| "Moon,Mercury,Mars,Venus".into());
        let quality = std::env::var("EXOGEN_QUALITY").unwrap_or_else(|_| "preview".into());
        for name in names.split(',') {
            let g = preset_graph(name, &quality);
            let out = run_coarse_loop(&g, &NoProgress).unwrap();
            let thermal = &out.thermal;
            let age = g.f64("body.age");
            let gravity = g.f64("body.gravity");
            let lids: Vec<String> = [0.1, 0.25, 0.5, 0.75, 1.0]
                .iter()
                .map(|f| {
                    format!(
                        "{:.2} Gyr {:.0} km",
                        f * age,
                        thermal.at(f * age).elastic / 1e3
                    )
                })
                .collect();
            println!("{name}: elastic lid {}", lids.join(", "));
            let n = out.grid.cell_count();
            let mut thickness = vec![0.0f64; n];
            let mut weighted_time = vec![0.0f64; n];
            let mut volume_time = (0.0, 0.0);
            for rec in &out.floods {
                for (idx, c) in rec.cells() {
                    thickness[idx as usize] += c.thickness as f64;
                    weighted_time[idx as usize] += c.thickness as f64 * rec.time_gyr;
                }
                volume_time.0 += rec.volume * rec.time_gyr;
                volume_time.1 += rec.volume;
            }
            let mut lava: Vec<f64> = thickness.iter().copied().filter(|&t| t > 1.0).collect();
            lava.sort_by(f64::total_cmp);
            let q = |v: &[f64], f: f64| {
                if v.is_empty() {
                    0.0
                } else {
                    v[((v.len() - 1) as f64 * f) as usize]
                }
            };
            let mean_time = volume_time.0 / volume_time.1.max(1.0);
            let te = thermal.at(mean_time).elastic;
            let rigidity = ROCK_YOUNG_MODULUS * te.powi(3) / (12.0 * (1.0 - ROCK_POISSON.powi(2)));
            let l = (rigidity / (DENSITY_ROCK * gravity)).powf(0.25);
            println!(
                "  lava over {:.2} of cells: thickness m p50 {:.0} p90 {:.0} p99 {:.0} max {:.0}; volume-weighted time {:.2} Gyr, lid then {:.0} km, flexural parameter (D / rho_m g)^1/4 {:.0} km; loop cell {:.0} km",
                lava.len() as f64 / n as f64,
                q(&lava, 0.5),
                q(&lava, 0.9),
                q(&lava, 0.99),
                q(&lava, 1.0),
                mean_time,
                te / 1e3,
                l / 1e3,
                out.grid.delta() * g.f64("body.radius") / 1e3
            );
            if let Some(r) = &out.resurfacing {
                let mut fill: Vec<f64> = r
                    .fill
                    .iter()
                    .map(|&f| f as f64)
                    .filter(|&f| f > 1.0)
                    .collect();
                fill.sort_by(f64::total_cmp);
                let te = thermal.at(r.time_gyr).elastic;
                println!(
                    "  resurfacing at {:.2} Gyr over {:.2} of cells: fill m p10 {:.0} p50 {:.0} p90 {:.0} max {:.0}; lid then {:.1} km",
                    r.time_gyr,
                    fill.len() as f64 / n as f64,
                    q(&fill, 0.1),
                    q(&fill, 0.5),
                    q(&fill, 0.9),
                    q(&fill, 1.0),
                    te / 1e3
                );
            }
        }
    }

    /// Prints the thrusts and ridges a body's loop lays (F155, F156): their numbers and
    /// lengths, the share of the ridges' length below the mean height, and, about a rise,
    /// the share running round it.
    #[test]
    #[ignore = "diagnostic: the thrusts and ridges a body's loop lays"]
    fn print_thrusts() {
        let names = std::env::var("EXOGEN_BODY").unwrap_or_else(|_| "Venus,Mars,Mercury".into());
        let quality = std::env::var("EXOGEN_QUALITY").unwrap_or_else(|_| "preview".into());
        for name in names.split(',') {
            let g = preset_graph(name, &quality);
            let t0 = Instant::now();
            let out = run_coarse_loop(&g, &NoProgress).unwrap();
            let r = g.f64("body.radius");
            let grid = out.grid;
            let layers = &out.keyframes.last().unwrap().layers;
            let height = layers.bedrock.as_slice();
            let (mut sum, mut area) = (0.0, 0.0);
            for (k, &h) in height.iter().enumerate() {
                let a = grid.cell_area(grid.cell(k));
                sum += h as f64 * a;
                area += a;
            }
            let mean = sum / area;
            let carried = out.scarps.iter().filter(|s| s.carried).count();
            let taken: f64 = out
                .scarps
                .iter()
                .map(|s| s.length * tectonics::scarp_shortening(s.height))
                .sum();
            let mut lengths: Vec<f64> = out.scarps.iter().map(|s| s.length / 1e3).collect();
            lengths.sort_by(f64::total_cmp);
            let q = |v: &[f64], f: f64| {
                if v.is_empty() {
                    0.0
                } else {
                    v[((v.len() - 1) as f64 * f) as usize]
                }
            };
            println!(
                "{name}: loop {:.0} s; scarps {} ({} carried), length km p50 {:.0} p90 {:.0} max {:.0}, shortening {:.3e} m2 ({:.2e} of the area)",
                t0.elapsed().as_secs_f64(),
                out.scarps.len(),
                carried,
                q(&lengths, 0.5),
                q(&lengths, 0.9),
                q(&lengths, 1.0),
                taken,
                taken / (4.0 * std::f64::consts::PI * r * r)
            );
            let mut times: Vec<f64> = out.scarps.iter().map(|s| s.time_gyr).collect();
            times.sort_by(f64::total_cmp);
            let lids: Vec<String> = [0.0, 0.05, 0.1, 0.2, 0.4, 0.7, 1.0]
                .iter()
                .map(|&t| format!("{t} Gyr {:.0} km", out.thermal.at(t).elastic / 1e3))
                .collect();
            println!(
                "  scarp times Gyr p10 {:.2} p50 {:.2} p90 {:.2}; lid {}",
                q(&times, 0.1),
                q(&times, 0.5),
                q(&times, 0.9),
                lids.join(", ")
            );
            let thrust = |id: u64| id >= 1 << 32;
            for (label, pick) in [("plains ridges", false), ("thrust ridges", true)] {
                let rs: Vec<&ridges::RidgeRecord> =
                    out.ridges.iter().filter(|x| thrust(x.id) == pick).collect();
                let total: f64 = rs.iter().map(|x| x.length).sum();
                let mut below = 0.0f64;
                let mut all = 0.0f64;
                let mut round = 0.0f64;
                let mut near = 0.0f64;
                for x in &rs {
                    for w in x.points.windows(2) {
                        let k = grid.index(grid.locate(w[0]).0);
                        all += 1.0;
                        if (height[k] as f64) < mean {
                            below += 1.0;
                        }
                        // Round the nearest rise: the step across the direction out from it.
                        if let Some(rise) = out.rises.iter().min_by(|a, b| {
                            angle_between(a.center, w[0]).total_cmp(&angle_between(b.center, w[0]))
                        }) {
                            let d = angle_between(rise.center, w[0]) * r;
                            if d < 4.0 * rise.radius && d > 0.5 * rise.radius {
                                near += 1.0;
                                let along = dot(rise.center, w[0]);
                                let o = [
                                    w[0][0] * along - rise.center[0],
                                    w[0][1] * along - rise.center[1],
                                    w[0][2] * along - rise.center[2],
                                ];
                                let on = dot(o, o).sqrt().max(1e-12);
                                let step =
                                    [w[1][0] - w[0][0], w[1][1] - w[0][1], w[1][2] - w[0][2]];
                                let sn = dot(step, step).sqrt().max(1e-12);
                                // Within 30 degrees of round: the step's share outward under sin 30.
                                if (dot(o, step) / (on * sn)).abs() < 0.5 {
                                    round += 1.0;
                                }
                            }
                        }
                    }
                }
                println!(
                    "  {label}: {} of {:.0} km in all; {:.2} of their length below the mean height; within 0.5 to 4 rise radii: {:.2} of it, {:.2} of that within 30 degrees of round",
                    rs.len(),
                    total / 1e3,
                    below / all.max(1.0),
                    near / all.max(1.0),
                    round / near.max(1.0)
                );
            }
            // The graben the loads opened (F158): how many, how long, when, and how far
            // from the lava they lie.
            let load_graben: Vec<&plume::Trough> = out
                .graben
                .iter()
                .filter(|g| (60_000_000_000..61_000_000_000).contains(&g.id))
                .collect();
            if !load_graben.is_empty() {
                let mut lengths: Vec<f64> = load_graben.iter().map(|g| g.length / 1e3).collect();
                lengths.sort_by(f64::total_cmp);
                let mut times: Vec<f64> = load_graben.iter().map(|g| g.start_gyr).collect();
                times.sort_by(f64::total_cmp);
                let units = layers.unit_id.as_slice();
                let lava_cells: Vec<[f64; 3]> = (0..grid.cell_count())
                    .step_by(7)
                    .filter(|&k| volcanism::is_lava_unit(units[k]))
                    .map(|k| grid.center(grid.cell(k)))
                    .collect();
                let mut near: Vec<f64> = load_graben
                    .iter()
                    .map(|g| {
                        let m = g.points[g.points.len() / 2];
                        lava_cells
                            .iter()
                            .map(|&c| angle_between(c, m) * r / 1e3)
                            .fold(f64::INFINITY, f64::min)
                    })
                    .collect();
                near.sort_by(f64::total_cmp);
                println!(
                    "  load graben: {} , length km p50 {:.0} max {:.0}; width km p50 {:.1}; formed Gyr p10 {:.2} p50 {:.2}; midpoint from lava km p50 {:.0} p90 {:.0}",
                    load_graben.len(),
                    q(&lengths, 0.5),
                    q(&lengths, 1.0),
                    {
                        let mut w: Vec<f64> = load_graben.iter().map(|g| g.width / 1e3).collect();
                        w.sort_by(f64::total_cmp);
                        q(&w, 0.5)
                    },
                    q(&times, 0.1),
                    q(&times, 0.5),
                    q(&near, 0.5),
                    q(&near, 0.9)
                );
            }
            // The plains below the mean height, for comparison.
            let (mut lava, mut lava_below) = (0.0, 0.0);
            for (k, &h) in height.iter().enumerate() {
                if layers.unit_id.as_slice()[k] == crate::units::UNIT_LAVA {
                    let a = grid.cell_area(grid.cell(k));
                    lava += a;
                    if (h as f64) < mean {
                        lava_below += a;
                    }
                }
            }
            println!(
                "  plains: {:.2} of the surface, {:.2} of them below the mean height",
                lava / area,
                lava_below / lava.max(1e-12)
            );
        }
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
        assert!(rayed > 0.001 && rayed < 0.3, "rayed fraction {rayed}");
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
            // Inside the floor and the wall: a complex crater's rim, scalloped by the
            // slumps of its outermost terrace (F67), comes inside 0.45 of its diameter at
            // places, and beyond it the ejecta keep the age of the ground under them (F100).
            let inner = 0.35 * c.diameter / c.ellipticity.sqrt();
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
        // The drawn scarps carry most of the budget; the rest is the lid's elastic strain,
        // the thrusts beneath the plains and the faults too short to draw (F155).
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

    /// F174: a rise is a pile as well as a dome: the magma that stays below the eruptions,
    /// laid where the melt rises from, gathers beneath Mars's rise as Tharsis's some 3 × 10⁸
    /// km³ (Phillips et al. 2001), and all the loop lays is the melt's intrusive share and the
    /// magma the full edifices could not lift, the melt one budget.
    #[test]
    fn mars_s_rise_gathers_its_pile() {
        let g = preset_graph("Mars", "preview");
        let out = run_coarse_loop(&g, &NoProgress).unwrap();
        let r = g.f64("body.radius");
        let grid = out.grid;
        let rise = &out.rises[0];
        let (mut all, mut under) = (0.0f64, 0.0f64);
        for k in 0..grid.cell_count() {
            let c = grid.cell(k);
            let a = grid.cell_area(c) * r * r;
            let laid: f64 = out.intrusions.iter().map(|(_, f)| f[k] as f64).sum();
            all += laid * a;
            if angle_between(grid.center(c), rise.center) * r < rise.radius {
                under += laid * a;
            }
        }
        let thermal = ThermalHistory::from_graph(&g);
        let floods = volcanism::FloodProcess::new(&g, grid, &thermal, &out.rises).unwrap();
        let steps = out.steps;
        let dt = out.age_gyr / steps as f64;
        let intrusive: f64 = (0..steps)
            .map(|s| floods.step_intrusive(dt * s as f64, dt * (s + 1) as f64))
            .sum();
        let surplus: f64 = out
            .shields
            .iter()
            .map(|s| {
                let holds = s.volume / (1.0 - s.overflow).max(1e-9);
                (plume::shield_supplied(out.age_gyr, s) - holds).max(0.0)
            })
            .sum();
        let thickness = |k: usize| out.intrusions.iter().map(|(_, f)| f[k] as f64).sum::<f64>();
        let centre = thickness(grid.index(grid.locate(rise.center).0));
        let least = (0..grid.cell_count())
            .map(thickness)
            .fold(f64::MAX, f64::min);
        println!(
            "Mars: {:.3e} km3 intruded, {:.3e} km3 beneath its rise, {:.1} km at its centre and {:.1} km at the least; the melt's intrusive share {:.3e} km3 and the edifices' surplus {:.3e} km3",
            all / 1e9,
            under / 1e9,
            centre / 1e3,
            least / 1e3,
            intrusive / 1e9,
            surplus / 1e9
        );
        assert!(
            ((intrusive + surplus) - all).abs() < 0.01 * all,
            "{all} of {}",
            intrusive + surplus
        );
        assert!(under > 1e17 && under < 1e18, "{under}");
    }

    /// F173: the ground stands on the mantle's flow: Mercury's surface rises and sinks as the
    /// pattern beneath it changes, late in its life by as much as the ground the flow holds
    /// up changes, little else moving it then, and the flow's highs and lows reach some
    /// hundreds of metres.
    #[test]
    fn the_ground_stands_on_the_mantle_s_flow() {
        let g = preset_graph("Mercury", "preview");
        let out = run_coarse_loop(&g, &NoProgress).unwrap();
        let support = out.support.as_ref().unwrap();
        let grid = out.grid;
        let frames = &out.keyframes;
        let (a, b) = (&frames[frames.len() - 2], &frames[frames.len() - 1]);
        let (ha, hb) = (a.layers.bedrock.as_slice(), b.layers.bedrock.as_slice());
        let (mut cov, mut var, mut lo, mut hi) = (0.0, 0.0, f64::MAX, f64::MIN);
        for k in 0..grid.cell_count() {
            let c = grid.cell(k);
            let d = grid.center(c);
            let area = grid.cell_area(c);
            let now = support.at(b.time_gyr, d);
            let dl = now - support.at(a.time_gyr, d);
            cov += area * (hb[k] - ha[k]) as f64 * dl;
            var += area * dl * dl;
            lo = lo.min(now);
            hi = hi.max(now);
        }
        let slope = cov / var.max(1e-12);
        println!(
            "Mercury from {:.2} to {:.2} Gyr: the ground moved {slope:.3} of the flow's change; the flow holds it {lo:.0} to {hi:.0} m",
            a.time_gyr, b.time_gyr
        );
        assert!((slope - 1.0).abs() < 0.15, "{slope}");
        assert!(hi - lo > 500.0 && hi - lo < 6_000.0, "{lo} {hi}");
    }

    /// Mars's thrusts run round Tharsis, which squeezes its foot past the lid's strength
    /// (F155, F156), and Venus's ridges keep to its lows, where the relief squeezes the
    /// plains and its highs stretch them (Bilotti and Suppe 1999); its thrusts form soon after
    /// its last overturn (F185).
    #[test]
    fn mars_s_thrusts_run_round_tharsis_and_venus_s_ridges_keep_to_its_lows() {
        let g = preset_graph("Mars", "preview");
        let mars = run_coarse_loop(&g, &NoProgress).unwrap();
        let r = g.f64("body.radius");
        let rise = &mars.rises[0];
        assert!(!mars.scarps.is_empty());
        let thrust_ridges: Vec<&ridges::RidgeRecord> =
            mars.ridges.iter().filter(|x| x.id >= 1 << 32).collect();
        assert!(!thrust_ridges.is_empty());
        // Within four radii of the rise, most of the thrusts' length runs round it.
        let (mut near, mut round) = (0usize, 0usize);
        let points = mars
            .scarps
            .iter()
            .map(|s| &s.points)
            .chain(thrust_ridges.iter().map(|x| &x.points));
        for trace in points {
            for w in trace.windows(2) {
                let d = angle_between(rise.center, w[0]) * r;
                if d < 0.5 * rise.radius || d > 4.0 * rise.radius {
                    continue;
                }
                near += 1;
                let along = dot(rise.center, w[0]);
                let out = [
                    w[0][0] * along - rise.center[0],
                    w[0][1] * along - rise.center[1],
                    w[0][2] * along - rise.center[2],
                ];
                let step = [w[1][0] - w[0][0], w[1][1] - w[0][1], w[1][2] - w[0][2]];
                let cos = dot(out, step) / (dot(out, out) * dot(step, step)).sqrt().max(1e-30);
                if cos.abs() < 0.5 {
                    round += 1;
                }
            }
        }
        assert!(
            near > 100 && round * 10 > near * 6,
            "{round} of {near} round"
        );

        let g = preset_graph("Venus", "preview");
        let venus = run_coarse_loop(&g, &NoProgress).unwrap();
        // Its thrusts form on the lid its last overturn left (F185).
        let overturn = ThermalHistory::from_graph(&g)
            .mantle
            .unwrap()
            .last_overturn()
            .unwrap();
        assert!(!venus.scarps.is_empty());
        let length: f64 = venus.scarps.iter().map(|s| s.length).sum();
        let mut by_age = [0.0f64; 5];
        for s in &venus.scarps {
            assert!(s.time_gyr > overturn, "{} after {overturn}", s.time_gyr);
            let k = (((s.time_gyr - overturn) / 0.05) as usize).min(4);
            by_age[k] += s.length / length;
        }
        eprintln!("Venus's thrusts by 50 Myr since the overturn: {by_age:.2?}, {:.0} km", length / 1e3);
        let grid = venus.grid;
        let height = venus.keyframes.last().unwrap().layers.bedrock.as_slice();
        let (mut sum, mut area) = (0.0, 0.0);
        for (k, &h) in height.iter().enumerate() {
            let a = grid.cell_area(grid.cell(k));
            sum += h as f64 * a;
            area += a;
        }
        let mean = sum / area;
        let (mut below, mut all) = (0usize, 0usize);
        for x in &venus.ridges {
            for p in &x.points {
                all += 1;
                if (height[grid.index(grid.locate(*p).0)] as f64) < mean {
                    below += 1;
                }
            }
        }
        assert!(
            all > 1000 && below * 10 > all * 8,
            "{below} of {all} below the mean"
        );
    }

    /// The ground bends under its loads (F162), a shield's lava ponds past its foot (F163)
    /// and its root weighs on the lid (F166): Mars's largest shield's moat holds the lava its
    /// flows sent past its foot, the shields' own lava lying far thicker there than farther
    /// off, and still the ground just past its foot stands below the ground beyond; its
    /// summit still stands three quarters of the height it built above the ground beyond;
    /// and the loop keeps how far the ground sank for the finish.
    #[test]
    fn mars_s_largest_shield_s_moat_holds_its_lava() {
        let g = preset_graph("Mars", "preview");
        let mars = run_coarse_loop(&g, &NoProgress).unwrap();
        assert!(!mars.flexure.is_empty());
        let r = g.f64("body.radius");
        let grid = mars.grid;
        let height = mars.keyframes.last().unwrap().layers.bedrock.as_slice();
        // The lava the shields' own eruptions laid past their feet on each coarse cell, poured
        // out at a shield's rates where the floods' are a flood basalt's (F163).
        let mut lava = vec![0.0f64; grid.cell_count()];
        for f in mars
            .floods
            .iter()
            .filter(|f| f.effusion < FLOOD_EFFUSION_RATE.0)
        {
            for (&idx, c) in &f.cells() {
                lava[idx as usize] += c.thickness as f64;
            }
        }
        let shield = mars
            .shields
            .iter()
            .max_by(|a, b| a.radius.total_cmp(&b.radius))
            .expect("a shield");
        let (east, north) = local_basis(shield.center);
        // The mean of a field round the shield at a distance from its summit.
        let ring = |field: &dyn Fn(usize) -> f64, d: f64| {
            let mut sum = 0.0;
            for k in 0..72 {
                let phi = std::f64::consts::TAU * k as f64 / 72.0;
                let p = offset(shield.center, east, north, d / r, phi);
                sum += field(grid.index(grid.locate(p).0));
            }
            sum / 72.0
        };
        let ground = |k: usize| height[k] as f64;
        let laid = |k: usize| lava[k];
        let moat = ring(&ground, 1.15 * shield.radius);
        let beyond = ring(&ground, 2.0 * shield.radius);
        let near = [1.1, 1.3, 1.5]
            .iter()
            .map(|&f| ring(&laid, f * shield.radius))
            .sum::<f64>()
            / 3.0;
        let far = [2.5, 3.0, 3.5]
            .iter()
            .map(|&f| ring(&laid, f * shield.radius))
            .sum::<f64>()
            / 3.0;
        assert!(near > 3.0 * far + 300.0, "{near} {far}");
        // Its root deepens the moat past what its lava fills (F166).
        assert!(moat < beyond - 300.0, "{moat} {beyond}");
        // Near: its root sinks it (F166), and the ground beyond keeps the relief the
        // ejecta about it draped (F169).
        assert!(
            ring(&ground, 0.0) - beyond > 0.75 * shield.height,
            "summit {} beyond {} height {}",
            ring(&ground, 0.0),
            beyond,
            shield.height
        );
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
        // The maria carry wrinkle ridges, all of them on lava, but for those a later
        // mega-impact's ejecta buried (F168); a thrust's ridge folds the plains where most of
        // its trace crossed lava (F155) and may cross older ground besides.
        let ridges: Vec<&ridges::RidgeRecord> = out
            .loop_output
            .ridges
            .iter()
            .filter(|r| r.id < 1 << 32)
            .collect();
        assert!(ridges.len() > 20, "{} ridges", ridges.len());
        let r_body = g.f64("body.radius");
        let megas: Vec<&CraterRecord> = out
            .loop_output
            .craters
            .iter()
            .filter(|c| mega::is_mega(c.diameter, r_body))
            .collect();
        let g = out.loop_output.grid;
        for r in ridges {
            if megas
                .iter()
                .any(|m| m.time_gyr > r.time_gyr && mega::renews(m, r.midpoint(), r_body))
            {
                continue;
            }
            let idx = g.index(g.locate(r.midpoint()).0);
            assert!(volcanism::is_lava_unit(units[idx]));
        }
    }

    /// Mars raises one plume: a rise kilometres high with shields on it, a canyon
    /// kilometres deep on its flank, graben, and lava concentrated around it.
    /// Diagnostic: each of Mars's canyons with its window, carved depth and the relief it
    /// keeps below its rims at a quarter, half and three quarters along its trace, and the
    /// shields' starts; run with `--ignored --nocapture`.
    #[test]
    #[ignore]
    fn print_mars_canyon_relief() {
        let mut r = Runner::new();
        let quality = std::env::var("EXOGEN_QUALITY").unwrap_or_else(|_| "preview".into());
        let g = preset_graph("Mars", &quality);
        let out = r.run(&g, &NoProgress).unwrap();
        let lo = &out.loop_output;
        let grid = out.finish.grid;
        let radius = g.f64("body.radius");
        let height = out.finish.height.as_slice();
        let at = |d: [f64; 3]| height[grid.index(grid.locate(d).0)] as f64;
        let relief = |canyon: &plume::Trough, share: f64| {
            let k = ((canyon.points.len() - 2) as f64 * share) as usize;
            let (mid, next) = (canyon.points[k], canyon.points[k + 1]);
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
        println!(
            "shields start: {:?}",
            lo.shields
                .iter()
                .map(|s| (s.start_gyr * 100.0).round() / 100.0)
                .collect::<Vec<_>>()
        );
        println!("start end depth_m width_km | relief_m at 1/4 1/2 3/4");
        let mut canyons: Vec<&plume::Trough> = lo.canyons.iter().collect();
        canyons.sort_by(|a, b| a.start_gyr.total_cmp(&b.start_gyr));
        let mut deepest = f64::NEG_INFINITY;
        for c in canyons {
            let rel: Vec<f64> = [0.25, 0.5, 0.75].iter().map(|&f| relief(c, f)).collect();
            deepest = rel.iter().cloned().fold(deepest, f64::max);
            println!(
                "{:.2} {:.2} {:.0} {:.0} | {:.0} {:.0} {:.0}",
                c.start_gyr,
                c.end_gyr,
                c.depth,
                c.width / 1e3,
                rel[0],
                rel[1],
                rel[2]
            );
        }
        println!("deepest standing {deepest:.0} m");
    }

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
        // Some trough the loop carried still stands over a kilometre deep below the
        // plateau either side along most of its length, the rims being the highest ground
        // within a few widths past the en-echelon neighbours, its median taken over 19
        // points along it since a local fill decides any one of them; a later flood may
        // have filled another. The canyons open as the dome rises (F116), and the rise's
        // eruptions fill them in part until its shields gather its dikes as its stretching
        // wanes (F119).
        let relief_at = |canyon: &plume::Trough, share: f64| {
            let k = ((canyon.points.len() - 2) as f64 * share) as usize;
            let (mid, next) = (canyon.points[k], canyon.points[k + 1]);
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
        let relief = |canyon: &plume::Trough| {
            let mut along: Vec<f64> = (1..20)
                .map(|k| relief_at(canyon, k as f64 / 20.0))
                .collect();
            along.sort_by(|a, b| a.total_cmp(b));
            along[9]
        };
        let deepest = lo
            .canyons
            .iter()
            .map(relief)
            .fold(f64::NEG_INFINITY, f64::max);
        assert!(deepest > 1_000.0, "deepest standing trough {deepest} m");
        // The Moon raises none.
        let moon = r
            .run(&preset_graph("Moon", "preview"), &NoProgress)
            .unwrap();
        assert!(moon.loop_output.rises.is_empty() && moon.loop_output.shields.is_empty());
    }

    /// F197 diagnostic: the lid that holds an icy body's relief up through its life against
    /// the tectonic lid, and how much of a crater's relief its shell lets survive from its
    /// time to the present, at twice its diameter; run with `--ignored --nocapture`.
    #[test]
    #[ignore]
    fn print_crater_relaxation() {
        for name in ["Ganymede", "Callisto", "Europa", "Enceladus", "Mimas", "Miranda"] {
            let g = preset_graph(name, "standard");
            let thermal = ThermalHistory::from_graph(&g);
            let Some(shell_at) = ice_shell_history(&g, &thermal) else {
                continue;
            };
            let age = g.f64("body.age");
            let lids: Vec<String> = [0.5, 2.0, age]
                .iter()
                .map(|&t| {
                    format!(
                        "{t:.1} Gyr {:.1} km (tectonic {:.1})",
                        shell_at(t).elastic / 1e3,
                        thermal.at(t).elastic / 1e3
                    )
                })
                .collect();
            println!("{name}: load lid {}", lids.join(", "));
            for before in [4.0, 3.0, 1.0] {
                let kept: Vec<String> = [5e3, 10e3, 20e3, 40e3, 80e3]
                    .iter()
                    .map(|&d| {
                        let s = relaxation::survival_since(
                            2.0 * d,
                            (age - before).max(0.0),
                            age,
                            &*shell_at,
                            1.0,
                        );
                        format!("{:.0} km {s:.2}", d / 1e3)
                    })
                    .collect();
                println!("  formed {before} Ga, kept at twice the diameter: {}", kept.join(", "));
            }
        }
    }

    /// T6 diagnostic: how much of a lane's grooves' relief the shell's relaxation lets
    /// survive from the lane's time to the present, at the groove spacing; run with
    /// `--ignored --nocapture`.
    #[test]
    #[ignore]
    fn print_lane_relaxation() {
        for name in ["Ganymede", "Enceladus", "Miranda"] {
            let g = preset_graph(name, "standard");
            let thermal = ThermalHistory::from_graph(&g);
            let Some(shell_at) = ice_shell_history(&g, &thermal) else {
                continue;
            };
            let grid = CubeSphere::new(64);
            let Some(p) = ice_shell::IceShellProcess::new(&g, grid) else {
                continue;
            };
            let age = g.f64("body.age");
            let mut kept: Vec<(f64, f64, f64)> = p
                .lanes
                .iter()
                .map(|l| {
                    let s = relaxation::survival_since(
                        l.groove_spacing,
                        l.time_gyr,
                        age,
                        &*shell_at,
                        1.0,
                    );
                    (s, l.groove_spacing, l.time_gyr)
                })
                .collect();
            kept.sort_by(|a, b| a.0.total_cmp(&b.0));
            let n = kept.len();
            if n == 0 {
                println!("{name}: no lanes");
                continue;
            }
            let q = |f: f64| kept[((n - 1) as f64 * f) as usize];
            println!(
                "{name}: {n} lanes, gravity {:.3}; survival p10 {:.3} p50 {:.3} p90 {:.3}; spacing p50 {:.0} m; time p50 {:.2} Gyr of {age:.2}",
                g.f64("body.gravity"),
                q(0.1).0,
                q(0.5).0,
                q(0.9).0,
                q(0.5).1,
                q(0.5).2
            );
        }
    }

    /// F187 diagnostic: the strain the tesserae's lid took from the flow, set by set, over
    /// the tesserae; run with `--ignored --nocapture`.
    #[test]
    #[ignore]
    fn print_tessera_strain() {
        let mut r = Runner::new();
        let g = preset_graph("Venus", "standard");
        let out = r.run(&g, &NoProgress).unwrap();
        let rec = out.loop_output.resurfacing.as_ref().unwrap();
        let n = rec.score.len();
        let tessera: Vec<usize> = (0..n).filter(|&k| rec.score[k] > rec.threshold).collect();
        println!("event {:.3} Gyr, {} tessera cells of {n}", rec.time_gyr, tessera.len());
        let q = |v: &mut Vec<f32>, f: f64| {
            if v.is_empty() {
                return 0.0;
            }
            v.sort_by(|a, b| a.total_cmp(b));
            v[((v.len() - 1) as f64 * f) as usize]
        };
        for s in &rec.strain {
            let (mut short, mut ext, mut iso) = (Vec::new(), Vec::new(), Vec::new());
            let mut any = 0usize;
            for k in 0..n {
                let [a, b, c] = s.strain[k];
                if a == 0.0 && b == 0.0 && c == 0.0 {
                    continue;
                }
                any += 1;
            }
            for &k in &tessera {
                let [a, b, c] = s.strain[k];
                if a == 0.0 && b == 0.0 && c == 0.0 {
                    continue;
                }
                let m = 0.5 * (a + c);
                let r = (0.25 * (a - c) * (a - c) + b * b).sqrt();
                short.push(m - r);
                ext.push(m + r);
                iso.push(2.0 * m);
            }
            let strained = short.len();
            let past = |v: &[f32]| -> Vec<f32> {
                tessera
                    .iter()
                    .map(|&k| v[k] / 1e6)
                    .filter(|&x| x > 0.0)
                    .collect()
            };
            let (mut pressed, mut pulled) = (past(&s.pressed), past(&s.pulled));
            println!(
                "  pressed past strength on {} tessera cells, p50 {:.1} p90 {:.1} MPa; pulled on {}, p50 {:.1} p90 {:.1} MPa",
                pressed.len(),
                q(&mut pressed, 0.5),
                q(&mut pressed, 0.9),
                pulled.len(),
                q(&mut pulled, 0.5),
                q(&mut pulled, 0.9)
            );
            println!(
                "set at {:.4} Gyr, Te {:.1} km: {any} cells strained, {strained} of the tesserae; least principal p10 {:.3} p50 {:.3} p90 {:.3}; greatest p10 {:.3} p50 {:.3} p90 {:.3}; areal p50 {:.3}",
                s.time,
                s.elastic / 1e3,
                q(&mut short, 0.1),
                q(&mut short, 0.5),
                q(&mut short, 0.9),
                q(&mut ext, 0.1),
                q(&mut ext, 0.5),
                q(&mut ext, 0.9),
                q(&mut iso, 0.5),
            );
        }
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
        // The tessera meets the plains by the share of each cell its margin covers (F131),
        // not along whole cells.
        {
            let grid = layers.grid;
            let units = layers.unit_id.as_slice();
            let cover = layers.cover.as_slice();
            let (mut contact, mut partial) = (0usize, 0usize);
            for k in 0..units.len() {
                if units[k] != crate::units::UNIT_TESSERA && units[k] != crate::units::UNIT_LAVA {
                    continue;
                }
                let cell = grid.cell(k);
                let other = [(1, 0), (-1, 0), (0, 1), (0, -1)].iter().any(|&(di, dj)| {
                    let u = units[grid.index(grid.neighbor(cell, di, dj))];
                    (u == crate::units::UNIT_TESSERA || u == crate::units::UNIT_LAVA)
                        && u != units[k]
                });
                if other {
                    contact += 1;
                    if cover[k] > 0.02 && cover[k] < 0.98 {
                        partial += 1;
                    }
                }
            }
            assert!(
                contact > 0 && partial * 4 > contact,
                "{partial} of {contact} cells on the tessera's contact partly covered"
            );
            // The plains onlap the tesserae (F186): the fabric's ridges stand as inliers
            // beyond the edge of the old ground's smoothed height, and the lava runs up its
            // valleys within it, as far as the grid draws the fabric (F187).
            let rec = lo.resurfacing.as_ref().unwrap();
            let coarse = crate::grid::CubeSphere::new(((rec.score.len() / 6) as f64).sqrt() as u32);
            let score = crate::grid::projection::PaddedField::new(&coarse, &rec.score);
            let (mut inliers, mut embayed, mut near) = (0usize, 0usize, 0usize);
            for (k, &unit) in units.iter().enumerate() {
                let d = grid.center(grid.cell(k));
                let s =
                    crate::grid::projection::sample_bilinear(&coarse, &score, d) - rec.threshold;
                if s.abs() > 300.0 {
                    continue;
                }
                near += 1;
                if s < 0.0 && unit == crate::units::UNIT_TESSERA {
                    inliers += 1;
                } else if s > 0.0 && unit == crate::units::UNIT_LAVA {
                    embayed += 1;
                }
            }
            assert!(
                inliers > 0 && embayed > 0 && (inliers + embayed) * 20 > near,
                "{inliers} inliers and {embayed} embayed of {near} cells about the edge"
            );
        }
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
        // The salts' margins are drawn by the share of each cell they cover (F127), not
        // as whole cells.
        let partial = units
            .iter()
            .zip(cover)
            .filter(|(u, c)| {
                (**u == crate::units::UNIT_CRATER_FLOOR || **u == crate::units::UNIT_DEEP_FLOOR)
                    && **c > 0.02
                    && **c < 0.98
            })
            .count();
        assert!(
            partial * 10 > salted,
            "{partial} of {salted} salted floor cells partly covered"
        );
        // The domes' salts end by the share of each cell their edifice covers (F131).
        let salt = CRYODOME_SALT as f32;
        let domed = units
            .iter()
            .filter(|u| **u == crate::units::UNIT_CRYODOME)
            .count();
        let edged = units
            .iter()
            .zip(cover)
            .filter(|(u, c)| **u == crate::units::UNIT_CRYODOME && **c > 0.01 && **c < 0.95 * salt)
            .count();
        assert!(
            domed > 0 && edged * 20 > domed,
            "{edged} of {domed} dome cells partly covered"
        );
    }

    /// A ray reaches no further than a quarter of the way round the body (F127): fading
    /// from RAY_COHERENT_K of the circular orbital speed, it is gone where its fragments
    /// would need that speed.
    #[test]
    fn rays_end_a_quarter_of_the_way_round() {
        assert_eq!(ray_coherence(0.1), 1.0);
        assert_eq!(ray_coherence(std::f64::consts::FRAC_PI_2), 0.0);
        let mut last = 1.0;
        for i in 0..=90 {
            let c = ray_coherence(i as f64 * std::f64::consts::FRAC_PI_2 / 90.0);
            assert!((0.0..=1.0).contains(&c) && c <= last + 1e-12, "{i}: {c}");
            last = c;
        }
        // The Moon's longest rays, Tycho's 1,500 km, are whole.
        assert!(ray_coherence(1_500e3 / 1_737.4e3) > 0.99);
    }

    /// F129: a basin's secondaries are as many as Orientale's census gives for its largest,
    /// lie in its rays' clusters beyond its continuous blanket, the larger nearer, and the
    /// same seed gives the same ones.
    #[test]
    fn a_basin_s_secondaries_lie_in_its_rays() {
        // The census's draw for a basin of its own, whatever the preset's seed.
        let mut g = preset_graph("Moon", "preview");
        g.set("seed.master", 13.0).unwrap();
        g.recompute();
        let inputs = impacts::CraterInputs::from_graph(&g);
        let body = inputs.body_radius;
        let cell = 2.0 * std::f64::consts::PI * body / (4.0 * 2048.0);
        let d_min = FINISHING_MIN_CELLS * cell;
        let basin = Primary {
            id: 75,
            center: [0.0, 0.0, 1.0],
            diameter: 378e3,
            time_gyr: 4.0,
            azimuth: 0.3,
            angle: 0.8,
        };
        let out = sample_secondaries(std::slice::from_ref(&basin), &inputs, d_min, 4.5, cell);
        let again = sample_secondaries(std::slice::from_ref(&basin), &inputs, d_min, 4.5, cell);
        assert_eq!(out.len(), again.len());
        assert!(
            out.iter()
                .zip(&again)
                .all(|(a, b)| a.0.center == b.0.center)
        );
        let largest = largest_secondary(basin.diameter, CONTINUOUS_EJECTA_MULTIPLE);
        let expected = inputs.secondary_fraction / SECONDARY_FRACTION_REFERENCE
            * secondaries_above(d_min, largest);
        let n = out.len() as f64;
        assert!(
            (n - expected).abs() < 4.0 * expected.sqrt() + 1.0,
            "{n} against {expected:.0}"
        );
        let radius = 0.5 * basin.diameter;
        let thrown = throw(
            basin.id,
            basin.diameter,
            basin.azimuth,
            basin.angle,
            &inputs,
            cell,
        );
        let (east, north) = local_basis(basin.center);
        let sites: Vec<([f64; 3], f64)> = thrown
            .clusters
            .iter()
            .filter(|c| c.r >= CONTINUOUS_EJECTA_MULTIPLE * radius)
            .map(|c| {
                (
                    offset(basin.center, east, north, c.r / body, c.phi),
                    c.length,
                )
            })
            .collect();
        let mut by_size: Vec<(f64, f64)> = Vec::new();
        for (s, parent) in &out {
            assert_eq!(*parent, basin.id);
            assert!(s.diameter >= d_min);
            let r = angle_between(basin.center, s.center) * body;
            assert!(r > 1.5 * radius, "a secondary {:.0} km out", r / 1e3);
            let held = sites
                .iter()
                .any(|(at, length)| angle_between(*at, s.center) * body < 3.0 * length);
            assert!(held, "a secondary {:.0} km out in no cluster", r / 1e3);
            by_size.push((s.diameter, r));
        }
        by_size.sort_by(|a, b| b.0.total_cmp(&a.0));
        let median = |v: &[(f64, f64)]| {
            let mut r: Vec<f64> = v.iter().map(|x| x.1).collect();
            r.sort_by(f64::total_cmp);
            r[r.len() / 2]
        };
        let tenth = by_size.len() / 10;
        let (near, rest) = (median(&by_size[..tenth]), median(&by_size[tenth..]));
        assert!(
            near < rest,
            "the largest at {near:.0} m, the rest at {rest:.0} m"
        );
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
        // body (Schenk 2002), a depth ratio near 0.012: relaxed, not flattened. The median
        // of some sixty craters moves by a tenth with which craters the history leaves:
        // 0.0085 in the fifty-sixth milestone and 0.0072 in the fifty-seventh, whose colder
        // surface left 54 of the same craters at the same median depth.
        assert!(small > 0.006, "small craters keep their depth: {small}");
        assert!(large < 0.25 * small, "large {large} vs small {small}");
    }

    /// Io keeps no craters, carries hundreds of paterae and dozens of mountains, its flow
    /// fields cover about the share Williams et al. (2011) map (27.8 percent; F201), and
    /// its surface is as bright and as pale a yellow as Io's measured colour (F204: green
    /// albedo 0.62 to 0.74, blue some 0.72 of green).
    #[test]
    fn io_is_craterless_with_paterae_and_mountains() {
        let mut r = Runner::new();
        let g = preset_graph("Io", "preview");
        let out = r.run(&g, &NoProgress).unwrap();
        let lo = &out.loop_output;
        assert!(lo.paterae.len() > 100, "{} paterae", lo.paterae.len());
        assert!(lo.mountains.len() > 20, "{} mountains", lo.mountains.len());
        assert!(!out.finish.floods.is_empty() && !out.finish.plume_episodes.is_empty());
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
        assert!(flows > 0.15 && flows < 0.4, "flow share {flows}");
        let albedo = out.finish.albedo.as_slice();
        let mean: [f64; 3] = albedo.iter().fold([0.0; 3], |a, c| {
            [a[0] + c[0] as f64, a[1] + c[1] as f64, a[2] + c[2] as f64]
        });
        let n = albedo.len() as f64;
        let (rr, gg, bb) = (mean[0] / n, mean[1] / n, mean[2] / n);
        assert!(
            rr > gg && gg > bb && gg > 0.5 && gg < 0.8,
            "mean colour {rr:.2} {gg:.2} {bb:.2}"
        );
        assert!(
            bb / gg > 0.62 && bb / gg < 0.82,
            "mean colour {rr:.2} {gg:.2} {bb:.2}"
        );
        // Frost does not hold on the dark lava still showing (F206).
        let layers = &out.finish.layers;
        let (deep, cover) = (layers.deep.as_slice(), layers.cover.as_slice());
        let dark: Vec<usize> = (0..deep.len()).filter(|&k| deep[k] > 0.9).collect();
        assert!(!dark.is_empty());
        let frosted = dark.iter().map(|&k| cover[k] as f64).sum::<f64>() / dark.len() as f64;
        assert!(frosted < 0.1, "frost over the dark lava {frosted:.3}");
    }

    /// Diagnostic (F170): the share of each wavelength's relief a rocky lid's lower crust
    /// lets survive from a time to the present, for each rocky preset, a crust 1.5 times the
    /// starting one and on Mars one 100 km thick; run with `--ignored --nocapture`.
    #[test]
    #[ignore]
    fn print_rock_relaxation() {
        for name in ["Moon", "Mercury", "Mars", "Venus", "Earth", "Io"] {
            let g = preset_graph(name, "preview");
            let thermal = ThermalHistory::from_graph(&g);
            let Some(crust_at) = rock_crust_history(&g, &thermal) else {
                println!("{name}: none");
                continue;
            };
            let age = g.f64("body.age");
            // The crust 1.5 times the starting one, and on Mars a crust 100 km thick too.
            let crusts: Vec<f64> = if name == "Mars" {
                vec![1.5 * g.f64("tectonics.initial_crust_thickness"), 100e3]
            } else if name == "Venus" {
                vec![g.f64("tectonics.initial_crust_thickness"), 65e3]
            } else {
                vec![1.5 * g.f64("tectonics.initial_crust_thickness")]
            };
            let overturn = thermal.mantle.as_ref().and_then(|m| m.last_overturn());
            let after: Vec<f64> = overturn
                .map(|o| [0.005, 0.02, 0.05, 0.1, 0.3].iter().map(|d| o + d).collect())
                .unwrap_or_default();
            for crust in crusts {
            for from in [0.1, 0.5, 1.0, 3.5, 4.0].into_iter().chain(after.clone()) {
                let steps = 200;
                let dt = (age - from) / steps as f64;
                let kept: Vec<String> = [300e3, 1_000e3, 3_000e3]
                    .iter()
                    .map(|&w| {
                        let mut k = 1.0;
                        for i in 0..steps {
                            let t = from + (i as f64 + 0.5) * dt;
                            k *= relaxation::survival(w, dt, &crust_at(t, crust));
                        }
                        format!("{:.0} km {:.3}", w / 1e3, k)
                    })
                    .collect();
                let s = crust_at(from, crust);
                println!("{name} crust {:.0} km from {from} Gyr (base {:.0} K, lid {:.0} km): {}", crust / 1e3, s.basal_temperature, s.elastic / 1e3, kept.join(", "));
            }
            }
        }
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
        // The rows within 40 degrees of the equator, every fourth.
        let rows: Vec<Vec<f64>> = (width / 4 - 113..=width / 4 + 113)
            .step_by(4)
            .map(|y| row[y * width..(y + 1) * width].iter().map(|v| *v as f64).collect())
            .collect();
        let slope = mean_spectral_slope(&rows, 4, 128);
        assert!(slope < -1.3 && slope > -3.8, "spectral slope {slope}");
    }

    /// Log-log slope of the periodogram averaged over periodic rows, between two
    /// wavenumbers: steadier than one row's, whose power scatters by a factor of several at
    /// each wavenumber.
    fn mean_spectral_slope(rows: &[Vec<f64>], k_lo: usize, k_hi: usize) -> f64 {
        let n = rows[0].len();
        let mut xs = Vec::new();
        let mut ys = Vec::new();
        let mut k = k_lo;
        while k <= k_hi {
            let mut power = 0.0;
            for row in rows {
                let mean = row.iter().sum::<f64>() / n as f64;
                let (mut re, mut im) = (0.0, 0.0);
                for (i, v) in row.iter().enumerate() {
                    let ang = 2.0 * std::f64::consts::PI * k as f64 * i as f64 / n as f64;
                    re += (v - mean) * ang.cos();
                    im -= (v - mean) * ang.sin();
                }
                power += re * re + im * im;
            }
            xs.push((k as f64).ln());
            ys.push((power / rows.len() as f64).max(1e-30).ln());
            k = (k as f64 * 1.25).ceil() as usize;
        }
        let mx = xs.iter().sum::<f64>() / xs.len() as f64;
        let my = ys.iter().sum::<f64>() / ys.len() as f64;
        let num: f64 = xs.iter().zip(&ys).map(|(x, y)| (x - mx) * (y - my)).sum();
        let den: f64 = xs.iter().map(|x| (x - mx).powi(2)).sum();
        num / den
    }
}
