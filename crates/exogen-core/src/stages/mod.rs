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

pub mod impacts;
pub mod materials;
pub mod primordial_crust;

use std::collections::HashMap;
use std::sync::Arc;
use std::time::{Duration, Instant};

use rayon::prelude::*;

use crate::constants::*;
use crate::crater::production::poisson;
use crate::crater::{
    CraterRecord, PopulationConfig, ProfileTable, RestoreParams, StampParams, Stamper, dimensions,
    restore, sample_population, stamp,
};
use crate::fields::{Field, LayerSet};
use crate::geom::{local_basis, offset};
use crate::grid::CubeSphere;
use crate::params::ParamGraph;
use crate::params::derive::physics::volatile;
use crate::rng::{Rng, hash_words};
use crate::units::{UNIT_BASIN_FLOOR, UNIT_CRATER_FLOOR, UNIT_EJECTA, UNIT_SECONDARY_FIELD};

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
            ]),
            primordial_crust::VERSION,
            impacts::VERSION,
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
pub const FINISH_VERSION: u64 = 3;

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

    let crust = primordial_crust::PrimordialCrust::new(graph, grid);
    let mut layers = crust.init();
    let mut impacts = impacts::ImpactProcess::new(graph, grid);
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
        impacts.step(&mut layers, t0, t1, s as u64);
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
    Ok(LoopOutput {
        grid,
        keyframes,
        steps,
        age_gyr: age,
        craters: impacts.records,
        loop_threshold,
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
    let mut stamper = Stamper::new(&grid);
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
    let base: Vec<f32> = layers.bedrock.as_slice().to_vec();
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
        Replay(usize),
    }
    let mut events: Vec<(f64, u8, u64, Event)> =
        Vec::with_capacity(loop_craters.len() + replay.len());
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
    for (i, (_, _, _, event)) in events.iter().enumerate() {
        match event {
            Event::Loop(li) => {
                let c = &loop_craters[*li];
                let dim = dimensions(c.diameter, &inputs.scaling);
                // The detail the coarse grid could not carry: the resolved profile minus
                // its coarse-cell-smoothed version.
                let fresh = ProfileTable::new(dim, 0.0);
                let smooth = ProfileTable::new(dim, 0.5 * smooth_sigma * smooth_sigma);
                let mut detail = fresh.clone();
                for (v, s) in detail.values.iter_mut().zip(&smooth.values) {
                    *v -= s;
                }
                let params = RestoreParams {
                    time_gyr: c.time_gyr as f32,
                    unit_floor: if dim.class == crate::crater::CraterClass::MultiRing {
                        UNIT_BASIN_FLOOR
                    } else {
                        UNIT_CRATER_FLOOR
                    },
                    unit_ejecta: UNIT_EJECTA,
                    azimuth: c.azimuth,
                    ellipticity: c.ellipticity,
                    detail_weight: 1.0,
                };
                restore(
                    &mut layers,
                    &grid,
                    body_radius,
                    c.center,
                    &base,
                    &fresh,
                    &detail,
                    &params,
                    &mut stamper,
                );
            }
            Event::Replay(j) => {
                let (c, secondary) = &replay[*j];
                let dim = dimensions(c.diameter, &inputs.scaling);
                let kappa_t = inputs.diffusivity_m2_per_myr * (present - c.time_gyr).max(0.0) * 1e3;
                let table = ProfileTable::new(dim, kappa_t);
                let params = StampParams {
                    time_gyr: c.time_gyr as f32,
                    unit_floor: if *secondary {
                        UNIT_SECONDARY_FIELD
                    } else {
                        UNIT_CRATER_FLOOR
                    },
                    unit_ejecta: UNIT_EJECTA,
                    azimuth: c.azimuth,
                    ellipticity: c.ellipticity,
                    amplitude: if *secondary {
                        SECONDARY_DEPTH_FACTOR
                    } else {
                        1.0
                    },
                    replace_floor: true,
                    set_age: true,
                };
                stamp(
                    &mut layers,
                    &grid,
                    body_radius,
                    c.center,
                    &table,
                    &params,
                    &mut stamper,
                );
                if c.diameter >= record_min && !*secondary {
                    records.push(CraterRecord {
                        center: c.center,
                        diameter: c.diameter,
                        time_gyr: c.time_gyr,
                        azimuth: c.azimuth,
                        ellipticity: c.ellipticity,
                        class: dim.class,
                        from_loop: false,
                        id: c.key,
                    });
                }
            }
        }
        if i % 20_000 == 0 {
            if sink.cancelled() {
                return Err(RunError::Cancelled);
            }
            sink.report(Progress {
                stage: "finishing",
                fraction: 0.3 + 0.4 * i as f32 / total as f32,
                message: format!("stamping craters {i} of {total}"),
            });
        }
    }
    let replayed_count = replay.len();

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
    let rays = ray_field(
        &grid,
        body_radius,
        &records,
        &inputs,
        present,
        cell,
        &mut stamper,
    );

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
        let clusters = (count / 6).max(1);
        for j in 0..count as u64 {
            let cluster = j % clusters as u64;
            let base = 8 + cluster * 4;
            let dist = radius * (r_lo + (r_hi - r_lo) * rng.uniform(p.id, base));
            let phi = 2.0 * std::f64::consts::PI * rng.uniform(p.id, base + 1);
            let draw = 1_000 + j * 4;
            let jitter_r = 0.15 * radius * (rng.uniform(p.id, draw) - 0.5);
            let jitter_phi = 0.3 * (rng.uniform(p.id, draw + 1) - 0.5) * radius / dist.max(1.0);
            let center = offset(
                p.center,
                east,
                north,
                (dist + jitter_r) / inputs.body_radius,
                phi + jitter_phi,
            );
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
                ellipticity: 1.2,
                key: (p.id << 16) | j,
            });
        }
    }
    out
}

/// Ray brightening from craters younger than the ray lifetime: a few bright radial sectors
/// out to several radii, plus the fresh continuous ejecta, fading with age.
fn ray_field(
    grid: &CubeSphere,
    body_radius: f64,
    craters: &[CraterRecord],
    inputs: &impacts::CraterInputs,
    present: f64,
    cell: f64,
    stamper: &mut Stamper,
) -> Field<f32> {
    let mut rays = Field::new(*grid);
    let lifetime = inputs.ray_lifetime.max(1e-6);
    let rng = Rng::new(inputs.seed ^ 0x9A7F_3C11);
    let (n_lo, n_hi) = RAYS_PER_CRATER;
    let (l_lo, l_hi) = RAY_LENGTH_MULTIPLE;
    for c in craters {
        let age = present - c.time_gyr;
        if age < 0.0 || age >= lifetime || c.diameter < RAY_MIN_CELLS * cell {
            continue;
        }
        let freshness = (1.0 - age / lifetime).powi(2);
        let radius = 0.5 * c.diameter;
        let n = n_lo + ((n_hi - n_lo) as f64 * rng.uniform(c.id, 0)) as u32;
        // Rays keep a roughly constant linear width, so their angular width shrinks with
        // distance, and large craters' rays are shorter in radii than small craters' rays.
        let (w_lo, w_hi) = RAY_WIDTH_FRACTION;
        let length_scale = (RAY_REFERENCE_DIAMETER / c.diameter).sqrt().min(1.0);
        let ray_params: Vec<(f64, f64, f64)> = (0..n)
            .map(|i| {
                let d = 4 + 3 * i as u64;
                let phi = 2.0 * std::f64::consts::PI * rng.uniform(c.id, d);
                let width = radius * (w_lo + (w_hi - w_lo) * rng.uniform(c.id, d + 1));
                let length =
                    radius * (l_lo + (l_hi - l_lo) * rng.uniform(c.id, d + 2)) * length_scale;
                (phi, width, length)
            })
            .collect();
        let max_len = ray_params
            .iter()
            .map(|p| p.2)
            .fold(0.0, f64::max)
            .max(CONTINUOUS_EJECTA_MULTIPLE * radius);
        stamper.gather(grid, body_radius, c.center, max_len);
        let out = rays.as_mut_slice();
        let ejecta_edge = CONTINUOUS_EJECTA_MULTIPLE * radius;
        for &(idx, r, phi) in &stamper.footprint {
            // Fresh ejecta: full inside the rim, fading to half at the continuous edge.
            let mut v: f64 = if r <= radius {
                1.0
            } else if r <= ejecta_edge {
                1.0 - 0.5 * (r - radius) / (ejecta_edge - radius)
            } else {
                0.0
            };
            if r > ejecta_edge {
                let falloff = (ejecta_edge / r).powf(1.2);
                for &(p, w, l) in &ray_params {
                    if r > l {
                        continue;
                    }
                    let dphi = (phi - p + std::f64::consts::PI)
                        .rem_euclid(2.0 * std::f64::consts::PI)
                        - std::f64::consts::PI;
                    let taper = ((l - r) / (0.3 * l)).clamp(0.0, 1.0);
                    let across = dphi * r / w;
                    v = v.max((-(across * across)).exp() * taper * falloff);
                }
            }
            let value = (v * freshness) as f32;
            if value > out[idx] {
                out[idx] = value;
            }
        }
    }
    rays
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
