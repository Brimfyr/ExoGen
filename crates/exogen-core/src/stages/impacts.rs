// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! Impact cratering inside the coarse loop (PLAN.md 4.5, ADR 0002): craters above the
//! loop threshold are sampled from the production function for each time step and
//! stamped onto the coarse layers at the step they form, so every later process acts on
//! them. Each one is recorded for the finishing pass.

use crate::constants::*;
use crate::crater::production::{expected_count_above, poisson};
use crate::crater::stamp::elliptical_radius;
use crate::crater::{
    CraterRecord, FluxHistory, ProductionFunction, ProfileTable, Scaling, StampParams, Stamper,
    dimensions, ellipticity_for, impact_angle, stamp,
};
use crate::fields::LayerSet;
use crate::grid::CubeSphere;
use crate::params::ParamGraph;
use crate::rng::Rng;
use crate::units::{
    UNIT_BASIN_FLOOR, UNIT_CRATER_FLOOR, UNIT_DEEP_EJECTA, UNIT_DEEP_FLOOR, UNIT_EJECTA,
    UNIT_LAVA_CRATER_FLOOR, UNIT_LAVA_EJECTA,
};

pub const VERSION: u64 = 13;

/// Body-level cratering inputs read once from the graph.
#[derive(Clone, Debug)]
pub struct CraterInputs {
    pub pf: ProductionFunction,
    /// Flux history in time since formation, with the record starting at solidification.
    pub flux: FluxHistory,
    pub scaling: Scaling,
    pub flux_scale: f64,
    pub body_radius: f64,
    pub age: f64,
    pub oblique_fraction: f64,
    pub secondary_fraction: f64,
    pub saturation_fraction: f64,
    pub retention_age: f64,
    /// Retention applied to loop craters. Regimes that resurface wholesale (mobile lid,
    /// episodic lid, a young ice shell) do not keep basins older than the retention age;
    /// stagnant lids keep theirs, since erosion removes small craters long before large
    /// ones. The tectonic and erosion stages replace this rule when they exist.
    pub loop_retention_age: f64,
    pub diffusivity_m2_per_myr: f64,
    pub ray_lifetime: f64,
    /// Thickness of the surface layer above the deep bedrock, metres.
    pub deep_bedrock_depth: f64,
    pub seed: u64,
}

/// Whether a crater's excavation reaches below the surface layer into the deep bedrock.
pub fn reaches_deep_bedrock(diameter: f64, deep_bedrock_depth: f64) -> bool {
    EXCAVATION_DEPTH_FRACTION * diameter > deep_bedrock_depth
}

/// Unit of a crater floor by class, excavation and target.
pub fn floor_unit(class: crate::crater::CraterClass, deep: bool, on_lava: bool) -> u16 {
    if deep {
        UNIT_DEEP_FLOOR
    } else if on_lava {
        UNIT_LAVA_CRATER_FLOOR
    } else if class == crate::crater::CraterClass::MultiRing {
        UNIT_BASIN_FLOOR
    } else {
        UNIT_CRATER_FLOOR
    }
}

/// Unit of continuous ejecta by excavation and target.
pub fn ejecta_unit(deep: bool, on_lava: bool) -> u16 {
    if deep {
        UNIT_DEEP_EJECTA
    } else if on_lava {
        UNIT_LAVA_EJECTA
    } else {
        UNIT_EJECTA
    }
}

/// Whether the ground at a point is lava plains, so a crater there is made of them.
pub fn target_is_lava(layers: &LayerSet, grid: &CubeSphere, center: [f64; 3]) -> bool {
    let idx = grid.index(grid.locate(center).0);
    crate::stages::volcanism::is_lava_unit(layers.unit_id.as_slice()[idx])
}

impl CraterInputs {
    pub fn from_graph(graph: &ParamGraph) -> Self {
        let pf = ProductionFunction::by_name(graph.str("cratering.production_function"));
        Self {
            pf,
            flux: FluxHistory {
                steady: pf.chronology.2,
                early_ratio: graph.f64("cratering.early_flux_ratio"),
                early_decay: graph.f64("cratering.early_flux_decay"),
                record_start: graph.f64("surface.solidification_time"),
            },
            scaling: Scaling {
                simple_complex: graph.f64("cratering.simple_complex_transition"),
                peak_ring: graph.f64("cratering.peak_ring_onset"),
                multi_ring: graph.f64("cratering.multiring_onset"),
                gravity: graph.f64("body.gravity"),
            },
            flux_scale: graph.f64("cratering.flux_relative_to_moon")
                * graph.f64("multiplier.crater_density"),
            body_radius: graph.f64("body.radius"),
            age: graph.f64("body.age"),
            oblique_fraction: graph.f64("cratering.oblique_fraction"),
            secondary_fraction: graph.f64("cratering.secondary_fraction"),
            saturation_fraction: graph.f64("cratering.empirical_saturation"),
            retention_age: graph.f64("cratering.retention_age"),
            loop_retention_age: if wholesale_resurfacing(graph) {
                graph.f64("cratering.retention_age")
            } else {
                graph.f64("body.age")
            },
            diffusivity_m2_per_myr: graph.f64("cratering.diffusivity"),
            ray_lifetime: graph.f64("cratering.ray_fading_timescale"),
            deep_bedrock_depth: graph.f64("materials.deep_bedrock_depth"),
            seed: graph.f64("seed.cratering") as u64,
        }
    }
}

/// Regimes that renew the whole surface within the retention age: plates, global
/// resurfacing, or an ice shell kept young by tidal heating.
pub fn wholesale_resurfacing(graph: &ParamGraph) -> bool {
    match graph.str("interior.tectonic_regime") {
        "mobile_lid" | "episodic_lid" => true,
        "ice_shell_ocean" => graph.f64("interior.tidal_flux") > ICE_RESURFACING_TIDAL_FLUX,
        _ => false,
    }
}

pub struct ImpactProcess {
    inputs: CraterInputs,
    grid: CubeSphere,
    d_loop: f64,
    d_max: f64,
    rng: Rng,
    stamper: Stamper,
    pub records: Vec<CraterRecord>,
}

impl ImpactProcess {
    pub fn new(graph: &ParamGraph, grid: CubeSphere) -> Self {
        let inputs = CraterInputs::from_graph(graph);
        let d_loop = graph.f64("cratering.loop_threshold");
        let d_max = 2.0 * MAX_CRATER_RADIUS_FRACTION * inputs.body_radius;
        Self {
            rng: Rng::new(inputs.seed ^ 0x1D2A_9F0B),
            stamper: Stamper::new(&grid),
            inputs,
            grid,
            d_loop,
            d_max,
            records: Vec::new(),
        }
    }

    pub fn loop_threshold(&self) -> f64 {
        self.d_loop
    }

    /// Measure the relief every emplaced crater still has on the final coarse surface,
    /// mark the recognisable ones, and return the records.
    pub fn into_records(mut self, layers: &LayerSet) -> Vec<CraterRecord> {
        let height = layers.bedrock.as_slice();
        for rec in self.records.iter_mut() {
            let dim = dimensions(rec.diameter, &self.inputs.scaling);
            self.stamper.gather(
                &self.grid,
                self.inputs.body_radius,
                rec.center,
                1.15 * dim.radius,
            );
            let (mut floor, mut n_floor, mut rim, mut n_rim) = (0.0, 0usize, 0.0, 0usize);
            for &(idx, r, phi) in &self.stamper.footprint {
                let x = elliptical_radius(r, phi, rec.azimuth, rec.ellipticity) / dim.radius;
                if x < 0.5 {
                    floor += height[idx] as f64;
                    n_floor += 1;
                } else if (0.9..=1.15).contains(&x) {
                    rim += height[idx] as f64;
                    n_rim += 1;
                }
            }
            if n_floor == 0 || n_rim == 0 {
                continue;
            }
            let relief = rim / n_rim as f64 - floor / n_floor as f64;
            rec.relief = relief as f32;
            rec.visible = relief > VISIBLE_RELIEF_FRACTION * (dim.depth + dim.rim_height);
        }
        self.records
    }

    /// Emplace the craters of one time step, oldest first.
    /// Craters formed between two times, in order of formation, not yet emplaced.
    pub fn sample(&self, t_start: f64, t_end: f64, step_index: u64) -> Vec<CraterRecord> {
        let inp = &self.inputs;
        if inp.flux_scale <= 0.0 || self.d_max <= self.d_loop {
            return Vec::new();
        }
        // Craters that the regime would have removed by the present are not emplaced.
        let cutoff = inp.age - inp.loop_retention_age;
        if t_end <= cutoff {
            return Vec::new();
        }
        let t_start = t_start.max(cutoff);
        let expected = expected_count_above(
            &inp.pf,
            &inp.flux,
            inp.flux_scale,
            inp.body_radius,
            self.d_loop,
            t_start,
            t_end,
        );
        let (count, _) = poisson(&self.rng, step_index, 0, expected);
        let mut batch: Vec<CraterRecord> = (0..count as u64)
            .map(|j| {
                let draw = 16 + j * 8;
                let u_d = self.rng.uniform(step_index, draw);
                let center = self.rng.unit_direction(step_index, draw + 1);
                let time_gyr = inp.flux.time_at_fraction(
                    t_start,
                    t_end,
                    self.rng.uniform(step_index, draw + 3),
                );
                // From the early, basin-forming population or the steady asteroid flux,
                // by their shares at the time.
                let early = self.rng.uniform(step_index, draw + 6) < inp.flux.early_share(time_gyr);
                let diameter = if early {
                    inp.pf
                        .sample_diameter_above(u_d, self.d_loop / 1e3, self.d_max / 1e3)
                } else {
                    inp.pf
                        .sample_steady_diameter_above(u_d, self.d_loop / 1e3, self.d_max / 1e3)
                } * 1e3;
                let angle = impact_angle(self.rng.uniform(step_index, draw + 4));
                let azimuth = self.rng.uniform(step_index, draw + 5) * 2.0 * std::f64::consts::PI;
                let ellipticity = ellipticity_for(angle, inp.oblique_fraction);
                CraterRecord {
                    center,
                    diameter,
                    time_gyr,
                    azimuth,
                    angle,
                    ellipticity,
                    class: inp.scaling.classify(diameter),
                    from_loop: true,
                    relief: 0.0,
                    visible: true,
                    on_lava: false,
                    id: (step_index << 32) | j,
                }
            })
            .collect();
        batch.sort_by(|a, b| a.time_gyr.total_cmp(&b.time_gyr).then(a.id.cmp(&b.id)));
        batch
    }

    /// Emplace one crater on the surface as it is now, and record it.
    pub fn apply(&mut self, layers: &mut LayerSet, mut rec: CraterRecord) {
        let inp = &self.inputs;
        {
            let dim = dimensions(rec.diameter, &inp.scaling);
            let deep = reaches_deep_bedrock(rec.diameter, inp.deep_bedrock_depth);
            rec.on_lava = target_is_lava(layers, &self.grid, rec.center);
            let table = ProfileTable::new(dim, 0.0);
            let params = StampParams {
                time_gyr: rec.time_gyr as f32,
                unit_floor: floor_unit(dim.class, deep, rec.on_lava),
                unit_ejecta: ejecta_unit(deep, rec.on_lava),
                azimuth: rec.azimuth,
                angle: rec.angle,
                ellipticity: rec.ellipticity,
                amplitude: 1.0,
                replace_floor: true,
                set_age: true,
                bury: true,
                seed: rec.id,
                deep,
            };
            stamp(
                layers,
                &self.grid,
                inp.body_radius,
                rec.center,
                &table,
                &params,
                &mut self.stamper,
            );
            // Large impacts thin the crust under them, by excavation and the mantle's
            // rise beneath the cavity. Later magma favours these thin spots, so basins
            // flood.
            let thinning = BASIN_CRUST_THINNING * EXCAVATION_DEPTH_FRACTION * rec.diameter;
            let radius = dim.radius * rec.ellipticity.sqrt();
            let crust = layers.crust_thickness.as_mut_slice();
            for &(idx, r, _) in &self.stamper.footprint {
                if r < radius {
                    let w = 1.0 - (r / radius).powi(2);
                    crust[idx] = (crust[idx] as f64 - thinning * w).max(CRUST_THICKNESS_MIN) as f32;
                }
            }
        }
        self.records.push(rec);
    }

    /// Sample and emplace a step's craters.
    pub fn step(&mut self, layers: &mut LayerSet, t_start: f64, t_end: f64, step_index: u64) {
        for rec in self.sample(t_start, t_end, step_index) {
            self.apply(layers, rec);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::presets::preset;
    use crate::stages::primordial_crust::PrimordialCrust;

    /// Over many seeds the loop population must follow the production function: the share
    /// of craters above twice the threshold, the number formed in the last 1 Gyr, and the
    /// rarity of young giant basins.
    #[test]
    fn loop_population_follows_the_production_function() {
        let (mut g, _) = preset("Moon").unwrap().to_graph();
        g.set("detail.quality", "preview").unwrap();
        g.set("detail.loop_resolution", 64.0).unwrap();
        let grid = CubeSphere::new(64);
        let age = g.f64("body.age");
        let steps = 30usize;
        let dt = age / steps as f64;
        let seeds = 60u64;
        let (mut total, mut above_twice, mut young, mut young_giant) =
            (0usize, 0usize, 0usize, 0usize);
        let mut d_loop = 0.0;
        for s in 0..seeds {
            g.set("seed.master", 1_000.0 + s as f64).unwrap();
            let mut layers = PrimordialCrust::new(&g, grid).init();
            let mut process = ImpactProcess::new(&g, grid);
            d_loop = process.loop_threshold();
            for k in 1..=steps {
                process.step(&mut layers, dt * (k - 1) as f64, dt * k as f64, k as u64);
            }
            total += process.records.len();
            above_twice += process
                .records
                .iter()
                .filter(|c| c.diameter >= 2.0 * d_loop)
                .count();
            young += process
                .records
                .iter()
                .filter(|c| age - c.time_gyr < 1.0)
                .count();
            young_giant += process
                .records
                .iter()
                .filter(|c| age - c.time_gyr < 1.0 && c.diameter >= 700_000.0)
                .count();
        }
        let inp = CraterInputs::from_graph(&g);
        let frac_expected =
            inp.pf.cumulative_1gyr(2.0 * d_loop / 1e3) / inp.pf.cumulative_1gyr(d_loop / 1e3);
        let frac = above_twice as f64 / total as f64;
        assert!(
            (frac - frac_expected).abs() < 0.2 * frac_expected + 3.0 / (total as f64).sqrt(),
            "share above twice the threshold {frac} expected {frac_expected} from {total} craters"
        );
        let young_expected = seeds as f64
            * expected_count_above(
                &inp.pf,
                &inp.flux,
                inp.flux_scale,
                inp.body_radius,
                d_loop,
                age - 1.0,
                age,
            );
        assert!(
            (young as f64 - young_expected).abs()
                < 0.25 * young_expected + 3.0 * young_expected.sqrt(),
            "young craters {young} expected {young_expected}"
        );
        let giant_expected = seeds as f64
            * expected_count_above(
                &inp.pf,
                &inp.flux,
                inp.flux_scale,
                inp.body_radius,
                700_000.0,
                age - 1.0,
                age,
            );
        assert!(
            (young_giant as f64) <= 2.0 * giant_expected + 4.0,
            "young giant basins {young_giant} expected {giant_expected}"
        );
        // Sanity: the whole-age basin count of a Moon-like body is Moon-like.
        let basins = seeds as f64
            * expected_count_above(
                &inp.pf,
                &inp.flux,
                inp.flux_scale,
                inp.body_radius,
                300_000.0,
                0.0,
                age,
            );
        assert!(
            (basins / seeds as f64 - 43.0).abs() < 12.0,
            "basins above 300 km per Moon: {}",
            basins / seeds as f64
        );
    }
}
