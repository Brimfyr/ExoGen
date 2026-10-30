// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! Impact cratering inside the coarse loop (PLAN.md 4.5, ADR 0002): craters above the
//! loop threshold are sampled from the production function for each time step and
//! stamped onto the coarse layers at the step they form, so every later process acts on
//! them. Each one is recorded for the finishing pass.

use crate::constants::*;
use crate::crater::Levelling;
use crate::crater::morphology::transient_diameter;
use crate::crater::production::{Oblique, Target, expected_count_above, poisson};
use crate::crater::stamp::elliptical_radius;
use crate::crater::{
    Bombardment, BombardmentParams, CraterRecord, DeepShare, ProfileTable, Scaling, StampParams,
    Stamper, dimensions, excavated_below, hold_diameter, impact_angle, stamp,
};
use crate::fields::LayerSet;
use crate::grid::CubeSphere;
use crate::params::ParamGraph;
use crate::rng::Rng;
use crate::units::{
    UNIT_BASIN_FLOOR, UNIT_CRATER_FLOOR, UNIT_DEEP_EJECTA, UNIT_DEEP_FLOOR, UNIT_EJECTA,
    UNIT_LAVA_CRATER_FLOOR, UNIT_LAVA_EJECTA,
};

pub const VERSION: u64 = 35;

/// Body-level cratering inputs read once from the graph.
#[derive(Clone, Debug)]
pub struct CraterInputs {
    /// The populations striking the body, each with its production and its history in time
    /// since formation, the record starting at `cratering.record_start` (F135).
    pub bombardment: Bombardment,
    pub scaling: Scaling,
    pub body_radius: f64,
    pub age: f64,
    /// The threshold angles of elliptical craters (F167).
    pub oblique: Oblique,
    /// The crust a mega-impact's melt sheet makes on its floor, m, and the height the
    /// ground stands per metre of crust (F168).
    pub mega_floor_crust: f64,
    pub airy: f64,
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
    /// Thickness of the surface layer above the deep bedrock, metres; and, where the deep
    /// bedrock is a layer of the crust, the crust's thickness at which it lies that deep, m
    /// (F215).
    pub deep_bedrock_depth: f64,
    pub deep_crust_reference: Option<f64>,
    /// The exposure on which the surface's maturity saturates, Gyr, for mixing thin
    /// ejecta with the ground it covers (F51).
    pub maturity_timescale: f64,
    /// The slopes, radians, from which the ground sheds what lies loose on it and at which
    /// it sheds it all (F98): an icy surface's lag from gentler slopes than a silicate
    /// regolith.
    pub shed: (f64, f64),
    pub seed: u64,
    /// The width over which an ejecta deposit levels the ground, by its thickness (F169).
    pub levelling: Levelling,
}

impl CraterInputs {
    /// What a crater of this diameter, m, brings up from the deep bedrock (F215, F216),
    /// where the crust under its centre is this thick, m: its ejecta by the share of its
    /// excavation cavity below it, the cavity reaching EXCAVATION_DEPTH_FRACTION of its
    /// transient cavity's diameter, and its floor by the share of the transient cavity's
    /// whole depth below it, TRANSIENT_DEPTH_FRACTION (Melosh 1989). A deep bedrock that is
    /// the crust's own lower layer lies as deep as the crust there is thick against the
    /// reference (the lunar central peaks more mafic where the crust is thinner, Cahill et
    /// al. 2009); an ice shell's dirty skin is its own depth.
    pub fn deep_share(&self, diameter: f64, crust: f64) -> DeepShare {
        let depth = match self.deep_crust_reference {
            Some(reference) => self.deep_bedrock_depth * crust.max(0.0) / reference.max(1.0),
            None => self.deep_bedrock_depth,
        };
        let (excavation, cavity) = self.cavity(diameter);
        DeepShare {
            floor: excavated_below(cavity, depth),
            ejecta: excavated_below(excavation, depth),
        }
    }

    /// The depths a crater of this diameter, m, excavates to and its transient cavity
    /// reaches, m (Melosh 1989).
    pub fn cavity(&self, diameter: f64) -> (f64, f64) {
        let transient = transient_diameter(diameter, &self.scaling);
        (
            EXCAVATION_DEPTH_FRACTION * transient,
            TRANSIENT_DEPTH_FRACTION * transient,
        )
    }
}

/// The crust's thickness under a point, m.
pub fn crust_at(layers: &LayerSet, grid: &CubeSphere, center: [f64; 3]) -> f64 {
    layers.crust_thickness.as_slice()[grid.index(grid.locate(center).0)] as f64
}

/// Unit of a crater floor by class, excavation and target.
pub fn floor_unit(class: crate::crater::CraterClass, deep: bool, on_lava: bool) -> u16 {
    // A crater in lava is a lava crater whatever it brings up from the deep bedrock, which
    // shows by its share (F222).
    if on_lava {
        UNIT_LAVA_CRATER_FLOOR
    } else if deep {
        UNIT_DEEP_FLOOR
    } else if class == crate::crater::CraterClass::MultiRing {
        UNIT_BASIN_FLOOR
    } else {
        UNIT_CRATER_FLOOR
    }
}

/// Unit of continuous ejecta by excavation and target.
pub fn ejecta_unit(deep: bool, on_lava: bool) -> u16 {
    if on_lava {
        UNIT_LAVA_EJECTA
    } else if deep {
        UNIT_DEEP_EJECTA
    } else {
        UNIT_EJECTA
    }
}

/// Whether the ground at a point is lava plains, so a crater there is made of them.
pub fn target_is_lava(layers: &LayerSet, grid: &CubeSphere, center: [f64; 3]) -> bool {
    let idx = grid.index(grid.locate(center).0);
    crate::stages::volcanism::is_lava_unit(layers.unit_id.as_slice()[idx])
}

/// The threshold angles of the body's elliptical craters (F167): its crater scaling at its
/// mean impact speed, from the graph.
pub fn oblique(graph: &ParamGraph) -> Oblique {
    Oblique::new(&Target {
        gravity: graph.f64("body.gravity"),
        speed: graph.f64("cratering.impact_velocity"),
        strength: if crate::stages::materials::is_ice(graph.str("materials.bedrock")) {
            TARGET_STRENGTH_ICE
        } else {
            TARGET_STRENGTH_ROCK
        },
        density: graph.f64("cratering.target_density"),
        simple_complex: graph.f64("cratering.simple_complex_transition"),
    })
}

/// The populations striking the body (F135, F136), from the graph.
pub fn bombardment(graph: &ParamGraph) -> Bombardment {
    BombardmentParams::from_nodes(
        |id| graph.f64(id),
        graph.str("cratering.production_function"),
        crate::stages::materials::is_ice(graph.str("materials.bedrock")),
        graph.f64("cratering.record_start"),
        true,
    )
    .build()
}

impl CraterInputs {
    pub fn from_graph(graph: &ParamGraph) -> Self {
        let mut inputs = Self {
            bombardment: bombardment(graph),
            scaling: Scaling {
                simple_complex: graph.f64("cratering.simple_complex_transition"),
                peak_ring: graph.f64("cratering.peak_ring_onset"),
                multi_ring: graph.f64("cratering.multiring_onset"),
                gravity: graph.f64("body.gravity"),
                icy: crate::stages::materials::is_ice(graph.str("materials.bedrock")),
            },
            body_radius: graph.f64("body.radius"),
            age: graph.f64("body.age"),
            oblique: oblique(graph),
            mega_floor_crust: MEGA_FLOOR_CRUST_SHARE
                * graph.f64("tectonics.initial_crust_thickness"),
            airy: crate::stages::primordial_crust::airy_factor(graph),
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
            deep_crust_reference: (!crate::stages::materials::is_ice(
                graph.str("materials.deep_bedrock"),
            ))
            .then(|| graph.f64("tectonics.initial_crust_thickness")),
            maturity_timescale: crate::stages::ice::maturity_timescale(graph),
            shed: if crate::stages::ice::is_ice_surface(graph) {
                (
                    ICE_LAG_SHED_START_DEG.to_radians(),
                    ICE_LAG_SHED_FULL_DEG.to_radians(),
                )
            } else {
                (
                    SLOPE_SHED_START_DEG.to_radians(),
                    REPOSE_SLOPE_DEG.to_radians(),
                )
            },
            seed: graph.f64("seed.cratering") as u64,
            levelling: Levelling::rings(),
        };
        inputs.levelling = Levelling::new(&inputs.scaling);
        inputs
    }
}

/// Regimes that renew the whole surface within the retention age: plates, global
/// resurfacing, or an ice shell kept young by tidal heating.
pub fn wholesale_resurfacing(graph: &ParamGraph) -> bool {
    match graph.str("interior.tectonic_regime") {
        "mobile_lid" | "episodic_lid" => true,
        "ice_shell_ocean" => {
            graph.f64("interior.tidal_flux") > ICE_RESURFACING_TIDAL_FLUX
                && graph.f64("interior.ice_shell_thickness") / graph.f64("body.radius").max(1.0)
                    < WHOLESALE_SHELL_RATIO
        }
        _ => crate::stages::io::is_io_regime(graph),
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
        // The largest crater the body holds (F135), any larger impact having come before its
        // record; those wider than the stamp draws are mega-impacts (F168).
        let d_max = hold_diameter(inputs.body_radius, inputs.scaling.icy);
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

    /// Craters formed between two times, in order of formation, not yet emplaced: drawn
    /// per fixed time bin and population with the bin as the random index, so the
    /// population does not depend on the step length.
    pub fn sample(&self, t_start: f64, t_end: f64) -> Vec<CraterRecord> {
        let inp = &self.inputs;
        if inp.bombardment.sources.is_empty() || self.d_max <= self.d_loop {
            return Vec::new();
        }
        // Craters that the regime would have removed by the present are not emplaced.
        let cutoff = inp.age - inp.loop_retention_age;
        let mut batch: Vec<CraterRecord> = Vec::new();
        for (bin, b0, b1) in crate::stages::time_bins(inp.age, t_start, t_end) {
            if b1 <= cutoff {
                continue;
            }
            let b0 = b0.max(cutoff);
            for source in &inp.bombardment.sources {
                // Each population draws from its own stream, the rocky leftovers' as before.
                let stream = match source.pf.impactors.map(|i| i.sizes) {
                    None => 0u64,
                    Some(crate::crater::ImpactorSizes::Comets) => 1,
                    Some(crate::crater::ImpactorSizes::Debris) => 2,
                };
                let index = bin | (stream << 40);
                let expected = expected_count_above(
                    &source.pf,
                    &source.flux,
                    source.scale,
                    inp.body_radius,
                    self.d_loop,
                    b0,
                    b1,
                );
                let (count, _) = poisson(&self.rng, index, 0, expected);
                for j in 0..count as u64 {
                    let draw = 16 + j * 8;
                    let u_d = self.rng.uniform(index, draw);
                    let center = self.rng.unit_direction(index, draw + 1);
                    let time_gyr =
                        source
                            .flux
                            .time_at_fraction(b0, b1, self.rng.uniform(index, draw + 3));
                    // From the early, basin-forming population or the steady flux, by their
                    // shares at the time.
                    let early =
                        self.rng.uniform(index, draw + 6) < source.flux.early_share(time_gyr);
                    let (lo, hi) = (self.d_loop / 1e3, self.d_max / 1e3);
                    let diameter = if early {
                        source.pf.sample_diameter_above(u_d, lo, hi)
                    } else {
                        source.pf.sample_steady_diameter_above(u_d, lo, hi)
                    } * 1e3;
                    let angle = impact_angle(self.rng.uniform(index, draw + 4));
                    let azimuth = self.rng.uniform(index, draw + 5) * 2.0 * std::f64::consts::PI;
                    let ellipticity = inp.oblique.ellipticity(diameter, angle);
                    batch.push(CraterRecord {
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
                        id: (bin << 32) | (stream << 28) | j,
                    });
                }
            }
        }
        batch.sort_by(|a, b| a.time_gyr.total_cmp(&b.time_gyr).then(a.id.cmp(&b.id)));
        batch
    }

    /// Emplace one crater on the surface as it is now, and record it.
    pub fn apply(&mut self, layers: &mut LayerSet, mut rec: CraterRecord) {
        let inp = &self.inputs;
        if crate::stages::mega::is_mega(rec.diameter, inp.body_radius) {
            crate::stages::mega::emplace(
                layers,
                &self.grid,
                inp.body_radius,
                &rec,
                inp.mega_floor_crust,
                inp.airy,
                &inp.levelling,
            );
            self.records.push(rec);
            return;
        }
        {
            let dim = dimensions(rec.diameter, &inp.scaling);
            let deep = inp.deep_share(rec.diameter, crust_at(layers, &self.grid, rec.center));
            rec.on_lava = target_is_lava(layers, &self.grid, rec.center);
            let table = ProfileTable::new(dim, 0.0);
            let params = StampParams {
                time_gyr: rec.time_gyr as f32,
                unit_floor: floor_unit(dim.class, deep.floor > 0.0, rec.on_lava),
                unit_ejecta: ejecta_unit(deep.ejecta > 0.0, rec.on_lava),
                azimuth: rec.azimuth,
                angle: rec.angle,
                ellipticity: rec.ellipticity,
                amplitude: 1.0,
                feature_amplitude: 1.0,
                replace_floor: true,
                set_age: true,
                bury: true,
                levelling: inp.levelling,
                seed: rec.id,
                deep,
                cavity: inp.cavity(rec.diameter),
                exhume: false,
                maturity_timescale: inp.maturity_timescale,
                renewal_guard: true,
                terraced: false,
                shed: inp.shed,
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
            // flood. The crust they excavate and push aside lands beyond the rim, thinning
            // outward as their ejecta does, so the crust is kept, thickest in an annulus
            // about the basin (F168; the lunar basins' annuli, Wieczorek et al. 2013).
            let thinning = BASIN_CRUST_THINNING * EXCAVATION_DEPTH_FRACTION * rec.diameter;
            let radius = dim.radius * rec.ellipticity.sqrt();
            let area = |idx: usize| self.grid.cell_area(self.grid.cell(idx));
            let crust = layers.crust_thickness.as_mut_slice();
            let mut volume = 0.0;
            for &(idx, r, _) in &self.stamper.footprint {
                if r < radius {
                    let w = 1.0 - (r / radius).powi(2);
                    let old = crust[idx] as f64;
                    let new = (old - thinning * w).max(CRUST_THICKNESS_MIN);
                    crust[idx] = new as f32;
                    volume += (old - new) * area(idx);
                }
            }
            let spread: f64 = self
                .stamper
                .footprint
                .iter()
                .filter(|&&(_, r, _)| r >= radius)
                .map(|&(idx, r, _)| (r / radius).powf(dim.ejecta_exponent) * area(idx))
                .sum();
            if volume > 0.0 && spread > 0.0 {
                for &(idx, r, _) in &self.stamper.footprint {
                    if r >= radius {
                        crust[idx] += (volume * (r / radius).powf(dim.ejecta_exponent) / spread)
                            as f32;
                    }
                }
            }
        }
        self.records.push(rec);
    }

    /// Sample and emplace a step's craters.
    pub fn step(&mut self, layers: &mut LayerSet, t_start: f64, t_end: f64, _step_index: u64) {
        for rec in self.sample(t_start, t_end) {
            self.apply(layers, rec);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::presets::preset;
    use crate::stages::primordial_crust::PrimordialCrust;

    /// F215, F216: a crater brings up the deep bedrock by the share of its excavation cavity
    /// below it, none until it reaches it, and shows it in its floor by the share of its
    /// transient cavity below it, sooner. The Moon's basins throw out mostly its upper crust,
    /// Orientale's size a quarter or so of its lower, and more where the crust is thinner;
    /// Mercury's craters of 150 km none of its deep graphite crust but floor themselves in
    /// it, and Caloris's size throws out most; an ice shell's craters of 20 km about two
    /// fifths of their clean ice. The deep share lands in a ring nearest the rim.
    #[test]
    fn a_crater_brings_up_the_deep_bedrock_by_the_cavity_below_it() {
        assert_eq!(excavated_below(10.0, 20.0), 0.0);
        assert!((excavated_below(20.0, 10.0) - 0.25).abs() < 1e-12);
        let inputs = |name: &str| {
            let (g, _) = preset(name).unwrap().to_graph();
            (
                CraterInputs::from_graph(&g),
                g.f64("tectonics.initial_crust_thickness"),
            )
        };
        let (moon, crust) = inputs("Moon");
        assert_eq!(moon.deep_share(300e3, crust).ejecta, 0.0);
        assert!(moon.deep_share(300e3, crust).floor > 0.2);
        let orientale = moon.deep_share(930e3, crust);
        assert!(
            orientale.ejecta > 0.2 && orientale.ejecta < 0.5,
            "Orientale {orientale:?}"
        );
        assert!(orientale.floor > orientale.ejecta);
        assert!(moon.deep_share(930e3, 0.75 * crust).ejecta > orientale.ejecta);
        let (mercury, crust) = inputs("Mercury");
        let crater = mercury.deep_share(150e3, crust);
        assert!(
            crater.ejecta == 0.0 && crater.floor > 0.3,
            "150 km {crater:?}"
        );
        let caloris = mercury.deep_share(1550e3, crust).ejecta;
        assert!(caloris > 0.6, "Caloris {caloris}");
        let (ganymede, crust) = inputs("Ganymede");
        let twenty = ganymede.deep_share(20e3, crust);
        assert!(
            twenty.ejecta > 0.3 && twenty.ejecta < 0.55,
            "a 20 km crater {twenty:?}"
        );
        assert_eq!(ganymede.deep_share(20e3, 0.5 * crust), twenty);
        // Half the excavated volume deep: it lands within two radii, all of it near the rim
        // and none past the ring's soft edge.
        let half = DeepShare {
            floor: 0.8,
            ejecta: 0.5,
        };
        assert_eq!(half.ejecta_at(1.2), 1.0);
        assert!((half.ejecta_at(2.0) - 0.5).abs() < 1e-9);
        assert_eq!(half.ejecta_at(2.5), 0.0);
        assert_eq!(DeepShare::NONE.ejecta_at(1.0), 0.0);
        // The floor's share inside its edge grades up the wall to the ring's at the rim,
        // with no step there, and the ring follows the ejecta's lobes.
        assert_eq!(half.at(0.3, 0.5, 1.0), 0.8);
        assert!((half.at(1.0, 0.5, 1.0) - half.at(1.0 + 1e-9, 0.5, 1.0)).abs() < 1e-6);
        let wall = half.at(0.75, 0.5, 1.0);
        assert!(wall > 0.8 && wall < 1.0, "the wall {wall}");
        assert!(half.at(2.0, 0.5, 1.2) > half.at(2.0, 0.5, 1.0));
        assert!(half.at(2.0, 0.5, 0.8) < half.at(2.0, 0.5, 1.0));
    }

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
        let b = &inp.bombardment;
        let frac_expected = b.cumulative_1gyr(2.0 * d_loop / 1e3) / b.cumulative_1gyr(d_loop / 1e3);
        let frac = above_twice as f64 / total as f64;
        assert!(
            (frac - frac_expected).abs() < 0.2 * frac_expected + 3.0 / (total as f64).sqrt(),
            "share above twice the threshold {frac} expected {frac_expected} from {total} craters"
        );
        let young_expected =
            seeds as f64 * b.expected_count_above(inp.body_radius, d_loop, age - 1.0, age);
        assert!(
            (young as f64 - young_expected).abs()
                < 0.25 * young_expected + 3.0 * young_expected.sqrt(),
            "young craters {young} expected {young_expected}"
        );
        let giant_expected =
            seeds as f64 * b.expected_count_above(inp.body_radius, 700_000.0, age - 1.0, age);
        assert!(
            (young_giant as f64) <= 2.0 * giant_expected + 4.0,
            "young giant basins {young_giant} expected {giant_expected}"
        );
        // Sanity: the whole-age basin count of a Moon-like body is Moon-like.
        let basins = seeds as f64 * b.expected_count_above(inp.body_radius, 300_000.0, 0.0, age);
        assert!(
            (basins / seeds as f64 - 43.0).abs() < 12.0,
            "basins above 300 km per Moon: {}",
            basins / seeds as f64
        );
    }

    /// Mimas, from Enceladus's preset with Mimas's radius, density and orbit.
    fn mimas() -> crate::params::ParamGraph {
        saturnian(1.982e5, 1149.0, 1.8552e8, 0.0196)
    }

    /// A moon of Saturn, from Enceladus's preset with its own radius, density and orbit.
    fn saturnian(radius: f64, density: f64, a: f64, e: f64) -> crate::params::ParamGraph {
        let (mut g, _) = preset("Enceladus").unwrap().to_graph();
        g.set("body.radius", radius).unwrap();
        g.set("body.density", density).unwrap();
        g.set("orbit.semi_major_axis", a).unwrap();
        g.set("orbit.eccentricity", e).unwrap();
        g
    }

    /// The comets the outer system's moons meet (F135): their speeds and numbers per unit
    /// area today against Nesvorny et al. 2023's, whose rates are Jupiter's 1 km impactor
    /// every 230 years times each moon's share.
    #[test]
    fn the_outer_moons_meet_their_comets_at_their_speeds_and_rates() {
        let jupiter_per_year = 1.0 / 230.0;
        for (name, speed_kms, share) in [
            ("Io", 31.6, 1.1e-4),
            ("Europa", 25.4, 5.9e-5),
            ("Ganymede", 20.3, 9.7e-5),
            ("Callisto", 16.0, 5.7e-5),
            ("Enceladus", 23.3, 1.5e-6),
            ("Miranda", 12.7, 3.5e-6),
        ] {
            let (g, _) = preset(name).unwrap().to_graph();
            let speed = g.f64("cratering.impact_velocity") / 1e3;
            assert!(
                (speed / speed_kms - 1.0).abs() < 0.1,
                "{name}: impact speed {speed:.1} km/s, Nesvorny et al. 2023 {speed_kms}"
            );
            let r = g.f64("body.radius") / 1e3;
            let rate = jupiter_per_year * share / (4.0 * std::f64::consts::PI * r * r) * 1e9;
            let flux = g.f64("cratering.icy_flux");
            assert!(
                (flux / rate).ln().abs() < 1.5f64.ln(),
                "{name}: comets of 1 km {flux:.2e} per km² per Gyr, Nesvorny et al. 2023 {rate:.2e}"
            );
            assert_eq!(g.f64("cratering.flux_relative_to_moon"), 0.0, "{name}");
        }
    }

    /// The moons keep the craters they were counted with (F136): Kirchoff and Schenk's counts
    /// as Bottke et al. 2024 replot them (Fig. 7 for Saturn's moons, Fig. 9 for Miranda),
    /// read off the figures, per km² from the record's start. Mimas, the least altered of
    /// the small moons, calibrates the debris its last unholdable impact rains on it, which
    /// makes its craters under 20 km; Tethys, Dione and Miranda check it; Rhea, whose record
    /// starts with its crust, takes no debris, and keeps the comets' craters at 0.6 of its
    /// counts.
    #[test]
    fn the_moons_keep_the_craters_they_were_counted_with() {
        let (miranda, _) = preset("Miranda").unwrap().to_graph();
        let cases = [
            (
                "Mimas",
                mimas(),
                vec![
                    (5.0, 4.47e-3),
                    (10.0, 1.37e-3),
                    (20.0, 2.79e-4),
                    (50.0, 2.3e-5),
                ],
                1.3,
            ),
            (
                "Tethys",
                saturnian(5.311e5, 984.0, 2.94619e8, 0.0001),
                vec![
                    (10.0, 7.4e-4),
                    (20.0, 1.6e-4),
                    (50.0, 1.45e-5),
                    (100.0, 3.9e-6),
                ],
                1.5,
            ),
            (
                "Dione",
                saturnian(5.614e5, 1478.0, 3.77396e8, 0.0022),
                vec![
                    (10.0, 6.9e-4),
                    (20.0, 1.29e-4),
                    (50.0, 1.17e-5),
                    (100.0, 2.3e-6),
                ],
                1.5,
            ),
            (
                "Miranda",
                miranda,
                vec![(5.0, 1.75e-3), (10.0, 7.16e-4), (20.0, 1.74e-4)],
                1.5,
            ),
            (
                "Rhea",
                saturnian(7.638e5, 1236.0, 5.27108e8, 0.001),
                vec![(10.0, 6.0e-4), (20.0, 1.48e-4), (50.0, 1.86e-5)],
                2.0,
            ),
        ];
        for (name, g, counts, within) in cases {
            let inp = CraterInputs::from_graph(&g);
            let age = g.f64("body.age");
            let start = g.f64("cratering.record_start");
            for (d, counted) in counts {
                let n = inp.bombardment.density_above(d, start, age);
                assert!(
                    (n / counted).ln().abs() < f64::ln(within),
                    "{name}: N({d} km) {n:.2e} per km², counted {counted:.2e}"
                );
            }
            let debris = g.f64("cratering.debris_fragments");
            assert_eq!(debris > 0.0, name != "Rhea", "{name}: debris {debris}");
        }
        let g = mimas();
        let formed_ga = g.f64("body.age") - g.f64("cratering.record_start");
        assert!(
            (4.0..=4.35).contains(&formed_ga),
            "Mimas's record from {formed_ga:.2} Ga"
        );
    }

    /// A small ice moon holds no crater wider than itself (F135): its record starts at the
    /// last such impact's median time, when ln 2 are still to come, and the loop draws none;
    /// the large moons keep their record from when their crust solidified, the Moon from
    /// within a few tens of millions of years of it.
    #[test]
    fn a_small_ice_moon_holds_no_crater_wider_than_itself() {
        let mut small = vec![("Mimas".to_string(), mimas())];
        for name in ["Enceladus", "Miranda"] {
            small.push((name.to_string(), preset(name).unwrap().to_graph().0));
        }
        for (name, g) in &small {
            let process = ImpactProcess::new(g, CubeSphere::new(16));
            let r = g.f64("body.radius");
            assert!((process.d_max / r - 1.0).abs() < 1e-9, "{name}");
            let start = g.f64("cratering.record_start");
            assert!(start > g.f64("surface.solidification_time"), "{name}");
            let after = process.inputs.bombardment.expected_count_above(
                r,
                process.d_max,
                start,
                g.f64("body.age"),
            );
            assert!(
                (after - std::f64::consts::LN_2).abs() < 0.01,
                "{name}: {after} impacts it could not hold"
            );
        }
        for name in ["Mercury", "Ganymede", "Callisto"] {
            let (g, _) = preset(name).unwrap().to_graph();
            assert_eq!(
                g.f64("cratering.record_start"),
                g.f64("surface.solidification_time"),
                "{name}"
            );
        }
        let (g, _) = preset("Moon").unwrap().to_graph();
        let late = g.f64("cratering.record_start") - g.f64("surface.solidification_time");
        assert!(
            (0.0..0.05).contains(&late),
            "the Moon's record {late} Gyr after solidifying"
        );
    }

    /// Around a red dwarf the leftovers clear on its short orbital times (F135): a planet
    /// where the Earth's insolation falls two thirds as strongly, as TRAPPIST-1e, meets them
    /// fast and keeps a record of the late flux alone.
    #[test]
    fn around_a_red_dwarf_the_early_bombardment_is_brief() {
        let (mut g, _) = preset("Earth").unwrap().to_graph();
        g.set("star.mass", 0.09).unwrap();
        g.set("star.luminosity", 5.5e-4).unwrap();
        g.set("orbit.star_distance", 0.0289 * crate::constants::AU)
            .unwrap();
        g.set("body.radius", 5.85e6).unwrap();
        let x = g.f64("cratering.snow_line_position");
        assert!((0.4..0.5).contains(&x), "place against the snow line {x}");
        let speed = g.f64("cratering.impact_velocity") / 1e3;
        // Kral et al. 2018: 15 to 110 km/s on TRAPPIST-1's planets.
        assert!((15.0..60.0).contains(&speed), "impact speed {speed} km/s");
        assert!(g.f64("cratering.early_flux_decay") < 0.005);
        let start = g.f64("cratering.record_start");
        assert_eq!(start, g.f64("surface.solidification_time"));
        let age = g.f64("body.age");
        let b = CraterInputs::from_graph(&g).bombardment;
        let equivalent = b.equivalent_gyr(start, age);
        assert!(
            equivalent < 1.05 * (age - start),
            "{equivalent} Gyr of steady flux over {} Gyr",
            age - start
        );
    }

    /// F167: on Mars, a few percent of the craters of tens of kilometres come out elliptical,
    /// and some three in ten of its giant basins (Collins et al. 2011, Figure 2: the model's
    /// share some 3 to 4 percent from 10 to 100 km, near 30 percent at 2,500 km; the observed
    /// share some 6.5 percent from 10 to 40 km), the share rising with the size on every
    /// body; the Moon's craters of 10 km a few percent (Bottke et al. 2000: some 5 percent).
    #[test]
    fn a_few_craters_and_many_giant_basins_are_elliptical() {
        let share = |o: &Oblique, d: f64| o.threshold(d).sin().powi(2);
        let (g, _) = preset("Mars").unwrap().to_graph();
        let mars = oblique(&g);
        for d in [1e4, 3e4] {
            let s = share(&mars, d);
            assert!(s > 0.03 && s < 0.1, "{d}: {s}");
        }
        let s = share(&mars, 2.5e6);
        assert!(s > 0.25 && s < 0.45, "{s}");
        let (g, _) = preset("Moon").unwrap().to_graph();
        let moon = oblique(&g);
        let s = share(&moon, 1e4);
        assert!(s > 0.02 && s < 0.08, "{s}");
        for o in [mars, moon] {
            for k in 0..24 {
                let d = 1e3 * 10f64.powf(0.2 * k as f64);
                assert!(o.threshold(1.6 * d) > o.threshold(d));
            }
        }
    }

}
