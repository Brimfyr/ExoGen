// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! Primordial crust (PLAN.md 4.1): the initial layered state at formation.
//!
//! Crustal thickness is a low-degree spherical-harmonic field, or a two-level crust for
//! mobile-lid bodies until plate tectonics exists; a hemispheric dichotomy comes of the
//! impacts the body takes (F168), not of the crust it starts with. Elevation follows from Airy isostasy, with a seeded self-affine relief from the
//! scale of the body's hemispheres down to the grid, as rough at its largest as a cratered
//! surface is at the scale of the largest crater the body holds (F137), standing in for
//! the relief its history before the record left.

use rayon::prelude::*;

use crate::constants::{
    DENSITY_CRUST_ICE, DENSITY_CRUST_ROCK, FINE_TEXTURE_SHARE, MAX_CRATER_RADIUS_FRACTION,
    ROUGHNESS_HURST, ROUGHNESS_RMS_SHARE,
};
use crate::fields::{Field, LayerSet};
use crate::grid::CubeSphere;
use crate::noise::GradientNoise;
use crate::params::ParamGraph;
use crate::rng::Rng;
use crate::spharm::Harmonics;
use crate::units::{UNIT_CONTINENTAL, UNIT_OCEANIC, UNIT_PRIMORDIAL};

/// Bumped when the behaviour changes, so caches are not reused across versions.
pub const VERSION: u64 = 7;

/// Two-level crust for mobile-lid bodies: thin oceanic crust and thick continental crust
/// covering a fraction of the area, which gives the bimodal hypsometry that lets an ocean
/// fill the basins and leave continents standing. Plate tectonics replaces it in Phase 6.
struct Continents {
    oceanic_thickness: f64,
    continental_thickness: f64,
    area_fraction: f64,
}

/// The height the ground stands per metre of crust, the crust's buoyancy in the mantle
/// beneath it (Airy): rock on rock, or ice on its ocean.
pub fn airy_factor(graph: &ParamGraph) -> f64 {
    let icy = matches!(
        graph.str("body.composition_class"),
        "rock_ice" | "ice_dominated"
    );
    let (rho_mantle, rho_crust) = if icy {
        (1000.0, DENSITY_CRUST_ICE)
    } else {
        (3300.0, DENSITY_CRUST_ROCK)
    };
    (rho_mantle - rho_crust) / rho_mantle
}

pub struct PrimordialCrust {
    grid: CubeSphere,
    noise: GradientNoise,
    harmonics: Harmonics,
    crust_mean: f64,
    crust_sigma: f64,
    continents: Option<Continents>,
    /// The primordial relief's largest octave, rms m, and its features, lattice cells per
    /// radian: the largest crater the body holds and its depth's share (F137).
    relief: f64,
    relief_frequency: f64,
    /// The strength's limit on the relief's whole rms, m.
    relief_cap: f64,
    isostatic_factor: f64,
    /// The relief's scale on its normalized field, m, measured on the loop's grid when
    /// the crust is made (F141): computed eagerly, never under a lock, since the texture
    /// is read from inside the parallel passes.
    relief_scale: f64,
}

impl PrimordialCrust {
    pub fn new(graph: &ParamGraph, grid: CubeSphere) -> Self {
        let seed = graph.f64("seed.terrain") as u64;
        let rng = Rng::new(seed);
        // F137: the relief of a cratered surface at the scale of the largest crater the
        // body holds, a share of its depth, as the history before the record left it;
        // never more than a quarter of the relief the crust's strength supports.
        let scaling = crate::crater::Scaling {
            simple_complex: graph.f64("cratering.simple_complex_transition"),
            peak_ring: graph.f64("cratering.peak_ring_onset"),
            multi_ring: graph.f64("cratering.multiring_onset"),
            gravity: graph.f64("body.gravity"),
            ice: graph.f64("cratering.target_ice"),
        };
        let radius = graph.f64("body.radius");
        let largest = crate::crater::hold_diameter(radius, scaling.ice)
            .min(2.0 * MAX_CRATER_RADIUS_FRACTION * radius);
        let relief = ROUGHNESS_RMS_SHARE
            * crate::crater::dimensions(largest, &scaling).depth
            * graph.f64("multiplier.relief");
        let relief_frequency = radius / largest.max(1.0);
        let relief_cap = 0.25 * graph.f64("interior.max_relief") * graph.f64("multiplier.relief");
        let continents = if graph.str("interior.tectonic_regime") == "mobile_lid" {
            Some(Continents {
                oceanic_thickness: 8_000.0,
                continental_thickness: 38_000.0,
                area_fraction: 0.30,
            })
        } else {
            None
        };
        let mut crust = Self {
            grid,
            noise: GradientNoise::new(seed ^ 0x9E37_79B9),
            harmonics: Harmonics::random(&rng, 2, 8, 2.0),
            crust_mean: graph.f64("tectonics.initial_crust_thickness"),
            crust_sigma: 0.12 * graph.f64("tectonics.initial_crust_thickness"),
            continents,
            relief,
            relief_frequency,
            relief_cap,
            isostatic_factor: airy_factor(graph),
            relief_scale: 0.0,
        };
        let raw = crust.raw_relief();
        let rms = (raw.iter().map(|v| v * v).sum::<f64>() / raw.len().max(1) as f64).sqrt();
        crust.relief_scale = crust.relief_scale_for(rms);
        crust
    }

    /// The relief's octaves on the loop's grid: from the largest crater the body holds
    /// down to two of its cells, or the octave past that (F137).
    fn relief_octaves(&self, grid: &CubeSphere) -> u32 {
        let finest = grid.face_size() as f64 / std::f64::consts::PI;
        ((finest / self.relief_frequency).log2().ceil() as u32 + 1).max(1)
    }

    /// The relief's normalized field on the loop's grid, as `fbm` sums it.
    fn raw_relief(&self) -> Vec<f64> {
        let grid = self.grid;
        let gain = 2f64.powf(-ROUGHNESS_HURST);
        let octaves = self.relief_octaves(&grid);
        (0..grid.cell_count())
            .into_par_iter()
            .map(|k| {
                self.noise.fbm(
                    grid.center(grid.cell(k)),
                    self.relief_frequency,
                    octaves,
                    2.0,
                    gain,
                )
            })
            .collect()
    }

    /// The scale on the normalized field, m: the largest octave's rms is the relief's and
    /// the others fall by the gain, the whole no more than the strength's limit, over the
    /// field's rms on the loop's grid.
    fn relief_scale_for(&self, raw_rms: f64) -> f64 {
        let gain = 2f64.powf(-ROUGHNESS_HURST);
        let total = self.relief
            * (0..self.relief_octaves(&self.grid))
                .map(|o| gain.powi(2 * o as i32))
                .sum::<f64>()
                .sqrt();
        total.min(self.relief_cap) / raw_rms.max(1e-12)
    }

    fn relief_scale(&self) -> f64 {
        self.relief_scale
    }

    /// The relief's finest octave on the loop's grid, the last init laid.
    pub fn loop_finest(&self) -> u32 {
        self.relief_octaves(&self.grid) - 1
    }

    /// The relief's finest octave a grid draws cleanly: the last whose wavelength is at
    /// least two of its cells.
    pub fn finest_on(&self, grid: &CubeSphere) -> u32 {
        let finest = grid.face_size() as f64 / std::f64::consts::PI;
        (finest / self.relief_frequency).log2().floor().max(0.0) as u32
    }

    /// The relief's octaves `first` to `last` at a point, m, each at the amplitude the
    /// loop's grid gave it (F141), weighted by how far below `max_wavelength` m it lies:
    /// in full at half that wavelength and below, nothing at it, between by the logarithm
    /// of the wavelength, so that a terrain levelling over a length keeps the relief below
    /// that length and two terrains of nearly the same length keep nearly the same; and
    /// below `fabric` m, the spacing of the terrain's own grooves, faults or folds, toned
    /// down to FINE_TEXTURE_SHARE over an octave, the fabric overprinting what is finer
    /// than itself (0 for no fabric). Octave o has the wavelength R / (f 2^o), f the
    /// relief's frequency in lattice cells per radian.
    pub fn texture_at(
        &self,
        d: [f64; 3],
        max_wavelength: f64,
        fabric: f64,
        body_radius: f64,
        first: u32,
        last: u32,
    ) -> f64 {
        let gain = 2f64.powf(-ROUGHNESS_HURST);
        let norm: f64 = (0..self.relief_octaves(&self.grid))
            .map(|o| gain.powi(o as i32))
            .sum();
        let scale = self.relief_scale() / norm.max(1e-12);
        let mut sum = 0.0;
        for o in first..=last {
            let wavelength = body_radius / (self.relief_frequency * 2f64.powi(o as i32));
            let mut weight = (max_wavelength / wavelength).log2().clamp(0.0, 1.0);
            if weight <= 0.0 {
                continue;
            }
            if fabric > 0.0 {
                let kept = ((wavelength / fabric).log2() + 1.0).clamp(0.0, 1.0);
                weight *= FINE_TEXTURE_SHARE + (1.0 - FINE_TEXTURE_SHARE) * kept;
            }
            sum += weight
                * self
                    .noise
                    .octaves(d, self.relief_frequency, o, 1, 2.0, gain);
        }
        scale * sum
    }

    /// `texture_at` over a grid.
    pub fn texture_field(
        &self,
        grid: &CubeSphere,
        max_wavelength: f64,
        fabric: f64,
        body_radius: f64,
        first: u32,
        last: u32,
    ) -> Vec<f32> {
        if last < first {
            return vec![0.0; grid.cell_count()];
        }
        (0..grid.cell_count())
            .into_par_iter()
            .map(|k| {
                self.texture_at(
                    grid.center(grid.cell(k)),
                    max_wavelength,
                    fabric,
                    body_radius,
                    first,
                    last,
                ) as f32
            })
            .collect()
    }

    /// The relief below `max_wavelength` m on a grid, down to two of its cells (F141).
    pub fn fine_relief(
        &self,
        grid: &CubeSphere,
        max_wavelength: f64,
        body_radius: f64,
    ) -> Vec<f32> {
        self.texture_field(
            grid,
            max_wavelength,
            0.0,
            body_radius,
            0,
            self.finest_on(grid),
        )
    }

    /// The relief's octaves below the loop's cells on a finer grid (F141), under a
    /// terrain whose fabric is spaced `fabric` m: what the loop's grid could not hold,
    /// which the floors it laid lack.
    pub fn below_loop(&self, grid: &CubeSphere, fabric: f64, body_radius: f64) -> Vec<f32> {
        self.texture_field(
            grid,
            f64::INFINITY,
            fabric,
            body_radius,
            self.loop_finest() + 1,
            self.finest_on(grid),
        )
    }

    /// Initial layered state at formation.
    pub fn init(&self) -> LayerSet {
        let grid = self.grid;
        let mut layers = LayerSet::new(grid);

        let raw: Vec<f64> = (0..grid.cell_count())
            .into_par_iter()
            .map(|k| self.harmonics.eval(grid.center(grid.cell(k))))
            .collect();
        let rms = (raw.iter().map(|v| v * v).sum::<f64>() / raw.len() as f64)
            .sqrt()
            .max(1e-12);

        let crust: Vec<f32> = match &self.continents {
            Some(cont) => {
                let mut sorted = raw.clone();
                sorted.sort_by(|a, b| a.total_cmp(b));
                let thr = sorted[((1.0 - cont.area_fraction) * (sorted.len() - 1) as f64) as usize];
                let width = 0.35 * rms;
                (0..grid.cell_count())
                    .into_par_iter()
                    .map(|k| {
                        let s = ((raw[k] - thr) / width + 0.5).clamp(0.0, 1.0);
                        let smooth = s * s * (3.0 - 2.0 * s);
                        (cont.oceanic_thickness
                            + (cont.continental_thickness - cont.oceanic_thickness) * smooth)
                            as f32
                    })
                    .collect()
            }
            None => (0..grid.cell_count())
                .into_par_iter()
                .map(|k| (self.crust_mean + self.crust_sigma * raw[k] / rms) as f32)
                .collect(),
        };
        let mean_crust = crust.iter().map(|v| *v as f64).sum::<f64>() / crust.len() as f64;

        // The relief under everything is self-affine, falling toward finer scales with the
        // Hurst exponent of cratered terrain (F54), so that what floods and renews it meets
        // rough ground at every scale and embays it; from the largest crater the body holds
        // down to two cells of the grid (F137).
        let raw = self.raw_relief();
        let raw_rms = (raw.iter().map(|v| v * v).sum::<f64>() / raw.len().max(1) as f64)
            .sqrt()
            .max(1e-12);
        let scale = self.relief_scale_for(raw_rms);

        let bedrock: Vec<f32> = (0..grid.cell_count())
            .into_par_iter()
            .map(|k| {
                let isostatic = (crust[k] as f64 - mean_crust) * self.isostatic_factor;
                (isostatic + scale * raw[k]) as f32
            })
            .collect();

        let units: Vec<u16> = match &self.continents {
            Some(cont) => {
                let mid = 0.5 * (cont.oceanic_thickness + cont.continental_thickness) as f32;
                crust
                    .iter()
                    .map(|&c| {
                        if c > mid {
                            UNIT_CONTINENTAL
                        } else {
                            UNIT_OCEANIC
                        }
                    })
                    .collect()
            }
            None => vec![UNIT_PRIMORDIAL; grid.cell_count()],
        };

        layers
            .crust_thickness
            .as_mut_slice()
            .copy_from_slice(&crust);
        layers.bedrock.as_mut_slice().copy_from_slice(&bedrock);
        layers.unit_id = Field::new(grid);
        layers.unit_id.as_mut_slice().copy_from_slice(&units);
        layers
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// F141: the relief's octaves regenerated on the loop's grid, all of them, are the
    /// relief init laid over the isostatic crust; those below a length carry little
    /// beyond it; and on a finer grid they run on below the loop's cells.
    #[test]
    fn the_fine_relief_is_the_relief_laid() {
        let (mut g, _) = crate::presets::preset("Mimas").unwrap().to_graph();
        let grid = crate::grid::CubeSphere::new(32);
        let r = g.f64("body.radius");
        let crust = PrimordialCrust::new(&g, grid);
        let bedrock = crust.init().bedrock;
        g.set("multiplier.relief", 0.0).unwrap();
        let flat = PrimordialCrust::new(&g, grid).init().bedrock;
        let all = crust.texture_field(&grid, f64::INFINITY, 0.0, r, 0, crust.loop_finest());
        let rms = |v: &[f32]| {
            (v.iter().map(|x| (*x as f64).powi(2)).sum::<f64>() / v.len() as f64).sqrt()
        };
        let diff: Vec<f32> = (0..all.len())
            .map(|k| bedrock.as_slice()[k] - flat.as_slice()[k] - all[k])
            .collect();
        assert!(rms(&all) > 100.0, "{}", rms(&all));
        assert!(
            rms(&diff) < 1e-3 * rms(&all),
            "{} of {}",
            rms(&diff),
            rms(&all)
        );
        // Below three tenths of the radius: little left after a low-pass over four times
        // that.
        let short = crust.fine_relief(&grid, 0.3 * r, r);
        let cell = grid.delta() * r;
        let mut field = crate::fields::Field::<f32>::new(grid);
        field.as_mut_slice().copy_from_slice(&short);
        let low = crate::stages::relaxation::low_pass(&field, cell, 1.2 * r);
        assert!(
            rms(&low) < 0.25 * rms(&short),
            "{} of {}",
            rms(&low),
            rms(&short)
        );
        // A grid four times as fine holds more octaves, with their share of the
        // amplitude, and the octaves below the loop's cells are that difference.
        let fine_grid = crate::grid::CubeSphere::new(128);
        let fine = crust.fine_relief(&fine_grid, 0.3 * r, r);
        assert!(
            rms(&fine) > rms(&short) && rms(&fine) < 1.5 * rms(&short),
            "{} vs {}",
            rms(&fine),
            rms(&short)
        );
        let below = crust.below_loop(&fine_grid, 0.0, r);
        assert!(
            rms(&below) > 0.0 && rms(&below) < rms(&short),
            "{}",
            rms(&below)
        );
        // Under a fabric spaced at the loop's finest octave, the octaves below it are
        // toned down to the share, the first of them in part.
        let toned = crust.below_loop(&fine_grid, 0.05 * r, r);
        assert!(
            rms(&toned) < 0.8 * rms(&below) && rms(&toned) > FINE_TEXTURE_SHARE * rms(&below),
            "{} of {}",
            rms(&toned),
            rms(&below)
        );
    }
}
