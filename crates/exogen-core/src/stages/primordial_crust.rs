// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! Primordial crust (PLAN.md 4.1): the initial layered state at formation.
//!
//! Crustal thickness is a low-degree spherical-harmonic field with an optional hemispheric
//! dichotomy term, or a two-level crust for mobile-lid bodies until plate tectonics
//! exists. Elevation follows from Airy isostasy, with a seeded self-affine relief from the
//! scale of the body's hemispheres down to the grid, as rough at its largest as a cratered
//! surface is at the scale of the largest crater the body holds (F137), standing in for
//! the relief its history before the record left.

use rayon::prelude::*;

use crate::constants::{
    DENSITY_CRUST_ICE, DENSITY_CRUST_ROCK, MAX_CRATER_RADIUS_FRACTION, ROUGHNESS_HURST,
    ROUGHNESS_RMS_SHARE,
};
use crate::fields::{Field, LayerSet};
use crate::grid::CubeSphere;
use crate::noise::GradientNoise;
use crate::params::ParamGraph;
use crate::rng::Rng;
use crate::spharm::Harmonics;
use crate::units::{UNIT_CONTINENTAL, UNIT_OCEANIC, UNIT_PRIMORDIAL};

/// Bumped when the behaviour changes, so caches are not reused across versions.
pub const VERSION: u64 = 5;

/// Two-level crust for mobile-lid bodies: thin oceanic crust and thick continental crust
/// covering a fraction of the area, which gives the bimodal hypsometry that lets an ocean
/// fill the basins and leave continents standing. Plate tectonics replaces it in Phase 6.
struct Continents {
    oceanic_thickness: f64,
    continental_thickness: f64,
    area_fraction: f64,
}

pub struct PrimordialCrust {
    grid: CubeSphere,
    noise: GradientNoise,
    harmonics: Harmonics,
    dichotomy_axis: [f64; 3],
    dichotomy_amplitude: f64,
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
    /// The relief's scale on its normalized field, m, once the loop's grid has measured
    /// it (F141).
    relief_scale: std::sync::OnceLock<f64>,
}

impl PrimordialCrust {
    pub fn new(graph: &ParamGraph, grid: CubeSphere) -> Self {
        let seed = graph.f64("seed.terrain") as u64;
        let rng = Rng::new(seed);
        let icy = matches!(
            graph.str("body.composition_class"),
            "rock_ice" | "ice_dominated"
        );
        let (rho_mantle, rho_crust) = if icy {
            (1000.0, DENSITY_CRUST_ICE)
        } else {
            (3300.0, DENSITY_CRUST_ROCK)
        };
        // F137: the relief of a cratered surface at the scale of the largest crater the
        // body holds, a share of its depth, as the history before the record left it;
        // never more than a quarter of the relief the crust's strength supports.
        let scaling = crate::crater::Scaling {
            simple_complex: graph.f64("cratering.simple_complex_transition"),
            peak_ring: graph.f64("cratering.peak_ring_onset"),
            multi_ring: graph.f64("cratering.multiring_onset"),
            gravity: graph.f64("body.gravity"),
            icy: crate::stages::materials::is_ice(graph.str("materials.bedrock")),
        };
        let radius = graph.f64("body.radius");
        let largest = crate::crater::hold_diameter(radius, scaling.icy)
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
        Self {
            grid,
            noise: GradientNoise::new(seed ^ 0x9E37_79B9),
            harmonics: Harmonics::random(&rng, 2, 8, 2.0),
            dichotomy_axis: rng.unit_direction(1, 0),
            dichotomy_amplitude: graph.f64("tectonics.dichotomy_amplitude"),
            crust_mean: graph.f64("tectonics.initial_crust_thickness"),
            crust_sigma: 0.12 * graph.f64("tectonics.initial_crust_thickness"),
            continents,
            relief,
            relief_frequency,
            relief_cap,
            isostatic_factor: (rho_mantle - rho_crust) / rho_mantle,
            relief_scale: std::sync::OnceLock::new(),
        }
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
        *self.relief_scale.get_or_init(|| {
            let raw = self.raw_relief();
            let rms = (raw.iter().map(|v| v * v).sum::<f64>() / raw.len().max(1) as f64).sqrt();
            self.relief_scale_for(rms)
        })
    }

    /// The first of the relief's octaves shorter than `max_wavelength` m: octave o has the
    /// wavelength R / (f 2^o), f the relief's frequency in lattice cells per radian.
    pub fn first_octave(&self, max_wavelength: f64, body_radius: f64) -> u32 {
        let x = body_radius / (self.relief_frequency * max_wavelength.max(1.0));
        if x < 1.0 {
            0
        } else {
            x.log2().floor() as u32 + 1
        }
    }

    /// The relief's octaves from `first` on, on any grid, down to two of its cells, each at
    /// the amplitude the loop's grid gave it (F141): the texture the base relief has under
    /// a terrain that renewed the ground, below the scale of the terrain's own fabric,
    /// which the renewal's smoothing took with the craters.
    pub fn fine_relief_from(&self, grid: &CubeSphere, first: u32) -> Vec<f32> {
        let gain = 2f64.powf(-ROUGHNESS_HURST);
        let norm: f64 = (0..self.relief_octaves(&self.grid))
            .map(|o| gain.powi(o as i32))
            .sum();
        let last = self.relief_octaves(grid) as i64 - 1;
        if last < first as i64 {
            return vec![0.0; grid.cell_count()];
        }
        let scale = self.relief_scale() / norm.max(1e-12);
        let count = (last - first as i64 + 1) as u32;
        (0..grid.cell_count())
            .into_par_iter()
            .map(|k| {
                (scale
                    * self.noise.octaves(
                        grid.center(grid.cell(k)),
                        self.relief_frequency,
                        first,
                        count,
                        2.0,
                        gain,
                    )) as f32
            })
            .collect()
    }

    /// The relief's octaves shorter than `max_wavelength` m on any grid (F141).
    pub fn fine_relief(
        &self,
        grid: &CubeSphere,
        max_wavelength: f64,
        body_radius: f64,
    ) -> Vec<f32> {
        self.fine_relief_from(grid, self.first_octave(max_wavelength, body_radius))
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

        let axis = self.dichotomy_axis;
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
                .map(|k| {
                    let d = grid.center(grid.cell(k));
                    let dot = d[0] * axis[0] + d[1] * axis[1] + d[2] * axis[2];
                    let dichotomy = 0.5 * self.dichotomy_amplitude * (3.0 * dot).tanh();
                    (self.crust_mean + dichotomy + self.crust_sigma * raw[k] / rms) as f32
                })
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
        let scale = *self
            .relief_scale
            .get_or_init(|| self.relief_scale_for(raw_rms));

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
    /// relief init laid over the isostatic crust; those shorter than a length carry
    /// little beyond it; and on a finer grid they run on below the loop's cells.
    #[test]
    fn the_fine_relief_is_the_relief_laid() {
        let (mut g, _) = crate::presets::preset("Mimas").unwrap().to_graph();
        let grid = crate::grid::CubeSphere::new(32);
        let r = g.f64("body.radius");
        let crust = PrimordialCrust::new(&g, grid);
        let bedrock = crust.init().bedrock;
        g.set("multiplier.relief", 0.0).unwrap();
        let flat = PrimordialCrust::new(&g, grid).init().bedrock;
        let all = crust.fine_relief(&grid, f64::INFINITY, r);
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
        // Shorter than a tenth of the radius: little left after a low-pass over twice that.
        let short = crust.fine_relief(&grid, 0.1 * r, r);
        let cell = grid.delta() * r;
        let mut field = crate::fields::Field::<f32>::new(grid);
        field.as_mut_slice().copy_from_slice(&short);
        let low = crate::stages::relaxation::low_pass(&field, cell, 0.4 * r);
        assert!(
            rms(&low) < 0.25 * rms(&short),
            "{} of {}",
            rms(&low),
            rms(&short)
        );
        // A grid twice as fine holds one more octave, with its share of the amplitude.
        let fine = crust.fine_relief(&crate::grid::CubeSphere::new(64), 0.1 * r, r);
        assert!(
            rms(&fine) > rms(&short) && rms(&fine) < 1.5 * rms(&short),
            "{} vs {}",
            rms(&fine),
            rms(&short)
        );
    }
}
