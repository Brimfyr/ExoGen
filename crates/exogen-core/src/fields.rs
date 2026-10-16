// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! Per-cell field storage on a cubesphere, the layered surface representation
//! (PLAN.md 3.2), resampling between resolutions and content hashing (ADR 0003).

use rayon::prelude::*;

use crate::grid::projection::PaddedField;
use crate::grid::{Cell, CubeSphere};
use crate::params::fnv1a64;

/// A value per cell of a cubesphere.
#[derive(Clone, Debug)]
pub struct Field<T> {
    grid: CubeSphere,
    data: Vec<T>,
}

impl<T: Copy + Default + Send + Sync> Field<T> {
    pub fn new(grid: CubeSphere) -> Self {
        Self {
            grid,
            data: vec![T::default(); grid.cell_count()],
        }
    }

    pub fn filled(grid: CubeSphere, v: T) -> Self {
        Self {
            grid,
            data: vec![v; grid.cell_count()],
        }
    }

    /// Build a field in parallel from a function of the cell and its unit direction.
    pub fn from_fn<F>(grid: CubeSphere, f: F) -> Self
    where
        F: Fn(Cell, [f64; 3]) -> T + Sync,
    {
        let data: Vec<T> = (0..grid.cell_count())
            .into_par_iter()
            .map(|k| {
                let c = grid.cell(k);
                f(c, grid.center(c))
            })
            .collect();
        Self { grid, data }
    }

    pub fn grid(&self) -> CubeSphere {
        self.grid
    }

    pub fn len(&self) -> usize {
        self.data.len()
    }

    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    #[inline]
    pub fn get(&self, c: Cell) -> T {
        self.data[self.grid.index(c)]
    }

    #[inline]
    pub fn set(&mut self, c: Cell, v: T) {
        let i = self.grid.index(c);
        self.data[i] = v;
    }

    pub fn as_slice(&self) -> &[T] {
        &self.data
    }

    pub fn as_mut_slice(&mut self) -> &mut [T] {
        &mut self.data
    }

    pub fn into_vec(self) -> Vec<T> {
        self.data
    }

    /// Apply a function to every value in parallel.
    pub fn map_in_place<F>(&mut self, f: F)
    where
        F: Fn(T) -> T + Sync,
    {
        self.data.par_iter_mut().for_each(|v| *v = f(*v));
    }
}

/// Deterministic parallel sum (ADR 0004): fixed-size chunks are summed in parallel, each
/// sequentially, and the chunk sums are added in order. The result does not depend on the
/// number of threads.
pub fn det_sum<F>(n: usize, f: F) -> f64
where
    F: Fn(usize) -> f64 + Sync,
{
    const CHUNK: usize = 4096;
    let chunks: Vec<f64> = (0..n.div_ceil(CHUNK))
        .into_par_iter()
        .map(|c| {
            let start = c * CHUNK;
            let end = (start + CHUNK).min(n);
            let mut s = 0.0;
            for k in start..end {
                s += f(k);
            }
            s
        })
        .collect();
    chunks.iter().sum()
}

impl Field<f32> {
    pub fn min_max(&self) -> (f32, f32) {
        self.data
            .par_iter()
            .fold(
                || (f32::INFINITY, f32::NEG_INFINITY),
                |(lo, hi), &v| (lo.min(v), hi.max(v)),
            )
            .reduce(
                || (f32::INFINITY, f32::NEG_INFINITY),
                |a, b| (a.0.min(b.0), a.1.max(b.1)),
            )
    }

    /// Area-weighted mean over the sphere.
    pub fn area_mean(&self) -> f64 {
        let g = self.grid;
        let n = g.cell_count();
        let sum = det_sum(n, |k| self.data[k] as f64 * g.cell_area(g.cell(k)));
        let area = det_sum(n, |k| g.cell_area(g.cell(k)));
        sum / area
    }

    /// Bilinear resample to another resolution of the same sphere.
    pub fn resample(&self, target: CubeSphere) -> Field<f32> {
        if target == self.grid {
            return self.clone();
        }
        let padded = PaddedField::new(&self.grid, &self.data);
        let src = self.grid;
        Field::from_fn(target, |_, d| {
            crate::grid::projection::sample_bilinear(&src, &padded, d)
        })
    }

    /// Resample by the nearest source cell: for a field of event times, which must not
    /// blend across a cell's edge into times no event had.
    pub fn resample_nearest(&self, target: CubeSphere) -> Field<f32> {
        if target == self.grid {
            return self.clone();
        }
        let src = self.grid;
        Field::from_fn(target, |_, d| self.data[src.index(src.locate(d).0)])
    }

    /// Stable content hash over the exact bit patterns.
    pub fn hash(&self) -> u64 {
        let mut bytes = Vec::with_capacity(self.data.len() * 4 + 8);
        bytes.extend_from_slice(&(self.grid.face_size() as u64).to_le_bytes());
        for v in &self.data {
            bytes.extend_from_slice(&v.to_bits().to_le_bytes());
        }
        fnv1a64(&bytes)
    }

    /// Area-weighted hypsometric fill (docs/parameters.md F26): the level `z` at which the
    /// volume below `z` equals `volume` on a sphere of radius `radius`. Heights and the
    /// result are in metres relative to the same datum.
    pub fn fill_level(&self, radius: f64, volume: f64) -> f64 {
        if volume <= 0.0 {
            let (lo, _) = self.min_max();
            return lo as f64 - 1.0;
        }
        let g = self.grid;
        let areas: Vec<f64> = (0..g.cell_count())
            .map(|k| g.cell_area(g.cell(k)) * radius * radius)
            .collect();
        let volume_below = |z: f64| -> f64 {
            det_sum(self.data.len(), |k| {
                (z - self.data[k] as f64).max(0.0) * areas[k]
            })
        };
        let (lo, hi) = self.min_max();
        let (mut lo, mut hi) = (
            lo as f64,
            hi as f64 + volume / areas.iter().sum::<f64>() + 1.0,
        );
        for _ in 0..80 {
            let mid = 0.5 * (lo + hi);
            if volume_below(mid) < volume {
                lo = mid
            } else {
                hi = mid
            }
        }
        0.5 * (lo + hi)
    }
}

impl Field<u16> {
    /// Nearest-cell resample for label fields.
    pub fn resample_nearest(&self, target: CubeSphere) -> Field<u16> {
        if target == self.grid {
            return self.clone();
        }
        let src = self.grid;
        Field::from_fn(target, |_, d| self.data[src.index(src.locate(d).0)])
    }

    pub fn hash(&self) -> u64 {
        let mut bytes = Vec::with_capacity(self.data.len() * 2);
        for v in &self.data {
            bytes.extend_from_slice(&v.to_le_bytes());
        }
        fnv1a64(&bytes)
    }
}

/// The layered surface state (PLAN.md 3.2). Heights are metres relative to the mean
/// radius. Ages are Gyr since formation of the body (time of formation for crust, time of
/// last resurfacing for the surface).
#[derive(Clone, Debug)]
pub struct LayerSet {
    pub grid: CubeSphere,
    pub bedrock: Field<f32>,
    pub regolith: Field<f32>,
    pub dust: Field<f32>,
    pub ice: Field<f32>,
    pub liquid: Field<f32>,
    pub crust_thickness: Field<f32>,
    pub crust_age: Field<f32>,
    pub surface_age: Field<f32>,
    pub unit_id: Field<u16>,
    /// Fraction of the surface material that came from the deep bedrock, 0 to 1.
    pub deep: Field<f32>,
    /// Impact melt cover, 0 to 1.
    pub melt: Field<f32>,
    /// Coverage of a deposit laid over the substrate, 0 to 1: lava by its thickness,
    /// plume deposits by their density. What is not covered shows the substrate.
    pub cover: Field<f32>,
    /// Share of the substrate that is clean ice renewed from below, 0 to 1: the lanes of
    /// grooved terrain and the polar terrain of an ice shell, which stay bright under
    /// whatever craters cut them.
    pub clean: Field<f32>,
    /// The composition of a lava surface as the drift of its eruption's lavas from the
    /// material's mean, -1 to 1 (F55): a flow unit's tone, which the ejecta of the same
    /// basalt keeps; 0 where no eruption set it.
    pub composition: Field<f32>,
    /// When the lava at the surface erupted, Gyr since formation (F95): the age of the
    /// contacts it makes with its neighbours, which gardening has mixed since, kept under
    /// the craters and deposits later laid on it; 0 where no eruption set it.
    pub lava_age: Field<f32>,
    /// The share of the ground on slopes steep enough to shed its regolith, by how far past
    /// the start they are, at the scale the slopes are steep (F98): the walls of the
    /// craters from their own profiles, finer than the grid; a renewed surface carries
    /// none of the old.
    pub steep: Field<f32>,
    /// The age whose exposure gives the maturity of the fresh ejecta that craters have laid
    /// thin over older ground beyond their rims, Gyr since formation (F100): the craters
    /// under such a deposit survive, so the surface age, which the crater population
    /// reads, keeps the ground's own, and the maturity reads the younger of the two.
    pub ejecta_age: Field<f32>,
    /// The share of the surface stained by the red-brown non-ice material of a crack's
    /// margins (F78; the triple bands of Geissler et al. 1998), 0 to 1, as it stands at
    /// the present: laid by each crack as it forms and faded since; a renewed surface
    /// carries none.
    pub stain: Field<f32>,
}

impl LayerSet {
    pub fn new(grid: CubeSphere) -> Self {
        Self {
            grid,
            bedrock: Field::new(grid),
            regolith: Field::new(grid),
            dust: Field::new(grid),
            ice: Field::new(grid),
            liquid: Field::new(grid),
            crust_thickness: Field::new(grid),
            crust_age: Field::new(grid),
            surface_age: Field::new(grid),
            unit_id: Field::new(grid),
            deep: Field::new(grid),
            melt: Field::new(grid),
            cover: Field::new(grid),
            clean: Field::new(grid),
            composition: Field::new(grid),
            lava_age: Field::new(grid),
            steep: Field::new(grid),
            ejecta_age: Field::new(grid),
            stain: Field::new(grid),
        }
    }

    /// Solid surface height: bedrock plus loose and frozen layers.
    pub fn height(&self) -> Field<f32> {
        let mut h = self.bedrock.clone();
        let (r, d, i) = (
            self.regolith.as_slice(),
            self.dust.as_slice(),
            self.ice.as_slice(),
        );
        h.as_mut_slice()
            .par_iter_mut()
            .enumerate()
            .for_each(|(k, v)| *v += r[k] + d[k] + i[k]);
        h
    }

    pub fn resample(&self, target: CubeSphere) -> LayerSet {
        LayerSet {
            grid: target,
            bedrock: self.bedrock.resample(target),
            regolith: self.regolith.resample(target),
            dust: self.dust.resample(target),
            ice: self.ice.resample(target),
            liquid: self.liquid.resample(target),
            crust_thickness: self.crust_thickness.resample(target),
            crust_age: self.crust_age.resample(target),
            // A surface age is the time of the event that made the surface: a cell takes
            // its nearest coarse cell's, never a blend of two events' times, which the
            // finish's passes would read as ramps of age along every coarse edge.
            surface_age: self.surface_age.resample_nearest(target),
            unit_id: self.unit_id.resample_nearest(target),
            deep: self.deep.resample(target),
            melt: self.melt.resample(target),
            cover: self.cover.resample(target),
            clean: self.clean.resample(target),
            composition: self.composition.resample_nearest(target),
            lava_age: self.lava_age.resample_nearest(target),
            steep: self.steep.resample(target),
            ejecta_age: self.ejecta_age.resample_nearest(target),
            stain: self.stain.resample(target),
        }
    }

    /// The age whose exposure sets each cell's maturity (F100): its surface's, or the thin
    /// ejecta's laid over it since, whichever is younger.
    pub fn exposure_age(&self, k: usize) -> f32 {
        self.surface_age.as_slice()[k].max(self.ejecta_age.as_slice()[k])
    }

    pub fn hash(&self) -> u64 {
        let parts = [
            self.bedrock.hash(),
            self.regolith.hash(),
            self.dust.hash(),
            self.ice.hash(),
            self.liquid.hash(),
            self.crust_thickness.hash(),
            self.crust_age.hash(),
            self.surface_age.hash(),
            self.unit_id.hash(),
            self.deep.hash(),
            self.melt.hash(),
            self.cover.hash(),
            self.clean.hash(),
            self.composition.hash(),
            self.lava_age.hash(),
            self.steep.hash(),
            self.ejecta_age.hash(),
            self.stain.hash(),
        ];
        crate::rng::hash_words(&parts)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::grid::projection::test_field;

    #[test]
    fn fill_level_matches_the_requested_volume_and_a_flat_sphere_is_exact() {
        let g = CubeSphere::new(16);
        let flat = Field::filled(g, 0.0f32);
        let r = 1.0e6;
        let area = 4.0 * std::f64::consts::PI * r * r;
        let z = flat.fill_level(r, 100.0 * area);
        assert!((z - 100.0).abs() < 1e-6, "{z}");

        let h = Field::from_fn(g, |_, d| (1000.0 * test_field(d)) as f32);
        let target = 300.0 * area;
        let z = h.fill_level(r, target);
        let got: f64 = (0..g.cell_count())
            .map(|k| (z - h.as_slice()[k] as f64).max(0.0) * g.cell_area(g.cell(k)) * r * r)
            .sum();
        assert!((got - target).abs() / target < 1e-9);
    }

    #[test]
    fn resampling_preserves_a_smooth_field() {
        let coarse = CubeSphere::new(32);
        let fine = CubeSphere::new(64);
        let f = Field::from_fn(coarse, |_, d| test_field(d) as f32);
        let up = f.resample(fine);
        let mut worst = 0.0f32;
        for c in fine.cells() {
            worst = worst.max((up.get(c) - test_field(fine.center(c)) as f32).abs());
        }
        assert!(worst < 0.02, "{worst}");
        assert_ne!(f.hash(), up.hash());
        assert_eq!(f.hash(), f.clone().hash());
    }

    #[test]
    fn area_mean_of_a_constant_is_the_constant() {
        let g = CubeSphere::new(8);
        assert!((Field::filled(g, 3.5f32).area_mean() - 3.5).abs() < 1e-9);
    }
}
