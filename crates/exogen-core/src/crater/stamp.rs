// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! Stamping craters onto cubesphere layers.
//!
//! A footprint is gathered by a breadth-first walk over edge neighbours from the centre
//! cell while the great-circle distance stays within the ejecta radius. The crater floor
//! replaces the pre-impact topography relative to the mean height inside the rim; the rim
//! and ejecta are added to whatever is there, with a short blend across the rim. Cells under
//! the continuous ejecta have their surface age reset to the impact time.

use super::morphology::ProfileTable;
use crate::constants::CONTINUOUS_EJECTA_MULTIPLE;
use crate::fields::LayerSet;
use crate::geom::{angle_between, azimuth_from, local_basis};
use crate::grid::CubeSphere;

/// Reusable scratch space for footprints.
pub struct Stamper {
    visited: Vec<u32>,
    generation: u32,
    queue: Vec<usize>,
    /// (cell index, distance in metres, azimuth)
    pub footprint: Vec<(usize, f64, f64)>,
}

impl Stamper {
    pub fn new(grid: &CubeSphere) -> Self {
        Self {
            visited: vec![0; grid.cell_count()],
            generation: 0,
            queue: Vec::new(),
            footprint: Vec::new(),
        }
    }

    /// Collect the cells within `radius_m` of `center` on a body of radius `body_radius`.
    pub fn gather(&mut self, grid: &CubeSphere, body_radius: f64, center: [f64; 3], radius_m: f64) {
        self.generation = self.generation.wrapping_add(1);
        if self.generation == 0 {
            self.visited.iter_mut().for_each(|v| *v = 0);
            self.generation = 1;
        }
        self.footprint.clear();
        self.queue.clear();
        let theta_max = (radius_m / body_radius).min(std::f64::consts::PI);
        let (east, north) = local_basis(center);
        let start = grid.index(grid.locate(center).0);
        self.visited[start] = self.generation;
        self.queue.push(start);
        let mut head = 0;
        while head < self.queue.len() {
            let idx = self.queue[head];
            head += 1;
            let cell = grid.cell(idx);
            let d = grid.center(cell);
            let theta = angle_between(center, d);
            // Accept cells whose centre is within the radius plus half a cell so the edge
            // is not under-covered; the profile returns zero beyond the ejecta anyway.
            if theta > theta_max + 0.5 * grid.delta() {
                continue;
            }
            let phi = azimuth_from(center, east, north, d);
            self.footprint.push((idx, theta * body_radius, phi));
            for (di, dj) in [(1, 0), (-1, 0), (0, 1), (0, -1)] {
                let n = grid.index(grid.neighbor(cell, di, dj));
                if self.visited[n] != self.generation {
                    self.visited[n] = self.generation;
                    self.queue.push(n);
                }
            }
        }
    }
}

/// Effective radius for an elliptical crater: scaled so the long axis lies along `azimuth`
/// and the area is preserved.
#[inline]
pub fn elliptical_radius(r: f64, phi: f64, azimuth: f64, ellipticity: f64) -> f64 {
    if ellipticity <= 1.0 {
        return r;
    }
    let (s, c) = (phi - azimuth).sin_cos();
    r * ((c * c) / ellipticity + (s * s) * ellipticity).sqrt()
}

pub struct StampParams {
    pub time_gyr: f32,
    pub unit_floor: u16,
    pub unit_ejecta: u16,
    pub azimuth: f64,
    pub ellipticity: f64,
    /// Multiply the profile height before applying (1 for a full stamp).
    pub amplitude: f64,
    /// When false the floor does not replace the old topography; the profile is purely
    /// additive (used to add resolved detail on top of a coarse stamp).
    pub replace_floor: bool,
    pub set_age: bool,
}

/// Stamp one crater. Returns the number of cells touched.
pub fn stamp(
    layers: &mut LayerSet,
    grid: &CubeSphere,
    body_radius: f64,
    center: [f64; 3],
    profile: &ProfileTable,
    p: &StampParams,
    stamper: &mut Stamper,
) -> usize {
    let dim = &profile.dim;
    stamper.gather(
        grid,
        body_radius,
        center,
        dim.ejecta_radius * p.ellipticity.sqrt(),
    );
    if stamper.footprint.is_empty() {
        return 0;
    }
    let bedrock = layers.bedrock.as_mut_slice();
    // Reference height: area-weighted mean of the pre-impact surface inside the rim.
    let mut sum = 0.0;
    let mut wsum = 0.0;
    for &(idx, r, phi) in &stamper.footprint {
        let re = elliptical_radius(r, phi, p.azimuth, p.ellipticity);
        if re <= dim.radius {
            let a = grid.cell_area(grid.cell(idx));
            sum += bedrock[idx] as f64 * a;
            wsum += a;
        }
    }
    let h_ref = if wsum > 0.0 {
        sum / wsum
    } else {
        bedrock[stamper.footprint[0].0] as f64
    };
    let continuous = CONTINUOUS_EJECTA_MULTIPLE * dim.radius;
    let ages = layers.surface_age.as_mut_slice();
    let units = layers.unit_id.as_mut_slice();
    for &(idx, r, phi) in &stamper.footprint {
        let re = elliptical_radius(r, phi, p.azimuth, p.ellipticity);
        let v = profile.height(re) * p.amplitude;
        let old = bedrock[idx] as f64;
        let w = if p.replace_floor {
            // 1 inside 0.95 R, 0 outside 1.05 R.
            let t = (1.05 - re / dim.radius) / 0.10;
            t.clamp(0.0, 1.0)
        } else {
            0.0
        };
        bedrock[idx] = (v + w * h_ref + (1.0 - w) * old) as f32;
        if p.set_age && re <= continuous {
            ages[idx] = p.time_gyr;
            units[idx] = if re <= dim.radius {
                p.unit_floor
            } else {
                p.unit_ejecta
            };
        }
    }
    stamper.footprint.len()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::crater::morphology::{Scaling, dimensions};
    use crate::fields::Field;

    #[test]
    fn a_stamped_fresh_crater_has_the_right_depth_and_resets_ages_under_its_ejecta() {
        let grid = CubeSphere::new(256);
        let body_radius = 1.7374e6;
        let mut layers = LayerSet::new(grid);
        layers.bedrock = Field::filled(grid, 100.0);
        layers.surface_age = Field::filled(grid, 0.0);
        let s = Scaling {
            simple_complex: 16_500.0,
            peak_ring: 140_000.0,
            multi_ring: 330_000.0,
            gravity: 1.62,
        };
        let dim = dimensions(200_000.0, &s);
        let table = ProfileTable::new(dim, 0.0);
        let mut stamper = Stamper::new(&grid);
        let center = crate::geom::normalize([0.3, 0.4, 0.86]);
        let params = StampParams {
            time_gyr: 3.0,
            unit_floor: 3,
            unit_ejecta: 4,
            azimuth: 0.0,
            ellipticity: 1.0,
            amplitude: 1.0,
            replace_floor: true,
            set_age: true,
        };
        let touched = stamp(
            &mut layers,
            &grid,
            body_radius,
            center,
            &table,
            &params,
            &mut stamper,
        );
        assert!(touched > 100);
        let (lo, hi) = layers.bedrock.min_max();
        let expected_depth = dim.depth + dim.rim_height;
        assert!(
            ((hi - lo) as f64 - expected_depth).abs() < 0.1 * expected_depth,
            "measured {} expected {expected_depth}",
            hi - lo
        );
        // The floor sits at the reference height minus the depth.
        assert!(
            (lo as f64 - (100.0 - dim.depth)).abs() < 0.05 * dim.depth,
            "floor {lo}"
        );
        // Ages reset under the continuous ejecta but not far away.
        let centre_cell = grid.locate(center).0;
        assert_eq!(layers.surface_age.get(centre_cell), 3.0);
        assert_eq!(layers.unit_id.get(centre_cell), 3);
        let far = grid.locate(crate::geom::normalize([-0.9, 0.1, -0.3])).0;
        assert_eq!(layers.surface_age.get(far), 0.0);
        // Cells reset: area of a disc of 2 R on the sphere.
        let reset: f64 = grid
            .cells()
            .filter(|c| layers.surface_age.get(*c) == 3.0)
            .map(|c| grid.cell_area(c))
            .sum();
        let cap = 2.0 * std::f64::consts::PI * (1.0 - (2.0 * dim.radius / body_radius).cos());
        assert!((reset - cap).abs() < 0.08 * cap, "reset {reset} cap {cap}");
    }

    #[test]
    fn footprint_handles_face_edges_and_poles() {
        let grid = CubeSphere::new(32);
        let mut stamper = Stamper::new(&grid);
        for center in [
            [0.0, 0.0, 1.0],
            crate::geom::normalize([1.0, 1.0, 0.0]),
            crate::geom::normalize([1.0, 1.0, 1.0]),
        ] {
            stamper.gather(&grid, 1.0e6, center, 0.3e6);
            let n = stamper.footprint.len() as f64;
            let expected = 2.0 * std::f64::consts::PI * (1.0 - 0.3f64.cos())
                / (4.0 * std::f64::consts::PI)
                * grid.cell_count() as f64;
            assert!(
                (n - expected).abs() < 0.25 * expected + 20.0,
                "{n} vs {expected}"
            );
        }
    }
}
