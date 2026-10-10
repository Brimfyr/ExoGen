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
use crate::constants::{
    BURIAL_FACTOR_HIGH, BURIAL_FACTOR_LOW, CONTINUOUS_EJECTA_MULTIPLE, OBLIQUE_ASYMMETRY_DEG,
    OBLIQUE_BUTTERFLY_DEG, OBLIQUE_BUTTERFLY_MAX, OBLIQUE_DOWNRANGE_MAX, OBLIQUE_FORBIDDEN_DEG,
    OBLIQUE_FORBIDDEN_DEPTH, OBLIQUE_FORBIDDEN_HALF_WIDTH_DEG,
};
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
    /// Downrange direction, radians in the local east-north frame.
    pub azimuth: f64,
    /// Impact angle from horizontal, radians; π/2 is vertical.
    pub angle: f64,
    pub ellipticity: f64,
    /// Multiply the profile height before applying (1 for a full stamp).
    pub amplitude: f64,
    /// When false the floor does not replace the old topography; the profile is purely
    /// additive (used to add resolved detail on top of a coarse stamp).
    pub replace_floor: bool,
    pub set_age: bool,
    /// Bury the relief already there under the ejecta according to its thickness.
    pub bury: bool,
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
    // Reference height: area-weighted mean of the pre-impact surface inside the rim; and
    // the azimuthal mean by radius under the ejecta, against which burial damps relief.
    const BINS: usize = 16;
    let span = (dim.ejecta_radius / dim.radius - 1.0).max(1e-6);
    let bin_of = |re: f64| -> usize {
        (((re / dim.radius - 1.0) / span) * BINS as f64)
            .floor()
            .clamp(0.0, (BINS - 1) as f64) as usize
    };
    let mut sum = 0.0;
    let mut wsum = 0.0;
    // Ring means are plain cell averages: cell areas differ by at most 40% on the grid,
    // and the mean only serves as a local reference.
    let mut ring_sum = [0.0f64; BINS];
    let mut ring_area = [0.0f64; BINS];
    for &(idx, r, phi) in &stamper.footprint {
        let re = elliptical_radius(r, phi, p.azimuth, p.ellipticity);
        if re <= dim.radius {
            let a = grid.cell_area(grid.cell(idx));
            sum += bedrock[idx] as f64 * a;
            wsum += a;
        } else if p.bury && re < dim.ejecta_radius {
            let b = bin_of(re);
            ring_sum[b] += bedrock[idx] as f64;
            ring_area[b] += 1.0;
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
        let v = profile.height(re)
            * p.amplitude
            * ejecta_factor(re, dim.radius, phi, p.azimuth, p.angle);
        let mut old = bedrock[idx] as f64;
        if p.bury && re > dim.radius && v > 0.0 {
            let b = bin_of(re);
            let mean = if ring_area[b] > 0.0 {
                ring_sum[b] / ring_area[b]
            } else {
                h_ref
            };
            old = mean + (old - mean) * burial_factor(old - mean, v);
        }
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

/// Azimuthal modulation of ejecta thickness for an oblique impact, with `phi` the azimuth
/// from the centre and `azimuth` the downrange direction: ejecta thickens downrange and
/// thins uprange below 45°, an uprange zone of avoidance opens below 30°, and side lobes
/// grow below 15° (Gault and Wedekind 1978; Herrick and Forsberg-Taylor 2003).
#[inline]
pub fn ejecta_modulation(phi: f64, azimuth: f64, angle: f64) -> f64 {
    let dphi = phi - azimuth;
    let (s, c) = dphi.sin_cos();
    let deg = angle.to_degrees();
    let down = (1.0 - deg / OBLIQUE_ASYMMETRY_DEG).clamp(0.0, 1.0);
    let mut m = 1.0 + OBLIQUE_DOWNRANGE_MAX * down * c;
    if deg < OBLIQUE_FORBIDDEN_DEG {
        let strength = 1.0 - deg / OBLIQUE_FORBIDDEN_DEG;
        let half_width = OBLIQUE_FORBIDDEN_HALF_WIDTH_DEG.to_radians() * strength;
        // Angle from the uprange direction.
        let from_up = dphi.rem_euclid(2.0 * std::f64::consts::PI) - std::f64::consts::PI;
        m *= 1.0 - OBLIQUE_FORBIDDEN_DEPTH * strength * (-(from_up / half_width).powi(2)).exp();
    }
    if deg < OBLIQUE_BUTTERFLY_DEG {
        let strength = 1.0 - deg / OBLIQUE_BUTTERFLY_DEG;
        m *= 1.0 + OBLIQUE_BUTTERFLY_MAX * strength * s * s;
    }
    m.max(0.0)
}

/// The modulation faded in across the rim so the cavity side stays symmetric.
#[inline]
fn ejecta_factor(re: f64, radius: f64, phi: f64, azimuth: f64, angle: f64) -> f64 {
    let ramp = ((re / radius - 1.0) / 0.3).clamp(0.0, 1.0);
    if ramp <= 0.0 {
        return 1.0;
    }
    1.0 + (ejecta_modulation(phi, azimuth, angle) - 1.0) * ramp
}

/// Parameters of a loop-crater restore at the finishing resolution.
pub struct RestoreParams {
    pub time_gyr: f32,
    pub unit_floor: u16,
    pub unit_ejecta: u16,
    pub azimuth: f64,
    pub angle: f64,
    pub ellipticity: f64,
    /// Weight of the resolved detail: the fraction of the crater that survived the loop.
    pub detail_weight: f64,
}

/// Relief under ejecta of thickness `thickness` is damped by this factor: depressions fill
/// and highs are shaken down, with the thickness measured against the relief itself so
/// thin ejecta on deep relief changes little and thick ejecta on shallow relief erases it.
#[inline]
pub fn burial_factor(relief: f64, thickness: f64) -> f64 {
    if thickness <= 0.0 {
        return 1.0;
    }
    let k = if relief < 0.0 {
        BURIAL_FACTOR_LOW
    } else {
        BURIAL_FACTOR_HIGH
    };
    (-k * thickness / relief.abs().max(1.0)).exp()
}

/// Re-stamp a loop crater at the finishing resolution, in its place in the age order.
///
/// `base` is the resampled coarse height before any finishing-pass stamp; it already holds
/// the crater's coarse form. Inside the rim the height returns to the base plus the
/// resolved detail, which removes whatever the pass stamped there earlier: those craters
/// are older and the impact destroyed them. Under the ejecta the relief added earlier in
/// the pass is buried according to the ejecta thickness. Ages and units within the
/// continuous ejecta are reset to the impact time.
#[allow(clippy::too_many_arguments)]
pub fn restore(
    layers: &mut LayerSet,
    grid: &CubeSphere,
    body_radius: f64,
    center: [f64; 3],
    base: &[f32],
    fresh: &ProfileTable,
    detail: &ProfileTable,
    p: &RestoreParams,
    stamper: &mut Stamper,
) -> usize {
    let dim = &fresh.dim;
    stamper.gather(
        grid,
        body_radius,
        center,
        dim.ejecta_radius * p.ellipticity.sqrt(),
    );
    if stamper.footprint.is_empty() {
        return 0;
    }
    let continuous = CONTINUOUS_EJECTA_MULTIPLE * dim.radius;
    let bedrock = layers.bedrock.as_mut_slice();
    let ages = layers.surface_age.as_mut_slice();
    let units = layers.unit_id.as_mut_slice();
    for &(idx, r, phi) in &stamper.footprint {
        let re = elliptical_radius(r, phi, p.azimuth, p.ellipticity);
        let fine = bedrock[idx] as f64 - base[idx] as f64;
        // 1 inside 0.95 R, 0 outside 1.05 R: the cavity removes the fine relief.
        let cut = ((1.05 - re / dim.radius) / 0.10).clamp(0.0, 1.0);
        let thickness = if re > dim.radius {
            (fresh.height(re) * ejecta_factor(re, dim.radius, phi, p.azimuth, p.angle)).max(0.0)
        } else {
            0.0
        };
        let kept = fine * (1.0 - cut) * burial_factor(fine, thickness);
        let d = detail.height(re)
            * p.detail_weight
            * ejecta_factor(re, dim.radius, phi, p.azimuth, p.angle);
        bedrock[idx] = (base[idx] as f64 + kept + d) as f32;
        if re <= continuous {
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
            angle: std::f64::consts::FRAC_PI_2,
            ellipticity: 1.0,
            amplitude: 1.0,
            replace_floor: true,
            set_age: true,
            bury: true,
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

    /// Ejecta buries older relief: a small crater under a larger crater's ejecta loses a
    /// good part of its depth; without burial the ejecta only drapes it.
    #[test]
    fn ejecta_buries_older_relief() {
        let grid = CubeSphere::new(256);
        let body_radius = 1.7374e6;
        let s = Scaling {
            simple_complex: 16_500.0,
            peak_ring: 140_000.0,
            multi_ring: 330_000.0,
            gravity: 1.62,
        };
        let small = dimensions(20_000.0, &s);
        let big = dimensions(200_000.0, &s);
        let big_center = [0.0, 0.0, 1.0];
        let (east, north) = crate::geom::local_basis(big_center);
        // The small crater sits at 1.3 big radii, under the continuous ejecta.
        let small_center =
            crate::geom::offset(big_center, east, north, 1.3 * big.radius / body_radius, 0.3);
        let relief = |layers: &LayerSet, stamper: &mut Stamper| -> f64 {
            stamper.gather(&grid, body_radius, small_center, 1.1 * small.radius);
            let h = layers.bedrock.as_slice();
            let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
            for &(idx, _, _) in &stamper.footprint {
                lo = lo.min(h[idx] as f64);
                hi = hi.max(h[idx] as f64);
            }
            hi - lo
        };
        let mut results = Vec::new();
        for bury in [true, false] {
            let mut layers = LayerSet::new(grid);
            let mut stamper = Stamper::new(&grid);
            let params = |t: f32| StampParams {
                time_gyr: t,
                unit_floor: 4,
                unit_ejecta: 5,
                azimuth: 0.0,
                angle: std::f64::consts::FRAC_PI_2,
                ellipticity: 1.0,
                amplitude: 1.0,
                replace_floor: true,
                set_age: true,
                bury,
            };
            stamp(
                &mut layers,
                &grid,
                body_radius,
                small_center,
                &ProfileTable::new(small, 0.0),
                &params(1.0),
                &mut stamper,
            );
            let before = relief(&layers, &mut stamper);
            stamp(
                &mut layers,
                &grid,
                body_radius,
                big_center,
                &ProfileTable::new(big, 0.0),
                &params(2.0),
                &mut stamper,
            );
            let after = relief(&layers, &mut stamper);
            results.push((before, after));
        }
        let (buried_before, buried_after) = results[0];
        let (draped_before, draped_after) = results[1];
        assert!((buried_before - draped_before).abs() < 1.0);
        // Draping adds the ejecta slope across the span and nothing else is lost.
        assert!(
            draped_after > 0.9 * draped_before,
            "draped {draped_before} -> {draped_after}"
        );
        // Burial removes a good part of the relief on top of the same slope.
        let lost = draped_after - buried_after;
        assert!(
            lost > 0.2 * draped_before && lost < 0.8 * draped_before,
            "buried {buried_before} -> {buried_after}, draped -> {draped_after}"
        );
    }

    /// Oblique ejecta is thicker downrange, nearly absent uprange and raised to the sides;
    /// a vertical impact is symmetric. Measured against the analytic thickness at each
    /// sampled cell's own distance, so the grid does not enter.
    #[test]
    fn oblique_ejecta_is_thicker_downrange_and_absent_uprange() {
        let grid = CubeSphere::new(512);
        let body_radius = 1.7374e6;
        let s = Scaling {
            simple_complex: 16_500.0,
            peak_ring: 140_000.0,
            multi_ring: 330_000.0,
            gravity: 1.62,
        };
        let dim = dimensions(200_000.0, &s);
        let table = ProfileTable::new(dim, 0.0);
        let center = [0.0, 0.0, 1.0];
        let (east, north) = crate::geom::local_basis(center);
        let azimuth = 1.0;
        let factor = |layers: &LayerSet, phi: f64| -> f64 {
            let d = crate::geom::offset(center, east, north, 1.5 * dim.radius / body_radius, phi);
            let cell = grid.locate(d).0;
            let r = crate::geom::angle_between(center, grid.center(cell)) * body_radius;
            layers.bedrock.get(cell) as f64 / table.height(r)
        };
        let mut out = Vec::new();
        for angle in [std::f64::consts::FRAC_PI_2, 10f64.to_radians()] {
            let mut layers = LayerSet::new(grid);
            let mut stamper = Stamper::new(&grid);
            let params = StampParams {
                time_gyr: 1.0,
                unit_floor: 4,
                unit_ejecta: 5,
                azimuth,
                angle,
                ellipticity: 1.0,
                amplitude: 1.0,
                replace_floor: true,
                set_age: true,
                bury: true,
            };
            stamp(
                &mut layers,
                &grid,
                body_radius,
                center,
                &table,
                &params,
                &mut stamper,
            );
            out.push((
                factor(&layers, azimuth),
                factor(&layers, azimuth + std::f64::consts::PI),
                factor(&layers, azimuth + std::f64::consts::FRAC_PI_2),
            ));
        }
        let (vd, vu, vs) = out[0];
        for v in [vd, vu, vs] {
            assert!((v - 1.0).abs() < 0.05, "vertical impact: {vd} {vu} {vs}");
        }
        let (od, ou, os) = out[1];
        assert!(
            od > 1.3 && ou < 0.4 && os > 1.1,
            "oblique: down {od} up {ou} side {os}"
        );
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
