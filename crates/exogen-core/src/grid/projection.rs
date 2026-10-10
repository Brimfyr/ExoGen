// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! Seamless sampling of cubesphere fields and resampling to equirectangular images.
//!
//! Bilinear interpolation in face-angle space needs the four cell centres around a sample
//! point. Within half a cell of a face edge some of those lie on a neighbouring face, so a
//! field is first padded with a one-cell ghost ring per face ([`PaddedField`]). An edge ghost
//! centre lies exactly on the neighbouring face's edge column or row (the equiangular mapping
//! guarantees this), so its value is a one-dimensional linear interpolation along that line
//! and is second-order accurate. A corner ghost lies on the third face at the cube corner and
//! is interpolated there with clamping; its error is first order but confined to the cells
//! touching the eight cube corners.
//!
//! Pixel `(px, py)` of a `width x height` equirectangular image covers longitude
//! `-pi + 2 pi (px + 1/2) / width` and latitude `pi/2 - pi (py + 1/2) / height`.

use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI};

use rayon::prelude::*;

use super::cubesphere::{Cell, CubeSphere};
use super::direction_from_lon_lat;

/// A cubesphere field with a one-cell ghost ring around each face.
pub struct PaddedField {
    f: usize,
    data: Vec<f32>,
}

impl PaddedField {
    pub fn new(grid: &CubeSphere, values: &[f32]) -> Self {
        assert_eq!(values.len(), grid.cell_count());
        let f = grid.face_size() as usize;
        let p = f + 2;
        let mut data = vec![0.0f32; 6 * p * p];
        for face in 0..6u8 {
            for j in 0..f {
                for i in 0..f {
                    let c = Cell {
                        face,
                        i: i as u32,
                        j: j as u32,
                    };
                    data[(face as usize * p + j + 1) * p + i + 1] = values[grid.index(c)];
                }
            }
        }
        let delta = grid.delta();
        for face in 0..6u8 {
            for k in -1..=(f as i64) {
                for (gi, gj) in [(-1i64, k), (f as i64, k), (k, -1i64), (k, f as i64)] {
                    let alpha = -FRAC_PI_4 + (gi as f64 + 0.5) * delta;
                    let beta = -FRAC_PI_4 + (gj as f64 + 0.5) * delta;
                    let d = CubeSphere::direction(face, alpha, beta);
                    let idx = (face as usize * p + (gj + 1) as usize) * p + (gi + 1) as usize;
                    data[idx] = interior_bilinear(grid, values, d);
                }
            }
        }
        Self { f, data }
    }

    /// Value at in-face indices in `-1..=f`.
    #[inline]
    pub fn get(&self, face: u8, i: i64, j: i64) -> f32 {
        let p = self.f + 2;
        self.data[(face as usize * p + (j + 1) as usize) * p + (i + 1) as usize]
    }

    pub fn face_size(&self) -> usize {
        self.f
    }
}

/// Bilinear interpolation using only the cells of the face that contains `d`, clamping the
/// stencil to the face. Exact whenever the four stencil cells lie inside that face.
fn interior_bilinear(grid: &CubeSphere, values: &[f32], d: [f64; 3]) -> f32 {
    let (c, x, y) = grid.locate(d);
    let f = grid.face_size() as i64;
    let i0 = ((x - 0.5).floor() as i64).clamp(0, f - 2);
    let j0 = ((y - 0.5).floor() as i64).clamp(0, f - 2);
    let fx = (x - 0.5 - i0 as f64).clamp(0.0, 1.0);
    let fy = (y - 0.5 - j0 as f64).clamp(0.0, 1.0);
    let v = |i: i64, j: i64| {
        values[grid.index(Cell {
            face: c.face,
            i: i as u32,
            j: j as u32,
        })] as f64
    };
    (v(i0, j0) * (1.0 - fx) * (1.0 - fy)
        + v(i0 + 1, j0) * fx * (1.0 - fy)
        + v(i0, j0 + 1) * (1.0 - fx) * fy
        + v(i0 + 1, j0 + 1) * fx * fy) as f32
}

/// Seamless bilinear sample of a padded field at a unit direction.
pub fn sample_bilinear(grid: &CubeSphere, field: &PaddedField, d: [f64; 3]) -> f32 {
    let (c, x, y) = grid.locate(d);
    let i0 = (x - 0.5).floor() as i64;
    let j0 = (y - 0.5).floor() as i64;
    let fx = x - 0.5 - i0 as f64;
    let fy = y - 0.5 - j0 as f64;
    let v = |i: i64, j: i64| field.get(c.face, i, j) as f64;
    (v(i0, j0) * (1.0 - fx) * (1.0 - fy)
        + v(i0 + 1, j0) * fx * (1.0 - fy)
        + v(i0, j0 + 1) * (1.0 - fx) * fy
        + v(i0 + 1, j0 + 1) * fx * fy) as f32
}

/// Resample a field to an equirectangular image of `width x height` samples.
pub fn to_equirectangular(
    grid: &CubeSphere,
    values: &[f32],
    width: usize,
    height: usize,
) -> Vec<f32> {
    let field = PaddedField::new(grid, values);
    let mut out = vec![0.0f32; width * height];
    out.par_chunks_mut(width).enumerate().for_each(|(py, row)| {
        let lat = FRAC_PI_2 - PI * (py as f64 + 0.5) / height as f64;
        for (px, o) in row.iter_mut().enumerate() {
            let lon = -PI + 2.0 * PI * (px as f64 + 0.5) / width as f64;
            *o = sample_bilinear(grid, &field, direction_from_lon_lat(lon, lat));
        }
    });
    out
}

/// Angular distance in radians from a direction to the nearest cube corner.
pub fn angle_to_nearest_cube_corner(d: [f64; 3]) -> f64 {
    let s = 1.0 / 3f64.sqrt();
    let dot = (d[0].abs() + d[1].abs() + d[2].abs()) * s;
    dot.clamp(-1.0, 1.0).acos()
}

/// A smooth analytic field on the sphere used by tests and spikes. It has structure at
/// several scales and no symmetry aligned with the cube faces.
pub fn test_field(d: [f64; 3]) -> f64 {
    let [x, y, z] = d;
    (3.0 * x).sin() * (2.0 * y).cos() + z * z * z + 0.5 * (5.0 * (x + y) + 1.0).sin()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Max error over an equirectangular resample, split into pixels near a cube corner
    /// (within 1.5 cells) and the rest.
    fn resample_errors(f: u32, w: usize, h: usize) -> (f64, f64) {
        let grid = CubeSphere::new(f);
        let values: Vec<f32> = grid
            .cells()
            .map(|c| test_field(grid.center(c)) as f32)
            .collect();
        let img = to_equirectangular(&grid, &values, w, h);
        let (mut max_err, mut max_corner) = (0.0f64, 0.0f64);
        for py in 0..h {
            let lat = FRAC_PI_2 - PI * (py as f64 + 0.5) / h as f64;
            for px in 0..w {
                let lon = -PI + 2.0 * PI * (px as f64 + 0.5) / w as f64;
                let d = direction_from_lon_lat(lon, lat);
                let err = (img[py * w + px] as f64 - test_field(d)).abs();
                if angle_to_nearest_cube_corner(d) < 1.5 * grid.delta() {
                    max_corner = max_corner.max(err);
                } else {
                    max_err = max_err.max(err);
                }
            }
        }
        (max_err, max_corner)
    }

    #[test]
    fn resampling_is_second_order_and_seamless() {
        let mut previous = f64::INFINITY;
        for f in [32u32, 64, 128] {
            let (max_err, max_corner) = resample_errors(f, 512, 256);
            // Second order: halving the cell size divides the error by about four.
            assert!(
                max_err < previous * 0.35,
                "f={f} max error {max_err} did not fall from {previous}"
            );
            // Corner cells are first order but must stay small.
            assert!(max_corner < 0.05, "f={f} corner error {max_corner}");
            previous = max_err;
        }
        assert!(previous < 2e-3, "max error at f=128 is {previous}");
    }

    /// Worst error of the edge ghosts (corners excluded) against the analytic field.
    fn worst_edge_ghost_error(f: u32) -> f64 {
        let grid = CubeSphere::new(f);
        let values: Vec<f32> = grid
            .cells()
            .map(|c| test_field(grid.center(c)) as f32)
            .collect();
        let field = PaddedField::new(&grid, &values);
        let delta = grid.delta();
        let n = grid.face_size() as i64;
        let mut worst = 0.0f64;
        for face in 0..6u8 {
            for k in 0..n {
                for (gi, gj) in [(-1i64, k), (n, k), (k, -1i64), (k, n)] {
                    let alpha = -FRAC_PI_4 + (gi as f64 + 0.5) * delta;
                    let beta = -FRAC_PI_4 + (gj as f64 + 0.5) * delta;
                    let d = CubeSphere::direction(face, alpha, beta);
                    worst = worst.max((field.get(face, gi, gj) as f64 - test_field(d)).abs());
                }
            }
        }
        worst
    }

    #[test]
    fn edge_ghosts_are_second_order_accurate() {
        // Edge ghost centres lie on the neighbouring face's edge line, so their values are
        // linear interpolations there: the error must fall by about four per halving of the
        // cell size, and be small in absolute terms at f=128.
        let e32 = worst_edge_ghost_error(32);
        let e64 = worst_edge_ghost_error(64);
        let e128 = worst_edge_ghost_error(128);
        assert!(
            e64 < e32 * 0.35,
            "ghost error {e32} -> {e64} is not second order"
        );
        assert!(
            e128 < e64 * 0.35,
            "ghost error {e64} -> {e128} is not second order"
        );
        assert!(e128 < 1e-3, "ghost error at f=128 is {e128}");
    }
}
