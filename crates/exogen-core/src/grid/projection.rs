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
#[derive(Clone, Debug)]
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
    interior_bilinear_with(grid, &|c: Cell| values[grid.index(c)], d)
}

/// The same, reading cells through `value`.
fn interior_bilinear_with(grid: &CubeSphere, value: &dyn Fn(Cell) -> f32, d: [f64; 3]) -> f32 {
    let (c, x, y) = grid.locate(d);
    let f = grid.face_size() as i64;
    let i0 = ((x - 0.5).floor() as i64).clamp(0, f - 2);
    let j0 = ((y - 0.5).floor() as i64).clamp(0, f - 2);
    let fx = (x - 0.5 - i0 as f64).clamp(0.0, 1.0);
    let fy = (y - 0.5 - j0 as f64).clamp(0.0, 1.0);
    let v = |i: i64, j: i64| {
        value(Cell {
            face: c.face,
            i: i as u32,
            j: j as u32,
        }) as f64
    };
    (v(i0, j0) * (1.0 - fx) * (1.0 - fy)
        + v(i0 + 1, j0) * fx * (1.0 - fy)
        + v(i0, j0 + 1) * (1.0 - fx) * fy
        + v(i0 + 1, j0 + 1) * fx * fy) as f32
}

/// Seamless bilinear sample of a padded field at a unit direction.
pub fn sample_bilinear(grid: &CubeSphere, field: &PaddedField, d: [f64; 3]) -> f32 {
    sample_bilinear_with(grid, &|face, i, j| field.get(face, i, j), d)
}

/// The same, reading padded values through `get` (in-face indices in `-1..=f`).
fn sample_bilinear_with(grid: &CubeSphere, get: &dyn Fn(u8, i64, i64) -> f32, d: [f64; 3]) -> f32 {
    let (c, x, y) = grid.locate(d);
    let i0 = (x - 0.5).floor() as i64;
    let j0 = (y - 0.5).floor() as i64;
    let fx = x - 0.5 - i0 as f64;
    let fy = y - 0.5 - j0 as f64;
    let v = |i: i64, j: i64| get(c.face, i, j) as f64;
    (v(i0, j0) * (1.0 - fx) * (1.0 - fy)
        + v(i0 + 1, j0) * fx * (1.0 - fy)
        + v(i0, j0 + 1) * (1.0 - fx) * fy
        + v(i0 + 1, j0 + 1) * fx * fy) as f32
}

/// A field on a cubesphere read at a direction by Steffen's monotone cubic (F291): along the
/// face's rows over the four columns about the direction, then across the four rows, a cell
/// past the face's edge read where its centre falls, by the bilinear reading the ghost ring
/// takes; so that, read so, a field shows none of its cells' edges and passes none of the
/// values about it.
pub fn sample_monotone(grid: &CubeSphere, values: &[f32], d: [f64; 3]) -> f32 {
    let (c, x, y) = grid.locate(d);
    let f = grid.face_size() as i64;
    let (i0, j0) = ((x - 0.5).floor() as i64, (y - 0.5).floor() as i64);
    let (tx, ty) = (x - 0.5 - i0 as f64, y - 0.5 - j0 as f64);
    let delta = grid.delta();
    let cell = |i: i64, j: i64| -> f64 {
        if (0..f).contains(&i) && (0..f).contains(&j) {
            values[grid.index(Cell {
                face: c.face,
                i: i as u32,
                j: j as u32,
            })] as f64
        } else {
            let alpha = -FRAC_PI_4 + (i as f64 + 0.5) * delta;
            let beta = -FRAC_PI_4 + (j as f64 + 0.5) * delta;
            interior_bilinear(grid, values, CubeSphere::direction(c.face, alpha, beta)) as f64
        }
    };
    let rows: [f64; 4] = std::array::from_fn(|b| {
        let j = j0 + b as i64 - 1;
        crate::geom::monotone(std::array::from_fn(|a| cell(i0 + a as i64 - 1, j)), tx)
    });
    crate::geom::monotone(rows, ty) as f32
}

/// A field that is zero except on a few cells, sampled exactly as a `PaddedField` of the
/// same values would be: ghost cells across face edges are built on demand by the same
/// in-face interpolation the padding uses.
pub struct SparseField {
    grid: CubeSphere,
    values: std::collections::HashMap<usize, f32>,
}

impl SparseField {
    pub fn new(grid: CubeSphere, entries: impl IntoIterator<Item = (usize, f32)>) -> Self {
        Self {
            grid,
            values: entries.into_iter().collect(),
        }
    }

    fn value(&self, c: Cell) -> f32 {
        self.values.get(&self.grid.index(c)).copied().unwrap_or(0.0)
    }

    /// Value at in-face indices in `-1..=f`, ghosts included.
    pub fn get(&self, face: u8, i: i64, j: i64) -> f32 {
        let f = self.grid.face_size() as i64;
        if (0..f).contains(&i) && (0..f).contains(&j) {
            return self.value(Cell {
                face,
                i: i as u32,
                j: j as u32,
            });
        }
        let delta = self.grid.delta();
        let alpha = -FRAC_PI_4 + (i as f64 + 0.5) * delta;
        let beta = -FRAC_PI_4 + (j as f64 + 0.5) * delta;
        let d = CubeSphere::direction(face, alpha, beta);
        interior_bilinear_with(&self.grid, &|c: Cell| self.value(c), d)
    }
}

/// Seamless bilinear sample of a sparse field at a unit direction.
pub fn sample_bilinear_sparse(grid: &CubeSphere, field: &SparseField, d: [f64; 3]) -> f32 {
    sample_bilinear_with(grid, &|face, i, j| field.get(face, i, j), d)
}

/// Resample a field to an equirectangular image of `width x height` pixels, each the mean of
/// the field over its footprint (F424): read between the cells at as many points across the
/// pixel, each way, as cells cross it, so that a pixel wider than the cells holds their mean
/// and not the few nearest its centre, the cells' own variation aliased into the pixels as
/// speckle; a pixel a cell wide or less read at its centre alone.
pub fn to_equirectangular(
    grid: &CubeSphere,
    values: &[f32],
    width: usize,
    height: usize,
) -> Vec<f32> {
    let field = PaddedField::new(grid, values);
    let cell = FRAC_PI_2 / grid.face_size() as f64;
    let (dlon, dlat) = (2.0 * PI / width as f64, PI / height as f64);
    let across = |angle: f64| ((angle / cell) - 1e-9).ceil().max(1.0) as usize;
    let mut out = vec![0.0f32; width * height];
    out.par_chunks_mut(width).enumerate().for_each(|(py, row)| {
        let top = FRAC_PI_2 - dlat * py as f64;
        let middle = top - 0.5 * dlat;
        let (rows, columns) = (across(dlat), across(dlon * middle.cos()));
        for (px, o) in row.iter_mut().enumerate() {
            let left = -PI + dlon * px as f64;
            let mut sum = 0.0f64;
            for i in 0..rows {
                let lat = top - dlat * (i as f64 + 0.5) / rows as f64;
                for j in 0..columns {
                    let lon = left + dlon * (j as f64 + 0.5) / columns as f64;
                    sum += sample_bilinear(grid, &field, direction_from_lon_lat(lon, lat)) as f64;
                }
            }
            *o = (sum / (rows * columns) as f64) as f32;
        }
    });
    out
}

/// Resample a map of each cell's own that crosses a threshold where its kinds' shares are a half
/// (F419), the height crossing the sea's level, to an equirectangular image of `width x height`
/// pixels, each the mean over its footprint (F424) of the cells on its greater side (F426, ADR
/// 0009). `shares`: the kinds' shares of each cell in their order, the sea's and then the
/// lakes', the rest the land's. A cell is on the first kind's side whose share of it is over a
/// half, else on the rest's; a pixel on the first kind's whose share of it, over its footprint,
/// is a half or more, as the masks draw it, else on the rest's. So that a pixel the sea covers
/// over half holds its shelf's floor and not that floor's mean with the land beside it, and the
/// map crosses its threshold where the masks' shares are a half and not a pixel off. A pixel
/// whose side has no cell within its reach holds the mean of all of them.
pub fn to_equirectangular_by_sides(
    grid: &CubeSphere,
    values: &[f32],
    shares: &[&[f32]],
    width: usize,
    height: usize,
) -> Vec<f32> {
    let side: Vec<usize> = (0..values.len())
        .map(|c| {
            shares
                .iter()
                .position(|s| s[c] > 0.5)
                .unwrap_or(shares.len())
        })
        .collect();
    let drawn: Vec<Vec<f32>> = shares
        .iter()
        .map(|s| to_equirectangular(grid, s, width, height))
        .collect();
    let pixel_side: Vec<usize> = (0..width * height)
        .map(|p| {
            drawn
                .iter()
                .position(|s| s[p] >= 0.5)
                .unwrap_or(shares.len())
        })
        .collect();
    let mut out = to_equirectangular(grid, values, width, height);
    for kind in 0..=shares.len() {
        if !pixel_side.contains(&kind) {
            continue;
        }
        let on: Vec<f32> = side
            .iter()
            .map(|&s| if s == kind { 1.0 } else { 0.0 })
            .collect();
        let weighted: Vec<f32> = values.iter().zip(&on).map(|(v, w)| v * w).collect();
        let sum = to_equirectangular(grid, &weighted, width, height);
        let weight = to_equirectangular(grid, &on, width, height);
        for p in (0..width * height).filter(|&p| pixel_side[p] == kind && weight[p] > 0.0) {
            out[p] = sum[p] / weight[p];
        }
    }
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

    /// F424: a pixel wider than the cells holds their mean: a field alternating +1 and -1
    /// from cell to cell over a 128-cell face, resampled to 64 by 32 pixels, some eight cells
    /// across each at the equator. Read at each pixel's centre alone, the pixels speckle near
    /// the cells' own +1 and -1; read over its footprint, each holds their mean, near none.
    /// And a smooth field so resampled stands within its own variation over a pixel.
    #[test]
    fn a_pixel_wider_than_the_cells_holds_their_mean() {
        let grid = CubeSphere::new(128);
        let checker: Vec<f32> = grid
            .cells()
            .map(|c| if (c.i + c.j) % 2 == 0 { 1.0 } else { -1.0 })
            .collect();
        let (w, h) = (64, 32);
        let img = to_equirectangular(&grid, &checker, w, h);
        let field = PaddedField::new(&grid, &checker);
        let rms = |v: &mut dyn Iterator<Item = f64>| {
            let (sum, count) = v.fold((0.0, 0), |(s, n), x| (s + x * x, n + 1));
            (sum / count as f64).sqrt()
        };
        let point = rms(&mut (0..w * h).map(|k| {
            let (px, py) = (k % w, k / w);
            let lat = FRAC_PI_2 - PI * (py as f64 + 0.5) / h as f64;
            let lon = -PI + 2.0 * PI * (px as f64 + 0.5) / w as f64;
            sample_bilinear(&grid, &field, direction_from_lon_lat(lon, lat)) as f64
        }));
        let mean = rms(&mut img.iter().map(|&v| v as f64));
        let smooth: Vec<f32> = grid
            .cells()
            .map(|c| test_field(grid.center(c)) as f32)
            .collect();
        let read = to_equirectangular(&grid, &smooth, w, h);
        let off = (0..w * h)
            .map(|k| {
                let (px, py) = (k % w, k / w);
                let lat = FRAC_PI_2 - PI * (py as f64 + 0.5) / h as f64;
                let lon = -PI + 2.0 * PI * (px as f64 + 0.5) / w as f64;
                (read[k] as f64 - test_field(direction_from_lon_lat(lon, lat))).abs()
            })
            .fold(0.0f64, f64::max);
        println!(
            "cells alternating by one: the pixels' spread {point:.3} read at their centres, {mean:.3} over their footprints; a smooth field off its value at the pixels' centres by {off:.3} at most"
        );
        assert!(point > 0.15, "{point}");
        assert!(mean < 0.25 * point, "{mean} {point}");
        assert!(off < 0.1, "{off}");
    }

    /// F426: each pixel its greater side's mean. A coast along a great circle off the grid's
    /// axes, its sea a shelf 10 m deep against land 200 m high, and within the land a lake 5 m
    /// deep over a floor at 150 m, each cell's share of each its distance from its margin
    /// through a logistic over a cell, each cell its own side's height (F419); resampled to
    /// pixels some four cells across. Over its footprint's mean a pixel the sea covers over
    /// half stands over the level where its shelf's floor meets the land's height, a pixel the
    /// lake covers over half over the lake's; over its greater side's, every pixel stands
    /// under the sea's level exactly where the sea's share of it is a half or more, and under
    /// the lake's where the lake's is; a pixel wholly of one kind as before.
    #[test]
    fn a_pixel_holds_its_greater_side_s_mean() {
        use crate::geom::{angle_between, cross, dot, normalize};
        let grid = CubeSphere::new(64);
        let cell = FRAC_PI_2 / 64.0;
        let coast = normalize([0.31, 0.17, 0.93]);
        let along = normalize(cross(coast, [0.0, 0.0, 1.0]));
        let inland = 12.0 * cell;
        let centre = normalize([
            coast[0] * inland.sin() + along[0] * inland.cos(),
            coast[1] * inland.sin() + along[1] * inland.cos(),
            coast[2] * inland.sin() + along[2] * inland.cos(),
        ]);
        let logistic = |x: f64| (1.0 / (1.0 + (2.0 * x).exp())) as f32;
        let (mut sea, mut lake, mut ground) = (vec![], vec![], vec![]);
        for c in grid.cells() {
            let d = grid.center(c);
            let s = logistic(dot(d, coast).asin() / cell);
            let l = logistic((angle_between(d, centre) - 4.0 * cell) / cell);
            sea.push(s);
            lake.push(l);
            ground.push(if s > 0.5 {
                -10.0
            } else if l > 0.5 {
                150.0
            } else {
                200.0
            });
        }
        let (w, h) = (64, 32);
        let shares: [&[f32]; 2] = [&sea, &lake];
        let mean = to_equirectangular(&grid, &ground, w, h);
        let sides = to_equirectangular_by_sides(&grid, &ground, &shares, w, h);
        let sea_px = to_equirectangular(&grid, &sea, w, h);
        let lake_px = to_equirectangular(&grid, &lake, w, h);
        let apart = |img: &[f32]| {
            let at_sea = (0..w * h)
                .filter(|&p| (img[p] < 0.0) != (sea_px[p] >= 0.5))
                .count();
            let at_lake = (0..w * h)
                .filter(|&p| sea_px[p] < 0.5 && (img[p] < 155.0) != (lake_px[p] >= 0.5))
                .count();
            (at_sea, at_lake)
        };
        let whole = (0..w * h)
            .filter(|&p| {
                let one = |s: f32| s <= 1e-6 || s >= 1.0 - 1e-6;
                one(sea_px[p]) && one(lake_px[p])
            })
            .map(|p| (sides[p] - mean[p]).abs())
            .fold(0.0f32, f32::max);
        println!(
            "pixels apart from their shares' half, the sea's and the lake's: {:?} over the footprint's mean, {:?} over the greater side's; a pixel wholly of one kind off by {whole}",
            apart(&mean),
            apart(&sides)
        );
        let (sea_mean, lake_mean) = apart(&mean);
        assert!(sea_mean >= 5 && lake_mean >= 2, "{sea_mean} {lake_mean}");
        assert_eq!(apart(&sides), (0, 0));
        assert!(whole < 1e-3, "{whole}");
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

    /// F291: a field read by the monotone cubic passes none of the values about it, a high
    /// cell among cells of none reading its own value at its centre and nothing below none
    /// about it, across the face's edges as within them.
    #[test]
    fn a_monotone_read_passes_no_cell_s_value() {
        let grid = CubeSphere::new(16);
        let mut values = vec![0.0f32; grid.cell_count()];
        let at = grid.index(grid.locate([1.0, 0.02, 0.0]).0);
        values[at] = 100.0;
        let edge = grid.index(grid.locate(crate::geom::normalize([1.0, 0.98, 0.0])).0);
        values[edge] = 50.0;
        let (mut least, mut most) = (f32::MAX, f32::MIN);
        for a in 0..=400 {
            for b in 0..=40 {
                let lon = -1.2 + 2.4 * a as f64 / 400.0;
                let lat = -0.3 + 0.6 * b as f64 / 40.0;
                let v = sample_monotone(&grid, &values, direction_from_lon_lat(lon, lat));
                (least, most) = (least.min(v), most.max(v));
            }
        }
        let centre = sample_monotone(&grid, &values, grid.center(grid.cell(at)));
        assert!((centre - 100.0).abs() < 1e-3, "{centre}");
        assert!(least >= 0.0 && most <= 100.0 + 1e-3, "{least} {most}");
    }
}
