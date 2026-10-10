// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! Spherical grids (ADR 0001): the equiangular cubesphere for all gridded fields and
//! the projection to equirectangular maps for export.

pub mod cubesphere;
pub mod projection;

pub use cubesphere::{Cell, CubeSphere};

/// Unit direction from longitude and latitude in radians. Longitude increases eastward
/// from the +X axis; +Z is the north pole.
pub fn direction_from_lon_lat(lon: f64, lat: f64) -> [f64; 3] {
    let (sin_lat, cos_lat) = lat.sin_cos();
    let (sin_lon, cos_lon) = lon.sin_cos();
    [cos_lat * cos_lon, cos_lat * sin_lon, sin_lat]
}

/// Longitude and latitude in radians of a unit direction.
pub fn lon_lat_from_direction(d: [f64; 3]) -> (f64, f64) {
    (d[1].atan2(d[0]), d[2].clamp(-1.0, 1.0).asin())
}

pub(crate) fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

pub(crate) fn normalize(v: [f64; 3]) -> [f64; 3] {
    let n = dot(v, v).sqrt();
    [v[0] / n, v[1] / n, v[2] / n]
}
