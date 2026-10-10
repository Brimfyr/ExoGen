// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! Shaded relief for judging maps: Lambert lighting of an equirectangular height map from
//! one fixed sun, applied to the albedo. It is a viewing aid, not the renderer's lighting,
//! so the sun is the same everywhere on the map.

/// Sun azimuth clockwise from north and elevation, degrees, and the ambient share of the
/// light, the cartographic convention of light from the upper left.
pub const SHADE_SUN_AZIMUTH_DEG: f64 = 315.0;
pub const SHADE_SUN_ELEVATION_DEG: f64 = 30.0;
pub const SHADE_AMBIENT: f32 = 0.3;
/// Mean albedo the preview is exposed to, so dark worlds read as well as bright ones.
pub const SHADE_PREVIEW_MEAN: f32 = 0.5;

/// Exposure gain that brings the mean albedo up to the preview mean; never darkens.
pub fn preview_gain(albedo: impl Iterator<Item = f32>) -> f32 {
    let (mut sum, mut n) = (0.0f64, 0usize);
    for a in albedo {
        if a.is_finite() {
            sum += a as f64;
            n += 1;
        }
    }
    if n == 0 || sum <= 0.0 {
        return 1.0;
    }
    (SHADE_PREVIEW_MEAN as f64 / (sum / n as f64)).max(1.0) as f32
}

/// Lambert shading, 0 to 1, of an equirectangular height map in metres on a body of
/// radius `body_radius`, with `width` columns and `rows` rows.
pub fn hillshade(height: &[f32], width: usize, rows: usize, body_radius: f64) -> Vec<f32> {
    let az = SHADE_SUN_AZIMUTH_DEG.to_radians();
    let el = SHADE_SUN_ELEVATION_DEG.to_radians();
    let sun = [az.sin() * el.cos(), az.cos() * el.cos(), el.sin()];
    let px_y = std::f64::consts::PI * body_radius / rows as f64;
    let mut out = vec![1.0f32; width * rows];
    for y in 0..rows {
        let lat =
            std::f64::consts::FRAC_PI_2 - std::f64::consts::PI * (y as f64 + 0.5) / rows as f64;
        let px_x = 2.0 * std::f64::consts::PI * body_radius * lat.cos().max(0.05) / width as f64;
        let (yn, ys) = (y.saturating_sub(1), (y + 1).min(rows - 1));
        for x in 0..width {
            let (xw, xe) = ((x + width - 1) % width, (x + 1) % width);
            let sx = (height[y * width + xe] as f64 - height[y * width + xw] as f64) / (2.0 * px_x);
            let sy = (height[yn * width + x] as f64 - height[ys * width + x] as f64) / (2.0 * px_y);
            let len = (sx * sx + sy * sy + 1.0).sqrt();
            let n = [-sx / len, -sy / len, 1.0 / len];
            let dot = n[0] * sun[0] + n[1] * sun[1] + n[2] * sun[2];
            out[y * width + x] = dot.max(0.0) as f32;
        }
    }
    out
}

/// Albedo, exposed by `gain`, lit by the shade: ambient plus the direct share.
pub fn lit(albedo: [f32; 3], shade: f32, gain: f32) -> [f32; 3] {
    let k = gain * (SHADE_AMBIENT + (1.0 - SHADE_AMBIENT) * shade);
    [
        (albedo[0] * k).min(1.0),
        (albedo[1] * k).min(1.0),
        (albedo[2] * k).min(1.0),
    ]
}
