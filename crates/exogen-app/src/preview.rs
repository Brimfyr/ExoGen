// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! Preview images: equirectangular channel renders and a software-rendered globe.
//!
//! The globe is ray-cast on the CPU into an image, which is enough for a preview and keeps
//! the renderer on the deterministic CPU path. A wgpu globe with atmosphere scattering is
//! Phase 7 work.

use eframe::egui;
use exogen_core::grid::projection::to_equirectangular;
use exogen_core::stages::FinishOutput;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Channel {
    Albedo,
    Shaded,
    Height,
    SurfaceAge,
    Liquid,
    Crust,
    Unit,
    Maturity,
    Rays,
}

impl Channel {
    pub const ALL: [Channel; 9] = [
        Channel::Albedo,
        Channel::Shaded,
        Channel::Height,
        Channel::SurfaceAge,
        Channel::Liquid,
        Channel::Crust,
        Channel::Unit,
        Channel::Maturity,
        Channel::Rays,
    ];

    pub fn label(self) -> &'static str {
        match self {
            Channel::Albedo => "albedo",
            Channel::Shaded => "shaded relief",
            Channel::Height => "height",
            Channel::SurfaceAge => "surface age",
            Channel::Liquid => "liquid depth",
            Channel::Crust => "crust thickness",
            Channel::Unit => "unit id",
            Channel::Maturity => "maturity",
            Channel::Rays => "rays",
        }
    }
}

/// A linear-light RGB equirectangular image.
pub struct EquirectRgb {
    pub width: usize,
    pub height: usize,
    pub rgb: Vec<[f32; 3]>,
    /// Range of the scalar shown, for the readout.
    pub range: (f32, f32),
}

fn ramp(t: f32, stops: &[(f32, [f32; 3])]) -> [f32; 3] {
    let t = t.clamp(0.0, 1.0);
    for w in stops.windows(2) {
        let (t0, c0) = w[0];
        let (t1, c1) = w[1];
        if t <= t1 {
            let u = if t1 > t0 { (t - t0) / (t1 - t0) } else { 0.0 };
            return [
                c0[0] + (c1[0] - c0[0]) * u,
                c0[1] + (c1[1] - c0[1]) * u,
                c0[2] + (c1[2] - c0[2]) * u,
            ];
        }
    }
    stops[stops.len() - 1].1
}

const HEIGHT_LAND: [(f32, [f32; 3]); 5] = [
    (0.0, [0.05, 0.12, 0.04]),
    (0.3, [0.25, 0.30, 0.10]),
    (0.6, [0.40, 0.28, 0.16]),
    (0.85, [0.45, 0.42, 0.40]),
    (1.0, [0.90, 0.90, 0.92]),
];
const HEIGHT_SEA: [(f32, [f32; 3]); 2] = [(0.0, [0.01, 0.02, 0.08]), (1.0, [0.05, 0.15, 0.30])];
const AGE_RAMP: [(f32, [f32; 3]); 4] = [
    (0.0, [0.02, 0.02, 0.08]),
    (0.35, [0.10, 0.20, 0.45]),
    (0.7, [0.15, 0.55, 0.35]),
    (1.0, [0.95, 0.90, 0.20]),
];
const GREY: [(f32, [f32; 3]); 2] = [(0.0, [0.02, 0.02, 0.02]), (1.0, [0.95, 0.95, 0.95])];

fn unit_colour(id: u16) -> [f32; 3] {
    match id {
        0 => [0.05, 0.05, 0.05],
        1 => [0.30, 0.25, 0.20],
        2 => [0.20, 0.45, 0.25],
        _ => {
            let h = (id as f32 * 0.618) % 1.0;
            ramp(h, &AGE_RAMP)
        }
    }
}

/// Render a channel of a finished state to an equirectangular image.
pub fn build_equirect(
    out: &FinishOutput,
    channel: Channel,
    width: usize,
    body_radius: f64,
) -> EquirectRgb {
    let height = width / 2;
    let grid = out.grid;
    let n = width * height;
    let mut rgb = vec![[0.0f32; 3]; n];
    let mut range = (0.0f32, 1.0f32);
    match channel {
        Channel::Albedo | Channel::Shaded => {
            for ch in 0..3 {
                let plane: Vec<f32> = out.albedo.as_slice().iter().map(|c| c[ch]).collect();
                let eq = to_equirectangular(&grid, &plane, width, height);
                for (k, v) in eq.iter().enumerate() {
                    rgb[k][ch] = *v;
                }
            }
            if channel == Channel::Shaded {
                let eq = to_equirectangular(&grid, out.height.as_slice(), width, height);
                let shade = exogen_core::shade::hillshade(&eq, width, height, body_radius);
                let gain = exogen_core::shade::preview_gain(
                    rgb.iter().map(|c| (c[0] + c[1] + c[2]) / 3.0),
                );
                for (k, s) in shade.iter().enumerate() {
                    rgb[k] = exogen_core::shade::lit(rgb[k], *s, gain);
                }
            }
        }
        Channel::Height => {
            let eq = to_equirectangular(&grid, out.height.as_slice(), width, height);
            let (lo, hi) = eq
                .iter()
                .fold((f32::INFINITY, f32::NEG_INFINITY), |(a, b), &v| {
                    (a.min(v), b.max(v))
                });
            range = (lo, hi);
            let sea = out.sea_level as f32;
            for (k, &h) in eq.iter().enumerate() {
                rgb[k] = if out.sea_level.is_finite() && h < sea {
                    ramp((h - lo) / (sea - lo).max(1.0), &HEIGHT_SEA)
                } else if out.sea_level.is_finite() {
                    ramp((h - sea) / (hi - sea).max(1.0), &HEIGHT_LAND)
                } else {
                    ramp((h - lo) / (hi - lo).max(1.0), &HEIGHT_LAND)
                };
            }
        }
        Channel::SurfaceAge => {
            let eq = to_equirectangular(&grid, out.layers.surface_age.as_slice(), width, height);
            let hi = eq.iter().cloned().fold(0.0f32, f32::max).max(1e-6);
            range = (0.0, hi);
            for (k, &a) in eq.iter().enumerate() {
                rgb[k] = ramp(a / hi, &AGE_RAMP);
            }
        }
        Channel::Liquid => {
            let eq = to_equirectangular(&grid, out.layers.liquid.as_slice(), width, height);
            let hi = eq.iter().cloned().fold(0.0f32, f32::max).max(1e-6);
            range = (0.0, hi);
            for (k, &d) in eq.iter().enumerate() {
                rgb[k] = if d > 0.0 {
                    ramp(1.0 - d / hi, &HEIGHT_SEA)
                } else {
                    [0.12, 0.10, 0.08]
                };
            }
        }
        Channel::Crust => {
            let eq =
                to_equirectangular(&grid, out.layers.crust_thickness.as_slice(), width, height);
            let (lo, hi) = eq
                .iter()
                .fold((f32::INFINITY, f32::NEG_INFINITY), |(a, b), &v| {
                    (a.min(v), b.max(v))
                });
            range = (lo, hi);
            for (k, &c) in eq.iter().enumerate() {
                rgb[k] = ramp((c - lo) / (hi - lo).max(1.0), &GREY);
            }
        }
        Channel::Maturity => {
            let eq = to_equirectangular(&grid, out.maturity.as_slice(), width, height);
            for (k, &m) in eq.iter().enumerate() {
                rgb[k] = ramp(m, &GREY);
            }
        }
        Channel::Rays => {
            let eq = to_equirectangular(&grid, out.rays.as_slice(), width, height);
            for (k, &v) in eq.iter().enumerate() {
                rgb[k] = ramp(v, &AGE_RAMP);
            }
        }
        Channel::Unit => {
            // Nearest-cell lookup for labels.
            for py in 0..height {
                let lat = std::f64::consts::FRAC_PI_2
                    - std::f64::consts::PI * (py as f64 + 0.5) / height as f64;
                for px in 0..width {
                    let lon = -std::f64::consts::PI
                        + 2.0 * std::f64::consts::PI * (px as f64 + 0.5) / width as f64;
                    let d = exogen_core::grid::direction_from_lon_lat(lon, lat);
                    let cell = grid.locate(d).0;
                    rgb[py * width + px] = unit_colour(out.layers.unit_id.get(cell));
                }
            }
        }
    }
    EquirectRgb {
        width,
        height,
        rgb,
        range,
    }
}

pub fn srgb8(v: f32) -> u8 {
    exogen_core::export::srgb_encode(v)
}

/// Convert to an egui image, sRGB encoded.
pub fn to_color_image(img: &EquirectRgb, exposure: f32) -> egui::ColorImage {
    let mut pixels = Vec::with_capacity(img.rgb.len());
    for c in &img.rgb {
        pixels.push(egui::Color32::from_rgb(
            srgb8(c[0] * exposure),
            srgb8(c[1] * exposure),
            srgb8(c[2] * exposure),
        ));
    }
    egui::ColorImage {
        size: [img.width, img.height],
        source_size: egui::vec2(img.width as f32, img.height as f32),
        pixels,
    }
}

/// Bilinear sample of the equirectangular image at a unit direction.
fn sample(img: &EquirectRgb, d: [f64; 3]) -> [f32; 3] {
    let (lon, lat) = exogen_core::grid::lon_lat_from_direction(d);
    let x = (lon + std::f64::consts::PI) / (2.0 * std::f64::consts::PI) * img.width as f64 - 0.5;
    let y = (std::f64::consts::FRAC_PI_2 - lat) / std::f64::consts::PI * img.height as f64 - 0.5;
    let x0 = x.floor();
    let y0 = y.floor();
    let fx = (x - x0) as f32;
    let fy = (y - y0) as f32;
    let w = img.width as i64;
    let h = img.height as i64;
    let px = |xi: i64, yi: i64| -> [f32; 3] {
        let xi = xi.rem_euclid(w) as usize;
        let yi = yi.clamp(0, h - 1) as usize;
        img.rgb[yi * img.width + xi]
    };
    let (xi, yi) = (x0 as i64, y0 as i64);
    let a = px(xi, yi);
    let b = px(xi + 1, yi);
    let c = px(xi, yi + 1);
    let dd = px(xi + 1, yi + 1);
    let mut out = [0.0f32; 3];
    for k in 0..3 {
        out[k] =
            (a[k] * (1.0 - fx) + b[k] * fx) * (1.0 - fy) + (c[k] * (1.0 - fx) + dd[k] * fx) * fy;
    }
    out
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GlobeView {
    /// Rotation of the body about its axis, radians.
    pub yaw: f32,
    /// Tilt of the view, radians, positive looks down on the north pole.
    pub pitch: f32,
    pub zoom: f32,
    pub sun_azimuth: f32,
    pub sun_elevation: f32,
    pub exposure: f32,
}

impl Default for GlobeView {
    fn default() -> Self {
        Self {
            yaw: 0.0,
            pitch: 0.3,
            zoom: 1.0,
            sun_azimuth: -0.8,
            sun_elevation: 0.35,
            exposure: 1.0,
        }
    }
}

/// Orthographic software render of the globe into a square image.
pub fn render_globe(img: &EquirectRgb, view: &GlobeView, size: usize) -> egui::ColorImage {
    let mut pixels = vec![egui::Color32::from_rgb(8, 9, 12); size * size];
    let (sy, cy) = view.yaw.sin_cos();
    let (sp, cp) = view.pitch.sin_cos();
    let sun = {
        let (sa, ca) = view.sun_azimuth.sin_cos();
        let (se, ce) = view.sun_elevation.sin_cos();
        [ce * sa, se, ce * ca]
    };
    let scale = 1.0 / view.zoom.max(0.1);
    for py in 0..size {
        for px in 0..size {
            let u = ((px as f32 + 0.5) / size as f32 * 2.0 - 1.0) * scale;
            let v = (1.0 - (py as f32 + 0.5) / size as f32 * 2.0) * scale;
            let r2 = u * u + v * v;
            if r2 > 1.0 {
                continue;
            }
            let w = (1.0 - r2).sqrt();
            // View-space point on the unit sphere: x right, y up, z toward the viewer.
            let n = [u, v, w];
            // Undo pitch (about x) then yaw (about the body's axis, which is view y before pitch).
            let y1 = n[1] * cp + n[2] * sp;
            let z1 = -n[1] * sp + n[2] * cp;
            let x1 = n[0];
            let xb = x1 * cy - z1 * sy;
            let zb = x1 * sy + z1 * cy;
            // Body frame: +Z is the pole, so map view y to body z.
            let d = [zb as f64, xb as f64, y1 as f64];
            let colour = sample(img, d);
            let lambert = (n[0] * sun[0] + n[1] * sun[1] + n[2] * sun[2]).max(0.0);
            let light = 0.06 + 0.94 * lambert;
            let e = view.exposure * light;
            pixels[py * size + px] = egui::Color32::from_rgb(
                srgb8(colour[0] * e),
                srgb8(colour[1] * e),
                srgb8(colour[2] * e),
            );
        }
    }
    egui::ColorImage {
        size: [size, size],
        source_size: egui::vec2(size as f32, size as f32),
        pixels,
    }
}
