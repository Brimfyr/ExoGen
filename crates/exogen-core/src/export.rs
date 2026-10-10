// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! Export of equirectangular maps and the sidecar (PLAN.md 0 and 4.13).
//!
//! Height goes out as float32 EXR in metres relative to the datum, plus a 16-bit PNG
//! normalised copy. Albedo goes out as linear float16 EXR and an sRGB-encoded 8-bit PNG.
//! Masks are 8-bit PNG. The sidecar JSON records units, ranges, the datum, the epoch,
//! every parameter with its state, and the plausibility warnings.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use exr::prelude::{
    AnyChannel, AnyChannels, Encoding, FlatSamples, Image, Layer, LayerAttributes, WritableImage,
};
use serde::Serialize;

use crate::grid::projection::to_equirectangular;
use crate::params::{Kind, ParamGraph};
use crate::stages::FinishOutput;

#[derive(Debug, thiserror::Error)]
pub enum ExportError {
    #[error("io error writing {path}: {source}")]
    Io {
        path: String,
        source: std::io::Error,
    },
    #[error("image encoding error: {0}")]
    Image(#[from] image::ImageError),
    #[error("exr encoding error: {0}")]
    Exr(String),
    #[error("json error: {0}")]
    Json(#[from] serde_json::Error),
}

#[derive(Serialize)]
pub struct MapInfo {
    pub file: String,
    pub encoding: String,
    pub units: String,
    pub min: f64,
    pub max: f64,
}

#[derive(Serialize)]
pub struct ParameterInfo {
    pub value: serde_json::Value,
    pub state: String,
    pub unit: String,
}

#[derive(Serialize)]
pub struct Sidecar {
    pub exogen_version: String,
    pub project: String,
    pub width: usize,
    pub height: usize,
    pub radius_m: f64,
    pub datum: String,
    pub datum_offset_m: f64,
    pub epoch_gyr_before_present: f64,
    pub epoch_time_since_formation_gyr: f64,
    pub maps: BTreeMap<String, MapInfo>,
    pub stage_outputs: BTreeMap<String, f64>,
    pub parameters: BTreeMap<String, ParameterInfo>,
    pub warnings: Vec<String>,
}

pub struct ExportSettings {
    pub width: usize,
    pub dir: PathBuf,
    pub project_name: String,
}

/// Write the Phase 1 map set and return the sidecar that describes it.
pub fn export(
    graph: &ParamGraph,
    out: &FinishOutput,
    settings: &ExportSettings,
) -> Result<Sidecar, ExportError> {
    let width = settings.width.max(4);
    let height_px = width / 2;
    std::fs::create_dir_all(&settings.dir).map_err(|source| ExportError::Io {
        path: settings.dir.display().to_string(),
        source,
    })?;
    let mut maps = BTreeMap::new();

    // Datum.
    let datum = graph.str("output.datum").to_string();
    let datum_offset = if datum == "sea_level" && out.sea_level.is_finite() {
        out.sea_level
    } else {
        0.0
    };

    // Height.
    let grid = out.grid;
    let height_eq = to_equirectangular(&grid, out.height.as_slice(), width, height_px);
    let height_rel: Vec<f32> = height_eq.iter().map(|h| h - datum_offset as f32).collect();
    let (hmin, hmax) = min_max(&height_rel);
    write_exr_f32(
        &settings.dir.join("height.exr"),
        width,
        height_px,
        "height",
        &height_rel,
    )?;
    maps.insert(
        "height".to_string(),
        MapInfo {
            file: "height.exr".into(),
            encoding: "float32 EXR, channel height".into(),
            units: "m above datum".into(),
            min: hmin as f64,
            max: hmax as f64,
        },
    );
    write_png16(
        &settings.dir.join("height.png"),
        width,
        height_px,
        &height_rel,
        hmin,
        hmax,
    )?;
    maps.insert(
        "height_png".to_string(),
        MapInfo {
            file: "height.png".into(),
            encoding: "16-bit grey PNG, linear from min to max".into(),
            units: "m above datum".into(),
            min: hmin as f64,
            max: hmax as f64,
        },
    );

    // Albedo.
    let channels: Vec<Vec<f32>> = (0..3)
        .map(|ch| {
            let plane: Vec<f32> = out.albedo.as_slice().iter().map(|c| c[ch]).collect();
            to_equirectangular(&grid, &plane, width, height_px)
        })
        .collect();
    write_exr_rgb(
        &settings.dir.join("albedo.exr"),
        width,
        height_px,
        &channels,
    )?;
    maps.insert(
        "albedo".to_string(),
        MapInfo {
            file: "albedo.exr".into(),
            encoding: "float16 EXR RGB, linear reflectance".into(),
            units: "reflectance".into(),
            min: 0.0,
            max: 1.0,
        },
    );
    write_png_srgb(
        &settings.dir.join("albedo.png"),
        width,
        height_px,
        &channels,
    )?;
    maps.insert(
        "albedo_png".to_string(),
        MapInfo {
            file: "albedo.png".into(),
            encoding: "8-bit sRGB PNG".into(),
            units: "reflectance".into(),
            min: 0.0,
            max: 1.0,
        },
    );

    // Ocean mask.
    let liquid_eq = to_equirectangular(&grid, out.layers.liquid.as_slice(), width, height_px);
    let mask: Vec<u8> = liquid_eq
        .iter()
        .map(|d| if *d > 0.0 { 255 } else { 0 })
        .collect();
    write_png8(
        &settings.dir.join("mask_liquid.png"),
        width,
        height_px,
        &mask,
    )?;
    maps.insert(
        "mask_liquid".to_string(),
        MapInfo {
            file: "mask_liquid.png".into(),
            encoding: "8-bit grey PNG".into(),
            units: "mask".into(),
            min: 0.0,
            max: 1.0,
        },
    );

    // Surface age.
    let age_eq = to_equirectangular(&grid, out.layers.surface_age.as_slice(), width, height_px);
    let (amin, amax) = min_max(&age_eq);
    write_exr_f32(
        &settings.dir.join("surface_age.exr"),
        width,
        height_px,
        "age",
        &age_eq,
    )?;
    maps.insert(
        "surface_age".to_string(),
        MapInfo {
            file: "surface_age.exr".into(),
            encoding: "float32 EXR".into(),
            units: "Gyr since formation of the current surface".into(),
            min: amin as f64,
            max: amax as f64,
        },
    );

    // Sidecar.
    let mut parameters = BTreeMap::new();
    for def in graph.defs() {
        let state = match def.kind {
            Kind::Root => {
                if graph.is_modified(def.id) {
                    "set"
                } else {
                    "default"
                }
            }
            Kind::Derived => {
                if graph.is_pinned(def.id) {
                    "pinned"
                } else {
                    "derived"
                }
            }
            Kind::Multiplier => "multiplier",
            Kind::StageOutput => "stage_output",
        };
        parameters.insert(
            def.id.to_string(),
            ParameterInfo {
                value: serde_json::to_value(graph.get(def.id))?,
                state: state.to_string(),
                unit: def.unit.to_string(),
            },
        );
    }
    let mut stage_outputs = BTreeMap::new();
    stage_outputs.insert("sea_level_m".to_string(), out.sea_level);
    stage_outputs.insert("ocean_fraction".to_string(), out.ocean_fraction);
    let sidecar = Sidecar {
        exogen_version: crate::VERSION.to_string(),
        project: settings.project_name.clone(),
        width,
        height: height_px,
        radius_m: graph.f64("body.radius"),
        datum,
        datum_offset_m: datum_offset,
        epoch_gyr_before_present: graph.f64("output.epoch"),
        epoch_time_since_formation_gyr: out.epoch_time_gyr,
        maps,
        stage_outputs,
        parameters,
        warnings: graph
            .warnings()
            .into_iter()
            .map(|w| format!("{:?}: {}", w.severity, w.text))
            .collect(),
    };
    let path = settings.dir.join("sidecar.json");
    let text = serde_json::to_string_pretty(&sidecar)?;
    std::fs::write(&path, text).map_err(|source| ExportError::Io {
        path: path.display().to_string(),
        source,
    })?;
    Ok(sidecar)
}

fn min_max(v: &[f32]) -> (f32, f32) {
    v.iter()
        .fold((f32::INFINITY, f32::NEG_INFINITY), |(lo, hi), &x| {
            (lo.min(x), hi.max(x))
        })
}

fn write_exr_f32(
    path: &Path,
    width: usize,
    height: usize,
    channel: &str,
    data: &[f32],
) -> Result<(), ExportError> {
    let layer = Layer::new(
        (width, height),
        LayerAttributes::named(channel),
        Encoding::FAST_LOSSLESS,
        AnyChannels::sort(smallvec::smallvec![AnyChannel::new(
            channel,
            FlatSamples::F32(data.to_vec())
        )]),
    );
    Image::from_layer(layer)
        .write()
        .to_file(path)
        .map_err(|e| ExportError::Exr(e.to_string()))
}

fn write_exr_rgb(
    path: &Path,
    width: usize,
    height: usize,
    channels: &[Vec<f32>],
) -> Result<(), ExportError> {
    let to_half =
        |v: &Vec<f32>| FlatSamples::F16(v.iter().map(|x| half::f16::from_f32(*x)).collect());
    let layer = Layer::new(
        (width, height),
        LayerAttributes::named("albedo"),
        Encoding::FAST_LOSSLESS,
        AnyChannels::sort(smallvec::smallvec![
            AnyChannel::new("R", to_half(&channels[0])),
            AnyChannel::new("G", to_half(&channels[1])),
            AnyChannel::new("B", to_half(&channels[2])),
        ]),
    );
    Image::from_layer(layer)
        .write()
        .to_file(path)
        .map_err(|e| ExportError::Exr(e.to_string()))
}

fn write_png16(
    path: &Path,
    width: usize,
    height: usize,
    data: &[f32],
    lo: f32,
    hi: f32,
) -> Result<(), ExportError> {
    let scale = if hi > lo { 65535.0 / (hi - lo) } else { 0.0 };
    let pixels: Vec<u16> = data
        .iter()
        .map(|&v| ((v - lo) * scale).round().clamp(0.0, 65535.0) as u16)
        .collect();
    let img = image::ImageBuffer::<image::Luma<u16>, Vec<u16>>::from_raw(
        width as u32,
        height as u32,
        pixels,
    )
    .expect("buffer size");
    img.save(path)?;
    Ok(())
}

fn write_png8(path: &Path, width: usize, height: usize, data: &[u8]) -> Result<(), ExportError> {
    let img = image::ImageBuffer::<image::Luma<u8>, Vec<u8>>::from_raw(
        width as u32,
        height as u32,
        data.to_vec(),
    )
    .expect("buffer size");
    img.save(path)?;
    Ok(())
}

/// Linear reflectance to 8-bit sRGB.
pub fn srgb_encode(v: f32) -> u8 {
    let v = v.clamp(0.0, 1.0);
    let s = if v <= 0.003_130_8 {
        12.92 * v
    } else {
        1.055 * v.powf(1.0 / 2.4) - 0.055
    };
    (s * 255.0).round() as u8
}

fn write_png_srgb(
    path: &Path,
    width: usize,
    height: usize,
    channels: &[Vec<f32>],
) -> Result<(), ExportError> {
    let mut pixels = Vec::with_capacity(width * height * 3);
    for ((r, g), b) in channels[0].iter().zip(&channels[1]).zip(&channels[2]) {
        pixels.push(srgb_encode(*r));
        pixels.push(srgb_encode(*g));
        pixels.push(srgb_encode(*b));
    }
    let img = image::ImageBuffer::<image::Rgb<u8>, Vec<u8>>::from_raw(
        width as u32,
        height as u32,
        pixels,
    )
    .expect("buffer size");
    img.save(path)?;
    Ok(())
}
