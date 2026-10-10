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
pub struct CraterInfo {
    pub lat_deg: f64,
    pub lon_deg: f64,
    pub diameter_km: f64,
    pub age_gyr: f64,
    pub class: String,
    pub impact_angle_deg: f64,
    /// Rim-to-floor relief left at the end of the coarse loop for loop craters, metres.
    pub relief_m: f64,
    /// Whether the crater is still recognisable on the surface.
    pub visible: bool,
    pub from_loop: bool,
}

#[derive(Serialize)]
pub struct FloodInfo {
    pub lat_deg: f64,
    pub lon_deg: f64,
    pub volume_km3: f64,
    pub area_km2: f64,
    pub age_gyr: f64,
}

#[derive(Serialize)]
pub struct RidgeInfo {
    pub lat_deg: f64,
    pub lon_deg: f64,
    pub length_km: f64,
    pub height_m: f64,
    pub age_gyr: f64,
}

#[derive(Serialize)]
pub struct ScarpInfo {
    pub lat_deg: f64,
    pub lon_deg: f64,
    /// Strike at the midpoint, degrees clockwise from north.
    pub strike_deg: f64,
    pub length_km: f64,
    pub height_m: f64,
    pub age_gyr: f64,
}

#[derive(Serialize)]
pub struct UnitInfo {
    pub id: u16,
    pub name: String,
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
    /// Craters emplaced in the loop and the larger replayed ones.
    pub craters: Vec<CraterInfo>,
    pub replayed_crater_count: usize,
    /// Contraction scarps present at the epoch.
    pub scarps: Vec<ScarpInfo>,
    /// Lava floods present at the epoch.
    pub floods: Vec<FloodInfo>,
    /// Wrinkle ridges present at the epoch.
    pub ridges: Vec<RidgeInfo>,
    /// Area fraction of the surface that is lava plains or craters formed in them.
    pub lava_area_fraction: f64,
    pub unit_legend: Vec<UnitInfo>,
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
    // All or nothing: make sure every output can be written before writing any, so that
    // a file held open by a viewer leaves the previous set intact instead of a mixed one.
    for name in [
        "height.exr",
        "height.png",
        "albedo.exr",
        "albedo.png",
        "shaded.png",
        "mask_liquid.png",
        "surface_age.exr",
        "surface_age.png",
        "units.png",
        "sidecar.json",
    ] {
        let path = settings.dir.join(name);
        retry(|| {
            std::fs::OpenOptions::new()
                .write(true)
                .create(true)
                .truncate(false)
                .open(&path)
        })
        .map_err(|source| ExportError::Io {
            path: path.display().to_string(),
            source,
        })?;
    }
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

    // Shaded relief preview: the albedo, exposed to a common mean, lit from one fixed sun
    // over the height map.
    let shade = crate::shade::hillshade(&height_eq, width, height_px, graph.f64("body.radius"));
    let gain = crate::shade::preview_gain(
        (0..width * height_px).map(|k| (channels[0][k] + channels[1][k] + channels[2][k]) / 3.0),
    );
    let lit: Vec<Vec<f32>> = (0..3)
        .map(|ch| {
            (0..width * height_px)
                .map(|k| {
                    crate::shade::lit(
                        [channels[0][k], channels[1][k], channels[2][k]],
                        shade[k],
                        gain,
                    )[ch]
                })
                .collect()
        })
        .collect();
    write_png_srgb(&settings.dir.join("shaded.png"), width, height_px, &lit)?;
    maps.insert(
        "shaded_png".to_string(),
        MapInfo {
            file: "shaded.png".into(),
            encoding: "8-bit sRGB PNG".into(),
            units: "preview: albedo exposed to a common mean, lit from the upper left".into(),
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
    write_png16(
        &settings.dir.join("surface_age.png"),
        width,
        height_px,
        &age_eq,
        amin,
        amax,
    )?;
    maps.insert(
        "surface_age_png".to_string(),
        MapInfo {
            file: "surface_age.png".into(),
            encoding: "16-bit grey PNG, linear from min to max".into(),
            units: "Gyr since formation of the current surface".into(),
            min: amin as f64,
            max: amax as f64,
        },
    );

    // Units, nearest cell per pixel, as an index map: the id times 16.
    let unit_ids = out.layers.unit_id.as_slice();
    let unit_eq: Vec<u8> = (0..width * height_px)
        .map(|k| {
            let (x, y) = (k % width, k / width);
            let lon = ((x as f64 + 0.5) / width as f64 - 0.5) * 2.0 * std::f64::consts::PI;
            let lat = (0.5 - (y as f64 + 0.5) / height_px as f64) * std::f64::consts::PI;
            let d = crate::grid::direction_from_lon_lat(lon, lat);
            (unit_ids[grid.index(grid.locate(d).0)] as u8).saturating_mul(16)
        })
        .collect();
    write_png8(&settings.dir.join("units.png"), width, height_px, &unit_eq)?;
    maps.insert(
        "units_png".to_string(),
        MapInfo {
            file: "units.png".into(),
            encoding: "8-bit grey PNG, unit id times 16".into(),
            units: "unit id, see unit_legend".into(),
            min: 0.0,
            max: 255.0,
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
        craters: out
            .craters
            .iter()
            .map(|c| {
                let (lon, lat) = crate::grid::lon_lat_from_direction(c.center);
                CraterInfo {
                    lat_deg: lat.to_degrees(),
                    lon_deg: lon.to_degrees(),
                    diameter_km: c.diameter / 1e3,
                    age_gyr: out.epoch_time_gyr - c.time_gyr,
                    class: c.class.name().to_string(),
                    impact_angle_deg: c.angle.to_degrees(),
                    relief_m: c.relief as f64,
                    visible: c.visible,
                    from_loop: c.from_loop,
                }
            })
            .collect(),
        replayed_crater_count: out.replayed_count,
        scarps: out
            .scarps
            .iter()
            .map(|s| {
                let (lon, lat) = crate::grid::lon_lat_from_direction(s.midpoint());
                ScarpInfo {
                    lat_deg: lat.to_degrees(),
                    lon_deg: lon.to_degrees(),
                    strike_deg: (90.0 - s.midpoint_heading().to_degrees()).rem_euclid(360.0),
                    length_km: s.length / 1e3,
                    height_m: s.height,
                    age_gyr: out.epoch_time_gyr - s.time_gyr,
                }
            })
            .collect(),
        floods: out
            .floods
            .iter()
            .map(|f| {
                let (lon, lat) = crate::grid::lon_lat_from_direction(f.vent);
                FloodInfo {
                    lat_deg: lat.to_degrees(),
                    lon_deg: lon.to_degrees(),
                    volume_km3: f.volume / 1e9,
                    area_km2: f.area / 1e6,
                    age_gyr: out.epoch_time_gyr - f.time_gyr,
                }
            })
            .collect(),
        ridges: out
            .ridges
            .iter()
            .map(|r| {
                let (lon, lat) = crate::grid::lon_lat_from_direction(r.midpoint());
                RidgeInfo {
                    lat_deg: lat.to_degrees(),
                    lon_deg: lon.to_degrees(),
                    length_km: r.length / 1e3,
                    height_m: r.height + r.crest,
                    age_gyr: out.epoch_time_gyr - r.time_gyr,
                }
            })
            .collect(),
        lava_area_fraction: crate::stages::volcanism::lava_fraction(&out.layers),
        unit_legend: crate::units::legend()
            .into_iter()
            .map(|(id, name)| UnitInfo {
                id,
                name: name.to_string(),
            })
            .collect(),
    };
    let path = settings.dir.join("sidecar.json");
    let text = serde_json::to_string_pretty(&sidecar)?;
    retry(|| std::fs::write(&path, &text)).map_err(|source| ExportError::Io {
        path: path.display().to_string(),
        source,
    })?;
    Ok(sidecar)
}

/// Run a file write, retrying for a couple of seconds: on Windows a viewer or a scanner
/// can hold a mapping on the file being replaced, which fails the write for a moment.
fn retry<T, E>(mut f: impl FnMut() -> Result<T, E>) -> Result<T, E> {
    let mut attempts = 0;
    loop {
        match f() {
            Err(_) if attempts < 10 => {
                attempts += 1;
                std::thread::sleep(std::time::Duration::from_millis(250));
            }
            result => return result,
        }
    }
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
    let image = Image::from_layer(layer);
    retry(|| image.write().non_parallel().to_file(path))
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
    let image = Image::from_layer(layer);
    retry(|| image.write().non_parallel().to_file(path))
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
    retry(|| img.save(path))?;
    Ok(())
}

fn write_png8(path: &Path, width: usize, height: usize, data: &[u8]) -> Result<(), ExportError> {
    let img = image::ImageBuffer::<image::Luma<u8>, Vec<u8>>::from_raw(
        width as u32,
        height as u32,
        data.to_vec(),
    )
    .expect("buffer size");
    retry(|| img.save(path))?;
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
    retry(|| img.save(path))?;
    Ok(())
}
