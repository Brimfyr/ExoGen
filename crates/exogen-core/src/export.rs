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

/// A map in the sidecar: its file, encoding and units, its range, how it may be filtered,
/// "linear" or "nearest" where its values are codes or phases that must not be blended
/// (F253), and its channels where it has several or they need saying.
#[derive(Serialize, Default)]
pub struct MapInfo {
    pub file: String,
    pub encoding: String,
    pub units: String,
    pub min: f64,
    pub max: f64,
    pub filter: String,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub channels: Vec<ChannelInfo>,
}

/// A channel of a map (F253): its name in the file, what it holds, its unit and its range.
#[derive(Serialize)]
pub struct ChannelInfo {
    pub name: String,
    pub holds: String,
    pub unit: String,
    pub min: f64,
    pub max: f64,
}

/// The ices that come and go through the year, as the cryosphere's maps draw them (F253).
#[derive(Serialize)]
pub struct CryosphereInfo {
    pub phase_origin: String,
    pub evaluation: String,
    pub snapshot_phase: Option<f64>,
    pub thickest_m: f64,
    pub ices: Vec<IceInfo>,
    pub clock: Option<ClockInfo>,
    pub year: Vec<YearRow>,
}

/// The year's water (F268), where the body has a working liquid under an atmosphere: its
/// liquid's vapour; kg/m² a year over the body, precipitation and evaporation; over its land,
/// precipitation, evaporation and runoff; and the shares of the land arid and hyper-arid, its
/// precipitation under 0.2 and 0.05 of its potential evaporation (UNEP 1992).
#[derive(Serialize)]
pub struct WaterInfo {
    pub vapour: String,
    pub precipitation_kg_m2_yr: f64,
    pub evaporation_kg_m2_yr: f64,
    pub land_precipitation_kg_m2_yr: f64,
    pub land_evaporation_kg_m2_yr: f64,
    pub land_runoff_kg_m2_yr: f64,
    pub arid_land_fraction: f64,
    pub hyper_arid_land_fraction: f64,
}

/// The year's water clouds (M66), a diagnostic the climate does not read until their means come
/// from the clouds themselves (M67): over the body by area, the share of it they cover and the
/// low cloud's, the layers', the anvils' (with the storms' high cloud at their level, M70) and
/// the towers', which overlap at random; the albedo
/// they add over the ground; and what they change of what it lets out to space and of what its
/// air sends down to its ground, W/m², each with them less without.
#[derive(Serialize)]
pub struct CloudsInfo {
    pub cover: f64,
    pub low: f64,
    pub layers: f64,
    pub anvils: f64,
    pub towers: f64,
    pub albedo_added: f64,
    pub outgoing_change_w_m2: f64,
    pub downward_change_w_m2: f64,
}

/// The rivers and the lakes of the year's runoff (F277): the runoff of the land about the lakes,
/// what the lakes' surfaces gain over what they lose to the air, and what reaches the sea,
/// m³/s; the share of the land whose water ends in a lake that does not spill to the sea
/// and the share the lakes cover; and the most any one river brings the sea, m³/s, and where,
/// deg; and what the rivers' erosion moved through today's climate's span (F278): the ground
/// they cut, what of it they laid in lakes and pits and what they brought the sea, km³, and the
/// rate they began at over the land, mm a year.
#[derive(Serialize)]
pub struct RiversInfo {
    pub runoff_m3_s: f64,
    pub lakes_gain_m3_s: f64,
    pub to_sea_m3_s: f64,
    pub closed_land_fraction: f64,
    pub lake_land_fraction: f64,
    pub largest_river_m3_s: f64,
    pub largest_mouth_lat_lon_deg: [f64; 2],
    pub eroded_km3: f64,
    pub laid_km3: f64,
    pub brought_to_sea_km3: f64,
    pub erosion_mm_yr: f64,
}

/// The year's surface wind (F267): its mean speed over the body, m/s, the fastest of its
/// places' year-mean speeds, m/s, the work the atmosphere's heat engine does, W/m² over the
/// body, the thermal equator's southernmost and northernmost through the year, deg, and where
/// the trades end north and south in the year's mean, deg.
#[derive(Serialize)]
pub struct WindInfo {
    pub mean_speed_m_s: f64,
    pub fastest_m_s: f64,
    pub work_w_m2: f64,
    pub thermal_equator_deg: [f64; 2],
    pub trades_end_deg: [f64; 2],
}

/// An ice of the cryosphere's maps (F253): its index in cryosphere_ice.png, its name and its
/// species, where it lies, its colour, linear reflectance, its density, the mass that hides
/// 1 - 1/e of what lies under it, the mass it lies at, past which half is hidden, which
/// opens and closes its window, and its year curve's exponents, none where no place's ice
/// comes and goes thickly enough to fit them.
#[derive(Serialize)]
pub struct IceInfo {
    pub index: usize,
    pub name: String,
    pub species: String,
    pub on: String,
    pub colour_linear: [f32; 3],
    pub density_kg_m3: f64,
    pub cover_mass_kg_m2: f64,
    pub lying_mass_kg_m2: f64,
    pub curve: Option<CurveInfo>,
}

#[derive(Serialize)]
pub struct CurveInfo {
    pub rise: f64,
    pub fall: f64,
}

/// The clock the phases run on (F253).
#[derive(Serialize)]
pub struct ClockInfo {
    pub cycle_s: f64,
    pub orbits: u32,
    pub steps: usize,
    pub solar_day_s: Option<f64>,
    pub obliquity_deg: f64,
    pub eccentricity: f64,
    pub perihelion_ls_deg: f64,
    pub northern_spring_phase: f64,
    pub northern_summer_phase: f64,
    pub northern_autumn_phase: f64,
    pub northern_winter_phase: f64,
}

/// A step of the year (F253): its end's phase, the star's longitude seen from the body, its
/// declination and flux, the surface pressure, and each ice's mean over each hemisphere,
/// kg/m² of the hemisphere.
#[derive(Serialize)]
pub struct YearRow {
    pub phase: f64,
    pub ls_deg: f64,
    pub declination_deg: f64,
    pub flux_w_m2: f64,
    pub pressure_pa: f64,
    pub frost_north_kg_m2: f64,
    pub frost_south_kg_m2: f64,
    pub sea_ice_north_kg_m2: f64,
    pub sea_ice_south_kg_m2: f64,
    pub snow_north_kg_m2: f64,
    pub snow_south_kg_m2: f64,
}

/// A parameter in the sidecar: a number in its display `unit`, its stored SI value `value`
/// times `scale`; any other value as it is.
#[derive(Serialize)]
pub struct ParameterInfo {
    pub value: serde_json::Value,
    pub state: String,
    pub unit: String,
    pub scale: f64,
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
pub struct ProvinceInfo {
    pub lat_deg: f64,
    pub lon_deg: f64,
    pub radius_deg: f64,
}

#[derive(Serialize)]
pub struct RiseInfo {
    pub lat_deg: f64,
    pub lon_deg: f64,
    pub diameter_km: f64,
    pub height_m: f64,
    /// Its latitude when it began to rise, before the body turned under it (F104).
    pub paleo_lat_deg: f64,
    /// The azimuth of the axis it grows along, degrees east of north, 0 to 180 (F114).
    pub axis_azimuth_deg: f64,
}

#[derive(Serialize)]
pub struct ShieldInfo {
    pub lat_deg: f64,
    pub lon_deg: f64,
    pub diameter_km: f64,
    pub height_m: f64,
    pub calderas: usize,
    pub age_gyr: f64,
}

#[derive(Serialize)]
pub struct TroughInfo {
    pub lat_deg: f64,
    pub lon_deg: f64,
    pub length_km: f64,
    pub width_km: f64,
    pub depth_m: f64,
    pub age_gyr: f64,
}

#[derive(Serialize)]
pub struct CoronaInfo {
    pub lat_deg: f64,
    pub lon_deg: f64,
    pub diameter_km: f64,
    pub centre_relief_m: f64,
    pub rim_height_m: f64,
    pub age_gyr: f64,
}

#[derive(Serialize)]
pub struct DomeInfo {
    pub lat_deg: f64,
    pub lon_deg: f64,
    pub diameter_km: f64,
    pub height_m: f64,
    pub age_gyr: f64,
}

#[derive(Serialize)]
pub struct PateraInfo {
    pub lat_deg: f64,
    pub lon_deg: f64,
    pub diameter_km: f64,
    pub depth_m: f64,
    pub active: bool,
}

#[derive(Serialize)]
pub struct LineaInfo {
    pub lat_deg: f64,
    pub lon_deg: f64,
    pub length_km: f64,
    pub width_km: f64,
    pub kind: String,
    pub age_gyr: f64,
}

#[derive(Serialize)]
pub struct MountainInfo {
    pub lat_deg: f64,
    pub lon_deg: f64,
    pub length_km: f64,
    pub width_km: f64,
    pub height_m: f64,
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

/// How many floods the sidecar lists; the rest are counted.
const SIDECAR_FLOOD_LIMIT: usize = 4_000;

/// The sidecar's lists of features, which go apart to features.json (F253).
const FEATURE_LISTS: [&str; 14] = [
    "craters",
    "scarps",
    "floods",
    "magma_provinces",
    "ridges",
    "rises",
    "shields",
    "canyons",
    "graben",
    "coronae",
    "domes",
    "paterae",
    "mountains",
    "lineae",
];

fn largest_floods(
    floods: &[crate::stages::volcanism::FloodRecord],
    limit: usize,
) -> Vec<&crate::stages::volcanism::FloodRecord> {
    let mut refs: Vec<&crate::stages::volcanism::FloodRecord> = floods.iter().collect();
    refs.sort_by(|a, b| b.area.total_cmp(&a.area).then(a.id.cmp(&b.id)));
    refs.truncate(limit);
    refs.sort_by(|a, b| a.time_gyr.total_cmp(&b.time_gyr).then(a.id.cmp(&b.id)));
    refs
}

fn trough_info(t: &crate::stages::plume::Trough, epoch: f64) -> TroughInfo {
    let (lon, lat) = crate::grid::lon_lat_from_direction(t.midpoint());
    TroughInfo {
        lat_deg: lat.to_degrees(),
        lon_deg: lon.to_degrees(),
        length_km: t.length / 1e3,
        width_km: t.width / 1e3,
        depth_m: t.depth,
        age_gyr: epoch - t.end_gyr.min(epoch),
    }
}

#[derive(Serialize)]
pub struct UnitInfo {
    pub id: u16,
    pub name: String,
}

/// The working liquid (ADR 0007): what `liquid.exr` holds and how the ground shows through it.
/// Where its depth H is above 0, ground of albedo A shows through it as, in each channel,
/// R = R_deep + (A - R_deep) exp(-2 K H) (Maritorena, Morel and Gentili 1994).
#[derive(Serialize)]
pub struct LiquidInfo {
    pub working_liquid: String,
    /// The sea's surface above the datum, m, none where only lakes lie.
    pub sea_level_m: Option<f64>,
    pub deep_reflectance_linear: [f32; 3],
    pub attenuation_per_m: [f32; 3],
    /// The wavelengths the channels' attenuation is read at, nm.
    pub channel_wavelengths_nm: [f64; 3],
    pub evaluation: String,
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
    /// Vertical exaggeration of the slopes in the shaded preview.
    pub preview_exaggeration: f64,
    pub epoch_gyr_before_present: f64,
    pub epoch_time_since_formation_gyr: f64,
    pub maps: BTreeMap<String, MapInfo>,
    /// The ices that come and go, how the cryosphere's maps hold them (F253).
    pub cryosphere: Option<CryosphereInfo>,
    /// The year's surface wind (F267), none without an atmosphere.
    pub wind: Option<WindInfo>,
    /// The year's water (F268), none without a working liquid under an atmosphere.
    pub water: Option<WaterInfo>,
    /// The year's water clouds (M66), none where its working liquid is not water.
    pub clouds: Option<CloudsInfo>,
    /// The rivers and the lakes (F277), where the body has water.
    pub rivers: Option<RiversInfo>,
    /// The working liquid's depth and optics (ADR 0007), none where no liquid lies.
    pub liquid: Option<LiquidInfo>,
    pub stage_outputs: BTreeMap<String, f64>,
    pub parameters: BTreeMap<String, ParameterInfo>,
    pub warnings: Vec<String>,
    /// Craters emplaced in the loop and the larger replayed ones.
    pub craters: Vec<CraterInfo>,
    pub replayed_crater_count: usize,
    /// Contraction scarps present at the epoch.
    pub scarps: Vec<ScarpInfo>,
    /// Lava floods present at the epoch: how many, and the largest by area; and the
    /// magmatic provinces the vents favoured.
    pub flood_count: usize,
    pub floods: Vec<FloodInfo>,
    pub magma_provinces: Vec<ProvinceInfo>,
    /// Wrinkle ridges present at the epoch.
    pub ridges: Vec<RidgeInfo>,
    /// Plume rises, and the shields, canyons and graben present at the epoch.
    pub rises: Vec<RiseInfo>,
    pub shields: Vec<ShieldInfo>,
    pub canyons: Vec<TroughInfo>,
    pub graben: Vec<TroughInfo>,
    /// An episodic lid: how long before the epoch the surface was renewed, the tessera
    /// share, and the coronae and domes present.
    pub resurfacing_age_gyr: Option<f64>,
    pub tessera_fraction: f64,
    pub coronae: Vec<CoronaInfo>,
    pub domes: Vec<DomeInfo>,
    /// The Io regime: paterae, mountains, how many flow fields, how many plumes erupt at
    /// the epoch and how many eruptions the surface remembers (F84).
    pub paterae: Vec<PateraInfo>,
    pub mountains: Vec<MountainInfo>,
    pub flow_field_count: usize,
    pub plume_count: usize,
    pub plume_episode_count: usize,
    pub small_eruption_count: usize,
    /// An ice shell's cracks present at the epoch: how many, and the longest.
    pub linea_count: usize,
    pub lineae: Vec<LineaInfo>,
    /// An ice shell's lanes of grooved terrain, and the shares of the surface that are
    /// grooved terrain and renewed polar terrain.
    pub lane_count: usize,
    pub lane_coverage: f64,
    pub grooved_fraction: f64,
    /// Chaos generations and the share of the surface that is chaos matrix; cryovolcanic
    /// domes built over the body's life.
    pub chaos_count: usize,
    pub chaos_fraction: f64,
    pub cryodome_count: usize,
    pub polar_terrain_fraction: f64,
    /// Past episodes of focused heating whose terrains the surface carries, and the share
    /// of the surface the lanes and the polar terrains, present and past, renewed.
    pub fossil_terrain_count: usize,
    pub renewed_fraction: f64,
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
    let species: Vec<String> = out
        .caps
        .iter()
        .map(|cap| cap.ice.species.to_lowercase())
        .collect();
    let mut names: Vec<String> = [
        "height.exr",
        "height.png",
        "albedo.exr",
        "albedo.png",
        "liquid.exr",
        "liquid.png",
        "snapshot.exr",
        "snapshot.png",
        "shaded.png",
        "mask_liquid.png",
        "surface_age.exr",
        "surface_age.png",
        "temperature.exr",
        "temperature.png",
        "temperature_min.exr",
        "temperature_min.png",
        "temperature_max.exr",
        "temperature_max.png",
        "frost.exr",
        "frost.png",
        "ice.exr",
        "ice.png",
        "sea_ice.exr",
        "sea_ice.png",
        "cryosphere_window.exr",
        "cryosphere_window.png",
        "cryosphere_thickness.exr",
        "cryosphere_thickness.png",
        "cryosphere_ice.png",
        "wind.exr",
        "wind.png",
        "clouds.exr",
        "clouds.png",
        "precipitation.exr",
        "precipitation.png",
        "aridity.exr",
        "aridity.png",
        "units.png",
        "sidecar.json",
        "features.json",
    ]
    .iter()
    .map(|n| n.to_string())
    .collect();
    for sp in &species {
        names.push(format!("ice_{sp}.exr"));
        names.push(format!("ice_{sp}.png"));
    }
    for name in &names {
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
            ..Default::default()
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
            ..Default::default()
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
            units: "reflectance of the ground itself, under its liquid too (ADR 0007)".into(),
            min: 0.0,
            max: 1.0,
            ..Default::default()
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
            units: "reflectance of the ground itself, under its liquid too (ADR 0007)".into(),
            min: 0.0,
            max: 1.0,
            ..Default::default()
        },
    );

    // The liquid (ADR 0007): its depth over the ground, m, seas and lakes alike, its surface
    // flat; an engine draws the ground through it by the optics the sidecar gives, at the sea
    // level it sets.
    let liquid_eq = to_equirectangular(&grid, out.layers.liquid.as_slice(), width, height_px);
    let depth: Vec<f32> = liquid_eq.iter().map(|d| d.max(0.0)).collect();
    let deepest = depth.iter().fold(0.0f32, |a, &d| a.max(d));
    let depth_units =
        "m of the working liquid over the ground, its surface flat; 0 where none lies";
    write_exr_f32(
        &settings.dir.join("liquid.exr"),
        width,
        height_px,
        "depth",
        &depth,
    )?;
    maps.insert(
        "liquid".to_string(),
        MapInfo {
            file: "liquid.exr".into(),
            encoding: "float32 EXR, channel depth".into(),
            units: depth_units.into(),
            min: 0.0,
            max: deepest as f64,
            ..Default::default()
        },
    );
    write_png16(
        &settings.dir.join("liquid.png"),
        width,
        height_px,
        &depth,
        0.0,
        deepest.max(1e-6),
    )?;
    maps.insert(
        "liquid_png".to_string(),
        MapInfo {
            file: "liquid.png".into(),
            encoding: "16-bit grey PNG, linear from min to max".into(),
            units: depth_units.into(),
            min: 0.0,
            max: deepest as f64,
            ..Default::default()
        },
    );
    let working = graph.str("volatiles.working_liquid");
    let optics = crate::stages::materials::LiquidOptics::of(working);
    let liquid = (deepest > 0.0).then(|| LiquidInfo {
        working_liquid: working.to_string(),
        sea_level_m: out.sea_level.is_finite().then_some(out.sea_level - datum_offset),
        deep_reflectance_linear: optics.deep,
        attenuation_per_m: optics.attenuation,
        channel_wavelengths_nm: [612.0, 549.0, 464.0],
        evaluation: "Where liquid.exr's depth H is above 0, the liquid's surface stands H over the ground, and the ground, of albedo A in albedo.exr, shows through it as R = R_deep + (A - R_deep) exp(-2 K H) in each channel, R_deep its deep reflectance and K its attenuation (Maritorena, Morel and Gentili 1994). At another sea level H is that level less height.exr's ground. snapshot.exr is this at ExoGen's own level, with the ices at the snapshot's phase over it.".into(),
    });

    // The snapshot (ADR 0007): the surface as it looks at the snapshot's season, the albedo
    // under its liquid and the ices over it at the snapshot's phase as the cryosphere's maps
    // draw them, composed on the grid and read between its cells as the albedo is.
    let composed = crate::snapshot::snapshot(out, working);
    let snap: Vec<Vec<f32>> = (0..3)
        .map(|ch| {
            let plane: Vec<f32> = composed.as_slice().iter().map(|c| c[ch]).collect();
            to_equirectangular(&grid, &plane, width, height_px)
        })
        .collect();
    let snapshot_units = "reflectance of the surface at the snapshot's season: the albedo under its liquid, the ices over it at the snapshot's phase (ADR 0007)";
    write_exr_rgb(&settings.dir.join("snapshot.exr"), width, height_px, &snap)?;
    maps.insert(
        "snapshot".to_string(),
        MapInfo {
            file: "snapshot.exr".into(),
            encoding: "float16 EXR RGB, linear reflectance".into(),
            units: snapshot_units.into(),
            min: 0.0,
            max: 1.0,
            ..Default::default()
        },
    );
    write_png_srgb(&settings.dir.join("snapshot.png"), width, height_px, &snap)?;
    maps.insert(
        "snapshot_png".to_string(),
        MapInfo {
            file: "snapshot.png".into(),
            encoding: "8-bit sRGB PNG".into(),
            units: snapshot_units.into(),
            min: 0.0,
            max: 1.0,
            ..Default::default()
        },
    );

    // Shaded relief preview: the snapshot, exposed to a common mean, lit from one fixed sun
    // over its surface, a liquid's flat (ADR 0007), the slopes exaggerated where the ground is
    // gentle.
    let exaggeration =
        crate::shade::preview_exaggeration(&height_eq, width, height_px, graph.f64("body.radius"));
    let surface_eq = to_equirectangular(
        &grid,
        crate::snapshot::surface(out).as_slice(),
        width,
        height_px,
    );
    let shade = crate::shade::hillshade(
        &surface_eq,
        width,
        height_px,
        graph.f64("body.radius"),
        exaggeration,
    );
    let gain = crate::shade::preview_gain(
        (0..width * height_px).map(|k| (snap[0][k] + snap[1][k] + snap[2][k]) / 3.0),
    );
    let lit: Vec<Vec<f32>> = (0..3)
        .map(|ch| {
            (0..width * height_px)
                .map(|k| {
                    crate::shade::lit([snap[0][k], snap[1][k], snap[2][k]], shade[k], gain)[ch]
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
            units: "preview: the snapshot exposed to a common mean, lit from the upper left over its surface, a liquid's flat".into(),
            min: 0.0,
            max: 1.0,
            ..Default::default()
        },
    );

    // Ocean mask.
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
            ..Default::default()
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
            ..Default::default()
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
            ..Default::default()
        },
    );

    // The year's mean surface temperature (F235).
    let temperature_eq = to_equirectangular(&grid, out.temperature.as_slice(), width, height_px);
    let (tmin, tmax) = min_max(&temperature_eq);
    write_exr_f32(
        &settings.dir.join("temperature.exr"),
        width,
        height_px,
        "temperature",
        &temperature_eq,
    )?;
    maps.insert(
        "temperature".to_string(),
        MapInfo {
            file: "temperature.exr".into(),
            encoding: "float32 EXR".into(),
            units: "K, the year's mean surface temperature through its days and seasons".into(),
            min: tmin as f64,
            max: tmax as f64,
            ..Default::default()
        },
    );
    write_png16(
        &settings.dir.join("temperature.png"),
        width,
        height_px,
        &temperature_eq,
        tmin,
        tmax,
    )?;
    maps.insert(
        "temperature_png".to_string(),
        MapInfo {
            file: "temperature.png".into(),
            encoding: "16-bit grey PNG, linear from min to max".into(),
            units: "K, the year's mean surface temperature through its days and seasons".into(),
            min: tmin as f64,
            max: tmax as f64,
            ..Default::default()
        },
    );

    // The year's coldest and warmest surface temperature at each cell (F253).
    for (key, field, what) in [
        (
            "temperature_min",
            &out.temperature_min,
            "K, the coldest the surface stands through its days and seasons",
        ),
        (
            "temperature_max",
            &out.temperature_max,
            "K, the warmest the surface stands through its days and seasons",
        ),
    ] {
        let eq = to_equirectangular(&grid, field.as_slice(), width, height_px);
        let (lo, hi) = min_max(&eq);
        write_exr_f32(
            &settings.dir.join(format!("{key}.exr")),
            width,
            height_px,
            "temperature",
            &eq,
        )?;
        write_png16(
            &settings.dir.join(format!("{key}.png")),
            width,
            height_px,
            &eq,
            lo,
            hi,
        )?;
        for (suffix, encoding) in [
            ("exr", "float32 EXR"),
            ("png", "16-bit grey PNG, linear from min to max"),
        ] {
            let name = if suffix == "exr" {
                key.to_string()
            } else {
                format!("{key}_png")
            };
            maps.insert(
                name,
                MapInfo {
                    file: format!("{key}.{suffix}"),
                    encoding: encoding.into(),
                    units: what.into(),
                    min: lo as f64,
                    max: hi as f64,
                    ..Default::default()
                },
            );
        }
    }

    // The frost on the ground at the snapshot's season (F239).
    let frost_eq = to_equirectangular(&grid, out.frost.as_slice(), width, height_px);
    let (fmin, fmax) = min_max(&frost_eq);
    let frost_units = match out.frost_gas {
        Some(gas) => format!("kg/m² of {gas} frost on the ground at the snapshot's season"),
        None => "kg/m² of frost on the ground at the snapshot's season: none".to_string(),
    };
    write_exr_f32(
        &settings.dir.join("frost.exr"),
        width,
        height_px,
        "frost",
        &frost_eq,
    )?;
    maps.insert(
        "frost".to_string(),
        MapInfo {
            file: "frost.exr".into(),
            encoding: "float32 EXR".into(),
            units: frost_units.clone(),
            min: fmin as f64,
            max: fmax as f64,
            ..Default::default()
        },
    );
    write_png16(
        &settings.dir.join("frost.png"),
        width,
        height_px,
        &frost_eq,
        fmin,
        fmax,
    )?;
    maps.insert(
        "frost_png".to_string(),
        MapInfo {
            file: "frost.png".into(),
            encoding: "16-bit grey PNG, linear from min to max".into(),
            units: frost_units,
            min: fmin as f64,
            max: fmax as f64,
            ..Default::default()
        },
    );

    // The caps' ice, its thickness (F241).
    let ice_eq = to_equirectangular(&grid, out.layers.ice.as_slice(), width, height_px);
    let (imin, imax) = min_max(&ice_eq);
    write_exr_f32(
        &settings.dir.join("ice.exr"),
        width,
        height_px,
        "ice",
        &ice_eq,
    )?;
    maps.insert(
        "ice".to_string(),
        MapInfo {
            file: "ice.exr".into(),
            encoding: "float32 EXR".into(),
            units: "m of the caps' ice the inventory laid in its cold traps".into(),
            min: imin as f64,
            max: imax as f64,
            ..Default::default()
        },
    );
    write_png16(
        &settings.dir.join("ice.png"),
        width,
        height_px,
        &ice_eq,
        imin,
        imax,
    )?;
    maps.insert(
        "ice_png".to_string(),
        MapInfo {
            file: "ice.png".into(),
            encoding: "16-bit grey PNG, linear from min to max".into(),
            units: "m of the caps' ice the inventory laid in its cold traps".into(),
            min: imin as f64,
            max: imax as f64,
            ..Default::default()
        },
    );

    // Each ice's caps apart (F253): water's, carbon dioxide's and the rest.
    for (cap, sp) in out.caps.iter().zip(&species) {
        let eq = to_equirectangular(&grid, cap.thickness.as_slice(), width, height_px);
        let (lo, hi) = min_max(&eq);
        let what = format!("m of the caps' {} ice", cap.ice.species);
        write_exr_f32(
            &settings.dir.join(format!("ice_{sp}.exr")),
            width,
            height_px,
            "ice",
            &eq,
        )?;
        write_png16(
            &settings.dir.join(format!("ice_{sp}.png")),
            width,
            height_px,
            &eq,
            lo,
            hi,
        )?;
        maps.insert(
            format!("ice_{sp}"),
            MapInfo {
                file: format!("ice_{sp}.exr"),
                encoding: "float32 EXR".into(),
                units: what.clone(),
                min: lo as f64,
                max: hi as f64,
                ..Default::default()
            },
        );
        maps.insert(
            format!("ice_{sp}_png"),
            MapInfo {
                file: format!("ice_{sp}.png"),
                encoding: "16-bit grey PNG, linear from min to max".into(),
                units: what,
                min: lo as f64,
                max: hi as f64,
                ..Default::default()
            },
        );
    }

    // The sea's ice at the snapshot's season (F244).
    let sea_eq = to_equirectangular(&grid, out.sea_ice.as_slice(), width, height_px);
    let (smin, smax) = min_max(&sea_eq);
    let sea_units = "kg/m² of the sea's ice at the snapshot's season";
    write_exr_f32(
        &settings.dir.join("sea_ice.exr"),
        width,
        height_px,
        "sea_ice",
        &sea_eq,
    )?;
    maps.insert(
        "sea_ice".to_string(),
        MapInfo {
            file: "sea_ice.exr".into(),
            encoding: "float32 EXR".into(),
            units: sea_units.into(),
            min: smin as f64,
            max: smax as f64,
            ..Default::default()
        },
    );
    write_png16(
        &settings.dir.join("sea_ice.png"),
        width,
        height_px,
        &sea_eq,
        smin,
        smax,
    )?;
    maps.insert(
        "sea_ice_png".to_string(),
        MapInfo {
            file: "sea_ice.png".into(),
            encoding: "16-bit grey PNG, linear from min to max".into(),
            units: sea_units.into(),
            min: smin as f64,
            max: smax as f64,
            ..Default::default()
        },
    );

    // The ices that come and go through the year (F246, F253), three maps of one meaning
    // each and no alpha. The window: R the phase of the year an ice freezes at, G the phase it
    // melts out at, from the northern spring equinox, 0 and 1 where it lies all year, B the
    // phase of its peak; nearest cell per pixel, so that a window round the year's end is not
    // blended through its middle. Its thickness at its peak, bilinear, but only where the
    // nearest cell's window lies, so that it does not bleed past the window at the ice's
    // margins (F248). Which ice lies, nearest cell, its index among the sidecar's ices.
    let nearest = |field: &crate::fields::Field<f32>| {
        to_equirectangular_nearest(&grid, field.as_slice(), width, height_px)
    };
    let window: Vec<Vec<f32>> = [0, 1, 3]
        .iter()
        .map(|&i| nearest(&out.cryosphere[i]))
        .collect();
    let mut thickness = nearest(&out.cryosphere[2]);
    let smooth = to_equirectangular(&grid, out.cryosphere[2].as_slice(), width, height_px);
    for (v, s) in thickness.iter_mut().zip(smooth) {
        if *v > 0.0 {
            *v = s.max(0.0);
        }
    }
    let thickest = thickness.iter().fold(0.0f32, |a, &v| a.max(v));
    let ice_index: Vec<u8> = nearest(&out.cryosphere_ice)
        .iter()
        .map(|&v| v.round().clamp(0.0, 255.0) as u8)
        .collect();
    let phase_unit = "share of the year's cycle from the northern spring equinox";
    let window_channels = || {
        [
            (
                "R",
                "the phase the ice freezes at, 0 where it lies all year",
            ),
            ("G", "the phase it melts out at, 1 where it lies all year"),
            ("B", "the phase of its peak"),
        ]
        .map(|(name, holds)| ChannelInfo {
            name: name.into(),
            holds: holds.into(),
            unit: phase_unit.into(),
            min: 0.0,
            max: 1.0,
        })
        .into_iter()
        .collect::<Vec<_>>()
    };
    let window_units =
        "the ice's window through the year: it lies from R round to G; all zero where no ice lies";
    write_exr_channels(
        &settings.dir.join("cryosphere_window.exr"),
        (width, height_px),
        "cryosphere_window",
        &["R", "G", "B"],
        &window,
    )?;
    maps.insert(
        "cryosphere_window".to_string(),
        MapInfo {
            file: "cryosphere_window.exr".into(),
            encoding: "float32 EXR, channels R, G, B".into(),
            units: window_units.into(),
            min: 0.0,
            max: 1.0,
            filter: "nearest".into(),
            channels: window_channels(),
        },
    );
    write_png_rgb16(
        &settings.dir.join("cryosphere_window.png"),
        (width, height_px),
        &window,
    )?;
    maps.insert(
        "cryosphere_window_png".to_string(),
        MapInfo {
            file: "cryosphere_window.png".into(),
            encoding: "16-bit RGB PNG, each channel from 0 to 1".into(),
            units: window_units.into(),
            min: 0.0,
            max: 1.0,
            filter: "nearest".into(),
            channels: window_channels(),
        },
    );
    let thickness_units = "m of the ice at its peak, over the height map, which holds none";
    write_exr_f32(
        &settings.dir.join("cryosphere_thickness.exr"),
        width,
        height_px,
        "thickness",
        &thickness,
    )?;
    maps.insert(
        "cryosphere_thickness".to_string(),
        MapInfo {
            file: "cryosphere_thickness.exr".into(),
            encoding: "float32 EXR, channel thickness".into(),
            units: thickness_units.into(),
            min: 0.0,
            max: thickest as f64,
            ..Default::default()
        },
    );
    write_png16(
        &settings.dir.join("cryosphere_thickness.png"),
        width,
        height_px,
        &thickness,
        0.0,
        thickest,
    )?;
    maps.insert(
        "cryosphere_thickness_png".to_string(),
        MapInfo {
            file: "cryosphere_thickness.png".into(),
            encoding: "16-bit grey PNG, linear from min to max".into(),
            units: thickness_units.into(),
            min: 0.0,
            max: thickest as f64,
            ..Default::default()
        },
    );
    write_png8(
        &settings.dir.join("cryosphere_ice.png"),
        width,
        height_px,
        &ice_index,
    )?;
    maps.insert(
        "cryosphere_ice_png".to_string(),
        MapInfo {
            file: "cryosphere_ice.png".into(),
            encoding: "8-bit grey PNG".into(),
            units: "the ice that lies, its index among the sidecar's cryosphere ices, 0 none"
                .into(),
            min: 0.0,
            max: out.cryosphere_ices.len() as f64,
            filter: "nearest".into(),
            ..Default::default()
        },
    );
    // What an engine draws them by: each ice, the clock and the year step by step.
    let cryosphere = (!out.cryosphere_ices.is_empty()).then(|| CryosphereInfo {
        phase_origin: "the northern spring equinox: phases are shares of the cycle's time after it".into(),
        evaluation: "Evaluate each texel at the map's resolution, then filter the result: the window and the ice index are nearest-cell and must not be filtered. An ice lies at phase p where p falls in its window, from R round to G; its mass, kg/m², is its thickness at its peak times its density, times (t / rise)^a from R to its peak's phase B, t the time since R and rise the time from R to B, then times 1 - (g / fall)^b from B to G, g the time since B and fall the time from B to G, (a, b) its curve; all year its peak where R is 0 and G 1. It hides 1 - exp(-m / cover_mass) of what lies under it and stands m / density thick.".into(),
        snapshot_phase: out.cryosphere_snapshot.is_finite().then_some(out.cryosphere_snapshot),
        thickest_m: thickest as f64,
        ices: out
            .cryosphere_ices
            .iter()
            .enumerate()
            .map(|(i, ice)| IceInfo {
                index: i + 1,
                name: ice.name.clone(),
                species: ice.species.to_string(),
                on: if ice.sea { "sea" } else { "land" }.into(),
                colour_linear: ice.colour,
                density_kg_m3: ice.density,
                cover_mass_kg_m2: ice.cover_mass,
                lying_mass_kg_m2: ice.cover_mass * std::f64::consts::LN_2,
                curve: ice.curve.map(|(rise, fall)| CurveInfo { rise, fall }),
            })
            .collect(),
        clock: out.cryosphere_clock.map(|c| ClockInfo {
            cycle_s: c.cycle_s,
            orbits: c.orbits,
            steps: c.steps,
            solar_day_s: c.solar_day_s,
            obliquity_deg: c.obliquity_deg,
            eccentricity: c.eccentricity,
            perihelion_ls_deg: c.perihelion_ls_deg,
            northern_spring_phase: c.seasons[0],
            northern_summer_phase: c.seasons[1],
            northern_autumn_phase: c.seasons[2],
            northern_winter_phase: c.seasons[3],
        }),
        year: out
            .cryosphere_series
            .as_ref()
            .map(|y| {
                (0..y.phase.len())
                    .map(|k| YearRow {
                        phase: y.phase[k],
                        ls_deg: y.ls_deg[k],
                        declination_deg: y.declination_deg[k],
                        flux_w_m2: y.flux[k],
                        pressure_pa: y.pressure[k],
                        frost_north_kg_m2: y.frost[0][k],
                        frost_south_kg_m2: y.frost[1][k],
                        sea_ice_north_kg_m2: y.sea_ice[0][k],
                        sea_ice_south_kg_m2: y.sea_ice[1][k],
                        snow_north_kg_m2: y.snow[0].get(k).copied().unwrap_or(0.0),
                        snow_south_kg_m2: y.snow[1].get(k).copied().unwrap_or(0.0),
                    })
                    .collect()
            })
            .unwrap_or_default(),
    });

    // The year's mean surface wind (F267): R its eastward and G its northward wind, B its
    // mean speed, m/s, bilinear; all zero without an atmosphere. The PNG holds R and G from
    // minus to plus the fastest component, B from 0 to the fastest mean speed.
    let calm = || vec![0.0f32; width * height_px];
    let wind_maps: Vec<Vec<f32>> = match &out.wind {
        Some(fields) => fields
            .iter()
            .map(|f| to_equirectangular(&grid, f.as_slice(), width, height_px))
            .collect(),
        None => vec![calm(), calm(), calm()],
    };
    let strongest = wind_maps[..2]
        .iter()
        .flatten()
        .fold(0.0f32, |a, &v| a.max(v.abs()));
    let fastest = wind_maps[2].iter().fold(0.0f32, |a, &v| a.max(v));
    let wind_channels = |png: bool| {
        let (lo, hi) = if png {
            (-strongest as f64, strongest as f64)
        } else {
            (
                wind_maps[0]
                    .iter()
                    .chain(&wind_maps[1])
                    .fold(0.0f32, |a, &v| a.min(v)) as f64,
                wind_maps[0]
                    .iter()
                    .chain(&wind_maps[1])
                    .fold(0.0f32, |a, &v| a.max(v)) as f64,
            )
        };
        [
            ("R", "the year's mean eastward wind", lo, hi),
            ("G", "the year's mean northward wind", lo, hi),
            ("B", "the year's mean speed", 0.0, fastest as f64),
        ]
        .map(|(name, holds, min, max)| ChannelInfo {
            name: name.into(),
            holds: holds.into(),
            unit: "m/s".into(),
            min,
            max,
        })
        .into_iter()
        .collect::<Vec<_>>()
    };
    let wind_units = "m/s at 10 m above the ground: the year's mean wind, east and north, and its mean speed; all zero without an atmosphere";
    write_exr_channels(
        &settings.dir.join("wind.exr"),
        (width, height_px),
        "wind",
        &["R", "G", "B"],
        &wind_maps,
    )?;
    maps.insert(
        "wind".to_string(),
        MapInfo {
            file: "wind.exr".into(),
            encoding: "float32 EXR, channels R, G, B".into(),
            units: wind_units.into(),
            min: 0.0,
            max: fastest as f64,
            channels: wind_channels(false),
            ..Default::default()
        },
    );
    let span = strongest.max(1e-6);
    let wind_png: Vec<Vec<f32>> = vec![
        wind_maps[0].iter().map(|&v| 0.5 + 0.5 * v / span).collect(),
        wind_maps[1].iter().map(|&v| 0.5 + 0.5 * v / span).collect(),
        wind_maps[2]
            .iter()
            .map(|&v| v / fastest.max(1e-6))
            .collect(),
    ];
    write_png_rgb16(
        &settings.dir.join("wind.png"),
        (width, height_px),
        &wind_png,
    )?;
    maps.insert(
        "wind_png".to_string(),
        MapInfo {
            file: "wind.png".into(),
            encoding:
                "16-bit RGB PNG: R and G linear from the channels' min to max, B from 0 to its max"
                    .into(),
            units: wind_units.into(),
            min: 0.0,
            max: fastest as f64,
            channels: wind_channels(true),
            ..Default::default()
        },
    );
    // The year's water clouds (M66): R the low cloud's year-mean cover, G the layers', B the
    // anvils', A the towers', shares of the place, bilinear; they overlap at random, so that
    // the share any covers is 1 - (1 - R)(1 - G)(1 - B)(1 - A). All zero where the working
    // liquid is not water.
    let cloud_maps: Vec<Vec<f32>> = match &out.clouds {
        Some(fields) => fields[1..]
            .iter()
            .map(|f| to_equirectangular(&grid, f.as_slice(), width, height_px))
            .collect(),
        None => vec![calm(), calm(), calm(), calm()],
    };
    let cloud_channels = || {
        [
            ("R", "the low cloud's year-mean cover"),
            ("G", "the layers' year-mean cover"),
            (
                "B",
                "the anvils' and the storms' high cloud's year-mean cover",
            ),
            ("A", "the towers' year-mean cover"),
        ]
        .map(|(name, holds)| ChannelInfo {
            name: name.into(),
            holds: holds.into(),
            unit: "share of the place".into(),
            min: 0.0,
            max: 1.0,
        })
        .into_iter()
        .collect::<Vec<_>>()
    };
    let cloud_units = "shares of the place each kind of the water's clouds covers through the year, overlapping at random: the share any covers is 1 - (1 - R)(1 - G)(1 - B)(1 - A); a diagnostic the climate does not yet read; all zero where the working liquid is not water";
    write_exr_channels(
        &settings.dir.join("clouds.exr"),
        (width, height_px),
        "clouds",
        &["R", "G", "B", "A"],
        &cloud_maps,
    )?;
    maps.insert(
        "clouds".to_string(),
        MapInfo {
            file: "clouds.exr".into(),
            encoding: "float32 EXR, channels R, G, B, A".into(),
            units: cloud_units.into(),
            min: 0.0,
            max: 1.0,
            channels: cloud_channels(),
            ..Default::default()
        },
    );
    write_png_rgba16(
        &settings.dir.join("clouds.png"),
        (width, height_px),
        &cloud_maps,
    )?;
    maps.insert(
        "clouds_png".to_string(),
        MapInfo {
            file: "clouds.png".into(),
            encoding: "16-bit RGBA PNG, each channel linear from 0 to 1".into(),
            units: cloud_units.into(),
            min: 0.0,
            max: 1.0,
            channels: cloud_channels(),
            ..Default::default()
        },
    );
    // The year's water (F268): its precipitation, kg/m² a year, and its aridity index, its
    // precipitation over its potential evaporation, held under ARIDITY_MAP_MOST, bilinear; all
    // zero without a working liquid under an atmosphere.
    let water_maps: Vec<Vec<f32>> = match &out.water {
        Some(fields) => fields
            .iter()
            .map(|f| to_equirectangular(&grid, f.as_slice(), width, height_px))
            .collect(),
        None => vec![calm(), calm()],
    };
    let wettest = water_maps[0].iter().fold(0.0f32, |a, &v| a.max(v));
    for (k, (name, unit, most)) in [
        (
            "precipitation",
            "kg/m² a year of the working liquid (mm a year of water), the year's mean; zero without a working liquid under an atmosphere",
            wettest,
        ),
        (
            "aridity",
            "the year's precipitation over its potential evaporation, the aridity index (UNEP 1992: under 0.05 hyper-arid, 0.2 arid, 0.5 semi-arid), held under 4; over the sea its open water's; zero without a working liquid under an atmosphere",
            crate::constants::ARIDITY_MAP_MOST as f32,
        ),
    ]
    .into_iter()
    .enumerate()
    {
        write_exr_f32(
            &settings.dir.join(format!("{name}.exr")),
            width,
            height_px,
            name,
            &water_maps[k],
        )?;
        maps.insert(
            name.to_string(),
            MapInfo {
                file: format!("{name}.exr"),
                encoding: format!("float32 EXR, channel {name}"),
                units: unit.into(),
                min: 0.0,
                max: most as f64,
                ..Default::default()
            },
        );
        write_png16(
            &settings.dir.join(format!("{name}.png")),
            width,
            height_px,
            &water_maps[k],
            0.0,
            most.max(1e-6),
        )?;
        maps.insert(
            format!("{name}_png"),
            MapInfo {
                file: format!("{name}.png"),
                encoding: "16-bit grey PNG, linear from min to max".into(),
                units: unit.into(),
                min: 0.0,
                max: most as f64,
                ..Default::default()
            },
        );
    }
    // The rivers (F277): the water each cell's river carries down from it, m³/s, the most of
    // the cells in each pixel so that a river a cell wide keeps its course; and the lakes'
    // depth, m, nearest-cell; both zero without the year's water.
    let (rivers_map, lakes_map) = match &out.rivers {
        Some(r) => (
            to_equirectangular_most(&grid, r.discharge.as_slice(), width, height_px),
            to_equirectangular_nearest(&grid, r.lakes.as_slice(), width, height_px),
        ),
        None => (calm(), calm()),
    };
    let widest = rivers_map.iter().fold(0.0f32, |a, &v| a.max(v));
    let deepest = lakes_map.iter().fold(0.0f32, |a, &v| a.max(v));
    let rivers_unit = "m³/s of the working liquid each cell's river carries down from it, the most of the cells in the pixel; zero at sea, in lakes and without the year's water";
    let lakes_unit = "m of the lakes' water over the ground, nearest-cell; zero where none stands and without the year's water";
    write_exr_f32(
        &settings.dir.join("rivers.exr"),
        width,
        height_px,
        "rivers",
        &rivers_map,
    )?;
    maps.insert(
        "rivers".into(),
        MapInfo {
            file: "rivers.exr".into(),
            encoding: "float32 EXR, channel rivers".into(),
            units: rivers_unit.into(),
            min: 0.0,
            max: widest as f64,
            filter: "nearest".into(),
            ..Default::default()
        },
    );
    let logged: Vec<f32> = rivers_map.iter().map(|&q| (1.0 + q).log10()).collect();
    let top = (1.0 + widest).log10().max(1e-6);
    write_png16(
        &settings.dir.join("rivers.png"),
        width,
        height_px,
        &logged,
        0.0,
        top,
    )?;
    maps.insert(
        "rivers_png".into(),
        MapInfo {
            file: "rivers.png".into(),
            encoding: "16-bit grey PNG, log10(1 + discharge) linear from min to max".into(),
            units: rivers_unit.into(),
            min: 0.0,
            max: top as f64,
            filter: "nearest".into(),
            ..Default::default()
        },
    );
    write_exr_f32(
        &settings.dir.join("lakes.exr"),
        width,
        height_px,
        "lakes",
        &lakes_map,
    )?;
    maps.insert(
        "lakes".into(),
        MapInfo {
            file: "lakes.exr".into(),
            encoding: "float32 EXR, channel lakes".into(),
            units: lakes_unit.into(),
            min: 0.0,
            max: deepest as f64,
            filter: "nearest".into(),
            ..Default::default()
        },
    );
    write_png16(
        &settings.dir.join("lakes.png"),
        width,
        height_px,
        &lakes_map,
        0.0,
        deepest.max(1e-6),
    )?;
    maps.insert(
        "lakes_png".into(),
        MapInfo {
            file: "lakes.png".into(),
            encoding: "16-bit grey PNG, linear from min to max".into(),
            units: lakes_unit.into(),
            min: 0.0,
            max: deepest as f64,
            filter: "nearest".into(),
            ..Default::default()
        },
    );
    let rivers = out.rivers.as_ref().map(|r| {
        let s = &r.summary;
        let (lon, lat) = crate::grid::lon_lat_from_direction(s.mouth);
        RiversInfo {
            runoff_m3_s: s.runoff,
            lakes_gain_m3_s: s.lakes_gain,
            to_sea_m3_s: s.to_sea,
            closed_land_fraction: s.closed,
            lake_land_fraction: s.lakes,
            largest_river_m3_s: s.largest,
            largest_mouth_lat_lon_deg: [lat.to_degrees(), lon.to_degrees()],
            eroded_km3: out.erosion.map_or(0.0, |e| e.cut / 1e9),
            laid_km3: out.erosion.map_or(0.0, |e| e.laid / 1e9),
            brought_to_sea_km3: out.erosion.map_or(0.0, |e| e.to_sea / 1e9),
            erosion_mm_yr: out.erosion.map_or(0.0, |e| 1e3 * e.rate),
        }
    });
    let water = out.water_summary.map(|w| WaterInfo {
        vapour: out.water_vapour.unwrap_or("").to_string(),
        precipitation_kg_m2_yr: w.precipitation,
        evaporation_kg_m2_yr: w.evaporation,
        land_precipitation_kg_m2_yr: w.land_precipitation,
        land_evaporation_kg_m2_yr: w.land_evaporation,
        land_runoff_kg_m2_yr: w.runoff,
        arid_land_fraction: w.arid,
        hyper_arid_land_fraction: w.hyper_arid,
    });
    let clouds = out.cloud_summary.map(|c| CloudsInfo {
        cover: c.cover,
        low: c.kinds[0],
        layers: c.kinds[1],
        anvils: c.kinds[2],
        towers: c.kinds[3],
        albedo_added: c.albedo,
        outgoing_change_w_m2: c.olr,
        downward_change_w_m2: c.down,
    });
    let wind = out.wind_summary.map(|w| WindInfo {
        mean_speed_m_s: w.speed,
        fastest_m_s: w.fastest,
        work_w_m2: w.work,
        thermal_equator_deg: [w.equator.0, w.equator.1],
        trades_end_deg: [w.ridges.0, w.ridges.1],
    });

    // Units, nearest cell per pixel, as an index map: the id times 8.
    let unit_ids = out.layers.unit_id.as_slice();
    let unit_eq: Vec<u8> = (0..width * height_px)
        .map(|k| {
            let (x, y) = (k % width, k / width);
            let lon = ((x as f64 + 0.5) / width as f64 - 0.5) * 2.0 * std::f64::consts::PI;
            let lat = (0.5 - (y as f64 + 0.5) / height_px as f64) * std::f64::consts::PI;
            let d = crate::grid::direction_from_lon_lat(lon, lat);
            (unit_ids[grid.index(grid.locate(d).0)] as u8).saturating_mul(8)
        })
        .collect();
    write_png8(&settings.dir.join("units.png"), width, height_px, &unit_eq)?;
    // The hidden layers, nearest cell per pixel, for inspecting what the finish wrote:
    // only when asked for, since they are not part of the map set.
    if std::env::var_os("EXOGEN_DUMP_LAYERS").is_some() {
        let dump = |name: &str, data: &[f32], lo: f32, hi: f32| -> Result<(), ExportError> {
            let eq: Vec<f32> = (0..width * height_px)
                .map(|k| {
                    let (x, y) = (k % width, k / width);
                    let lon = ((x as f64 + 0.5) / width as f64 - 0.5) * 2.0 * std::f64::consts::PI;
                    let lat = (0.5 - (y as f64 + 0.5) / height_px as f64) * std::f64::consts::PI;
                    let d = crate::grid::direction_from_lon_lat(lon, lat);
                    data[grid.index(grid.locate(d).0)]
                })
                .collect();
            write_png16(
                &settings.dir.join(format!("dbg_{name}.png")),
                width,
                height_px,
                &eq,
                lo,
                hi,
            )
        };
        dump("cover", out.layers.cover.as_slice(), 0.0, 1.0)?;
        dump("deep", out.layers.deep.as_slice(), 0.0, 1.0)?;
        dump("melt", out.layers.melt.as_slice(), 0.0, 1.0)?;
        dump("clean", out.layers.clean.as_slice(), 0.0, 1.0)?;
        dump("composition", out.layers.composition.as_slice(), -1.0, 1.0)?;
        dump("stain", out.layers.stain.as_slice(), 0.0, 1.0)?;
        dump("steep", out.layers.steep.as_slice(), 0.0, 1.0)?;
        dump("maturity", out.maturity.as_slice(), 0.0, 1.0)?;
        dump("rays", out.rays.as_slice(), 0.0, 1.0)?;
        // The lanes of grooved terrain as polylines in degrees, for plotting.
        #[derive(Serialize)]
        struct LaneDump {
            id: u64,
            time: f64,
            half_width: f64,
            half_length: f64,
            grooved: bool,
            pts: Vec<[f64; 2]>,
        }
        let lanes: Vec<LaneDump> = out
            .lanes
            .iter()
            .map(|l| LaneDump {
                id: l.id,
                time: l.time_gyr,
                half_width: l.half_width,
                half_length: l.half_length,
                grooved: l.grooved,
                pts: l
                    .points
                    .iter()
                    .map(|p| {
                        let (lon, lat) = crate::grid::lon_lat_from_direction(*p);
                        [lon.to_degrees(), lat.to_degrees()]
                    })
                    .collect(),
            })
            .collect();
        let path = settings.dir.join("dbg_lanes.json");
        std::fs::write(&path, serde_json::to_vec(&lanes)?).map_err(|source| ExportError::Io {
            path: path.display().to_string(),
            source,
        })?;
    }
    maps.insert(
        "units_png".to_string(),
        MapInfo {
            file: "units.png".into(),
            encoding: "8-bit grey PNG, unit id times 8".into(),
            units: "unit id, see unit_legend".into(),
            min: 0.0,
            max: 255.0,
            filter: "nearest".into(),
            ..Default::default()
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
                if !graph.is_pinned(def.id) {
                    "derived"
                } else if def.input {
                    "set"
                } else {
                    "pinned"
                }
            }
            Kind::Multiplier => "multiplier",
            Kind::StageOutput => "stage_output",
        };
        // A number in the unit it is labelled with (the review of the cryosphere maps,
        // 2026-11-01: Mars's year, stored in seconds, read as 59 million days).
        let value = match graph.get(def.id) {
            crate::params::Value::F64(v) => serde_json::to_value(def.display(*v))?,
            other => serde_json::to_value(other)?,
        };
        parameters.insert(
            def.id.to_string(),
            ParameterInfo {
                value,
                state: state.to_string(),
                unit: def.unit.to_string(),
                scale: def.scale,
            },
        );
    }
    let mut stage_outputs = BTreeMap::new();
    stage_outputs.insert("sea_level_m".to_string(), out.sea_level);
    stage_outputs.insert("ocean_fraction".to_string(), out.ocean_fraction);
    let (mean_t, equator_t, poles_t) = crate::climate::bands(&out.temperature);
    stage_outputs.insert("mean_temperature_k".to_string(), mean_t);
    stage_outputs.insert("equator_temperature_k".to_string(), equator_t);
    stage_outputs.insert("polar_temperature_k".to_string(), poles_t);
    // The share of its star's light the body reflected through its year (M67).
    if out.planetary_albedo.is_finite() {
        stage_outputs.insert("planetary_albedo".to_string(), out.planetary_albedo);
    }
    stage_outputs.insert("coldest_temperature_k".to_string(), out.temperature_range.0);
    // The poles' shadows (F243): the share of the body in them and the coldest the map holds.
    stage_outputs.insert("shadow_fraction".to_string(), out.shadow_fraction);
    stage_outputs.insert("sea_ice_fraction".to_string(), out.sea_ice_fraction);
    stage_outputs.insert("cryosphere_thickest_m".to_string(), thickest as f64);
    if out.cryosphere_snapshot.is_finite() {
        stage_outputs.insert(
            "cryosphere_snapshot_phase".to_string(),
            out.cryosphere_snapshot,
        );
    }
    let mapped = out
        .temperature
        .as_slice()
        .iter()
        .cloned()
        .fold(f32::MAX, f32::min);
    stage_outputs.insert("coldest_mapped_temperature_k".to_string(), mapped as f64);
    stage_outputs.insert("warmest_temperature_k".to_string(), out.temperature_range.1);
    let (frost_mean, ..) = crate::climate::bands(&out.frost);
    let cover = crate::fields::Field::from_fn(out.frost.grid(), |c, _| {
        crate::climate::frost_cover(out.frost.get(c)) as f32
    });
    stage_outputs.insert("frost_mean_kg_m2".to_string(), frost_mean);
    stage_outputs.insert(
        "frost_cover_fraction".to_string(),
        crate::climate::bands(&cover).0,
    );
    if let Some((lo, hi)) = out.pressure_range {
        stage_outputs.insert("pressure_lowest_pa".to_string(), lo);
        stage_outputs.insert("pressure_highest_pa".to_string(), hi);
    }
    if let Some(growth) = out.cap_growth {
        stage_outputs.insert("cap_growth_kg_m2_yr".to_string(), growth);
    }
    if let Some(gain) = out.reservoir_gain {
        stage_outputs.insert("reservoir_gain_kg_m2_yr".to_string(), gain);
    }
    // The cold traps (F240): each ice's coldest and the vapour it holds the air at; the
    // frost's own gas's ice's best growth and the share of the body it would grow on.
    for (ice, _) in &out.traps.vapour {
        if let Some((lat, _, vapour)) = out.traps.coldest(ice.species) {
            let species = ice.species.to_lowercase();
            stage_outputs.insert(format!("trap_{species}_vapour_pa"), vapour);
            stage_outputs.insert(format!("trap_{species}_latitude_deg"), lat);
        }
    }
    if let Some((ice, growth)) = &out.traps.growth {
        let (most, share) = crate::climate::growth_reach(growth);
        let species = ice.species.to_lowercase();
        stage_outputs.insert(format!("trap_{species}_growth_kg_m2_yr"), most);
        stage_outputs.insert(format!("trap_{species}_growth_fraction"), share);
    }
    // The caps (F241): each ice's laid, m of global layer of its liquid, and of that the
    // sea's; the share of the body it covers and its thickest, m.
    for cap in &out.caps {
        let species = cap.ice.species.to_lowercase();
        let grid = cap.thickness.grid();
        let (mut covered, mut whole) = (0.0, 0.0);
        for (k, &t) in cap.thickness.as_slice().iter().enumerate() {
            let a = grid.cell_area(grid.cell(k));
            whole += a;
            if t > 0.0 {
                covered += a;
            }
        }
        let thickest = cap
            .thickness
            .as_slice()
            .iter()
            .cloned()
            .fold(0.0f32, f32::max);
        stage_outputs.insert(format!("cap_{species}_laid_m"), cap.laid);
        stage_outputs.insert(format!("cap_{species}_from_sea_m"), cap.from_sea);
        stage_outputs.insert(format!("cap_{species}_fraction"), covered / whole);
        stage_outputs.insert(format!("cap_{species}_thickest_m"), thickest as f64);
    }
    let sidecar = Sidecar {
        exogen_version: crate::VERSION.to_string(),
        project: settings.project_name.clone(),
        width,
        height: height_px,
        radius_m: graph.f64("body.radius"),
        datum,
        datum_offset_m: datum_offset,
        preview_exaggeration: exaggeration,
        epoch_gyr_before_present: graph.f64("output.epoch"),
        epoch_time_since_formation_gyr: out.epoch_time_gyr,
        maps: maps
            .into_iter()
            .map(|(name, mut info)| {
                if info.filter.is_empty() {
                    info.filter = "linear".into();
                }
                (name, info)
            })
            .collect(),
        cryosphere,
        wind,
        water,
        clouds,
        rivers,
        liquid,
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
        flood_count: out.floods.len(),
        magma_provinces: out
            .provinces
            .iter()
            .map(|p| {
                let (lon, lat) = crate::grid::lon_lat_from_direction(p.center);
                ProvinceInfo {
                    lat_deg: lat.to_degrees(),
                    lon_deg: lon.to_degrees(),
                    radius_deg: p.sigma.to_degrees(),
                }
            })
            .collect(),
        floods: largest_floods(&out.floods, SIDECAR_FLOOD_LIMIT)
            .into_iter()
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
        rises: out
            .rises
            .iter()
            .map(|r| {
                let (lon, lat) = crate::grid::lon_lat_from_direction(r.center);
                let (east, north) = crate::geom::local_basis(r.center);
                let azimuth = crate::geom::dot(r.axis, east)
                    .atan2(crate::geom::dot(r.axis, north))
                    .to_degrees()
                    .rem_euclid(180.0);
                RiseInfo {
                    lat_deg: lat.to_degrees(),
                    lon_deg: lon.to_degrees(),
                    diameter_km: 2.0 * r.radius / 1e3,
                    height_m: r.height,
                    paleo_lat_deg: r.paleo_latitude.to_degrees(),
                    axis_azimuth_deg: azimuth,
                }
            })
            .collect(),
        shields: out
            .shields
            .iter()
            // A volcano that never grew stands nowhere (F165).
            .filter(|sh| crate::stages::plume::shield_progress(out.epoch_time_gyr, sh) > 0.0)
            .map(|sh| {
                let (lon, lat) = crate::grid::lon_lat_from_direction(sh.center);
                ShieldInfo {
                    lat_deg: lat.to_degrees(),
                    lon_deg: lon.to_degrees(),
                    diameter_km: 2.0 * sh.radius / 1e3,
                    height_m: sh.height
                        * crate::stages::plume::shield_grown(out.epoch_time_gyr, sh),
                    calderas: sh.calderas.len(),
                    age_gyr: out.epoch_time_gyr - sh.start_gyr,
                }
            })
            .collect(),
        canyons: out
            .canyons
            .iter()
            .map(|t| trough_info(t, out.epoch_time_gyr))
            .collect(),
        graben: out
            .graben
            .iter()
            .map(|t| trough_info(t, out.epoch_time_gyr))
            .collect(),
        resurfacing_age_gyr: out
            .resurfacing
            .as_ref()
            .map(|r| out.epoch_time_gyr - r.time_gyr),
        tessera_fraction: crate::stages::episodic::tessera_fraction(&out.layers),
        // Every corona takes the corona process (F133, F145): an ice mantle's its outline's
        // mean diameter; a rock lid's its profile at the epoch too, at its centre and its
        // rim's crest.
        coronae: out
            .ice_coronae
            .iter()
            .map(|c| {
                let (lon, lat) = crate::grid::lon_lat_from_direction(c.center);
                let (centre, rim) = match &c.head {
                    Some(h) => {
                        let a = c.outline_radius();
                        let at = |r: f64| h.profile(r, c.time_gyr, out.epoch_time_gyr);
                        (
                            at(0.0),
                            (0..=20)
                                .map(|i| at(a * (0.8 + 0.02 * i as f64)))
                                .fold(f64::NEG_INFINITY, f64::max),
                        )
                    }
                    None => (0.0, 0.0),
                };
                CoronaInfo {
                    lat_deg: lat.to_degrees(),
                    lon_deg: lon.to_degrees(),
                    diameter_km: 2.0 * c.outline_radius() / 1e3,
                    centre_relief_m: centre,
                    rim_height_m: rim,
                    age_gyr: out.epoch_time_gyr - c.time_gyr,
                }
            })
            .collect(),
        domes: out
            .domes
            .iter()
            .map(|d| {
                let (lon, lat) = crate::grid::lon_lat_from_direction(d.center);
                DomeInfo {
                    lat_deg: lat.to_degrees(),
                    lon_deg: lon.to_degrees(),
                    diameter_km: 2.0 * d.radius / 1e3,
                    height_m: d.height,
                    age_gyr: out.epoch_time_gyr - d.time_gyr,
                }
            })
            .collect(),
        paterae: out
            .paterae
            .iter()
            .map(|p| {
                let (lon, lat) = crate::grid::lon_lat_from_direction(p.center);
                PateraInfo {
                    lat_deg: lat.to_degrees(),
                    lon_deg: lon.to_degrees(),
                    diameter_km: 2.0 * p.extent() / 1e3,
                    depth_m: p.depth,
                    active: p.active,
                }
            })
            .collect(),
        mountains: out
            .mountains
            .iter()
            .map(|m| {
                let (lon, lat) = crate::grid::lon_lat_from_direction(m.center);
                MountainInfo {
                    lat_deg: lat.to_degrees(),
                    lon_deg: lon.to_degrees(),
                    length_km: m.length / 1e3,
                    width_km: m.width / 1e3,
                    height_m: m.height,
                }
            })
            .collect(),
        // Io's lava erupts as floods from its paterae (F117).
        flow_field_count: if out.paterae.is_empty() {
            0
        } else {
            out.floods.len()
        },
        plume_count: out
            .plume_episodes
            .iter()
            .filter(|e| e.end_yr == 0.0)
            .count(),
        plume_episode_count: out.plume_episodes.len(),
        small_eruption_count: out.small_eruption_count,
        linea_count: out.lineae.len(),
        lane_count: out.lanes.len(),
        chaos_count: out.chaos.len(),
        chaos_fraction: out.chaos_fraction,
        cryodome_count: out.cryodomes.len(),
        lane_coverage: out.lane_coverage,
        grooved_fraction: out.grooved_fraction,
        polar_terrain_fraction: out.polar_fraction,
        fossil_terrain_count: out.episodes.len(),
        renewed_fraction: out.renewed_fraction,
        lineae: {
            let mut sorted: Vec<&crate::stages::lineae::Linea> = out.lineae.iter().collect();
            sorted.sort_by(|a, b| b.length.total_cmp(&a.length));
            sorted
                .into_iter()
                .take(SIDECAR_FLOOD_LIMIT)
                .map(|l| {
                    let (lon, lat) = crate::grid::lon_lat_from_direction(l.midpoint());
                    LineaInfo {
                        lat_deg: lat.to_degrees(),
                        lon_deg: lon.to_degrees(),
                        length_km: l.length / 1e3,
                        width_km: l.width / 1e3,
                        kind: match l.kind {
                            crate::stages::lineae::LineaKind::Ridge => "ridge",
                            crate::stages::lineae::LineaKind::Cycloid => "cycloid",
                            crate::stages::lineae::LineaKind::Band => "band",
                            crate::stages::lineae::LineaKind::Vent => "vent",
                            crate::stages::lineae::LineaKind::Fracture => "fracture",
                        }
                        .to_string(),
                        age_gyr: out.epoch_time_gyr - l.time_gyr,
                    }
                })
                .collect()
        },
        lava_area_fraction: crate::stages::volcanism::lava_fraction(&out.layers),
        unit_legend: crate::units::legend()
            .into_iter()
            .map(|(id, name)| UnitInfo {
                id,
                name: name.to_string(),
            })
            .collect(),
    };
    // The features, most of the sidecar's size, apart in features.json, compact (F253).
    let mut whole = serde_json::to_value(&sidecar)?;
    let mut features = serde_json::Map::new();
    if let Some(object) = whole.as_object_mut() {
        for key in FEATURE_LISTS {
            if let Some(list) = object.remove(key) {
                features.insert(key.to_string(), list);
            }
        }
        object.insert("features".to_string(), "features.json".into());
    }
    for (name, text) in [
        ("sidecar.json", serde_json::to_string_pretty(&whole)?),
        (
            "features.json",
            serde_json::to_string(&serde_json::Value::Object(features))?,
        ),
    ] {
        let path = settings.dir.join(name);
        retry(|| std::fs::write(&path, &text)).map_err(|source| ExportError::Io {
            path: path.display().to_string(),
            source,
        })?;
    }
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

fn write_exr_channels(
    path: &Path,
    (width, height): (usize, usize),
    layer: &str,
    names: &[&str],
    data: &[Vec<f32>],
) -> Result<(), ExportError> {
    let channels = names
        .iter()
        .zip(data)
        .map(|(name, v)| AnyChannel::new(*name, FlatSamples::F32(v.clone())))
        .collect();
    let layer = Layer::new(
        (width, height),
        LayerAttributes::named(layer),
        Encoding::FAST_LOSSLESS,
        AnyChannels::sort(channels),
    );
    let image = Image::from_layer(layer);
    retry(|| image.write().non_parallel().to_file(path))
        .map_err(|e| ExportError::Exr(e.to_string()))
}

/// Three channels from 0 to 1 to a 16-bit RGB PNG (F253).
fn write_png_rgb16(
    path: &Path,
    (width, height): (usize, usize),
    channels: &[Vec<f32>],
) -> Result<(), ExportError> {
    let mut pixels = Vec::with_capacity(width * height * 3);
    for k in 0..width * height {
        for channel in channels {
            pixels.push((channel[k].clamp(0.0, 1.0) * 65535.0).round() as u16);
        }
    }
    let img = image::ImageBuffer::<image::Rgb<u16>, Vec<u16>>::from_raw(
        width as u32,
        height as u32,
        pixels,
    )
    .expect("buffer size");
    retry(|| img.save(path))?;
    Ok(())
}

/// Four channels from 0 to 1 to a 16-bit RGBA PNG (M66).
fn write_png_rgba16(
    path: &Path,
    (width, height): (usize, usize),
    channels: &[Vec<f32>],
) -> Result<(), ExportError> {
    let mut pixels = Vec::with_capacity(width * height * 4);
    for k in 0..width * height {
        for channel in channels {
            pixels.push((channel[k].clamp(0.0, 1.0) * 65535.0).round() as u16);
        }
    }
    let img = image::ImageBuffer::<image::Rgba<u16>, Vec<u16>>::from_raw(
        width as u32,
        height as u32,
        pixels,
    )
    .expect("buffer size");
    retry(|| img.save(path))?;
    Ok(())
}

/// A field on the cube to equirectangular pixels (F277), each the most of the cells that fall in
/// it, or where none does the cell its centre falls in: so that a line a cell wide stays
/// unbroken where the pixels are the coarser.
fn to_equirectangular_most(
    grid: &crate::grid::CubeSphere,
    values: &[f32],
    width: usize,
    height: usize,
) -> Vec<f32> {
    let mut out = to_equirectangular_nearest(grid, values, width, height);
    let mut most = vec![f32::NEG_INFINITY; width * height];
    for (k, &v) in values.iter().enumerate() {
        let (lon, lat) = crate::grid::lon_lat_from_direction(grid.center(grid.cell(k)));
        let x = ((lon / (2.0 * std::f64::consts::PI) + 0.5) * width as f64).floor() as usize;
        let y = ((0.5 - lat / std::f64::consts::PI) * height as f64).floor() as usize;
        let p = y.min(height - 1) * width + x.min(width - 1);
        most[p] = most[p].max(v);
    }
    for (o, m) in out.iter_mut().zip(most) {
        if m.is_finite() {
            *o = m;
        }
    }
    out
}

/// A field on the cube to equirectangular pixels, each the value of the cell its centre
/// falls in.
fn to_equirectangular_nearest(
    grid: &crate::grid::CubeSphere,
    values: &[f32],
    width: usize,
    height: usize,
) -> Vec<f32> {
    (0..width * height)
        .map(|k| {
            let (x, y) = (k % width, k / width);
            let lon = ((x as f64 + 0.5) / width as f64 - 0.5) * 2.0 * std::f64::consts::PI;
            let lat = (0.5 - (y as f64 + 0.5) / height as f64) * std::f64::consts::PI;
            let d = crate::grid::direction_from_lon_lat(lon, lat);
            values[grid.index(grid.locate(d).0)]
        })
        .collect()
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
