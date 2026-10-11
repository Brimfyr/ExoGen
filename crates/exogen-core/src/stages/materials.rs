// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! Material model (PLAN.md 4.10): per-cell mixtures of end-members, each with a fresh and
//! a space-weathered colour, blended by optical maturity and brightened by rays. Colours
//! are linear reflectance. The spectral library and the integration against colour
//! matching functions arrive with the colour-management work in Phase 7; the colours here
//! are the integrated values of typical published spectra.
//!
//! Two controls shape a surface: the bedrock is what impacts excavate, so a fresh surface
//! shows its fresh colour; the regolith mix is the weathered surface, so a mature surface
//! shows the mixture's mature colour. Maturity blends between the two.

use std::collections::BTreeMap;

use rayon::prelude::*;

use crate::constants::*;
use crate::fields::{Field, LayerSet};
use crate::geom::dot;
use crate::grid::CubeSphere;
use crate::params::ParamGraph;
use crate::units::{
    UNIT_BAND, UNIT_CHAOS, UNIT_CONTINENTAL, UNIT_CRATER_FLOOR, UNIT_CRYODOME, UNIT_DEEP_FLOOR,
    UNIT_IO_FLOW, UNIT_LINEA, UNIT_MOUNTAIN, UNIT_OCEANIC, UNIT_PATERA_FLOOR,
};

pub const VERSION: u64 = 18;

#[derive(Clone, Copy, Debug)]
pub struct EndMember {
    pub id: &'static str,
    pub fresh: [f32; 3],
    pub mature: [f32; 3],
    pub roughness: f32,
}

/// The initial library. Fresh colours are immature regolith or clean deposits; mature
/// colours are after space weathering (darkened and reddened on silicates).
pub const END_MEMBERS: &[EndMember] = &[
    EndMember {
        id: "anorthosite",
        fresh: [0.34, 0.34, 0.35],
        mature: [0.150, 0.138, 0.125],
        roughness: 0.85,
    },
    EndMember {
        id: "basalt",
        fresh: [0.100, 0.090, 0.085],
        mature: [0.065, 0.058, 0.052],
        roughness: 0.80,
    },
    EndMember {
        id: "felsic",
        fresh: [0.26, 0.23, 0.19],
        mature: [0.24, 0.21, 0.17],
        roughness: 0.75,
    },
    EndMember {
        id: "mercurian_regolith",
        fresh: [0.17, 0.17, 0.17],
        mature: [0.085, 0.085, 0.085],
        roughness: 0.85,
    },
    EndMember {
        id: "water_ice",
        fresh: [0.80, 0.82, 0.86],
        mature: [0.64, 0.65, 0.68],
        roughness: 0.55,
    },
    EndMember {
        id: "water_ice_dirty",
        fresh: [0.36, 0.35, 0.34],
        mature: [0.26, 0.25, 0.24],
        roughness: 0.70,
    },
    // Ceres's crust: hydrated silicates, carbonates and clathrates with a third of ice,
    // dark even when fresh (Haulani's ejecta are 1.6 times the surroundings).
    EndMember {
        id: "carbonaceous_ice",
        fresh: [0.160, 0.155, 0.150],
        mature: [0.095, 0.090, 0.085],
        roughness: 0.80,
    },
    EndMember {
        id: "carbonaceous",
        fresh: [0.120, 0.115, 0.110],
        mature: [0.090, 0.085, 0.080],
        roughness: 0.85,
    },
    EndMember {
        id: "ammonia_water_ice",
        fresh: [0.70, 0.72, 0.75],
        mature: [0.50, 0.50, 0.52],
        roughness: 0.55,
    },
    EndMember {
        id: "rock_ice_mixture",
        fresh: [0.45, 0.44, 0.43],
        mature: [0.30, 0.29, 0.28],
        roughness: 0.75,
    },
    EndMember {
        id: "sulfur_coated",
        fresh: [0.65, 0.55, 0.25],
        mature: [0.50, 0.42, 0.20],
        roughness: 0.70,
    },
    // Lunar lower crust: noritic, darker than the anorthositic highlands even when fresh.
    EndMember {
        id: "norite",
        fresh: [0.24, 0.23, 0.22],
        mature: [0.11, 0.10, 0.09],
        roughness: 0.85,
    },
    // Mercury's low-reflectance material, excavated by basins and large craters.
    EndMember {
        id: "low_reflectance_material",
        fresh: [0.11, 0.11, 0.12],
        mature: [0.065, 0.065, 0.07],
        roughness: 0.85,
    },
    // Mare basalt: iron- and titanium-rich flood lavas, dark and faintly blue-grey when
    // mature, brighter where fresh craters turn them over.
    EndMember {
        id: "mare_basalt",
        fresh: [0.19, 0.18, 0.165],
        mature: [0.075, 0.072, 0.070],
        roughness: 0.55,
    },
    // Mercury's smooth plains: low-iron lavas a shade brighter and redder than the
    // intercrater plains around them (Denevi et al. 2013).
    EndMember {
        id: "smooth_plains",
        fresh: [0.27, 0.255, 0.235],
        mature: [0.125, 0.112, 0.098],
        roughness: 0.55,
    },
    // Glassy impact melt on silicate bodies: distinctly dark when fresh, buried by regolith
    // with age.
    // Hydrated salts and sulfuric acid hydrates from an ocean, brown-red when fresh,
    // paling as frost settles on them.
    EndMember {
        id: "hydrated_salt",
        fresh: [0.40, 0.26, 0.17],
        mature: [0.52, 0.42, 0.34],
        roughness: 0.6,
    },
    // Sodium carbonates and chlorides from an ocean where no magnetosphere irradiates
    // them: white, greying as they take dust (Ceres's faculae, Enceladus's stripes).
    EndMember {
        id: "carbonate_salt",
        fresh: [0.86, 0.85, 0.80],
        mature: [0.72, 0.70, 0.66],
        roughness: 0.6,
    },
    EndMember {
        id: "impact_melt",
        fresh: [0.10, 0.095, 0.09],
        mature: [0.07, 0.065, 0.06],
        roughness: 0.60,
    },
];

/// Ice end-members, whose excavated bright material darkens back to the lag it came from.
pub fn is_ice(id: &str) -> bool {
    id.starts_with("water_ice")
        || id == "ammonia_water_ice"
        || id == "rock_ice_mixture"
        || id == "carbonaceous_ice"
}

/// The share of an end-member that is not ice and stays behind as lag when the ice
/// sublimes: clean ices a trace, dirty ice a third, a rock-ice mixture half, rock all.
pub fn non_ice_fraction(id: &str) -> f64 {
    match id {
        "water_ice" | "ammonia_water_ice" => 0.02,
        "water_ice_dirty" => 0.3,
        "carbonaceous_ice" => 0.4,
        "rock_ice_mixture" => 0.5,
        _ => 1.0,
    }
}

pub fn end_member(id: &str) -> Option<&'static EndMember> {
    END_MEMBERS.iter().find(|e| e.id == id)
}

pub fn end_member_ids() -> Vec<&'static str> {
    END_MEMBERS.iter().map(|e| e.id).collect()
}

fn lerp3(a: [f32; 3], b: [f32; 3], t: f32) -> [f32; 3] {
    [
        a[0] + (b[0] - a[0]) * t,
        a[1] + (b[1] - a[1]) * t,
        a[2] + (b[2] - a[2]) * t,
    ]
}

/// Colour of a mixture at a given maturity.
fn mixture_colour(mix: &[(&'static EndMember, f32)], maturity: f32) -> [f32; 3] {
    let mut out = [0.0f32; 3];
    for (e, f) in mix {
        let c = lerp3(e.fresh, e.mature, maturity);
        out[0] += c[0] * f;
        out[1] += c[1] * f;
        out[2] += c[2] * f;
    }
    out
}

/// Resolve a map of end-member fractions into references, normalised to sum to one.
fn resolve_mix(map: &BTreeMap<String, f64>) -> Vec<(&'static EndMember, f32)> {
    let mut out: Vec<(&'static EndMember, f32)> = map
        .iter()
        .filter_map(|(k, v)| end_member(k).map(|e| (e, *v as f32)))
        .filter(|(_, f)| *f > 0.0)
        .collect();
    let total: f32 = out.iter().map(|(_, f)| f).sum();
    if out.is_empty() || total <= 0.0 {
        return vec![(end_member("basalt").unwrap(), 1.0)];
    }
    for (_, f) in &mut out {
        *f /= total;
    }
    out
}

/// A fresh colour and a mature colour between which maturity blends.
#[derive(Clone, Copy)]
struct Pair {
    fresh: [f32; 3],
    mature: [f32; 3],
}

/// Lava coverage mixed sideways by impact gardening: the coverage blurred over the
/// mixing zone of the oldest surfaces, blended in by how far each cell's own age has
/// let the zone grow, so fresh flow margins stay sharp and old mare margins soften.
fn mixed_cover(
    grid: CubeSphere,
    body_radius: f64,
    cover: &[f32],
    ages: &[f32],
    present: f64,
) -> Vec<f32> {
    let cell = grid.delta() * body_radius;
    let sigma_full =
        REGOLITH_MIXING_LENGTH * (present / REGOLITH_MIXING_REFERENCE_GYR).sqrt().max(1e-3);
    // Each pass of the 1/8, 1/2, 1/8 kernel adds a quarter of a cell of variance.
    let passes = (4.0 * (sigma_full / cell).powi(2)).round() as usize;
    if passes == 0 {
        return cover.to_vec();
    }
    let passes = passes.min(REGOLITH_MIXING_MAX_PASSES);
    let sigma_blur = 0.5 * (passes as f64).sqrt() * cell;
    let blurred = crate::stages::episodic::smooth(&grid, cover, passes);
    cover
        .iter()
        .zip(blurred.iter())
        .zip(ages.iter())
        .map(|((&c, &b), &a)| {
            let exposure = (present - a as f64).max(0.0);
            let sigma = REGOLITH_MIXING_LENGTH * (exposure / REGOLITH_MIXING_REFERENCE_GYR).sqrt();
            let f = (sigma / sigma_blur).clamp(0.0, 1.0) as f32;
            c + f * (b - c)
        })
        .collect()
}

/// Linear-reflectance albedo of the surface.
pub fn albedo(
    graph: &ParamGraph,
    layers: &LayerSet,
    maturity: &Field<f32>,
    rays: &Field<f32>,
    present: f64,
    slope: Option<&Field<f32>>,
) -> Field<[f32; 3]> {
    let bedrock =
        end_member(graph.str("materials.bedrock")).unwrap_or_else(|| end_member("basalt").unwrap());
    let regolith = graph
        .get("materials.regolith_mix")
        .as_map()
        .map(resolve_mix)
        .unwrap_or_else(|| vec![(bedrock, 1.0)]);
    let base = Pair {
        fresh: bedrock.fresh,
        mature: mixture_colour(&regolith, 1.0),
    };
    // What deep craters excavate. Silicates keep their contrast when mature; ice brought
    // up from under a lag darkens back to the lag.
    let deep = end_member(graph.str("materials.deep_bedrock")).unwrap_or(bedrock);
    let deep_pair = Pair {
        fresh: deep.fresh,
        mature: if is_ice(deep.id) {
            base.mature
        } else {
            deep.mature
        },
    };
    // Melt: dark glass on silicates, refrozen clean ice on ice shells; either converges on
    // the regolith as gardening buries it.
    let melt = end_member("impact_melt").unwrap();
    let melt_pair = Pair {
        fresh: if is_ice(deep.id) {
            deep.fresh
        } else {
            melt.fresh
        },
        mature: base.mature,
    };
    let felsic = end_member("felsic").unwrap();
    let basalt = end_member("basalt").unwrap();
    let continental = Pair {
        fresh: felsic.fresh,
        mature: felsic.mature,
    };
    let oceanic = Pair {
        fresh: basalt.fresh,
        mature: basalt.mature,
    };
    let lava_member = end_member(graph.str("materials.lava")).unwrap_or(bedrock);
    let lava = Pair {
        fresh: lava_member.fresh,
        mature: lava_member.mature,
    };
    let (deep, shallow): ([f32; 3], [f32; 3]) = match graph.str("volatiles.working_liquid") {
        "water" | "ammonia_water" => ([0.012, 0.035, 0.085], [0.06, 0.22, 0.26]),
        "methane_ethane" => ([0.02, 0.02, 0.015], [0.07, 0.07, 0.05]),
        _ => ([0.03, 0.03, 0.03], [0.08, 0.08, 0.08]),
    };
    let grid = layers.grid;
    let liquid = layers.liquid.as_slice();
    let units = layers.unit_id.as_slice();
    let deep_f = layers.deep.as_slice();
    let melt_f = layers.melt.as_slice();
    let cover_f = layers.cover.as_slice();
    let clean_f = layers.clean.as_slice();
    let m = maturity.as_slice();
    // On an airless surface small impacts mix the regolith sideways, and a contact
    // between lava and its substrate blurs in colour as it ages.
    let io = graph.str("materials.bedrock") == "sulfur_coated";
    let gardened = !io && graph.f64("atmosphere.surface_pressure") < SPACE_WEATHERING_MAX_PRESSURE;
    // An icy surface wears frost where it is cold enough to keep what the warm
    // latitudes sublimed.
    let frost = if is_ice(graph.str("materials.bedrock")) && !io {
        Some(crate::stages::ice::FrostCover::new(graph))
    } else {
        None
    };
    let frost_colour = end_member("water_ice").unwrap().fresh;
    let salt_colour = end_member(graph.str("materials.salt"))
        .unwrap_or_else(|| end_member("hydrated_salt").unwrap())
        .fresh;
    // A plume's fallout, and the plasma on a locked moon's trailing hemisphere.
    let fallout = frost
        .as_ref()
        .and_then(|_| crate::stages::cryo::Fallout::new(graph));
    let trailing = graph.f64("surface.trailing_irradiation");
    // Ice renewed from below, the grooved lanes and the polar terrain, is clean whatever
    // lies around it.
    let clean_ice = {
        let w = end_member("water_ice").unwrap();
        Pair {
            fresh: w.fresh,
            mature: w.mature,
        }
    };
    let ages_f = layers.surface_age.as_slice();
    let mixed = if gardened {
        mixed_cover(
            grid,
            graph.f64("body.radius"),
            cover_f,
            layers.surface_age.as_slice(),
            present,
        )
    } else {
        cover_f.to_vec()
    };
    let r = rays.as_slice();
    // The Io regime colours itself from what the eruptions left and the temperature the
    // latitude's insolation gives: the day's peak temperature decides where the global
    // frost fall survives and where radiation keeps sulfur red.
    let io_subsolar_k = {
        let absorbed = graph.f64("orbit.insolation") * (1.0 - graph.f64("atmosphere.bond_albedo"));
        (absorbed.max(0.0) / SIGMA_SB).powf(0.25)
    };
    let bedrock_f = layers.bedrock.as_slice();
    let io_cell = grid.delta() * graph.f64("body.radius");
    let data: Vec<[f32; 3]> = (0..grid.cell_count())
        .into_par_iter()
        .map(|k| {
            if io {
                let d = grid.center(grid.cell(k));
                let cos_lat = (1.0 - d[2] * d[2]).max(0.0).sqrt();
                let t_peak = io_subsolar_k * cos_lat.max(0.02).powf(0.25);
                let logistic = |x: f64| 1.0 / (1.0 + (-x).exp());
                let persist = logistic((IO_FROST_HALF_K - t_peak) / IO_FROST_WIDTH_K) as f32;
                let red = logistic((IO_RED_HALF_K - t_peak) / IO_RED_WIDTH_K) as f32;
                let plains = lerp3(IO_SULFUR, IO_FROST, IO_FROST_BACKGROUND as f32 * persist);
                let plains = lerp3(plains, IO_POLAR, red);
                // What lies on the cell: frost, or on lava the lava's own coverage; the
                // fresh silicate exposure; the red fallout over all but hot lava.
                let coat = cover_f[k].clamp(0.0, 1.0);
                let fresh = deep_f[k].clamp(0.0, 1.0);
                let c = match units[k] {
                    UNIT_IO_FLOW => lerp3(lerp3(plains, IO_FLOW_OLD, coat), IO_FLOW_DARK, fresh),
                    UNIT_PATERA_FLOOR => lerp3(
                        lerp3(plains, IO_PATERA_FLOOR, 0.8 * coat),
                        IO_FLOW_DARK,
                        fresh,
                    ),
                    UNIT_MOUNTAIN => {
                        // Steep flanks shed their coat and show the rock.
                        let cell = grid.cell(k);
                        let h = |di: i32, dj: i32| {
                            bedrock_f[grid.index(grid.neighbor(cell, di, dj))] as f64
                        };
                        let slope = ((h(1, 0) - h(-1, 0)).powi(2) + (h(0, 1) - h(0, -1)).powi(2))
                            .sqrt()
                            / (2.0 * io_cell);
                        let (lo, hi) = IO_MOUNTAIN_BARE_SLOPE_DEG;
                        let bare =
                            ((slope.atan().to_degrees() - lo) / (hi - lo)).clamp(0.0, 1.0) as f32;
                        lerp3(
                            lerp3(plains, IO_FROST, coat),
                            IO_MOUNTAIN_ROCK,
                            0.25 + 0.75 * bare,
                        )
                    }
                    _ => lerp3(plains, IO_FROST, coat),
                };
                return lerp3(c, IO_PLUME_RED, melt_f[k].clamp(0.0, 1.0) * (1.0 - fresh));
            }
            if liquid[k] > 0.0 {
                let t = (-liquid[k] / 60.0).exp();
                lerp3(deep, shallow, t)
            } else {
                let pair = match units[k] {
                    UNIT_CONTINENTAL => continental,
                    UNIT_OCEANIC => oceanic,
                    _ => base,
                };
                // Ice renewed from below stays clean under whatever cuts it.
                let pair = if clean_f[k] > 0.0 {
                    Pair {
                        fresh: lerp3(pair.fresh, clean_ice.fresh, clean_f[k]),
                        mature: lerp3(pair.mature, clean_ice.mature, clean_f[k]),
                    }
                } else {
                    pair
                };
                // The substrate: deep-sourced material by its fraction. Lava lies over it
                // by its coverage, mixed sideways into its surroundings, so a thin flow
                // shows the substrate through and an old contact is blurred; then melt
                // over the fresh surface.
                let fresh = lerp3(pair.fresh, deep_pair.fresh, deep_f[k]);
                let mature = lerp3(pair.mature, deep_pair.mature, deep_f[k]);
                let c = mixed[k].clamp(0.0, 1.0);
                let salted = matches!(
                    units[k],
                    UNIT_LINEA | UNIT_BAND | UNIT_CHAOS | UNIT_CRYODOME
                );
                let (fresh, mature) = if c > 0.0
                    && !salted
                    && units[k] != UNIT_CONTINENTAL
                    && units[k] != UNIT_OCEANIC
                {
                    (lerp3(fresh, lava.fresh, c), lerp3(mature, lava.mature, c))
                } else {
                    (fresh, mature)
                };
                let fresh = lerp3(fresh, melt_pair.fresh, melt_f[k]);
                let effective = (m[k] * (1.0 - r[k])).clamp(0.0, 1.0);
                let colour = lerp3(fresh, mature, effective);
                match &frost {
                    Some(f) => {
                        let d = grid.center(grid.cell(k));
                        // Frost: what the cold latitudes keep, and what a plume lays down.
                        let tilt = slope.map_or(0.0, |s| s.as_slice()[k] as f64);
                        let mut frost_share = f.at_slope(d, tilt);
                        let mut burial = f64::INFINITY;
                        if let Some(fall) = &fallout {
                            frost_share = frost_share.max(fall.frost_at(d));
                            burial = fall.burial_gyr(d);
                        }
                        let colour = lerp3(colour, frost_colour, frost_share as f32);
                        // The ocean's salts: in a crack or a chaos matrix fading under
                        // frost, on a dome or a deep crater's floor darkening with age,
                        // and under a plume buried as fast as they fall.
                        let fade = match units[k] {
                            UNIT_LINEA | UNIT_BAND | UNIT_CHAOS => Some(SALT_FADE_GYR),
                            UNIT_CRYODOME | UNIT_CRATER_FLOOR | UNIT_DEEP_FLOOR => {
                                Some(FACULA_FADE_GYR)
                            }
                            _ => None,
                        };
                        let colour = match fade {
                            Some(fade) if cover_f[k] > 0.0 => {
                                let age = (present - ages_f[k] as f64).max(0.0);
                                let salt = cover_f[k] * (-age / fade.min(burial)).exp() as f32;
                                lerp3(colour, salt_colour, salt.clamp(0.0, 1.0))
                            }
                            _ => colour,
                        };
                        // A plume's warm vents sinter the ice in the fissures to coarse
                        // grains, darker and bluer than the frost round them.
                        let colour = match &fallout {
                            Some(fall) if units[k] == UNIT_LINEA || units[k] == UNIT_BAND => {
                                let w = (cover_f[k] as f64 * fall.vent_weight(d)) as f32;
                                let coarse = [
                                    colour[0] * VENT_ICE_TINT[0],
                                    colour[1] * VENT_ICE_TINT[1],
                                    colour[2] * VENT_ICE_TINT[2],
                                ];
                                lerp3(colour, coarse, w.clamp(0.0, 1.0))
                            }
                            _ => colour,
                        };
                        // The trailing hemisphere's plasma darkens and reddens the ice.
                        if trailing > 0.0 {
                            let w = (trailing
                                * dot(d, crate::stages::cryo::TRAILING_APEX).max(0.0).powi(2))
                                as f32;
                            let tinted = [
                                colour[0] * TRAILING_TINT[0],
                                colour[1] * TRAILING_TINT[1],
                                colour[2] * TRAILING_TINT[2],
                            ];
                            lerp3(colour, tinted, w)
                        } else {
                            colour
                        }
                    }
                    None => colour,
                }
            }
        })
        .collect();
    let mut f = Field::new(grid);
    f.as_mut_slice().copy_from_slice(&data);
    f
}
