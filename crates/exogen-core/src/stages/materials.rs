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
    UNIT_LINEA, UNIT_MOUNTAIN, UNIT_OCEANIC,
};

pub const VERSION: u64 = 23;

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

/// What gardening mixes sideways (F95): the lava's coverage, and its coverage times its
/// tint, over the mixing zone of the contacts each cell lies by. Small impacts move the
/// regolith about as a diffusion does, so a contact's zone grows as the square root of
/// the gardening since it formed, the body's diffusivity times the impacts' equivalent
/// time (the law that degrades its craters), from the zone Li and Mustard measure across
/// the lunar maria's contacts. A contact is as old as its younger side, so a cell of
/// lava mixes over no more than its own lava's zone, and any cell over no more than the
/// youngest lava within reach allows, a power mean weighted to the least gardened: a
/// young flow on old mare keeps a sharp contact on both sides, and old highlands mix
/// with the mare beside them only for as long as the mare has lain there. `dose` is each
/// lava cell's gardening relative to the reference's, its zone the reference's times its
/// square root. The zones are the levels of a Gaussian pyramid: a cell takes the widest
/// its contacts allow, between two levels by its share of the step.
fn garden_mix(
    grid: CubeSphere,
    cell: f64,
    cover: &[f32],
    tinted: &[f32],
    dose: &[Option<f32>],
) -> (Vec<f32>, Vec<f32>) {
    use crate::stages::relaxation::{pyramid_coarse, pyramid_expand, pyramid_spread};
    let zone = |d: f64| REGOLITH_MIXING_LENGTH * d.max(0.0).sqrt();
    let spread = |l: usize| cell * pyramid_spread(l);
    let widest = zone(dose.iter().flatten().fold(0.0f32, |a, &d| a.max(d)) as f64);
    if widest < 0.1 * spread(0) {
        return (cover.to_vec(), tinted.to_vec());
    }
    let wanted = (0..24).find(|&l| spread(l) >= widest).unwrap_or(23) + 1;
    let p = REGOLITH_MIXING_YOUNGEST_POWER;
    let own: Vec<f64> = dose
        .iter()
        .map(|d| d.map_or(f64::INFINITY, |d| zone(d as f64)))
        .collect();
    // Each lava cell's weight toward the youngest, by its coverage.
    let youth: Vec<f32> = cover
        .iter()
        .zip(dose)
        .map(|(&c, d)| d.map_or(0.0, |d| c * ((d as f64).max(1e-6).powf(-p)) as f32))
        .collect();
    let field = |v: &[f32]| {
        let mut f = Field::<f32>::new(grid);
        f.as_mut_slice().copy_from_slice(v);
        f
    };
    let pyramids: Vec<_> = [cover, tinted, youth.as_slice()]
        .iter()
        .map(|v| pyramid_coarse(&field(v), wanted))
        .collect();
    // Per cell: the channels at the widest zone allowed so far, and whether a wider one
    // may still be.
    let mut state: Vec<(f32, f32, bool)> = cover
        .iter()
        .zip(tinted)
        .map(|(&c, &t)| (c, t, true))
        .collect();
    let mut below = 0.0;
    for l in 0..pyramids[0].1.len() {
        let here = spread(l);
        let level: Vec<Field<f32>> = pyramids
            .iter()
            .map(|(grids, coarse)| pyramid_expand(grids, &coarse[l], l))
            .collect();
        let (c_l, t_l, y_l) = (
            level[0].as_slice(),
            level[1].as_slice(),
            level[2].as_slice(),
        );
        state.par_iter_mut().enumerate().for_each(|(k, s)| {
            if !s.2 {
                return;
            }
            let reach = if c_l[k] > 1e-6 && y_l[k] > 0.0 {
                zone(((y_l[k] / c_l[k]) as f64).powf(-1.0 / p))
            } else {
                f64::INFINITY
            };
            let allowed = reach.min(own[k]);
            if allowed >= here {
                s.0 = c_l[k];
                s.1 = t_l[k];
            } else {
                let f = ((allowed - below) / (here - below)).clamp(0.0, 1.0) as f32;
                s.0 += f * (c_l[k] - s.0);
                s.1 += f * (t_l[k] - s.1);
                s.2 = false;
            }
        });
        below = here;
    }
    state.into_iter().map(|(c, t, _)| (c, t)).unzip()
}

/// The luminance of a regolith mixture when fresh and when mature, for the body's own
/// albedo (F52).
pub fn mixture_luminance(map: &BTreeMap<String, f64>) -> (f64, f64) {
    let mix = resolve_mix(map);
    let lum = |c: [f32; 3]| (c[0] + c[1] + c[2]) as f64 / 3.0;
    (
        lum(mixture_colour(&mix, 0.0)),
        lum(mixture_colour(&mix, 1.0)),
    )
}

/// The reflectance a surface has when fresh and when mature, as the lag law needs them:
/// the body's regolith and, for ice renewed from below, clean ice.
pub struct SurfaceAlbedos {
    pub fresh: f64,
    pub mature: f64,
    pub clean_fresh: f64,
    pub clean_mature: f64,
}

pub fn surface_albedos(graph: &ParamGraph) -> SurfaceAlbedos {
    let luminance = |c: [f32; 3]| (c[0] + c[1] + c[2]) as f64 / 3.0;
    let bedrock =
        end_member(graph.str("materials.bedrock")).unwrap_or_else(|| end_member("basalt").unwrap());
    let regolith = graph
        .get("materials.regolith_mix")
        .as_map()
        .map(resolve_mix)
        .unwrap_or_else(|| vec![(bedrock, 1.0)]);
    let ice = end_member("water_ice").unwrap();
    SurfaceAlbedos {
        fresh: luminance(bedrock.fresh),
        mature: luminance(mixture_colour(&regolith, 1.0)),
        clean_fresh: luminance(ice.fresh),
        clean_mature: luminance(ice.mature),
    }
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
    let comp_f = layers.composition.as_slice();
    let stain_f = layers.stain.as_slice();
    let salted = |u: u16| matches!(u, UNIT_LINEA | UNIT_BAND | UNIT_CHAOS | UNIT_CRYODOME);
    // A flow unit's reflectance departs from the mean by its eruption's composition (F55),
    // which the unit's own ejecta keeps.
    let tint = |k: usize| {
        if salted(units[k]) || units[k] == UNIT_CONTINENTAL || units[k] == UNIT_OCEANIC {
            1.0
        } else {
            (1.0 + LAVA_COMPOSITION_RANGE as f32 * comp_f[k]).max(0.1)
        }
    };
    let tinted: Vec<f32> = (0..grid.cell_count())
        .into_par_iter()
        .map(|k| cover_f[k] * tint(k))
        .collect();
    let (mixed, mixed_tinted) = if gardened {
        // Each lava cell's gardening since its lava erupted, relative to the lunar maria's
        // at the reference age (F95).
        let flux = crate::crater::FluxHistory {
            steady: crate::crater::ProductionFunction::by_name(
                graph.str("cratering.production_function"),
            )
            .chronology
            .2,
            early_ratio: graph.f64("cratering.early_flux_ratio"),
            early_decay: graph.f64("cratering.early_flux_decay"),
            record_start: graph.f64("surface.solidification_time"),
        };
        let lunar = crate::crater::FluxHistory::from_chronology(
            crate::crater::ProductionFunction::lunar().chronology,
            0.0,
        );
        let reference = DIFFUSIVITY_LUNAR
            * lunar.equivalent_gyr(
                NEUKUM_REFERENCE_AGE_GYR - REGOLITH_MIXING_REFERENCE_GYR,
                NEUKUM_REFERENCE_AGE_GYR,
            );
        let kappa = graph.f64("cratering.diffusivity");
        let lava_age = layers.lava_age.as_slice();
        let dose: Vec<Option<f32>> = (0..grid.cell_count())
            .into_par_iter()
            .map(|k| {
                (cover_f[k] > 0.0).then(|| {
                    let erupted = if lava_age[k] > 0.0 {
                        lava_age[k]
                    } else {
                        ages_f[k]
                    };
                    (kappa * flux.equivalent_gyr(erupted as f64, present) / reference) as f32
                })
            })
            .collect();
        let clock = std::time::Instant::now();
        let mixed = garden_mix(
            grid,
            grid.delta() * graph.f64("body.radius"),
            cover_f,
            &tinted,
            &dose,
        );
        if std::env::var_os("EXOGEN_DEBUG").is_some() {
            eprintln!("regolith mixing: {:.1} s", clock.elapsed().as_secs_f64());
        }
        mixed
    } else {
        (cover_f.to_vec(), tinted.clone())
    };
    let r = rays.as_slice();
    let bedrock_f = layers.bedrock.as_slice();
    let io_cell = grid.delta() * graph.f64("body.radius");
    let data: Vec<[f32; 3]> = (0..grid.cell_count())
        .into_par_iter()
        .map(|k| {
            if io {
                // What the plumes have left on the cell (F84): sulfur, altered by the
                // plasma where the atmosphere does not shield it and red where a plume's
                // sulfur is fresh; dark silicate, bare lava or pyroclasts; and frost
                // over all.
                let sulfur = lerp3(
                    lerp3(IO_SULFUR, IO_POLAR, stain_f[k].clamp(0.0, 1.0)),
                    IO_PLUME_RED,
                    melt_f[k].clamp(0.0, 1.0),
                );
                let ground = lerp3(sulfur, IO_FLOW_DARK, deep_f[k].clamp(0.0, 1.0));
                let frost = cover_f[k].clamp(0.0, 1.0);
                return match units[k] {
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
                        lerp3(lerp3(ground, IO_FROST, frost), IO_MOUNTAIN_ROCK, bare)
                    }
                    _ => lerp3(ground, IO_FROST, frost),
                };
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
                let (fresh, mature) = if c > 0.0
                    && !salted(units[k])
                    && units[k] != UNIT_CONTINENTAL
                    && units[k] != UNIT_OCEANIC
                {
                    // The lava's tint, mixed sideways with its coverage (F95).
                    let t = if mixed[k] > 1e-4 {
                        mixed_tinted[k] / mixed[k]
                    } else {
                        tint(k)
                    };
                    let tinted = |v: [f32; 3]| [v[0] * t, v[1] * t, v[2] * t];
                    (
                        lerp3(fresh, tinted(lava.fresh), c),
                        lerp3(mature, tinted(lava.mature), c),
                    )
                } else {
                    (fresh, mature)
                };
                // A ray freshens the ground its secondaries dug up, not the melt of an
                // older crater it crosses, which keeps the maturity of its own age (F91).
                let effective = (m[k] * (1.0 - r[k])).clamp(0.0, 1.0);
                let colour = lerp3(
                    lerp3(fresh, mature, effective),
                    lerp3(melt_pair.fresh, mature, m[k]),
                    melt_f[k],
                );
                // A crack's margins, stained red-brown by the non-ice material it brought
                // up, under whatever frost lies on them (F78).
                let colour = if stain_f[k] > 0.0 {
                    lerp3(colour, salt_colour, stain_f[k].clamp(0.0, 1.0))
                } else {
                    colour
                };
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
                        // The trailing hemisphere's plasma darkens and reddens the ice
                        // in the time it takes to process the optical surface (F58); a
                        // plume's fallout that buries that surface first hides the
                        // darkening by the ratio of the two times.
                        if trailing > 0.0 {
                            let tau = TRAILING_DARKENING_GYR / trailing;
                            let exposed = if burial.is_finite() {
                                burial / (burial + tau)
                            } else {
                                1.0
                            };
                            let w = (trailing
                                * exposed
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::grid::Cell;

    #[test]
    fn gardening_blurs_an_old_contact_and_keeps_a_young_one_sharp() {
        // F95: two lava units of different tints meet along a line. At the reference
        // gardening the contact blurs over REGOLITH_MIXING_LENGTH; where one unit has had
        // a hundredth of it, the contact stays sharp on both sides, the old unit mixing
        // with the young only as long as the young has lain there.
        let grid = CubeSphere::new(128);
        let body_radius = 100e3;
        let cell = grid.delta() * body_radius;
        let n = grid.cell_count();
        let side = |k: usize| {
            let c = grid.cell(k);
            c.face == 0 && c.i >= 64
        };
        let cover = vec![1.0f32; n];
        let tinted: Vec<f32> = (0..n).map(|k| if side(k) { 0.8 } else { 1.2 }).collect();
        let at = |v: &[f32], i: u32| v[grid.index(Cell { face: 0, i, j: 64 })];
        let old = vec![Some(1.0f32); n];
        let (_, mixed) = garden_mix(grid, cell, &cover, &tinted, &old);
        assert!(
            at(&mixed, 61) < 1.15 && at(&mixed, 66) > 0.85,
            "an old contact blurs: {} and {}",
            at(&mixed, 61),
            at(&mixed, 66)
        );
        let young: Vec<Option<f32>> = (0..n)
            .map(|k| Some(if side(k) { 0.01 } else { 1.0 }))
            .collect();
        let (_, mixed) = garden_mix(grid, cell, &cover, &tinted, &young);
        assert!(
            at(&mixed, 61) > 1.19 && at(&mixed, 65) < 0.85,
            "a young contact stays sharp: {} and {}",
            at(&mixed, 61),
            at(&mixed, 65)
        );
    }
}
