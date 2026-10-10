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

use crate::fields::{Field, LayerSet};
use crate::params::ParamGraph;
use crate::units::{UNIT_CONTINENTAL, UNIT_OCEANIC};

pub const VERSION: u64 = 3;

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
        mature: [0.55, 0.56, 0.58],
        roughness: 0.55,
    },
    EndMember {
        id: "water_ice_dirty",
        fresh: [0.36, 0.35, 0.34],
        mature: [0.26, 0.25, 0.24],
        roughness: 0.70,
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
];

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

/// Linear-reflectance albedo of the surface.
pub fn albedo(
    graph: &ParamGraph,
    layers: &LayerSet,
    maturity: &Field<f32>,
    rays: &Field<f32>,
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
    let (deep, shallow): ([f32; 3], [f32; 3]) = match graph.str("volatiles.working_liquid") {
        "water" | "ammonia_water" => ([0.012, 0.035, 0.085], [0.06, 0.22, 0.26]),
        "methane_ethane" => ([0.02, 0.02, 0.015], [0.07, 0.07, 0.05]),
        _ => ([0.03, 0.03, 0.03], [0.08, 0.08, 0.08]),
    };
    let grid = layers.grid;
    let liquid = layers.liquid.as_slice();
    let units = layers.unit_id.as_slice();
    let m = maturity.as_slice();
    let r = rays.as_slice();
    let data: Vec<[f32; 3]> = (0..grid.cell_count())
        .into_par_iter()
        .map(|k| {
            if liquid[k] > 0.0 {
                let t = (-liquid[k] / 60.0).exp();
                lerp3(deep, shallow, t)
            } else {
                let pair = match units[k] {
                    UNIT_CONTINENTAL => continental,
                    UNIT_OCEANIC => oceanic,
                    _ => base,
                };
                let effective = (m[k] * (1.0 - r[k])).clamp(0.0, 1.0);
                lerp3(pair.fresh, pair.mature, effective)
            }
        })
        .collect();
    let mut f = Field::new(grid);
    f.as_mut_slice().copy_from_slice(&data);
    f
}
