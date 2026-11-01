// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! The caps the inventory lays in its cold traps (PLAN.md Phase 5B, F241).
//!
//! Each ice of the inventory fills its traps in their order (F240), from the coldest, each
//! place's trap read at its own height. Its surface stands as the plastic profile of an ice
//! that holds a yield stress (Nye 1951), rising as the square root of the distance in from
//! the margin of the ground it covers, so that the ice stands in domes on open ground and
//! fills a basin toward its rim; the ground it covers grows until its caps hold the
//! inventory. An ice that is the body's working liquid stands only where it never melts and
//! leaves the rest to the sea.

use std::cmp::Ordering;
use std::collections::BinaryHeap;

use rayon::prelude::*;

use crate::climate::{Ice, Traps};
use crate::constants::*;
use crate::fields::{Field, LayerSet};
use crate::geom::angle_between;
use crate::grid::CubeSphere;
use crate::grid::projection::{PaddedField, sample_bilinear};
use crate::params::ParamGraph;

/// A cap the inventory laid (F241): its ice, its thickness on the finish's grid, m, the ice
/// it holds, m of the body's global layer of its liquid, and of that what the sea gave.
#[derive(Clone, Debug)]
pub struct Cap {
    pub ice: Ice,
    pub thickness: Field<f32>,
    pub laid: f64,
    pub from_sea: f64,
}

impl Cap {
    /// The share of the ground's colour the cap hides at a cell, its mass hiding it as the
    /// frost's does (F239).
    pub fn cover(&self, k: usize) -> f64 {
        let mass = self.thickness.as_slice()[k] as f64 * ice_density(self.ice.species);
        crate::climate::frost_cover(mass as f32)
    }
}

/// An ice's density, kg/m³.
pub fn ice_density(species: &str) -> f64 {
    ICE_DENSITIES
        .iter()
        .find(|(s, _)| *s == species)
        .map_or(ICE_DENSITIES[0].1, |(_, rho)| *rho)
}

/// The species a working liquid is made of.
pub fn liquid_species(working: &str) -> &'static [&'static str] {
    match working {
        "water" | "ammonia_water" => &["H2O"],
        "methane_ethane" => &["CH4", "C2H6"],
        "nitrogen" => &["N2"],
        "sulfur_dioxide" => &["SO2"],
        "carbon_dioxide" => &["CO2"],
        _ => &[],
    }
}

/// Lays each ice of the inventory in its traps (F241), the least volatile first, each on the
/// ground the last left: the ice covers the places from its coldest trap on, a sea's floor
/// holding none; its surface stands as the plastic profile of an ice that holds
/// CAP_YIELD_STRESS (Nye 1951) from the margin of the ground it covers, until its caps hold
/// its inventory or, a liquid's, until it covers all the ground it can. Each place's trap is
/// its region's on the climate's grid, read smoothly between its cells: a cap kilometres
/// thick buries the relief under it, which decides where it stands no more. An ice's
/// coldest traps are those where it needs the least vapour in the air (F240), the frost's
/// own gas's those where its ice grows the most (F242); an ice that is the body's working
/// liquid stands only where it never melts, the frost's own gas's where its ice grows.
pub fn lay(graph: &ParamGraph, layers: &mut LayerSet, traps: &Traps) -> Vec<Cap> {
    let radius = graph.f64("body.radius");
    let gravity = graph.f64("body.gravity");
    let empty = std::collections::BTreeMap::new();
    let inventory = graph
        .get("volatiles.inventory")
        .as_map()
        .unwrap_or(&empty)
        .clone();
    let liquid = liquid_species(graph.str("volatiles.working_liquid"));
    // Each ice's traps on the climate's grid, the lower the colder, and where a liquid's ice
    // may stand, at or above none.
    let mut ices: Vec<(Ice, Vec<f32>, Vec<f32>)> = traps
        .vapour
        .iter()
        .zip(&traps.warmest)
        .map(|((ice, vapour), warmest)| {
            let rank = vapour
                .as_slice()
                .iter()
                .map(|&v| (v as f64).max(1e-300).ln() as f32)
                .collect();
            let melt = ice.volatile.t_triple as f32;
            let stands = warmest.as_slice().iter().map(|&w| melt - w).collect();
            (*ice, rank, stands)
        })
        .collect();
    if let Some((ice, growth)) = &traps.growth {
        let rank = growth.as_slice().iter().map(|&g| -g).collect();
        ices.push((*ice, rank, growth.as_slice().to_vec()));
    }
    // The least volatile first, its triple point the warmest.
    ices.sort_by(|a, b| b.0.volatile.t_triple.total_cmp(&a.0.volatile.t_triple));
    let grid = layers.grid;
    let whole = 4.0 * std::f64::consts::PI * radius * radius;
    let areas: Vec<f64> = (0..grid.cell_count())
        .map(|k| grid.cell_area(grid.cell(k)) * radius * radius)
        .collect();
    let coarse = traps.grid();
    let mut caps = Vec::new();
    for (ice, rank, stands) in ices {
        let gel = inventory.get(ice.species).copied().unwrap_or(0.0);
        if gel <= 0.0 {
            continue;
        }
        let rho = ice_density(ice.species);
        let wanted = gel * ice.volatile.liquid_density / rho * whole;
        let melts = liquid.contains(&ice.species);
        let ground = layers.height();
        // Each place's region's trap, read smoothly on the climate's grid; where its ice may
        // stand, coldest first.
        let trap = PaddedField::new(&coarse, &rank);
        let stands = PaddedField::new(&coarse, &stands);
        let mut ranks: Vec<(usize, f64)> = (0..grid.cell_count())
            .into_par_iter()
            .filter_map(|k| {
                if layers.liquid.as_slice()[k] > 0.0 {
                    return None;
                }
                let d = grid.center(grid.cell(k));
                if melts && sample_bilinear(&coarse, &stands, d) < 0.0 {
                    return None;
                }
                let rank = sample_bilinear(&coarse, &trap, d) as f64;
                rank.is_finite().then_some((k, rank))
            })
            .collect();
        ranks.sort_by(|a, b| a.1.total_cmp(&b.1).then(a.0.cmp(&b.0)));
        if ranks.is_empty() {
            continue;
        }
        let rise = (2.0 * CAP_YIELD_STRESS / (rho * gravity)).sqrt();
        let held = |t: &[f32]| -> f64 { t.iter().zip(&areas).map(|(&t, a)| t as f64 * a).sum() };
        let shape = |n: usize| profile(grid, &ground, &ranks[..n], rise, radius, wanted);
        // The fewest coldest places whose caps hold the inventory: a liquid's ice that covers
        // all the ground it may holds what it can, and ice that may cover every place
        // spreads over all of them only if no fewer hold it.
        let everywhere = ranks.len() == grid.cell_count();
        let all = shape(ranks.len());
        let (mut thickness, n) = if (!everywhere && held(&all) <= wanted)
            || (everywhere && held(&shape(ranks.len() - 1)) < wanted)
        {
            (all, ranks.len())
        } else {
            let (mut lo, mut hi) = (0, ranks.len());
            while hi - lo > 1 {
                let mid = (lo + hi) / 2;
                if held(&shape(mid)) >= wanted {
                    hi = mid;
                } else {
                    lo = mid;
                }
            }
            (shape(hi), hi)
        };
        let volume = held(&thickness);
        if volume <= 0.0 {
            continue;
        }
        // Just the inventory: the profile's shape, its thickness scaled to it.
        let kept = wanted.min(volume);
        let scale = (kept / volume) as f32;
        for t in thickness.iter_mut() {
            *t *= scale;
        }
        for (ice_layer, t) in layers.ice.as_mut_slice().iter_mut().zip(&thickness) {
            *ice_layer += t;
        }
        let laid = kept * rho / ice.volatile.liquid_density / whole;
        if std::env::var_os("EXOGEN_DEBUG").is_some() {
            let thickest = thickness.iter().cloned().fold(0.0f32, f32::max);
            eprintln!(
                "caps: {} lays {laid:.3} m of its {gel:.3} m over {n} places, {:.2}% of the body, at most {thickest:.0} m thick",
                ice.species,
                100.0 * ranks[..n].iter().map(|r| areas[r.0]).sum::<f64>() / whole
            );
        }
        let mut field = Field::new(grid);
        field.as_mut_slice().copy_from_slice(&thickness);
        caps.push(Cap {
            ice,
            thickness: field,
            laid,
            from_sea: if melts { laid } else { 0.0 },
        });
    }
    caps
}

/// A place on the margin's front, by its surface: the lowest first.
#[derive(PartialEq)]
struct Front(f64, usize);

impl Eq for Front {}

impl PartialOrd for Front {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Front {
    fn cmp(&self, other: &Self) -> Ordering {
        other.0.total_cmp(&self.0).then(other.1.cmp(&self.1))
    }
}

/// The thickness of the caps an ice covering `cells` stands in (F241), m: its surface the
/// lowest of the plastic profiles from every point of the margin of the ground it covers,
/// the margin's ground plus `rise` times the square root of the distance in from it, found
/// outward from the margin as distances are; where no ground is left uncovered, `whole`
/// m³ spread evenly.
fn profile(
    grid: CubeSphere,
    ground: &Field<f32>,
    cells: &[(usize, f64)],
    rise: f64,
    radius: f64,
    whole: f64,
) -> Vec<f32> {
    let n = grid.cell_count();
    let mut thickness = vec![0.0f32; n];
    if cells.len() == n {
        let layer = whole / (4.0 * std::f64::consts::PI * radius * radius);
        thickness.fill(layer as f32);
        return thickness;
    }
    let mut inside = vec![false; n];
    for &(k, _) in cells {
        inside[k] = true;
    }
    let level = |k: usize| {
        let h = ground.as_slice()[k] as f64;
        if h.is_finite() { h } else { 0.0 }
    };
    // Each covered place's surface, the margin's ground it rises from and its distance in.
    let mut best: Vec<(f64, f64, f64)> = vec![(f64::INFINITY, 0.0, 0.0); n];
    let mut front = BinaryHeap::new();
    let steps = [(1, 0), (-1, 0), (0, 1), (0, -1)];
    for &(k, _) in cells {
        let c = grid.cell(k);
        let here = grid.center(c);
        for (di, dj) in steps {
            let nb = grid.neighbor(c, di, dj);
            let j = grid.index(nb);
            if inside[j] {
                continue;
            }
            // The margin lies half way to the uncovered ground beside it.
            let d = 0.5 * angle_between(here, grid.center(nb)) * radius;
            let surface = level(j) + rise * d.sqrt();
            if surface < best[k].0 {
                best[k] = (surface, level(j), d);
            }
        }
        if best[k].0.is_finite() {
            front.push(Front(best[k].0, k));
        }
    }
    while let Some(Front(surface, k)) = front.pop() {
        if surface > best[k].0 {
            continue;
        }
        let (_, from, d) = best[k];
        let c = grid.cell(k);
        let here = grid.center(c);
        for (di, dj) in steps {
            let nb = grid.neighbor(c, di, dj);
            let j = grid.index(nb);
            if !inside[j] {
                continue;
            }
            let dj = d + angle_between(here, grid.center(nb)) * radius;
            let s = from + rise * dj.sqrt();
            if s < best[j].0 {
                best[j] = (s, from, dj);
                front.push(Front(s, j));
            }
        }
    }
    for &(k, _) in cells {
        thickness[k] = (best[k].0 - level(k)).max(0.0) as f32;
    }
    thickness
}

/// The ground's colour under the caps (F241), each in the order laid hiding what is under
/// it by its cover.
pub fn over(albedo: &Field<[f32; 3]>, caps: &[Cap]) -> Field<[f32; 3]> {
    if caps.is_empty() {
        return albedo.clone();
    }
    let grid = albedo.grid();
    Field::from_fn(grid, |c, _| {
        let k = grid.index(c);
        let mut a = albedo.get(c);
        for cap in caps {
            let cover = cap.cover(k) as f32;
            for (v, &ice) in a.iter_mut().zip(&cap.ice.colour) {
                *v += (ice - *v) * cover;
            }
        }
        a
    })
}

/// The caps on the climate's grid (F241): the share of each cell they cover and their
/// albedo there under the star's light (F250), the uppermost's showing.
pub fn on_climate(
    caps: &[Cap],
    climate: CubeSphere,
    balance: &crate::climate::EnergyBalance,
) -> (Field<f32>, Field<f32>) {
    let grid = caps[0].thickness.grid();
    let cover = Field::from_fn(grid, |c, _| {
        let k = grid.index(c);
        let mut hidden = 0.0;
        for cap in caps {
            let s = cap.cover(k);
            hidden += (1.0 - hidden) * s;
        }
        hidden as f32
    });
    let reflected = Field::from_fn(grid, |c, _| {
        let k = grid.index(c);
        let (mut shown, mut sum) = (1.0, 0.0);
        for cap in caps.iter().rev() {
            let s = cap.cover(k);
            sum += shown * s * balance.ice_albedo(&cap.ice);
            shown *= 1.0 - s;
        }
        sum as f32
    });
    let share = crate::climate::coarsen(&cover, climate);
    let reflected = crate::climate::coarsen(&reflected, climate);
    let albedo = Field::from_fn(climate, |c, _| {
        let s = share.get(c);
        if s > 0.0 { reflected.get(c) / s } else { 0.0 }
    });
    (share, albedo)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::climate::{Downscale, EnergyBalance, Ground};
    use crate::grid::lon_lat_from_direction;

    /// A preset's graph and its cold traps on ground at the datum.
    fn traps_of(name: &str) -> (ParamGraph, Traps) {
        let graph = crate::presets::preset(name).unwrap().to_graph().0;
        let eb = EnergyBalance::from_graph(&graph);
        let climate = eb.solve(CubeSphere::new(CLIMATE_FACE_SIZE));
        let ground = Ground::from_graph(&graph);
        let cycle = eb.cycle(&climate, &ground, None);
        let ices = Ice::of_inventory(&graph, [0.80, 0.82, 0.86]);
        let heights = Field::new(climate.temperature.grid());
        let downscale = Downscale::from_graph(&graph);
        let traps = eb.traps(&climate, &ground, &cycle, &ices, &heights, &downscale);
        (graph, traps)
    }

    /// A cap's share of the sphere, its thickest and the lowest latitude it reaches, deg.
    fn reach(cap: &Cap) -> (f64, f64, f64) {
        let grid = cap.thickness.grid();
        let (mut covered, mut thickest, mut lowest) = (0.0, 0.0f64, 90.0f64);
        for (k, &t) in cap.thickness.as_slice().iter().enumerate() {
            if t > 0.0 {
                let c = grid.cell(k);
                covered += grid.cell_area(c);
                thickest = thickest.max(t as f64);
                let (_, lat) = lon_lat_from_direction(grid.center(c));
                lowest = lowest.min(lat.abs().to_degrees());
            }
        }
        (covered / (4.0 * std::f64::consts::PI), thickest, lowest)
    }

    /// F241: Mars's 25 m of water lays all of itself in caps at its poles, a few kilometres
    /// thick at most over a few percent of the body, as its polar layered deposits stand:
    /// some 3 km thick over some 2.6 percent of it (Zuber et al. 1998; Plaut et al. 2007).
    #[test]
    fn mars_s_water_lays_its_polar_caps() {
        let (graph, traps) = traps_of("Mars");
        let mut layers = LayerSet::new(CubeSphere::new(64));
        let caps = lay(&graph, &mut layers, &traps);
        let cap = caps
            .iter()
            .find(|c| c.ice.species == "H2O")
            .expect("water caps");
        assert!((cap.laid - 25.0).abs() < 1e-6, "laid {}", cap.laid);
        let (share, thickest, lowest) = reach(cap);
        assert!((0.01..0.05).contains(&share), "share {share}");
        assert!((2000.0..6000.0).contains(&thickest), "thickest {thickest}");
        assert!(lowest > 60.0, "reaches {lowest}");
        let ice: f64 = layers.ice.as_slice().iter().map(|&t| t as f64).sum();
        assert!(ice > 0.0 && cap.from_sea == 0.0);
    }

    /// F241: an ice spreads over every place only when no fewer hold it: Pluto's 5 m of
    /// methane gathers in caps at its coldest traps rather than lying over the whole body.
    #[test]
    fn an_ice_gathers_in_its_coldest_traps() {
        let (graph, traps) = traps_of("Pluto");
        let mut layers = LayerSet::new(CubeSphere::new(64));
        let caps = lay(&graph, &mut layers, &traps);
        let cap = caps
            .iter()
            .find(|c| c.ice.species == "CH4")
            .expect("methane caps");
        let (share, thickest, _) = reach(cap);
        assert!(
            share < 0.2 && thickest > 100.0,
            "share {share}, thickest {thickest}"
        );
    }

    /// F242: the frost's own gas lays its caps where its ice grows the most, last, on the
    /// caps of the ices less volatile: Mars's carbon dioxide on ground at the datum at its
    /// poles, its whole inventory, the 0.12 m of its residual cap and the deposits buried
    /// under its south polar layered deposits (Phillips et al. 2011).
    #[test]
    fn the_frost_s_own_gas_lays_its_caps_where_its_ice_grows() {
        let (graph, traps) = traps_of("Mars");
        let mut layers = LayerSet::new(CubeSphere::new(64));
        let caps = lay(&graph, &mut layers, &traps);
        let species: Vec<&str> = caps.iter().map(|c| c.ice.species).collect();
        assert_eq!(species, ["H2O", "CO2"]);
        let cap = &caps[1];
        assert!((cap.laid - 0.12).abs() < 1e-9, "laid {}", cap.laid);
        let (share, thickest, lowest) = reach(cap);
        assert!(
            share > 0.0 && thickest > 0.0 && lowest > 75.0,
            "{share} {thickest} {lowest}"
        );
    }

    /// F242: a working liquid that is the frost's own gas lays its ice only where the ice
    /// grows: Titan's methane, its air over 90 K holding its ice above its point, lays none.
    #[test]
    fn titan_s_methane_lays_no_cap() {
        let (graph, traps) = traps_of("Titan");
        let mut layers = LayerSet::new(CubeSphere::new(32));
        let caps = lay(&graph, &mut layers, &traps);
        assert!(
            caps.is_empty(),
            "{:?}",
            caps.iter().map(|c| c.ice.species).collect::<Vec<_>>()
        );
    }

    /// F241: ice covering a disc on flat ground stands as Nye's plastic dome, its centre
    /// `rise` times the square root of the disc's radius above the ground.
    #[test]
    fn a_disc_of_ice_stands_as_a_plastic_dome() {
        let grid = CubeSphere::new(128);
        let ground = Field::new(grid);
        let radius = 1.0e6;
        let cells: Vec<(usize, f64)> = (0..grid.cell_count())
            .filter(|&k| {
                let (_, lat) = lon_lat_from_direction(grid.center(grid.cell(k)));
                lat > 80f64.to_radians()
            })
            .map(|k| (k, 0.0))
            .collect();
        let rise = 5.0;
        let t = profile(grid, &ground, &cells, rise, radius, 0.0);
        let pole = grid.index(grid.locate([0.0, 0.0, 1.0]).0);
        let expect = rise * (10f64.to_radians() * radius).sqrt();
        let centre = t[pole] as f64;
        assert!(
            (centre / expect - 1.0).abs() < 0.03,
            "{centre} against {expect}"
        );
    }
}
