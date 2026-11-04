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

use crate::climate::{Ice, Traps, read_table};
use crate::constants::*;
use crate::fields::{Field, LayerSet};
use crate::geom::{angle_between, normalize};
use crate::grid::projection::{PaddedField, sample_bilinear};
use crate::grid::{CubeSphere, lon_lat_from_direction};
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
/// liquid stands only where it never melts, the frost's own gas's where its ice grows. An ice
/// `gathered` through the body's age (F259) covers the places it gathered the most at first,
/// its masses on the climate's grid, then its traps' coldest.
pub fn lay(
    graph: &ParamGraph,
    layers: &mut LayerSet,
    traps: &Traps,
    gathered: &[(&'static str, Field<f32>)],
) -> Vec<Cap> {
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
    // An ice gathered covers the places it gathered at, the most first, then, where its caps
    // need more ground than they gathered on, the coldest of its traps.
    for (ice, rank, _) in ices.iter_mut() {
        if let Some((_, mass)) = gathered.iter().find(|(s, _)| *s == ice.species) {
            let least = rank
                .iter()
                .filter(|r| r.is_finite())
                .fold(f32::MAX, |a, &r| a.min(r));
            *rank = mass
                .as_slice()
                .iter()
                .zip(rank.iter())
                .map(|(&m, &r)| if m > 0.0 { -m } else { r - least })
                .collect();
        }
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
/// outward from the margin as distances are, each place carrying the point of the margin it
/// rises from and its distance the true one from it (F280): summed along the grid's steps,
/// the distance ran long across them, its domes diamonds, its surface ruled along the grid
/// and stepped where one margin's point gave way to another's. Where no ground is left
/// uncovered, `whole` m³ spread evenly.
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
    // Each covered place's surface, the margin's ground it rises from and that margin's point.
    let mut best: Vec<(f64, f64, [f64; 3])> = vec![(f64::INFINITY, 0.0, [0.0; 3]); n];
    let mut front = BinaryHeap::new();
    let steps = [(1, 0), (-1, 0), (0, 1), (0, -1)];
    let around = [
        (1, 0),
        (-1, 0),
        (0, 1),
        (0, -1),
        (1, 1),
        (1, -1),
        (-1, 1),
        (-1, -1),
    ];
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
            let there = grid.center(nb);
            let point = normalize([here[0] + there[0], here[1] + there[1], here[2] + there[2]]);
            let d = angle_between(here, point) * radius;
            let surface = level(j) + rise * d.sqrt();
            if surface < best[k].0 {
                best[k] = (surface, level(j), point);
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
        let (_, from, point) = best[k];
        let c = grid.cell(k);
        for (di, dj) in around {
            let nb = grid.neighbor(c, di, dj);
            let j = grid.index(nb);
            if !inside[j] {
                continue;
            }
            let d = angle_between(point, grid.center(nb)) * radius;
            let s = from + rise * d.sqrt();
            if s < best[j].0 {
                best[j] = (s, from, point);
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
    shown(&caps.iter().collect::<Vec<_>>(), climate, balance)
}

/// The frost's own gas's caps on the climate's grid as the reservoir its air draws on
/// (F258): the share of each cell they cover, by the cover their ice hides the ground with,
/// the ice they hold over that share, kg/m², and the albedo of the ground about them under
/// the other ices' caps; none where the gas laid no cap.
pub fn reservoir(
    caps: &[Cap],
    climate: CubeSphere,
    balance: &crate::climate::EnergyBalance,
    species: &str,
) -> Option<crate::climate::Reservoir> {
    let cap = caps.iter().find(|c| c.ice.species == species)?;
    let grid = cap.thickness.grid();
    let rho = ice_density(species);
    let cover = Field::from_fn(grid, |c, _| cap.cover(grid.index(c)) as f32);
    let mass = Field::from_fn(grid, |c, _| (cap.thickness.get(c) as f64 * rho) as f32);
    let share = crate::climate::coarsen(&cover, climate);
    let held = crate::climate::coarsen(&mass, climate);
    let mass = Field::from_fn(climate, |c, _| {
        let s = share.get(c);
        if s > 0.0 { held.get(c) / s } else { 0.0 }
    });
    let others: Vec<&Cap> = caps.iter().filter(|c| c.ice.species != species).collect();
    let over = (!others.is_empty()).then(|| shown(&others, climate, balance));
    Some(crate::climate::Reservoir {
        share,
        mass,
        ground: balance.ground_albedo(climate, over.as_ref()),
    })
}

/// Caps on the climate's grid, as `on_climate`.
fn shown(
    caps: &[&Cap],
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

/// An ice's traps on the cycle's latitudes' rows as the seasons turn (F259): at each latitude
/// of the rows and each level of the table, latitude after latitude, at the present and their
/// mean through the turning. The frost's own gas's ice's growth, kg/m² a year; every other
/// ice's year-mean saturation pressure over it, Pa.
#[derive(Clone, Debug)]
pub struct Cycled {
    pub ice: Ice,
    pub present: Vec<f64>,
    pub mean: Vec<f64>,
}

/// What the ices gather under (F259): the traps' table's levels, shares of the datum's
/// pressure; the pressure's scale height, m, none on an airless body; the period over which
/// the seasons turn, s, none where they hold; the body's year and age, s; its mean surface
/// temperature, K; and whether it has air.
#[derive(Clone, Copy, Debug)]
pub struct Gathering<'a> {
    pub levels: &'a [f64],
    pub scale_height: f64,
    pub cycle: f64,
    pub year: f64,
    pub age: f64,
    pub temperature: f64,
    pub airless: bool,
}

/// Each ice of the inventory gathered through the body's age (F259), from a layer over the
/// whole body, on the climate's grid of `heights`, the ground's above the datum, m: what each
/// place holds at the end, kg/m².
///
/// An ice exchanges where it lies uppermost, the most volatile on top, and lays itself where
/// it would grow, at the rate its traps give at its surface's height, the ground's and the
/// ices' over it. The frost's own gas's ice gains its growth less what holds its gas's whole
/// mass, the air buffered by its caps (F258); every other ice sublimes into the air or takes
/// from it by the saturation pressure over it against the vapour the air holds there, at the
/// air's turbulent exchange, C_H U / (R_s T) per pascal, or into the vacuum at the
/// Hertz-Knudsen rate on an airless body, the air's vapour the one at which what all of it
/// gives and takes balances (Ingersoll 1970).
///
/// As the seasons turn, an ice follows its traps' departure from their mean through the
/// turning by the share 1 / sqrt(1 + (2 pi tau / P)^2) a store of exchange time tau answers a
/// cycle of period P by (the lag left out), its exchange time the time its layer over the
/// body takes to move at its rates there: the frost's own gas, moving in years, the present
/// season's; a slow ice, moving in millions of years, the turning's mean, so that a pole the
/// present favours does not take what the other gathered.
pub fn gather(
    graph: &ParamGraph,
    heights: &Field<f32>,
    traps: &[Cycled],
    frost: Option<&'static str>,
    g: &Gathering,
) -> Vec<(&'static str, Field<f32>)> {
    let grid = heights.grid();
    let n = grid.cell_count();
    let empty = std::collections::BTreeMap::new();
    let inventory = graph
        .get("volatiles.inventory")
        .as_map()
        .unwrap_or(&empty)
        .clone();
    let area: Vec<f64> = (0..n).map(|k| grid.cell_area(grid.cell(k))).collect();
    let whole: f64 = area.iter().sum();
    let lat: Vec<f64> = (0..n)
        .map(|k| lon_lat_from_direction(grid.center(grid.cell(k))).1)
        .collect();
    let ground: Vec<f64> = heights
        .as_slice()
        .iter()
        .map(|&h| if h.is_finite() { h as f64 } else { 0.0 })
        .collect();
    let thin = |h: f64| {
        if g.scale_height > 0.0 {
            (-h / g.scale_height).exp()
        } else {
            1.0
        }
    };
    let count = traps.len();
    let own: Vec<bool> = traps.iter().map(|t| Some(t.ice.species) == frost).collect();
    let rho: Vec<f64> = traps.iter().map(|t| ice_density(t.ice.species)).collect();
    // Each ice's exchange per pascal of saturation pressure over the air's, kg/m²/s/Pa.
    let exchange: Vec<f64> = traps
        .iter()
        .map(|t| {
            let rt = t.ice.volatile.r_specific * g.temperature.max(1.0);
            if g.airless {
                (1.0 / (std::f64::consts::TAU * rt)).sqrt()
            } else {
                AIR_EXCHANGE_COEFFICIENT * SURFACE_WIND_SPEED / rt
            }
        })
        .collect();
    // The most volatile uppermost.
    let mut order: Vec<usize> = (0..count).collect();
    order.sort_by(|&a, &b| {
        let (ta, tb) = (
            traps[a].ice.volatile.t_triple,
            traps[b].ice.volatile.t_triple,
        );
        ta.total_cmp(&tb)
    });
    let mut mass: Vec<Vec<f64>> = traps
        .iter()
        .map(|t| {
            let gel = inventory.get(t.ice.species).copied().unwrap_or(0.0);
            vec![gel.max(0.0) * t.ice.volatile.liquid_density; n]
        })
        .collect();
    // The rates of an ice's table over the places, kg/m²/s, each lying uppermost or not, and
    // whether each takes part.
    let rates = |i: usize, table: &[f64], surface: &[f64], exposed: &dyn Fn(usize) -> bool| {
        let value = |k: usize| read_table(table, g.levels, lat[k], thin(surface[k]), !own[i]);
        if own[i] {
            let r: Vec<f64> = (0..n).map(|k| value(k) / g.year).collect();
            let c = buffer(&r, &area, exposed);
            (0..n)
                .map(|k| {
                    if exposed(k) || r[k] > c {
                        r[k] - c
                    } else {
                        0.0
                    }
                })
                .collect::<Vec<f64>>()
        } else {
            let v: Vec<f64> = (0..n).map(value).collect();
            let th: Vec<f64> = surface.iter().map(|&h| thin(h)).collect();
            let air = vapour(&v, &th, &area, exposed);
            (0..n)
                .map(|k| {
                    let flux = exchange[i] * (air * th[k] - v[k]);
                    if exposed(k) || flux > 0.0 { flux } else { 0.0 }
                })
                .collect()
        }
    };
    // Each ice's table, its present's departure from its mean followed by the share its
    // exchange time answers the cycles by.
    let tables: Vec<Vec<f64>> = (0..count)
        .map(|i| {
            let t = &traps[i];
            if g.cycle <= 0.0 {
                return t.present.clone();
            }
            let layer = mass[i][0];
            let moved: f64 = rates(i, &t.present, &ground, &|_| true)
                .iter()
                .zip(&area)
                .map(|(r, a)| r.abs() * a)
                .sum::<f64>()
                / whole;
            let tau = if moved > 0.0 {
                layer / moved
            } else {
                f64::INFINITY
            };
            let follows = 1.0 / (1.0 + (std::f64::consts::TAU * tau / g.cycle).powi(2)).sqrt();
            t.mean
                .iter()
                .zip(&t.present)
                .map(|(m, p)| m + follows * (p - m))
                .collect()
        })
        .collect();
    let mut elapsed = 0.0;
    let mut steps = 0;
    for _ in 0..GATHER_MAX_STEPS {
        if elapsed >= g.age {
            break;
        }
        steps += 1;
        let surface: Vec<f64> = (0..n)
            .map(|k| ground[k] + (0..count).map(|i| mass[i][k] / rho[i]).sum::<f64>())
            .collect();
        let top: Vec<Option<usize>> = (0..n)
            .map(|k| order.iter().copied().find(|&i| mass[i][k] > 0.0))
            .collect();
        let flows: Vec<Vec<f64>> = (0..count)
            .map(|i| rates(i, &tables[i], &surface, &|k| top[k] == Some(i)))
            .collect();
        // A step moves at most GATHER_STEP_SHARE of what lies.
        let lying: f64 = (0..count)
            .map(|i| mass[i].iter().zip(&area).map(|(m, a)| m * a).sum::<f64>())
            .sum();
        let moving: f64 = flows
            .iter()
            .map(|f| f.iter().zip(&area).map(|(r, a)| r.abs() * a).sum::<f64>())
            .sum();
        if moving <= 0.0 || lying <= 0.0 {
            break;
        }
        let dt = (GATHER_STEP_SHARE * lying / moving).min(g.age - elapsed);
        for i in 0..count {
            // What each place loses, no more than it holds, and its gains scaled to match.
            let lost: f64 = (0..n)
                .map(|k| (-flows[i][k] * dt).clamp(0.0, mass[i][k]) * area[k])
                .sum();
            let won: f64 = (0..n).map(|k| (flows[i][k] * dt).max(0.0) * area[k]).sum();
            let scale = if won > 0.0 { lost / won } else { 0.0 };
            for k in 0..n {
                let d = flows[i][k] * dt;
                if d < 0.0 {
                    mass[i][k] -= (-d).min(mass[i][k]);
                } else {
                    mass[i][k] += d * scale;
                }
            }
        }
        elapsed += dt;
        // Settled once what moves through the rest of the age is a sliver of what lies.
        if moving * (g.age - elapsed) < GATHER_SETTLED * lying {
            break;
        }
    }
    if std::env::var_os("EXOGEN_DEBUG").is_some() {
        let bands = [-90.0, -60.0, -40.0, -20.0, 20.0, 40.0, 60.0, 90.0];
        for (i, t) in traps.iter().enumerate() {
            let total: f64 = (0..n).map(|k| mass[i][k] * area[k]).sum();
            let shares: Vec<String> = bands
                .windows(2)
                .map(|b| {
                    let part: f64 = (0..n)
                        .filter(|&k| (b[0]..b[1]).contains(&lat[k].to_degrees()))
                        .map(|k| mass[i][k] * area[k])
                        .sum();
                    format!("{}..{}: {:.2}", b[0], b[1], part / total.max(1e-300))
                })
                .collect();
            eprintln!(
                "gather: {} after {:.3e} years in {steps} steps, by band {}",
                t.ice.species,
                elapsed / 3.155_76e7,
                shares.join(" ")
            );
        }
    }
    traps
        .iter()
        .zip(mass)
        .map(|(t, m)| {
            let field = Field::from_fn(grid, |c, _| m[grid.index(c)] as f32);
            (t.ice.species, field)
        })
        .collect()
}

/// What holds the frost's gas's whole mass while its caps exchange with its air (F259): the
/// growth `c` the places' take less of so that what they gain and lose balances, the places
/// that lie under its ice taking part and the others where they would grow. Exact, the places
/// that would grow added from the fastest.
fn buffer(r: &[f64], area: &[f64], exposed: &dyn Fn(usize) -> bool) -> f64 {
    let (mut sum, mut weight) = (0.0, 0.0);
    let mut others: Vec<usize> = Vec::new();
    for (k, (&rk, &a)) in r.iter().zip(area).enumerate() {
        if exposed(k) {
            sum += a * rk;
            weight += a;
        } else {
            others.push(k);
        }
    }
    others.sort_by(|&a, &b| r[b].total_cmp(&r[a]));
    let mut c = if weight > 0.0 {
        sum / weight
    } else {
        f64::NEG_INFINITY
    };
    for &k in &others {
        if r[k] <= c {
            break;
        }
        sum += area[k] * r[k];
        weight += area[k];
        c = sum / weight;
    }
    c
}

/// The vapour of the air at the datum, Pa, at which what an ice's places give and take
/// balances (F259): the places under it taking part, and the others where it would lay
/// itself, the saturation pressure over it `v` below the air's there, its share `th` of the
/// datum's. Exact, the places added from the coldest.
fn vapour(v: &[f64], th: &[f64], area: &[f64], exposed: &dyn Fn(usize) -> bool) -> f64 {
    let (mut given, mut held) = (0.0, 0.0);
    let mut others: Vec<usize> = Vec::new();
    for k in 0..v.len() {
        if exposed(k) {
            given += area[k] * v[k];
            held += area[k] * th[k];
        } else {
            others.push(k);
        }
    }
    let level = |k: usize| v[k] / th[k].max(1e-300);
    others.sort_by(|&a, &b| level(a).total_cmp(&level(b)));
    let mut air = if held > 0.0 { given / held } else { 0.0 };
    for &k in &others {
        if level(k) >= air {
            break;
        }
        given += area[k] * v[k];
        held += area[k] * th[k];
        air = given / held;
    }
    air
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
        let cycle = eb.cycle(&climate, &ground, None, None);
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
        let caps = lay(&graph, &mut layers, &traps, &[]);
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
        let caps = lay(&graph, &mut layers, &traps, &[]);
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
        let caps = lay(&graph, &mut layers, &traps, &[]);
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

    /// F258: Pluto's nitrogen lays its caps in the bands some 30° from its equator that its
    /// tilt of 120° and its ground's seasons leave the coldest, where volatile transport
    /// models gather it (Bertrand et al. 2018), Sputnik Planitia's 25° north among them.
    #[test]
    fn pluto_s_nitrogen_lays_its_caps_in_its_coldest_bands() {
        let (graph, traps) = traps_of("Pluto");
        let mut layers = LayerSet::new(CubeSphere::new(64));
        let caps = lay(&graph, &mut layers, &traps, &[]);
        let cap = caps
            .iter()
            .find(|c| c.ice.species == "N2")
            .expect("nitrogen caps");
        let grid = cap.thickness.grid();
        let (mut banded, mut covered) = (0.0, 0.0);
        for (k, &t) in cap.thickness.as_slice().iter().enumerate() {
            if t > 0.0 {
                let c = grid.cell(k);
                let (_, lat) = lon_lat_from_direction(grid.center(c));
                let a = grid.cell_area(c);
                covered += a;
                if (10.0..40.0).contains(&lat.to_degrees().abs()) {
                    banded += a;
                }
            }
        }
        let share = banded / covered;
        assert!(share > 0.8, "{share} of the caps within 10° to 40°");
    }

    /// F242: a working liquid that is the frost's own gas lays its ice only where the ice
    /// grows: Titan's methane, its air over 90 K holding its ice above its point, lays none.
    #[test]
    fn titan_s_methane_lays_no_cap() {
        let (graph, traps) = traps_of("Titan");
        let mut layers = LayerSet::new(CubeSphere::new(32));
        let caps = lay(&graph, &mut layers, &traps, &[]);
        assert!(
            caps.is_empty(),
            "{:?}",
            caps.iter().map(|c| c.ice.species).collect::<Vec<_>>()
        );
    }

    /// A table of an ice's traps on the cycle's rows, one level, by latitude in degrees.
    fn rows_table(value: impl Fn(f64) -> f64) -> Vec<f64> {
        (0..CYCLE_LATITUDES)
            .map(|i| value(-90.0 + 180.0 * i as f64 / (CYCLE_LATITUDES - 1) as f64))
            .collect()
    }

    /// An ice of a preset's inventory.
    fn ice_of(graph: &ParamGraph, species: &str) -> Ice {
        *Ice::of_inventory(graph, [0.80, 0.82, 0.86])
            .iter()
            .find(|i| i.species == species)
            .unwrap()
    }

    /// The share of what an ice gathered on the places a test picks, and its whole over the
    /// body, kg/m².
    fn held(mass: &Field<f32>, pick: impl Fn(f64) -> bool) -> (f64, f64) {
        let grid = mass.grid();
        let (mut part, mut all, mut area) = (0.0, 0.0, 0.0);
        for k in 0..grid.cell_count() {
            let c = grid.cell(k);
            let a = grid.cell_area(c);
            let m = mass.as_slice()[k] as f64 * a;
            let (_, lat) = lon_lat_from_direction(grid.center(c));
            if pick(lat.to_degrees()) {
                part += m;
            }
            all += m;
            area += a;
        }
        (part / all, all / area)
    }

    const YEAR: f64 = 3.155_76e7;

    /// F259: the frost's own gas, gathered from a layer over the whole body, lies where its
    /// ice grows, its whole mass kept: Pluto's nitrogen, its traps growing within 20° of the
    /// equator and losing beyond.
    #[test]
    fn the_frost_s_gas_gathers_where_its_ice_grows() {
        let graph = crate::presets::preset("Pluto").unwrap().to_graph().0;
        let ice = ice_of(&graph, "N2");
        let table = rows_table(|lat| 2.0 - 0.1 * lat.abs());
        let cycled = [Cycled {
            ice,
            present: table.clone(),
            mean: table,
        }];
        let gathering = Gathering {
            levels: &[1.0],
            scale_height: 0.0,
            cycle: 0.0,
            year: YEAR,
            age: 1e7 * YEAR,
            temperature: 40.0,
            airless: false,
        };
        let heights = Field::new(CubeSphere::new(16));
        let gathered = gather(&graph, &heights, &cycled, Some("N2"), &gathering);
        let (share, whole) = held(&gathered[0].1, |lat| lat.abs() < 20.0);
        let layer = 30.0 * ice.volatile.liquid_density;
        assert!(share > 0.95, "{share} within 20°");
        assert!(
            (whole / layer - 1.0).abs() < 1e-4,
            "{whole} against {layer}"
        );
    }

    /// F259: a slow ice, gathered from a layer over the whole body, leaves its warm places and
    /// keeps both its poles where one is barely the colder, the little it would move between
    /// them in the body's age; a fast one, as the frost's gas moves, gathers all of it at the
    /// colder (Mars's polar layered deposits stand at both its poles).
    #[test]
    fn a_slow_ice_keeps_both_its_poles() {
        let graph = crate::presets::preset("Mars").unwrap().to_graph().0;
        let ice = ice_of(&graph, "H2O");
        let table = rows_table(|lat| {
            let polar = if lat > 0.0 { 1.0e-4 } else { 1.05e-4 };
            polar + 1.0e-2 * (1.0 - (lat.abs() / 90.0).powi(2))
        });
        let cycled = [Cycled {
            ice,
            present: table.clone(),
            mean: table,
        }];
        let heights = Field::new(CubeSphere::new(16));
        let run = |temperature: f64, age: f64| {
            let gathering = Gathering {
                levels: &[1.0],
                scale_height: 0.0,
                cycle: 0.0,
                year: YEAR,
                age: age * YEAR,
                temperature,
                airless: false,
            };
            let gathered = gather(&graph, &heights, &cycled, None, &gathering);
            let north = held(&gathered[0].1, |lat| lat > 60.0).0;
            let south = held(&gathered[0].1, |lat| lat < -60.0).0;
            (north, south)
        };
        let (north, south) = run(200.0, 4.5e9);
        println!("slow: north {north:.3}, south {south:.3}");
        assert!(north + south > 0.9, "{north} {south}");
        assert!(north > 0.2 && south > 0.2, "{north} {south}");
        // The same ice a million times as quick: all of it at the colder pole.
        let (north, south) = run(200.0e-6, 4.5e9);
        println!("fast: north {north:.3}, south {south:.3}");
        assert!(north > 0.9, "{north} {south}");
    }

    /// A disc of ice 10 degrees across from the north pole on flat ground on a grid of
    /// `cells` a face, its profile's thickness at each cell, m, `rise` 5 on a body 1,000 km
    /// in radius; and the thickness Nye's plastic dome gives at a cell, the square root of its
    /// centre's distance from the disc's edge.
    fn disc(cells: u32) -> (CubeSphere, Vec<f32>, impl Fn(usize) -> f64) {
        let grid = CubeSphere::new(cells);
        let ground = Field::new(grid);
        let (radius, rise, edge) = (1.0e6, 5.0, 10f64.to_radians());
        let covered: Vec<(usize, f64)> = (0..grid.cell_count())
            .filter(|&k| {
                let (_, lat) = lon_lat_from_direction(grid.center(grid.cell(k)));
                lat > 90f64.to_radians() - edge
            })
            .map(|k| (k, 0.0))
            .collect();
        let t = profile(grid, &ground, &covered, rise, radius, 0.0);
        let nye = move |k: usize| {
            let r = angle_between([0.0, 0.0, 1.0], grid.center(grid.cell(k)));
            rise * ((edge - r) * radius).sqrt()
        };
        (grid, t, nye)
    }

    /// F241: ice covering a disc on flat ground stands as Nye's plastic dome, its centre
    /// `rise` times the square root of its distance from the disc's edge above the ground:
    /// the cell by the pole, whose centre lies off it by half a cell's diagonal on a grid of
    /// an even count of cells a face, the distance from its own centre (F280).
    #[test]
    fn a_disc_of_ice_stands_as_a_plastic_dome() {
        let (grid, t, nye) = disc(128);
        let pole = grid.index(grid.locate([0.0, 0.0, 1.0]).0);
        let (centre, expect) = (t[pole] as f64, nye(pole));
        assert!(
            (centre / expect - 1.0).abs() < 0.03,
            "{centre} against {expect}"
        );
    }

    /// F280: the dome is round: half way in from its margin, the ice stands as thick along the
    /// grid's rows as across them, its distance in the true one from the margin, not summed
    /// along the grid's steps, which ran long across them by up to the square root of two.
    #[test]
    fn a_dome_of_ice_stands_round_across_the_grid() {
        let (grid, t, nye) = disc(128);
        let (s, c) = (5f64.to_radians().sin(), 5f64.to_radians().cos());
        let h = s / 2f64.sqrt();
        let along = [[s, 0.0, c], [0.0, s, c], [-s, 0.0, c], [0.0, -s, c]];
        let across = [[h, h, c], [h, -h, c], [-h, h, c], [-h, -h, c]];
        // Each place's thickness over what Nye's dome gives at its cell's centre.
        let share = |ds: &[[f64; 3]]| {
            ds.iter()
                .map(|&d| {
                    let k = grid.index(grid.locate(d).0);
                    t[k] as f64 / nye(k)
                })
                .sum::<f64>()
                / ds.len() as f64
        };
        let (rows, diagonals) = (share(&along), share(&across));
        println!("half way in, of Nye's dome: {rows:.3} along the grid, {diagonals:.3} across it");
        assert!((rows - 1.0).abs() < 0.05, "{rows}");
        assert!((diagonals / rows - 1.0).abs() < 0.02, "{diagonals} {rows}");
    }
}
