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
use crate::grid::projection::sample_monotone;
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

/// An ice's traps as its exchange answers the turning of the seasons (F259): its gas and its
/// traps on the latitudes' rows at `levels`, `read_table`'s values (the frost's own gas's
/// growth, kg/m² a year, the others' saturation pressure over the ice, Pa), which its caps are
/// laid by on the finer grid (F427).
#[derive(Clone, Debug)]
pub struct Turned {
    pub species: &'static str,
    pub table: Vec<f64>,
    pub levels: Vec<f64>,
}

/// What a cap reads of the ground's own climate on the finer grid: how much warmer each
/// cell's year peaks at its own height and slope than at its place's, K (F282); the climate's
/// year-mean temperature on its own grid, K, and the heat the air exchanges with the ground
/// per kelvin, W/m²/K, which set the share of the sunlight a slope faces its temperature
/// keeps (F237, F427); and the sea's full level the traps' heights stand above, m, not
/// finite without a sea, its ground then above the datum.
pub struct Fine<'a> {
    pub warmer: Option<&'a Field<f32>>,
    pub facing: Option<(&'a Field<f32>, f64)>,
    pub datum: f64,
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
        "ammonia" => &["NH3"],
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
/// its inventory or, a liquid's, until it covers all the ground it can, over the ground under
/// it rising or falling from its margin's (F418). An ice's coldest traps are those where it
/// needs the least vapour in the air (F240), the frost's own gas's those where its ice grows
/// the most (F242); an ice that is the body's working liquid stands only where it never
/// melts, the frost's own gas's where its ice grows.
///
/// Each cell's trap is its own on the finer grid (F427): its ice's table on the latitudes'
/// rows, an ice the orbit's cycles move the one its exchange answers them by (`turned`,
/// F259), the traps' otherwise, read at the height above the datum of the ground's shape the
/// grid resolves (F423), the caps already laid on it, and at the latitude that shape's slope
/// faces, by the share of the sunlight it faces its temperature keeps (F237); where each of
/// the climate's cells runs its own traps, theirs read smoothly at the point it faces. Read
/// at their places, smoothly between the climate's cells some 170 km apart on Mars, its
/// caps' outlines were the contours of their places' traps, blobs about the climate's cells;
/// and an ice gathered through the body's age covered the places it gathered at the most
/// first, piled without flowing some 16 km deep on a dozen of Mars's places, cut along the
/// grid's rows where their traps read alike, the rest of its caps its present season's
/// traps', not its own. A working liquid's ice stands where its own ground keeps it frozen
/// (F282): its place's warmest, read smoothly between the climate's places, `fine`'s
/// `warmer` at each cell by its own height and slope than at its place's, K, so that where
/// the ice thins to its margin the ground under it decides where it melts.
pub fn lay(
    graph: &ParamGraph,
    layers: &mut LayerSet,
    traps: &Traps,
    turned: &[Turned],
    fine: Option<&Fine>,
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
    let warmer = fine.and_then(|f| f.warmer);
    let datum = fine.map_or(f64::NAN, |f| f.datum);
    let scale_height = if graph.f64("atmosphere.surface_pressure") > 0.0 {
        graph.f64("atmosphere.scale_height")
    } else {
        0.0
    };
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
        // Where an ice stands by how far its warmest keeps under its melting, its own ground's
        // warmest (F282); the frost's own gas's stands where its ice grows.
        let by_heat = traps
            .growth
            .as_ref()
            .is_none_or(|(g, _)| g.species != ice.species);
        let ground = standing(layers);
        // The ice's own traps (F427): the table its exchange answers the orbit's cycles by or
        // the traps', the frost's own gas's its growth; none where the climate's cells run
        // their own.
        let grows = !by_heat;
        let table = turned
            .iter()
            .find(|t| t.species == ice.species)
            .map(|t| (t.table.as_slice(), t.levels.as_slice()))
            .or_else(|| {
                let rows = traps.rows.as_ref()?;
                let values = rows.values.iter().find(|(s, _)| *s == ice.species)?;
                Some((values.1.as_slice(), rows.levels.as_slice()))
            });
        // The ground's shape the grid resolves above the datum (F423), the caps already laid
        // on it: a cell's own relief at its scale decides neither its height nor its slope.
        let land: Vec<f32> = layers
            .liquid
            .as_slice()
            .iter()
            .map(|&l| if l > 0.0 { 0.0 } else { 1.0 })
            .collect();
        let above: Vec<f32> = ground
            .as_slice()
            .iter()
            .map(|&h| {
                let h = if h.is_finite() { h } else { 0.0 };
                if datum.is_finite() {
                    (h - datum as f32).max(0.0)
                } else {
                    h
                }
            })
            .collect();
        let resolved = crate::shares::shape_of(&grid, &above, &land);
        let resolved = Field::from_fn(grid, |c, _| resolved[grid.index(c)]);
        let own_trap = OwnTrap {
            table,
            grows,
            coarse,
            rank: &rank,
            resolved: &resolved,
            scale_height,
            radius,
            facing: fine.and_then(|f| f.facing),
        };
        // Each cell's own trap on the finer grid; where its ice may stand, coldest first.
        // Its place's warmest read by the monotone cubic (F291): read bilinearly, the edge where
        // the ice may stand kinked at the climate's cells where their warmest hovers at its
        // melting.
        let mut ranks: Vec<(usize, f64)> = (0..grid.cell_count())
            .into_par_iter()
            .filter_map(|k| {
                if layers.liquid.as_slice()[k] > 0.0 {
                    return None;
                }
                let c = grid.cell(k);
                let d = grid.center(c);
                let own = match warmer {
                    Some(w) if by_heat => w.as_slice()[k],
                    _ => 0.0,
                };
                if melts && sample_monotone(&coarse, &stands, d) - own < 0.0 {
                    return None;
                }
                let rank = own_trap.at(c, d);
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

/// The surface the caps stand on and their margins rise from (F287): the ground, and over the
/// sea the sea's surface, where an ice that ends at the sea meets the water; from the sea's
/// floor beside it, its plastic profile stood below the sea a cell and more in, and left the
/// low ground along its coasts bare.
fn standing(layers: &LayerSet) -> Field<f32> {
    let mut surface = layers.height();
    for (h, &l) in surface
        .as_mut_slice()
        .iter_mut()
        .zip(layers.liquid.as_slice())
    {
        if l > 0.0 {
            *h += l;
        }
    }
    surface
}

/// How an ice reads its own trap at a cell of the finer grid (F427), the lower the colder:
/// its `table` on the latitudes' rows at its levels where it has one, read at the height of
/// the ground's `resolved` shape, the share of the datum's pressure that height keeps by the
/// `scale_height`, and at the latitudes the cell and its slope face, by the share of the
/// sunlight it faces that its temperature keeps (F237) at the climate's year-mean temperature
/// with the air's exchange, `facing`: the saturation pressure over it, over that share, in the
/// logarithm, or its growth, `grows`, less; else its places' `rank` on the `coarse` grid, read
/// smoothly at the cell and at the point it faces.
struct OwnTrap<'a> {
    table: Option<(&'a [f64], &'a [f64])>,
    grows: bool,
    coarse: CubeSphere,
    rank: &'a [f32],
    resolved: &'a Field<f32>,
    scale_height: f64,
    radius: f64,
    facing: Option<(&'a Field<f32>, f64)>,
}

impl OwnTrap<'_> {
    fn at(&self, c: crate::grid::Cell, d: [f64; 3]) -> f64 {
        let (faced, keep) = match self.facing {
            Some((mean, exchange)) => {
                let t = sample_monotone(&mean.grid(), mean.as_slice(), d) as f64;
                let normal = crate::climate::surface_normal(self.resolved, c, d, self.radius);
                (normal, crate::climate::kept_share(t, exchange))
            }
            None => (d, 0.0),
        };
        match self.table {
            Some((values, levels)) => {
                let h = self.resolved.get(c) as f64;
                let thin = if self.scale_height > 0.0 {
                    (-h / self.scale_height).exp()
                } else {
                    1.0
                };
                // Read within its table's levels: beyond them its traps stood at the last
                // level's, all alike.
                assert!(
                    levels.len() < 2
                        || (levels[0] * (1.0 - 1e-6)..=levels[levels.len() - 1] * (1.0 + 1e-6))
                            .contains(&thin),
                    "a trap read at {h:.0} m, beyond its table's levels"
                );
                let read = |p: [f64; 3]| {
                    let lat = lon_lat_from_direction(p).1;
                    let v = read_table(values, levels, lat, thin, !self.grows);
                    if self.grows { -v } else { v.max(1e-300).ln() }
                };
                let r = (1.0 - keep) * read(d) + keep * read(faced);
                if self.grows { r } else { r - thin.ln() }
            }
            None => {
                let read = |p: [f64; 3]| sample_monotone(&self.coarse, self.rank, p) as f64;
                (1.0 - keep) * read(d) + keep * read(faced)
            }
        }
    }
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
/// each the ground plus the plastic thickness the distance in from that point gives over the
/// ground rising or falling from the margin's to its own (F418, `plastic`), on flat ground
/// `rise` times the square root of the distance: measured from the margin's ground alone, the
/// profile ran under ground rising inland as fast as its square root, a coast's slope bare
/// between the ice at its foot and the ice above it, its edges where each cell's own relief
/// decided; found outward from the margin as distances are, the margin's ground the higher
/// of its edge's own and the surface beyond it (F287): so that an ice fills a basin toward its
/// rim, and on a
/// plateau ending in a cliff stands from the plateau's edge, not from the ground or the sea
/// below it, which left a high coast bare; each place carrying the point of the margin it
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
            let base = level(j).max(level(k));
            let surface = level(k) + plastic(d, (level(k) - base) / d, rise);
            if surface < best[k].0 {
                best[k] = (surface, base, point);
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
            let s = level(j) + plastic(d, (level(j) - from) / d, rise);
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

/// The thickness of a perfectly plastic ice `d` m in from its margin over ground rising
/// `slope` m a metre from the margin's (F418), m: its surface's slope times its thickness the
/// yield stress over its weight, H (dH/dx + slope) = c, c = rise² / 2 = CAP_YIELD_STRESS /
/// (ρ g) (Nye 1951), from none at the margin. On flat ground `rise` √d; over ground rising
/// it thickens toward c / slope and is never none, its surface rising with the ground; over
/// ground falling it is the thicker. For a uniform slope, with u = slope H / c and y = d
/// slope² / c, -u - ln(1 - u) = y, solved by Newton's method in w = -ln(1 - u) over rising
/// ground, w - 1 + e^-w = y, and in v = -u over falling ground, v - ln(1 + v) = y, each
/// convex and rising from its root's side.
fn plastic(d: f64, slope: f64, rise: f64) -> f64 {
    let c = 0.5 * rise * rise;
    if d <= 0.0 || c <= 0.0 {
        return 0.0;
    }
    let y = d * slope * slope / c;
    if !y.is_finite() || y < 1e-12 {
        return rise * d.sqrt();
    }
    // Each start from the series near none, s (1 + s / 6) and s (1 + s / 3), s = √(2y), or
    // the root's growth far from it.
    let s = (2.0 * y).sqrt();
    if slope > 0.0 {
        let mut w = if y < 1.0 {
            s * (1.0 + s / 6.0)
        } else {
            y + 1.0
        };
        for _ in 0..CAP_PROFILE_NEWTON {
            let f = w + (-w).exp_m1() - y;
            let df = -(-w).exp_m1();
            if df <= 0.0 {
                break;
            }
            w = (w - f / df).max(f64::MIN_POSITIVE);
        }
        c * -(-w).exp_m1() / slope
    } else {
        let mut v = if y < 1.0 {
            s * (1.0 + s / 3.0)
        } else {
            y + (1.0 + y).ln()
        };
        for _ in 0..CAP_PROFILE_NEWTON {
            let f = v - v.ln_1p() - y;
            let df = v / (1.0 + v);
            if df <= 0.0 {
                break;
            }
            v = (v - f / df).max(f64::MIN_POSITIVE);
        }
        c * v / -slope
    }
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

/// What an ice's traps answer the turning of the seasons by (F259): the traps' table's
/// levels, shares of the datum's pressure; the pressure's scale height, m, none on an airless
/// body; the period over which the seasons turn, s, none where they hold; the body's year, s;
/// its mean surface temperature, K; and whether it has air.
#[derive(Clone, Copy, Debug)]
pub struct Turning<'a> {
    pub levels: &'a [f64],
    pub scale_height: f64,
    pub cycle: f64,
    pub year: f64,
    pub temperature: f64,
    pub airless: bool,
}

/// Each ice's traps as its exchange answers the turning of the seasons (F259), on the
/// latitudes' rows at the levels of `turning`, which its caps are laid by on the finer grid
/// (F427). An ice follows its traps' departure from their mean through the turning by the
/// share 1 / sqrt(1 + (2 pi tau / P)^2) a store of exchange time tau answers a cycle of period
/// P by (the lag left out), its exchange time the time its inventory, a layer over the whole
/// body on the climate's grid of `heights` (the ground's above the datum, m), takes to move at
/// its rates there: the frost's own gas, moving in years, the present season's; a slow ice,
/// moving in millions of years, the turning's mean, so that a pole the present favours does
/// not take what the other holds. Its rates: the frost's own gas's its growth less its mean
/// over the body, the air buffered by its caps (F258); every other ice's the saturation
/// pressure over it against the vapour the air holds, at the air's turbulent exchange,
/// C_H U / (R_s T) per pascal, or into the vacuum at the Hertz-Knudsen rate on an airless
/// body, the air's vapour the one at which what all of it gives and takes balances
/// (Ingersoll 1970).
///
/// The ice was gathered through the body's age besides, from that layer, without flowing, and
/// its caps covered the places it gathered at the most first: on Mars it piled 9 to 16 km deep
/// on a dozen of the climate's places, and its caps were the blocks of those places. Its caps
/// are now laid by its table on the finer grid, and the gathering is gone (F427).
pub fn turned(
    graph: &ParamGraph,
    heights: &Field<f32>,
    traps: &[Cycled],
    frost: Option<&'static str>,
    turning: &Turning,
) -> Vec<Turned> {
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
    let mean = |values: &[f64]| -> f64 {
        values.iter().zip(&area).map(|(v, a)| v * a).sum::<f64>() / whole
    };
    let lat: Vec<f64> = (0..n)
        .map(|k| lon_lat_from_direction(grid.center(grid.cell(k))).1)
        .collect();
    let thin: Vec<f64> = heights
        .as_slice()
        .iter()
        .map(|&h| {
            let h = if h.is_finite() { h as f64 } else { 0.0 };
            if turning.scale_height > 0.0 {
                (-h / turning.scale_height).exp()
            } else {
                1.0
            }
        })
        .collect();
    traps
        .iter()
        .map(|t| {
            let own = Some(t.ice.species) == frost;
            let table = if turning.cycle <= 0.0 {
                t.present.clone()
            } else {
                // Its rates at the present season's traps over its layer, kg/m²/s.
                let value =
                    |k: usize| read_table(&t.present, turning.levels, lat[k], thin[k], !own);
                let rates: Vec<f64> = if own {
                    let growth: Vec<f64> = (0..n).map(|k| value(k) / turning.year).collect();
                    let held = mean(&growth);
                    growth.iter().map(|g| g - held).collect()
                } else {
                    let rt = t.ice.volatile.r_specific * turning.temperature.max(1.0);
                    let exchange = if turning.airless {
                        (1.0 / (std::f64::consts::TAU * rt)).sqrt()
                    } else {
                        AIR_EXCHANGE_COEFFICIENT * SURFACE_WIND_SPEED / rt
                    };
                    let v: Vec<f64> = (0..n).map(value).collect();
                    let air = mean(&v) / mean(&thin);
                    (0..n).map(|k| exchange * (air * thin[k] - v[k])).collect()
                };
                let moved = mean(&rates.iter().map(|r| r.abs()).collect::<Vec<f64>>());
                let gel = inventory.get(t.ice.species).copied().unwrap_or(0.0);
                let layer = gel.max(0.0) * t.ice.volatile.liquid_density;
                let tau = if moved > 0.0 {
                    layer / moved
                } else {
                    f64::INFINITY
                };
                let follows =
                    1.0 / (1.0 + (std::f64::consts::TAU * tau / turning.cycle).powi(2)).sqrt();
                t.mean
                    .iter()
                    .zip(&t.present)
                    .map(|(m, p)| m + follows * (p - m))
                    .collect()
            };
            Turned {
                species: t.ice.species,
                table,
                levels: turning.levels.to_vec(),
            }
        })
        .collect()
}

#[cfg(test)]
pub(crate) mod tests {
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

    /// A preset's caps laid on its traps over a 64-cell face, laid once for the whole-body
    /// checks that read them, each of which laid them on its own before (the tests' tiers).
    pub(crate) struct Laid {
        pub(crate) layers: LayerSet,
        pub(crate) caps: Vec<Cap>,
    }

    impl Laid {
        fn of(name: &str) -> Self {
            let (graph, traps) = traps_of(name);
            let mut layers = LayerSet::new(CubeSphere::new(64));
            let caps = lay(&graph, &mut layers, &traps, &[], None);
            Self { layers, caps }
        }
    }

    /// Mars's caps, read by each check of them; Pluto's are laid on the year its climate's
    /// checks read (climate.rs, `pluto_s_and_triton_s_years`).
    #[test]
    #[cfg_attr(not(feature = "whole-body"), ignore = "whole body")]
    fn mars_s_caps() {
        let mars = Laid::of("Mars");
        crate::testkit::checks(&[
            ("mars_s_water_lays_its_polar_caps", &|| {
                mars_s_water_lays_its_polar_caps(&mars)
            }),
            (
                "the_frost_s_own_gas_lays_its_caps_where_its_ice_grows",
                &|| the_frost_s_own_gas_lays_its_caps_where_its_ice_grows(&mars),
            ),
        ]);
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
    fn mars_s_water_lays_its_polar_caps(laid: &Laid) {
        let Laid { layers, caps, .. } = laid;
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
    pub(crate) fn an_ice_gathers_in_its_coldest_traps(laid: &Laid) {
        let Laid { caps, .. } = laid;
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
    fn the_frost_s_own_gas_lays_its_caps_where_its_ice_grows(laid: &Laid) {
        let Laid { caps, .. } = laid;
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
    pub(crate) fn pluto_s_nitrogen_lays_its_caps_in_its_coldest_bands(laid: &Laid) {
        let Laid { caps, .. } = laid;
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
    #[cfg_attr(not(feature = "whole-body"), ignore = "whole body")]
    fn titan_s_methane_lays_no_cap() {
        let (graph, traps) = traps_of("Titan");
        let mut layers = LayerSet::new(CubeSphere::new(32));
        let caps = lay(&graph, &mut layers, &traps, &[], None);
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

    const YEAR: f64 = 3.155_76e7;

    /// F259: an ice follows its traps' departure from their mean through the turning of the
    /// seasons by how quickly it moves: Mars's water, its inventory moving in millions of
    /// years, by the turning's mean, so that a pole the present season favours does not take
    /// what the other holds; the same ice under a turning a million times as slow by the
    /// present season's; its carbon dioxide, the frost's own gas, moving in years, by the
    /// present season's.
    #[test]
    fn an_ice_answers_the_turning_by_how_quickly_it_moves() {
        let mars = crate::presets::preset("Mars").unwrap().to_graph().0;
        let follows =
            |graph: &ParamGraph, species: &str, frost: Option<&'static str>, slower: f64| {
                let ice = ice_of(graph, species);
                let (present, mean) = if frost.is_some() {
                    (rows_table(|lat| 100.0 * lat / 90.0), rows_table(|_| 0.0))
                } else {
                    (
                        rows_table(|lat| 1.0e-3 * (1.0 + lat / 90.0) + 1.0e-4),
                        rows_table(|_| 1.1e-3),
                    )
                };
                let cycled = [Cycled {
                    ice,
                    present: present.clone(),
                    mean: mean.clone(),
                }];
                let turning = Turning {
                    levels: &[1.0],
                    scale_height: 0.0,
                    cycle: slower * graph.f64("orbit.precession_period"),
                    year: YEAR,
                    temperature: graph.f64("atmosphere.surface_temperature_mean"),
                    airless: false,
                };
                let heights = Field::new(CubeSphere::new(16));
                let t = &turned(graph, &heights, &cycled, frost, &turning)[0];
                let pole = CYCLE_LATITUDES - 1;
                (t.table[pole] - mean[pole]) / (present[pole] - mean[pole])
            };
        let slow = follows(&mars, "H2O", None, 1.0);
        let quick = follows(&mars, "H2O", None, 1.0e6);
        let frost = follows(&mars, "CO2", Some("CO2"), 1.0);
        println!(
            "the share of the present's departure followed: Mars's water {slow:.4}, under a turning a million times as slow {quick:.4}, its carbon dioxide {frost:.4}"
        );
        assert!(slow < 0.05, "{slow}");
        assert!(quick > 0.95, "{quick}");
        assert!(frost > 0.95, "{frost}");
    }

    /// F259, F427: Mars's water, laid by the traps its exchange answers the turning of the
    /// seasons by, the turning's mean holding each pole alike, stands at both its poles on flat
    /// ground, where its present season's traps alone lay all of it at the north (Mars's polar
    /// layered deposits stand at both its poles).
    #[test]
    #[cfg_attr(not(feature = "whole-body"), ignore = "whole body")]
    fn mars_s_water_stands_at_both_its_poles() {
        let (graph, traps) = traps_of("Mars");
        let balance = EnergyBalance::from_graph(&graph);
        let grid = CubeSphere::new(CLIMATE_FACE_SIZE);
        let ground = Ground::from_graph(&graph);
        let downscale = Downscale::from_graph(&graph);
        let ices = Ice::of_inventory(&graph, [0.80, 0.82, 0.86]);
        let heights = Field::new(grid);
        let levels = [1.0];
        let world = (grid, None, &ground, None, None, None);
        let cycled =
            super::super::cycled_traps(&graph, &balance, world, &ices, &levels, &downscale);
        let turning = Turning {
            levels: &levels,
            scale_height: downscale.scale_height,
            cycle: graph.f64("orbit.precession_period"),
            year: graph.f64("orbit.year"),
            temperature: graph.f64("atmosphere.surface_temperature_mean"),
            airless: false,
        };
        let frost = balance.frost.map(|f| f.species);
        let turned = turned(&graph, &heights, &cycled, frost, &turning);
        let poles = |turned: &[Turned]| {
            let mut layers = LayerSet::new(CubeSphere::new(64));
            let laid = lay(&graph, &mut layers, &traps, turned, None);
            let cap = laid.iter().find(|c| c.ice.species == "H2O").expect("water");
            let grid = cap.thickness.grid();
            let (mut north, mut south, mut all) = (0.0, 0.0, 0.0);
            for (k, &t) in cap.thickness.as_slice().iter().enumerate() {
                let v = t as f64 * grid.cell_area(grid.cell(k));
                let lat = lat_lon(grid, k).0;
                if lat > 60.0 {
                    north += v;
                } else if lat < -60.0 {
                    south += v;
                }
                all += v;
            }
            (north / all, south / all)
        };
        let (north, south) = poles(&turned);
        let (now_north, now_south) = poles(&[]);
        println!(
            "Mars's water: north {north:.3}, south {south:.3}; by the present season's traps {now_north:.3}, {now_south:.3}"
        );
        assert!(north > 0.25 && south > 0.25, "north {north}, south {south}");
        assert!(north + south > 0.9, "north {north}, south {south}");
        assert!(now_north > 0.9, "{now_north}");
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

    /// F282: a working liquid's ice stands where its own ground keeps it frozen: on the
    /// Earth-like, whose water's ice stands nowhere in its tropics on its places' warmest,
    /// ground standing 40 K colder than its places within 20 degrees of the equator, as high
    /// ground stands, holds it there.
    #[test]
    #[cfg_attr(not(feature = "whole-body"), ignore = "whole body")]
    fn a_working_liquid_s_ice_stands_where_its_own_ground_keeps_it_frozen() {
        let (graph, traps) = traps_of("Earth");
        let grid = CubeSphere::new(32);
        let tropics = |k: usize| {
            lon_lat_from_direction(grid.center(grid.cell(k))).1.abs() < 20f64.to_radians()
        };
        let covered = |warmer: Option<&Field<f32>>| {
            let mut layers = LayerSet::new(grid);
            let fine = warmer.map(|w| Fine {
                warmer: Some(w),
                facing: None,
                datum: f64::NAN,
            });
            let caps = lay(&graph, &mut layers, &traps, &[], fine.as_ref());
            let Some(water) = caps.iter().find(|c| c.ice.species == "H2O") else {
                return 0;
            };
            let t = water.thickness.as_slice();
            (0..grid.cell_count())
                .filter(|&k| tropics(k) && t[k] > 0.0)
                .count()
        };
        let cold = Field::from_fn(
            grid,
            |c, _| if tropics(grid.index(c)) { -40.0 } else { 0.0 },
        );
        let (flat, high) = (covered(None), covered(Some(&cold)));
        println!(
            "the tropics' cells under the water's ice: {flat} on their places' warmest, {high} on their own"
        );
        assert_eq!(flat, 0);
        assert!(high > 0);
    }

    /// F287: an ice cap over an island in a sea 3 km deep stands from the higher of its edge's
    /// own ground and the sea's surface, every cell it covers under ice, a low island's and a
    /// high plateau's ending in cliffs alike; from the sea's floor its plastic profile stood
    /// below the coast's ground near its margin and left it bare, and from the sea's surface
    /// below the plateau's edge.
    #[test]
    fn an_island_s_cap_stands_from_its_own_edge() {
        let grid = CubeSphere::new(64);
        for high in [200.0f32, 1_500.0] {
            let mut layers = LayerSet::new(grid);
            let mut cells = Vec::new();
            for k in 0..grid.cell_count() {
                let (_, lat) = lon_lat_from_direction(grid.center(grid.cell(k)));
                if lat > 70f64.to_radians() {
                    layers.bedrock.as_mut_slice()[k] = high;
                    cells.push((k, 0.0));
                } else {
                    layers.bedrock.as_mut_slice()[k] = -3_000.0;
                    layers.liquid.as_mut_slice()[k] = 3_000.0;
                }
            }
            let t = profile(grid, &standing(&layers), &cells, 5.0, 1.0e6, 0.0);
            let bare = cells.iter().filter(|&&(k, _)| t[k] <= 0.0).count();
            println!(
                "an island {high} m high: {bare} of {} cells bare",
                cells.len()
            );
            assert_eq!(bare, 0, "{high}");
        }
    }

    /// F418: the plastic thickness over ground rising or falling at a uniform slope is the
    /// profile Nye's law gives, H (dH/dx + slope) = c, found independently: the distance in at
    /// which it reaches each thickness, x(H) = ∫ h / (c - slope h) dh from none, by Simpson's
    /// rule; on flat ground rise √d.
    #[test]
    fn a_plastic_ice_thickens_by_nye_s_law_over_sloping_ground() {
        let rise = 3.34f64;
        let c = 0.5 * rise * rise;
        let mut worst = 0.0f64;
        for slope in [-0.1, -0.01, -1e-4, 1e-4, 0.01, 0.1] {
            // Thicknesses up to nine tenths of the most a rising slope holds.
            let top = if slope > 0.0 {
                0.9 * c / slope
            } else {
                3_000.0
            };
            for f in [0.05, 0.3, 0.6, 1.0] {
                let h = f * top;
                let n = 2_000;
                let step = h / n as f64;
                let g = |x: f64| x / (c - slope * x);
                let x = (0..n)
                    .map(|i| {
                        let (a, b) = (i as f64 * step, (i + 1) as f64 * step);
                        step / 6.0 * (g(a) + 4.0 * g(0.5 * (a + b)) + g(b))
                    })
                    .sum::<f64>();
                let off = (plastic(x, slope, rise) / h - 1.0).abs();
                worst = worst.max(off);
                assert!(
                    off < 1e-6,
                    "slope {slope}, {h} m at {x} m: {}",
                    plastic(x, slope, rise)
                );
            }
        }
        for d in [1.0, 1e3, 1e5, 1e6] {
            assert!((plastic(d, 0.0, rise) / (rise * d.sqrt()) - 1.0).abs() < 1e-12);
        }
        println!("the plastic thickness within {worst:.1e} of Nye's law integrated");
    }

    /// F418: an ice cap over ground rising inland at a hundredth from the sea about it, 3.5 km
    /// at its pole 20 degrees in, covers it all, its thickness toward c / slope inland, as
    /// Nye's law over that ground gives: measured from the margin's ground alone its profile
    /// ran under the ground from some 250 km in, the inner cap bare.
    #[test]
    fn a_cap_covers_ground_rising_inland() {
        let grid = CubeSphere::new(64);
        let (radius, rise, slope) = (1.0e6, 5.0, 0.01);
        let edge = 70f64.to_radians();
        let mut layers = LayerSet::new(grid);
        let mut cells = Vec::new();
        for k in 0..grid.cell_count() {
            let (_, lat) = lon_lat_from_direction(grid.center(grid.cell(k)));
            if lat > edge {
                layers.bedrock.as_mut_slice()[k] = (slope * (lat - edge) * radius) as f32;
                cells.push((k, 0.0));
            } else {
                layers.bedrock.as_mut_slice()[k] = -3_000.0;
                layers.liquid.as_mut_slice()[k] = 3_000.0;
            }
        }
        let t = profile(grid, &standing(&layers), &cells, rise, radius, 0.0);
        let bare = cells.iter().filter(|&&(k, _)| t[k] <= 0.0).count();
        let pole = grid.index(grid.locate([0.0, 0.0, 1.0]).0);
        let (_, lat) = lon_lat_from_direction(grid.center(grid.cell(pole)));
        let nye = plastic((lat - edge) * radius, slope, rise);
        let most = 0.5 * rise * rise / slope;
        println!(
            "{bare} of {} cells bare; the pole's ice {:.0} m, Nye's law {nye:.0}, at most {most:.0}",
            cells.len(),
            t[pole]
        );
        assert_eq!(bare, 0);
        assert!(
            (t[pole] as f64 / nye - 1.0).abs() < 0.05,
            "{} {nye}",
            t[pole]
        );
        assert!(nye > 0.9 * most);
    }

    /// Traps on the climate's grid from a table on the latitudes' rows at one level, the
    /// saturation pressure over the ice `ln_v` in the logarithm by latitude in degrees.
    fn rows_traps(ice: Ice, ln_v: impl Fn(f64) -> f64) -> Traps {
        let coarse = CubeSphere::new(CLIMATE_FACE_SIZE);
        let table = rows_table(|lat| ln_v(lat).exp());
        let vapour = Field::from_fn(coarse, |_, d| {
            read_table(&table, &[1.0], lon_lat_from_direction(d).1, 1.0, true) as f32
        });
        Traps {
            growth: None,
            vapour: vec![(ice, vapour)],
            warmest: vec![Field::from_fn(coarse, |_, _| 150.0)],
            rows: Some(crate::climate::Rows {
                levels: vec![1.0],
                values: vec![(ice.species, table)],
            }),
        }
    }

    /// The thickness of a preset's water laid on `layers` under `traps`.
    fn water_laid(
        graph: &ParamGraph,
        layers: &mut LayerSet,
        traps: &Traps,
        turned: &[Turned],
    ) -> Vec<f32> {
        let caps = lay(graph, layers, traps, turned, None);
        let cap = caps.iter().find(|c| c.ice.species == "H2O").expect("water");
        cap.thickness.as_slice().to_vec()
    }

    /// Latitude and longitude of a cell, degrees.
    fn lat_lon(grid: CubeSphere, k: usize) -> (f64, f64) {
        let (lon, lat) = lon_lat_from_direction(grid.center(grid.cell(k)));
        (lat.to_degrees(), lon.to_degrees())
    }

    /// F427: a cap's margin follows the ground's own traps on the finer grid: Mars's water,
    /// under traps colder toward its north pole, covers the floor of a basin 2 km deep beyond
    /// the margin it reaches over flat ground, the floor's denser air its colder trap, and
    /// leaves the flat ground at the basin's latitude bare; read at their places' traps on
    /// the climate's grid, which stand at their heights' means, the basin's cells stood as
    /// the flat ground's.
    #[test]
    fn a_cap_s_margin_follows_its_own_trap_on_the_finer_grid() {
        let graph = crate::presets::preset("Mars").unwrap().to_graph().0;
        let traps = rows_traps(ice_of(&graph, "H2O"), |lat| -8.0 + 0.02 * (90.0 - lat));
        let grid = CubeSphere::new(128);
        let centre = crate::grid::direction_from_lon_lat(0.0, 70f64.to_radians());
        let reach = 6f64.to_radians();
        let mut layers = LayerSet::new(grid);
        for k in 0..grid.cell_count() {
            let r = angle_between(centre, grid.center(grid.cell(k)));
            if r < reach {
                layers.bedrock.as_mut_slice()[k] = (-2_000.0 * (1.0 - (r / reach).powi(2))) as f32;
            }
        }
        let t = water_laid(&graph, &mut layers, &traps, &[]);
        let at = |lon: f64, lat: f64| {
            let d = crate::grid::direction_from_lon_lat(lon.to_radians(), lat.to_radians());
            t[grid.index(grid.locate(d).0)]
        };
        // The margin over flat ground: the lowest latitude covered on the far side.
        let margin = (0..grid.cell_count())
            .filter(|&k| t[k] > 0.0 && lat_lon(grid, k).1.abs() > 90.0)
            .map(|k| lat_lon(grid, k).0)
            .fold(90.0f64, f64::min);
        println!(
            "the margin over flat ground at {margin:.1}°; the basin's floor at 70° under {:.0} m",
            at(0.0, 70.0)
        );
        assert!(margin > 72.0, "{margin}");
        assert!(at(0.0, 70.0) > 0.0);
        assert_eq!(at(180.0, 70.0), 0.0);
    }

    /// F427: a slope reads its trap at the latitude it faces, by the share of the sunlight it
    /// faces its temperature keeps: ground at 75° falling 2° toward the pole, its air's
    /// exchange all but none, reads the trap of 77°; with the air holding it at its
    /// surroundings' its own.
    #[test]
    fn a_slope_reads_its_trap_where_it_faces() {
        let grid = CubeSphere::new(128);
        let radius = 3.39e6;
        let tilt = 2f64.to_radians();
        let resolved = Field::from_fn(grid, |_, d| {
            let lat = lon_lat_from_direction(d).1;
            (-(lat - 75f64.to_radians()) * radius * tilt.tan()) as f32
        });
        let ln_v = |lat: f64| -8.0 + 0.2 * (90.0 - lat);
        let table = rows_table(|lat| ln_v(lat).exp());
        let mean = Field::from_fn(CubeSphere::new(CLIMATE_FACE_SIZE), |_, _| 150.0);
        let reader = |exchange: f64| OwnTrap {
            table: Some((table.as_slice(), &[1.0][..])),
            grows: false,
            coarse: mean.grid(),
            rank: &[],
            resolved: &resolved,
            scale_height: 0.0,
            radius,
            facing: Some((&mean, exchange)),
        };
        let c = grid
            .locate(crate::grid::direction_from_lon_lat(0.3, 75f64.to_radians()))
            .0;
        let d = grid.center(c);
        let lat = lon_lat_from_direction(d).1.to_degrees();
        let (radiating, held) = (reader(1e-9).at(c, d), reader(1e9).at(c, d));
        println!(
            "at {lat:.2}°: {radiating:.3} where it radiates, {held:.3} held by its air; the trap of {:.2}° {:.3}",
            lat + 2.0,
            ln_v(lat + 2.0)
        );
        assert!((radiating - ln_v(lat + 2.0)).abs() < 0.02, "{radiating}");
        assert!((held - ln_v(lat)).abs() < 0.02, "{held}");
    }

    /// F427: a cap's margin is not speckled by each cell's own relief: Mars's water over flat
    /// ground roughened cell by cell by 400 m, its trap read at the ground's shape the grid
    /// resolves (F423), keeps no lone covered cell and no lone bare one; read at each cell's
    /// own height, the roughness moved its trap by as much as a degree and a half of
    /// latitude and scattered its margin.
    #[test]
    fn a_cap_s_margin_is_not_speckled_by_its_cells_own_relief() {
        let graph = crate::presets::preset("Mars").unwrap().to_graph().0;
        let traps = rows_traps(ice_of(&graph, "H2O"), |lat| -8.0 + 0.02 * (90.0 - lat));
        let grid = CubeSphere::new(128);
        let mut layers = LayerSet::new(grid);
        for (k, h) in layers.bedrock.as_mut_slice().iter_mut().enumerate() {
            let u =
                ((k as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15) >> 40) as f64 / (1u64 << 24) as f64;
            *h = (800.0 * (u - 0.5)) as f32;
        }
        let t = water_laid(&graph, &mut layers, &traps, &[]);
        let (mut lone, mut holes) = (0, 0);
        for k in 0..grid.cell_count() {
            let c = grid.cell(k);
            let covered = [(1, 0), (-1, 0), (0, 1), (0, -1)]
                .iter()
                .filter(|&&(di, dj)| t[grid.index(grid.neighbor(c, di, dj))] > 0.0)
                .count();
            if t[k] > 0.0 && covered == 0 {
                lone += 1;
            }
            if t[k] <= 0.0 && covered == 4 {
                holes += 1;
            }
        }
        println!("{lone} lone covered cells and {holes} lone bare ones");
        assert_eq!(lone + holes, 0);
    }

    /// F427: an ice the orbit's cycles move lays its caps by the traps its exchange answers
    /// them by, not its present season's: traps colder toward the north, its turned table
    /// colder toward the south, Mars's water lies in the south.
    #[test]
    fn an_ice_lays_its_caps_by_its_turned_traps() {
        let graph = crate::presets::preset("Mars").unwrap().to_graph().0;
        let traps = rows_traps(ice_of(&graph, "H2O"), |lat| -8.0 + 0.02 * (90.0 - lat));
        let grid = CubeSphere::new(64);
        let turned = [Turned {
            species: "H2O",
            table: rows_table(|lat| (-8.0 + 0.02 * (90.0 + lat)).exp()),
            levels: vec![1.0],
        }];
        let mut layers = LayerSet::new(grid);
        let t = water_laid(&graph, &mut layers, &traps, &turned);
        let south = (0..grid.cell_count())
            .filter(|&k| t[k] > 0.0)
            .filter(|&k| lat_lon(grid, k).0 < 0.0)
            .count();
        let all = t.iter().filter(|&&t| t > 0.0).count();
        println!("{south} of {all} covered cells in the south");
        assert!(all > 0 && south == all);
    }
}
