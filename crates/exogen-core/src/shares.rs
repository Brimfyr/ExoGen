// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! The maps by shares, not cells (F416, ADR 0008). A cell's value is its mean; the ground
//! within it and about it at the grid's own scale is more than its cell places. Where the
//! model decides a threshold cell by cell (its sea over its ground, a lake, an ice), drawing
//! each cell as all or none lets the relief at the cell's scale decide which side each cell
//! falls, and a margin where the ground hovers near the threshold speckles cell by cell. Drawn
//! instead with the cells about each cell read as samples of the ground within it, each kind
//! by the share of those samples it holds, the margin follows the ground's shape the grid
//! resolves, softened over as wide a band as the samples are read over: the physics places
//! it, the relief at the cell's scale realises its edge (PLAN principle 16). Each kind is drawn
//! over its own share of a cell: the sea's depth and its ice over the sea's share, the land's
//! ground and its ice over the rest (F417).

use rayon::prelude::*;

use crate::constants::RELIEF_SPREAD_CELLS;
use crate::grid::CubeSphere;
use crate::stages::FinishOutput;
use crate::stages::episodic::smooth;

/// A quantity's mean over the cells about each cell (F416), the cells read as samples of the
/// ground within it: over a Gaussian of RELIEF_SPREAD_CELLS cells' spread, by passes of the
/// grid's neighbour average, each spreading by a quarter of a cell squared. A kind's share of
/// each cell is its cells' mean as ones and noughts; a quantity of one kind's, over that
/// kind's share of a cell of the other kind, is its mean with the other kind's cells as
/// noughts, over that share (F417).
pub fn over_cells(grid: &CubeSphere, values: &[f32]) -> Vec<f32> {
    let passes = (RELIEF_SPREAD_CELLS * RELIEF_SPREAD_CELLS / 0.25).round() as usize;
    smooth(grid, values, passes)
}

/// Each cell under the finish's sea, 1, or not, 0: the model's own decision, cell by cell.
pub fn wet(out: &FinishOutput) -> Vec<f32> {
    wet_of(out.layers.liquid.as_slice())
}

/// Each cell its `liquid`, m, covers, 1, or not, 0.
fn wet_of(liquid: &[f32]) -> Vec<f32> {
    liquid
        .par_iter()
        .map(|&l| if l > 0.0 { 1.0 } else { 0.0 })
        .collect()
}

/// The least share a map draws (F417): half the 8-bit mask's least step, so that a quantity
/// drawn by its share stands above none exactly where its share does in the mask.
const LEAST_SHARE: f32 = 0.5 / 255.0;

/// The cells about a cell.
const AROUND: [(i32, i32); 8] = [
    (1, 0),
    (-1, 0),
    (0, 1),
    (0, -1),
    (1, 1),
    (1, -1),
    (-1, 1),
    (-1, -1),
];

/// Values a margin drawn by shares takes from the cells it is drawn from (F417): a cell `on`
/// keeps its own `fields`; one off them within `reach` takes those of the nearest cell on, ring
/// by ring of the cells about it, the first such in each ring; the rest none. So that the
/// timing of an ice's window, which no mean can blend round the year's end, and which ice
/// lies, reach the margin its share draws.
pub fn carried(grid: &CubeSphere, on: &[bool], reach: &[bool], fields: &[&[f32]]) -> Vec<Vec<f32>> {
    const NONE: u32 = u32::MAX;
    let n = on.len();
    let mut from: Vec<u32> = (0..n)
        .map(|k| if on[k] { k as u32 } else { NONE })
        .collect();
    loop {
        let next: Vec<u32> = (0..n)
            .into_par_iter()
            .map(|k| {
                if from[k] != NONE || !reach[k] {
                    return from[k];
                }
                let c = grid.cell(k);
                AROUND
                    .iter()
                    .map(|&(di, dj)| from[grid.index(grid.neighbor(c, di, dj))])
                    .find(|&f| f != NONE)
                    .unwrap_or(NONE)
            })
            .collect();
        let grew = next != from;
        from = next;
        if !grew {
            break;
        }
    }
    fields
        .iter()
        .map(|f| {
            from.iter()
                .map(|&s| if s == NONE { 0.0 } else { f[s as usize] })
                .collect()
        })
        .collect()
}

/// A quantity that lies on some cells and not others (a sea's depth, a lake's, a cap's
/// thickness), `values` above none where it lies, drawn by its share of each cell (F416, F417):
/// its share the mean over the cells about each cell of those it lies on as ones and the rest
/// as noughts, none under LEAST_SHARE; over that share, a cell's own where it lies on the cell
/// and the mean of the cells it lies on about it where not. Returns the share of each cell and
/// the quantity's mean over each cell, the share times that: a cell it covers whole keeps its
/// own, its depth or its thickness unblurred, and only a margin's cells are drawn by shares.
pub fn by_share(grid: &CubeSphere, values: &[f32]) -> (Vec<f32>, Vec<f32>) {
    let on = wet_of(values);
    let share = over_cells(grid, &on);
    let sum = over_cells(grid, values);
    (0..values.len())
        .into_par_iter()
        .map(|k| {
            if share[k] < LEAST_SHARE {
                return (0.0, 0.0);
            }
            let over = if on[k] > 0.0 {
                values[k]
            } else {
                sum[k] / share[k]
            };
            (share[k], share[k] * over)
        })
        .unzip()
}

/// A map of each cell's own value drawn so that its kinds agree with their margins drawn by
/// shares (F419, ADR 0009): of two kinds of cell, `kind` and `other` as ones and noughts, a
/// cell of one whose kind holds the lesser share of the cells about it takes the other's mean
/// of `values` about it; the rest, and cells of neither kind, keep their own. So that a map
/// drawn cell by cell, the ground's height, divides its kinds where their shares are even, as
/// the maps drawn by shares do, and only the cells the relief at the cells' scale set apart
/// from the ground about them change.
pub fn by_majority(grid: &CubeSphere, kind: &[f32], other: &[f32], values: &[f32]) -> Vec<f32> {
    let n = values.len();
    let mean = |flags: &[f32]| {
        let share = over_cells(grid, flags);
        let weighted: Vec<f32> = (0..n).map(|k| values[k] * flags[k]).collect();
        (share, over_cells(grid, &weighted))
    };
    let (kind_share, kind_sum) = mean(kind);
    let (other_share, other_sum) = mean(other);
    (0..n)
        .into_par_iter()
        .map(|k| {
            if kind[k] > 0.0 && kind_share[k] < other_share[k] {
                other_sum[k] / other_share[k]
            } else if other[k] > 0.0 && other_share[k] < kind_share[k] {
                kind_sum[k] / kind_share[k]
            } else {
                values[k]
            }
        })
        .collect()
}

/// Values of kinds that cannot be blended, a unit's id, drawn as `by_majority` draws a cell's
/// own (F419): a cell of one kind holding the lesser share of the cells about it takes the
/// nearest cell's of the other kind (`carried`), its own where none is reached.
pub fn by_majority_nearest(
    grid: &CubeSphere,
    kind: &[f32],
    other: &[f32],
    values: &[f32],
) -> Vec<f32> {
    let n = values.len();
    let (kind_share, other_share) = (over_cells(grid, kind), over_cells(grid, other));
    let lesser = |of: &[f32], own: &[f32], theirs: &[f32]| -> Vec<bool> {
        (0..n).map(|k| of[k] > 0.0 && own[k] < theirs[k]).collect()
    };
    let kind_lesser = lesser(kind, &kind_share, &other_share);
    let other_lesser = lesser(other, &other_share, &kind_share);
    let ones = vec![1.0f32; n];
    // From the other kind's cells on their greater side to the kind's on its lesser, and back.
    let from = |of: &[f32], skip: &[bool], to: &[bool]| {
        let on: Vec<bool> = (0..n).map(|k| of[k] > 0.0 && !skip[k]).collect();
        carried(grid, &on, to, &[values, &ones])
    };
    let to_kind = from(other, &other_lesser, &kind_lesser);
    let to_other = from(kind, &kind_lesser, &other_lesser);
    (0..n)
        .map(|k| {
            let taken = if kind_lesser[k] {
                &to_kind
            } else if other_lesser[k] {
                &to_other
            } else {
                return values[k];
            };
            if taken[1][k] > 0.0 {
                taken[0][k]
            } else {
                values[k]
            }
        })
        .collect()
}

/// The rivers as the maps draw them (F421, ADR 0009), each cell's `discharge`, m³/s, none under
/// the water, drawn where the land is drawn over more than half of each cell, `share` the sea's
/// and the lakes' shares of each cell together and `water` 1 on the cells under them: a river
/// cell the water is drawn over half of or more is cut, the river reaching the water before it;
/// and a river cell beside the water carries its river on through the water's cells, to the
/// one beside it the water is drawn over most of, until it reaches one the water is drawn over
/// half of, so that a river meets the water where the water is drawn and not a cell or two
/// short of it, at a narrow inlet's head or a lake's outlet.
pub fn rivers(grid: &CubeSphere, discharge: &[f32], water: &[f32], share: &[f32]) -> Vec<f32> {
    let n = discharge.len();
    let mut drawn: Vec<f32> = (0..n)
        .map(|k| {
            if water[k] <= 0.0 && share[k] < 0.5 {
                discharge[k]
            } else {
                0.0
            }
        })
        .collect();
    for k in 0..n {
        let q = discharge[k];
        if q <= 0.0 || water[k] > 0.0 || share[k] >= 0.5 {
            continue;
        }
        let mut at = k;
        for _ in 0..RIVER_CARRY_CELLS {
            let c = grid.cell(at);
            let next = AROUND
                .iter()
                .map(|&(di, dj)| grid.index(grid.neighbor(c, di, dj)))
                .filter(|&j| water[j] > 0.0 && share[j] > share[at])
                .max_by(|&a, &b| share[a].total_cmp(&share[b]).then(b.cmp(&a)));
            let Some(j) = next else { break };
            if share[j] >= 0.5 {
                break;
            }
            drawn[j] = drawn[j].max(q);
            at = j;
        }
    }
    drawn
}

/// The most cells a river is carried on through the water's cells at a margin (F421): eight of
/// the shares' Gaussian's spreads, RELIEF_SPREAD_CELLS cells each, past which it holds nothing.
const RIVER_CARRY_CELLS: usize = 16;

/// The ices that come and go through the year as the maps draw them (F417), from each cell's
/// window, `window` (the phases it freezes at, melts out at and peaks at), its thickness at its
/// peak, `thickness`, m, and which ice lies, `ice`: the share of each cell the cells an ice
/// lies on cover; the window and the ice of a cell an ice lies on its own, and of a cell at
/// their margin the nearest such cell's (`carried`); the thickness over the share a cell's own
/// where an ice lies on it and the mean of those about it where not (`by_share`); none past
/// the share.
pub struct Seasons {
    pub share: Vec<f32>,
    pub window: [Vec<f32>; 3],
    pub thickness: Vec<f32>,
    pub ice: Vec<f32>,
}

impl Seasons {
    pub fn of(grid: &CubeSphere, window: [&[f32]; 3], thickness: &[f32], ice: &[f32]) -> Self {
        let (share, mean) = by_share(grid, thickness);
        let on: Vec<bool> = thickness.iter().map(|&t| t > 0.0).collect();
        let reach: Vec<bool> = share.iter().map(|&s| s > 0.0).collect();
        let [freeze, melt, peak] = window;
        let mut kept = carried(grid, &on, &reach, &[freeze, melt, peak, ice]);
        let ice = kept.pop().expect("the ices' index");
        let peak = kept.pop().expect("the peak");
        let melt = kept.pop().expect("the melt");
        let freeze = kept.pop().expect("the freeze");
        let thickness = (0..thickness.len())
            .into_par_iter()
            .map(|k| {
                if on[k] {
                    thickness[k]
                } else if share[k] > 0.0 {
                    mean[k] / share[k]
                } else {
                    0.0
                }
            })
            .collect();
        Seasons {
            share,
            window: [freeze, melt, peak],
            thickness,
            ice,
        }
    }
}

/// A quantity of one kind of cell over that kind's share of each cell (F417, F420), `kind` 1 on
/// the kind's cells and 0 on the rest: a cell's own `values` where the cell is of the kind, and
/// where not the mean of the kind's cells' about it, each read as a sample of the ground within
/// it; its own where none of the kind lies about it.
pub fn of_kind(grid: &CubeSphere, values: &[f32], kind: &[f32]) -> Vec<f32> {
    let share = over_cells(grid, kind);
    let weighted: Vec<f32> = values.iter().zip(kind).map(|(v, w)| v * w).collect();
    let sum = over_cells(grid, &weighted);
    (0..values.len())
        .into_par_iter()
        .map(|k| {
            if kind[k] > 0.0 || share[k] < LEAST_SHARE {
                values[k]
            } else {
                sum[k] / share[k]
            }
        })
        .collect()
}

/// The shape the grid resolves of a quantity over one kind of cell (F423), `kind` 1 on the
/// kind's cells and 0 on the rest: at a cell of the kind, the mean of the kind's cells about
/// it, each read as a sample of the ground within it (the land's ground at a coast not drawn
/// down by the sea's floor beside it); at another, its own.
pub fn shape_of(grid: &CubeSphere, values: &[f32], kind: &[f32]) -> Vec<f32> {
    let share = over_cells(grid, kind);
    let weighted: Vec<f32> = values.iter().zip(kind).map(|(v, w)| v * w).collect();
    let sum = over_cells(grid, &weighted);
    (0..values.len())
        .into_par_iter()
        .map(|k| {
            if kind[k] > 0.0 && share[k] > 0.0 {
                sum[k] / share[k]
            } else {
                values[k]
            }
        })
        .collect()
}

/// A finish's liquids and ground as its maps draw them (F416, F417, F419, F420), each liquid by
/// `by_share`. Each cell's kind, 1 or 0: under the sea (`wet`), under a lake (`lake`), or land,
/// neither (`land`). The sea's share of each cell and its mean depth over each cell, m, its
/// depth over the share it covers this over that share, a wet cell's own; the lakes' share and
/// mean depth likewise. The surface the maps are lit over, m: the sea's flat level over the
/// sea's share of each cell, each lake's flat surface over the lakes' share, a lake cell's own
/// and the mean of the lake cells' about it, and the land's ground over the rest, a land cell's
/// own, another's the land cells' about it. And the ground's height as the maps draw it, m, a
/// cell's own but for a cell on the lesser side of the cells about it, between the sea and the
/// land and then between a lake and the land, drawn at the greater side's mean (F419,
/// `by_majority`): so that the ground crosses the sea's level where the sea's share is a half,
/// as the mask does, and a lake's where the lakes' is (`agreeing` draws other such maps alike).
/// And the rivers, m³/s, drawn where the land is drawn over more than half of each cell (F421,
/// `rivers`), none without them.
pub struct Drawn {
    pub grid: CubeSphere,
    pub rivers: Vec<f32>,
    pub wet: Vec<f32>,
    pub lake: Vec<f32>,
    pub land: Vec<f32>,
    pub sea_share: Vec<f32>,
    pub sea: Vec<f32>,
    pub lake_share: Vec<f32>,
    pub lakes: Vec<f32>,
    pub surface: Vec<f32>,
    pub height: Vec<f32>,
}

impl Drawn {
    /// A map of each cell's own `values` of the kinds of cell, the sea's floor and the land
    /// (the temperature, the surface's age), drawn as the height is (F419): a cell on the
    /// lesser side of the cells about it, between the sea and the land and then between a lake
    /// and the land, at the greater side's mean; the rest their own.
    pub fn agreeing(&self, values: &[f32]) -> Vec<f32> {
        let dry: Vec<f32> = self.wet.iter().map(|&w| 1.0 - w).collect();
        let values = by_majority(&self.grid, &self.wet, &dry, values);
        by_majority(&self.grid, &self.lake, &self.land, &values)
    }

    /// As `agreeing`, for values that cannot be blended, a unit's id: a cell on the lesser
    /// side takes the nearest cell's of the greater (`by_majority_nearest`).
    pub fn agreeing_nearest(&self, values: &[f32]) -> Vec<f32> {
        let dry: Vec<f32> = self.wet.iter().map(|&w| 1.0 - w).collect();
        let values = by_majority_nearest(&self.grid, &self.wet, &dry, values);
        by_majority_nearest(&self.grid, &self.lake, &self.land, &values)
    }

    pub fn of(out: &FinishOutput) -> Self {
        let grid = out.grid;
        let n = grid.cell_count();
        let ground = out.height.as_slice();
        let wet = wet(out);
        let (sea_share, sea) = by_share(&grid, out.layers.liquid.as_slice());
        let depth = match &out.rivers {
            Some(r) => r.lakes.as_slice().to_vec(),
            None => vec![0.0; n],
        };
        let lake = wet_of(&depth);
        let (lake_share, lakes) = by_share(&grid, &depth);
        let land: Vec<f32> = (0..n)
            .map(|k| {
                if wet[k] > 0.0 || lake[k] > 0.0 {
                    0.0
                } else {
                    1.0
                }
            })
            .collect();
        let dry: Vec<f32> = wet.iter().map(|&w| 1.0 - w).collect();
        let height = by_majority(&grid, &wet, &dry, ground);
        let height = by_majority(&grid, &lake, &land, &height);
        let on_land = of_kind(&grid, ground, &land);
        let lake_level: Vec<f32> = ground.iter().zip(&depth).map(|(g, d)| g + d).collect();
        let lake_level = of_kind(&grid, &lake_level, &lake);
        let level = out.sea_level as f32;
        let surface = (0..n)
            .into_par_iter()
            .map(|k| {
                let (s, l) = (sea_share[k], lake_share[k]);
                let rest = (1.0 - s - l).max(0.0);
                let mut h = rest * on_land[k];
                if s > 0.0 {
                    h += s * level;
                }
                if l > 0.0 {
                    h += l * lake_level[k];
                }
                h / (s + l + rest).max(f32::MIN_POSITIVE)
            })
            .collect();
        let rivers = match &out.rivers {
            Some(r) => {
                let water: Vec<f32> = (0..n).map(|k| wet[k].max(lake[k])).collect();
                let share: Vec<f32> = (0..n).map(|k| sea_share[k] + lake_share[k]).collect();
                rivers(&grid, r.discharge.as_slice(), &water, &share)
            }
            None => vec![0.0; n],
        };
        Drawn {
            grid,
            rivers,
            wet,
            lake,
            land,
            sea_share,
            sea,
            lake_share,
            lakes,
            surface,
            height,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fields::Field;

    /// A sea's depth over each cell, m, its level at none over the `ground`.
    fn under(ground: &Field<f32>) -> Vec<f32> {
        ground.as_slice().iter().map(|&g| (-g).max(0.0)).collect()
    }

    /// The cells of face 0 whose Gaussian keeps to the face, and each one's distance from the
    /// face's middle row, cells.
    fn inner(grid: &CubeSphere) -> impl Fn(usize) -> Option<f64> + '_ {
        let f = grid.face_size();
        let cell = grid.delta();
        move |k| {
            let c = grid.cell(k);
            (c.face == 0 && (20..f - 20).contains(&c.i) && (20..f - 20).contains(&c.j))
                .then(|| grid.center(c)[1] / cell)
        }
    }

    /// F416: a coast the ground's shape places and its relief at the cells' scale speckles:
    /// over a plain rising through the sea's level, each cell's ground spread at random by
    /// 100 m, the cells drawn all or none speckle the coast with cells cut off from their kind
    /// on either side; drawn by their shares, half the cells' share lies where the plain itself
    /// crosses the level, within half a cell on the mean and a cell and a half about it, and
    /// hardly a cell stands apart from its kind.
    #[test]
    fn a_coast_is_placed_by_the_ground_s_shape_not_its_relief() {
        let grid = CubeSphere::new(128);
        let n = grid.cell_count();
        let rng = crate::rng::Rng::new(5);
        let rise = 1_500.0;
        let ground = Field::from_fn(grid, |c, d| {
            let noise = 100.0 * (rng.uniform(grid.index(c) as u64, 0) - 0.5) * 12f64.sqrt();
            (rise * d[1] + noise) as f32
        });
        let liquid = under(&ground);
        let wet = wet_of(&liquid);
        let (share, _) = by_share(&grid, &liquid);
        let inner = inner(&grid);
        let apart = |wet: &dyn Fn(usize) -> bool| {
            (0..n)
                .filter(|&k| inner(k).is_some())
                .filter(|&k| {
                    let c = grid.cell(k);
                    let alike = AROUND
                        .iter()
                        .filter(|&&(di, dj)| wet(grid.index(grid.neighbor(c, di, dj))) == wet(k))
                        .count();
                    alike <= 1
                })
                .count()
        };
        let (cells, shares) = (apart(&|k| wet[k] > 0.5), apart(&|k| share[k] >= 0.5));
        // Where half the share lies against where the plain crosses the level, in cells: on
        // the mean, and how far about it.
        let (mut off, mut wide, mut count) = (0.0, 0.0, 0);
        for (k, s) in share.iter().enumerate() {
            let Some(y) = inner(k) else { continue };
            if (s - 0.5).abs() < 0.15 {
                off += y;
                wide += y.abs();
                count += 1;
            }
        }
        let (off, wide) = (off / count.max(1) as f64, wide / count.max(1) as f64);
        println!(
            "cells apart from their kind: {cells} drawn all or none, {shares} by shares; half the share {off:.2} cells from where the plain crosses the level on the mean, {wide:.2} about it ({count} cells)"
        );
        assert!(cells > 50, "{cells}");
        assert!(shares <= cells / 20, "{shares} {cells}");
        assert!(count > 0 && off.abs() < 0.5 && wide < 1.5, "{off} {wide}");
    }

    /// F417: a steep coast the grid resolves is drawn where its cells meet: a deep sea, 3,000 m,
    /// against land whose ridges stand 2,500 m over valleys 100 m over the level every four
    /// cells, its first ridge on the coast. Half the sea's share lies within half a cell of the
    /// line between its last wet cells and its first dry ones, not drawn up the land where the
    /// sea's depth outweighs the land's height; and the valleys five cells and more inland,
    /// the ridges' relief about them none of the sea's, carry under a hundredth of its share.
    #[test]
    fn a_steep_coast_is_drawn_where_its_cells_meet() {
        let grid = CubeSphere::new(128);
        let inner = inner(&grid);
        let cell = grid.delta();
        let ground = Field::from_fn(grid, |c, _| {
            let y = grid.center(c)[1] / cell;
            if y < 0.0 {
                -3_000.0
            } else {
                (1_300.0 + 1_200.0 * (std::f64::consts::TAU * y / 4.0).cos()) as f32
            }
        });
        let (share, _) = by_share(&grid, &under(&ground));
        let (mut off, mut count, mut inland) = (0.0, 0, 0.0f32);
        for (k, &s) in share.iter().enumerate() {
            let Some(y) = inner(k) else { continue };
            if (s - 0.5).abs() < 0.15 {
                off += y;
                count += 1;
            }
            if y >= 5.0 {
                inland = inland.max(s);
            }
        }
        let off = off / count.max(1) as f64;
        println!(
            "half the share {off:.2} cells from where the cells meet ({count} cells); the most share five cells inland and more {inland:.4}"
        );
        assert!(count > 0 && off.abs() < 0.5, "{off}");
        assert!(inland < 0.01, "{inland}");
    }

    /// F419: the ground's height drawn so that it crosses the sea's level where the sea's share
    /// is a half: over the plain the relief at the cells' scale speckles, each cell's drawn
    /// ground stands under the level exactly where more than half the cells about it are under
    /// the sea; the cells that change are those the relief set on the lesser side, each drawn
    /// at the greater side's mean, and every other keeps its own height.
    #[test]
    fn the_ground_crosses_the_level_where_the_sea_s_share_is_a_half() {
        let grid = CubeSphere::new(128);
        let n = grid.cell_count();
        let rng = crate::rng::Rng::new(5);
        let ground = Field::from_fn(grid, |c, d| {
            let noise = 100.0 * (rng.uniform(grid.index(c) as u64, 0) - 0.5) * 12f64.sqrt();
            (1_500.0 * d[1] + noise) as f32
        });
        let ground = ground.as_slice();
        let wet: Vec<f32> = ground
            .iter()
            .map(|&g| if g < 0.0 { 1.0 } else { 0.0 })
            .collect();
        let land: Vec<f32> = wet.iter().map(|&w| 1.0 - w).collect();
        let drawn = by_majority(&grid, &wet, &land, ground);
        let share = over_cells(&grid, &wet);
        let inner = inner(&grid);
        let (mut disagree, mut changed, mut stray, mut kept) = (0, 0, 0, 0);
        for k in (0..n).filter(|&k| inner(k).is_some()) {
            if (share[k] - 0.5).abs() < 1e-6 {
                continue;
            }
            let sea = share[k] > 0.5;
            disagree += usize::from((drawn[k] < 0.0) != sea);
            if drawn[k] != ground[k] {
                changed += 1;
                // Only a cell its own ground set on the lesser side changes.
                stray += usize::from((ground[k] < 0.0) == sea);
            } else {
                kept += usize::from((ground[k] < 0.0) == sea);
            }
        }
        println!(
            "{disagree} cells' drawn ground on the other side of the level from the sea's half share; {changed} cells drawn at the greater side's mean, {stray} of them on it already; {kept} kept their own"
        );
        assert_eq!(disagree, 0);
        assert!(changed > 50, "{changed}");
        assert_eq!(stray, 0);
    }

    /// F419: a unit's id, which no mean can blend, drawn as the height is: over the plain the
    /// relief at the cells' scale speckles, the floor's unit 7 and the land's 3, each cell's
    /// drawn unit is the floor's exactly where more than half the cells about it are under the
    /// sea, a cell set on the lesser side taking its nearest cell's on the greater.
    #[test]
    fn a_unit_is_drawn_where_its_kind_holds_the_greater_share() {
        let grid = CubeSphere::new(128);
        let n = grid.cell_count();
        let rng = crate::rng::Rng::new(5);
        let wet: Vec<f32> = (0..n)
            .map(|k| {
                let d = grid.center(grid.cell(k));
                let noise = 100.0 * (rng.uniform(k as u64, 0) - 0.5) * 12f64.sqrt();
                if 1_500.0 * d[1] + noise < 0.0 {
                    1.0
                } else {
                    0.0
                }
            })
            .collect();
        let land: Vec<f32> = wet.iter().map(|&w| 1.0 - w).collect();
        let ids: Vec<f32> = wet
            .iter()
            .map(|&w| if w > 0.0 { 7.0 } else { 3.0 })
            .collect();
        let drawn = by_majority_nearest(&grid, &wet, &land, &ids);
        let share = over_cells(&grid, &wet);
        let inner = inner(&grid);
        let (mut disagree, mut changed) = (0, 0);
        for k in (0..n).filter(|&k| inner(k).is_some()) {
            if (share[k] - 0.5).abs() < 1e-6 {
                continue;
            }
            disagree += usize::from((drawn[k] == 7.0) != (share[k] > 0.5));
            changed += usize::from(drawn[k] != ids[k]);
        }
        println!(
            "{disagree} cells' units on the lesser side; {changed} cells drawn with their nearest's"
        );
        assert_eq!(disagree, 0);
        assert!(changed > 50, "{changed}");
    }

    /// F421: a river meets the water where the water is drawn. A river down the land to the
    /// head of an inlet two cells wide and eight long, the sea beyond, the inlet drawn under
    /// half its cells: the river is carried down it to where the sea is drawn over half, where
    /// as the model leaves it it stops at the inlet's head, cells short of the sea. And a river
    /// along a spit a cell wide out into the sea, drawn as the sea, is cut: no river is drawn
    /// where the water is drawn over half of a cell.
    #[test]
    fn a_river_meets_the_water_where_it_is_drawn() {
        let grid = CubeSphere::new(128);
        let n = grid.cell_count();
        let cell = grid.delta();
        let at = |k: usize| {
            let c = grid.cell(k);
            let d = grid.center(c);
            (c.face == 0, d[1] / cell, d[2] / cell)
        };
        let (inlet, river, spit) = (
            |k: usize| {
                let (face, y, x) = at(k);
                face && (0.0..8.0).contains(&y) && x.abs() < 1.0
            },
            |k: usize| {
                let (face, y, x) = at(k);
                face && (8.0..30.0).contains(&y) && (0.0..1.0).contains(&x)
            },
            |k: usize| {
                let (face, y, x) = at(k);
                face && (-6.0..0.0).contains(&y) && (10.0..11.0).contains(&x)
            },
        );
        let wet: Vec<f32> = (0..n)
            .map(|k| {
                let (_, y, _) = at(k);
                if (y < 0.0 && !spit(k)) || inlet(k) {
                    1.0
                } else {
                    0.0
                }
            })
            .collect();
        let discharge: Vec<f32> = (0..n)
            .map(|k| if river(k) || spit(k) { 1_000.0 } else { 0.0 })
            .collect();
        let share = over_cells(&grid, &wet);
        let drawn = rivers(&grid, &discharge, &wet, &share);
        let meets = |q: &[f32]| {
            (0..n).any(|k| {
                q[k] > 0.0
                    && AROUND.iter().any(|&(di, dj)| {
                        share[grid.index(grid.neighbor(grid.cell(k), di, dj))] >= 0.5
                    })
            })
        };
        let under = (0..n)
            .filter(|&k| drawn[k] > 0.0 && share[k] >= 0.5)
            .count();
        let carried = (0..n).filter(|&k| drawn[k] > 0.0 && wet[k] > 0.0).count();
        let cut = (0..n).filter(|&k| spit(k) && drawn[k] <= 0.0).count();
        // The inland river alone, as the model leaves it and as it is drawn.
        let inland =
            |q: &[f32]| -> Vec<f32> { (0..n).map(|k| if spit(k) { 0.0 } else { q[k] }).collect() };
        let (model, shown) = (meets(&inland(&discharge)), meets(&inland(&drawn)));
        println!(
            "the model's river meets the drawn sea: {model}; drawn: {shown}, carried through {carried} of the inlet's cells; {under} cells drawn under the water; {cut} of the spit's cells cut"
        );
        assert!(!model);
        assert!(shown);
        assert!(carried >= 3, "{carried}");
        assert_eq!(under, 0);
        assert!(cut > 0, "{cut}");
    }

    /// F420: a quantity of the land read over the land's share: flats wet and dry at random, the
    /// land's aridity 0.3 and the open water's 2. Over the land's share every cell reads the
    /// land's 0.3, a wet cell its land cells' about it; read cell by cell, half the flats' cells
    /// gave the land's share of a pixel the water's.
    #[test]
    fn the_land_s_share_reads_the_land_s_values() {
        let grid = CubeSphere::new(64);
        let n = grid.cell_count();
        let rng = crate::rng::Rng::new(11);
        let land: Vec<f32> = (0..n)
            .map(|k| {
                if rng.uniform(k as u64, 0) < 0.5 {
                    1.0
                } else {
                    0.0
                }
            })
            .collect();
        let values: Vec<f32> = land
            .iter()
            .map(|&l| if l > 0.0 { 0.3 } else { 2.0 })
            .collect();
        let read = of_kind(&grid, &values, &land);
        let worst = read.iter().fold(0.0f32, |a, &v| a.max((v - 0.3).abs()));
        let water = values.iter().filter(|&&v| v > 1.0).count();
        println!(
            "the land's share reads off the land's by {worst:.1e} at most; {water} of {n} cells hold the water's own"
        );
        assert!(worst < 1e-5, "{worst}");
        assert!(water > n / 3, "{water}");
    }

    /// F423: the land's shape the grid resolves at a coast: land 100 m over the level, each
    /// cell's ground spread at random by 100 m, against a sea 3,000 m deep. The land's shape is
    /// the mean of the land cells about each, within a few tens of metres of the land's 100 m
    /// to the coast, not drawn down toward the sea's floor as the mean of every cell about it
    /// is; the sea's cells keep their own.
    #[test]
    fn the_land_s_shape_is_its_own_cells_mean() {
        let grid = CubeSphere::new(128);
        let n = grid.cell_count();
        let rng = crate::rng::Rng::new(13);
        let cell = grid.delta();
        let ground: Vec<f32> = (0..n)
            .map(|k| {
                if grid.center(grid.cell(k))[1] / cell < 0.0 {
                    -3_000.0
                } else {
                    (100.0 + 100.0 * (rng.uniform(k as u64, 0) - 0.5) * 12f64.sqrt()) as f32
                }
            })
            .collect();
        let land: Vec<f32> = ground
            .iter()
            .map(|&g| if g > -1_000.0 { 1.0 } else { 0.0 })
            .collect();
        let shape = shape_of(&grid, &ground, &land);
        let all = over_cells(&grid, &ground);
        let inner = inner(&grid);
        let (mut off, mut pulled, mut sea) = (0.0f32, 0.0f32, 0.0f32);
        for k in (0..n).filter(|&k| inner(k).is_some()) {
            if land[k] > 0.0 {
                off = off.max((shape[k] - 100.0).abs());
                pulled = pulled.max(100.0 - all[k]);
            } else {
                sea = sea.max((shape[k] - ground[k]).abs());
            }
        }
        println!(
            "the land's shape off its 100 m by {off:.0} m at most, the mean of every cell about it by {pulled:.0}; the sea's cells off their own by {sea}"
        );
        assert!(off < 60.0, "{off}");
        assert!(pulled > 500.0, "{pulled}");
        assert_eq!(sea, 0.0);
    }

    /// F419: a steep coast's ground is drawn as it stands: where the cells meet along a coast
    /// the grid resolves, no cell stands on the lesser side of the cells about it, and none
    /// changes.
    #[test]
    fn a_steep_coast_s_ground_is_drawn_as_it_stands() {
        let grid = CubeSphere::new(128);
        let cell = grid.delta();
        let ground = Field::from_fn(grid, |c, _| {
            let y = grid.center(c)[1] / cell;
            if y < 0.0 {
                -3_000.0
            } else {
                (1_300.0 + 1_200.0 * (std::f64::consts::TAU * y / 4.0).cos()) as f32
            }
        });
        let wet = wet_of(&under(&ground));
        let land: Vec<f32> = wet.iter().map(|&w| 1.0 - w).collect();
        let drawn = by_majority(&grid, &wet, &land, ground.as_slice());
        let inner = inner(&grid);
        let changed = (0..grid.cell_count())
            .filter(|&k| inner(k).is_some() && drawn[k] != ground.as_slice()[k])
            .count();
        println!("{changed} cells of a steep coast drawn other than they stand");
        assert_eq!(changed, 0);
    }

    /// F417: a sea's shallows keep their depth: a sea 2,000 m deep off a coast, a shelf 10 m
    /// deep along it two cells wide, and a bank 5 m deep a cell wide ten cells out. Each cell
    /// the sea covers is drawn at its own depth over the share it covers, the bank at its 5 m
    /// and the shelf at its 10 m, where the depth's mean over the cells about them drew the bank
    /// hundreds of metres deep and the shelf's water dark. And the sea's depth stands above
    /// none exactly where its share does, none eight cells and more inland.
    #[test]
    fn a_sea_s_shallows_keep_their_depth() {
        let grid = CubeSphere::new(128);
        let inner = inner(&grid);
        let cell = grid.delta();
        let ground = Field::from_fn(grid, |c, _| {
            let y = grid.center(c)[1] / cell;
            if y >= 0.0 {
                500.0
            } else if y > -2.0 {
                -10.0
            } else if (-11.0..-10.0).contains(&y) {
                -5.0
            } else {
                -2_000.0
            }
        });
        let liquid = under(&ground);
        let (share, mean) = by_share(&grid, &liquid);
        let blurred = over_cells(&grid, &liquid);
        let (mut bank, mut shelf, mut was, mut cells) = (0.0f32, 0.0f32, f32::INFINITY, 0);
        let mut halo = 0;
        for (k, (&s, &m)) in share.iter().zip(&mean).enumerate() {
            if (m > 0.0) != (s > 0.0) {
                halo += 1;
            }
            let Some(y) = inner(k) else { continue };
            if (-11.0..-10.0).contains(&y) {
                bank = bank.max((m / s - 5.0).abs());
                was = was.min(blurred[k]);
                cells += 1;
            } else if (-2.0..0.0).contains(&y) {
                shelf = shelf.max((m / s - 10.0).abs());
            }
            if y > 8.0 {
                halo += usize::from(s > 0.0);
            }
        }
        println!(
            "the bank off its 5 m by {bank:.1e} m, the shelf off its 10 m by {shelf:.1e} m ({cells} bank cells), the bank by the depth's mean {was:.0} m at the shallowest; {halo} cells' depth and share apart or drawn far inland"
        );
        assert!(cells > 0 && bank < 1e-3 && shelf < 1e-3, "{bank} {shelf}");
        assert!(was > 100.0, "{was}");
        assert_eq!(halo, 0);
    }

    /// F417: a margin takes the values of the nearest cell it is drawn from: half a face's
    /// cells on, each holding its own column's value, the three rows beyond them within reach.
    /// Each cell on keeps its own; each in reach takes its own column's, the nearest cell on's,
    /// ring by ring; those out of reach none.
    #[test]
    fn a_margin_takes_its_nearest_cell_s_values() {
        let grid = CubeSphere::new(64);
        let n = grid.cell_count();
        let f = grid.face_size();
        let rows = |k: usize, to: u32| {
            let c = grid.cell(k);
            c.face == 0 && c.i < to
        };
        let on: Vec<bool> = (0..n).map(|k| rows(k, f / 2)).collect();
        let reach: Vec<bool> = (0..n).map(|k| rows(k, f / 2 + 3)).collect();
        let column = |k: usize| grid.cell(k).j as f32 + 1.0;
        let field: Vec<f32> = (0..n)
            .map(|k| if on[k] { column(k) } else { 0.0 })
            .collect();
        let got = &carried(&grid, &on, &reach, &[&field])[0];
        let (mut off, mut margin, mut beyond) = (0.0f32, 0, 0.0f32);
        for (k, &g) in got.iter().enumerate() {
            let c = grid.cell(k);
            let inside = c.face == 0 && (2..f - 2).contains(&c.j);
            if on[k] || (reach[k] && inside) {
                off = off.max((g - column(k)).abs());
                margin += usize::from(!on[k]);
            } else if !reach[k] {
                beyond = beyond.max(g.abs());
            }
        }
        println!(
            "off their own column's value by {off} ({margin} margin cells); beyond reach {beyond}"
        );
        assert!(margin > 0 && off == 0.0, "{off}");
        assert_eq!(beyond, 0.0);
    }

    /// F417: a snowfield that comes and goes, drawn by its share: a disc of cells whose snow
    /// freezes at 0.2 of the year, peaks at 0.4 and melts out at 0.7, 1 to 2 m at its peak.
    /// Its cells keep their own; the cells at its margin carry its window and its ice, and the
    /// snow's mean thickness over the share they hold, between 1 and 2 m, where the share times
    /// that would thin it to nothing as the share falls; past the share, none.
    #[test]
    fn a_snowfield_s_margin_carries_its_seasons() {
        let grid = CubeSphere::new(64);
        let n = grid.cell_count();
        let f = grid.face_size() as f64;
        let distance = |k: usize| {
            let c = grid.cell(k);
            let (i, j) = (c.i as f64 - f / 2.0, c.j as f64 - f / 2.0);
            if c.face == 0 {
                (i * i + j * j).sqrt()
            } else {
                f64::INFINITY
            }
        };
        let on: Vec<bool> = (0..n).map(|k| distance(k) < 10.0).collect();
        let at = |v: f32| -> Vec<f32> { on.iter().map(|&o| if o { v } else { 0.0 }).collect() };
        let (freeze, melt, peak, ice) = (at(0.2), at(0.7), at(0.4), at(1.0));
        let thickness: Vec<f32> = (0..n)
            .map(|k| {
                if on[k] {
                    1.0 + grid.cell(k).j as f32 / f as f32
                } else {
                    0.0
                }
            })
            .collect();
        let s = Seasons::of(&grid, [&freeze, &melt, &peak], &thickness, &ice);
        let (mut own, mut margin, mut thin, mut thick, mut beyond) =
            (0.0f32, 0, 2.0f32, 0.0f32, 0.0f32);
        for k in 0..n {
            let got = [s.window[0][k], s.window[1][k], s.window[2][k], s.ice[k]];
            if on[k] {
                own = own.max((s.thickness[k] - thickness[k]).abs());
            } else if s.share[k] > 0.0 {
                margin += 1;
                own = own.max(
                    (got[0] - 0.2).abs()
                        + (got[1] - 0.7).abs()
                        + (got[2] - 0.4).abs()
                        + (got[3] - 1.0).abs(),
                );
                thin = thin.min(s.thickness[k]);
                thick = thick.max(s.thickness[k]);
            } else {
                beyond = beyond.max(got.iter().fold(s.thickness[k], |a, v| a.max(v.abs())));
            }
        }
        println!(
            "its own and its margin's window off by {own}; {margin} margin cells {thin:.2} to {thick:.2} m; past the share {beyond}"
        );
        assert!(margin > 0 && own == 0.0, "{own}");
        assert!(thin >= 1.0 && thick <= 2.0, "{thin} {thick}");
        assert_eq!(beyond, 0.0);
    }
}
