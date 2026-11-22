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

/// A finish's liquids as its maps draw them (F416, F417), each by `by_share`: the sea's share
/// of each cell, its cells' about it under the sea; the sea's mean depth over each cell, m, its
/// depth over the share it covers this over that share, a wet cell's own; the lakes' mean depth
/// over each cell, m; and the surface the maps are lit over, m, the sea's flat level over the
/// sea's share of each cell and the land's ground over the rest, a dry cell's own, a wet cell's
/// its dry cells' about it.
pub struct Drawn {
    pub sea_share: Vec<f32>,
    pub sea: Vec<f32>,
    pub lakes: Vec<f32>,
    pub surface: Vec<f32>,
}

impl Drawn {
    pub fn of(out: &FinishOutput) -> Self {
        let grid = out.grid;
        let n = grid.cell_count();
        let ground = out.height.as_slice();
        let wet = wet(out);
        let (sea_share, sea) = by_share(&grid, out.layers.liquid.as_slice());
        let lakes = match &out.rivers {
            Some(r) => by_share(&grid, r.lakes.as_slice()).1,
            None => vec![0.0; n],
        };
        let dry: Vec<f32> = (0..n)
            .map(|k| if wet[k] > 0.0 { 0.0 } else { ground[k] })
            .collect();
        let dry = over_cells(&grid, &dry);
        let level = out.sea_level as f32;
        let surface = (0..n)
            .into_par_iter()
            .map(|k| {
                let share = sea_share[k];
                let land = if wet[k] <= 0.0 {
                    ground[k]
                } else if share < 1.0 {
                    dry[k] / (1.0 - share)
                } else {
                    level
                };
                if share > 0.0 {
                    share * level + (1.0 - share) * land
                } else {
                    land
                }
            })
            .collect();
        Drawn {
            sea_share,
            sea,
            lakes,
            surface,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fields::Field;

    /// The cells about a cell.
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
