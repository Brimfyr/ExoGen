// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! The maps by shares, not cells (F416, ADR 0008). A cell's ground is its mean; the relief
//! within it and about it at the grid's own scale is more than its cell places. Where a map
//! draws a threshold (a liquid's level over the ground, a lake's), drawing each cell as all
//! or none lets that relief decide which side of the threshold each cell falls, and a margin
//! where the ground hovers near the threshold speckles cell by cell. Drawn instead by the
//! share of the cell past the threshold, its ground spread about the relief's mean as the
//! relief at the cell's scale spreads it, the margin follows the ground's shape the grid
//! resolves, softened over as wide a band as the relief at the cell's scale makes uncertain:
//! the physics places it, the relief at the cell's scale realises its edge (PLAN principle 16).

use rayon::prelude::*;

use crate::constants::RELIEF_SPREAD_CELLS;
use crate::fields::Field;
use crate::grid::CubeSphere;
use crate::params::derive::physics::normal_cdf;
use crate::stages::FinishOutput;
use crate::stages::episodic::smooth;

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

/// The least spread a cell's ground is read with, m: a numerical floor, so that a ground the
/// grid places exactly draws its threshold as a step.
const LEAST_SPREAD: f64 = 1e-3;

/// The ground's relief at its cells' scale (F416): at each cell its mean over the cells about
/// it, m, and the spread of the ground about that mean, m, the root mean square of the ground
/// less its mean over the same cells: both read over a Gaussian of RELIEF_SPREAD_CELLS cells'
/// spread, by passes of the grid's neighbour average, each spreading by a quarter of a cell
/// squared. An even slope is its own mean and spreads by none.
pub struct Relief {
    pub mean: Vec<f32>,
    pub spread: Vec<f32>,
}

impl Relief {
    pub fn of(ground: &Field<f32>) -> Self {
        let grid = ground.grid();
        let h = ground.as_slice();
        let passes = (RELIEF_SPREAD_CELLS * RELIEF_SPREAD_CELLS / 0.25).round() as usize;
        let mean = smooth(&grid, h, passes);
        let square: Vec<f32> = h
            .par_iter()
            .zip(&mean)
            .map(|(&g, &m)| (g - m) * (g - m))
            .collect();
        let spread = smooth(&grid, &square, passes)
            .into_par_iter()
            .map(|v| v.max(0.0).sqrt())
            .collect();
        Relief { mean, spread }
    }

    /// The share of cell `k` under a liquid whose surface stands at `level`, m, and the
    /// liquid's mean depth over the cell, m: its ground spread about its relief's mean by its
    /// spread as a normal distribution, the share the distribution's below the level and the
    /// depth its mean of the level less the ground where that is above none.
    pub fn under(&self, k: usize, level: f64) -> (f64, f64) {
        let spread = (self.spread[k] as f64).max(LEAST_SPREAD);
        let z = (level - self.mean[k] as f64) / spread;
        let share = normal_cdf(z);
        let density = (-0.5 * z * z).exp() / (2.0 * std::f64::consts::PI).sqrt();
        (share, (spread * (z * share + density)).max(0.0))
    }
}

/// A quantity's mean over the cells about each cell (F416), the cells read as samples of the
/// ground within it: a cover's share of each cell, where what the cover lies over is decided
/// cell by cell, over the same Gaussian the relief at the cell's scale is read over.
pub fn over_cells(grid: &CubeSphere, values: &[f32]) -> Vec<f32> {
    let passes = (RELIEF_SPREAD_CELLS * RELIEF_SPREAD_CELLS / 0.25).round() as usize;
    smooth(grid, values, passes)
}

/// A finish's liquids as its maps draw them (F416): the share of each cell under its sea and
/// the sea's mean depth over it, m; the lakes' mean depth over each cell, m; and the surface
/// the maps are lit over, m, the sea's flat level over the share of each cell it covers and
/// the ground over the rest.
pub struct Drawn {
    pub sea_share: Vec<f32>,
    pub sea: Vec<f32>,
    pub lakes: Vec<f32>,
    pub surface: Vec<f32>,
}

impl Drawn {
    pub fn of(out: &FinishOutput) -> Self {
        let ground = out.height.as_slice();
        let n = ground.len();
        let grid = out.grid;
        let relief = Relief::of(&out.height);
        let level = out.sea_level;
        let sea: Vec<(f64, f64)> = (0..n)
            .into_par_iter()
            .map(|k| {
                if level.is_finite() {
                    relief.under(k, level)
                } else {
                    (0.0, 0.0)
                }
            })
            .collect();
        // A lake's level is its water's over the ground where it stands (F277), carried to
        // the cells about it as far as the relief at the cells' scale is read, each beside a
        // lake taking the highest lake's beside it: its share under that level its own relief's.
        let lakes: Vec<f32> = match &out.rivers {
            Some(r) => {
                let depth = r.lakes.as_slice();
                let mut lake_level: Vec<f32> = (0..n)
                    .map(|k| {
                        if depth[k] > 0.0 {
                            ground[k] + depth[k]
                        } else {
                            f32::NEG_INFINITY
                        }
                    })
                    .collect();
                for _ in 0..RELIEF_SPREAD_CELLS.ceil() as usize {
                    lake_level = (0..n)
                        .into_par_iter()
                        .map(|k| {
                            if depth[k] > 0.0 {
                                return lake_level[k];
                            }
                            let c = grid.cell(k);
                            AROUND.iter().fold(lake_level[k], |a, &(di, dj)| {
                                a.max(lake_level[grid.index(grid.neighbor(c, di, dj))])
                            })
                        })
                        .collect();
                }
                (0..n)
                    .into_par_iter()
                    .map(|k| {
                        if lake_level[k].is_finite() {
                            relief.under(k, lake_level[k] as f64).1 as f32
                        } else {
                            0.0
                        }
                    })
                    .collect()
            }
            None => vec![0.0; n],
        };
        let surface = (0..n)
            .map(|k| {
                let (share, _) = sea[k];
                (share * level.max(ground[k] as f64) + (1.0 - share) * ground[k] as f64) as f32
            })
            .collect();
        Drawn {
            sea_share: sea.iter().map(|s| s.0 as f32).collect(),
            sea: sea.iter().map(|s| s.1 as f32).collect(),
            lakes,
            surface,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// F416: a coast the ground's shape places and its relief at the cells' scale speckles:
    /// over a plain rising through the sea's level, each cell's ground spread at random by
    /// 100 m, the cells drawn all or none speckle the coast with cells cut off from their kind
    /// on either side; drawn by their shares, half the cells' share lies where the plain itself
    /// crosses the level, within a cell, and hardly a cell stands apart from its kind.
    #[test]
    fn a_coast_is_placed_by_the_ground_s_shape_not_its_relief() {
        let grid = CubeSphere::new(128);
        let n = grid.cell_count();
        let f = grid.face_size();
        let rng = crate::rng::Rng::new(5);
        let rise = 1_500.0;
        let ground = Field::from_fn(grid, |c, d| {
            let noise = 100.0 * (rng.uniform(grid.index(c) as u64, 0) - 0.5) * 12f64.sqrt();
            (rise * d[1] + noise) as f32
        });
        let relief = Relief::of(&ground);
        let inner = |k: usize| {
            let c = grid.cell(k);
            c.face == 0 && (20..f - 20).contains(&c.i) && (20..f - 20).contains(&c.j)
        };
        let wet_all = |k: usize| ground.as_slice()[k] < 0.0;
        let wet_share = |k: usize| relief.under(k, 0.0).0 >= 0.5;
        let apart = |wet: &dyn Fn(usize) -> bool| {
            (0..n)
                .filter(|&k| inner(k))
                .filter(|&k| {
                    let c = grid.cell(k);
                    let alike = [
                        (1, 0),
                        (-1, 0),
                        (0, 1),
                        (0, -1),
                        (1, 1),
                        (1, -1),
                        (-1, 1),
                        (-1, -1),
                    ]
                    .iter()
                    .filter(|&&(di, dj)| wet(grid.index(grid.neighbor(c, di, dj))) == wet(k))
                    .count();
                    alike <= 1
                })
                .count()
        };
        let (cells, shares) = (apart(&wet_all), apart(&wet_share));
        // Where half the share lies against where the plain crosses the level, in cells.
        let cell = grid.delta();
        let (mut off, mut count) = (0.0, 0);
        for k in (0..n).filter(|&k| inner(k)) {
            let (share, _) = relief.under(k, 0.0);
            if (share - 0.5).abs() < 0.1 {
                off += (grid.center(grid.cell(k))[1] / cell).abs();
                count += 1;
            }
        }
        let off = off / count.max(1) as f64;
        println!(
            "cells apart from their kind: {cells} drawn all or none, {shares} by shares; half the share {off:.2} cells from where the plain crosses the level ({count} cells)"
        );
        assert!(cells > 50, "{cells}");
        assert!(shares <= cells / 20, "{shares} {cells}");
        assert!(count > 0 && off < 1.0, "{off}");
    }

    /// F416: an even slope has no relief at its cells' scale, its every cell its
    /// neighbourhood's mean; a ground spread at random about a level by an amount spreads by
    /// near that amount; at its mean a ground is half under a liquid there, its mean depth the
    /// spread's 0.399 (the normal's density at its mean); far under, all of it under, its
    /// depth the level less the ground.
    #[test]
    fn a_cell_s_relief_is_what_its_grid_does_not_place() {
        let grid = CubeSphere::new(128);
        let n = grid.cell_count();
        // A plane along each face's own rows, read on cells whose Gaussian keeps to the face.
        let slope = Field::from_fn(grid, |c, _| (10.0 * c.i as f64 + 3.0 * c.j as f64) as f32);
        let even = Relief::of(&slope);
        let f = grid.face_size();
        let worst = (0..n)
            .filter(|&k| {
                let c = grid.cell(k);
                (20..f - 20).contains(&c.i) && (20..f - 20).contains(&c.j)
            })
            .fold(0.0f32, |a, k| a.max(even.spread[k]));
        let rng = crate::rng::Rng::new(3);
        let noisy = Field::from_fn(grid, |c, _| {
            (100.0 * (rng.uniform(grid.index(c) as u64, 0) - 0.5) * 12f64.sqrt()) as f32
        });
        let rough = Relief::of(&noisy);
        let mean_spread = rough.spread.iter().map(|&s| s as f64).sum::<f64>() / n as f64;
        println!(
            "an even slope rising 10 m a cell spreads by {worst:.2e} m at most; a ground spread by 100 m at random spreads by {mean_spread:.1} m"
        );
        assert!(worst < 1e-3, "{worst}");
        assert!((mean_spread / 100.0 - 1.0).abs() < 0.1, "{mean_spread}");
        let k = 0;
        let relief = Relief {
            mean: vec![-50.0; n],
            spread: vec![100.0; n],
        };
        let (half, depth) = relief.under(k, -50.0);
        assert!(
            (half - 0.5).abs() < 1e-6 && (depth - 39.894).abs() < 1e-2,
            "{half} {depth}"
        );
        let (all, deep) = relief.under(k, 1_000.0);
        assert!(
            all > 1.0 - 1e-9 && (deep - 1_050.0).abs() < 1e-6,
            "{all} {deep}"
        );
        let (none, dry) = relief.under(k, -1_000.0);
        assert!(none < 1e-9 && dry < 1e-9, "{none} {dry}");
    }
}
