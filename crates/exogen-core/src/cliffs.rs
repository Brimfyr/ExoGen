// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! The cliffs the waves cut and the platforms at their feet (F430).
//!
//! The waves cut a cliff's toe while their force there, rho g H_f, over its rock's uniaxial
//! compressive strength exceeds a threshold (Sunamura 2015, eq. 7). The cliff retreats and
//! leaves a platform, across which the waves fall as they cross it from where they break at
//! its edge (Matsumoto et al. 2016; Shadrick et al. 2021): the wider the platform, the less
//! reaches the cliff, and under a sea at one level its retreat slows until it stops (Trenhaile
//! 2000; Walkden and Hall 2005; Shadrick et al. 2021). Under the waves exceeded 12 hours a year
//! (F428), a platform widens to where they no longer cut: within some centuries for ground
//! the waves cut at tens of centimetres a year (Premaillon et al. 2018), well within the span
//! the finish holds its sea at its level (F278). A beach the rivers laid before the cliff
//! (F414, F429) is width the waves cross too; Everts found the cliffs behind a beach 60 m
//! wide no longer cut (Sunamura 2015). The ground rises from the shore at its cell's slope, the
//! rivers' loose deposits over the body's rock: the toe stands in the cover until the rock's
//! top rises above the level. Platforms of some hundreds of metres lie far below the cells
//! the finish resolves: the coast's kind, its platform, its cliff and its beach are context
//! for the detail pass (PLAN section 4.11; the user, 2026-11-24), and what the waves cut is
//! laid at the cliffs' feet.

use crate::constants::*;
use crate::fields::LayerSet;
use crate::grid::CubeSphere;
use crate::waves::Waves;

/// No coast here, or one its waves do not reach: the land's inland cells, the coasts under
/// ice, and those without waves (F430).
pub const COAST_NONE: u8 = 0;
/// A coast whose beach holds its waves off its toe: the rivers' load laid before it.
pub const COAST_BEACH: u8 = 1;
/// A coast whose waves cut nothing: its ground stands as the terrain left it.
pub const COAST_STANDING: u8 = 2;
/// A cliff cut in the loose ground the rivers laid, its platform at its foot.
pub const COAST_LOOSE_CLIFF: u8 = 3;
/// A cliff cut in the body's rock, its platform at its foot.
pub const COAST_ROCK_CLIFF: u8 = 4;

/// The coasts the waves work (F430): whether each cell of the finish's grid is on the land
/// beside the liquid along an edge; at each, the width of the platform its waves cut, m; the
/// height of the cliff behind it, m; the width of the beach before it, m; and its kind, one of
/// COAST_*; none elsewhere. What they cut, m³, all laid at the cliffs' feet.
#[derive(Clone, Debug)]
pub struct Coast {
    pub coastal: Vec<bool>,
    pub platform: Vec<f32>,
    pub cliff: Vec<f32>,
    pub beach: Vec<f32>,
    pub kind: Vec<u8>,
    pub cut: f64,
}

/// The coasts in brief (F430): the cells on the land beside the liquid; the share of them of
/// each kind, COAST_* (the first those under ice or without waves); the median and the most
/// of the platforms and the cliffs where the waves cut, m; and what they cut, m³.
#[derive(Clone, Copy, Debug, Default)]
pub struct Summary {
    pub cells: usize,
    pub shares: [f64; 5],
    pub platform_median: f64,
    pub platform_most: f64,
    pub cliff_median: f64,
    pub cliff_most: f64,
    pub cut: f64,
}

impl Coast {
    pub fn summary(&self) -> Summary {
        let coast: Vec<usize> = (0..self.kind.len()).filter(|&k| self.coastal[k]).collect();
        let mut shares = [0.0f64; 5];
        for &k in &coast {
            shares[self.kind[k] as usize] += 1.0;
        }
        if !coast.is_empty() {
            shares.iter_mut().for_each(|s| *s /= coast.len() as f64);
        }
        let cut: Vec<usize> = coast
            .iter()
            .copied()
            .filter(|&k| self.platform[k] > 0.0)
            .collect();
        let median_most = |values: &[f32]| -> (f64, f64) {
            let mut v: Vec<f64> = cut.iter().map(|&k| values[k] as f64).collect();
            v.sort_by(f64::total_cmp);
            (
                v.get(v.len() / 2).copied().unwrap_or(0.0),
                v.last().copied().unwrap_or(0.0),
            )
        };
        let (platform_median, platform_most) = median_most(&self.platform);
        let (cliff_median, cliff_most) = median_most(&self.cliff);
        Summary {
            cells: coast.len(),
            shares,
            platform_median,
            platform_most,
            cliff_median,
            cliff_most,
            cut: self.cut,
        }
    }
}

/// The compressive strength of the body's rock, Pa, by its bedrock's kind `bedrock` (F430):
/// ice, salt, or the silicate crusts'.
pub fn rock_strength(bedrock: &str) -> f64 {
    if crate::stages::materials::is_ice(bedrock) {
        CLIFF_ICE_STRENGTH
    } else if bedrock.ends_with("_salt") {
        CLIFF_SALT_STRENGTH
    } else {
        CLIFF_ROCK_STRENGTH
    }
}

/// The width the platform at a cliff's foot widens to (F430), m: where waves `height` m high
/// where they break at its edge, falling as exp(-CLIFF_WAVE_DECAY x) across it, no longer cut
/// its toe of compressive strength `strength`, Pa, in a liquid of `density` kg/m³ under gravity
/// `g`, rho g H exp(-y x) / S_c = CLIFF_THRESHOLD (Sunamura 2015, eq. 7); none where they
/// never could.
pub fn still_width(height: f64, strength: f64, (density, g): (f64, f64)) -> f64 {
    let force = density * g * height / (CLIFF_THRESHOLD * strength);
    if force > 1.0 {
        force.ln() / CLIFF_WAVE_DECAY
    } else {
        0.0
    }
}

/// One coast's cut (F430): under waves exceeded 12 hours a year `height` m high, working the
/// shore down to their depth of closure `closure` m; a deposit `laid` m deep over the liquid's
/// cell before it, `cell` m across; the land rising `rise` m over its level at its cell's
/// centre, half a cell from its shore, its loose cover `cover` m deep over a rock of
/// compressive strength `rock` Pa; in a liquid of density and gravity `liquid`. The width of
/// the platform its waves cut, m, the height of the cliff behind it, m, the width of the beach
/// before it, m, and its kind.
pub fn shore(
    (height, closure): (f64, f64),
    (laid, cell): (f64, f64),
    (rise, cover): (f64, f64),
    rock: f64,
    liquid: (f64, f64),
) -> (f64, f64, f64, u8) {
    if height <= 0.0 {
        return (0.0, 0.0, 0.0, COAST_NONE);
    }
    // Waves that cut not even loose ground work nothing on the coast: it stands, and what its
    // rivers laid before it is their own, no beach the waves built.
    let loose = still_width(height, CLIFF_LOOSE_STRENGTH, liquid);
    if loose <= 0.0 {
        return (0.0, 0.0, 0.0, COAST_STANDING);
    }
    // The beach: the deposit over a cell offshore, spread over the shoreface its waves work,
    // as deep as their closure (F429), no wider than the cell it lies on: wider, it is the
    // coast's own advance, which the cells draw.
    let beach = if laid > 0.0 {
        (laid * cell / closure.max(height)).min(cell)
    } else {
        0.0
    };
    let loose = (loose - beach).max(0.0);
    let hard = (still_width(height, rock, liquid) - beach).max(0.0);
    // The ground rises from its shore at its cell's slope; the toe, at the level, stands in its
    // cover until the rock's top rises above it, the cover's depth over the slope inland.
    let slope = rise.max(0.0) / (0.5 * cell);
    let reach = if cover <= 0.0 {
        0.0
    } else if slope > 0.0 {
        cover / slope
    } else {
        f64::INFINITY
    };
    let (platform, kind) = if loose <= reach {
        (loose, COAST_LOOSE_CLIFF)
    } else if hard > reach {
        (hard, COAST_ROCK_CLIFF)
    } else {
        (reach, COAST_LOOSE_CLIFF)
    };
    if platform <= 0.0 {
        let kind = if beach > 0.0 {
            COAST_BEACH
        } else {
            COAST_STANDING
        };
        return (0.0, 0.0, beach, kind);
    }
    (platform, (slope * platform).min(rise.max(0.0)), beach, kind)
}

/// The coasts the waves work on the finish's `grid` (F430), its liquid where `layers` hold it
/// at `sea_level`, m, on a body of `radius` m: each cell of the land beside the liquid along an
/// edge cut by `waves` as `shore` gives, the deposits `laid` over the span at each cell, m,
/// before it, its rock of compressive strength `rock` Pa, in a liquid of density and gravity
/// `liquid`. What each cliff's retreat cut, the wedge of ground between its shore and its
/// line, is taken from its cell, its cover first, and laid at its foot on the liquid's cells
/// beside it, the liquid giving way to it.
pub fn cut(
    grid: CubeSphere,
    layers: &mut LayerSet,
    laid: &[f32],
    waves: &Waves,
    (sea_level, radius): (f64, f64),
    rock: f64,
    liquid: (f64, f64),
) -> Coast {
    let n = grid.cell_count();
    let sea: Vec<bool> = layers.liquid.as_slice().iter().map(|&l| l > 0.0).collect();
    let height = layers.height();
    let mut coast = Coast {
        coastal: vec![false; n],
        platform: vec![0.0; n],
        cliff: vec![0.0; n],
        beach: vec![0.0; n],
        kind: vec![COAST_NONE; n],
        cut: 0.0,
    };
    let mut moves: Vec<(usize, f64, Vec<usize>)> = Vec::new();
    for k in 0..n {
        if sea[k] {
            continue;
        }
        let c = grid.cell(k);
        let beside: Vec<usize> = [(1, 0), (-1, 0), (0, 1), (0, -1)]
            .iter()
            .map(|&(di, dj)| grid.index(grid.neighbor(c, di, dj)))
            .filter(|&m| sea[m])
            .collect();
        if beside.is_empty() {
            continue;
        }
        coast.coastal[k] = true;
        let w = if layers.ice.as_slice()[k] > 0.0 {
            None
        } else {
            waves.cell_at(grid.center(c))
        };
        let Some(w) = w else {
            continue;
        };
        let cell = grid.cell_area(c).sqrt() * radius;
        let front =
            beside.iter().map(|&m| laid[m].max(0.0) as f64).sum::<f64>() / beside.len() as f64;
        let rise = height.as_slice()[k] as f64 - sea_level;
        let (platform, cliff, beach, kind) = shore(
            (waves.exceeded[w] as f64, waves.closure(w, liquid.1)),
            (front, cell),
            (rise, layers.regolith.as_slice()[k].max(0.0) as f64),
            rock,
            liquid,
        );
        coast.platform[k] = platform as f32;
        coast.cliff[k] = cliff as f32;
        coast.beach[k] = beach as f32;
        coast.kind[k] = kind;
        // The wedge between the shore and the cliff's line, along the cell's shore.
        let volume = 0.5 * platform * cliff * cell * beside.len() as f64;
        if volume > 0.0 {
            moves.push((k, volume, beside));
        }
    }
    let area = |k: usize| grid.cell_area(grid.cell(k)) * radius * radius;
    for (k, volume, beside) in moves {
        coast.cut += volume;
        let mut take = volume / area(k);
        let cover = layers.regolith.as_slice()[k] as f64;
        let from_cover = cover.min(take);
        layers.regolith.as_mut_slice()[k] = (cover - from_cover) as f32;
        take -= from_cover;
        layers.bedrock.as_mut_slice()[k] -= take as f32;
        for &m in &beside {
            let laid = volume / beside.len() as f64 / area(m);
            layers.regolith.as_mut_slice()[m] += laid as f32;
            let l = layers.liquid.as_slice()[m] as f64;
            layers.liquid.as_mut_slice()[m] = (l - laid).max(0.0) as f32;
        }
    }
    coast
}

#[cfg(test)]
mod tests {
    use super::*;

    const WATER: (f64, f64) = (1025.0, 9.81);

    /// F430: a cliff's platform widens until the waves crossing it no longer cut its toe, the
    /// wider the weaker its ground and the higher its waves; the Earth's storm waves cut
    /// nothing from its silicate rock, and Titan's calm methane nothing from its loose ground.
    #[test]
    fn a_cliff_is_cut_until_its_platform_stills_the_waves() {
        let width = still_width(4.0, CLIFF_LOOSE_STRENGTH, WATER);
        assert!(width > 50.0 && width < 500.0, "platform {width:.0} m");
        // At its edge the waves at its toe stand at the threshold.
        let at_toe = 4.0 * (-CLIFF_WAVE_DECAY * width).exp();
        let force = WATER.0 * WATER.1 * at_toe / CLIFF_LOOSE_STRENGTH;
        assert!(
            (force / CLIFF_THRESHOLD - 1.0).abs() < 1e-9,
            "force {force:.6} at the toe"
        );
        assert!(still_width(2.0, CLIFF_LOOSE_STRENGTH, WATER) < width);
        assert!(still_width(4.0, 1.0e6, WATER) < width);
        assert!(still_width(4.0, CLIFF_SALT_STRENGTH, WATER) == 0.0);
        assert!(still_width(20.0, CLIFF_ROCK_STRENGTH, WATER) == 0.0);
        assert!(still_width(0.07, CLIFF_LOOSE_STRENGTH, (450.0, 1.352)) == 0.0);
        assert!(still_width(40.0, CLIFF_SALT_STRENGTH, WATER) > 0.0);
    }

    /// F430: a beach before a cliff is width its waves cross: the platform they cut narrows by
    /// it, and a beach as wide as the platform would be holds them off the cliff altogether, a
    /// beach no wider than the cell it lies on; a sea whose waves cut not even loose ground
    /// builds no beach, its coast standing.
    #[test]
    fn a_beach_holds_the_waves_off_its_cliff() {
        let (cell, closure) = (10_000.0, 8.0);
        let width = still_width(4.0, CLIFF_LOOSE_STRENGTH, WATER);
        let bare = shore(
            (4.0, closure),
            (0.0, cell),
            (10.0, 100.0),
            CLIFF_ROCK_STRENGTH,
            WATER,
        );
        assert_eq!(bare.3, COAST_LOOSE_CLIFF);
        assert!((bare.0 - width).abs() < 1e-9);
        // A deposit making a beach half as wide as the platform.
        let half = 0.5 * width * closure / cell;
        let fronted = shore(
            (4.0, closure),
            (half, cell),
            (10.0, 100.0),
            CLIFF_ROCK_STRENGTH,
            WATER,
        );
        assert!(
            (fronted.2 - 0.5 * width).abs() < 1e-6,
            "beach {}",
            fronted.2
        );
        assert!(
            (fronted.0 - 0.5 * width).abs() < 1e-6,
            "platform {}",
            fronted.0
        );
        assert_eq!(fronted.3, COAST_LOOSE_CLIFF);
        let held = shore(
            (4.0, closure),
            (2.0 * half, cell),
            (10.0, 100.0),
            CLIFF_ROCK_STRENGTH,
            WATER,
        );
        assert_eq!(held.0, 0.0);
        assert_eq!(held.1, 0.0);
        assert_eq!(held.3, COAST_BEACH);
        // A deposit deep over its cell makes a beach no wider than the cell.
        let deep = shore(
            (4.0, closure),
            (50.0, cell),
            (10.0, 100.0),
            CLIFF_ROCK_STRENGTH,
            WATER,
        );
        assert_eq!((deep.2, deep.3), (cell, COAST_BEACH));
        // Titan's calm methane cuts not even loose ground: its coast stands behind what its
        // rivers laid, their own, no beach its waves built.
        let calm = shore(
            (0.07, 0.16),
            (1.0, cell),
            (10.0, 100.0),
            CLIFF_ICE_STRENGTH,
            (450.0, 1.352),
        );
        assert_eq!(calm, (0.0, 0.0, 0.0, COAST_STANDING));
    }

    /// F430: the waves cut a loose cover down to the rock beneath it, and stop where the rock's
    /// top rises above their level, the cliff then as high as the cover is deep; a cover deep
    /// enough is cut to its platform's own width; bare silicate rock stands, a salt rock under
    /// great waves is cut; a cliff stands no higher than its ground.
    #[test]
    fn the_waves_cut_the_cover_down_to_its_rock() {
        let cell = 10_000.0;
        let width = still_width(4.0, CLIFF_LOOSE_STRENGTH, WATER);
        // Ground 100 m up half a cell inland: a slope of 0.02, the cover's 1 m reaching 50 m.
        let thin = shore(
            (4.0, 8.0),
            (0.0, cell),
            (100.0, 1.0),
            CLIFF_ROCK_STRENGTH,
            WATER,
        );
        assert_eq!(thin.3, COAST_LOOSE_CLIFF);
        assert!((thin.0 - 50.0).abs() < 1e-9, "platform {}", thin.0);
        assert!((thin.1 - 1.0).abs() < 1e-9, "cliff {}", thin.1);
        let deep = shore(
            (4.0, 8.0),
            (0.0, cell),
            (100.0, 50.0),
            CLIFF_ROCK_STRENGTH,
            WATER,
        );
        assert!((deep.0 - width).abs() < 1e-9);
        assert!((deep.1 - 0.02 * width).abs() < 1e-9);
        let bare = shore(
            (4.0, 8.0),
            (0.0, cell),
            (100.0, 0.0),
            CLIFF_ROCK_STRENGTH,
            WATER,
        );
        assert_eq!((bare.0, bare.3), (0.0, COAST_STANDING));
        let salt = shore(
            (40.0, 60.0),
            (0.0, cell),
            (100.0, 0.0),
            CLIFF_SALT_STRENGTH,
            WATER,
        );
        assert_eq!(salt.3, COAST_ROCK_CLIFF);
        assert!(salt.0 > 0.0);
        let low = shore(
            (4.0, 8.0),
            (0.0, 400.0),
            (1.0, 50.0),
            CLIFF_ROCK_STRENGTH,
            WATER,
        );
        assert!(low.1 <= 1.0, "cliff {}", low.1);
        assert_eq!(
            shore((0.0, 0.0), (0.0, cell), (100.0, 50.0), 1.0, WATER).3,
            COAST_NONE
        );
    }

    /// F430: on a straight coast of loose ground the waves cut every cell beside the sea, none
    /// inland, a cliff on each as its slope gives, and lay all they cut at the cliffs' feet,
    /// the sea's surface where it stood.
    #[test]
    fn what_the_waves_cut_lies_at_the_cliffs_feet() {
        let (radius, level) = (6.371e6, 0.0);
        let fine = CubeSphere::new(256);
        let n = fine.cell_count();
        let mut layers = LayerSet::new(fine);
        for k in 0..n {
            let d = fine.center(fine.cell(k));
            if d[0] < 0.0 {
                layers.bedrock.as_mut_slice()[k] = -100.0;
                layers.liquid.as_mut_slice()[k] = 100.0;
            } else {
                layers.bedrock.as_mut_slice()[k] = 20.0;
                layers.regolith.as_mut_slice()[k] = 40.0;
            }
        }
        let grid = CubeSphere::new(16);
        let m = grid.cell_count();
        let waves = Waves {
            grid,
            sea: (0..m).map(|k| grid.center(grid.cell(k))[0] < 0.0).collect(),
            mean: vec![1.5; m],
            exceeded: vec![4.0; m],
            exceeded_period: vec![10.0; m],
            net_height: vec![2.0; m],
            net_period: vec![8.0; m],
            energy: vec![[0.0; WAVE_DIRECTIONS]; m],
            fetch: vec![[1.0e6; WAVE_DIRECTIONS]; m],
        };
        let area = |k: usize| fine.cell_area(fine.cell(k)) * radius * radius;
        let ground = |l: &LayerSet| -> Vec<f64> {
            (0..n)
                .map(|k| l.bedrock.as_slice()[k] as f64 + l.regolith.as_slice()[k] as f64)
                .collect()
        };
        let before = ground(&layers);
        let surface: Vec<f64> = layers
            .height()
            .as_slice()
            .iter()
            .zip(layers.liquid.as_slice())
            .map(|(&h, &l)| h as f64 + l as f64)
            .collect();
        let laid = vec![0.0f32; n];
        let coast = cut(
            fine,
            &mut layers,
            &laid,
            &waves,
            (level, radius),
            CLIFF_ROCK_STRENGTH,
            WATER,
        );
        let width = still_width(4.0, CLIFF_LOOSE_STRENGTH, WATER);
        let (mut shore_cells, mut inland) = (0, 0);
        for k in 0..n {
            let c = fine.cell(k);
            if fine.center(c)[0] < 0.0 {
                continue;
            }
            let beside = [(1, 0), (-1, 0), (0, 1), (0, -1)]
                .iter()
                .any(|&(di, dj)| fine.center(fine.neighbor(c, di, dj))[0] < 0.0);
            if beside {
                shore_cells += 1;
                assert_eq!(coast.kind[k], COAST_LOOSE_CLIFF);
                assert!((coast.platform[k] as f64 - width).abs() < 1e-3);
                let cell = fine.cell_area(c).sqrt() * radius;
                let rise = 20.0 + 40.0;
                let expected = rise / (0.5 * cell) * width;
                assert!(
                    (coast.cliff[k] as f64 - expected).abs() < 1e-3 * expected,
                    "cliff {} where {expected}",
                    coast.cliff[k]
                );
            } else {
                inland += 1;
                assert_eq!(coast.kind[k], COAST_NONE);
                assert_eq!(coast.platform[k], 0.0);
            }
        }
        assert!(shore_cells > 0 && inland > 0);
        assert!(coast.cut > 0.0);
        // What the land lost is what the sea's floor took, the sea's surface where it stood.
        let after = ground(&layers);
        let (mut lost, mut taken) = (0.0, 0.0);
        let height = layers.height();
        for k in 0..n {
            let change = (after[k] - before[k]) * area(k);
            if fine.center(fine.cell(k))[0] < 0.0 {
                taken += change;
                let now = height.as_slice()[k] as f64 + layers.liquid.as_slice()[k] as f64;
                assert!((now - surface[k]).abs() < 1e-3, "the sea's surface moved");
            } else {
                lost -= change;
            }
        }
        for (what, v) in [("lost", lost), ("taken", taken)] {
            assert!(
                (v / coast.cut - 1.0).abs() < 1e-2,
                "{what} {v:.4e} m³ of {:.4e} cut",
                coast.cut
            );
        }
        let s = coast.summary();
        assert_eq!(s.cells, shore_cells);
        assert!((s.shares[COAST_LOOSE_CLIFF as usize] - 1.0).abs() < 1e-12);
    }
}
