// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! The body as it looks at the snapshot's season (ADR 0007): the ground's own colour, under
//! its liquid by the liquid's optics, under the ices that come and go as the cryosphere's maps
//! draw them at the snapshot's phase; and the surface it is lit over, its liquid's own where a
//! liquid lies. A composite for viewing and for an engine that wants it baked; the albedo, the
//! liquid's depth and the cryosphere's maps are the parts an engine composes itself.

use crate::climate::Window;
use crate::fields::Field;
use crate::grid::CubeSphere;
use crate::shares::{Drawn, over_cells, wet};
use crate::stages::FinishOutput;
use crate::stages::materials::LiquidOptics;

/// The surface's height, m: the ground's, and where a sea lies its flat surface, over the
/// share of each cell it covers (F416, F417).
pub fn surface(out: &FinishOutput, drawn: &Drawn) -> Field<f32> {
    let grid = out.grid;
    Field::from_fn(grid, |c, _| drawn.surface[grid.index(c)])
}

/// The colour of the surface at the snapshot's season, linear reflectance: each cell the sea's
/// share of it and the land's (F416), each its own kind's colour under its own ice (F417). The
/// sea's share its water as it shows (`waters`), each wet cell's floor under its working
/// liquid, `working`, at its own depth (`LiquidOptics::over`), under the sea's ice; the land's
/// its ground, its albedo (`own_kind`), under the land's ice. Each ice is the one that lies on
/// its cells at the snapshot's phase as `Window::at_phase` reads the cryosphere's maps, hiding
/// 1 - exp(-m / cover mass) of what lies under it: its cover and its colour each a mean over
/// the cells of its kind about the cell, so that an ice's edge the ground's relief at the
/// cell's scale decides cell by cell is drawn by its share of the cell, and a sea's ice covers
/// the sea's share of a cell and a land's snow the land's.
pub fn snapshot(out: &FinishOutput, working: &str, drawn: &Drawn) -> Field<[f32; 3]> {
    let optics = LiquidOptics::of(working);
    let grid = out.grid;
    let ices = ices(out);
    let wet = wet(out);
    let dry: Vec<f32> = wet.iter().map(|w| 1.0 - w).collect();
    let albedo = out.albedo.as_slice();
    let water = waters(&grid, &optics, albedo, out.layers.liquid.as_slice(), &wet);
    let land = own_kind(&grid, albedo, &dry);
    Field::from_fn(grid, |c, _| {
        let k = grid.index(c);
        by_shares(
            (water[k], land[k]),
            drawn.sea_share[k],
            |under, share| ices.sea.over(k, under, share),
            |under, share| ices.land.over(k, under, share),
        )
    })
}

/// A colour of one kind over its share of each cell (F417), `kind` 1 on the kind's cells and 0
/// on the rest: a cell's own `values` where the cell is of the kind, and where not the mean of
/// the kind's cells' about it, each read as a sample of the ground within it; its own where
/// none of the kind lies about it. The materials the model lays on each cell by its side of
/// the coast (a floor's sediment, a land's rock) are drawn over their own kind's share, not
/// over the other's.
fn own_kind(grid: &CubeSphere, values: &[[f32; 3]], kind: &[f32]) -> Vec<[f32; 3]> {
    let share = over_cells(grid, kind);
    let sums: [Vec<f32>; 3] = std::array::from_fn(|ch| {
        let plane: Vec<f32> = values.iter().zip(kind).map(|(v, w)| v[ch] * w).collect();
        over_cells(grid, &plane)
    });
    (0..values.len())
        .map(|k| {
            if kind[k] > 0.0 || share[k] <= 0.0 {
                values[k]
            } else {
                std::array::from_fn(|ch| (sums[ch][k] / share[k]).clamp(0.0, 1.0))
            }
        })
        .collect()
}

/// The sea's water as it shows over the sea's share of each cell (F417), each cell under it,
/// 1, or not, 0, by `wet`: a wet cell's own floor, its `albedo`, under its own depth of the
/// liquid, `liquid`, m, by the liquid's `optics`; a dry cell's the mean of the wet cells'
/// colours about it, each a sample of the water within it, so that a shore beside its shallows
/// shows them and not the colour of the mean depth of the shallows and the deep beyond.
fn waters(
    grid: &CubeSphere,
    optics: &LiquidOptics,
    albedo: &[[f32; 3]],
    liquid: &[f32],
    wet: &[f32],
) -> Vec<[f32; 3]> {
    let shown: Vec<[f32; 3]> = albedo
        .iter()
        .zip(liquid)
        .map(|(&a, &d)| optics.over(a, d))
        .collect();
    own_kind(grid, &shown, wet)
}

/// A cell's colour, linear reflectance, `share` of it under the sea (F417): the sea's share
/// its water as it shows under the sea's ice, `sea_ice`, and the rest the land's ground under
/// the land's, `land_ice`, `(water, ground)` the colours under them, each ice given the colour
/// under it and its kind's share.
fn by_shares(
    (water, ground): ([f32; 3], [f32; 3]),
    share: f32,
    sea_ice: impl Fn([f32; 3], f32) -> [f32; 3],
    land_ice: impl Fn([f32; 3], f32) -> [f32; 3],
) -> [f32; 3] {
    let land = land_ice(ground, 1.0 - share);
    if share <= 0.0 {
        return land;
    }
    let sea = sea_ice(water, share);
    std::array::from_fn(|ch| share * sea[ch] + (1.0 - share) * land[ch])
}

/// An ice over one kind's share of each cell (F417): the share of each cell it hides and its
/// colour by it in each channel, each summed over the cells of that kind about the cell, the
/// other kind's as none.
pub struct Cover {
    pub cover: Vec<f32>,
    pub toned: [Vec<f32>; 3],
}

impl Cover {
    /// Over its own cells about each cell, each cell's own `cover` and its colour by it,
    /// `toned`, where `kind` holds 1 and none where it holds 0.
    fn of(grid: &CubeSphere, kind: &[f32], cover: &[f32], toned: &[Vec<f32>; 3]) -> Self {
        let only = |v: &[f32]| -> Vec<f32> { v.iter().zip(kind).map(|(a, b)| a * b).collect() };
        Cover {
            cover: over_cells(grid, &only(cover)),
            toned: std::array::from_fn(|ch| over_cells(grid, &only(&toned[ch]))),
        }
    }

    /// `under`, the colour of the kind holding `share` of cell `k`, under the ice.
    fn over(&self, k: usize, under: [f32; 3], share: f32) -> [f32; 3] {
        let cover = self.cover[k];
        if cover <= 0.0 || share <= 0.0 {
            return under;
        }
        let hidden = (cover / share).min(1.0);
        std::array::from_fn(|ch| under[ch] + (self.toned[ch][k] / cover - under[ch]) * hidden)
    }
}

/// The ices that lie at the snapshot's phase as the snapshot draws them (F416, F417).
pub struct Ices {
    /// The share of each cell its own ice hides, 1 - exp(-m / cover mass).
    pub own: Vec<f32>,
    /// That share drawn, its mean over the cells about each cell: the sea's ice and the
    /// land's together.
    pub drawn: Vec<f32>,
    /// The ice over the sea's share of each cell.
    pub sea: Cover,
    /// The ice over the land's share of each cell.
    pub land: Cover,
}

/// The ices that lie at the snapshot's phase as the snapshot draws them (F416, F417): each
/// cell's own, and each kind's over its share of each cell.
pub fn ices(out: &FinishOutput) -> Ices {
    let grid = out.grid;
    let n = grid.cell_count();
    let ice = out.cryosphere_ice.as_slice();
    let [freeze, melt, thickness, at] = &out.cryosphere;
    let phase = out.cryosphere_snapshot;
    // Each cell's own ice: its cover, and its colour by its cover.
    let mut cover = vec![0.0f32; n];
    let mut toned = [vec![0.0f32; n], vec![0.0f32; n], vec![0.0f32; n]];
    for k in 0..n {
        let index = ice[k].round() as usize;
        let lying = (index > 0 && phase.is_finite())
            .then(|| out.cryosphere_ices.get(index - 1))
            .flatten();
        if let Some(ice) = lying {
            let window = Window {
                freeze: freeze.as_slice()[k] as f64,
                melt: melt.as_slice()[k] as f64,
                peak: thickness.as_slice()[k] as f64,
                at: at.as_slice()[k] as f64,
            };
            let mass = window.at_phase(phase, ice.curve.unwrap_or((1.0, 1.0))) * ice.density;
            cover[k] = (1.0 - (-mass / ice.cover_mass.max(1e-30)).exp()) as f32;
            for (t, &tone) in toned.iter_mut().zip(&ice.colour) {
                t[k] = cover[k] * tone;
            }
        }
    }
    let wet = wet(out);
    let dry: Vec<f32> = wet.iter().map(|w| 1.0 - w).collect();
    Ices {
        drawn: over_cells(&grid, &cover),
        sea: Cover::of(&grid, &wet, &cover, &toned),
        land: Cover::of(&grid, &dry, &cover, &toned),
        own: cover,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// F417: a cell a fifth under a deep sea is that fifth the sea's colour and four fifths its
    /// ground bare: drawn as the sea over the whole cell, a pale ground under a sea a fifth wet
    /// shows as dark as the sea.
    #[test]
    fn a_cell_part_under_the_sea_is_its_shares_colours() {
        let optics = LiquidOptics::of("water");
        let bare = |under: [f32; 3], _: f32| under;
        let ground = [0.4f32, 0.35, 0.3];
        let water = optics.over(ground, 500.0);
        let colour = by_shares((water, ground), 0.2, bare, bare);
        let want: [f32; 3] = std::array::from_fn(|ch| 0.2 * water[ch] + 0.8 * ground[ch]);
        println!("a fifth under 500 m: {colour:?}, its shares {want:?}; the sea {water:?}");
        for ch in 0..3 {
            assert!((colour[ch] - want[ch]).abs() < 1e-6, "{colour:?} {want:?}");
        }
        assert!(colour[0] > 0.75 * ground[0], "{colour:?}");
    }

    /// F417: a shore shows its shallows: land beside a shelf 10 m deep two cells wide, a sea
    /// 2,000 m deep beyond, over a pale floor. Each wet cell's water is its own floor under its
    /// own depth; a dry cell's within a cell of the shore the mean of its wet cells' colours,
    /// the shelf's green showing in it, where the colour of the shallows' and the deep's mean
    /// depth, some hundreds of metres, showed the deep sea's.
    #[test]
    fn a_shore_shows_its_shallows() {
        let grid = CubeSphere::new(128);
        let n = grid.cell_count();
        let f = grid.face_size();
        let cell = grid.delta();
        let optics = LiquidOptics::of("water");
        let y = |k: usize| grid.center(grid.cell(k))[1] / cell;
        let inner = |k: usize| {
            let c = grid.cell(k);
            c.face == 0 && (20..f - 20).contains(&c.i) && (20..f - 20).contains(&c.j)
        };
        let liquid: Vec<f32> = (0..n)
            .map(|k| match y(k) {
                y if y >= 0.0 => 0.0,
                y if y > -2.0 => 10.0,
                _ => 2_000.0,
            })
            .collect();
        let wet: Vec<f32> = liquid
            .iter()
            .map(|&d| if d > 0.0 { 1.0 } else { 0.0 })
            .collect();
        let albedo = vec![[0.3f32; 3]; n];
        let water = waters(&grid, &optics, &albedo, &liquid, &wet);
        let (share, depth) = (over_cells(&grid, &wet), over_cells(&grid, &liquid));
        let (mut own, mut shore, mut mean, mut cells) = (0.0f32, 1.0f32, 0.0f32, 0);
        for k in (0..n).filter(|&k| inner(k)) {
            if wet[k] > 0.0 {
                own = own.max((water[k][1] - optics.over(albedo[k], liquid[k])[1]).abs());
            } else if y(k) < 1.0 {
                shore = shore.min(water[k][1]);
                mean = mean.max(optics.over(albedo[k], depth[k] / share[k])[1]);
                cells += 1;
            }
        }
        let shelf = optics.over([0.3; 3], 10.0)[1];
        println!(
            "wet cells off their own colour by {own:.1e}; the shore's water green {shore:.3} at least ({cells} cells), the shelf's {shelf:.3}, the mean depth's {mean:.3} at most"
        );
        assert!(own < 1e-6, "{own}");
        assert!(cells > 0 && shore > mean + 0.02, "{shore} {mean}");
    }

    /// F417: flats where the cells fall wet and dry at random, the model's floor pale under its
    /// sediment, 0.6, and its land dark, 0.1. Each kind's ground drawn from its own cells, the
    /// land's share shows the land's 0.1 and the floor the floor's 0.6 on every cell, and a
    /// cell's colour follows its share alone; each cell's own ground drawn under both, a wet
    /// cell's land shows the floor's pale sediment and the flats speckle cell by cell.
    #[test]
    fn each_kind_shows_its_own_ground() {
        let grid = CubeSphere::new(64);
        let n = grid.cell_count();
        let rng = crate::rng::Rng::new(9);
        let wet: Vec<f32> = (0..n)
            .map(|k| {
                if rng.uniform(k as u64, 0) < 0.5 {
                    1.0
                } else {
                    0.0
                }
            })
            .collect();
        let albedo: Vec<[f32; 3]> = wet
            .iter()
            .map(|&w| if w > 0.0 { [0.6; 3] } else { [0.1; 3] })
            .collect();
        let dry: Vec<f32> = wet.iter().map(|w| 1.0 - w).collect();
        let (floor, land) = (
            own_kind(&grid, &albedo, &wet),
            own_kind(&grid, &albedo, &dry),
        );
        let (mut worst, mut own) = (0.0f32, 0.0f32);
        for ((l, f), a) in land.iter().zip(&floor).zip(&albedo) {
            worst = worst.max((l[0] - 0.1).abs()).max((f[0] - 0.6).abs());
            // Each cell's own ground under both kinds instead.
            own = own.max((a[0] - 0.1).abs());
        }
        println!(
            "the kinds' grounds off theirs by {worst:.2e} at most; each cell's own under the land off it by {own:.2}"
        );
        assert!(worst < 1e-5, "{worst}");
        assert!(own > 0.4, "{own}");
    }

    /// F417: a sea frozen whole beside bare land, and a land under snow beside an open sea.
    /// Each ice drawn over its own kind's share of each cell, no cell shows the frozen sea's
    /// water nor the snowy land's ground; drawn over the whole of each cell, the cells along
    /// the coast show both, the ice's mean over them thinned by the other kind's cells.
    #[test]
    fn each_ice_covers_its_own_kind_s_share() {
        let grid = CubeSphere::new(64);
        let n = grid.cell_count();
        let wet: Vec<f32> = (0..n)
            .map(|k| {
                if grid.center(grid.cell(k))[1] < 0.0 {
                    1.0
                } else {
                    0.0
                }
            })
            .collect();
        let dry: Vec<f32> = wet.iter().map(|w| 1.0 - w).collect();
        let share = over_cells(&grid, &wet);
        let white = |own: &[f32]| -> [Vec<f32>; 3] { std::array::from_fn(|_| own.to_vec()) };
        for (frozen, sea, name) in [(&wet, true, "a frozen sea"), (&dry, false, "a snowy land")] {
            let ice = Cover::of(&grid, frozen, frozen, &white(frozen));
            let whole = over_cells(&grid, frozen);
            let (mut shown, mut thinned) = (0.0f32, 0.0f32);
            for k in 0..n {
                let s = if sea { share[k] } else { 1.0 - share[k] };
                if s <= 0.0 {
                    continue;
                }
                // The frozen kind's own colour, none, under its ice: the share of the cell it
                // shows over.
                let drawn = ice.over(k, [0.0; 3], s)[0];
                shown = shown.max(s * (1.0 - drawn));
                // Its ice drawn over the whole cell instead, its mean over all the cells.
                thinned = thinned.max(s * (1.0 - whole[k]));
            }
            println!(
                "{name}: the most of a cell its own ground shows over {shown:.2e}; its ice drawn over the whole cell {thinned:.2}"
            );
            assert!(shown < 1e-5, "{name} {shown}");
            assert!(thinned > 0.2, "{name} {thinned}");
        }
    }
}
