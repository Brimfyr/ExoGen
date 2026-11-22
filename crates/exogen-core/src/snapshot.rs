// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! The body as it looks at the snapshot's season (ADR 0007): the ground's own colour, under
//! its liquid by the liquid's optics, under the ices that come and go as the cryosphere's maps
//! draw them at the snapshot's phase; and the surface it is lit over, its liquid's own where a
//! liquid lies. A composite for viewing and for an engine that wants it baked; the albedo, the
//! liquid's depth and the cryosphere's maps are the parts an engine composes itself.

use crate::climate::Window;
use crate::fields::Field;
use crate::shares::{Drawn, over_cells};
use crate::stages::FinishOutput;
use crate::stages::materials::LiquidOptics;

/// The surface's height, m: the ground's, and where a sea lies its flat surface, over the
/// share of each cell it covers (F416).
pub fn surface(out: &FinishOutput, drawn: &Drawn) -> Field<f32> {
    let grid = out.grid;
    Field::from_fn(grid, |c, _| drawn.surface[grid.index(c)])
}

/// The colour of the surface at the snapshot's season, linear reflectance: the albedo under
/// its working liquid, `working`, at the sea's mean depth over each cell as the maps draw it
/// (`LiquidOptics::over`, F416), then the ice that lies on each cell at the snapshot's phase as
/// `Window::at_phase` reads the cryosphere's maps, hiding 1 - exp(-m / cover mass) of what lies
/// under it: its cover and its colour each a mean over the cells about the cell (F416), so
/// that an ice's edge the ground's relief at the cell's scale decides cell by cell is drawn by
/// its share of the cell.
pub fn snapshot(out: &FinishOutput, working: &str, drawn: &Drawn) -> Field<[f32; 3]> {
    let optics = LiquidOptics::of(working);
    let grid = out.grid;
    let liquid = &drawn.sea;
    let (_, cover, toned) = ices(out);
    Field::from_fn(grid, |c, _| {
        let k = grid.index(c);
        let mut colour = optics.over(out.albedo.get(c), liquid[k]);
        if cover[k] > 0.0 {
            for (ch, v) in colour.iter_mut().enumerate() {
                let tone = toned[ch][k] / cover[k];
                *v += (tone - *v) * cover[k];
            }
        }
        colour
    })
}

/// The ices that lie at the snapshot's phase as the snapshot draws them (F416): the share of
/// each cell its own ice hides, 1 - exp(-m / cover mass); that share drawn, and its colour by
/// it in each channel, each a mean over the cells about the cell.
pub fn ices(out: &FinishOutput) -> (Vec<f32>, Vec<f32>, [Vec<f32>; 3]) {
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
    let drawn = over_cells(&grid, &cover);
    (cover, drawn, toned.map(|t| over_cells(&grid, &t)))
}
