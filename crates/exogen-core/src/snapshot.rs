// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! The body as it looks at the snapshot's season (ADR 0007): the ground's own colour, under
//! its liquid by the liquid's optics, under the ices that come and go as the cryosphere's maps
//! draw them at the snapshot's phase; and the surface it is lit over, its liquid's own where a
//! liquid lies. A composite for viewing and for an engine that wants it baked; the albedo, the
//! liquid's depth and the cryosphere's maps are the parts an engine composes itself.

use crate::climate::Window;
use crate::fields::Field;
use crate::stages::FinishOutput;
use crate::stages::materials::LiquidOptics;

/// The surface's height, m: the ground's, and where a liquid lies its liquid's surface, flat
/// over a sea or a lake.
pub fn surface(out: &FinishOutput) -> Field<f32> {
    let grid = out.grid;
    let (height, liquid) = (out.height.as_slice(), out.layers.liquid.as_slice());
    Field::from_fn(grid, |c, _| {
        let k = grid.index(c);
        height[k] + liquid[k].max(0.0)
    })
}

/// The colour of the surface at the snapshot's season, linear reflectance: the albedo under
/// the depth of its working liquid, `working` (`LiquidOptics::over`), then the ice that lies
/// on each cell at the snapshot's phase as `Window::at_phase` reads the cryosphere's maps,
/// hiding 1 - exp(-m / cover mass) of what lies under it.
pub fn snapshot(out: &FinishOutput, working: &str) -> Field<[f32; 3]> {
    let optics = LiquidOptics::of(working);
    let grid = out.grid;
    let liquid = out.layers.liquid.as_slice();
    let ice = out.cryosphere_ice.as_slice();
    let [freeze, melt, thickness, at] = &out.cryosphere;
    let phase = out.cryosphere_snapshot;
    Field::from_fn(grid, |c, _| {
        let k = grid.index(c);
        let mut colour = optics.over(out.albedo.get(c), liquid[k]);
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
            let cover = 1.0 - (-mass / ice.cover_mass.max(1e-30)).exp();
            for (v, &tone) in colour.iter_mut().zip(&ice.colour) {
                *v += (tone - *v) * cover as f32;
            }
        }
        colour
    })
}
