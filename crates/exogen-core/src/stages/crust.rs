// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! The crust the mantle's flow carries (F184).
//!
//! Where a rocky lid's crust reaches below the lid's stiff base, into the mantle convecting
//! beneath (F173), that part of it is as hot as the mantle and moves with its flow: toward
//! the downwellings, where it piles up and the crust thickens, and away from the upwellings,
//! where it thins (the crustal plateaus over the mantle's downwellings of Bindschadler,
//! Schubert and Kaula 1992). The flow along the lid's base runs from an upwelling to the
//! next downwelling, half a cell, in the time the pattern takes to renew itself, so that
//! its convergence is that speed times the pattern's wavenumber times the pattern's depth
//! there, and the crust thickens by what it carries times that, as far as the depth at which
//! its base turns to eclogite and founders (CRUST_MAX_THICKNESS). What thickens stands
//! compensated on its root, (ρ_m − ρ_c) / ρ_m of the thickening above the ground. The lid
//! holds the crust within it rigid, so a stagnant lid thicker than its crust carries none
//! of it; the lower crust's own flow (F170) spreads what piles up on the time its
//! temperature sets, and what is left when the lid has grown past the crust's base stands.

use rayon::prelude::*;

use crate::constants::*;
use crate::convection::{Support, cells};
use crate::grid::CubeSphere;
use crate::mantle::MantleHistory;

/// What the flow moved over a step: how many cell-substeps carried crust, the most a cell's
/// crust changed by, m, the thickest and the mean crust after, m, and the crust that
/// foundered past the cap, m of the mean.
#[derive(Clone, Copy, Debug, Default)]
pub struct Moved {
    pub carried: usize,
    pub most: f64,
    pub thickest: f64,
    pub mean: f64,
    pub foundered: f64,
}

/// Carries the crust beneath the stiff lid with the mantle's flow from `t0` to `t1`, Gyr,
/// thickening it over the downwellings and thinning it over the upwellings, and raises or
/// lowers the ground by the share its root compensates; what founders past the cap lowers
/// the body's mean radius, the heights' datum, with it (F177), so the ground keeps its
/// height above that.
pub fn drag(
    grid: CubeSphere,
    crust: &mut [f32],
    bedrock: &mut [f32],
    support: &Support,
    mantle: &MantleHistory,
    t0: f64,
    t1: f64,
) -> Moved {
    let airy = ((DENSITY_ROCK - DENSITY_CRUST_ROCK) / DENSITY_ROCK) as f32;
    let n = ((t1 - t0) / CRUST_DRAG_SUBSTEP_GYR).ceil().max(1.0) as usize;
    let dt = (t1 - t0) / n as f64;
    let mut moved = Moved::default();
    for i in 0..n {
        let t = t0 + (i as f64 + 0.5) * dt;
        let stiff = mantle.stiff_lid_at(t) as f32;
        let thickest = crust.iter().fold(0.0f32, |m, &c| m.max(c));
        if thickest <= stiff {
            continue;
        }
        let c = cells(mantle, t);
        if c.renewal <= 0.0 || c.support <= 0.0 {
            continue;
        }
        // The flow's speed along the lid's base, m/s, the pattern's wavenumber, 1/m, and
        // the strain their convergence puts on the carried crust over the substep per unit
        // of the pattern.
        let speed = 0.5 * CONVECTION_CELL_WAVELENGTH * c.depth / (c.renewal * SECONDS_PER_GYR);
        let wavenumber = (c.degree * (c.degree + 1.0)).sqrt() / c.middle;
        let strain = speed * wavenumber * dt * SECONDS_PER_GYR;
        // What the flow takes from under the upwellings, no deeper than the lid's base, is
        // what it brings to the downwellings: the thickening is scaled to the thinning that
        // supplies it, and what piles past the depth at which the crust's base turns to
        // eclogite founders.
        let changes: Vec<f32> = (0..grid.cell_count())
            .into_par_iter()
            .map(|k| {
                let carried = (crust[k] - stiff) as f64;
                if carried <= 0.0 {
                    return 0.0;
                }
                let p = support.pattern(t, grid.center(grid.cell(k)));
                (-carried * strain * p).max(-carried) as f32
            })
            .collect();
        let (mut taken, mut brought, mut whole) = (0.0f64, 0.0f64, 0.0f64);
        for (k, &dc) in changes.iter().enumerate() {
            let a = grid.cell_area(grid.cell(k));
            whole += a;
            if dc < 0.0 {
                taken -= dc as f64 * a;
            } else {
                brought += dc as f64 * a;
            }
        }
        // The flow moves the lesser of what the upwellings give up and what the downwellings
        // take, each side scaled to it, so that no crust is lost but what founders.
        let moving = taken.min(brought);
        let give = if taken > 0.0 { moving / taken } else { 0.0 };
        let take = if brought > 0.0 { moving / brought } else { 0.0 };
        let mut lost = 0.0f64;
        for (k, &dc) in changes.iter().enumerate() {
            if dc == 0.0 {
                continue;
            }
            let now = crust[k] as f64;
            let next = if dc > 0.0 {
                let piled = now + dc as f64 * take;
                let kept = piled.min(CRUST_MAX_THICKNESS.max(now));
                lost += (piled - kept) * grid.cell_area(grid.cell(k));
                kept
            } else {
                now + dc as f64 * give
            };
            let dc = (next - now) as f32;
            crust[k] += dc;
            bedrock[k] += airy * dc;
            moved.carried += 1;
            moved.most = moved.most.max(dc.abs() as f64);
        }
        if lost > 0.0 {
            let fell = (airy as f64 * lost / whole) as f32;
            for b in bedrock.iter_mut() {
                *b += fell;
            }
            moved.foundered += lost / whole;
        }
    }
    let (mut sum, mut whole) = (0.0f64, 0.0f64);
    for (k, &c) in crust.iter().enumerate() {
        let a = grid.cell_area(grid.cell(k));
        sum += c as f64 * a;
        whole += a;
    }
    moved.thickest = crust.iter().fold(0.0f32, |m, &c| m.max(c)) as f64;
    moved.mean = sum / whole;
    moved
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::convection::Convection;
    use crate::geom::Turn;
    use crate::thermal::ThermalHistory;

    /// F184: after an overturn Venus's crust reaches below its regrowing lid, and the flow
    /// piles it over the downwellings, as deep as the eclogite's cap at the strongest, and
    /// thins it over the upwellings, the whole kept; the ground rises over the thickened
    /// crust by its root's share. Mars's stagnant lid, thicker than its crust, carries none.
    #[test]
    fn the_flow_piles_the_crust_over_the_downwellings_beneath_a_thin_lid() {
        let grid = CubeSphere::new(16);
        let run = |name: &str, from: Option<f64>| {
            let (g, _) = crate::presets::preset(name).unwrap().to_graph();
            let thermal = ThermalHistory::from_graph(&g);
            let mantle = thermal.mantle.clone().unwrap();
            let convection = Convection::from_graph(&g, &thermal).unwrap();
            let support = Support {
                convection,
                turn: Turn::default(),
            };
            let t0 = from.unwrap_or(2.0);
            let c0 = g.f64("tectonics.initial_crust_thickness") as f32;
            let mut crust = vec![c0; grid.cell_count()];
            let mut bedrock = vec![0.0f32; grid.cell_count()];
            let moved = drag(
                grid,
                &mut crust,
                &mut bedrock,
                &support,
                &mantle,
                t0,
                t0 + 0.05,
            );
            (support, mantle, crust, bedrock, moved, c0, t0)
        };
        let (_, mantle, crust, bedrock, moved, c0, _) = run("Mars", None);
        assert!(mantle.stiff_lid_at(2.0) > c0 as f64);
        assert_eq!(moved.carried, 0);
        assert!(crust.iter().all(|&c| c == c0) && bedrock.iter().all(|&b| b == 0.0));

        let (g, _) = crate::presets::preset("Venus").unwrap().to_graph();
        let overturn = ThermalHistory::from_graph(&g)
            .mantle
            .unwrap()
            .last_overturn()
            .unwrap();
        let (support, _, crust, bedrock, moved, c0, t0) = run("Venus", Some(overturn));
        assert!(moved.carried > 0 && moved.most > 100.0, "{moved:?}");
        assert!(
            (moved.thickest - CRUST_MAX_THICKNESS).abs() < 1.0,
            "{moved:?}"
        );
        // Thickened where the pattern stands low and thinned where it stands high, the
        // pattern read just after the overturn, when the crust moved most.
        let (mut low, mut high, mut sum, mut area) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
        for k in 0..grid.cell_count() {
            let cell = grid.cell(k);
            let p = support.pattern(t0 + 0.005, grid.center(cell));
            let a = grid.cell_area(cell);
            let dc = (crust[k] - c0) as f64;
            if p < -1.0 {
                low += dc;
            } else if p > 1.0 {
                high += dc;
            }
            sum += dc * a;
            area += a;
            // The ground stands by its root's share above the datum, which fell with the
            // crust that foundered.
            let airy = (DENSITY_ROCK - DENSITY_CRUST_ROCK) / DENSITY_ROCK;
            assert!((bedrock[k] as f64 - airy * (dc + moved.foundered)).abs() < 1.0);
        }
        assert!(low > 0.0 && high < 0.0, "{low} {high}");
        // Kept but for what foundered past the cap.
        assert!(
            (sum / area + moved.foundered).abs() < 1.0,
            "{} {}",
            sum / area,
            moved.foundered
        );
        assert!(
            moved.foundered > 0.0 && moved.foundered < 0.1 * c0 as f64,
            "{}",
            moved.foundered
        );
    }
}
