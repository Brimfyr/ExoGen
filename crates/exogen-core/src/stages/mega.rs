// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! Mega-impacts (F168): the impacts too wide for the crater stamp, from three quarters of
//! the body's radius in crater radius up to a hemisphere, take a form of their own. The
//! basin excavates the crust inside its rim, the crust its melt sheet makes taking its
//! place, and throws the excavated crust beyond the rim, where it thickens the crust; the
//! ground stands where the crust's buoyancy holds it (Airy), the floor level where its
//! melt ponded, the older relief under the ejecta buried by its thickness. Where the
//! bombardment brings such an impact, the body has a hemispheric dichotomy (T5).

use rayon::prelude::*;

use crate::constants::*;
use crate::crater::stamp::{burial_factor, elliptical_radius};
use crate::crater::{CellWrite, CraterRecord, Levelling, StampOutput};
use crate::fields::LayerSet;
use crate::geom::{angle_between, azimuth_from, local_basis};
use crate::grid::CubeSphere;
use crate::units::{UNIT_BASIN_FLOOR, UNIT_DEEP_EJECTA};

/// Whether an impact of `diameter` m is a mega-impact on a body of `body_radius` m: wider
/// than the crater stamp draws.
pub fn is_mega(diameter: f64, body_radius: f64) -> bool {
    diameter > 2.0 * MAX_CRATER_RADIUS_FRACTION * body_radius
}

/// Whether a mega-impact's continuous ejecta renewed the ground at `at`: within
/// CONTINUOUS_EJECTA_MULTIPLE of its radii, as its form renews it.
pub fn renews(rec: &CraterRecord, at: [f64; 3], body_radius: f64) -> bool {
    let (east, north) = local_basis(rec.center);
    let r = angle_between(rec.center, at) * body_radius;
    let phi = azimuth_from(rec.center, east, north, at);
    elliptical_radius(r, phi, rec.azimuth, rec.ellipticity)
        <= CONTINUOUS_EJECTA_MULTIPLE * 0.5 * rec.diameter
}

/// Smooth step from 0 at `a` to 1 at `b`.
fn smooth(a: f64, b: f64, x: f64) -> f64 {
    let t = ((x - a) / (b - a)).clamp(0.0, 1.0);
    t * t * (3.0 - 2.0 * t)
}

/// The share of a cell's crust a mega-basin excavates, by its elliptical distance from the
/// centre in basin radii: all of it inside the foot of the basin's wall, none beyond its
/// rim.
pub fn excavated(x: f64) -> f64 {
    1.0 - smooth(1.0 - MEGA_WALL_SHARE, 1.0, x)
}

/// Emplace a mega-impact on the layers as they are now. `floor_crust` is the thickness of
/// the crust the basin's melt sheet makes, m, and `airy` the height the ground stands per
/// metre of crust, the crust's buoyancy in the mantle.
pub fn emplace(
    layers: &mut LayerSet,
    grid: &CubeSphere,
    body_radius: f64,
    rec: &CraterRecord,
    floor_crust: f64,
    airy: f64,
    levelling: &Levelling,
) {
    let n = grid.cell_count();
    let radius = 0.5 * rec.diameter;
    let (east, north) = local_basis(rec.center);
    // Each cell's elliptical distance from the centre, in basin radii, and its area, m².
    let geo: Vec<(f64, f64)> = (0..n)
        .into_par_iter()
        .map(|k| {
            let c = grid.cell(k);
            let d = grid.center(c);
            let r = angle_between(rec.center, d) * body_radius;
            let phi = azimuth_from(rec.center, east, north, d);
            (
                elliptical_radius(r, phi, rec.azimuth, rec.ellipticity) / radius,
                grid.cell_area(c) * body_radius * body_radius,
            )
        })
        .collect();
    // The body's mean ground and crust before the impact: the column the floor's new crust
    // floats against.
    let total: f64 = geo.iter().map(|(_, a)| a).sum();
    let mean_of = |v: &[f32]| -> f64 {
        v.iter().zip(&geo).map(|(&h, (_, a))| h as f64 * a).sum::<f64>() / total
    };
    let mean_ground = mean_of(layers.bedrock.as_slice());
    let mean_crust = mean_of(layers.crust_thickness.as_slice());
    // The crust excavated, m per cell, down to what the melt sheet makes.
    let removed: Vec<f64> = geo
        .iter()
        .zip(layers.crust_thickness.as_slice())
        .map(|(&(x, _), &h)| (excavated(x) * (h as f64 - floor_crust)).max(0.0))
        .collect();
    let volume: f64 = removed.iter().zip(&geo).map(|(r, (_, a))| r * a).sum();
    // The excavated crust lands beyond the rim, thinning with distance as a basin's ejecta
    // does.
    let spread: f64 = geo
        .iter()
        .filter(|(x, _)| *x > 1.0)
        .map(|(x, a)| x.powf(BASIN_EJECTA_EXPONENT) * a)
        .sum();
    let at_rim = if spread > 0.0 { volume / spread } else { 0.0 };
    let ejecta: Vec<f64> = geo
        .iter()
        .map(|&(x, _)| {
            if x > 1.0 {
                at_rim * x.powf(BASIN_EJECTA_EXPONENT)
            } else {
                0.0
            }
        })
        .collect();
    // The crust, and the ground its buoyancy holds up.
    {
        let time = rec.time_gyr as f32;
        let crust = layers.crust_thickness.as_mut_slice();
        let crust_age = layers.crust_age.as_mut_slice();
        let bedrock = layers.bedrock.as_mut_slice();
        for k in 0..n {
            let change = ejecta[k] - removed[k];
            crust[k] = (crust[k] as f64 + change) as f32;
            bedrock[k] = (bedrock[k] as f64 + airy * change) as f32;
            if excavated(geo[k].0) > 0.5 {
                crust_age[k] = time;
            }
        }
    }
    // The floor's melt ponds level and its new crust settles where its buoyancy holds it
    // against the body's mean column, whatever stood there before.
    let bed = layers.bedrock.as_slice();
    let level = mean_ground + airy * (floor_crust - mean_crust);
    // The ejecta levels the ground only over the hollows it can fill (F169): the relief
    // shorter than that width it buries by its thickness, the longer it drapes.
    let longest = levelling.length(at_rim).min(radius);
    let lows = crate::stages::relaxation::low_passes(
        &layers.bedrock,
        grid.delta() * body_radius,
        longest,
    );
    let time = rec.time_gyr as f32;
    let cells: Vec<CellWrite> = (0..n)
        .map(|k| {
            let x = geo[k].0;
            let b = bed[k] as f64;
            let inside = x <= 1.0;
            let bedrock = if inside {
                b + (level - b) * excavated(x)
            } else if lows.is_empty() {
                b
            } else {
                let l = crate::stages::relaxation::blend_levels(
                    &lows,
                    levelling.length(ejecta[k]),
                    k,
                ) as f64;
                let relief = b - l;
                l + relief * burial_factor(relief, ejecta[k])
            };
            CellWrite {
                idx: k as u32,
                bedrock: bedrock as f32,
                // The continuous ejecta renews the surface, as a crater's does.
                surface: x <= CONTINUOUS_EJECTA_MULTIPLE,
                age: time,
                unit: if inside {
                    UNIT_BASIN_FLOOR
                } else {
                    UNIT_DEEP_EJECTA
                },
                deep: 1.0,
                melt: if inside { excavated(x) as f32 } else { 0.0 },
                cover: 0.0,
            }
        })
        .collect();
    StampOutput {
        deposit: None,
        cells,
        ..Default::default()
    }
    .apply(layers);
}

/// A mega-impact in its place in the finish's age order (F168): its form is in the
/// resampled `base`, so inside its rim the ground is the base again, and beyond it what the
/// finish drew before it is buried by the ejecta's thickness against the base, its
/// continuous ejecta renewing the surface; whatever older the finish drew over the ground
/// the impact renewed is gone, as a crater's restore removes what its impact destroyed.
/// `excavated_crust` is the crust the basin excavated, m, from the starting crust to the
/// melt sheet's.
pub fn renew_compute(
    layers: &LayerSet,
    grid: &CubeSphere,
    body_radius: f64,
    rec: &CraterRecord,
    base: &[f32],
    excavated_crust: f64,
) -> StampOutput {
    let n = grid.cell_count();
    let radius = 0.5 * rec.diameter;
    let (east, north) = local_basis(rec.center);
    let geo: Vec<(f64, f64)> = (0..n)
        .into_par_iter()
        .map(|k| {
            let c = grid.cell(k);
            let d = grid.center(c);
            let r = angle_between(rec.center, d) * body_radius;
            let phi = azimuth_from(rec.center, east, north, d);
            (
                elliptical_radius(r, phi, rec.azimuth, rec.ellipticity) / radius,
                grid.cell_area(c) * body_radius * body_radius,
            )
        })
        .collect();
    let volume: f64 = geo
        .iter()
        .map(|&(x, a)| excavated(x) * excavated_crust * a)
        .sum();
    let spread: f64 = geo
        .iter()
        .filter(|(x, _)| *x > 1.0)
        .map(|(x, a)| x.powf(BASIN_EJECTA_EXPONENT) * a)
        .sum();
    let at_rim = if spread > 0.0 { volume / spread } else { 0.0 };
    let bed = layers.bedrock.as_slice();
    let time = rec.time_gyr as f32;
    let cells = (0..n)
        .into_par_iter()
        .map(|k| {
            let x = geo[k].0;
            let inside = x <= 1.0;
            let bedrock = if inside {
                base[k]
            } else {
                let relief = (bed[k] - base[k]) as f64;
                let thickness = at_rim * x.powf(BASIN_EJECTA_EXPONENT);
                (base[k] as f64 + relief * burial_factor(relief, thickness)) as f32
            };
            CellWrite {
                idx: k as u32,
                bedrock,
                surface: x <= CONTINUOUS_EJECTA_MULTIPLE,
                age: time,
                unit: if inside {
                    UNIT_BASIN_FLOOR
                } else {
                    UNIT_DEEP_EJECTA
                },
                deep: 1.0,
                melt: if inside { excavated(x) as f32 } else { 0.0 },
                cover: 0.0,
            }
        })
        .collect();
    StampOutput {
        deposit: None,
        cells,
        ..Default::default()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::crater::CraterClass;
    use crate::geom::offset;
    use crate::units::UNIT_PRIMORDIAL;

    fn record(center: [f64; 3], diameter: f64, ellipticity: f64) -> CraterRecord {
        CraterRecord {
            center,
            diameter,
            time_gyr: 0.6,
            azimuth: 0.3,
            angle: 0.4,
            ellipticity,
            class: CraterClass::MultiRing,
            from_loop: true,
            relief: 0.0,
            visible: true,
            on_lava: false,
            id: 7,
        }
    }

    /// F168: a mega-impact on a uniform crust 40 km thick, Mars's size: the crust inside the
    /// basin is the melt sheet's, the excavated crust lands beyond the rim, the whole kept;
    /// the floor stands as far below the old ground as the crust's buoyancy lets it, level
    /// however rough the old ground was; the rim's crust is the thickest, thinning away from
    /// it; the basin is longer along the impact's azimuth by its elongation; the floor is a
    /// basin floor of the impact's age and the ejecta deep ejecta.
    #[test]
    fn a_mega_impact_excavates_its_basin_and_throws_the_crust_beyond_its_rim() {
        let grid = CubeSphere::new(64);
        let body_radius = 3.39e6;
        let (crust0, floor, airy) = (40e3, 20e3, 0.12);
        let mut layers = LayerSet::new(grid);
        for k in 0..grid.cell_count() {
            let d = grid.center(grid.cell(k));
            layers.crust_thickness.as_mut_slice()[k] = crust0 as f32;
            // Rough old ground, some 1 km either way.
            layers.bedrock.as_mut_slice()[k] = (1e3 * (7.0 * d[0]).sin() * (5.0 * d[1]).cos()) as f32;
            layers.unit_id.as_mut_slice()[k] = UNIT_PRIMORDIAL;
        }
        let area = |k: usize| grid.cell_area(grid.cell(k)) * body_radius * body_radius;
        let volume = |l: &LayerSet| -> f64 {
            l.crust_thickness
                .as_slice()
                .iter()
                .enumerate()
                .map(|(k, &h)| h as f64 * area(k))
                .sum()
        };
        let before = volume(&layers);
        let center = crate::geom::normalize([0.2, 0.5, 0.8]);
        let rec = record(center, 1.8 * body_radius, 1.4);
        emplace(&mut layers, &grid, body_radius, &rec, floor, airy, &Levelling::rings());
        let after = volume(&layers);
        assert!((after - before).abs() < 1e-3 * before, "{before:.4e} {after:.4e}");

        let at = |theta: f64, phi: f64| {
            let (east, north) = local_basis(center);
            grid.index(grid.locate(offset(center, east, north, theta, phi)).0)
        };
        let radius = 0.5 * rec.diameter / body_radius;
        let c = at(0.0, 0.0);
        let crust = layers.crust_thickness.as_slice();
        let bed = layers.bedrock.as_slice();
        assert!((crust[c] as f64 - floor).abs() < 1.0, "{}", crust[c]);
        // The floor: level, and as low as the crust it lost holds it.
        let floor_cells: Vec<f64> = bed
            .iter()
            .enumerate()
            .filter(|(k, _)| angle_between(center, grid.center(grid.cell(*k))) < 0.6 * radius)
            .map(|(_, &b)| b as f64)
            .collect();
        let mean = floor_cells.iter().sum::<f64>() / floor_cells.len() as f64;
        let spread = floor_cells.iter().map(|b| (b - mean).abs()).fold(0.0, f64::max);
        assert!(spread < 1.0, "{spread}");
        assert!((mean + airy * (crust0 - floor)).abs() < 150.0, "{mean}");
        // The rim's crust the thickest, thinning outward, never below the old.
        let (az, e) = (rec.azimuth, rec.ellipticity);
        let across = az + std::f64::consts::FRAC_PI_2;
        let rim = crust[at(1.08 * radius / e.sqrt(), across)] as f64;
        let far = crust[at(2.5 * radius / e.sqrt(), across)] as f64;
        assert!(rim > far && far > crust0, "{rim} {far}");
        // Longer along the azimuth: the crust is the floor's out to near the rim along it,
        // and the old ground's well short of that across it.
        let along_in = crust[at(0.9 * radius * e.sqrt(), az)] as f64;
        let across_out = crust[at(1.05 * radius / e.sqrt() * 1.1, across)] as f64;
        assert!((along_in - floor).abs() < 1.0, "{along_in}");
        assert!(across_out > crust0, "{across_out}");
        let units = layers.unit_id.as_slice();
        let ages = layers.surface_age.as_slice();
        assert_eq!(units[c], UNIT_BASIN_FLOOR);
        assert_eq!(ages[c], 0.6);
        let out = at(1.5 * radius / e.sqrt(), across);
        assert_eq!(units[out], UNIT_DEEP_EJECTA);
        assert_eq!(units[at(3.0, 0.0)], UNIT_PRIMORDIAL);
    }
}
