// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! Ring faults of a collapsing basin (PLAN.md Phase 4F, F69).
//!
//! A large crater collapses by drawing the weak substrate under the lithosphere into its
//! cavity, and the inflow drags the lithosphere toward the basin. Where the drag beats
//! the lithosphere's strength, the lithosphere fails in concentric normal faults: ring
//! tectonics (Melosh and McKinnon 1978; McKinnon and Melosh 1980). A thick, strong
//! lithosphere holds and the basin keeps only the rings of its own collapse; a thin one
//! over a substrate that flows breaks into rings of graben far beyond the rim, as around
//! Valhalla and Asgard on Callisto and Tyre and Callanish on Europa, and in Ganymede's
//! furrows, the ring system of a single impact across most of the moon (Hirata et al.
//! 2020). The process reads the lithosphere and the substrate at the impact's time, so a
//! body's ancient basins, struck while its lid was thin, carry ring systems its younger
//! ones do not.

use crate::constants::*;
use crate::crater::morphology::{Scaling, transient_diameter};
use crate::crater::{CellWrite, StampOutput, with_stamper};
use crate::fields::LayerSet;
use crate::geom::{angle_between, local_basis, offset};
use crate::grid::CubeSphere;
use crate::noise::GradientNoise;
use crate::rng::Rng;
use crate::stages::plume::{Trough, trough_depths};
use crate::stages::relaxation::Shell;
use rayon::prelude::*;
use std::collections::HashMap;

/// Melosh's (1982) measure of the drag of the inflow against the strength of a thin
/// plastic lithosphere, gamma = sqrt(3) S D_tr / (4 Y z): the basal drag S a share of
/// rho g H, H the transient cavity's depth, and the strength of a brittle lithosphere
/// Y z = (1/2) delta-sigma(z) z, its frictional strength rising with depth to
/// delta-sigma(z) = (1 - 1/k) rho g z for Andersonian normal faulting, k the ratio of the
/// principal stresses at failure for the friction of ice (Singer et al. 2023). The
/// density and the gravity cancel: gamma goes as the square of the transient diameter
/// over the lithosphere's thickness.
pub fn melosh_gamma(transient: f64, lithosphere: f64) -> f64 {
    if lithosphere <= 0.0 {
        return f64::INFINITY;
    }
    let mu = RING_FRICTION;
    let k = ((1.0 + mu * mu).sqrt() + mu).powi(2);
    let drag = RING_DRAG * RING_CAVITY_DEPTH;
    3f64.sqrt() * drag / (4.0 * 0.5 * (1.0 - 1.0 / k)) * (transient / lithosphere).powi(2)
}

/// The lithosphere through a basin's collapse, m: the ice colder than RING_WEAK_HOMOLOGOUS
/// of its melting point, on the shell's conductive profile as the thermal history lays
/// its elastic lid. On Callisto it was 19 km four billion years ago, where McKinnon and
/// Melosh (1980) read 15 to 20 km from Valhalla's rings.
pub fn collapse_lithosphere(shell: &Shell) -> f64 {
    let ts = shell.surface_temperature;
    let tm = shell.basal_temperature.max(ts + 1.0);
    let weak = RING_WEAK_HOMOLOGOUS * tm;
    shell.thickness.max(0.0) * ((weak - ts) / (tm - ts)).clamp(0.0, 1.0)
}

/// The dip of a normal fault in a lithosphere of the friction of ice: 45 degrees plus
/// half its friction angle (Anderson's theory; Singer et al. 2023), radians.
fn fault_dip() -> f64 {
    std::f64::consts::FRAC_PI_4 + 0.5 * RING_FRICTION.atan()
}

/// One basin's ring faults: the arcs of graben its collapse broke the lithosphere into.
#[derive(Clone, Debug)]
pub struct RingSystem {
    pub center: [f64; 3],
    pub time_gyr: f64,
    pub id: u64,
    pub transient: f64,
    pub lithosphere: f64,
    pub gamma: f64,
    /// The graben annulus's inner and outer distances from the centre, m.
    pub inner: f64,
    pub reach: f64,
    pub rings: usize,
    pub faults: Vec<Trough>,
}

/// What a basin needs for its ring faults: where and when it struck, how thick the
/// lithosphere was there then, and what the substrate did.
pub struct Basin<'a> {
    pub center: [f64; 3],
    pub diameter: f64,
    pub time_gyr: f64,
    pub seed: u64,
    /// The lithosphere through the collapse at the impact's time and place, m.
    pub lithosphere: f64,
    /// The share of the basin's own relief the substrate relaxed by flowing in during the
    /// collapse and its aftermath: what drags the lithosphere (0 where nothing flows).
    pub inflow: f64,
    /// The share of a wavelength's relief, m, that survives from the impact to the
    /// present, for the graben's own relaxation.
    pub kept: &'a (dyn Fn(f64) -> f64 + Sync),
}

/// The ring faults of a basin, or none when its lithosphere held (F69). The lithosphere
/// fails above RING_FAILURE of Melosh's gamma and breaks into an annulus of graben above
/// RING_GRABEN; between, the annulus reaches out by degrees. The graben run from
/// RING_INNER transient radii out to RING_REACH of them, as far as the inflow drags,
/// each a pair of normal faults meeting at the base of the lithosphere, spaced by the
/// stress each one relieves, deep by the strain it takes up, broken into arcs; inside
/// them, from RING_RIDGE_INNER transient radii, the blocks' raised edges stand as ridges.
pub fn ring_system(
    basin: &Basin,
    scaling: &Scaling,
    body_radius: f64,
    cell: f64,
) -> Option<RingSystem> {
    let transient = transient_diameter(basin.diameter, scaling);
    let gamma = melosh_gamma(transient, basin.lithosphere);
    let ramp = ((gamma - RING_FAILURE) / (RING_GRABEN - RING_FAILURE)).clamp(0.0, 1.0)
        * basin.inflow.clamp(0.0, 1.0);
    if ramp <= 0.0 {
        return None;
    }
    let a = 0.5 * transient;
    let inner = RING_INNER * a;
    // Nothing beyond the antipode's last spacing.
    let farthest = std::f64::consts::PI * body_radius - inner.min(0.1 * body_radius);
    let reach = (a * (1.0 + (RING_REACH - 1.0) * ramp)).min(farthest);
    // The graben's width at the annulus's inner edge: two faults meeting at the base of
    // the lithosphere.
    let width0 = 2.0 * basin.lithosphere / fault_dip().tan();
    let spacing0 = RING_SPACING * width0;
    if reach <= inner + 0.5 * spacing0 || width0 < 0.5 * cell {
        return None;
    }
    let rng = Rng::new(basin.seed ^ 0x2196_0000_0000);
    let wander = GradientNoise::new(basin.seed ^ 0x2196_0A00);
    let (east, north) = local_basis(basin.center);
    let tau = 2.0 * std::f64::consts::PI;
    let mut faults = Vec::new();
    // One ring of arcs at a distance: a graben where `depth` is positive, a ridge where it
    // is negative.
    let lay =
        |faults: &mut Vec<Trough>, ring: u64, r: f64, spacing: f64, width: f64, depth: f64| {
            let theta = r / body_radius;
            let circumference = tau * body_radius * theta.sin();
            if depth == 0.0 || circumference <= spacing {
                return;
            }
            let mut along = circumference * rng.uniform(ring, 0);
            let end = along + circumference;
            let mut arc = 0u64;
            while along < end - spacing && faults.len() < RING_MAX_FAULTS {
                let id = ring * 4096 + arc;
                let length = spacing * rng.range(id, 1, RING_ARC.0, RING_ARC.1);
                let length = length.min(end - along);
                let shift = RING_OFFSET * spacing * (2.0 * rng.uniform(id, 2) - 1.0);
                // A frame every half width: the arcs bend over hundreds of kilometres.
                let step = (0.5 * width).max(cell);
                let n = ((length / step).ceil() as usize).max(2);
                let points: Vec<[f64; 3]> = (0..=n)
                    .map(|i| {
                        let s = along + length * i as f64 / n as f64;
                        let phi = s / circumference * tau;
                        let nominal = offset(basin.center, east, north, theta, phi);
                        let swing = RING_WANDER
                            * spacing
                            * wander.fbm(nominal, body_radius / (2.0 * spacing), 2, 2.0, 0.5);
                        offset(
                            basin.center,
                            east,
                            north,
                            (r + shift + swing) / body_radius,
                            phi,
                        )
                    })
                    .collect();
                faults.push(Trough {
                    points,
                    spacing: length / n as f64,
                    length,
                    width: width * rng.range(id, 3, 0.8, 1.2),
                    depth: depth * rng.range(id, 4, 0.6, 1.4),
                    flank: if depth > 0.0 { RING_FLANK } else { 0.0 },
                    start_gyr: basin.time_gyr,
                    end_gyr: basin.time_gyr,
                    id: basin.seed ^ (0x2196_0000_0000_0000 | (ring << 20) | arc),
                });
                let gap = length * RING_GAP / (1.0 - RING_GAP) * rng.range(id, 5, 0.5, 1.5);
                along += length + gap;
                arc += 1;
            }
        };
    let mut ring = 0u64;
    // Inside the graben, where the strain is greatest and the faults form last and dip
    // toward the basin, the lithosphere breaks into blocks whose raised edges stand as
    // closely spaced ridges, no farther apart outward (Valhalla's sinuous ridge zone).
    let ridge_spacing = RING_RIDGE_SPACING * spacing0;
    let ridge_width = RING_RIDGE_WIDTH * width0;
    let ridge_height = RING_STRAIN * ridge_spacing * basin.inflow * (basin.kept)(2.0 * ridge_width);
    let mut r = RING_RIDGE_INNER * a;
    while r < inner - 0.5 * ridge_spacing {
        lay(
            &mut faults,
            ring,
            r,
            ridge_spacing,
            ridge_width,
            -ridge_height,
        );
        r += ridge_spacing * rng.range(ring, 6, 0.75, 1.25);
        ring += 1;
    }
    let mut r = inner;
    while r < reach && faults.len() < RING_MAX_FAULTS {
        let x = r / inner;
        let spacing = spacing0 * x.powf(RING_SPACING_GROWTH);
        let width = width0 * x.powf(-RING_WIDTH_DECAY);
        // The throw takes up the radial strain the inflow imposes here: each graben's two
        // walls shift the ground by about its depth.
        let depth = RING_STRAIN * spacing * x.powf(-RING_STRAIN_DECAY) * basin.inflow;
        let depth = depth * (basin.kept)(2.0 * width);
        if depth > 0.0 {
            lay(&mut faults, ring, r, spacing, width, depth);
        }
        r += spacing * rng.range(ring, 6, 0.75, 1.25);
        ring += 1;
    }
    if faults.is_empty() {
        return None;
    }
    Some(RingSystem {
        center: basin.center,
        time_gyr: basin.time_gyr,
        id: basin.seed,
        transient,
        lithosphere: basin.lithosphere,
        gamma,
        inner,
        reach,
        rings: ring as usize,
        faults,
    })
}

impl RingSystem {
    /// The cap the system's arcs, their walls and raised flanks lie within: its centre and
    /// angular radius.
    pub fn cap(&self, body_radius: f64) -> ([f64; 3], f64) {
        let mut reach: f64 = 0.0;
        for f in &self.faults {
            let margin = 2.5 * f.width / body_radius;
            for p in [
                f.points[0],
                f.points[f.points.len() / 2],
                f.points[f.points.len() - 1],
            ] {
                reach = reach.max(angle_between(self.center, p) + margin);
            }
        }
        (self.center, reach.min(std::f64::consts::PI))
    }
}

/// A basin's ring faults at the finishing resolution, its arcs drawn in parallel: the
/// ground dropped by each graben's profile and raised by each ridge's, whatever stood on it
/// lowered or raised with it (a graben's floor is the old surface let down between its
/// faults, not a new one). Where two arcs touch, their reliefs add.
pub fn system_compute(
    layers: &LayerSet,
    grid: &CubeSphere,
    body_radius: f64,
    system: &RingSystem,
    blur: f64,
) -> StampOutput {
    let bedrock = layers.bedrock.as_slice();
    let cell = grid.delta() * body_radius;
    let parts: Vec<Vec<(usize, f64)>> = system
        .faults
        .par_iter()
        .map(|f| {
            with_stamper(grid, |stamper| {
                trough_depths(grid, body_radius, f, 1.0, cell, blur, stamper)
            })
        })
        .collect();
    let mut total: HashMap<usize, f64> = HashMap::new();
    for part in parts {
        for (idx, depth) in part {
            *total.entry(idx).or_insert(0.0) += depth;
        }
    }
    let mut cells: Vec<CellWrite> = total
        .into_iter()
        .map(|(idx, depth)| CellWrite {
            idx: idx as u32,
            bedrock: (bedrock[idx] as f64 - depth) as f32,
            surface: false,
            age: 0.0,
            unit: 0,
            deep: 0.0,
            melt: 0.0,
            cover: 0.0,
        })
        .collect();
    cells.sort_by_key(|w| w.idx);
    StampOutput {
        cells,
        composition: None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::angle_between;

    fn callisto() -> Scaling {
        Scaling {
            simple_complex: 3_244.0,
            peak_ring: 27_577.0,
            multi_ring: 64_887.0,
            gravity: 1.235,
        }
    }

    fn basin<'a>(
        diameter: f64,
        lithosphere: f64,
        kept: &'a (dyn Fn(f64) -> f64 + Sync),
    ) -> Basin<'a> {
        Basin {
            center: [0.0, 0.0, 1.0],
            diameter,
            time_gyr: 0.3,
            seed: 7,
            lithosphere,
            inflow: 1.0,
            kept,
        }
    }

    #[test]
    fn melosh_gamma_matches_singer_s_evaluation() {
        // Singer et al. 2023 evaluate gamma for Tyre and Callanish at about 2: marginal
        // ring graben. Their numbers, a 20 km transient cavity on a 2.5 km brittle layer,
        // with this model's frictional strength give a few.
        let g = melosh_gamma(20e3, 2.5e3);
        assert!((2.0..6.0).contains(&g), "{g}");
        // Valhalla: a 450 km transient cavity on Callisto's lithosphere of 19 km four
        // billion years ago fails far beyond the graben threshold.
        assert!(melosh_gamma(450e3, 19e3) > 5.0 * RING_GRABEN);
    }

    #[test]
    fn a_thick_lithosphere_holds_and_a_thin_one_breaks_into_rings() {
        let keep = |_: f64| 1.0;
        let s = callisto();
        let r = 2_410e3;
        // A 1,000 km basin, Valhalla's size, on a lithosphere of 19 km: an annulus of graben from
        // about three transient radii out to about eight, in some tens of arcs.
        let v = ring_system(&basin(1_000e3, 19e3, &keep), &s, r, 2e3).unwrap();
        let a = 0.5 * v.transient;
        assert!((400e3..500e3).contains(&v.transient), "{}", v.transient);
        assert!((v.reach / a - RING_REACH).abs() < 0.01, "{}", v.reach / a);
        assert!(v.rings >= 6, "{} rings", v.rings);
        // Graben outside RING_INNER transient radii, ridges inside it from
        // RING_RIDGE_INNER.
        for f in &v.faults {
            let d = angle_between([0.0, 0.0, 1.0], f.points[0]) * r;
            if f.depth > 0.0 {
                assert!(d > 0.8 * v.inner && d < v.reach + 0.5 * a, "{d}");
            } else {
                assert!(d > 0.7 * RING_RIDGE_INNER * a && d < 1.2 * v.inner, "{d}");
            }
        }
        assert!(v.faults.iter().any(|f| f.depth < 0.0));
        // The same basin on a lithosphere of 150 km holds, gamma about half: no ring faults.
        assert!(ring_system(&basin(1_000e3, 150e3, &keep), &s, r, 2e3).is_none());
        // Nor does anything form where the substrate does not flow in.
        let mut b = basin(1_000e3, 19e3, &keep);
        b.inflow = 0.0;
        assert!(ring_system(&b, &s, r, 2e3).is_none());
    }

    #[test]
    fn graben_widen_apart_and_shallow_outward() {
        let keep = |_: f64| 1.0;
        let v = ring_system(&basin(1_000e3, 19e3, &keep), &callisto(), 2_410e3, 2e3).unwrap();
        let near = |f: &Trough| angle_between([0.0, 0.0, 1.0], f.points[0]) * 2_410e3;
        let inner: Vec<&Trough> = v
            .faults
            .iter()
            .filter(|f| f.depth > 0.0 && near(f) < 1.3 * v.inner)
            .collect();
        let outer: Vec<&Trough> = v
            .faults
            .iter()
            .filter(|f| near(f) > 0.8 * v.reach)
            .collect();
        let mean = |fs: &[&Trough], g: &dyn Fn(&Trough) -> f64| {
            fs.iter().map(|f| g(f)).sum::<f64>() / fs.len() as f64
        };
        assert!(!inner.is_empty() && !outer.is_empty());
        assert!(mean(&inner, &|f| f.depth) > 2.0 * mean(&outer, &|f| f.depth));
        assert!(mean(&inner, &|f| f.width) > mean(&outer, &|f| f.width));
        // The graben are about as wide as the faults are deep: twice the lithosphere over
        // the tangent of a 62 degree dip at the annulus's inner edge.
        let w = mean(&inner, &|f| f.width);
        assert!((17e3..24e3).contains(&w), "{w}");
    }
}
