// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! The tesserae's fabric from the strain their lid took (docs/parameters.md F187).
//!
//! Where the mantle's flow dragged the young lid past its strength, the lid broke and went
//! with the flow (crust::drag), shortening where the flow converged and stretching where
//! it spread, each way by the flow's own curvature. A strong layer shortened buckles into
//! folds across the shortening, some five times as long as the layer is thick, and one
//! stretched necks into graben across the stretching, three and a half times; the lid has
//! two such layers, the brittle layer its faults cut and the elastic lid beneath, so each
//! strain makes folds or graben at two scales, the brittle layer's a few kilometres apart
//! and the lid's some tens (the short and long folds and the ribbons of Ghent and Hansen).
//! A fold grows as the layer's length takes up the shortening, until its limbs dip as
//! steeply as the layer's friction lets them slip, or its relief weighs on the layer as
//! much as the layer's strength holds; a graben deepens with the stretching until its
//! relief weighs as much as the layer pulled apart holds. Each state of the flow made its
//! own sets, the later laid over the earlier.
//!
//! The fabric is drawn as a sum of oriented wave packets on the coarse cells, each with the
//! cell's strain: the packets' places, phases and a scatter in their spacing and relief are
//! drawn, so that the folds come as en echelon ridges of a few wavelengths each rather than
//! a grating, while where they stand, which way they run, how far apart and how high is the
//! strain's.

use rayon::prelude::*;

use crate::constants::*;
use crate::geom::{cross, dot, local_basis, normalize, offset};
use crate::grid::CubeSphere;
use crate::rng::Rng;
use crate::stages::episodic::LidStrain;
use crate::stages::stress::brittle_yields;

/// One wave of a packet: across `dir`, a fold or a graben.
#[derive(Clone, Copy, Debug)]
struct Wave {
    /// The direction the strain ran, which the wave's crests cross.
    dir: [f64; 3],
    wavenumber: f64,
    /// The fold's amplitude or the graben's depth, m.
    relief: f64,
    phase: f64,
    graben: bool,
    /// The packet's reach across its crests and along them, m.
    across: f64,
    along: f64,
}

/// A wave packet: the strain of one coarse cell in one layer of one state, its two waves
/// across the two principal directions.
#[derive(Clone, Debug)]
struct Packet {
    center: [f64; 3],
    family: usize,
    waves: Vec<Wave>,
}

/// The fabric of the tesserae: the packets by the coarse cell they stand in.
pub struct TesseraFabric {
    coarse: CubeSphere,
    body_radius: f64,
    /// Per family, the packets by coarse cell and how many coarse cells about a point they
    /// reach.
    packets: Vec<(Vec<Vec<Packet>>, i32)>,
    families: usize,
    /// Per family, the layer it folds or necks: the brittle layer or the elastic lid.
    layers: Vec<usize>,
    /// The relief a layer's strength holds, m.
    hold: f64,
    /// The highest relief any one wave reaches, m, for the fabric's bound.
    pub bound: f64,
}

/// The mean of a graben's profile over its spacing, which the profile is taken about so
/// that the packets sum to no net change.
const GRABEN_MEAN: f64 = 0.2180;

fn profile(phase: f64, graben: bool) -> f64 {
    if graben {
        // A trough across a third of the spacing, flat between, taken about its mean.
        -((phase.cos() - 0.5).max(0.0) / 0.5 - GRABEN_MEAN)
    } else {
        // A narrow crest and a broad trough, as the tesserae's ridges are narrower than the
        // troughs between (Hansen and Willis 1996), taken about its mean.
        4.0 / std::f64::consts::PI - 2.0 * (0.5 * phase).sin().abs()
    }
}

impl TesseraFabric {
    /// The fabric from the lid's strain sets on the coarse grid, drawn on a grid of cells
    /// `cell` m wide: the waves the grid cannot resolve are left out. `near` says which
    /// coarse cells the fabric is wanted about.
    pub fn new(
        coarse: &CubeSphere,
        sets: &[LidStrain],
        body_radius: f64,
        gravity: f64,
        cell: f64,
        seed: u64,
        near: &(dyn Fn(usize) -> bool + Sync),
    ) -> Self {
        let (yield_tension, yield_compression) = brittle_yields(gravity);
        let weight = DENSITY_CRUST_ROCK * gravity;
        let friction = ROCK_FRICTION;
        // Anderson's normal faults dip at 45 degrees plus half the friction angle.
        let dip = std::f64::consts::FRAC_PI_4 + 0.5 * friction.atan();
        let fold_cap = yield_compression / (weight * TESSERA_BUCKLING_RATIO);
        let graben_cap = yield_tension / (weight * GROOVE_NECKING_WAVELENGTH);
        let coarse_cell = coarse.delta() * body_radius;
        let rng = Rng::new(seed ^ 0x7E55_E2A0);
        let n = coarse.cell_count();
        let mut packets: Vec<(Vec<Vec<Packet>>, i32)> = Vec::new();
        let mut families = 0usize;
        let mut layers = Vec::new();
        let mut bound = 0.0f64;
        // The layers: the brittle layer the faults cut and the elastic lid.
        for (s, set) in sets.iter().enumerate() {
            for (l, &share) in [GROOVE_RIFT_LID_SHARE, 1.0].iter().enumerate() {
                let layer = share * set.elastic;
                let fold_length = TESSERA_BUCKLING_RATIO * layer;
                let graben_length = GROOVE_NECKING_WAVELENGTH * layer;
                // The grid draws a wave where it spans this many cells, fading in from two.
                let drawn = |length: f64| {
                    let x = ((length / cell - 2.0) / (GROOVE_MIN_CELLS - 2.0)).clamp(0.0, 1.0);
                    x * x * (3.0 - 2.0 * x)
                };
                let (fold_drawn, graben_drawn) = (drawn(fold_length), drawn(graben_length));
                if fold_drawn <= 0.0 && graben_drawn <= 0.0 {
                    continue;
                }
                let family = families;
                families += 1;
                layers.push(l);
                // As many packets as cover each place this many times over, so that each
                // shows as a segment of its own rather than merging into a grating.
                let shortest = fold_length.min(graben_length);
                let area = TESSERA_PACKET_WIDTH * TESSERA_PACKET_LENGTH * shortest * shortest;
                let per_cell = ((TESSERA_PACKET_DENSITY * coarse_cell * coarse_cell / area).ceil()
                    as usize)
                    .clamp(1, 16);
                let made: Vec<(usize, Vec<Packet>)> = (0..n)
                    .into_par_iter()
                    .filter(|&k| near(k))
                    .filter_map(|k| {
                        let [a, b, c] = set.strain[k];
                        if a == 0.0 && b == 0.0 && c == 0.0 {
                            return None;
                        }
                        let (a, b, c) = (a as f64, b as f64, c as f64);
                        let mean = 0.5 * (a + c);
                        let radius = (0.25 * (a - c) * (a - c) + b * b).sqrt();
                        let theta = 0.5 * (2.0 * b).atan2(a - c);
                        let center = coarse.center(coarse.cell(k));
                        let (east, north) = local_basis(center);
                        let axis = |angle: f64| {
                            let (sn, cs) = angle.sin_cos();
                            normalize([
                                east[0] * cs + north[0] * sn,
                                east[1] * cs + north[1] * sn,
                                east[2] * cs + north[2] * sn,
                            ])
                        };
                        // The greatest principal strain along `theta`, the least across it.
                        let id = ((s as u64) << 40) ^ ((l as u64) << 36) ^ k as u64;
                        let mut out = Vec::with_capacity(per_cell);
                        for j in 0..per_cell as u64 {
                            let draw = |i: u64| rng.uniform(id, 16 * j + i);
                            // The strike scattered about the strain's, as a fold set's
                            // segments step and turn.
                            let turn = TESSERA_STRIKE_SCATTER * (2.0 * draw(8) - 1.0);
                            let principal = [
                                (mean + radius, axis(theta + turn)),
                                (
                                    mean - radius,
                                    axis(theta + turn + std::f64::consts::FRAC_PI_2),
                                ),
                            ];
                            let place = offset(
                                center,
                                east,
                                north,
                                0.5 * coarse.delta() * draw(0).sqrt(),
                                std::f64::consts::TAU * draw(1),
                            );
                            let mut waves = Vec::with_capacity(2);
                            for (w, &(strain, dir)) in principal.iter().enumerate() {
                                let graben = strain > 0.0;
                                let (length, drawn) = if graben {
                                    (graben_length, graben_drawn)
                                } else {
                                    (fold_length, fold_drawn)
                                };
                                if drawn <= 0.0 || strain == 0.0 {
                                    continue;
                                }
                                let length = length
                                    * 2f64.powf(
                                        TESSERA_SPACING_SCATTER * (2.0 * draw(2 + w as u64) - 1.0),
                                    );
                                let relief = if graben {
                                    // Each graben takes the stretching across its spacing
                                    // on two faults, sinking as they slip.
                                    (0.5 * strain * length * dip.tan()).min(graben_cap)
                                } else {
                                    // A fold takes up the shortening in its layer's length,
                                    // to its limbs' friction or its layer's strength.
                                    ((length / std::f64::consts::PI) * (-strain).sqrt())
                                        .min(friction * length / std::f64::consts::TAU)
                                        .min(fold_cap)
                                };
                                let relief = relief
                                    * drawn
                                    * (1.0
                                        + TESSERA_RELIEF_SCATTER
                                            * (2.0 * draw(4 + w as u64) - 1.0));
                                // The wave's direction moved to the packet's place.
                                let dir = normalize(cross(cross(place, dir), place));
                                waves.push(Wave {
                                    dir,
                                    wavenumber: std::f64::consts::TAU / length,
                                    relief,
                                    phase: std::f64::consts::TAU * draw(6 + w as u64),
                                    graben,
                                    across: TESSERA_PACKET_WIDTH * length,
                                    along: TESSERA_PACKET_LENGTH * length,
                                });
                            }
                            // The folds of one layer together weigh on it no more than
                            // its strength holds.
                            let folded: f64 =
                                waves.iter().filter(|w| !w.graben).map(|w| w.relief).sum();
                            let most = fold_cap * (1.0 + TESSERA_RELIEF_SCATTER);
                            if folded > most {
                                for w in waves.iter_mut().filter(|w| !w.graben) {
                                    w.relief *= most / folded;
                                }
                            }
                            if !waves.is_empty() {
                                out.push(Packet {
                                    center: place,
                                    family,
                                    waves,
                                });
                            }
                        }
                        Some((k, out))
                    })
                    .collect();
                let mut cells: Vec<Vec<Packet>> = vec![Vec::new(); n];
                let mut widest = 0.0f64;
                for (k, list) in made {
                    for p in &list {
                        for w in &p.waves {
                            widest = widest.max(w.along).max(w.across);
                            bound = bound.max(w.relief);
                        }
                    }
                    cells[k] = list;
                }
                let reach = ((PACKET_CUTOFF * widest + coarse_cell) / coarse_cell).ceil() as i32;
                packets.push((cells, reach));
            }
        }
        Self {
            coarse: *coarse,
            body_radius,
            packets,
            families,
            layers,
            hold: fold_cap,
            bound,
        }
    }

    pub fn is_empty(&self) -> bool {
        self.families == 0
    }

    /// The fabric's relief at a direction, m.
    pub fn at(&self, d: [f64; 3]) -> f64 {
        if self.families == 0 {
            return 0.0;
        }
        // Per family and wave, the packets' weighted sum and their weights' squares.
        let mut sums = vec![[0.0f64; 4]; self.families];
        let (home, _, _) = self.coarse.locate(d);
        let r = self.body_radius;
        for (cells, reach) in &self.packets {
            for di in -reach..=*reach {
                for dj in -reach..=*reach {
                    let c = self.coarse.neighbor(home, di, dj);
                    for p in &cells[self.coarse.index(c)] {
                        let v = [
                            (d[0] - p.center[0]) * r,
                            (d[1] - p.center[1]) * r,
                            (d[2] - p.center[2]) * r,
                        ];
                        let normal = p.center;
                        for (w, wave) in p.waves.iter().enumerate() {
                            let u = dot(v, wave.dir);
                            let side = cross(normal, wave.dir);
                            let s = dot(v, side);
                            let q = (u / wave.across).powi(2) + (s / wave.along).powi(2);
                            if q > PACKET_CUTOFF * PACKET_CUTOFF {
                                continue;
                            }
                            let weight = (-0.5 * q).exp();
                            let e = &mut sums[p.family];
                            e[2 * w] += weight
                                * wave.relief
                                * profile(wave.wavenumber * u + wave.phase, wave.graben);
                            e[2 * w + 1] += weight * weight;
                        }
                    }
                }
            }
        }
        // Each layer's relief, every set it took laid together, weighs on it no more than
        // its strength holds: past that the relief its sets would add gives way.
        let mut layer = [0.0f64; 2];
        for (f, e) in sums.iter().enumerate() {
            layer[self.layers[f]] += e[0] / (e[1] + 1.0).sqrt() + e[2] / (e[3] + 1.0).sqrt();
        }
        layer
            .iter()
            .map(|&v| self.hold * (v / self.hold).tanh())
            .sum()
    }
}

/// A packet's reach in its widths: its weight is negligible beyond.
const PACKET_CUTOFF: f64 = 2.5;

#[cfg(test)]
mod tests {
    use super::*;

    /// F187: a lid shortened one way folds across that way, its crests running along the
    /// other, as far apart as five times its layer, no higher than its strength holds; one
    /// stretched one way has graben across it.
    #[test]
    fn a_strained_lid_folds_across_its_shortening_and_necks_across_its_stretching() {
        let coarse = CubeSphere::new(64);
        let r = 6.05e6;
        let n = coarse.cell_count();
        let elastic = 8e3;
        let fold = |strain: [f32; 3]| LidStrain {
            time: 0.0,
            elastic,
            strain: vec![strain; n],
            pressed: vec![0.0; n],
            pulled: vec![0.0; n],
        };
        let at = coarse.center(coarse.cell(coarse.index(crate::grid::Cell {
            face: 0,
            i: 32,
            j: 32,
        })));
        let (east, north) = local_basis(at);
        let cell = 2e3;
        let near = |k: usize| crate::geom::angle_between(coarse.center(coarse.cell(k)), at) < 0.05;
        // Shortened east-west: the relief varies east-west and hardly north-south.
        let fabric =
            TesseraFabric::new(&coarse, &[fold([-0.2, 0.0, 0.0])], r, 8.87, cell, 7, &near);
        // The relief's total variation along a line from the centre at an azimuth.
        let walk = |fabric: &TesseraFabric, phi: f64| -> f64 {
            let mut sum = 0.0;
            let mut last = fabric.at(at);
            for i in 1..400 {
                let h = fabric.at(offset(at, east, north, i as f64 * 250.0 / r, phi));
                sum += (h - last).abs();
                last = h;
            }
            sum
        };
        let across = walk(&fabric, 0.0);
        let along = walk(&fabric, std::f64::consts::FRAC_PI_2);
        assert!(across > 2.0 * along, "across {across}, along {along}");
        let (_, yc) = brittle_yields(8.87);
        let cap = yc / (DENSITY_CRUST_ROCK * 8.87 * TESSERA_BUCKLING_RATIO);
        assert!(
            fabric.bound <= cap * (1.0 + TESSERA_RELIEF_SCATTER) + 1.0,
            "{}",
            fabric.bound
        );
        // Three sets shortened every way, laid together: each layer's relief no more than
        // its strength holds (F189).
        let sets = [
            fold([-0.3, 0.05, -0.2]),
            fold([-0.25, -0.1, -0.3]),
            fold([-0.2, 0.0, -0.2]),
        ];
        let fabric = TesseraFabric::new(&coarse, &sets, r, 8.87, cell, 7, &near);
        let highest = (0..400)
            .map(|i| {
                fabric
                    .at(offset(at, east, north, i as f64 * 250.0 / r, 0.7))
                    .abs()
            })
            .fold(0.0f64, f64::max);
        assert!(
            highest <= 2.0 * cap + 1.0,
            "{highest} against {cap} a layer"
        );
        // Stretched: graben, the relief's low side deeper than its high side is high.
        let fabric = TesseraFabric::new(&coarse, &[fold([0.3, 0.0, 0.0])], r, 8.87, cell, 7, &near);
        let mut hs: Vec<f64> = (0..400)
            .map(|i| fabric.at(offset(at, east, north, i as f64 * 250.0 / r, 0.0)))
            .collect();
        hs.sort_by(f64::total_cmp);
        let median = hs[200];
        assert!(
            median - hs[4] > 1.5 * (hs[395] - median),
            "{} {} {}",
            hs[4],
            median,
            hs[395]
        );
    }
}
