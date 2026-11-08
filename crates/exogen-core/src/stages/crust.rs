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
use crate::convection::{Support, cells, drag_stress};
use crate::geom::{local_basis, offset};
use crate::grid::CubeSphere;
use crate::mantle::MantleHistory;
use crate::stages::stress::brittle_yields;

/// The strain one state of the mantle's flow put on the lid where its drag broke the lid and
/// carried it along (F187).
#[derive(Clone, Debug)]
pub struct StrainSet {
    /// The convection's state whose flow strained the lid.
    pub state: usize,
    /// The time the strain was taken at, Gyr, weighted by it.
    pub time: f64,
    weight: f64,
    /// Per cell, the horizontal strain in the cell's east and north basis, [ee, en, nn],
    /// extension positive.
    pub strain: Vec<[f32; 3]>,
    /// Per cell, the most the drag pressed and pulled the lid past its strength, Pa.
    pub pressed: Vec<f32>,
    pub pulled: Vec<f32>,
}

/// The strain the lid took from the flow between two times, Gyr, a set for each state of the
/// flow (F187): the lid an episodic body's last overturn keeps, from the overturn before.
#[derive(Clone, Debug)]
pub struct StrainRecord {
    pub from: f64,
    pub to: f64,
    pub sets: Vec<StrainSet>,
    /// The strain the lower crust's flow put on the lid as the relief relaxed (F188).
    pub collapse: Option<StrainSet>,
}

impl StrainRecord {
    pub fn new(from: f64, to: f64) -> Self {
        Self {
            from,
            to,
            sets: Vec::new(),
            collapse: None,
        }
    }

    /// Whether a step from `t0` to `t1`, Gyr, falls in the record's times, by its middle.
    pub fn covers(&self, t0: f64, t1: f64) -> bool {
        let t = 0.5 * (t0 + t1);
        t >= self.from && t < self.to
    }

    /// Adds the strain the lower crust's flow put on the lid over a step ending at `t`, Gyr,
    /// as the relief relaxed by `relaxed`, m (F188). The flow runs down the slope of its
    /// potential, `flow` (relaxation::Bands), so its strain is the potential's curvature over
    /// the crust's root share and thickness, `crust`, m, its trace the crust's thinning: the
    /// lid, whose crust thins as its lower part flows out from under the highs, stretches
    /// with it there and is pressed where it flows in. The potential is scaled so that its
    /// Laplacian is the relief relaxed, which the bands' spread of wavelengths about the
    /// ones they stand for blurs.
    pub fn add_collapse(
        &mut self,
        grid: CubeSphere,
        body_radius: f64,
        flow: &[f32],
        relaxed: &[f32],
        crust: &[f32],
        t: f64,
    ) {
        let airy = (DENSITY_ROCK - DENSITY_CRUST_ROCK) / DENSITY_ROCK;
        let padded = crate::grid::projection::PaddedField::new(&grid, flow);
        let h = grid.delta();
        let n = grid.cell_count();
        let curvatures: Vec<[f64; 3]> = (0..n)
            .into_par_iter()
            .map(|k| {
                curvature(
                    &|e| crate::grid::projection::sample_bilinear(&grid, &padded, e) as f64,
                    grid.center(grid.cell(k)),
                    h,
                )
            })
            .collect();
        let r2 = body_radius * body_radius;
        let (mut cross_sum, mut square) = (0.0f64, 0.0f64);
        for (c, &dh) in curvatures.iter().zip(relaxed) {
            let laplacian = (c[0] + c[2]) / r2;
            cross_sum -= laplacian * dh as f64;
            square += laplacian * laplacian;
        }
        if square <= 0.0 {
            return;
        }
        let scale = cross_sum / square / (airy * r2);
        let strain: Vec<[f32; 3]> = curvatures
            .iter()
            .zip(crust)
            .map(|(&[ee, en, nn], &c)| {
                if c <= 0.0 {
                    return [0.0; 3];
                }
                let f = scale / c as f64;
                [(ee * f) as f32, (en * f) as f32, (nn * f) as f32]
            })
            .collect();
        let set = self.collapse.get_or_insert_with(|| StrainSet {
            state: usize::MAX,
            time: 0.0,
            weight: 0.0,
            strain: vec![[0.0; 3]; n],
            pressed: vec![0.0; n],
            pulled: vec![0.0; n],
        });
        let mut weight = 0.0f64;
        for (s, inc) in set.strain.iter_mut().zip(&strain) {
            s[0] += inc[0];
            s[1] += inc[1];
            s[2] += inc[2];
            weight += (inc[0].abs() + inc[2].abs()) as f64;
        }
        set.time = (set.time * set.weight + t * weight) / (set.weight + weight).max(1e-300);
        set.weight += weight;
    }

    fn set(&mut self, state: usize, cells: usize) -> &mut StrainSet {
        match self.sets.iter().position(|s| s.state == state) {
            Some(i) => &mut self.sets[i],
            None => {
                self.sets.push(StrainSet {
                    state,
                    time: 0.0,
                    weight: 0.0,
                    strain: vec![[0.0; 3]; cells],
                    pressed: vec![0.0; cells],
                    pulled: vec![0.0; cells],
                });
                self.sets.last_mut().unwrap()
            }
        }
    }
}

/// The pattern's second derivatives at a direction, per radian², in the direction's east
/// and north basis: [ee, en, nn].
fn hessian(support: &Support, t: f64, d: [f64; 3], h: f64) -> [f64; 3] {
    curvature(&|e| support.pattern(t, e), d, h)
}

/// A field's second derivatives at a direction, per radian², read `h` radians about it, in
/// the direction's east and north basis: [ee, en, nn].
fn curvature(field: &dyn Fn([f64; 3]) -> f64, d: [f64; 3], h: f64) -> [f64; 3] {
    let (east, north) = local_basis(d);
    let at = |theta: f64, phi: f64| field(offset(d, east, north, theta, phi));
    let p0 = field(d);
    let q = std::f64::consts::FRAC_PI_4;
    let diagonal = h * std::f64::consts::SQRT_2;
    let ee = (at(h, 0.0) - 2.0 * p0 + at(h, 4.0 * q)) / (h * h);
    let nn = (at(h, 2.0 * q) - 2.0 * p0 + at(h, -2.0 * q)) / (h * h);
    let en = (at(diagonal, q) - at(diagonal, 3.0 * q) - at(diagonal, -q) + at(diagonal, -3.0 * q))
        / (4.0 * h * h);
    [ee, en, nn]
}

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
///
/// Where the flow's drag on the lid passes the lid's brittle strength, in compression over
/// the downwellings or in tension over the upwellings (F185), the lid breaks and goes with
/// the flow beneath it, and the strain it takes is written into `record` within its times
/// (F187): the flow's velocity runs down the pattern's slope, so its strain is the pattern's
/// curvature, of which the drag's thickening is the trace; its part alike every way is the
/// thickening the drag gives the crust, and the rest the curvature's.
#[allow(clippy::too_many_arguments)]
pub fn drag(
    grid: CubeSphere,
    crust: &mut [f32],
    bedrock: &mut [f32],
    support: &Support,
    mantle: &MantleHistory,
    t0: f64,
    t1: f64,
    mut record: Option<&mut StrainRecord>,
) -> Moved {
    let airy = ((DENSITY_ROCK - DENSITY_CRUST_ROCK) / DENSITY_ROCK) as f32;
    let (yield_tension, yield_compression) = brittle_yields(mantle.gravity);
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
        if let Some(rec) = record.as_deref_mut().filter(|r| t >= r.from && t < r.to) {
            let drag = drag_stress(mantle, t);
            let curvature = c.degree * (c.degree + 1.0);
            // The curvature read over twice the pattern's sampling, which its kinks would
            // otherwise make.
            let h = 2.0 * support.convection.spacing_at(t);
            let increments: Vec<Option<([f32; 3], f32)>> = (0..grid.cell_count())
                .into_par_iter()
                .map(|k| {
                    let dc = changes[k] as f64;
                    let carried = (crust[k] - stiff) as f64;
                    if dc == 0.0 || carried <= 0.0 {
                        return None;
                    }
                    let d = grid.center(grid.cell(k));
                    let p = support.pattern(t, d);
                    let stress = drag * p;
                    if stress < yield_tension && -stress < yield_compression {
                        return None;
                    }
                    let scale = if dc > 0.0 { take } else { give };
                    let areal = -dc * scale / carried;
                    let [ee, en, nn] = hessian(support, t, d, h);
                    let mean = 0.5 * (ee + nn);
                    let shear = -strain * scale / curvature;
                    Some((
                        [
                            (0.5 * areal + shear * (ee - mean)) as f32,
                            (shear * en) as f32,
                            (0.5 * areal + shear * (nn - mean)) as f32,
                        ],
                        stress as f32,
                    ))
                })
                .collect();
            if increments.iter().any(Option::is_some) {
                let set = rec.set(support.convection.state_at(t), grid.cell_count());
                let mut weight = 0.0f64;
                for (k, inc) in increments.iter().enumerate() {
                    if let Some(([a, b, cc], stress)) = inc {
                        let s = &mut set.strain[k];
                        s[0] += a;
                        s[1] += b;
                        s[2] += cc;
                        weight += (a.abs() + cc.abs()) as f64;
                        let (t_, c_) = (yield_tension as f32, yield_compression as f32);
                        set.pressed[k] = set.pressed[k].max(-stress - c_);
                        set.pulled[k] = set.pulled[k].max(stress - t_);
                    }
                }
                set.time = (set.time * set.weight + t * weight) / (set.weight + weight).max(1e-300);
                set.weight += weight;
            }
        }
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

    /// F187: where the flow's drag passes the young lid's strength, the lid takes the flow's
    /// strain: shortening every way at the median over the downwellings it was pressed
    /// past its strength, stretching over the upwellings it was pulled past it; nowhere
    /// else; and the strain is more than its part alike every way, the flow's curvature
    /// shaping it.
    #[test]
    fn the_lid_takes_the_flow_s_strain_where_the_drag_breaks_it() {
        let grid = CubeSphere::new(16);
        let (g, _) = crate::presets::preset("Venus").unwrap().to_graph();
        let thermal = ThermalHistory::from_graph(&g);
        let mantle = thermal.mantle.clone().unwrap();
        let support = Support {
            convection: Convection::from_graph(&g, &thermal).unwrap(),
            turn: Turn::default(),
        };
        let overturn = mantle.last_overturn().unwrap();
        let c0 = g.f64("tectonics.initial_crust_thickness") as f32;
        let mut crust = vec![c0; grid.cell_count()];
        let mut bedrock = vec![0.0f32; grid.cell_count()];
        let mut record = StrainRecord::new(overturn, overturn + 0.02);
        drag(
            grid,
            &mut crust,
            &mut bedrock,
            &support,
            &mantle,
            overturn,
            overturn + 0.05,
            Some(&mut record),
        );
        assert!(!record.sets.is_empty());
        let (mut pressed, mut pulled, mut deviatoric) = (Vec::new(), Vec::new(), 0usize);
        let mut strained = 0usize;
        for set in &record.sets {
            for k in 0..grid.cell_count() {
                let [a, b, c] = set.strain[k];
                if a == 0.0 && b == 0.0 && c == 0.0 {
                    assert!(set.pressed[k] <= 0.0 && set.pulled[k] <= 0.0);
                    continue;
                }
                strained += 1;
                assert!(set.pressed[k] > 0.0 || set.pulled[k] > 0.0, "{k}");
                let radius = (0.25 * (a - c) * (a - c) + b * b).sqrt();
                if radius > 0.5 * (a + c).abs() {
                    deviatoric += 1;
                }
                if set.pressed[k] > 0.0 && set.pulled[k] <= 0.0 {
                    pressed.push(a + c);
                } else if set.pulled[k] > 0.0 && set.pressed[k] <= 0.0 {
                    pulled.push(a + c);
                }
            }
        }
        let median = |v: &mut Vec<f32>| {
            v.sort_by(|a, b| a.total_cmp(b));
            v[v.len() / 2]
        };
        assert!(!pressed.is_empty() && !pulled.is_empty());
        assert!(median(&mut pressed) < 0.0 && median(&mut pulled) > 0.0);
        assert!(strained < grid.cell_count() * record.sets.len());
        assert!(deviatoric * 4 > strained, "{deviatoric} of {strained}");
    }

    /// F188: a high relaxing stretches the lid every way on its crest, and at its foot,
    /// where the lower crust flows in, stretches it around the high and presses it toward
    /// the high: graben radial there, folds concentric. Its trace is the crust's thinning.
    #[test]
    fn a_relaxing_high_stretches_its_crest_and_presses_its_foot() {
        let grid = CubeSphere::new(64);
        let r = 6.05e6;
        let center = grid.center(crate::grid::Cell {
            face: 0,
            i: 32,
            j: 32,
        });
        let sigma = 5.0 * grid.delta();
        let n = grid.cell_count();
        // A relaxed high of 100 m at its crest, and the potential whose Laplacian it is.
        let x2 = |k: usize| {
            let a = crate::geom::angle_between(grid.center(grid.cell(k)), center) / sigma;
            a * a
        };
        let flow: Vec<f32> = (0..n).map(|k| (-(-0.5 * x2(k)).exp()) as f32).collect();
        let relaxed: Vec<f32> = (0..n)
            .map(|k| {
                let q = x2(k);
                (-100.0 * (2.0 - q) * (-0.5 * q).exp()) as f32
            })
            .collect();
        let crust = vec![30e3f32; n];
        let mut record = StrainRecord::new(0.0, 1.0);
        record.add_collapse(grid, r, &flow, &relaxed, &crust, 0.5);
        let set = record.collapse.as_ref().unwrap();
        let principal = |i: u32, j: u32| {
            let k = grid.index(crate::grid::Cell { face: 0, i, j });
            let [a, b, c] = set.strain[k];
            let m = 0.5 * (a + c);
            let rad = (0.25 * (a - c) * (a - c) + b * b).sqrt();
            (m - rad, m + rad, a + c)
        };
        let airy = ((DENSITY_ROCK - DENSITY_CRUST_ROCK) / DENSITY_ROCK) as f32;
        let (least, greatest, trace) = principal(32, 32);
        assert!(least > 0.0 && greatest > 0.0, "{least} {greatest}");
        // The crest thins by its relaxed relief over the root's share.
        let thinning = 200.0 / airy / 30e3;
        assert!(
            (trace - thinning).abs() < 0.3 * thinning,
            "{trace} {thinning}"
        );
        let (least, greatest, _) = principal(32 + 9, 32);
        assert!(least < 0.0 && greatest > 0.0, "{least} {greatest}");
    }

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
                None,
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
