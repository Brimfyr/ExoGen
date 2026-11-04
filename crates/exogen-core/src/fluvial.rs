// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! The rivers' erosion (F278).
//!
//! The rivers cut their beds by stream power, E = K A^m S (Howard 1994; Whipple and Tucker
//! 1999): A the area draining through a cell, weighted by its runoff against Earth's mean, S
//! its bed's slope to the cell its water runs to, solved implicitly down each river from its
//! mouth or its lake up (Braun and Willett 2013), so that no step cuts a bed below the one it
//! runs to. A full lake's spill cuts the pass it spills over. What the rivers cut lies in the
//! basins their water ends in, filling them from their deepest ground by the flood's basins as
//! the water's lakes settle (Fill-Spill-Merge by volume), or is lost to the sea. The ground
//! rebounds as the lithosphere's flexure gives (Gilchrist and Summerfield 1990): it lifts by
//! the crust's density over the mantle's of the ground the rivers took, averaged over its
//! flexural wavelength, the sea's floor losing none, so that a range wider than the wavelength
//! rebounds by the crust's share and lowers by the rest, while the plate holds up a loss
//! narrower than it: a coast's outlets lift the land beside them by the share of their loss the
//! wavelength spreads to it, the rest spread under the sea, whose floor the finish leaves where
//! it lies. Pits fill and their rims are cut until the land drains, the faster the wetter its
//! climate, over the span of today's climate the finish erodes by (until the loop carries the
//! erosion through the body's history, with its tectonics' uplift: the user, 2026-11-04).

use crate::constants::*;
use crate::fields::{Field, LayerSet};
use crate::geom::dot;
use crate::grid::CubeSphere;
use crate::rivers::{Drainage, NONE, drain};

/// What the rivers' erosion moved through its span (F278): the ground they cut, what of it they
/// laid in lakes and what they brought the sea, m³; and the rate they began at, m a year over
/// the land, the land's mean.
#[derive(Clone, Copy, Debug, Default)]
pub struct Erosion {
    pub cut: f64,
    pub laid: f64,
    pub to_sea: f64,
    pub rate: f64,
}

/// The rivers' erosion of the ground in `layers` on a grid of a body of `radius` m (F278), its
/// sea where `sea` at `sea_level`, m, under a year whose land sheds `runoff` and whose open
/// water gains `gain`, kg/m²/s at each cell, of a liquid of `density` kg/m³, the ground
/// rebounding by `rebound`, the crust's density over the mantle's, of what it lost over `hold`,
/// the lithosphere's flexural wavelength, m.
pub fn erode(
    grid: CubeSphere,
    layers: &mut LayerSet,
    sea: &[bool],
    water: (&[f32], &[f32]),
    (radius, density, sea_level): (f64, f64, f64),
    (hold, rebound): (f64, f64),
) -> Erosion {
    let n = grid.cell_count();
    let area: Vec<f64> = (0..n)
        .map(|k| grid.cell_area(grid.cell(k)) * radius * radius)
        .collect();
    let year = SECONDS_PER_DAY * 365.25;
    // The runoff the stream power's erodibility is taken at, kg/m²/s: the water a cell
    // carries over it the area that would drain through it there, m².
    let reference = FLUVIAL_REFERENCE_RUNOFF * density / year;
    let body = (radius, reference, sea_level);
    let cell = grid.delta() * radius;
    // The plate's flexure under a load: the load's mean over the flexural wavelength, the sea's
    // floor bearing none of what the land lost. Taken over the land alone, a coast's outlets,
    // cut to the sea's level in a step, lifted the dry land beside them by their whole loss
    // each step, by kilometres through the span.
    let held = |values: &[f64]| -> Vec<f64> {
        let mut f = Field::new(grid);
        for (x, &v) in f.as_mut_slice().iter_mut().zip(values) {
            *x = v as f32;
        }
        crate::stages::relaxation::low_pass(&f, cell, hold)
            .into_iter()
            .map(f64::from)
            .collect()
    };
    let mut out = Erosion::default();
    let dt = FLUVIAL_SPAN_YR / FLUVIAL_STEPS as f64;
    for step in 0..FLUVIAL_STEPS {
        let h = layers.height();
        let d = drain(grid, h.as_slice(), sea, water, (radius, density));
        let z = incise(grid, h.as_slice(), sea, &d, dt, body);
        // What each land cell lost, carried with its water to its pit's basin, or the sea; and
        // laid in the basins by their tree.
        let mut cut = vec![0.0f64; n];
        for k in (0..n).filter(|&k| !sea[k] && d.lake_of[k] == NONE) {
            cut[k] = (h.as_slice()[k] as f64 - z[k]).max(0.0);
        }
        let (laid, lost) = lay(h.as_slice(), sea, &d, &cut, &area);
        out.to_sea += lost;
        if step == 0 {
            let (lost, land) = (0..n)
                .filter(|&k| !sea[k])
                .fold((0.0, 0.0), |(s, w), k| (s + cut[k] * area[k], w + area[k]));
            out.rate = if land > 0.0 { lost / land / dt } else { 0.0 };
        }
        // The ground's rebound from its net loss over the flexural wavelength; then what it
        // lost, from its loose cover first, and what it was laid.
        let loss: Vec<f64> = (0..n).map(|k| cut[k] - laid[k]).collect();
        let lift: Vec<f64> = held(&loss).into_iter().map(|v| rebound * v).collect();
        for k in (0..n).filter(|&k| !sea[k]) {
            out.cut += cut[k] * area[k];
            out.laid += laid[k] * area[k];
            let regolith = layers.regolith.as_slice()[k] as f64;
            let take = regolith.min(cut[k]);
            layers.regolith.as_mut_slice()[k] = (regolith - take + laid[k]) as f32;
            layers.bedrock.as_mut_slice()[k] += (lift[k] - (cut[k] - take)) as f32;
        }
    }
    out
}

/// The rivers' beds after a step of `dt` years (F278), each cell cut by stream power down to the
/// bed its water runs to, implicitly from the rivers' mouths and lakes up: the bed at its mouth
/// the sea's level `sea_level`, at a lake the lake's, each full lake's spill cutting the pass it
/// spills over; `reference` the runoff the erodibility is taken at, kg/m²/s, on a body of
/// `radius` m.
fn incise(
    grid: CubeSphere,
    height: &[f32],
    sea: &[bool],
    d: &Drainage,
    dt: f64,
    (radius, reference, sea_level): (f64, f64, f64),
) -> Vec<f64> {
    let mut z: Vec<f64> = height.iter().map(|&h| h as f64).collect();
    let mut outlet: std::collections::HashMap<u32, (u32, f64)> = Default::default();
    for &(own, beyond, spills) in &d.outlets {
        if d.lake_of[own as usize] == NONE {
            outlet.insert(own, (beyond, spills));
        }
    }
    let reach = |a: usize, b: usize| {
        let cos = dot(grid.center(grid.cell(a)), grid.center(grid.cell(b)));
        (radius * cos.clamp(-1.0, 1.0).acos()).max(1.0)
    };
    for &k in &d.order {
        let k = k as usize;
        if sea[k] || d.lake_of[k] != NONE {
            continue;
        }
        let lifted = height[k] as f64;
        let (r, carried) = match outlet.get(&(k as u32)) {
            Some(&(beyond, spills)) => (beyond as usize, d.carried[k] + spills),
            None => (d.receiver[k] as usize, d.carried[k]),
        };
        if r == k {
            z[k] = lifted;
            continue;
        }
        let base = if sea[r] {
            sea_level
        } else if d.lake_of[r] != NONE {
            d.levels[d.lake_of[r] as usize] as f64
        } else {
            z[r]
        };
        let drained = carried / reference;
        let f = FLUVIAL_ERODIBILITY * drained.powf(FLUVIAL_AREA_EXPONENT) * dt / reach(k, r);
        z[k] = if lifted > base {
            (lifted + f * base) / (1.0 + f)
        } else {
            lifted
        };
    }
    z
}

/// What each basin holds of ground laid in it (F278), m³, at a level `level`, m, at or past
/// where its halves met: the ground below its halves' meeting, all of which lies under that
/// level, `area` m² of it standing `heights` m·m² over the datum, and its own cells' below the
/// level.
fn holds(
    basin: &crate::rivers::Basin,
    (below_area, below_heights): (f64, f64),
    level: f64,
    (height, area): (&[f32], &[f64]),
) -> f64 {
    let mut room = level * below_area - below_heights;
    for &c in &basin.own {
        let h = height[c as usize] as f64;
        if h >= level {
            break;
        }
        room += (level - h) * area[c as usize];
    }
    room.max(0.0)
}

/// The level the ground `volume` m³ laid in a basin's own range fills to (F278), from where its
/// halves met, over the ground below that meeting, `below_area` m², and its own cells by height,
/// at its pass at most.
fn fill_level(
    basin: &crate::rivers::Basin,
    below_area: f64,
    volume: f64,
    (height, area): (&[f32], &[f64]),
) -> f64 {
    let cap = basin.spill as f64;
    let mut level = if basin.joined.is_finite() {
        basin.joined as f64
    } else {
        basin
            .own
            .first()
            .map_or(cap, |&c| height[c as usize] as f64)
    };
    let (mut under, mut left, mut i) = (below_area, volume, 0usize);
    loop {
        while i < basin.own.len() && height[basin.own[i] as usize] as f64 <= level {
            under += area[basin.own[i] as usize];
            i += 1;
        }
        let next = basin
            .own
            .get(i)
            .map_or(cap, |&c| (height[c as usize] as f64).min(cap));
        let room = (next - level) * under;
        if left <= room || under <= 0.0 || next >= cap {
            let rise = if under > 0.0 { left / under } else { 0.0 };
            return (level + rise).min(cap);
        }
        left -= room;
        level = next;
    }
}

/// The ground each catchment cut laid in the flood's basins (F278), `cut` m at each cell: each
/// pit's basin brought its catchment's, the basins settled as the water's lakes are but by
/// volume, from the tops down in the order they settle. Returns the depth laid at each cell, m,
/// and what reached the sea, m³.
fn lay(height: &[f32], sea: &[bool], d: &Drainage, cut: &[f64], area: &[f64]) -> (Vec<f64>, f64) {
    let n = height.len();
    let basins = &d.basins;
    let count = basins.len();
    // Each basin's catchment's ground, and the ground below where its halves met: its area and
    // its heights' sum over it.
    let mut brought = vec![0.0f64; count];
    let mut to_sea = 0.0;
    for k in (0..n).filter(|&k| !sea[k]) {
        let v = cut[k] * area[k];
        match d.label[k] {
            0 | NONE => to_sea += v,
            l => brought[l as usize] += v,
        }
    }
    let mut below = vec![(0.0f64, 0.0f64); count];
    let mut full = vec![(0.0f64, 0.0f64); count];
    for b in 1..count {
        let own = basins[b].own.iter().fold((0.0, 0.0), |(a, h), &c| {
            (
                a + area[c as usize],
                h + height[c as usize] as f64 * area[c as usize],
            )
        });
        if let Some((x, y)) = basins[b].halves {
            let (x, y) = (x as usize, y as usize);
            brought[b] = brought[x] + brought[y];
            below[b] = (full[x].0 + full[y].0, full[x].1 + full[y].1);
        }
        full[b] = (below[b].0 + own.0, below[b].1 + own.1);
    }
    let holds_pit = |b: &crate::rivers::Basin, cell: u32| {
        (b.span.0..b.span.1).contains(&d.pit_place[d.label[cell as usize] as usize])
    };
    let mut levels: Vec<(u32, f64)> = Vec::new();
    let mut passed: Vec<Vec<(u32, f64)>> = vec![Vec::new(); count];
    let state = (height, area);
    for &root in &d.roots {
        let mut stack = vec![(root, std::mem::take(&mut passed[root as usize]))];
        while let Some((b, inflows)) = stack.pop() {
            let basin = &basins[b as usize];
            let given: f64 = brought[b as usize] + inflows.iter().map(|i| i.1).sum::<f64>();
            let at_join = if basin.joined.is_finite() {
                basin.joined as f64 * below[b as usize].0 - below[b as usize].1
            } else {
                0.0
            };
            match basin.halves {
                Some((x, y)) if given < at_join => {
                    let (mut into_x, mut into_y): (Vec<_>, Vec<_>) = inflows
                        .into_iter()
                        .partition(|&(cell, _)| holds_pit(&basins[x as usize], cell));
                    let over = |h: u32, into: &[(u32, f64)]| {
                        let half = &basins[h as usize];
                        let room = half.spill as f64 * full[h as usize].0 - full[h as usize].1;
                        brought[h as usize] + into.iter().map(|i| i.1).sum::<f64>() - room
                    };
                    let (ox, oy) = (over(x, &into_x), over(y, &into_y));
                    if ox > 0.0 {
                        into_y.push((basins[x as usize].pass.1, ox));
                    }
                    if oy > 0.0 {
                        into_x.push((basins[y as usize].pass.1, oy));
                    }
                    stack.push((x, into_x));
                    stack.push((y, into_y));
                }
                _ => {
                    let below_area = below[b as usize].0;
                    let in_range = given - at_join;
                    let cap = holds(basin, below[b as usize], basin.spill as f64, state);
                    if given >= cap && basin.spill.is_finite() {
                        levels.push((b, basin.spill as f64));
                        let rest = given - cap;
                        let beyond = basin.pass.1;
                        if b == root && rest > 0.0 {
                            let next = if beyond == NONE || sea[beyond as usize] {
                                NONE
                            } else {
                                match d.label[beyond as usize] {
                                    0 | NONE => NONE,
                                    l => d.top_of[l as usize],
                                }
                            };
                            if next == NONE || next == root {
                                to_sea += rest;
                            } else {
                                passed[next as usize].push((beyond, rest));
                            }
                        }
                    } else {
                        let level = fill_level(basin, below_area, in_range.max(0.0), state);
                        levels.push((b, level));
                    }
                }
            }
        }
    }
    // The depth laid: each filled basin's own cells below its level and all its halves'.
    let mut laid = vec![0.0f64; n];
    for &(b, level) in &levels {
        let mut stack = vec![b];
        while let Some(x) = stack.pop() {
            let basin = &basins[x as usize];
            for &c in &basin.own {
                let h = height[c as usize] as f64;
                if h < level {
                    laid[c as usize] = laid[c as usize].max(level - h);
                }
            }
            if let Some((p, q)) = basin.halves {
                stack.extend([p, q]);
            }
        }
    }
    (laid, to_sea)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A grid's layers whose bedrock is `f` of each cell's centre, the sea below none.
    fn ground(grid: CubeSphere, f: impl Fn([f64; 3]) -> f64) -> (LayerSet, Vec<bool>) {
        let mut layers = LayerSet::new(grid);
        for (k, b) in layers.bedrock.as_mut_slice().iter_mut().enumerate() {
            *b = f(grid.center(grid.cell(k))) as f32;
        }
        let sea = layers.bedrock.as_slice().iter().map(|&h| h < 0.0).collect();
        (layers, sea)
    }

    /// F278: on a wet island whose slope is pitted, the rivers fill the pits and cut their
    /// rims until the land drains: the lakes the pits held are gone; what they cut is laid or
    /// brought to the sea. The island lowers by its net loss less its rebound: by no less than
    /// the crust's share leaves, near that under a plate whose flexural wavelength is short
    /// beside the island, and the more under one whose wavelength spans it, the plate holding up
    /// a loss narrower than its wavelength.
    #[test]
    fn wet_land_drains_its_pits() {
        let grid = CubeSphere::new(48);
        let top = [0.0, 0.0, 1.0];
        let noise = crate::noise::GradientNoise::new(7);
        let (start, sea) = ground(grid, |d| {
            let r = dot(d, top).clamp(-1.0, 1.0).acos();
            1_000.0 - 2_000.0 * r + 200.0 * noise.fbm(d, 60.0, 3, 2.0, 0.5)
        });
        let n = grid.cell_count();
        let radius = 2e6;
        let (runoff, gain) = (vec![2e-5f32; n], vec![-1e-6f32; n]);
        let route = |l: &LayerSet| {
            crate::rivers::route(
                grid,
                l.height().as_slice(),
                &sea,
                (&runoff, &gain),
                (radius, 1_000.0),
            )
        };
        let area = |k: usize| grid.cell_area(grid.cell(k)) * radius * radius;
        let land: f64 = (0..n).filter(|&k| !sea[k]).map(area).sum();
        let land_mean = |l: &LayerSet| {
            let h = l.height();
            (0..n)
                .filter(|&k| !sea[k])
                .map(|k| h.as_slice()[k] as f64 * area(k))
                .sum::<f64>()
                / land
        };
        let rebound = DENSITY_CRUST_ROCK / DENSITY_ROCK;
        // The island's fall over its net loss, under a plate of flexural wavelength `hold`, m.
        let fall = |hold: f64| {
            let mut layers = start.clone();
            let moved = erode(
                grid,
                &mut layers,
                &sea,
                (&runoff, &gain),
                (radius, 1_000.0, 0.0),
                (hold, rebound),
            );
            let fell = land_mean(&start) - land_mean(&layers);
            let share = fell / ((moved.cut - moved.laid) / land);
            (layers, moved, share)
        };
        let (layers, moved, narrow) = fall(2e5);
        let (_, _, wide) = fall(2e6);
        let (before, after) = (route(&start), route(&layers));
        println!(
            "lakes {:.4} of the land before, {:.4} after; the land's mean {:.0} m, {:.0} after; cut {:.3e} m³, laid {:.3e}, to the sea {:.3e}; begun at {:.2} mm a year",
            before.summary.lakes,
            after.summary.lakes,
            land_mean(&start),
            land_mean(&layers),
            moved.cut,
            moved.laid,
            moved.to_sea,
            1e3 * moved.rate
        );
        assert!(before.summary.lakes > 0.005, "{:?}", before.summary);
        assert!(
            after.summary.lakes < 0.25 * before.summary.lakes,
            "{:?}",
            after.summary
        );
        println!(
            "the island fell by {narrow:.3} of its net loss under a plate whose wavelength is short beside it, {wide:.3} under one that spans it; the crust's share leaves {:.3}",
            1.0 - rebound
        );
        assert!(
            narrow > 0.99 * (1.0 - rebound) && narrow < 1.0 - 0.5 * rebound,
            "{narrow}"
        );
        assert!(wide > narrow && wide < 1.0, "{narrow} {wide}");
        assert!(
            ((moved.laid + moved.to_sea - moved.cut) / moved.cut).abs() < 0.05,
            "{moved:?}"
        );
    }

    /// F278: on a dry plain sloping to the sea, a broad basin 600 m deep stays closed through
    /// the span, its rivers bringing too little to fill it, while the small pits about it fill
    /// with what theirs bring: the land closed falls to the basin's catchment, no further.
    #[test]
    fn dry_land_keeps_its_broad_closed_basin() {
        let grid = CubeSphere::new(48);
        let top = [0.0, 0.0, 1.0];
        let noise = crate::noise::GradientNoise::new(7);
        let (mut layers, sea) = ground(grid, |d| {
            let r = dot(d, top).clamp(-1.0, 1.0).acos();
            let basin = if r < 0.2 {
                -600.0 * (1.0 - (r / 0.2).powi(2))
            } else {
                0.0
            };
            1_000.0 - 2_000.0 * r + basin + 200.0 * noise.fbm(d, 60.0, 3, 2.0, 0.5)
        });
        let n = grid.cell_count();
        let radius = 2e6;
        let land: f64 = (0..n)
            .filter(|&k| !sea[k])
            .map(|k| grid.cell_area(grid.cell(k)))
            .sum();
        let inside: f64 = (0..n)
            .filter(|&k| dot(grid.center(grid.cell(k)), top).acos() < 0.2)
            .map(|k| grid.cell_area(grid.cell(k)))
            .sum();
        let (runoff, gain) = (vec![2e-8f32; n], vec![-3e-5f32; n]);
        let route = |l: &LayerSet| {
            crate::rivers::route(
                grid,
                l.height().as_slice(),
                &sea,
                (&runoff, &gain),
                (radius, 1_000.0),
            )
        };
        let before = route(&layers);
        erode(
            grid,
            &mut layers,
            &sea,
            (&runoff, &gain),
            (radius, 1_000.0, 0.0),
            (2e6, DENSITY_CRUST_ROCK / DENSITY_ROCK),
        );
        let after = route(&layers);
        println!(
            "closed {:.3} of the land before, {:.3} after; the basin {:.3} of it",
            before.summary.closed,
            after.summary.closed,
            inside / land
        );
        assert!(
            before.summary.closed > inside / land,
            "{:?}",
            before.summary
        );
        assert!(
            after.summary.closed > 0.9 * inside / land
                && after.summary.closed < before.summary.closed,
            "{:?}",
            after.summary
        );
    }
}
