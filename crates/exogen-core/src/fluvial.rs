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
/// the land, the land's mean. What they brought the sea laid on its floor under its level
/// (F414), m³; what full seas held back and the rivers carried on, m³, each time they carried
/// it; what was laid over the level where at the last step no sea had room left, m³; the
/// most the floor sank under its load (F415), m; and the mouths whose load their waves laid
/// along their coasts and what of it they moved from them, m³, through the span (F429).
#[derive(Clone, Copy, Debug, Default)]
pub struct Erosion {
    pub cut: f64,
    pub laid: f64,
    pub to_sea: f64,
    pub rate: f64,
    pub floor: f64,
    pub carried_on: f64,
    pub over: f64,
    pub sank: f64,
    pub mouths: usize,
    pub waved: usize,
    pub along: f64,
}

/// The rivers' erosion of the ground in `layers` on a grid of a body of `radius` m (F278), its
/// sea where `sea` at `sea_level`, m, under a year whose land sheds `runoff` and whose open
/// water gains `gain`, kg/m²/s at each cell, of a liquid of `density` kg/m³, over `hold`, the
/// lithosphere's flexural wavelength, m, a crust of `crust` kg/m³ over a mantle of `mantle`:
/// the ground rebounding by the crust's density over the mantle's of what it lost. What each
/// step's rivers bring the sea is laid on its floor from their mouths at the step's end (F413,
/// F414), the floor sinking under it (F415); the sea below the liquid's level as each step leaves
/// it, so that where the load fills a sea to its level the next step's rivers run over it, and
/// carry on with them what a full sea held back, into the basins and the seas beyond, steps of
/// carrying alone after the span's until none is left, the last spreading what no sea has room
/// for over all the sea it reached. Where the sea has `waves` under gravity `g`, each step's
/// mouths' load is laid along their coasts as their waves carry it before it is built out
/// (F429). Through today's climate's span, FLUVIAL_SPAN_YR in FLUVIAL_STEPS steps (`erode_over`).
pub fn erode(
    grid: CubeSphere,
    layers: &mut LayerSet,
    sea: &[bool],
    water: (&[f32], &[f32]),
    liquid: (f64, f64, f64),
    lithosphere: (f64, f64, f64),
    waves: Option<(&crate::waves::Waves, f64)>,
) -> Erosion {
    let span = (FLUVIAL_SPAN_YR, FLUVIAL_STEPS);
    erode_over(
        grid,
        layers,
        sea,
        water,
        liquid,
        lithosphere,
        waves,
        span,
        None,
    )
}

/// The rivers' erosion as `erode`'s through `span` years in `steps` steps, the drainage routed
/// again at each (F278): the finish's through today's climate, the loop's through each of its
/// steps (F437). Where `carry` is given, what full seas held back at their mouths, m³ at each
/// cell, is carried on from it and back into it, held for the rivers of the span after, as the
/// liquid stands again over the floor the load raised, where the finish's last step spreads it
/// over the sea it reached: through a step of the loop's history a small sea took too little of
/// what its rivers brought, and the rest stood over its few cells by kilometres.
#[allow(clippy::too_many_arguments)]
pub fn erode_over(
    grid: CubeSphere,
    layers: &mut LayerSet,
    sea: &[bool],
    water: (&[f32], &[f32]),
    (radius, density, sea_level): (f64, f64, f64),
    (hold, crust, mantle): (f64, f64, f64),
    waves: Option<(&crate::waves::Waves, f64)>,
    (span, steps): (f64, usize),
    carry: Option<&mut Vec<f64>>,
) -> Erosion {
    let rebound = crust / mantle;
    let n = grid.cell_count();
    let areas = grid.face_areas();
    let area: Vec<f64> = (0..n)
        .map(|k| areas[k % areas.len()] * radius * radius)
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
    let laying = sea_level.is_finite();
    let mut out = Erosion::default();
    let steps = steps.max(1);
    let dt = span / steps as f64;
    // The sea as each step leaves it; what full seas held back at their mouths, m³; and what
    // reached the sea, carried on again or not.
    let mut wet = sea.to_vec();
    // The seas a load filled to their level: land for the steps after, whatever they sink.
    let mut filled = vec![false; n];
    let mut pending = match carry.as_deref() {
        Some(c) if c.len() == n => c.clone(),
        _ => vec![0.0f64; n],
    };
    let holding = carry.is_some();
    let mut reached = 0.0;
    let (mut routing_time, mut laying_time, mut passes) = (0.0, 0.0, 0);
    for pass in 0..2 * steps {
        let eroding = pass < steps;
        let carried: f64 = pending.iter().sum();
        if !eroding && carried <= 0.0 {
            break;
        }
        passes += 1;
        let clock = std::time::Instant::now();
        let h = layers.height();
        let d = drain(grid, h.as_slice(), &wet, water, (radius, density));
        // What each land cell lost, carried with its water to its pit's basin, or the sea, with
        // what full seas held back; laid in the basins by their tree, and on the sea's floor
        // from its mouths.
        let mut cut = vec![0.0f64; n];
        if eroding {
            let z = incise(grid, h.as_slice(), &wet, &d, dt, body);
            for k in (0..n).filter(|&k| !wet[k] && d.lake_of[k] == NONE) {
                cut[k] = (h.as_slice()[k] as f64 - z[k]).max(0.0);
            }
        }
        // What a full sea held back runs on from where it was held: down the rivers from the
        // land it filled, into the sea's own floor where that is still under its level.
        let volumes: Vec<f64> = (0..n)
            .map(|k| cut[k] * area[k] + if wet[k] { 0.0 } else { pending[k] })
            .collect();
        out.carried_on += carried;
        let (mut mouths, mut astray) = (vec![0.0f64; n], 0.0);
        for (k, p) in pending.iter_mut().enumerate() {
            if wet[k] {
                mouths[k] += *p;
                reached += *p;
            }
            *p = 0.0;
        }
        let (laid, lost) = lay(
            h.as_slice(),
            &wet,
            &d,
            (&volumes, &area),
            (&mut mouths, &mut astray),
        );
        reached += lost;
        routing_time += clock.elapsed().as_secs_f64();
        if pass == 0 {
            let (lost, land) = (0..n)
                .filter(|&k| !wet[k])
                .fold((0.0, 0.0), |(s, w), k| (s + cut[k] * area[k], w + area[k]));
            out.rate = if land > 0.0 { lost / land / dt } else { 0.0 };
        }
        // What found no way to the sea, given to its mouths by their shares (none is lost); to
        // the sea's every cell by its area where none took any.
        let taken: f64 = mouths.iter().sum();
        if astray > 0.0 && taken > 0.0 {
            mouths.iter_mut().for_each(|m| *m *= 1.0 + astray / taken);
        } else if astray > 0.0 {
            let open: f64 = (0..n).filter(|&k| wet[k]).map(|k| area[k]).sum();
            for k in (0..n).filter(|&k| wet[k]) {
                mouths[k] = astray * area[k] / open.max(1e-30);
            }
        }
        // What a full sea cannot hold is held for the next step's rivers. At the last step of
        // carrying, it is given to the mouths whose seas had room, by their shares, as what
        // found no way is, and laid from them; what then finds no room is spread over all its
        // sea reached, over its level.
        let last = !holding && pass + 1 == 2 * steps;
        let clock = std::time::Instant::now();
        // The waves lay each mouth's load along its coast (F429).
        if let (Some((w, g)), true) = (waves, laying) {
            let (all, waved, along) = w.spread(
                grid,
                &wet,
                (h.as_slice(), sea_level),
                &mut mouths,
                (dt * year, steps),
                (radius, g),
            );
            out.mouths += all;
            out.waved += waved;
            out.along += along;
        }
        let floor = if laying && mouths.iter().any(|&m| m > 0.0) {
            let (mut floor, rests, full) = crate::stages::seafloor::wedges(
                grid,
                h.as_slice(),
                &wet,
                &mouths,
                (sea_level, radius),
                true,
            );
            for &c in &full {
                filled[c as usize] = true;
            }
            if last && !rests.is_empty() {
                let back: f64 = rests.iter().map(|r| r.1).sum();
                let open_sea: Vec<bool> = (0..n).map(|k| wet[k] && !filled[k]).collect();
                let mut open: Vec<f64> = (0..n)
                    .map(|k| if open_sea[k] { mouths[k] } else { 0.0 })
                    .collect();
                let room: f64 = open.iter().sum();
                if room > 0.0 {
                    open.iter_mut().for_each(|m| *m *= back / room);
                } else {
                    let span: f64 = (0..n).filter(|&k| open_sea[k]).map(|k| area[k]).sum();
                    for k in (0..n).filter(|&k| open_sea[k]) {
                        open[k] = back * area[k] / span.max(1e-30);
                    }
                }
                let now: Vec<f32> = (0..n)
                    .map(|k| (h.as_slice()[k] as f64 + floor[k]) as f32)
                    .collect();
                let (more, left, _) = crate::stages::seafloor::wedges(
                    grid,
                    &now,
                    &open_sea,
                    &open,
                    (sea_level, radius),
                    false,
                );
                out.over += left.iter().map(|r| r.1).sum::<f64>();
                for (f, m) in floor.iter_mut().zip(more) {
                    *f += m;
                }
            } else {
                for (m, rest) in rests {
                    pending[m] += rest;
                }
            }
            floor
        } else {
            vec![0.0f64; n]
        };
        laying_time += clock.elapsed().as_secs_f64();
        // The ground's rebound from its net loss over the flexural wavelength, and its sinking
        // under what the sea's floor took; then what it lost, from its loose cover first, and
        // what it was laid.
        let loss: Vec<f64> = (0..n).map(|k| cut[k] - laid[k]).collect();
        let lift: Vec<f64> = held(&loss).into_iter().map(|v| rebound * v).collect();
        let sink = sinking(&floor, &wet, (density, crust, mantle), &held);
        for k in 0..n {
            out.sank = out.sank.max(sink[k]);
            if wet[k] {
                out.floor += floor[k] * area[k];
                layers.regolith.as_mut_slice()[k] += floor[k] as f32;
                layers.bedrock.as_mut_slice()[k] -= sink[k] as f32;
                continue;
            }
            out.cut += cut[k] * area[k];
            out.laid += laid[k] * area[k];
            let regolith = layers.regolith.as_slice()[k] as f64;
            let take = regolith.min(cut[k]);
            layers.regolith.as_mut_slice()[k] = (regolith - take + laid[k]) as f32;
            layers.bedrock.as_mut_slice()[k] += (lift[k] - sink[k] - (cut[k] - take)) as f32;
        }
        // The sea below the liquid's level as the step leaves it (F26), ground within
        // LEVEL_TOLERANCE of it standing at it; a sea its load filled land from then on, a
        // basin that sinks as it fills taking the rest as a lake's would (F278), up to its pass.
        if laying {
            let h = layers.height();
            wet = h
                .as_slice()
                .iter()
                .zip(&filled)
                .map(|(&g, &f)| (g as f64) < sea_level - LEVEL_TOLERANCE && !f)
                .collect();
        }
    }
    if let Some(c) = carry {
        *c = pending;
    }
    // Of what the floor took, what lies over the level where no sea had room left.
    out.floor -= out.over;
    out.to_sea = reached - out.carried_on;
    if std::env::var_os("EXOGEN_DEBUG").is_some() {
        eprintln!(
            "seafloor: the sea took {:.6e} m³ net, {:.6e} laid on its floor, {:.4e} carried on from full seas, {:.4e} left over their level; {passes} steps, {routing_time:.1} s routing and cutting, {laying_time:.1} s laying on the floor",
            out.to_sea, out.floor, out.carried_on, out.over
        );
    }
    out
}

/// How far the ground sinks under what a step laid on the sea's floor (F415), m: `floor` the
/// depth laid at each cell, m, the sea where `wet`, under a liquid of `density` kg/m³, a crust
/// of `crust` kg/m³ over a mantle of `mantle` holding it, `held` its load spread over its
/// flexural wavelength. The load of what is laid where liquid stood is the crust's density less
/// the liquid's, laid as the ground it was cut from; it sinks the floor by that over the
/// mantle's less the liquid's, the liquid filling what it sinks, and the land beside it by that
/// over the mantle's alone.
fn sinking(
    floor: &[f64],
    wet: &[bool],
    (density, crust, mantle): (f64, f64, f64),
    held: &dyn Fn(&[f64]) -> Vec<f64>,
) -> Vec<f64> {
    if !floor.iter().any(|&f| f > 0.0) {
        return vec![0.0; floor.len()];
    }
    // The load as the crust's own thickness of it.
    let load: Vec<f64> = floor
        .iter()
        .map(|&f| f * (crust - density) / crust)
        .collect();
    held(&load)
        .into_iter()
        .zip(wet)
        .map(|(spread, &w)| spread * crust / if w { mantle - density } else { mantle })
        .collect()
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

/// The ground `volume` m³ the sea takes from the land at `from` (F413), carried down its
/// receivers to the cell of the sea it enters, added to `mouths`; to `astray` where its way
/// ends short of the sea.
fn to_mouth(
    d: &Drainage,
    sea: &[bool],
    from: u32,
    volume: f64,
    (mouths, astray): (&mut [f64], &mut f64),
) {
    let mut c = from;
    for _ in 0..sea.len() {
        if c == NONE {
            break;
        }
        if sea[c as usize] {
            mouths[c as usize] += volume;
            return;
        }
        let r = d.receiver[c as usize];
        if r == c {
            break;
        }
        c = r;
    }
    *astray += volume;
}

/// The ground each catchment moved laid in the flood's basins (F278), `moved` m³ at each cell
/// of cells of `area` m²: each pit's basin brought its catchment's, the basins settled as the water's lakes are but by
/// volume, from the tops down in the order they settle. Returns the depth laid at each cell, m,
/// and what reached the sea, m³; what reached it carried down the rivers to the cell of the sea
/// each enters, added to `mouths`, m³, and to `astray` where its way ends short of it (F413).
fn lay(
    height: &[f32],
    sea: &[bool],
    d: &Drainage,
    (moved, area): (&[f64], &[f64]),
    (mouths, astray): (&mut [f64], &mut f64),
) -> (Vec<f64>, f64) {
    let n = height.len();
    let basins = &d.basins;
    let count = basins.len();
    // Each basin's catchment's ground, and the ground below where its halves met: its area and
    // its heights' sum over it. What drains to the sea carried down its rivers, each cell
    // before the one it drains to.
    let mut brought = vec![0.0f64; count];
    let mut to_sea = 0.0;
    let mut carried = vec![0.0f64; n];
    for &k in d.order.iter().rev() {
        let k = k as usize;
        if sea[k] {
            continue;
        }
        let v = moved[k];
        match d.label[k] {
            0 | NONE => {
                to_sea += v;
                carried[k] += v;
                let r = d.receiver[k];
                if r == NONE || r as usize == k {
                    *astray += carried[k];
                } else if sea[r as usize] {
                    mouths[r as usize] += carried[k];
                } else {
                    carried[r as usize] += carried[k];
                }
            }
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
                                to_mouth(d, sea, beyond, rest, (&mut *mouths, &mut *astray));
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
                (hold, DENSITY_CRUST_ROCK, DENSITY_ROCK),
                None,
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

    /// F437: the rivers through a span of the loop's history erode alike whatever their steps:
    /// on the pitted island under a dry year, 200 Myr in steps of LOOP_FLUVIAL_STEP_YR cut
    /// within 15 percent of what steps four times finer cut, and more than twice what today's
    /// span cuts; what they cut laid or brought to the sea within a twentieth.
    #[test]
    fn the_history_s_rivers_wear_alike_whatever_their_steps() {
        let grid = CubeSphere::new(48);
        let top = [0.0, 0.0, 1.0];
        let noise = crate::noise::GradientNoise::new(7);
        let (start, sea) = ground(grid, |d| {
            let r = dot(d, top).clamp(-1.0, 1.0).acos();
            1_000.0 - 2_000.0 * r + 200.0 * noise.fbm(d, 60.0, 3, 2.0, 0.5)
        });
        let n = grid.cell_count();
        let (runoff, gain) = (vec![2e-7f32; n], vec![-1e-7f32; n]);
        let wear = |span: f64, steps: usize| {
            let mut layers = start.clone();
            erode_over(
                grid,
                &mut layers,
                &sea,
                (&runoff, &gain),
                (2e6, 1_000.0, 0.0),
                (2e5, DENSITY_CRUST_ROCK, DENSITY_ROCK),
                None,
                (span, steps),
                None,
            )
        };
        let steps = (2e8 / LOOP_FLUVIAL_STEP_YR).ceil() as usize;
        let (few, many) = (wear(2e8, steps), wear(2e8, 4 * steps));
        let today = wear(FLUVIAL_SPAN_YR, FLUVIAL_STEPS);
        println!(
            "200 Myr: cut {:.3e} m³ in {steps} steps, {:.3e} in four times as many; today's span {:.3e}; laid {:.3e} and {:.3e}, to the sea {:.3e} and {:.3e}",
            few.cut, many.cut, today.cut, few.laid, many.laid, few.to_sea, many.to_sea
        );
        assert!(
            today.cut > 0.0 && many.cut > 2.0 * today.cut,
            "{today:?} {many:?}"
        );
        assert!(
            ((few.cut - many.cut) / many.cut).abs() < 0.15,
            "{few:?} {many:?}"
        );
        for m in [few, many] {
            assert!(((m.laid + m.to_sea - m.cut) / m.cut).abs() < 0.05, "{m:?}");
        }
    }

    /// F437: a small sea its rivers fill through a long step of the loop's history, what it has no
    /// room for held at its mouths for the step after (`erode_over`'s `carry`): on a smooth
    /// plateau draining to a sea of a few cells, through 50 Myr in one step, what the rivers cut
    /// is laid, brought to the sea or held, within a twentieth; the ground rises nowhere over the
    /// plateau's highest.
    #[test]
    fn a_full_sea_holds_its_rivers_load_for_the_step_after() {
        let grid = CubeSphere::new(48);
        let top = [0.0, 0.0, 1.0];
        let (start, _) = ground(grid, |d| {
            let r = dot(d, top).clamp(-1.0, 1.0).acos();
            let pit = if r < 0.08 {
                -800.0 * (1.0 - (r / 0.08).powi(2))
            } else {
                0.0
            };
            500.0 + 200.0 * r + pit
        });
        let n = grid.cell_count();
        let height = start.height();
        let sea: Vec<bool> = height.as_slice().iter().map(|&h| h < 0.0).collect();
        let highest = height.as_slice().iter().copied().fold(f32::MIN, f32::max);
        let (runoff, gain) = (vec![2e-5f32; n], vec![-1e-6f32; n]);
        let mut layers = start.clone();
        let mut held = vec![0.0f64; n];
        let m = erode_over(
            grid,
            &mut layers,
            &sea,
            (&runoff, &gain),
            (2e6, 1_000.0, 0.0),
            (2e5, DENSITY_CRUST_ROCK, DENSITY_ROCK),
            None,
            (5e7, 1),
            Some(&mut held),
        );
        let kept: f64 = held.iter().sum();
        let most = layers
            .height()
            .as_slice()
            .iter()
            .copied()
            .fold(f32::MIN, f32::max);
        println!(
            "a sea of {} cells: cut {:.3e} m³, laid {:.3e}, to the sea {:.3e}, held {kept:.3e}, carried on {:.3e}; the ground's highest {most:.0} m where {highest:.0}",
            sea.iter().filter(|&&s| s).count(),
            m.cut,
            m.laid,
            m.to_sea,
            m.carried_on
        );
        assert!(
            ((m.laid + m.to_sea + kept - m.cut) / m.cut).abs() < 0.05,
            "{m:?} {kept}"
        );
        assert!(most <= highest, "{most} {highest}");
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
            (2e6, DENSITY_CRUST_ROCK, DENSITY_ROCK),
            None,
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

    /// F413: what the sea takes from a pitted island reaches it at its mouths: the cells of the
    /// sea its rivers and its full pits' spills enter, each beside the land, all of it, none
    /// astray.
    #[test]
    fn a_catchment_s_ground_reaches_the_sea_at_its_mouths() {
        let grid = CubeSphere::new(48);
        let top = [0.0, 0.0, 1.0];
        let noise = crate::noise::GradientNoise::new(7);
        let (layers, sea) = ground(grid, |d| {
            let r = dot(d, top).clamp(-1.0, 1.0).acos();
            1_000.0 - 2_000.0 * r + 200.0 * noise.fbm(d, 60.0, 3, 2.0, 0.5)
        });
        let n = grid.cell_count();
        let radius = 2e6;
        let (runoff, gain) = (vec![2e-5f32; n], vec![-1e-6f32; n]);
        let h = layers.height();
        let d = drain(
            grid,
            h.as_slice(),
            &sea,
            (&runoff, &gain),
            (radius, 1_000.0),
        );
        let area: Vec<f64> = (0..n)
            .map(|k| grid.cell_area(grid.cell(k)) * radius * radius)
            .collect();
        // A metre cut from every land cell no lake covers.
        let moved: Vec<f64> = (0..n)
            .map(|k| {
                if sea[k] || d.lake_of[k] != NONE {
                    0.0
                } else {
                    area[k]
                }
            })
            .collect();
        let (mut mouths, mut astray) = (vec![0.0f64; n], 0.0);
        let (_, to_sea) = lay(
            h.as_slice(),
            &sea,
            &d,
            (&moved, &area),
            (&mut mouths, &mut astray),
        );
        let taken: f64 = mouths.iter().sum();
        let beside_land = |k: usize| {
            let c = grid.cell(k);
            [
                (1, 0),
                (-1, 0),
                (0, 1),
                (0, -1),
                (1, 1),
                (1, -1),
                (-1, 1),
                (-1, -1),
            ]
            .iter()
            .any(|&(di, dj)| !sea[grid.index(grid.neighbor(c, di, dj))])
        };
        let count = (0..n).filter(|&k| mouths[k] > 0.0).count();
        let spills = d.outlets.len();
        println!(
            "the sea takes {to_sea:.4e} m³, {taken:.4e} at {count} mouths, {astray:.3e} astray; {spills} full lakes spill"
        );
        assert!(to_sea > 0.0 && spills > 0);
        assert!(astray == 0.0, "{astray}");
        assert!(((taken - to_sea) / to_sea).abs() < 1e-9, "{taken} {to_sea}");
        assert!((0..n).all(|k| mouths[k] == 0.0 || (sea[k] && beside_land(k))));
    }

    /// F415: the floor sinks under what is laid on it by its densities: where the plate holds a
    /// load where it lies, by the load's density less the liquid's over the mantle's less the
    /// liquid's, the liquid filling what it sinks; where the plate spreads it over a wavelength
    /// that spans the body, the sea's floor by that spread load over the mantle's less the
    /// liquid's and the land beside it by it over the mantle's alone.
    #[test]
    fn the_floor_sinks_under_its_load_by_its_densities() {
        let grid = CubeSphere::new(32);
        let n = grid.cell_count();
        let (crust, mantle, liquid) = (DENSITY_CRUST_ROCK, DENSITY_ROCK, 1_000.0);
        let densities = (liquid, crust, mantle);
        let wet: Vec<bool> = (0..n).map(|k| grid.center(grid.cell(k))[2] < 0.3).collect();
        let floor: Vec<f64> = wet.iter().map(|&w| if w { 10.0 } else { 0.0 }).collect();
        // A plate that holds each load where it lies, and one that spreads it over the body.
        let here = sinking(&floor, &wet, densities, &|v: &[f64]| v.to_vec());
        let spread = |v: &[f64]| vec![v.iter().sum::<f64>() / v.len() as f64; v.len()];
        let wide = sinking(&floor, &wet, densities, &spread);
        let local = 10.0 * (crust - liquid) / (mantle - liquid);
        let share = wet.iter().filter(|&&w| w).count() as f64 / n as f64;
        let mean = 10.0 * (crust - liquid) / crust * share;
        let (under, beside) = (mean * crust / (mantle - liquid), mean * crust / mantle);
        let at = |f: &[f64], w: bool| f[(0..n).find(|&k| wet[k] == w).unwrap()];
        println!(
            "held where it lies the floor sank {:.3} m, {local:.3} by the densities; spread, the floor {:.3} m where {under:.3}, the land {:.3} where {beside:.3}",
            at(&here, true),
            at(&wide, true),
            at(&wide, false)
        );
        for (k, &w) in wet.iter().enumerate() {
            let want = if w { local } else { 0.0 };
            assert!((here[k] - want).abs() < 1e-9, "{k}: {} {want}", here[k]);
            let want = if w { under } else { beside };
            assert!((wide[k] - want).abs() < 1e-9, "{k}: {} {want}", wide[k]);
        }
    }

    /// F414, F415: an island's rivers lay what they bring its sea on its floor, step by step,
    /// all of it: what the sea took, net of what full seas held back and the rivers carried on,
    /// is what its floor took and what was left over the level.
    #[test]
    fn the_sea_s_floor_takes_what_its_rivers_bring() {
        let grid = CubeSphere::new(48);
        let top = [0.0, 0.0, 1.0];
        let noise = crate::noise::GradientNoise::new(7);
        let (mut layers, sea) = ground(grid, |d| {
            let r = dot(d, top).clamp(-1.0, 1.0).acos();
            1_000.0 - 2_000.0 * r + 200.0 * noise.fbm(d, 60.0, 3, 2.0, 0.5)
        });
        let n = grid.cell_count();
        let (runoff, gain) = (vec![2e-5f32; n], vec![-1e-6f32; n]);
        let moved = erode(
            grid,
            &mut layers,
            &sea,
            (&runoff, &gain),
            (2e6, 1_000.0, 0.0),
            (2e5, DENSITY_CRUST_ROCK, DENSITY_ROCK),
            None,
        );
        println!("{moved:?}");
        assert!(moved.to_sea > 0.0 && moved.floor > 0.0);
        let kept = moved.floor + moved.over;
        assert!(
            ((kept - moved.to_sea) / moved.to_sea).abs() < 1e-9,
            "{moved:?}"
        );
        assert!(moved.sank > 0.0);
    }

    /// A pitted island whose slope holds a hollow `wide` radians across, 300 m below the sea's
    /// level, filled by the hypsometric liquid (F26), eroded under a plate of flexural wavelength
    /// `hold`, m: what the rivers moved, the hollow's lowest ground after, m, and what of the
    /// ground they cut is in neither the lakes nor the sea, m³.
    fn hollow_island(wide: f64, hold: f64) -> (Erosion, f32, f64) {
        let grid = CubeSphere::new(48);
        let top = [0.0, 0.0, 1.0];
        let hollow = crate::geom::normalize([0.25, 0.0, 1.0]);
        let noise = crate::noise::GradientNoise::new(7);
        let (mut layers, sea) = ground(grid, |d| {
            let r = dot(d, top).clamp(-1.0, 1.0).acos();
            let off = dot(d, hollow).clamp(-1.0, 1.0).acos();
            if off < wide {
                -300.0 * (1.0 - (off / wide).powi(2))
            } else {
                1_000.0 - 2_000.0 * r + 200.0 * noise.fbm(d, 60.0, 3, 2.0, 0.5)
            }
        });
        let n = grid.cell_count();
        let inside: Vec<usize> = (0..n)
            .filter(|&k| dot(grid.center(grid.cell(k)), hollow).acos() < wide)
            .collect();
        assert!(!inside.is_empty() && inside.iter().all(|&k| sea[k]));
        let (runoff, gain) = (vec![2e-5f32; n], vec![-1e-6f32; n]);
        let moved = erode(
            grid,
            &mut layers,
            &sea,
            (&runoff, &gain),
            (2e6, 1_000.0, 0.0),
            (hold, DENSITY_CRUST_ROCK, DENSITY_ROCK),
            None,
        );
        let h = layers.height();
        let lowest = inside
            .iter()
            .map(|&k| h.as_slice()[k])
            .fold(f32::INFINITY, f32::min);
        let lost = moved.cut - moved.laid - moved.to_sea;
        println!("{moved:?}; the hollow's lowest ground {lowest:.4} m; lost {lost:.3e} m³");
        (moved, lowest, lost)
    }

    /// F414: a hollow below the sea's level in an island's slope, filled by the hypsometric
    /// liquid (F26), takes its rivers' load until it stands at the level; then the rivers run
    /// over it, carry on what it held back into the basin it is above the level and the sea
    /// beyond, and none is left piled over the level: what left the land is what the floor and
    /// the basins took.
    #[test]
    fn a_full_sea_s_load_is_carried_on() {
        let (moved, lowest, lost) = hollow_island(0.04, 2e5);
        assert!(moved.carried_on > 0.0 && moved.over == 0.0, "{moved:?}");
        assert!(lowest as f64 >= -LEVEL_TOLERANCE, "{lowest}");
        assert!(lost.abs() < 1e-6 * moved.cut, "{lost}");
        let kept = moved.floor + moved.over;
        assert!(
            ((kept - moved.to_sea) / moved.to_sea).abs() < 1e-9,
            "{moved:?}"
        );
    }

    /// F414: a small hollow under a plate that holds its load where it lies sinks by most of
    /// what fills it; filled to the level once, it stays land, a basin that takes the rest as
    /// a lake's would up to its pass and spills the rest to the sea, none left over the level.
    /// The island's only small sea, it holds back once, the step it fills: what the rivers carry
    /// on is no more than a step brings the sea. Read as sea again each time it sank, it took its
    /// rivers' load back step after step.
    #[test]
    fn a_full_sea_that_sinks_is_filled_over_its_rim() {
        let (moved, lowest, lost) = hollow_island(0.03, 1.0);
        assert!(moved.carried_on > 0.0 && moved.over == 0.0, "{moved:?}");
        assert!(
            moved.carried_on < moved.to_sea / FLUVIAL_STEPS as f64,
            "{moved:?}"
        );
        assert!(lowest as f64 >= -LEVEL_TOLERANCE, "{lowest}");
        assert!(lost.abs() < 1e-6 * moved.cut, "{lost}");
    }
}
