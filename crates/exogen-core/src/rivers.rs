// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! The rivers and the lakes (F277).
//!
//! The water the land sheds through its year, its rain and its snow's melt past what its ground
//! gives back to the air, runs down the terrain to the sea or to the lowest ground of a closed
//! basin. The basins are found as the priority flood finds them (Barnes et al. 2014, 2020): a
//! flood from the sea and from every pit, each cell taken in the order of its height by the
//! basin that reached it first and draining to the cell it was reached from, the lowest pass
//! between two basins where they meet; the basins merge at their passes in the order of the
//! passes' heights, each merger the basin of the two together above its pass, so that the
//! basins nest within one another to the one that spills to the sea. In each a lake stands
//! where what it is brought, its catchment's runoff and what spills into it, balances what its
//! surface loses to the air over what falls on it, at the year's steady state (Fill-Spill-Merge,
//! Barnes, Callaghan and Wickert 2021, its volumes the year's fluxes): filled to its pass, it
//! spills into the basin beyond, and two basins both filled to the pass between them fill
//! together past it. The rivers carry their catchments' runoff and the lakes' spill down to the
//! lakes and the sea.

use crate::fields::Field;
use crate::grid::CubeSphere;
use std::cmp::{Ordering, Reverse};
use std::collections::{BinaryHeap, HashMap};

/// A cell no basin has taken yet; a cell no lake covers.
pub const NONE: u32 = u32::MAX;
/// The sea's basin.
const SEA: u32 = 0;
/// The cells about a cell the water may run to.
const AROUND: [(i32, i32); 8] = [
    (1, 0),
    (-1, 0),
    (0, 1),
    (0, -1),
    (1, 1),
    (1, -1),
    (-1, 1),
    (-1, -1),
];

/// The rivers and the lakes of a year's runoff (F277): the water each cell's river carries
/// down from it, m³/s, none at sea or in a lake; the lakes' depth, m; and the year in brief.
#[derive(Clone, Debug)]
pub struct Rivers {
    pub discharge: Field<f32>,
    pub lakes: Field<f32>,
    pub summary: Summary,
}

/// The rivers' year in brief (F277): the runoff of the land about the lakes, what the lakes'
/// surfaces gain over what they lose to the air, and what reaches the sea, m³/s; the share of
/// the land, its lakes among it, whose water ends in a lake that does not spill to the sea; the
/// share of the land the lakes cover; and the most any one river brings the sea, m³/s, and
/// where.
#[derive(Clone, Copy, Debug, Default)]
pub struct Summary {
    pub runoff: f64,
    pub lakes_gain: f64,
    pub to_sea: f64,
    pub closed: f64,
    pub lakes: f64,
    pub largest: f64,
    pub mouth: [f64; 3],
}

/// A cell's height and its index, ordered by height, then index.
#[derive(Clone, Copy, PartialEq)]
struct Key(f32, u32);

impl Eq for Key {}

impl PartialOrd for Key {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Key {
    fn cmp(&self, other: &Self) -> Ordering {
        self.0.total_cmp(&other.0).then(self.1.cmp(&other.1))
    }
}

/// A basin of the flood, a pit's or the merger of two (F277): the level at which it meets
/// another, m, and the cells either side of the pass there, its own first, its level infinite
/// where it never does; whether beyond that pass lies the sea's region, the sea or the basins
/// that already spill toward it, rather than its other half, and at which of the passes it
/// joined that region; the level at which its halves
/// met and its halves, none for a pit; its catchment's runoff, kg/s; what its lake gains, less
/// what the land it covers would shed, kg/s, flooded through its halves' ranges and through its
/// own and theirs, at most none; its own cells, whose heights lie between its halves' meeting
/// and its pass, by height; and the pits it holds, as a range of the pits in the order of a walk
/// down the basins.
pub(crate) struct Basin {
    pub(crate) spill: f32,
    pub(crate) pass: (u32, u32),
    pub(crate) sea: bool,
    pub(crate) joined_sea: u32,
    pub(crate) joined: f32,
    pub(crate) halves: Option<(u32, u32)>,
    catchment: f64,
    below: f64,
    flooded: f64,
    pub(crate) own: Vec<u32>,
    pub(crate) span: (u32, u32),
}

impl Basin {
    fn new() -> Self {
        Basin {
            spill: f32::INFINITY,
            pass: (NONE, NONE),
            sea: false,
            joined_sea: NONE,
            joined: f32::NEG_INFINITY,
            halves: None,
            catchment: 0.0,
            below: 0.0,
            flooded: 0.0,
            own: Vec::new(),
            span: (0, 0),
        }
    }
}

/// A lake (F277): its basin, its level, m, its own cells flooded and the share of the next
/// that is, whether it stands full at its pass, and what it spills there, kg/s.
struct Lake {
    basin: u32,
    level: f32,
    count: usize,
    share: f64,
    full: bool,
    spills: f64,
}

/// The drainage of a year's runoff (F277, F278): each cell's receiver, the cell its water runs
/// to, itself at sea and at a dry pit; the cells in the order the flood took them, each after
/// its receiver; the lake that covers each cell, NONE where none does, and each lake's level, m;
/// the lakes' depth at each cell, m; the water each land cell carries down, kg/s; each full
/// lake's outlet, the cells either side of its pass, its own first, and what it spills over it,
/// kg/s (F278); and the year in brief.
pub struct Drainage {
    pub receiver: Vec<u32>,
    pub order: Vec<u32>,
    pub lake_of: Vec<u32>,
    pub levels: Vec<f32>,
    pub depth: Vec<f32>,
    pub carried: Vec<f64>,
    pub outlets: Vec<(u32, u32, f64)>,
    pub summary: Summary,
    /// The flood's basins, the basins at the top in the order they settle, each spilling into
    /// one settled after it; each cell's pit; and each pit's basin at the top and its place in
    /// the walk down the basins (F278).
    pub(crate) basins: Vec<Basin>,
    pub(crate) roots: Vec<u32>,
    pub(crate) label: Vec<u32>,
    pub(crate) top_of: Vec<u32>,
    pub(crate) pit_place: Vec<u32>,
}

/// The rivers and the lakes on a grid of a body of `radius` m (F277), the ground at `height`,
/// m, the sea where `sea`, under a year whose land sheds `runoff` and whose open water gains
/// `gain` over its surface, kg/m²/s at each cell, of a liquid of `density` kg/m³.
pub fn route(
    grid: CubeSphere,
    height: &[f32],
    sea: &[bool],
    water: (&[f32], &[f32]),
    (radius, density): (f64, f64),
) -> Rivers {
    let drainage = drain(grid, height, sea, water, (radius, density));
    let discharge = Field::from_fn(grid, |c, _| {
        let k = grid.index(c);
        if sea[k] || drainage.lake_of[k] != NONE {
            0.0
        } else {
            (drainage.carried[k] / density) as f32
        }
    });
    let mut lakes = Field::new(grid);
    lakes.as_mut_slice().copy_from_slice(&drainage.depth);
    Rivers {
        discharge,
        lakes,
        summary: drainage.summary,
    }
}

/// The drainage of a year's runoff on a grid (F277), as `route` reads it: the basins of the
/// flood, their lakes and the rivers between them.
pub fn drain(
    grid: CubeSphere,
    height: &[f32],
    sea: &[bool],
    (runoff, gain): (&[f32], &[f32]),
    (radius, density): (f64, f64),
) -> Drainage {
    let n = grid.cell_count();
    let area = |k: usize| grid.cell_area(grid.cell(k)) * radius * radius;
    let around = |k: usize| {
        let c = grid.cell(k);
        AROUND.map(|(di, dj)| grid.index(grid.neighbor(c, di, dj)))
    };
    // What each cell sheds as land, kg/s, and what it gains as a lake over that, at most none.
    let shed: Vec<f32> = (0..n)
        .map(|k| {
            if sea[k] {
                0.0
            } else {
                (runoff[k].max(0.0) as f64 * area(k)) as f32
            }
        })
        .collect();
    let rise: Vec<f32> = (0..n)
        .map(|k| ((gain[k] as f64 * area(k) - shed[k] as f64).min(0.0)) as f32)
        .collect();

    // The flood: from the sea and every pit, each cell in the order of its height taken by the
    // basin that reached it first, draining to the cell it was reached from; where two basins
    // meet, the lowest pass between them.
    let mut label = vec![NONE; n];
    let mut receiver = vec![NONE; n];
    let mut order: Vec<u32> = Vec::with_capacity(n);
    let mut heap = BinaryHeap::new();
    for k in (0..n).filter(|&k| sea[k]) {
        (label[k], receiver[k]) = (SEA, k as u32);
        order.push(k as u32);
    }
    for k in (0..n).filter(|&k| sea[k]) {
        for m in around(k) {
            if label[m] == NONE {
                (label[m], receiver[m]) = (SEA, k as u32);
                heap.push(Reverse(Key(height[m], m as u32)));
            }
        }
    }
    let mut pits = 0u32;
    for k in 0..n {
        if label[k] == NONE && around(k).iter().all(|&m| height[m] >= height[k]) {
            pits += 1;
            (label[k], receiver[k]) = (pits, k as u32);
            heap.push(Reverse(Key(height[k], k as u32)));
        }
    }
    let mut passes: HashMap<u64, (f32, u32, u32)> = HashMap::new();
    while let Some(Reverse(Key(h, k))) = heap.pop() {
        let k = k as usize;
        order.push(k as u32);
        let l = label[k];
        for m in around(k) {
            if label[m] == NONE {
                (label[m], receiver[m]) = (l, k as u32);
                heap.push(Reverse(Key(height[m], m as u32)));
            } else if label[m] != l {
                let level = h.max(height[m]);
                let (lo, hi, a, b) = if l < label[m] {
                    (l, label[m], k as u32, m as u32)
                } else {
                    (label[m], l, m as u32, k as u32)
                };
                let key = ((lo as u64) << 32) | hi as u64;
                let found = passes.entry(key).or_insert((level, a, b));
                if (level, a, b) < *found {
                    *found = (level, a, b);
                }
            }
        }
    }
    let mut passes: Vec<(f32, u32, u32, u32, u32)> = passes
        .into_iter()
        .map(|(key, (level, a, b))| (level, (key >> 32) as u32, key as u32, a, b))
        .collect();
    passes.sort_by(|x, y| {
        x.0.total_cmp(&y.0)
            .then((x.1, x.2, x.3, x.4).cmp(&(y.1, y.2, y.3, y.4)))
    });

    // The basins merge at their passes in the order of the passes' heights, each cell given to
    // the basin whose range holds its height as the level rises past it.
    let mut basins: Vec<Basin> = (0..=pits).map(|_| Basin::new()).collect();
    let mut parent: Vec<u32> = (0..=pits).collect();
    let mut holds: Vec<u32> = (0..=pits).collect();
    let mut next = 0;
    for &k in &order {
        let k = k as usize;
        while next < passes.len() && passes[next].0 <= height[k] {
            merge((passes[next], next), &mut basins, &mut parent, &mut holds);
            next += 1;
        }
        if label[k] == SEA {
            continue;
        }
        let root = find(&mut parent, label[k]);
        if root == SEA {
            continue;
        }
        basins[holds[root as usize] as usize].own.push(k as u32);
    }
    while next < passes.len() {
        merge((passes[next], next), &mut basins, &mut parent, &mut holds);
        next += 1;
    }

    // Each basin's catchment, and what flooding it gives: its pits' cells all shed into it;
    // its halves come before it.
    for k in 0..n {
        if label[k] != SEA && label[k] != NONE {
            basins[label[k] as usize].catchment += shed[k] as f64;
        }
    }
    let mut is_half = vec![false; basins.len()];
    for b in 1..basins.len() {
        let own: f64 = basins[b].own.iter().map(|&c| rise[c as usize] as f64).sum();
        if let Some((x, y)) = basins[b].halves {
            let (x, y) = (x as usize, y as usize);
            (is_half[x], is_half[y]) = (true, true);
            basins[b].catchment = basins[x].catchment + basins[y].catchment;
            basins[b].below = basins[x].flooded + basins[y].flooded;
        }
        basins[b].flooded = basins[b].below + own;
    }
    let mut roots: Vec<u32> = (1..basins.len() as u32)
        .filter(|&b| !is_half[b as usize])
        .collect();
    // The pits each basin holds, as a range of a walk down the basins, and the basin each pit
    // lies in at the top.
    let mut walked = 0u32;
    let mut pit_place = vec![0u32; pits as usize + 1];
    let mut top_of = vec![NONE; pits as usize + 1];
    for &r in &roots {
        let mut stack = vec![(r, false)];
        while let Some((b, done)) = stack.pop() {
            match (basins[b as usize].halves, done) {
                (Some((x, y)), false) => {
                    basins[b as usize].span.0 = walked;
                    stack.push((b, true));
                    stack.push((y, false));
                    stack.push((x, false));
                }
                (Some(_), true) => basins[b as usize].span.1 = walked,
                (None, _) => {
                    (pit_place[b as usize], top_of[b as usize]) = (walked, r);
                    basins[b as usize].span = (walked, walked + 1);
                    walked += 1;
                }
            }
        }
    }
    // The lakes, from the basins at the top down: where a basin's water, with all its halves'
    // ranges flooded, keeps it filled past their meeting, its lake stands in its own range; else
    // each half apart, a full one spilling into the other at their pass. The basins at the top
    // in the reverse of the order they joined the sea's region, each full one's spill brought
    // to the basin beyond its pass, which joined before it.
    roots.sort_by_key(|&r| Reverse((basins[r as usize].joined_sea, r)));
    let mut brought_to: Vec<Vec<(u32, f64)>> = vec![Vec::new(); basins.len()];
    let mut lakes: Vec<Lake> = Vec::new();
    for &root in &roots {
        let mut stack = vec![(root, std::mem::take(&mut brought_to[root as usize]))];
        let first = lakes.len();
        settle(
            &mut stack,
            &basins,
            (&rise, height),
            (&mut lakes, &pit_place, &label),
        );
        for lake in &lakes[first..] {
            let basin = &basins[lake.basin as usize];
            if lake.basin == root && lake.full && lake.spills > 0.0 && basin.pass.1 != NONE {
                let beyond = label[basin.pass.1 as usize];
                if beyond != SEA && beyond != NONE && top_of[beyond as usize] != NONE {
                    let to = top_of[beyond as usize] as usize;
                    brought_to[to].push((basin.pass.1, lake.spills));
                }
            }
        }
    }
    // The lakes' cells: their own flooded, and all their halves'.
    let mut depth = vec![0.0f32; n];
    let mut lake_of = vec![NONE; n];
    let (mut lake_area, mut lakes_gain) = (0.0, 0.0);
    for (i, lake) in lakes.iter().enumerate() {
        let basin = &basins[lake.basin as usize];
        let mut cells: Vec<u32> = basin.own[..lake.count].to_vec();
        let mut below: Vec<u32> = basin.halves.map_or_else(Vec::new, |(x, y)| vec![x, y]);
        while let Some(h) = below.pop() {
            let half = &basins[h as usize];
            cells.extend_from_slice(&half.own);
            if let Some((x, y)) = half.halves {
                below.extend([x, y]);
            }
        }
        for &c in &cells {
            let c = c as usize;
            depth[c] = (lake.level - height[c]).max(0.0);
            lake_of[c] = i as u32;
            lake_area += area(c);
            lakes_gain += (shed[c] + rise[c]) as f64;
        }
        if let Some(&c) = basin.own.get(lake.count).filter(|_| lake.share > 0.0) {
            let c = c as usize;
            lake_of[c] = i as u32;
            lake_area += lake.share * area(c);
            lakes_gain += shed[c] as f64 + lake.share * rise[c] as f64;
        }
    }

    // The rivers: each land cell's runoff carried down the cells it drains to, to a lake or
    // the sea; and what each full lake spills, from the cell beyond its pass down, the land's
    // runoff counted outside the lakes.
    let mut carried = vec![0.0f64; n];
    let mut mouth = (0.0f64, NONE);
    let mut to_sea = 0.0;
    for &k in order.iter().rev() {
        let k = k as usize;
        if sea[k] || lake_of[k] != NONE {
            continue;
        }
        carried[k] += shed[k] as f64;
        let r = receiver[k] as usize;
        if r == k {
            continue;
        }
        carried[r] += carried[k];
        if sea[r] {
            to_sea += carried[k];
            if carried[k] > mouth.0 {
                mouth = (carried[k], k as u32);
            }
        }
    }
    let mut outlets = Vec::new();
    for lake in lakes.iter().filter(|l| l.full && l.spills > 0.0) {
        let (own, beyond) = basins[lake.basin as usize].pass;
        if own != NONE && beyond != NONE {
            outlets.push((own, beyond, lake.spills));
        }
        let mut c = beyond;
        let mut last = NONE;
        while c != NONE && !sea[c as usize] && lake_of[c as usize] == NONE {
            carried[c as usize] += lake.spills;
            last = c;
            let r = receiver[c as usize];
            if r == c {
                break;
            }
            c = r;
        }
        if c != NONE && sea[c as usize] {
            to_sea += lake.spills;
            if last != NONE && carried[last as usize] > mouth.0 {
                mouth = (carried[last as usize], last);
            }
        }
    }

    // Which lakes' water reaches the sea: a full one at the top whose spill runs to the sea or
    // into a lake whose water does, from the lowest pass up; a half's lake holds its water.
    let mut open = vec![false; lakes.len()];
    let mut by_pass: Vec<usize> = (0..lakes.len()).collect();
    by_pass.sort_by(|&x, &y| {
        let (a, b) = (
            &basins[lakes[x].basin as usize],
            &basins[lakes[y].basin as usize],
        );
        a.spill.total_cmp(&b.spill).then(x.cmp(&y))
    });
    for i in by_pass {
        let basin = &basins[lakes[i].basin as usize];
        if !(lakes[i].full && basin.sea && basin.pass.1 != NONE) {
            continue;
        }
        let mut c = basin.pass.1 as usize;
        loop {
            if sea[c] {
                open[i] = true;
                break;
            }
            if lake_of[c] != NONE {
                open[i] = open[lake_of[c] as usize];
                break;
            }
            let r = receiver[c] as usize;
            if r == c {
                break;
            }
            c = r;
        }
    }
    // Where each land cell's water ends: the sea, or a lake whose water reaches it, or one that
    // holds it, or a dry pit.
    let mut reaches = vec![true; n];
    let (mut land, mut closed) = (0.0, 0.0);
    for &k in &order {
        let k = k as usize;
        if sea[k] {
            continue;
        }
        reaches[k] = if lake_of[k] != NONE {
            open[lake_of[k] as usize]
        } else if receiver[k] as usize == k {
            false
        } else {
            reaches[receiver[k] as usize]
        };
        let a = area(k);
        land += a;
        if !reaches[k] {
            closed += a;
        }
    }
    let volume = |kg: f64| kg / density;
    let total: f64 = (0..n)
        .filter(|&k| !sea[k] && lake_of[k] == NONE)
        .map(|k| shed[k] as f64)
        .sum();
    Drainage {
        receiver,
        order,
        levels: lakes.iter().map(|l| l.level).collect(),
        lake_of,
        depth,
        carried,
        outlets,
        basins,
        roots,
        label,
        top_of,
        pit_place,
        summary: Summary {
            runoff: volume(total),
            lakes_gain: volume(lakes_gain),
            to_sea: volume(to_sea),
            closed: if land > 0.0 { closed / land } else { 0.0 },
            lakes: if land > 0.0 { lake_area / land } else { 0.0 },
            largest: volume(mouth.0),
            mouth: if mouth.1 == NONE {
                [0.0; 3]
            } else {
                grid.center(grid.cell(mouth.1 as usize))
            },
        },
    }
}

/// The basin a pit's basin has merged into, the sea's own for those that spill to it.
fn find(parent: &mut [u32], mut x: u32) -> u32 {
    while parent[x as usize] != x {
        let up = parent[parent[x as usize] as usize];
        parent[x as usize] = up;
        x = up;
    }
    x
}

/// The basins either side of a pass at its level (F277): two merge into one above it; one
/// that meets the sea's region spills into it there.
fn merge(
    ((level, la, lb, ca, cb), at): ((f32, u32, u32, u32, u32), usize),
    basins: &mut Vec<Basin>,
    parent: &mut [u32],
    holds: &mut [u32],
) {
    let (ra, rb) = (find(parent, la), find(parent, lb));
    if ra == rb {
        return;
    }
    if ra == SEA || rb == SEA {
        let (r, own, other) = if ra == SEA {
            (rb, cb, ca)
        } else {
            (ra, ca, cb)
        };
        let b = &mut basins[holds[r as usize] as usize];
        (b.spill, b.pass, b.sea, b.joined_sea) = (level, (own, other), true, at as u32);
        parent[r as usize] = SEA;
        return;
    }
    let (na, nb) = (holds[ra as usize], holds[rb as usize]);
    (basins[na as usize].spill, basins[na as usize].pass) = (level, (ca, cb));
    (basins[nb as usize].spill, basins[nb as usize].pass) = (level, (cb, ca));
    let mut joined = Basin::new();
    (joined.joined, joined.halves) = (level, Some((na, nb)));
    basins.push(joined);
    parent[rb as usize] = ra;
    holds[ra as usize] = (basins.len() - 1) as u32;
}

/// The lakes in the basins on `stack` and below (F277), each with what it is brought at the
/// cells it enters by: where a basin's water, with all its halves' ranges flooded, keeps it
/// filled past their meeting, its lake stands in its own range, filling it from the bottom
/// while what floods each cell leaves its balance at least none, full at its pass if all of
/// them; else each half apart, a full one spilling into the other at their pass.
fn settle(
    stack: &mut Vec<(u32, Vec<(u32, f64)>)>,
    basins: &[Basin],
    (rise, height): (&[f32], &[f32]),
    (lakes, pit_place, label): (&mut Vec<Lake>, &[u32], &[u32]),
) {
    let holds_pit = |b: &Basin, cell: u32| {
        (b.span.0..b.span.1).contains(&pit_place[label[cell as usize] as usize])
    };
    while let Some((b, inflows)) = stack.pop() {
        let basin = &basins[b as usize];
        let brought: f64 = inflows.iter().map(|i| i.1).sum();
        let start = basin.catchment + basin.below + brought;
        match basin.halves {
            Some((x, y)) if start < 0.0 => {
                let (mut into_x, mut into_y): (Vec<_>, Vec<_>) = inflows
                    .into_iter()
                    .partition(|&(cell, _)| holds_pit(&basins[x as usize], cell));
                let full = |h: u32, into: &[(u32, f64)]| {
                    let half = &basins[h as usize];
                    half.catchment + half.flooded + into.iter().map(|i| i.1).sum::<f64>()
                };
                let (fx, fy) = (full(x, &into_x), full(y, &into_y));
                if fx >= 0.0 {
                    into_y.push((basins[x as usize].pass.1, fx));
                }
                if fy >= 0.0 {
                    into_x.push((basins[y as usize].pass.1, fy));
                }
                stack.push((x, into_x));
                stack.push((y, into_y));
            }
            _ => {
                let mut balance = start;
                let mut lake = Lake {
                    basin: b,
                    level: basin.spill,
                    count: basin.own.len(),
                    share: 0.0,
                    full: true,
                    spills: 0.0,
                };
                for (i, &c) in basin.own.iter().enumerate() {
                    let d = rise[c as usize] as f64;
                    if balance + d < 0.0 {
                        lake.level = height[c as usize];
                        (lake.count, lake.full) = (i, false);
                        lake.share = if d < 0.0 { balance / -d } else { 0.0 };
                        break;
                    }
                    balance += d;
                }
                if lake.full {
                    lake.spills = balance;
                    if !lake.level.is_finite() {
                        lake.level = basin
                            .own
                            .last()
                            .map_or(basin.joined, |&c| height[c as usize]);
                    }
                }
                lakes.push(lake);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::dot;

    /// A grid's heights and sea from a function of each cell's centre.
    fn terrain(grid: CubeSphere, f: impl Fn([f64; 3]) -> f64) -> (Vec<f32>, Vec<bool>) {
        let heights: Vec<f32> = (0..grid.cell_count())
            .map(|k| f(grid.center(grid.cell(k))) as f32)
            .collect();
        let sea = heights.iter().map(|&h| h < 0.0).collect();
        (heights, sea)
    }

    /// F277: a cone of an island under an even runoff sheds all of it to the sea, no lake on
    /// it, none of it closed, the rivers carrying more the lower they run.
    #[test]
    fn an_island_sheds_all_its_runoff_to_the_sea() {
        let grid = CubeSphere::new(48);
        let top = [0.0, 0.0, 1.0];
        let (heights, sea) = terrain(grid, |d| 2_000.0 - 8_000.0 * dot(d, top).acos());
        let n = grid.cell_count();
        let runoff = vec![1e-5f32; n];
        let gain = vec![-1e-5f32; n];
        let rivers = route(grid, &heights, &sea, (&runoff, &gain), (1e6, 1_000.0));
        let s = rivers.summary;
        println!(
            "runoff {:.1} m³/s, to the sea {:.1}, closed {:.4}, lakes {:.4}, the largest river {:.2}",
            s.runoff, s.to_sea, s.closed, s.lakes, s.largest
        );
        assert!(((s.to_sea - s.runoff) / s.runoff).abs() < 1e-9, "{s:?}");
        assert!(s.lakes_gain == 0.0, "{s:?}");
        assert!(s.closed == 0.0 && s.lakes == 0.0, "{s:?}");
        assert!(s.largest > 0.0 && s.largest < s.runoff, "{s:?}");
    }

    /// F277: a bowl on a high plain that drains to the sea holds a lake its climate balances:
    /// dry, the lake covers what the bowl's catchment brings it over what its surface loses,
    /// and none of it spills; wet, it fills to its rim and spills the rest to the sea.
    #[test]
    fn a_bowl_holds_the_lake_its_climate_balances() {
        let grid = CubeSphere::new(64);
        let centre = [0.0, 0.0, 1.0];
        // A plain falling gently from the bowl to a sea 60° off; the bowl 0.15 rad across,
        // 300 m deep, its rim standing 50 m over the plain about it, wide enough for the grid to
        // hold its crest.
        let shape = |d: [f64; 3]| {
            let r = dot(d, centre).clamp(-1.0, 1.0).acos();
            let plain = 1_000.0 - 1_200.0 * r;
            let bowl = if r < 0.15 {
                -300.0 * (1.0 - (r / 0.15).powi(2))
            } else {
                50.0 * (-((r - 0.15) / 0.06).powi(2)).exp()
            };
            plain + bowl
        };
        let (heights, sea) = terrain(grid, shape);
        let n = grid.cell_count();
        let radius = 1e6;
        let inside: Vec<bool> = (0..n)
            .map(|k| dot(grid.center(grid.cell(k)), centre).acos() < 0.15)
            .collect();
        let bowl_area: f64 = (0..n)
            .filter(|&k| inside[k])
            .map(|k| grid.cell_area(grid.cell(k)) * radius * radius)
            .sum();
        for (name, shed, loss) in [("dry", 1e-6f32, 1e-5f32), ("wet", 2e-5, 1e-6)] {
            let runoff = vec![shed; n];
            let gain = vec![-loss; n];
            let rivers = route(grid, &heights, &sea, (&runoff, &gain), (radius, 1_000.0));
            let lake: f64 = (0..n)
                .filter(|&k| rivers.lakes.as_slice()[k] > 0.0)
                .map(|k| grid.cell_area(grid.cell(k)) * radius * radius)
                .sum();
            let s = rivers.summary;
            println!(
                "{name}: the lake {:.3} of the bowl, closed {:.4}, to the sea {:.1} of {:.1} m³/s",
                lake / bowl_area,
                s.closed,
                s.to_sea,
                s.runoff
            );
            if name == "dry" {
                // Its catchment, the bowl and its rim's inner side, brings its runoff; the lake
                // loses its own the less: area = catchment x shed / (loss + shed).
                let catchment_share = 1.1;
                let expected = catchment_share * (shed / (loss + shed)) as f64;
                assert!(
                    lake / bowl_area > 0.5 * expected && lake / bowl_area < 2.0 * expected,
                    "{} {expected}",
                    lake / bowl_area
                );
                assert!(s.closed > 0.0 && s.to_sea < s.runoff, "{s:?}");
            } else {
                assert!(lake / bowl_area > 0.95, "{}", lake / bowl_area);
                assert!(s.closed == 0.0, "{s:?}");
                // The lake's cells no longer shed, the land's runoff counted outside it; the sea
                // takes what the lake does not lose to the air.
                let lost = lake * loss as f64 / 1_000.0;
                assert!(((s.lakes_gain + lost) / lost).abs() < 0.02, "{s:?}");
                assert!(
                    ((s.to_sea + lost - s.runoff) / s.runoff).abs() < 0.02,
                    "{s:?}"
                );
            }
        }
    }

    /// F277: two bowls side by side under a closed plain's floor, one wet and one dry: the wet
    /// one fills and spills over the low pass between them into the dry one, whose lake grows
    /// by what it is brought; the basin closed but for its rim's seaward face.
    #[test]
    fn a_full_lake_spills_into_its_neighbour() {
        let grid = CubeSphere::new(64);
        let (a, b) = (
            [0.0, 0.12f64.sin(), 0.12f64.cos()],
            [0.0, -0.12f64.sin(), 0.12f64.cos()],
        );
        let shape = |d: [f64; 3]| {
            let (ra, rb) = (dot(d, a).acos(), dot(d, b).acos());
            let well = |r: f64| -400.0 * (-(r / 0.08).powi(2)).exp();
            // A closed basin's floor, its rim far over the bowls, the sea beyond.
            let r = dot(d, [0.0, 0.0, 1.0]).acos();
            if r > 1.2 {
                return -2_000.0;
            }
            let floor = if r < 0.5 {
                500.0
            } else {
                500.0 + 4_000.0 * (r - 0.5)
            };
            floor + well(ra) + well(rb)
        };
        let (heights, sea) = terrain(grid, shape);
        let n = grid.cell_count();
        let near_a: Vec<bool> = (0..n)
            .map(|k| {
                let d = grid.center(grid.cell(k));
                dot(d, a) > dot(d, b)
            })
            .collect();
        let lake_in = |rivers: &Rivers, side: bool| -> usize {
            (0..n)
                .filter(|&k| near_a[k] == side && rivers.lakes.as_slice()[k] > 0.0)
                .count()
        };
        let loss = vec![-2e-5f32; n];
        let alone = {
            let runoff: Vec<f32> = (0..n).map(|k| if near_a[k] { 0.0 } else { 1e-7 }).collect();
            route(grid, &heights, &sea, (&runoff, &loss), (1e6, 1_000.0))
        };
        let together = {
            let runoff: Vec<f32> = (0..n)
                .map(|k| if near_a[k] { 5e-5 } else { 1e-7 })
                .collect();
            route(grid, &heights, &sea, (&runoff, &loss), (1e6, 1_000.0))
        };
        println!(
            "the dry bowl's lake {} cells alone, {} with its wet neighbour's spill; the wet one's {}",
            lake_in(&alone, false),
            lake_in(&together, false),
            lake_in(&together, true)
        );
        assert!(lake_in(&together, false) > lake_in(&alone, false));
        assert!(lake_in(&together, true) > 0);
        assert!(together.summary.closed > 0.9, "{:?}", together.summary);
    }

    /// The rivers keep the water: on rough ground of many nested basins, under a wet year, a
    /// dry one and one wet in places and dry in others, the sea takes the land's runoff and
    /// what the lakes' surfaces gain, each closed lake losing to the air all its catchment
    /// brings it; under the wet year every lake spills and all the land drains to the sea;
    /// and each land cell's water runs to the sea, a lake or a pit without turning back (the
    /// tests' tiers).
    #[test]
    fn the_rivers_keep_the_water() {
        let grid = CubeSphere::new(48);
        let n = grid.cell_count();
        for seed in 0..3u64 {
            let ground = crate::noise::GradientNoise::new(seed);
            let weather = crate::noise::GradientNoise::new(seed + 100);
            let (heights, sea) =
                terrain(grid, |d| 3_000.0 * ground.fbm(d, 4.0, 6, 2.0, 0.55) + 200.0);
            let at = |k: usize, f: f64| weather.fbm(grid.center(grid.cell(k)), f, 3, 2.0, 0.5);
            let runoff: Vec<f32> = (0..n)
                .map(|k| (1e-5 * (1.0 + 0.5 * at(k, 3.0))) as f32)
                .collect();
            for year in ["wet", "dry", "mixed"] {
                let gain: Vec<f32> = (0..n)
                    .map(|k| match year {
                        "wet" => 2e-5,
                        "dry" => -3e-5,
                        _ => (3e-5 * at(k, 7.0)) as f32,
                    })
                    .collect();
                let drainage = drain(grid, &heights, &sea, (&runoff, &gain), (1e6, 1_000.0));
                let s = drainage.summary;
                let off = (s.to_sea - s.runoff - s.lakes_gain).abs() / s.runoff;
                println!(
                    "seed {seed}, {year}: {} lakes, runoff {:.1} m³/s, lakes gain {:.1}, to the sea {:.1}, closed {:.3}, off {off:.1e}",
                    drainage.levels.len(),
                    s.runoff,
                    s.lakes_gain,
                    s.to_sea,
                    s.closed
                );
                assert!(off < 1e-6, "seed {seed}, {year}: {s:?}");
                if year == "wet" {
                    assert!(s.closed == 0.0, "seed {seed}: {s:?}");
                }
                for k in (0..n).filter(|&k| !sea[k]) {
                    let (mut c, mut steps) = (k, 0);
                    while !sea[c] && drainage.lake_of[c] == NONE {
                        let r = drainage.receiver[c] as usize;
                        if r == c {
                            break;
                        }
                        (c, steps) = (r, steps + 1);
                        assert!(
                            steps <= n,
                            "seed {seed}, {year}: cell {k}'s water turns back"
                        );
                    }
                }
            }
        }
    }
}
