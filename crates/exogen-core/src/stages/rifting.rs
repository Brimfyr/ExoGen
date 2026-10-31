// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! A rock lid's tension faults (F178), opened through the loop by what the lid's field adds
//! to its tension past the yield, each sized by what drives it.
//!
//! The lid's field (F150, F157) holds the rises' stresses, the upwelling limbs between them,
//! the relief's potential energy, the domes' among it, and the loads of the lava and of the
//! magma that stayed below (F162, F174), each bending the lid of its time. Where the field's
//! most tensile stress passes what breaks the faulted layer, the lid's normal faults open.
//! The relief alone, which the lid has borne since it formed, opens none; what the rest add
//! beyond it is shared among three drivers, each read beyond those beneath it: the rises'
//! domes beyond the bare relief, the limbs beyond the rises, the loads beyond both.
//!
//! A rise's dome and a limb concentrate the strain on a crest or along a limb, and stretch
//! the lid by their growth, the dikes and faults taking that up whatever the stress past the
//! yield (F38, F150): their extension is released as they grow, and their faults cut as deep
//! as a rise's canyon and run as long as a rise is wide or a limb long. The loads spread a
//! small strain broadly, the elastic strain past the yield, their share of the excess over
//! the lid's stiffness, taken up by graben within the faulted layer (F158). Each driver's
//! faults start where its share grew, run across the most tensile stress of the whole field
//! while their driver still stretches the lid past the yield there, and form within the step
//! its share grew in, a rise's or a limb's while the rises rose; a rift wide enough for the
//! loop's grid is carved as the mantle beneath makes its melt (F120).
//!
//! A rise's and a limb's faults grow as their driver stretches the lid (F179): its
//! extension is spent in increments of its whole, and in each the faults standing take up
//! the growth within their shadows, as far as a neighbour stands beside them, before what
//! falls outside every shadow opens new faults. A fault runs its free tips on as the field
//! leads them, its throw growing with its length, until its driver no longer stretches the
//! lid there or it meets a fault at least half as wide, which it links to; what its length
//! cannot take up deepens it, as far as the length of the system it linked into allows.

use std::collections::HashMap;

use rayon::prelude::*;

use crate::constants::*;
use crate::geom::{angle_between, local_basis};
use crate::grid::CubeSphere;
use crate::params::ParamGraph;
use crate::rng::Rng;
use crate::stages::plume::{
    FaultLayer, FaultRules, Segment, Trough, TroughKind, populate, rise_progress, trace_in,
};
use crate::stages::stress::FieldHistory;
use crate::thermal::{MeltClock, ThermalHistory};

/// What drives a share of the lid's tension past the yield.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Source {
    /// The rises' domes, beyond the bare relief's.
    Rises,
    /// The upwelling limbs between them, beyond the rises'.
    Limbs,
    /// The loads laid by then, beyond the rest.
    Loads,
    /// An icy lid stretched alike every way by its body's change of radius (F147, F233).
    Lid,
}

/// The deepest an icy lid's faults reach at a time (F147): within its elastic part, its
/// ice and a thin shell's share of the radius (`strain::brittle_bound`).
struct LidBound {
    thermal: ThermalHistory,
    shell: f64,
}

impl LidBound {
    fn at(&self, t: f64) -> f64 {
        self.thermal.at(t).elastic.min(self.shell).max(1e3)
    }
}

/// A share of the lid's tension and the faults it opens.
struct Driver {
    source: Source,
    /// What it adds to the tension past the yield, Pa, as last read, and the most it has
    /// reached at each cell, which its faults have taken up.
    added: Vec<f64>,
    high: Vec<f64>,
    rules: FaultRules,
    /// The keys its populations draw from.
    keys: (u64, u64),
    /// The extension its growth released and its faults have yet to take up, m², for a
    /// rise or a limb, and the share of its whole extension it is spent in at a time.
    pending: f64,
    increment: f64,
}

/// A fault the loop grows (F179): its trough as it stands, the ends of its segments along it
/// but the last, its driver, the draw its shape was chosen with, whether each tip, the
/// start's and the end's, may still run on, and whether it is carved by the loop and has
/// changed since the loop last heard of it.
struct Fault {
    trough: Trough,
    breaks: Vec<f64>,
    driver: usize,
    draw: u64,
    open: [bool; 2],
    carried: bool,
    changed: bool,
}

pub struct RiftProcess {
    /// The grid the field is read on and its cells' areas, m².
    coarse: CubeSphere,
    area: Vec<f64>,
    drivers: Vec<Driver>,
    /// The extension each rise takes up across its rift zone at full growth, m² (F38), and
    /// each limb along its length between the rises it joins, m² (F150).
    extension: Vec<f64>,
    limb_extension: Vec<(usize, usize, f64)>,
    /// Every fault opened so far, and their points by the coarse cell they lie in.
    faults: Vec<Fault>,
    index: HashMap<usize, Vec<(usize, [f64; 3])>>,
    /// The systems the faults linked into: each fault's parent toward its system's root, and
    /// at a root the system's summed length, m.
    parent: Vec<usize>,
    system: Vec<f64>,
    /// The width a trough the loop's grid resolves stands at, m.
    carried: f64,
    /// The rises' growth and the epochs read when last read, when, and whether yet.
    grown: Vec<f64>,
    epochs: usize,
    last: f64,
    read: bool,
    rng: Rng,
    clock: MeltClock,
    record_start: f64,
    age: f64,
    /// The populations drawn so far, for their keys, and the span of the last step a rise
    /// grew in.
    drawn: u64,
    rising: (f64, f64),
    /// The rounds the rises' and limbs' extension has been spent in so far, for their keys.
    rounds: u64,
    /// An icy lid's faults' bound through time, where the process is an icy lid's.
    lid: Option<LidBound>,
}

impl RiftProcess {
    /// A stagnant or episodic rock lid's tension faults; none on other lids.
    pub fn new(
        graph: &ParamGraph,
        grid: CubeSphere,
        thermal: &ThermalHistory,
        history: &FieldHistory,
    ) -> Option<Self> {
        if !matches!(
            graph.str("interior.tectonic_regime"),
            "stagnant_lid" | "episodic_lid"
        ) || crate::stages::io::is_io_regime(graph)
        {
            return None;
        }
        let body_radius = history.body_radius;
        let coarse = CubeSphere::new(COMPRESSION_FACE);
        let area = (0..coarse.cell_count())
            .map(|i| coarse.cell_area(coarse.cell(i)) * body_radius * body_radius)
            .collect();
        let age = graph.f64("body.age");
        let n = coarse.cell_count();
        let rules = |l_max: f64, depth_cap: f64| FaultRules {
            l_min: RIFT_MIN_LENGTH,
            l_max: l_max.max(2.0 * RIFT_MIN_LENGTH),
            depth_per_length: RIFT_DEPTH_PER_LENGTH,
            depth_cap,
            compliance: 1.0,
            layer: FaultLayer::Surface,
            cell: grid.delta() * body_radius,
            body_radius,
        };
        let driver = |source: Source, rules: FaultRules, keys: (u64, u64)| Driver {
            source,
            added: vec![0.0; n],
            high: vec![0.0; n],
            rules,
            keys,
            pending: 0.0,
            increment: f64::INFINITY,
        };
        let rises = &history.rises;
        let limbs = &history.limbs;
        let lengths = limbs.lengths(rises, body_radius);
        let mut drivers = Vec::new();
        // Where a rise's dome or a limb concentrates the strain, a fault may cut as deep as a
        // rise's canyon and run as long as the rise is wide or the limb long (F105, F150).
        let canyon = (CANYON_DEPTH_FRACTION * rises.iter().map(|r| r.height).fold(0.0, f64::max))
            .max(RIFT_MIN_DEPTH);
        if !rises.is_empty() {
            let widest = rises.iter().map(|r| r.radius).fold(0.0, f64::max);
            drivers.push(driver(
                Source::Rises,
                rules(RIFT_MAX_LENGTH_FRACTION * 2.0 * widest, canyon),
                (0x41F7_0000_0000, 70_000_000_000),
            ));
        }
        if !limbs.pairs.is_empty() {
            let longest = lengths.iter().copied().fold(0.0, f64::max);
            drivers.push(driver(
                Source::Limbs,
                rules(longest + 2.0 * limbs.width, canyon),
                (9_000_000_000, 8_000_000_000),
            ));
        }
        // The loads' strain is small and broad: graben within the faulted layer, as the lid's
        // own faults, up to a quarter of the circumference long (F158, a stated choice).
        drivers.push(driver(
            Source::Loads,
            rules(
                0.5 * std::f64::consts::PI * body_radius,
                FAULT_LAYER_RANGE.0,
            ),
            (0x10AD_0000_0000, 60_000_000_000),
        ));
        // A rise's dome stretches the lid by its growth, its dikes and faults taking that up
        // whatever the stress past the yield (F38; rifting assisted by magma, Buck 2004):
        // Valles Marineris some 30 km across Tharsis, Devana Chasma 10 to 20 km across Beta
        // Regio; a limb as much across it along its length (F150).
        let across = graph.f64("tectonics.rift_extension");
        let nominal = graph.f64("tectonics.plume_rise_height").max(1.0);
        let extension = rises
            .iter()
            .map(|r| {
                across
                    * (r.height / nominal)
                    * 2.0
                    * std::f64::consts::PI
                    * RIFT_TRANSECT_FRACTION
                    * r.radius
            })
            .collect();
        let limb_extension: Vec<(usize, usize, f64)> = limbs
            .pairs
            .iter()
            .zip(&lengths)
            .map(|(&(a, b), &length)| (a, b, across * length))
            .collect();
        let extension: Vec<f64> = extension;
        for d in &mut drivers {
            let all = match d.source {
                Source::Rises => extension.iter().sum(),
                Source::Limbs => limb_extension.iter().map(|l| l.2).sum(),
                Source::Loads | Source::Lid => f64::INFINITY,
            };
            d.increment = all / RIFT_GROWTH_INCREMENTS;
        }
        Some(Self {
            coarse,
            area,
            drivers,
            extension,
            limb_extension,
            faults: Vec::new(),
            index: HashMap::new(),
            parent: Vec::new(),
            system: Vec::new(),
            carried: RIFT_LOOP_WIDTH_CELLS * grid.delta() * body_radius,
            grown: Vec::new(),
            epochs: 0,
            last: 0.0,
            read: false,
            rng: Rng::new(graph.f64("seed.tectonics") as u64 ^ 0x41F7_E5A7),
            clock: thermal.melt_clock(4.0 * std::f64::consts::PI * body_radius * body_radius),
            record_start: graph
                .f64("cratering.record_start")
                .max(age - graph.f64("cratering.retention_age"))
                .clamp(0.0, age),
            age,
            drawn: 0,
            rising: (0.0, 0.0),
            rounds: 0,
            lid: None,
        })
    }

    /// An icy lid's faults (F147), grown as its body's change of radius stretches it alike
    /// every way (F233): one driver over the whole lid, its faults' shadows read on a grid
    /// `face` cells a side, within the lid's brittle part at each time.
    pub(crate) fn for_lid(
        graph: &ParamGraph,
        thermal: &ThermalHistory,
        rules: FaultRules,
        face: u32,
    ) -> Self {
        let body_radius = rules.body_radius;
        let coarse = CubeSphere::new(face);
        let area: Vec<f64> = (0..coarse.cell_count())
            .map(|i| coarse.cell_area(coarse.cell(i)) * body_radius * body_radius)
            .collect();
        let n = coarse.cell_count();
        let age = graph.f64("body.age");
        let carried = RIFT_LOOP_WIDTH_CELLS * rules.cell;
        Self {
            coarse,
            drivers: vec![Driver {
                source: Source::Lid,
                added: vec![1.0; n],
                high: vec![0.0; n],
                rules,
                keys: (0x5714_0000_0000, 50_000_000_000),
                pending: 0.0,
                increment: f64::INFINITY,
            }],
            clock: thermal.melt_clock(area.iter().sum()),
            area,
            extension: Vec::new(),
            limb_extension: Vec::new(),
            faults: Vec::new(),
            index: HashMap::new(),
            parent: Vec::new(),
            system: Vec::new(),
            carried,
            grown: Vec::new(),
            epochs: 0,
            last: 0.0,
            read: false,
            rng: Rng::new(graph.f64("seed.tectonics") as u64 ^ 0x5714_1D00),
            record_start: graph
                .f64("cratering.record_start")
                .max(age - graph.f64("cratering.retention_age"))
                .clamp(0.0, age),
            age,
            drawn: 0,
            rising: (0.0, 0.0),
            rounds: 0,
            lid: Some(LidBound {
                thermal: thermal.clone(),
                shell: crate::stages::strain::shell_bound(graph),
            }),
        }
    }

    /// Spends `budget`, m², of an icy lid's extension over the span of time it was taken up
    /// in, Gyr (F233): the faults standing take up the share that falls within their
    /// shadows and grow by it, the rest opens new faults outside them. What the faults took
    /// up beyond it, a new fault's being whole, the next increment's budget makes up.
    pub(crate) fn stretch(&mut self, budget: f64, span: (f64, f64)) {
        let target = self.drivers[0].pending + budget;
        if target <= 0.0 {
            self.drivers[0].pending = target;
            return;
        }
        let body_radius = self.drivers[0].rules.body_radius;
        let growth = self.area.clone();
        let whole: f64 = growth.iter().sum();
        // With the lid stretched alike every way, a fault keeps its own heading.
        let isotropic = |p: [f64; 3]| {
            let (east, north) = local_basis(p);
            (1.0, 1.0, east, east, north)
        };
        let before: f64 = self.faults.iter().map(|f| taken(&f.trough)).sum();
        self.round(
            0,
            &growth,
            whole,
            target,
            &isotropic,
            0.0,
            body_radius,
            span,
        );
        let after: f64 = self.faults.iter().map(|f| taken(&f.trough)).sum();
        self.drivers[0].pending = target - (after - before);
    }

    /// Reads the field at a time if a rise has grown by FIELD_REREAD_GROWTH or the loop has
    /// read a new epoch since it last did, grows and opens the faults it drives, and gives
    /// the troughs wide enough for the loop's grid that opened or changed since (F179).
    pub fn step(&mut self, history: &FieldHistory, t: f64) -> Vec<Trough> {
        let grown: Vec<f64> = history.rises.iter().map(|r| rise_progress(r, t)).collect();
        let moved = grown
            .iter()
            .zip(&self.grown)
            .any(|(a, b)| (a - b).abs() > FIELD_REREAD_GROWTH);
        let epochs = history.read_by(t);
        if self.read && !moved && epochs == self.epochs {
            return Vec::new();
        }
        let since = if self.read { self.last } else { 0.0 };
        // The rises grew only while they rose: their faults and their limbs' formed then.
        let before = |i: usize| self.grown.get(i).copied().unwrap_or(0.0);
        let (from, until) = history
            .rises
            .iter()
            .enumerate()
            .filter(|&(i, _)| grown[i] > before(i))
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), (_, r)| {
                (a.min(r.start_gyr), b.max(r.end_gyr))
            });
        // Their extension waiting to be spent is spent as of the last step they grew in.
        if from < until {
            self.rising = (since.max(from).min(t), t.min(until).max(since));
        }
        let rising = self.rising;
        // The extension the rises and their limbs released by their growth since.
        let by_rises: f64 = self
            .extension
            .iter()
            .enumerate()
            .map(|(i, e)| e * (grown[i] - before(i)).max(0.0))
            .sum();
        let by_limbs: f64 = self
            .limb_extension
            .iter()
            .map(|&(a, b, e)| e * (grown[a].min(grown[b]) - before(a).min(before(b))).max(0.0))
            .sum();
        for d in &mut self.drivers {
            match d.source {
                Source::Rises => d.pending += by_rises,
                Source::Limbs => d.pending += by_limbs,
                Source::Loads | Source::Lid => {}
            }
        }
        self.grown = grown;
        self.epochs = epochs;
        self.last = t;
        self.read = true;
        // What each driver adds to the excess past the yield beyond the rest beneath it: the
        // rises beyond the bare relief, the limbs beyond the rises, the loads beyond both.
        let field = history.field(t);
        let unloaded = history.field_unloaded(t);
        let domes = history.field_rises(t);
        let bare = history.field_bare(t);
        let coarse = self.coarse;
        let shares: Vec<[f64; 3]> = (0..coarse.cell_count())
            .into_par_iter()
            .map(|i| {
                let p = coarse.center(coarse.cell(i));
                let (f, u, d, b) = (
                    field.excess(p, t),
                    unloaded.excess(p, t),
                    domes.excess(p, t),
                    bare.excess(p, t),
                );
                [(d - b).max(0.0), (u - d).max(0.0), (f - u).max(0.0)]
            })
            .collect();
        for d in &mut self.drivers {
            let j = d.source as usize;
            d.added = shares.iter().map(|s| s[j]).collect();
        }
        let body_radius = history.body_radius;
        let yield_tension = field.yield_tension;
        let at = |q: [f64; 3]| field.at(q, t);
        self.drawn += 1;
        for di in 0..self.drivers.len() {
            let span = if self.drivers[di].source == Source::Loads {
                (since, t)
            } else {
                rising
            };
            let driver = &mut self.drivers[di];
            // The excess's growth past what the faults have taken up.
            let growth: Vec<f64> = driver
                .added
                .iter()
                .zip(&driver.high)
                .zip(&self.area)
                .map(|((a, h), area)| (a - h).max(0.0) * area)
                .collect();
            for (h, a) in driver.high.iter_mut().zip(&driver.added) {
                *h = h.max(*a);
            }
            let whole: f64 = growth.iter().sum();
            if whole <= 0.0 {
                continue;
            }
            let mut budget = match driver.source {
                // The loads' strain is their share of the excess over the lid's stiffness.
                Source::Loads => whole / (ROCK_YOUNG_MODULUS / (1.0 - ROCK_POISSON)),
                // A rise's or a limb's extension is spent where its share of the field grew
                // past the yield, and waits for that.
                _ => std::mem::take(&mut driver.pending),
            };
            // A rise's or a limb's in increments of its whole, each taken up by the faults
            // standing before the next opens more, so that its faults grow as it stretches
            // the lid, whatever the loop's step (F179).
            let increment = match driver.source {
                Source::Loads => f64::INFINITY,
                _ => driver.increment,
            };
            let (mut grew, mut linked, mut opened, mut rounds) = (0, 0, 0, 0);
            while budget > 0.0 {
                let chunk = budget.min(increment);
                budget -= chunk;
                let (unspent, g, l, o) = self.round(
                    di,
                    &growth,
                    whole,
                    chunk,
                    &at,
                    yield_tension,
                    body_radius,
                    span,
                );
                (grew, linked, opened, rounds) = (grew + g, linked + l, opened + o, rounds + 1);
                if unspent > 0.0 {
                    self.drivers[di].pending += unspent + budget;
                    break;
                }
            }
            if std::env::var_os("EXOGEN_DEBUG").is_some() {
                let mine = self.faults.iter().filter(|f| f.driver == di);
                let rifts = mine.clone().filter(|f| f.trough.kind == TroughKind::Rift);
                let longest = rifts.clone().map(|f| f.trough.length).fold(0.0, f64::max);
                let deepest = rifts.clone().map(|f| f.trough.depth).fold(0.0, f64::max);
                eprintln!(
                    "rifts {since:.3} to {t:.3}, {:?}: {rounds} rounds, {grew} grew ({linked} linked), {opened} opened; {} rifts ({} carried), longest {:.0} km, deepest {:.0} m",
                    self.drivers[di].source,
                    rifts.count(),
                    mine.filter(|f| f.trough.width >= self.carried).count(),
                    longest / 1e3,
                    deepest
                );
            }
        }
        // The troughs the loop's grid resolves that opened or changed in the step, carved
        // over the melt their mantle makes about the time they did.
        let now = 0.5 * (rising.0 + rising.1);
        let mut out = Vec::new();
        for f in &mut self.faults {
            if !f.changed {
                continue;
            }
            f.changed = false;
            if f.trough.width < self.carried {
                continue;
            }
            let (_, end) = carving(&self.clock, self.age, f.trough.start_gyr, now, true);
            f.trough.end_gyr = f.trough.end_gyr.max(end);
            f.carried = true;
            out.push(f.trough.clone());
        }
        out
    }

    /// The faults the loop did not carve, the graben formed before the surface began to
    /// keep a record gone with it, a fault a step long none (F150); and those it carved if
    /// it carves none (`carves` false).
    pub fn finish(self, carves: bool) -> Vec<Trough> {
        let record_start = self.record_start;
        self.faults
            .into_iter()
            .filter(|f| f.trough.points.len() >= 3)
            .filter(|f| {
                if f.carried {
                    !carves
                } else {
                    f.trough.width >= self.carried || f.trough.start_gyr >= record_start
                }
            })
            .map(|f| f.trough)
            .collect()
    }

    /// Spends `budget`, m², of a driver's extension over the growth of its share of the field
    /// in a step (F179): a rise's or a limb's faults take up the growth within their shadows
    /// and grow by it, the rest, with what they could not take up, opens new faults where it
    /// grew outside them; the loads' opens new faults wherever it grew. Gives what it could
    /// not spend, there being nowhere to open faults, and how many faults grew, linked and
    /// opened.
    #[allow(clippy::too_many_arguments)]
    fn round(
        &mut self,
        di: usize,
        growth: &[f64],
        whole: f64,
        budget: f64,
        at: &dyn Fn([f64; 3]) -> crate::stages::plume::FieldAt,
        yield_tension: f64,
        body_radius: f64,
        (begin, end): (f64, f64),
    ) -> (f64, usize, usize, usize) {
        // The loads draw once a step, the rises and limbs once a round.
        let k = if self.drivers[di].source == Source::Loads {
            self.drawn - 1
        } else {
            self.rounds += 1;
            self.rounds - 1
        };
        let coarse = self.coarse;
        let mut growth = growth.to_vec();
        let mut budget = budget;
        let (mut grew, mut linked) = (0, 0);
        if self.drivers[di].source != Source::Loads {
            let owner = self.shadows(di, body_radius);
            let mut captured: HashMap<usize, f64> = HashMap::new();
            for (i, &o) in owner.iter().enumerate() {
                if let Some(f) = o {
                    *captured.entry(f).or_default() += growth[i];
                    growth[i] = 0.0;
                }
            }
            let mut takers: Vec<(usize, f64)> = captured.into_iter().collect();
            takers.sort_by_key(|&(f, _)| f);
            let share = budget / whole;
            budget *= growth.iter().sum::<f64>() / whole;
            for (f, g) in takers {
                let (left, link) = self.grow(
                    f,
                    share * g,
                    at,
                    yield_tension,
                    body_radius,
                    (k, 0.5 * (begin + end)),
                );
                budget += left;
                grew += 1;
                linked += usize::from(link);
            }
        }
        let sum: f64 = growth.iter().sum();
        if sum <= 0.0 {
            return (budget, grew, linked, 0);
        }
        let mut cumulative = Vec::with_capacity(growth.len());
        let mut running = 0.0;
        for g in &growth {
            running += g;
            cumulative.push(running);
        }
        let n = growth.len();
        let driver = &self.drivers[di];
        let rng = &self.rng;
        let start_of = |draw: u64| -> [f64; 3] {
            let target = rng.uniform(draw, 20) * sum;
            let i = cumulative.partition_point(|&c| c <= target).min(n - 1);
            let p = coarse.center(coarse.cell(i));
            let (east, north) = local_basis(p);
            let (u, v) = (rng.uniform(draw, 22) - 0.5, rng.uniform(draw, 23) - 0.5);
            let step = coarse.delta();
            crate::geom::normalize([
                p[0] + step * (u * east[0] + v * north[0]),
                p[1] + step * (u * east[1] + v * north[1]),
                p[2] + step * (u * east[2] + v * north[2]),
            ])
        };
        let added = &driver.added;
        let unstretched = |q: [f64; 3]| added[coarse.index(coarse.locate(q).0)] <= 0.0;
        // An icy lid's new fault runs its own way until it meets a fault standing at least
        // half as wide, which takes up the stress it would have run on (F147).
        let lid = self.lid.as_ref();
        let (index, faults) = (&self.index, &self.faults);
        let meets = |q: [f64; 3], width: f64| {
            let c = coarse.locate(q).0;
            (-1..=1).any(|di| {
                (-1..=1).any(|dj| {
                    let cell = coarse.index(coarse.neighbor(c, di, dj));
                    index.get(&cell).is_some_and(|list| {
                        list.iter().any(|&(other, y)| {
                            faults[other].trough.width >= 0.5 * width
                                && angle_between(q, y) * body_radius < 0.5 * width
                        })
                    })
                })
            })
        };
        let fresh = |draw: u64, length: f64, spacing: f64, width: f64| {
            let stop = |q: [f64; 3]| unstretched(q) || (lid.is_some() && meets(q, width));
            trace_in(
                at,
                start_of(draw),
                length,
                spacing,
                body_radius,
                rng,
                draw,
                width,
                None,
                yield_tension,
                Some(&stop),
            )
        };
        let timing = |u: f64, _spent: f64, _draw: u64| begin + (end - begin) * u;
        let (clock, age) = (&self.clock, self.age);
        // A lid's fault drops its floor at once; the loop carves a wide one within the step
        // it forms (F147).
        let window = |time: f64, carried: bool| match lid {
            Some(_) => (time, time + 1e-6),
            None => carving(clock, age, begin, time, carried),
        };
        let layer_at = |t: f64| lid.map_or(0.0, |l| l.at(t));
        let keys = (
            driver.keys.0 + k * 10_000_000,
            driver.keys.1 + k * 1_000_000,
        );
        let opened = populate(
            &driver.rules,
            budget,
            rng,
            keys,
            &fresh,
            &timing,
            &window,
            &layer_at,
            None,
        );
        let count = opened.len();
        for trough in opened {
            let j = trough.id - keys.1;
            let breaks = trough
                .segments
                .iter()
                .map(|g| g.end)
                .take(trough.segments.len().saturating_sub(1))
                .collect();
            self.add(Fault {
                trough,
                breaks,
                driver: di,
                draw: keys.0 + (j - 1) * 16,
                open: [true; 2],
                carried: false,
                changed: true,
            });
        }
        (0.0, grew, linked, count)
    }

    fn add(&mut self, fault: Fault) {
        let id = self.faults.len();
        for &p in &fault.trough.points {
            let c = self.coarse.index(self.coarse.locate(p).0);
            self.index.entry(c).or_default().push((id, p));
        }
        self.parent.push(id);
        self.system.push(fault.trough.length);
        self.faults.push(fault);
    }

    /// The root of the system a fault linked into.
    fn root(&mut self, mut i: usize) -> usize {
        while self.parent[i] != i {
            self.parent[i] = self.parent[self.parent[i]];
            i = self.parent[i];
        }
        i
    }

    /// Links two faults' systems into one.
    fn link(&mut self, a: usize, b: usize) {
        let (a, b) = (self.root(a), self.root(b));
        if a != b {
            self.parent[b] = a;
            self.system[a] += self.system[b];
        }
    }

    /// The fault whose shadow each coarse cell lies in, of the driver given: within
    /// RIFT_CLUSTER_OFFSET's farthest of its widths, as far as a neighbour stands beside it
    /// (F112), the nearest by the share of its shadow where shadows meet.
    fn shadows(&self, driver: usize, body_radius: f64) -> Vec<Option<usize>> {
        let coarse = self.coarse;
        let mut owner = vec![(f64::INFINITY, None); coarse.cell_count()];
        for (fi, f) in self.faults.iter().enumerate() {
            if f.driver != driver || f.trough.points.len() < 2 {
                continue;
            }
            let reach = RIFT_CLUSTER_OFFSET.1 * f.trough.width / body_radius;
            let k = (reach / coarse.delta()).ceil() as i32;
            let stride = ((0.5 * coarse.delta() * body_radius / f.trough.spacing) as usize).max(1);
            for &p in f.trough.points.iter().step_by(stride) {
                let c = coarse.locate(p).0;
                for di in -k..=k {
                    for dj in -k..=k {
                        let nb = coarse.neighbor(c, di, dj);
                        let d = angle_between(coarse.center(nb), p) / reach;
                        let i = coarse.index(nb);
                        if d <= 1.0 && d < owner[i].0 {
                            owner[i] = (d, Some(fi));
                        }
                    }
                }
            }
        }
        owner.into_iter().map(|(_, f)| f).collect()
    }

    /// Grows a fault by the extension `budget`, m², it took up (F179): its free tips run on
    /// as the field leads them, displacement growing with length (Cowie and Scholz 1992),
    /// each until its driver no longer stretches the lid past the yield there, until the
    /// fault reaches its rules' longest, or until it comes within half its width of another
    /// fault, which it links to and runs no farther; what its length could not take up
    /// deepens it, within its rules' cap. Gives what it could not take up, and whether it
    /// linked.
    fn grow(
        &mut self,
        fi: usize,
        budget: f64,
        at: &dyn Fn([f64; 3]) -> crate::stages::plume::FieldAt,
        yield_tension: f64,
        body_radius: f64,
        (k, now): (u64, f64),
    ) -> (f64, bool) {
        let coarse = self.coarse;
        let rules = &self.drivers[self.faults[fi].driver].rules;
        let added = &self.drivers[self.faults[fi].driver].added;
        let f = &self.faults[fi];
        let before = taken(&f.trough);
        let length = f.trough.length;
        let spacing = f.trough.spacing;
        let reach = 0.5 * f.trough.width;
        let free = f.open.iter().filter(|&&o| o).count();
        let mut points = f.trough.points.clone();
        let mut breaks = f.breaks.clone();
        let mut open = f.open;
        let mut fresh: Vec<[f64; 3]> = Vec::new();
        let mut met_with: Vec<usize> = Vec::new();
        if free > 0 && length < rules.l_max && points.len() >= 2 {
            // A rift's throw grows with its length, so what it takes up as the square of its
            // length; a graben's stays within its layer and a rift at its cap stays there, so
            // theirs as their length.
            let ratio = (before + budget) / before.max(1.0);
            // An icy lid's graben deepens and widens with its length until its fault cuts
            // its layer (F233), so what it takes up grows as the square of its length until
            // then.
            let confined = match (rules.layer, &self.lid) {
                (FaultLayer::Lid(deepest), Some(lid)) => {
                    length / LID_FAULT_ASPECT >= lid.at(now).min(deepest)
                }
                _ => f.trough.kind == TroughKind::Graben || f.trough.depth >= rules.depth_cap,
            };
            let grown = if confined { ratio } else { ratio.sqrt() };
            let target = (length * grown).min(rules.l_max);
            let run = (target - length) / free as f64;
            // Each tip, the start's and the end's; the loop closes the one it runs from.
            #[allow(clippy::needless_range_loop)]
            for tip in 0..2 {
                if !open[tip] || run < spacing {
                    continue;
                }
                let n = points.len();
                let (p, q) = if tip == 0 {
                    (points[0], points[1])
                } else {
                    (points[n - 1], points[n - 2])
                };
                let heading = [p[0] - q[0], p[1] - q[1], p[2] - q[2]];
                let met = std::cell::Cell::new(None);
                let index = &self.index;
                // A tip stops at a fault at least half as wide as its own, which it links to;
                // narrower ones it cuts through (a stated choice).
                let faults = &self.faults;
                let stop = |x: [f64; 3]| {
                    if added[coarse.index(coarse.locate(x).0)] <= 0.0 {
                        return true;
                    }
                    let c = coarse.locate(x).0;
                    for di in -1..=1 {
                        for dj in -1..=1 {
                            let cell = coarse.index(coarse.neighbor(c, di, dj));
                            let near = index.get(&cell).and_then(|list| {
                                list.iter().find(|&&(other, y)| {
                                    other != fi
                                        && faults[other].trough.width >= reach
                                        && angle_between(x, y) * body_radius < reach
                                })
                            });
                            if let Some(&(other, _)) = near {
                                met.set(Some(other));
                                return true;
                            }
                        }
                    }
                    false
                };
                let draw = self.faults[fi].draw ^ (0x6E77_0000_0000 + 2 * k + tip as u64) << 8;
                let (trace, mut more) = trace_in(
                    at,
                    p,
                    run,
                    spacing,
                    body_radius,
                    &self.rng,
                    draw,
                    self.faults[fi].trough.width,
                    Some(heading),
                    yield_tension,
                    Some(&stop),
                );
                if let Some(other) = met.get() {
                    open[tip] = false;
                    met_with.push(other);
                }
                let added_points: Vec<[f64; 3]> = trace.into_iter().skip(1).collect();
                if added_points.is_empty() {
                    continue;
                }
                let run_length = added_points.len() as f64 * spacing;
                more.retain(|&b| b > 0.0 && b < run_length);
                if tip == 0 {
                    let mut ahead: Vec<[f64; 3]> = added_points.iter().rev().copied().collect();
                    let mut shifted: Vec<f64> = more.iter().rev().map(|b| run_length - b).collect();
                    shifted.extend(breaks.iter().map(|b| b + run_length));
                    breaks = shifted;
                    fresh.extend(&added_points);
                    ahead.extend(points);
                    points = ahead;
                } else {
                    let old = (points.len() - 1) as f64 * spacing;
                    breaks.extend(more.iter().map(|b| b + old));
                    fresh.extend(&added_points);
                    points.extend(added_points);
                }
            }
            if (points.len() - 1) as f64 * spacing >= rules.l_max {
                open = [false; 2];
            }
        }
        let length = (points.len().max(2) - 1) as f64 * spacing;
        let mean = 1.0 - SEGMENT_LINK_DEFICIT * (1.0 - std::f64::consts::FRAC_PI_4);
        let least = (before + budget) / (RIFT_HEAVE_PER_DEPTH * mean * length.max(1.0));
        // The system it belongs to, its length as it now stands.
        let old = self.faults[fi].trough.length;
        let root = self.root(fi);
        self.system[root] += length - old;
        for &other in &met_with {
            self.link(fi, other);
        }
        let root = self.root(fi);
        let linked = self.system[root].max(length);
        let link = !met_with.is_empty();
        let rules = &self.drivers[self.faults[fi].driver].rules;
        let mut trough = self.faults[fi].trough.clone();
        trough.points = points;
        trough.length = length;
        let bound = self.lid.as_ref().map_or(0.0, |l| l.at(now));
        shape(
            rules,
            &self.rng,
            self.faults[fi].draw,
            &mut trough,
            &breaks,
            least,
            linked,
            bound,
        );
        let left = (before + budget - taken(&trough)).max(0.0);
        // Its time: when it first broke; a graben's when it last grew, a lid's dropping its
        // floor at once (F147).
        if self.lid.is_some() {
            trough.start_gyr = trough.start_gyr.max(now.min(self.age));
            trough.end_gyr = trough.start_gyr + 1e-6;
        } else if trough.start_gyr == trough.end_gyr && !self.faults[fi].carried {
            trough.end_gyr = trough.start_gyr.max(now.min(self.age));
            trough.start_gyr = trough.end_gyr;
        }
        for p in fresh {
            let c = coarse.index(coarse.locate(p).0);
            self.index.entry(c).or_default().push((fi, p));
        }
        let f = &mut self.faults[fi];
        f.trough = trough;
        f.breaks = breaks;
        f.open = open;
        f.changed = true;
        (left, link)
    }
}

/// The times a trough is carved over (F120): those the loop carries, over the melt their
/// mantle makes about their time, from no earlier than `begin`; the rest, at once.
fn carving(clock: &MeltClock, age: f64, begin: f64, time: f64, carried: bool) -> (f64, f64) {
    if carried {
        let m = clock.at(time);
        (
            clock.time_of(m - 0.5 * RIFT_CARVING_MELT).max(begin),
            clock.time_of(m + 0.5 * RIFT_CARVING_MELT).min(age),
        )
    } else {
        (time, time)
    }
}

/// The extension a fault has taken up, m²: its heave along its length, each segment's mean
/// throw over it (F111), as the faults' generator spends it.
fn taken(trough: &Trough) -> f64 {
    let mean = 1.0 - SEGMENT_LINK_DEFICIT * (1.0 - std::f64::consts::FRAC_PI_4);
    RIFT_HEAVE_PER_DEPTH
        * mean
        * trough
            .segments
            .iter()
            .map(|g| g.throw * (g.end - g.start))
            .sum::<f64>()
}

/// A fault's trough shaped as the faults' generator shapes it (F111, F112) from its draw,
/// its trace and length and the ends of its segments, its depth at least `least`, what the
/// extension it has taken up asks of it (F179), as far as the length of the system it
/// linked into, `linked`, m, lets a fault's throw grow (Cartwright et al. 1995: linked
/// faults' throw catches up with their summed length), within its rules' cap; `bound` the
/// deepest an icy lid's faults reach at its time (F233).
#[allow(clippy::too_many_arguments)]
fn shape(
    rules: &FaultRules,
    rng: &Rng,
    draw: u64,
    trough: &mut Trough,
    breaks: &[f64],
    least: f64,
    linked: f64,
    bound: f64,
) {
    let length = trough.length;
    let scatter = (RIFT_DEPTH_SCATTER * rng.normal(draw, 2)).exp();
    let (layer_width, layer) = rules.layer.draw(rng, draw, bound, length);
    let most = (rules.depth_per_length * linked.max(length) * scatter).min(rules.depth_cap);
    let zone = (rules.depth_per_length * length * scatter)
        .clamp(RIFT_MIN_DEPTH, rules.depth_cap)
        .max(least.min(most));
    let rift = zone >= layer && matches!(rules.layer, FaultLayer::Surface);
    let width = if rift {
        (RIFT_WIDTH_PER_DEPTH * zone).max(RIFT_MIN_WIDTH)
    } else {
        layer_width
    };
    let mut bounds = vec![0.0];
    bounds.extend(breaks.iter().copied().filter(|&b| b > 0.0 && b < length));
    bounds.push(length);
    let mut polarity = if rng.uniform(draw ^ 0x5E67_0000_0000, 1) < 0.5 {
        -1.0
    } else {
        1.0
    };
    let segments: Vec<Segment> = bounds
        .windows(2)
        .enumerate()
        .map(|(j, w)| {
            let throw = if rift {
                zone
            } else {
                let cap = layer
                    .min(GRABEN_DEPTH_PER_WIDTH * width)
                    .max(RIFT_MIN_DEPTH);
                let mean_segment = length / (bounds.len() - 1).max(1) as f64;
                let plateau =
                    (SEGMENT_THROW_PER_LENGTH * rules.compliance * mean_segment * scatter)
                        .clamp(RIFT_MIN_DEPTH, cap);
                let own = (SEGMENT_THROW_SCATTER
                    * rng.normal(draw ^ 0x5E67_0000_0000, 2 * j as u64))
                .exp();
                (plateau * own).clamp(RIFT_MIN_DEPTH, cap)
            };
            if j > 0 && rng.uniform(draw ^ 0x5E67_0000_0000, 2 * j as u64 + 1) < RIFT_POLARITY_FLIP
            {
                polarity = -polarity;
            }
            Segment {
                start: w[0],
                end: w[1],
                throw,
                polarity,
            }
        })
        .collect();
    let half_graben = rng.uniform(draw ^ 0x4A1F_0000_0000, 0) < GRABEN_HALF_SHARE;
    trough.segments = if rift || half_graben {
        segments
    } else {
        segments
            .into_iter()
            .map(|g| Segment { polarity: 0.0, ..g })
            .collect()
    };
    trough.depth = trough.segments.iter().map(|g| g.throw).fold(0.0, f64::max);
    trough.width = width;
    trough.kind = if rift {
        TroughKind::Rift
    } else {
        TroughKind::Graben
    };
}

/// The faults the rises and their limbs alone open through a body's history, read every
/// hundredth of a gigayear, for the tests of what they drive.
#[cfg(test)]
pub(crate) fn rise_faults(
    graph: &ParamGraph,
    grid: CubeSphere,
    thermal: &ThermalHistory,
    rises: &[crate::stages::plume::Rise],
    limbs: &crate::stages::stress::Limbs,
) -> Vec<Trough> {
    let history = FieldHistory::new(graph, rises.to_vec(), limbs.clone());
    let Some(mut rifts) = RiftProcess::new(graph, grid, thermal, &history) else {
        return Vec::new();
    };
    let steps = (graph.f64("body.age") / 0.01).ceil() as usize;
    for s in 1..=steps {
        rifts.step(&history, (0.01 * s as f64).min(rifts.age));
    }
    rifts.finish(false)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::{angle_between, dot};
    use crate::stages::load::{flexural_parameter, lay};
    use crate::stages::stress::{FieldEpoch, Relief};

    /// A fault grown by the loop is shaped as the faults' generator shaped it (F179): given
    /// nothing more to take up, its trough comes back as drawn.
    #[test]
    fn a_fault_reshaped_as_drawn_is_unchanged() {
        let r = 3.39e6;
        let rules = FaultRules {
            l_min: RIFT_MIN_LENGTH,
            l_max: 3.0e6,
            depth_per_length: RIFT_DEPTH_PER_LENGTH,
            depth_cap: 8e3,
            compliance: 1.0,
            layer: FaultLayer::Surface,
            cell: 20e3,
            body_radius: r,
        };
        let rng = Rng::new(7);
        let fresh = |draw: u64, length: f64, spacing: f64, _width: f64| {
            let start = rng.uniform(draw, 20) * std::f64::consts::TAU;
            let n = (length / spacing).ceil() as usize;
            let points = (0..=n)
                .map(|i| {
                    let a = start + i as f64 * spacing / r;
                    [a.cos(), a.sin(), 0.0]
                })
                .collect();
            let breaks = (1..)
                .map(|j| j as f64 * 60e3)
                .take_while(|&b| b < length)
                .collect();
            (points, breaks)
        };
        let keys = (0x1234_0000_0000, 5_000_000);
        let drawn = populate(
            &rules,
            2e10,
            &rng,
            keys,
            &fresh,
            &|u, _, _| u,
            &|t, _| (t, t),
            &|_| 0.0,
            None,
        );
        assert!(drawn.len() > 20);
        let mut rifts = 0;
        for t in &drawn {
            let breaks: Vec<f64> = t.segments[..t.segments.len() - 1]
                .iter()
                .map(|g| g.end)
                .collect();
            let mut again = t.clone();
            shape(
                &rules,
                &rng,
                keys.0 + (t.id - keys.1 - 1) * 16,
                &mut again,
                &breaks,
                0.0,
                t.length,
                0.0,
            );
            assert_eq!(again.kind, t.kind, "{}", t.id);
            assert!((again.width - t.width).abs() < 1e-6 && (again.depth - t.depth).abs() < 1e-6);
            for (a, b) in again.segments.iter().zip(&t.segments) {
                assert!((a.throw - b.throw).abs() < 1e-6 && a.polarity == b.polarity);
            }
            assert!((taken(&again) - taken(t)).abs() < 1e-3);
            rifts += usize::from(t.kind == TroughKind::Rift);
        }
        assert!(rifts > 0);
    }

    /// A basin filled with lava opens graben beyond its margin (F158, F178), where its load
    /// stretches the lid outward, running round it as the arcuate graben about the maria
    /// do, as the load was laid, and none before the record began.
    #[test]
    fn a_filled_basin_opens_graben_round_its_margin() {
        let grid = CubeSphere::new(128);
        let (r, g, te) = (1_737e3, 1.62, 60e3);
        let n = grid.cell_count();
        let centre = crate::geom::normalize([1.0, 0.3, 0.2]);
        let radius = 0.25;
        let thickness: Vec<f32> = (0..n)
            .map(|k| {
                if angle_between(centre, grid.center(grid.cell(k))) < radius {
                    2_000.0
                } else {
                    0.0
                }
            })
            .collect();
        let flat = Relief {
            grid,
            height: vec![0.0; n],
        };
        let history = FieldHistory {
            rises: Vec::new(),
            limbs: Default::default(),
            body_radius: r,
            gravity: g,
            epochs: vec![
                FieldEpoch {
                    flow: None,
                    time: 0.0,
                    elastic: te,
                    relief: flat.clone(),
                    load: None,
                },
                FieldEpoch {
                    flow: None,
                    time: 3.0,
                    elastic: te,
                    relief: flat,
                    load: lay(grid, &thickness, r, g, te, 3_000.0, None).map(|l| l.stress),
                },
            ],
        };
        let (graph, _) = crate::presets::preset("Moon").unwrap().to_graph();
        let thermal = ThermalHistory::from_graph(&graph);
        // Over a few draws of the faults, as their population is a draw.
        let run = |record_start: f64| -> Vec<Trough> {
            (0..6u64)
                .flat_map(|seed| {
                    let mut rifts = RiftProcess::new(&graph, grid, &thermal, &history).unwrap();
                    rifts.record_start = record_start;
                    rifts.rng = Rng::new(seed);
                    for s in 1..=90 {
                        rifts.step(&history, 0.05 * s as f64);
                    }
                    rifts.finish(false)
                })
                .collect()
        };
        let faults = run(0.0);
        assert!(faults.len() >= 5, "{}", faults.len());
        let flexural = flexural_parameter(te, g) / r;
        // Near the margin, within a flexural parameter of it, and beyond.
        let (mut steps, mut round, mut inside) = ([0usize; 2], [0usize; 2], 0usize);
        for f in &faults {
            assert!(
                f.start_gyr > 0.0 && f.start_gyr <= 3.0 + 1e-9,
                "{}",
                f.start_gyr
            );
            let m = f.points[f.points.len() / 2];
            let d = angle_between(centre, m);
            if d > radius - flexural && d < radius + 6.0 * flexural {
                inside += 1;
            }
            for w in f.points.windows(2) {
                let band = usize::from(angle_between(centre, w[0]) > radius + flexural);
                let along = dot(centre, w[0]);
                let out = crate::geom::normalize([
                    w[0][0] * along - centre[0],
                    w[0][1] * along - centre[1],
                    w[0][2] * along - centre[2],
                ]);
                let step = crate::geom::normalize([
                    w[1][0] - w[0][0],
                    w[1][1] - w[0][1],
                    w[1][2] - w[0][2],
                ]);
                steps[band] += 1;
                if dot(out, step).abs() < 0.5 {
                    round[band] += 1;
                }
            }
        }
        eprintln!(
            "near: {} of {} round; beyond: {} of {} round",
            round[0], steps[0], round[1], steps[1]
        );
        // Nine in ten lie in the ring the load stretches, from a flexural parameter within
        // the margin to six beyond it; the rest beside one that does.
        assert!(
            inside * 10 >= faults.len() * 9,
            "{inside} of {}",
            faults.len()
        );
        assert!(round[0] * 10 > steps[0] * 7, "{round:?} of {steps:?}");
        // Formed before the record began, they are gone with it.
        assert!(run(3.1).is_empty());
    }
}
