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

use rayon::prelude::*;

use crate::constants::*;
use crate::geom::local_basis;
use crate::grid::CubeSphere;
use crate::params::ParamGraph;
use crate::rng::Rng;
use crate::stages::plume::{
    FaultLayer, FaultRules, Trough, TroughKind, populate, rise_progress, trace_in,
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
    /// rise or a limb.
    pending: f64,
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
    /// The rises' growth and the epochs read when last read, when, and whether yet.
    grown: Vec<f64>,
    epochs: usize,
    last: f64,
    read: bool,
    rng: Rng,
    clock: MeltClock,
    record_start: f64,
    age: f64,
    /// The populations drawn so far, for their keys.
    drawn: u64,
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
        let limb_extension = limbs
            .pairs
            .iter()
            .zip(&lengths)
            .map(|(&(a, b), &length)| (a, b, across * length))
            .collect();
        Some(Self {
            coarse,
            area,
            drivers,
            extension,
            limb_extension,
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
        })
    }

    /// The faults the field opens by a time since it was last read: none unless a rise has
    /// grown by FIELD_REREAD_GROWTH or the loop has read a new epoch since; the graben formed
    /// before the surface began to keep a record are gone with it.
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
        let rising = if from < until {
            (since.max(from).min(t), t.min(until).max(since))
        } else {
            (since, t)
        };
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
                Source::Loads => {}
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
        let k = self.drawn;
        self.drawn += 1;
        let mut faults = Vec::new();
        let body_radius = history.body_radius;
        let yield_tension = field.yield_tension;
        let rng = &self.rng;
        let (clock, age, area, record_start) =
            (&self.clock, self.age, &self.area, self.record_start);
        for driver in &mut self.drivers {
            let (begin, end) = if driver.source == Source::Loads {
                (since, t)
            } else {
                rising
            };
            let timing = |u: f64, _spent: f64, _draw: u64| begin + (end - begin) * u;
            let window = |time: f64, carried: bool| {
                if carried {
                    // Carved over the melt its mantle makes about its time (F120).
                    let m = clock.at(time);
                    (
                        clock.time_of(m - 0.5 * RIFT_CARVING_MELT).max(begin),
                        clock.time_of(m + 0.5 * RIFT_CARVING_MELT).min(age),
                    )
                } else {
                    (time, time)
                }
            };
            // The excess's growth past what the faults have taken up.
            let growth: Vec<f64> = driver
                .added
                .iter()
                .zip(&driver.high)
                .map(|(a, h)| (a - h).max(0.0))
                .collect();
            for (h, a) in driver.high.iter_mut().zip(&driver.added) {
                *h = h.max(*a);
            }
            let mut cumulative = Vec::with_capacity(growth.len());
            let mut sum = 0.0;
            for (g, a) in growth.iter().zip(area) {
                sum += g * a;
                cumulative.push(sum);
            }
            if sum <= 0.0 {
                continue;
            }
            // A rise's or a limb's extension is spent where its share of the field grew past
            // the yield, and waits for that; the loads' strain is their share over the lid's
            // stiffness.
            let budget = match driver.source {
                Source::Loads => sum / (ROCK_YOUNG_MODULUS / (1.0 - ROCK_POISSON)),
                _ => std::mem::take(&mut driver.pending),
            };
            if budget <= 0.0 {
                continue;
            }
            let n = growth.len();
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
            let fresh = |draw: u64, length: f64, spacing: f64, width: f64| {
                trace_in(
                    &|q| field.at(q, t),
                    start_of(draw),
                    length,
                    spacing,
                    body_radius,
                    rng,
                    draw,
                    width,
                    None,
                    yield_tension,
                    Some(&unstretched),
                )
            };
            let mut drawn = populate(
                &driver.rules,
                budget,
                rng,
                (
                    driver.keys.0 + k * 10_000_000,
                    driver.keys.1 + k * 1_000_000,
                ),
                &fresh,
                &timing,
                &window,
                &|_| 0.0,
                None,
            );
            // A fault that ran out of its driver's stretch within a step is no fault (F150). A
            // graben the finish draws, formed before the surface began to keep a record, is
            // gone with it; a rift the loop carves stands, degraded by what came after.
            let count = drawn.len();
            let carried = RIFT_LOOP_WIDTH_CELLS * driver.rules.cell;
            drawn.retain(|f| {
                f.points.len() >= 3 && (f.width >= carried || f.start_gyr >= record_start)
            });
            if std::env::var_os("EXOGEN_DEBUG").is_some() {
                let rifts = drawn.iter().filter(|f| f.kind == TroughKind::Rift);
                let longest = rifts.clone().map(|f| f.length).fold(0.0, f64::max);
                eprintln!(
                    "rifts {since:.3} to {t:.3}, {:?}: budget {budget:.3e} m2, {count} drawn, {} kept, {} rifts ({} carried), longest {:.0} km",
                    driver.source,
                    drawn.len(),
                    rifts.count(),
                    drawn.iter().filter(|f| f.width >= carried).count(),
                    longest / 1e3
                );
            }
            faults.extend(drawn);
        }
        faults
    }
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
    (1..=steps)
        .flat_map(|s| rifts.step(&history, (0.01 * s as f64).min(rifts.age)))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::{angle_between, dot};
    use crate::stages::load::{flexural_parameter, lay};
    use crate::stages::stress::{FieldEpoch, Relief};

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
                    time: 0.0,
                    elastic: te,
                    relief: flat.clone(),
                    load: None,
                },
                FieldEpoch {
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
                    (1..=90)
                        .flat_map(|s| rifts.step(&history, 0.05 * s as f64))
                        .collect::<Vec<Trough>>()
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
