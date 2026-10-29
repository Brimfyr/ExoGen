// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! The ridged plains' fabric (docs/parameters.md F81). The cracks narrower than a map
//! draws are the tail of the same population as the lineae (F78): placed by the same
//! tension, traced along the same field at the times they opened, and built into ridges
//! of the same shape. None is resolved: each renews the share of the cells it crosses
//! that its width covers, adds its cross-section to their relief and stains them by its
//! margins as they stand at the present. A surface renewed throughout renews itself one
//! surface's worth in each retention age; the fabric renews what the drawn lineae and the
//! chaos do not within that age, and all of it before, where what they made has been
//! overprinted into the plains. Its history is laid in slices, each one event over the
//! body, so the surface's ages come from what renewed it.

use rayon::prelude::*;

use crate::constants::*;
use crate::crater::{CellWrite, StampOutput, with_stamper};
use crate::fields::LayerSet;
use crate::geom::angle_between;
use crate::grid::CubeSphere;
use crate::params::ParamGraph;
use crate::rng::Rng;
use crate::stages::lineae::{Cracks, Linea, LineaKind, mean_axis_at};

pub const VERSION: u64 = 1;

/// A slice of the fabric's history, laid as one event over the body at its middle: the
/// narrow ridges that opened in it, as many as cover its share of the renewal.
#[derive(Clone, Debug)]
pub struct FabricSlice {
    pub start: f64,
    pub end: f64,
    /// The area its ridges cover, m².
    pub budget: f64,
    pub id: u64,
    /// The slice's place in the history, which seeds its ridges.
    index: u64,
    /// For the first slice, the time, Gyr since formation, that the surface older than
    /// the history takes: the renewal keeps no memory, so what it has not reached since
    /// the history began is as old as the process makes it on average, its start less one
    /// retention age.
    pub tail: Option<f64>,
}

impl FabricSlice {
    pub fn time_gyr(&self) -> f64 {
        0.5 * (self.start + self.end)
    }
}

/// A narrow ridge of the fabric.
#[derive(Clone, Debug)]
pub struct FabricRidge {
    pub points: Vec<[f64; 3]>,
    pub width: f64,
    pub relief: f64,
    pub time_gyr: f64,
}

/// The area of a ridge's cross-section over its height and width: the two crests of the
/// double ridge less the trough between them (`lineae::linea_profile`), each a Gaussian
/// of a sixth of the width.
fn cross_section_share() -> f64 {
    (2.0 - LINEA_TROUGH_FRACTION) * std::f64::consts::PI.sqrt() / 6.0
}

/// The mean over a crack's length of its opening, which closes toward the tips as an
/// elastic crack's does, (4 s (1 - s))^(1/2): pi / 4, the share of its width along its
/// length a ridge renews.
pub const OPENING_MEAN: f64 = std::f64::consts::FRAC_PI_4;

/// The length of a trace, m.
fn trace_length(points: &[[f64; 3]], body_radius: f64) -> f64 {
    points
        .windows(2)
        .map(|w| angle_between(w[0], w[1]))
        .sum::<f64>()
        * body_radius
}

/// The share of the surface the drawn lineae renew: each ridge and cycloid its width
/// along its length, a band, a lens, two thirds of its widest.
pub fn drawn_coverage(lineae: &[Linea], body_radius: f64) -> f64 {
    let area = 4.0 * std::f64::consts::PI * body_radius * body_radius;
    lineae
        .iter()
        .map(|l| {
            let width = match l.kind {
                LineaKind::Band => 2.0 / 3.0 * l.width,
                LineaKind::Vent | LineaKind::Fracture => l.width,
                _ => OPENING_MEAN * l.width,
            };
            width * trace_length(&l.points, body_radius)
        })
        .sum::<f64>()
        / area
}

/// The fabric's history at the epoch (F81), on a surface renewed throughout: slices a
/// FABRIC_WINDOW_SLICES-th of the retention age long on the body's own timeline, back to
/// where FABRIC_UNRENEWED of the surface is left unrenewed, which takes the age the
/// renewal gives what it has not reached. Within the retention age the fabric covers what
/// the drawn lineae and the chaos matrix, each a share of the surface in that age, leave
/// of one surface's worth; before it, all of it.
pub fn fabric_slices(graph: &ParamGraph, present: f64, drawn: f64, chaos: f64) -> Vec<FabricSlice> {
    if !crate::stages::impacts::wholesale_resurfacing(graph)
        || graph.f64("tectonics.ice_crack_density") <= 0.0
    {
        return Vec::new();
    }
    let tau = graph.f64("cratering.retention_age");
    if tau <= 0.0 || present <= 0.0 {
        return Vec::new();
    }
    let radius = graph.f64("body.radius");
    let area = 4.0 * std::f64::consts::PI * radius * radius;
    let within = (1.0 - drawn - chaos).max(0.0);
    let window = present - tau;
    let start = (present - tau * (1.0 / FABRIC_UNRENEWED).ln()).max(0.0);
    let slice = tau / FABRIC_WINDOW_SLICES as f64;
    let first = (start / slice).floor() as u64;
    let last = (present / slice).ceil() as u64;
    (first..last)
        .filter_map(|k| {
            let (s, e) = (
                (k as f64 * slice).max(start),
                ((k + 1) as f64 * slice).min(present),
            );
            if e <= s {
                return None;
            }
            // The share of the slice before the retention age, and within it.
            let before = (window.min(e) - s).max(0.0);
            let after = (e - window.max(s)).max(0.0);
            let budget = (before + within * after) / tau * area;
            (budget > 0.0).then_some(FabricSlice {
                start: s,
                end: e,
                budget,
                id: (0xFAB << 40) | k,
                index: k,
                tail: (k == first && start > 0.0).then_some((start - tau).max(0.0)),
            })
        })
        .collect()
}

impl FabricSlice {
    /// The ridges that opened in the slice: drawn one by one as the lineae are, in
    /// batches, until their widths along their lengths cover the slice's budget.
    pub fn ridges(&self, cracks: &Cracks, present: f64, cell: f64, seed: u64) -> Vec<FabricRidge> {
        let rng = Rng::new(seed ^ 0xFAB0_0000 ^ self.index.wrapping_mul(0x9E37_79B9_7F4A_7C15));
        let radius = cracks.body_radius;
        let step = FABRIC_STEP_CELLS * cell;
        let (l_lo, l_hi) = (LINEA_LENGTH_RAD.0 * radius, LINEA_LENGTH_RAD.1 * radius);
        let scale = cracks.amplitude - ICE_CRACK_STRESS;
        let thinnest = 1.0 - cracks.shell_map.thinning;
        let draw = |j: u64| -> Option<FabricRidge> {
            let time = self.start + (self.end - self.start) * rng.uniform(j, 0);
            let axis = mean_axis_at(time, present, cracks.nsr_period_gyr);
            // Where the tension above the threshold and the thin shell favour, by
            // rejection, as a drawn crack starts.
            let mut start = None;
            let freezing = cracks.freezing_at(Some(time));
            for attempt in 0..24u64 {
                let p = rng.unit_direction(j, 2 + 2 * attempt);
                let (peak, phase) = cracks.shell.peak_tension(axis, p);
                let excess = peak + freezing - ICE_CRACK_STRESS;
                if excess <= 0.0 {
                    continue;
                }
                let thin = thinnest / cracks.shell_map.factor(p).max(1e-3);
                if rng.uniform(j ^ crate::rng::ACCEPT_STREAM, attempt)
                    < excess / scale.max(excess) * thin
                {
                    start = Some((p, phase));
                    break;
                }
            }
            let (p0, phase) = start?;
            // The narrow tail of the lineae's heavy tail of widths.
            let (w_lo, w_hi) = (FABRIC_WIDTH_MIN, LINEA_WIDTH_RANGE.0);
            let k = 1.0 - LINEA_WIDTH_EXPONENT;
            let u = rng.uniform(j, 60);
            let width = (w_lo.powf(k) + u * (w_hi.powf(k) - w_lo.powf(k))).powf(1.0 / k);
            let relief = LINEA_RELIEF_SHARE * width * (0.6 + 0.8 * rng.uniform(j, 61));
            let length = l_lo * (l_hi / l_lo).powf(rng.uniform(j, 62));
            let id = (self.index << 32) ^ j ^ 0xFAB0_0000_0000_0000;
            let points = cracks.trace_at(axis, p0, phase, step, length, id, time);
            (points.len() >= 3).then_some(FabricRidge {
                points,
                width,
                relief,
                time_gyr: time,
            })
        };
        let mut out = Vec::new();
        let mut covered = 0.0;
        let mut next = 0u64;
        const BATCH: u64 = 512;
        // Never more than a few times the expected count: a field too weak to crack
        // anywhere draws nothing.
        let expected = self.budget / (FABRIC_WIDTH_MIN * l_lo);
        while covered < self.budget && (next as f64) < 4.0 * expected + BATCH as f64 {
            let batch: Vec<Option<FabricRidge>> =
                (next..next + BATCH).into_par_iter().map(draw).collect();
            next += BATCH;
            for r in batch.into_iter().flatten() {
                if covered >= self.budget {
                    break;
                }
                covered += OPENING_MEAN * r.width * trace_length(&r.points, radius);
                out.push(r);
            }
        }
        out
    }
}

/// The stain of a ridge's margins as it stands at an age, Gyr (F78): laid after the
/// crack opens, faded since to the ancient band's.
pub fn margin_stain(age: f64) -> f64 {
    LINEA_STAIN
        * (1.0 - (-age / LINEA_MARGIN_ONSET_GYR).exp())
        * (LINEA_STAIN_RESIDUAL + (1.0 - LINEA_STAIN_RESIDUAL) * (-age / SALT_FADE_GYR).exp())
}

/// The share of a cell, `cell` wide, that a strip `width` wide renews whose axis passes
/// `distance` from the cell's centre, measured across the strip. A strip narrower than
/// LINEA_RESOLVED_CELLS is drawn that wide by its share, as its relief is (F78): a
/// narrower one would fall on one cell where its axis meets a cell's centre and on two
/// halves between, and bead along its length.
pub fn strip_share(distance: f64, width: f64, cell: f64) -> f64 {
    let x = distance.abs();
    let drawn = width.max(LINEA_RESOLVED_CELLS * cell);
    width / drawn
        * ((x + 0.5 * cell).min(0.5 * drawn) - (x - 0.5 * cell).max(-0.5 * drawn)).max(0.0)
        / cell
}

/// The stain a ridge narrower than the grid resolves lays at a distance from its axis,
/// for a full stain (F78, F81): its margins, each a Gaussian of LINEA_MARGIN_SHARE of
/// its width, and the dimmer medial stripe over a third of its width, each spread over
/// at least LINEA_STAIN_MIN_CELLS by its share, so that neither beads along the line.
pub fn narrow_margin(distance: f64, width: f64, cell: f64) -> f64 {
    let spread = |true_width: f64| {
        let drawn = true_width
            .hypot(0.5 * cell)
            .max(LINEA_STAIN_MIN_CELLS * cell);
        (true_width / drawn, drawn)
    };
    let (margin_share, margin) = spread(LINEA_MARGIN_SHARE * width);
    let (stripe_share, stripe) = spread(0.3 * width);
    margin_share
        * (-(distance / margin).powi(2)).exp()
        * (1.0 - LINEA_STRIPE_DIM * stripe_share * (-(distance / stripe).powi(2)).exp())
}

/// A ridge's mark on a cell: the cell, the share of it renewed, the relief added and the
/// stain laid.
type Mark = (u32, f32, f32, f32);

/// What a slice's ridges leave on a cell: the share of it they leave unrenewed, the sum
/// of their times by the shares they renew, the sum of those shares, the relief they add
/// and the strongest stain they lay.
#[derive(Clone, Copy, Default)]
struct Acc {
    kept: f32,
    timed: f32,
    shared: f32,
    relief: f32,
    laid: f32,
}

/// One ridge's marks on the cells it crosses: the share of each it renews, the relief it
/// adds there and the stain its margins lay.
fn ridge_marks(
    ridge: &FabricRidge,
    grid: &CubeSphere,
    body_radius: f64,
    cell: f64,
    present: f64,
) -> Vec<Mark> {
    let w = ridge.width;
    // The dark material a ridge brings up goes as its cross-section, spread over margins
    // as wide as it is: the stain it lays goes as its width, full from the narrowest
    // ridge a map draws (a stated choice).
    let stain =
        margin_stain((present - ridge.time_gyr).max(0.0)) * (w / LINEA_WIDTH_RANGE.0).min(1.0);
    let section = cross_section_share() * ridge.relief;
    let margin = (LINEA_MARGIN_SHARE * w)
        .hypot(0.5 * cell)
        .max(LINEA_STAIN_MIN_CELLS * cell);
    let reach = (0.5 * w.max(LINEA_RESOLVED_CELLS * cell) + cell).max(2.0 * margin);
    with_stamper(grid, |stamper| {
        stamper.gather_corridor(grid, body_radius, &ridge.points, reach);
        stamper
            .footprint
            .iter()
            .filter_map(|&(idx, x, along)| {
                // A crack closes toward its tips as an elastic crack opens, an ellipse
                // along it: its width, its relief and its stain with it.
                let s = along / (ridge.points.len() - 1).max(1) as f64;
                let tip = (4.0 * s * (1.0 - s)).max(0.0).sqrt();
                let share = strip_share(x, w * tip, cell);
                let laid = stain * tip * narrow_margin(x, w, cell);
                (share > 0.0 || laid > 1e-4).then_some((
                    idx as u32,
                    share as f32,
                    (section * share) as f32,
                    laid as f32,
                ))
            })
            .collect()
    })
}

/// A slice of the fabric on the resolved grid (F81). Each cell takes from the ridges that
/// cross it the share of its surface they renew, their mean time for the surface it
/// becomes, their cross-sections for its relief, and the strongest stain of their
/// margins, cutting away the older stain under them.
pub fn fabric_compute(
    layers: &LayerSet,
    grid: &CubeSphere,
    cell: f64,
    present: f64,
    slice: &FabricSlice,
    cracks: &Cracks,
    seed: u64,
) -> StampOutput {
    let body_radius = cracks.body_radius;
    let ridges = slice.ridges(cracks, present, cell, seed);
    let n = grid.cell_count();
    let mut acc = vec![
        Acc {
            kept: 1.0,
            ..Acc::default()
        };
        n
    ];
    // The marks of each batch of ridges are summed in parallel over disjoint ranges of
    // cells, each range reading every mark and taking its own.
    let span = n.div_ceil(rayon::current_num_threads().max(1)).max(1);
    for chunk in ridges.chunks(2048) {
        let marks: Vec<(f32, Vec<Mark>)> = chunk
            .par_iter()
            .map(|r| {
                (
                    r.time_gyr as f32,
                    ridge_marks(r, grid, body_radius, cell, present),
                )
            })
            .collect();
        acc.par_chunks_mut(span)
            .enumerate()
            .for_each(|(range, part)| {
                let lo = range * span;
                let hi = lo + part.len();
                for (time, cells) in &marks {
                    for &(idx, share, rel, stain) in cells {
                        let k = idx as usize;
                        if k < lo || k >= hi {
                            continue;
                        }
                        let a = &mut part[k - lo];
                        a.kept *= 1.0 - share;
                        a.timed += share * time;
                        a.shared += share;
                        a.relief += rel;
                        a.laid = a.laid.max(stain);
                    }
                }
            });
    }
    let bedrock = layers.bedrock.as_slice();
    let ages = layers.surface_age.as_slice();
    let units = layers.unit_id.as_slice();
    let deep = layers.deep.as_slice();
    let melt = layers.melt.as_slice();
    let cover = layers.cover.as_slice();
    let tail = slice.tail.map(|t| t as f32);
    let writes: Vec<(CellWrite, [f32; 2])> = (0..n)
        .into_par_iter()
        .filter_map(|k| {
            let a = acc[k];
            // The surface older than the history takes the age of the renewal's tail.
            let aged = tail.is_some_and(|t| ages[k] < t);
            if a.shared <= 0.0 && a.laid <= 0.0 && !aged {
                return None;
            }
            let age = match tail {
                Some(t) if aged => t,
                _ => ages[k],
            };
            let renewed = 1.0 - a.kept;
            let mut w = CellWrite {
                idx: k as u32,
                bedrock: bedrock[k] + a.relief,
                surface: aged,
                age,
                unit: units[k],
                deep: deep[k],
                melt: melt[k],
                cover: cover[k],
            };
            if renewed > 0.0 {
                let time = a.timed / a.shared;
                w.surface = true;
                if time > age {
                    w.age = (1.0 - renewed) * age + renewed * time;
                }
                w.deep = deep[k] * (1.0 - renewed);
                w.melt = melt[k] * (1.0 - renewed);
                w.cover = cover[k] * (1.0 - renewed);
            }
            Some((w, [a.laid, renewed]))
        })
        .collect();
    let (cells, stain): (Vec<CellWrite>, Vec<[f32; 2]>) = writes.into_iter().unzip();
    StampOutput {
        deposit: None,
        cells,
        composition: None,
        stain: Some(stain),
        steep: None,
        ejecta_age: None,
        base: None,
        compositions: None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A narrow ridge renews the share of each cell its width covers, and adds the area
    /// of its cross-section, wherever its axis falls between cell centres (F81).
    #[test]
    fn a_fabric_ridge_renews_its_share_of_the_cells_it_crosses() {
        let radius = 300e3;
        let grid = CubeSphere::new(256);
        let cell = grid.delta() * radius;
        // An arc of 200 km near the middle of a face, crossing the grid obliquely.
        let a = crate::geom::normalize([1.0, -0.05, -0.2]);
        let b = crate::geom::normalize([1.0, 0.25, 0.3]);
        let n = 40;
        let points: Vec<[f64; 3]> = (0..=n)
            .map(|i| {
                let t = i as f64 / n as f64;
                crate::geom::normalize([
                    a[0] + (b[0] - a[0]) * t,
                    a[1] + (b[1] - a[1]) * t,
                    a[2] + (b[2] - a[2]) * t,
                ])
            })
            .collect();
        let length = trace_length(&points, radius);
        let ridge = FabricRidge {
            points,
            width: 800.0,
            relief: 40.0,
            time_gyr: 1.0,
        };
        let marks = ridge_marks(&ridge, &grid, radius, cell, 1.0);
        let area = |idx: u32| grid.cell_area(grid.cell(idx as usize)) * radius * radius;
        let renewed: f64 = marks.iter().map(|m| m.1 as f64 * area(m.0)).sum();
        let raised: f64 = marks.iter().map(|m| m.2 as f64 * area(m.0)).sum();
        // The ridge closes toward its tips, renewing OPENING_MEAN of its width along it.
        let expected = OPENING_MEAN * ridge.width * length;
        assert!(
            (renewed / expected - 1.0).abs() < 0.1,
            "renewed {renewed:.3e} m² of {expected:.3e}"
        );
        let section = cross_section_share() * ridge.relief * expected;
        assert!(
            (raised / section - 1.0).abs() < 0.1,
            "raised {raised:.3e} m³ of {section:.3e}"
        );
        let most = marks.iter().map(|m| m.1).fold(0.0f32, f32::max) as f64;
        assert!(most <= ridge.width / cell + 1e-6, "{most} of a cell");
        // Fresh, it has laid no stain yet.
        assert!(marks.iter().all(|m| m.3 < 1e-3));
    }

    /// Within the retention age the fabric covers what the drawn lineae and the chaos
    /// matrix leave of one surface's worth, and all of it before, back to where a
    /// hundredth of the surface is left unrenewed; a slice's ridges cover its budget
    /// with widths from the lineae's narrow tail (F81).
    #[test]
    fn the_fabric_closes_the_renewal_budget() {
        let (g, _) = crate::presets::preset("Europa").unwrap().to_graph();
        let present = g.f64("body.age");
        let tau = g.f64("cratering.retention_age");
        let radius = g.f64("body.radius");
        let area = 4.0 * std::f64::consts::PI * radius * radius;
        let slices = fabric_slices(&g, present, 0.2, 0.15);
        let total: f64 = slices.iter().map(|s| s.budget).sum();
        let history = (1.0 / FABRIC_UNRENEWED).ln();
        let expected = history - 1.0 + (1.0 - 0.2 - 0.15);
        assert!(
            (total / area - expected).abs() < 0.01,
            "{} surfaces against {expected}",
            total / area
        );
        let span = slices.last().unwrap().end - slices[0].start;
        assert!((span - history * tau).abs() < 1e-9, "{span}");
        let cracks = Cracks::new(&g).unwrap();
        let grid = CubeSphere::new(256);
        let cell = grid.delta() * radius;
        let last = slices.last().unwrap();
        let ridges = last.ridges(&cracks, present, cell, 5);
        let covered: f64 = ridges
            .iter()
            .map(|r| OPENING_MEAN * r.width * trace_length(&r.points, radius))
            .sum();
        assert!(
            covered >= last.budget && covered < 1.01 * last.budget,
            "{covered:.3e} of {:.3e}",
            last.budget
        );
        assert!(
            ridges
                .iter()
                .all(|r| (FABRIC_WIDTH_MIN..=LINEA_WIDTH_RANGE.0).contains(&r.width))
        );
        assert!(
            ridges
                .iter()
                .all(|r| r.time_gyr >= last.start && r.time_gyr <= last.end)
        );
    }
}
