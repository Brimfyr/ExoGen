// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! Global contraction expressed as lobate thrust scarps (PLAN.md 4.3, stagnant lid).
//!
//! The thermal history gives the cumulative radial contraction, which squeezes the lid
//! alike every way; the lid's stress field adds its own compression and tension (F150).
//! The lid holds the squeeze elastically until it passes what breaks its faulted layer in
//! compression, and its thrusts take up whatever lies beyond (F155): twice that overstress
//! over the lid's stiffness, as a share of its area, place by place, so that a thrust forms
//! where the overstress grows and a place shortens only past the most it has borne. A
//! scarp of length L and height h shortens the crust by h / tan(dip) across its length. As
//! the budget grows through the loop, scarps are sampled from a power-law length
//! distribution, the longest growing with the strain the place takes up, and stamped at
//! the step they form, so they cut what is older and are cut by what is younger. Each
//! strikes across the most compressive stress as far as the field decides it and otherwise
//! runs its own way, scalloped into lobes, and stops where the lid is no longer squeezed
//! past its strength. Each has a steep face on the hanging-wall side and a gentle back
//! slope, tapering to nothing at its ends. A thrust beneath the layered lava of the plains
//! breaks no scarp: it folds them into a wrinkle ridge.

use std::collections::HashMap;

use crate::constants::*;
use crate::crater::{CellWrite, StampOutput, Stamper};
use crate::fields::LayerSet;
use crate::geom::{angle_between, azimuth_from, cross, dot, local_basis, normalize, offset};
use crate::grid::CubeSphere;
use crate::params::ParamGraph;
use crate::rng::Rng;
use crate::stages::ridges::RidgeRecord;
use crate::stages::stress::{Compression, FieldHistory, Squeeze, shrinking};
use crate::stages::volcanism::is_lava_unit;
use crate::thermal::ThermalHistory;

pub const VERSION: u64 = 5;

/// A scarp that was emplaced, kept for the sidecar and for re-stamping its face at the
/// finishing resolution.
#[derive(Clone, Debug)]
pub struct ScarpRecord {
    /// The trace, evenly spaced points along the fault.
    pub points: Vec<[f64; 3]>,
    /// Spacing of the trace points, metres.
    pub spacing: f64,
    pub length: f64,
    pub height: f64,
    pub time_gyr: f64,
    pub id: u64,
    /// Whether the loop's grid carried it; the finish draws the rest whole (F155).
    pub carried: bool,
}

/// How the thrusts' strikes lie (F182), for the debug log: by length, in bands of latitude
/// below 30°, to 60° and beyond, the share running within 30° of north and south, and how
/// alike they run within each sixty degrees of longitude, the mean resultant length of
/// their doubled strikes, 0 at random and 1 all alike.
pub fn strike_summary(records: &[ScarpRecord]) -> String {
    let mut lines = Vec::new();
    for (lo, hi) in [(0.0, 30.0), (30.0, 60.0), (60.0, 90.0)] {
        let (mut north, mut whole) = (0.0, 0.0);
        let mut bands = [[0.0f64; 3]; 6];
        for r in records {
            let n = r.points.len();
            if n < 2 {
                continue;
            }
            let m = r.midpoint();
            let lat = m[2].clamp(-1.0, 1.0).asin().to_degrees();
            if lat.abs() < lo || lat.abs() >= hi {
                continue;
            }
            let (east, up) = crate::geom::local_basis(m);
            let d = [
                r.points[n - 1][0] - r.points[0][0],
                r.points[n - 1][1] - r.points[0][1],
                r.points[n - 1][2] - r.points[0][2],
            ];
            let theta = crate::geom::dot(d, east).atan2(crate::geom::dot(d, up));
            let w = r.length;
            whole += w;
            let off = (theta.rem_euclid(std::f64::consts::PI))
                .min(std::f64::consts::PI - theta.rem_euclid(std::f64::consts::PI));
            if off < 30f64.to_radians() {
                north += w;
            }
            let lon = m[1].atan2(m[0]).to_degrees() + 180.0;
            let b = ((lon / 60.0) as usize).min(5);
            bands[b][0] += w * (2.0 * theta).cos();
            bands[b][1] += w * (2.0 * theta).sin();
            bands[b][2] += w;
        }
        let alike: Vec<String> = bands
            .iter()
            .map(|b| {
                if b[2] > 0.0 {
                    format!("{:.2}", (b[0] * b[0] + b[1] * b[1]).sqrt() / b[2])
                } else {
                    "-".into()
                }
            })
            .collect();
        lines.push(format!(
            "lat {lo:.0}-{hi:.0}: {:.0} km, north-south {:.2} (random 0.33), alike by longitude {}",
            whole / 1e3,
            north / whole.max(1.0),
            alike.join(" ")
        ));
    }
    lines.join("; ")
}

impl ScarpRecord {
    pub fn midpoint(&self) -> [f64; 3] {
        self.points[self.points.len() / 2]
    }

    /// Angular radius of the cap holding the scarp and its back slope.
    pub fn cap_radius(&self, body_radius: f64) -> f64 {
        (0.5 * self.length + influence_width(self.height)) / body_radius
    }

    /// Strike at the midpoint, the trace heading there, as `geom::offset` measures it:
    /// counterclockwise from east.
    pub fn midpoint_heading(&self) -> f64 {
        let m = self.points.len() / 2;
        let p = self.points[m];
        let (east, north) = local_basis(p);
        azimuth_from(
            p,
            east,
            north,
            self.points[(m + 1).min(self.points.len() - 1)],
        )
    }
}

/// Height of a scarp from its length (Enterprise Rupes: 1,000 km long, 3 km high).
pub fn scarp_height(length: f64) -> f64 {
    SCARP_HEIGHT_COEFF_KM * (length / 1e3).powf(SCARP_HEIGHT_EXPONENT) * 1e3
}

/// Horizontal shortening across a thrust of this height.
pub fn scarp_shortening(height: f64) -> f64 {
    height / SCARP_DIP_DEG.to_radians().tan()
}

/// How far from the trace a scarp changes the surface.
fn influence_width(height: f64) -> f64 {
    3.0 * SCARP_BACK_SLOPE_FACTOR * height
}

/// Uplift at signed across-strike distance `d`, positive on the hanging wall: a step up
/// the face, then a back slope decaying away from the fault. The step is a tanh whose
/// central slope is the face slope, so its half-width is h / (2 tan(face slope)). `blur`
/// broadens the face, for degradation and for the coarse grid's smoothing.
pub fn scarp_profile(d: f64, height: f64, blur: f64) -> f64 {
    let face = (0.5 * height / SCARP_FACE_SLOPE_DEG.to_radians().tan()).hypot(blur);
    let step = 0.5 * (1.0 + (d / face).tanh());
    let back = SCARP_BACK_SLOPE_FACTOR * height;
    let decay = if d > 0.0 { (-d / back).exp() } else { 1.0 };
    height * step * decay
}

/// A length from the truncated power law N(>L) ∝ L^-a between the limits.
pub(crate) fn sample_length(u: f64, l_min: f64, l_max: f64) -> f64 {
    let a = SCARP_LENGTH_EXPONENT;
    let lo = l_min.powf(-a);
    let hi = l_max.powf(-a);
    (lo - u.clamp(0.0, 1.0) * (lo - hi)).powf(-1.0 / a)
}

/// A thrust's trace (F155): from `start` along `heading`, a unit vector in the tangent
/// plane, turning toward the strike the lid's compression gives as far as the field decides
/// it, and where it does not, bending its own way by `turn` over its length; lobes from a
/// few harmonics of heading wobble at the given `phases` scallop it as lobate scarps are.
/// It stops where the lid is no longer squeezed past its strength: where the uniform
/// compression's excess over the strength, `drive`, Pa, and the field's compression there
/// sum to nothing, or where the field's least compressive stress outpulls the uniform
/// compression, `squeeze`, and the lid would slip sideways instead (F157), or where `stop`
/// says it has met a thrust drawn before (F160).
#[allow(clippy::too_many_arguments)]
fn thrust_trace(
    start: [f64; 3],
    heading: [f64; 3],
    turn: f64,
    phases: [f64; SCARP_LOBE_HARMONICS],
    length: f64,
    spacing: f64,
    body_radius: f64,
    field: &dyn Fn([f64; 3]) -> Squeeze,
    drive: f64,
    squeeze: f64,
    stop: Option<&dyn Fn([f64; 3]) -> bool>,
) -> Vec<[f64; 3]> {
    let n = (length / spacing).ceil().max(2.0) as usize;
    let mut points = Vec::with_capacity(n + 1);
    let tangent = |p: [f64; 3], v: [f64; 3]| {
        let r = dot(v, p);
        normalize([v[0] - r * p[0], v[1] - r * p[1], v[2] - r * p[2]])
    };
    let mut p = start;
    let mut base = tangent(p, heading);
    for i in 0..=n {
        let here = field(p);
        if i > 0
            && (drive + here.stress <= 0.0 || here.least > squeeze || stop.is_some_and(|f| f(p)))
        {
            break;
        }
        points.push(p);
        if i == n {
            break;
        }
        let (east, north) = local_basis(p);
        // Toward the field's strike, as far as it decides it, the way the trace runs.
        let strike = if dot(here.strike, base) < 0.0 {
            [-here.strike[0], -here.strike[1], -here.strike[2]]
        } else {
            here.strike
        };
        let pull = here.decided * THRUST_FIELD_FOLLOWING;
        base = tangent(
            p,
            [
                (1.0 - pull) * base[0] + pull * strike[0],
                (1.0 - pull) * base[1] + pull * strike[1],
                (1.0 - pull) * base[2] + pull * strike[2],
            ],
        );
        // Its own bend, where the field leaves it free.
        let azimuth =
            dot(base, north).atan2(dot(base, east)) + (1.0 - here.decided) * turn / n as f64;
        base = [
            east[0] * azimuth.cos() + north[0] * azimuth.sin(),
            east[1] * azimuth.cos() + north[1] * azimuth.sin(),
            east[2] * azimuth.cos() + north[2] * azimuth.sin(),
        ];
        let s = (i as f64 + 0.5) / n as f64;
        let mut head = azimuth;
        for (k, phase) in phases.iter().enumerate() {
            let k = (k + 1) as f64;
            head += SCARP_LOBE_DEG.to_radians() / k
                * (2.0 * std::f64::consts::PI * k * s + phase).sin();
        }
        let q = offset(p, east, north, spacing / body_radius, head);
        base = tangent(q, base);
        p = q;
    }
    points
}

/// The height change one scarp makes: the profile with `blur_fine`, minus the profile
/// with `blur_coarse` when restoring detail over a coarse stamp.
pub fn scarp_compute(
    layers: &LayerSet,
    grid: &CubeSphere,
    body_radius: f64,
    rec: &ScarpRecord,
    blur_fine: f64,
    blur_coarse: Option<f64>,
    stamper: &mut Stamper,
) -> StampOutput {
    let bedrock = layers.bedrock.as_slice();
    let cell = grid.delta() * body_radius;
    let width = influence_width(rec.height) + 2.0 * cell;
    // Bands a cell wider than the spacing, so the fan on the outside of a bend leaves no
    // gaps; the nearest point still wins where bands overlap.
    let half_band = 0.5 * rec.spacing + cell;
    // Each cell takes the value from the nearest trace point along strike.
    let mut best: HashMap<usize, (f64, f64)> = HashMap::new();
    let n = rec.points.len();
    for k in 0..n {
        let p = rec.points[k];
        let prev = rec.points[k.saturating_sub(1)];
        let next = rec.points[(k + 1).min(n - 1)];
        let chord = [next[0] - prev[0], next[1] - prev[1], next[2] - prev[2]];
        let radial = dot(chord, p);
        let t = normalize([
            chord[0] - radial * p[0],
            chord[1] - radial * p[1],
            chord[2] - radial * p[2],
        ]);
        if dot(t, t) < 0.5 {
            continue;
        }
        let normal = normalize(cross(p, t));
        let along_trace = k as f64 * rec.spacing;
        let taper = (std::f64::consts::PI * along_trace / rec.length)
            .sin()
            .max(0.0)
            .sqrt();
        if taper <= 0.0 {
            continue;
        }
        stamper.gather(grid, body_radius, p, width);
        for &(idx, _, _) in &stamper.footprint {
            let d = grid.center(grid.cell(idx));
            let along = dot(d, t) * body_radius;
            if along.abs() > half_band {
                continue;
            }
            let across = dot(d, normal) * body_radius;
            let mut u = scarp_profile(across, rec.height, blur_fine);
            if let Some(c) = blur_coarse {
                u -= scarp_profile(across, rec.height, c);
            }
            let u = taper * u;
            let entry = best.entry(idx).or_insert((f64::INFINITY, 0.0));
            if along.abs() < entry.0 {
                *entry = (along.abs(), u);
            }
        }
    }
    let mut cells: Vec<CellWrite> = best
        .into_iter()
        .map(|(idx, (_, u))| CellWrite {
            idx: idx as u32,
            bedrock: bedrock[idx] + u as f32,
            surface: false,
            age: 0.0,
            unit: 0,
            deep: 0.0,
            melt: 0.0,
            cover: 0.0,
        })
        .collect();
    cells.sort_by_key(|c| c.idx);
    StampOutput {
        cells,
        composition: None,
        stain: None,
        steep: None,
        ejecta_age: None,
        compositions: None,
    }
}

/// A thrust or a ridge drawn before (F160): its trace, its trace's spacing, its length and
/// its reach, the width over which it changes the ground, m.
pub(crate) struct Drawn {
    pub points: Vec<[f64; 3]>,
    pub spacing: f64,
    pub length: f64,
    pub reach: f64,
}

/// A drawn trace's point as the systems index it: where, the half-width it holds off others,
/// its trace's length, and which trace and which point of it.
type Indexed = ([f64; 3], f64, f64, u32, u32);

/// Where the thrusts or the ridges drawn so far run (F160), as the faults' generator keeps
/// them (F112, F147): their traces, for the half that run beside one drawn before near where
/// they start, en-echelon along a span of it, as thrusts gather in systems and the maria's
/// ridges in rings broken by gaps; and their points on a coarse grid, for the arrest of the
/// rest, which stop where they meet one drawn before at least as long, whose own slip took up
/// the squeeze they would have run on.
pub(crate) struct Systems {
    coarse: CubeSphere,
    body_radius: f64,
    /// By coarse cell: each point, the half-width it holds off others, its trace's length,
    /// and which trace and which point of it.
    index: HashMap<usize, Vec<Indexed>>,
    pub drawn: Vec<Drawn>,
}

impl Systems {
    /// On a coarse grid with cells as wide as `widest`, m, the most a thrust and another can
    /// stand apart and still meet.
    pub fn new(widest: f64, body_radius: f64) -> Self {
        let face = ((std::f64::consts::TAU * body_radius / (4.0 * widest.max(1.0))).round() as u32)
            .clamp(8, 256);
        Self {
            coarse: CubeSphere::new(face),
            body_radius,
            index: HashMap::new(),
            drawn: Vec::new(),
        }
    }

    /// Note one drawn.
    pub fn add(&mut self, points: Vec<[f64; 3]>, spacing: f64, length: f64, reach: f64) {
        let half = SYSTEM_ARREST_SHARE * reach;
        let which = self.drawn.len() as u32;
        for (k, &p) in points.iter().enumerate() {
            let c = self.coarse.index(self.coarse.locate(p).0);
            self.index
                .entry(c)
                .or_default()
                .push((p, half, length, which, k as u32));
        }
        self.drawn.push(Drawn {
            points,
            spacing,
            length,
            reach,
        });
    }

    /// Whether a point of one `length` long, holding off others by `half`, meets one drawn
    /// before at least as long.
    pub fn meets(&self, p: [f64; 3], half: f64, length: f64) -> bool {
        let (c, _, _) = self.coarse.locate(p);
        for di in -1..=1 {
            for dj in -1..=1 {
                let k = self.coarse.index(self.coarse.neighbor(c, di, dj));
                if let Some(points) = self.index.get(&k)
                    && points.iter().any(|&(q, h, l, _, _)| {
                        l >= length && angle_between(p, q) * self.body_radius < half + h
                    })
                {
                    return true;
                }
            }
        }
        false
    }

    /// A trace beside one drawn before near `start`, or none: whether to run beside one
    /// (`join` below RIFT_CLUSTER_PROBABILITY), which of those whose traces pass within a
    /// coarse cell or so of the start, by length (`pick`), on which side (`side`) and how
    /// far, SYSTEM_SPACING of its reach (`gap`), along a span as long as `length` about its
    /// point nearest the start; walked at the host's spacing, which it returns with the
    /// trace.
    pub fn beside(
        &self,
        start: [f64; 3],
        length: f64,
        join: f64,
        pick: f64,
        side: f64,
        gap: f64,
    ) -> Option<(Vec<[f64; 3]>, f64)> {
        if join >= RIFT_CLUSTER_PROBABILITY {
            return None;
        }
        // Each nearby trace, with its point nearest the start.
        let (c, _, _) = self.coarse.locate(start);
        let mut near: Vec<(u32, u32, f64)> = Vec::new();
        for di in -1..=1 {
            for dj in -1..=1 {
                let k = self.coarse.index(self.coarse.neighbor(c, di, dj));
                if let Some(points) = self.index.get(&k) {
                    for &(q, _, _, which, at) in points {
                        let d = angle_between(start, q);
                        match near.iter_mut().find(|e| e.0 == which) {
                            Some(e) if d < e.2 => *e = (which, at, d),
                            Some(_) => {}
                            None => near.push((which, at, d)),
                        }
                    }
                }
            }
        }
        if near.is_empty() {
            return None;
        }
        near.sort_by_key(|e| e.0);
        let total: f64 = near.iter().map(|e| self.drawn[e.0 as usize].length).sum();
        let mut target = pick * total;
        let mut chosen = near[near.len() - 1];
        for e in &near {
            let l = self.drawn[e.0 as usize].length;
            if target < l {
                chosen = *e;
                break;
            }
            target -= l;
        }
        let host = &self.drawn[chosen.0 as usize];
        let n = host.points.len();
        if n < 3 {
            return None;
        }
        let span = ((length / host.spacing).round() as usize).clamp(2, n);
        let i0 = (chosen.1 as usize).saturating_sub(span / 2).min(n - span);
        let (lo, hi) = SYSTEM_SPACING;
        let distance = host.reach * (lo + (hi - lo) * gap);
        let sign = if side < 0.5 { -1.0 } else { 1.0 };
        let trace = crate::stages::plume::offset_trace(
            &host.points[i0..i0 + span],
            sign * distance,
            self.body_radius,
            host.spacing,
        );
        (trace.len() >= 3).then_some((trace, host.spacing))
    }
}

/// Contraction scarps forming through the coarse loop.
pub struct ContractionProcess {
    grid: CubeSphere,
    body_radius: f64,
    thermal: ThermalHistory,
    age: f64,
    /// When the surface began to keep a record, Gyr: the thrusts before it are gone with
    /// it, as the lid's faults are (F147), though the lid failed and bore them.
    record_start: f64,
    rng: Rng,
    stamper: Stamper,
    /// Shortening area still to be taken up, carried between steps.
    carry: f64,
    l_min: f64,
    stamp_min: f64,
    /// The lid's compression the thrusts read (F154), its stiffness alike every way, Pa,
    /// and each coarse cell's overstress at its highest so far, Pa: the lid there has
    /// failed up to it, and only what passes it shortens it anew (F155).
    compression: Compression,
    stiffness: f64,
    high: Vec<f64>,
    pub records: Vec<ScarpRecord>,
    /// The thrusts beneath the plains, which fold them into ridges (F155).
    pub ridges: Vec<RidgeRecord>,
    /// The thrusts drawn so far, for their systems (F160).
    systems: Systems,
}

/// The longest thrust a strain makes: deeper faulting under higher stress and the linking
/// of segments into longer systems. Mercury's 0.4 per cent beyond its lid's elastic strain
/// gives thousand-kilometre scarps; a hundredth of that, kilometres.
fn longest_thrust(strain: f64, body_radius: f64) -> f64 {
    (SCARP_REFERENCE_MAX_LENGTH
        * (strain / SCARP_REFERENCE_STRAIN).powf(SCARP_LENGTH_STRAIN_EXPONENT))
    .min(SCARP_MAX_LENGTH_FRACTION * body_radius)
}

impl ContractionProcess {
    /// None for regimes that do not express contraction as scarps: a whole rock lid, not
    /// Io's, whose surface its eruptions lay afresh. `history` is the lid's field through
    /// the loop (F150, F154, F157).
    pub fn new(
        graph: &ParamGraph,
        grid: CubeSphere,
        thermal: &ThermalHistory,
        history: &FieldHistory,
    ) -> Option<Self> {
        if graph.str("interior.tectonic_regime") != "stagnant_lid"
            || crate::stages::io::is_io_regime(graph)
        {
            return None;
        }
        let body_radius = graph.f64("body.radius");
        let age = graph.f64("body.age");
        let cell = 2.0 * std::f64::consts::PI * body_radius / (4.0 * grid.face_size() as f64);
        let stiffness = ROCK_YOUNG_MODULUS / (1.0 - ROCK_POISSON);
        let compression = Compression::new(history, COMPRESSION_FACE);
        let cells = compression.grid.cell_count();
        Some(Self {
            grid,
            body_radius,
            thermal: thermal.clone(),
            age,
            record_start: graph
                .f64("cratering.record_start")
                .max(age - graph.f64("cratering.retention_age"))
                .clamp(0.0, age),
            rng: Rng::new(graph.f64("seed.tectonics") as u64 ^ 0x5CA2_9F11),
            stamper: Stamper::new(&grid),
            carry: 0.0,
            l_min: SCARP_MIN_LENGTH,
            stamp_min: SCARP_LOOP_MIN_CELLS * cell,
            compression,
            stiffness,
            high: vec![0.0; cells],
            records: Vec::new(),
            ridges: Vec::new(),
            systems: Systems::new(SYSTEM_INDEX_WIDTH, body_radius),
        })
    }

    /// The uniform compression's excess over the lid's strength at a time, Pa: what the
    /// body's shrinking puts on the lid alike every way, over its stiffness.
    fn drive(&self, t: f64) -> f64 {
        shrinking(&self.thermal, self.body_radius, t) - self.compression.yield_compression
    }

    /// Scarps formed between two times, in order of formation, not yet emplaced. The
    /// budget is spent on the whole population; those the loop grid can carry are marked
    /// so, and the rest are left to the finish.
    pub fn sample(&mut self, t0: f64, t1: f64, history: &FieldHistory) -> Vec<ScarpRecord> {
        let mut out = Vec::new();
        for (bin, b0, b1) in crate::stages::time_bins(self.age, t0, t1) {
            self.sample_bin(bin, b0, b1, history, &mut out);
        }
        out.sort_by(|a, b| a.time_gyr.total_cmp(&b.time_gyr).then(a.id.cmp(&b.id)));
        out
    }

    /// The scarps of one time bin, drawn with the bin as the random index.
    fn sample_bin(
        &mut self,
        bin: u64,
        t0: f64,
        t1: f64,
        history: &FieldHistory,
        out: &mut Vec<ScarpRecord>,
    ) {
        let step_index = bin;
        // Where the lid's overstress grew in the bin, and the shortening it asks.
        self.compression.at(t1, history);
        let drive = self.drive(t1);
        let shrink = drive + self.compression.yield_compression;
        // The shrinking still to come, which every place's strain will take up as well.
        let to_come = (self.drive(self.age) - drive).max(0.0);
        let mut weights = Vec::with_capacity(self.high.len());
        let mut sum = 0.0;
        let mut budget = self.carry;
        for (i, q) in self.compression.squeeze.iter().enumerate() {
            // A thrust where both horizontal stresses press harder than the overburden.
            let over = if q.least <= shrink {
                (drive + q.stress).max(0.0)
            } else {
                0.0
            };
            let grown = (over - self.high[i]).max(0.0);
            self.high[i] = self.high[i].max(over);
            let w = grown * self.compression.area[i];
            budget += 2.0 * w / self.stiffness;
            sum += w;
            weights.push(sum);
        }
        if t1 <= self.record_start {
            // The lid failed and bore them, but the surface keeps nothing of them.
            self.carry = 0.0;
            return;
        }
        if sum <= 0.0 {
            self.carry = budget.min(0.0);
            return;
        }
        let grid = self.compression.grid;
        let half_cell = 0.5 * grid.delta();
        let mut j = 0u64;
        while budget > 0.0 && j < 200_000 {
            let draw = 16 + j * 16;
            j += 1;
            // Where it starts: as the overstress grew, anywhere in its cell.
            let u = self.rng.uniform(step_index, draw + 1) * sum;
            let i = weights.partition_point(|&w| w < u).min(weights.len() - 1);
            let center = grid.center(grid.cell(i));
            let (east, north) = local_basis(center);
            let start = offset(
                center,
                east,
                north,
                half_cell * self.rng.uniform(step_index, draw + 9).sqrt(),
                2.0 * std::f64::consts::PI * self.rng.uniform(step_index, draw + 2),
            );
            let reach = longest_thrust(
                2.0 * (self.high[i] + to_come) / self.stiffness,
                self.body_radius,
            );
            if reach <= self.l_min {
                // Its strain goes into faults too short to draw, as the shortest drawn.
                budget -= self.l_min * scarp_shortening(scarp_height(self.l_min));
                continue;
            }
            let length = sample_length(self.rng.uniform(step_index, draw), self.l_min, reach);
            let time_gyr = t0 + (t1 - t0) * self.rng.uniform(step_index, draw + 5);
            // Half beside a thrust drawn before near where it starts, en-echelon along a
            // span of it, so that systems grow where the lid is squeezed (F160).
            let beside = self.systems.beside(
                start,
                length,
                self.rng.uniform(step_index, draw + 10),
                self.rng.uniform(step_index, draw + 11),
                self.rng.uniform(step_index, draw + 13),
                self.rng.uniform(step_index, draw + 14),
            );
            let (points, spacing) = match beside {
                Some(found) => found,
                None => {
                    let spacing = (0.5 * self.stamp_min / SCARP_LOOP_MIN_CELLS).min(length / 8.0);
                    // Its strike: across the most compressive stress as far as the field
                    // decides it, its own way otherwise.
                    let here = self.compression.reading(start);
                    let phi = 2.0 * std::f64::consts::PI * self.rng.uniform(step_index, draw + 3);
                    let free = [
                        east[0] * phi.cos() + north[0] * phi.sin(),
                        east[1] * phi.cos() + north[1] * phi.sin(),
                        east[2] * phi.cos() + north[2] * phi.sin(),
                    ];
                    let strike = if dot(here.strike, free) < 0.0 {
                        [-here.strike[0], -here.strike[1], -here.strike[2]]
                    } else {
                        here.strike
                    };
                    let a = here.decided;
                    let heading = [
                        a * strike[0] + (1.0 - a) * free[0],
                        a * strike[1] + (1.0 - a) * free[1],
                        a * strike[2] + (1.0 - a) * free[2],
                    ];
                    let turn = SCARP_TURN_DEG.to_radians()
                        * (2.0 * self.rng.uniform(step_index, draw + 4) - 1.0);
                    let mut phases = [0.0; SCARP_LOBE_HARMONICS];
                    for (k, phase) in phases.iter_mut().enumerate() {
                        *phase = 2.0
                            * std::f64::consts::PI
                            * self.rng.uniform(step_index, draw + 6 + k as u64);
                    }
                    // It stops where it meets a thrust drawn before at least as long.
                    let half = SYSTEM_ARREST_SHARE * influence_width(scarp_height(length));
                    let systems = &self.systems;
                    let stop = |p: [f64; 3]| systems.meets(p, half, length);
                    let compression = &self.compression;
                    let points = thrust_trace(
                        start,
                        heading,
                        turn,
                        phases,
                        length,
                        spacing,
                        self.body_radius,
                        &|p| compression.reading(p),
                        drive,
                        shrink,
                        Some(&stop),
                    );
                    (points, spacing)
                }
            };
            // What it takes up, as far as it ran.
            let actual = (points.len() - 1) as f64 * spacing;
            let height = scarp_height(actual.max(1.0));
            budget -= actual * scarp_shortening(height);
            if actual < self.l_min || points.len() < 3 {
                continue;
            }
            self.systems
                .add(points.clone(), spacing, actual, influence_width(height));
            out.push(ScarpRecord {
                points,
                spacing,
                length: actual,
                height,
                time_gyr,
                id: (step_index << 32) | (j - 1),
                carried: actual >= self.stamp_min,
            });
        }
        self.carry = budget.min(0.0);
    }

    /// Cut one scarp into the surface as it is now, if the loop's grid carries it, and
    /// record it; beneath the plains, record the ridge it folds them into instead.
    pub fn apply(&mut self, layers: &mut LayerSet, rec: ScarpRecord) {
        let units = layers.unit_id.as_slice();
        let on_lava = rec
            .points
            .iter()
            .filter(|p| is_lava_unit(units[self.grid.index(self.grid.locate(**p).0)]))
            .count();
        if 2 * on_lava > rec.points.len() {
            // A thrust beneath layered lava folds it into a ridge as the plains' own are,
            // its slip taken up in the fold (F155).
            let height =
                crate::stages::ridges::ridge_height(rec.length, self.rng.uniform(rec.id, 0xA1));
            let crest = RIDGE_CREST_HEIGHT * height / RIDGE_HEIGHT
                * (0.6 + 0.8 * self.rng.uniform(rec.id, 0xA2));
            self.ridges.push(RidgeRecord {
                points: rec.points,
                spacing: rec.spacing,
                length: rec.length,
                height,
                crest,
                time_gyr: rec.time_gyr,
                id: rec.id ^ 0x7A12_0000_0000,
            });
            return;
        }
        if rec.carried {
            let out = scarp_compute(
                layers,
                &self.grid,
                self.body_radius,
                &rec,
                0.0,
                None,
                &mut self.stamper,
            );
            out.apply(layers);
        }
        self.records.push(rec);
    }

    /// Sample and emplace a step's scarps.
    pub fn step(&mut self, layers: &mut LayerSet, t0: f64, t1: f64, history: &FieldHistory) {
        for rec in self.sample(t0, t1, history) {
            self.apply(layers, rec);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_profile_steps_up_the_face_and_decays_down_the_back_slope() {
        let h = 1_000.0;
        assert!(scarp_profile(-50_000.0, h, 0.0).abs() < 1.0);
        let face = scarp_profile(2_000.0, h, 0.0);
        assert!(face > 0.5 * h, "{face}");
        let far = scarp_profile(100_000.0, h, 0.0);
        assert!(far < 0.1 * h, "{far}");
        // The central slope is the face slope.
        let face = 0.5 * h / SCARP_FACE_SLOPE_DEG.to_radians().tan();
        let slope = (scarp_profile(1.0, h, 0.0) - scarp_profile(-1.0, h, 0.0)) / 2.0;
        assert!(
            (slope - SCARP_FACE_SLOPE_DEG.to_radians().tan()).abs() < 0.02,
            "{slope}"
        );
        assert!(face > 1_000.0 && face < 3_000.0, "{face}");
        // Blur broadens the face without changing its height far from it.
        let blurred = scarp_profile(0.0, h, 10_000.0);
        assert!(blurred > 0.4 * h && blurred < 0.6 * h);
        assert!((scarp_shortening(1_000.0) - 2_144.5).abs() < 5.0);
        assert!((scarp_height(1.0e6) - 3_000.0).abs() < 50.0);
    }

    /// A scarp stamped on a flat grid: nothing on the footwall, most of the height just
    /// behind the face, little far behind it, and little near the ends of the trace.
    #[test]
    fn a_stamped_scarp_steps_up_behind_its_face_and_tapers_at_its_ends() {
        let grid = CubeSphere::new(256);
        let body_radius = 300e3;
        let mut layers = LayerSet::new(grid);
        let length = 120e3;
        let height = scarp_height(length);
        let spacing = 0.5 * grid.delta() * body_radius;
        let start = normalize([1.0, 0.2, 0.1]);
        let (east, north) = local_basis(start);
        let heading = [
            east[0] * 0.3f64.cos() + north[0] * 0.3f64.sin(),
            east[1] * 0.3f64.cos() + north[1] * 0.3f64.sin(),
            east[2] * 0.3f64.cos() + north[2] * 0.3f64.sin(),
        ];
        let points = thrust_trace(
            start,
            heading,
            0.2,
            [0.0; SCARP_LOBE_HARMONICS],
            length,
            spacing,
            body_radius,
            &|_| Squeeze::default(),
            1e9,
            1e9,
            None,
        );
        let rec = ScarpRecord {
            points,
            spacing,
            length,
            height,
            time_gyr: 1.0,
            id: 1,
            carried: true,
        };
        let mut stamper = Stamper::new(&grid);
        let out = scarp_compute(&layers, &grid, body_radius, &rec, 0.0, None, &mut stamper);
        assert!(out.cells.len() > 100);
        out.apply(&mut layers);
        let face = 0.5 * height / SCARP_FACE_SLOPE_DEG.to_radians().tan();
        let back = SCARP_BACK_SLOPE_FACTOR * height;
        let at = |k: usize, across: f64| -> f64 {
            let p = rec.points[k];
            let prev = rec.points[k - 1];
            let next = rec.points[k + 1];
            let chord = [next[0] - prev[0], next[1] - prev[1], next[2] - prev[2]];
            let radial = dot(chord, p);
            let t = normalize([
                chord[0] - radial * p[0],
                chord[1] - radial * p[1],
                chord[2] - radial * p[2],
            ]);
            let n = normalize(cross(p, t));
            let s = across / body_radius;
            let q = normalize([p[0] + s * n[0], p[1] + s * n[1], p[2] + s * n[2]]);
            layers.bedrock.as_slice()[grid.index(grid.locate(q).0)] as f64
        };
        let m = rec.points.len() / 2;
        let footwall = at(m, -4.0 * face);
        let crest = at(m, 4.0 * face);
        let far = at(m, 3.0 * back);
        assert!(footwall.abs() < 0.05 * height, "{footwall} of {height}");
        assert!(
            crest > 0.5 * height && crest < 0.75 * height,
            "{crest} of {height}"
        );
        assert!(far < 0.1 * height, "{far} of {height}");
        // Near the ends the taper has brought it most of the way down.
        let end = at(2, 4.0 * face);
        assert!(end < 0.5 * crest, "{end} vs {crest}");
    }

    /// A thrust's trace turns across the squeeze where the field decides its strike, and
    /// stops where the lid is no longer squeezed past its strength (F155).
    #[test]
    fn a_thrust_turns_across_the_squeeze_and_stops_where_it_ends() {
        let body_radius = 1.0e6;
        let spacing = 5e3;
        // The strike runs east everywhere, wholly decided; the lid is squeezed past its
        // strength west of 20 degrees east and not beyond.
        let field = |p: [f64; 3]| Squeeze {
            stress: if p[1].atan2(p[0]).to_degrees() < 20.0 {
                1e6
            } else {
                -1e6
            },
            strike: local_basis(p).0,
            decided: 1.0,
            least: -2e6,
        };
        let start = [1.0, 0.0, 0.0];
        let (east, north) = local_basis(start);
        let heading = [
            0.2 * east[0] + north[0],
            0.2 * east[1] + north[1],
            0.2 * east[2] + north[2],
        ];
        let points = thrust_trace(
            start,
            heading,
            0.0,
            [0.0; SCARP_LOBE_HARMONICS],
            1.0e6,
            spacing,
            body_radius,
            &field,
            0.0,
            0.0,
            None,
        );
        // It stopped short of its length, at the squeeze's end.
        let last = *points.last().unwrap();
        assert!(points.len() < 150, "{}", points.len());
        let lon = last[1].atan2(last[0]).to_degrees();
        assert!(lon < 20.0 && lon > 15.0, "{lon}");
        // Its last steps run east, across the squeeze.
        let n = points.len();
        let (a, b) = (points[n - 2], points[n - 1]);
        let step = normalize([b[0] - a[0], b[1] - a[1], b[2] - a[2]]);
        assert!(
            dot(step, local_basis(a).0) > 0.95,
            "{}",
            dot(step, local_basis(a).0)
        );
    }

    /// Thrusts in systems (F160): one beside a thrust drawn before runs along a span of it
    /// about where it starts, as far off as the host's reach sets, and none joins a thrust
    /// that is not near; one running its own way meets a host at least as long within their
    /// reaches, and not a shorter one, nor one farther off.
    #[test]
    fn thrusts_gather_beside_their_hosts_and_stop_at_longer_ones() {
        let r = 2.44e6;
        let mut systems = Systems::new(SYSTEM_INDEX_WIDTH, r);
        let spacing = 5e3;
        let along = |d: f64, off: f64| {
            let (a, b) = (d / r, off / r);
            [a.cos() * b.cos(), a.sin() * b.cos(), b.sin()]
        };
        // A host along the equator, 600 km long, reaching 40 km.
        let host: Vec<[f64; 3]> = (0..=120).map(|k| along(k as f64 * spacing, 0.0)).collect();
        systems.add(host, spacing, 600e3, 40e3);
        let start = along(300e3, 0.0);
        let (trace, sp) = systems
            .beside(start, 200e3, 0.0, 0.5, 0.9, 0.5)
            .expect("a neighbour");
        assert_eq!(sp, spacing);
        let (lo, hi) = SYSTEM_SPACING;
        let expected = 40e3 * (lo + (hi - lo) * 0.5);
        for p in &trace {
            let off = p[2].asin() * r;
            assert!(
                (off - expected).abs() < 0.1 * expected,
                "{off} vs {expected}"
            );
        }
        let mid = trace[trace.len() / 2];
        assert!(angle_between(mid, along(300e3, expected)) * r < 30e3);
        assert!((trace.len() - 1) as f64 * spacing > 150e3);
        // Not when the draw says to run its own way, nor beside nothing near.
        assert!(systems.beside(start, 200e3, 0.9, 0.5, 0.9, 0.5).is_none());
        assert!(
            systems
                .beside([0.0, 0.0, 1.0], 200e3, 0.0, 0.5, 0.9, 0.5)
                .is_none()
        );
        // The arrest: within a quarter of their reaches of a host at least as long.
        assert!(systems.meets(along(300e3, 12e3), 5e3, 100e3));
        assert!(!systems.meets(along(300e3, 12e3), 5e3, 800e3));
        assert!(!systems.meets(along(300e3, 60e3), 5e3, 100e3));
    }

    #[test]
    fn lengths_follow_the_power_law_within_the_limits() {
        let (lo, hi) = (10_000.0, 1_000_000.0);
        let n = 20_000;
        let mut above_100 = 0;
        for i in 0..n {
            let l = sample_length((i as f64 + 0.5) / n as f64, lo, hi);
            assert!(l >= lo - 1e-6 && l <= hi + 1e-6);
            if l >= 100_000.0 {
                above_100 += 1;
            }
        }
        // N(>100 km) / N(>10 km) for a -1.5 slope, truncated at 1,000 km.
        let a = SCARP_LENGTH_EXPONENT;
        let expected = (100_000f64.powf(-a) - hi.powf(-a)) / (lo.powf(-a) - hi.powf(-a));
        let got = above_100 as f64 / n as f64;
        assert!((got - expected).abs() < 0.01, "{got} vs {expected}");
    }
}
