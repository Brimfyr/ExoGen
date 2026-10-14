// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! The cracked shell (PLAN.md 4.3 and Phase 4): lineae, cycloids and bands on an ice
//! shell over an ocean.
//!
//! Cracks open where the tidal tension exceeds what the ice can hold, perpendicular to
//! the most tensile stress of the moment. Where the slowly turning field of
//! non-synchronous rotation dominates, a crack runs for hundreds of kilometres along the
//! field and freezes into a double ridge; where the diurnal field dominates, the crack
//! propagates only while the tension lasts and turns with the field through the orbit,
//! so it arcs, stops, and starts again at a cusp in the next orbit: a cycloid. A share
//! of the cracks pull apart into bands of new ice. The population is drawn against the
//! crack density the tidal stress and the shell thickness give, per fixed time bin so
//! every quality sees the same cracks, within the surface's retention age; the cracks
//! are cut at the finish, reset the surface along their width, and carry the salts of
//! the ocean that fade under frost as they age. Cracks answer one another: a share of
//! them start at the tip of an older one, where its stress concentrates, and a crack
//! that comes within the stress shadow of an open fissure turns to meet it at a right
//! angle and stops there, as at a free surface, while healed ridges are crossed.

use std::collections::HashMap;

use crate::constants::*;
use crate::crater::{CellWrite, StampOutput, Stamper};
use crate::fields::LayerSet;
use crate::geom::{angle_between, cross, dot, local_basis, normalize};
use crate::grid::CubeSphere;
use crate::noise::GradientNoise;
use crate::params::ParamGraph;
use crate::rng::Rng;
use crate::stages::tidal::TidalShell;
use crate::units::*;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LineaKind {
    /// A double ridge along a crack.
    Ridge,
    /// A chain of arcs with cusps, the crack having turned with the diurnal field.
    Cycloid,
    /// A dilational band of new ice between the halves of a crack.
    Band,
    /// A fissure under a plume: a trough widened and deepened by the venting, with low
    /// flanks.
    Vent,
}

/// A crack in the shell and what grew along it.
#[derive(Clone, Debug)]
pub struct Linea {
    pub points: Vec<[f64; 3]>,
    pub spacing: f64,
    pub length: f64,
    /// Full width of the feature, m, and its relief: the ridge height for a ridge or
    /// cycloid, the depth of the floor for a band.
    pub width: f64,
    pub relief: f64,
    pub kind: LineaKind,
    /// The orbital phase the crack opened at, for the record.
    pub phase: f64,
    pub time_gyr: f64,
    pub id: u64,
}

impl Linea {
    pub fn midpoint(&self) -> [f64; 3] {
        self.points[self.points.len() / 2]
    }

    pub fn cap_radius(&self, body_radius: f64) -> f64 {
        // The trace may wander: take the farthest point from the midpoint.
        let mid = self.midpoint();
        let far = self
            .points
            .iter()
            .map(|p| crate::geom::angle_between(mid, *p))
            .fold(0.0f64, f64::max);
        far + self.width / body_radius
    }
}

/// The unit tangent of a trace at its midpoint, from the chord of its neighbours.
pub fn midpoint_tangent(points: &[[f64; 3]]) -> [f64; 3] {
    let n = points.len();
    let m = n / 2;
    let mid = points[m];
    let (a, b) = (points[m.saturating_sub(1)], points[(m + 1).min(n - 1)]);
    let chord = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let radial = dot(chord, mid);
    normalize([
        chord[0] - radial * mid[0],
        chord[1] - radial * mid[1],
        chord[2] - radial * mid[2],
    ])
}

/// The tidal axis in the shell's frame at a time: with non-synchronous rotation the
/// shell has turned against it since.
pub fn mean_axis_at(time_gyr: f64, present: f64, nsr_period_gyr: f64) -> [f64; 3] {
    if nsr_period_gyr <= 0.0 {
        return [1.0, 0.0, 0.0];
    }
    let turns = (present - time_gyr) / nsr_period_gyr;
    let lon = -2.0 * std::f64::consts::PI * turns;
    [lon.cos(), lon.sin(), 0.0]
}

/// The direction a crack runs in at a point under a stress: perpendicular to the most
/// tensile principal stress, as a unit tangent vector, in the sense nearest `along`.
pub fn crack_direction(
    p: [f64; 3],
    stress: &crate::stages::tidal::Stress,
    along: Option<[f64; 3]>,
) -> [f64; 3] {
    let (east, north) = local_basis(p);
    let (_, _, psi) = stress.principal();
    let run = psi + std::f64::consts::FRAC_PI_2;
    let mut d = normalize([
        east[0] * run.cos() + north[0] * run.sin(),
        east[1] * run.cos() + north[1] * run.sin(),
        east[2] * run.cos() + north[2] * run.sin(),
    ]);
    if let Some(a) = along
        && dot(d, a) < 0.0
    {
        d = [-d[0], -d[1], -d[2]];
    }
    d
}

/// Step a point along a tangent direction by a distance on the sphere.
fn step(p: [f64; 3], d: [f64; 3], distance: f64, body_radius: f64) -> [f64; 3] {
    let s = distance / body_radius;
    normalize([p[0] + s * d[0], p[1] + s * d[1], p[2] + s * d[2]])
}

/// The cracks already open, bucketed on the sphere, for the ones that follow.
struct CrackIndex {
    /// Buckets per unit of direction: one over the shadow range in radians.
    scale: f64,
    /// Per bucket: the crack's index, a point of it, the tangent there, and whether it
    /// is an open fissure.
    buckets: HashMap<(i32, i32, i32), Bucket>,
}

/// A bucket's entries: the crack's index, a point of it, the tangent there, open or not.
type Bucket = Vec<(usize, [f64; 3], [f64; 3], bool)>;

impl CrackIndex {
    fn new(range_rad: f64) -> Self {
        Self {
            scale: 1.0 / range_rad.max(1e-6),
            buckets: HashMap::new(),
        }
    }

    fn key(&self, p: [f64; 3]) -> (i32, i32, i32) {
        (
            (p[0] * self.scale).floor() as i32,
            (p[1] * self.scale).floor() as i32,
            (p[2] * self.scale).floor() as i32,
        )
    }

    fn insert(&mut self, crack: usize, points: &[[f64; 3]], open: bool) {
        let n = points.len();
        for i in 0..n {
            let a = points[i.saturating_sub(1)];
            let b = points[(i + 1).min(n - 1)];
            let p = points[i];
            let chord = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
            let radial = dot(chord, p);
            let t = normalize([
                chord[0] - radial * p[0],
                chord[1] - radial * p[1],
                chord[2] - radial * p[2],
            ]);
            self.buckets
                .entry(self.key(p))
                .or_default()
                .push((crack, p, t, open));
        }
    }

    /// The nearest point of another crack to a direction, within a bucket's reach: its
    /// distance, rad, its tangent and whether it is open.
    fn nearest(&self, p: [f64; 3], exclude: Option<usize>) -> Option<(f64, [f64; 3], bool)> {
        let (i, j, k) = self.key(p);
        let mut best: Option<(f64, [f64; 3], bool)> = None;
        for di in -1..=1 {
            for dj in -1..=1 {
                for dk in -1..=1 {
                    let Some(v) = self.buckets.get(&(i + di, j + dj, k + dk)) else {
                        continue;
                    };
                    for &(c, q, t, open) in v {
                        if Some(c) == exclude {
                            continue;
                        }
                        let dist = angle_between(p, q);
                        if best.is_none_or(|b| dist < b.0) {
                            best = Some((dist, t, open));
                        }
                    }
                }
            }
        }
        best
    }
}

/// What a trace answers to besides the field.
struct TraceRules<'a> {
    /// The crack belongs to the thinned polar zone: it refracts round the margin and
    /// stops in the thick shell.
    polar: bool,
    /// The cracks already open, if the trace is to answer them.
    index: Option<&'a CrackIndex>,
    /// A crack the trace does not answer: the parent it branched from.
    exclude: Option<usize>,
    /// A branch runs one way only, outward from the tip it starts at.
    one_way: Option<[f64; 3]>,
    /// Within this distance of an open fissure, m, the trace ends.
    junction: f64,
    /// The trace's wander: radians of turn over a length, m.
    wobble: (f64, f64),
    /// A branch's bias: the turn it leaves its parent with, radians, fading over a
    /// length, m, as it bends back to the field.
    bias: (f64, f64),
}

impl TraceRules<'_> {
    /// A trace that answers nothing but the field: the lanes' centre lines.
    fn free() -> Self {
        Self {
            polar: false,
            index: None,
            exclude: None,
            one_way: None,
            junction: 0.0,
            wobble: (LINEA_WOBBLE_RAD, LINEA_WOBBLE_LENGTH),
            bias: (0.0, 1.0),
        }
    }
}

/// The population of cracks on the present surface: those within the retention age,
/// drawn per time bin.
pub struct Cracks {
    pub shell: TidalShell,
    /// The shell's thickness across the body: cracks favour where it is thin.
    pub shell_map: crate::stages::ice::ShellMap,
    /// The field's largest peak tension, Pa.
    pub amplitude: f64,
    pub nsr_period_gyr: f64,
    pub density: f64,
    pub age: f64,
    pub retention: f64,
    pub body_radius: f64,
    pub orbit_period_s: f64,
    /// Whether a plume vents from the fissures of the polar terrain.
    pub vent: bool,
    /// The shell's mean thickness, m: the scale of a crack's stress shadow.
    pub mean_shell: f64,
    /// The surface gravity, m/s^2, and the elastic thickness of the shell, m, for the
    /// flexure of the thinned zone under its vents' ridges (F57).
    pub gravity: f64,
    pub elastic: f64,
    /// The past episodes of focused heating, whose fissures the body has carried from
    /// the pole.
    pub episodes: Vec<crate::stages::ice::FocusEpisode>,
}

impl Cracks {
    pub fn new(graph: &ParamGraph) -> Option<Self> {
        let density = graph.f64("tectonics.ice_crack_density");
        if density <= 0.0 {
            return None;
        }
        let mut shell = TidalShell::new(graph)?;
        let nsr_period_gyr = graph.f64("tectonics.nsr_period");
        if nsr_period_gyr > 0.0 {
            shell.nsr_shift = NSR_SHIFT_DEG.to_radians();
        }
        let amplitude = shell.amplitude();
        if amplitude <= ICE_CRACK_STRESS {
            return None;
        }
        let shell_map = crate::stages::ice::ShellMap::new(graph);
        let episodes = crate::stages::ice::focusing_episodes(graph, &shell_map);
        Some(Self {
            shell,
            shell_map,
            episodes,
            amplitude,
            nsr_period_gyr,
            density,
            age: graph.f64("body.age"),
            retention: graph.f64("cratering.retention_age"),
            body_radius: graph.f64("body.radius"),
            orbit_period_s: graph.f64("orbit.period"),
            vent: graph.f64("volcanism.plume_fallout_rate") > 0.0,
            mean_shell: graph.f64("interior.ice_shell_thickness"),
            gravity: graph.f64("body.gravity"),
            elastic: graph.f64("interior.elastic_thickness"),
        })
    }

    /// The direction a crack takes near a point in the present field, for the folds of
    /// the polar terrain, which lie along the cracks between them.
    pub fn crack_direction_near(&self, pole: [f64; 3], present: f64) -> [f64; 3] {
        let p = normalize([pole[0] + 0.3, pole[1] + 0.1, pole[2]]);
        let axis = mean_axis_at(present, present, self.nsr_period_gyr);
        let (_, phase) = self.shell.peak_tension(axis, p);
        crack_direction(p, &self.shell.total(axis, p, phase), None)
    }

    /// The expected total length of cracks on the surface, m: a transect crosses
    /// `density` cracks per metre, so the length over the area is π density / 2.
    pub fn total_length(&self) -> f64 {
        let area = 4.0 * std::f64::consts::PI * self.body_radius * self.body_radius;
        0.5 * std::f64::consts::PI * self.density * area
    }

    /// The cracks alive at the present: formed since the retention age, drawn per fixed
    /// time bin with the bin as the random index.
    pub fn generate(&self, present: f64, grid: &CubeSphere, seed: u64) -> Vec<Linea> {
        let rng = Rng::new(seed ^ 0x11CE_A000);
        let cell = grid.delta() * self.body_radius;
        let spacing = 0.5 * cell;
        let (l_lo, l_hi) = (
            LINEA_LENGTH_RAD.0 * self.body_radius,
            LINEA_LENGTH_RAD.1 * self.body_radius,
        );
        let mean_length = (l_hi - l_lo) / (l_hi / l_lo).ln();
        let expected = self.total_length() / mean_length;
        let window_start = (present - self.retention).max(0.0);
        let mut out: Vec<Linea> = Vec::new();
        let mut index = CrackIndex::new(CRACK_SHADOW_SHELLS * self.mean_shell / self.body_radius);
        // Under a plume the vents of the thinned zone come first, in their cascade (F57),
        // for the present episode and for each past one; the zone's other cracks are its
        // lesser fractures.
        if self.vent
            && self.shell_map.has_active_terrain()
            && let Some(pole) = self.shell_map.pole
        {
            self.cascade(
                pole, present, None, 0, spacing, cell, &rng, &mut out, &mut index,
            );
            for (vi, e) in self.episodes.iter().enumerate() {
                self.cascade(
                    pole,
                    present,
                    Some(e),
                    vi as u64 + 1,
                    spacing,
                    cell,
                    &rng,
                    &mut out,
                    &mut index,
                );
            }
        }
        for (bin, b0, b1) in crate::stages::time_bins(self.age, 0.0, present) {
            if b1 <= window_start {
                continue;
            }
            let b0 = b0.max(window_start);
            let share = ((b1 - b0) / self.retention.max(1e-9)).min(1.0);
            let lambda = expected * share;
            let (count, _) = crate::crater::production::poisson(&rng, bin, 0, lambda);
            for j in 0..count as u64 {
                let draw = 16 + j * 64;
                let time = b0 + (b1 - b0) * rng.uniform(bin, draw);
                let axis = mean_axis_at(time, present, self.nsr_period_gyr);
                // A share of cracks start at the tip of an older one, where its stress
                // concentrates, and run on from it; the rest where the tension above the
                // threshold and the thin shell favour (the crack density goes as 1 / d),
                // by rejection.
                let branch = (!out.is_empty() && rng.uniform(bin, draw + 60) < CRACK_BRANCH_SHARE)
                    .then(|| {
                        let i = ((rng.uniform(bin, draw + 61) * out.len() as f64) as usize)
                            .min(out.len() - 1);
                        let pts = &out[i].points;
                        let (tip, prev) = if rng.uniform(bin, draw + 62) < 0.5 {
                            (pts[0], pts[1])
                        } else {
                            (pts[pts.len() - 1], pts[pts.len() - 2])
                        };
                        let chord = [tip[0] - prev[0], tip[1] - prev[1], tip[2] - prev[2]];
                        let radial = dot(chord, tip);
                        let outward = normalize([
                            chord[0] - radial * tip[0],
                            chord[1] - radial * tip[1],
                            chord[2] - radial * tip[2],
                        ]);
                        (i, tip, outward)
                    });
                let scale = self.amplitude - ICE_CRACK_STRESS;
                let thinnest = 1.0 - self.shell_map.thinning;
                let mut start = None;
                if let Some((_, tip, _)) = branch {
                    let (peak, phase) = self.shell.peak_tension(axis, tip);
                    if peak > ICE_CRACK_STRESS {
                        start = Some((tip, phase));
                    }
                } else {
                    for attempt in 0..24u64 {
                        let p = rng.unit_direction(bin, draw + 1 + 2 * attempt);
                        let (peak, phase) = self.shell.peak_tension(axis, p);
                        let excess = peak - ICE_CRACK_STRESS;
                        if excess <= 0.0 {
                            continue;
                        }
                        let thin = thinnest / self.shell_map.factor(p).max(1e-3);
                        if rng.uniform(bin, draw + 2 + 2 * attempt)
                            < excess / scale.max(excess) * thin
                        {
                            start = Some((p, phase));
                            break;
                        }
                    }
                }
                let Some((p0, phase)) = start else {
                    continue;
                };
                // The polar terrain renews itself: a crack there is as young as its
                // surface, and opened in the field of its own time.
                let polar = self.shell_map.thin(p0);
                // The thinned zone has a history (F53): it sits at the pole now and sat
                // there through each past episode, between which the body reoriented.
                // Each episode cracked the zone alike, so a crack of the zone is drawn
                // once for the present episode, as young as the renewed surface, and once
                // for each past one, at a time within it and carried to where the body
                // left that terrain.
                let mut variants: Vec<(f64, Option<&crate::stages::ice::FocusEpisode>)> =
                    vec![(time, None)];
                if polar {
                    variants[0].0 =
                        (present - POLAR_RENEWAL_GYR * rng.uniform(bin, draw + 54)).max(0.0);
                    for e in &self.episodes {
                        variants.push((
                            e.start + (e.end - e.start) * rng.uniform(bin, draw + 54),
                            Some(e),
                        ));
                    }
                }
                for (vi, &(time, carry)) in variants.iter().enumerate() {
                    let mut axis = axis;
                    let mut phase = phase;
                    if polar {
                        axis = mean_axis_at(time, present, self.nsr_period_gyr);
                        phase = self.shell.peak_tension(axis, p0).1;
                    }
                    // What kind: the diurnal share of the tension decides between a cycloid
                    // and a ridge; a share of the ridges pull apart into bands.
                    let (d1, _, _) = self.shell.diurnal(axis, p0, phase).principal();
                    let (t1, _, _) = self.shell.total(axis, p0, phase).principal();
                    let diurnal_share = if t1 > 0.0 {
                        (d1 / t1).clamp(0.0, 1.0)
                    } else {
                        0.0
                    };
                    let kind = if diurnal_share >= LINEA_DIURNAL_SHARE {
                        LineaKind::Cycloid
                    } else if rng.uniform(bin, draw + 50) < BAND_FRACTION {
                        LineaKind::Band
                    } else {
                        LineaKind::Ridge
                    };
                    let length = l_lo * (l_hi / l_lo).powf(rng.uniform(bin, draw + 51));
                    let id = (bin << 32) | j | ((vi as u64) << 56);
                    let (width, relief) = match kind {
                        LineaKind::Band => {
                            let (w_lo, w_hi) = (
                                BAND_WIDTH_RAD.0 * self.body_radius,
                                BAND_WIDTH_RAD.1 * self.body_radius,
                            );
                            (
                                w_lo * (w_hi / w_lo).powf(rng.uniform(bin, draw + 52)),
                                BAND_DEPTH * (0.6 + 0.8 * rng.uniform(bin, draw + 53)),
                            )
                        }
                        LineaKind::Vent => (
                            VENT_WIDTH,
                            VENT_TROUGH_DEPTH * (0.7 + 0.6 * rng.uniform(bin, draw + 53)),
                        ),
                        _ => (
                            LINEA_WIDTH,
                            LINEA_RIDGE_HEIGHT * (0.6 + 0.8 * rng.uniform(bin, draw + 53)),
                        ),
                    };
                    let rules = TraceRules {
                        polar,
                        // A carried crack traces in the pole's frame, where the present
                        // cracks are not its neighbours.
                        index: if carry.is_some() { None } else { Some(&index) },
                        exclude: branch.map(|b| b.0),
                        one_way: branch.map(|b| b.2),
                        junction: (0.5 * width).max(cell),
                        wobble: (LINEA_WOBBLE_RAD, LINEA_WOBBLE_LENGTH),
                        bias: (0.0, 1.0),
                    };
                    let points = match kind {
                        LineaKind::Cycloid => {
                            self.cycloid(axis, p0, phase, spacing, length, &rules)
                        }
                        _ => self.ridge_trace(axis, p0, phase, spacing, length, id, &rules),
                    };
                    if points.len() < 3 {
                        continue;
                    }
                    // A past episode's crack is carried to where the body left its terrain,
                    // and is long healed.
                    let points: Vec<[f64; 3]> = match carry {
                        Some(e) => points.iter().map(|p| e.carry(*p)).collect(),
                        None => points,
                    };
                    index.insert(
                        out.len(),
                        &points,
                        kind == LineaKind::Vent && carry.is_none(),
                    );
                    out.push(Linea {
                        length: (points.len() - 1) as f64 * spacing,
                        points,
                        spacing,
                        width,
                        relief,
                        kind,
                        phase,
                        time_gyr: time,
                        id,
                    });
                }
            }
        }
        out
    }

    /// The fissures of a plume's thinned zone form in a cascade (F57; Hemingway, Rudolph
    /// and Manga 2020): the first opens through the thinnest shell, at the pole; the ice
    /// its eruptions pile along its flanks loads a plate broken at the fissure, whose
    /// bending tension at the surface peaks a quarter flexural parameter away on either
    /// side, where the next fissure opens, and so on outward while the site is still in
    /// the thinned core and the bending and tidal tensions together exceed the strength.
    /// Each fissure runs along the field while the thin shell lets it, as a polar crack
    /// does, and refracts round the margin. The fissures are the vents of the present
    /// episode, or the healed vents of a past one carried to where the body left its
    /// terrain.
    #[allow(clippy::too_many_arguments)]
    fn cascade(
        &self,
        pole: [f64; 3],
        present: f64,
        carry: Option<&crate::stages::ice::FocusEpisode>,
        vi: u64,
        spacing: f64,
        cell: f64,
        rng: &Rng,
        out: &mut Vec<Linea>,
        index: &mut CrackIndex,
    ) {
        let bin = u64::MAX - vi;
        let time = match carry {
            None => present,
            Some(e) => 0.5 * (e.start + e.end),
        };
        let axis = mean_axis_at(time, present, self.nsr_period_gyr);
        let l_max = LINEA_LENGTH_RAD.1 * self.body_radius;
        let junction = (0.5 * VENT_WIDTH).max(cell);
        // The load of a fissure's flanking ridges, N/m: the two Gaussian flanks of the
        // vent profile, each VENT_RIDGE_SHARE of the trough deep over an eighth of the
        // width.
        let load = DENSITY_ICE
            * self.gravity
            * 2.0
            * VENT_RIDGE_SHARE
            * VENT_TROUGH_DEPTH
            * std::f64::consts::PI.sqrt()
            * (VENT_WIDTH / 8.0);
        // A line load V on a plate broken at the load (Turcotte and Schubert): the
        // deflection V alpha^3 / (2 D) e^(-x / alpha) cos(x / alpha), whose surface
        // tension peaks at x = alpha pi / 4 with the moment 0.3225 V alpha and the stress
        // 6 M / te^2, with the elastic thickness the shell's local thinning leaves.
        let site = |p: [f64; 3]| -> (f64, f64) {
            let te = (self.elastic * self.shell_map.factor(p)).max(1.0);
            let alpha = crate::stages::ice_shell::flexural_parameter(te, self.gravity);
            (
                std::f64::consts::FRAC_PI_4 * alpha,
                6.0 * 0.3225 * load * alpha / (te * te),
            )
        };
        // The first fissure: through the thinnest shell, a little off the pole.
        let (east, north) = local_basis(pole);
        let (d0, _) = site(pole);
        let off = VENT_CASCADE_FIRST_OFFSET * d0 * rng.uniform(bin, 1) / self.body_radius;
        let ang = std::f64::consts::TAU * rng.uniform(bin, 2);
        let p0 = normalize([
            pole[0] + off * (east[0] * ang.cos() + north[0] * ang.sin()),
            pole[1] + off * (east[1] * ang.cos() + north[1] * ang.sin()),
            pole[2] + off * (east[2] * ang.cos() + north[2] * ang.sin()),
        ]);
        let (peak, phase) = self.shell.peak_tension(axis, p0);
        if peak <= ICE_CRACK_STRESS {
            return;
        }
        let id_of = |k: u64| (VENT_CASCADE_BIN << 32) | k | (vi << 56);
        let first = {
            let rules = TraceRules {
                polar: true,
                index: if carry.is_none() { Some(&*index) } else { None },
                exclude: None,
                one_way: None,
                junction,
                wobble: (VENT_WOBBLE_RAD, VENT_WOBBLE_LENGTH),
                bias: (0.0, 1.0),
            };
            self.ridge_trace(axis, p0, phase, spacing, l_max, id_of(0), &rules)
        };
        if first.len() < 3 {
            return;
        }
        let push = |pts: &[[f64; 3]],
                    phase: f64,
                    k: u64,
                    scale: f64,
                    out: &mut Vec<Linea>,
                    index: &mut CrackIndex| {
            let time_gyr = match carry {
                None => (present - POLAR_RENEWAL_GYR * rng.uniform(bin, 10 + 2 * k)).max(0.0),
                // A past episode's fissures vented to its end and stand on the last
                // surface its terrain renewed.
                Some(e) => e.end + FOSSIL_VENT_LAG_GYR * (1.0 + rng.uniform(bin, 10 + 2 * k)),
            };
            let relief = scale * VENT_TROUGH_DEPTH * (0.7 + 0.6 * rng.uniform(bin, 11 + 2 * k));
            let points: Vec<[f64; 3]> = match carry {
                Some(e) => pts.iter().map(|p| e.carry(*p)).collect(),
                None => pts.to_vec(),
            };
            index.insert(out.len(), &points, carry.is_none());
            out.push(Linea {
                length: (points.len() - 1) as f64 * spacing,
                points,
                spacing,
                width: scale * VENT_WIDTH,
                relief,
                kind: LineaKind::Vent,
                phase,
                time_gyr,
                id: id_of(k),
            });
        };
        // Each fissure ends in a fork: from each tip a branch leaves at an angle and bends
        // back to the field, the hooked termination, as young as its parent and narrower.
        let fork = |pts: &[[f64; 3]], k: u64, out: &mut Vec<Linea>, index: &mut CrackIndex| {
            let n = pts.len();
            let parent = out.len() - 1;
            let length = VENT_FORK_LENGTH_SHARE * (n as f64 - 1.0) * spacing;
            for (end, (tip, prev)) in [(pts[0], pts[1]), (pts[n - 1], pts[n - 2])]
                .into_iter()
                .enumerate()
            {
                let chord = [tip[0] - prev[0], tip[1] - prev[1], tip[2] - prev[2]];
                let radial = dot(chord, tip);
                let outward = normalize([
                    chord[0] - radial * tip[0],
                    chord[1] - radial * tip[1],
                    chord[2] - radial * tip[2],
                ]);
                let (peak, phase) = self.shell.peak_tension(axis, tip);
                if peak <= ICE_CRACK_STRESS {
                    continue;
                }
                let kf = 64 + 2 * k + end as u64;
                let sign = if rng.uniform(bin, 40 + kf) < 0.5 {
                    1.0
                } else {
                    -1.0
                };
                let branch = {
                    let rules = TraceRules {
                        polar: true,
                        index: if carry.is_none() { Some(&*index) } else { None },
                        exclude: Some(parent),
                        one_way: Some(outward),
                        junction,
                        wobble: (VENT_WOBBLE_RAD, VENT_WOBBLE_LENGTH),
                        bias: (
                            sign * VENT_FORK_ANGLE_DEG.to_radians(),
                            VENT_FORK_BIAS_LENGTH,
                        ),
                    };
                    self.ridge_trace(axis, tip, phase, spacing, length, id_of(kf), &rules)
                };
                if branch.len() >= 3 {
                    push(&branch, phase, kf, VENT_FORK_WIDTH_SHARE, out, index);
                }
            }
        };
        let mut k = 0u64;
        push(&first, phase, k, 1.0, out, index);
        fork(&first, k, out, index);
        // Outward from the first on either side, each where the last one's load bends
        // the plate most, while the shell there is still the thinned core's and the
        // tension opens it.
        for sign in [1.0f64, -1.0] {
            let mut last = first.clone();
            for _ in 0..VENT_CASCADE_MAX {
                let mid = last[last.len() / 2];
                let tangent = midpoint_tangent(&last);
                let side = normalize(cross(mid, tangent));
                let side = [sign * side[0], sign * side[1], sign * side[2]];
                let (distance, _) = site(mid);
                let theta = distance / self.body_radius;
                let p = normalize([
                    mid[0] * theta.cos() + side[0] * theta.sin(),
                    mid[1] * theta.cos() + side[1] * theta.sin(),
                    mid[2] * theta.cos() + side[2] * theta.sin(),
                ]);
                if self.shell_map.profile(p) <= POLAR_CORE_PROFILE {
                    break;
                }
                let (_, bending) = site(p);
                let (peak, phase) = self.shell.peak_tension(axis, p);
                if peak + bending <= ICE_CRACK_STRESS {
                    break;
                }
                k += 1;
                let pts = {
                    let rules = TraceRules {
                        polar: true,
                        index: if carry.is_none() { Some(&*index) } else { None },
                        exclude: None,
                        one_way: None,
                        junction,
                        wobble: (VENT_WOBBLE_RAD, VENT_WOBBLE_LENGTH),
                        bias: (0.0, 1.0),
                    };
                    self.ridge_trace(axis, p, phase, spacing, l_max, id_of(k), &rules)
                };
                if pts.len() < 3 {
                    break;
                }
                push(&pts, phase, k, 1.0, out, index);
                fork(&pts, k, out, index);
                last = pts;
            }
        }
    }

    /// How a propagating crack answers the cracks already open near it (Pollard and
    /// Aydin 1988): within the stress shadow of an open fissure, half a shell thickness,
    /// it turns toward meeting it at a right angle, and at the junction range it stops,
    /// having reached a free surface; healed ridges and bands it crosses. None: stop.
    fn steer(&self, p: [f64; 3], d: [f64; 3], rules: &TraceRules) -> Option<[f64; 3]> {
        let Some(index) = rules.index else {
            return Some(d);
        };
        let Some((dist_rad, tangent, open)) = index.nearest(p, rules.exclude) else {
            return Some(d);
        };
        if !open {
            return Some(d);
        }
        let dist = dist_rad * self.body_radius;
        if dist <= rules.junction {
            return None;
        }
        let range = CRACK_SHADOW_SHELLS * self.mean_shell * self.shell_map.factor(p).max(1e-3);
        if dist >= range {
            return Some(d);
        }
        let mut perp = normalize(cross(p, tangent));
        if dot(perp, d) < 0.0 {
            perp = [-perp[0], -perp[1], -perp[2]];
        }
        let w = CRACK_DEFLECTION * (1.0 - dist / range).powi(2);
        Some(normalize([
            d[0] + (perp[0] - d[0]) * w,
            d[1] + (perp[1] - d[1]) * w,
            d[2] + (perp[2] - d[2]) * w,
        ]))
    }

    /// A centre line for a lane of grooved terrain, traced along the field from a start
    /// in both directions as a ridge is, with its wander and without the polar rules.
    pub fn trace(
        &self,
        axis: [f64; 3],
        p0: [f64; 3],
        phase: f64,
        spacing: f64,
        length: f64,
        id: u64,
    ) -> Vec<[f64; 3]> {
        self.ridge_trace(axis, p0, phase, spacing, length, id, &TraceRules::free())
    }

    /// A crack along the slowly turning field: from the start in both directions,
    /// perpendicular to the most tensile stress at the phase it opened, while the
    /// tension holds, with a little wander. A crack of the thinned polar zone refracts
    /// round the zone's margin and stops in the thick shell beyond; a branch runs one
    /// way, outward from the tip it starts at, for the whole length.
    #[allow(clippy::too_many_arguments)]
    fn ridge_trace(
        &self,
        axis: [f64; 3],
        p0: [f64; 3],
        phase: f64,
        spacing: f64,
        length: f64,
        id: u64,
        rules: &TraceRules,
    ) -> Vec<[f64; 3]> {
        let polar = rules.polar;
        let noise = GradientNoise::new(id ^ 0x7E1A_0000);
        let steps = (length / spacing).ceil() as usize;
        let half = if rules.one_way.is_some() {
            steps
        } else {
            steps.div_ceil(2)
        };
        let mut sides: [Vec<[f64; 3]>; 2] = [Vec::new(), Vec::new()];
        for (side, sign) in [1.0f64, -1.0].iter().enumerate() {
            if side == 1 && rules.one_way.is_some() {
                break;
            }
            let mut p = p0;
            let mut along: Option<[f64; 3]> = rules.one_way;
            for i in 0..half {
                let stress = self.shell.total(axis, p, phase);
                let (s1, _, _) = stress.principal();
                if s1 < ICE_CRACK_STRESS && i > 2 {
                    break;
                }
                let mut d = crack_direction(p, &stress, along);
                if along.is_none() {
                    d = [sign * d[0], sign * d[1], sign * d[2]];
                }
                // Wander: a slow turn by noise along the trace; and for a branch the
                // bias it leaves its parent with, fading as it bends back to the field.
                let (wobble, wobble_length) = rules.wobble;
                let turn = wobble
                    * noise.fbm(
                        [i as f64 * spacing / wobble_length, 0.3, 0.7],
                        1.0,
                        2,
                        2.0,
                        0.5,
                    )
                    / 0.6
                    + rules.bias.0 * (-(i as f64 * spacing) / rules.bias.1.max(1.0)).exp();
                let side_dir = normalize(cross(p, d));
                let mut d = normalize([
                    d[0] * turn.cos() + side_dir[0] * turn.sin(),
                    d[1] * turn.cos() + side_dir[1] * turn.sin(),
                    d[2] * turn.cos() + side_dir[2] * turn.sin(),
                ]);
                if polar {
                    // Leaving the thinned core, the crack meets a thicker, stiffer shell
                    // and refracts along the contrast: it turns to run round the margin,
                    // and in the thick shell beyond the tension no longer opens it.
                    let profile = self.shell_map.profile(p);
                    if profile < CRACK_STOP_PROFILE && i > 2 {
                        break;
                    }
                    let w = ((POLAR_CORE_PROFILE - profile)
                        / (POLAR_CORE_PROFILE - POLAR_OUTER_PROFILE))
                        .clamp(0.0, 1.0);
                    if w > 0.0
                        && let Some(pole) = self.shell_map.pole
                    {
                        let mut t = normalize(cross(pole, p));
                        if dot(t, d) < 0.0 {
                            t = [-t[0], -t[1], -t[2]];
                        }
                        d = normalize([
                            d[0] + (t[0] - d[0]) * w,
                            d[1] + (t[1] - d[1]) * w,
                            d[2] + (t[2] - d[2]) * w,
                        ]);
                    }
                }
                // The cracks already open: a turn toward an open fissure's perpendicular,
                // or the end against it.
                let Some(d) = self.steer(p, d, rules) else {
                    break;
                };
                p = step(p, d, spacing, self.body_radius);
                along = Some(d);
                sides[side].push(p);
            }
        }
        let mut points: Vec<[f64; 3]> = sides[1].iter().rev().copied().collect();
        points.push(p0);
        points.extend(sides[0].iter().copied());
        points
    }

    /// A cycloid: the crack runs while the diurnal tension exceeds the threshold,
    /// turning with the field as the orbit advances, stops when it falls, and starts
    /// again at a cusp when it rises in the next orbit.
    fn cycloid(
        &self,
        axis: [f64; 3],
        p0: [f64; 3],
        phase0: f64,
        spacing: f64,
        length: f64,
        rules: &TraceRules,
    ) -> Vec<[f64; 3]> {
        let dt = self.orbit_period_s / CYCLOID_PHASE_STEPS as f64;
        let d_phase = 2.0 * std::f64::consts::PI / CYCLOID_PHASE_STEPS as f64;
        let advance = CYCLOID_SPEED * dt;
        let mut points = vec![p0];
        let mut p = p0;
        let mut phase = phase0;
        let mut along: Option<[f64; 3]> = rules.one_way;
        let mut run = 0.0;
        let mut since_point = 0.0;
        let mut idle = 0usize;
        let max_steps = CYCLOID_MAX_ORBITS * CYCLOID_PHASE_STEPS;
        for _ in 0..max_steps {
            let stress = self.shell.total(axis, p, phase);
            let (s1, _, _) = stress.principal();
            if s1 >= ICE_CRACK_STRESS {
                idle = 0;
                let d = crack_direction(p, &stress, along);
                let Some(d) = self.steer(p, d, rules) else {
                    break;
                };
                p = step(p, d, advance, self.body_radius);
                along = Some(d);
                run += advance;
                since_point += advance;
                if since_point >= spacing {
                    points.push(p);
                    since_point = 0.0;
                }
                if run >= length {
                    break;
                }
            } else {
                // The crack waits for the tension to return; it gives up after an orbit
                // without it.
                idle += 1;
                if idle > CYCLOID_PHASE_STEPS {
                    break;
                }
            }
            phase += d_phase;
        }
        points
    }
}

/// The cross-profile of a crack's relief at a distance from its trace, m: a double
/// ridge with a trough between for a ridge or cycloid, a floor with soft walls for a
/// band.
pub fn linea_profile(kind: LineaKind, d: f64, width: f64, relief: f64, blur: f64) -> f64 {
    let x = d.abs();
    match kind {
        LineaKind::Band => {
            let wall = (0.1 * width).hypot(blur).max(1.0);
            -relief * 0.5 * (1.0 - ((x - 0.5 * width) / wall).tanh())
        }
        LineaKind::Vent => {
            let w = width.hypot(2.0 * blur).max(1.0);
            let sigma = (w / 4.0).max(1.0);
            let trough = -relief * (-(x / sigma).powi(2)).exp();
            let crest = 0.6 * w;
            let flanks = VENT_RIDGE_SHARE * relief * (-((x - crest) / (0.5 * sigma)).powi(2)).exp();
            trough + flanks
        }
        _ => {
            let w = width.hypot(2.0 * blur).max(1.0);
            let crest = 0.25 * w;
            let sigma = (w / 6.0).max(1.0);
            let ridges = relief * (-((x - crest) / sigma).powi(2)).exp();
            let trough = LINEA_TROUGH_FRACTION * relief * (-(x / sigma).powi(2)).exp();
            ridges - trough
        }
    }
}

/// Stamp a crack on the resolved grid: its relief across the trace, and its unit, age
/// and salts along it. A ridge or a cycloid that opened under loose cover deeper than
/// its own relief built no ridges at the surface: the cover drained into it, and it
/// shows as a chain of pits (F64), as wide and as deep as the cover that lies on it at
/// `present`.
#[allow(clippy::too_many_arguments)]
pub fn linea_compute(
    layers: &LayerSet,
    grid: &CubeSphere,
    body_radius: f64,
    rec: &Linea,
    blur: f64,
    loose: Option<&crate::stages::cryo::LooseCover>,
    present: f64,
    stamper: &mut Stamper,
) -> StampOutput {
    let bedrock = layers.bedrock.as_slice();
    let ages = layers.surface_age.as_slice();
    let cell = grid.delta() * body_radius;
    let width = rec.width.max(2.0 * cell);
    // A ridge or a cycloid is a chain of segments (F62): where each ends along the
    // crack, as a share of its length, and how far each stands aside from the trace, to
    // the one side from each to the next. A band and a vent are one opening.
    let rng = Rng::new(rec.id ^ 0x11CE_5E61);
    let chained = !matches!(rec.kind, LineaKind::Band | LineaKind::Vent);
    let segments = if chained {
        (rec.length / (CRACK_SEGMENT_WIDTHS * width))
            .round()
            .max(1.0) as usize
    } else {
        1
    };
    let bounds: Vec<f64> = (0..=segments)
        .map(|i| {
            if i == 0 {
                0.0
            } else if i == segments {
                1.0
            } else {
                (i as f64 + 0.6 * (rng.uniform(i as u64, 1) - 0.5)) / segments as f64
            }
        })
        .collect();
    let sense = if rng.uniform(0, 0) < 0.5 { -1.0 } else { 1.0 };
    let asides: Vec<f64> = (0..segments)
        .map(|i| {
            sense
                * CRACK_STEP_SHARE
                * width
                * (i as f64 - 0.5 * (segments - 1) as f64)
                * (0.5 + rng.uniform(i as u64, 2))
        })
        .collect();
    let aside_reach = asides.iter().fold(0.0f64, |m, a| m.max(a.abs()));
    let reach = match rec.kind {
        LineaKind::Band => 0.5 * width + 0.3 * width + 2.0 * cell,
        _ => width + 2.0 * cell + aside_reach,
    };
    let half_band = 0.5 * rec.spacing + cell;
    let n = rec.points.len();
    // Per cell: the nearest point along the crack, the signed distance across it, and
    // how far along the crack's length that point lies.
    let mut best: HashMap<usize, (f64, f64, f64)> = HashMap::new();
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
        stamper.gather(grid, body_radius, p, reach);
        for &(idx, _, _) in &stamper.footprint {
            let d = grid.center(grid.cell(idx));
            let along = dot(d, t) * body_radius;
            if along.abs() > half_band {
                continue;
            }
            let across = dot(d, normal) * body_radius;
            let entry = best.entry(idx).or_insert((f64::INFINITY, 0.0, 0.0));
            if along.abs() < entry.0 {
                *entry = (along.abs(), across, k as f64 / (n - 1).max(1) as f64);
            }
        }
    }
    let time = rec.time_gyr as f32;
    let lineation_spacing = BAND_LINEATION_SPACING.max(3.0 * cell);
    let variation = GradientNoise::new(rec.id ^ 0x11CE_7A47);
    let mut cells: Vec<CellWrite> = best
        .into_iter()
        .map(|(idx, (_, across, s))| {
            // The segment of the chain the cell lies along: how far it stands aside, and
            // how open it is there, an ellipse along the segment and along the crack.
            let segment = (0..segments)
                .find(|&i| s < bounds[i + 1])
                .unwrap_or(segments - 1);
            let u = ((s - bounds[segment]) / (bounds[segment + 1] - bounds[segment]).max(1e-9))
                .clamp(0.0, 1.0);
            let across = across - asides[segment];
            let open = if chained {
                (CRACK_TIP_RELIEF
                    + (1.0 - CRACK_TIP_RELIEF) * (4.0 * u * (1.0 - u)).max(0.0).sqrt())
                    * (4.0 * s * (1.0 - s)).max(0.0).powf(0.25)
            } else {
                1.0
            };
            // A band is a lens, widest in the middle and closing at its ends, and is
            // lineated inside along its axis.
            // A vent is an elliptical opening, narrowing and shallowing to its tips; and
            // every crack's relief and width rise and fall along it, in segments.
            let lens = (4.0 * s * (1.0 - s)).max(0.0).sqrt();
            let vary = (1.0
                + CRACK_RELIEF_VARIATION
                    * variation.fbm(
                        [s * rec.length / CRACK_VARIATION_LENGTH, 0.3, 0.7],
                        1.0,
                        2,
                        2.0,
                        0.5,
                    )
                    / 0.6)
                .max(0.3);
            let width_here = match rec.kind {
                LineaKind::Band => {
                    width
                        * (4.0 * s * (1.0 - s))
                            .max(0.0)
                            .powf(BAND_TAPER_POWER)
                            .max(0.05)
                }
                LineaKind::Vent => width * (VENT_TIP_WIDTH + (1.0 - VENT_TIP_WIDTH) * lens),
                _ => width,
            };
            let width_here = width_here * vary;
            let relief_here = vary
                * match rec.kind {
                    LineaKind::Vent => {
                        rec.relief * (VENT_TIP_DEPTH + (1.0 - VENT_TIP_DEPTH) * lens)
                    }
                    _ => rec.relief * open,
                };
            let mut h = linea_profile(rec.kind, across, width_here, relief_here, blur);
            // Under loose cover (F64): by how much the crack shows as pits and not as
            // ridges, from the cover that lay there when it opened against its relief;
            // and the pits, from the cover that lies there now.
            let mut buried = 0.0;
            if chained && let Some(cover) = loose {
                let d = grid.center(grid.cell(idx));
                let since = ages[idx] as f64;
                let x = (cover.thickness_between(d, since, rec.time_gyr) / rec.relief.max(1.0))
                    .clamp(0.0, 1.0);
                buried = x * x * (3.0 - 2.0 * x);
                if buried > 0.0 {
                    let snow = cover.drained(cover.thickness_between(d, since, present));
                    let pit_width = PIT_WIDTH_COVERS * snow;
                    let mut pits = 0.0f64;
                    if pit_width >= PIT_MIN_CELLS * cell && across.abs() < pit_width {
                        let step = PIT_SPACING * pit_width;
                        let along = s * rec.length;
                        let offset = rng.uniform(0, 60);
                        let nearest = (along / step + offset).floor() as i64;
                        for m in nearest - 1..=nearest + 1 {
                            let pit = (m + (1 << 30)) as u64;
                            if rng.uniform(pit, 61) > open {
                                continue;
                            }
                            let centre =
                                (m as f64 + 0.5 + PIT_JITTER * (rng.uniform(pit, 62) - 0.5)
                                    - offset)
                                    * step;
                            let radius = 0.5 * pit_width * (0.7 + 0.6 * rng.uniform(pit, 63));
                            let rho = across.hypot(along - centre) / radius;
                            if rho < 1.0 {
                                pits = pits.max((1.0 - rho * rho).powf(1.5));
                            }
                        }
                    }
                    h = h * (1.0 - buried) - snow * pits * buried;
                }
            }
            let inside = across.abs() < 0.5 * width_here && buried < 0.5;
            if rec.kind == LineaKind::Band && inside {
                let edge = 1.0 - (across / (0.5 * width_here)).powi(2);
                h += BAND_LINEATION
                    * edge
                    * (2.0 * std::f64::consts::PI * across / lineation_spacing).sin();
            }
            let mut write = CellWrite {
                idx: idx as u32,
                bedrock: (bedrock[idx] as f64 + h) as f32,
                surface: false,
                age: 0.0,
                unit: 0,
                deep: 0.0,
                melt: 0.0,
                cover: 0.0,
            };
            if inside && ages[idx] <= time + 1e-6 {
                write.surface = true;
                write.age = time;
                write.unit = match rec.kind {
                    LineaKind::Band => UNIT_BAND,
                    _ => UNIT_LINEA,
                };
                write.deep = 0.0;
                write.melt = 0.0;
                // The ocean's salts, full strength when the crack is fresh.
                write.cover = 1.0;
            }
            write
        })
        .collect();
    cells.sort_by_key(|c| c.idx);
    StampOutput {
        cells,
        composition: None,
    }
}

/// Where the lineae cross a cell as a unit: the share of the surface that is crack.
pub fn crack_fraction(layers: &LayerSet) -> f64 {
    let units = layers.unit_id.as_slice();
    let n = units
        .iter()
        .filter(|&&u| u == UNIT_LINEA || u == UNIT_BAND)
        .count();
    n as f64 / units.len().max(1) as f64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn profiles_have_their_shapes() {
        // A double ridge: higher at a quarter width than at the centre or the edge.
        let centre = linea_profile(LineaKind::Ridge, 0.0, 2e3, 150.0, 0.0);
        let crest = linea_profile(LineaKind::Ridge, 500.0, 2e3, 150.0, 0.0);
        let edge = linea_profile(LineaKind::Ridge, 1.5e3, 2e3, 150.0, 0.0);
        assert!(crest > centre && crest > edge, "{centre} {crest} {edge}");
        // A band: a flat floor below the ground with soft walls.
        let floor = linea_profile(LineaKind::Band, 0.0, 20e3, 100.0, 0.0);
        let outside = linea_profile(LineaKind::Band, 15e3, 20e3, 100.0, 0.0);
        assert!(floor < -90.0 && outside > -5.0, "{floor} {outside}");
    }

    /// Under a plume the vents of the thinned zone come in a cascade (F57): a parallel
    /// set a quarter flexural parameter apart, each long enough to cross the core.
    #[test]
    fn the_tiger_stripes_form_a_parallel_cascade() {
        let (g, _) = crate::presets::preset("Enceladus").unwrap().to_graph();
        let cracks = Cracks::new(&g).unwrap();
        let present = g.f64("body.age");
        let all = cracks.generate(present, &CubeSphere::new(128), 7);
        let stripes: Vec<&Linea> = all
            .iter()
            .filter(|l| l.kind == LineaKind::Vent && l.time_gyr > present - 0.01)
            .collect();
        let r = cracks.body_radius;
        let longest = stripes.iter().map(|l| l.length).fold(0.0, f64::max);
        let first = stripes
            .iter()
            .find(|l| l.id == VENT_CASCADE_BIN << 32)
            .expect("the first stripe");
        let (mut nearest, mut parallel) = (f64::INFINITY, 0.0);
        for l in &stripes {
            if l.id == first.id {
                continue;
            }
            let d = angle_between(first.midpoint(), l.midpoint()) * r;
            if d < nearest {
                nearest = d;
                parallel = dot(midpoint_tangent(&first.points), midpoint_tangent(&l.points)).abs();
            }
        }
        eprintln!(
            "enceladus: {} stripes, longest {:.0} km, the first's neighbour {:.0} km away, parallel {parallel:.2}; {} cracks in all",
            stripes.len(),
            longest / 1e3,
            nearest / 1e3,
            all.len()
        );
        assert!(stripes.len() >= 2, "{} stripes", stripes.len());
        assert!(longest > 80e3, "longest stripe {:.0} km", longest / 1e3);
        assert!(
            (15e3..=90e3).contains(&nearest),
            "neighbour {:.0} km away",
            nearest / 1e3
        );
        assert!(parallel > 0.6, "parallel {parallel:.2}");
    }

    /// A ridge is a chain of segments (F62): its relief is greatest in the middle of a
    /// segment, falls where one hands over to the next, and closes at the crack's tips.
    #[test]
    fn a_crack_opens_as_a_chain_of_segments() {
        // A grid of 500 m cells, and a ridge 96 km long and 2 km wide along the equator.
        let grid = CubeSphere::new(256);
        let r = 500.0 / grid.delta();
        let spacing = 250.0;
        let n = (96e3 / spacing) as usize;
        let points: Vec<[f64; 3]> = (0..=n)
            .map(|i| {
                let a = (i as f64 * spacing - 48e3) / r;
                [a.cos(), a.sin(), 0.0]
            })
            .collect();
        let rec = Linea {
            points,
            spacing,
            length: 96e3,
            width: 2e3,
            relief: 150.0,
            kind: LineaKind::Ridge,
            phase: 0.0,
            time_gyr: 1.0,
            id: 7,
        };
        let layers = LayerSet::new(grid);
        let out = crate::crater::with_stamper(&grid, |stamper| {
            linea_compute(&layers, &grid, r, &rec, 0.0, None, 1.0, stamper)
        });
        // The crest along the crack: the highest cell in each kilometre of its length.
        let mut crest = vec![0.0f32; 96];
        for c in &out.cells {
            let d = grid.center(grid.cell(c.idx as usize));
            let along = d[1].atan2(d[0]) * r + 48e3;
            if (0.0..96e3).contains(&along) {
                let b = (along / 1e3) as usize;
                crest[b] = crest[b].max(c.bedrock);
            }
        }
        let top = crest.iter().cloned().fold(0.0f32, f32::max);
        assert!(top > 60.0 && top < 260.0, "{top}");
        assert!(crest[0] < 0.5 * top && crest[95] < 0.5 * top);
        // Four segments in 96 km of a 2 km crack: somewhere away from the tips the
        // relief falls where one hands over to the next.
        let low = crest[12..84].iter().cloned().fold(f32::INFINITY, f32::min);
        assert!(low < 0.6 * top, "{low} against {top}");
    }

    /// A crack that opens under loose cover deeper than its relief builds no ridges at
    /// the surface (F64): the cover drains into it, and pits follow its line.
    #[test]
    fn a_crack_under_loose_cover_shows_as_pits() {
        let grid = CubeSphere::new(256);
        let r = 500.0 / grid.delta();
        let spacing = 250.0;
        let n = (96e3 / spacing) as usize;
        let points: Vec<[f64; 3]> = (0..=n)
            .map(|i| {
                let a = (i as f64 * spacing - 48e3) / r;
                [a.cos(), a.sin(), 0.0]
            })
            .collect();
        let rec = Linea {
            points,
            spacing,
            length: 96e3,
            width: 2e3,
            relief: 150.0,
            kind: LineaKind::Ridge,
            phase: 0.0,
            time_gyr: 1.0,
            id: 7,
        };
        let layers = LayerSet::new(grid);
        // A plume at the crack's middle that laid 300 m there before the crack opened.
        let cover = crate::stages::cryo::LooseCover::of_plumes(
            vec![([1.0, 0.0, 0.0], 0.0, 0.5)],
            600.0,
            1e3,
        );
        let out = crate::crater::with_stamper(&grid, |stamper| {
            linea_compute(&layers, &grid, r, &rec, 0.0, Some(&cover), 1.0, stamper)
        });
        // The middle of the crack, where the cover is deepest.
        let middle: Vec<&CellWrite> = out
            .cells
            .iter()
            .filter(|c| {
                let d = grid.center(grid.cell(c.idx as usize));
                (d[1].atan2(d[0]) * r).abs() < 8e3
            })
            .collect();
        assert!(middle.len() > 100);
        let high = middle.iter().map(|c| c.bedrock).fold(f32::MIN, f32::max);
        let low = middle.iter().map(|c| c.bedrock).fold(f32::MAX, f32::min);
        assert!(high < 20.0, "a ridge of {high} m stands under the cover");
        assert!(low < -100.0 && low > -301.0, "{low}");
        // The surface there is the cover's, not the crack's.
        assert!(middle.iter().all(|c| !c.surface));
        // Without the cover the same crack builds its ridges.
        let bare = crate::crater::with_stamper(&grid, |stamper| {
            linea_compute(&layers, &grid, r, &rec, 0.0, None, 1.0, stamper)
        });
        assert!(
            bare.cells
                .iter()
                .map(|c| c.bedrock)
                .fold(f32::MIN, f32::max)
                > 60.0
        );
    }
}
