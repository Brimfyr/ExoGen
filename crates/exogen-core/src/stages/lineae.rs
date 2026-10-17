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
use crate::geom::{angle_between, cross, dot, local_basis, normalize, rotate};
use crate::grid::CubeSphere;
use crate::grid::projection::{PaddedField, sample_bilinear};
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

    /// The cap a stamp of the crack touches on a grid of `cell`, radians about its
    /// midpoint: its farthest point, the trace wandering, and the stamp's reach beyond
    /// the trace, so that stamps whose caps are apart never touch the same cell.
    pub fn cap_radius(&self, body_radius: f64, cell: f64) -> f64 {
        let mid = self.midpoint();
        let far = self
            .points
            .iter()
            .map(|p| crate::geom::angle_between(mid, *p))
            .fold(0.0f64, f64::max);
        far + self.reach(cell) / body_radius
    }

    /// How far the crack's stamp reaches from its trace on a grid of `cell`, m: a band
    /// its floor and walls, a vent its trough and flanks, a ridge or a cycloid its
    /// margins, offset by its segments.
    pub fn reach(&self, cell: f64) -> f64 {
        let width = self.width.max(2.0 * cell);
        match self.kind {
            // A band pushes the plates it opens between over half its length (F83).
            LineaKind::Band => (0.5 * self.length).max(width + 2.0 * cell),
            LineaKind::Vent => width + 2.0 * cell,
            _ => {
                let (_, asides) = chain(self, width);
                let aside = asides.iter().fold(0.0f64, |m, a| m.max(a.abs()));
                let margin = (LINEA_MARGIN_SHARE * self.width)
                    .hypot(0.5 * cell)
                    .max(LINEA_STAIN_MIN_CELLS * cell);
                (width + 2.0 * cell + aside).max(2.5 * margin + aside)
            }
        }
    }
}

/// A ridge or a cycloid is a chain of segments (F62): where each ends along the crack, as
/// a share of its length, and how far each stands aside from the trace, to the one side
/// from each to the next, for a crack drawn `width` wide. A band and a vent are one
/// opening.
fn chain(rec: &Linea, width: f64) -> (Vec<f64>, Vec<f64>) {
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
    (bounds, asides)
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

/// The heading a running crack takes at a point where the field would have it run
/// `toward`, coming in along `prev`, under principal stresses `s1` and `s2`, most tensile
/// first. It kinks by the maximum tangential stress criterion (Erdogan and Sih 1963):
/// misaligned by an angle, it bears a normal stress and a shear on its own plane, and turns
/// by 2 atan((K_I - (K_I² + 8 K_II²)^1/2) / (4 K_II)), the two loadings standing for K_I
/// and K_II, never past the field's heading. Under near-equal tension in all directions the
/// shear vanishes and it runs on straight; misaligned in a strongly anisotropic field it
/// kinks sharply, as at a cycloid's cusp; a crack held closed, in pure shear, kinks by
/// the criterion's 70.5 degrees, or to the field if it is nearer. Near a right angle the
/// nearer of the field's two headings is ambiguous, and the crack keeps kinking the way it
/// last kinked (`sense`, updated), as a crack that has begun to turn does.
pub fn kinked(
    prev: [f64; 3],
    toward: [f64; 3],
    p: [f64; 3],
    s1: f64,
    s2: f64,
    sense: &mut f64,
) -> [f64; 3] {
    // The heading the crack came in on, laid in the tangent plane here: a step on the
    // sphere tilts it out of the plane it was tangent to.
    let prev = in_plane(prev, p);
    let misalign = dot(prev, toward).clamp(-1.0, 1.0).acos();
    if misalign < 1e-9 {
        return toward;
    }
    let (normal, shear) = plane_loading(prev, toward, p, s1, s2);
    let turn = if normal <= 0.0 {
        misalign.min((1.0f64 / 3.0).acos())
    } else if shear <= 0.0 {
        0.0
    } else {
        let k = (normal - (normal * normal + 8.0 * shear * shear).sqrt()) / (4.0 * shear);
        (2.0 * k.atan()).abs().min(misalign)
    };
    let mut way = dot(cross(prev, toward), p).signum();
    if misalign > KINK_AMBIGUOUS_DEG.to_radians() && *sense != 0.0 {
        way = *sense;
    }
    if turn > 0.0 {
        *sense = way;
    }
    normalize(rotate(prev, p, way * turn))
}

/// The normal stress across a crack running along `prev` and the shear along it, where
/// the field would have a crack run `toward`, under principal stresses `s1` and `s2`: the
/// crack opens, and runs on, while the normal stress exceeds the cracking stress.
pub fn plane_loading(
    prev: [f64; 3],
    toward: [f64; 3],
    p: [f64; 3],
    s1: f64,
    s2: f64,
) -> (f64, f64) {
    let misalign = dot(in_plane(prev, p), toward).clamp(-1.0, 1.0).acos();
    let (sn, cs) = misalign.sin_cos();
    (s1 * cs * cs + s2 * sn * sn, (s1 - s2) * sn * cs)
}

/// A direction laid in the tangent plane at a point.
fn in_plane(v: [f64; 3], p: [f64; 3]) -> [f64; 3] {
    let r = dot(v, p);
    normalize([v[0] - r * p[0], v[1] - r * p[1], v[2] - r * p[2]])
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
    /// The fissure is held open by more than the tide: a plume's fissure by the bending
    /// of its neighbours' flanking ridges, which place it (F57), and by its own venting.
    /// It runs while its zone lets it, whatever the tide across its own plane.
    held_open: bool,
    /// When the crack opened, Gyr, for the shell's stresses besides the tide (F82); none
    /// for a trace of the tide alone, and for a plume's fissure, which the bending of its
    /// neighbours' flanking ridges steers along them far more strongly than the shell's
    /// heterogeneity could (some 135 kPa on Enceladus against 16).
    time: Option<f64>,
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
            held_open: false,
            time: None,
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
    /// The shell's stresses besides the tide (F82): its freezing, once its thermal history
    /// is known (`with_freezing`), and its heterogeneity.
    pub freezing: Option<crate::stages::tidal::Freezing>,
    pub heterogeneity: Option<crate::stages::tidal::Heterogeneity>,
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
        let shell_map = crate::stages::ice::ShellMap::new(graph);
        // A floating shell is stressed as its local thickness has it (F76).
        if graph.bool("interior.has_subsurface_ocean") {
            shell.local = Some(crate::stages::tidal::LocalShell {
                map: shell_map.clone(),
                love: crate::stages::tidal::ShellLove::new(
                    graph.f64("body.radius"),
                    graph.f64("body.gravity"),
                    graph.f64("body.density"),
                ),
                mean: graph.f64("interior.ice_shell_thickness"),
            });
        }
        let amplitude = shell.amplitude();
        if amplitude <= ICE_CRACK_STRESS {
            return None;
        }
        let episodes = crate::stages::ice::focusing_episodes(graph, &shell_map);
        let heterogeneity = Some(crate::stages::tidal::Heterogeneity::new(
            graph.f64("seed.tectonics") as u64,
            graph.f64("body.radius"),
            crate::stages::ice_shell::flexural_parameter(
                graph.f64("interior.elastic_thickness"),
                graph.f64("body.gravity"),
            ),
        ));
        Some(Self {
            shell,
            shell_map,
            episodes,
            freezing: None,
            heterogeneity,
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

    /// The cracks of a shell whose thermal history is known, which freezes and thaws as
    /// it thickens and thins (F82).
    pub fn with_freezing(mut self, thermal: &crate::thermal::ThermalHistory) -> Self {
        self.freezing = crate::stages::tidal::Freezing::new(thermal, self.body_radius);
        self
    }

    /// The isotropic stress the shell's freezing adds at a time, Pa (F82).
    pub fn freezing_at(&self, time: Option<f64>) -> f64 {
        match (&self.freezing, time) {
            (Some(f), Some(t)) => f.at(t),
            _ => 0.0,
        }
    }

    /// The whole stress a crack opening at a time meets at a point and phase: the tide's,
    /// the shell's freezing and its heterogeneity (F82); the tide's alone without a time.
    pub fn stress_at(
        &self,
        axis: [f64; 3],
        p: [f64; 3],
        phase: f64,
        time: Option<f64>,
    ) -> crate::stages::tidal::Stress {
        let (field, heterogeneity) = self.stress_parts(axis, p, phase, time);
        field + heterogeneity
    }

    /// The stress at a point and phase in its two parts: the field, the tide's and the
    /// shell's freezing, and the heterogeneity (F82), which a crack feels along its length.
    fn stress_parts(
        &self,
        axis: [f64; 3],
        p: [f64; 3],
        phase: f64,
        time: Option<f64>,
    ) -> (crate::stages::tidal::Stress, crate::stages::tidal::Stress) {
        let tide = self.shell.total(axis, p, phase);
        let none = crate::stages::tidal::Stress {
            ee: 0.0,
            nn: 0.0,
            en: 0.0,
        };
        let Some(t) = time else {
            return (tide, none);
        };
        let iso = self.freezing_at(time);
        let field = tide
            + crate::stages::tidal::Stress {
                ee: iso,
                nn: iso,
                en: 0.0,
            };
        let heterogeneity = match &self.heterogeneity {
            Some(h) => {
                let (s1, s2, _) = field.principal();
                h.at(p, t, s1.abs().max(s2.abs()))
            }
            None => none,
        };
        (field, heterogeneity)
    }

    /// The heterogeneity a crack's tip feels, `felt` so far, as the crack grows by `step`
    /// to `run` long, meeting `local` at its tip (F82): the tip's stress intensity weighs
    /// the stress along the crack by the inverse square root of the distance back from
    /// the tip, whose mean is a third of the crack's length, so the crack answers a
    /// heterogeneity finer than itself less than a coarser one.
    fn felt(
        felt: Option<crate::stages::tidal::Stress>,
        local: crate::stages::tidal::Stress,
        step: f64,
        run: f64,
    ) -> crate::stages::tidal::Stress {
        match felt {
            None => local,
            Some(f) => f + (local - f).scaled((step / (run / 3.0).max(step)).min(1.0)),
        }
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
                    if peak + self.freezing_at(Some(time)) > ICE_CRACK_STRESS {
                        start = Some((tip, phase));
                    }
                } else {
                    let freezing = self.freezing_at(Some(time));
                    for attempt in 0..24u64 {
                        let p = rng.unit_direction(bin, draw + 1 + 2 * attempt);
                        let (peak, phase) = self.shell.peak_tension(axis, p);
                        let excess = peak + freezing - ICE_CRACK_STRESS;
                        if excess <= 0.0 {
                            continue;
                        }
                        let thin = thinnest / self.shell_map.factor(p).max(1e-3);
                        if rng.uniform(bin ^ crate::rng::ACCEPT_STREAM, draw + attempt)
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
                    let (t1, _, _) = self.stress_at(axis, p0, phase, Some(time)).principal();
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
                                BAND_HEIGHT * (0.6 + 0.8 * rng.uniform(bin, draw + 53)),
                            )
                        }
                        LineaKind::Vent => (
                            VENT_WIDTH,
                            VENT_TROUGH_DEPTH * (0.7 + 0.6 * rng.uniform(bin, draw + 53)),
                        ),
                        _ => {
                            // A heavy tail of widths (F78): the inverse of the cumulative
                            // share of a truncated power law.
                            let (w_lo, w_hi) = LINEA_WIDTH_RANGE;
                            let k = 1.0 - LINEA_WIDTH_EXPONENT;
                            let u = rng.uniform(bin, draw + 52);
                            let w =
                                (w_lo.powf(k) + u * (w_hi.powf(k) - w_lo.powf(k))).powf(1.0 / k);
                            (
                                w,
                                LINEA_RELIEF_SHARE * w * (0.6 + 0.8 * rng.uniform(bin, draw + 53)),
                            )
                        }
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
                        held_open: false,
                        time: Some(time),
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
                held_open: true,
                time: None,
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
                // A fissure of the present episode vents to the present: its trough is
                // renewed by its own eruptions until now, the last surface of its zone,
                // which no lesser crack of the zone crosses.
                None => present,
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
                        held_open: true,
                        time: None,
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
                        held_open: true,
                        time: None,
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

    /// A narrow crack of the ridged plains' fabric (F81), traced as a ridge is but
    /// answering no other crack, in the stresses of the time it opened (F82).
    #[allow(clippy::too_many_arguments)]
    pub fn trace_at(
        &self,
        axis: [f64; 3],
        p0: [f64; 3],
        phase: f64,
        spacing: f64,
        length: f64,
        id: u64,
        time: f64,
    ) -> Vec<[f64; 3]> {
        let rules = TraceRules {
            time: Some(time),
            ..TraceRules::free()
        };
        self.ridge_trace(axis, p0, phase, spacing, length, id, &rules)
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
            let mut kink_sense = 0.0;
            let mut felt = None;
            for i in 0..half {
                let (field, local) = self.stress_parts(axis, p, phase, rules.time);
                // Both tips grow together: the crack is twice as long as either side.
                let tip = Self::felt(felt, local, spacing, 2.0 * i as f64 * spacing);
                felt = Some(tip);
                let stress = field + tip;
                let (s1, s2, _) = stress.principal();
                let mut d = crack_direction(p, &stress, along);
                // The crack runs on while the stress across its own plane opens it.
                let opening = along.map_or(s1, |prev| plane_loading(prev, d, p, s1, s2).0);
                if !rules.held_open && opening < ICE_CRACK_STRESS && i > 2 {
                    break;
                }
                if along.is_none() {
                    d = [sign * d[0], sign * d[1], sign * d[2]];
                }
                // Wander: a slow turn by noise along the trace; and for a branch the
                // bias it leaves its parent with, fading as it bends back to the field.
                // Both turn the heading the crack seeks, not the one it has, so they do
                // not add up step by step.
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
                let sought = normalize([
                    d[0] * turn.cos() + side_dir[0] * turn.sin(),
                    d[1] * turn.cos() + side_dir[1] * turn.sin(),
                    d[2] * turn.cos() + side_dir[2] * turn.sin(),
                ]);
                // A running crack turns toward the heading it seeks only as the field
                // shears it (`kinked`), so it does not turn at right angles where the
                // field's principal directions swap.
                let mut d = match along {
                    Some(prev) => kinked(prev, sought, p, s1, s2, &mut kink_sense),
                    None => sought,
                };
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
        let mut active = 0usize;
        let mut kink_sense = 0.0;
        let max_steps = CYCLOID_MAX_ORBITS * CYCLOID_PHASE_STEPS;
        let mut felt = None;
        for _ in 0..max_steps {
            let (field, local) = self.stress_parts(axis, p, phase, rules.time);
            let tip = Self::felt(felt, local, advance, run);
            felt = Some(tip);
            let stress = field + tip;
            let (s1, s2, _) = stress.principal();
            let toward = crack_direction(p, &stress, along);
            // The crack runs while the stress across its own plane opens it (Hoppa et
            // al.'s timing, judged on the crack's plane): where the field's principal
            // directions swap through the orbit, as they do on the equator, it waits for
            // the next orbit instead of turning.
            let opening = along.map_or(s1, |prev| plane_loading(prev, toward, p, s1, s2).0);
            if opening >= ICE_CRACK_STRESS {
                idle = 0;
                // A crack the tension holds open through half an orbit no longer stops and
                // starts again with the orbit, as a cycloid's does (F82: a shell's freezing
                // can hold it so): it would follow the field's turn into a loop, and ends.
                active += 1;
                if active > CYCLOID_PHASE_STEPS / 2 {
                    break;
                }
                let d = match along {
                    Some(prev) => kinked(prev, toward, p, s1, s2, &mut kink_sense),
                    None => toward,
                };
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
                active = 0;
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
            // Degraded or drawn wider than it is, a ridge keeps the area of its
            // cross-section: its height falls as its width grows.
            let w = width.hypot(2.0 * blur).max(1.0);
            let relief = relief * width / w;
            let crest = 0.25 * w;
            let sigma = (w / 6.0).max(1.0);
            let ridges = relief * (-((x - crest) / sigma).powi(2)).exp();
            let trough = LINEA_TROUGH_FRACTION * relief * (-(x / sigma).powi(2)).exp();
            ridges - trough
        }
    }
}

/// The relief of a band across it, m, at `x` from its axis, `half` its half-width, standing
/// `height` high (F83): a plateau of buoyant ice with a V-shaped trough along the spreading
/// axis down to the old surface's level, a ridge along each margin, and soft walls; its
/// lineations parallel to the axis where the grid resolves them; `blur` the spreading its
/// relaxation since has brought.
fn band_profile(x: f64, half: f64, height: f64, blur: f64, lineation: f64) -> f64 {
    let a = x.abs();
    let wall = (0.1 * half).hypot(blur).max(1.0);
    let plateau = 0.5 * (1.0 - ((a - half) / wall).tanh());
    let trough_half = (BAND_TROUGH_SHARE * half).hypot(blur).max(1.0);
    let trough = (1.0 - a / trough_half).max(0.0);
    let margin = (-((a - 0.9 * half) / (0.06 * half).hypot(blur).max(1.0)).powi(2)).exp();
    let lines = if a < half {
        BAND_LINEATION
            * (1.0 - (a / half).powi(2))
            * (2.0 * std::f64::consts::PI * a / lineation).sin()
    } else {
        0.0
    };
    height * (plateau * (1.0 - trough) + BAND_MARGIN_RIDGE * margin) + lines
}

/// A band on the resolved grid (F83). The crack pulled apart, and buoyant ice spread from
/// its axis between the halves: inside, the band's own surface, the axis youngest, standing
/// above the old surface's level at the crack; outside, the plates it opened between,
/// pushed apart by its half-width at each margin, the push falling off to nothing at half
/// the band's length, so that what crossed the crack before it opened is parted by it.
fn band_compute(
    layers: &LayerSet,
    grid: &CubeSphere,
    body_radius: f64,
    rec: &Linea,
    blur: f64,
    present: f64,
    stamper: &mut Stamper,
) -> StampOutput {
    let cell = grid.delta() * body_radius;
    let corridor = rec.reach(cell);
    let n = rec.points.len();
    if n < 2 {
        return StampOutput::default();
    }
    stamper.gather_corridor(grid, body_radius, &rec.points, corridor);
    let footprint = std::mem::take(&mut stamper.footprint);
    let ages = layers.surface_age.as_slice();
    let units = layers.unit_id.as_slice();
    let deeps = layers.deep.as_slice();
    let melts = layers.melt.as_slice();
    let covers = layers.cover.as_slice();
    let bedrock = PaddedField::new(grid, layers.bedrock.as_slice());
    let stain = PaddedField::new(grid, layers.stain.as_slice());
    let lineation = BAND_LINEATION_SPACING.max(3.0 * cell);
    let mut cells: Vec<(CellWrite, [f32; 2])> = footprint
        .iter()
        .filter_map(|&(idx, _, along)| {
            let d = grid.center(grid.cell(idx));
            let k = (along.floor().max(0.0) as usize).min(n - 2);
            let s = (along / (n - 1) as f64).clamp(0.0, 1.0);
            let (a, b) = (rec.points[k], rec.points[k + 1]);
            let normal = normalize(cross(a, b));
            // The signed distance across the crack, m, and the band's half-width here: a
            // lens, widest in the middle and closing to its tips.
            let x = dot(d, normal) * body_radius;
            let half = 0.5 * rec.width * (4.0 * s * (1.0 - s)).max(0.0).powf(BAND_TAPER_POWER);
            let side = x.signum();
            // The share of the cell the band covers, across its edge.
            let covered = if half > 0.5 * cell {
                (0.5 + (half - x.abs()) / cell).clamp(0.0, 1.0)
            } else {
                0.0
            };
            let push = half
                * (1.0 - (x.abs() - half).max(0.0) / (corridor - half).max(1.0)).clamp(0.0, 1.0);
            if covered <= 0.0 && push < 0.05 * cell {
                return None;
            }
            // The plate, pushed away from the axis: the old surface from nearer the crack.
            let source = normalize([
                d[0] - normal[0] * side * push / body_radius,
                d[1] - normal[1] * side * push / body_radius,
                d[2] - normal[2] * side * push / body_radius,
            ]);
            let from = grid.index(grid.locate(source).0);
            let plate_bedrock = sample_bilinear(grid, &bedrock, source) as f64;
            let plate_stain = sample_bilinear(grid, &stain, source) as f64;
            let mut write = CellWrite {
                idx: idx as u32,
                bedrock: plate_bedrock as f32,
                surface: true,
                age: ages[from],
                unit: units[from],
                deep: deeps[from],
                melt: melts[from],
                cover: covers[from],
            };
            if covered <= 0.0 {
                return Some((write, [plate_stain as f32, 1.0]));
            }
            // The band's own ice, at the old surface's level where the crack was, its age
            // spreading from the axis.
            let foot = normalize([
                d[0] - normal[0] * x / body_radius,
                d[1] - normal[1] * x / body_radius,
                d[2] - normal[2] * x / body_radius,
            ]);
            let base = sample_bilinear(grid, &bedrock, foot) as f64;
            let time = rec.time_gyr - BAND_SPREAD_GYR * (x.abs() / half.max(1.0)).min(1.0);
            let band_bedrock = base + band_profile(x, half, rec.relief, blur, lineation);
            let band_stain =
                BAND_STAIN_SHARE * crate::stages::cryo::matrix_stain((present - time).max(0.0));
            write.bedrock = (covered * band_bedrock + (1.0 - covered) * plate_bedrock) as f32;
            if covered >= 0.5 {
                write.age = time as f32;
                write.unit = UNIT_BAND;
                write.deep = 0.0;
                write.melt = 0.0;
                write.cover = 0.0;
            }
            let laid = covered * band_stain + (1.0 - covered) * plate_stain;
            Some((write, [laid as f32, 1.0]))
        })
        .collect();
    stamper.footprint = footprint;
    cells.sort_by_key(|c| c.0.idx);
    let (cells, stain): (Vec<CellWrite>, Vec<[f32; 2]>) = cells.into_iter().unzip();
    StampOutput {
        cells,
        composition: None,
        stain: Some(stain),
        steep: None,
        ejecta_age: None,
        compositions: None,
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
    if rec.kind == LineaKind::Band {
        return band_compute(layers, grid, body_radius, rec, blur, present, stamper);
    }
    let bedrock = layers.bedrock.as_slice();
    let ages = layers.surface_age.as_slice();
    let units = layers.unit_id.as_slice();
    let deeps = layers.deep.as_slice();
    let melts = layers.melt.as_slice();
    let covers = layers.cover.as_slice();
    let cell = grid.delta() * body_radius;
    let width = rec.width.max(2.0 * cell);
    // A ridge narrower than the grid draws is drawn at the grid's width by the area of
    // its cross-section, and makes no surface of its own (F78).
    let ridged = matches!(rec.kind, LineaKind::Ridge | LineaKind::Cycloid);
    let area_share = if ridged {
        (rec.width / width).min(1.0)
    } else {
        1.0
    };
    let resolved = !ridged || rec.width >= LINEA_RESOLVED_CELLS * cell;
    // Its margins' stain as it stands at the present (F78): laid after the crack opens,
    // faded since, to the ancient band's.
    let age = (present - rec.time_gyr).max(0.0);
    let stain_now = if rec.kind == LineaKind::Vent {
        0.0
    } else {
        LINEA_STAIN
            * (1.0 - (-age / LINEA_MARGIN_ONSET_GYR).exp())
            * (LINEA_STAIN_RESIDUAL + (1.0 - LINEA_STAIN_RESIDUAL) * (-age / SALT_FADE_GYR).exp())
    };
    let margin = (LINEA_MARGIN_SHARE * rec.width).hypot(0.5 * cell);
    // A ridge or a cycloid is a chain of segments (F62); the stamp reaches as its cap
    // says.
    let rng = Rng::new(rec.id ^ 0x11CE_5E61);
    let chained = !matches!(rec.kind, LineaKind::Band | LineaKind::Vent);
    let (bounds, asides) = chain(rec, width);
    let segments = asides.len();
    let reach = rec.reach(cell);
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
    let mut cells: Vec<(CellWrite, [f32; 2])> = best
        .into_iter()
        .map(|(idx, (_, across, s))| {
            // The segment of the chain the cell lies along: how far it stands aside, and
            // how open it is there, an ellipse along the segment and along the crack.
            let segment = (0..segments)
                .find(|&i| s < bounds[i + 1])
                .unwrap_or(segments - 1);
            let u = ((s - bounds[segment]) / (bounds[segment + 1] - bounds[segment]).max(1e-9))
                .clamp(0.0, 1.0);
            // The margins' stain follows the crack's trace, a deposit spread over the
            // whole lineament; the relief follows its segments.
            let from_trace = across;
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
            // A crack closes toward its tips: an elastic crack under a uniform load opens
            // as an ellipse, so a ridge's or a cycloid's width, its stain and the surface
            // it renews follow the crack's opening, (4 s (1 - s))^(1/2) along it.
            let tip = if chained {
                (4.0 * s * (1.0 - s)).max(0.0).sqrt()
            } else {
                1.0
            };
            let width_here = width_here * vary * tip;
            let relief_here = vary
                * area_share
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
            // Where the grid resolves it: a ridge by its width, a band by its width here,
            // which closes to nothing toward its tips.
            let drawn = resolved
                && (rec.kind != LineaKind::Band || width_here >= LINEA_RESOLVED_CELLS * cell);
            let inside = drawn && across.abs() < 0.5 * width_here && buried < 0.5;
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
            // A ridge narrower than the grid renews the share of the cell its width covers
            // (F81): the cell's surface takes the crack's time by that share.
            let renewed = if ridged && !drawn {
                crate::stages::fabric::strip_share(across, rec.width * vary * open * tip, cell)
                    * (1.0 - buried)
            } else {
                0.0
            };
            if renewed > 0.0 && ages[idx] < time {
                let r = renewed as f32;
                write.surface = true;
                write.age = (1.0 - r) * ages[idx] + r * time;
                write.unit = units[idx];
                write.deep = deeps[idx] * (1.0 - r);
                write.melt = melts[idx] * (1.0 - r);
                write.cover = covers[idx] * (1.0 - r);
            }
            if inside && ages[idx] <= time + 1e-6 {
                write.surface = true;
                write.age = time;
                write.unit = match rec.kind {
                    LineaKind::Band => UNIT_BAND,
                    _ => UNIT_LINEA,
                };
                write.deep = 0.0;
                write.melt = 0.0;
                // A vent's warm ice sinters coarse (the plume's tint); a crack's colour
                // is its stain.
                write.cover = if rec.kind == LineaKind::Vent {
                    1.0
                } else {
                    0.0
                };
            }
            // The stain (F78): a ridge's margins, a Gaussian about its axis dimmer over
            // its medial stripe, cutting away the older stain under its own width; a
            // band's floor, stained across it. A crack buried under loose cover laid it
            // under the snow.
            let core = (0.5 * width_here).hypot(0.5 * cell);
            let cut = if ridged && !drawn {
                renewed
            } else {
                (1.0 - (from_trace / core).powi(2)).clamp(0.0, 1.0)
            };
            let laid = match rec.kind {
                LineaKind::Vent => 0.0,
                LineaKind::Band => stain_now * cut,
                // A ridge narrower than the grid lays its stain by its share (F81).
                _ if ridged && !drawn => {
                    stain_now * crate::stages::fabric::narrow_margin(from_trace, rec.width, cell)
                }
                _ => {
                    stain_now
                        * (-(from_trace / margin).powi(2)).exp()
                        * (1.0 - LINEA_STRIPE_DIM * (-(from_trace / (0.6 * core)).powi(2)).exp())
                }
            } * (1.0 - buried)
                * tip;
            (write, [laid as f32, (cut * tip) as f32])
        })
        .collect();
    cells.sort_by_key(|c| c.0.idx);
    let stain = (rec.kind != LineaKind::Vent).then(|| cells.iter().map(|c| c.1).collect());
    StampOutput {
        cells: cells.into_iter().map(|c| c.0).collect(),
        composition: None,
        stain,
        steep: None,
        ejecta_age: None,
        compositions: None,
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
        // A band (F83): it stands high, a trough along its axis down to the old level and
        // a ridge along each margin, and nothing beyond it.
        let at = |x: f64| band_profile(x, 10e3, 150.0, 0.0, 1e9);
        assert!(at(0.0).abs() < 1.0, "{}", at(0.0));
        assert!(at(5e3) > 140.0, "{}", at(5e3));
        assert!(at(9e3) > at(6e3) + 20.0, "{} {}", at(9e3), at(6e3));
        assert!(at(15e3).abs() < 5.0, "{}", at(15e3));
    }

    /// A band pulls the plates it opens between apart (F83): a bump beside it is pushed
    /// outward by about its half-width, less as the push spreads into the plate, and the
    /// band's ice spreads from its axis, the axis youngest.
    #[test]
    fn a_band_pulls_the_plates_apart() {
        let r = 300e3;
        let grid = CubeSphere::new(256);
        let spacing = 1e3;
        let n = (300e3 / spacing) as usize;
        let points: Vec<[f64; 3]> = (0..=n)
            .map(|i| {
                let a = (i as f64 * spacing - 150e3) / r;
                [a.cos(), a.sin(), 0.0]
            })
            .collect();
        let rec = Linea {
            points,
            spacing,
            length: 300e3,
            width: 20e3,
            relief: 150.0,
            kind: LineaKind::Band,
            phase: 0.0,
            time_gyr: 1.0,
            id: 5,
        };
        let mut layers = LayerSet::new(grid);
        let bump = crate::grid::direction_from_lon_lat(0.0, 15e3 / r);
        for k in 0..grid.cell_count() {
            let d = grid.center(grid.cell(k));
            let dist = angle_between(d, bump) * r;
            layers.bedrock.as_mut_slice()[k] = (500.0 * (-(dist / 3e3).powi(2)).exp()) as f32;
        }
        let out = crate::crater::with_stamper(&grid, |stamper| {
            linea_compute(&layers, &grid, r, &rec, 0.0, None, 1.0, stamper)
        });
        // Where the bump stands now: the highest cell north of the band near its meridian.
        let (mut peak, mut lat) = (f32::MIN, 0.0);
        let (mut axis, mut margin) = (0.0f32, 1.0f32);
        for w in &out.cells {
            let d = grid.center(grid.cell(w.idx as usize));
            let (lon, la) = crate::grid::lon_lat_from_direction(d);
            let north = la * r;
            if lon.abs() * r < 3e3 && north > 12e3 && w.bedrock > peak {
                peak = w.bedrock;
                lat = north;
            }
            if w.unit == UNIT_BAND && lon.abs() * r < 20e3 {
                if north.abs() < 2e3 {
                    axis = axis.max(w.age);
                } else if north.abs() > 7e3 {
                    margin = margin.min(w.age);
                }
            }
        }
        assert!((22e3..26e3).contains(&lat), "bump at {:.1} km", lat / 1e3);
        assert!(peak > 400.0, "{peak}");
        assert!(axis > margin && axis <= 1.0, "axis {axis}, margin {margin}");
    }

    /// A crack's heading is carried from step to step in the plane tangent where it is:
    /// Europa's ridges and cycloids almost never turn back and forth step after step (a
    /// heading left tilted out of the plane zigzagged at every step). A cycloid's cusp
    /// reverses its curvature once; a zigzag reverses it at every step.
    #[test]
    fn europa_s_cracks_do_not_zigzag() {
        let (g, _) = crate::presets::preset("Europa").unwrap().to_graph();
        let cracks = Cracks::new(&g).unwrap();
        let present = g.f64("body.age");
        let all = cracks.generate(present, &CubeSphere::new(256), 3);
        // Per kind, ridges then cycloids: runs of three turns, those alternating in sense
        // by over half a degree each, and the worst of them.
        let mut tally = [(0usize, 0usize, 0.0f64); 2];
        for l in all.iter().filter(|l| l.kind != LineaKind::Band) {
            let (steps, reversals, worst) = &mut tally[usize::from(l.kind == LineaKind::Cycloid)];
            let turns: Vec<f64> = l
                .points
                .windows(3)
                .map(|w| {
                    let a = normalize([w[1][0] - w[0][0], w[1][1] - w[0][1], w[1][2] - w[0][2]]);
                    let b = normalize([w[2][0] - w[1][0], w[2][1] - w[1][1], w[2][2] - w[1][2]]);
                    dot(cross(a, b), w[1]).signum() * dot(a, b).clamp(-1.0, 1.0).acos().to_degrees()
                })
                .collect();
            for t in turns.windows(3) {
                *steps += 1;
                if t[0] * t[1] < 0.0 && t[1] * t[2] < 0.0 && t.iter().all(|x| x.abs() > 0.5) {
                    *reversals += 1;
                    *worst = worst.max(t.iter().fold(f64::INFINITY, |m, x| m.min(x.abs())));
                }
            }
        }
        let [(rs, rr, rw), (cs, cr, cw)] = tally;
        eprintln!("ridges {rr} of {rs}, worst {rw:.1}; cycloids {cr} of {cs}, worst {cw:.1}");
        assert!(
            (rr as f64) < 0.005 * rs as f64 && rw < 5.0,
            "ridges {rr} of {rs}, {rw:.1}"
        );
        // A cycloid's arcs curve steadily and its cusps turn sharply one way; the few
        // alternations sit at the cusps, small.
        assert!(
            (cr as f64) < 0.06 * cs as f64 && cw < 8.0,
            "cycloids {cr} of {cs}, {cw:.1}"
        );
    }

    /// The shell's heterogeneity (F82) bends Europa's cracks: traced from the same starts
    /// in the stresses of their time, they turn more along their length than in the
    /// tide's alone, and none turns back on itself.
    #[test]
    fn the_shell_s_heterogeneity_bends_the_cracks() {
        let (g, _) = crate::presets::preset("Europa").unwrap().to_graph();
        let cracks = Cracks::new(&g).unwrap();
        let present = g.f64("body.age");
        let r = cracks.body_radius;
        let rng = Rng::new(17);
        // Turning per 100 km, in degrees, of a trace.
        let turning = |pts: &[[f64; 3]]| {
            let (mut turn, mut length) = (0.0, 0.0);
            for w in pts.windows(3) {
                let a = normalize([w[1][0] - w[0][0], w[1][1] - w[0][1], w[1][2] - w[0][2]]);
                let b = normalize([w[2][0] - w[1][0], w[2][1] - w[1][1], w[2][2] - w[1][2]]);
                turn += dot(a, b).clamp(-1.0, 1.0).acos().to_degrees();
                length += angle_between(w[0], w[1]) * r;
            }
            100e3 * turn / length.max(1.0)
        };
        let (mut tide, mut both, mut n) = (0.0, 0.0, 0);
        for j in 0..200u64 {
            let time = present - 0.1 * rng.uniform(j, 0);
            let axis = mean_axis_at(time, present, cracks.nsr_period_gyr);
            let p0 = rng.unit_direction(j, 1);
            let (peak, phase) = cracks.shell.peak_tension(axis, p0);
            if peak <= ICE_CRACK_STRESS {
                continue;
            }
            let plain = cracks.trace(axis, p0, phase, 5e3, 800e3, j);
            let bent = cracks.trace_at(axis, p0, phase, 5e3, 800e3, j, time);
            if plain.len() < 20 || bent.len() < 20 {
                continue;
            }
            tide += turning(&plain);
            both += turning(&bent);
            n += 1;
        }
        let (tide, both) = (tide / n as f64, both / n as f64);
        eprintln!(
            "{n} cracks turn {tide:.1} degrees per 100 km in the tide, {both:.1} with the shell's heterogeneity"
        );
        assert!(n > 50, "{n}");
        assert!(both > 1.5 * tide, "{both:.1} against {tide:.1}");
    }

    /// A crack runs on straight where the tension is near equal in all directions, so it
    /// does not turn at right angles where the field's principal directions swap: traced
    /// from the equator, where the field's symmetry sets them east and north, no step of
    /// Europa's cracks turns by as much as 45 degrees.
    #[test]
    fn a_crack_does_not_stair_step_where_the_field_swaps() {
        let (g, _) = crate::presets::preset("Europa").unwrap().to_graph();
        let cracks = Cracks::new(&g).unwrap();
        let present = g.f64("body.age");
        let rules = TraceRules {
            polar: false,
            index: None,
            exclude: None,
            one_way: None,
            junction: 1e3,
            wobble: (0.0, 1.0),
            bias: (0.0, 1.0),
            held_open: false,
            time: None,
        };
        let mut sharpest = 0.0f64;
        let mut sharpest_cycloid = 0.0f64;
        let turn = |pts: &[[f64; 3]]| {
            let mut most = 0.0f64;
            for w in pts.windows(3) {
                let a = normalize([w[1][0] - w[0][0], w[1][1] - w[0][1], w[1][2] - w[0][2]]);
                let b = normalize([w[2][0] - w[1][0], w[2][1] - w[1][1], w[2][2] - w[1][2]]);
                most = most.max(dot(a, b).clamp(-1.0, 1.0).acos().to_degrees());
            }
            most
        };
        for k in 0..72 {
            let lon = (k as f64 * 5.0).to_radians();
            for lat in [-0.03f64, 0.0, 0.04] {
                let p0 = [lat.cos() * lon.cos(), lat.cos() * lon.sin(), lat.sin()];
                for time in [present - 0.02, present - 0.06] {
                    let axis = mean_axis_at(time, present, cracks.nsr_period_gyr);
                    let (peak, phase) = cracks.shell.peak_tension(axis, p0);
                    if peak <= ICE_CRACK_STRESS {
                        continue;
                    }
                    let pts = cracks.ridge_trace(axis, p0, phase, 1.2e3, 600e3, k, &rules);
                    sharpest = sharpest.max(turn(&pts));
                    // A cycloid kinks at its cusps, by no more than the criterion's
                    // 70.5 degrees: never the right angles of a stair.
                    let pts = cracks.cycloid(axis, p0, phase, 1.2e3, 600e3, &rules);
                    sharpest_cycloid = sharpest_cycloid.max(turn(&pts));
                }
            }
        }
        assert!(sharpest < 45.0, "sharpest step {sharpest:.1} degrees");
        assert!(
            sharpest_cycloid < 75.0,
            "sharpest cycloid step {sharpest_cycloid:.1} degrees"
        );
    }

    /// Europa's cracks nucleate over every longitude: no sector of the body is short of
    /// them, as a test drawn from the stream that placed the candidate would make it.
    #[test]
    fn europa_s_cracks_reach_every_longitude() {
        let (g, _) = crate::presets::preset("Europa").unwrap().to_graph();
        let cracks = Cracks::new(&g).unwrap();
        let present = g.f64("body.age");
        let all = cracks.generate(present, &CubeSphere::new(64), 7);
        let mut sectors = [0usize; 6];
        for l in &all {
            let m = l.midpoint();
            let lon = m[1].atan2(m[0]).to_degrees() + 180.0;
            sectors[((lon / 60.0) as usize).min(5)] += 1;
        }
        let mean = all.len() as f64 / 6.0;
        assert!(
            sectors.iter().all(|&n| n as f64 > 0.5 * mean),
            "{sectors:?} of {}",
            all.len()
        );
    }

    /// Under a plume the vents of the thinned zone come in a cascade (F57): a parallel
    /// set a quarter flexural parameter apart, each long enough to cross the core, and
    /// each running its whole length.
    #[test]
    fn the_tiger_stripes_form_a_parallel_cascade() {
        let (g, _) = crate::presets::preset("Enceladus").unwrap().to_graph();
        let cracks = Cracks::new(&g).unwrap();
        let present = g.f64("body.age");
        // As the preset renders: its own seed, at the standard finish's grid.
        let seed = g.f64("seed.cratering") as u64;
        let all = cracks.generate(present, &CubeSphere::new(FINISHING_FACE_SIZE[1]), seed);
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
        // Held open by their neighbours' bending, the stripes run their whole length: a
        // stop by the tide across their own plane cut half of them short.
        let traced = |l: &Linea| {
            l.points
                .windows(2)
                .map(|w| angle_between(w[0], w[1]))
                .sum::<f64>()
                * r
        };
        let whole: Vec<&&Linea> = stripes
            .iter()
            .filter(|l| l.length > 150e3 && traced(l) > 0.9 * l.length)
            .collect();
        assert!(whole.len() >= 5, "{} stripes run their length", whole.len());
        // And they keep parallel along their length, steered by their neighbours'
        // bending: each within 50 degrees of the first at each tenth of its length (the
        // shell's heterogeneity, had it bent them, let two converge at 63).
        let tangent = |l: &Linea, k: usize| {
            let i = (l.points.len() * k / 10).clamp(1, l.points.len() - 2);
            let (a, b) = (l.points[i - 1], l.points[i + 1]);
            normalize([b[0] - a[0], b[1] - a[1], b[2] - a[2]])
        };
        let worst = whole
            .iter()
            .flat_map(|l| (1..10).map(move |k| dot(tangent(l, k), tangent(first, k)).abs()))
            .fold(1.0f64, f64::min);
        assert!(
            worst > 0.65,
            "stripes as far apart in heading as {worst:.2}"
        );
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

    /// A ridge's margins (F78): a red-brown stain about twice its width across, dimmer
    /// over its medial stripe, absent while the crack is fresh and faded to the ancient
    /// band's when it is old; on a grid coarser than the ridge, its relief keeps the area
    /// of its cross-section and it renews its share of the cells it crosses.
    #[test]
    fn a_ridge_stains_its_margins_as_it_ages() {
        // A grid of 2 km cells and a ridge 3 km wide along the equator.
        let grid = CubeSphere::new(256);
        let r = 2e3 / grid.delta();
        let spacing = 1e3;
        let n = (400e3 / spacing) as usize;
        let points: Vec<[f64; 3]> = (0..=n)
            .map(|i| {
                let a = (i as f64 * spacing - 200e3) / r;
                [a.cos(), a.sin(), 0.0]
            })
            .collect();
        let stamp = |present: f64| {
            let rec = Linea {
                points: points.clone(),
                spacing,
                length: 400e3,
                width: 3e3,
                relief: 180.0,
                kind: LineaKind::Ridge,
                phase: 0.0,
                time_gyr: 1.0,
                id: 9,
            };
            let layers = LayerSet::new(grid);
            crate::crater::with_stamper(&grid, |stamper| {
                linea_compute(&layers, &grid, r, &rec, 0.0, None, present, stamper)
            })
        };
        // The stain across the middle of the ridge, by distance from its axis in 2 km
        // bins, a cell each.
        let across = |out: &StampOutput| {
            let mut by_km = vec![0.0f32; 12];
            let stain = out.stain.as_ref().expect("a ridge lays a stain");
            for (c, s) in out.cells.iter().zip(stain) {
                let d = grid.center(grid.cell(c.idx as usize));
                if (d[1].atan2(d[0]) * r).abs() < 50e3 {
                    let x = (d[2].asin() * r / 2e3).abs() as usize;
                    if x < 12 {
                        by_km[x] = by_km[x].max(s[0]);
                    }
                }
            }
            by_km
        };
        let fresh = across(&stamp(1.0 + 1e-5));
        let young = across(&stamp(1.003));
        let old = across(&stamp(1.2));
        assert!(fresh.iter().all(|&v| v < 0.05), "{fresh:?}");
        // A triple band: margins some 9 km across on a 3 km ridge (Geissler et al.: 8 to
        // 12 km on 3 to 6 km ridges).
        assert!(
            young[0] > 0.3 && young[1] > 0.3 && young[2] > 0.1,
            "{young:?}"
        );
        assert!(young[3] < 0.1 && young[5] < 0.02, "{young:?}");
        // The ancient band: faded to a sixth or so of the triple band's.
        let ratio = old[1] / young[1];
        assert!((0.1..0.3).contains(&ratio), "{ratio}");
        // A 3 km ridge on 2 km cells: drawn by its area, it renews the share of each
        // cell its width covers, the cell taking its time by that share, and makes no
        // unit of its own (F81).
        let out = stamp(1.003);
        assert!(out.cells.iter().all(|c| c.unit == UNIT_NONE));
        let most = out
            .cells
            .iter()
            .filter(|c| c.surface)
            .map(|c| c.age)
            .fold(0.0f32, f32::max);
        assert!(most > 0.5 && most < 1.0, "{most}");
        let top = out.cells.iter().map(|c| c.bedrock).fold(0.0f32, f32::max);
        assert!(top > 50.0 && top < 160.0, "{top}");
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
