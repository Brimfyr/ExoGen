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
//! the ocean that fade under frost as they age.

use std::collections::HashMap;

use crate::constants::*;
use crate::crater::{CellWrite, StampOutput, Stamper};
use crate::fields::LayerSet;
use crate::geom::{cross, dot, local_basis, normalize};
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

/// The population of cracks on the present surface: those within the retention age,
/// drawn per time bin.
pub struct Cracks {
    pub shell: TidalShell,
    /// The field's largest peak tension, Pa.
    pub amplitude: f64,
    pub nsr_period_gyr: f64,
    pub density: f64,
    pub age: f64,
    pub retention: f64,
    pub body_radius: f64,
    pub orbit_period_s: f64,
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
        Some(Self {
            shell,
            amplitude,
            nsr_period_gyr,
            density,
            age: graph.f64("body.age"),
            retention: graph.f64("cratering.retention_age"),
            body_radius: graph.f64("body.radius"),
            orbit_period_s: graph.f64("orbit.period"),
        })
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
        let (l_lo, l_hi) = LINEA_LENGTH;
        let mean_length = (l_hi - l_lo) / (l_hi / l_lo).ln();
        let expected = self.total_length() / mean_length;
        let window_start = (present - self.retention).max(0.0);
        let mut out = Vec::new();
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
                // Where: in proportion to the tension above the threshold, by rejection.
                let scale = self.amplitude - ICE_CRACK_STRESS;
                let mut start = None;
                for attempt in 0..24u64 {
                    let p = rng.unit_direction(bin, draw + 1 + 2 * attempt);
                    let (peak, phase) = self.shell.peak_tension(axis, p);
                    let excess = peak - ICE_CRACK_STRESS;
                    if excess <= 0.0 {
                        continue;
                    }
                    if rng.uniform(bin, draw + 2 + 2 * attempt) < excess / scale.max(excess) {
                        start = Some((p, phase));
                        break;
                    }
                }
                let Some((p0, phase)) = start else { continue };
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
                let id = (bin << 32) | j;
                let points = match kind {
                    LineaKind::Cycloid => self.cycloid(axis, p0, phase, spacing, length, id),
                    _ => self.ridge_trace(axis, p0, phase, spacing, length, id),
                };
                if points.len() < 3 {
                    continue;
                }
                let (width, relief) = match kind {
                    LineaKind::Band => {
                        let (w_lo, w_hi) = BAND_WIDTH;
                        (
                            w_lo * (w_hi / w_lo).powf(rng.uniform(bin, draw + 52)),
                            BAND_DEPTH * (0.6 + 0.8 * rng.uniform(bin, draw + 53)),
                        )
                    }
                    _ => (
                        LINEA_WIDTH,
                        LINEA_RIDGE_HEIGHT * (0.6 + 0.8 * rng.uniform(bin, draw + 53)),
                    ),
                };
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
        out
    }

    /// A crack along the slowly turning field: from the start in both directions,
    /// perpendicular to the most tensile stress at the phase it opened, while the
    /// tension holds, with a little wander.
    fn ridge_trace(
        &self,
        axis: [f64; 3],
        p0: [f64; 3],
        phase: f64,
        spacing: f64,
        length: f64,
        id: u64,
    ) -> Vec<[f64; 3]> {
        let noise = GradientNoise::new(id ^ 0x7E1A_0000);
        let half = (0.5 * length / spacing).ceil() as usize;
        let mut sides: [Vec<[f64; 3]>; 2] = [Vec::new(), Vec::new()];
        for (side, sign) in [1.0f64, -1.0].iter().enumerate() {
            let mut p = p0;
            let mut along: Option<[f64; 3]> = None;
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
                // Wander: a slow turn by noise along the trace.
                let turn = LINEA_WOBBLE_RAD
                    * noise.fbm(
                        [i as f64 * spacing / LINEA_WOBBLE_LENGTH, 0.3, 0.7],
                        1.0,
                        2,
                        2.0,
                        0.5,
                    )
                    / 0.6;
                let side_dir = normalize(cross(p, d));
                let d = normalize([
                    d[0] * turn.cos() + side_dir[0] * turn.sin(),
                    d[1] * turn.cos() + side_dir[1] * turn.sin(),
                    d[2] * turn.cos() + side_dir[2] * turn.sin(),
                ]);
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
        _id: u64,
    ) -> Vec<[f64; 3]> {
        let dt = self.orbit_period_s / CYCLOID_PHASE_STEPS as f64;
        let d_phase = 2.0 * std::f64::consts::PI / CYCLOID_PHASE_STEPS as f64;
        let advance = CYCLOID_SPEED * dt;
        let mut points = vec![p0];
        let mut p = p0;
        let mut phase = phase0;
        let mut along: Option<[f64; 3]> = None;
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
/// and salts along it.
pub fn linea_compute(
    layers: &LayerSet,
    grid: &CubeSphere,
    body_radius: f64,
    rec: &Linea,
    blur: f64,
    stamper: &mut Stamper,
) -> StampOutput {
    let bedrock = layers.bedrock.as_slice();
    let ages = layers.surface_age.as_slice();
    let cell = grid.delta() * body_radius;
    let width = rec.width.max(2.0 * cell);
    let reach = match rec.kind {
        LineaKind::Band => 0.5 * width + 0.3 * width + 2.0 * cell,
        _ => width + 2.0 * cell,
    };
    let half_band = 0.5 * rec.spacing + cell;
    let n = rec.points.len();
    let mut best: HashMap<usize, (f64, f64)> = HashMap::new();
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
            let entry = best.entry(idx).or_insert((f64::INFINITY, 0.0));
            if along.abs() < entry.0 {
                *entry = (along.abs(), across);
            }
        }
    }
    let time = rec.time_gyr as f32;
    let mut cells: Vec<CellWrite> = best
        .into_iter()
        .map(|(idx, (_, across))| {
            let h = linea_profile(rec.kind, across, width, rec.relief, blur);
            let inside = across.abs() < 0.5 * width;
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
    StampOutput { cells }
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
}
