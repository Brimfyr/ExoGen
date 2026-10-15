// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! Stamping craters onto cubesphere layers.
//!
//! A footprint is gathered by a breadth-first walk over edge neighbours from the centre
//! cell while the great-circle distance stays within the ejecta radius. The crater floor
//! replaces the pre-impact topography relative to the mean height inside the rim; the rim
//! and ejecta are added to whatever is there, with a short blend across the rim. Cells under
//! the continuous ejecta have their surface age reset to the impact time.

use super::morphology::{CraterClass, Dimensions, ProfileTable};
use crate::constants::{
    BASIN_SCULPTURE, BURIAL_FACTOR_HIGH, BURIAL_FACTOR_LOW, CONTINUOUS_EJECTA_MULTIPLE,
    EJECTA_LOBE_MIN_CELLS, EJECTA_TEXTURE, FLOOR_TEXTURE, OBLIQUE_ASYMMETRY_DEG,
    OBLIQUE_BUTTERFLY_DEG, OBLIQUE_BUTTERFLY_MAX, OBLIQUE_DOWNRANGE_MAX, OBLIQUE_FORBIDDEN_DEG,
    OBLIQUE_FORBIDDEN_DEPTH, OBLIQUE_FORBIDDEN_HALF_WIDTH_DEG, RIM_CRENULATION_COMPLEX,
    RIM_CRENULATION_SIMPLE, SCULPTURE_AZIMUTH_CELLS, SCULPTURE_RADIAL_CELLS, TEXTURE_AZIMUTH_CELLS,
    TEXTURE_FADE_RADII, TEXTURE_MIN_CELLS, TEXTURE_RADIAL_CELLS,
};
use crate::constants::{
    EJECTA_COVERAGE_TAPER, EJECTA_EDGE_RAGGEDNESS, EJECTA_EXPONENT, EJECTA_FULL_COVER,
    EJECTA_PURE_RADII, MELT_EDGE_WIDTH, MELT_POND_CELLS, MELT_POND_THRESHOLD,
    OUTER_RING_HEIGHT_FRACTION, RING_CRENULATION, RING_PRESENCE, RING_WIDTH_FRACTION,
};
use crate::constants::{
    FLOOR_TEXTURE_CELLS, PEAK_MAX_MASSIFS, PEAK_SIGMA_FRACTION, PEAK_SPREAD, PIT_RIM_PRESENCE,
    PIT_WANDER, PLANAR_TEXTURE_RADII, RIM_HEIGHT_VARIATION, RIM_SHARE_BAND,
};
use crate::constants::{
    TERRACE_ABSENT, TERRACE_BENCH_SLOPE, TERRACE_BLOCK_LENGTH, TERRACE_BOW, TERRACE_HUMMOCK_SCALE,
    TERRACE_HUMMOCKS, TERRACE_INWARD, TERRACE_LEAST, TERRACE_OFFSET, TERRACE_RIM_SCALLOP,
    TERRACE_SCARP, TERRACE_SCARP_SPREAD, TERRACE_TOE, TERRACE_WALL_SHARE, TERRACE_WOBBLE,
};
use crate::fields::LayerSet;
use crate::geom::{angle_between, azimuth_from, cross, dot, local_basis, normalize};
use crate::grid::CubeSphere;
use crate::rng::Rng;
use rayon::prelude::*;

/// Reusable scratch space for footprints.
pub struct Stamper {
    visited: Vec<u32>,
    generation: u32,
    queue: Vec<usize>,
    /// For a corridor, the segment each queued cell came from.
    segments: Vec<u32>,
    /// (cell index, distance in metres, azimuth)
    pub footprint: Vec<(usize, f64, f64)>,
}

impl Stamper {
    pub fn new(grid: &CubeSphere) -> Self {
        Self {
            visited: vec![0; grid.cell_count()],
            generation: 0,
            queue: Vec::new(),
            segments: Vec::new(),
            footprint: Vec::new(),
        }
    }

    /// Collect the cells within `radius_m` of `center` on a body of radius `body_radius`.
    pub fn gather(&mut self, grid: &CubeSphere, body_radius: f64, center: [f64; 3], radius_m: f64) {
        self.generation = self.generation.wrapping_add(1);
        if self.generation == 0 {
            self.visited.iter_mut().for_each(|v| *v = 0);
            self.generation = 1;
        }
        self.footprint.clear();
        self.queue.clear();
        let theta_max = (radius_m / body_radius).min(std::f64::consts::PI);
        let (east, north) = local_basis(center);
        let start = grid.index(grid.locate(center).0);
        self.visited[start] = self.generation;
        self.queue.push(start);
        let mut head = 0;
        while head < self.queue.len() {
            let idx = self.queue[head];
            head += 1;
            let cell = grid.cell(idx);
            let d = grid.center(cell);
            let theta = angle_between(center, d);
            // Accept cells whose centre is within the radius plus half a cell so the edge
            // is not under-covered; the profile returns zero beyond the ejecta anyway.
            if theta > theta_max + 0.5 * grid.delta() {
                continue;
            }
            let phi = azimuth_from(center, east, north, d);
            self.footprint.push((idx, theta * body_radius, phi));
            for (di, dj) in [(1, 0), (-1, 0), (0, 1), (0, -1)] {
                let n = grid.index(grid.neighbor(cell, di, dj));
                if self.visited[n] != self.generation {
                    self.visited[n] = self.generation;
                    self.queue.push(n);
                }
            }
        }
    }

    /// Collect the cells within `radius_m` of a polyline of unit directions: the walk
    /// starts from the cells along each segment and spreads across it, each cell taking
    /// its distance from the nearest of the segments about the one it came from. The
    /// footprint holds the distance in metres and the position along the line, as the
    /// segment's index plus the share of it.
    pub fn gather_corridor(
        &mut self,
        grid: &CubeSphere,
        body_radius: f64,
        points: &[[f64; 3]],
        radius_m: f64,
    ) {
        self.generation = self.generation.wrapping_add(1);
        if self.generation == 0 {
            self.visited.iter_mut().for_each(|v| *v = 0);
            self.generation = 1;
        }
        self.footprint.clear();
        self.queue.clear();
        self.segments.clear();
        let n = points.len();
        if n < 2 {
            return;
        }
        let reach = radius_m / body_radius + 0.5 * grid.delta();
        for k in 0..n - 1 {
            let (a, b) = (points[k], points[k + 1]);
            let steps = (angle_between(a, b) / (0.5 * grid.delta())).ceil().max(1.0) as usize;
            for i in 0..=steps {
                let t = i as f64 / steps as f64;
                let p = normalize([
                    a[0] + (b[0] - a[0]) * t,
                    a[1] + (b[1] - a[1]) * t,
                    a[2] + (b[2] - a[2]) * t,
                ]);
                let idx = grid.index(grid.locate(p).0);
                if self.visited[idx] != self.generation {
                    self.visited[idx] = self.generation;
                    self.queue.push(idx);
                    self.segments.push(k as u32);
                }
            }
        }
        let mut head = 0;
        while head < self.queue.len() {
            let idx = self.queue[head];
            let from = self.segments[head] as usize;
            head += 1;
            let cell = grid.cell(idx);
            let d = grid.center(cell);
            let mut best = (f64::INFINITY, 0.0, from);
            for k in from.saturating_sub(1)..(from + 2).min(n - 1) {
                let (dist, share) = segment_distance(points[k], points[k + 1], d);
                if dist < best.0 {
                    best = (dist, k as f64 + share, k);
                }
            }
            if best.0 > reach {
                continue;
            }
            self.footprint.push((idx, best.0 * body_radius, best.1));
            for (di, dj) in [(1, 0), (-1, 0), (0, 1), (0, -1)] {
                let m = grid.index(grid.neighbor(cell, di, dj));
                if self.visited[m] != self.generation {
                    self.visited[m] = self.generation;
                    self.queue.push(m);
                    self.segments.push(best.2 as u32);
                }
            }
        }
    }
}

/// The angular distance from a direction to a short great-circle arc between two others,
/// and the share of the arc at which its nearest point lies: across the arc the sine of
/// the distance to its great circle, beyond its ends the chord to the nearer, both the
/// angle itself to the accuracy a corridor of a few cells needs.
fn segment_distance(a: [f64; 3], b: [f64; 3], d: [f64; 3]) -> (f64, f64) {
    let ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let ad = [d[0] - a[0], d[1] - a[1], d[2] - a[2]];
    let len2 = dot(ab, ab);
    let chord = |p: [f64; 3]| {
        let r = [d[0] - p[0], d[1] - p[1], d[2] - p[2]];
        dot(r, r).sqrt()
    };
    if len2 < 1e-30 {
        return (chord(a), 0.0);
    }
    let t = dot(ad, ab) / len2;
    if t <= 0.0 {
        return (chord(a), 0.0);
    }
    if t >= 1.0 {
        return (chord(b), 1.0);
    }
    let n = cross(a, b);
    (dot(d, n).abs() / dot(n, n).sqrt(), t)
}

/// One cell's new values from a stamp. Stamps are computed first and applied afterwards,
/// so stamps with disjoint footprints can be computed in parallel and applied in any order.
#[derive(Clone, Copy)]
pub struct CellWrite {
    pub idx: u32,
    pub bedrock: f32,
    /// Whether the surface fields below are written.
    pub surface: bool,
    pub age: f32,
    pub unit: u16,
    pub deep: f32,
    pub melt: f32,
    pub cover: f32,
}

/// The cells one stamp changes.
#[derive(Default)]
pub struct StampOutput {
    pub cells: Vec<CellWrite>,
    /// The composition an eruption's lava cells take (F55); none for every other stamp,
    /// whose cells keep the composition under them.
    pub composition: Option<f32>,
    /// Per cell, the stain a crack lays (F78) and the share of the older stain it cuts
    /// away; none for a stamp that lays none, whose renewed cells lose the stain under
    /// them.
    pub stain: Option<Vec<[f32; 2]>>,
}

impl StampOutput {
    /// Write the cells into the layers.
    pub fn apply(&self, layers: &mut LayerSet) {
        let bedrock = layers.bedrock.as_mut_slice();
        let ages = layers.surface_age.as_mut_slice();
        let units = layers.unit_id.as_mut_slice();
        let deep = layers.deep.as_mut_slice();
        let melt = layers.melt.as_mut_slice();
        let cover = layers.cover.as_mut_slice();
        let composition = layers.composition.as_mut_slice();
        let stain = layers.stain.as_mut_slice();
        for (j, c) in self.cells.iter().enumerate() {
            let i = c.idx as usize;
            bedrock[i] = c.bedrock;
            if c.surface {
                ages[i] = c.age;
                units[i] = c.unit;
                deep[i] = c.deep;
                melt[i] = c.melt;
                cover[i] = c.cover;
                if let Some(x) = self.composition
                    && c.unit == crate::units::UNIT_LAVA
                {
                    composition[i] = x;
                }
                if self.stain.is_none() {
                    stain[i] = 0.0;
                }
            }
            if let Some(s) = &self.stain {
                let [laid, cut] = s[j];
                stain[i] = (stain[i] * (1.0 - cut)).max(laid);
            }
        }
    }
}

thread_local! {
    static STAMPER: std::cell::RefCell<Option<Stamper>> = const { std::cell::RefCell::new(None) };
}

/// Run `f` with this thread's stamper for `grid`, creating or resizing it as needed, so
/// parallel stamping does not allocate a visited array per stamp.
pub fn with_stamper<R>(grid: &CubeSphere, f: impl FnOnce(&mut Stamper) -> R) -> R {
    // The stamper leaves the slot while in use: `f` may run parallel work, and rayon can
    // run another stamp on this thread while it waits, which then makes its own.
    let taken = STAMPER.with(|cell| cell.borrow_mut().take());
    let mut stamper = match taken {
        Some(s) if s.visited.len() == grid.cell_count() => s,
        _ => Stamper::new(grid),
    };
    let result = f(&mut stamper);
    STAMPER.with(|cell| {
        let mut slot = cell.borrow_mut();
        if slot.is_none() {
            *slot = Some(stamper);
        }
    });
    result
}

/// Footprints at least this large are computed in parallel chunks; a basin's restore at
/// the finishing resolution covers millions of cells.
const PARALLEL_FOOTPRINT: usize = 32_768;

/// Compute the writes of a footprint, in parallel when it is large. The result is the
/// same either way: cells are independent.
fn map_footprint<F>(footprint: &[(usize, f64, f64)], f: F) -> Vec<CellWrite>
where
    F: Fn(&(usize, f64, f64)) -> CellWrite + Sync + Send,
{
    if footprint.len() >= PARALLEL_FOOTPRINT {
        footprint.par_iter().map(f).collect()
    } else {
        footprint.iter().map(f).collect()
    }
}

/// Effective radius for an elliptical crater: scaled so the long axis lies along `azimuth`
/// and the area is preserved.
#[inline]
pub fn elliptical_radius(r: f64, phi: f64, azimuth: f64, ellipticity: f64) -> f64 {
    if ellipticity <= 1.0 {
        return r;
    }
    let (s, c) = (phi - azimuth).sin_cos();
    r * ((c * c) / ellipticity + (s * s) * ellipticity).sqrt()
}

pub struct StampParams {
    pub time_gyr: f32,
    pub unit_floor: u16,
    pub unit_ejecta: u16,
    /// Downrange direction, radians in the local east-north frame.
    pub azimuth: f64,
    /// Impact angle from horizontal, radians; π/2 is vertical.
    pub angle: f64,
    pub ellipticity: f64,
    /// Multiply the profile height before applying (1 for a full stamp).
    pub amplitude: f64,
    /// Multiply the wall's terraces, the ridge rings and the central peak, the crater's
    /// short wavelengths, which an ice shell relaxes less than its bowl (F70); equal to
    /// `amplitude` where the profile's table carries no relaxation.
    pub feature_amplitude: f64,
    /// When false the floor does not replace the old topography; the profile is purely
    /// additive (used to add resolved detail on top of a coarse stamp).
    pub replace_floor: bool,
    pub set_age: bool,
    /// Bury the relief already there under the ejecta according to its thickness.
    pub bury: bool,
    /// Seed of the crater's rim scalloping and texture; zero for a smooth stamp.
    pub seed: u64,
    /// The crater excavated the deep bedrock: its floor and near ejecta are deep material.
    pub deep: bool,
    /// The exposure on which the surface's maturity saturates, Gyr, for mixing thin
    /// ejecta with the ground it covers (F51).
    pub maturity_timescale: f64,
    /// Leave a cell alone where clean ice renewed it after this crater. The coarse loop
    /// lays a step's lanes before the step's craters, so an older crater of the same step
    /// must not land on them. The finish runs in the order things formed and needs no
    /// such guard: there the renewal's own turn takes what stood before it, by its
    /// weight, and a crater across its edge is softened rather than cut.
    pub renewal_guard: bool,
    /// Draw a complex crater's wall terraced as it collapsed (F67): in the finish, where
    /// the grid resolves it. The coarse loop keeps the plain wall, and the finish restores
    /// its terraces with the rest of its detail.
    pub terraced: bool,
}

/// Stamp one crater. Returns the number of cells touched.
pub fn stamp(
    layers: &mut LayerSet,
    grid: &CubeSphere,
    body_radius: f64,
    center: [f64; 3],
    profile: &ProfileTable,
    p: &StampParams,
    stamper: &mut Stamper,
) -> usize {
    let out = stamp_compute(layers, grid, body_radius, center, profile, p, stamper);
    out.apply(layers);
    out.cells.len()
}

/// Compute one crater's stamp without writing it.
pub fn stamp_compute(
    layers: &LayerSet,
    grid: &CubeSphere,
    body_radius: f64,
    center: [f64; 3],
    profile: &ProfileTable,
    p: &StampParams,
    stamper: &mut Stamper,
) -> StampOutput {
    let dim = &profile.dim;
    stamper.gather(
        grid,
        body_radius,
        center,
        dim.ejecta_radius * p.ellipticity.sqrt(),
    );
    if stamper.footprint.is_empty() {
        return StampOutput::default();
    }
    let shape = RimShape::new(p.seed, dim, stamper.footprint.len() >= TEXTURE_MIN_CELLS);
    let texture = (p.seed != 0 && stamper.footprint.len() >= TEXTURE_MIN_CELLS)
        .then(|| CraterTexture::new(p.seed, dim));
    // The wall as it collapsed (F67), where the grid draws it.
    let terraces = (p.terraced && p.seed != 0 && stamper.footprint.len() >= TEXTURE_MIN_CELLS)
        .then(|| Terraces::new(p.seed, dim))
        .flatten();
    let rings = Rings::new(p.seed, dim);
    let peak = Peak::new(p.seed, dim);
    let pit = Pit::new(p.seed, dim);
    let rise = RimRise::new(p.seed, stamper.footprint.len() >= TEXTURE_MIN_CELLS);
    let edge = EdgeShape::new(p.seed, dim.ejecta_radius, grid.delta() * body_radius);
    let bedrock = layers.bedrock.as_slice();
    let deep = layers.deep.as_slice();
    // Reference height: area-weighted mean of the pre-impact surface inside the rim; and
    // the azimuthal mean by radius under the ejecta, against which burial damps relief.
    const BINS: usize = 16;
    let span = (dim.ejecta_radius / dim.radius - 1.0).max(1e-6);
    let bin_of = |re: f64| -> usize {
        (((re / dim.radius - 1.0) / span) * BINS as f64)
            .floor()
            .clamp(0.0, (BINS - 1) as f64) as usize
    };
    let mut sum = 0.0;
    let mut wsum = 0.0;
    let mut deep_sum = 0.0;
    // Ring means are plain cell averages: cell areas differ by at most 40% on the grid,
    // and the mean only serves as a local reference.
    let mut ring_sum = [0.0f64; BINS];
    let mut ring_area = [0.0f64; BINS];
    for &(idx, r, phi) in &stamper.footprint {
        let re = elliptical_radius(r, phi, p.azimuth, p.ellipticity) * shape.factor(phi);
        if re <= dim.radius {
            let a = grid.cell_area(grid.cell(idx));
            sum += bedrock[idx] as f64 * a;
            deep_sum += deep[idx] as f64 * a;
            wsum += a;
        } else if p.bury && re < dim.ejecta_radius {
            let b = bin_of(re);
            ring_sum[b] += bedrock[idx] as f64;
            ring_area[b] += 1.0;
        }
    }
    let h_ref = if wsum > 0.0 {
        sum / wsum
    } else {
        bedrock[stamper.footprint[0].0] as f64
    };
    // The material the crater excavated: the deep bedrock, or what was already there.
    let source_deep = if p.deep {
        1.0
    } else if wsum > 0.0 {
        deep_sum / wsum
    } else {
        0.0
    };
    let continuous = CONTINUOUS_EJECTA_MULTIPLE * dim.radius;
    let ages = layers.surface_age.as_slice();
    let units = layers.unit_id.as_slice();
    let melts = layers.melt.as_slice();
    let covers = layers.cover.as_slice();
    let cleans = layers.clean.as_slice();
    let compute = |&(idx, r, phi): &(usize, f64, f64)| -> CellWrite {
        // A cell whose substrate was renewed from below after this crater, the clean ice
        // of a lane or a polar terrain, holds nothing of its deposit; nor of its relief,
        // where the renewal is not to take that at its own turn.
        let renewed_later = cleans[idx] > 0.0 && ages[idx] > p.time_gyr + 1e-6;
        if p.renewal_guard && renewed_later {
            return CellWrite {
                idx: idx as u32,
                bedrock: bedrock[idx],
                surface: false,
                age: 0.0,
                unit: 0,
                deep: 0.0,
                melt: 0.0,
                cover: 0.0,
            };
        }
        let re0 = elliptical_radius(r, phi, p.azimuth, p.ellipticity);
        let re = re0 * shape.stretch(re0 / dim.radius, phi, dim);
        let rw = terraces
            .as_ref()
            .map_or(re, |t| t.radius(re, phi, dim, profile.sigma));
        let mut v = profile.height(rw)
            * p.amplitude
            * ejecta_factor(re, dim.radius, phi, p.azimuth, p.angle);
        if let Some(t) = &texture {
            v = t.height(rw, phi, v);
        }
        if let Some(rr) = &rise
            && v > 0.0
        {
            v *= rr.at(re0 / dim.radius, phi);
        }
        if let Some(t) = &terraces {
            v += t.relief(re, phi, dim) * p.feature_amplitude;
        }
        if let Some(rg) = &rings {
            // Each ring keeps an outline of its own, as irregular as the rim's (F73).
            v += rg.relief(re0, phi, profile.sigma) * p.feature_amplitude;
        }
        if let Some(pk) = &peak {
            v += pk.relief(re0 / dim.radius, phi, profile.sigma / dim.radius) * p.feature_amplitude;
        }
        if let Some(pt) = &pit {
            v += pt.relief(re0 / dim.radius, phi, profile.sigma / dim.radius) * p.feature_amplitude;
        }
        let mut old = bedrock[idx] as f64;
        if p.bury && re > dim.radius && v > 0.0 {
            let b = bin_of(re);
            let mean = if ring_area[b] > 0.0 {
                ring_sum[b] / ring_area[b]
            } else {
                h_ref
            };
            old = mean + (old - mean) * burial_factor(old - mean, v);
        }
        let w = if p.replace_floor {
            // 1 inside 0.95 R, 0 outside 1.05 R.
            let t = (1.05 - re / dim.radius) / 0.10;
            t.clamp(0.0, 1.0)
        } else {
            0.0
        };
        let mut write = CellWrite {
            idx: idx as u32,
            bedrock: (v + w * h_ref + (1.0 - w) * old) as f32,
            surface: false,
            age: 0.0,
            unit: 0,
            deep: 0.0,
            melt: 0.0,
            cover: 0.0,
        };
        let cont = continuous * edge.factor_at(phi, re / dim.radius);
        let cover = ejecta_coverage(re, dim.radius, cont, dim.ejecta_radius);
        if p.set_age && cover > 0.0 && !renewed_later {
            // The deposit covers the ground by its cover (F51); the cell keeps the one
            // age whose exposure gives the mixture's maturity. Where the cover is under
            // a half the ground keeps its unit, its melt and its salts or lava.
            let age = mixed_age(p.time_gyr, ages[idx], cover, p.maturity_timescale);
            let thick = re <= dim.radius || cover >= 0.5;
            // Deep material by ballistic mixing: all ejecta near the rim, all local at the
            // continuous edge. Fresh ejecta covers older melt.
            let edge_x = cont / dim.radius;
            let mix = ((edge_x - re / dim.radius) / (edge_x - EJECTA_PURE_RADII)).clamp(0.0, 1.0);
            write.surface = true;
            write.age = age;
            write.unit = if re <= dim.radius {
                p.unit_floor
            } else if thick {
                p.unit_ejecta
            } else {
                units[idx]
            };
            write.deep = deep[idx] + ((source_deep - deep[idx] as f64) * mix) as f32;
            write.melt = if thick {
                texture.as_ref().map_or(0.0, |t| t.melt_fraction(rw, phi)) as f32
            } else {
                melts[idx]
            };
            // A crater in lava is made of lava through and through.
            write.cover = if !thick {
                covers[idx]
            } else if crate::stages::volcanism::is_lava_unit(write.unit) {
                1.0
            } else {
                0.0
            };
        }
        write
    };
    StampOutput {
        cells: map_footprint(&stamper.footprint, compute),
        composition: None,
        stain: None,
    }
}

/// Azimuthal modulation of ejecta thickness for an oblique impact, with `phi` the azimuth
/// from the centre and `azimuth` the downrange direction: ejecta thickens downrange and
/// thins uprange below 45°, an uprange zone of avoidance opens below 30°, and side lobes
/// grow below 15° (Gault and Wedekind 1978; Herrick and Forsberg-Taylor 2003).
#[inline]
pub fn ejecta_modulation(phi: f64, azimuth: f64, angle: f64) -> f64 {
    let dphi = phi - azimuth;
    let (s, c) = dphi.sin_cos();
    let deg = angle.to_degrees();
    let down = (1.0 - deg / OBLIQUE_ASYMMETRY_DEG).clamp(0.0, 1.0);
    let mut m = 1.0 + OBLIQUE_DOWNRANGE_MAX * down * c;
    if deg < OBLIQUE_FORBIDDEN_DEG {
        let strength = 1.0 - deg / OBLIQUE_FORBIDDEN_DEG;
        let half_width = OBLIQUE_FORBIDDEN_HALF_WIDTH_DEG.to_radians() * strength;
        // Angle from the uprange direction.
        let from_up = dphi.rem_euclid(2.0 * std::f64::consts::PI) - std::f64::consts::PI;
        m *= 1.0 - OBLIQUE_FORBIDDEN_DEPTH * strength * (-(from_up / half_width).powi(2)).exp();
    }
    if deg < OBLIQUE_BUTTERFLY_DEG {
        let strength = 1.0 - deg / OBLIQUE_BUTTERFLY_DEG;
        m *= 1.0 + OBLIQUE_BUTTERFLY_MAX * strength * s * s;
    }
    m.max(0.0)
}

/// The modulation faded in across the rim so the cavity side stays symmetric.
#[inline]
fn ejecta_factor(re: f64, radius: f64, phi: f64, azimuth: f64, angle: f64) -> f64 {
    let ramp = ((re / radius - 1.0) / 0.3).clamp(0.0, 1.0);
    if ramp <= 0.0 {
        return 1.0;
    }
    1.0 + (ejecta_modulation(phi, azimuth, angle) - 1.0) * ramp
}

/// Value-noise lattice of one crater on (azimuth, ln r/R), periodic around the crater,
/// hashed once per lattice point and interpolated per cell.
struct TextureLattice {
    period: usize,
    n_r: f64,
    v0: i64,
    rows: usize,
    values: Vec<f32>,
}

impl TextureLattice {
    fn new(rng: &Rng, octave: u64, period: u32, n_r: u32, x_min: f64, x_max: f64) -> Self {
        let v0 = (x_min.ln() * n_r as f64).floor() as i64 - 1;
        let v1 = (x_max.ln() * n_r as f64).ceil() as i64 + 1;
        let rows = (v1 - v0 + 1).max(2) as usize;
        let mut values = Vec::with_capacity(rows * period as usize);
        for row in 0..rows {
            let cv = ((v0 + row as i64) as u64) & 0xF_FFFF;
            for cu in 0..period as u64 {
                let key = (octave << 40) | (cu << 20) | cv;
                values.push((2.0 * rng.uniform(key, 0) - 1.0) as f32);
            }
        }
        Self {
            period: period as usize,
            n_r: n_r as f64,
            v0,
            rows,
            values,
        }
    }

    /// Sample at an azimuth in turns and a log radius.
    #[inline]
    fn sample(&self, turns: f64, log_x: f64) -> f64 {
        let u = turns * self.period as f64;
        let v = log_x * self.n_r;
        let (iu, iv) = (u.floor(), v.floor());
        let smooth = |t: f64| t * t * (3.0 - 2.0 * t);
        let (su, sv) = (smooth(u - iu), smooth(v - iv));
        let cu0 = (iu as i64).rem_euclid(self.period as i64) as usize;
        let cu1 = (cu0 + 1) % self.period;
        let row = ((iv as i64) - self.v0).clamp(0, self.rows as i64 - 2) as usize;
        let p = self.period;
        let (a, b) = (self.values[row * p + cu0], self.values[row * p + cu1]);
        let (c, d) = (
            self.values[(row + 1) * p + cu0],
            self.values[(row + 1) * p + cu1],
        );
        let lerp = |a: f32, b: f32, t: f64| a as f64 + (b as f64 - a as f64) * t;
        lerp(a, b, su) + (lerp(c, d, su) - lerp(a, b, su)) * sv
    }
}

/// Value-noise lattice on the crater's tangent plane in units of its radius, for the
/// floor, where a polar lattice would shrink its cells without limit toward the centre.
struct PlanarLattice {
    size: usize,
    half: f64,
    values: Vec<f32>,
}

impl PlanarLattice {
    fn new(rng: &Rng, octave: u64, cells_per_radius: u32, half: f64) -> Self {
        let n = (cells_per_radius as f64 * half).ceil() as usize;
        let size = 2 * n + 1;
        let mut values = Vec::with_capacity(size * size);
        for j in 0..size as u64 {
            for i in 0..size as u64 {
                let key = (1 << 60) | (octave << 40) | (i << 20) | j;
                values.push((2.0 * rng.uniform(key, 0) - 1.0) as f32);
            }
        }
        Self { size, half, values }
    }

    #[inline]
    fn sample(&self, u: f64, v: f64) -> f64 {
        let last = (self.size - 1) as f64;
        let fu = ((u / self.half + 1.0) * 0.5 * last).clamp(0.0, last - 1e-6);
        let fv = ((v / self.half + 1.0) * 0.5 * last).clamp(0.0, last - 1e-6);
        let (iu, iv) = (fu.floor(), fv.floor());
        let smooth = |t: f64| t * t * (3.0 - 2.0 * t);
        let (su, sv) = (smooth(fu - iu), smooth(fv - iv));
        let (i, j) = (iu as usize, iv as usize);
        let at = |i: usize, j: usize| self.values[j * self.size + i] as f64;
        let a = at(i, j) + (at(i + 1, j) - at(i, j)) * su;
        let b = at(i, j + 1) + (at(i + 1, j + 1) - at(i, j + 1)) * su;
        a + (b - a) * sv
    }
}

/// The texture of one crater: two octaves of hummocks on the ejecta and the floor, and
/// radial sculpture for basins and peak-ring craters.
pub struct CraterTexture {
    dim: Dimensions,
    hummocks: [TextureLattice; 2],
    sculpture: Option<[TextureLattice; 2]>,
    floor: Option<[PlanarLattice; 2]>,
    ponds: Option<PlanarLattice>,
}

impl CraterTexture {
    pub fn new(seed: u64, dim: &Dimensions) -> Self {
        let rng = Rng::new(seed ^ 0x7E57_A11E);
        // The lattice covers the whole interior: the floor texture runs from the centre
        // out, and a lattice that began at the floor radius would clamp every floor cell
        // to one row and turn the floor into spokes.
        let x_min = if dim.class == CraterClass::Simple {
            1.0
        } else {
            0.05
        };
        let x_max = TEXTURE_FADE_RADII.1 + 0.1;
        let lattice = |rng: &Rng, octave: u64, n_phi: u32, n_r: u32| {
            TextureLattice::new(rng, octave, n_phi << octave, n_r << octave, x_min, x_max)
        };
        let hummocks = [
            lattice(&rng, 0, TEXTURE_AZIMUTH_CELLS, TEXTURE_RADIAL_CELLS),
            lattice(&rng, 1, TEXTURE_AZIMUTH_CELLS, TEXTURE_RADIAL_CELLS),
        ];
        let sculpture = (dim.class == CraterClass::MultiRing).then(|| {
            let rng = Rng::new(seed ^ 0x5C0_1F7E);
            [
                lattice(&rng, 0, SCULPTURE_AZIMUTH_CELLS, SCULPTURE_RADIAL_CELLS),
                lattice(&rng, 1, SCULPTURE_AZIMUTH_CELLS, SCULPTURE_RADIAL_CELLS),
            ]
        });
        let floor = (dim.class != CraterClass::Simple).then(|| {
            let rng = Rng::new(seed ^ 0xF100_12A7);
            [
                PlanarLattice::new(&rng, 0, FLOOR_TEXTURE_CELLS, PLANAR_TEXTURE_RADII),
                PlanarLattice::new(&rng, 1, 2 * FLOOR_TEXTURE_CELLS, PLANAR_TEXTURE_RADII),
            ]
        });
        let ponds = (dim.class != CraterClass::Simple)
            .then(|| PlanarLattice::new(&Rng::new(seed ^ 0x90ED_5EED), 0, MELT_POND_CELLS, 1.1));
        Self {
            dim: *dim,
            hummocks,
            sculpture,
            floor,
            ponds,
        }
    }

    /// Floor noise in -1 to 1 on the tangent plane at a radius in crater radii.
    #[inline]
    pub fn noise_floor(&self, phi: f64, x: f64) -> f64 {
        match &self.floor {
            Some(f) => {
                let (s, c) = phi.sin_cos();
                (f[0].sample(x * c, x * s) + 0.5 * f[1].sample(x * c, x * s)) / 1.5
            }
            None => 0.0,
        }
    }

    /// Azimuth in turns and log radius, the lattice coordinates.
    #[inline]
    fn coords(phi: f64, x: f64) -> (f64, f64) {
        let tau = 2.0 * std::f64::consts::PI;
        (phi.rem_euclid(tau) / tau, x.max(1e-3).ln())
    }

    /// Two-octave noise in -1 to 1 at an azimuth and a radius in crater radii.
    #[inline]
    pub fn noise(&self, phi: f64, x: f64) -> f64 {
        let (t, l) = Self::coords(phi, x);
        self.noise_at(t, l)
    }

    #[inline]
    fn noise_at(&self, turns: f64, log_x: f64) -> f64 {
        (self.hummocks[0].sample(turns, log_x) + 0.5 * self.hummocks[1].sample(turns, log_x)) / 1.5
    }

    /// Impact melt cover, 0 to 1: ponds on the floor of a complex crater or basin and flows
    /// just outside the rim, where the coarse texture octave is high, with soft edges.
    #[inline]
    pub fn melt_fraction(&self, re: f64, phi: f64) -> f64 {
        let dim = &self.dim;
        if dim.class == CraterClass::Simple {
            return 0.0;
        }
        if re >= dim.floor_radius {
            return 0.0;
        }
        // Ponds on a coarse planar lattice, so they are patches rather than specks.
        let x = re / dim.radius;
        let (s, c) = phi.sin_cos();
        let n = self.ponds.as_ref().map_or(0.0, |p| p.sample(x * c, x * s));
        let s = ((n - MELT_POND_THRESHOLD) / MELT_EDGE_WIDTH + 0.5).clamp(0.0, 1.0);
        s * s * (3.0 - 2.0 * s)
    }

    /// Apply the texture to a profile height at radius `re` and azimuth `phi`: hummocky
    /// ejecta beyond the rim fading out past the continuous ejecta, radial sculpture for
    /// basins, and a rough floor inside.
    #[inline]
    pub fn height(&self, re: f64, phi: f64, v: f64) -> f64 {
        let dim = &self.dim;
        let x = re / dim.radius;
        if re > dim.radius {
            let (fade_in, fade_out) = TEXTURE_FADE_RADII;
            let taper = ((fade_out - x) / (fade_out - fade_in)).clamp(0.0, 1.0);
            if taper <= 0.0 {
                return v;
            }
            let (t, l) = Self::coords(phi, x);
            let mut f = 1.0 + EJECTA_TEXTURE * taper * self.noise_at(t, l);
            if let Some(s) = &self.sculpture {
                let n = (s[0].sample(t, l) + 0.5 * s[1].sample(t, l)) / 1.5;
                f *= 1.0 + BASIN_SCULPTURE * taper * n;
            }
            v * f.max(0.0)
        } else if re < dim.floor_radius && dim.class != CraterClass::Simple {
            v + FLOOR_TEXTURE * dim.depth * self.noise_floor(phi, x)
        } else {
            v
        }
    }
}

/// The central peak of a complex crater as a cluster of massifs on the tangent plane, in
/// units of the crater radius. Degradation widens and lowers each massif.
struct Peak {
    /// Per massif: offset u, offset v, sigma, height in metres.
    massifs: Vec<(f64, f64, f64, f64)>,
}

impl Peak {
    fn new(seed: u64, dim: &Dimensions) -> Option<Self> {
        if dim.peak_height <= 0.0 || dim.radius <= 0.0 {
            return None;
        }
        let radius = dim.peak_radius / dim.radius;
        let rng = Rng::new(seed ^ 0x9EA4_0C11);
        let count = if seed == 0 {
            1
        } else {
            1 + (rng.uniform(0, 0) * PEAK_MAX_MASSIFS as f64) as usize
        };
        let massifs = (0..count as u64)
            .map(|i| {
                let (offset, angle) = if count == 1 {
                    (0.0, 0.0)
                } else {
                    (
                        PEAK_SPREAD * radius * rng.uniform(i, 1).sqrt(),
                        2.0 * std::f64::consts::PI * rng.uniform(i, 2),
                    )
                };
                let sigma = PEAK_SIGMA_FRACTION * radius * (0.7 + 0.6 * rng.uniform(i, 3))
                    / (count as f64).sqrt();
                let height = dim.peak_height * (0.7 + 0.5 * rng.uniform(i, 4));
                (offset * angle.cos(), offset * angle.sin(), sigma, height)
            })
            .collect();
        Some(Self { massifs })
    }

    /// Peak relief at a radius in crater radii and an azimuth, after a blur in radii.
    #[inline]
    fn relief(&self, x: f64, phi: f64, blur: f64) -> f64 {
        let (s, c) = phi.sin_cos();
        let (u, v) = (x * c, x * s);
        let mut h = 0.0;
        for &(mu, mv, sigma, height) in &self.massifs {
            let s2 = sigma * sigma + blur * blur;
            let d2 = (u - mu) * (u - mu) + (v - mv) * (v - mv);
            if d2 < 16.0 * s2 {
                h += height * (sigma * sigma / s2) * (-d2 / s2).exp();
            }
        }
        h
    }
}

/// The central pit of a crater on ice (F72): a depression its uplift collapsed into,
/// flat-floored over its inner half, ringed by a raised rim that breaks up around the
/// larger pits, and in the largest craters a dome risen in it; the pit's outline, its rim
/// and its dome each of their own shape.
struct Pit {
    /// Radii in crater radii; depth, rim height and dome height in metres.
    radius: f64,
    rim_radius: f64,
    depth: f64,
    rim_height: f64,
    dome_radius: f64,
    dome_height: f64,
    seed: u64,
}

impl Pit {
    fn new(seed: u64, dim: &Dimensions) -> Option<Self> {
        if dim.pit_radius <= 0.0 || dim.radius <= 0.0 {
            return None;
        }
        let r = dim.radius;
        Some(Self {
            radius: dim.pit_radius / r,
            rim_radius: dim.pit_rim_radius.max(1.2 * dim.pit_radius) / r,
            depth: dim.pit_depth,
            rim_height: dim.pit_rim_height,
            dome_radius: dim.dome_radius / r,
            dome_height: dim.dome_height,
            seed,
        })
    }

    /// The outline factor of one of the pit's parts at an azimuth.
    fn wander(&self, salt: u64, phi: f64) -> f64 {
        if self.seed == 0 {
            return 1.0;
        }
        1.0 + PIT_WANDER * harmonics(self.seed ^ salt, phi, &[2.0, 3.0, 5.0]).clamp(-2.0, 2.0)
    }

    /// Relief at a radius in crater radii and an azimuth, after a blur in radii, m.
    fn relief(&self, x: f64, phi: f64, blur: f64) -> f64 {
        let reach = self.rim_radius * (1.0 + 2.0 * PIT_WANDER) + 3.0 * blur;
        if x > reach {
            return 0.0;
        }
        let soft = |t: f64| {
            let t = t.clamp(0.0, 1.0);
            t * t * (3.0 - 2.0 * t)
        };
        let rp = self.radius * self.wander(0x9171_0001, phi);
        let rr = (self.rim_radius * self.wander(0x9171_0002, phi)).max(1.1 * rp);
        // The pit, flat-floored over its inner half.
        let wall = (0.55 * rp).hypot(blur);
        let pit = 1.0 - soft((x - 0.45 * rp) / wall);
        // Its rim: a crest at the pit's edge, falling into the pit and outward to the
        // rim's outer radius, broken around its circuit.
        let presence = if self.seed == 0 {
            1.0
        } else {
            (PIT_RIM_PRESENCE.0
                + PIT_RIM_PRESENCE.1 * harmonics(self.seed ^ 0x9171_0003, phi, &[3.0, 5.0, 8.0]))
            .clamp(0.0, 1.0)
        };
        let inner = soft((x - 0.7 * rp) / (0.3 * rp).hypot(blur));
        let outer = 1.0 - soft((x - rp) / (rr - rp).hypot(blur));
        let mut h = self.rim_height * presence * inner * outer - self.depth * pit;
        if self.dome_radius > 0.0 {
            let rd = self.dome_radius * self.wander(0x9171_0004, phi);
            let t = (x / rd.hypot(blur)).min(1.0);
            h += self.dome_height * (1.0 - t * t).max(0.0).powf(1.5);
        }
        h
    }
}

/// A crater's rim height around it (F73): a factor by azimuth, applied in full on the rim
/// and its flanks and fading out over the ejecta beyond.
struct RimRise {
    table: [f64; 64],
}

impl RimRise {
    fn new(seed: u64, enabled: bool) -> Option<Self> {
        if seed == 0 || !enabled {
            return None;
        }
        let mut table = [1.0; 64];
        for (b, t) in table.iter_mut().enumerate() {
            let phi = b as f64 / 64.0 * 2.0 * std::f64::consts::PI;
            *t = (1.0
                + RIM_HEIGHT_VARIATION * harmonics(seed ^ 0x81E5_0000, phi, &[1.0, 2.0, 3.0, 5.0]))
            .max(0.3);
        }
        Some(Self { table })
    }

    /// The factor at an elliptical radius in crater radii and an azimuth.
    fn at(&self, x: f64, phi: f64) -> f64 {
        let tau = 2.0 * std::f64::consts::PI;
        let u = phi.rem_euclid(tau) / tau * 64.0;
        let i = (u.floor() as usize) % 64;
        let f = u - u.floor();
        let t = self.table[i] * (1.0 - f) + self.table[(i + 1) % 64] * f;
        let fade = 1.0 - ((x - 1.3) / 0.7).clamp(0.0, 1.0);
        1.0 + (t - 1.0) * fade
    }
}

/// The share of the rim's scalloped outline a radius carries (F73): all of it at the rim
/// and beyond, none inward of the outer RIM_SHARE_BAND of the wall, so the inner terraces,
/// the floor and the central features keep shapes of their own.
fn rim_share(x: f64, dim: &Dimensions) -> f64 {
    let floor = if dim.class == CraterClass::Simple || dim.radius <= 0.0 {
        0.0
    } else {
        (dim.floor_radius / dim.radius).clamp(0.0, 0.95)
    };
    let x_in = 1.0 - RIM_SHARE_BAND * (1.0 - floor);
    let t = ((x - x_in) / (1.0 - x_in).max(1e-3)).clamp(0.0, 1.0);
    t * t * (3.0 - 2.0 * t)
}

/// A sum of azimuthal harmonics with seeded phases and weights, scaled to unit
/// root-mean-square.
fn harmonics(seed: u64, phi: f64, orders: &[f64]) -> f64 {
    let rng = Rng::new(seed);
    let mut s = 0.0;
    let mut power = 0.0;
    for (i, m) in orders.iter().enumerate() {
        let phase = rng.uniform(i as u64, 0) * 2.0 * std::f64::consts::PI;
        let w = (0.5 + rng.uniform(i as u64, 1)) / m.sqrt();
        s += w * (m * phi + phase).cos();
        power += 0.5 * w * w;
    }
    s / power.sqrt()
}

/// Scalloping of a rim as a function of azimuth: a few low harmonics.
pub fn crenulation(seed: u64, phi: f64) -> f64 {
    harmonics(seed ^ 0xC2E4_0A17, phi, &[3.0, 5.0, 8.0, 13.0])
}

/// Lobes of a continuous-ejecta edge as a function of azimuth: high harmonics at
/// irregular orders, so the edge is ragged in many small lobes, a fraction of a radius
/// wide, rather than petalled in a few.
pub fn lobes(seed: u64, phi: f64, max_order: f64) -> f64 {
    let orders: Vec<f64> = [11.0, 17.0, 26.0, 40.0, 61.0]
        .into_iter()
        .filter(|&m| m <= max_order)
        .collect();
    if orders.is_empty() {
        return 0.0;
    }
    harmonics(seed ^ 0xED6E_5A1F, phi, &orders)
}

/// The lobate outline of a crater's continuous ejecta, tabulated around the azimuth as a
/// radius factor. The stamp's age reset and the ray field share it, so the maturity disc
/// and the bright ejecta end on the same line.
pub struct EdgeShape {
    table: [f64; 64],
}

impl EdgeShape {
    /// The edge of a blanket of the given radius, m, drawn on a grid of the given cell,
    /// m: lobes narrower than a few cells along the edge are left out.
    pub fn new(seed: u64, ejecta_radius: f64, cell: f64) -> Self {
        let mut table = [1.0; 64];
        if seed != 0 {
            let max_order = 2.0 * std::f64::consts::PI * ejecta_radius
                / (EJECTA_LOBE_MIN_CELLS * cell.max(1.0));
            for (b, t) in table.iter_mut().enumerate() {
                let phi = b as f64 / 64.0 * 2.0 * std::f64::consts::PI;
                *t = 1.0 + EJECTA_EDGE_RAGGEDNESS * lobes(seed, phi, max_order);
            }
        }
        Self { table }
    }

    #[inline]
    pub fn factor(&self, phi: f64) -> f64 {
        let x = phi.rem_euclid(2.0 * std::f64::consts::PI) / (2.0 * std::f64::consts::PI) * 64.0;
        let i = (x.floor() as usize) % 64;
        let f = x - x.floor();
        self.table[i] * (1.0 - f) + self.table[(i + 1) % 64] * f
    }

    /// The factor at a distance `x` in crater radii: the blanket is continuous round the
    /// rim and its lobes develop outward, in full from two radii on.
    #[inline]
    pub fn factor_at(&self, phi: f64, x: f64) -> f64 {
        1.0 + (self.factor(phi) - 1.0) * (x - 1.0).clamp(0.0, 1.0)
    }
}

/// The ridge rings of one crater, each with an outline of its own as irregular as the
/// rim's and a few gaps, tabulated around the azimuth (F73); a multi-ring basin also
/// carries a subdued outer ring at sqrt(2) radii. Degradation widens and lowers the
/// ridges.
struct Rings {
    /// Per ring: radius by azimuth, presence by azimuth, height.
    items: Vec<([f64; 64], [f64; 64], f64)>,
    sigma: f64,
}

impl Rings {
    fn new(seed: u64, dim: &Dimensions) -> Option<Self> {
        let mut specs: Vec<(f64, f64)> = dim.rings[..dim.ring_count].to_vec();
        if dim.class == CraterClass::MultiRing {
            specs.push((
                std::f64::consts::SQRT_2 * dim.radius,
                OUTER_RING_HEIGHT_FRACTION * dim.rim_height,
            ));
        }
        if specs.is_empty() {
            return None;
        }
        let tau = 2.0 * std::f64::consts::PI;
        let items = specs
            .iter()
            .enumerate()
            .map(|(k, &(r0, height))| {
                let mut radius = [r0; 64];
                let mut presence = [1.0; 64];
                if seed != 0 {
                    let shape_seed = seed ^ (0x9A5E_0000 + k as u64);
                    let gap_seed = seed ^ (0xA6B7_0000 + k as u64);
                    for b in 0..64 {
                        let phi = b as f64 / 64.0 * tau;
                        radius[b] = r0 * (1.0 + RING_CRENULATION * crenulation(shape_seed, phi));
                        presence[b] = (RING_PRESENCE.0
                            + RING_PRESENCE.1 * crenulation(gap_seed, phi))
                        .clamp(0.0, 1.0);
                    }
                }
                (radius, presence, height)
            })
            .collect();
        Some(Self {
            items,
            sigma: RING_WIDTH_FRACTION * dim.radius,
        })
    }

    /// Ring relief at a radius and azimuth after a degradation blur of `blur` metres.
    #[inline]
    fn relief(&self, re: f64, phi: f64, blur: f64) -> f64 {
        let tau = 2.0 * std::f64::consts::PI;
        let x = phi.rem_euclid(tau) / tau * 64.0;
        let i = (x.floor() as usize) % 64;
        let j = (i + 1) % 64;
        let f = x - x.floor();
        let sigma = (self.sigma * self.sigma + blur * blur).sqrt();
        let scale = self.sigma / sigma;
        let mut h = 0.0;
        for (radius, presence, height) in &self.items {
            let r = radius[i] * (1.0 - f) + radius[j] * f;
            let p = presence[i] * (1.0 - f) + presence[j] * f;
            let z = (re - r) / sigma;
            if z.abs() < 4.0 {
                h += height * p * scale * (-z * z).exp();
            }
        }
        h
    }
}

/// Azimuth bins of a crater's terrace tables.
const TERRACE_BINS: usize = 256;

/// The wall of a complex crater as it collapsed (F67): the transient cavity's rim slumped
/// inward on inward-facing normal faults in rings of blocks, terrace after terrace from
/// the rim inward, each narrower than the one outside it (Pearce and Melosh 1986). A ring
/// is broken into blocks of their own length, throw and bow, some of which slumped no
/// terrace at all, so that the terraces come and go around the crater and the wall's width
/// changes from block to block; the outermost blocks scallop the rim crest and the
/// innermost blocks' toes the floor's edge.
struct Terraces {
    /// Per terrace, from the rim inward, by azimuth: its inner scarp's radius, whether the
    /// terrace stands there and its scarp's width, in crater radii.
    levels: Vec<(Vec<f32>, Vec<f32>, Vec<f32>)>,
    /// The broken surface of the slumped blocks, and its relief, m.
    hummocks: PlanarLattice,
    hummock_relief: f64,
    /// The floor's edge by azimuth, crater radii.
    toe: Vec<f32>,
    /// The rim crest's radius factor by azimuth, as `RimShape` keeps it.
    rim: [f64; 64],
    /// The width of the rim's headscarp, crater radii.
    head: f64,
}

impl Terraces {
    fn new(seed: u64, dim: &Dimensions) -> Option<Self> {
        if !matches!(
            dim.class,
            CraterClass::Complex
                | CraterClass::PeakRing
                | CraterClass::CentralPit
                | CraterClass::CentralDome
        ) || dim.terrace_width <= 0.0
            || dim.radius <= 0.0
        {
            return None;
        }
        let r = dim.radius;
        let wall = (r - dim.floor_radius).max(0.0);
        // The terraces' widths from the rim inward.
        let mut widths = Vec::new();
        let (mut w, mut sum) = (dim.terrace_width, 0.0);
        while w >= TERRACE_LEAST * dim.terrace_width && sum + w <= TERRACE_WALL_SHARE * wall {
            widths.push(w);
            sum += w;
            w *= TERRACE_INWARD;
        }
        if widths.is_empty() {
            return None;
        }
        let tau = 2.0 * std::f64::consts::PI;
        let rng = Rng::new(seed ^ 0x7E44_ACE5);
        // A ring of blocks round a circumference of `length` m, each about `block` m long:
        // by bin, the position within the block (0 to 1), the block's offset draw (-1 to
        // 1), whether it stands and its scarp's width draw (-1 to 1).
        let ring = |level: u64, length: f64, block: f64| -> Vec<(f64, f64, bool, f64)> {
            let n = ((length / block).round() as usize).max(3);
            let key = |i: usize| level * 4096 + i as u64;
            let lengths: Vec<f64> = (0..n)
                .map(|i| 2f64.powf(2.0 * rng.uniform(key(i), 0) - 1.0))
                .collect();
            let total: f64 = lengths.iter().sum();
            let mut cuts = Vec::with_capacity(n);
            let mut acc = 0.0;
            for l in &lengths {
                cuts.push(acc);
                acc += l / total;
            }
            let start = rng.uniform(key(4095), 1);
            (0..TERRACE_BINS)
                .map(|b| {
                    let turn = (b as f64 / TERRACE_BINS as f64 - start).rem_euclid(1.0);
                    let i = cuts.partition_point(|&c| c <= turn).max(1) - 1;
                    let end = if i + 1 < n { cuts[i + 1] } else { 1.0 };
                    let u = ((turn - cuts[i]) / (end - cuts[i]).max(1e-12)).clamp(0.0, 1.0);
                    let offset = 2.0 * rng.uniform(key(i), 2) - 1.0;
                    let stands = rng.uniform(key(i), 3) >= TERRACE_ABSENT;
                    let width = 2.0 * rng.uniform(key(i), 4) - 1.0;
                    (u, offset, stands, width)
                })
                .collect()
        };
        let bow = |u: f64| (std::f64::consts::PI * u).sin() - 2.0 / std::f64::consts::PI;
        let mut levels = Vec::with_capacity(widths.len());
        let mut outer = r;
        let mut first = Vec::new();
        let mut last = Vec::new();
        for (k, &w) in widths.iter().enumerate() {
            let inner = outer - w;
            let blocks = ring(k as u64, tau * (inner + 0.5 * w), TERRACE_BLOCK_LENGTH * w);
            // The scarp's trace wanders at the scale of a few terrace widths.
            let wander_seed = seed ^ (0x7E44_0000 + k as u64);
            let orders = tau * (inner + 0.5 * w) / (2.0 * w);
            let scarp: Vec<f32> = blocks
                .iter()
                .enumerate()
                .map(|(b, &(u, offset, _, _))| {
                    let phi = b as f64 / TERRACE_BINS as f64 * tau;
                    let at = inner
                        + w * (TERRACE_OFFSET * offset
                            + TERRACE_BOW * (std::f64::consts::PI * u).sin()
                            + TERRACE_WOBBLE * lobes(wander_seed, phi, orders));
                    (at / r) as f32
                })
                .collect();
            let stands: Vec<f32> = blocks
                .iter()
                .map(|&(_, _, s, _)| if s { 1.0 } else { 0.0 })
                .collect();
            let scarp_width: Vec<f32> = blocks
                .iter()
                .map(|&(_, _, _, d)| {
                    (TERRACE_SCARP * w * (1.0 + TERRACE_SCARP_SPREAD * d) / r) as f32
                })
                .collect();
            levels.push((scarp, stands, scarp_width));
            if k == 0 {
                first = blocks.clone();
            }
            last = blocks;
            outer = inner;
        }
        // The rim crest, scalloped by the outermost slumps: out in each one's middle, in
        // at the cusps between them.
        let w1 = widths[0];
        let mut rim = [1.0; 64];
        for (b, t) in rim.iter_mut().enumerate() {
            let (u, _, _, _) = first[b * TERRACE_BINS / 64];
            *t = 1.0 / (1.0 + TERRACE_RIM_SCALLOP * w1 / r * bow(u));
        }
        // The floor's edge at the innermost blocks' toes, bulging in at each one's middle.
        let wk = widths[widths.len() - 1];
        let toe: Vec<f32> = last
            .iter()
            .map(|&(u, _, _, _)| ((dim.floor_radius - TERRACE_TOE * wk * bow(u)) / r) as f32)
            .collect();
        // The slumped blocks are broken: hummocks at half the outermost terrace's width,
        // standing a share of a terrace's drop.
        let cells = ((r / (TERRACE_HUMMOCK_SCALE * w1)).round() as u32).clamp(2, 64);
        let hummocks = PlanarLattice::new(&Rng::new(seed ^ 0x40CC_0C55), 0, cells, 1.1);
        let hummock_relief =
            TERRACE_HUMMOCKS * (dim.depth + dim.rim_height) / (widths.len() + 1) as f64;
        Some(Self {
            levels,
            hummocks,
            hummock_relief,
            toe,
            rim,
            head: TERRACE_SCARP * w1 / r,
        })
    }

    /// The radius, m, at which to read a crater's plain profile for the terraced wall at
    /// the elliptical radius `re`, m, and azimuth `phi`, with the scarps softened by the
    /// crater's degradation `blur`, m: the floor reaches the blocks' toes, each terrace's
    /// bench keeps a little of the wall's slope and its scarp drops to the next, and the
    /// wall below the terraces is a plain slope to the floor.
    fn radius(&self, re: f64, phi: f64, dim: &Dimensions, blur: f64) -> f64 {
        let r = dim.radius;
        let x = re / r;
        if x >= 1.0 {
            return re;
        }
        let xf = dim.floor_radius / r;
        let blur = blur / r;
        let turn = (phi / (2.0 * std::f64::consts::PI)).rem_euclid(1.0) * TERRACE_BINS as f64;
        let i = (turn.floor() as usize) % TERRACE_BINS;
        let j = (i + 1) % TERRACE_BINS;
        let f = turn - turn.floor();
        let at = |v: &[f32]| v[i] as f64 * (1.0 - f) + v[j] as f64 * f;
        let toe = at(&self.toe).clamp(0.05, 0.98);
        if x <= toe {
            return x * xf / toe * r;
        }
        let lin = |y: f64| xf + (y - toe) * (1.0 - xf) / (1.0 - toe);
        let soft = |t: f64| {
            let t = t.clamp(0.0, 1.0);
            t * t * (3.0 - 2.0 * t)
        };
        let mut value = lin(x);
        let (mut outer, mut outer_w) = (1.0, self.head.hypot(blur));
        for (k, (scarp, stands, width)) in self.levels.iter().enumerate() {
            let w = at(width).hypot(blur);
            let inner = at(scarp).min(outer - 0.5 * w).max(toe + 0.5 * w);
            let mid = 0.5 * (inner + outer);
            let dev = (1.0 - TERRACE_BENCH_SLOPE) * (lin(mid) - lin(x));
            // The bench rises across its own scarp and gives way across the one above it;
            // the outermost gives way to the rim's headscarp before the crest.
            let rise = soft((x - inner) / w + 0.5);
            let fall = if k == 0 {
                soft((outer - x) / outer_w - 0.5)
            } else {
                soft((outer - x) / outer_w + 0.5)
            };
            value += at(stands) * rise * fall * dev;
            outer = inner;
            outer_w = w;
        }
        value * r
    }

    /// The broken surface of the slumped blocks at the elliptical radius `re`, m, and
    /// azimuth `phi`: hummocks across the terraced wall, fading at the floor's edge and the
    /// rim's headscarp, m.
    fn relief(&self, re: f64, phi: f64, dim: &Dimensions) -> f64 {
        let r = dim.radius;
        let x = re / r;
        if x >= 1.0 {
            return 0.0;
        }
        let xf = dim.floor_radius / r;
        let soft = |t: f64| {
            let t = t.clamp(0.0, 1.0);
            t * t * (3.0 - 2.0 * t)
        };
        let window = soft((x - xf) / (0.5 * (1.0 - xf))) * soft((1.0 - x) / self.head - 0.5);
        if window <= 0.0 {
            return 0.0;
        }
        let (s, c) = phi.sin_cos();
        self.hummock_relief * window * self.hummocks.sample(x * c, x * s)
    }
}

/// Rim scalloping of one crater, tabulated around the azimuth as a radius factor.
struct RimShape {
    table: [f64; 64],
    /// The outline of the floor and the bowl where no terraces shape them, a crenulation
    /// of their own (F73).
    inner: [f64; 64],
    flat: bool,
}

impl RimShape {
    fn new(seed: u64, dim: &Dimensions, enabled: bool) -> Self {
        let mut table = [1.0; 64];
        let mut inner = [1.0; 64];
        let flat = seed == 0 || !enabled;
        // A complex crater's rim is scalloped by the slumps of its outermost terrace (F67),
        // and its terraces and their toes shape the wall and the floor's edge within.
        if !flat && let Some(t) = Terraces::new(seed, dim) {
            return Self {
                table: t.rim,
                inner,
                flat: false,
            };
        }
        if !flat {
            let amplitude = if dim.class == CraterClass::Simple {
                RIM_CRENULATION_SIMPLE
            } else {
                RIM_CRENULATION_COMPLEX
            };
            for (b, (t, u)) in table.iter_mut().zip(inner.iter_mut()).enumerate() {
                let phi = b as f64 / 64.0 * 2.0 * std::f64::consts::PI;
                *t = 1.0 / (1.0 + amplitude * crenulation(seed, phi));
                *u = 1.0 / (1.0 + amplitude * crenulation(seed ^ 0x1DE5_0000, phi));
            }
        }
        Self { table, inner, flat }
    }

    #[inline]
    fn lookup(table: &[f64; 64], phi: f64) -> f64 {
        let x = phi.rem_euclid(2.0 * std::f64::consts::PI) / (2.0 * std::f64::consts::PI) * 64.0;
        let i = (x.floor() as usize) % 64;
        let f = x - x.floor();
        table[i] * (1.0 - f) + table[(i + 1) % 64] * f
    }

    #[inline]
    fn factor(&self, phi: f64) -> f64 {
        if self.flat {
            return 1.0;
        }
        Self::lookup(&self.table, phi)
    }

    /// The stretch of an elliptical radius in crater radii at an azimuth (F73): the rim's
    /// outline at the rim and beyond, the floor's own within.
    #[inline]
    fn stretch(&self, x: f64, phi: f64, dim: &Dimensions) -> f64 {
        if self.flat {
            return 1.0;
        }
        let share = rim_share(x, dim);
        1.0 + (Self::lookup(&self.table, phi) - 1.0) * share
            + (Self::lookup(&self.inner, phi) - 1.0) * (1.0 - share)
    }
}

/// Parameters of a loop-crater restore at the finishing resolution.
pub struct RestoreParams {
    pub time_gyr: f32,
    pub unit_floor: u16,
    pub unit_ejecta: u16,
    pub azimuth: f64,
    pub angle: f64,
    pub ellipticity: f64,
    /// The share of the wall's terraces, the ridge rings and the peak an ice shell has
    /// left, the crater's short wavelengths (F70); the resolved profile's detail comes
    /// relaxed by band in its table.
    pub detail_weight: f64,
    /// Seed of the crater's rim scalloping and texture; zero for a smooth stamp.
    pub seed: u64,
    /// The crater excavated the deep bedrock.
    pub deep: bool,
    /// Degradation blur of the resolved profile and of its coarse-smoothed version, metres.
    pub fine_blur: f64,
    pub coarse_blur: f64,
    /// The exposure on which the surface's maturity saturates, Gyr (F51).
    pub maturity_timescale: f64,
}

/// The share of the ground an ejecta blanket covers at the elliptical radius `re` of a
/// crater of radius `radius` with its continuous edge at `cont` (F51): the deposit thins
/// as r^EJECTA_EXPONENT and the cover falls as 1 - 2^-(cont / r)^3, complete out to 1.3
/// radii, half at the edge, tapering to nothing over the outer part of the footprint
/// `r_max`, where the deposit is too thin to leave a record.
pub fn ejecta_coverage(re: f64, radius: f64, cont: f64, r_max: f64) -> f64 {
    if re <= radius {
        return 1.0;
    }
    let raw = 1.0 - (-std::f64::consts::LN_2 * (cont / re).powf(-EJECTA_EXPONENT)).exp();
    let cover = (raw / EJECTA_FULL_COVER).min(1.0);
    let taper = ((r_max - re) / (EJECTA_COVERAGE_TAPER * r_max)).clamp(0.0, 1.0);
    cover * taper
}

/// The surface age a cell carries after a deposit of the event's time covered the share
/// `cover` of it (F51): the one age whose exposure gives the mixture's maturity under a
/// law that saturates on `tau`, the covered ground keeping its own exposure; exact for
/// an exponential law at every later time. All covered: the event's time; none: the old
/// age.
pub fn mixed_age(event: f32, old: f32, cover: f64, tau: f64) -> f32 {
    if cover >= 1.0 {
        return event;
    }
    if cover <= 0.0 || event <= old {
        return old;
    }
    let tau = tau.max(1e-6);
    let gap = (event - old) as f64;
    (event as f64 + tau * (cover + (1.0 - cover) * (-gap / tau).exp()).ln()) as f32
}

/// Relief under ejecta of thickness `thickness` is damped by this factor: depressions fill
/// and highs are shaken down, with the thickness measured against the relief itself so
/// thin ejecta on deep relief changes little and thick ejecta on shallow relief erases it.
#[inline]
pub fn burial_factor(relief: f64, thickness: f64) -> f64 {
    if thickness <= 0.0 {
        return 1.0;
    }
    let k = if relief < 0.0 {
        BURIAL_FACTOR_LOW
    } else {
        BURIAL_FACTOR_HIGH
    };
    (-k * thickness / relief.abs().max(1.0)).exp()
}

/// Re-stamp a loop crater at the finishing resolution, in its place in the age order.
///
/// `base` is the resampled coarse height before any finishing-pass stamp; it already holds
/// the crater's coarse form. Inside the rim the height returns to the base plus the
/// resolved detail, which removes whatever the pass stamped there earlier: those craters
/// are older and the impact destroyed them. Under the ejecta the relief added earlier in
/// the pass is buried according to the ejecta thickness. Ages and units within the
/// continuous ejecta are reset to the impact time.
#[allow(clippy::too_many_arguments)]
pub fn restore(
    layers: &mut LayerSet,
    grid: &CubeSphere,
    body_radius: f64,
    center: [f64; 3],
    base: &[f32],
    base_deep: &[f32],
    fresh: &ProfileTable,
    detail: &ProfileTable,
    p: &RestoreParams,
    stamper: &mut Stamper,
) -> usize {
    let out = restore_compute(
        layers,
        grid,
        body_radius,
        center,
        base,
        base_deep,
        fresh,
        detail,
        p,
        stamper,
    );
    out.apply(layers);
    out.cells.len()
}

/// Compute a loop crater's restore without writing it.
#[allow(clippy::too_many_arguments)]
pub fn restore_compute(
    layers: &LayerSet,
    grid: &CubeSphere,
    body_radius: f64,
    center: [f64; 3],
    base: &[f32],
    base_deep: &[f32],
    fresh: &ProfileTable,
    detail: &ProfileTable,
    p: &RestoreParams,
    stamper: &mut Stamper,
) -> StampOutput {
    let dim = &fresh.dim;
    stamper.gather(
        grid,
        body_radius,
        center,
        dim.ejecta_radius * p.ellipticity.sqrt(),
    );
    if stamper.footprint.is_empty() {
        return StampOutput::default();
    }
    let shape = RimShape::new(p.seed, dim, stamper.footprint.len() >= TEXTURE_MIN_CELLS);
    let texture = (p.seed != 0 && stamper.footprint.len() >= TEXTURE_MIN_CELLS)
        .then(|| CraterTexture::new(p.seed, dim));
    let terraces = (p.seed != 0 && stamper.footprint.len() >= TEXTURE_MIN_CELLS)
        .then(|| Terraces::new(p.seed, dim))
        .flatten();
    let rings = Rings::new(p.seed, dim);
    let peak = Peak::new(p.seed, dim);
    let pit = Pit::new(p.seed, dim);
    let rise = RimRise::new(p.seed, stamper.footprint.len() >= TEXTURE_MIN_CELLS);
    let edge = EdgeShape::new(p.seed, dim.ejecta_radius, grid.delta() * body_radius);
    let continuous = CONTINUOUS_EJECTA_MULTIPLE * dim.radius;
    let bedrock = layers.bedrock.as_slice();
    let ages = layers.surface_age.as_slice();
    let units = layers.unit_id.as_slice();
    let melts = layers.melt.as_slice();
    let covers = layers.cover.as_slice();
    let compute = |&(idx, r, phi): &(usize, f64, f64)| -> CellWrite {
        let re0 = elliptical_radius(r, phi, p.azimuth, p.ellipticity);
        let re = re0 * shape.stretch(re0 / dim.radius, phi, dim);
        let fine = bedrock[idx] as f64 - base[idx] as f64;
        // 1 inside 0.95 R, 0 outside 1.05 R: the cavity removes the fine relief.
        let cut = ((1.05 - re / dim.radius) / 0.10).clamp(0.0, 1.0);
        let thickness = if re > dim.radius {
            let t = fresh.height(re) * ejecta_factor(re, dim.radius, phi, p.azimuth, p.angle);
            match &texture {
                Some(tx) => tx.height(re, phi, t).max(0.0),
                None => t.max(0.0),
            }
        } else {
            0.0
        };
        let kept = fine * (1.0 - cut) * burial_factor(fine, thickness);
        let rw = terraces
            .as_ref()
            .map_or(re, |t| t.radius(re, phi, dim, p.fine_blur));
        let mut d = detail.height(re) * ejecta_factor(re, dim.radius, phi, p.azimuth, p.angle);
        if let Some(t) = &texture {
            d = t.height(rw, phi, d);
        }
        if let Some(rr) = &rise
            && d > 0.0
        {
            d *= rr.at(re0 / dim.radius, phi);
        }
        // The terraces the coarse grid could not carry: the wall as it collapsed, less the
        // plain wall the loop drew (F67).
        if let Some(t) = &terraces
            && re < dim.radius
        {
            d += (fresh.height(rw) - fresh.height(re) + t.relief(re, phi, dim)) * p.detail_weight;
        }
        if let Some(rg) = &rings {
            // The ring detail the coarse grid could not carry.
            d += (rg.relief(re0, phi, p.fine_blur) - rg.relief(re0, phi, p.coarse_blur))
                * p.detail_weight;
        }
        if let Some(pk) = &peak {
            let x = re0 / dim.radius;
            d += (pk.relief(x, phi, p.fine_blur / dim.radius)
                - pk.relief(x, phi, p.coarse_blur / dim.radius))
                * p.detail_weight;
        }
        if let Some(pt) = &pit {
            let x = re0 / dim.radius;
            d += (pt.relief(x, phi, p.fine_blur / dim.radius)
                - pt.relief(x, phi, p.coarse_blur / dim.radius))
                * p.detail_weight;
        }
        let mut write = CellWrite {
            idx: idx as u32,
            bedrock: (base[idx] as f64 + kept + d) as f32,
            surface: false,
            age: 0.0,
            unit: 0,
            deep: 0.0,
            melt: 0.0,
            cover: 0.0,
        };
        let cont = continuous * edge.factor_at(phi, re / dim.radius);
        let cover = ejecta_coverage(re, dim.radius, cont, dim.ejecta_radius);
        // The surface is the crater's only where nothing younger, lava or a later
        // crater, has resurfaced the cell since; its relief is restored regardless.
        if cover > 0.0 && ages[idx] <= p.time_gyr + 1e-6 {
            // The loop mixed this crater's deposit into the ground already; here the
            // cell is made at least as fresh as that cover on mature ground leaves it
            // (F51), which redraws the edge at this resolution without mixing twice.
            let floor_age =
                p.time_gyr as f64 + p.maturity_timescale.max(1e-6) * cover.min(1.0).ln();
            let thick = re <= dim.radius || cover >= 0.5;
            write.surface = true;
            write.age = ages[idx].max(floor_age as f32);
            write.unit = if re <= dim.radius {
                p.unit_floor
            } else if thick {
                p.unit_ejecta
            } else {
                units[idx]
            };
            // The loop already laid down this crater's deep material; put it back over
            // whatever older craters left, and start the melt afresh.
            write.deep = base_deep[idx];
            write.melt = if thick {
                texture.as_ref().map_or(0.0, |t| t.melt_fraction(rw, phi)) as f32
            } else {
                melts[idx]
            };
            write.cover = if !thick {
                covers[idx]
            } else if crate::stages::volcanism::is_lava_unit(write.unit) {
                1.0
            } else {
                0.0
            };
        }
        write
    };
    StampOutput {
        cells: map_footprint(&stamper.footprint, compute),
        composition: None,
        stain: None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::crater::morphology::{Scaling, dimensions};
    use crate::fields::Field;

    #[test]
    fn a_stamped_fresh_crater_has_the_right_depth_and_resets_ages_under_its_ejecta() {
        let grid = CubeSphere::new(256);
        let body_radius = 1.7374e6;
        let mut layers = LayerSet::new(grid);
        layers.bedrock = Field::filled(grid, 100.0);
        layers.surface_age = Field::filled(grid, 0.0);
        let s = Scaling {
            simple_complex: 16_500.0,
            peak_ring: 140_000.0,
            multi_ring: 330_000.0,
            gravity: 1.62,
            icy: false,
        };
        let dim = dimensions(200_000.0, &s);
        let table = ProfileTable::new(dim, 0.0);
        let mut stamper = Stamper::new(&grid);
        let center = crate::geom::normalize([0.3, 0.4, 0.86]);
        let params = StampParams {
            time_gyr: 3.0,
            unit_floor: 3,
            unit_ejecta: 4,
            azimuth: 0.0,
            angle: std::f64::consts::FRAC_PI_2,
            ellipticity: 1.0,
            amplitude: 1.0,
            feature_amplitude: 1.0,
            replace_floor: true,
            set_age: true,
            bury: true,
            seed: 0,
            deep: false,
            maturity_timescale: 0.5,
            renewal_guard: true,
            terraced: false,
        };
        let touched = stamp(
            &mut layers,
            &grid,
            body_radius,
            center,
            &table,
            &params,
            &mut stamper,
        );
        assert!(touched > 100);
        let (lo, hi) = layers.bedrock.min_max();
        let expected_depth = dim.depth + dim.rim_height;
        assert!(
            ((hi - lo) as f64 - expected_depth).abs() < 0.1 * expected_depth,
            "measured {} expected {expected_depth}",
            hi - lo
        );
        // The floor sits at the reference height minus the depth.
        assert!(
            (lo as f64 - (100.0 - dim.depth)).abs() < 0.05 * dim.depth,
            "floor {lo}"
        );
        // Ages reset under the continuous ejecta but not far away.
        let centre_cell = grid.locate(center).0;
        assert_eq!(layers.surface_age.get(centre_cell), 3.0);
        assert_eq!(layers.unit_id.get(centre_cell), 3);
        let far = grid.locate(crate::geom::normalize([-0.9, 0.1, -0.3])).0;
        assert_eq!(layers.surface_age.get(far), 0.0);
        // Cells covered by at least half (F51): a disc of 2 R, the continuous edge, by
        // the age that cover leaves on ground of age 0 with the test's timescale.
        let threshold = mixed_age(3.0, 0.0, 0.5 / EJECTA_FULL_COVER, 0.5);
        let reset: f64 = grid
            .cells()
            .filter(|c| layers.surface_age.get(*c) >= threshold - 1e-4)
            .map(|c| grid.cell_area(c))
            .sum();
        let half = 2.0 * dim.radius / body_radius;
        let cap = 2.0 * std::f64::consts::PI * (1.0 - half.cos());
        assert!((reset - cap).abs() < 0.08 * cap, "reset {reset} cap {cap}");
    }

    /// Ejecta buries older relief: a small crater under a larger crater's ejecta loses a
    /// good part of its depth; without burial the ejecta only drapes it.
    #[test]
    fn ejecta_buries_older_relief() {
        let grid = CubeSphere::new(256);
        let body_radius = 1.7374e6;
        let s = Scaling {
            simple_complex: 16_500.0,
            peak_ring: 140_000.0,
            multi_ring: 330_000.0,
            gravity: 1.62,
            icy: false,
        };
        let small = dimensions(20_000.0, &s);
        let big = dimensions(200_000.0, &s);
        let big_center = [0.0, 0.0, 1.0];
        let (east, north) = crate::geom::local_basis(big_center);
        // The small crater sits at 1.3 big radii, under the continuous ejecta.
        let small_center =
            crate::geom::offset(big_center, east, north, 1.3 * big.radius / body_radius, 0.3);
        let relief = |layers: &LayerSet, stamper: &mut Stamper| -> f64 {
            stamper.gather(&grid, body_radius, small_center, 1.1 * small.radius);
            let h = layers.bedrock.as_slice();
            let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
            for &(idx, _, _) in &stamper.footprint {
                lo = lo.min(h[idx] as f64);
                hi = hi.max(h[idx] as f64);
            }
            hi - lo
        };
        let mut results = Vec::new();
        for bury in [true, false] {
            let mut layers = LayerSet::new(grid);
            let mut stamper = Stamper::new(&grid);
            let params = |t: f32| StampParams {
                time_gyr: t,
                unit_floor: 4,
                unit_ejecta: 5,
                azimuth: 0.0,
                angle: std::f64::consts::FRAC_PI_2,
                ellipticity: 1.0,
                amplitude: 1.0,
                feature_amplitude: 1.0,
                replace_floor: true,
                set_age: true,
                bury,
                seed: 0,
                deep: false,
                maturity_timescale: 0.5,
                renewal_guard: true,
                terraced: false,
            };
            stamp(
                &mut layers,
                &grid,
                body_radius,
                small_center,
                &ProfileTable::new(small, 0.0),
                &params(1.0),
                &mut stamper,
            );
            let before = relief(&layers, &mut stamper);
            stamp(
                &mut layers,
                &grid,
                body_radius,
                big_center,
                &ProfileTable::new(big, 0.0),
                &params(2.0),
                &mut stamper,
            );
            let after = relief(&layers, &mut stamper);
            results.push((before, after));
        }
        let (buried_before, buried_after) = results[0];
        let (draped_before, draped_after) = results[1];
        assert!((buried_before - draped_before).abs() < 1.0);
        // Draping adds the ejecta slope across the span and nothing else is lost.
        assert!(
            draped_after > 0.9 * draped_before,
            "draped {draped_before} -> {draped_after}"
        );
        // Burial removes a good part of the relief on top of the same slope.
        let lost = draped_after - buried_after;
        assert!(
            lost > 0.2 * draped_before && lost < 0.8 * draped_before,
            "buried {buried_before} -> {buried_after}, draped -> {draped_after}"
        );
    }

    /// Oblique ejecta is thicker downrange, nearly absent uprange and raised to the sides;
    /// a vertical impact is symmetric. Measured against the analytic thickness at each
    /// sampled cell's own distance, so the grid does not enter.
    #[test]
    fn oblique_ejecta_is_thicker_downrange_and_absent_uprange() {
        let grid = CubeSphere::new(512);
        let body_radius = 1.7374e6;
        let s = Scaling {
            simple_complex: 16_500.0,
            peak_ring: 140_000.0,
            multi_ring: 330_000.0,
            gravity: 1.62,
            icy: false,
        };
        let dim = dimensions(200_000.0, &s);
        let table = ProfileTable::new(dim, 0.0);
        let center = [0.0, 0.0, 1.0];
        let (east, north) = crate::geom::local_basis(center);
        let azimuth = 1.0;
        let factor = |layers: &LayerSet, phi: f64| -> f64 {
            let d = crate::geom::offset(center, east, north, 1.5 * dim.radius / body_radius, phi);
            let cell = grid.locate(d).0;
            let r = crate::geom::angle_between(center, grid.center(cell)) * body_radius;
            layers.bedrock.get(cell) as f64 / table.height(r)
        };
        let mut out = Vec::new();
        for angle in [std::f64::consts::FRAC_PI_2, 10f64.to_radians()] {
            let mut layers = LayerSet::new(grid);
            let mut stamper = Stamper::new(&grid);
            let params = StampParams {
                time_gyr: 1.0,
                unit_floor: 4,
                unit_ejecta: 5,
                azimuth,
                angle,
                ellipticity: 1.0,
                amplitude: 1.0,
                feature_amplitude: 1.0,
                replace_floor: true,
                set_age: true,
                bury: true,
                seed: 0,
                deep: false,
                maturity_timescale: 0.5,
                renewal_guard: true,
                terraced: false,
            };
            stamp(
                &mut layers,
                &grid,
                body_radius,
                center,
                &table,
                &params,
                &mut stamper,
            );
            out.push((
                factor(&layers, azimuth),
                factor(&layers, azimuth + std::f64::consts::PI),
                factor(&layers, azimuth + std::f64::consts::FRAC_PI_2),
            ));
        }
        let (vd, vu, vs) = out[0];
        for v in [vd, vu, vs] {
            assert!((v - 1.0).abs() < 0.05, "vertical impact: {vd} {vu} {vs}");
        }
        let (od, ou, os) = out[1];
        assert!(
            od > 1.3 && ou < 0.4 && os > 1.1,
            "oblique: down {od} up {ou} side {os}"
        );
    }

    /// The floor texture must vary with radius as well as azimuth: a lattice that clamps
    /// the floor to one row turns it into spokes from the centre.
    #[test]
    fn floor_texture_varies_with_radius_and_melt_is_soft() {
        let s = Scaling {
            simple_complex: 16_500.0,
            peak_ring: 140_000.0,
            multi_ring: 330_000.0,
            gravity: 1.62,
            icy: false,
        };
        let dim = dimensions(90_000.0, &s);
        let tex = CraterTexture::new(0x5EED, &dim);
        let mut same = 0;
        let n = 64;
        for b in 0..n {
            let phi = b as f64 / n as f64 * 2.0 * std::f64::consts::PI;
            let a = tex.noise(phi, 0.1);
            let c = tex.noise(phi, 0.35);
            if (a - c).abs() < 0.05 {
                same += 1;
            }
        }
        assert!(
            same < n / 4,
            "{same} of {n} azimuths have the same texture at two radii"
        );
        // Melt is a soft fraction: intermediate values exist, and none outside the floor
        // and near-rim bands.
        let mut partial = 0;
        for b in 0..n {
            let phi = b as f64 / n as f64 * 2.0 * std::f64::consts::PI;
            for x in [0.1, 0.25, 0.4] {
                let m = tex.melt_fraction(x * dim.radius, phi);
                assert!((0.0..=1.0).contains(&m));
                if m > 0.05 && m < 0.95 {
                    partial += 1;
                }
            }
            assert_eq!(tex.melt_fraction(2.0 * dim.radius, phi), 0.0);
        }
        assert!(partial > 0, "no soft pond edges");
        // The floor texture is continuous through the centre: all azimuths agree there.
        let centre: Vec<f64> = (0..n)
            .map(|b| tex.noise_floor(b as f64 / n as f64 * 2.0 * std::f64::consts::PI, 0.01))
            .collect();
        let spread = centre.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
            - centre.iter().cloned().fold(f64::INFINITY, f64::min);
        assert!(spread < 0.1, "floor texture spread at the centre {spread}");
    }

    /// A complex crater's central peak rises from its floor by most of the Pike height.
    #[test]
    fn central_peak_rises_from_the_floor() {
        let grid = CubeSphere::new(1024);
        let body_radius = 1.7374e6;
        let s = Scaling {
            simple_complex: 16_500.0,
            peak_ring: 140_000.0,
            multi_ring: 330_000.0,
            gravity: 1.62,
            icy: false,
        };
        let dim = dimensions(60_000.0, &s);
        assert!(dim.peak_height > 500.0);
        let mut layers = LayerSet::new(grid);
        let mut stamper = Stamper::new(&grid);
        let center = [0.0, 0.0, 1.0];
        let params = StampParams {
            time_gyr: 1.0,
            unit_floor: 4,
            unit_ejecta: 5,
            azimuth: 0.0,
            angle: std::f64::consts::FRAC_PI_2,
            ellipticity: 1.0,
            amplitude: 1.0,
            feature_amplitude: 1.0,
            replace_floor: true,
            set_age: true,
            bury: true,
            seed: 0,
            deep: false,
            maturity_timescale: 0.5,
            renewal_guard: true,
            terraced: false,
        };
        stamp(
            &mut layers,
            &grid,
            body_radius,
            center,
            &ProfileTable::new(dim, 0.0),
            &params,
            &mut stamper,
        );
        stamper.gather(&grid, body_radius, center, 0.4 * dim.radius);
        let h = layers.bedrock.as_slice();
        let top = stamper
            .footprint
            .iter()
            .map(|&(idx, _, _)| h[idx] as f64)
            .fold(f64::NEG_INFINITY, f64::max);
        let rise = top + dim.depth;
        assert!(
            rise > 0.4 * dim.peak_height,
            "peak rises {rise} m of {} m",
            dim.peak_height
        );
    }

    #[test]
    fn footprint_handles_face_edges_and_poles() {
        let grid = CubeSphere::new(32);
        let mut stamper = Stamper::new(&grid);
        for center in [
            [0.0, 0.0, 1.0],
            crate::geom::normalize([1.0, 1.0, 0.0]),
            crate::geom::normalize([1.0, 1.0, 1.0]),
        ] {
            stamper.gather(&grid, 1.0e6, center, 0.3e6);
            let n = stamper.footprint.len() as f64;
            let expected = 2.0 * std::f64::consts::PI * (1.0 - 0.3f64.cos())
                / (4.0 * std::f64::consts::PI)
                * grid.cell_count() as f64;
            assert!(
                (n - expected).abs() < 0.25 * expected + 20.0,
                "{n} vs {expected}"
            );
        }
    }

    #[test]
    fn a_complex_crater_wall_is_terraced_and_differs_around_it() {
        // F67: a lunar crater the size of Copernicus.
        let moon = Scaling {
            simple_complex: 16_457.0,
            peak_ring: 140_000.0,
            multi_ring: 330_000.0,
            gravity: 1.62,
            icy: false,
        };
        let dim = dimensions(93e3, &moon);
        // The outermost terrace is Pearce and Melosh's 0.09 D^0.87 km.
        let expected = 0.09 * 93f64.powf(0.87);
        assert!(
            (dim.terrace_width / 1e3 - expected).abs() < 0.01,
            "{}",
            dim.terrace_width
        );
        // Mercury's, at its smaller strength length, is a tenth narrower.
        let mercury = Scaling {
            simple_complex: 7_221.0,
            ..moon
        };
        let ratio = dimensions(93e3, &mercury).terrace_width / dim.terrace_width;
        assert!((0.88..0.92).contains(&ratio), "{ratio}");
        let t = Terraces::new(7, &dim).expect("a terraced wall");
        assert!(t.levels.len() >= 3, "{}", t.levels.len());
        let table = ProfileTable::new(dim, 0.0);
        let r = dim.radius;
        let relief = dim.depth + dim.rim_height;
        let n = 800;
        let (mut edges, mut firsts, mut benches) = (Vec::new(), Vec::new(), 0usize);
        for a in 0..16 {
            let phi = a as f64 * std::f64::consts::TAU / 16.0;
            let heights: Vec<f64> = (0..=n)
                .map(|i| table.height(t.radius(i as f64 / n as f64 * r, phi, &dim, 0.0)))
                .collect();
            // The wall never falls on the way out to the crest (the profile's table, not the
            // terraces, rounds the crest itself between its samples).
            for w in heights[..(0.97 * n as f64) as usize].windows(2) {
                assert!(w[1] >= w[0] - 1e-6 * relief, "{a}");
            }
            let edge = heights
                .iter()
                .position(|&h| h > -dim.depth + 0.01 * relief)
                .expect("a wall");
            edges.push(edge as f64 / n as f64);
            // Benches: runs of the wall that climb at under 0.6 of its mean slope.
            let mean = relief / (n - edge) as f64;
            let mut run = 0;
            let mut first = None;
            for (i, w) in heights.windows(2).enumerate().skip(edge) {
                if w[1] - w[0] < 0.6 * mean {
                    run += 1;
                    if run == 4 {
                        benches += 1;
                        first.get_or_insert(i as f64 / n as f64);
                    }
                } else {
                    run = 0;
                }
            }
            firsts.push(first.unwrap_or(0.0));
        }
        let spread = |v: &[f64]| {
            v.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
                - v.iter().cloned().fold(f64::INFINITY, f64::min)
        };
        // Two benches or more on most radii, and neither the floor's edge nor the benches
        // copy one another round the crater.
        assert!(benches >= 2 * 16, "{benches}");
        assert!(spread(&edges) > 0.01, "{edges:?}");
        assert!(spread(&firsts) > 0.01, "{firsts:?}");
        // The rim crest is scalloped by the outermost slumps.
        let rim = RimShape::new(7, &dim, true);
        let (lo, hi) = rim
            .table
            .iter()
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), &v| {
                (lo.min(v), hi.max(v))
            });
        assert!(hi - lo > 0.02, "{lo} {hi}");
    }
}
