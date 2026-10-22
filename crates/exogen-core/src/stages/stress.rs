// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! The lid's stress field (F150): one field, in pascals, that every fault on a rock lid
//! reads. It is the sum of the rises' stresses, their domes' bending and the membrane
//! stress of their uplift on the sphere (F105, F114), of the upwelling between neighbouring
//! plumes, which pulls the lid apart across the ridge of the mantle's convection joining
//! them (F173), and of the
//! potential energy of the relief: a high stands in tension, a low in compression, each
//! alike every way within and across the slope at its margin, where a plateau spreads down
//! its slope (Bilotti and Suppe 1999). A fault forms where the most tensile stress passes
//! what breaks the brittle layer, and runs as long as it does. A thrust reads the same
//! field from its other side (F154): it strikes across the most compressive stress, and it
//! forms where the field's compression, added to what the body's shrinking puts on the lid
//! alike every way, passes what breaks the layer in compression. Through the loop the field
//! is read at epochs (F157): the relief as the ground stands then, carried by the lid as
//! thick as it is then, and the stresses of the loads the floods have laid, each bending the
//! lid of its own time.

use rayon::prelude::*;

use crate::constants::*;
use crate::convection::Support;
use crate::geom::{angle_between, cross, dot, local_basis, normalize};
use crate::grid::CubeSphere;
use crate::params::ParamGraph;
use crate::spharm::LatLon;
use crate::stages::load::LoadStress;
use crate::stages::plume::{Axial, FieldAt, Rise, dome, rise_progress, total_stress};

/// The upwelling limbs between neighbouring plumes (F150): the pairs of rises they join, by
/// index, the width the limb's pull spreads over, m, and where the mantle convects each
/// limb's way along the ridge of its pattern (F173), evenly spaced points from one plume's
/// centre to the other's; without them a limb runs along the great circle.
#[derive(Clone, Debug, Default)]
pub struct Limbs {
    pub pairs: Vec<(usize, usize)>,
    pub width: f64,
    pub paths: Vec<Vec<[f64; 3]>>,
}

impl Limbs {
    /// The limbs where the body's mantle convects (F173): the ridges of its pattern as the
    /// plumes rose from it, each joining two plumes whose ridge stands above
    /// UPWELLING_LIMB_SADDLE of the pattern all the way between them, along the way that
    /// keeps highest; where it keeps no convection, the great circles between plumes that
    /// near (`between`).
    pub fn along(rises: &[Rise], support: Option<&Support>, body_radius: f64) -> Self {
        let Some(support) = support else {
            return Self::between(rises, body_radius);
        };
        let width = rises.iter().map(|r| r.radius).fold(0.0, f64::max) * UPWELLING_LIMB_WIDTH;
        let mut limbs = Self {
            pairs: Vec::new(),
            width,
            paths: Vec::new(),
        };
        if rises.len() < 2 {
            return limbs;
        }
        let start = rises
            .iter()
            .map(|r| r.start_gyr)
            .fold(f64::INFINITY, f64::min);
        let (nlat, nlon) = (90usize, 180usize);
        let grid = LatLon {
            nlat,
            nlon,
            values: Vec::new(),
        };
        let nodes: Vec<[f64; 3]> = (0..nlat * nlon)
            .map(|n| grid.direction(n / nlon, n % nlon))
            .collect();
        let values: Vec<f64> = nodes
            .iter()
            .map(|&d| support.convection.pattern(start, d))
            .collect();
        let nearest = |d: [f64; 3]| {
            (0..nodes.len())
                .max_by(|&a, &b| dot(nodes[a], d).total_cmp(&dot(nodes[b], d)))
                .unwrap_or(0)
        };
        let spacing = (0.25 * width / body_radius).max(1e-3);
        for a in 0..rises.len() {
            for b in a + 1..rises.len() {
                let from = nearest(support.turn.undo(rises[a].center));
                let to = nearest(support.turn.undo(rises[b].center));
                let Some(way) = ridge_path(&nodes, &values, nlon, from, to) else {
                    continue;
                };
                if way.iter().any(|&n| values[n] < UPWELLING_LIMB_SADDLE) {
                    continue;
                }
                // The way in the map's frame from one plume's centre to the other's,
                // smoothed off the grid's steps and evenly spaced.
                let mut path: Vec<[f64; 3]> =
                    way.iter().map(|&n| support.turn.apply(nodes[n])).collect();
                path[0] = rises[a].center;
                let last = path.len() - 1;
                path[last] = rises[b].center;
                for _ in 0..4 {
                    let before = path.clone();
                    for k in 1..last {
                        let (p, q, r) = (before[k - 1], before[k], before[k + 1]);
                        path[k] = normalize([
                            p[0] + 2.0 * q[0] + r[0],
                            p[1] + 2.0 * q[1] + r[1],
                            p[2] + 2.0 * q[2] + r[2],
                        ]);
                    }
                }
                limbs.pairs.push((a, b));
                limbs.paths.push(resample(&path, spacing));
            }
        }
        limbs
    }

    /// Each limb's length, m: along its way, or along the great circle between its plumes.
    pub fn lengths(&self, rises: &[Rise], body_radius: f64) -> Vec<f64> {
        self.pairs
            .iter()
            .enumerate()
            .map(|(k, &(a, b))| match self.paths.get(k) {
                Some(path) => {
                    path.windows(2)
                        .map(|w| angle_between(w[0], w[1]))
                        .sum::<f64>()
                        * body_radius
                }
                None => angle_between(rises[a].center, rises[b].center) * body_radius,
            })
            .collect()
    }

    pub fn between(rises: &[Rise], body_radius: f64) -> Self {
        let mut pairs = Vec::new();
        for a in 0..rises.len() {
            for b in a + 1..rises.len() {
                let apart = angle_between(rises[a].center, rises[b].center) * body_radius;
                if apart <= UPWELLING_LIMB_REACH * (rises[a].radius + rises[b].radius) {
                    pairs.push((a, b));
                }
            }
        }
        let width = rises.iter().map(|r| r.radius).fold(0.0, f64::max) * UPWELLING_LIMB_WIDTH;
        Self {
            pairs,
            width,
            paths: Vec::new(),
        }
    }

    /// How much of a limb's pull a point feels at a time, 0 to 1, and the direction across
    /// the limb there: a Gaussian of the distance from the great circle between the two
    /// rises, over the arc between them and fading beyond their centres over the width,
    /// growing as the lesser of the two plumes rises.
    pub fn pull(&self, rises: &[Rise], p: [f64; 3], t: f64, body_radius: f64) -> (f64, [f64; 3]) {
        let mut best = (0.0, [0.0; 3]);
        if !self.paths.is_empty() {
            // Along each limb's way: a Gaussian of the distance from it, the ends its
            // plumes' centres, and across it the way's pole there.
            let w = self.width.max(1.0);
            for (k, &(a, b)) in self.pairs.iter().enumerate() {
                let grown = rise_progress(&rises[a], t).min(rise_progress(&rises[b], t));
                if grown <= best.0 {
                    continue;
                }
                let (off, across) = nearest_on(&self.paths[k], p);
                let pull = (-(off * body_radius / w).powi(2)).exp() * grown;
                if pull > best.0 {
                    best = (pull, across);
                }
            }
            return best;
        }
        for &(a, b) in &self.pairs {
            let (ca, cb) = (rises[a].center, rises[b].center);
            let pole = normalize(cross(ca, cb));
            // The distance off the great circle, and where along the arc the point lies.
            let off = dot(p, pole).clamp(-1.0, 1.0).asin() * body_radius;
            let foot = normalize([
                p[0] - dot(p, pole) * pole[0],
                p[1] - dot(p, pole) * pole[1],
                p[2] - dot(p, pole) * pole[2],
            ]);
            // The arc runs from a to b counterclockwise about the pole; the foot's angle
            // from a along it, and how far it lies beyond either end.
            let arc = angle_between(ca, cb);
            let phi = dot(cross(ca, foot), pole).atan2(dot(ca, foot));
            let beyond = if phi < 0.0 {
                -phi * body_radius
            } else if phi > arc {
                (phi - arc) * body_radius
            } else {
                0.0
            };
            let w = self.width.max(1.0);
            let share = (-(off / w).powi(2)).exp() * (-(beyond / w).powi(2)).exp();
            let grown = rise_progress(&rises[a], t).min(rise_progress(&rises[b], t));
            let pull = share * grown;
            if pull > best.0 {
                // Across the limb: the pole's direction in the tangent plane.
                let r = dot(pole, p);
                let across = [pole[0] - r * p[0], pole[1] - r * p[1], pole[2] - r * p[2]];
                let across = if dot(across, across) < 1e-20 {
                    local_basis(p).0
                } else {
                    normalize(across)
                };
                best = (pull, across);
            }
        }
        best
    }
}

/// The way across a grid of the convection's pattern between two nodes that keeps highest
/// (F173): the least cost by Dijkstra's search over each node's eight neighbours, each
/// step's arc weighed by exp(−UPWELLING_LIMB_RIDGE × the pattern where it ends).
fn ridge_path(
    nodes: &[[f64; 3]],
    values: &[f64],
    nlon: usize,
    from: usize,
    to: usize,
) -> Option<Vec<usize>> {
    use std::cmp::Reverse;
    use std::collections::BinaryHeap;
    let n = nodes.len();
    let nlat = n / nlon;
    let mut cost = vec![f64::INFINITY; n];
    let mut prev = vec![usize::MAX; n];
    let mut heap = BinaryHeap::new();
    cost[from] = 0.0;
    // Costs are never negative, so their bits order as they do.
    heap.push(Reverse((0.0f64.to_bits(), from)));
    while let Some(Reverse((bits, u))) = heap.pop() {
        let c = f64::from_bits(bits);
        if c > cost[u] {
            continue;
        }
        if u == to {
            break;
        }
        let (i, j) = (u / nlon, u % nlon);
        for di in -1i64..=1 {
            let ii = i as i64 + di;
            if ii < 0 || ii >= nlat as i64 {
                continue;
            }
            for dj in -1i64..=1 {
                if di == 0 && dj == 0 {
                    continue;
                }
                let jj = (j as i64 + dj).rem_euclid(nlon as i64) as usize;
                let v = ii as usize * nlon + jj;
                let step =
                    angle_between(nodes[u], nodes[v]) * (-UPWELLING_LIMB_RIDGE * values[v]).exp();
                if c + step < cost[v] {
                    cost[v] = c + step;
                    prev[v] = u;
                    heap.push(Reverse(((c + step).to_bits(), v)));
                }
            }
        }
    }
    if !cost[to].is_finite() {
        return None;
    }
    let mut way = vec![to];
    let mut k = to;
    while k != from {
        k = prev[k];
        if k == usize::MAX {
            return None;
        }
        way.push(k);
    }
    way.reverse();
    Some(way)
}

/// A path of unit directions evenly spaced, its points `spacing` rad apart along it, its
/// ends kept.
fn resample(path: &[[f64; 3]], spacing: f64) -> Vec<[f64; 3]> {
    let lengths: Vec<f64> = path.windows(2).map(|w| angle_between(w[0], w[1])).collect();
    let total: f64 = lengths.iter().sum();
    let n = ((total / spacing).ceil() as usize).max(1);
    let mut out = Vec::with_capacity(n + 1);
    let (mut k, mut walked) = (0usize, 0.0);
    for s in 0..=n {
        let at = total * s as f64 / n as f64;
        while k + 1 < lengths.len() && walked + lengths[k] < at {
            walked += lengths[k];
            k += 1;
        }
        let f = if lengths.is_empty() || lengths[k] <= 0.0 {
            0.0
        } else {
            ((at - walked) / lengths[k]).clamp(0.0, 1.0)
        };
        let (a, b) = (path[k], path[(k + 1).min(path.len() - 1)]);
        out.push(normalize([
            a[0] + f * (b[0] - a[0]),
            a[1] + f * (b[1] - a[1]),
            a[2] + f * (b[2] - a[2]),
        ]));
    }
    out
}

/// A point's angular distance from a path of unit directions, rad, and the direction across
/// the path there, the pole of its nearest piece in the point's tangent plane.
fn nearest_on(path: &[[f64; 3]], p: [f64; 3]) -> (f64, [f64; 3]) {
    let mut best = (f64::INFINITY, local_basis(p).0);
    for w in path.windows(2) {
        let (a, b) = (w[0], w[1]);
        let c = cross(a, b);
        let size = dot(c, c).sqrt();
        if size < 1e-12 {
            continue;
        }
        let pole = [c[0] / size, c[1] / size, c[2] / size];
        let off = dot(p, pole);
        let foot = normalize([
            p[0] - off * pole[0],
            p[1] - off * pole[1],
            p[2] - off * pole[2],
        ]);
        let phi = dot(cross(a, foot), pole).atan2(dot(a, foot));
        let d = if phi >= 0.0 && phi <= angle_between(a, b) {
            off.abs().min(1.0).asin()
        } else {
            angle_between(p, a).min(angle_between(p, b))
        };
        if d < best.0 {
            let across = [
                pole[0] - off * p[0],
                pole[1] - off * p[1],
                pole[2] - off * p[2],
            ];
            let across = if dot(across, across) < 1e-20 {
                local_basis(p).0
            } else {
                normalize(across)
            };
            best = (d, across);
        }
    }
    best
}

/// The relief the potential-energy stress reads: the ground before the rises, low-passed
/// over the lid's flexural wavelength, on the loop's grid.
#[derive(Clone, Debug)]
pub struct Relief {
    pub grid: CubeSphere,
    pub height: Vec<f32>,
}

impl Relief {
    /// The relief at an epoch (F157): the ground as it stands less the rises' domes, which
    /// the field holds itself, low-passed over the flexural wavelength of a lid `elastic` m
    /// thick.
    pub fn read(
        bedrock: &crate::fields::Field<f32>,
        domes: Option<&[f32]>,
        cell: f64,
        elastic: f64,
        gravity: f64,
    ) -> Self {
        let mut ground = bedrock.clone();
        if let Some(d) = domes {
            for (g, h) in ground.as_mut_slice().iter_mut().zip(d) {
                *g -= h;
            }
        }
        let wavelength = StressField::flexural_wavelength(elastic, gravity);
        Self {
            grid: ground.grid(),
            height: crate::stages::relaxation::low_pass(&ground, cell, wavelength),
        }
    }
}

/// The lid's stress field on a body.
pub struct StressField<'a> {
    pub rises: &'a [Rise],
    axes: Vec<Axial>,
    /// Each rise's tension at its crest at full growth, Pa.
    crest: Vec<f64>,
    pub limbs: Limbs,
    /// The limb's pull at full growth, Pa, and the stress along it, a share of that.
    limb_tension: f64,
    relief: Option<&'a Relief>,
    /// The stresses the loads laid by then put on the lid (F157).
    loads: Vec<&'a LoadStress>,
    body_radius: f64,
    gravity: f64,
    /// The lid's elastic thickness, m, which carries the relief's force, and the length its
    /// margins spread over, its flexural wavelength, m.
    elastic: f64,
    flexural: f64,
    /// The crust's density and the mantle's it floats on, kg/m³, and Poisson's ratio.
    density: f64,
    mantle: f64,
    poisson: f64,
    /// The tension that breaks the brittle layer, Pa, and the compression (F154).
    pub yield_tension: f64,
    pub yield_compression: f64,
}

/// The lid's compression at a place and time as a thrust reads it (F154).
#[derive(Clone, Copy, Debug, Default)]
pub struct Squeeze {
    /// The field's most compressive stress, Pa, compression positive: negative where the
    /// field stands in tension every way.
    pub stress: f64,
    /// The field's least compressive stress, Pa, tension positive as the field's stresses
    /// are: a thrust breaks only where this too, with the shrinking's squeeze, presses
    /// harder than the overburden, which is then the least of the three, as Anderson's
    /// thrusts need; where it pulls the lid would slip sideways instead (F157).
    pub least: f64,
    /// The strike a thrust takes there: across the most compressive direction.
    pub strike: [f64; 3],
    /// How far the field decides that strike, 0 to 1: its differential stress against the
    /// scatter of the lid's own strength, which a thrust in a field alike every way follows
    /// instead.
    pub decided: f64,
}

impl<'a> StressField<'a> {
    /// The field over a rock lid of elastic thickness `elastic`, m, with the relief given.
    pub(crate) fn new(
        rises: &'a [Rise],
        axes: Vec<Axial>,
        body_radius: f64,
        gravity: f64,
        elastic: f64,
        relief: Option<&'a Relief>,
        loads: Vec<&'a LoadStress>,
    ) -> Self {
        let crest = rises
            .iter()
            .map(|r| RISE_MEMBRANE_SHARE * ROCK_YOUNG_MODULUS * r.height / body_radius)
            .collect::<Vec<f64>>();
        let limbs = Limbs::between(rises, body_radius);
        let least = crest.iter().copied().fold(f64::INFINITY, f64::min);
        let limb_tension = if crest.is_empty() {
            0.0
        } else {
            UPWELLING_LIMB_TENSION * least
        };
        let (density, mantle, poisson) = (DENSITY_CRUST_ROCK, DENSITY_ROCK, ROCK_POISSON);
        let rigidity =
            ROCK_YOUNG_MODULUS * elastic.max(1e3).powi(3) / (12.0 * (1.0 - poisson * poisson));
        let flexural = std::f64::consts::TAU * (rigidity / (mantle * gravity.max(1e-3))).powf(0.25);
        let k = ((1.0 + ROCK_FRICTION.powi(2)).sqrt() + ROCK_FRICTION).powi(2);
        let yield_tension = (1.0 - 1.0 / k) * density * gravity * FAULT_LAYER_DEPTH;
        let yield_compression = (k - 1.0) * density * gravity * FAULT_LAYER_DEPTH;
        Self {
            rises,
            axes,
            crest,
            limbs,
            limb_tension,
            relief,
            loads,
            body_radius,
            gravity,
            elastic: elastic.max(1e3),
            flexural,
            density,
            mantle,
            poisson,
            yield_tension,
            yield_compression,
        }
    }

    /// The field with these limbs between its plumes (F173) in place of the great circles.
    pub(crate) fn with_limbs(mut self, limbs: Limbs) -> Self {
        self.limbs = limbs;
        self
    }

    /// The flexural wavelength of a lid of elastic thickness `elastic`, m, the length the
    /// relief's potential energy reads it over.
    pub fn flexural_wavelength(elastic: f64, gravity: f64) -> f64 {
        let rigidity = ROCK_YOUNG_MODULUS * elastic.max(1e3).powi(3)
            / (12.0 * (1.0 - ROCK_POISSON * ROCK_POISSON));
        std::f64::consts::TAU * (rigidity / (DENSITY_ROCK * gravity.max(1e-3))).powf(0.25)
    }

    /// The relief at a point and time, m: the ground before the rises and the rises' domes
    /// as grown by then.
    fn height(&self, p: [f64; 3], t: f64) -> f64 {
        let base = match &self.relief {
            Some(r) => r.height[r.grid.index(r.grid.locate(p).0)] as f64,
            None => 0.0,
        };
        base + self
            .rises
            .iter()
            .map(|r| {
                let x = angle_between(r.center, p) * self.body_radius / r.radius;
                if x < 1.0 {
                    r.height * dome(x) * rise_progress(r, t)
                } else {
                    0.0
                }
            })
            .sum::<f64>()
    }

    /// The horizontal force per length the relief's potential energy puts on the lid at a
    /// point, N/m, tension positive: a compensated high of height h stands on a root, and
    /// its column holds ρ g h² / 2 × ρ_m / (ρ_m − ρ_c) more energy per area than the
    /// reference column (Turcotte and Schubert 2002, the force on a compensated plateau),
    /// which it spends spreading; a low holds less and is squeezed.
    fn potential_force(&self, p: [f64; 3], t: f64) -> f64 {
        let h = self.height(p, t);
        let factor = self.mantle / (self.mantle - self.density).max(1.0);
        h.signum() * 0.5 * self.density * self.gravity * h * h * factor
    }

    /// The stress tensor at a point and time in the local east-north frame, Pa, tension
    /// positive: (σ_ee, σ_nn, σ_en).
    fn tensor(&self, p: [f64; 3], t: f64) -> (f64, f64, f64) {
        let (east, north) = local_basis(p);
        let mut s = (0.0, 0.0, 0.0);
        // A tensor given in a frame (u, v), added in the east-north frame.
        let mut add = |u: [f64; 3], v: [f64; 3], suu: f64, svv: f64, suv: f64| {
            let (ue, un, ve, vn) = (dot(u, east), dot(u, north), dot(v, east), dot(v, north));
            s.0 += suu * ue * ue + 2.0 * suv * ue * ve + svv * ve * ve;
            s.1 += suu * un * un + 2.0 * suv * un * vn + svv * vn * vn;
            s.2 += suu * ue * un + suv * (ue * vn + ve * un) + svv * ve * vn;
        };
        for (i, r) in self.rises.iter().enumerate() {
            let scale = self.crest[i] * rise_progress(r, t);
            if scale <= 0.0 {
                continue;
            }
            let (radial, hoop, shear, out, round) =
                total_stress(r, Some(self.axes[i]), p, self.body_radius);
            add(out, round, radial * scale, hoop * scale, shear * scale);
        }
        if self.limb_tension > 0.0 {
            let (pull, across) = self.limbs.pull(self.rises, p, t, self.body_radius);
            if pull > 0.0 {
                let along = normalize(cross(p, across));
                let tension = self.limb_tension * pull;
                add(across, along, tension, self.poisson * tension, 0.0);
            }
        }
        // The relief's potential energy: a plateau's interior is stressed alike every way;
        // at its margin it spreads down its slope, so the stress across the slope stands
        // and the stress along it falls toward Poisson's share of it, as a plateau's margin
        // in plane strain, over the lid's flexural wavelength (a stated choice).
        let force = self.potential_force(p, t);
        let step = 0.5 * self.flexural / self.body_radius;
        let at = |dir: [f64; 3], sign: f64| {
            let q = normalize([
                p[0] + sign * step * dir[0],
                p[1] + sign * step * dir[1],
                p[2] + sign * step * dir[2],
            ]);
            self.potential_force(q, t)
        };
        let ge = (at(east, 1.0) - at(east, -1.0)) / self.flexural;
        let gn = (at(north, 1.0) - at(north, -1.0)) / self.flexural;
        let grad = (ge * ge + gn * gn).sqrt();
        let iso = force / self.elastic;
        let anisotropy = if force != 0.0 {
            (grad * self.flexural / force.abs()).min(1.0)
        } else {
            0.0
        };
        if anisotropy > 0.0 && grad > 0.0 {
            let down = normalize([
                east[0] * ge / grad + north[0] * gn / grad,
                east[1] * ge / grad + north[1] * gn / grad,
                east[2] * ge / grad + north[2] * gn / grad,
            ]);
            let along = normalize(cross(p, down));
            add(
                down,
                along,
                iso,
                iso * (1.0 - (1.0 - self.poisson) * anisotropy),
                0.0,
            );
        } else {
            s.0 += iso;
            s.1 += iso;
        }
        // The loads' bending, as each left the lid (F157).
        for load in &self.loads {
            let (a, b, c) = load.at(p, east, north);
            s.0 += a;
            s.1 += b;
            s.2 += c;
        }
        s
    }

    pub fn body_radius(&self) -> f64 {
        self.body_radius
    }

    /// The stress a fault reads at a point and time: the most and least tensile principal
    /// stresses, Pa, the strike across the most tensile direction, and the frame of the
    /// direction outward from the nearest rise and round it.
    pub fn at(&self, p: [f64; 3], t: f64) -> FieldAt {
        let (see, snn, sen) = self.tensor(p, t);
        let (east, north) = local_basis(p);
        let mean = 0.5 * (see + snn);
        let half = (0.25 * (see - snn).powi(2) + sen * sen).sqrt();
        // The most tensile direction lies at this angle from east toward north; the strike
        // runs across it.
        let theta = 0.5 * (2.0 * sen).atan2(see - snn);
        let strike = normalize([
            -east[0] * theta.sin() + north[0] * theta.cos(),
            -east[1] * theta.sin() + north[1] * theta.cos(),
            -east[2] * theta.sin() + north[2] * theta.cos(),
        ]);
        let (out, round) = match self
            .rises
            .iter()
            .min_by(|a, b| angle_between(a.center, p).total_cmp(&angle_between(b.center, p)))
        {
            Some(r) => {
                let along = dot(r.center, p);
                let out = [
                    p[0] * along - r.center[0],
                    p[1] * along - r.center[1],
                    p[2] * along - r.center[2],
                ];
                let out = if dot(out, out) < 1e-20 {
                    east
                } else {
                    normalize(out)
                };
                (out, normalize(cross(p, out)))
            }
            None => (east, north),
        };
        (mean + half, mean - half, strike, out, round)
    }

    /// The most tensile stress beyond what breaks the brittle layer at a point and time,
    /// Pa, or nothing.
    pub fn excess(&self, p: [f64; 3], t: f64) -> f64 {
        (self.at(p, t).0 - self.yield_tension).max(0.0)
    }

    /// The compression a thrust reads at a point and time (F154): the most compressive
    /// stress, the strike across it, and how far the field decides that strike, its
    /// differential stress against THRUST_STRIKE_SCATTER of the lid's strength in
    /// compression or of the mean squeeze, whichever is greater (F157): a weakness lowers a
    /// plane's strength by a share of the stress that presses it shut, so a squeeze nearly
    /// alike every way decides little however strong.
    pub fn squeeze(&self, p: [f64; 3], t: f64) -> Squeeze {
        let (s1, s2, strike, _, _) = self.at(p, t);
        // The normal faults' strike runs along the most compressive direction; a thrust's
        // runs across it.
        let across = cross(p, strike);
        let strike = if dot(across, across) < 1e-20 {
            local_basis(p).0
        } else {
            normalize(across)
        };
        let differential = (s1 - s2).max(0.0);
        let mean = (-0.5 * (s1 + s2)).max(self.yield_compression);
        let scatter = THRUST_STRIKE_SCATTER * mean;
        Squeeze {
            stress: -s2,
            least: s1,
            strike,
            decided: differential / (differential + scatter.max(1.0)),
        }
    }
}

/// The compression the body's shrinking puts on a whole rock lid alike every way at a time,
/// Pa (F154): its contraction over its radius, times the lid's stiffness alike every way.
pub fn shrinking(thermal: &crate::thermal::ThermalHistory, body_radius: f64, t: f64) -> f64 {
    ROCK_YOUNG_MODULUS / (1.0 - ROCK_POISSON) * thermal.at(t).contraction / body_radius
}

/// One reading of the ground through the loop (F157): its time, the lid's elastic thickness
/// then, the relief as the field reads it then, and the stress of the load laid since the
/// epoch before, bending the lid of its time.
#[derive(Clone, Debug)]
pub struct FieldEpoch {
    pub time: f64,
    pub elastic: f64,
    pub relief: Relief,
    pub load: Option<LoadStress>,
}

/// What the lid's field is built from through the loop (F154, F157): the rises with their
/// axes and the epochs at which the loop read the ground.
#[derive(Clone, Debug)]
pub struct FieldHistory {
    pub rises: Vec<Rise>,
    /// The upwelling limbs between them (F150, F173).
    pub limbs: Limbs,
    pub body_radius: f64,
    pub gravity: f64,
    pub epochs: Vec<FieldEpoch>,
}

impl FieldHistory {
    pub fn new(graph: &ParamGraph, rises: Vec<Rise>, limbs: Limbs) -> Self {
        Self {
            rises,
            limbs,
            body_radius: graph.f64("body.radius"),
            gravity: graph.f64("body.gravity"),
            epochs: Vec::new(),
        }
    }

    pub fn push(&mut self, epoch: FieldEpoch) {
        self.epochs.push(epoch);
    }

    /// How many epochs the loop had read by a time.
    pub fn read_by(&self, t: f64) -> usize {
        self.epochs.partition_point(|e| e.time <= t)
    }

    /// The field at a time: the ground and the lid of the last epoch read by then, or of the
    /// first, and the stresses of every load laid by then; each rise growing along its axis
    /// (F114).
    pub fn field(&self, t: f64) -> StressField<'_> {
        self.field_with(t, true)
    }

    /// The field at a time without the loads' stresses, for what the loads alone add.
    pub fn field_unloaded(&self, t: f64) -> StressField<'_> {
        self.field_with(t, false)
    }

    fn field_with(&self, t: f64, loaded: bool) -> StressField<'_> {
        let k = self.read_by(t);
        let current = self.epochs.get(k.saturating_sub(1));
        let axes = self
            .rises
            .iter()
            .map(|rise| Axial {
                axis: rise.axis,
                stress: RISE_AXIAL_STRESS,
            })
            .collect();
        StressField::new(
            &self.rises,
            axes,
            self.body_radius,
            self.gravity,
            current.map_or(1e3, |e| e.elastic),
            current.map(|e| &e.relief),
            if loaded {
                self.epochs[..k]
                    .iter()
                    .filter_map(|e| e.load.as_ref())
                    .collect()
            } else {
                Vec::new()
            },
        )
        .with_limbs(self.limbs.clone())
    }

    /// Whether any load has been laid.
    pub fn loaded(&self) -> bool {
        self.epochs.iter().any(|e| e.load.is_some())
    }
}

/// The lid's compression on a coarse grid through the loop (F154): each cell's reading at
/// the time last asked for, read again as the rises grow and as the loop reads the ground
/// anew (F157).
pub struct Compression {
    pub grid: CubeSphere,
    /// Each cell's area, m².
    pub area: Vec<f64>,
    pub squeeze: Vec<Squeeze>,
    pub yield_compression: f64,
    /// The rises' growth and the epochs read when the cells were last read, and whether
    /// they have been.
    grown: Vec<f64>,
    epochs: usize,
    read: bool,
}

impl Compression {
    pub fn new(history: &FieldHistory, face: u32) -> Self {
        let grid = CubeSphere::new(face);
        let r = history.body_radius;
        let area = (0..grid.cell_count())
            .map(|i| grid.cell_area(grid.cell(i)) * r * r)
            .collect();
        let yield_compression = history.field(0.0).yield_compression;
        Self {
            grid,
            area,
            squeeze: vec![Squeeze::default(); grid.cell_count()],
            yield_compression,
            grown: Vec::new(),
            epochs: 0,
            read: false,
        }
    }

    /// Read the cells at a time, unless no rise has grown by FIELD_REREAD_GROWTH and the
    /// loop has read no new epoch since they were last read.
    pub fn at(&mut self, t: f64, history: &FieldHistory) {
        let grown: Vec<f64> = history.rises.iter().map(|r| rise_progress(r, t)).collect();
        let moved = grown
            .iter()
            .zip(&self.grown)
            .any(|(a, b)| (a - b).abs() > FIELD_REREAD_GROWTH);
        let epochs = history.read_by(t);
        if self.read && !moved && epochs == self.epochs {
            return;
        }
        let field = history.field(t);
        let grid = self.grid;
        self.squeeze = (0..grid.cell_count())
            .into_par_iter()
            .map(|i| field.squeeze(grid.center(grid.cell(i)), t))
            .collect();
        self.grown = grown;
        self.epochs = epochs;
        self.read = true;
        if std::env::var_os("EXOGEN_DEBUG_FIELD").is_some() {
            let mut decided: Vec<f64> = self.squeeze.iter().map(|q| q.decided).collect();
            let mut diff: Vec<f64> = self
                .squeeze
                .iter()
                .map(|q| (q.stress + q.least) / 1e6)
                .collect();
            decided.sort_by(f64::total_cmp);
            diff.sort_by(f64::total_cmp);
            let q = |v: &[f64], f: f64| v[((v.len() - 1) as f64 * f) as usize];
            eprintln!(
                "field at {t:.2} Gyr ({epochs} epochs): decided p50 {:.2} p90 {:.2}; differential MPa p50 {:.1} p90 {:.1} p99 {:.1}",
                q(&decided, 0.5),
                q(&decided, 0.9),
                q(&diff, 0.5),
                q(&diff, 0.9),
                q(&diff, 0.99)
            );
        }
    }

    /// The reading at the cell that holds a point.
    pub fn reading(&self, p: [f64; 3]) -> Squeeze {
        self.squeeze[self.grid.index(self.grid.locate(p).0)]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A thrust reads the field from its other side (F154): a low squeezes the lid and a
    /// high stretches it; at the low's margin the squeeze is greatest across the slope, so a
    /// thrust strikes along the contour, and the field decides it there, while in the low's
    /// middle, squeezed alike every way, it decides little.
    #[test]
    fn a_thrust_strikes_along_a_low_s_margin_as_the_field_decides() {
        let (r, g, te) = (6.0e6, 8.9, 30e3);
        let grid = CubeSphere::new(128);
        let ramp = |x: f64| {
            let u = x.clamp(0.0, 1.0);
            u * u * (3.0 - 2.0 * u)
        };
        let cap = |p: [f64; 3]| {
            let lat = p[2].asin().to_degrees();
            3000.0 * (ramp((lat - 55.0) / 10.0) - ramp((-lat - 55.0) / 10.0))
        };
        let height: Vec<f32> = (0..grid.cell_count())
            .map(|i| cap(grid.center(grid.cell(i))) as f32)
            .collect();
        let relief = Relief { grid, height };
        let field = StressField::new(&[], Vec::new(), r, g, te, Some(&relief), Vec::new());
        let at = |lat: f64, lon: f64| {
            let (la, lo) = (lat.to_radians(), lon.to_radians());
            [la.cos() * lo.cos(), la.cos() * lo.sin(), la.sin()]
        };
        assert!(field.yield_compression > 3.0 * field.yield_tension);
        let margin = field.squeeze(at(-60.0, 20.0), 1.0);
        let middle = field.squeeze(at(-85.0, 20.0), 1.0);
        let high = field.squeeze(at(85.0, 20.0), 1.0);
        let (east, _) = local_basis(at(-60.0, 20.0));
        assert!(
            margin.stress > 0.0 && middle.stress > margin.stress,
            "{} {}",
            margin.stress,
            middle.stress
        );
        assert!(high.stress < 0.0, "{}", high.stress);
        assert!(
            dot(margin.strike, east).abs() > 0.9,
            "{}",
            dot(margin.strike, east)
        );
        assert!(
            margin.decided > 0.2 && margin.decided > 5.0 * middle.decided,
            "{} {}",
            margin.decided,
            middle.decided
        );
    }

    /// Prints the compression the field puts on each rocky lid: the most compressive
    /// stress's spread, the differential stress and the lid's strength in compression.
    #[test]
    #[ignore = "diagnostic: the field's compression on the rocky lids"]
    fn print_compression() {
        let names =
            std::env::var("EXOGEN_BODY").unwrap_or_else(|_| "Mercury,Moon,Mars,Venus".into());
        for name in names.split(',') {
            let (g, _) = crate::presets::preset(name).unwrap().to_graph();
            let thermal = crate::thermal::ThermalHistory::from_graph(&g);
            let grid = CubeSphere::new(g.f64("detail.loop_resolution") as u32);
            let r = g.f64("body.radius");
            let gravity = g.f64("body.gravity");
            let age = g.f64("body.age");
            let elastic = g.f64("interior.elastic_thickness");
            let layers = crate::stages::primordial_crust::PrimordialCrust::new(&g, grid).init();
            let wavelength = StressField::flexural_wavelength(elastic, gravity);
            let cell = grid.delta() * r;
            let relief = Relief {
                grid,
                height: crate::stages::relaxation::low_pass(&layers.bedrock, cell, wavelength),
            };
            let (hmin, hmax) = relief
                .height
                .iter()
                .fold((f32::MAX, f32::MIN), |a, &h| (a.0.min(h), a.1.max(h)));
            let plume =
                crate::stages::plume::PlumeProcess::new(&g, grid, &thermal, Some(relief.clone()));
            let rises = plume.as_ref().map(|p| p.rises.clone()).unwrap_or_default();
            let axes: Vec<Axial> = rises
                .iter()
                .map(|rise| Axial {
                    axis: rise.axis,
                    stress: RISE_AXIAL_STRESS,
                })
                .collect();
            let field =
                StressField::new(&rises, axes, r, gravity, elastic, Some(&relief), Vec::new());
            let k = ((1.0 + ROCK_FRICTION.powi(2)).sqrt() + ROCK_FRICTION).powi(2);
            let yc = (k - 1.0) * DENSITY_CRUST_ROCK * gravity * FAULT_LAYER_DEPTH;
            println!(
                "{name}: Te {:.0} km, flexural {:.0} km, relief {:.1}..{:.1} km, rises {}, yield {:.1} MPa tension, {:.1} MPa compression",
                elastic / 1e3,
                wavelength / 1e3,
                hmin / 1e3,
                hmax / 1e3,
                rises.len(),
                field.yield_tension / 1e6,
                yc / 1e6
            );
            let probe = CubeSphere::new(48);
            for t in [0.25 * age, 0.5 * age, 0.75 * age, age] {
                let mut s2s = Vec::new();
                let mut diffs = Vec::new();
                for i in 0..probe.cell_count() {
                    let p = probe.center(probe.cell(i));
                    let (s1, s2, _, _, _) = field.at(p, t);
                    s2s.push(s2 / 1e6);
                    diffs.push((s1 - s2) / 1e6);
                }
                s2s.sort_by(f64::total_cmp);
                diffs.sort_by(f64::total_cmp);
                let q = |v: &[f64], f: f64| v[((v.len() - 1) as f64 * f) as usize];
                let beyond =
                    s2s.iter().filter(|&&x| x < -yc / 1e6).count() as f64 / s2s.len() as f64;
                println!(
                    "  t {t:.2}: s2 MPa p1 {:.1} p10 {:.1} p50 {:.1} p90 {:.1} p99 {:.1}; s1-s2 p50 {:.1} p90 {:.1} p99 {:.1}; s2 beyond -Yc {:.3}",
                    q(&s2s, 0.01),
                    q(&s2s, 0.1),
                    q(&s2s, 0.5),
                    q(&s2s, 0.9),
                    q(&s2s, 0.99),
                    q(&diffs, 0.5),
                    q(&diffs, 0.9),
                    q(&diffs, 0.99),
                    beyond
                );
            }
        }
    }

    /// The relief's potential energy stresses a high in tension, alike every way at its
    /// middle, across its slope at its margin, where the faults would run along the
    /// contour; a low stands in compression; level ground is unstressed.
    #[test]
    fn a_high_stands_in_tension_and_spreads_down_its_margin() {
        let (r, g, te) = (6.0e6, 8.9, 30e3);
        let grid = CubeSphere::new(128);
        // A high plateau north of 60° and a low south of it, each ramping over ten degrees.
        let ramp = |x: f64| {
            let u = x.clamp(0.0, 1.0);
            u * u * (3.0 - 2.0 * u)
        };
        let cap = |p: [f64; 3]| {
            let lat = p[2].asin().to_degrees();
            3000.0 * (ramp((lat - 55.0) / 10.0) - ramp((-lat - 55.0) / 10.0))
        };
        let height: Vec<f32> = (0..grid.cell_count())
            .map(|i| cap(grid.center(grid.cell(i))) as f32)
            .collect();
        let relief = Relief { grid, height };
        let field = StressField::new(&[], Vec::new(), r, g, te, Some(&relief), Vec::new());
        let at = |lat: f64, lon: f64| {
            let (la, lo) = (lat.to_radians(), lon.to_radians());
            [la.cos() * lo.cos(), la.cos() * lo.sin(), la.sin()]
        };
        // The middle of the high: tension, nearly alike every way, tens of MPa.
        let (s1, s2, _, _, _) = field.at(at(85.0, 20.0), 1.0);
        assert!(s1 > 10e6 && s2 > 0.8 * s1, "{s1} {s2}");
        // Its margin: the most tensile stress across the slope, the strike along the
        // contour, east-west.
        let (m1, m2, strike, _, _) = field.at(at(60.0, 20.0), 1.0);
        let (east, _) = local_basis(at(60.0, 20.0));
        assert!(m1 > 0.0 && m2 < 0.8 * m1, "{m1} {m2}");
        assert!(dot(strike, east).abs() > 0.9, "{}", dot(strike, east));
        // The low: compression; the plains: nothing.
        let (l1, _, _, _, _) = field.at(at(-85.0, 20.0), 1.0);
        assert!(l1 < -10e6, "{l1}");
        let (p1, p2, _, _, _) = field.at(at(0.0, 20.0), 1.0);
        assert!(p1.abs() < 1.0 && p2.abs() < 1.0, "{p1} {p2}");
    }
}
