// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! The lid's stress field (F150): one field, in pascals, that every fault on a rock lid
//! reads. It is the sum of the rises' stresses, their domes' bending and the membrane
//! stress of their uplift on the sphere (F105, F114), of the upwelling between neighbouring
//! plumes, which pulls the lid apart across the great circle joining them, and of the
//! potential energy of the relief: a high stands in tension, a low in compression, each
//! alike every way within and across the slope at its margin, where a plateau spreads down
//! its slope (Bilotti and Suppe 1999). A fault forms where the most tensile stress passes
//! what breaks the brittle layer, and runs as long as it does.

use crate::constants::*;
use crate::geom::{angle_between, cross, dot, local_basis, normalize};
use crate::grid::CubeSphere;
use crate::stages::plume::{Axial, FieldAt, Rise, dome, rise_progress, total_stress};

/// The upwelling limbs between neighbouring plumes: pairs of rises, by index, closer than
/// UPWELLING_LIMB_REACH of their radii together, and the width the limb's pull spreads
/// over, m.
#[derive(Clone, Debug, Default)]
pub struct Limbs {
    pub pairs: Vec<(usize, usize)>,
    pub width: f64,
}

impl Limbs {
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
        Self { pairs, width }
    }

    /// How much of a limb's pull a point feels at a time, 0 to 1, and the direction across
    /// the limb there: a Gaussian of the distance from the great circle between the two
    /// rises, over the arc between them and fading beyond their centres over the width,
    /// growing as the lesser of the two plumes rises.
    pub fn pull(&self, rises: &[Rise], p: [f64; 3], t: f64, body_radius: f64) -> (f64, [f64; 3]) {
        let mut best = (0.0, [0.0; 3]);
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

/// The relief the potential-energy stress reads: the ground before the rises, low-passed
/// over the lid's flexural wavelength, on the loop's grid.
#[derive(Clone, Debug)]
pub struct Relief {
    pub grid: CubeSphere,
    pub height: Vec<f32>,
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
    relief: Option<Relief>,
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
    /// The tension that breaks the brittle layer, Pa.
    pub yield_tension: f64,
}

impl<'a> StressField<'a> {
    /// The field over a rock lid of elastic thickness `elastic`, m, with the relief given.
    pub(crate) fn new(
        rises: &'a [Rise],
        axes: Vec<Axial>,
        body_radius: f64,
        gravity: f64,
        elastic: f64,
        relief: Option<Relief>,
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
        Self {
            rises,
            axes,
            crest,
            limbs,
            limb_tension,
            relief,
            body_radius,
            gravity,
            elastic: elastic.max(1e3),
            flexural,
            density,
            mantle,
            poisson,
            yield_tension,
        }
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
}

#[cfg(test)]
mod tests {
    use super::*;

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
        let field = StressField::new(&[], Vec::new(), r, g, te, Some(Relief { grid, height }));
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
