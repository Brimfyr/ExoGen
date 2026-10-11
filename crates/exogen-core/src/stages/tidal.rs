// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! The tidal stress field of a thin ice shell over an ocean (PLAN.md 4.3, Greenberg et
//! al. 1998, Hoppa et al. 1999).
//!
//! A moon on an eccentric orbit, locked to its planet, carries a permanent tidal bulge of
//! height H = h₂ (M_p / M)(R / a)³ R towards the planet, and the bulge breathes through
//! each orbit: by 3e in height as the distance changes and by 2e in longitude as the
//! moon's constant spin runs ahead of and behind its uneven orbital motion. The shell has
//! relaxed to the permanent bulge; the breathing stresses it, by some hundred kilopascals
//! on Europa. Where the shell also rotates slowly against the tidal axis (non-synchronous
//! rotation), the bulge's slow march round the body adds a stress of the same form for
//! the angle the shell lags before it relaxes. Cracks open perpendicular to the most
//! tensile stress of the moment: along the slowly turning field for the long lineae, and
//! along the field as it turns through an orbit for the cycloids.
//!
//! The stresses are the membrane stresses of a thin elastic shell under a degree-2
//! deformation: with u_r = H P₂(cos θ) and u_θ = L dP₂/dθ, L = (l₂ / h₂) H, the strains
//! are ε_θθ = (L P₂'' + H P₂) / R and ε_φφ = (L cot θ P₂' + H P₂) / R about the bulge axis,
//! and plane stress gives σ_θθ and σ_φφ (Vening Meinesz 1947, Melosh 1977).

use crate::constants::*;
use crate::geom::{cross, dot, local_basis, normalize, rotate};
use crate::params::ParamGraph;

/// A symmetric stress tensor in the tangent plane, in the local east-north frame,
/// tension positive.
#[derive(Clone, Copy, Debug, Default)]
pub struct Stress {
    pub ee: f64,
    pub nn: f64,
    pub en: f64,
}

impl std::ops::Add for Stress {
    type Output = Stress;
    fn add(self, o: Stress) -> Stress {
        Stress {
            ee: self.ee + o.ee,
            nn: self.nn + o.nn,
            en: self.en + o.en,
        }
    }
}

impl std::ops::Sub for Stress {
    type Output = Stress;
    fn sub(self, o: Stress) -> Stress {
        Stress {
            ee: self.ee - o.ee,
            nn: self.nn - o.nn,
            en: self.en - o.en,
        }
    }
}

impl Stress {
    /// The principal stresses, most tensile first, and the azimuth of the most tensile
    /// one from east towards north, radians.
    pub fn principal(&self) -> (f64, f64, f64) {
        let mean = 0.5 * (self.ee + self.nn);
        let r = (0.25 * (self.ee - self.nn).powi(2) + self.en * self.en).sqrt();
        let psi = 0.5 * (2.0 * self.en).atan2(self.ee - self.nn);
        (mean + r, mean - r, psi)
    }
}

/// The shell and orbit that set the field.
#[derive(Clone, Debug)]
pub struct TidalShell {
    pub radius: f64,
    /// Height of the permanent bulge at the sub-planet point, m.
    pub bulge: f64,
    pub l2_over_h2: f64,
    pub eccentricity: f64,
    /// E / (1 − ν²) and ν of the shell.
    pub stiffness: f64,
    pub poisson: f64,
    /// The angle the shell lags the tidal axis before it relaxes, rad; zero without
    /// non-synchronous rotation.
    pub nsr_shift: f64,
}

impl TidalShell {
    /// The shell of a locked moon from the graph, or none for a body without a planet
    /// to raise a tide or one that is not locked to it.
    pub fn new(graph: &ParamGraph) -> Option<Self> {
        Self::build(
            graph.str("orbit.primary"),
            graph.str("rotation.mode"),
            graph.f64("body.radius"),
            graph.f64("body.mass"),
            graph.f64("parent.mass"),
            graph.f64("orbit.semi_major_axis"),
            graph.f64("orbit.eccentricity"),
            graph.f64("interior.love_number_h2"),
        )
    }

    /// The same from a derivation context.
    pub fn from_ctx(c: &crate::params::Ctx) -> Option<Self> {
        Self::build(
            c.str("orbit.primary"),
            c.str("rotation.mode"),
            c.f64("body.radius"),
            c.f64("body.mass"),
            c.f64("parent.mass"),
            c.f64("orbit.semi_major_axis"),
            c.f64("orbit.eccentricity"),
            c.f64("interior.love_number_h2"),
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn build(
        primary: &str,
        rotation: &str,
        radius: f64,
        mass: f64,
        parent: f64,
        a: f64,
        eccentricity: f64,
        h2: f64,
    ) -> Option<Self> {
        if primary != "planet" || rotation != "locked" {
            return None;
        }
        if mass <= 0.0 || parent <= 0.0 || a <= 0.0 || radius <= 0.0 {
            return None;
        }
        let bulge = h2 * (parent / mass) * (radius / a).powi(3) * radius;
        Some(Self {
            radius,
            bulge,
            l2_over_h2: L2_OVER_H2,
            eccentricity,
            stiffness: ICE_YOUNG_MODULUS / (1.0 - ICE_POISSON * ICE_POISSON),
            poisson: ICE_POISSON,
            nsr_shift: 0.0,
        })
    }

    /// The membrane stresses, meridional and azimuthal about the bulge axis, at the
    /// cosine of the angle from the axis, for a bulge `amplitude` times the permanent one.
    fn components(&self, cos_theta: f64, amplitude: f64) -> (f64, f64) {
        let c2 = cos_theta * cos_theta;
        let h = amplitude * self.bulge;
        let l = self.l2_over_h2 * h;
        let p2 = 0.5 * (3.0 * c2 - 1.0);
        let p2_dd = -3.0 * (2.0 * c2 - 1.0);
        let cot_p2_d = -3.0 * c2;
        let e_tt = (l * p2_dd + h * p2) / self.radius;
        let e_pp = (l * cot_p2_d + h * p2) / self.radius;
        (
            self.stiffness * (e_tt + self.poisson * e_pp),
            self.stiffness * (e_pp + self.poisson * e_tt),
        )
    }

    /// The stress at a point for a bulge of the given amplitude about an axis, in the
    /// point's east-north frame.
    pub fn stress_about(&self, axis: [f64; 3], p: [f64; 3], amplitude: f64) -> Stress {
        let (east, north) = local_basis(p);
        let c = dot(axis, p).clamp(-1.0, 1.0);
        let (s_tt, s_pp) = self.components(c, amplitude);
        // The meridional direction points away from the axis along the great circle.
        let away = [axis[0] - c * p[0], axis[1] - c * p[1], axis[2] - c * p[2]];
        let len = dot(away, away).sqrt();
        if len < 1e-9 {
            // On the axis the field is isotropic.
            return Stress {
                ee: s_tt,
                nn: s_tt,
                en: 0.0,
            };
        }
        let t = normalize([-away[0] / len, -away[1] / len, -away[2] / len]);
        let f = normalize(cross(p, t));
        let (te, tn) = (dot(t, east), dot(t, north));
        let (fe, fn_) = (dot(f, east), dot(f, north));
        Stress {
            ee: s_tt * te * te + s_pp * fe * fe,
            nn: s_tt * tn * tn + s_pp * fn_ * fn_,
            en: s_tt * te * tn + s_pp * fe * fn_,
        }
    }

    /// The tidal axis at an orbital phase: the sub-planet direction librates in
    /// longitude by 2e sin M about the mean axis.
    pub fn axis_at(mean_axis: [f64; 3], eccentricity: f64, phase: f64) -> [f64; 3] {
        rotate(mean_axis, [0.0, 0.0, 1.0], 2.0 * eccentricity * phase.sin())
    }

    /// The diurnal stress at a point at an orbital phase (mean anomaly from pericentre):
    /// the bulge of the moment less the permanent one the shell has relaxed to.
    pub fn diurnal(&self, mean_axis: [f64; 3], p: [f64; 3], phase: f64) -> Stress {
        let amplitude = 1.0 + 3.0 * self.eccentricity * phase.cos();
        let axis = Self::axis_at(mean_axis, self.eccentricity, phase);
        self.stress_about(axis, p, amplitude) - self.stress_about(mean_axis, p, 1.0)
    }

    /// The stress of non-synchronous rotation: the bulge has marched on by the lag
    /// angle since the shell last relaxed to it.
    pub fn nsr(&self, mean_axis: [f64; 3], p: [f64; 3]) -> Stress {
        if self.nsr_shift == 0.0 {
            return Stress::default();
        }
        let shifted = rotate(mean_axis, [0.0, 0.0, 1.0], self.nsr_shift);
        self.stress_about(shifted, p, 1.0) - self.stress_about(mean_axis, p, 1.0)
    }

    /// The whole field at a point and phase.
    pub fn total(&self, mean_axis: [f64; 3], p: [f64; 3], phase: f64) -> Stress {
        self.diurnal(mean_axis, p, phase) + self.nsr(mean_axis, p)
    }

    /// The most tensile stress a point sees through an orbit, and the phase of it.
    pub fn peak_tension(&self, mean_axis: [f64; 3], p: [f64; 3]) -> (f64, f64) {
        let mut best = (f64::NEG_INFINITY, 0.0);
        for k in 0..TIDAL_PHASE_SAMPLES {
            let phase = 2.0 * std::f64::consts::PI * k as f64 / TIDAL_PHASE_SAMPLES as f64;
            let (s1, _, _) = self.total(mean_axis, p, phase).principal();
            if s1 > best.0 {
                best = (s1, phase);
            }
        }
        best
    }

    /// The largest peak tension anywhere on the body, Pa: the amplitude of the field.
    pub fn amplitude(&self) -> f64 {
        let axis = [1.0, 0.0, 0.0];
        let mut best: f64 = 0.0;
        for i in 0..24 {
            let lat = -1.5 + 3.0 * (i as f64 + 0.5) / 24.0;
            for j in 0..48 {
                let lon = std::f64::consts::PI * (j as f64 + 0.5) / 24.0;
                let p = [lat.cos() * lon.cos(), lat.cos() * lon.sin(), lat.sin()];
                best = best.max(self.peak_tension(axis, p).0);
            }
        }
        best
    }
}

/// The Love number h₂ of a body: 1 + k₂ for a fluid-like interior, reduced by the
/// rigidity of an ice shell in proportion to its thickness over the radius.
pub fn love_number_h2(k2: f64, shell_thickness: f64, radius: f64, icy: bool) -> f64 {
    let base = 1.0 + k2;
    if icy && radius > 0.0 {
        base / (1.0 + SHELL_RIGIDITY_FACTOR * shell_thickness / radius)
    } else {
        base
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn europa() -> TidalShell {
        TidalShell {
            radius: 1.5608e6,
            bulge: 1.1 * (1.898e27 / 4.8e22) * (1.5608e6f64 / 6.711e8).powi(3) * 1.5608e6,
            l2_over_h2: L2_OVER_H2,
            eccentricity: 0.0094,
            stiffness: ICE_YOUNG_MODULUS / (1.0 - ICE_POISSON * ICE_POISSON),
            poisson: ICE_POISSON,
            nsr_shift: 1f64.to_radians(),
        }
    }

    #[test]
    fn europa_breathes_by_a_hundred_kilopascals() {
        let e = europa();
        assert!(e.bulge > 500.0 && e.bulge < 2_000.0, "bulge {}", e.bulge);
        let axis = [1.0, 0.0, 0.0];
        // On the axis the static field is isotropic; off it the principal axes lie
        // along and across the meridian from the axis.
        let on = e.stress_about(axis, axis, 1.0);
        assert!((on.ee - on.nn).abs() < 1e-6 && on.en.abs() < 1e-6);
        let p = normalize([0.5, 0.5, 0.0]);
        let s = e.stress_about(axis, p, 1.0);
        assert!(s.en.abs() < 1e-6 * s.ee.abs().max(s.nn.abs()));
        // The diurnal tension peaks between 30 and 400 kPa and changes sign through
        // the orbit.
        let amp = e.amplitude();
        assert!(amp > 3e4 && amp < 4e5, "amplitude {amp}");
        let (peri, _, _) = e.diurnal(axis, axis, 0.0).principal();
        let (apo, _, _) = e.diurnal(axis, axis, std::f64::consts::PI).principal();
        assert!(peri * apo < 0.0, "{peri} {apo}");
        // Non-synchronous rotation by a degree stresses the shell as much again.
        let (n1, _, _) = e.nsr(axis, p).principal();
        assert!(n1.abs() > 1e4, "{n1}");
    }
}
