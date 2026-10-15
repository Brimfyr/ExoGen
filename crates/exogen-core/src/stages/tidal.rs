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
    pub fn scaled(self, k: f64) -> Stress {
        Stress {
            ee: k * self.ee,
            nn: k * self.nn,
            en: k * self.en,
        }
    }

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
    /// Where the shell is thinner than its mean and floats on an ocean, it deforms, and
    /// is stressed, as its local Love number has it (F76); none for a uniform shell.
    pub local: Option<LocalShell>,
}

/// A floating shell of varying thickness: its thinning map, its Love numbers and its mean
/// thickness, m (F76).
#[derive(Clone, Debug)]
pub struct LocalShell {
    pub map: crate::stages::ice::ShellMap,
    pub love: ShellLove,
    pub mean: f64,
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
            local: None,
        })
    }

    /// The membrane strains, meridional and azimuthal about the bulge axis, at the cosine
    /// of the angle from the axis, for a bulge `amplitude` times the permanent one.
    fn strains(&self, cos_theta: f64, amplitude: f64) -> (f64, f64) {
        let c2 = cos_theta * cos_theta;
        let h = amplitude * self.bulge;
        let l = self.l2_over_h2 * h;
        let p2 = 0.5 * (3.0 * c2 - 1.0);
        let p2_dd = -3.0 * (2.0 * c2 - 1.0);
        let cot_p2_d = -3.0 * c2;
        (
            (l * p2_dd + h * p2) / self.radius,
            (l * cot_p2_d + h * p2) / self.radius,
        )
    }

    /// The membrane stresses from those strains by plane stress.
    fn components(&self, cos_theta: f64, amplitude: f64) -> (f64, f64) {
        let (e_tt, e_pp) = self.strains(cos_theta, amplitude);
        (
            self.stiffness * (e_tt + self.poisson * e_pp),
            self.stiffness * (e_pp + self.poisson * e_tt),
        )
    }

    /// How much more the shell at a point deforms than the mean shell would (F76).
    pub fn amplification(&self, p: [f64; 3]) -> f64 {
        self.local
            .as_ref()
            .map_or(1.0, |s| s.love.local_amplification(s.mean, s.map.factor(p)))
    }

    /// The strain at a point for a bulge of the given amplitude about an axis, in the
    /// point's east-north frame, as a symmetric tensor.
    pub fn strain_about(&self, axis: [f64; 3], p: [f64; 3], amplitude: f64) -> Stress {
        let c = dot(axis, p).clamp(-1.0, 1.0);
        let (e_tt, e_pp) = self.strains(c, amplitude);
        about_axis(axis, p, e_tt, e_pp)
    }

    /// The stress at a point for a bulge of the given amplitude about an axis, in the
    /// point's east-north frame.
    pub fn stress_about(&self, axis: [f64; 3], p: [f64; 3], amplitude: f64) -> Stress {
        self.uniform_about(axis, p, amplitude)
            .scaled(self.amplification(p))
    }

    /// The same for the mean shell, before its local thickness amplifies it.
    fn uniform_about(&self, axis: [f64; 3], p: [f64; 3], amplitude: f64) -> Stress {
        let c = dot(axis, p).clamp(-1.0, 1.0);
        let (s_tt, s_pp) = self.components(c, amplitude);
        about_axis(axis, p, s_tt, s_pp)
    }

    /// The tidal axis at an orbital phase: the sub-planet direction librates in
    /// longitude by 2e sin M about the mean axis.
    pub fn axis_at(mean_axis: [f64; 3], eccentricity: f64, phase: f64) -> [f64; 3] {
        rotate(mean_axis, [0.0, 0.0, 1.0], 2.0 * eccentricity * phase.sin())
    }

    /// The diurnal stress at a point at an orbital phase (mean anomaly from pericentre):
    /// the bulge of the moment less the permanent one the shell has relaxed to.
    pub fn diurnal(&self, mean_axis: [f64; 3], p: [f64; 3], phase: f64) -> Stress {
        self.uniform_diurnal(mean_axis, p, phase)
            .scaled(self.amplification(p))
    }

    fn uniform_diurnal(&self, mean_axis: [f64; 3], p: [f64; 3], phase: f64) -> Stress {
        let amplitude = 1.0 + 3.0 * self.eccentricity * phase.cos();
        let axis = Self::axis_at(mean_axis, self.eccentricity, phase);
        self.uniform_about(axis, p, amplitude) - self.uniform_about(mean_axis, p, 1.0)
    }

    /// The stress of non-synchronous rotation: the bulge has marched on by the lag
    /// angle since the shell last relaxed to it.
    pub fn nsr(&self, mean_axis: [f64; 3], p: [f64; 3]) -> Stress {
        self.uniform_nsr(mean_axis, p).scaled(self.amplification(p))
    }

    fn uniform_nsr(&self, mean_axis: [f64; 3], p: [f64; 3]) -> Stress {
        if self.nsr_shift == 0.0 {
            return Stress::default();
        }
        let shifted = rotate(mean_axis, [0.0, 0.0, 1.0], self.nsr_shift);
        self.uniform_about(shifted, p, 1.0) - self.uniform_about(mean_axis, p, 1.0)
    }

    /// The whole field at a point and phase.
    pub fn total(&self, mean_axis: [f64; 3], p: [f64; 3], phase: f64) -> Stress {
        (self.uniform_diurnal(mean_axis, p, phase) + self.uniform_nsr(mean_axis, p))
            .scaled(self.amplification(p))
    }

    /// The most tensile stress a point sees through an orbit, and the phase of it.
    pub fn peak_tension(&self, mean_axis: [f64; 3], p: [f64; 3]) -> (f64, f64) {
        let nsr = self.uniform_nsr(mean_axis, p);
        let mut best = (f64::NEG_INFINITY, 0.0);
        for k in 0..TIDAL_PHASE_SAMPLES {
            let phase = 2.0 * std::f64::consts::PI * k as f64 / TIDAL_PHASE_SAMPLES as f64;
            let (s1, _, _) = (self.uniform_diurnal(mean_axis, p, phase) + nsr).principal();
            if s1 > best.0 {
                best = (s1, phase);
            }
        }
        (best.0 * self.amplification(p), best.1)
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

/// A symmetric tangent tensor with components `tt` along the great circle from an axis
/// and `pp` across it, in a point's east-north frame.
fn about_axis(axis: [f64; 3], p: [f64; 3], tt: f64, pp: f64) -> Stress {
    let (east, north) = local_basis(p);
    let c = dot(axis, p).clamp(-1.0, 1.0);
    // The meridional direction points away from the axis along the great circle.
    let away = [axis[0] - c * p[0], axis[1] - c * p[1], axis[2] - c * p[2]];
    let len = dot(away, away).sqrt();
    if len < 1e-9 {
        // On the axis the field is isotropic.
        return Stress {
            ee: tt,
            nn: tt,
            en: 0.0,
        };
    }
    let t = normalize([-away[0] / len, -away[1] / len, -away[2] / len]);
    let f = normalize(cross(p, t));
    let (te, tn) = (dot(t, east), dot(t, north));
    let (fe, fn_) = (dot(f, east), dot(f, north));
    Stress {
        ee: tt * te * te + pp * fe * fe,
        nn: tt * tn * tn + pp * fn_ * fn_,
        en: tt * te * tn + pp * fe * fn_,
    }
}

/// The shear modulus of ice, Pa.
pub fn ice_shear_modulus() -> f64 {
    ICE_YOUNG_MODULUS / (2.0 * (1.0 + ICE_POISSON))
}

/// The Love numbers of an ice shell (F75; Beuthe 2018, Icarus 302, Eqs. 53, 58, 65 and
/// 66). Over an ocean, with the core taken as rigid (a good approximation for periodic
/// tides), the shell is a membrane spring whose degree-2 stiffness against the ocean's
/// buoyancy is Λ = 8 (1 + ν) / (5 + ν) μ̂ d / R, μ̂ = μ / (ρ g R) its shear modulus over
/// the ocean's weight, so that h₂ = 1 / (1 − ξ₂ + Λ) and k₂ = (1 + Λ) h₂ − 1, with
/// ξ₂ = 3 ρ / (5 ρ̄) the ocean's self-gravity against the body's. A large moon's shell is
/// soft and follows the ocean (Europa: Λ about 0.06, h₂ 1.16); a small moon's is hard,
/// and its own stiffness, not gravity, sets its tide (Enceladus: Λ about 23 for 23 km, h₂
/// 0.04, against Beuthe's 22 and 0.043). Ice frozen to the rock under it deforms as a
/// homogeneous elastic sphere: h₂ = (5/2) / (1 + 19 μ̂ / 2), μ̂ over the bulk density
/// (Love 1909).
#[derive(Clone, Copy, Debug)]
pub struct ShellLove {
    /// 1 − ξ₂.
    pub base: f64,
    /// Λ per metre of shell.
    pub stiffness_per_m: f64,
    /// h₂ of ice frozen to the rock.
    pub grounded: f64,
}

impl ShellLove {
    pub fn new(radius: f64, gravity: f64, bulk_density: f64) -> Self {
        let (r, g, rho) = (radius.max(1.0), gravity.max(1e-9), bulk_density.max(1.0));
        let mu = ice_shear_modulus();
        let membrane = 8.0 * (1.0 + ICE_POISSON) / (5.0 + ICE_POISSON);
        Self {
            base: (1.0 - 3.0 * DENSITY_OCEAN / (5.0 * rho)).max(0.05),
            stiffness_per_m: membrane * mu / (DENSITY_OCEAN * g * r) / r,
            grounded: 2.5 / (1.0 + 9.5 * mu / (rho * g * r)),
        }
    }

    /// Λ for a shell `d` thick.
    pub fn stiffness(&self, d: f64) -> f64 {
        self.stiffness_per_m * d.max(0.0)
    }

    pub fn h2(&self, d: f64) -> f64 {
        1.0 / (self.base + self.stiffness(d))
    }

    pub fn k2(&self, d: f64) -> f64 {
        (1.0 + self.stiffness(d)) * self.h2(d) - 1.0
    }

    /// How much more a shell thinned to `factor` of a mean `d` deforms than the mean
    /// shell, as its h₂ against the mean's (F76): a hard shell's strain, and so its
    /// stress, goes nearly as the inverse of its local thickness (Beuthe 2018: "surface
    /// stresses are approximately inversely proportional to the local shell thickness"),
    /// a soft shell's barely changes.
    pub fn local_amplification(&self, d: f64, factor: f64) -> f64 {
        self.h2(d * factor.max(1e-3)) / self.h2(d)
    }
}

/// The diurnal tide's strain factor κ (F75): the orbit- and surface-averaged square of the
/// deviatoric strain of the breathing bulge (F41), with the radial strain plane stress
/// gives the membrane, per unit of e H / R; so ⟨ė²⟩ = κ (e ω H / R)². With l₂ / h₂ = 0.25
/// and ν = 0.33 it is about 1.92.
pub fn diurnal_strain_factor() -> f64 {
    static FACTOR: std::sync::OnceLock<f64> = std::sync::OnceLock::new();
    *FACTOR.get_or_init(|| {
        let e = 1e-4;
        let shell = TidalShell {
            radius: 1.0,
            bulge: 1.0,
            l2_over_h2: L2_OVER_H2,
            eccentricity: e,
            stiffness: 1.0,
            poisson: ICE_POISSON,
            nsr_shift: 0.0,
            local: None,
        };
        let axis = [1.0, 0.0, 0.0];
        let (rows, cols, phases) = (48, 96, 32);
        let radial = -ICE_POISSON / (1.0 - ICE_POISSON);
        let (mut sum, mut weight) = (0.0, 0.0);
        for i in 0..rows {
            let lat = std::f64::consts::PI * ((i as f64 + 0.5) / rows as f64 - 0.5);
            let w = lat.cos();
            for j in 0..cols {
                let lon = std::f64::consts::TAU * (j as f64 + 0.5) / cols as f64;
                let p = [lat.cos() * lon.cos(), lat.cos() * lon.sin(), lat.sin()];
                let still = shell.strain_about(axis, p, 1.0);
                for k in 0..phases {
                    let m = std::f64::consts::TAU * k as f64 / phases as f64;
                    let now = shell.strain_about(
                        TidalShell::axis_at(axis, e, m),
                        p,
                        1.0 + 3.0 * e * m.cos(),
                    );
                    let s = now - still;
                    let rr = radial * (s.ee + s.nn);
                    let mean = (s.ee + s.nn + rr) / 3.0;
                    sum += w
                        * ((s.ee - mean).powi(2)
                            + (s.nn - mean).powi(2)
                            + (rr - mean).powi(2)
                            + 2.0 * s.en * s.en);
                }
                weight += w * phases as f64;
            }
        }
        sum / weight / (e * e)
    })
}

/// The tides an ice shell over an ocean dissipates in itself (F75; Ojakangas and
/// Stevenson 1989; the formula of Nimmo et al. 2007, as Lemasquerier et al. 2023 give it).
/// A Maxwell solid strained at the tidal frequency ω turns 2 μ ⟨ė²⟩ / ω · ωτ / (1 + (ωτ)²)
/// of work per unit volume into heat, τ = η / μ its Maxwell time and ⟨ė²⟩ the
/// orbit-averaged square of its deviatoric strain rate. The cold top of a shell is elastic
/// and dissipates nothing; the warm base, where the viscosity falls through μ / ω, takes
/// nearly all of it ("nearly all of the tidal dissipation is concentrated in the lowermost
/// few kilometers of the shell"). Across a conductive shell the viscosity falls by a factor
/// e in each R T_m² / (E ΔT) of the thickness d, so the shell dissipates, per unit area,
/// H = 2 μ ⟨ė²⟩ / ω · arctan(μ / (η_m ω)) · R T_m² d / (E ΔT), with the viscosity at the
/// melting point and the activation energy of F39. The strain rate is the diurnal tide's
/// (F41): ⟨ė²⟩ = κ (e ω h₂ (M / m)(R / a)³)², Europa's 2.4e-20 s⁻² against Ojakangas and
/// Stevenson's 3.3e-20. A soft shell dissipates in proportion to its thickness, so its
/// heating and its conduction meet at one thickness; a hard shell's h₂, and so its
/// heating, falls as it thickens.
#[derive(Clone, Copy, Debug)]
pub struct ShellTides {
    pub love: ShellLove,
    /// (M / m)(R / a)³: the static tide's height over the radius per unit h₂.
    pub tide: f64,
    /// The tidal frequency, the mean motion of a locked moon, rad/s.
    pub frequency: f64,
    /// k ΔT, W/m: the heat a shell a metre thick conducts, per square metre.
    pub conduction: f64,
    /// 2 μ / ω · arctan(μ / (η_m ω)) · R T_m² / (E ΔT), Pa s per metre of shell: the
    /// heating is this times ⟨ė²⟩ times the thickness.
    pub dissipation: f64,
    /// The H2O layer, m: a floating shell is thinner.
    pub layer: f64,
    pub radius: f64,
}

/// A shell in equilibrium: its thickness, m, the tides it dissipates, W/m², and whether
/// it floats on an ocean.
#[derive(Clone, Copy, Debug)]
pub struct ShellState {
    pub thickness: f64,
    pub heating: f64,
    pub floating: bool,
}

impl ShellTides {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        love: ShellLove,
        primary_mass: f64,
        mass: f64,
        radius: f64,
        semi_major_axis: f64,
        mean_motion: f64,
        surface_temperature: f64,
        melting_point: f64,
        layer: f64,
    ) -> Self {
        let dt = (melting_point - surface_temperature).max(1.0);
        let mu = ice_shear_modulus();
        let w = mean_motion.max(1e-30);
        Self {
            love,
            tide: primary_mass / mass.max(1.0) * (radius / semi_major_axis.max(1.0)).powi(3),
            frequency: mean_motion,
            conduction: ICE_K_NUM / (0.5 * (melting_point + surface_temperature)).max(20.0) * dt,
            dissipation: 2.0 * mu / w
                * (mu / (ICE_VISCOSITY_MELT * w)).atan()
                * GAS_CONSTANT
                * melting_point
                * melting_point
                / (ICE_CREEP_ACTIVATION * dt),
            layer,
            radius,
        }
    }

    /// The shell of an icy moon locked to its planet, from the graph, mirroring
    /// `derive::shell_tides`; none for any other body.
    pub fn from_graph(graph: &ParamGraph) -> Option<Self> {
        let icy = matches!(
            graph.str("body.composition_class"),
            "rock_ice" | "ice_dominated"
        );
        if !icy || graph.str("orbit.primary") != "planet" || graph.str("rotation.mode") != "locked"
        {
            return None;
        }
        let ammonia = graph
            .get("volatiles.inventory")
            .as_map()
            .and_then(|m| m.get("NH3").copied())
            .unwrap_or(0.0)
            > 0.0;
        let radius = graph.f64("body.radius");
        Some(Self::new(
            ShellLove::new(radius, graph.f64("body.gravity"), graph.f64("body.density")),
            graph.f64("orbit.primary_mass"),
            graph.f64("body.mass"),
            radius,
            graph.f64("orbit.semi_major_axis"),
            graph.f64("orbit.mean_motion"),
            graph.f64("atmosphere.surface_temperature_mean"),
            if ammonia {
                T_EUTECTIC_AMMONIA_WATER
            } else {
                T_MELT_WATER
            },
            graph.f64("body.ice_layer_thickness"),
        ))
    }

    /// ⟨ė²⟩, s⁻², of a shell with this h₂ at this eccentricity.
    pub fn strain_rate_squared(&self, h2: f64, eccentricity: f64) -> f64 {
        diurnal_strain_factor() * (eccentricity * self.frequency * h2 * self.tide).powi(2)
    }

    /// The tides a shell `d` thick dissipates, W/m².
    pub fn heating(&self, d: f64, eccentricity: f64, floating: bool) -> f64 {
        let h2 = if floating {
            self.love.h2(d)
        } else {
            self.love.grounded
        };
        self.dissipation * self.strain_rate_squared(h2, eccentricity) * d.max(0.0)
    }

    /// The shell whose conduction carries the heat from below, `basal` W/m², and its own:
    /// the thickest such shell thinner than the water, or ice frozen to the rock where
    /// even the whole layer conducts more than both.
    pub fn equilibrium(&self, basal: f64, eccentricity: f64) -> ShellState {
        let k = self.conduction;
        let top = self.layer.min(self.radius);
        let gap = |d: f64| k / d - basal - self.heating(d, eccentricity, true);
        if top <= 1.0 || gap(top) > 0.0 {
            let heating = self.heating(top.max(1.0), eccentricity, false);
            let q = basal + heating;
            let thickness = if q > 0.0 {
                (k / q).min(self.radius)
            } else {
                self.radius
            };
            return ShellState {
                thickness,
                heating,
                floating: false,
            };
        }
        // Down from the water in steps until the conduction outruns the heat, then
        // between the last two.
        let mut hi = top;
        let mut lo = top / 1.25;
        while lo > 1.0 && gap(lo) <= 0.0 {
            hi = lo;
            lo /= 1.25;
        }
        for _ in 0..50 {
            let mid = (lo * hi).sqrt();
            if gap(mid) > 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let thickness = (lo * hi).sqrt();
        ShellState {
            thickness,
            heating: self.heating(thickness, eccentricity, true),
            floating: true,
        }
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
            local: None,
        }
    }

    /// A floating shell's Love numbers follow Beuthe (2018): a large moon's shell is soft
    /// and follows its ocean, a small one's is hard (F75).
    #[test]
    fn shell_love_numbers_follow_the_membrane_spring() {
        // Enceladus, a 23 km shell: Λ 22.2 and h₂ near 0.043 (Beuthe's Table 3 and Eq. 65).
        let enceladus = ShellLove::new(252.1e3, 0.1134, 1609.0);
        let lam = enceladus.stiffness(23e3);
        assert!((20.0..24.5).contains(&lam), "Λ {lam}");
        let h2 = enceladus.h2(23e3);
        assert!((0.038..0.048).contains(&h2), "h2 {h2}");
        // Europa, a 25 km shell: soft, h₂ near the fluid-crust 1.25 and k₂ near the 0.21
        // to 0.23 Bagheri et al. (2026) predict.
        let europa = ShellLove::new(1560.8e3, 1.315, 3013.0);
        assert!(europa.stiffness(25e3) < 0.1);
        let (h2, k2) = (europa.h2(25e3), europa.k2(25e3));
        assert!(
            (1.1..1.25).contains(&h2) && (0.2..0.26).contains(&k2),
            "{h2} {k2}"
        );
        // Thinned to a quarter, the hard shell deforms some three and a half times as
        // much; the soft one barely more.
        let hard = enceladus.local_amplification(23e3, 0.25);
        assert!((3.0..4.2).contains(&hard), "{hard}");
        assert!(europa.local_amplification(25e3, 0.25) < 1.06);
    }

    /// The diurnal strain of the breathing bulge, averaged over the sphere and the orbit
    /// (F75), and Europa's strain rate against Ojakangas and Stevenson's (1989) 3.3e-20
    /// s⁻² (Lemasquerier et al. 2023).
    #[test]
    fn the_diurnal_strain_rate_matches_europa_s() {
        let kappa = diurnal_strain_factor();
        assert!((1.85..2.0).contains(&kappa), "κ {kappa}");
        let love = ShellLove::new(1560.8e3, 1.315, 3013.0);
        let tides = ShellTides::new(
            love, 1.89813e27, 4.7998e22, 1560.8e3, 6.711e8, 2.0473e-5, 96.0, 273.15, 120e3,
        );
        let rate = tides.strain_rate_squared(love.h2(30e3), 0.0094);
        assert!((1.5e-20..4.5e-20).contains(&rate), "{rate}");
        // Its shell meets its heat at the 29 ± 10 km Juno measured (Levin et al. 2025),
        // dissipating some half a terawatt (Lemasquerier et al. 2023: 0.56 TW for 25 km).
        let s = tides.equilibrium(0.0072, 0.0094);
        assert!(
            s.floating && (19e3..39e3).contains(&s.thickness),
            "{}",
            s.thickness
        );
        let power = s.heating * 4.0 * std::f64::consts::PI * 1560.8e3f64.powi(2);
        assert!((2e11..9e11).contains(&power), "{power}");
        // A fivefold eccentricity thins it some fivefold (Kihoulou et al. 2025: 35 km to
        // 5 km).
        let pulse = tides.equilibrium(0.0072, 0.047);
        assert!((3e3..9e3).contains(&pulse.thickness), "{}", pulse.thickness);
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
