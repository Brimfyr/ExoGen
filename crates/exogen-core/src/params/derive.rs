// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! Derivation functions for the parameter graph.
//!
//! `physics` holds pure functions of plain numbers, tested on their own. The node
//! functions below read upstream values through the context and call them. Formula
//! numbers refer to docs/parameters.md.

use std::collections::BTreeMap;
use std::f64::consts::PI;

use super::value::Curve;
use super::{Ctx, Value};
use crate::constants::*;

pub const ATMOSPHERE_SPECIES: &[&str] = &[
    "N2", "O2", "CO2", "CH4", "H2", "He", "Ar", "SO2", "NH3", "H2O",
];
pub const VOLATILE_SPECIES: &[&str] = &["H2O", "CO2", "CH4", "N2", "NH3", "SO2", "C2H6"];

/// Phase state codes used in the `volatiles.phase_state` map.
pub const PHASE_ABSENT: f64 = 0.0;
pub const PHASE_GAS: f64 = 1.0;
pub const PHASE_LIQUID: f64 = 2.0;
pub const PHASE_SOLID: f64 = 3.0;
pub const PHASE_LIQUID_AND_SOLID: f64 = 4.0;

pub fn phase_name(code: f64) -> &'static str {
    match code as i32 {
        1 => "gas",
        2 => "liquid",
        3 => "solid",
        4 => "liquid and solid",
        _ => "absent",
    }
}

pub mod physics {
    use super::*;

    pub fn molar_mass(species: &str) -> f64 {
        MOLAR_MASS_BY_GAS
            .iter()
            .find(|(n, _)| *n == species)
            .map(|(_, m)| *m)
            .unwrap_or(0.028)
    }

    pub fn heat_capacity(species: &str) -> f64 {
        CP_BY_GAS
            .iter()
            .find(|(n, _)| *n == species)
            .map(|(_, c)| *c)
            .unwrap_or(1_000.0)
    }

    pub fn volatile(species: &str) -> Option<&'static Volatile> {
        VOLATILES.iter().find(|v| v.name == species)
    }

    // ---------------------------------------------------------------- star

    /// F1: main-sequence luminosity at the reference age, solar units (Eker et al. 2018).
    pub fn mass_luminosity(m: f64) -> f64 {
        let last = MASS_LUMINOSITY[MASS_LUMINOSITY.len() - 1];
        let (mut a, mut b) = (last.1, last.2);
        for (upper, aa, bb) in MASS_LUMINOSITY {
            if m <= upper {
                a = aa;
                b = bb;
                break;
            }
        }
        10f64.powf(a * m.max(0.05).log10() + b) * MASS_LUMINOSITY_SOLAR_NORMALISATION
    }

    /// Main-sequence lifetime in Gyr from the fuel-over-burn-rate scaling.
    pub fn main_sequence_lifetime_gyr(m: f64, l_ref: f64) -> f64 {
        MAIN_SEQUENCE_LIFETIME_SUN_GYR * m / l_ref
    }

    /// F3: luminosity at a given age from the reference luminosity (Gough 1981 form).
    pub fn luminosity_at_age(l_ref: f64, age_gyr: f64, t_ms_gyr: f64) -> f64 {
        let t_ref = GOUGH_REFERENCE_FRACTION * t_ms_gyr;
        let age = age_gyr.clamp(0.0, t_ms_gyr);
        l_ref / (1.0 + GOUGH_COEFFICIENT * (1.0 - age / t_ref))
    }

    /// F2: main-sequence radius in solar units (Eker et al. 2018 quadratic, extended above
    /// 1.5 solar masses with a power law continuous at the break).
    pub fn stellar_radius(m: f64) -> f64 {
        let (a, b, c) = MASS_RADIUS_QUADRATIC;
        let quad = |x: f64| a * x * x + b * x + c;
        if m <= MASS_RADIUS_BREAK {
            quad(m.max(0.1))
        } else {
            quad(MASS_RADIUS_BREAK) * (m / MASS_RADIUS_BREAK).powf(MASS_RADIUS_UPPER_EXPONENT)
        }
    }

    pub fn effective_temperature(l_sun: f64, r_sun: f64) -> f64 {
        T_SUN * (l_sun / (r_sun * r_sun)).powf(0.25)
    }

    pub fn spectral_type(teff: f64) -> String {
        const CLASSES: [(&str, f64, f64); 7] = [
            ("O", 30_000.0, 60_000.0),
            ("B", 10_000.0, 30_000.0),
            ("A", 7_500.0, 10_000.0),
            ("F", 6_000.0, 7_500.0),
            ("G", 5_200.0, 6_000.0),
            ("K", 3_700.0, 5_200.0),
            ("M", 2_400.0, 3_700.0),
        ];
        for (c, lo, hi) in CLASSES {
            if teff >= lo {
                let sub = (((hi - teff) / (hi - lo)) * 10.0).floor().clamp(0.0, 9.0) as i32;
                return format!("{c}{sub}V");
            }
        }
        "M9V".to_string()
    }

    /// F4: XUV flux relative to the present Sun at 1 AU, at distance `d_au` from the star.
    pub fn xuv_relative(age_gyr: f64, d_au: f64) -> f64 {
        (age_gyr.max(XUV_SATURATION_AGE_GYR) / AGE_SUN_GYR).powf(XUV_AGE_EXPONENT) / (d_au * d_au)
    }

    // ---------------------------------------------------------------- orbit and body

    pub fn kepler_period(a: f64, m_total: f64) -> f64 {
        2.0 * PI * (a * a * a / (G * m_total)).sqrt()
    }

    pub fn orbital_speed(m_central: f64, d: f64) -> f64 {
        (G * m_central / d).sqrt()
    }

    pub fn gravity(m: f64, r: f64) -> f64 {
        G * m / (r * r)
    }

    pub fn escape_velocity(m: f64, r: f64) -> f64 {
        (2.0 * G * m / r).sqrt()
    }

    /// Bulk-to-uncompressed density factor (F8 helper).
    pub fn compression_factor(mass: f64) -> f64 {
        1.0 + COMPRESSION_A * (mass / M_EARTH).powf(COMPRESSION_B)
    }

    /// F8: (core, rock, ice) mass fractions from the uncompressed density.
    pub fn composition_fractions(rho_u: f64) -> (f64, f64, f64) {
        let core_share = CHONDRITIC_CORE_TO_ROCK / (1.0 + CHONDRITIC_CORE_TO_ROCK);
        let rho_rock_iron = 1.0 / (core_share / DENSITY_IRON + (1.0 - core_share) / DENSITY_ROCK);
        if rho_u >= rho_rock_iron {
            let f = ((1.0 / DENSITY_ROCK - 1.0 / rho_u)
                / (1.0 / DENSITY_ROCK - 1.0 / DENSITY_IRON))
                .clamp(0.0, 1.0);
            (f, 1.0 - f, 0.0)
        } else {
            let fi = ((1.0 / rho_u - 1.0 / rho_rock_iron)
                / (1.0 / DENSITY_ICE - 1.0 / rho_rock_iron))
                .clamp(0.0, 1.0);
            let core = (1.0 - fi) * core_share;
            (core, 1.0 - fi - core, fi)
        }
    }

    pub fn composition_class(core: f64, ice: f64) -> &'static str {
        if core > IRON_RICH_CORE_FRACTION {
            "iron_rich"
        } else if ice < ROCK_MAX_ICE_FRACTION {
            "rock"
        } else if ice < ROCK_ICE_MAX_ICE_FRACTION {
            "rock_ice"
        } else {
            "ice_dominated"
        }
    }

    /// Thickness of the outer ice layer for a given ice mass fraction, metres.
    pub fn ice_layer_thickness(r: f64, mass: f64, ice_fraction: f64) -> f64 {
        if ice_fraction <= 0.0 {
            return 0.0;
        }
        let v_ice = ice_fraction * mass / DENSITY_ICE;
        let inner3 = r * r * r - 3.0 * v_ice / (4.0 * PI);
        if inner3 <= 0.0 { r } else { r - inner3.cbrt() }
    }

    /// F9: Darwin–Radau flattening.
    pub fn oblateness(omega: f64, r: f64, m: f64, c: f64) -> f64 {
        let q = omega * omega * r * r * r / (G * m);
        2.5 * q / (1.0 + 6.25 * (1.0 - 1.5 * c).powi(2))
    }

    /// F10: tidal locking timescale in seconds.
    pub fn locking_timescale(
        omega: f64,
        a: f64,
        moment: f64,
        q: f64,
        m_primary: f64,
        k2: f64,
        r: f64,
    ) -> f64 {
        omega * a.powi(6) * moment * q / (3.0 * G * m_primary * m_primary * k2 * r.powi(5))
    }

    /// F11: rotation regime from the deformation radius against the body radius.
    pub fn rotation_regime(g: f64, scale_height: f64, omega: f64, r: f64) -> &'static str {
        if omega <= 0.0 {
            return "slow";
        }
        let c = (g * scale_height.max(1.0)).sqrt();
        let l_ro = c / (2.0 * omega);
        let ratio = l_ro / r;
        if ratio < 0.5 {
            "rapid"
        } else if ratio > 2.0 {
            "slow"
        } else {
            "rhines"
        }
    }

    // ---------------------------------------------------------------- interior

    /// F12: radiogenic heat per kilogram of silicate for a body `age_gyr` after formation,
    /// assuming solar-system initial abundances.
    pub fn radiogenic_per_kg(age_gyr: f64) -> f64 {
        let dt = SOLAR_SYSTEM_AGE_GYR - age_gyr;
        let total: f64 = RADIOGENIC_ISOTOPES.iter().map(|(s, _)| s).sum();
        let mut f = 0.0;
        for (share, half_life) in RADIOGENIC_ISOTOPES {
            f += share / total * (std::f64::consts::LN_2 / half_life * dt).exp();
        }
        RADIOGENIC_NOW_W_PER_KG * f
    }

    /// F13: tidal power in watts.
    pub fn tidal_power(k2: f64, q: f64, m_primary: f64, r: f64, n: f64, e: f64, a: f64) -> f64 {
        10.5 * (k2 / q) * G * m_primary * m_primary * r.powi(5) * n * e * e / a.powi(6)
    }

    pub fn secular_fraction(r: f64) -> f64 {
        let t = ((r - SECULAR_RADIUS_SMALL) / (SECULAR_RADIUS_LARGE - SECULAR_RADIUS_SMALL))
            .clamp(0.0, 1.0);
        SECULAR_FRACTION_SMALL + (SECULAR_FRACTION_LARGE - SECULAR_FRACTION_SMALL) * t
    }

    pub fn ice_conductivity(t: f64) -> f64 {
        ICE_K_NUM / t.max(20.0)
    }

    /// F14 and F18: conductive thickness for a given temperature drop and flux.
    pub fn conductive_thickness(k: f64, dt: f64, q: f64) -> f64 {
        if q <= 0.0 {
            f64::INFINITY
        } else {
            (k * dt / q).max(0.0)
        }
    }

    /// F15: strength-limited relief, metres.
    pub fn relief_limit(strength: f64, rho: f64, g: f64) -> f64 {
        RELIEF_COEFFICIENT * strength / (rho * g)
    }

    /// F16: cumulative radial contraction, metres.
    pub fn contraction(r: f64, uncompressed_density: f64, age_gyr: f64) -> f64 {
        // Mantle cooling shrinks every body a little; a large core, cooling and
        // solidifying, is what shrinks Mercury by kilometres. The iron content comes
        // from the uncompressed density, which the composition model gets right where
        // its core fraction for small bodies is still rough.
        let iron = ((uncompressed_density - CONTRACTION_DENSITY_ROCK)
            / (CONTRACTION_DENSITY_IRON_RICH - CONTRACTION_DENSITY_ROCK))
            .clamp(0.0, 1.0);
        let hot_start = (r / CONTRACTION_HOT_START_RADIUS).powi(2).min(1.0);
        let dt_total = CONTRACTION_DT_MANTLE * hot_start + CONTRACTION_DT_CORE * iron * iron;
        let since_onset = (age_gyr - CONTRACTION_ONSET_GYR).max(0.0);
        let dt = dt_total * (1.0 - (-since_onset / CONTRACTION_TIMESCALE_GYR).exp());
        THERMAL_EXPANSIVITY / 3.0 * dt * r
    }

    /// F17: magma production rate, m³/s.
    pub fn magma_rate(q: f64, area: f64, efficiency: f64, heat_per_m3: f64) -> f64 {
        (efficiency * q * area / heat_per_m3).max(0.0)
    }

    // ---------------------------------------------------------------- atmosphere

    pub fn mean_molecular_mass(composition: &BTreeMap<String, f64>) -> f64 {
        let total: f64 = composition.values().sum();
        if total <= 0.0 {
            return 0.028;
        }
        composition
            .iter()
            .map(|(s, f)| f * molar_mass(s))
            .sum::<f64>()
            / total
    }

    pub fn mean_heat_capacity(composition: &BTreeMap<String, f64>) -> f64 {
        let total: f64 = composition.values().sum();
        if total <= 0.0 {
            return 1_000.0;
        }
        // mass-weighted
        let mass_total: f64 = composition.iter().map(|(s, f)| f * molar_mass(s)).sum();
        composition
            .iter()
            .map(|(s, f)| f * molar_mass(s) * heat_capacity(s))
            .sum::<f64>()
            / mass_total.max(1e-12)
    }

    /// Column mass of each species, kg/m², from volume fractions.
    pub fn species_columns(
        composition: &BTreeMap<String, f64>,
        column_total: f64,
    ) -> BTreeMap<String, f64> {
        let mu = mean_molecular_mass(composition);
        let total: f64 = composition.values().sum::<f64>().max(1e-12);
        composition
            .iter()
            .map(|(s, f)| (s.clone(), column_total * (f / total) * molar_mass(s) / mu))
            .collect()
    }

    pub fn scale_height(t: f64, mu_kg_per_mol: f64, g: f64) -> f64 {
        const R_GAS: f64 = 8.314_462_618;
        R_GAS * t / (mu_kg_per_mol * g)
    }

    /// F22.
    pub fn equilibrium_temperature(l_watts: f64, d: f64, albedo: f64) -> f64 {
        (l_watts * (1.0 - albedo) / (16.0 * PI * SIGMA_SB * d * d)).powf(0.25)
    }

    /// Grey optical depth from the dry columns (F21).
    pub fn tau_dry(columns: &BTreeMap<String, f64>) -> f64 {
        let u = |s: &str| columns.get(s).copied().unwrap_or(0.0);
        TAU_CO2.0 * u("CO2").powf(TAU_CO2.1)
            + TAU_CH4.0 * u("CH4").powf(TAU_CH4.1)
            + TAU_N2_CIA.0 * u("N2").powf(TAU_N2_CIA.1)
    }

    /// Water vapour column at temperature `t`, kg/m², capped by the available inventory.
    pub fn water_column(t: f64, inventory_gel_m: f64) -> f64 {
        let u = WATER_COLUMN_REF * (WATER_L_OVER_RV * (1.0 / WATER_COLUMN_REF_T - 1.0 / t)).exp();
        u.min(WATER_COLUMN_MAX).min(inventory_gel_m * 1_000.0)
    }

    pub fn grey_surface_temperature(t_eq: f64, tau: f64) -> f64 {
        t_eq * (1.0 + 0.75 * tau).powf(0.25)
    }

    /// Surface temperature with the water vapour feedback: the smallest root of
    /// T = T_eq (1 + 3 tau(T) / 4)^(1/4) above T_eq, by scan and bisection.
    pub fn surface_temperature(t_eq: f64, tau_dry: f64, water_inventory_gel_m: f64) -> f64 {
        let tau = |t: f64| {
            tau_dry
                + if water_inventory_gel_m > 0.0 {
                    TAU_H2O.0 * water_column(t, water_inventory_gel_m).powf(TAU_H2O.1)
                } else {
                    0.0
                }
        };
        let f = |t: f64| t - grey_surface_temperature(t_eq, tau(t));
        let t0 = t_eq.max(1.0);
        if f(t0) >= 0.0 {
            return t0;
        }
        let mut lo = t0;
        let mut hi = t0;
        let step = (t0 * 0.02).max(0.5);
        for _ in 0..4000 {
            hi += step;
            if f(hi) >= 0.0 {
                break;
            }
            lo = hi;
        }
        if f(hi) < 0.0 {
            return hi;
        }
        for _ in 0..60 {
            let mid = 0.5 * (lo + hi);
            if f(mid) < 0.0 { lo = mid } else { hi = mid }
        }
        0.5 * (lo + hi)
    }

    /// F23: Jeans parameter of a species at temperature `t`.
    pub fn jeans_parameter(m: f64, r: f64, t: f64, molar_mass_kg: f64) -> f64 {
        const N_A: f64 = 6.022_140_76e23;
        let m_mol = molar_mass_kg / N_A;
        G * m * m_mol / (K_B * t * r)
    }

    pub fn exobase_temperature(t_s: f64, xuv: f64, f_co2: f64) -> f64 {
        t_s + EXOBASE_HEATING * xuv.powf(EXOBASE_XUV_EXPONENT) * (1.0 - EXOBASE_CO2_COOLING * f_co2)
    }

    /// F24: atmospheric crater cut-off diameter, metres.
    pub fn crater_cutoff(p: f64, g: f64) -> f64 {
        if p <= 0.0 {
            return 0.0;
        }
        let c = CUTOFF_EARTH_DIAMETER / (101_325.0 / 9.81f64).powf(CUTOFF_EXPONENT);
        c * (p / g).powf(CUTOFF_EXPONENT)
    }

    // ---------------------------------------------------------------- volatiles

    /// F25: saturation vapour pressure, Pa.
    pub fn vapour_pressure(v: &Volatile, t: f64) -> f64 {
        let l = if t >= v.t_triple { v.l_vap } else { v.l_sub };
        v.p_triple * ((l / v.r_specific) * (1.0 / v.t_triple - 1.0 / t)).exp()
    }

    /// Temperature at which the saturation pressure equals `p` (boiling above the triple
    /// point, frost point below).
    pub fn saturation_temperature(v: &Volatile, p: f64) -> f64 {
        let l = if p >= v.p_triple { v.l_vap } else { v.l_sub };
        1.0 / (1.0 / v.t_triple - (v.r_specific / l) * (p / v.p_triple).ln())
    }

    /// Phase state code of a species with a condensed inventory at surface conditions.
    pub fn phase_state(v: &Volatile, t: f64, p_total: f64, has_inventory: bool) -> f64 {
        if !has_inventory {
            return PHASE_GAS;
        }
        let t_melt = v.t_triple;
        if t < t_melt {
            return PHASE_SOLID;
        }
        if p_total < v.p_triple {
            return PHASE_GAS;
        }
        let t_boil = saturation_temperature(v, p_total);
        if t < t_boil {
            if t < t_melt + 0.05 * t_melt {
                PHASE_LIQUID_AND_SOLID
            } else {
                PHASE_LIQUID
            }
        } else {
            PHASE_GAS
        }
    }

    // ---------------------------------------------------------------- climate

    /// F32: Held–Hou cell edge in degrees.
    pub fn held_hou_edge_deg(g: f64, h_trop: f64, delta_h: f64, omega: f64, r: f64) -> f64 {
        if omega <= 0.0 {
            return 90.0;
        }
        let phi2 = 5.0 * g * h_trop * delta_h / (3.0 * omega * omega * r * r);
        phi2.max(0.0).sqrt().min(PI / 2.0).to_degrees()
    }

    /// F33: transport coefficient, W/m²/K.
    pub fn transport_coefficient(p_pa: f64, cp: f64, mu: f64, omega: f64) -> f64 {
        if p_pa <= 0.0 {
            return 0.0;
        }
        let d = D_EARTH
            * (p_pa / 101_325.0)
            * (cp / 1_004.0)
            * (0.028_97 / mu)
            * (7.2921e-5 / omega.max(1e-12)).powi(2);
        d.min(D_CAP_MULTIPLE * D_EARTH)
    }

    // ---------------------------------------------------------------- cratering

    /// F30: simple-to-complex transition diameter, metres.
    pub fn transition_diameter(g: f64, icy: bool) -> f64 {
        TRANSITION_LUNAR_KM
            * 1e3
            * (GRAVITY_MOON / g)
            * if icy { TRANSITION_ICE_FACTOR } else { 1.0 }
    }
}

use physics::*;

// ------------------------------------------------------------------ node functions: star

pub fn star_luminosity_reference(c: &Ctx) -> Value {
    mass_luminosity(c.f64("star.mass")).into()
}

pub fn star_main_sequence_lifetime(c: &Ctx) -> Value {
    main_sequence_lifetime_gyr(c.f64("star.mass"), c.f64("star.luminosity_reference")).into()
}

pub fn star_luminosity(c: &Ctx) -> Value {
    luminosity_at_age(
        c.f64("star.luminosity_reference"),
        c.f64("star.age"),
        c.f64("star.main_sequence_lifetime"),
    )
    .into()
}

pub fn star_radius(c: &Ctx) -> Value {
    stellar_radius(c.f64("star.mass")).into()
}

pub fn star_teff(c: &Ctx) -> Value {
    effective_temperature(c.f64("star.luminosity"), c.f64("star.radius")).into()
}

pub fn star_spectral_type(c: &Ctx) -> Value {
    spectral_type(c.f64("star.teff")).into()
}

pub fn star_xuv_relative(c: &Ctx) -> Value {
    xuv_relative(c.f64("star.age"), c.f64("orbit.star_distance") / AU).into()
}

pub fn star_luminosity_history(c: &Ctx) -> Value {
    let l_ref = c.f64("star.luminosity_reference");
    let t_ms = c.f64("star.main_sequence_lifetime");
    let star_age = c.f64("star.age");
    let body_age = c.f64("body.age");
    let offset = star_age - body_age;
    let n = 24;
    let points = (0..=n)
        .map(|k| {
            let t = body_age * k as f64 / n as f64;
            (t, luminosity_at_age(l_ref, offset + t, t_ms))
        })
        .collect();
    Curve::new(points).into()
}

// ------------------------------------------------------------------ orbit

pub fn orbit_primary_mass(c: &Ctx) -> Value {
    match c.str("orbit.primary") {
        "planet" => c.f64("parent.mass"),
        _ => c.f64("star.mass") * M_SUN,
    }
    .into()
}

pub fn orbit_period(c: &Ctx) -> Value {
    kepler_period(
        c.f64("orbit.semi_major_axis"),
        c.f64("orbit.primary_mass") + c.f64("body.mass"),
    )
    .into()
}

pub fn orbit_mean_motion(c: &Ctx) -> Value {
    (2.0 * PI / c.f64("orbit.period")).into()
}

pub fn orbit_star_distance(c: &Ctx) -> Value {
    match c.str("orbit.primary") {
        "planet" => c.f64("parent.semi_major_axis"),
        _ => c.f64("orbit.semi_major_axis"),
    }
    .into()
}

pub fn orbit_insolation(c: &Ctx) -> Value {
    let d = c.f64("orbit.star_distance");
    (c.f64("star.luminosity") * L_SUN / (4.0 * PI * d * d)).into()
}

// ------------------------------------------------------------------ body

pub fn body_mass(c: &Ctx) -> Value {
    let r = c.f64("body.radius");
    (4.0 / 3.0 * PI * r * r * r * c.f64("body.density")).into()
}

pub fn body_gravity(c: &Ctx) -> Value {
    gravity(c.f64("body.mass"), c.f64("body.radius")).into()
}

pub fn body_escape_velocity(c: &Ctx) -> Value {
    escape_velocity(c.f64("body.mass"), c.f64("body.radius")).into()
}

pub fn body_age(c: &Ctx) -> Value {
    (c.f64("star.age") - 0.05).max(0.01).into()
}

pub fn body_compression_factor(c: &Ctx) -> Value {
    compression_factor(c.f64("body.mass")).into()
}

pub fn body_uncompressed_density(c: &Ctx) -> Value {
    (c.f64("body.density") / c.f64("body.compression_factor")).into()
}

pub fn body_core_mass_fraction(c: &Ctx) -> Value {
    composition_fractions(c.f64("body.uncompressed_density"))
        .0
        .into()
}

pub fn body_ice_mass_fraction(c: &Ctx) -> Value {
    composition_fractions(c.f64("body.uncompressed_density"))
        .2
        .into()
}

pub fn body_rock_mass_fraction(c: &Ctx) -> Value {
    (1.0 - c.f64("body.core_mass_fraction") - c.f64("body.ice_mass_fraction"))
        .max(0.0)
        .into()
}

pub fn body_composition_class(c: &Ctx) -> Value {
    composition_class(
        c.f64("body.core_mass_fraction"),
        c.f64("body.ice_mass_fraction"),
    )
    .into()
}

pub fn body_moment_of_inertia_factor(c: &Ctx) -> Value {
    match c.str("body.composition_class") {
        "iron_rich" | "rock" => 0.33,
        "rock_ice" => 0.35,
        _ => 0.40,
    }
    .into()
}

pub fn body_ice_layer_thickness(c: &Ctx) -> Value {
    ice_layer_thickness(
        c.f64("body.radius"),
        c.f64("body.mass"),
        c.f64("body.ice_mass_fraction"),
    )
    .into()
}

pub fn body_oblateness(c: &Ctx) -> Value {
    oblateness(
        c.f64("rotation.angular_velocity"),
        c.f64("body.radius"),
        c.f64("body.mass"),
        c.f64("body.moment_of_inertia_factor"),
    )
    .into()
}

// ------------------------------------------------------------------ rotation

pub fn rotation_locking_timescale(c: &Ctx) -> Value {
    let m = c.f64("body.mass");
    let r = c.f64("body.radius");
    let moment = c.f64("body.moment_of_inertia_factor") * m * r * r;
    let omega = 2.0 * PI / c.f64("rotation.period");
    let t = locking_timescale(
        omega,
        c.f64("orbit.semi_major_axis"),
        moment,
        c.f64("interior.dissipation_q"),
        c.f64("orbit.primary_mass"),
        c.f64("interior.love_number_k2"),
        r,
    );
    (t / GYR).min(1e12).into()
}

pub fn rotation_mode(c: &Ctx) -> Value {
    if c.f64("rotation.locking_timescale") < c.f64("body.age") {
        "locked"
    } else {
        "free"
    }
    .into()
}

pub fn rotation_effective_period(c: &Ctx) -> Value {
    match c.str("rotation.mode") {
        "locked" => c.f64("orbit.period"),
        _ => c.f64("rotation.period"),
    }
    .into()
}

pub fn rotation_effective_obliquity(c: &Ctx) -> Value {
    match c.str("rotation.mode") {
        "locked" => 0.0,
        _ => c.f64("rotation.obliquity"),
    }
    .into()
}

pub fn rotation_angular_velocity(c: &Ctx) -> Value {
    (2.0 * PI / c.f64("rotation.effective_period")).into()
}

pub fn rotation_regime(c: &Ctx) -> Value {
    physics::rotation_regime(
        c.f64("body.gravity"),
        c.f64("atmosphere.scale_height"),
        c.f64("rotation.angular_velocity"),
        c.f64("body.radius"),
    )
    .into()
}

// ------------------------------------------------------------------ interior

fn is_icy(class: &str) -> bool {
    matches!(class, "rock_ice" | "ice_dominated")
}

pub fn interior_love_number_k2(c: &Ctx) -> Value {
    match c.str("body.composition_class") {
        "rock_ice" => K2_ROCK_ICE,
        "ice_dominated" => K2_ICE,
        _ => K2_ROCK,
    }
    .into()
}

pub fn interior_love_number_h2(c: &Ctx) -> Value {
    crate::stages::tidal::love_number_h2(
        c.f64("interior.love_number_k2"),
        c.f64("interior.ice_shell_thickness"),
        c.f64("body.radius"),
        is_icy(c.str("body.composition_class")),
    )
    .into()
}

pub fn tectonics_tidal_stress(c: &Ctx) -> Value {
    match crate::stages::tidal::TidalShell::from_ctx(c) {
        Some(shell) => shell.amplitude().into(),
        None => 0.0.into(),
    }
}

pub fn tectonics_ice_crack_density(c: &Ctx) -> Value {
    let stress = c.f64("tectonics.tidal_stress");
    if c.str("interior.tectonic_regime") != "ice_shell_ocean" || stress <= ICE_CRACK_STRESS {
        return 0.0.into();
    }
    (CRACK_DENSITY_FACTOR * (stress / ICE_TENSILE_STRENGTH)
        / c.f64("interior.ice_shell_thickness").max(1.0))
    .into()
}

pub fn tectonics_polar_thinning(c: &Ctx) -> Value {
    if c.str("interior.tectonic_regime") != "ice_shell_ocean" {
        return 0.0.into();
    }
    let ratio = c.f64("interior.ice_shell_thickness") / c.f64("body.radius").max(1.0);
    let focus = 1.0 - (-ratio / POLAR_FOCUS_SHELL_RATIO).exp();
    (POLAR_THINNING_MAX * (c.f64("interior.tidal_flux") / POLAR_FOCUS_FLUX).min(1.0) * focus).into()
}

pub fn tectonics_ice_resurfaced_fraction(c: &Ctx) -> Value {
    if c.str("interior.tectonic_regime") != "ice_shell_ocean" {
        return 0.0.into();
    }
    let ratio = c.f64("interior.ice_shell_thickness") / c.f64("body.radius").max(1.0);
    if ratio < WHOLESALE_SHELL_RATIO {
        // A thin shell renews itself by cracking everywhere instead.
        return 0.0.into();
    }
    // What the focusing takes to the pole does not spread as lanes.
    let diffuse = (1.0 - c.f64("tectonics.polar_thinning") / POLAR_THINNING_MAX).max(0.0);
    // The peak tidal flux of the body's history, F13 over the eccentricity history: the
    // pulse of a resonance episode (F50), or the steady flux.
    let e_now = c.f64("orbit.eccentricity");
    let e_max = c
        .curve("orbit.eccentricity_history")
        .points
        .iter()
        .map(|p| p.1)
        .fold(e_now, f64::max);
    let peak = if e_now > 0.0 {
        c.f64("interior.tidal_flux") * (e_max / e_now).powi(2)
    } else {
        c.f64("interior.tidal_flux")
    };
    // A shell renews itself at most once over: the renewal share saturates at one, and
    // the focusing's share is taken from that.
    let renewal = ((peak - ICE_RESURFACING_FLUX_FLOOR) / ICE_RESURFACING_FLUX).clamp(0.0, 1.0);
    (renewal * diffuse).clamp(0.0, 0.95).into()
}

pub fn tectonics_chaos_fraction(c: &Ctx) -> Value {
    if c.str("interior.tectonic_regime") != "ice_shell_ocean" {
        return 0.0.into();
    }
    let ratio = c.f64("interior.ice_shell_thickness") / c.f64("body.radius").max(1.0);
    if ratio >= WHOLESALE_SHELL_RATIO {
        return 0.0.into();
    }
    (CHAOS_FRACTION_MAX * (c.f64("interior.tidal_flux") / CHAOS_FLUX).min(1.0)).into()
}

pub fn volcanism_cryodome_count(c: &Ctx) -> Value {
    if c.str("interior.tectonic_regime") != "ice_shell_rock" {
        return 0.0.into();
    }
    let r_km3 = (c.f64("body.radius") / 1e3).powi(3);
    let (lo, hi) = (
        CRYODOME_VOLUME_SHARE.0 * r_km3,
        CRYODOME_VOLUME_SHARE.1 * r_km3,
    );
    let mean = (hi - lo) / (hi / lo).ln();
    // The magma rate is kept in m³/s; the volumes are km³.
    (CRYODOME_SHARE * c.f64("volcanism.rate") * GYR / 1e9 * c.f64("body.age") / mean)
        .round()
        .into()
}

pub fn volcanism_plume_fallout_rate(c: &Ctx) -> Value {
    if c.f64("tectonics.polar_thinning") <= POLAR_ACTIVE_THINNING
        || c.f64("tectonics.ice_crack_density") <= 0.0
    {
        return 0.0.into();
    }
    (PLUME_FALLOUT_RATE * (c.f64("interior.tidal_flux") / POLAR_FOCUS_FLUX).min(1.0)).into()
}

pub fn materials_salt(c: &Ctx) -> Value {
    if c.str("orbit.primary") == "planet"
        && c.f64("orbit.primary_mass") > MAGNETOSPHERE_JUPITER_MASS
    {
        "hydrated_salt".into()
    } else {
        "carbonate_salt".into()
    }
}

pub fn surface_trailing_irradiation(c: &Ctx) -> Value {
    if c.str("rotation.mode") != "locked" || c.str("orbit.primary") != "planet" {
        return 0.0.into();
    }
    let m = c.f64("orbit.primary_mass");
    if m > MAGNETOSPHERE_JUPITER_MASS {
        1.0.into()
    } else if m > MAGNETOSPHERE_SATURN_MASS {
        MAGNETOSPHERE_SATURN_SHARE.into()
    } else {
        0.0.into()
    }
}

/// When the lanes of grooved terrain form (F44): the body's history weighted by the
/// tidal flux above the resurfacing floor, as the mean time since formation and the
/// width of the uniform episode with the same variance (sqrt(12) sigma). The history's
/// samples are uniform in time, so plain sums serve.
fn resurfacing_episode(c: &Ctx) -> Option<(f64, f64)> {
    let curve = c.curve("interior.tidal_flux_history");
    let (mut w, mut wt, mut wtt) = (0.0, 0.0, 0.0);
    for p in &curve.points {
        let x = (p.1 - ICE_RESURFACING_FLUX_FLOOR).max(0.0);
        w += x;
        wt += x * p.0;
        wtt += x * p.0 * p.0;
    }
    if w <= 0.0 {
        return None;
    }
    let mean = wt / w;
    let variance = (wtt / w - mean * mean).max(0.0);
    Some((mean, (12.0 * variance).sqrt()))
}

pub fn tectonics_ice_resurfacing_duration(c: &Ctx) -> Value {
    if c.f64("tectonics.ice_resurfaced_fraction") <= 0.0 {
        return 0.0.into();
    }
    resurfacing_episode(c)
        .map_or(0.0, |(_, width)| width)
        .into()
}

pub fn tectonics_ice_resurfacing_age(c: &Ctx) -> Value {
    if c.f64("tectonics.ice_resurfaced_fraction") <= 0.0 {
        return 0.0.into();
    }
    let age = c.f64("body.age");
    resurfacing_episode(c)
        .map_or(0.0, |(mean, _)| (age - mean).max(0.0))
        .into()
}

pub fn tectonics_nsr_period(c: &Ctx) -> Value {
    if c.f64("tectonics.ice_crack_density") > 0.0 {
        (NSR_PERIOD_RETENTION_FACTOR * c.f64("cratering.retention_age")).into()
    } else {
        0.0.into()
    }
}

pub fn interior_dissipation_q(c: &Ctx) -> Value {
    match c.str("body.composition_class") {
        "rock_ice" => Q_ROCK_ICE,
        "ice_dominated" => Q_ICE,
        _ => {
            if c.str("volatiles.working_liquid") == "water" {
                Q_ROCK_SURFACE_OCEAN
            } else {
                Q_ROCK
            }
        }
    }
    .into()
}

pub fn interior_radiogenic_heat_per_kg(c: &Ctx) -> Value {
    radiogenic_per_kg(c.f64("body.age")).into()
}

pub fn interior_radiogenic_power(c: &Ctx) -> Value {
    (c.f64("interior.radiogenic_heat_per_kg")
        * c.f64("body.mass")
        * c.f64("body.rock_mass_fraction"))
    .into()
}

pub fn interior_radiogenic_flux(c: &Ctx) -> Value {
    let r = c.f64("body.radius");
    (c.f64("interior.radiogenic_power") / (4.0 * PI * r * r)).into()
}

pub fn interior_tidal_power(c: &Ctx) -> Value {
    tidal_power(
        c.f64("interior.love_number_k2"),
        c.f64("interior.dissipation_q"),
        c.f64("orbit.primary_mass"),
        c.f64("body.radius"),
        c.f64("orbit.mean_motion"),
        c.f64("orbit.eccentricity"),
        c.f64("orbit.semi_major_axis"),
    )
    .into()
}

pub fn interior_tidal_flux(c: &Ctx) -> Value {
    let r = c.f64("body.radius");
    (c.f64("interior.tidal_power") / (4.0 * PI * r * r)).into()
}

pub fn interior_secular_cooling_fraction(c: &Ctx) -> Value {
    secular_fraction(c.f64("body.radius")).into()
}

pub fn interior_surface_heat_flux(c: &Ctx) -> Value {
    (c.f64("interior.radiogenic_flux") * (1.0 + c.f64("interior.secular_cooling_fraction"))
        + c.f64("interior.tidal_flux"))
    .into()
}

/// The thermal evolution's inputs from the graph context, mirroring
/// `ThermalInputs::from_graph`.
fn thermal_inputs(c: &Ctx) -> crate::thermal::ThermalInputs {
    let regime = c.str("interior.tectonic_regime");
    crate::thermal::ThermalInputs {
        age_gyr: c.f64("body.age"),
        radius: c.f64("body.radius"),
        gravity: c.f64("body.gravity"),
        uncompressed_density: c.f64("body.uncompressed_density"),
        rock_mass: c.f64("body.mass") * c.f64("body.rock_mass_fraction"),
        secular_fraction: c.f64("interior.secular_cooling_fraction"),
        tidal_flux_now: c.f64("interior.tidal_flux"),
        eccentricity_now: c.f64("orbit.eccentricity"),
        eccentricity_history: c.curve("orbit.eccentricity_history").clone(),
        icy: is_icy(c.str("body.composition_class")),
        surface_temperature: c.f64("atmosphere.surface_temperature_mean"),
        melting_point: ice_melting_point(c),
        crust_strength: c.f64("interior.crust_strength"),
        magma_efficiency_base: match regime {
            "mobile_lid" => MAGMA_EFFICIENCY_MOBILE,
            "episodic_lid" => MAGMA_EFFICIENCY_EPISODIC,
            "ice_shell_ocean" | "ice_shell_rock" => MAGMA_EFFICIENCY_ICE,
            _ => MAGMA_EFFICIENCY_STAGNANT,
        },
    }
}

/// F50: the present eccentricity, with the pulse of a resonance episode centred
/// `orbit.resonance_age` before the present and `orbit.resonance_duration` wide at its
/// base, rising over the first fifth, holding at the forced value, falling over the
/// last fifth; constant when there is none.
pub fn orbit_eccentricity_history(c: &Ctx) -> Value {
    let e = c.f64("orbit.eccentricity");
    let age = c.f64("body.age");
    let peak = c.f64("orbit.resonance_eccentricity");
    let before = c.f64("orbit.resonance_age");
    let width = c.f64("orbit.resonance_duration");
    if peak <= e || width <= 0.0 || before <= 0.0 || before > age {
        return Curve::constant(e).into();
    }
    let t_peak = age - before;
    let half = 0.5 * width;
    let pulse = |t: f64| {
        let x = (t - t_peak).abs() / half;
        e + (peak - e) * ((1.0 - x) / RESONANCE_RAMP_SHARE).clamp(0.0, 1.0)
    };
    let hold = half * (1.0 - RESONANCE_RAMP_SHARE);
    let mut times = vec![
        0.0,
        t_peak - half,
        t_peak - hold,
        t_peak + hold,
        t_peak + half,
        age,
    ];
    times.retain(|t| (0.0..=age).contains(t));
    times.sort_by(|a, b| a.total_cmp(b));
    times.dedup();
    Curve::new(times.into_iter().map(|t| (t, pulse(t))).collect()).into()
}

pub fn interior_tidal_flux_history(c: &Ctx) -> Value {
    crate::thermal::ThermalHistory::new(&thermal_inputs(c))
        .curve(|s| s.tidal_flux)
        .into()
}

pub fn interior_heat_flux_history(c: &Ctx) -> Value {
    crate::thermal::ThermalHistory::new(&thermal_inputs(c))
        .curve(|s| s.heat_flux)
        .into()
}

pub fn interior_mantle_temperature_history(c: &Ctx) -> Value {
    crate::thermal::ThermalHistory::new(&thermal_inputs(c))
        .curve(|s| s.mantle_temperature)
        .into()
}

pub fn interior_lithosphere_history(c: &Ctx) -> Value {
    crate::thermal::ThermalHistory::new(&thermal_inputs(c))
        .curve(|s| s.lithosphere)
        .into()
}

pub fn interior_relief_limit_history(c: &Ctx) -> Value {
    crate::thermal::ThermalHistory::new(&thermal_inputs(c))
        .curve(|s| s.relief_limit)
        .into()
}

pub fn interior_contraction_history(c: &Ctx) -> Value {
    crate::thermal::ThermalHistory::new(&thermal_inputs(c))
        .curve(|s| s.contraction)
        .into()
}

pub fn volcanism_rate(c: &Ctx) -> Value {
    c.f64("interior.magma_production_rate").into()
}

pub fn volcanism_resurfacing_rate(c: &Ctx) -> Value {
    let r = c.f64("body.radius");
    let area = 4.0 * PI * r * r;
    (c.f64("volcanism.rate") / area * MYR).into()
}

pub fn volcanism_history(c: &Ctx) -> Value {
    crate::thermal::ThermalHistory::new(&thermal_inputs(c))
        .curve(|s| s.magma_rate)
        .into()
}

fn ice_melting_point(c: &Ctx) -> f64 {
    if c.map_get("volatiles.inventory", "NH3") > 0.0 {
        T_EUTECTIC_AMMONIA_WATER
    } else {
        T_MELT_WATER
    }
}

pub fn interior_mantle_temperature(c: &Ctx) -> Value {
    if is_icy(c.str("body.composition_class")) {
        ice_melting_point(c)
    } else {
        MANTLE_T_BASE
            + MANTLE_T_SCALE
                * (c.f64("interior.surface_heat_flux") / EARTH_HEAT_FLUX)
                    .max(0.0)
                    .sqrt()
    }
    .into()
}

pub fn interior_lithosphere_thickness(c: &Ctx) -> Value {
    let t_s = c.f64("atmosphere.surface_temperature_mean");
    let t_m = c.f64("interior.mantle_temperature");
    let q = c.f64("interior.surface_heat_flux");
    let dt = (t_m - t_s).max(1.0);
    let z = if is_icy(c.str("body.composition_class")) {
        conductive_thickness(ice_conductivity(0.5 * (t_m + t_s)), dt, q)
    } else {
        conductive_thickness(K_ROCK, dt, MANTLE_FLUX_FRACTION * q)
    };
    z.min(c.f64("body.radius")).into()
}

pub fn interior_elastic_thickness(c: &Ctx) -> Value {
    let t_s = c.f64("atmosphere.surface_temperature_mean");
    let t_m = c.f64("interior.mantle_temperature");
    let iso = if is_icy(c.str("body.composition_class")) {
        ELASTIC_ISOTHERM_ICE
    } else {
        ELASTIC_ISOTHERM_ROCK
    };
    let frac = ((iso - t_s) / (t_m - t_s).max(1.0)).clamp(0.0, 1.0);
    (c.f64("interior.lithosphere_thickness") * frac).into()
}

pub fn interior_crust_strength(c: &Ctx) -> Value {
    if is_icy(c.str("body.composition_class")) {
        STRENGTH_ICE
    } else {
        STRENGTH_ROCK
    }
    .into()
}

pub fn interior_max_relief(c: &Ctx) -> Value {
    let icy = is_icy(c.str("body.composition_class"));
    let rho = if icy {
        DENSITY_CRUST_ICE
    } else {
        DENSITY_CRUST_ROCK
    };
    relief_limit(c.f64("interior.crust_strength"), rho, c.f64("body.gravity")).into()
}

pub fn interior_contraction_total(c: &Ctx) -> Value {
    contraction(
        c.f64("body.radius"),
        c.f64("body.uncompressed_density"),
        c.f64("body.age"),
    )
    .into()
}

pub fn interior_ice_shell_thickness(c: &Ctx) -> Value {
    if is_icy(c.str("body.composition_class")) {
        c.f64("interior.lithosphere_thickness")
    } else {
        0.0
    }
    .into()
}

pub fn interior_has_subsurface_ocean(c: &Ctx) -> Value {
    (is_icy(c.str("body.composition_class"))
        && c.f64("body.ice_layer_thickness") > c.f64("interior.ice_shell_thickness"))
    .into()
}

pub fn interior_tectonic_regime(c: &Ctx) -> Value {
    let class = c.str("body.composition_class");
    let q = c.f64("interior.surface_heat_flux");
    let r = c.f64("body.radius");
    let regime = if is_icy(class) {
        if c.bool("interior.has_subsurface_ocean") {
            "ice_shell_ocean"
        } else {
            "ice_shell_rock"
        }
    } else if c.str("volatiles.working_liquid") == "water"
        && q >= MOBILE_LID_MIN_FLUX
        && r >= MOBILE_LID_MIN_RADIUS
    {
        "mobile_lid"
    } else if q >= EPISODIC_LID_MIN_FLUX && r >= EPISODIC_LID_MIN_RADIUS {
        "episodic_lid"
    } else {
        "stagnant_lid"
    };
    regime.into()
}

pub fn interior_magma_production_rate(c: &Ctx) -> Value {
    let regime = c.str("interior.tectonic_regime");
    let icy = is_icy(c.str("body.composition_class"));
    let base = match regime {
        "mobile_lid" => MAGMA_EFFICIENCY_MOBILE,
        "episodic_lid" => MAGMA_EFFICIENCY_EPISODIC,
        "ice_shell_ocean" | "ice_shell_rock" => MAGMA_EFFICIENCY_ICE,
        _ => MAGMA_EFFICIENCY_STAGNANT,
    };
    let tidal = (MAGMA_EFFICIENCY_TIDAL_MAX
        * (c.f64("interior.tidal_flux") / MAGMA_TIDAL_FLUX_SCALE).min(1.0))
    .max(0.0);
    let efficiency = base.max(tidal);
    let r = c.f64("body.radius");
    let heat = if icy {
        MAGMA_HEAT_PER_M3_ICE
    } else {
        MAGMA_HEAT_PER_M3_ROCK
    };
    magma_rate(
        c.f64("interior.surface_heat_flux"),
        4.0 * PI * r * r,
        efficiency,
        heat,
    )
    .into()
}

// ------------------------------------------------------------------ atmosphere

pub fn atmosphere_mean_molecular_mass(c: &Ctx) -> Value {
    mean_molecular_mass(c.map("atmosphere.composition")).into()
}

pub fn atmosphere_column_mass(c: &Ctx) -> Value {
    (c.f64("atmosphere.surface_pressure") / c.f64("body.gravity")).into()
}

pub fn atmosphere_bond_albedo(c: &Ctx) -> Value {
    let p = c.f64("atmosphere.surface_pressure");
    let class = c.str("body.composition_class");
    let a = if p < 10.0 {
        match class {
            "rock_ice" => ALBEDO_AIRLESS_ROCK_ICE,
            "ice_dominated" => ALBEDO_AIRLESS_ICE,
            _ => ALBEDO_AIRLESS_ROCK,
        }
    } else if p < THIN_ATMOSPHERE_PRESSURE {
        ALBEDO_THIN_ATMOSPHERE
    } else if p < THICK_ATMOSPHERE_PRESSURE {
        ALBEDO_EARTHLIKE
    } else {
        ALBEDO_THICK_ATMOSPHERE
    };
    a.into()
}

pub fn atmosphere_equilibrium_temperature(c: &Ctx) -> Value {
    equilibrium_temperature(
        c.f64("star.luminosity") * L_SUN,
        c.f64("orbit.star_distance"),
        c.f64("atmosphere.bond_albedo"),
    )
    .into()
}

pub fn atmosphere_surface_temperature_mean(c: &Ctx) -> Value {
    let t_eq = c.f64("atmosphere.equilibrium_temperature");
    let p = c.f64("atmosphere.surface_pressure");
    if p <= 0.0 {
        return t_eq.into();
    }
    let columns = species_columns(
        c.map("atmosphere.composition"),
        c.f64("atmosphere.column_mass"),
    );
    let water = c.map_get("volatiles.inventory", "H2O")
        + if c.map_get("atmosphere.composition", "H2O") > 0.0 {
            1.0
        } else {
            0.0
        };
    surface_temperature(t_eq, tau_dry(&columns), water).into()
}

pub fn atmosphere_greenhouse_optical_depth(c: &Ctx) -> Value {
    let t_eq = c.f64("atmosphere.equilibrium_temperature");
    let t_s = c.f64("atmosphere.surface_temperature_mean");
    ((4.0 / 3.0) * ((t_s / t_eq).powi(4) - 1.0)).max(0.0).into()
}

pub fn atmosphere_scale_height(c: &Ctx) -> Value {
    scale_height(
        c.f64("atmosphere.surface_temperature_mean"),
        c.f64("atmosphere.mean_molecular_mass"),
        c.f64("body.gravity"),
    )
    .into()
}

pub fn atmosphere_exobase_temperature(c: &Ctx) -> Value {
    exobase_temperature(
        c.f64("atmosphere.surface_temperature_mean"),
        c.f64("star.xuv_relative"),
        c.map_get("atmosphere.composition", "CO2"),
    )
    .into()
}

pub fn atmosphere_jeans_parameter(c: &Ctx) -> Value {
    let m = c.f64("body.mass");
    let r = c.f64("body.radius");
    let t = c.f64("atmosphere.exobase_temperature");
    let mut out = BTreeMap::new();
    for (s, f) in c.map("atmosphere.composition") {
        if *f > 0.0 {
            out.insert(s.clone(), jeans_parameter(m, r, t, molar_mass(s)));
        }
    }
    Value::Map(out)
}

pub fn atmosphere_retention_warning(c: &Ctx) -> Value {
    let escaping: Vec<String> = c
        .map("atmosphere.jeans_parameter")
        .iter()
        .filter(|(_, l)| **l < JEANS_RETENTION_THRESHOLD)
        .map(|(s, l)| format!("{s} ({l:.0})"))
        .collect();
    if escaping.is_empty() {
        String::new()
    } else {
        escaping.join(", ")
    }
    .into()
}

pub fn atmosphere_oxidant_availability(c: &Ctx) -> Value {
    let p_bar = c.f64("atmosphere.surface_pressure") / 1e5;
    ((c.map_get("atmosphere.composition", "O2") + 0.1 * c.map_get("atmosphere.composition", "H2O"))
        * p_bar)
        .clamp(0.0, 1.0)
        .into()
}

pub fn atmosphere_haze_production(c: &Ctx) -> Value {
    let xuv = c.f64("star.xuv_relative");
    let p_bar = c.f64("atmosphere.surface_pressure") / 1e5;
    let tholin = c.map_get("atmosphere.composition", "CH4") * xuv * 50.0 * p_bar.min(2.0);
    let sulfuric = c.map_get("atmosphere.composition", "SO2") * xuv * 50.0 * p_bar.min(2.0);
    (tholin + sulfuric).into()
}

pub fn atmosphere_crater_cutoff_diameter(c: &Ctx) -> Value {
    crater_cutoff(c.f64("atmosphere.surface_pressure"), c.f64("body.gravity")).into()
}

pub fn atmosphere_condensables(c: &Ctx) -> Value {
    let comp = c.map("atmosphere.composition");
    let inv = c.map("volatiles.inventory");
    let mut out = Vec::new();
    for s in VOLATILE_SPECIES {
        let present = comp.get(*s).copied().unwrap_or(0.0) > 1e-4
            || inv.get(*s).copied().unwrap_or(0.0) > 0.0;
        if present {
            out.push((*s).to_string());
        }
    }
    Value::List(out)
}

// ------------------------------------------------------------------ volatiles

pub fn volatiles_inventory(c: &Ctx) -> Value {
    let mut m = BTreeMap::new();
    if is_icy(c.str("body.composition_class")) {
        m.insert("H2O".to_string(), c.f64("body.ice_layer_thickness"));
    }
    Value::Map(m)
}

pub fn volatiles_phase_state(c: &Ctx) -> Value {
    let t = c.f64("atmosphere.surface_temperature_mean");
    let p = c.f64("atmosphere.surface_pressure");
    let inv = c.map("volatiles.inventory");
    let mut out = BTreeMap::new();
    for s in VOLATILE_SPECIES {
        let Some(v) = volatile(s) else { continue };
        let has = inv.get(*s).copied().unwrap_or(0.0) > 0.0;
        let in_air = c.map_get("atmosphere.composition", s) > 0.0;
        if !has && !in_air {
            out.insert((*s).to_string(), PHASE_ABSENT);
        } else {
            out.insert((*s).to_string(), phase_state(v, t, p, has));
        }
    }
    Value::Map(out)
}

pub fn volatiles_phase_summary(c: &Ctx) -> Value {
    let parts: Vec<String> = c
        .map("volatiles.phase_state")
        .iter()
        .filter(|(_, code)| **code != PHASE_ABSENT)
        .map(|(s, code)| format!("{s}: {}", phase_name(*code)))
        .collect();
    parts.join("; ").into()
}

pub fn volatiles_working_liquid(c: &Ctx) -> Value {
    let states = c.map("volatiles.phase_state");
    let inv = c.map("volatiles.inventory");
    let liquid = |s: &str| {
        let code = states.get(s).copied().unwrap_or(PHASE_ABSENT);
        code == PHASE_LIQUID || code == PHASE_LIQUID_AND_SOLID
    };
    let amount = |s: &str| inv.get(s).copied().unwrap_or(0.0);
    let mut best: Option<(&str, f64)> = None;
    for s in VOLATILE_SPECIES {
        if liquid(s) {
            let a = amount(s);
            if best.is_none_or(|(_, b)| a > b) {
                best = Some((s, a));
            }
        }
    }
    let name = match best.map(|b| b.0) {
        Some("H2O") => {
            if amount("NH3") > 0.0 && liquid("NH3") {
                "ammonia_water"
            } else {
                "water"
            }
        }
        Some("CH4") | Some("C2H6") => "methane_ethane",
        Some("N2") => "nitrogen",
        Some("NH3") => "ammonia_water",
        Some("SO2") => "sulfur_dioxide",
        Some("CO2") => "carbon_dioxide",
        _ => "none",
    };
    name.into()
}

pub fn volatiles_frost_points(c: &Ctx) -> Value {
    let p = c.f64("atmosphere.surface_pressure");
    let comp = c.map("atmosphere.composition");
    let total: f64 = comp.values().sum::<f64>().max(1e-12);
    let mut out = BTreeMap::new();
    for s in VOLATILE_SPECIES {
        let Some(v) = volatile(s) else { continue };
        let partial = p * comp.get(*s).copied().unwrap_or(0.0) / total;
        let pp = if partial > 0.0 {
            partial.min(v.p_triple * 0.999)
        } else {
            1.0e-3
        };
        out.insert((*s).to_string(), saturation_temperature(v, pp));
    }
    Value::Map(out)
}

pub fn volatiles_cap_mode(c: &Ctx) -> Value {
    let any_solid = c
        .map("volatiles.phase_state")
        .values()
        .any(|code| *code == PHASE_SOLID || *code == PHASE_LIQUID_AND_SOLID);
    let mode = if !any_solid {
        "none"
    } else if c.str("rotation.mode") == "locked" {
        "hemispheric"
    } else if c.f64("rotation.effective_obliquity") > HIGH_OBLIQUITY_DEG {
        "equatorial"
    } else {
        "polar"
    };
    mode.into()
}

// ------------------------------------------------------------------ tectonics (subset)

pub fn tectonics_dichotomy_amplitude(c: &Ctx) -> Value {
    if c.str("tectonics.dichotomy_mode") == "none" {
        return 0.0.into();
    }
    match c.str("interior.tectonic_regime") {
        "mobile_lid" => 0.0,
        "episodic_lid" => 5_000.0,
        "stagnant_lid" => 20_000.0,
        _ => 2_000.0,
    }
    .into()
}

pub fn tectonics_initial_crust_thickness(c: &Ctx) -> Value {
    if is_icy(c.str("body.composition_class")) {
        20_000.0
    } else {
        40_000.0
    }
    .into()
}

// ------------------------------------------------------------------ cratering (scaling subset)

fn outer_system(c: &Ctx) -> bool {
    c.f64("orbit.star_distance") / AU > OUTER_SYSTEM_AU
}

pub fn cratering_production_function(c: &Ctx) -> Value {
    let d = c.f64("orbit.star_distance") / AU;
    if outer_system(c) {
        "outer_system"
    } else if d > 1.2 {
        "neukum_mars"
    } else {
        "neukum_lunar"
    }
    .into()
}

pub fn cratering_impactor_density(c: &Ctx) -> Value {
    if outer_system(c) {
        IMPACTOR_DENSITY_OUTER
    } else {
        IMPACTOR_DENSITY_INNER
    }
    .into()
}

pub fn cratering_target_density(c: &Ctx) -> Value {
    if is_icy(c.str("body.composition_class")) {
        DENSITY_CRUST_ICE
    } else {
        2_700.0
    }
    .into()
}

pub fn cratering_impact_velocity(c: &Ctx) -> Value {
    let d = c.f64("orbit.star_distance");
    let v_orb = orbital_speed(c.f64("star.mass") * M_SUN, d);
    let factor = (VINF_FACTOR_INNER * (AU / d).sqrt()).max(VINF_FACTOR_MIN);
    let v_inf = factor * v_orb;
    let v_esc = c.f64("body.escape_velocity");
    let v2 = if c.str("orbit.primary") == "planet" {
        let v_esc_parent = (2.0 * G * c.f64("parent.mass") / c.f64("orbit.semi_major_axis")).sqrt();
        v_inf * v_inf + v_esc_parent * v_esc_parent + v_esc * v_esc
    } else {
        v_inf * v_inf + v_esc * v_esc
    };
    v2.sqrt().into()
}

pub fn cratering_flux_relative_to_moon(c: &Ctx) -> Value {
    let d = c.f64("orbit.star_distance") / AU;
    let flux = if outer_system(c) {
        if c.str("orbit.primary") == "planet" {
            let v_orb = orbital_speed(c.f64("star.mass") * M_SUN, c.f64("orbit.star_distance"));
            let v_inf = VINF_FACTOR_MIN * v_orb;
            let v_esc_parent2 = 2.0 * G * c.f64("parent.mass") / c.f64("orbit.semi_major_axis");
            0.1 * (1.0 + v_esc_parent2 / (v_inf * v_inf))
        } else {
            0.1
        }
    } else if d > 1.2 {
        2.0
    } else {
        1.0
    };
    flux.into()
}

pub fn cratering_simple_complex_transition(c: &Ctx) -> Value {
    // The transition is a property of the target the crater forms in: the surface
    // material, not the bulk composition.
    transition_diameter(
        c.f64("body.gravity"),
        crate::stages::materials::is_ice(c.str("materials.bedrock")),
    )
    .into()
}

pub fn cratering_peak_ring_onset(c: &Ctx) -> Value {
    (PEAK_RING_MULTIPLE * c.f64("cratering.simple_complex_transition")).into()
}

pub fn cratering_multiring_onset(c: &Ctx) -> Value {
    (MULTIRING_MULTIPLE * c.f64("cratering.simple_complex_transition")).into()
}

pub fn cratering_diffusivity(c: &Ctx) -> Value {
    let mut k = DIFFUSIVITY_LUNAR * (c.f64("body.gravity") / GRAVITY_MOON).sqrt();
    if c.f64("atmosphere.surface_pressure") > DIFFUSIVITY_ATMOSPHERE_PRESSURE {
        k *= DIFFUSIVITY_ATMOSPHERE_MULTIPLE;
    }
    k.into()
}

pub fn cratering_ray_fading_timescale(c: &Ctx) -> Value {
    if c.f64("atmosphere.surface_pressure") > DIFFUSIVITY_ATMOSPHERE_PRESSURE {
        RAY_FADING_GYR_ATMOSPHERE
    } else {
        let d = c.f64("orbit.star_distance") / AU;
        (RAY_FADING_GYR_1AU * d * d).min(RAY_FADING_MAX_GYR)
    }
    .into()
}

fn production_function(c: &Ctx) -> crate::crater::ProductionFunction {
    crate::crater::ProductionFunction::by_name(c.str("cratering.production_function"))
}

/// F31: the early rate at formation relative to the steady rate, from the chronology.
pub fn cratering_early_flux_ratio(c: &Ctx) -> Value {
    crate::crater::production::early_flux_ratio(production_function(c).chronology).into()
}

/// F31: e-folding time of the early bombardment, from the chronology.
pub fn cratering_early_flux_decay(c: &Ctx) -> Value {
    (1.0 / production_function(c).chronology.1).into()
}

pub fn cratering_flux_history(c: &Ctx) -> Value {
    // Production rate over the body's life relative to the present rate.
    let age = c.f64("body.age");
    let flux = crate::crater::FluxHistory {
        steady: production_function(c).chronology.2,
        early_ratio: c.f64("cratering.early_flux_ratio"),
        early_decay: c.f64("cratering.early_flux_decay"),
        record_start: c.f64("surface.solidification_time"),
    };
    let now = flux.rate1(age).max(1e-30);
    let n = 48;
    let points = (0..=n)
        .map(|k| {
            let t = age * k as f64 / n as f64;
            (t, flux.rate1(t) / now)
        })
        .collect();
    Curve::new(points).into()
}

/// Time after formation when the crust is solid and rigid enough to keep large relief.
/// Provisional scaling until the Phase 3 heat budget derives it.
pub fn surface_solidification_time(c: &Ctx) -> Value {
    let r = c.f64("body.radius");
    let ice = if is_icy(c.str("body.composition_class")) {
        SOLIDIFICATION_ICE_FACTOR
    } else {
        1.0
    };
    let t = SOLIDIFICATION_REFERENCE_GYR * (r / SOLIDIFICATION_REFERENCE_RADIUS).sqrt() * ice;
    t.clamp(SOLIDIFICATION_RANGE_GYR.0, SOLIDIFICATION_RANGE_GYR.1)
        .into()
}

// ------------------------------------------------------------------ climate (zero-dimensional subset)

pub fn climate_greenhouse_offset(c: &Ctx) -> Value {
    (c.f64("atmosphere.surface_temperature_mean") - c.f64("atmosphere.equilibrium_temperature"))
        .into()
}

pub fn climate_hadley_edge_latitude(c: &Ctx) -> Value {
    held_hou_edge_deg(
        c.f64("body.gravity"),
        TROPOPAUSE_SCALE_HEIGHTS * c.f64("atmosphere.scale_height"),
        HELD_HOU_DELTA_H,
        c.f64("rotation.angular_velocity"),
        c.f64("body.radius"),
    )
    .into()
}

pub fn climate_transport_coefficient(c: &Ctx) -> Value {
    let comp = c.map("atmosphere.composition");
    transport_coefficient(
        c.f64("atmosphere.surface_pressure"),
        mean_heat_capacity(comp),
        mean_molecular_mass(comp),
        c.f64("rotation.angular_velocity"),
    )
    .into()
}

pub fn climate_lapse_rate(c: &Ctx) -> Value {
    let cp = mean_heat_capacity(c.map("atmosphere.composition"));
    let moist = if c.str("volatiles.working_liquid") != "none" {
        MOIST_LAPSE_FACTOR
    } else {
        1.0
    };
    (c.f64("body.gravity") / cp * moist).into()
}

// ------------------------------------------------------------------ tectonics (contraction)

fn io_regime(c: &Ctx) -> bool {
    c.str("body.composition_class") == "rock" && c.f64("interior.tidal_flux") > IO_REGIME_TIDAL_FLUX
}

pub fn volcanism_patera_count(c: &Ctx) -> Value {
    if !io_regime(c) {
        return 0.0.into();
    }
    (IO_PATERAE_PER_TW * c.f64("interior.tidal_power") / 1e12)
        .round()
        .into()
}

pub fn tectonics_mountain_count(c: &Ctx) -> Value {
    if !io_regime(c) {
        return 0.0.into();
    }
    (IO_MOUNTAINS_PER_TW * c.f64("interior.tidal_power") / 1e12)
        .round()
        .into()
}

pub fn tectonics_resurfacing_event_age(c: &Ctx) -> Value {
    if c.str("interior.tectonic_regime") != "episodic_lid" {
        return 0.0.into();
    }
    c.f64("cratering.retention_age").into()
}

pub fn tectonics_tessera_fraction(c: &Ctx) -> Value {
    if c.str("interior.tectonic_regime") == "episodic_lid" {
        TESSERA_FRACTION.into()
    } else {
        0.0.into()
    }
}

pub fn volcanism_corona_count(c: &Ctx) -> Value {
    if c.str("interior.tectonic_regime") != "episodic_lid" {
        return 0.0.into();
    }
    let r = c.f64("body.radius");
    (CORONAE_PER_AREA * 4.0 * PI * r * r).round().into()
}

pub fn volcanism_dome_count(c: &Ctx) -> Value {
    if c.str("interior.tectonic_regime") != "episodic_lid" {
        return 0.0.into();
    }
    let r = c.f64("body.radius");
    (DOMES_PER_AREA * 4.0 * PI * r * r).round().into()
}

pub fn tectonics_plume_count(c: &Ctx) -> Value {
    let regime = c.str("interior.tectonic_regime");
    if !matches!(regime, "stagnant_lid" | "episodic_lid") || c.f64("body.radius") < PLUME_MIN_RADIUS
    {
        return 0.0.into();
    }
    (c.f64("interior.surface_heat_flux") / PLUME_FLUX_PER_PLUME)
        .floor()
        .clamp(0.0, 3.0)
        .into()
}

pub fn tectonics_rise_diameter(c: &Ctx) -> Value {
    if c.f64("tectonics.plume_rise_count") < 1.0 {
        return 0.0.into();
    }
    let fraction = if c.str("interior.tectonic_regime") == "episodic_lid" {
        PLUME_RISE_DIAMETER_FRACTION_EPISODIC
    } else {
        PLUME_RISE_DIAMETER_FRACTION
    };
    (fraction * c.f64("body.radius")).into()
}

pub fn tectonics_rise_height(c: &Ctx) -> Value {
    if c.f64("tectonics.plume_rise_count") < 1.0 {
        return 0.0.into();
    }
    (RISE_HEIGHT_FRACTION * c.f64("interior.max_relief")).into()
}

pub fn tectonics_rift_extension(c: &Ctx) -> Value {
    if c.f64("tectonics.plume_rise_count") >= 1.0 {
        (RIFT_EXTENSION_PER_UPLIFT * c.f64("tectonics.plume_rise_height")).into()
    } else {
        0.0.into()
    }
}

pub fn volcanism_province_count(c: &Ctx) -> Value {
    (c.f64("body.radius") / MAGMA_PROVINCE_RADIUS_PER)
        .round()
        .clamp(1.0, 3.0)
        .into()
}

pub fn volcanism_shield_count(c: &Ctx) -> Value {
    if c.f64("tectonics.plume_rise_count") >= 1.0 {
        SHIELDS_PER_RISE.into()
    } else {
        0.0.into()
    }
}

pub fn tectonics_contraction_strain(c: &Ctx) -> Value {
    (2.0 * c.f64("interior.contraction_total") / c.f64("body.radius")).into()
}

pub fn tectonics_scarp_total_shortening(c: &Ctx) -> Value {
    (c.f64("tectonics.contraction_strain") * 2.0 * PI * c.f64("body.radius")).into()
}

// ------------------------------------------------------------------ cratering (population)

/// Mean cell width of a cubesphere face size on this body, metres.
fn cell_width(radius: f64, face: f64) -> f64 {
    2.0 * PI * radius / (4.0 * face.max(1.0))
}

pub fn cratering_loop_threshold(c: &Ctx) -> Value {
    (LOOP_THRESHOLD_CELLS * cell_width(c.f64("body.radius"), LOOP_REFERENCE_FACE)).into()
}

pub fn cratering_min_diameter(c: &Ctx) -> Value {
    let by_resolution = FINISHING_MIN_CELLS
        * cell_width(c.f64("body.radius"), c.f64("detail.finishing_resolution"));
    by_resolution
        .max(c.f64("atmosphere.crater_cutoff_diameter"))
        .into()
}

pub fn cratering_retention_age(c: &Ctx) -> Value {
    let age = c.f64("body.age");
    let regime = c.str("interior.tectonic_regime");
    let p = c.f64("atmosphere.surface_pressure");
    let mut retention = match regime {
        "mobile_lid" => RETENTION_MOBILE_LID_GYR,
        "episodic_lid" => RETENTION_EPISODIC_GYR,
        "ice_shell_ocean" => {
            // A thin shell kept warm renews itself everywhere; a thick one keeps its
            // craters and renews only where it is thin.
            let ratio = c.f64("interior.ice_shell_thickness") / c.f64("body.radius").max(1.0);
            if c.f64("interior.tidal_flux") > ICE_RESURFACING_TIDAL_FLUX
                && ratio < WHOLESALE_SHELL_RATIO
            {
                RETENTION_ICE_OCEAN_GYR
            } else {
                age
            }
        }
        _ => {
            if p > SPACE_WEATHERING_MAX_PRESSURE {
                RETENTION_THIN_ATMOSPHERE_GYR
            } else {
                age
            }
        }
    };
    // Volcanic burial applies only where volcanism resurfaces the whole body, which the
    // tidal flux threshold of the Io regime marks; elsewhere flows are local and the
    // volcanism stage places them.
    if c.f64("interior.tidal_flux") > IO_REGIME_TIDAL_FLUX {
        let r = c.f64("body.radius");
        let rate_m_per_gyr = c.f64("interior.magma_production_rate") / (4.0 * PI * r * r) * GYR;
        if rate_m_per_gyr > 0.0 {
            retention = retention.min(RETENTION_BURIAL_DEPTH_M / rate_m_per_gyr);
        }
    }
    retention.clamp(1e-4, age).into()
}

pub fn surface_space_weathering_timescale(c: &Ctx) -> Value {
    if c.f64("atmosphere.surface_pressure") > SPACE_WEATHERING_MAX_PRESSURE {
        1.0e6
    } else {
        let d = c.f64("orbit.star_distance") / AU;
        (SPACE_WEATHERING_GYR_1AU * d * d).min(SPACE_WEATHERING_MAX_GYR)
    }
    .into()
}

// ------------------------------------------------------------------ materials

pub fn materials_bedrock(c: &Ctx) -> Value {
    let class = c.str("body.composition_class");
    let regime = c.str("interior.tectonic_regime");
    let airless = c.f64("atmosphere.surface_pressure") < SPACE_WEATHERING_MAX_PRESSURE;
    let bedrock = match class {
        "iron_rich" => "mercurian_regolith",
        "rock" => {
            if c.f64("interior.tidal_flux") > IO_REGIME_TIDAL_FLUX {
                "sulfur_coated"
            } else if regime == "mobile_lid" {
                "felsic"
            } else if airless && regime == "stagnant_lid" {
                "anorthosite"
            } else {
                "basalt"
            }
        }
        "rock_ice" => {
            if regime == "ice_shell_ocean" && airless {
                "water_ice"
            } else {
                "water_ice_dirty"
            }
        }
        _ => "water_ice",
    };
    bedrock.into()
}

pub fn materials_regolith_mix(c: &Ctx) -> Value {
    let mut m = BTreeMap::new();
    m.insert(c.str("materials.bedrock").to_string(), 1.0);
    Value::Map(m)
}

pub fn materials_lava(c: &Ctx) -> Value {
    let regime = c.str("interior.tectonic_regime");
    let airless = c.f64("atmosphere.surface_pressure") < SPACE_WEATHERING_MAX_PRESSURE;
    let lava = match c.str("body.composition_class") {
        "iron_rich" => "smooth_plains",
        "rock" => {
            if airless
                && regime == "stagnant_lid"
                && c.f64("interior.tidal_flux") <= IO_REGIME_TIDAL_FLUX
            {
                "mare_basalt"
            } else {
                "basalt"
            }
        }
        "rock_ice" => "water_ice_dirty",
        _ => "water_ice",
    };
    lava.into()
}

pub fn materials_deep_bedrock(c: &Ctx) -> Value {
    let regime = c.str("interior.tectonic_regime");
    let airless = c.f64("atmosphere.surface_pressure") < SPACE_WEATHERING_MAX_PRESSURE;
    let deep = match c.str("body.composition_class") {
        "iron_rich" => "low_reflectance_material",
        "rock" => {
            if airless && regime == "stagnant_lid" {
                "norite"
            } else {
                "basalt"
            }
        }
        // Under a shell that floats on an ocean the deep ice is clean; a crust over rock
        // is the same mixture all the way down.
        "rock_ice" | "ice_dominated" => {
            if regime == "ice_shell_ocean" {
                "water_ice"
            } else {
                c.str("materials.bedrock")
            }
        }
        _ => c.str("materials.bedrock"),
    };
    deep.into()
}

pub fn materials_deep_bedrock_depth(c: &Ctx) -> Value {
    match c.str("body.composition_class") {
        "iron_rich" => DEEP_BEDROCK_DEPTH_IRON_RICH,
        "rock" => DEEP_BEDROCK_DEPTH_ROCK,
        _ => DEEP_BEDROCK_DEPTH_ICE,
    }
    .into()
}

// ------------------------------------------------------------------ detail, output, seeds

fn quality_index(c: &Ctx) -> usize {
    match c.str("detail.quality") {
        "preview" => 0,
        "high" => 2,
        _ => 1,
    }
}

pub fn detail_loop_resolution(c: &Ctx) -> Value {
    (LOOP_FACE_SIZE[quality_index(c)] as f64).into()
}

pub fn detail_finishing_resolution(c: &Ctx) -> Value {
    (FINISHING_FACE_SIZE[quality_index(c)] as f64).into()
}

pub fn output_datum(c: &Ctx) -> Value {
    if c.str("volatiles.working_liquid") != "none" {
        "sea_level"
    } else {
        "mean_radius"
    }
    .into()
}

/// Sub-seed for a stage: a 53-bit value so it is exact in an f64.
pub fn stage_seed(master: f64, stage: &str) -> f64 {
    (crate::rng::stream_seed(master as u64, stage) >> 11) as f64
}

macro_rules! seed_fn {
    ($name:ident, $stage:literal) => {
        pub fn $name(c: &Ctx) -> Value {
            stage_seed(c.f64("seed.master"), $stage).into()
        }
    };
}

seed_fn!(seed_terrain, "terrain");
seed_fn!(seed_tectonics, "tectonics");
seed_fn!(seed_volcanism, "volcanism");
seed_fn!(seed_cratering, "cratering");
seed_fn!(seed_climate, "climate");
seed_fn!(seed_weather, "weather");
seed_fn!(seed_detail, "detail");

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sun_is_a_g2_star_at_present_luminosity() {
        let l_ref = mass_luminosity(1.0);
        assert!((l_ref - 1.0).abs() < 1e-3, "{l_ref}");
        let t_ms = main_sequence_lifetime_gyr(1.0, l_ref);
        assert!((t_ms - 10.0).abs() < 1e-2);
        assert!((luminosity_at_age(l_ref, 4.57, t_ms) - 1.0).abs() < 1e-3);
        assert!((luminosity_at_age(l_ref, 0.0, t_ms) - 1.0 / 1.4).abs() < 1e-3);
        assert_eq!(
            spectral_type(effective_temperature(1.0, stellar_radius(1.0))),
            "G2V"
        );
        assert_eq!(spectral_type(5772.0), "G2V");
        assert!((stellar_radius(1.0) - 1.0).abs() < 0.02);
    }

    #[test]
    fn earth_derived_quantities() {
        let m = 5.9722e24;
        let r = 6.371e6;
        assert!((gravity(m, r) - 9.82).abs() < 0.02);
        assert!((escape_velocity(m, r) - 11_186.0).abs() < 20.0);
        assert!((equilibrium_temperature(L_SUN, AU, 0.3) - 254.6).abs() < 1.0);
        assert!((kepler_period(AU, M_SUN + m) / 86_400.0 - 365.25).abs() < 0.5);
        let f = oblateness(7.2921e-5, r, m, 0.3307);
        assert!((f - 1.0 / 298.0).abs() < 2e-4, "flattening {f}");
    }

    #[test]
    fn composition_classes_match_the_reference_bodies() {
        let class = |rho: f64, mass_earth: f64| {
            let rho_u = rho / compression_factor(mass_earth * M_EARTH);
            let (core, _, ice) = composition_fractions(rho_u);
            composition_class(core, ice)
        };
        assert_eq!(class(5_427.0, 0.0553), "iron_rich");
        assert_eq!(class(5_514.0, 1.0), "rock");
        assert_eq!(class(3_934.0, 0.107), "rock");
        assert_eq!(class(3_344.0, 0.0123), "rock");
        assert_eq!(class(3_528.0, 0.015), "rock");
        assert_eq!(class(3_013.0, 0.008), "rock_ice");
        assert_eq!(class(1_936.0, 0.0248), "rock_ice");
        assert_eq!(class(1_880.0, 0.0225), "rock_ice");
        assert_eq!(class(1_609.0, 1.8e-5), "ice_dominated");
    }

    #[test]
    fn radiogenic_heat_is_present_value_now_and_higher_in_the_past() {
        assert!((radiogenic_per_kg(SOLAR_SYSTEM_AGE_GYR) - RADIOGENIC_NOW_W_PER_KG).abs() < 1e-18);
        let young = radiogenic_per_kg(1.0);
        assert!(
            young > 2.5 * RADIOGENIC_NOW_W_PER_KG && young < 3.5 * RADIOGENIC_NOW_W_PER_KG,
            "{young}"
        );
    }

    #[test]
    fn io_tidal_heating_is_within_a_factor_of_two() {
        let n = 2.0 * PI / (1.769 * 86_400.0);
        let p = tidal_power(K2_ROCK, Q_ROCK, M_JUP, 1.8216e6, n, 0.0041, 4.217e8);
        assert!(p > 0.5e14 && p < 2.0e14, "Io tidal power {p:e}");
    }

    #[test]
    fn surface_temperatures_of_the_reference_atmospheres() {
        // Venus: CO2 column from 92 bar.
        let mut venus = BTreeMap::new();
        venus.insert("CO2".to_string(), 1.0e6);
        let t = surface_temperature(226.6, tau_dry(&venus), 0.0);
        assert!((t - 737.0).abs() < 15.0, "Venus {t}");
        // Earth: trace CO2 plus the water feedback.
        let mut earth = BTreeMap::new();
        earth.insert("CO2".to_string(), 6.3);
        earth.insert("N2".to_string(), 7_800.0);
        let t = surface_temperature(254.6, tau_dry(&earth), 2_700.0);
        assert!((t - 288.0).abs() < 4.0, "Earth {t}");
        // Mars.
        let mut mars = BTreeMap::new();
        mars.insert("CO2".to_string(), 156.0);
        let t = surface_temperature(210.1, tau_dry(&mars), 25.0);
        assert!((t - 213.0).abs() < 4.0, "Mars {t}");
    }

    #[test]
    fn vapour_pressure_curves_pass_through_known_points() {
        let w = volatile("H2O").unwrap();
        assert!((vapour_pressure(w, 373.15) - 101_325.0).abs() / 101_325.0 < 0.08);
        assert!((saturation_temperature(w, 101_325.0) - 373.15).abs() < 3.0);
        let n2 = volatile("N2").unwrap();
        assert!((saturation_temperature(n2, 101_325.0) - 77.4).abs() < 1.5);
        let ch4 = volatile("CH4").unwrap();
        assert!((saturation_temperature(ch4, 101_325.0) - 111.7).abs() < 3.0);
        assert!((saturation_temperature(ch4, 1.47e5) - 117.0).abs() < 4.0);
    }

    #[test]
    fn phase_states_of_titan_and_mars() {
        let ch4 = volatile("CH4").unwrap();
        assert_eq!(
            phase_state(ch4, 94.0, 1.47e5, true),
            PHASE_LIQUID_AND_SOLID
                .max(PHASE_LIQUID)
                .min(phase_state(ch4, 94.0, 1.47e5, true))
        );
        assert!(matches!(phase_state(ch4, 94.0, 1.47e5, true) as i32, 2 | 4));
        let n2 = volatile("N2").unwrap();
        assert_eq!(phase_state(n2, 94.0, 1.47e5, true), PHASE_GAS);
        let w = volatile("H2O").unwrap();
        assert_eq!(phase_state(w, 213.0, 610.0, true), PHASE_SOLID);
        assert_eq!(phase_state(w, 288.0, 101_325.0, true), PHASE_LIQUID);
        assert_eq!(phase_state(n2, 44.0, 1.0, true), PHASE_SOLID);
    }

    #[test]
    fn held_hou_and_transport_for_earth() {
        let edge = held_hou_edge_deg(9.81, 2.0 * 8_500.0, 0.25, 7.2921e-5, 6.371e6);
        assert!(edge > 28.0 && edge < 36.0, "{edge}");
        let d = transport_coefficient(101_325.0, 1_004.0, 0.028_97, 7.2921e-5);
        assert!((d - 0.58).abs() < 0.01);
    }

    #[test]
    fn crater_cutoff_anchors() {
        assert!((crater_cutoff(101_325.0, 9.81) - 150.0).abs() < 1.0);
        let venus = crater_cutoff(9.2e6, 8.87);
        assert!(venus > 1_500.0 && venus < 3_000.0, "{venus}");
    }
}
