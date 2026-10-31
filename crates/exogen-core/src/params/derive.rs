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

    /// F8, F122: (core, rock, ice) mass fractions from the uncompressed density. A body
    /// formed where water ice condensed (`wet` 1) carries iron and rock in the chondritic
    /// ratio and ice for the rest of a deficit below their mix; one formed too warm for
    /// ice (`wet` 0) carries none, and a deficit below the chondritic mix is iron it lacks
    /// (the Moon's, Mars's); between them a blend. Above the chondritic mix both are the
    /// same: iron and rock.
    pub fn composition_fractions(rho_u: f64, wet: f64) -> (f64, f64, f64) {
        let core_share = CHONDRITIC_CORE_TO_ROCK / (1.0 + CHONDRITIC_CORE_TO_ROCK);
        let rho_rock_iron = 1.0 / (core_share / DENSITY_IRON + (1.0 - core_share) / DENSITY_ROCK);
        let f = ((1.0 / DENSITY_ROCK - 1.0 / rho_u) / (1.0 / DENSITY_ROCK - 1.0 / DENSITY_IRON))
            .clamp(0.0, 1.0);
        let dry = (f, 1.0 - f, 0.0);
        if rho_u >= rho_rock_iron {
            return dry;
        }
        let fi = ((1.0 / rho_u - 1.0 / rho_rock_iron) / (1.0 / DENSITY_ICE - 1.0 / rho_rock_iron))
            .clamp(0.0, 1.0);
        let core = (1.0 - fi) * core_share;
        let icy = (core, 1.0 - fi - core, fi);
        let w = wet.clamp(0.0, 1.0);
        (
            w * icy.0 + (1.0 - w) * dry.0,
            w * icy.1 + (1.0 - w) * dry.1,
            w * icy.2 + (1.0 - w) * dry.2,
        )
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

    /// F68: the mean porosity of a body's rocky interior under its H2O layer, by volume:
    /// each depth's porosity from the pressure there, the H2O layer's weight and the
    /// interior's own above it, as sediments compact (Athy's law).
    pub fn core_porosity(radius: f64, mass: f64, ice_fraction: f64, ice_layer: f64) -> f64 {
        let rc = (radius - ice_layer).max(0.0);
        if rc <= 0.0 || mass <= 0.0 {
            return 0.0;
        }
        let core_mass = mass * (1.0 - ice_fraction).clamp(0.0, 1.0);
        let rho = core_mass / (4.0 / 3.0 * PI * rc.powi(3));
        // The H2O layer's weight at the interior's surface, under the mean of the gravity
        // there and at the surface.
        let g_top = G * core_mass / (rc * rc);
        let g_surface = G * mass / (radius * radius);
        let top = DENSITY_ICE * 0.5 * (g_top + g_surface) * ice_layer;
        const SHELLS: usize = 64;
        let (mut sum, mut volume) = (0.0, 0.0);
        for i in 0..SHELLS {
            let r = (i as f64 + 0.5) / SHELLS as f64 * rc;
            let pressure = top + 2.0 / 3.0 * PI * G * rho * rho * (rc * rc - r * r);
            let dv = r * r;
            sum += POROSITY_GRANULAR * (-pressure / POROSITY_PRESSURE).exp() * dv;
            volume += dv;
        }
        sum / volume
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

    /// F30: simple-to-complex transition diameter, metres, on a target `ice` of the way
    /// from rock to ice (audit S9).
    pub fn transition_diameter(g: f64, ice: f64) -> f64 {
        TRANSITION_LUNAR_KM
            * 1e3
            * (GRAVITY_MOON / g)
            * crate::crater::by_ice(1.0, TRANSITION_ICE_FACTOR, ice)
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

/// F234: the year, the orbit about the star: the body's own, or its parent's for a moon.
pub fn orbit_year(c: &Ctx) -> Value {
    match c.str("orbit.primary") {
        "planet" => kepler_period(
            c.f64("parent.semi_major_axis"),
            c.f64("star.mass") * M_SUN + c.f64("parent.mass"),
        ),
        _ => c.f64("orbit.period"),
    }
    .into()
}

/// F234: the eccentricity of the orbit about the star: the body's own, or its parent's.
pub fn orbit_star_eccentricity(c: &Ctx) -> Value {
    match c.str("orbit.primary") {
        "planet" => c.f64("parent.eccentricity"),
        _ => c.f64("orbit.eccentricity"),
    }
    .into()
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

/// How far a body formed where water ice could condense (F122): 1 at or below the snow
/// line's temperature, 0 from SNOW_LINE_DRY_K, by its blackbody temperature at its distance
/// from its star (zero albedo, its whole surface sharing the heat).
fn formed_wet(c: &Ctx) -> f64 {
    ((SNOW_LINE_DRY_K - formation_temperature(c)) / (SNOW_LINE_DRY_K - SNOW_LINE_WET_K))
        .clamp(0.0, 1.0)
}

/// A body's blackbody temperature at its distance from its star, K, standing for where it
/// formed (F122, F125).
fn formation_temperature(c: &Ctx) -> f64 {
    (c.f64("orbit.insolation") / (4.0 * SIGMA_SB))
        .max(0.0)
        .powf(0.25)
}

/// F125: the mantle's FeO, mass fraction, from the warmth where the body formed.
pub fn body_mantle_feo(c: &Ctx) -> Value {
    let reduced = ((FEO_REDUCED_K - formation_temperature(c)) / (FEO_REDUCED_K - FEO_OXIDIZED_K))
        .clamp(0.0, 1.0);
    (FEO_REDUCED + (FEO_OXIDIZED - FEO_REDUCED) * reduced).into()
}

/// F125: what a rocky body's magma ocean floated, as shares of its primary crust: plagioclase
/// where its last liquid lay within the plagioclase field and was iron-rich enough to float
/// it, graphite where the ocean was too iron-poor to oxidize its carbon.
fn flotation(c: &Ctx) -> (f64, f64) {
    let r = c.f64("body.radius");
    let core = c.f64("body.mass") * c.f64("body.core_mass_fraction");
    let core_radius = (3.0 * core / (4.0 * PI * DENSITY_IRON)).cbrt();
    let field = PLAGIOCLASE_MAX_PRESSURE / (DENSITY_ROCK * c.f64("body.gravity")).max(1e-9);
    let share = ((r.powi(3) - (r - field.min(r)).powi(3))
        / (r.powi(3) - core_radius.powi(3)).max(1.0))
    .clamp(0.0, 1.0);
    let room = ((share / PLAGIOCLASE_RESIDUAL - 0.5) / 0.5).clamp(0.0, 1.0);
    let feo = c.f64("body.mantle_feo");
    let iron = ((feo - FLOTATION_IRON_FEO.0) / (FLOTATION_IRON_FEO.1 - FLOTATION_IRON_FEO.0))
        .clamp(0.0, 1.0);
    (room * iron, 1.0 - iron)
}

/// Whether a body's surface is rock under the materials' process rules (F125): a rocky or
/// iron-rich body neither renewed by a tidally heated body's plumes nor recycled by plates.
fn rocky_crust(c: &Ctx) -> bool {
    matches!(c.str("body.composition_class"), "rock" | "iron_rich")
        && c.f64("interior.tidal_flux") <= IO_REGIME_TIDAL_FLUX
        && c.str("interior.tectonic_regime") != "mobile_lid"
}

pub fn body_core_mass_fraction(c: &Ctx) -> Value {
    composition_fractions(c.f64("body.uncompressed_density"), formed_wet(c))
        .0
        .into()
}

pub fn body_ice_mass_fraction(c: &Ctx) -> Value {
    composition_fractions(c.f64("body.uncompressed_density"), formed_wet(c))
        .2
        .into()
}

pub fn body_rock_mass_fraction(c: &Ctx) -> Value {
    (1.0 - c.f64("body.core_mass_fraction") - c.f64("body.ice_mass_fraction"))
        .max(0.0)
        .into()
}

/// F126: the share of an icy crust's volume that is not ice, as it accreted: the body's
/// rock and metal against its ice. With grains of one size it is the share of their
/// cross-sections, which an intimate mixture's reflectance weighs (Hapke 1981).
pub fn body_crust_impurity(c: &Ctx) -> Value {
    let rock = c.f64("body.core_mass_fraction") / DENSITY_IRON
        + c.f64("body.rock_mass_fraction") / DENSITY_ROCK;
    let ice = c.f64("body.ice_mass_fraction") / DENSITY_ICE;
    if rock + ice <= 0.0 {
        return 1.0.into();
    }
    (rock / (rock + ice)).into()
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

/// F234: tides damp the tilt of a body locked to its star; a moon locked to its planet keeps
/// its tilt to its planet's orbit about the star.
pub fn rotation_effective_obliquity(c: &Ctx) -> Value {
    if c.str("rotation.mode") == "locked" && c.str("orbit.primary") != "planet" {
        0.0
    } else {
        c.f64("rotation.obliquity")
    }
    .into()
}

/// F234: sidereal turns per orbit about the star.
pub fn rotation_spins_per_year(c: &Ctx) -> Value {
    (c.f64("orbit.year") / c.f64("rotation.effective_period")).into()
}

/// F234: the spins per orbit of an exact low-order resonance with the star, 0 when there is
/// none; a moon's spin is held by its planet, not by the star.
pub fn rotation_resonance(c: &Ctx) -> Value {
    if c.str("orbit.primary") == "planet" {
        return 0.0.into();
    }
    crate::insolation::spin_orbit_resonance(c.f64("rotation.spins_per_year"))
        .map_or(0.0, |(p, q)| p as f64 / q as f64)
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

/// F75: an ice shell's h₂ from its stiffness against the ocean's buoyancy (Beuthe 2018),
/// or the homogeneous elastic sphere's for ice frozen to the rock; 1 + k₂ for a
/// fluid-like rocky interior.
pub fn interior_love_number_h2(c: &Ctx) -> Value {
    if !is_icy(c.str("body.composition_class")) {
        return (1.0 + c.f64("interior.love_number_k2")).into();
    }
    let love = crate::stages::tidal::ShellLove::new(
        c.f64("body.radius"),
        c.f64("body.gravity"),
        c.f64("body.density"),
    );
    if c.bool("interior.has_subsurface_ocean") {
        love.h2(c.f64("interior.ice_shell_thickness"))
    } else {
        love.grounded
    }
    .into()
}

/// The tides an icy moon's shell dissipates (F75), from its orbit and its interior, for a
/// locked moon of a planet; none otherwise.
pub fn shell_tides(c: &Ctx) -> Option<crate::stages::tidal::ShellTides> {
    if !is_icy(c.str("body.composition_class"))
        || c.str("orbit.primary") != "planet"
        || c.str("rotation.mode") != "locked"
    {
        return None;
    }
    let radius = c.f64("body.radius");
    Some(crate::stages::tidal::ShellTides::new(
        crate::stages::tidal::ShellLove::new(radius, c.f64("body.gravity"), c.f64("body.density")),
        c.f64("orbit.primary_mass"),
        c.f64("body.mass"),
        radius,
        c.f64("orbit.semi_major_axis"),
        c.f64("orbit.mean_motion"),
        c.f64("atmosphere.surface_temperature_mean"),
        ice_melting_point(c),
        c.f64("body.ice_layer_thickness"),
    ))
}

/// F138: an icy moon's shell and porous interior through its life under its eccentricity
/// history.
fn shell_history(c: &Ctx) -> Option<(crate::orbital::OrbitalInputs, crate::orbital::ShellHistory)> {
    let o = crate::orbital::OrbitalInputs::from_ctx(c)?;
    let h = o.shell_history(c.curve("orbit.eccentricity_history"));
    Some((o, h))
}

/// F75 and F138: the tides an icy moon's shell dissipates in itself at the thickness its
/// history leaves it, W/m².
pub fn interior_shell_tidal_flux(c: &Ctx) -> Value {
    shell_history(c)
        .map_or(0.0, |(o, h)| {
            let dry = crate::orbital::ShellState {
                wet: false,
                ..h.now()
            };
            o.tidal_flux(&dry, c.f64("orbit.eccentricity"))
        })
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
    if c.str("interior.tectonic_regime") != "ice_shell_ocean" {
        return 0.0.into();
    }
    // The shell cracks where it is stressed past the cracking stress: a hard shell most
    // where the focusing thins it (F76).
    let love = crate::stages::tidal::ShellLove::new(
        c.f64("body.radius"),
        c.f64("body.gravity"),
        c.f64("body.density"),
    );
    let thinnest = 1.0 - c.f64("tectonics.polar_thinning");
    let peak = stress * love.local_amplification(c.f64("interior.ice_shell_thickness"), thinnest);
    if peak <= ICE_CRACK_STRESS {
        return 0.0.into();
    }
    (CRACK_DENSITY_FACTOR * (stress / ICE_TENSILE_STRENGTH)
        / c.f64("interior.ice_shell_thickness").max(1.0))
    .into()
}

/// F138 and F140: how long the shell has floated on an ocean up to the present, Gyr: an icy
/// moon's from its history, any other body's its whole age when it has an ocean now.
pub fn interior_ocean_age(c: &Ctx) -> Value {
    if let Some((o, h)) = shell_history(c) {
        return h.ocean_age(o.layer()).into();
    }
    if c.bool("interior.has_subsurface_ocean") {
        c.f64("body.age").into()
    } else {
        0.0.into()
    }
}

pub fn tectonics_polar_thinning(c: &Ctx) -> Value {
    if c.str("interior.tectonic_regime") != "ice_shell_ocean" {
        return 0.0.into();
    }
    let d = c.f64("interior.ice_shell_thickness");
    let ratio = d / c.f64("body.radius").max(1.0);
    let focus = 1.0 - (-ratio / POLAR_FOCUS_SHELL_RATIO).exp();
    // F140: the focusing thins the shell where its tides meet its conduction, in a steady
    // state; a shell that conducts more than is made in it is freezing, and thins by the
    // share made. The thinning grows over the shell's conductive time from its ocean's
    // forming.
    let made = c.f64("interior.radiogenic_flux")
        * (1.0 + c.f64("interior.secular_cooling_fraction"))
        + c.f64("interior.tidal_flux");
    let balance = (made / c.f64("interior.surface_heat_flux").max(1e-9)).min(1.0);
    let grown = 1.0
        - (-c.f64("interior.ocean_age") * SECONDS_PER_GYR * ICE_THERMAL_DIFFUSIVITY
            / (d * d).max(1.0))
        .exp();
    (POLAR_THINNING_MAX
        * (c.f64("interior.tidal_flux") / POLAR_FOCUS_FLUX).min(1.0)
        * focus
        * balance
        * grown)
        .into()
}

/// Audit S5: how far the cracks the tides open from above at the thinnest shell reach
/// toward the water a crack from the ocean holds there, which stands below the surface by
/// the share of the shell the ice's buoyancy leaves dry (CRACK_REACH_WIDTH): ocean water
/// reaches the surface past one.
pub fn tectonics_crack_reach(c: &Ctx) -> Value {
    if c.str("interior.tectonic_regime") != "ice_shell_ocean" {
        return 0.0.into();
    }
    let d = c.f64("interior.ice_shell_thickness");
    let thinnest = 1.0 - c.f64("tectonics.polar_thinning");
    let love = crate::stages::tidal::ShellLove::new(
        c.f64("body.radius"),
        c.f64("body.gravity"),
        c.f64("body.density"),
    );
    let stress = c.f64("tectonics.tidal_stress") * love.local_amplification(d, thinnest);
    let opened = stress / (DENSITY_CRUST_ICE * c.f64("body.gravity").max(1e-3));
    let standing = d * thinnest * (1.0 - DENSITY_CRUST_ICE / DENSITY_OCEAN);
    (opened / standing.max(1.0)).into()
}

/// F53: the focused heating returns every so many conductive times of the shell, where
/// it thins the shell enough for the ocean to reach the surface (audit S5).
pub fn tectonics_focusing_period(c: &Ctx) -> Value {
    if c.f64("tectonics.crack_reach") < 1.0 {
        return 0.0.into();
    }
    let d = c.f64("interior.ice_shell_thickness");
    (FOCUSING_PERIOD_FACTOR * d * d / ICE_THERMAL_DIFFUSIVITY / SECONDS_PER_GYR).into()
}

pub fn tectonics_focusing_episode_count(c: &Ctx) -> Value {
    let period = c.f64("tectonics.focusing_period");
    if period <= 0.0 {
        return 0.0.into();
    }
    ((c.f64("body.age") - FOCUSING_EPISODE_DURATION_GYR) / period)
        .floor()
        .max(0.0)
        .into()
}

pub fn tectonics_ice_resurfaced_fraction(c: &Ctx) -> Value {
    if c.str("interior.tectonic_regime") != "ice_shell_ocean" {
        return 0.0.into();
    }
    // A thin shell renews itself by chaos and cracking instead, by its degree (audit S4).
    let thick = 1.0 - thin_shell(c);
    if thick <= 0.0 {
        return 0.0.into();
    }
    // What the focusing takes to the pole does not spread as lanes.
    let diffuse = (1.0 - c.f64("tectonics.polar_thinning") / POLAR_THINNING_MAX).max(0.0);
    // The tidal heat its lid carried above the floor over its life (F140): the pulse of a
    // resonance episode (F50), or the steady flux, for as long as it lasted, against the
    // heat F44 was calibrated on, Ganymede's two thirds renewed over its episode. A shell
    // renews itself at most once over: the renewal share saturates at one, and the
    // focusing's share is taken from that.
    // A lane is a rift, which the lid breaks in tension: while the body's change of radius
    // has the lid yielding in compression, its thrusts taking up the shortening (F147), as
    // a melting ocean's shrinking does, the heat opens no lane.
    let excess = renewing_excess(c);
    let lid = lid_strain(c);
    let squeezed = |t: f64| -> bool {
        let i = lid.times.partition_point(|&x| x < t);
        i > 0 && i < lid.times.len() && lid.shortening[i] > lid.shortening[i - 1]
    };
    let carried: f64 = excess
        .windows(2)
        .filter(|w| !squeezed(0.5 * (w[0].0 + w[1].0)))
        .map(|w| 0.5 * (w[0].1 + w[1].1) * (w[1].0 - w[0].0))
        .sum();
    let renewal = (carried / (ICE_RESURFACING_FLUX * LANE_RENEWAL_GYR)).clamp(0.0, 1.0);
    (renewal * diffuse * thick).clamp(0.0, 0.95).into()
}

/// How far a shell over an ocean is thin against its body (audit S4): wholly below
/// 1 - THIN_SHELL_WIDTH of WHOLESALE_SHELL_RATIO of the radius, not at all above
/// 1 + THIN_SHELL_WIDTH of it, smoothly between.
fn thin_shell(c: &Ctx) -> f64 {
    if c.str("interior.tectonic_regime") != "ice_shell_ocean" {
        return 0.0;
    }
    let ratio = c.f64("interior.ice_shell_thickness") / c.f64("body.radius").max(1.0);
    let x = ((1.0 + THIN_SHELL_WIDTH) * WHOLESALE_SHELL_RATIO - ratio)
        / (2.0 * THIN_SHELL_WIDTH * WHOLESALE_SHELL_RATIO);
    let x = x.clamp(0.0, 1.0);
    x * x * (3.0 - 2.0 * x)
}

/// How far a shell renews its whole surface (audit S6): as far as it is thin, at a rate in
/// proportion to its tides up to ICE_RESURFACING_TIDAL_FLUX.
fn wholesale_renewal(c: &Ctx) -> f64 {
    thin_shell(c) * (c.f64("interior.tidal_flux") / ICE_RESURFACING_TIDAL_FLUX).clamp(0.0, 1.0)
}

/// The lid's strain from the body's change of radius (F147), from the nodes a derivation
/// declares.
fn lid_strain(c: &Ctx) -> crate::stages::strain::LidStrain {
    let age = c.f64("body.age");
    let inputs = crate::stages::strain::LidInputs {
        radius: c.f64("body.radius"),
        ice: c.f64("body.ice_layer_thickness"),
        age,
        start: c
            .f64("cratering.record_start")
            .max(age - c.f64("cratering.retention_age"))
            .clamp(0.0, age),
        gravity: c.f64("body.gravity"),
    };
    let thermal = crate::thermal::ThermalHistory::new(&thermal_inputs(c));
    crate::stages::strain::LidStrain::from_inputs(&inputs, &thermal)
}

pub fn tectonics_chaos_fraction(c: &Ctx) -> Value {
    if c.str("interior.tectonic_regime") != "ice_shell_ocean" {
        return 0.0.into();
    }
    (CHAOS_FRACTION_MAX * (c.f64("interior.tidal_flux") / CHAOS_FLUX).min(1.0) * thin_shell(c))
        .into()
}

pub fn volcanism_cryodome_count(c: &Ctx) -> Value {
    // A shell over rock builds them, as far as no ocean floats it (audit S1).
    let open = 1.0 - c.f64("interior.ocean_degree");
    if !c.str("interior.tectonic_regime").starts_with("ice_shell") || open <= 0.0 {
        return 0.0.into();
    }
    let r_km3 = (c.f64("body.radius") / 1e3).powi(3);
    let (lo, hi) = (
        CRYODOME_VOLUME_SHARE.0 * r_km3,
        CRYODOME_VOLUME_SHARE.1 * r_km3,
    );
    let mean = (hi - lo) / (hi / lo).ln();
    // The magma rate is kept in m³/s; the volumes are km³.
    (CRYODOME_SHARE * c.f64("volcanism.rate") * GYR / 1e9 * c.f64("body.age") / mean * open)
        .round()
        .into()
}

pub fn volcanism_plume_fallout_rate(c: &Ctx) -> Value {
    // The plumes erupt where the ocean reaches the surface, the more the farther the
    // cracks reach past it (audit S5).
    let reach = c.f64("tectonics.crack_reach");
    if reach < 1.0 || c.f64("tectonics.ice_crack_density") <= 0.0 {
        return 0.0.into();
    }
    (PLUME_FALLOUT_RATE
        * (c.f64("interior.tidal_flux") / POLAR_FOCUS_FLUX).min(1.0)
        * ((reach - 1.0) / CRACK_REACH_WIDTH).min(1.0))
    .into()
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
/// F44 and F140: the tidal heat an icy shell's lid carries above the floor that renews its
/// surface, at the thermal history's samples: the lesser of its tides and the heat its lid
/// conducts, so that a shell melting from below under a lid still thick, or freezing from
/// a warm start, renews nothing.
fn renewing_excess(c: &Ctx) -> Vec<(f64, f64)> {
    let heat = c.curve("interior.heat_flux_history");
    c.curve("interior.tidal_flux_history")
        .points
        .iter()
        .map(|&(t, q)| {
            (
                t,
                (q.min(heat.eval(t)) - ICE_RESURFACING_FLUX_FLOOR).max(0.0),
            )
        })
        .collect()
}

fn resurfacing_episode(c: &Ctx) -> Option<(f64, f64)> {
    let (mut w, mut wt, mut wtt) = (0.0, 0.0, 0.0);
    for (t, x) in renewing_excess(c) {
        w += x;
        wt += x * t;
        wtt += x * t * t;
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

/// F13 for a rocky body; for an icy one the tides its shell dissipates (F75), and for
/// either those of a porous interior (F68).
pub fn interior_tidal_power(c: &Ctx) -> Value {
    let r = c.f64("body.radius");
    let own = if is_icy(c.str("body.composition_class")) {
        c.f64("interior.shell_tidal_flux") * 4.0 * PI * r * r
    } else {
        tidal_power(
            c.f64("interior.love_number_k2"),
            c.f64("interior.dissipation_q"),
            c.f64("orbit.primary_mass"),
            r,
            c.f64("orbit.mean_motion"),
            c.f64("orbit.eccentricity"),
            c.f64("orbit.semi_major_axis"),
        )
    };
    (own + c.f64("interior.core_tidal_power")).into()
}

pub fn interior_core_porosity(c: &Ctx) -> Value {
    core_porosity(
        c.f64("body.radius"),
        c.f64("body.mass"),
        c.f64("body.ice_mass_fraction"),
        c.f64("body.ice_layer_thickness"),
    )
    .into()
}

/// F68: the tides a porous rocky interior dissipates: the eccentricity tides of F13 on the
/// interior's radius, with an imaginary Love number from its porosity above that of a
/// consolidated rock.
pub fn interior_core_tidal_power(c: &Ctx) -> Value {
    // F138: an icy moon's porous interior dissipates only while its pores hold water.
    if let Some((o, h)) = shell_history(c) {
        let e = c.f64("orbit.eccentricity");
        return if h.now().wet {
            o.core_per_e2 * e * e
        } else {
            0.0
        }
        .into();
    }
    let excess = (c.f64("interior.core_porosity") - POROSITY_CONSOLIDATED).max(0.0);
    if excess <= 0.0 {
        return 0.0.into();
    }
    tidal_power(
        POROUS_DISSIPATION * excess,
        1.0,
        c.f64("orbit.primary_mass"),
        (c.f64("body.radius") - c.f64("body.ice_layer_thickness")).max(0.0),
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

/// The heat flux of a mantle in equilibrium with its heat: radiogenic plus secular plus
/// tidal.
fn equilibrium_heat_flux(c: &Ctx) -> f64 {
    c.f64("interior.radiogenic_flux") * (1.0 + c.f64("interior.secular_cooling_fraction"))
        + c.f64("interior.tidal_flux")
}

/// The mantle's inputs from the graph context (F121).
fn mantle_inputs(c: &Ctx) -> crate::mantle::MantleInputs {
    crate::mantle::MantleInputs {
        age_gyr: c.f64("body.age"),
        radius: c.f64("body.radius"),
        gravity: c.f64("body.gravity"),
        mass: c.f64("body.mass"),
        rock_fraction: c.f64("body.rock_mass_fraction"),
        core_fraction: c.f64("body.core_mass_fraction"),
        surface_temperature: c.f64("atmosphere.surface_temperature_mean"),
        tidal_flux_now: c.f64("interior.tidal_flux"),
        eccentricity_now: c.f64("orbit.eccentricity"),
        eccentricity_history: c.curve("orbit.eccentricity_history").clone(),
        crust_thickness: c.f64("tectonics.initial_crust_thickness"),
    }
}

/// The history of a rocky mantle whose lid may break (F121).
pub fn mantle_history(c: &Ctx) -> crate::mantle::MantleHistory {
    crate::mantle::evolve(&mantle_inputs(c))
}

/// Whether the plates' rule makes a mobile lid: surface water, enough heat, a large body.
fn plates(c: &Ctx) -> bool {
    c.str("volatiles.working_liquid") == "water"
        && equilibrium_heat_flux(c) >= MOBILE_LID_MIN_FLUX
        && c.f64("body.radius") >= MOBILE_LID_MIN_RADIUS
}

/// Whether a body's mantle keeps its own heat's history (F121, F122): a rocky body with no
/// plates, not heated past the Io regime.
fn mantle_keeps_heat(c: &Ctx) -> bool {
    !is_icy(c.str("body.composition_class"))
        && !plates(c)
        && c.f64("interior.tidal_flux") <= IO_REGIME_TIDAL_FLUX
}

/// Whether a dry rocky lid has broken (F121): its mantle's convective stress has passed the
/// lid's friction at least once, however long ago (audit S1).
fn lid_has_broken(c: &Ctx) -> bool {
    mantle_keeps_heat(c) && mantle_history(c).last_overturn().is_some()
}

/// Audit S1: how far a lid that has broken is episodic now, 1 while it last overturned
/// within EPISODIC_RECENT_SHARE's first share of the body's age, none past its second or
/// on a lid that never broke, by degrees between.
pub fn interior_episodic_degree(c: &Ctx) -> Value {
    if c.str("interior.tectonic_regime") != "episodic_lid" {
        return 0.0.into();
    }
    let since = c.f64("tectonics.resurfacing_event_age") / c.f64("body.age").max(1e-9);
    let (a, b) = EPISODIC_RECENT_SHARE;
    let x = ((b - since) / (b - a)).clamp(0.0, 1.0);
    (x * x * (3.0 - 2.0 * x)).into()
}

pub fn interior_surface_heat_flux(c: &Ctx) -> Value {
    // A rocky mantle that keeps its heat: the flux its history gives now (F121, F122).
    if mantle_keeps_heat(c) {
        return mantle_history(c)
            .heat_flux
            .last()
            .copied()
            .unwrap_or_else(|| equilibrium_heat_flux(c))
            .into();
    }
    // An icy moon's shell conducts what its history leaves it to (F138): the heat made
    // below it in a steady state, more while it freezes and less while it melts.
    if let Some((o, h)) = shell_history(c) {
        return o
            .surface_flux(&h.now(), c.f64("orbit.eccentricity"), c.f64("body.age"))
            .into();
    }
    equilibrium_heat_flux(c).into()
}

/// The thermal evolution's inputs from the graph context, mirroring
/// `ThermalInputs::from_graph`.
fn thermal_inputs(c: &Ctx) -> crate::thermal::ThermalInputs {
    let regime = c.str("interior.tectonic_regime");
    let mantle = mantle_keeps_heat(c).then(|| mantle_history(c));
    crate::thermal::ThermalInputs {
        age_gyr: c.f64("body.age"),
        radius: c.f64("body.radius"),
        gravity: c.f64("body.gravity"),
        uncompressed_density: c.f64("body.uncompressed_density"),
        rock_mass: c.f64("body.mass") * c.f64("body.rock_mass_fraction"),
        secular_fraction: c.f64("interior.secular_cooling_fraction"),
        mantle: mantle.clone(),
        plumes: c.f64("tectonics.plume_rise_count") >= 1.0,
        tidal_flux_now: c.f64("interior.tidal_flux"),
        eccentricity_now: c.f64("orbit.eccentricity"),
        eccentricity_history: c.curve("orbit.eccentricity_history").clone(),
        orbital: shell_history(c),
        icy: is_icy(c.str("body.composition_class")),
        surface_temperature: c.f64("atmosphere.surface_temperature_mean"),
        melting_point: ice_melting_point(c),
        crust_strength: c.f64("interior.crust_strength"),
        magma_efficiency_base: crate::thermal::magma_efficiency(
            regime,
            c.f64("interior.episodic_degree"),
        ),
    }
}

/// F138: an icy moon locked to its planet evolves its eccentricity with its interior
/// (`orbital`). For any other body, F50 and F77: the eccentricity about the mean the
/// present value stands at, with the pulse of a resonance episode centred
/// `orbit.resonance_age` before the present and `orbit.resonance_duration` wide at its
/// base, rising over the first fifth, holding at the forced value, falling over the last
/// fifth, and a resonance's cycle from its onset; constant when there is neither.
pub fn orbit_eccentricity_history(c: &Ctx) -> Value {
    if let Some(o) = crate::orbital::OrbitalInputs::from_ctx(c) {
        let forcing = crate::orbital::Forcing::from_ctx(c, &o);
        return o.eccentricity_history(&forcing).into();
    }
    let e_now = c.f64("orbit.eccentricity");
    let age = c.f64("body.age");
    let cycle = crate::orbital::ResonanceCycle::new(
        c.f64("orbit.cycle_period"),
        c.f64("orbit.cycle_onset"),
        c.f64("orbit.cycle_amplitude"),
        c.f64("orbit.cycle_phase"),
        age,
    );
    // The mean that the cycle swings about, so that the present takes its own value.
    let e = cycle
        .as_ref()
        .map_or(e_now, |k| e_now / k.factor(age).max(1e-6));
    let peak = c.f64("orbit.resonance_eccentricity");
    let before = c.f64("orbit.resonance_age");
    let width = c.f64("orbit.resonance_duration");
    let pulsed = peak > e && width > 0.0 && before > 0.0 && before <= age;
    if !pulsed && cycle.is_none() {
        return Curve::constant(e_now).into();
    }
    let t_peak = age - before;
    let half = 0.5 * width;
    let pulse = |t: f64| {
        if !pulsed {
            return e;
        }
        let x = (t - t_peak).abs() / half;
        e + (peak - e) * ((1.0 - x) / RESONANCE_RAMP_SHARE).clamp(0.0, 1.0)
    };
    let mut times = vec![0.0, age];
    if pulsed {
        let hold = half * (1.0 - RESONANCE_RAMP_SHARE);
        times.extend([t_peak - half, t_peak - hold, t_peak + hold, t_peak + half]);
    }
    if let Some(k) = &cycle {
        let step = k.period / CYCLE_SAMPLES_PER_PERIOD as f64;
        let count = ((age - k.onset) / step).ceil() as usize;
        times.extend((0..=count).map(|i| k.onset + i as f64 * step));
    }
    times.retain(|t| (0.0..=age).contains(t));
    times.sort_by(|a, b| a.total_cmp(b));
    times.dedup_by(|a, b| (*a - *b).abs() < 1e-9);
    let at = |t: f64| pulse(t) * cycle.as_ref().map_or(1.0, |k| k.factor(t));
    Curve::new(times.into_iter().map(|t| (t, at(t))).collect()).into()
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

/// Audit S1: how far an ice shell floats on its ocean, by the ocean's share of the ice
/// layer: none on a shell frozen to its rock, wholly from OCEAN_SHARE_WIDTH up, smoothly
/// between.
pub fn interior_ocean_degree(c: &Ctx) -> Value {
    if !c.bool("interior.has_subsurface_ocean") {
        return 0.0.into();
    }
    let layer = c.f64("body.ice_layer_thickness");
    let share = (layer - c.f64("interior.ice_shell_thickness")) / layer.max(1.0);
    let x = (share / OCEAN_SHARE_WIDTH).clamp(0.0, 1.0);
    (x * x * (3.0 - 2.0 * x)).into()
}

pub fn interior_tectonic_regime(c: &Ctx) -> Value {
    let regime = if is_icy(c.str("body.composition_class")) {
        if c.bool("interior.has_subsurface_ocean") {
            "ice_shell_ocean"
        } else {
            "ice_shell_rock"
        }
    } else if plates(c) {
        "mobile_lid"
    } else if lid_has_broken(c) {
        // A dry rocky lid that its mantle has broken (F121, audit S1).
        "episodic_lid"
    } else {
        "stagnant_lid"
    };
    regime.into()
}

pub fn interior_magma_production_rate(c: &Ctx) -> Value {
    let r = c.f64("body.radius");
    let regime = c.str("interior.tectonic_regime");
    let icy = is_icy(c.str("body.composition_class"));
    let base = crate::thermal::magma_efficiency(regime, c.f64("interior.episodic_degree"));
    // The tidal boost is a silicate mantle's heat pipe (Io), its melt lighter than its
    // lid; an ice mantle's melt is denser than the ice over it (F133).
    let tidal = if icy {
        0.0
    } else {
        (MAGMA_EFFICIENCY_TIDAL_MAX
            * (c.f64("interior.tidal_flux") / MAGMA_TIDAL_FLUX_SCALE).min(1.0))
        .max(0.0)
    };
    let efficiency = base.max(tidal);
    let heat = if icy {
        MAGMA_HEAT_PER_M3_ICE
    } else {
        MAGMA_HEAT_PER_M3_ROCK
    };
    // A dry rocky mantle's melt comes of its own part of the heat flux, not the crust's
    // radioactive heat, and only while it or its plumes stand above the solidus at its stiff
    // lid's base (F123).
    if mantle_keeps_heat(c) {
        let h = mantle_history(c);
        let age = c.f64("body.age");
        let crust = h.crust_share
            * radiogenic_per_kg(age)
            * c.f64("body.mass")
            * c.f64("body.rock_mass_fraction")
            / (4.0 * PI * r * r);
        let gate = h.melt_gate_at(age, c.f64("tectonics.plume_rise_count") >= 1.0);
        return (magma_rate(
            (c.f64("interior.surface_heat_flux") - crust).max(0.0),
            4.0 * PI * r * r,
            efficiency,
            heat,
        ) * gate)
            .into();
    }
    let rate = magma_rate(
        c.f64("interior.surface_heat_flux"),
        4.0 * PI * r * r,
        efficiency,
        heat,
    );
    if !icy {
        return rate.into();
    }
    // An ice mantle melts only where its conductive profile reaches the melting point
    // inside the body (F133).
    let t_s = c.f64("atmosphere.surface_temperature_mean");
    let t_m = c.f64("interior.mantle_temperature");
    let conductive = conductive_thickness(
        ice_conductivity(0.5 * (t_m + t_s)),
        (t_m - t_s).max(1.0),
        c.f64("interior.surface_heat_flux"),
    );
    (rate * crate::thermal::ice_melt_gate(r, conductive, t_s, t_m)).into()
}

// ------------------------------------------------------------------ atmosphere

pub fn atmosphere_mean_molecular_mass(c: &Ctx) -> Value {
    mean_molecular_mass(c.map("atmosphere.composition")).into()
}

pub fn atmosphere_column_mass(c: &Ctx) -> Value {
    (c.f64("atmosphere.surface_pressure") / c.f64("body.gravity")).into()
}

/// F52: an airless body's Bond albedo is its surface's own reflectance, the regolith
/// mixture at maturity and fresh over the share a plume's fallout keeps under frost;
/// both read from the previous pass of the graph, since the surface depends on the
/// temperature the albedo sets. Under an atmosphere the clouds and haze set it, which
/// the climate stage will compute (Phase 5); until then a guess by pressure class.
pub fn atmosphere_bond_albedo(c: &Ctx) -> Value {
    let p = c.f64("atmosphere.surface_pressure");
    if p >= AIRLESS_PRESSURE {
        return if p < THIN_ATMOSPHERE_PRESSURE {
            ALBEDO_THIN_ATMOSPHERE
        } else if p < THICK_ATMOSPHERE_PRESSURE {
            ALBEDO_EARTHLIKE
        } else {
            ALBEDO_THICK_ATMOSPHERE
        }
        .into();
    }
    // An airless surface is gardened, its basalt with the mantle's iron (F125) and its dirty
    // ice with its crust's rock (F126).
    use crate::stages::materials::{
        Library, end_member, is_ice, mixture_luminance, weathering_intensity,
    };
    let library = Library {
        feo: c.f64("body.mantle_feo"),
        gardened: true,
        impurity: c.f64("body.crust_impurity"),
        lag_depth: 0.0,
        weathering: weathering_intensity(
            c.f64("orbit.star_distance"),
            c.f64("cratering.diffusivity"),
        ),
    };
    let mature = mixture_luminance(&library, c.map("materials.regolith_mix")).1;
    let ice = end_member("water_ice").unwrap();
    let lum = |x: [f32; 3]| (x[0] + x[1] + x[2]) as f64 / 3.0;
    // F126: on an icy crust, the share its shell renewed in lanes or at a pole is clean ice,
    // with what impacts have brought it of the old crust's regolith.
    let mature = if is_ice(c.str("materials.bedrock")) {
        let renewed = c.f64("tectonics.ice_resurfaced_fraction").clamp(0.0, 1.0);
        let lanes =
            lum(ice.mature) + (mature - lum(ice.mature)) * REGOLITH_EXCHANGE * (1.0 - renewed);
        mature + (lanes - mature) * renewed
    } else {
        mature
    };
    // A plume's fallout is fresh frost.
    let frosted = crate::stages::cryo::mean_frost_share(c.f64("volcanism.plume_fallout_rate"));
    (mature + (lum(ice.fresh) - mature) * frosted).into()
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

pub fn tectonics_initial_crust_thickness(c: &Ctx) -> Value {
    if is_icy(c.str("body.composition_class")) {
        20_000.0
    } else {
        40_000.0
    }
    .into()
}

// ------------------------------------------------------------------ cratering (scaling subset)

/// F135: a body's place against its star's snow line, its distance over the distance
/// where its blackbody temperature is the snow line's, from its insolation: 1 AU is 0.372
/// of the Sun's.
fn snow_line_position(c: &Ctx) -> f64 {
    c.f64("cratering.snow_line_position")
}

pub fn cratering_snow_line_position(c: &Ctx) -> Value {
    (SNOW_LINE_WET_K / formation_temperature(c).max(1.0))
        .powi(2)
        .into()
}

fn ramp(x: f64, (a, b): (f64, f64)) -> f64 {
    let t = ((x - a) / (b - a)).clamp(0.0, 1.0);
    t * t * (3.0 - 2.0 * t)
}

/// The rocky leftovers' production function by place against the snow line (F135): the
/// lunar one inside ROCKY_MARS_X, the Mars one beyond.
pub fn cratering_production_function(c: &Ctx) -> Value {
    if snow_line_position(c) < ROCKY_MARS_X {
        "neukum_lunar"
    } else {
        "neukum_mars"
    }
    .into()
}

/// F135: the density of the impactors that make the body's craters today, each population's
/// by its share of the craters of 1 km.
pub fn cratering_impactor_density(c: &Ctx) -> Value {
    let rocky =
        c.f64("cratering.flux_relative_to_moon") * production_function(c).cumulative_1gyr(1.0);
    let icy = crate::crater::ImpactorCraters::new(
        c.f64("cratering.icy_flux"),
        crate::crater::ImpactorSizes::Comets,
        &crater_target(c, c.f64("cratering.impact_velocity")),
    )
    .cumulative(1.0);
    if rocky + icy <= 0.0 {
        return IMPACTOR_DENSITY_INNER.into();
    }
    ((rocky * IMPACTOR_DENSITY_INNER + icy * IMPACTOR_DENSITY_OUTER) / (rocky + icy)).into()
}

pub fn cratering_target_density(c: &Ctx) -> Value {
    if is_icy(c.str("body.composition_class")) {
        DENSITY_CRUST_ICE
    } else {
        2_700.0
    }
    .into()
}

/// F29, F135: the impactors' speed far from the body or its planet, a share of the orbital
/// speed about the star that is the larger the farther inside the snow line.
fn encounter_speed(c: &Ctx) -> f64 {
    let d = c.f64("orbit.star_distance");
    let v_orb = orbital_speed(c.f64("star.mass") * M_SUN, d);
    let x = snow_line_position(c).max(1e-6);
    (VINF_FACTOR_INNER * (VINF_REFERENCE_X / x).sqrt()).max(VINF_FACTOR_MIN) * v_orb
}

pub fn cratering_impact_velocity(c: &Ctx) -> Value {
    let v_inf = encounter_speed(c);
    let v_esc = c.f64("body.escape_velocity");
    let v2 = if c.str("orbit.primary") == "planet" {
        // A moon meets the impactors its planet focused, at its own orbital speed besides
        // (Zahnle et al. 2003): three times the square of that speed.
        let v_moon2 = G * c.f64("parent.mass") / c.f64("orbit.semi_major_axis");
        v_inf * v_inf + 3.0 * v_moon2 + v_esc * v_esc
    } else {
        v_inf * v_inf + v_esc * v_esc
    };
    v2.sqrt().into()
}

/// F135: the rocky leftovers' flux per unit area relative to the Moon's, by place against
/// the snow line, times the late flux's multiplier.
pub fn cratering_flux_relative_to_moon(c: &Ctx) -> Value {
    let x = snow_line_position(c);
    let belt = 1.0 + (ROCKY_BELT_FLUX - 1.0) * ramp(x, ROCKY_BELT_RAMP);
    let fade = 1.0 - ramp(x, ROCKY_OUTER_FADE);
    (belt * fade * c.f64("cratering.late_flux_multiplier")).into()
}

/// F135: the outer system's impactors at or above 1 km per km² per Gyr at the body today,
/// focused by its planet (or by itself), from the snow line out, times the late flux's
/// multiplier.
pub fn cratering_icy_flux(c: &Ctx) -> Value {
    let x = snow_line_position(c);
    let share = ramp(x, ICY_SNOW_RAMP);
    if share <= 0.0 {
        return 0.0.into();
    }
    let (x0, rise) = ICY_FLUX_RISE;
    let unfocused = ICY_FLUX * (x / x0).max(1.0).powf(rise);
    let v_inf = encounter_speed(c).max(1.0);
    let v_esc2 = if c.str("orbit.primary") == "planet" {
        2.0 * G * c.f64("parent.mass") / c.f64("orbit.semi_major_axis")
    } else {
        c.f64("body.escape_velocity").powi(2)
    };
    (unfocused * (1.0 + v_esc2 / (v_inf * v_inf)) * share * c.f64("cratering.late_flux_multiplier"))
        .into()
}

/// F136: what a crater forms in on the body, struck at `speed`, m/s: its gravity, the
/// strength of its bedrock, cold ice's or rock's, its density and its transition.
fn crater_target(c: &Ctx, speed: f64) -> crate::crater::Target {
    crate::crater::Target {
        gravity: c.f64("body.gravity"),
        speed,
        strength: crate::crater::by_ice(
            TARGET_STRENGTH_ROCK,
            TARGET_STRENGTH_ICE,
            c.f64("cratering.target_ice"),
        ),
        density: c.f64("cratering.target_density"),
        simple_complex: c.f64("cratering.simple_complex_transition"),
    }
}

/// F135, F136: the crater a 1 km impactor of the outer system makes on the body, km.
pub fn cratering_icy_crater_1km(c: &Ctx) -> Value {
    (crater_target(c, c.f64("cratering.impact_velocity")).crater(1e3) / 1e3).into()
}

/// F136: the speed a moon's debris strikes it at, m/s: DEBRIS_ECCENTRICITY times its
/// orbital speed about its planet; nothing for a body that circles its star.
pub fn cratering_debris_velocity(c: &Ctx) -> Value {
    if c.str("orbit.primary") != "planet" {
        return 0.0.into();
    }
    (DEBRIS_ECCENTRICITY * (G * c.f64("parent.mass") / c.f64("orbit.semi_major_axis")).sqrt())
        .into()
}

/// F136: the fragments at or above DEBRIS_BREAK_KM per km² that strike a moon after the last
/// impact it could not hold, as the square of the comets' impact speed (stated choice);
/// none where the record starts when the crust solidifies, or on a body that circles its
/// star, whose debris disperses.
pub fn cratering_debris_fragments(c: &Ctx) -> Value {
    let reset = c.f64("cratering.record_start") > c.f64("surface.solidification_time") + 1e-9;
    if c.str("orbit.primary") != "planet" || !reset {
        return 0.0.into();
    }
    (DEBRIS_FRAGMENTS * (c.f64("cratering.impact_velocity") / DEBRIS_REFERENCE_SPEED).powi(2))
        .into()
}

/// F135: around another star, the early bombardment's timescale against the Sun's and its
/// amount per unit area against the Sun's.
fn early_scalings(c: &Ctx) -> (f64, f64) {
    let l = c.f64("star.luminosity").max(1e-9);
    let m = c.f64("star.mass").max(1e-6);
    let (a, b) = EARLY_TIME_EXPONENTS;
    (l.powf(a) * m.powf(b), m / l)
}

pub fn cratering_icy_early_total(c: &Ctx) -> Value {
    let (_, solids) = early_scalings(c);
    (ICY_EARLY_TOTAL * solids / c.f64("cratering.late_flux_multiplier").max(1e-6)).into()
}

pub fn cratering_icy_early_timescale(c: &Ctx) -> Value {
    (ICY_EARLY_TAU_GYR * early_scalings(c).0).into()
}

/// F135, F136: the body's bombardment from the graph, its record starting at
/// `record_start`, with a moon's debris or without.
fn bombardment_from(c: &Ctx, record_start: f64, with_debris: bool) -> crate::crater::Bombardment {
    crate::crater::BombardmentParams::from_nodes(
        |id| c.f64(id),
        c.str("cratering.production_function"),
        c.f64("cratering.target_ice"),
        record_start,
        with_debris,
    )
    .build()
}

/// F135: when the surface begins its crater record: when it solidifies, or after the last
/// impact too large for the body to hold, which shattered and resurfaced it, whichever is
/// later. The last such impact is taken at its median time: impacts come at random, so
/// none comes after a time with the chance exp(-n), n the number still expected, and the
/// chance is one half where ln 2 are still to come.
pub fn cratering_record_start(c: &Ctx) -> Value {
    let start = c.f64("surface.solidification_time");
    let age = c.f64("body.age");
    let radius = c.f64("body.radius");
    let hold = crate::crater::hold_diameter(radius, c.f64("cratering.target_ice")) / 1e3;
    let area = 4.0 * PI * radius * radius / 1e6;
    let b = bombardment_from(c, 0.0, false);
    let after = |t: f64| b.density_above(hold, t, age) * area;
    let median = std::f64::consts::LN_2;
    if after(start) < median {
        return start.into();
    }
    let (mut lo, mut hi) = (start, age);
    for _ in 0..60 {
        let mid = 0.5 * (lo + hi);
        if after(mid) > median {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    (0.5 * (lo + hi)).into()
}

/// Audit S9: how far the target a crater forms in behaves as ice rather than rock, by the
/// share of its bedrock that is ice: rock below TARGET_ICE_SHARE's first share, ice above
/// its second, by degrees between (Durham et al. 1992). A body's dirty ice is its crust's
/// rock and ice (F126); its clean ice a trace of rock; anything else rock.
pub fn cratering_target_ice(c: &Ctx) -> Value {
    let bedrock = c.str("materials.bedrock");
    if !crate::stages::materials::is_ice(bedrock) {
        return 0.0.into();
    }
    let rock = if bedrock == "water_ice_dirty" {
        c.f64("body.crust_impurity")
    } else {
        crate::stages::materials::non_ice_fraction(bedrock)
    };
    let x = ((1.0 - rock - TARGET_ICE_SHARE.0) / (TARGET_ICE_SHARE.1 - TARGET_ICE_SHARE.0))
        .clamp(0.0, 1.0);
    (x * x * (3.0 - 2.0 * x)).into()
}

pub fn cratering_simple_complex_transition(c: &Ctx) -> Value {
    // The transition is a property of the target the crater forms in: the surface
    // material, not the bulk composition.
    transition_diameter(c.f64("body.gravity"), c.f64("cratering.target_ice")).into()
}

pub fn cratering_peak_ring_onset(c: &Ctx) -> Value {
    (PEAK_RING_MULTIPLE * c.f64("cratering.simple_complex_transition")).into()
}

pub fn cratering_multiring_onset(c: &Ctx) -> Value {
    (MULTIRING_MULTIPLE * c.f64("cratering.simple_complex_transition")).into()
}

/// F135: the gardening diffusivity, the Moon's scaled by gravity and by the body's present
/// production of craters of 1 km against the Moon's, the small impacts that do the gardening
/// coming as often as those.
pub fn cratering_diffusivity(c: &Ctx) -> Value {
    let flux = bombardment_from(c, 0.0, false).steady_1gyr(1.0)
        / crate::crater::ProductionFunction::lunar().cumulative_1gyr(1.0);
    let mut k = DIFFUSIVITY_LUNAR * (c.f64("body.gravity") / GRAVITY_MOON).sqrt() * flux;
    if c.f64("atmosphere.surface_pressure") > DIFFUSIVITY_ATMOSPHERE_PRESSURE {
        k *= DIFFUSIVITY_ATMOSPHERE_MULTIPLE;
    }
    k.into()
}

pub fn cratering_ray_fading_timescale(c: &Ctx) -> Value {
    if c.f64("atmosphere.surface_pressure") > DIFFUSIVITY_ATMOSPHERE_PRESSURE {
        RAY_FADING_GYR_ATMOSPHERE
    } else {
        (RAY_FADING_GYR_1AU / maturing_rate(c)).min(RAY_FADING_MAX_GYR)
    }
    .into()
}

/// F214: how fast a fresh surface matures against the Moon's, by its weathering and its
/// gardening.
fn maturing_rate(c: &Ctx) -> f64 {
    crate::stages::materials::maturing_rate(
        c.f64("orbit.star_distance"),
        c.f64("cratering.diffusivity"),
    )
    .max(1e-9)
}

fn production_function(c: &Ctx) -> crate::crater::ProductionFunction {
    crate::crater::ProductionFunction::by_name(c.str("cratering.production_function"))
}

/// F31, F135: the rocky leftovers' early rate at formation relative to the steady rate,
/// from the chronology, as much more per unit area as another star's disc leaves and
/// over as much less time as its orbits take, against a steady flux the late multiplier
/// scales.
pub fn cratering_early_flux_ratio(c: &Ctx) -> Value {
    let (time, solids) = early_scalings(c);
    (crate::crater::production::early_flux_ratio(production_function(c).chronology) * solids
        / time.max(1e-9)
        / c.f64("cratering.late_flux_multiplier").max(1e-6))
    .into()
}

/// F31, F135: e-folding time of the rocky leftovers' early bombardment, from the
/// chronology, on another star's orbital times.
pub fn cratering_early_flux_decay(c: &Ctx) -> Value {
    (early_scalings(c).0 / production_function(c).chronology.1).into()
}

pub fn cratering_flux_history(c: &Ctx) -> Value {
    // The production of craters of 1 km over the body's life relative to the present,
    // every population's (F135).
    let age = c.f64("body.age");
    let b = bombardment_from(c, c.f64("cratering.record_start"), true);
    let rate = |t: f64| -> f64 {
        b.sources
            .iter()
            .map(|s| {
                s.scale * s.pf.cumulative_1gyr(1.0) * s.flux.rate1(t) / s.flux.steady.max(1e-30)
            })
            .sum()
    };
    let now = rate(age).max(1e-30);
    let n = 48;
    let points = (0..=n)
        .map(|k| {
            let t = age * k as f64 / n as f64;
            (t, rate(t) / now)
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
    // When the mantle last broke the lid (F121).
    let age = c.f64("body.age");
    mantle_history(c)
        .last_overturn()
        .map_or(0.0, |t| (age - t).max(0.0))
        .into()
}

pub fn tectonics_tessera_fraction(c: &Ctx) -> Value {
    if c.str("interior.tectonic_regime") == "episodic_lid" {
        TESSERA_FRACTION.into()
    } else {
        0.0.into()
    }
}

pub fn volcanism_corona_count(c: &Ctx) -> Value {
    // A rock lid's coronae over its mantle's diapirs (F145), as the loop makes them.
    if matches!(
        c.str("interior.tectonic_regime"),
        "episodic_lid" | "stagnant_lid"
    ) {
        let mantle = mantle_history(c);
        let thermal = crate::thermal::ThermalHistory::new(&thermal_inputs(c));
        // The plumes and the upwelling limbs between them, as the loop makes them (F150,
        // F151).
        let (rises, turn) = crate::stages::plume::place_rises(c, &thermal).unwrap_or_default();
        let support = crate::convection::Convection::from_graph(c, &thermal)
            .map(|convection| crate::convection::Support { convection, turn });
        let limbs =
            crate::stages::stress::Limbs::along(&rises, support.as_ref(), c.f64("body.radius"));
        return (crate::stages::corona::rock_coronae(c, &mantle, &rises, &limbs).len() as f64)
            .into();
    }
    // One over each upwelling of an ice mantle's convection, where the lid breaks (F133).
    match ice_corona_episode(c) {
        Some(e) if e.strength > CORONA_STRENGTH_FLOOR => {
            e.convection.upwellings.round().max(1.0).into()
        }
        _ => 0.0.into(),
    }
}

/// F133: the convection a tidal heating episode drives in an ice mantle on rock, at the
/// episode's thinnest lid.
fn ice_corona_episode(c: &Ctx) -> Option<crate::stages::corona::IceEpisode> {
    if c.str("interior.tectonic_regime") != "ice_shell_rock" {
        return None;
    }
    let history = crate::thermal::ThermalHistory::new(&thermal_inputs(c));
    crate::stages::corona::IceEpisode::new(
        &history,
        c.f64("body.radius"),
        c.f64("body.ice_layer_thickness"),
        c.f64("body.gravity"),
    )
}

pub fn volcanism_corona_strength(c: &Ctx) -> Value {
    ice_corona_episode(c).map_or(0.0, |e| e.strength).into()
}

/// The coronae of an ice mantle take their time from the tidal heating as the lanes do
/// (F44, F133).
fn ice_corona_window(c: &Ctx) -> Option<(f64, f64)> {
    if c.str("interior.tectonic_regime") != "ice_shell_rock"
        || c.f64("volcanism.corona_count") <= 0.0
    {
        return None;
    }
    resurfacing_episode(c)
}

pub fn volcanism_corona_age(c: &Ctx) -> Value {
    let age = c.f64("body.age");
    ice_corona_window(c)
        .map_or(0.0, |(mean, _)| (age - mean).max(0.0))
        .into()
}

pub fn volcanism_corona_duration(c: &Ctx) -> Value {
    ice_corona_window(c).map_or(0.0, |(_, width)| width).into()
}

pub fn volcanism_dome_count(c: &Ctx) -> Value {
    if c.str("interior.tectonic_regime") != "episodic_lid" {
        return 0.0.into();
    }
    // The domes stand on an episodic lid's young plains, by its degree (audit S1).
    let r = c.f64("body.radius");
    (DOMES_PER_AREA * 4.0 * PI * r * r * c.f64("interior.episodic_degree"))
        .round()
        .into()
}

pub fn tectonics_plume_count(c: &Ctx) -> Value {
    let regime = c.str("interior.tectonic_regime");
    if !matches!(regime, "stagnant_lid" | "episodic_lid") || c.f64("body.radius") < PLUME_MIN_RADIUS
    {
        return 0.0.into();
    }
    // A plume carries a share of the heat its mantle loses while it rises, so their number
    // follows what the mantle has lost on average over its life (F122), not what it loses
    // today, long after its plumes rose.
    let r = c.f64("body.radius");
    let h = mantle_history(c);
    let mean = if h.heat_flux.is_empty() {
        c.f64("interior.surface_heat_flux")
    } else {
        h.mean_heat_flux()
    };
    (mean * 4.0 * PI * r * r / PLUME_POWER_PER_PLUME)
        .floor()
        .clamp(0.0, 3.0)
        .into()
}

pub fn tectonics_rise_diameter(c: &Ctx) -> Value {
    if c.f64("tectonics.plume_rise_count") < 1.0 {
        return 0.0.into();
    }
    // An episodic lid's narrower rises by its degree (audit S1).
    let episodic = c.f64("interior.episodic_degree");
    let fraction = if episodic <= 0.0 {
        PLUME_RISE_DIAMETER_FRACTION
    } else if episodic >= 1.0 {
        PLUME_RISE_DIAMETER_FRACTION_EPISODIC
    } else {
        PLUME_RISE_DIAMETER_FRACTION
            * (PLUME_RISE_DIAMETER_FRACTION_EPISODIC / PLUME_RISE_DIAMETER_FRACTION).powf(episodic)
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

/// Whether an ice shell's bedrock is its ocean's clean ice: where it renews its whole
/// surface more than half over (audit S6); its regolith takes the renewed ice by the
/// renewal's degree (`materials_regolith_mix`).
fn shell_renews_everywhere(c: &Ctx) -> bool {
    wholesale_renewal(c) > 0.5
}

pub fn cratering_retention_age(c: &Ctx) -> Value {
    let age = c.f64("body.age");
    let regime = c.str("interior.tectonic_regime");
    let p = c.f64("atmosphere.surface_pressure");
    let stagnant = if p > SPACE_WEATHERING_MAX_PRESSURE {
        RETENTION_THIN_ATMOSPHERE_GYR
    } else {
        age
    };
    let mut retention = match regime {
        "mobile_lid" => RETENTION_MOBILE_LID_GYR,
        // The surface the lid's last overturn laid (F121), within what the lid would keep
        // without it (audit S1).
        "episodic_lid" => c.f64("tectonics.resurfacing_event_age").min(stagnant),
        // A thin shell renews its surface within RETENTION_ICE_OCEAN_GYR over the degree
        // of its renewal (audit S6).
        "ice_shell_ocean" => (RETENTION_ICE_OCEAN_GYR / wholesale_renewal(c).max(1e-9)).min(age),
        _ => stagnant,
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
        (SPACE_WEATHERING_GYR_1AU / maturing_rate(c)).min(SPACE_WEATHERING_MAX_GYR)
    }
    .into()
}

// ------------------------------------------------------------------ materials

pub fn materials_bedrock(c: &Ctx) -> Value {
    let class = c.str("body.composition_class");
    let regime = c.str("interior.tectonic_regime");
    // F125: a tidally heated body's plumes coat it with sulfur; plates recycle a crust with
    // water into a felsic one; otherwise what the magma ocean floated, or the basalt of the
    // body's own volcanism.
    let bedrock = match class {
        "rock" | "iron_rich" => {
            if c.f64("interior.tidal_flux") > IO_REGIME_TIDAL_FLUX {
                "sulfur_deposit"
            } else if regime == "mobile_lid" {
                "felsic"
            } else {
                let (plagioclase, graphite) = flotation(c);
                if plagioclase >= 0.5 && plagioclase >= graphite {
                    "anorthosite"
                } else if graphite >= 0.5 {
                    "graphite_bearing_crust"
                } else {
                    "basalt"
                }
            }
        }
        // F126: an icy crust is the rock and ice the body accreted, unless its shell renews
        // itself everywhere and shows the ocean's clean ice.
        _ => {
            if shell_renews_everywhere(c) {
                "water_ice"
            } else {
                "water_ice_dirty"
            }
        }
    };
    bedrock.into()
}

pub fn materials_regolith_mix(c: &Ctx) -> Value {
    let mut m = BTreeMap::new();
    // F125: a rocky crust's surface mixes what its magma ocean floated with its basalt, by
    // their shares, so that a body between a flotation crust and none has a surface between.
    if rocky_crust(c) {
        let (plagioclase, graphite) = flotation(c);
        let basalt = (1.0 - plagioclase - graphite).max(0.0);
        for (id, share) in [
            ("anorthosite", plagioclase),
            ("graphite_bearing_crust", graphite),
            ("basalt", basalt),
        ] {
            if share > 1e-6 {
                m.insert(id.to_string(), share);
            }
        }
        if !m.is_empty() {
            return Value::Map(m);
        }
    }
    // F126: an airless icy crust's regolith carries the lag its warm latitudes leave, its
    // rock, which impact ejecta spread over the old crust: the crust and the lag by the
    // share of the surface whose ice sublimes into a lag over the surface's age. Impacts
    // also trade regolith between the old crust and the ice its shell renewed, the old
    // crust taking REGOLITH_EXCHANGE of the renewed ice by its share of the body.
    let bedrock = c.str("materials.bedrock");
    if crate::stages::materials::is_ice(bedrock)
        && c.f64("atmosphere.surface_pressure") < AIRLESS_PRESSURE
    {
        use crate::stages::materials::{LAG_MATERIAL, Library};
        let library = Library {
            feo: c.f64("body.mantle_feo"),
            gardened: true,
            impurity: c.f64("body.crust_impurity"),
            lag_depth: 0.0,
            weathering: 0.0,
        };
        if let (Some(crust), Some(lag)) = (library.member(bedrock), library.member(LAG_MATERIAL)) {
            let lum = |x: [f32; 3]| (x[0] + x[1] + x[2]) as f64 / 3.0;
            let cover = crate::stages::ice::lag_cover(
                c.f64("orbit.insolation"),
                library.non_ice(bedrock),
                c.f64("cratering.retention_age"),
                (lum(crust.fresh), lum(crust.mature), lum(lag.mature)),
            );
            // The ice the lanes renewed, traded onto the old crust, and the ice a thin
            // shell renews everywhere, by its degree (audit S6).
            let renewed = REGOLITH_EXCHANGE * c.f64("tectonics.ice_resurfaced_fraction")
                + wholesale_renewal(c);
            let cover = if cover > 1e-4 { cover } else { 0.0 };
            let renewed = if renewed > 1e-4 {
                renewed.min(1.0)
            } else {
                0.0
            };
            if cover > 0.0 || renewed > 0.0 {
                for (id, share) in [
                    (bedrock, (1.0 - cover) * (1.0 - renewed)),
                    (LAG_MATERIAL, cover * (1.0 - renewed)),
                    ("water_ice", renewed),
                ] {
                    if share > 1e-6 {
                        *m.entry(id.to_string()).or_insert(0.0) += share;
                    }
                }
                return Value::Map(m);
            }
        }
    }
    m.insert(bedrock.to_string(), 1.0);
    Value::Map(m)
}

pub fn materials_lava(c: &Ctx) -> Value {
    // F125: a rocky body's lavas are its basalt, the iron of its mantle in them. F126:
    // cryolava freezes to the ice of its liquid, the rock having settled out of it.
    let lava = match c.str("body.composition_class") {
        "rock" | "iron_rich" => "basalt",
        _ => {
            if c.str("volatiles.working_liquid") == "ammonia_water" {
                "ammonia_water_ice"
            } else {
                "water_ice"
            }
        }
    };
    lava.into()
}

pub fn materials_deep_bedrock(c: &Ctx) -> Value {
    let regime = c.str("interior.tectonic_regime");
    let deep = match c.str("body.composition_class") {
        // F125: under a flotation crust its mafic lower crust, under a reduced body's lavas
        // the graphite flotation crust, otherwise the basalt.
        "rock" | "iron_rich" if regime != "mobile_lid" => {
            let (plagioclase, graphite) = flotation(c);
            if plagioclase >= 0.5 && plagioclase >= graphite {
                "norite"
            } else if graphite >= 0.5 {
                "graphite_rich_crust"
            } else {
                "basalt"
            }
        }
        "rock" | "iron_rich" => "basalt",
        // Under a shell that floats on an ocean the deep ice is clean; a crust over rock
        // is the same mixture all the way down, and a shell partly afloat is between them
        // by its ocean's degree (audit S1, `materials::deep_clean_share`).
        "rock_ice" | "ice_dominated" => {
            if c.f64("interior.ocean_degree") >= 1.0 {
                "water_ice"
            } else {
                c.str("materials.bedrock")
            }
        }
        _ => c.str("materials.bedrock"),
    };
    deep.into()
}

/// F126: the share of the surface layer an airless icy crust's lag has taken over its age
/// (`ice::lag_depth_share`): where it is deep, fresh craters and slopes bring up lag.
pub fn materials_lag_depth_share(c: &Ctx) -> Value {
    use crate::stages::materials::{LAG_MATERIAL, Library, is_ice};
    let bedrock = c.str("materials.bedrock");
    if !is_ice(bedrock) || c.f64("atmosphere.surface_pressure") >= AIRLESS_PRESSURE {
        return 0.0.into();
    }
    let library = Library {
        feo: 0.0,
        gardened: true,
        impurity: c.f64("body.crust_impurity"),
        lag_depth: 0.0,
        weathering: 0.0,
    };
    let (Some(crust), Some(lag)) = (library.member(bedrock), library.member(LAG_MATERIAL)) else {
        return 0.0.into();
    };
    let lum = |x: [f32; 3]| (x[0] + x[1] + x[2]) as f64 / 3.0;
    crate::stages::ice::lag_depth_share(
        c.f64("orbit.insolation"),
        library.non_ice(bedrock),
        c.f64("cratering.retention_age"),
        (lum(crust.fresh), lum(crust.mature), lum(lag.mature)),
        c.f64("materials.deep_bedrock_depth"),
    )
    .into()
}

pub fn materials_deep_bedrock_depth(c: &Ctx) -> Value {
    // F125: the graphite flotation crust lies under a reduced body's lavas; a flotation
    // crust's mafic base, or a basaltic crust's, deeper.
    match c.str("body.composition_class") {
        "rock" | "iron_rich" => {
            if flotation(c).1 >= 0.5 {
                DEEP_BEDROCK_DEPTH_GRAPHITE
            } else {
                DEEP_BEDROCK_DEPTH_ROCK
            }
        }
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
    fn a_small_moon_s_porous_core_dissipates_its_tides() {
        // F68: Enceladus's rocky interior stays porous, about a quarter of it pores, and
        // dissipates 10 to 30 GW (Choblet et al. 2017), which brings its shell to its
        // measured thickness; the large moons' interiors are consolidated.
        let (g, _) = crate::presets::preset("Enceladus").unwrap().to_graph();
        let porosity = g.f64("interior.core_porosity");
        assert!((0.2..0.3).contains(&porosity), "{porosity}");
        let core = g.f64("interior.core_tidal_power");
        assert!((10e9..30e9).contains(&core), "{core}");
        let shell = g.f64("interior.ice_shell_thickness");
        assert!((15e3..30e3).contains(&shell), "{shell}");
        for name in ["Europa", "Ganymede", "Callisto", "Titan"] {
            let (g, _) = crate::presets::preset(name).unwrap().to_graph();
            assert!(
                g.f64("interior.core_porosity") < POROSITY_CONSOLIDATED,
                "{name}"
            );
            assert_eq!(g.f64("interior.core_tidal_power"), 0.0, "{name}");
        }
    }

    #[test]
    fn composition_classes_match_the_reference_bodies() {
        let class = |rho: f64, mass_earth: f64| {
            let rho_u = rho / compression_factor(mass_earth * M_EARTH);
            let (core, _, ice) = composition_fractions(rho_u, 1.0);
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
    fn a_body_formed_inside_the_snow_line_carries_no_ice() {
        // F122: the Moon and Mars formed too warm for water ice, so the deficit of their
        // densities below the chondritic mix is iron they lack: the Moon is all rock (its
        // core is one or two percent of its mass, Weber et al. 2011), Mars's core under a
        // quarter of its mass (Stähler et al. 2021). Bodies formed beyond the snow line keep
        // their ice, and between the two temperatures the shares blend.
        let fractions = |name: &str| {
            let (g, _) = crate::presets::preset(name).unwrap().to_graph();
            (
                g.f64("body.core_mass_fraction"),
                g.f64("body.ice_mass_fraction"),
            )
        };
        let (core, ice) = fractions("Moon");
        assert!(ice == 0.0 && core < 0.05, "Moon core {core} ice {ice}");
        let (core, ice) = fractions("Mars");
        assert!(
            ice == 0.0 && core > 0.1 && core < 0.25,
            "Mars core {core} ice {ice}"
        );
        for name in ["Europa", "Ganymede", "Ceres"] {
            let (_, ice) = fractions(name);
            assert!(ice > 0.1, "{name} ice {ice}");
        }
        // Between pure rock and the chondritic mix a dry body lacks iron where a wet one
        // carries ice beside the chondritic share.
        let (wet, half, dry) = (
            composition_fractions(3_600.0, 1.0),
            composition_fractions(3_600.0, 0.5),
            composition_fractions(3_600.0, 0.0),
        );
        assert!(
            dry.2 == 0.0 && wet.2 > 0.0 && dry.0 < wet.0,
            "{wet:?} {dry:?}"
        );
        assert!((half.0 - 0.5 * (wet.0 + dry.0)).abs() < 1e-12);
        assert!((half.2 - 0.5 * wet.2).abs() < 1e-12);
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
