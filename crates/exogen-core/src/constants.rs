// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! Physical constants and model constants.
//!
//! Physical constants are CODATA 2018 or IAU nominal values. Model constants carry their
//! source and, where the value is a provisional fit awaiting calibration in a later phase,
//! the word "provisional" and the phase that owns it. Nothing else in the engine hard-codes
//! a number that appears here.

// ---------------------------------------------------------------- physical constants (SI)

/// Newtonian constant of gravitation, CODATA 2018.
pub const G: f64 = 6.674_30e-11;
/// Stefan–Boltzmann constant, CODATA 2018 (exact in the 2019 SI).
pub const SIGMA_SB: f64 = 5.670_374_419e-8;
/// Boltzmann constant, exact in the 2019 SI.
pub const K_B: f64 = 1.380_649e-23;
/// Atomic mass unit, CODATA 2018.
pub const M_U: f64 = 1.660_539_066_60e-27;
/// Astronomical unit, IAU 2012, exact.
pub const AU: f64 = 1.495_978_707e11;
/// Julian year in seconds.
pub const YEAR: f64 = 3.155_76e7;
/// Gigayear in seconds.
pub const GYR: f64 = 3.155_76e16;
/// Megayear in seconds.
pub const MYR: f64 = 3.155_76e13;

/// Solar mass, IAU 2015 nominal (GM_sun / G).
pub const M_SUN: f64 = 1.988_47e30;
/// Solar luminosity, IAU 2015 nominal.
pub const L_SUN: f64 = 3.828e26;
/// Solar radius, IAU 2015 nominal.
pub const R_SUN: f64 = 6.957e8;
/// Solar effective temperature implied by the nominal L and R.
pub const T_SUN: f64 = 5772.0;
/// Present solar age in Gyr.
pub const AGE_SUN_GYR: f64 = 4.57;

/// Earth mass, IAU 2015 nominal (GM_earth / G).
pub const M_EARTH: f64 = 5.972_2e24;
/// Earth mean radius.
pub const R_EARTH: f64 = 6.371_0e6;
/// Jupiter mass, IAU 2015 nominal.
pub const M_JUP: f64 = 1.898_13e27;
/// Jupiter equatorial radius, IAU 2015 nominal.
pub const R_JUP: f64 = 7.149_2e7;

// ---------------------------------------------------------------- stellar relations

/// Eker et al. 2018 six-piece main-sequence mass–luminosity relation, log L = a log M + b
/// in solar units. Rows: (upper mass bound, a, b). Valid 0.179 to 31 solar masses.
pub const MASS_LUMINOSITY: [(f64, f64, f64); 6] = [
    (0.45, 2.028, -0.976),
    (0.72, 4.572, -0.102),
    (1.05, 5.743, -0.007),
    (2.40, 4.329, 0.010),
    (7.0, 3.967, 0.093),
    (31.0, 2.865, 1.105),
];
/// The fit gives 0.984 L_sun at one solar mass; this factor calibrates it to the Sun.
pub const MASS_LUMINOSITY_SOLAR_NORMALISATION: f64 = 1.0 / 0.984;
/// Eker et al. 2018 main-sequence mass–radius relation R = a M² + b M + c for 0.179 to
/// 1.5 solar masses. Above the break the radius continues as R(1.5) (M / 1.5)^0.555, the
/// slope of the upper branch of Demircan and Kahraman 1991.
pub const MASS_RADIUS_QUADRATIC: (f64, f64, f64) = (0.438, 0.479, 0.075);
pub const MASS_RADIUS_BREAK: f64 = 1.5;
pub const MASS_RADIUS_UPPER_EXPONENT: f64 = 0.555;

/// Gough 1981: L(t) = L_present / (1 + 0.4 (1 - t / t_present)) for the Sun, with the
/// present age 4.57 Gyr. Generalised to other masses by scaling the reference age with the
/// main-sequence lifetime (docs/parameters.md F3).
pub const GOUGH_COEFFICIENT: f64 = 0.4;
/// Reference age as a fraction of the main-sequence lifetime: 4.57 Gyr of a 10 Gyr lifetime.
pub const GOUGH_REFERENCE_FRACTION: f64 = 0.457;
/// Main-sequence lifetime of the Sun used by the lifetime scaling, Gyr.
pub const MAIN_SEQUENCE_LIFETIME_SUN_GYR: f64 = 10.0;

/// Ribas et al. 2005: XUV flux of solar analogues scales as age^-1.23 beyond 0.1 Gyr.
pub const XUV_AGE_EXPONENT: f64 = -1.23;
pub const XUV_SATURATION_AGE_GYR: f64 = 0.1;

// ---------------------------------------------------------------- bulk composition (F8)

/// Zero-pressure densities of the three end members, kg/m³. Iron alloy, silicate rock,
/// mixed ices. Provisional until the compressible three-layer model in Phase 4.
pub const DENSITY_IRON: f64 = 7_900.0;
pub const DENSITY_ROCK: f64 = 3_300.0;
pub const DENSITY_ICE: f64 = 1_000.0;
/// Chondritic core-to-silicate mass ratio assumed when the body also carries ice
/// (Earth: 0.325 / 0.675 = 0.48).
pub const CHONDRITIC_CORE_TO_ROCK: f64 = 0.48;
/// Compression correction: bulk density = uncompressed density x (1 + a (M / M_earth)^b).
/// Fit to Earth (1.30), Venus (1.30), Mars (1.05), Mercury (1.02). Provisional, Phase 4.
pub const COMPRESSION_A: f64 = 0.30;
pub const COMPRESSION_B: f64 = 0.6;
/// Class thresholds on the derived fractions.
pub const IRON_RICH_CORE_FRACTION: f64 = 0.5;
pub const ROCK_MAX_ICE_FRACTION: f64 = 0.10;
pub const ROCK_ICE_MAX_ICE_FRACTION: f64 = 0.45;

// ---------------------------------------------------------------- interior

/// Present-day radiogenic heat production per kilogram of silicate, W/kg, and the share of
/// each isotope with its half-life in Gyr. Bulk silicate Earth abundances (McDonough and
/// Sun 1995) with the per-isotope heat rates of Turcotte and Schubert (2002): U 20.8 ppb,
/// Th 79.5 ppb, K 240 ppm give 4.96e-12 W/kg; shares U238 0.389, U235 0.017, Th232 0.431,
/// K40 0.161; half-lives 4.468, 0.704, 14.05, 1.248 Gyr.
pub const RADIOGENIC_NOW_W_PER_KG: f64 = 4.96e-12;
pub const RADIOGENIC_ISOTOPES: [(f64, f64); 4] = [
    (0.389, 4.468),
    (0.017, 0.704),
    (0.431, 14.05),
    (0.161, 1.248),
];
/// Age of the solar system's rock inventory used to place a body's present on the decay
/// curves, Gyr.
pub const SOLAR_SYSTEM_AGE_GYR: f64 = 4.567;

/// Secular cooling contribution as a multiple of the radiogenic flux, interpolated linearly
/// in radius between (500 km, 0.2) and (3000 km, 1.0). Earth's total flux is about twice its
/// radiogenic flux. Provisional, Phase 3.
pub const SECULAR_FRACTION_SMALL: f64 = 0.2;
pub const SECULAR_FRACTION_LARGE: f64 = 1.0;
pub const SECULAR_RADIUS_SMALL: f64 = 5.0e5;
pub const SECULAR_RADIUS_LARGE: f64 = 3.0e6;

/// Love number k2 and dissipation factor Q by composition class. Provisional, Phase 3 and 4.
/// Rock: k2 0.3, Q 25 (Io's measured k2/Q is 0.015; Lainey et al. 2009). Rock with a
/// surface ocean: Q 12 (Earth). Rock-ice shell over an ocean: k2 0.25, Q 50. Ice-dominated:
/// k2 0.2, Q 20 (Enceladus needs about 0.01 to reach its observed 16 GW).
pub const K2_ROCK: f64 = 0.3;
pub const Q_ROCK: f64 = 25.0;
pub const Q_ROCK_SURFACE_OCEAN: f64 = 12.0;
pub const K2_ROCK_ICE: f64 = 0.25;
pub const Q_ROCK_ICE: f64 = 50.0;
pub const K2_ICE: f64 = 0.2;
pub const Q_ICE: f64 = 20.0;

/// Mantle temperature model: T_m = a + b (q / q_earth)^0.5 K. Earth 1600 K. Provisional, Phase 3.
pub const MANTLE_T_BASE: f64 = 1_400.0;
pub const MANTLE_T_SCALE: f64 = 200.0;
pub const EARTH_HEAT_FLUX: f64 = 0.087;
/// Thermal conductivity of rock, W/m/K, and the ice conductivity law k = ICE_K_NUM / T.
pub const K_ROCK: f64 = 3.0;
pub const ICE_K_NUM: f64 = 651.0;
/// Fraction of the surface heat flux that passes through the mantle lithosphere of a rocky
/// body (the rest is produced in the crust). Provisional, Phase 3.
pub const MANTLE_FLUX_FRACTION: f64 = 0.5;
/// Isotherms bounding elastic behaviour, K.
pub const ELASTIC_ISOTHERM_ROCK: f64 = 650.0;
pub const ELASTIC_ISOTHERM_ICE: f64 = 170.0;
/// Crustal strength, Pa, and crustal density, kg/m³, for the relief limit (F15).
pub const STRENGTH_ROCK: f64 = 150.0e6;
pub const STRENGTH_ICE: f64 = 5.0e6;
pub const DENSITY_CRUST_ROCK: f64 = 2_800.0;
pub const DENSITY_CRUST_ICE: f64 = 920.0;
/// Relief coefficient c in h_max = c strength / (rho g), set so that Earth gives 9.5 km,
/// Mars 25 km and Venus 10.5 km with the rock strength above.
pub const RELIEF_COEFFICIENT: f64 = 1.74;
/// Volumetric thermal expansivity for global contraction, per K, and the interior cooling
/// model dT = 300 K (1 - exp(-age / 2 Gyr)). Mercury gives 5 to 7 km. Provisional, Phase 3.
pub const THERMAL_EXPANSIVITY: f64 = 3.0e-5;
pub const CONTRACTION_DT_MAX: f64 = 300.0;
pub const CONTRACTION_TIMESCALE_GYR: f64 = 2.0;
/// Magma production (F17): fraction of the heat flux carried by magma, by regime, with a
/// tidal boost that reaches 0.9 at 1 W/m². Latent plus sensible heat per m³ of magma.
/// Provisional, Phase 3.
pub const MAGMA_EFFICIENCY_MOBILE: f64 = 0.15;
pub const MAGMA_EFFICIENCY_EPISODIC: f64 = 0.05;
pub const MAGMA_EFFICIENCY_STAGNANT: f64 = 0.01;
pub const MAGMA_EFFICIENCY_ICE: f64 = 0.02;
pub const MAGMA_EFFICIENCY_TIDAL_MAX: f64 = 0.9;
pub const MAGMA_TIDAL_FLUX_SCALE: f64 = 1.0;
pub const MAGMA_HEAT_PER_M3_ROCK: f64 = 2_800.0 * (4.0e5 + 1_200.0 * 1_000.0);
pub const MAGMA_HEAT_PER_M3_ICE: f64 = 920.0 * (3.34e5 + 2_100.0 * 100.0);
/// Melting point of water ice and the ammonia–water eutectic, K.
pub const T_MELT_WATER: f64 = 273.15;
pub const T_EUTECTIC_AMMONIA_WATER: f64 = 176.0;
/// Regime thresholds (F19). Provisional, Phase 3.
pub const MOBILE_LID_MIN_FLUX: f64 = 0.05;
pub const MOBILE_LID_MIN_RADIUS: f64 = 4.0e6;
pub const EPISODIC_LID_MIN_FLUX: f64 = 0.04;
pub const EPISODIC_LID_MIN_RADIUS: f64 = 5.0e6;

// ---------------------------------------------------------------- atmosphere and climate

/// Grey greenhouse optical depth per species: tau = a u^b with u the species column mass in
/// kg/m². CO2 fit to Venus (tau 148 at 1.0e6 kg/m²) and Mars (0.076 at 156 kg/m²); water
/// fit to Earth's 25 kg/m² column giving the balance of Earth's tau 0.85; methane and the
/// N2 collision-induced term fit to Titan (0.30 and 0.37). Provisional, Phase 5.
pub const TAU_CO2: (f64, f64) = (9.7e-4, 0.864);
pub const TAU_H2O: (f64, f64) = (0.169, 0.5);
pub const TAU_CH4: (f64, f64) = (5.3e-3, 0.5);
pub const TAU_N2_CIA: (f64, f64) = (1.09e-8, 1.5);
/// Water vapour column model: 25 kg/m² at 288 K and 70% relative humidity, scaled by
/// Clausius–Clapeyron with L / R_v = 5420 K and capped by the inventory.
pub const WATER_COLUMN_REF: f64 = 25.0;
pub const WATER_COLUMN_REF_T: f64 = 288.0;
pub const WATER_COLUMN_REF_RH: f64 = 0.7;
pub const WATER_L_OVER_RV: f64 = 5_420.0;
pub const WATER_COLUMN_MAX: f64 = 1.0e4;
/// Bond albedo guesses by class before the climate stage exists. Provisional, Phase 5.
pub const ALBEDO_AIRLESS_ROCK: f64 = 0.10;
pub const ALBEDO_AIRLESS_ROCK_ICE: f64 = 0.45;
pub const ALBEDO_AIRLESS_ICE: f64 = 0.70;
pub const ALBEDO_THIN_ATMOSPHERE: f64 = 0.25;
pub const ALBEDO_EARTHLIKE: f64 = 0.30;
pub const ALBEDO_THICK_ATMOSPHERE: f64 = 0.75;
pub const THICK_ATMOSPHERE_PRESSURE: f64 = 10.0e5;
pub const THIN_ATMOSPHERE_PRESSURE: f64 = 0.1e5;
/// Jeans parameter below which a species is flagged as escaping. The rule of thumb that a
/// gas survives when the escape speed exceeds six times its most probable thermal speed is
/// lambda = G M m / (k T R) above 36, evaluated at the exobase temperature.
pub const JEANS_RETENTION_THRESHOLD: f64 = 36.0;
/// Exobase temperature model: T_exo = T_s + EXOBASE_HEATING x xuv^0.3 x (1 - 0.8 f_CO2).
/// Earth gives about 1000 K, Mars about 250 K. Provisional, Phase 5.
pub const EXOBASE_HEATING: f64 = 700.0;
pub const EXOBASE_XUV_EXPONENT: f64 = 0.3;
pub const EXOBASE_CO2_COOLING: f64 = 0.8;
/// Crater cut-off (F24): D_min = c (P / g)^0.6 with Earth anchored at 150 m.
pub const CUTOFF_EXPONENT: f64 = 0.6;
pub const CUTOFF_EARTH_DIAMETER: f64 = 150.0;
/// Held–Hou (F32): tropopause height as a multiple of the scale height and the fractional
/// equator-to-pole radiative contrast. Earth gives about 32°.
pub const TROPOPAUSE_SCALE_HEIGHTS: f64 = 2.0;
pub const HELD_HOU_DELTA_H: f64 = 0.25;
/// Williams and Kasting transport coefficient (F33): Earth value and the cap found in the
/// Phase 0 spike to reproduce published locked-world day–night contrasts.
pub const D_EARTH: f64 = 0.58;
pub const D_CAP_MULTIPLE: f64 = 10.0;
/// Specific heat capacities at constant pressure, J/kg/K, by gas, for the lapse rate.
pub const CP_BY_GAS: [(&str, f64); 10] = [
    ("N2", 1_040.0),
    ("O2", 918.0),
    ("CO2", 844.0),
    ("CH4", 2_220.0),
    ("H2", 14_300.0),
    ("He", 5_193.0),
    ("Ar", 520.0),
    ("SO2", 622.0),
    ("NH3", 2_090.0),
    ("H2O", 1_996.0),
];
/// Molar masses, kg/mol.
pub const MOLAR_MASS_BY_GAS: [(&str, f64); 10] = [
    ("N2", 0.028_014),
    ("O2", 0.031_998),
    ("CO2", 0.044_01),
    ("CH4", 0.016_04),
    ("H2", 0.002_016),
    ("He", 0.004_003),
    ("Ar", 0.039_95),
    ("SO2", 0.064_07),
    ("NH3", 0.017_03),
    ("H2O", 0.018_015),
];
/// Moist lapse-rate factor applied when the working liquid is present.
pub const MOIST_LAPSE_FACTOR: f64 = 0.65;
/// Obliquity above which the annual insolation maximum moves to the poles, degrees.
pub const HIGH_OBLIQUITY_DEG: f64 = 54.0;

// ---------------------------------------------------------------- volatiles

/// Triple point temperature (K), triple point pressure (Pa), latent heat of vaporisation
/// (J/kg), latent heat of sublimation (J/kg), specific gas constant (J/kg/K), liquid
/// density (kg/m³), for each volatile species.
pub struct Volatile {
    pub name: &'static str,
    pub t_triple: f64,
    pub p_triple: f64,
    pub l_vap: f64,
    pub l_sub: f64,
    pub r_specific: f64,
    pub liquid_density: f64,
}

pub const VOLATILES: [Volatile; 7] = [
    // Water's latent heat of vaporisation is the effective mean over 273 to 373 K (2.50e6 at
    // the triple point, 2.26e6 at the normal boiling point), which puts the boiling point at
    // 101325 Pa within 0.3 K.
    Volatile {
        name: "H2O",
        t_triple: 273.16,
        p_triple: 611.657,
        l_vap: 2.40e6,
        l_sub: 2.83e6,
        r_specific: 461.5,
        liquid_density: 1_000.0,
    },
    Volatile {
        name: "CO2",
        t_triple: 216.55,
        p_triple: 5.185e5,
        l_vap: 3.5e5,
        l_sub: 5.71e5,
        r_specific: 188.9,
        liquid_density: 1_100.0,
    },
    Volatile {
        name: "CH4",
        t_triple: 90.69,
        p_triple: 1.17e4,
        l_vap: 5.1e5,
        l_sub: 6.0e5,
        r_specific: 518.3,
        liquid_density: 423.0,
    },
    Volatile {
        name: "N2",
        t_triple: 63.15,
        p_triple: 1.2523e4,
        l_vap: 1.99e5,
        l_sub: 2.25e5,
        r_specific: 296.8,
        liquid_density: 808.0,
    },
    Volatile {
        name: "NH3",
        t_triple: 195.4,
        p_triple: 6.06e3,
        l_vap: 1.37e6,
        l_sub: 1.76e6,
        r_specific: 488.2,
        liquid_density: 682.0,
    },
    Volatile {
        name: "SO2",
        t_triple: 197.6,
        p_triple: 1.67e3,
        l_vap: 3.89e5,
        l_sub: 4.6e5,
        r_specific: 129.8,
        liquid_density: 1_460.0,
    },
    Volatile {
        name: "C2H6",
        t_triple: 90.4,
        p_triple: 1.1,
        l_vap: 4.89e5,
        l_sub: 5.8e5,
        r_specific: 276.5,
        liquid_density: 544.0,
    },
];

// ---------------------------------------------------------------- cratering (scaling only)

/// Simple-to-complex transition (F30): lunar anchor, lunar gravity, ice factor.
pub const TRANSITION_LUNAR_KM: f64 = 16.5;
pub const GRAVITY_MOON: f64 = 1.62;
pub const TRANSITION_ICE_FACTOR: f64 = 0.15;
pub const PEAK_RING_MULTIPLE: f64 = 8.5;
pub const MULTIRING_MULTIPLE: f64 = 20.0;
/// Lunar topographic diffusivity, m²/Myr (Fassett and Thomson 2014), and the atmosphere
/// multiplier above 1 kPa. Provisional, Phase 2.
pub const DIFFUSIVITY_LUNAR: f64 = 5.5;
pub const DIFFUSIVITY_ATMOSPHERE_MULTIPLE: f64 = 10.0;
pub const DIFFUSIVITY_ATMOSPHERE_PRESSURE: f64 = 1.0e3;
/// Ray fading timescale at 1 AU, Gyr, airless; and under an atmosphere. Provisional, Phase 2.
pub const RAY_FADING_GYR_1AU: f64 = 1.0;
pub const RAY_FADING_GYR_ATMOSPHERE: f64 = 0.05;
/// Impactor densities by region, kg/m³.
pub const IMPACTOR_DENSITY_INNER: f64 = 2_700.0;
pub const IMPACTOR_DENSITY_OUTER: f64 = 600.0;
/// Heliocentric distance beyond which the outer-system impactor population applies, AU.
pub const OUTER_SYSTEM_AU: f64 = 4.0;
/// Impact velocity model (F29): v_inf = max(0.5 (1 AU / d)^0.5, 0.4) x local orbital speed.
pub const VINF_FACTOR_INNER: f64 = 0.5;
pub const VINF_FACTOR_MIN: f64 = 0.4;
/// Neukum chronology function coefficients: N(>=1 km) = A (exp(B T) - 1) + C T per km².
pub const NEUKUM_A: f64 = 5.44e-14;
pub const NEUKUM_B: f64 = 6.93;
pub const NEUKUM_C: f64 = 8.38e-4;

// ---------------------------------------------------------------- resolution by quality

pub const LOOP_FACE_SIZE: [u32; 3] = [128, 256, 512];
pub const FINISHING_FACE_SIZE: [u32; 3] = [512, 1024, 2048];
