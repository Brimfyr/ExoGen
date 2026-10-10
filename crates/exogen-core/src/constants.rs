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
/// The relief limit through time is scaled by the elastic thickness relative to the
/// present one, to no less than this fraction: a thin early lithosphere carries less.
/// Provisional until a flexural model in Phase 6.
pub const RELIEF_SUPPORT_MIN: f64 = 0.15;
/// Volumetric thermal expansivity for global contraction, per K, and the interior cooling
/// model dT = (50 K + 220 K i²) (1 - exp(-age / 2 Gyr)) with i the iron index from the
/// uncompressed density, 0 at 3.3 g/cm³ and 1 at the 5.15 the composition model gives Mercury: every body shrinks a
/// little from mantle cooling, and a large core, cooling and solidifying, is what shrinks
/// Mercury by 5 to 7 km (Byrne et al. 2014) while the Moon and Mars lose under a kilometre
/// and a few kilometres (Watters et al. 2010; Nahm and Schultz 2011). Provisional, Phase 3.
pub const THERMAL_EXPANSIVITY: f64 = 3.0e-5;
pub const CONTRACTION_DT_MANTLE: f64 = 50.0;
pub const CONTRACTION_DT_CORE: f64 = 250.0;
/// Radius above which a body accreted hot throughout. Smaller bodies started cool inside
/// and heated for a long time, offsetting the cooling of their outer shells, so their net
/// mantle contraction scales down with (R / this)²: the Moon's few hundred metres against
/// the kilometre thermal models give for a hot start (Watters et al. 2010).
pub const CONTRACTION_HOT_START_RADIUS: f64 = 3.39e6;
pub const CONTRACTION_DENSITY_ROCK: f64 = 3_300.0;
pub const CONTRACTION_DENSITY_IRON_RICH: f64 = 5_150.0;
pub const CONTRACTION_TIMESCALE_GYR: f64 = 2.0;
/// Time after formation before net contraction begins: early radiogenic heating and the
/// magma-ocean aftermath offset secular cooling at first (Tosi et al. 2013 for Mercury;
/// the Hesperian onset of wrinkle ridges on Mars), so the scarps that survive post-date
/// the heavy bombardment.
pub const CONTRACTION_ONSET_GYR: f64 = 0.7;
/// Lobate scarps (Byrne et al. 2014; Watters et al. 2009): lengths from a power law with
/// this cumulative exponent between 10 km and this fraction of the radius; height from
/// length as h_km = 0.03 L_km^0.667 (Enterprise Rupes: 1,000 km, 3 km); thrust dip for
/// the shortening; face slope; back slope length in heights; total heading change along
/// a trace; and the shortest scarp stamped in the loop, in loop cells. Provisional.
pub const SCARP_LENGTH_EXPONENT: f64 = 1.5;
pub const SCARP_MIN_LENGTH: f64 = 10_000.0;
pub const SCARP_MAX_LENGTH_FRACTION: f64 = 0.4;
pub const SCARP_HEIGHT_COEFF_KM: f64 = 0.03;
pub const SCARP_HEIGHT_EXPONENT: f64 = 0.667;
pub const SCARP_DIP_DEG: f64 = 25.0;
pub const SCARP_FACE_SLOPE_DEG: f64 = 12.0;
pub const SCARP_BACK_SLOPE_FACTOR: f64 = 20.0;
pub const SCARP_TURN_DEG: f64 = 40.0;
pub const SCARP_LOOP_MIN_CELLS: f64 = 2.0;
/// Lobes along a trace: harmonics of heading wobble and the amplitude of the first.
pub const SCARP_LOBE_HARMONICS: usize = 3;
pub const SCARP_LOBE_DEG: f64 = 18.0;
/// Longest scarp as a function of strain: the reference length at the reference strain
/// (Mercury's thousand-kilometre scarps at half a per cent), scaling with the strain to
/// this power as faulting deepens and segments link.
pub const SCARP_REFERENCE_MAX_LENGTH: f64 = 1.0e6;
pub const SCARP_REFERENCE_STRAIN: f64 = 0.005;
pub const SCARP_LENGTH_STRAIN_EXPONENT: f64 = 2.0;
/// Flood volcanism (PLAN.md 4.4). The extrusive share of magma production (intrusive to
/// extrusive 5:1 to 10:1, Crisp 1984); the lithosphere thickness over which ascent fails
/// (the Moon's maria waned as its lithosphere passed a few hundred kilometres, Head and
/// Wilson 1992); the contraction strain that closes the pathways (Mercury's volcanism
/// ended as it shrank, Byrne et al. 2016); vent weighting scales for crust thickness and
/// height; and the range of single eruption volumes, m³ (mare flow units of hundreds to
/// thousands of km³). Provisional.
pub const VOLCANISM_EXTRUSIVE_FRACTION: f64 = 0.10;
/// Lithosphere thickness past which magma no longer reaches the surface, and the width of
/// that cutoff, m: the Moon's maria waned once its lithosphere passed some 500 km.
pub const VOLCANISM_ASCENT_LITHOSPHERE: f64 = 400e3;
pub const VOLCANISM_ASCENT_WIDTH: f64 = 30e3;
pub const VOLCANISM_SHUTOFF_STRAIN: f64 = 1.2e-3;
pub const VOLCANISM_VENT_CRUST_SCALE: f64 = 10e3;
pub const VOLCANISM_VENT_HEIGHT_SCALE: f64 = 2e3;
pub const VOLCANISM_EVENT_VOLUME: (f64, f64) = (1e11, 1e13);
/// Cells one flood's lake may take in, and the thickness of a flow, m: a lake spreads no
/// thinner than this on average, and a flow drapes bumps that stand less than this above
/// its level (flood basalt flows are tens of metres thick).
pub const FLOOD_MAX_CELLS: usize = 40_000;
pub const FLOOD_SHEET_THICKNESS: f64 = 50.0;
/// Crust lost under a crater as a fraction of its excavation depth, after the mantle's
/// rise and the melt sheet, and the floor the crust cannot thin below, m. The lunar
/// basins thinned the crust to under 10 km (Wieczorek et al. 2013), which is where the
/// maria are.
pub const BASIN_CRUST_THINNING: f64 = 0.5;
pub const CRUST_THICKNESS_MIN: f64 = 3e3;
/// Thinner lava than this is not recorded, m.
pub const FLOOD_MIN_THICKNESS: f64 = 1.0;
/// How far from its vent a flow can travel, m: the longest flood-basalt flows run a
/// thousand kilometres or so (Tharsis, the Deccan), not around the body.
pub const FLOOD_MAX_REACH: f64 = 1.5e6;
/// Flow fronts at the finishing resolution: the smeared coarse footprint is kept where
/// it exceeds this threshold (below a half, so a single flooded cell keeps its area and
/// fronts sit a little outside the coarse edge), with lobes of noise this many loop cells
/// long moving the threshold by up to this much.
pub const FLOOD_MARGIN_THRESHOLD: f64 = 0.35;
pub const FLOOD_MARGIN_LOBE_CELLS: f64 = 1.5;
pub const FLOOD_MARGIN_WOBBLE: f64 = 0.2;
/// Subsidence of a lava fill as a fraction of its thickness, from compaction of the pile
/// and the sag of the floor under it, so a surface over a buried crater dips over the
/// deep fill and shows the rim as a ring (ghost craters).
pub const LAVA_SUBSIDENCE: f64 = 0.05;
/// Wrinkle ridges on plains (PLAN.md 4.4): ridge length per plains area, 1/m (lunar
/// maria carry about a kilometre of ridge per forty square kilometres); the range of
/// one ridge's length, m; arch height and width, crest height and width, m; how far the
/// crest wanders within the arch and over what length; the floor of the seed weight
/// relative to the mean thickness gradient; how strongly a trace follows the contours;
/// the heading wobble; and the delay after the last flow before the ridges stand.
pub const RIDGE_LENGTH_PER_AREA: f64 = 1.0 / 40e3;
pub const RIDGE_LENGTH: (f64, f64) = (20e3, 300e3);
pub const RIDGE_HEIGHT: f64 = 200.0;
pub const RIDGE_ARCH_WIDTH: f64 = 6e3;
pub const RIDGE_CREST_HEIGHT: f64 = 80.0;
pub const RIDGE_CREST_WIDTH: f64 = 1.2e3;
pub const RIDGE_CREST_WANDER: f64 = 0.6;
pub const RIDGE_CREST_WANDER_LENGTH: f64 = 25e3;
pub const RIDGE_SEED_FLOOR: f64 = 0.3;
pub const RIDGE_CONTOUR_FOLLOWING: f64 = 0.35;
pub const RIDGE_WOBBLE_DEG: f64 = 25.0;
pub const RIDGE_WOBBLE_LENGTH: f64 = 60e3;
pub const RIDGE_DELAY_GYR: f64 = 0.2;
/// Extra weight of a vent on a volcanic rise, per unit of dome height over the rise's
/// full height, and the dome's share in the load field the ridges follow.
pub const PLUME_VENT_BOOST: f64 = 20.0;
pub const RIDGE_LOAD_WEIGHT: f64 = 0.5;
/// Plume loading (PLAN.md 4.3). A stagnant lid at least this big raises one plume per
/// this much surface heat flux, up to three (Mars: one, Tharsis); the rise spans this
/// fraction of the radius and stands this fraction of the strength-limited relief
/// (Tharsis: 5,000 km across, 7 km high); it grows between these times.
pub const PLUME_MIN_RADIUS: f64 = 3.0e6;
pub const PLUME_FLUX_PER_PLUME: f64 = 0.02;
pub const PLUME_RISE_DIAMETER_FRACTION: f64 = 1.5;
/// A thin, hot episodic lid supports narrower rises (Beta Regio: 2,500 km on Venus).
pub const PLUME_RISE_DIAMETER_FRACTION_EPISODIC: f64 = 0.4;
pub const RISE_HEIGHT_FRACTION: f64 = 0.28;
pub const RISE_START_GYR: f64 = 0.3;
pub const RISE_END_GYR: f64 = 1.5;
/// Shields on a rise: how many, the largest as a fraction of the rise diameter (Olympus
/// Mons: 600 km on a 5,000 km rise), the flank slope, the flank shape, where on the rise
/// they sit (fraction of its radius), when they start, their growth timescale, and the
/// thickness of flows that renews a flank's surface. Calderas as fractions of the
/// shield's radius and height, with a cap on depth and the wall width as a fraction of
/// the caldera radius. Flank texture: radial flow ridges of this relative amplitude and
/// length.
pub const SHIELDS_PER_RISE: f64 = 5.0;
pub const SHIELD_MAX_DIAMETER_FRACTION: f64 = 0.12;
pub const SHIELD_SLOPE_DEG: f64 = 5.0;
pub const SHIELD_FLANK_EXPONENT: f64 = 1.6;
pub const SHIELD_FLANK_POWER: f64 = 0.9;
pub const SHIELD_PLACEMENT: (f64, f64) = (0.1, 0.6);
pub const SHIELD_START_GYR: (f64, f64) = (1.0, 2.5);
pub const SHIELD_GROWTH_GYR: f64 = 1.0;
pub const SHIELD_RENEWAL_THICKNESS: f64 = 5.0;
pub const CALDERA_RADIUS_FRACTION: f64 = 0.1;
pub const CALDERA_DEPTH_FRACTION: f64 = 0.12;
pub const CALDERA_MAX_DEPTH: f64 = 3e3;
pub const CALDERA_WALL_WIDTH: f64 = 0.15;
pub const SHIELD_FLOW_TEXTURE: f64 = 0.3;
pub const SHIELD_FLOW_TEXTURE_LENGTH: f64 = 15e3;
/// The canyon (Valles Marineris: 4,000 km long, 200 km wide, 7 km deep) on a rise at
/// least this wide: its length as a fraction of the rise diameter, where it starts
/// (fraction of the rise radius from the centre), width, depth, growth interval, and the
/// carving that renews the floor's surface in a step. Graben on the flanks: how many
/// per rise, where they start, their lengths, widths and depth. Trough traces wobble by
/// this much over this length; their walls take this fraction of the width; walls are
/// scalloped by this fraction over this length; the width pulses by this fraction over
/// this length (chasmata); and the floor tapers to the ends with this power.
pub const CANYON_MIN_RISE_DIAMETER: f64 = 2.0e6;
pub const CANYON_LENGTH_FRACTION: f64 = 0.8;
pub const CANYON_START_FRACTION: f64 = 0.6;
/// Width as a fraction of the rise diameter and depth as a fraction of the rise height
/// (Valles Marineris: 200 km wide and 7 km deep on a 5,000 km, 7 km rise; Devana
/// Chasma on Venus: 100 km wide and 2 to 3 km deep on a smaller, lower rise).
pub const CANYON_WIDTH_FRACTION: f64 = 0.04;
pub const CANYON_DEPTH_FRACTION: f64 = 1.0;
pub const CANYON_START_GYR: f64 = 1.0;
pub const CANYON_END_GYR: f64 = 2.0;
pub const CANYON_RENEWAL_THICKNESS: f64 = 50.0;
pub const GRABEN_PER_RISE: f64 = 60.0;
pub const GRABEN_START_FRACTION: (f64, f64) = (0.35, 1.0);
pub const GRABEN_LENGTH: (f64, f64) = (200e3, 1_500e3);
pub const GRABEN_WIDTH: (f64, f64) = (8e3, 30e3);
pub const GRABEN_DEPTH: f64 = 400.0;
pub const TROUGH_WOBBLE_DEG: f64 = 12.0;
pub const TROUGH_WOBBLE_LENGTH: f64 = 400e3;
pub const TROUGH_WALL_FRACTION: f64 = 0.05;
pub const TROUGH_SCALLOP: f64 = 0.12;
pub const TROUGH_SCALLOP_LENGTH: f64 = 60e3;
pub const TROUGH_WIDTH_PULSE: f64 = 0.4;
pub const TROUGH_PULSE_LENGTH: f64 = 700e3;
pub const TROUGH_TAPER_POWER: f64 = 0.3;
/// The episodic lid (PLAN.md 4.3, Venus). The event fills every depression this far
/// towards its spill point and lays this veneer everywhere else on the plains, m; the
/// tessera are the highest-standing fraction of the surface by a height smoothed this
/// many passes, keep a surface this much older than the event, Gyr, and carry a fabric
/// of this relief, m, at this wavelength, m. After the event, eruptions carry this share
/// of the extrusive budget, the rest having gone into the event.
pub const TESSERA_FRACTION: f64 = 0.08;
pub const TESSERA_SMOOTHING_PASSES: usize = 20;
pub const TESSERA_AGE_LEAD_GYR: f64 = 0.5;
pub const TESSERA_RELIEF: f64 = 1_000.0;
pub const TESSERA_FABRIC_LENGTH: f64 = 30e3;
pub const EPISODIC_FILL_FRACTION: f64 = 0.7;
pub const EPISODIC_VENEER: f64 = 20.0;
pub const EPISODIC_POST_EVENT_FRACTION: f64 = 0.02;
/// Coronae on an episodic lid: how many per m² (Venus: some 500), their diameters, m,
/// with the draw skewed to small ones by this power, the share whose centre sags, the
/// central relief and rim height, m, the rim's position and width as fractions of the
/// radius, the moat as a fraction of the rim height, the rim's height wobble, and the
/// concentric fractures' relief and spacing. Pancake domes: how many per m² (Venus:
/// some 150), diameters, height and top roughness.
pub const CORONAE_PER_AREA: f64 = 1.1e-12;
pub const CORONA_DIAMETER: (f64, f64) = (100e3, 1_000e3);
pub const CORONA_SIZE_SKEW: f64 = 2.0;
pub const CORONA_SAG_FRACTION: f64 = 0.4;
pub const CORONA_CENTRE_RELIEF: f64 = 400.0;
pub const CORONA_RIM_HEIGHT: f64 = 800.0;
pub const CORONA_RIM_POSITION: f64 = 0.85;
pub const CORONA_RIM_WIDTH: f64 = 0.05;
pub const CORONA_MOAT_FRACTION: f64 = 0.3;
pub const CORONA_RIM_WOBBLE: f64 = 0.3;
pub const CORONA_FRACTURE_RELIEF: f64 = 0.15;
pub const CORONA_FRACTURE_SPACING: f64 = 0.08;
pub const DOMES_PER_AREA: f64 = 3.3e-13;
pub const DOME_DIAMETER: (f64, f64) = (20e3, 60e3);
pub const DOME_HEIGHT: f64 = 800.0;
pub const DOME_TOP_ROUGHNESS: f64 = 0.08;
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

// ---------------------------------------------------------------- cratering (population)

/// Neukum lunar production function (Neukum, Ivanov and Hartmann 2001): log10 N(>= D) per
/// km² for a 1 Gyr surface as a polynomial in log10 D (km), valid 0.01 to 300 km.
pub const NPF_LUNAR: [f64; 12] = [
    -3.0876, -3.557528, 0.781027, 1.021521, -0.156012, -0.444058, 0.019977, 0.086850, -0.005874,
    -0.006809, 8.25e-4, 5.54e-5,
];
/// Ivanov 2001 Mars production function and the Hartmann and Neukum 2001 Mars chronology.
/// Provisional until the Mars stage in Phase 5 checks them against the source.
pub const NPF_MARS: [f64; 12] = [
    -3.384, -3.197, 1.257, 0.7915, -0.4861, -0.3630, 0.1016, 6.756e-2, -1.181e-2, -4.753e-3,
    6.233e-4, 5.805e-5,
];
pub const MARS_CHRONOLOGY: (f64, f64, f64) = (2.68e-14, 6.93, 4.13e-4);
/// Diameter range of the production polynomials, km; continued with the edge slope outside.
pub const NPF_RANGE_KM: (f64, f64) = (0.01, 300.0);
/// The Neukum chronology counts time before present, and the present is this long after
/// the Moon formed. Rewriting the chronology in time since formation fixes the early flux
/// at formation (F31), so bodies of any age get their heavy bombardment when they form.
pub const NEUKUM_REFERENCE_AGE_GYR: f64 = 4.5;
/// The steady flux after the heavy bombardment comes from the collisionally evolved
/// asteroid belt, whose large bodies are depleted next to the leftover planetesimals
/// that dug the basins: above this crater diameter, km, its production falls steeper
/// than the lunar function by this extra cumulative slope, and no stray body makes a
/// crater larger than this. The Moon and Mars have no crater over a few hundred
/// kilometres younger than the basin era.
pub const STEADY_FLUX_KNEE_KM: f64 = 150.0;
pub const STEADY_FLUX_EXTRA_SLOPE: f64 = 2.0;
pub const STEADY_FLUX_MAX_KM: f64 = 600.0;
/// Surface solidification: time after formation when the crust is solid and rigid enough
/// to keep large relief, so the crater record begins. 0.3 Gyr for a Moon-sized rocky body
/// (the oldest dated lunar surfaces are 4.2 Gyr old), scaling with the square root of the
/// radius since larger bodies cool more slowly, and shorter for ice shells. Provisional:
/// the Phase 3 heat budget derives it.
pub const SOLIDIFICATION_REFERENCE_GYR: f64 = 0.3;
pub const SOLIDIFICATION_REFERENCE_RADIUS: f64 = 1.7374e6;
pub const SOLIDIFICATION_ICE_FACTOR: f64 = 0.7;
pub const SOLIDIFICATION_RANGE_GYR: (f64, f64) = (0.02, 1.0);
/// Geometric saturation N(>= D) = 1.54 D^-2 per km² (Hartmann 1984) and the fraction of it
/// at which crater populations reach equilibrium (below 5%).
pub const GEOMETRIC_SATURATION: f64 = 1.54;
pub const EMPIRICAL_SATURATION_FRACTION: f64 = 0.05;
/// Pike 1977 lunar morphometry, (coefficient, exponent) in km: d = a D^b.
pub const PIKE_SIMPLE_DEPTH: (f64, f64) = (0.196, 1.010);
pub const PIKE_SIMPLE_RIM: (f64, f64) = (0.036, 1.014);
pub const PIKE_COMPLEX_DEPTH: (f64, f64) = (1.044, 0.301);
pub const PIKE_COMPLEX_RIM: (f64, f64) = (0.236, 0.399);
/// Central peak height; Pike's 0.032 D^0.9 over-predicts Copernicus, so the coefficient is
/// fitted to it. Provisional.
pub const CENTRAL_PEAK_HEIGHT: (f64, f64) = (0.02, 0.9);
/// Flat floor diameter of complex craters, km (Pike 1977 form).
pub const FLOOR_DIAMETER: (f64, f64) = (0.19, 1.25);
/// Caps for the largest basins, metres at lunar gravity; scaled as (g_moon / g)^0.25.
pub const MAX_BASIN_DEPTH_LUNAR: f64 = 9_000.0;
pub const MAX_RIM_HEIGHT_LUNAR: f64 = 3_500.0;
/// Ejecta extent in crater radii and the radial thickness exponent (McGetchin et al. 1973
/// give -3 for the thickness decay beyond the rim).
pub const EJECTA_RADIUS_MULTIPLE: f64 = 3.0;
pub const EJECTA_EXPONENT: f64 = -3.0;
pub const BASIN_EJECTA_RADIUS_MULTIPLE: f64 = 4.0;
pub const BASIN_EJECTA_EXPONENT: f64 = -2.5;
/// Extent of the continuous ejecta that resets the surface age, in radii.
pub const CONTINUOUS_EJECTA_MULTIPLE: f64 = 2.0;
/// Ejecta burial: relief under ejecta is damped by exp(-k t / |relief|) for thickness t.
/// Depressions fill (ejecta ponds in lows) faster than highs are shaken down. Provisional.
pub const BURIAL_FACTOR_LOW: f64 = 1.0;
pub const BURIAL_FACTOR_HIGH: f64 = 0.5;
/// Oblique impacts (Gault and Wedekind 1978; Bottke et al. 2000; Herrick and Forsberg-Taylor
/// 2003). Impact angles follow sin 2θ. Cavities stay circular down to the angle below which
/// the elliptical fraction of craters forms, with the elongation rising toward grazing
/// incidence by at most this excess. Ejecta thickens downrange and thins uprange below
/// 45°, an uprange zone of avoidance opens below 30°, and side lobes grow below 15°.
pub const ELLIPTICITY_MAX_EXCESS: f64 = 0.8;
pub const OBLIQUE_ASYMMETRY_DEG: f64 = 45.0;
pub const OBLIQUE_DOWNRANGE_MAX: f64 = 0.7;
pub const OBLIQUE_FORBIDDEN_DEG: f64 = 30.0;
pub const OBLIQUE_FORBIDDEN_HALF_WIDTH_DEG: f64 = 60.0;
pub const OBLIQUE_FORBIDDEN_DEPTH: f64 = 0.85;
pub const OBLIQUE_BUTTERFLY_DEG: f64 = 15.0;
pub const OBLIQUE_BUTTERFLY_MAX: f64 = 1.0;
/// Secondaries strike at low angles along the radial from their primary.
pub const SECONDARY_IMPACT_ANGLE: f64 = 0.5;
/// Secondaries come in groups of this size; a group is a radial chain with this fraction,
/// spaced along the radial by this fraction of the primary radius, otherwise a cluster.
pub const SECONDARY_GROUP_SIZE: u32 = 8;
pub const SECONDARY_CHAIN_FRACTION: f64 = 0.6;
pub const SECONDARY_CHAIN_SPACING: f64 = 0.12;
/// Multi-ring basins: rings at successive 1/sqrt(2) of the rim radius (Hartmann and Kuiper
/// 1962), one more for each sqrt(2) of diameter above the onset, at most this many; the
/// inner rings are lower by this factor each. Orientale has three inside its rim.
pub const MAX_RINGS: usize = 3;
pub const RING_HEIGHT_FRACTION: f64 = 0.3;
pub const RING_HEIGHT_DECAY: f64 = 0.8;
/// Ring ridges: width as a fraction of the rim radius, their own scalloping (larger than
/// the rim's so each ring has its own shape), presence around the azimuth as mean and
/// swing (a few gaps where it drops to zero), and the subdued outer ring of a multi-ring
/// basin at sqrt(2) radii as a fraction of the rim height.
pub const RING_WIDTH_FRACTION: f64 = 0.04;
pub const RING_CRENULATION: f64 = 0.02;
pub const RING_PRESENCE: (f64, f64) = (0.85, 0.25);
pub const OUTER_RING_HEIGHT_FRACTION: f64 = 0.3;
/// Central peaks: basal radius as a fraction of the crater radius (peak diameters are a
/// fifth to a quarter of the crater's; Hale and Head 1979), the Gaussian width of a massif
/// as a fraction of that, the most massifs in a peak complex, and their spread as a
/// fraction of the peak radius.
pub const PEAK_RADIUS_FRACTION: f64 = 0.22;
pub const PEAK_SIGMA_FRACTION: f64 = 0.5;
pub const PEAK_MAX_MASSIFS: usize = 3;
pub const PEAK_SPREAD: f64 = 0.5;
/// Planar texture lattice: cells per crater radius on the tangent plane, and its extent in
/// radii; it serves floors, melt and the brightness mottling of fresh ejecta.
pub const FLOOR_TEXTURE_CELLS: u32 = 10;
pub const PLANAR_TEXTURE_RADII: f64 = 2.2;
/// Smallest crater, in finishing cells of diameter, whose fresh ejecta is mottled and
/// whose continuous-ejecta edge is ragged in the ray field; smaller ones are smooth.
pub const RAY_TEXTURE_MIN_CELLS: f64 = 24.0;
/// Rim and ring scalloping as a fraction of the radius, and ejecta, basin-sculpture and
/// floor texture amplitudes. Provisional, to be tuned against imagery in Phase 7.
pub const RIM_CRENULATION_SIMPLE: f64 = 0.015;
pub const RIM_CRENULATION_COMPLEX: f64 = 0.035;
pub const EJECTA_TEXTURE: f64 = 0.5;
pub const BASIN_SCULPTURE: f64 = 0.25;
pub const FLOOR_TEXTURE: f64 = 0.04;
/// Texture lattice: cells around the azimuth and per e-fold of radius, per octave; the
/// basin sculpture lattice is stretched radially. Craters with fewer footprint cells than
/// this are stamped without texture.
pub const TEXTURE_AZIMUTH_CELLS: u32 = 24;
pub const TEXTURE_RADIAL_CELLS: u32 = 6;
pub const SCULPTURE_AZIMUTH_CELLS: u32 = 144;
pub const SCULPTURE_RADIAL_CELLS: u32 = 3;
pub const TEXTURE_MIN_CELLS: usize = 1800;
/// The ejecta texture fades out between these radii, in crater radii, where the ejecta is
/// too thin for it to show.
pub const TEXTURE_FADE_RADII: (f64, f64) = (1.3, 2.0);
/// Deep bedrock: what craters excavate from below the surface layer. A crater reaches it
/// when its excavation depth, this fraction of its diameter (Melosh 1989), exceeds the
/// layer thickness: the anorthositic upper crust of a Moon-like body, the depth of
/// Mercury's low-reflectance material, or the dirty surface layer of an ice shell.
/// Provisional: Phase 4 for ice, Phase 7 for the material library.
pub const EXCAVATION_DEPTH_FRACTION: f64 = 0.1;
pub const DEEP_BEDROCK_DEPTH_ROCK: f64 = 25_000.0;
pub const DEEP_BEDROCK_DEPTH_IRON_RICH: f64 = 10_000.0;
pub const DEEP_BEDROCK_DEPTH_ICE: f64 = 500.0;
/// Brightness mottling of fresh ejecta in the ray field, and the lobing of the continuous
/// ejecta edge as a fraction of its radius. The lobing is off: at map resolution a few
/// harmonics read as petals, and the real lumpy outline belongs to the detail pass.
pub const EJECTA_BRIGHTNESS_TEXTURE: f64 = 0.25;
pub const EJECTA_EDGE_RAGGEDNESS: f64 = 0.0;
/// The surface-age reset fades out over this outer fraction of the continuous ejecta zone,
/// where the ejecta is thin and renews the surface only partly.
pub const EJECTA_EDGE_FADE: f64 = 0.08;
/// A loop crater counts as recognisable while its rim-to-floor relief on the coarse surface
/// is above this fraction of the fresh relief.
pub const VISIBLE_RELIEF_FRACTION: f64 = 0.15;
/// Impact melt: ponds on the floors of complex craters and basins, placed where a coarse
/// planar lattice of this many cells per crater radius exceeds this threshold, about a
/// fifth of the floor.
pub const MELT_POND_CELLS: u32 = 4;
pub const MELT_POND_THRESHOLD: f64 = 0.4;
/// Width of the soft edge of a melt pond in texture units.
pub const MELT_EDGE_WIDTH: f64 = 0.3;
/// Deep-sourced ejecta is mixed with local material by ballistic sedimentation (Oberbeck
/// 1975): the deposit is all primary ejecta out to this many radii and all local material
/// at the continuous ejecta edge.
pub const EJECTA_PURE_RADII: f64 = 1.1;
/// Rays: length range in radii, count range, and the smallest crater (in finishing cells)
/// that gets rays.
pub const RAY_LENGTH_MULTIPLE: (f64, f64) = (6.0, 15.0);
pub const RAYS_PER_CRATER: (u32, u32) = (8, 24);
pub const RAY_MIN_CELLS: f64 = 6.0;
/// Rays of craters larger than this reference (Tycho, 85 km) are shortened in radii by the
/// square root of the diameter ratio, since ejecta speed grows slowly with crater size.
pub const RAY_REFERENCE_DIAMETER: f64 = 85_000.0;
/// Linear ray width range as a fraction of the crater radius.
pub const RAY_WIDTH_FRACTION: (f64, f64) = (0.10, 0.30);
/// Secondaries: count per 100 km of primary diameter at a secondary fraction of 0.2, size
/// range as a fraction of the primary, distance range in primary radii, and the depth
/// factor relative to a primary of the same size.
pub const SECONDARIES_PER_100KM: f64 = 400.0;
pub const SECONDARY_SIZE_FRACTION: (f64, f64) = (0.005, 0.04);
pub const SECONDARY_DISTANCE_MULTIPLE: (f64, f64) = (2.5, 8.0);
pub const SECONDARY_DEPTH_FACTOR: f64 = 0.5;
/// Split between simulated and replayed craters, in loop cells; the smallest replayed
/// crater in finishing cells; and the size bins per octave of the replay.
pub const LOOP_THRESHOLD_CELLS: f64 = 8.0;
pub const FINISHING_MIN_CELLS: f64 = 2.0;
pub const BINS_PER_OCTAVE: u32 = 4;
/// Largest basin radius as a fraction of the body radius (South Pole–Aitken is 0.72).
pub const MAX_CRATER_RADIUS_FRACTION: f64 = 0.75;
/// Crater retention ages by regime, Gyr, until the erosion stages exist. Provisional.
pub const RETENTION_MOBILE_LID_GYR: f64 = 0.1;
pub const RETENTION_EPISODIC_GYR: f64 = 0.5;
pub const RETENTION_THIN_ATMOSPHERE_GYR: f64 = 3.5;
pub const RETENTION_ICE_OCEAN_GYR: f64 = 0.1;
/// Tidal flux above which an ice shell over an ocean is treated as actively resurfacing
/// (Europa 0.075 W/m², Enceladus 0.02) rather than ancient (Callisto, Ganymede).
pub const ICE_RESURFACING_TIDAL_FLUX: f64 = 0.01;
/// Burial depth that erases the replayed craters, metres, against the volcanic
/// resurfacing rate, applied only above the Io-regime tidal flux where volcanism
/// resurfaces the whole body.
pub const RETENTION_BURIAL_DEPTH_M: f64 = 10.0;
pub const IO_REGIME_TIDAL_FLUX: f64 = 0.5;
/// Space weathering timescale at 1 AU, Gyr, and the pressure above which it is off.
pub const SPACE_WEATHERING_GYR_1AU: f64 = 0.3;
pub const SPACE_WEATHERING_MAX_PRESSURE: f64 = 100.0;
/// Cap on the optical maturation timescale: far from the Sun, radiolysis, sputtering and
/// lag accumulation still darken old surfaces relative to fresh ejecta within a few Gyr.
pub const SPACE_WEATHERING_MAX_GYR: f64 = 4.0;
/// Lunar topographic diffusivity at the present impact flux, m²/Myr (Fassett and Thomson
/// 2014, whose fit spans the last 3 Gyr, over which the flux varied little). The flux
/// history scales it back in time, so a crater's degradation follows the impacts that
/// came after it. The atmosphere multiplier applies above 1 kPa. Provisional, Phase 5.
pub const DIFFUSIVITY_LUNAR: f64 = 5.5;
pub const DIFFUSIVITY_ATMOSPHERE_MULTIPLE: f64 = 10.0;
pub const DIFFUSIVITY_ATMOSPHERE_PRESSURE: f64 = 1.0e3;
/// Ray fading timescale at 1 AU, Gyr, airless; and under an atmosphere. Provisional, Phase 2.
pub const RAY_FADING_GYR_1AU: f64 = 1.0;
pub const RAY_FADING_GYR_ATMOSPHERE: f64 = 0.05;
/// Cap on the ray lifetime: far from the Sun, sputtering and regolith gardening still erase
/// rays within about this time (rayed craters are rare on the Galilean satellites).
pub const RAY_FADING_MAX_GYR: f64 = 1.5;
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
