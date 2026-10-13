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
pub const SECONDS_PER_GYR: f64 = 3.155_76e16;
pub const WATER_MOLECULE_MASS: f64 = 2.991e-26;
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
/// Magmatic provinces: the mantle is not uniformly fertile. The cumulate pile a magma
/// ocean leaves concentrates the heat-producing elements in a few regions, and the
/// eruptions favour them: the Moon's Procellarum KREEP Terrane holds most of the maria,
/// Mercury's northern plains most of its smooth plains. One province per this much
/// radius, their Gaussian angular radius, and the weight of the rest of the body.
pub const MAGMA_PROVINCE_RADIUS_PER: f64 = 2.0e6;
pub const MAGMA_PROVINCE_SIGMA_RAD: f64 = 0.7;
pub const MAGMA_PROVINCE_FLOOR: f64 = 0.03;
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
/// Lava thinner than this is mixed with its substrate by gardening and shows it through:
/// a flow's coverage of the surface is its thickness over this, m. Gardening mixes to
/// about the regolith's depth, 4 to 5 m on the maria (Oberbeck and Quaide 1968;
/// Shkuratov and Bondarenko 2001), so only lava a few metres thin shows what is under it.
pub const LAVA_COVER_DEPTH: f64 = 6.0;
/// How far from its vent a flow can travel, m: the longest flood-basalt flows run a
/// thousand kilometres or so (Tharsis, the Deccan), not around the body.
pub const FLOOD_MAX_REACH: f64 = 1.5e6;
/// Flow fronts at the finishing resolution: the coarse footprint's indicator, each cell
/// half its own and a sixteenth each of its neighbours', is kept where it exceeds this
/// threshold (below a half, so a single flooded cell keeps its area and fronts sit a
/// little outside the coarse edge, and the contour rounds the coarse cells' corners),
/// with three octaves of lobes of noise this many loop cells long moving the threshold
/// by up to this much, the unresolved shape of the front.
pub const FLOOD_MARGIN_THRESHOLD: f64 = 0.45;
pub const FLOOD_MARGIN_LOBE_CELLS: f64 = 2.0;
pub const FLOOD_MARGIN_WOBBLE: f64 = 0.3;
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
/// Rifting on a rise: the uplift stretches the lithosphere, and the extension across
/// the rift zone comes to a few times the uplift (Valles Marineris some 30 km on
/// Tharsis's 8 km, Andrews-Hanna 2012; Devana Chasma 10 to 20 km on Beta Regio's 4 km).
/// It is spent along a transect round the rise as a population of normal-fault troughs:
/// lengths from N(>L) ∝ L^−2 between a lower bound and 1.5 rise diameters, depth a
/// thousandth and a half of the length with log-normal scatter, capped at the uplift;
/// width thirty depths; heave 1.2 depths, two faults dipping 60°. Half the faults join an
/// existing one as an en-echelon neighbour, a width or a few to its side, so the big
/// ones gather parallel chasmata and the small ones form swarms. Troughs wider than two
/// coarse cells are carved through the loop over a gigayear; the rest are cut at the
/// finish. On a stagnant lid the rifting follows the dome; on an episodic lid it
/// follows the resurfacing event.
pub const RIFT_EXTENSION_PER_UPLIFT: f64 = 4.0;
pub const RIFT_TRANSECT_FRACTION: f64 = 0.7;
pub const RIFT_LENGTH_SLOPE: f64 = 2.0;
pub const RIFT_MIN_LENGTH: f64 = 50e3;
pub const RIFT_MAX_LENGTH_FRACTION: f64 = 1.5;
pub const RIFT_DEPTH_PER_LENGTH: f64 = 1.5e-3;
pub const RIFT_DEPTH_SCATTER: f64 = 0.4;
pub const RIFT_MIN_DEPTH: f64 = 50.0;
pub const RIFT_WIDTH_PER_DEPTH: f64 = 30.0;
pub const RIFT_MIN_WIDTH: f64 = 2e3;
pub const RIFT_HEAVE_PER_DEPTH: f64 = 1.2;
pub const RIFT_START_FRACTION: (f64, f64) = (0.3, 1.1);
pub const RIFT_CLUSTER_PROBABILITY: f64 = 0.5;
pub const RIFT_CLUSTER_OFFSET: (f64, f64) = (0.8, 3.5);
pub const RIFT_LOOP_WIDTH_CELLS: f64 = 2.0;
pub const RIFT_CARVING_GYR: f64 = 1.0;
pub const RIFT_END_GYR: f64 = 2.0;
pub const RIFT_MAX_FAULTS: u64 = 20_000;
pub const CANYON_DEPTH_FRACTION: f64 = 1.0;
pub const CANYON_RENEWAL_THICKNESS: f64 = 50.0;
/// Trough traces wander by this much over this length; their walls take this fraction
/// of the width; the walls are scalloped by this fraction over this length; the width
/// swells and narrows by this fraction over this length (chasmata); and the floor tapers
/// to the ends with this power.
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
/// Broader folds under the fine fabric, as the tessera's ridge belts are.
pub const TESSERA_FOLD_RELIEF: f64 = 800.0;
pub const TESSERA_FOLD_LENGTH: f64 = 120e3;
/// Tessera stand as crustal plateaus: thickened crust held up isostatically, 1 to 4 km
/// above the plains (Ovda, Thetis, Alpha Regio). The uplift ramps in over a width of the
/// smoothed-height score about the tessera threshold, so the plateau edge is a slope.
pub const TESSERA_PLATEAU_UPLIFT: f64 = 2_000.0;
pub const TESSERA_EDGE_WIDTH: f64 = 300.0;
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
/// Novae: radial fracture systems round a share of the coronae, which with the rises'
/// graben give Venus its spidery fracture belts. Fractures per nova, their length as a
/// fraction of the corona radius, width and depth.
pub const NOVA_FRACTION: f64 = 0.3;
pub const NOVA_FRACTURES: (u64, u64) = (6, 14);
pub const NOVA_LENGTH: (f64, f64) = (1.0, 3.0);
pub const NOVA_WIDTH: (f64, f64) = (6e3, 15e3);
pub const NOVA_DEPTH: f64 = 250.0;
pub const DOME_DIAMETER: (f64, f64) = (20e3, 60e3);
pub const DOME_HEIGHT: f64 = 800.0;
pub const DOME_TOP_ROUGHNESS: f64 = 0.08;
/// The Io regime: paterae per terawatt of tidal power (Io: some 400), the latitude weight
/// of the paterae (1 at the equator, this at the poles, as the tidal heat comes out at low
/// latitudes), their diameters, m, log-uniform and skewed small, their depth, m, the
/// wall width as a fraction of the radius, the polygon share of a fault-bounded outline
/// and its noise, the sides of that polygon, and the share of paterae active at a time.
pub const IO_PATERAE_PER_TW: f64 = 5.4;
pub const IO_PATERA_POLAR_WEIGHT: f64 = 0.3;
pub const IO_PATERA_DIAMETER: (f64, f64) = (10e3, 200e3);
pub const IO_PATERA_SIZE_SKEW: f64 = 2.0;
pub const IO_PATERA_DEPTH: f64 = 1_200.0;
pub const IO_PATERA_WALL_FRACTION: f64 = 0.08;
pub const IO_PATERA_POLYGON: f64 = 0.6;
pub const IO_PATERA_OUTLINE: f64 = 0.12;
pub const IO_PATERA_SIDES: (u32, u32) = (5, 8);
pub const IO_PATERA_ACTIVE_FRACTION: f64 = 0.2;
/// Years since a quiet patera's floor last held a lake, log-uniform, and the time over
/// which a floor goes from black silicate to a sulfur coat.
pub const IO_FLOOR_AGE_YR: (f64, f64) = (1.0, 1e6);
pub const IO_FLOOR_COAT_YR: f64 = 2_000.0;
/// Mountains: tilted crustal blocks raised by the subsidence of the crust under the
/// burial, two per terawatt (Io: some 150), their diameters and heights, m, the tilt of a
/// block as a fraction of its height, elongation up to this aspect ratio, the noise on
/// the outline as a fraction of the radius, the ridged texture's relative amplitude and
/// length, and their lifetimes, Gyr: raised and buried in a few hundred million years.
pub const IO_MOUNTAINS_PER_TW: f64 = 2.0;
pub const IO_MOUNTAIN_DIAMETER: (f64, f64) = (80e3, 300e3);
pub const IO_MOUNTAIN_HEIGHT: (f64, f64) = (3e3, 15e3);
pub const IO_MOUNTAIN_TILT: f64 = 0.5;
pub const IO_MOUNTAIN_ASPECT: f64 = 2.5;
pub const IO_MOUNTAIN_OUTLINE: f64 = 0.3;
pub const IO_MOUNTAIN_TEXTURE: f64 = 0.08;
pub const IO_MOUNTAIN_TEXTURE_LENGTH: f64 = 12e3;
pub const IO_MOUNTAIN_LIFETIME_GYR: (f64, f64) = (0.1, 0.5);
/// Flow fields: from every active patera and this share of the quiet ones; their reach,
/// m, log-uniform; the vent's pool as a fraction of the reach; lobes per field, the
/// angular spread of a lobe, rad, the lobing of the margins, the softness of the edge as
/// a fraction of the reach, the thickness, m; lobe ages, years, log-uniform, and the
/// time over which a flow goes from black to brown under the sulfur rain. A field is a
/// connected body of this fraction of the reach with the lobes reaching out from it.
pub const IO_QUIET_FLOW_FRACTION: f64 = 0.15;
pub const IO_FLOW_REACH: (f64, f64) = (40e3, 300e3);
pub const IO_FLOW_VENT: f64 = 0.15;
pub const IO_FLOW_BODY: f64 = 0.25;
pub const IO_FLOW_LOBES: (u32, u32) = (6, 14);
pub const IO_FLOW_SPREAD: (f64, f64) = (0.12, 0.45);
pub const IO_FLOW_LOBING: f64 = 0.3;
pub const IO_FLOW_EDGE: f64 = 0.08;
pub const IO_FLOW_THICKNESS: f64 = 20.0;
pub const IO_FLOW_AGE_YR: (f64, f64) = (1.0, 3e5);
pub const IO_FLOW_DARK_YR: f64 = 150.0;
/// Io's deposits accumulate (F56): a flow's colour is coated over by the fallout of the
/// plains in this many years (Geissler et al. 2004 on the resurfacing at a centimetre a
/// year), so old flows fade into the plains; every patera carries this many generations
/// of flow fields and of plumes, the older ones faded by their age; and a flow's surface
/// varies by this share over this length, m, the overlapping flow units within a lobe.
pub const IO_FLOW_COAT_YR: f64 = 2.0e4;
pub const IO_FLOW_GENERATIONS: u32 = 3;
pub const IO_PLUME_GENERATIONS: u32 = 3;
pub const IO_FLOW_TEXTURE: f64 = 0.4;
pub const IO_FLOW_TEXTURE_LENGTH: f64 = 15e3;
/// Plumes: the share of centres that are gas-rich, Pele-class, and the share of active
/// paterae of each class with a plume at a time; the eruption speeds of each class, m/s,
/// whose ballistic range v²/g is the ring radius (Pele's ring 1,200 to 1,400 km across
/// from about a kilometre a second in Io's gravity; Prometheus-class rings a few hundred
/// kilometres from half that), with log-normal scatter; the ring's width and the inner
/// fall as fractions of the radius, a Pele-class plume leaving mostly the ring and a
/// Prometheus-class one filling its disc; the frost share of a Pele-class deposit and
/// the red share of a Prometheus-class one; the wobble of the ring and the mottling of
/// the fallout; years since a quiet centre's last plume, log-uniform; the burial time of
/// frost fallout, longer by this factor where the cold keeps the frost; and the fading
/// time of red sulfur, which Surt's and Tvashtar's rings lost within a couple of years.
pub const IO_PELE_CLASS_FRACTION: f64 = 0.05;
pub const IO_PLUME_FRACTION_PELE: f64 = 0.5;
pub const IO_PLUME_FRACTION_PROMETHEUS: f64 = 0.1;
pub const IO_PLUME_SPEED_PELE: f64 = 1_000.0;
pub const IO_PLUME_SPEED_PROMETHEUS: f64 = 500.0;
pub const IO_PLUME_SPEED_SCATTER: f64 = 0.25;
pub const IO_RING_WIDTH: f64 = 0.25;
pub const IO_RING_INNER_PELE: f64 = 0.3;
pub const IO_RING_INNER_PROMETHEUS: f64 = 0.8;
pub const IO_PELE_FROST: f64 = 0.25;
pub const IO_PROMETHEUS_RED: f64 = 0.15;
pub const IO_RING_EDGE_WOBBLE: f64 = 0.08;
pub const IO_DEPOSIT_MOTTLING: f64 = 0.35;
pub const IO_PLUME_AGE_YR: (f64, f64) = (10.0, 1e6);
pub const IO_FALLOUT_FADE_YR: f64 = 100.0;
pub const IO_FROST_COLD_FACTOR: f64 = 3.0;
pub const IO_RED_FADE_YR: f64 = 2.0;
pub const IO_SURFACE_AGE_GYR: f64 = 1e-4;
/// Slack on a feature's birth and death against the step times' rounding, Gyr.
pub const IO_ALIVE_TOLERANCE_GYR: f64 = 1e-6;
/// Colours, linear reflectance: sulfur plains, sulfur reddened by radiation in the
/// cold, SO2 frost, fresh silicate, an old flow under its sulfur coat, a quiet patera
/// floor, a plume's red sulfur and bare mountain rock. The global frost fall covers this
/// share where it survives the day; the day's peak temperature at which half of it
/// survives and the width of that transition, K (SO2 sublimes fast above about 120 K);
/// the temperature below which radiation keeps sulfur red and its width, K; and the
/// slopes between which a mountain sheds its coat, degrees.
pub const IO_SULFUR: [f32; 3] = [0.65, 0.55, 0.25];
pub const IO_POLAR: [f32; 3] = [0.38, 0.24, 0.13];
pub const IO_FROST: [f32; 3] = [0.86, 0.86, 0.88];
pub const IO_FLOW_DARK: [f32; 3] = [0.05, 0.05, 0.05];
pub const IO_FLOW_OLD: [f32; 3] = [0.20, 0.15, 0.09];
pub const IO_PATERA_FLOOR: [f32; 3] = [0.14, 0.09, 0.05];
pub const IO_PLUME_RED: [f32; 3] = [0.55, 0.20, 0.08];
pub const IO_MOUNTAIN_ROCK: [f32; 3] = [0.32, 0.28, 0.22];
pub const IO_FROST_BACKGROUND: f64 = 0.2;
pub const IO_FROST_HALF_K: f64 = 120.0;
pub const IO_FROST_WIDTH_K: f64 = 5.0;
pub const IO_RED_HALF_K: f64 = 118.0;
pub const IO_RED_WIDTH_K: f64 = 6.0;
pub const IO_MOUNTAIN_BARE_SLOPE_DEG: (f64, f64) = (5.0, 15.0);
/// Magma production (F17): fraction of the heat flux carried by magma, by regime, with a
/// tidal boost that reaches 0.9 at 1 W/m². Latent plus sensible heat per m³ of magma.
/// Provisional, Phase 3.
pub const MAGMA_EFFICIENCY_MOBILE: f64 = 0.15;
pub const MAGMA_EFFICIENCY_EPISODIC: f64 = 0.05;
pub const MAGMA_EFFICIENCY_STAGNANT: f64 = 0.01;
pub const MAGMA_EFFICIENCY_ICE: f64 = 3e-5;
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
/// Below this pressure, Pa, a body is airless for its albedo: the surface's own
/// reflectance stands for the Bond albedo (F52). Above it dust, haze and clouds take
/// part, and until the climate stage computes them the pressure classes below stand in.
pub const AIRLESS_PRESSURE: f64 = 10.0;
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
/// Sub-grid crater roughness (F54): the rms relief of a saturated field of craters as a
/// share of one crater's depth, and the Hurst exponent of cratered terrain, by which the
/// relief falls toward finer scales (Rosenburg et al. 2011: 0.76 to 0.96 on the Moon).
pub const ROUGHNESS_RMS_SHARE: f64 = 0.35;
pub const ROUGHNESS_HURST: f64 = 0.8;
/// The roughness's strength follows the surface age smoothed over this many passes of
/// the one-cell kernel (a sigma of half the square root of that in cells, about two
/// cells), so that the edge of a young deposit on old ground is a short ramp of roughness
/// and never a staircase of cells.
pub const ROUGHNESS_BLEND_PASSES: usize = 12;
/// Variation within units (F55): a body's lavas drift in composition over its eruptive
/// history, so a flow unit's reflectance departs from the mean by up to this share,
/// drifting over this period, Gyr, and written by the eruption into the surface's
/// composition layer, which its own ejecta keeps (the maria's titanium: albedos 0.06 to
/// 0.11, basins filled by flows of different ages and tones); slopes shed their regolith from
/// this steepness, degrees, the regolith's stable slope, and stay this share fresher than
/// flat ground at the angle of repose, degrees (optical maturity maps: the walls of young
/// craters bright, those of old craters near the stable slope only a little less mature
/// than their surroundings); and a lane of renewed ice is clean by a share drawn in this
/// range.
pub const LAVA_COMPOSITION_RANGE: f64 = 0.18;
pub const LAVA_COMPOSITION_PERIOD_GYR: f64 = 1.0;
pub const SLOPE_SHED_START_DEG: f64 = 20.0;
pub const SLOPE_FRESHNESS: f64 = 0.3;
pub const REPOSE_SLOPE_DEG: f64 = 30.0;
pub const GROOVE_CLEAN_RANGE: (f64, f64) = (0.55, 1.0);
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
/// ejecta edge as a fraction of its radius: the edge wanders by 15% in eleven to
/// sixty-one lobes round the crater, as the ragged margins of real blankets do; a few
/// low harmonics read as petals at map scale.
pub const EJECTA_BRIGHTNESS_TEXTURE: f64 = 0.25;
pub const EJECTA_EDGE_RAGGEDNESS: f64 = 0.15;
/// A lobe of the edge narrower than this many cells along it is not drawn: at the grid's
/// scale it would be a tooth, not a lobe.
pub const EJECTA_LOBE_MIN_CELLS: f64 = 4.0;
/// The blanket's cover of the ground (F51): the deposit thins as r^EJECTA_EXPONENT and
/// covers the old surface by 1 - 2^-(edge / r)^3, all of it at the rim, half at the
/// continuous edge, a fifth at three radii, tapering to nothing over this outer share of
/// the footprint, where the deposit is too thin to leave a record. The surface age and
/// the ray brightness both follow the cover, so a blanket grades into its rays and into
/// the ground instead of ending on a circle.
pub const EJECTA_COVERAGE_TAPER: f64 = 0.25;
/// A raw cover above this counts as complete, the few per cent of old ground showing
/// through a deposit that thick being below what the maturity can tell: the cover is
/// scaled by it, so the blanket renews the ground wholly out to 1.3 radii.
pub const EJECTA_FULL_COVER: f64 = 0.92;
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
/// The loop threshold is taken on this coarse grid at every quality, the standard one,
/// so the preview's loop population is the export's; a preview stamps its smallest loop
/// craters on half as many cells, and the finish restores them.
pub const LOOP_REFERENCE_FACE: f64 = 256.0;
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
/// Lateral mixing of the regolith by small impacts: a contact between two materials
/// blurs in colour over a zone this wide, m, after the reference age, Gyr. Li and
/// Mustard 2000 measure 4 to 5 km across lunar mare-highland contacts of 3.5 Gyr; the
/// zone widens as the square root of the age, as a diffusive process does.
pub const REGOLITH_MIXING_LENGTH: f64 = 4.5e3;
pub const REGOLITH_MIXING_REFERENCE_GYR: f64 = 3.5;
pub const REGOLITH_MIXING_MAX_PASSES: usize = 24;
/// Sublimation lag, the maturity of an icy surface: ice sublimes at the latitude's mean
/// temperature (Murphy and Koop 2005 for the vapour pressure, the Hertz-Knudsen rate) and
/// leaves its non-ice share behind; the lag is opaque past this mass per area, a couple
/// of metres of dust, since gardening and frost keep a thinner one from hiding the ice
/// (Moore et al. 1999 on Callisto's lag, metres thick after tens of metres of ice lost).
/// The poles hold no less than this temperature. The non-ice share of each ice
/// end-member is in the materials.
pub const LAG_OPAQUE_KG_M2: f64 = 2_000.0;
pub const ICE_POLAR_TEMPERATURE_FLOOR_K: f64 = 40.0;
/// Clean ice that an impact brings up is mixed back into the regolith by gardening over
/// this time, Gyr, whatever the temperature: Callisto's bright craters are its young ones.
pub const ICE_EJECTA_FADE_GYR: f64 = 0.5;
/// Viscous relaxation of an ice shell (F39): the viscosity of ice at its melting point,
/// Pa s, and the activation energy of its creep, J/mol (Goldsby and Kohlstedt 2001; the
/// basal viscosity of a convecting shell is 1e13 to 1e15), below a temperature at which
/// the law is capped; the depth a load of a wavelength reaches into the shell as a
/// fraction of the wavelength.
pub const ICE_VISCOSITY_MELT: f64 = 1.0e14;
pub const ICE_CREEP_ACTIVATION: f64 = 60e3;
pub const ICE_VISCOSITY_FLOOR_K: f64 = 60.0;
pub const RELAXATION_DEPTH_FRACTION: f64 = 0.16;
/// The cold elastic lid holds loads shorter than this many times its thickness up
/// without flow; longer ones ride on the ductile ice below (the flexural wavelength of a
/// thin plate is a few times its thickness).
pub const ELASTIC_SUPPORT_FACTOR: f64 = 4.0;
/// Rock in the ice stiffens it: the viscosity rises by this many decades per unit of the
/// crust's non-ice share, so a crust a third rock and salts is some four thousand times
/// stiffer than clean ice (Ceres's needs 1e25 Pa s near the surface against clean ice's
/// 3e21 at the same temperature: Fu et al. 2017).
pub const ICE_ROCK_VISCOSITY_DECADES: f64 = 12.0;
/// The tidal stress field (F41): the elastic moduli of ice (Young's 9 GPa, Poisson 0.33;
/// Gammon et al. 1983); the shell's reduction of the Love number h₂ per unit of its
/// thickness over the radius (a 10 km shell on Europa takes h₂ from 1.25 to 1.1, a 40 km
/// shell on Enceladus to a quarter); l₂ over h₂ for a shell over an ocean (Wahr et al.
/// 2006: Europa 0.3 and 1.2); phases sampled through an orbit; the tension at which the
/// shell cracks (Hoppa et al. 1999 fit cycloids with 40 kPa) and the tensile strength of
/// ice; the crack density per unit of tension over strength per metre of shell; the lag
/// of a non-synchronously rotating shell behind the tidal axis before it relaxes, degrees,
/// and its period as a multiple of the surface's retention age (the lineae's
/// orientations turn by half a turn over the oldest ones Europa shows, so the zones of tension sweep
/// every longitude).
pub const ICE_YOUNG_MODULUS: f64 = 9.0e9;
pub const ICE_POISSON: f64 = 0.33;
pub const SHELL_RIGIDITY_FACTOR: f64 = 20.0;
pub const L2_OVER_H2: f64 = 0.25;
pub const TIDAL_PHASE_SAMPLES: usize = 32;
pub const ICE_CRACK_STRESS: f64 = 4.0e4;
pub const ICE_TENSILE_STRENGTH: f64 = 1.0e6;
pub const CRACK_DENSITY_FACTOR: f64 = 3.0;
pub const NSR_SHIFT_DEG: f64 = 1.0;
pub const NSR_PERIOD_RETENTION_FACTOR: f64 = 2.0;
/// Cracks answer the cracks already open (Pollard and Aydin 1988): within this many
/// shell thicknesses of an open fissure a propagating crack feels its stress shadow and
/// turns toward meeting it at a right angle, by this share of the way at most; within
/// half a width of it (or a cell) it stops, having reached a free surface; healed ridges
/// and bands it crosses. This share of new cracks start at the tip of an older one,
/// where its stress concentrates, and run on from it at the field's heading of their
/// own time, so that sets branch and splay.
pub const CRACK_SHADOW_SHELLS: f64 = 0.5;
pub const CRACK_DEFLECTION: f64 = 0.6;
pub const CRACK_BRANCH_SHARE: f64 = 0.35;
/// Lineae (F42): lengths as a share of the radius, log-uniform (Europa 100 to 1,500 km; a
/// crack scales with the body it crosses); the width of a double ridge and its ridge
/// height, m, with the central trough as a fraction of it (Europa's ridges stand 100 to
/// 300 m over 1 to 2 km); the wander of a crack along the slow field; the share of the
/// tension that must be diurnal for a crack to be a cycloid; a cycloid's propagation
/// speed, m/s (Hoppa et al. 1999: about 3 km/h), the orbits it may run and the steps per
/// orbit; the share of cracks that pull apart into bands, their widths as a share of the
/// radius (Europa 10 to 30 km) and floor depth, m; and the time the ocean's salts take to fade under frost, Gyr.
pub const LINEA_LENGTH_RAD: (f64, f64) = (0.064, 0.96);
pub const LINEA_WIDTH: f64 = 2e3;
pub const LINEA_RIDGE_HEIGHT: f64 = 150.0;
pub const LINEA_TROUGH_FRACTION: f64 = 0.6;
pub const LINEA_WOBBLE_RAD: f64 = 0.08;
pub const LINEA_WOBBLE_LENGTH: f64 = 200e3;
/// A crack's relief and width vary along it by this rms share over segments this long,
/// m: ridges rise and fall, troughs pinch and widen, as the real ones do.
pub const CRACK_RELIEF_VARIATION: f64 = 0.35;
pub const CRACK_VARIATION_LENGTH: f64 = 15e3;
pub const LINEA_DIURNAL_SHARE: f64 = 0.6;
pub const CYCLOID_SPEED: f64 = 0.8;
pub const CYCLOID_MAX_ORBITS: usize = 12;
pub const CYCLOID_PHASE_STEPS: usize = 64;
pub const BAND_FRACTION: f64 = 0.15;
pub const BAND_WIDTH_RAD: (f64, f64) = (0.0064, 0.019);
pub const BAND_DEPTH: f64 = 100.0;
/// A band is a lens: its width scaled by (4 s (1 - s)) to this power along its length, s
/// the share of the length; and it is lineated inside, parallel to its axis, at this
/// spacing, m (or three cells), by this much, m.
pub const BAND_TAPER_POWER: f64 = 0.35;
pub const BAND_LINEATION_SPACING: f64 = 1.5e3;
pub const BAND_LINEATION: f64 = 40.0;
pub const SALT_FADE_GYR: f64 = 0.03;
/// Focused heating (F43): a shell over an ocean thins under one pole where the tidal
/// dissipation concentrates, by up to this share of its thickness, in proportion to the
/// tidal flux over this flux, W/m², rising with the shell's thickness over the radius
/// over this ratio (a thin shell like Europa's cracks everywhere instead); the thinning
/// tapers to nothing this far from the pole; the terrain within the inner half of the
/// taper is active once the thinning passes this share: renewed on this timescale, Gyr,
/// and lying this share of its isostatic depression low (Enceladus's south polar terrain
/// lies about a kilometre low over a shell a few kilometres thick under a mean of
/// twenty); its folds between the cracks, m and m; the margin ridge where the spreading
/// terrain meets the old, m and m; and the wobble of the thinning's edge round the pole
/// as a share of its radius, by a two-cycle field.
pub const POLAR_THINNING_MAX: f64 = 0.8;
pub const POLAR_FOCUS_FLUX: f64 = 0.02;
pub const POLAR_FOCUS_SHELL_RATIO: f64 = 0.03;
pub const POLAR_TERRAIN_RADIUS_DEG: f64 = 70.0;
pub const POLAR_ACTIVE_PROFILE: f64 = 0.5;
/// The renewed terrain is the thinned shell's core, the taper above this, and the ground
/// round the fissures that opened in the thinned zone: within this distance, m, of them
/// in the core, narrowing to nothing where the taper falls to this outer profile, so
/// that the terrain follows the fissures out as tongues; its folds fade over this
/// distance, m, from its edge. A fissure leaving the core meets a thicker, stiffer shell
/// and refracts along the contrast, turning to run round the margin between the core
/// and the outer profile, and stops in the thick shell where the taper falls below the
/// last value.
pub const POLAR_CORE_PROFILE: f64 = 0.75;
pub const POLAR_SPREAD_WIDTH: f64 = 30e3;
pub const POLAR_OUTER_PROFILE: f64 = 0.5;
pub const POLAR_FOLD_TAPER: f64 = 10e3;
pub const CRACK_STOP_PROFILE: f64 = 0.3;
pub const POLAR_ACTIVE_THINNING: f64 = 0.3;
/// The focusing comes and goes (F53; O'Neill and Nimmo 2010: overturns every 0.1 to 1
/// Gyr, each active about 10 Myr): episodes recur every this many conductive times of
/// the shell, d² / κ (a fit at the short end of their range, so that Enceladus's come
/// 0.53 Gyr apart and with the lanes renew 0.40 of it, against the 57% mapped as
/// tectonized), each this long, Gyr. Between them the thinned zone is thick again, and the
/// body reorients to put each new one at the pole (Nimmo and Pappalardo 2006), so the
/// terrain of each old one lies where the wander left it, its fissures turned with it.
pub const FOCUSING_PERIOD_FACTOR: f64 = 10.0;
pub const FOCUSING_EPISODE_DURATION_GYR: f64 = 0.01;
/// Thermal diffusivity of ice, m²/s.
pub const ICE_THERMAL_DIFFUSIVITY: f64 = 1.2e-6;
pub const POLAR_RENEWAL_GYR: f64 = 0.002;
/// A renewed polar terrain keeps no relief shorter than this share of the radius while
/// its episode lasts: in the thinned, hot zone every wavelength up to the zone's scale
/// relaxes in a geological instant, so no ghost of an older crater survives there.
pub const POLAR_RENEWAL_LENGTH_SHARE: f64 = 0.2;
pub const POLAR_DEPRESSION_FRACTION: f64 = 0.6;
pub const POLAR_FOLD_AMPLITUDE: f64 = 250.0;
pub const POLAR_FOLD_SPACING: f64 = 2.5e3;
/// The folds' spacing varies by this share with the set field, and their amplitude runs
/// from this share of the full in the slack of a set to the full at its heart (the set
/// field counted full at this value), so the fabric comes in fields of ridges of varied
/// spacing rather than one grating; a segment dies out where its field falls below this.
pub const POLAR_FOLD_SPACING_VARIATION: f64 = 0.3;
pub const POLAR_FOLD_SLACK: f64 = 0.5;
pub const POLAR_FOLD_SET_FULL: f64 = 0.3;
pub const POLAR_FOLD_SEGMENT_LEVEL: f64 = 0.75;
/// The fabric of a renewed plain: its ridges wander by this many radians of phase over
/// this length, m, so that they are sinuous rather than straight; they come in segments
/// about this long, m, that die out and relay en echelon; and under them lies the
/// tectonic roughness of a plain deformed at every scale, this rms, m, self-affine from
/// this longest wavelength, m, down to the grid (the Hurst exponent of F54), so that the
/// renewed plains are rough at every scale as the cratered terrain is.
pub const POLAR_FOLD_WANDER: f64 = 3.0;
/// The fold sets curve: the distance that sets a ridge's phase is warped by this much, m,
/// over this length, m, so that the sets bend as ridge belts do instead of running
/// straight across a terrain where no fissure steers them.
pub const POLAR_FOLD_WARP: f64 = 8e3;
pub const POLAR_FOLD_WARP_LENGTH: f64 = 60e3;
pub const POLAR_FOLD_WANDER_LENGTH: f64 = 6e3;
pub const POLAR_FOLD_SEGMENT: f64 = 10e3;
pub const POLAR_TECTONIC_RMS: f64 = 120.0;
pub const POLAR_TECTONIC_LENGTH: f64 = 20e3;
pub const POLAR_MARGIN_HEIGHT: f64 = 400.0;
pub const POLAR_MARGIN_WIDTH: f64 = 8e3;
pub const POLAR_MARGIN_WOBBLE: f64 = 0.15;
/// The renewed terrain's floor and folds blend into the old surface over this width, m,
/// centred on the edge: the thermal transition at the margin, about a thin shell's
/// thickness, across which the old craters sink part way and fade in rather than end on
/// a cut.
pub const POLAR_MARGIN_RAMP: f64 = 15e3;
/// The focused heating is a field, not a shape: the thinning varies across the zone by
/// this share between the upwellings and downwellings of the shell's convection, in
/// cells this many shell thicknesses wide (a layer's cells are about as wide as it is
/// deep), so the thinned zone's outline, its depression and its renewed terrain are
/// lobate rather than circular. The polar folds follow the nearest fissure within this
/// distance, m, and only beyond it the field's direction.
pub const CONVECTION_STRUCTURE: f64 = 0.3;
pub const CONVECTION_CELL_SHELLS: f64 = 2.0;
pub const POLAR_FABRIC_RANGE: f64 = 40e3;
/// The folds' phase is a soft minimum, over this length, m, of the distances across the
/// fissures near (to each point's local line, not to the point, so that the ridges run
/// parallel past a fissure's end instead of wrapping round it), so that the sets of
/// neighbouring fissures bend into one another instead of meeting on a seam.
pub const POLAR_FABRIC_SOFTNESS: f64 = 3e3;
/// The margin belt: the spreading is taken up in folds parallel to the renewed terrain's
/// edge, decaying inward over this distance, m, at this spacing, m, this high, m, broken
/// into segments this long, m, along the edge that shift against their neighbours by up
/// to this phase, rad, as the folds of a thrust belt do: en echelon, not concentric rings.
pub const POLAR_BELT_WIDTH: f64 = 15e3;
pub const POLAR_BELT_SPACING: f64 = 12e3;
pub const POLAR_BELT_HEIGHT: f64 = 300.0;
pub const POLAR_BELT_SEGMENT: f64 = 25e3;
pub const POLAR_BELT_OFFSET: f64 = 2.5;
/// A shell renews its whole surface within the retention age only when it is thin
/// against the body, under this share of the radius; a thick shell keeps its old craters
/// and renews only where it is thin.
pub const WHOLESALE_SHELL_RATIO: f64 = 0.02;
/// Grooved terrain (F44): lanes of new ice drawn as great-circle swaths, their lengths
/// log-uniform between these angles, rad, their widths between these shares of the
/// radius (Ganymede 30 to 250 km, Enceladus 3 to 24 km, its ridged belts), wobbling along
/// them by this share, gathered in belts of extension drawn by a field of this many
/// cycles round the body raised to this power (the extension localises, so the dark
/// terrain between the belts stays whole), this share of them laid beside an earlier
/// lane turned by up to
/// this angle and set this many widths to its side (so the lanes come in parallel sets,
/// as sulci do), until they cover the resurfaced share, at most this many; this share of
/// the lanes is smooth
/// bright terrain without grooves; the old relief inside is smoothed over this length,
/// m (or a twentieth of the radius), and dropped by this much, m; the grooves'
/// amplitude, m, and spacing, m, which the grid coarsens to at least this many cells;
/// the length over which groove sets turn, m, and the turn; the share of the surface
/// resurfaced per unit of tidal flux above this floor, W/m², over this flux, W/m², on a
/// shell thicker than the wholesale ratio, less what the focusing takes to the pole
/// (Callisto's 0.08 mW/m² never stretched its shell); and the age of the
/// resurfacing as a share of the body's age (Ganymede's bright terrain is about two
/// gigayears old).
/// A lane is a swath of grooved terrain, a rift zone of the shell: its width is this many
/// elastic thicknesses of the shell at the episode's time, log-uniform, as rift zones are
/// a few to a few tens of times the layer that breaks, and does not scale with the body
/// (Ganymede's swaths run 50 to 300 km on an 11 km lithosphere); its length is
/// this many widths, log-uniform, and at most this share of the radius in radians
/// (Ganymede's lanes run hundreds to a couple of thousand kilometres, Enceladus's sulci
/// to 380 km, a quarter of the way round it).
pub const GROOVE_LANE_WIDTH_ELASTIC: (f64, f64) = (4.0, 24.0);
pub const GROOVE_LANE_ASPECT: (f64, f64) = (6.0, 40.0);
/// No lane is wider than this share of the radius, whatever the shell: a geometric cap for
/// small bodies (Enceladus's sulci, 38 km).
pub const GROOVE_LANE_MAX_WIDTH_SHARE: f64 = 0.15;
pub const GROOVE_LANE_MAX_LENGTH_RAD: f64 = 1.5;
/// A lane's edges wander by this rms share of its width along it, in a self-affine
/// profile from three widths down to half a width (the Hurst exponent of F54).
pub const GROOVE_WIDTH_WOBBLE: f64 = 0.25;
pub const GROOVE_BELT_FREQUENCY: f64 = 1.2;
pub const GROOVE_BELT_POWER: f64 = 2.0;
pub const GROOVE_NEIGHBOUR_SHARE: f64 = 0.85;
pub const GROOVE_NEIGHBOUR_TURN_DEG: f64 = 3.0;
pub const GROOVE_NEIGHBOUR_OFFSET: (f64, f64) = (0.5, 1.2);
pub const GROOVE_MAX_LANES: u64 = 12000;
/// A lane laid beside an earlier one belongs to its generation: its time is the
/// parent's within this share of the episode, and it is no shorter than this share of
/// the parent. A lane that starts on its own runs along the strike of the shell's field
/// of extension, which the stress alone does not resolve (Ganymede's sets): a smooth
/// field of this many cells per radian about the body, read as an angle from east, which
/// drifts by this much through the episode so that generations cross-cut.
pub const GROOVE_SET_TIME_SPREAD: f64 = 0.05;
pub const GROOVE_SET_LENGTH_SHARE: f64 = 0.5;
/// A set grows beside one of this many lanes laid lately, while the extension is there.
pub const GROOVE_SET_RECENT: usize = 24;
pub const GROOVE_HEADING_FREQUENCY: f64 = 1.2;
pub const GROOVE_HEADING_DRIFT: f64 = 1.0;
pub const GROOVE_SMOOTH_LANE_SHARE: f64 = 0.2;
pub const GROOVE_SMOOTH_LENGTH: f64 = 100e3;
/// The lanes' drop spreads over the lithosphere's flexural parameter at their time, so a
/// lane narrower than it sinks little and none ends on a step.
pub const GROOVE_LANE_DROP: f64 = 200.0;
pub const GROOVE_AMPLITUDE: f64 = 250.0;
pub const GROOVE_SPACING: f64 = 7e3;
pub const GROOVE_MIN_CELLS: f64 = 3.0;
pub const GROOVE_SET_LENGTH: f64 = 300e3;
pub const GROOVE_SET_TURN: f64 = 0.4;
/// A lane's centre line is traced along the tidal field where the shell cracks and along
/// a great circle where nothing steers it, in points this share of the radius apart (or
/// the half width, if less); its width closes over this many widths at each end, as the
/// bounding faults die out, (the distance to the end over that) to this power, so the
/// ends are blunt rather than drawn to a point.
pub const GROOVE_LANE_POINT_SPACING_RAD: f64 = 0.02;
pub const GROOVE_END_TAPER_WIDTHS: f64 = 1.0;
pub const GROOVE_LANE_TAPER_POWER: f64 = 0.35;
/// A lane's edge is drawn over this many cells: the resolution's own ramp, so that the cut
/// through the old, rough terrain is not a staircase of cells.
pub const GROOVE_EDGE_CELLS: f64 = 2.5;
/// The grooves inside a lane (F44): their spacing varies by this share and their strength
/// from this share of the full to the full with the set field; they wander by this many
/// radians of phase over this length, m, and come in segments about this long, m, that
/// die out and relay (this level sets how much of the lane carries a groove at all); a
/// bounding trough this deep, m, runs just inside each edge; and under them lies the
/// tectonic roughness of a swath deformed at every scale, this rms, m, self-affine from
/// this longest wavelength, m, down to the grid.
pub const GROOVE_SPACING_VARIATION: f64 = 0.3;
pub const GROOVE_SLACK: f64 = 0.5;
pub const GROOVE_WANDER: f64 = 2.0;
pub const GROOVE_WANDER_LENGTH: f64 = 25e3;
pub const GROOVE_SEGMENT_LENGTH: f64 = 60e3;
pub const GROOVE_SEGMENT_LEVEL: f64 = 0.7;
pub const GROOVE_MARGIN_TROUGH: f64 = 150.0;
pub const GROOVE_ROUGHNESS: f64 = 80.0;
pub const GROOVE_ROUGH_LENGTH: f64 = 20e3;
/// The lanes' amount (F44): the share renewed is the peak tidal flux of the body's
/// history above this floor, W/m², over this flux, W/m²: a fit, so that the resonance
/// pulse that renews two thirds of Ganymede, 8 mW/m² at an eccentricity of 0.008 (the
/// low end of Showman and Malhotra 1997's range), leaves its shell thick enough, near
/// 40 km, to keep the dark terrain's 70 km craters through the episode while the 150 km
/// ones relax to palimpsests; a stronger pulse thins the shell under 30 km and flattens
/// the 70 km craters too, which the dark terrain does not show. The resurfacing is an
/// episode, not an instant: the lanes'
/// times are spread uniformly about the excess-weighted mean time of the flux above the
/// floor, over the width that matches its spread, so that younger lanes cut older ones;
/// a resonance episode (F50) gives a pulse, a steady flux the body's whole age.
pub const ICE_RESURFACING_FLUX_FLOOR: f64 = 0.005;
pub const ICE_RESURFACING_FLUX: f64 = 0.0045;
/// A resonance episode (F50): captured, the eccentricity rises over this share of each
/// half of the episode, holds at the forced value, and falls over the same share after
/// the escape; a trapezoid with a plateau of three fifths.
pub const RESONANCE_RAMP_SHARE: f64 = 0.4;
/// Chaos (F45): on a shell over an ocean thin against the body, the share of the surface
/// disrupted into chaos, up to this much, in proportion to the tidal flux over this flux,
/// W/m²; the patches' diameters, m, log-uniform, placed in a heat field of this many
/// cycles round the body whose troughs take this share fewer than its crests; the drop of the matrix, m, its hummocks, m, at this spacing, m;
/// the share of a patch left standing as rafts of the old surface and their scale, m;
/// and the wobble of a patch's outline.
pub const CHAOS_FRACTION_MAX: f64 = 0.3;
pub const CHAOS_FLUX: f64 = 0.05;
pub const CHAOS_DIAMETER: (f64, f64) = (5e3, 300e3);
pub const CHAOS_HEAT_FREQUENCY: f64 = 2.0;
pub const CHAOS_HEAT_CONTRAST: f64 = 0.5;
pub const CHAOS_DROP: f64 = 150.0;
pub const CHAOS_HUMMOCK: f64 = 60.0;
pub const CHAOS_HUMMOCK_SPACING: f64 = 1.5e3;
pub const CHAOS_RAFT_SHARE: f64 = 0.45;
pub const CHAOS_RAFT_SCALE: f64 = 15e3;
pub const CHAOS_EDGE_WOBBLE: f64 = 0.35;
/// The outline's wobble has about this many waves round a patch.
pub const CHAOS_EDGE_WAVES: f64 = 6.0;
/// Cryovolcanic domes (F46): on a shell over rock this share of the ice magma budget
/// extrudes as domes of volumes drawn log-uniform between these shares of the body's
/// radius cubed (Ceres 100 to 3,000 km³, Ahuna Mons about 1,500), with this height over
/// diameter (Ahuna Mons 4 km over 17); the streaks on their flanks as a share of the
/// height; the share of a dome's substrate that is clean ice; and the salts on it.
pub const CRYODOME_SHARE: f64 = 0.7;
pub const CRYODOME_VOLUME_SHARE: (f64, f64) = (1e-6, 3e-5);
pub const CRYODOME_ASPECT: f64 = 0.24;
pub const CRYODOME_FLANK_ROUGHNESS: f64 = 0.06;
pub const CRYODOME_CLEAN: f64 = 0.5;
pub const CRYODOME_SALT: f64 = 0.4;
/// Faculae (F46): on a shell over rock a crater deeper than this, m, brings brines up
/// into its floor, which leave salts over this share of the floor's radius in patches at
/// this scale, m, that darken over this time, Gyr (Occator's are a few million years old
/// in a crater of twenty million).
pub const FACULA_BRINE_DEPTH: f64 = 3e3;
pub const FACULA_FLOOR_SHARE: f64 = 0.35;
pub const FACULA_PATCH_SCALE: f64 = 4e3;
pub const FACULA_FADE_GYR: f64 = 0.15;
/// Plume fallout (F47): the fresh ice a polar plume lays down at its pole, m/Gyr, when
/// the tidal flux focuses fully (Enceladus about 0.1 mm/yr near the stripes), falling
/// off over this angle, deg, to this share of it far away; the depth that buries a lag or
/// a salt, m; and the rate that keeps a full frost cover, m/Gyr.
pub const PLUME_FALLOUT_RATE: f64 = 1e5;
pub const PLUME_FALLOUT_ANGLE_DEG: f64 = 25.0;
pub const PLUME_FALLOUT_FAR: f64 = 0.01;
pub const PLUME_BURIAL_DEPTH: f64 = 1e-3;
pub const PLUME_FROST_RATE: f64 = 1e3;
/// The warm fissures a plume vents from sinter the ice to coarse grains, darker and bluer
/// than the frost round them (Enceladus's tiger stripes): the tint at the pole.
pub const VENT_ICE_TINT: [f32; 3] = [0.66, 0.78, 0.92];
/// The fissures under a plume are widened and deepened by the venting: troughs this deep,
/// m, this wide, m, with flanking ridges this share of the depth (Enceladus's tiger
/// stripes: 500 m troughs 2 km wide with 100 m flanks).
pub const VENT_TROUGH_DEPTH: f64 = 450.0;
pub const VENT_WIDTH: f64 = 6e3;
pub const VENT_RIDGE_SHARE: f64 = 0.4;
/// A vent is an elliptical opening: its width and its depth fall to these shares of the
/// middle's at its tips.
pub const VENT_TIP_WIDTH: f64 = 0.35;
pub const VENT_TIP_DEPTH: f64 = 0.4;
/// The fissures of a plume's thinned zone form in a cascade (F57; Hemingway, Rudolph and
/// Manga 2020): the first through the thinnest shell, within this share of the cascade's
/// spacing of the pole, then each next one where the bending tension of the last one's
/// flanking ridges peaks on the plate broken at it, a quarter flexural parameter away
/// (35 km on Enceladus for an elastic thickness of 5.2 km), up to this many on either
/// side; and the time-bin tag the cascade's fissures carry in their ids.
pub const VENT_CASCADE_FIRST_OFFSET: f64 = 0.5;
pub const VENT_CASCADE_MAX: usize = 8;
pub const VENT_CASCADE_BIN: u64 = 0xFFFF;
/// A plume's fissures wander by this many radians over this length, m, more than the
/// lesser cracks do (the stripes' sinuous traces), and each ends in a fork: from each tip
/// a branch leaves at this angle, degrees, and bends back to the field over this length,
/// m, running this share of the fissure's length at this share of its width (the hooked,
/// bifurcating terminations of the real stripes).
pub const VENT_WOBBLE_RAD: f64 = 0.3;
pub const VENT_WOBBLE_LENGTH: f64 = 50e3;
pub const VENT_FORK_ANGLE_DEG: f64 = 60.0;
pub const VENT_FORK_BIAS_LENGTH: f64 = 15e3;
pub const VENT_FORK_LENGTH_SHARE: f64 = 0.25;
pub const VENT_FORK_WIDTH_SHARE: f64 = 0.6;
/// Hemispheric asymmetry (F48): a locked moon's trailing hemisphere takes its parent's
/// magnetospheric plasma, which darkens and reddens the ice there by this tint at the
/// apex in a Jupiter-class magnetosphere (a parent over this mass, kg; Europa's and
/// Ganymede's trailing sides are the darker and redder) and by this share of it in a
/// Saturn-class one (over this mass).
pub const TRAILING_TINT: [f32; 3] = [0.76, 0.69, 0.62];
pub const MAGNETOSPHERE_JUPITER_MASS: f64 = 1e27;
pub const MAGNETOSPHERE_SATURN_MASS: f64 = 1e26;
pub const MAGNETOSPHERE_SATURN_SHARE: f64 = 0.3;
/// The plasma processes a trailing hemisphere's optical surface, its top millimetre, in
/// this time, Gyr, in a Jupiter-class magnetosphere and proportionally longer in a weaker
/// one (F58; Cooper et al. 2001: micrometres on Europa in about ten years, centimetres in
/// under a million); a plume's fallout that buries the millimetre faster hides the
/// darkening by the ratio of the two times.
pub const TRAILING_DARKENING_GYR: f64 = 1e-5;
pub const GAS_CONSTANT: f64 = 8.314_462;
/// Thermal segregation: the ice the warm latitudes sublime condenses again where it is
/// cold, so a body that loses ice at its equator wears frost at its poles (Ganymede's
/// caps above about 40° latitude, Callisto's fainter ones). The frost covers up to this
/// share where the mean temperature is below the half point by the width, K, once the
/// equator has sublimed this much ice, kg/m².
pub const FROST_COVER_MAX: f64 = 0.45;
pub const FROST_HALF_K: f64 = 105.0;
pub const FROST_WIDTH_K: f64 = 3.0;
pub const FROST_SUPPLY_KG_M2: f64 = 500.0;
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
