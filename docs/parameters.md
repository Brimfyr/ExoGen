# Parameter graph specification

Revision 0.6, 2026-10-11. This document is the contract for the `params` module of `exogen-core`. Every node the user can see or set is listed here with its unit, kind, default or derivation, and the source of the relation. The formula appendix at the end gives the relations by number. Values marked *fit* are calibrated against the reference bodies during the phase that implements the stage; the derivation form is fixed here, the constants are not.

## Implementation status

Phase 1 implements the graph mechanism and the groups star, orbit, body, rotation, interior, atmosphere, volatiles, climate (zero-dimensional subset), cratering (scaling subset), detail, output, seeds and multipliers, plus `tectonics.dichotomy_mode`, `tectonics.dichotomy_amplitude` and `tectonics.initial_crust_thickness`. Phase 2 adds the cratering population nodes (`cratering.loop_threshold`, `cratering.min_diameter`, `cratering.retention_age`, `cratering.empirical_saturation`, `cratering.secondary_fraction`, `cratering.oblique_fraction`), `surface.space_weathering_timescale`, `surface.solidification_time`, the flux history nodes `cratering.early_flux_ratio` and `cratering.early_flux_decay`, and the material controls `materials.bedrock`, `materials.regolith_mix`, `materials.deep_bedrock`, `materials.deep_bedrock_depth` and `materials.colour_space`. The impact flux is written in time since formation, so the crater record begins at the solidification time whatever the body's age; craters carry an impact angle; and each preset carries its own master seed. The executable table is `crates/exogen-core/src/params/registry.rs`; every constant lives in `crates/exogen-core/src/constants.rs` with its source, and constants marked provisional there are calibrated in the phase named. Phase 3 adds the thermal history curves (`interior.heat_flux_history`, `interior.mantle_temperature_history`, `interior.lithosphere_history`, `interior.relief_limit_history`, `interior.contraction_history`, `orbit.eccentricity_history`), `volcanism.rate`, `volcanism.resurfacing_rate` and `volcanism.history`, the contraction nodes `tectonics.contraction_strain` and `tectonics.scarp_total_shortening`, the plume nodes `tectonics.plume_rise_count`, `tectonics.plume_rise_diameter`, `tectonics.plume_rise_height` and `tectonics.rift_extension`, the episodic-lid nodes `tectonics.resurfacing_event_age` and `tectonics.tessera_fraction`, the ice-shell nodes `interior.love_number_h2`, `tectonics.tidal_stress`, `tectonics.ice_crack_density` and `tectonics.nsr_period` (Phase 4), the Io-regime node `tectonics.mountain_count`, the volcanism counts `volcanism.shield_count`, `volcanism.corona_count`, `volcanism.dome_count` and `volcanism.paterae_count`, `volcanism.province_count`, and `materials.lava`. Phase 4 adds the partial-resurfacing nodes `tectonics.polar_thinning`, `tectonics.ice_resurfaced_fraction` and `tectonics.ice_resurfacing_age`. The remaining tectonics, volcanism, weather, surface-process and material nodes are added by the phases that implement those stages. A test asserts that every derivation reads only the nodes it declares, and `exogen validate` checks the derived values of the fourteen reference presets against measured values.

## Conventions

- Node ids are `group.name`, snake case.
- Units are SI unless the unit column says otherwise. Display units may differ (AU, bar, Gyr, km) and are converted at the UI boundary only.
- **Kinds**: `root` (an input with a default, no derivation), `derived` (computed from upstream nodes; shown dimmed; can be set, which turns it into an override), `expert` (derived or root, hidden unless expert mode is on), `multiplier` (an artistic scale factor, default 1, always visible, never derived), `curve` (a piecewise-linear history over simulated time; its default is a constant or a derived curve), `list` (a user-editable table of records).
- Every derived node records the ids it was derived from, so the UI can show them on hover.
- Setting a derived node pins it. Downstream nodes keep deriving from the pinned value. "Reset to derived" clears the pin.
- Ranges are soft: values outside them raise a plausibility warning, never an error.
- References point to PLAN.md section 12 by author and year.

## star

| id | unit | kind | default or derivation | source |
|---|---|---|---|---|
| star.mass | M_sun | root | 1.0; range 0.08 to 3 | |
| star.luminosity | L_sun | derived | F1 mass–luminosity relation, then F3 age brightening | Eker et al. 2018; Gough 1981 |
| star.radius | R_sun | derived | F2 | Eker et al. 2018 |
| star.teff | K | derived | Stefan–Boltzmann from luminosity and radius | |
| star.spectral_type | text | derived | Lookup from teff; display only | |
| star.age | Gyr | root | 4.57; warn if above star.main_sequence_lifetime | |
| star.main_sequence_lifetime | Gyr | derived | 10 x mass / luminosity_zams | |
| star.luminosity_history | L_sun | curve | F3 evaluated over body age | Gough 1981 |
| star.xuv_relative | ratio | derived | F4; relative to the present Sun at 1 AU, scaled by 1/a² at the body | Ribas et al. 2005 |

## orbit

| id | unit | kind | default or derivation | source |
|---|---|---|---|---|
| orbit.primary | enum star, planet | root | star | |
| orbit.semi_major_axis | m | root | 1 AU; display AU for a star primary, parent radii for a planet primary | |
| orbit.eccentricity | 1 | root; curve allowed | 0.0167 | |
| orbit.period | s | derived | F5 Kepler, primary mass plus body mass | |
| orbit.mean_motion | rad/s | derived | 2π / period | |
| orbit.insolation | W/m² | derived | L / (4π d²) where d is the distance to the star (the parent's orbit for a moon) | |
| parent.mass | kg | root (planet primary) | 1 M_jup | |
| parent.radius | m | root (planet primary) | 1 R_jup | |
| parent.semi_major_axis | m | root (planet primary) | 5.2 AU | |
| parent.eccentricity | 1 | root (planet primary) | 0.05 | |

## body

| id | unit | kind | default or derivation | source |
|---|---|---|---|---|
| body.radius | m | root | 6.371e6; warn below 4e5 (spheroid assumption) | |
| body.density | kg/m³ | root | 5514 | |
| body.mass | kg | derived | (4/3)π R³ ρ | |
| body.gravity | m/s² | derived | F6 | |
| body.escape_velocity | m/s | derived | F7 | |
| body.age | Gyr | derived | star.age minus 0.05 | |
| body.core_mass_fraction | 1 | derived | F8 three-layer mass–radius fit | Zeng et al. 2016; Fortney et al. 2007 |
| body.ice_mass_fraction | 1 | derived | F8 | |
| body.rock_mass_fraction | 1 | derived | 1 minus core minus ice | |
| body.composition_class | enum iron_rich, rock, rock_ice, ice_dominated | derived | From the three fractions and surface temperature (an ice fraction above 0.1 with a surface above the water triple point is an ocean world, class rock) | |
| body.uncompressed_density | kg/m³ | derived | From F8 layer densities at zero pressure; display only | |
| body.moment_of_inertia_factor | 1 | derived | 0.33 differentiated rock-iron, 0.35 rock-ice, 0.40 undifferentiated | |
| body.oblateness | 1 | derived | F9; display and export datum | |

## rotation

| id | unit | kind | default or derivation | source |
|---|---|---|---|---|
| rotation.locking_timescale | Gyr | derived | F10 | Gladman et al. 1996 |
| rotation.mode | enum free, locked | derived | locked if locking_timescale is below body.age | |
| rotation.period | s | root when free; derived equal to orbit.period when locked | 86164 | |
| rotation.obliquity | deg | root when free; derived 0 when locked | 23.44 | |
| rotation.angular_velocity | rad/s | derived | 2π / period | |
| rotation.regime | enum rapid, rhines, slow | derived | F11 from the equatorial Rossby deformation radius against the radius | Haqq-Misra et al. 2018 |

## interior

| id | unit | kind | default or derivation | source |
|---|---|---|---|---|
| interior.love_number_k2 | 1 | expert derived | 0.3 rock, 0.02 small cold rock, 0.25 ice shell over ocean, 0.05 solid ice | |
| interior.love_number_h2 | 1 | derived | 1 + k2 for a fluid-like interior, over 1 + 20 x shell thickness / radius for an ice shell (Europa 1.14, Enceladus 0.27, Ganymede 0.77) | Phase 4, Wahr et al. 2006 |
| interior.dissipation_q | 1 | expert derived | 100 rock, 30 ice, 12 with surface ocean | |
| interior.radiogenic_heat_per_kg | W/kg | derived | F12 at body.age, times rock mass fraction | Turcotte and Schubert |
| interior.radiogenic_flux | W/m² | derived | radiogenic power / (4π R²) | |
| interior.tidal_power | W | derived | F13 | Peale, Cassen and Reynolds 1979 |
| interior.tidal_flux | W/m² | derived | tidal_power / (4π R²) | |
| interior.secular_cooling_fraction | 1 | expert derived | 1.0 for R above 3000 km falling to 0.2 at 500 km, times radiogenic flux | fit |
| interior.surface_heat_flux | W/m² | derived | radiogenic plus tidal plus secular | |
| interior.heat_flux_history | W/m² | curve | F12 and F13 evaluated over body age with the eccentricity curve, plus secular cooling | Phase 3 |
| interior.mantle_temperature_history | K | curve | Thermal evolution over the body's life, ending at interior.mantle_temperature | Phase 3 |
| interior.lithosphere_history | m | curve | F14 along the heat-flux history | Phase 3 |
| interior.relief_limit_history | m | curve | F15 scaled by the elastic thickness relative to today's, floored at 0.15 | Phase 3 |
| interior.contraction_history | m | curve | F16 accumulated over the body's life | Phase 3 |
| interior.mantle_temperature | K | expert derived | Thermal evolution, calibrated so Earth is 1600 K now | fit |
| interior.lithosphere_thickness | m | derived | F14 conductive thickness | |
| interior.elastic_thickness | m | derived | Depth to the 650 K isotherm for rock, 170 K for ice, from the conductive profile | |
| interior.crust_strength | Pa | derived | 150 MPa rock, weakened linearly to 50 MPa as surface temperature approaches 750 K; 5 MPa ice | fit |
| interior.max_relief | m | derived | F15 | |
| interior.contraction_total | m | derived | F16 at the present | Byrne et al. 2014 |
| interior.magma_production_rate | m³/s | derived | F17 | fit |
| interior.ice_shell_thickness | m | derived | F18 for rock_ice and ice_dominated bodies | |
| interior.has_subsurface_ocean | bool | derived | Ice layer thicker than the conductive shell thickness | |
| interior.tectonic_regime | enum mobile_lid, episodic_lid, stagnant_lid, ice_shell_ocean, ice_shell_rock | derived | F19 | |

## atmosphere

| id | unit | kind | default or derivation | source |
|---|---|---|---|---|
| atmosphere.surface_pressure | Pa | root; curve allowed | 101325 | |
| atmosphere.composition | map species to volume fraction | root | N2 0.78, O2 0.21, Ar 0.01; species N2, O2, CO2, CH4, H2, He, Ar, SO2, NH3, H2O | |
| atmosphere.mean_molecular_mass | kg/mol | derived | Fraction-weighted | |
| atmosphere.scale_height | m | derived | F20 at the mean surface temperature | |
| atmosphere.column_mass | kg/m² | derived | pressure / gravity | |
| atmosphere.greenhouse_optical_depth | 1 | derived | F21 | fit |
| atmosphere.equilibrium_temperature | K | derived | F22 with the Bond albedo | |
| atmosphere.surface_temperature_mean | K | derived | Global mean from the climate stage; pinning it rescales the greenhouse term to match | |
| atmosphere.bond_albedo | 1 | derived | From the climate stage's surface and cloud state; initial guess by composition class | |
| atmosphere.jeans_parameter | map species to value | derived | F23 at the mean surface temperature | |
| atmosphere.exobase_temperature | K | expert derived | Surface temperature plus 700 K x xuv^0.3 x (1 − 0.8 f_CO2). Earth about 1000 K, Mars about 250 K | provisional, Phase 5 |
| atmosphere.retention_warning | text | derived | Species with Jeans parameter below 36 at the exobase temperature | |
| atmosphere.oxidant_availability | 1 | derived | O2 fraction plus 0.1 x water vapour fraction, times pressure in bar; clipped to 1 | fit |
| atmosphere.haze_production | kg/m²/Gyr | derived | Tholin: CH4 fraction x xuv_relative x constant; sulfuric: SO2 fraction x constant | fit |
| atmosphere.crater_cutoff_diameter | m | derived | F24 | fit |
| atmosphere.condensables | list species | derived | Species in composition or inventory whose saturation curve crosses the surface temperature range | |
| atmosphere.rayleigh_coefficients | per m at 440, 550, 680 nm | derived | Per-gas cross sections weighted by fraction and number density; sidecar only | |
| atmosphere.aerosols | list | derived | Per aerosol: optical depth, scale height, asymmetry, colour, latitude profile; sidecar only | |
| atmosphere.dust_background_optical_depth | 1 | derived | From aeolian dust lifting and settling; sidecar only | fit |

## volatiles

| id | unit | kind | default or derivation | source |
|---|---|---|---|---|
| volatiles.inventory | map species to global equivalent layer in m | root; curve allowed | H2O 2700; species H2O, CO2, CH4, N2, NH3, SO2, C2H6 | |
| volatiles.phase_state | map species to enum gas, liquid, solid, liquid_and_solid, absent | derived | F25 against the surface temperature range | |
| volatiles.working_liquid | enum none, water, methane_ethane, ammonia_water, nitrogen | derived | The liquid species with the largest inventory | |
| volatiles.sea_level | m | derived | F26 from inventory and hypsometry | |
| volatiles.ocean_fraction | 1 | derived | Area below sea level | |
| volatiles.cap_mode | enum none, polar, equatorial, hemispheric | derived | Polar below 54° obliquity, equatorial above, hemispheric when locked | |
| volatiles.frost_points | map species to K | derived | F25 inverted at the species partial pressure | |
| volatiles.seasonal_caps | bool | derived | True when a condensable's frost point lies inside the seasonal temperature swing | |

## tectonics

| id | unit | kind | default or derivation | source |
|---|---|---|---|---|
| tectonics.regime | enum as interior.tectonic_regime | derived | Mirror of interior.tectonic_regime; the user-facing control | |
| tectonics.dichotomy_mode | enum none, harmonic, giant_impact | root | harmonic | |
| tectonics.dichotomy_amplitude | m of crust | derived | 5000 rock, 30000 for a Mars-type stagnant lid, 0 for mobile lid default | fit |
| tectonics.initial_crust_thickness | m | derived | 40000 rock-iron and rock, 20000 ice shells | |
| tectonics.plate_count | 1 | derived | F27 | heuristic |
| tectonics.plate_speed_scale | m/s | derived | 5 cm/yr x (surface_heat_flux / 0.087)^0.5 | heuristic |
| tectonics.continental_fraction_initial | 1 | derived | 0.1 mobile lid; grows in the loop | |
| tectonics.hotspot_count | 1 | derived | 10 x (area / Earth area) | heuristic |
| tectonics.reorganisation_interval | Myr | expert | 150 | |
| tectonics.contraction_strain | 1 | derived | 2 x contraction_total / radius; the scarps take this fraction of the area up as shortening, the longest of them scaling with its square (F35) | Phase 3 |
| tectonics.scarp_total_shortening | m | derived | contraction_strain x 2π R | Phase 3 |
| tectonics.plume_rise_count | 1 | derived | Stagnant or episodic lid at least 3,000 km in radius: floor(surface_heat_flux / 0.02 W/m²), at most 3; else 0. Mars 1, Venus 3 | Phase 3, heuristic |
| tectonics.plume_rise_diameter | m | derived | 1.5 R on a stagnant lid (Tharsis), 0.4 R on an episodic lid (Beta Regio) | Phase 3, fit |
| tectonics.plume_rise_height | m | derived | 0.28 x max_relief (Tharsis 7 km, Beta Regio about 3 km) | Phase 3, fit |
| tectonics.rift_extension | m | derived | 4 x plume_rise_height when there is a rise: the extension across a rise's rift zone (Valles Marineris some 30 km on Tharsis's 8 km; Devana Chasma 10 to 20 km on Beta Regio's 4 km), taken up by a population of troughs from chasmata to graben drawn until it is spent (F38) | Phase 3, Andrews-Hanna 2012 |
| tectonics.resurfacing_event_age | Gyr | derived | The crater retention age for an episodic lid, 0.5 on Venus; else 0 | Phase 3 |
| tectonics.tessera_fraction | 1 | derived | 0.08 episodic lid: the highest-standing share of the surface, measured without the plume rises' domes, that survives the event; it stands as crustal plateaus 2 km above the plains | Phase 3 |
| tectonics.mountain_count | 1 | derived | Io regime: 2 per TW of tidal power, some 150 on Io | Phase 3, heuristic |
| tectonics.tidal_stress | Pa | derived | The largest tension the diurnal tide puts on a locked moon's shell through an orbit (F41); Europa 81 kPa, Enceladus 119 kPa, Ganymede 3 kPa | Phase 4 |
| tectonics.ice_crack_density | per m | derived | On a shell over an ocean whose tidal tension exceeds 40 kPa: 3 x (tidal_stress / 1 MPa) / shell thickness, the cracks a transect crosses per metre (Europa one per 32 km); else 0 | Phase 4, fit |
| tectonics.nsr_period | Gyr | derived | Period of a cracking shell's rotation against the tidal axis: 2 x cratering.retention_age, so the lineae turn half a turn over the surface's age and the zones of tension sweep every longitude; 0 without cracks | Phase 4, heuristic |
| tectonics.polar_thinning | 1 | derived | Share of the shell's thickness lost under one pole where the tidal dissipation focuses (F43): 0.8 x min(tidal_flux / 0.02 W/m², 1) x (1 − exp(−(d / R) / 0.03)) on a shell over an ocean, else 0; Enceladus about 0.75, Europa 0.3, Ganymede 0. Which pole is drawn from seed.tectonics | Phase 4, fit |
| tectonics.ice_resurfaced_fraction | 1 | derived | Share of a shell over an ocean, thicker than 0.02 R, renewed as lanes of grooved terrain (F44): clipped ((tidal_flux − 0.005 W/m²) / 0.05 W/m²) x (1 − polar_thinning / 0.8), so Callisto's 0.08 mW/m² gives none and Enceladus's 19 mW/m² a fiftieth; a thinner shell cracks everywhere instead and gets 0. Ganymede's preset pins 0.65, its heating episode being outside the eccentricity history | Phase 4, fit |
| tectonics.ice_resurfacing_age | Gyr | derived | 0.4 x body.age when there are lanes (Ganymede's bright terrain about 2 Gyr); else 0 | Phase 4, Zahnle et al. 2003 |
| tectonics.chaos_fraction | 1 | derived | Fraction of the shell thinner than 0.5 x mean, from the heat flux pattern | |
| tectonics.focused_heating_pole | superseded | | Replaced by tectonics.polar_thinning; which pole is a coin tossed from seed.tectonics | Phase 4 |
| tectonics.nonsynchronous_rotation | 1 | expert | 0 | |

## volcanism

| id | unit | kind | default or derivation | source |
|---|---|---|---|---|
| volcanism.rate | m³/s | derived | interior.magma_production_rate | Phase 3 |
| volcanism.resurfacing_rate | m/Myr | derived | rate / area | Phase 3 |
| volcanism.style | enum none, silicate, sulfur_silicate, cryo_water_ammonia, cryo_brine | derived | By composition class and heat flux; sulfur_silicate above 1 W/m² on rock | |
| volcanism.province_count | 1 | derived | Magmatic provinces the eruptions favour: the plume rises first, then one random region per 2,000 km of radius, Gaussian 40° wide, the rest of the body at 3% weight; Moon 1 (the Procellarum terrane), Mercury 1, Mars 2, Venus 3 | Phase 3, heuristic |
| volcanism.shield_count | 1 | derived | 5 on the main plume rise, half as many on each further rise; hotspot_count for a mobile lid in Phase 6 | Phase 3, heuristic |
| volcanism.corona_count | 1 | derived | Episodic lid: 1.1e-12 per m², some 500 on Venus | Phase 3, heuristic |
| volcanism.dome_count | 1 | derived | Episodic lid: 3.3e-13 per m², some 150 on Venus | Phase 3, heuristic |
| volcanism.arc_volcano_density | per m of arc | derived | 1 per 70 km of convergent boundary | |
| volcanism.flood_episode_count | 1 | derived | From the heat flux history: one per 10 Myr of rate above threshold, capped | fit |
| volcanism.paterae_count | 1 | derived | Io regime: 5.4 per TW of tidal power, some 400 on Io | Phase 3, heuristic |
| volcanism.fresh_lava_albedo | 1 | derived | 0.07 basalt, 0.4 cryolava, 0.1 sulfur-coated | |
| volcanism.active_vent_count | 1 | derived | Vents younger than 1 Myr at the output epoch | |
| volcanism.history | m³/s | curve | F17 along the thermal history; the extrusive share erupts as floods (F36) | Phase 3 |

## cratering

| id | unit | kind | default or derivation | source |
|---|---|---|---|---|
| cratering.production_function | enum neukum_lunar, neukum_mars, outer_system, custom | derived | By location: inner system by default, outer_system for a planet primary beyond 4 AU | Neukum et al. 2001; Zahnle et al. 2003 |
| cratering.flux_relative_to_moon | 1 | derived | F28 | Ivanov 2001 |
| cratering.impact_velocity | m/s | derived | F29 | |
| cratering.impactor_density | kg/m³ | derived | 2700 inner system, 600 outer | |
| cratering.target_density | kg/m³ | derived | 2700 rock crust, 920 ice | |
| cratering.simple_complex_transition | m | derived | F30 | Pike 1980; Schenk 2002 |
| cratering.peak_ring_onset | m | derived | 8.5 x simple_complex_transition | |
| cratering.multiring_onset | m | derived | 20 x simple_complex_transition | |
| cratering.morphometry_coefficients | table | expert | Pike lunar set by default; per-class overrides | Pike 1977, 1980 |
| cratering.min_diameter | m | derived | Larger of atmosphere.crater_cutoff_diameter and two export pixels | |
| cratering.early_flux_ratio | x steady | expert | Rate at formation over the steady rate: a b e^(4.5 b) / c of the chronology, about 1.57e4 (F31) | Neukum et al. 2001 |
| cratering.early_flux_decay | Gyr | expert | 1 / b of the chronology, 0.144 Gyr (F31) | Neukum et al. 2001 |
| cratering.flux_history | 1 | curve | Rate over the body's life relative to the present: zero before surface.solidification_time, then the early bombardment decaying into the steady flux (F31) | Neukum et al. 2001 |
| cratering.giant_basins | list of diameter, lat, lon, age | root | empty: sampled from the production function | |
| cratering.diffusivity | m²/Myr | derived | Gardening diffusivity at the present flux: 5.5 x (gravity / 1.62)^0.5, x 10 with an atmosphere above 1 kPa; the flux history scales it back in time, so a crater's degradation follows the impacts that came after it | Fassett and Thomson 2014; fit |
| cratering.ray_fading_timescale | Gyr | derived | 1.0 x (a / 1 AU)² airless; 0.05 with an atmosphere | fit |
| cratering.rampart_ejecta | bool | derived | Regolith volatile content above threshold | |
| cratering.secondary_fraction | 1 | expert | 0.2 | |
| cratering.oblique_fraction | 1 | expert | 0.05: the share of elliptical cavities, formed below about 12°. Impact angles follow sin 2θ; ejecta and rays thicken downrange below 45°, open an uprange zone of avoidance below 30° and grow side lobes below 15° | Bottke et al. 2000; Gault and Wedekind 1978 |

## climate

| id | unit | kind | default or derivation | source |
|---|---|---|---|---|
| climate.greenhouse_offset | K | derived | surface_temperature_mean minus equilibrium_temperature | |
| climate.equator_pole_contrast | K | derived | From the energy balance solution | |
| climate.day_night_contrast | K | derived | From the energy balance solution when locked | |
| climate.hadley_edge_latitude | deg | derived | F32 | Held and Hou 1980 |
| climate.transport_coefficient | W/m²/K | derived | F33 | Williams and Kasting 1997 |
| climate.lapse_rate | K/m | derived | gravity / c_p, times 0.65 when the working liquid is present | |
| climate.precipitation_scale | m/yr | derived | From evaporation over liquid area and transport | |
| climate.season | 1 | root | 0; orbital phase of the snapshot | |
| climate.subsolar_longitude | deg | root | 0 | |
| climate.ice_albedo_feedback | bool | expert | true | |

## weather

| id | unit | kind | default or derivation | source |
|---|---|---|---|---|
| weather.enabled | bool | derived | surface_pressure above 100 Pa | |
| weather.resolution | cells per face | derived | 128 standard, 256 high | |
| weather.spinup | model days | expert | 60 | |
| weather.diurnal_cycle | bool | derived | true when free-rotating and column_mass below 1e5 kg/m² | |
| weather.cloud_fraction_target | 1 | derived | From humidity statistics; reported, not imposed | |
| weather.tropical_cyclone_count | 1 | derived | From warm-liquid area between 8° and 30° and shear | fit |
| weather.dust_storm_activity | 1 | derived | From aeolian lifting rate | |
| weather.global_dust_storm | bool | root | false | |
| weather.global_dust_optical_depth | 1 | root | 2.0 when the storm is on | |
| weather.planetary_wave_contrast | 1 | derived | 0.02 for superrotating thick atmospheres, else 0 | |
| weather.plume_activity | 1 | derived | volcanism.active_vent_count | |
| weather.sequence_frames | 1 | root | 1 | |
| weather.sequence_interval | model hours | root | 6 | |

## surface processes

| id | unit | kind | default or derivation | source |
|---|---|---|---|---|
| surface.fluvial_enabled | bool | derived | working_liquid present and precipitation_scale above 0 | |
| surface.erodibility | 1/yr per unit | expert derived | By bedrock material; calibrated on Earth relief | fit |
| surface.stream_power_m | 1 | expert | 0.5 | Whipple and Tucker 1999 |
| surface.stream_power_n | 1 | expert | 1.0 | |
| surface.glacial_enabled | bool | derived | Any cell with positive ice mass balance | |
| surface.ice_softness | 1/(Pa³ s) | expert derived | Glen's law A(T) for the ice species | |
| surface.aeolian_enabled | bool | derived | surface_pressure above 100 Pa and loose material present | |
| surface.threshold_friction_speed | m/s | derived | F34 | Iversen and White 1982 |
| surface.mass_wasting_diffusivity | m²/Myr | derived | Same as cratering.diffusivity | |
| surface.angle_of_repose | deg | derived | 32 rock, 30 ice | |
| surface.solidification_time | Gyr after formation | derived | 0.3 x (R / 1737 km)^0.5, x 0.7 for ice shells: when the crust is solid and rigid enough to keep large relief, so the crater record begins. Pin to set | provisional, Phase 3 heat budget |
| surface.space_weathering_timescale | Gyr | derived | 1.0 x (a / 1 AU)² airless, capped at 4 Gyr; disabled above 100 Pa | fit |
| surface.oxidation_rate | 1/Gyr | derived | oxidant_availability x temperature factor | fit |
| surface.tholin_deposition_rate | m/Gyr | derived | haze_production / deposit density | fit |
| surface.sulfur_deposition | bool | derived | volcanism.style is sulfur_silicate | |
| surface.sublimation_rates | map species to m/Myr | derived | From vapour pressure and temperature | |
| surface.thermal_segregation | bool | derived | Airless, icy, locked, with a dark contaminant | |
| surface.hemispheric_asymmetry | bool | derived | Locked moon of a planet primary | |

## materials

Composition class sets only the defaults of this group. Surface colour comes from the per-cell mixture of end-members that the processes fill: the bedrock and regolith mix as the base, volatile ices from the inventory wherever the climate makes them stable (nitrogen and methane ice over water-ice bedrock on a cold body, for example), and dust, tholin, sulfur, salt and weathering deposits from their stages. Every derived node here can be pinned, so a surface can be composed by hand as well as derived.

| id | unit | kind | default or derivation | source |
|---|---|---|---|---|
| materials.bedrock | enum basalt, anorthosite, felsic, water_ice, ammonia_water_ice, rock_ice_mixture, carbonaceous, sulfur_coated | derived | Default from composition class, regime and volcanism style: basalt for rock, anorthosite for a thick primary crust on an airless rock body, felsic for continental crust on a mobile lid, water_ice for ice classes, carbonaceous for dark primitive bodies. Pin to choose | Phase 2 |
| materials.regolith_mix | map end-member to fraction | derived | Baseline surface mixture before deposits, from the bedrock and the class, for example water_ice 0.6 and rock 0.4 on a dark rock-ice terrain. Pin to set the surface composition by hand | Phase 2 |
| materials.deep_bedrock | enum as bedrock, plus norite, low_reflectance_material, impact_melt | derived | What large craters excavate from below the surface layer: norite on an airless rocky stagnant lid, low-reflectance material on an iron-rich body, clean water ice under the lag of an ice shell, otherwise the bedrock. Pin to choose | Phase 2, provisional |
| materials.lava | enum as bedrock, plus mare_basalt, smooth_plains | derived | What floods and flows leave: mare basalt on an airless rocky stagnant lid, the brighter smooth plains of an iron-rich body, basalt otherwise, dirty ice for cryolava. Pin to choose | Phase 3, provisional |
| materials.deep_bedrock_depth | m | derived | 25 km rocky, 10 km iron-rich, 0.5 km ice shells. A crater excavating a tenth of its diameter (Melosh 1989) deeper than this exposes the deep bedrock in its floor and ejecta | Phase 2, provisional |
| materials.end_members | list | derived | The end-members active for this world: the bedrock, the regolith mix, the volatile ices present in the inventory, and the deposits the processes produce | USGS Spectral Library v7; RELAB |
| materials.custom_end_members | list | root | User-defined materials with a name, a reflectance spectrum or an RGB, roughness and specular; usable in the regolith mix | Phase 2 |
| materials.colour_space | enum linear_srgb, acescg | root | linear_srgb | |
| materials.wetness_darkening | bool | expert | true | |
| materials.vegetation_enabled | bool | derived | Phase 8: liquid water, temperature band, atmosphere | |
| materials.vegetation_pigment | linear RGB | root | chlorophyll reflectance | Phase 8 |

## detail and output

| id | unit | kind | default or derivation | source |
|---|---|---|---|---|
| detail.quality | enum preview, standard, high | root | standard | |
| detail.loop_resolution | cells per face | derived | 256, 256, 512 by quality | |
| detail.finishing_resolution | cells per face | derived | 512, 1024, 2048 by quality | |
| detail.recipes_enabled | bool | root | true | |
| output.width | px | root | 8192; height is width / 2 | |
| output.maps | set | root | albedo, height, roughness, specular, normal, masks, cloud classes, unit_id | |
| output.cloud_classes | set | root | all six | |
| output.datum | enum mean_radius, sea_level, ellipsoid | root | sea_level when a liquid exists, else mean_radius | |
| output.normal_convention | enum opengl, directx | root | opengl | |
| output.tile_size | px | root | 2048 | |
| output.epoch | Gyr before present | root | 0 | |
| output.sequence_epochs | list | root | empty | |

## seeds and multipliers

| id | kind | default |
|---|---|---|
| seed.master | root | 1 |
| seed.<stage> | derived | hash(master, stage id); one per stage; any may be pinned |
| multiplier.crater_density | multiplier | 1 |
| multiplier.erosion_intensity | multiplier | 1 |
| multiplier.volcanic_output | multiplier | 1 |
| multiplier.dust_cover | multiplier | 1 |
| multiplier.frost_cover | multiplier | 1 (scales volatile frost and ice deposit extent) |
| multiplier.relief | multiplier | 1 |
| multiplier.cloud_fraction | multiplier | 1 |
| multiplier.storm_activity | multiplier | 1 |
| multiplier.plume_activity | multiplier | 1 |
| multiplier.haze_optical_depth | multiplier | 1 (sidecar value) |

## Plausibility rules

Each rule is a predicate with a severity (note, warning) and an explanation that names the quantities. Rules never block.

| Rule | Severity | Explanation template |
|---|---|---|
| Liquid requested but not stable | warning | "{species} liquid requested. Surface temperature {T} K is outside its liquid range at {P} bar ({T_melt} to {T_boil} K)." |
| Atmosphere retention | warning | "Jeans parameter below 36 at the exobase temperature {T} K for {species}. These species are unlikely to survive {age} Gyr." |
| Pressure below triple point | note | "Surface pressure {P} is below the {species} triple point; no liquid can exist." |
| Locking inconsistency | warning | "Locking timescale {t} Gyr is shorter than body age {age} Gyr but rotation is pinned free." |
| Relief above strength limit | warning | "Requested relief {h} km exceeds the crustal strength limit {h_max} km." |
| Body too small for spheroid | note | "Radius {R} km is below 400 km; the body is unlikely to be spheroidal." |
| Density outside physical range | warning | "Density {ρ} is outside 900 to 8000 kg/m³." |
| Star age beyond main sequence | warning | "Star age {age} Gyr exceeds the main-sequence lifetime {t_ms} Gyr." |
| Thick atmosphere on low gravity | note | "Column mass {m} kg/m² on gravity {g} m/s² has no solar-system analogue." |
| Io regime with surviving craters pinned | note | "Resurfacing rate {r} m/Myr erases craters faster than they form; the pinned crater density will be mostly buried." |

## Formula appendix

Symbols: G gravitational constant, σ Stefan–Boltzmann, k_B Boltzmann, m_u atomic mass unit. Solar units where stated.

- **F1 Mass–luminosity.** Six-piece fit of Eker et al. 2018, log L = a log M + b (solar units), with (a, b) = (2.028, −0.976) below 0.45 M_sun, (4.572, −0.102) to 0.72, (5.743, −0.007) to 1.05, (4.329, 0.010) to 2.4, (3.967, 0.093) to 7, (2.865, 1.105) to 31; normalised so that one solar mass gives one solar luminosity. This is the reference value at the age fraction of F3.
- **F2 Mass–radius (stars).** R = 0.438 M² + 0.479 M + 0.075 for 0.18 to 1.5 M_sun (Eker et al. 2018); above 1.5 the radius continues as R(1.5) (M / 1.5)^0.555.
- **F3 Luminosity history.** L(t) = L_zams / (1 − 0.286 (t / t_ref)) with t_ref = 0.457 x main-sequence lifetime, which reproduces the Sun at 71% of present luminosity at zero age (Gough 1981 form, rescaled by lifetime).
- **F4 XUV history.** F_xuv ∝ t^−1.23 for t above 0.1 Gyr, saturated before (Ribas et al. 2005), normalised to the present Sun.
- **F5 Kepler.** P = 2π sqrt(a³ / (G (M_primary + m))).
- **F6 Gravity.** g = G m / R².
- **F7 Escape velocity.** v_esc = sqrt(2 G m / R).
- **F8 Bulk composition.** Given (m, R), find the (core, rock, ice) mass fractions whose three-layer model reproduces R, using mass–radius curves for iron, silicate and water-ice end members (Zeng et al. 2016 fits). Two unknowns from one constraint: assume a chondritic core-to-rock ratio of 0.45 unless density exceeds pure rock, in which case ice is zero and the core fraction is solved for.
- **F9 Oblateness.** f = (5/4) Ω² R³ / (G m) x correction for the moment of inertia (Darwin–Radau).
- **F10 Locking timescale.** t_lock = ω a⁶ I Q / (3 G M_p² k₂ R⁵), I = α m R² (Gladman et al. 1996).
- **F11 Rotation regime.** L_Ro = sqrt(sqrt(g H) / (2 Ω / R))... computed as the equatorial Rossby deformation radius L = sqrt(c R / (2 Ω)) with c = sqrt(g H), H the scale height. rapid if L below R/3; slow if L above R; rhines between (Haqq-Misra et al. 2018, simplified).
- **F12 Radiogenic heating.** H(τ) = H_0 Σ_i f_i exp(λ_i τ) for τ = time before present, with H_0 = 5e−12 W per kg of silicate now and (f_i, half-life Gyr): U238 (0.37, 4.47), U235 (0.02, 0.70), Th232 (0.43, 14.0), K40 (0.18, 1.25) (Turcotte and Schubert).
- **F13 Tidal heating.** Ė = (21/2) (k₂/Q) G M_p² R⁵ n e² / a⁶ (Peale, Cassen and Reynolds 1979).
- **F14 Lithosphere thickness.** z = k (T_mantle − T_surface) / q, with k = 3 W/m/K for rock and k = 651 / T W/m/K for ice.
- **F15 Maximum relief.** h_max = c Δσ / (ρ_crust g), c calibrated so Earth gives 9 km, Mars 22 km and Venus 11 km with the strength law of interior.crust_strength. For ice, h_max is the smaller of the strength value and the relaxation-limited value for the body's age.
- **F16 Global contraction.** ΔR = (α_v / 3) ΔT R with α_v = 3e−5 per K and ΔT = 50 K x min(1, (R / R_Mars)²) of mantle cooling plus 250 K times the squared iron index from the uncompressed density (0 at 3.3 g/cm³, 1 at 5.15), accumulating as 1 − exp(−(t − 0.7 Gyr) / 2 Gyr) after a 0.7 Gyr onset during which early heating offsets the cooling. Mercury 5.7 km, Mars 1.7 km, the Moon 0.2 km (Byrne et al. 2014; Watters et al. 2010; Tosi et al. 2013).
- **F17 Magma production.** Rate = max(0, q − q_conductive) x 4π R² / (ρ L_fusion + ρ c_p ΔT) x efficiency, calibrated so Earth gives about 20 km³/yr and Io about 1 cm/yr of resurfacing.
- **F18 Ice shell thickness.** d = k_ice (T_melt − T_surface) / q with k_ice = 651 / T evaluated at the mean shell temperature; T_melt lowered by ammonia if present.
- **F19 Tectonic regime.** ice_shell_ocean if has_subsurface_ocean; ice_shell_rock if composition is ice_dominated or rock_ice without an ocean; mobile_lid if surface liquid water exists and surface_heat_flux is above 0.05 W/m² and radius above 4000 km; episodic_lid if rock with surface_heat_flux above 0.04 and no surface water and radius above 5000 km; stagnant_lid otherwise. Overridable.
- **F20 Scale height.** H = k_B T / (μ m_u g).
- **F21 Greenhouse optical depth.** τ = Σ_species a_s u_s^b_s with u_s the species column mass, a_s and b_s fit to Earth (+33 K), Mars (+5 K), Venus (+500 K) and Titan (+3 K net); water column derived from relative humidity and temperature by Clausius–Clapeyron; N2–N2 and CH4–N2 collision-induced terms scale as column squared. Surface temperature follows T_s⁴ = T_eq⁴ (1 + 3τ/4) in the zero-dimensional estimate.
- **F22 Equilibrium temperature.** T_eq = (L (1 − A_B) / (16 π σ d²))^(1/4); substellar temperature without redistribution is √2 T_eq.
- **F23 Jeans parameter.** λ = G m μ_s m_u / (k_B T R) per species, evaluated at the exobase temperature T_exo = T_s + 700 K x xuv^0.3 x (1 − 0.8 f_CO2) (provisional). A species is flagged below λ = 36, which is the rule of thumb that a gas survives when the escape speed exceeds six times its most probable thermal speed.
- **F24 Crater cut-off.** D_min = c (P/g)^0.6, c such that Earth gives 150 m; this yields about 2.4 km for Venus and 13 m for Mars, to be re-fit against Venus and Titan counts.
- **F25 Vapour pressure.** ln(p/p_t) = (L / R_s)(1/T_t − 1/T) per species about its triple point (T_t, p_t) with latent heat L of vaporisation above T_t and sublimation below.
- **F26 Sea level.** Find z such that ∫ max(0, z − h) dA over the surface equals the liquid inventory volume, by bisection on the hypsometric curve.
- **F27 Plate count.** N = round(7 (A / A_Earth)^0.5 (q / 0.087)^0.5), minimum 2.
- **F28 Impact flux.** Relative to the Moon: inner system by heliocentric distance and gravitational focusing per Ivanov 2001 (Mars about 2); outer system per Zahnle et al. 2003 by planet and satellite distance.
- **F29 Impact velocity.** v² = v_∞² + v_esc², with v_∞ by location: 0.4 to 1.0 x the local orbital speed for heliocentric impactors, plus parent-planet focusing for moons.
- **F30 Simple-to-complex transition.** D_sc = 16.5 km x (1.62 / g) x f_material, f = 1 rock, 0.15 ice.
- **F31 Flux history.** The Neukum chronology N(≥1 km) = 5.44e−14 (e^(6.93 T) − 1) + 8.38e−4 T per km², T in Gyr before present, rewritten in time since formation τ with the present 4.5 Gyr after the Moon formed: rate = c [1 + ρ e^(−τ/τ_d)] with ρ = a b e^(4.5 b) / c ≈ 1.57e4 and τ_d = 1 / b = 0.144 Gyr, zero before the surface solidifies. A body of any age gets its heavy bombardment at formation and the steady flux ever after; the production function is the Neukum 2001 polynomial.
- **F32 Hadley edge.** φ_H = sqrt(5 g H Δ_h / (3 Ω² R²)) in radians, capped at π/2, with H the tropopause height (taken as 2 scale heights, which gives Earth 32°) and Δ_h the fractional equator-to-pole radiative-equilibrium temperature contrast (Held and Hou 1980).
- **F33 Heat transport coefficient.** D = 0.58 W/m²/K x (P / 1 bar)(c_p / c_p,air)(μ_air / μ)(Ω_Earth / Ω)², capped at 100 x the Earth value (Williams and Kasting 1997 scaling as used in exoplanet energy-balance models).
- **F35 Scarp population.** Lengths from N(>L) ∝ L^−1.5 between 10 km and L_max = min(0.4 R, 1,000 km x (ε / 0.005)²) for areal strain ε; height h = 0.03 L_km^0.667 km (Enterprise Rupes), shortening h / tan 25° per unit length, the population drawn until the shortening area ε x 4π R² is spent; a face of 12° central slope and a back slope 20 h long.
- **F36 Flood volcanism.** Eruptive fraction of F17 = 0.1 / (1 + exp((z_lith − 400 km) / 30 km)) x exp(−(ε / 0.0012)²); eruptions of 100 to 10,000 km³ (log-uniform) from vents weighted by exp(−(crust − min) / 10 km) exp(−(h − median) / 2 km) x (1 + 20 x rise height fraction) x the magmatic province weight (volcanism.province_count); each fills one lake from the nearest low, taking cells in height order across passes, spreading no thinner than 50 m and no further than 1,500 km, subsiding by 5% of its thickness.
- **F41 Tidal stress of a thin shell.** The permanent bulge H = h₂ (M_p / M)(R / a)³ R towards the planet; membrane strains ε_θθ = (L P₂'' + H P₂) / R and ε_φφ = (L cot θ P₂' + H P₂) / R about the bulge axis with L = 0.25 H, plane stress with E = 9 GPa and ν = 0.33. The diurnal field at mean anomaly M is the field of a bulge (1 + 3e cos M) times the permanent one about an axis librated by 2e sin M, less the permanent field; non-synchronous rotation adds the field of the axis advanced by 1° less the permanent field. Cracks open perpendicular to the most tensile principal stress where it exceeds 40 kPa.
- **F42 Lineae.** The population carries π density / 2 of crack length per unit area, drawn per time bin within the retention age, lengths log-uniform 0.064 to 0.96 radii (Europa 100 to 1,500 km); a crack nucleates in proportion to the tension above the threshold, runs along the field of its time where the slow field dominates, and where the diurnal tension is more than 60% of the whole it propagates at 0.8 m/s only while the tension lasts and turns with the field through the orbit, arcing between cusps; 15% of the ridges pull apart into bands 0.0064 to 0.019 radii wide (Europa 10 to 30 km). The ocean's salts in a crack fade under frost over 0.03 Gyr.
- **F43 Focused heating.** The shell under one pole is thinner than the mean by the polar thinning times a cosine taper that reaches the mean 70° from the pole, the edge wobbling by 15% with a two-cycle field; the replayed craters' survival and the crack density follow the local thickness (the density as 1 / d); the loop's relaxation keeps the mean shell's rates, a soft zone given its own band factors beside a stiff one sinking many times its offset. Where the taper stands above a half, within 35° of the pole, the terrain is active once the thinning passes 0.3: it floats 0.6 of its isostatic depression low, which the relaxation treats as the shell's equilibrium surface rather than a load; at every step it keeps only the wavelengths longer than the renewal length (100 km or R / 20), so that it carries no craters, and its cracks are drawn within the last 2 Myr; and the finish lays folds of 100 m at 1.5 km along the cracks' direction and a 400 m margin ridge 8 km wide at its edge. A shell thinner than 0.02 R renews itself everywhere within the retention age instead; a thicker one keeps its old craters.
- **F44 Grooved terrain.** Lanes of new ice are parts of great circles, lengths log-uniform 0.2 to 2 rad, widths 20 to 200 km (at most 0.3 R) wobbling by a quarter along the lane, their midpoints gathered by rejection into belts of extension drawn by a 1.5-cycle field, three in five laid beside an earlier lane turned by up to 12° and one to three widths to its side so that they come in parallel sets, drawn until they cover the resurfaced fraction, at most 1,500; a fifth of them are smooth bright terrain. At the resurfacing age the relief inside is smoothed over 100 km (or R / 20) and dropped 200 m, the surface age resets and the unit becomes grooved terrain, clean ice whatever the regolith around it; the finish redraws the lanes' edges at its own resolution, lays the floor back and cuts grooves of 150 m along each lane at 7 km (or three cells), fading over the outer 15% of the width and turning in sets 300 km long. Craters formed since cut the grooves; craters formed before are gone.
- **F39 Viscous relaxation of an ice shell.** A load of wavelength λ relaxes on τ = 4π η / (ρ g λ) with η the viscosity of the ice at depth 0.16 λ into the shell, along the conductive profile T(z) = T_s (T_m / T_s)^(z/d), η = 1e14 Pa s x exp((60 kJ/mol / R)(1/T − 1/T_m)) x 10^(12 f) for a crust with non-ice share f; loads shorter than 4 times the elastic thickness ride on the lid and do not relax, and the share that relaxes rises to 1 at twice that. The band low-passes are area averages (a sampled pyramid let one cell stand for a patch). Where the shell is thinned (F43) the topography relaxes towards the thinned shell's equilibrium surface rather than treating the depression as a load, and a renewed terrain drops the wavelengths shorter than its renewal length and relaxes the rest with its neighbours. The loop applies the survival per wavelength band at every step with the shell as the thermal history has it; replayed craters relax by the same law over their lives.
- **F40 Sublimation lag and frost.** Ice sublimes at the latitude's mean temperature T = (S (1 − A) cos φ / π σ)^(1/4) at the Hertz-Knudsen rate of its vapour pressure (Murphy and Koop 2005); the non-ice share of the regolith stays as lag, which chokes its growth to the square root of time and is opaque past 2,000 kg/m²; clean excavated ice is gardened back in over 0.5 Gyr; the ice the equator lost condenses where the mean temperature is below 105 K, covering up to 45% there.
- **F38 Rift population.** On each plume rise the extension E = tectonics.rift_extension x (rise height / nominal height) is spent along a transect 0.7 x 2π a long as Σ heave_i x L_i, the troughs drawn from N(>L) ∝ L^−2 between 50 km and 1.5 rise diameters with depth 1.5e−3 L (log-normal scatter 0.4, capped at the uplift), width 30 D and heave 1.2 D; half of them join an existing trough chosen by length as an en-echelon neighbour 0.8 to 3.5 widths to its side. Troughs wider than two coarse cells are carved through the loop over 1 Gyr, the rest cut at the finish; the rifting follows the dome on a stagnant lid and the resurfacing event on an episodic lid.
- **F37 Steady-flux production.** After the heavy bombardment the lunar production function falls two cumulative slopes steeper above 150 km and makes nothing over 600 km; each crater draws from the early or steady population by the early share r e^(−t/τ) / (1 + r e^(−t/τ)) of the flux at its time (F31).
- **F34 Threshold friction speed.** u*_t = 0.1 sqrt(ρ_particle g d / ρ_air) for sand-size grains (Iversen and White 1982, large-particle limit); d = 200 µm default.
