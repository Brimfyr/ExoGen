# Parameter graph specification

Revision 0.2, 2026-10-10. This document is the contract for the `params` module of `exogen-core`. Every node the user can see or set is listed here with its unit, kind, default or derivation, and the source of the relation. The formula appendix at the end gives the relations by number. Values marked *fit* are calibrated against the reference bodies during the phase that implements the stage; the derivation form is fixed here, the constants are not.

## Implementation status

Phase 1 implements the graph mechanism and the groups star, orbit, body, rotation, interior, atmosphere, volatiles, climate (zero-dimensional subset), cratering (scaling subset), detail, output, seeds and multipliers, plus `tectonics.dichotomy_mode`, `tectonics.dichotomy_amplitude` and `tectonics.initial_crust_thickness`. The executable table is `crates/exogen-core/src/params/registry.rs`; every constant lives in `crates/exogen-core/src/constants.rs` with its source, and constants marked provisional there are calibrated in the phase named. The remaining tectonics, volcanism, weather, surface-process and material nodes are added by the phases that implement those stages. A test asserts that every derivation reads only the nodes it declares, and `exogen validate` checks the derived values of the fourteen reference presets against measured values.

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
| interior.dissipation_q | 1 | expert derived | 100 rock, 30 ice, 12 with surface ocean | |
| interior.radiogenic_heat_per_kg | W/kg | derived | F12 at body.age, times rock mass fraction | Turcotte and Schubert |
| interior.radiogenic_flux | W/m² | derived | radiogenic power / (4π R²) | |
| interior.tidal_power | W | derived | F13 | Peale, Cassen and Reynolds 1979 |
| interior.tidal_flux | W/m² | derived | tidal_power / (4π R²) | |
| interior.secular_cooling_fraction | 1 | expert derived | 1.0 for R above 3000 km falling to 0.2 at 500 km, times radiogenic flux | fit |
| interior.surface_heat_flux | W/m² | derived | radiogenic plus tidal plus secular | |
| interior.heat_flux_history | W/m² | curve | F12 and F13 evaluated over body age with the eccentricity curve | |
| interior.mantle_temperature | K | expert derived | Thermal evolution, calibrated so Earth is 1600 K now | fit |
| interior.lithosphere_thickness | m | derived | F14 conductive thickness | |
| interior.elastic_thickness | m | derived | Depth to the 650 K isotherm for rock, 170 K for ice, from the conductive profile | |
| interior.crust_strength | Pa | derived | 150 MPa rock, weakened linearly to 50 MPa as surface temperature approaches 750 K; 5 MPa ice | fit |
| interior.max_relief | m | derived | F15 | |
| interior.contraction_total | m | derived | F16 | |
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
| tectonics.contraction_strain | 1 | derived | 2 x contraction_total / radius | |
| tectonics.scarp_total_shortening | m | derived | contraction_strain x 2π R | |
| tectonics.plume_rise_count | 1 | derived | 0 mobile lid; stagnant lid: 2 when surface_heat_flux above 0.02, else 0 | heuristic |
| tectonics.plume_rise_diameter | m | derived | 0.75 x radius | fit |
| tectonics.plume_rise_height | m | derived | 0.3 x max_relief | fit |
| tectonics.canyon_system_count | 1 | derived | 1 per plume rise with diameter above 0.5 R | |
| tectonics.resurfacing_event_age | Gyr | derived | 0.5 for episodic lid | |
| tectonics.tessera_fraction | 1 | derived | 0.08 episodic lid | |
| tectonics.ice_resurfaced_fraction | 1 | derived | ice_shell_ocean: clipped (tidal_flux / 0.05); ice_shell_rock: 0 | fit |
| tectonics.ice_crack_density | per m | derived | Tidal stress amplitude / ice strength, scaled by shell thickness | fit |
| tectonics.chaos_fraction | 1 | derived | Fraction of the shell thinner than 0.5 x mean, from the heat flux pattern | |
| tectonics.focused_heating_pole | enum none, north, south | derived | south when tidal_flux is above 0.1 and radius is below 400 km | heuristic |
| tectonics.nonsynchronous_rotation | 1 | expert | 0 | |

## volcanism

| id | unit | kind | default or derivation | source |
|---|---|---|---|---|
| volcanism.rate | m³/s | derived | interior.magma_production_rate | |
| volcanism.resurfacing_rate | m/Myr | derived | rate / area | |
| volcanism.style | enum none, silicate, sulfur_silicate, cryo_water_ammonia, cryo_brine | derived | By composition class and heat flux; sulfur_silicate above 1 W/m² on rock | |
| volcanism.shield_count | 1 | derived | plume_rise_count x 4 for stagnant lid; hotspot_count for mobile lid | heuristic |
| volcanism.arc_volcano_density | per m of arc | derived | 1 per 70 km of convergent boundary | |
| volcanism.flood_episode_count | 1 | derived | From the heat flux history: one per 10 Myr of rate above threshold, capped | fit |
| volcanism.paterae_count | 1 | derived | 400 x (tidal_flux / 2) for sulfur_silicate | heuristic |
| volcanism.fresh_lava_albedo | 1 | derived | 0.07 basalt, 0.4 cryolava, 0.1 sulfur-coated | |
| volcanism.active_vent_count | 1 | derived | Vents younger than 1 Myr at the output epoch | |
| volcanism.history | m³/s | curve | Follows interior.heat_flux_history | |

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
| cratering.flux_history | 1 | curve | Chronology derivative F31 | Neukum et al. 2001 |
| cratering.late_spike | bool | root | false | |
| cratering.giant_basins | list of diameter, lat, lon, age | root | empty: sampled from the production function | |
| cratering.diffusivity | m²/Myr | derived | 5.5 x (gravity / 1.62)^0.5, x 10 with an atmosphere above 1 kPa | Fassett and Thomson 2014; fit |
| cratering.ray_fading_timescale | Gyr | derived | 1.0 x (a / 1 AU)² airless; 0.05 with an atmosphere | fit |
| cratering.rampart_ejecta | bool | derived | Regolith volatile content above threshold | |
| cratering.secondary_fraction | 1 | expert | 0.2 | |
| cratering.oblique_fraction | 1 | expert | 0.05 below 15° | |

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
| surface.space_weathering_timescale | Gyr | derived | 1.0 x (a / 1 AU)² airless; disabled above 100 Pa | fit |
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
- **F16 Global contraction.** ΔR = (α_v / 3) ΔT R with α_v = 3e−5 per K and ΔT the mean interior cooling from the thermal evolution (Mercury gives 5 to 7 km).
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
- **F31 Chronology.** N(≥1 km) = 5.44e−14 (e^(6.93 T) − 1) + 8.38e−4 T per km² with T in Gyr; the flux history is its derivative; the production function is the Neukum 2001 polynomial.
- **F32 Hadley edge.** φ_H = sqrt(5 g H Δ_h / (3 Ω² R²)) in radians, capped at π/2, with H the tropopause height (taken as 2 scale heights, which gives Earth 32°) and Δ_h the fractional equator-to-pole radiative-equilibrium temperature contrast (Held and Hou 1980).
- **F33 Heat transport coefficient.** D = 0.58 W/m²/K x (P / 1 bar)(c_p / c_p,air)(μ_air / μ)(Ω_Earth / Ω)², capped at 100 x the Earth value (Williams and Kasting 1997 scaling as used in exoplanet energy-balance models).
- **F34 Threshold friction speed.** u*_t = 0.1 sqrt(ρ_particle g d / ρ_air) for sand-size grains (Iversen and White 1982, large-particle limit); d = 200 µm default.
