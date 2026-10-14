# Audit against principles 13 to 15

Principles 13 to 15 went into PLAN.md on 2026-10-14, after most of the model was built. This note lists what does not meet them yet, so that the work can be planned and each entry closed by a milestone. The ice-world milestones close the entries they touch; the rest are closed when the earlier phases are revisited before Phase 5 (PLAN.md, Phase 4F). An entry closes with a test along a sweep of bodies that shows the feature coming in by degrees, or with a named physical threshold and its width along the quantity that controls it.

- **13, transitions are terrain.** A contact is made by the younger terrain's own process at part strength. Its width is in metres, set by that process, never in grid cells. It is sharp only where the process is sharp.
- **14, features from fields.** Where, when, how large and in what form a feature appears is read from physical fields at its place and time. Events (impacts, eruptions) stay objects. Noise stands only for the unresolved variation of a named quantity.
- **15, in-between worlds.** A calibration fixes a constant of a process, never one body's outcome. A class of feature comes or goes only where a named physical threshold is crossed. Validation carries interpolated bodies.

Status at commit 7088060: every entry is open.

Which reference bodies trigger which switch today (from `exogen derive` on the presets):

| Switch | Bodies |
|---|---|
| Mobile lid | Earth |
| Episodic lid | Venus |
| Io regime | Io |
| Ice shell over an ocean | Europa, Ganymede, Callisto, Titan, Enceladus, Triton |
| Ice shell over rock (cryodomes) | Pluto, Ceres |
| Chaos (thin-shell rule) | Europa |
| Lanes | Ganymede, Enceladus |
| Polar focusing, plumes, fallout | Enceladus |

A rule that one body alone triggers is the first place to look for a per-body tuning.

## Switches (principle 15)

### S1. The tectonic regime is a name

`interior_tectonic_regime` (params/derive.rs) picks one of five names:

- `ice_shell_ocean` or `ice_shell_rock` for an icy composition class, by whether there is an ocean;
- `mobile_lid` with liquid water, a surface heat flux of at least MOBILE_LID_MIN_FLUX (0.05 W/m²) and a radius of at least MOBILE_LID_MIN_RADIUS (4,000 km);
- `episodic_lid` from EPISODIC_LID_MIN_FLUX (0.04 W/m²) and EPISODIC_LID_MIN_RADIUS (5,000 km);
- `stagnant_lid` otherwise.

Eighteen derivations compare the name as a string: the retention age, the magma production rate, the surface heat flux, the bedrock, deep bedrock and lava materials, the shares of chaos, lanes and tesserae, the dichotomy amplitude, the crack density, the plume count, the polar thinning, the resurfacing age, the rise diameter, and the counts of coronae, cryodomes and domes. So do the cryo, episodic, impacts, plume, primordial crust, tectonics, volcanism and thermal stages.

A body that crosses a cut changes all of these at once. A Venus a little under 5,000 km has no tesserae and no coronae; a little over, it has both, its retention age falls from 3.5 to 0.5 Gyr and its magma production rises fivefold, all in one step. The cuts stand for real changes (an ocean, a lid that yields), but they are drawn on proxies (a radius, a heat flux) and have no width.

Direction: keep the name as a label for reports. Each derivation reads the continuous quantity the name stands for: the lid's convective stress against its yield stress for a mobile or episodic lid (Moresi and Solomatov 1998), the ocean's depth for an ice shell. A derivation that must switch names its threshold and ramps over a stated width of that quantity.

When: the earlier-phases revisit; the ice-shell names with Europa.

### S2. The composition class

`composition_class` (params/derive.rs) cuts the core and ice mass fractions into four names at IRON_RICH_CORE_FRACTION, ROCK_MAX_ICE_FRACTION and ROCK_ICE_MAX_ICE_FRACTION. The derivations read the name in 23 places. Among them are the moment of inertia factor (0.33, 0.35 or 0.40 by class), the bedrock and lava materials, the Io regime, and whether the body is icy at all (`is_icy`), which decides between an ice shell and a rock lid.

The fractions are continuous inputs. A body just under the ice cut has no ice shell however much ice it holds, and its moment of inertia factor jumps by a fifth across a cut.

Direction: the moment of inertia from the layered density profile the composition already implies; an ice shell wherever the ice layer has a thickness (`ice_layer_thickness`), not by class; materials by process (S7).

When: the earlier-phases revisit.

### S3. The Io regime

`is_io_regime` (stages/io.rs) and `io_regime` (params/derive.rs): a rock body with a tidal flux above IO_REGIME_TIDAL_FLUX (0.5 W/m²). It turns on paterae (IO_PATERAE_PER_TW), mountains, plume rings, the sulfur-coated bedrock, the burial term of the retention age and the Io stages, and turns off mare basalt. A rock body at 0.4 W/m² gets none of Io's features; one at 0.6 W/m² gets all of them.

Direction: paterae from the heat-pipe share of the heat flux, the part carried up by eruptions rather than conducted (Moore and Webb 2013), which grows with the flux; mountains from the compressive stress that burial and subsidence put on the lithosphere against its strength (Schenk and Bulmer 1998; McKinnon et al. 2001); the sulfur coating from the plume deposition rate against the weathering time. The burial term already scales with the magma rate and needs no threshold.

When: the Io milestone.

### S4. The thin-shell rule

WHOLESALE_SHELL_RATIO (a shell under 0.02 of the radius) is read in `tectonics_ice_resurfaced_fraction` (no lanes under it), `tectonics_chaos_fraction` (no chaos over it), `cratering_retention_age` (wholesale renewal under it, with S6) and the impacts stage.

Measured (2026-10-14, `lanes_come_in_by_degrees_as_the_lid_thins_against_the_body`): an Enceladus half again as large has a shell 0.035 of its radius and makes 81 lanes; one twice as large has 0.010 and makes none, and chaos takes over in the same step. Europa itself has both bands, which are lanes that spread, and chaos, so the reference body contradicts the switch.

Direction: the share of the renewal that goes to chaos from the processes that make chaos (melt-through and brine lenses in a thin or warm shell, Schmidt et al. 2011; diapirs in a convecting shell, Pappalardo et al. 1998), and the share that goes to lanes from the lid's extension. Both rise and fall with the shell's state, so a body between Ganymede and Europa has some of each.

When: the Europa milestone.

### S5. Active thinning

POLAR_ACTIVE_THINNING (0.3 of the shell) is read in `tectonics_focusing_period`, `volcanism_plume_fallout_rate` and the ice stage's activity test. Under it a body has no focusing episodes, no plumes and no fallout; over it, all three. Enceladus (0.74) is the only reference body over it; Europa is at 0.12.

The threshold is a share of the shell, not a physical event.

Direction: the event is a crack that reaches the ocean. Tidal stress opens a dry crack only so deep (ICE_CRACK_STRESS over the weight of the ice, about 380 m on Enceladus), and a water-filled crack rises by its own buoyancy (Crawford and Stevenson 1988). Episodes and plumes begin where the thinned shell lets a crack through, and their strength grows from there.

When: the Europa milestone, whose plumes are the test.

### S6. Wholesale renewal of a thin shell

`cratering_retention_age` gives a thin shell (S4) over ICE_RESURFACING_TIDAL_FLUX (0.01 W/m²) a retention of RETENTION_ICE_OCEAN_GYR, and any other shell the body's age: a step in the crater count.

Direction: the retention age from the area that lanes and chaos renew per unit time, so that it follows the same quantities they do.

When: the Europa milestone.

### S7. Materials by name

`materials_bedrock`, `materials_lava` and `materials_deep_bedrock` choose by class and regime: felsic under a mobile lid, anorthosite on an airless stagnant-lid rock body, sulfur-coated over the Io flux, water ice or dirty water ice by regime and airlessness.

Direction: materials from processes: anorthosite from flotation in a magma ocean deep enough for plagioclase to crystallize, felsic crust from recycling with water, the coating from deposition against weathering.

When: the earlier-phases revisit, with Phase 5's materials work.

### S8. The lanes' share along a sweep in radius

Measured (same test as S4): the share renewed as lanes is 0.073 for Enceladus, 0.058 at 1.2 times the radius, 0.310 at 1.5 times and 0 at twice (S4). The dip between the first two needs a physical reason or a fix. The share is the renewal of the peak tidal flux times the part the polar focusing does not take, and both change with the shell's ratio to the radius.

When: the ring-systems milestone, since it needs only a reading of `tectonics_ice_resurfaced_fraction`.

## Widths and scales in grid cells (principle 13)

### W1. Flood margins

FLOOD_MARGIN_LOBE_CELLS sets the lobes on a flood's outline as noise two loop cells long, so the lobes change with the quality and are noise where they should be flow.

Direction: lobes from the flow, a lobe's width from the lava's yield strength and the flow's thickness (Hulme 1974), as planned for Io's branching lobes.

When: the Io milestone for Io; the earlier-phases revisit for lava plains elsewhere.

### W2. The roughness blend

ROUGHNESS_BLEND_PASSES (12 passes of the one-cell kernel, about two cells) makes the edge of a young deposit on old ground a ramp of roughness two cells wide whatever the cell.

Direction: each terrain's roughness follows its own age, and the contact's width comes from the process that laid the younger terrain, as the regolith mixing already does with a width in metres (REGOLITH_MIXING_LENGTH).

When: the earlier-phases revisit.

### W3. Tessera outlines

`TESSERA_SMOOTHING_PASSES` (20 passes on the loop grid) smooths the old crust's height before the tesserae are chosen as its highest share. The smoothing is about two cells, so on Venus it is about 80 km at standard quality and 170 km at preview: the tesserae's outlines change with the quality.

Direction: smooth over a physical length, such as the lithosphere's flexural wavelength or the scale of the mantle's convection. The tessera texture (`tessera_texture`, stages/episodic.rs) is noise and belongs under T1.

When: the earlier-phases revisit.

### W4. The cap on regolith mixing

REGOLITH_MIXING_MAX_PASSES (24) caps the blur that carries the mixing zone's physical width. The mixing runs on the finish grid, and at high quality the cap holds the Moon's zone to about 3.3 km of the 5.1 km its age gives. Minor; the cap exists to bound the cost.

Direction: a blur whose cost does not grow with the pass count (a multigrid or recursive Gaussian), so that no cap is needed.

When: the earlier-phases revisit.

Floors that are not entries: GROOVE_EDGE_CELLS, GROOVE_MIN_CELLS, PIT_MIN_CELLS, EJECTA_LOBE_MIN_CELLS, FINISHING_MIN_CELLS, LOOP_THRESHOLD_CELLS, SCARP_LOOP_MIN_CELLS and RIFT_LOOP_WIDTH_CELLS decide where a feature is drawn (in the loop or the finish) or keep it from aliasing. None of them decides whether the feature forms or how wide it is. MELT_POND_CELLS counts cells of a lattice per crater radius, not grid cells.

## Stamped shapes (principle 14)

### T1. Shapes from a profile and noise

These features take a fixed profile, with noise for variety, rather than a form read from fields at their place and time:

| Feature | Where | First idea of the field |
|---|---|---|
| Chaos | `chaos_population`, `chaos_compute` (stages/cryo.rs) | outline from the thin or warm patch of the shell and a brine lens under it (Schmidt et al. 2011); matrix and rafts from the lens's collapse |
| Cryodomes | `cryodome_population`, `cryodome_compute` (stages/cryo.rs) | form from the extruded volume and the viscosity of the cryolava (Ahuna Mons; Ruesch et al. 2016) |
| Coronae | `corona_core` (stages/episodic.rs) | a plume head under a lid: uplift, then a moat and rim as it spreads and sinks |
| Tessera texture | `tessera_texture` (stages/episodic.rs) | fold and fault fabric oriented by the stress of the event |
| Rises | `dome` (stages/plume.rs) | a plume's buoyancy under the lid's flexure |
| Shields | `shield_height` (stages/plume.rs) | built from flows, each laid down where the last left the slope |
| Canyons | `trough_depth_at` (stages/plume.rs) | strike from the stress field, depth and width from the extension and the lid |
| Io's mountains | `block_radius`, `block_profile` (stages/io.rs) | thrust blocks on faults across a compressive field, their outline noise now |
| Io's plume rings | `ring_profile` (stages/io.rs) | ballistic deposit from the vent speed, which sets the ring's radius |
| Paterae | `patera_compute` (stages/io.rs) | outline from the faults of the block it sank into, a polygon share with noise now |
| Crater walls and rims | `fresh_height` (crater/morphology.rs), `RimShape` (crater/stamp.rs) | one radial profile read at a radius stretched by a crenulation of the outline, so the floor's edge, the wall and the rim are scaled copies of the outline and the wall is as wide all round; the wall should come from its collapse: slump blocks on arcuate faults, each with its own arc, width and drop, terraces that come and go around the crater, the floor's edge at their toes and a rim crest scalloped where large slumps took it back |

Impact craters are events and stay objects. Their sizes follow scaling laws, but their walls do not yet (the last row).

When: crater walls with the ring systems of large impacts, the next ice milestone; chaos with Europa; Io's features with Io; the rest in the earlier-phases revisit.

## Calibrations on one body (principle 15)

Each constant here was fitted so that one body comes out right. Each should become, or be shown to be, a constant of its process, with a check on what it gives a body between the references.

| Constant | Fitted to | Note |
|---|---|---|
| ICE_RESURFACING_FLUX, ICE_RESURFACING_FLUX_FLOOR | Ganymede: two thirds renewed, the shell left near 40 km | a constant of renewal per excess flux, or the strain an episode's extension opens |
| FOCUSING_PERIOD_FACTOR | Enceladus: episodes 0.53 Gyr apart, 0.40 renewed | O'Neill and Nimmo 2010 give 0.1 to 1 Gyr from the shell's convection |
| POLAR_FOCUS_FLUX, POLAR_THINNING_MAX | Enceladus | how strongly tidal heating focuses at a pole |
| PLUME_FALLOUT_RATE | Enceladus: about 0.1 mm/yr near the stripes | scale by the plume's mass flux, which the vent sets |
| CHAOS_FRACTION_MAX, CHAOS_FLUX | Europa's share of chaos | replaced under S4 |
| IO_PATERAE_PER_TW | Io's count over its tidal power | per heat-pipe power once S3 is done |
| GROOVE_RIFT_LID_SHARE | Ganymede's 7 km grooves | a process constant, calibrated once; check on a second body's measured spacing |
| CENTRAL_PEAK_HEIGHT | Copernicus | check how it scales with gravity |

Constants marked provisional for Phases 5 and 7 (atmosphere, albedo, imagery) are left to those phases.

## Validation (principle 15)

### V1. Interpolated bodies in validation

`exogen validate` checks the reference bodies only. Principle 15 asks for bodies between them as well. The first sweeps to carry:

- an Enceladus scaled in radius (in the test suite now, S4 and S8);
- a rock body from the Moon's tidal flux to Io's (S3);
- an ice body from Ganymede's shell to Europa's (S4 to S6);
- a body from Mars to Venus in radius and heat flux (S1);
- a surface from dry to wet at Venus's and Earth's size (S1).

Each sweep prints the count or share of each feature class along the path and fails on a step where no named threshold is crossed.

When: each milestone adds the sweep for what it touches.
