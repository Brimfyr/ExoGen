# Audit against principles 13 to 15

Principles 13 to 15 went into PLAN.md on 2026-10-14, after most of the model was built. This note lists what does not meet them yet, so that the work can be planned and each entry closed by a milestone. The ice-world milestones close the entries they touch; the rest are closed when the earlier phases are revisited before Phase 5 (PLAN.md, Phase 4F). An entry closes with a test along a sweep of bodies that shows the feature coming in by degrees, or with a named physical threshold and its width along the quantity that controls it.

- **13, transitions are terrain.** A contact is made by the younger terrain's own process at part strength. Its width is in metres, set by that process, never in grid cells. It is sharp only where the process is sharp.
- **14, features from fields.** Where, when, how large and in what form a feature appears is read from physical fields at its place and time. Events (impacts, eruptions) stay objects. Noise stands only for the unresolved variation of a named quantity.
- **15, in-between worlds.** A calibration fixes a constant of a process, never one body's outcome. A class of feature comes or goes only where a named physical threshold is crossed. Validation carries interpolated bodies.

Status: every entry is open unless marked closed, with the milestone that closed it.

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

Since F120 a rise keeps the melt clock of its mantle, which takes each body's pace from its magma production: Venus's rises grow within the half gigayear since its overturn because its mantle makes 4.3 times Mars's melt per unit area, its regime's efficiency (MAGMA_EFFICIENCY_EPISODIC, five times the stagnant lid's) being a constant of the name. The next milestone takes this entry's direction for the episodic lid: a mantle whose heat builds under its lid until the lid fails and overturns, its melting following its own temperature.

Partly closed in the fifty-second milestone (F121): the episodic lid is no longer a cut on proxies. A dry rocky body's mantle stores its heat under a lid that breaks where the convecting sublayer's stress along a cell as deep as the mantle exceeds the lid's friction over its thickness (after Moresi and Solomatov 1998): Venus's breaks 21 times and last 0.61 Ga; Mars's, the Moon's and Mercury's never. A lid broken within the last half of its body's life is episodic, its retention and resurfacing ages the time since. Still open: the plates' rule and the ice shells' names; the Io regime (S3); the stagnant bodies' histories, which keep the equilibrium heat flux until the composition model gives the Moon its rock (S2); the melt efficiency by name, so that Venus's melt clock still runs on the episodic lid's constant.

In the fifty-third milestone (F122) the stagnant bodies' mantles keep their own heat too, the composition model giving the Moon its rock: Mars, the Moon and Mercury cool under lids that grow and thin by the heat at their bases and never break, and Venus's breaks 42 times, last 0.49 Ga. Still open: the plates' rule and the ice shells' names; the Io regime (S3); the melt efficiency by name.

In the fifty-fourth (F123) a mantle melts only while it, or its plumes, stand above the solidus at its stiff lid's base, so that a mantle's volcanism ends as its lid thickens and it cools, whatever its regime: the Moon's by about 1 Ga. The share of its heat that goes into melt is still the regime's.

### S2. The composition class

`composition_class` (params/derive.rs) cuts the core and ice mass fractions into four names at IRON_RICH_CORE_FRACTION, ROCK_MAX_ICE_FRACTION and ROCK_ICE_MAX_ICE_FRACTION. The derivations read the name in 23 places. Among them are the moment of inertia factor (0.33, 0.35 or 0.40 by class), the bedrock and lava materials, the Io regime, and whether the body is icy at all (`is_icy`), which decides between an ice shell and a rock lid.

The fractions are continuous inputs. A body just under the ice cut has no ice shell however much ice it holds, and its moment of inertia factor jumps by a fifth across a cut.

Direction: the moment of inertia from the layered density profile the composition already implies; an ice shell wherever the ice layer has a thickness (`ice_layer_thickness`), not by class; materials by process (S7).

When: the earlier-phases revisit.

The fractions themselves were wrong for the bodies formed inside the snow line: a density below the chondritic mix was taken for ice wherever a body formed, so that the Moon carried 8 percent ice and a 30 percent core. Since the fifty-third milestone (F122) a body takes up ice by the temperature where it formed, and one formed too warm for it lacks iron instead: the Moon is all rock, Mars carries a 16 percent core and no ice. The cuts on the fractions are still open.

### S3. The Io regime

`is_io_regime` (stages/io.rs) and `io_regime` (params/derive.rs): a rock body with a tidal flux above IO_REGIME_TIDAL_FLUX (0.5 W/m²). It turns on paterae (IO_PATERAE_PER_TW), mountains, plume rings, the sulfur-coated bedrock, the burial term of the retention age and the Io stages, and turns off mare basalt. A rock body at 0.4 W/m² gets none of Io's features; one at 0.6 W/m² gets all of them.

Direction: paterae from the heat-pipe share of the heat flux, the part carried up by eruptions rather than conducted (Moore and Webb 2013), which grows with the flux; mountains from the compressive stress that burial and subsidence put on the lithosphere against its strength (Schenk and Bulmer 1998; McKinnon et al. 2001); the sulfur coating from the plume deposition rate against the weathering time. The burial term already scales with the magma rate and needs no threshold.

F84 made the coating emergent: the surface is each eruption's fall against the reversion, sublimation and alteration that age it, so a body with fewer plumes is altered further from its poles. The paterae, the plumes and the mountains still switch on at the threshold.

When: the Io milestone.

### S4. The thin-shell rule

WHOLESALE_SHELL_RATIO (a shell under 0.025 of the radius; 0.02 until Europa's shell came to its measured 30 km, F75, which put it at 0.019) is read in `tectonics_ice_resurfaced_fraction` (no lanes under it), `tectonics_chaos_fraction` (no chaos over it), `cratering_retention_age` (wholesale renewal under it, with S6) and the impacts stage.

Measured (2026-10-14, `lanes_come_in_by_degrees_as_the_lid_thins_against_the_body`): an Enceladus half again as large had a shell 0.035 of its radius and made 81 lanes; one twice as large had 0.010 and made none, and chaos took over in the same step. Since F75 (2026-10-15) a larger Enceladus's core is less porous and its shell heats in its own tides, so the shell stays thick against the body (0.094 of the radius at twice the size, 0.8 of the eccentricity) and the step is not reached in that sweep; it remains for a body whose shell is thin against it. Europa itself has both bands, which are lanes that spread, and chaos, so the reference body contradicts the switch.

Direction: the share of the renewal that goes to chaos from the processes that make chaos (melt-through and brine lenses in a thin or warm shell, Schmidt et al. 2011; diapirs in a convecting shell, Pappalardo et al. 1998), and the share that goes to lanes from the lid's extension. Both rise and fall with the shell's state, so a body between Ganymede and Europa has some of each.

When: the Europa milestone.

### S5. Active thinning

POLAR_ACTIVE_THINNING (0.3 of the shell) is read in `tectonics_focusing_period`, `volcanism_plume_fallout_rate` and the ice stage's activity test. Under it a body has no focusing episodes, no plumes and no fallout; over it, all three. Enceladus (0.76) is the only reference body over it; Europa is at 0.27 (0.12 until F75 brought its tidal flux to 14 mW/m² and its shell to 30 km, which focus more).

The threshold is a share of the shell, not a physical event.

Direction: the event is a crack that reaches the ocean. Tidal stress opens a dry crack only so deep (ICE_CRACK_STRESS over the weight of the ice, about 380 m on Enceladus), and a water-filled crack rises by its own buoyancy (Crawford and Stevenson 1988). Episodes and plumes begin where the thinned shell lets a crack through, and their strength grows from there.

When: the Europa milestone, whose plumes are the test.

### S6. Wholesale renewal of a thin shell

`cratering_retention_age` gives a thin shell (S4) over ICE_RESURFACING_TIDAL_FLUX (0.01 W/m²) a retention of RETENTION_ICE_OCEAN_GYR, and any other shell the body's age: a step in the crater count.

Measured (2026-10-15, the twenty-second milestone): the rule renews no surface. It sets the retention age of the crater record, but the plains between Europa's drawn lineae and chaos keep the coarse loop's first crust, so more than half of the moon's surface carries a surface age of 4.5 Gyr. The sub-grid crater roughness (F54) reads that age and roughens the plains as a saturated ancient surface: the grain the shaded maps show on Europa's plains, against which its young chaos matrix looks smooth.

Since F81 (2026-10-15, the twenty-third milestone) the plains' ages come from what renews them: the fabric of narrow ridges renews what the drawn lineae and the chaos matrix leave of one surface's worth in each retention age, so Europa's surface has a mean age of 97 Myr against its retention of 100 Myr. Still open: the retention age itself is RETENTION_ICE_OCEAN_GYR for a thin shell over the flux threshold and the body's age otherwise, a step, and the fabric's rate closes the renewal to it rather than following from the cracks' own rate of formation.

Direction: the retention age from the rate at which cracks, lanes and chaos form, so that it follows the same quantities they do. The freezing stresses (F82) now bend and cross the cracks but do not yet set how many form: a thickening shell's strain, accommodated by the cracks' opening, is the measure of that rate they offer.

When: the earlier-phases revisit, with the crack population's rate.

### S7. Materials by name

`materials_bedrock`, `materials_lava` and `materials_deep_bedrock` choose by class and regime: felsic under a mobile lid, anorthosite on an airless stagnant-lid rock body, sulfur-coated over the Io flux, water ice or dirty water ice by regime and airlessness.

Direction: materials from processes: anorthosite from flotation in a magma ocean deep enough for plagioclase to crystallize, felsic crust from recycling with water, the coating from deposition against weathering.

When: the earlier-phases revisit, with Phase 5's materials work.

The crust's composition also decides how much of a mantle's magma reaches the surface (F124): an anorthositic crust, lighter than basalt, holds most of it down, as the Moon's does (Head and Wilson 1992), while Mercury's reduced, iron-poor magma ocean floated no such crust (Vander Kaaden and McCubbin 2015). The model gives every body one extrusive share, which leaves Mercury's smooth plains at 18.5 percent of its surface against the 27 they cover; with the crust from the magma ocean, the extrusive share should follow it.

Partly closed in the fifty-sixth milestone (F125): a rocky crust's materials come from its processes, the list holds materials rather than bodies and units, and the crust's buoyancy decides how much magma erupts, Mercury's plains rising to 21.8 percent of its surface. Still open: the icy bodies' materials, and the presets that pin them (Callisto, Ceres, Ganymede).

Closed in the fifty-seventh milestone (F126): an icy crust is the rock and ice the body accreted, its colour the intimate mixture at the crust's rock, clean where its shell renews itself, and its regolith carries the lag its warm latitudes leave; Callisto, Ceres and Ganymede pin nothing (F127 since lets cryolava mature to the ground and the lag hide a crater's floor as it hides its rays). Two stated choices carry it, the grain sizes in the mixture and the share of regolith impacts trade between terrains (REGOLITH_EXCHANGE, from the lunar maria's), and one switch remains, the same as the retention age's (S6): a shell renews its surface everywhere or keeps it, by the tidal flux and the shell's ratio to the radius.

### S8. The lanes' share along a sweep in radius (closed)

Measured (same test as S4): the share renewed as lanes is 0.073 for Enceladus, 0.058 at 1.2 times the radius, 0.310 at 1.5 times and 0 at twice (S4). The dip between the first two needs a physical reason or a fix. The share is the renewal of the peak tidal flux times the part the polar focusing does not take, and both change with the shell's ratio to the radius.

When: the ring-systems milestone, since it needs only a reading of `tectonics_ice_resurfaced_fraction`.

Closed in the eighteenth milestone. The share is the renewal of the peak flux, saturated on Enceladus, times one minus the polar thinning over its maximum, and the thinning is the product of two factors that move apart along the sweep: the flux against POLAR_FOCUS_FLUX, capped at one, and the focus of a shell that is thick against the body, 1 - exp(-ratio / POLAR_FOCUS_SHELL_RATIO). Enceladus stood at 19 mW/m², just under the cap, so the first step of the sweep raised the flux term faster than the thinner shell lowered the focus, and the share dipped before it rose. Since F68 Enceladus stands at 37 mW/m², past the cap, and the share follows the shell's ratio alone: 0.002 at 0.8 of its radius, 0.053 at 1, 0.116 at 1.1, 0.194 at 1.2, 0.277 at 1.3, 0.358 at 1.4 and 0.433 at 1.5, rising with every step until the thin-shell rule's step at 1.7 (S4). A body whose flux is under the cap still has the dip: more heat focuses more until the focusing is full, and a thinner shell focuses less; both are the model's stated dependencies, not a fault in it.

### S10. Tides by composition class

F75 gives an icy moon's shell the tides its rheology dissipates, and F68 a porous interior its own; a rocky body's tides still come from `interior.love_number_k2` and `interior.dissipation_q` by composition class (F13): Io's k₂/Q is 0.012 against the measured 0.015 (Lainey et al. 2009), but a rocky body's dissipation does not change with its temperature, so Io's heat cannot run away or limit itself, and the Laplace resonance's cycles (Bagheri et al. 2026: driven by Io's melt-dependent dissipation) cannot arise from it. The same k₂ and Q set every body's locking timescale. An icy moon's rocky interior below its ocean dissipates nothing unless it is porous.

Direction: a rocky interior's imaginary Love number from its Maxwell or Andrade rheology at its temperature and melt fraction, as F75 does for a shell (Bagheri et al. 2026; Renaud and Henning 2018), so that its tides rise and fall with its heat.

When: the Io milestone.

### S9. Central pits by the target's material

A crater on a target `is_ice` reads has a central pit from 8.5 transitions and a dome in it from 20 (F72); on any other, a peak ring and a multi-ring basin. A crust of rock and ice, Ceres's or Callisto's dirty ice, gets one or the other whole. The same switch picks the transition's ice factor (S7).

Since F135 the same switch picks the largest crater a body holds, as wide as its radius on ice, as wide as itself on rock, and so when a body of rock and ice begins its record.

Direction: the pit from the impact melt's volume against the crater's, the melt-drainage model Alzate and Barlow (2011) and Elder et al. (2012) favour; ice melts under a far lower shock pressure than rock (Pierazzo et al. 1997), so a target's ice share would pit its craters by degrees. Since F126 the crust's share of rock is a node of its own (body.crust_impurity), which such a model would read.

When: the earlier-phases revisit, with S7.

## Widths and scales in grid cells (principle 13)

### W1. Flood margins (closed by F94)

FLOOD_MARGIN_LOBE_CELLS sets the lobes on a flood's outline as noise two loop cells long, so the lobes change with the quality and are noise where they should be flow.

Direction: lobes from the flow, a lobe's width from the lava's yield strength and the flow's thickness (Hulme 1974), as planned for Io's branching lobes.

When: the Io milestone for Io; the earlier-phases revisit for lava plains elsewhere.

On Io, F85 makes a field's margins its breakouts' lobes, each traced over the ground and laid at its own time, so the outline is the flow's. The lobes' sizes are stated shares of the field's reach, not yet the lava's yield strength and supply. Flood margins on other bodies are as before.

Closed in the thirty-first milestone (F94) for the floods that can show at the end: each is drawn as its flows, a trunk, branches and breakouts traced over the ground from where its lava came into its low, so its margin and the gaps it leaves are its lobes', at the finish's resolution, and a mare is a mosaic of flow units. A flood buried under four or more successors at every coarse cell of its low keeps the coarse outline, hidden. Still open: the lobes' sizes are stated shares of the low's reach, not the lava's yield strength and effusion rate (Hulme 1974), and a flood's extent is still the coarse low the loop filled, compact on flat ground, so Imbrium's 1,200 km flows cannot form; the loop's own spreading would have to follow the flows for that. Since F109 the loop's lava leaves its sheets down the slopes it runs over, a large flood's for hundreds of kilometres in a swath as wide as its volume spreads it, which is the loop's spreading following the flows in outline; on flat ground the sheet still ends where its low does.

### W2. The roughness blend

ROUGHNESS_BLEND_PASSES (12 passes of the one-cell kernel, about two cells) makes the edge of a young deposit on old ground a ramp of roughness two cells wide whatever the cell.

Direction: each terrain's roughness follows its own age, and the contact's width comes from the process that laid the younger terrain, as the regolith mixing already does with a width in metres (REGOLITH_MIXING_LENGTH).

When: the earlier-phases revisit.

### W3. Tessera outlines

`TESSERA_SMOOTHING_PASSES` (20 passes on the loop grid) smooths the old crust's height before the tesserae are chosen as its highest share. The smoothing is about two cells, so on Venus it is about 80 km at standard quality and 170 km at preview: the tesserae's outlines change with the quality.

Direction: smooth over a physical length, such as the lithosphere's flexural wavelength or the scale of the mantle's convection. The tessera texture (`tessera_texture`, stages/episodic.rs) is noise and belongs under T1.

When: the earlier-phases revisit.

### W4. The cap on regolith mixing (closed by F95)

REGOLITH_MIXING_MAX_PASSES (24) caps the blur that carries the mixing zone's physical width. The mixing runs on the finish grid, and at high quality the cap holds the Moon's zone to about 3.3 km of the 5.1 km its age gives. Minor; the cap exists to bound the cost.

Direction: a blur whose cost does not grow with the pass count (a multigrid or recursive Gaussian), so that no cap is needed.

When: the earlier-phases revisit.

Closed in the thirty-second milestone (F95): the zones are the levels of a Gaussian pyramid, each cell taking its share between the two levels its own zone lies between, so the zone is its physical width at any quality. Each contact's zone now follows its own age, the younger side's, and the body's gardening (its diffusivity times the impacts' equivalent time), and the lava's tint is mixed with its coverage. Still open: the diffusivity itself scales with gravity alone (5.5 m²/Myr times the square root of g over the Moon's), which gives Mercury 1.5 times the Moon's where Fassett et al. (2017) measure about twice from its simple craters; the impact velocity's part waits on the diffusivity's own revisit.

Floors that are not entries: GROOVE_EDGE_CELLS, GROOVE_MIN_CELLS, PIT_MIN_CELLS, EJECTA_LOBE_MIN_CELLS, FINISHING_MIN_CELLS, LOOP_THRESHOLD_CELLS, SCARP_LOOP_MIN_CELLS and RIFT_LOOP_WIDTH_CELLS decide where a feature is drawn (in the loop or the finish) or keep it from aliasing. None of them decides whether the feature forms or how wide it is. MELT_POND_CELLS counts cells of a lattice per crater radius, not grid cells.

### W5. Sub-grid crests of narrow ridges

A drawn ridge narrower than two cells is drawn two cells wide by the area of its cross-section (F78), but its double-ridge profile keeps two crests a quarter of that width apart, each narrower than a cell. On a diagonal they bead along the ridge in the shaded relief, as the narrow stains did before F81 spread them over at least 1.2 cells by their share.

Direction: a narrow ridge's relief drawn as the cell-averaged cross-section, a single bump of the same area, where its crests are finer than the grid.

When: the earlier-phases revisit.

### W6. Flows finer than the grid

A lobe sets its lava where its cover reaches half a cell (F94, F101), so a flow narrower than about a cell and a half does not show at all, and one that narrows and widens along its path shows only where it is wide. On a shield's flanks, where flows narrow on the steep ground (F101), the flows are finer than the grid at the standard quality and their toes show where they spread; at higher qualities more of the flank flows show.

Direction: a flow's lava blended into its cells by its cover, its composition at least, so that the flows finer than the grid show as a faint striation, as the Tharsis Montes' flanks do at a few kilometres a pixel.

When: the earlier-phases revisit.

### W7. Buried floods' margins in coarse steps (closed by F107)

A flood the finish draws as its flows (F94) has lobate margins, but one buried under four or more successors keeps the outline the loop gave it, cell by coarse cell, and the steps of its thickness stay in the ground under the thin sheets above. On Mars's volcanic plains the hillshade shows them as a staircase of straight edges a loop cell long (21 km), the plains' most visible artefact at a kilometre a pixel.

Direction: a buried flood's relief laid with the same lobate margins as a visible one's, or its coarse steps softened over the lobes' scale, so that the old flows under the young show as lobes too.

When: the earlier-phases revisit, before the Mars dichotomy's lowland burial, which buries with floods.

Diagnosed in the thirty-ninth milestone's session. The steps are not the coarse outlines: a cell's trace through the flood passes (EXOGEN_DEBUG_CELL) shows a flood filling a cell of the ring round its low some 800 m, to its level, while the cell next to it, outside the cells the pass visits, keeps its ground. F95's hollow fill takes any fine ground in the ring below the level for a crater the loop never saw; but where the ring stands high between the low and lower ground beyond, the coarse resampling's ramp toward the lower neighbour carries the fine ground below the level too, an open slope, and the fill ponds it out to the ring's edge, which runs in coarse cells. Sampling the level smoothly between coarse cells changed nothing; draining the open slopes, hollows that still run below the level past the reach of any unseen crater, took the cliffs away but bared the ground the finish pieces together from the coarse smears, blocky at a few fine cells, which the ponds had hidden. Both were taken back. The fix needs that reconstructed ground made smooth first, and then the open slopes drained.

Closed in the fortieth milestone (F107). Gating the ring's fill to the finer relief took the large cliffs away and left walls a fine cell wide. A trace of one wall showed the cause in the loop. A flood's lake took in the ground beyond its pass, down to where the lake would spread thinner than 50 m on average, and filled it all to one level, so the lake stood over the open ground where it stopped: by more than 100 m at a fifth of Mars's floods, and by up to 4.6 km. The finish resampled those steps, the later floods laid thin sheets over them, and the ground the finish pieced together from the coarse smears was blocky because the loop's lakes were. A full lake now spills the rest of its lava into the next low, so every lake ends in a flow's front, and the ring's fill takes only the finer relief. Mars's and Venus's plains are flow fields with lobate margins. The only blockiness left, at a fine cell or two, is graben a few kilometres wide on Mars's rise, partly buried by flows, a matter of their own shapes at the standard grid.

### W8. Deposits drawn by whole cells (closed by F131)

A deposit set where a field passes its level in each cell is a staircase of whole cells at any output resolution. Ceres's faculae were (`if patch > 0.45 { cover = 1 }`), and the user found them pixelated, as lava flows had been before.

Direction: every deposit a field bounds is drawn by the share of each cell its margin covers (`fields::level_share`, F127), the margin a straight edge at the cell's scale from the field's value and gradient. A test holds the faculae to it.

When: each with its process's revisit. Closed for the faculae in the fifty-eighth milestone (F127). Still drawn by whole cells: an overturn's lava and a corona's interior (episodic.rs), a shield's renewed lava (plume.rs), a patera's floor (io.rs), Enceladus's polar terrain (ice_shell.rs, under its plume's frost), and a flood's lobes, which set their lava where they cover half a cell (W6), the regolith's mixing blurring an old flow's margin but not a young one's.

Closed in the sixty-second milestone (F131): the overturn's tessera and plains, a corona's interior, a pancake dome, a shield's flows, pond and edifice, a flood's lobes and outline, a patera's floor and Ceres's cryodomes cover each cell by their share, and a deposit whose cover the albedo reads by its unit takes its unit wherever it lies. Enceladus's polar terrain was already drawn by share in the finish. Tests hold the tessera's contact and the domes' edges to it, with the faculae. The rule stands for every deposit to come.

## Stamped shapes (principle 14)

### T1. Shapes from a profile and noise

These features take a fixed profile, with noise for variety, rather than a form read from fields at their place and time:

| Feature | Where | First idea of the field |
|---|---|---|
| Chaos | `chaos_population`, `chaos_compute` (stages/cryo.rs) | closed by F80, the twenty-second milestone: generations where the warmth the ocean brings up stands above a threshold, structured at every scale from 1,300 km down to the rafts and strongest at low latitudes; the lid breaks into rafts grouped in plates that founder with the warmth's excess, or drift and turn carrying the old surface, in a matrix at the old surface's smoothed level. The warmth does not yet thin the shell that the cracks read (F76) |
| Cryodomes | `cryodome_population`, `cryodome_compute` (stages/cryo.rs) | form from the extruded volume and the viscosity of the cryolava (Ahuna Mons; Ruesch et al. 2016) |
| Coronae | `corona_core` (stages/episodic.rs) | a plume head under a lid: uplift, then a moat and rim as it spreads and sinks |
| Tessera texture | `tessera_texture` (stages/episodic.rs) | fold and fault fabric oriented by the stress of the event |
| Rises | `dome` (stages/plume.rs) | a plume's buoyancy under the lid's flexure |
| Shields | `shield_height`, `shield_flows` (stages/plume.rs) | partly by F101 and F102, the thirty-fifth milestone: the surface is the eruptions since the edifice last grew by three sheets, compound flow fields from fissures round the caldera, on the flanks and on the rift zone, whose eruptions run the farthest and spread as aprons, the flows narrowing on the steep flanks; the summit sinks in collapses that wander, shrink and deepen; the edifice stretches along its rift and its foot is lobed. Partly by F115, the forty-eighth: its lava buries the ground it grows over, craters and all, as a deposit buries what it falls on. Still drawn: the edifice's profile, one power law at five degrees on every body, where its flows' lengths and the lithosphere's support would set it (Venus's large volcanoes are broad and gentle); the foot's lobing; the vents' shares and the run-outs, where the supply and the effusion rates would set them; the edifice's outline its own rule, where its lava would pond in the depressions at its foot (audit M19) |
| Canyons | `trough_depth_at`, `stress_trace`, `trough_throw` (stages/plume.rs) | partly by F105, the thirty-eighth milestone: the faults start where the rise's stress is tensile, strike across its most tensile direction and die out where the tension fades, so the flanks carry radial graben, the crest a maze and the canyons run out from it; they form as the rise's load grows. Partly by F111, the forty-fourth: a fault grows in straight segments that link, turning toward the stress field and stepping en-echelon where they meet; a graben is as wide as its faulted layer sets it, its throw peaking in each segment and dipping where they link; a fault whose throw would exceed the layer is a rift, its segments half-graben whose master side flips. Partly by F112, the forty-fifth: a graben keeps its throw along its length, a confined fault's flat-topped profile; its walls are rough fault traces, so it pinches and swells, and every fault's trace wanders, self-affine; parallel graben stand apart by a few thicknesses of their faulted layer. Partly by F113, the forty-sixth: dike swarms radiate from the rise and its volcanoes, a graben above each dike. Partly by F114, the forty-seventh: a rise grows along an axis, its volcanoes on it and its stress a pure shear about it added to the load's, so the faults beyond its foot gather along the axis and spare the ground across it. Partly by F116, the forty-eighth: the faults form as the dome rises, the dike swarms with the magma that feeds them. Still drawn: the profile across a trough (a flat floor between tanh walls); the faulted layer drawn per fault rather than as a field over the rise; the rise's axis drawn (audit M18); a rift's depth and width from its length, not the strain where it runs, so a long fault far out along a rise's axis can open a canyon alone on the plains; and the stress field's shape, a dome's bending and a load's hoop tension rather than a thin shell solved for the rise's load and support |
| Io's mountains | `mountain_height`, `raise_mountains` (stages/io.rs) | partly by F87, the twenty-eighth milestone: each a block thrust up along a fault striking with the tide's compression, its uplift the throw tapering to the fault's tips times its fall across the block, buried by the plains so that its back and tips are where it dips under them, its trace bowed by the heave and its scarp cut back by landslides. Its sizes, its form's shares and its relief and ridges are still drawn, where the compression of the burial against the lithosphere's strength would set them (Schenk and Bulmer 1998; McKinnon et al. 2001) |
| Io's plume deposits | `deposit_profile`, `plume_history` (stages/io.rs) | partly by F84, the twenty-sixth milestone: each eruption's range from its speed (v²/g), its oval from its vent's elongation, its range scattered with azimuth, its gas condensing in a ring and its pyroclasts falling within; the ring's radial profile is still a fixed Gaussian each side and the fall within a fixed ramp, where the canopy shock and the particles' sizes would set them (Zhang et al. 2003's DSMC plumes) |
| Paterae | `patera_radius`, `foot_edges`, `populate` (stages/io.rs) | partly by F88 and F90: the outline is a chain of nested collapses along the crust's fabric, elongated as Radebaugh et al. measured, cut straight by its faults, and ending at the foot of a mountain it abuts. Partly by F118, the forty-ninth: it collapses piecemeal, more collapses of comparable size breaking out of the first, its margin wanders at every scale, and its floor lies in patches of different ages. The sizes, shapes and shares are Radebaugh et al.'s counts, not yet read from the magma reservoirs beneath; the faults that cut a patera show nowhere outside it |
| Crater walls and rims | `fresh_height` (crater/morphology.rs), `RimShape` (crater/stamp.rs) | closed by F67, the sixteenth milestone: the wall is made by its collapse, terraces of slump blocks on arcuate faults whose widths follow Pearce and Melosh (1986) and the strength length, a scalloped rim and the floor's edge at the blocks' toes |

Impact craters are events and stay objects. Their sizes follow scaling laws, and since F67 their walls follow their collapse (the last row).

When: Io's features with Io; the rest in the earlier-phases revisit.

### T2. Terraces as rings (closed by F92)

Each terrace level of F67 is a ring of slump blocks aligned round the crater, so a deep crater shows four or five evenly spaced scarps on its wall, like the rungs of a ladder. Real terraced walls are fewer, wider slump blocks that merge and split, with scarps of varying height, and near the simple-to-complex transition they are a few slumps rather than terraces. The user saw it on Enceladus (2026-10-15), whose largest craters are laddered; but every crater it shows there is larger than any on the real Enceladus, whose plains' craters are under about 40 km, relaxed by its past heat and muted by fallout (Bland et al. 2012; Kirchoff and Schenk 2009), which M3 and the crater population's high end own.

Direction: terrace levels whose blocks are not aligned round the crater, merging and splitting, with the count rising with the diameter over the transition.

When: the earlier-phases revisit, with the crater population.

Closed in the thirtieth milestone (F92). How much of the wall slumps follows how far the crater collapsed, the fall of its depth-to-diameter ratio from the transition's under Pike's complex law, so a crater just past the transition slumps a few blocks of its outermost terrace and the inner terraces come in as the collapse deepens: a lunar crater of 20 km stands 0.7 terraces at an azimuth, one of 93 km 3.4, where every size had about four. The terraces' spacing varies from crater to crater, and blocks stand off far enough to merge with the terrace inside or outside them. Still drawn: the collapse at which the wall is terraced in full (TERRACE_FULL_COLLAPSE, stated), and each level is still a ring of blocks, if a broken one.


### T3. Bands as lenses, apart from the lanes

A band (F83) opens as a symmetric lens along its crack, widest in the middle. Many of Europa's bands open as wedges, the plates on either side rotating about a pole as they part, as in Argadnel Regio. Bands are drawn by their own stamp, while Ganymede's grooved terrain uses the lane mechanism (F61) of a rifting lid: the two are one process, a lid pulling apart, in lids of different thickness. Lanes also still trace the tide alone, without the freezing and heterogeneity of the shell (F82).

Direction: a band's opening from the rotation of the plates it parts, about the pole their motion defines; bands and lanes as one rift process whose width and fabric follow the lid; the lanes' strikes in the shell's whole stress.

When: the earlier-phases revisit, with the lanes.

### T4. Rays as spokes (closed by F91)

A fresh crater's rays are drawn as straight, blurred spokes about a glow, a lens flare, on every airless body; the survey of 2026-10-16 saw them on the Moon and Mercury, and the fourteenth milestone's note saw Callisto's as a starfield. Real rays are discontinuous filaments of secondary craters and the bright ejecta they throw, in chains and herringbone patterns that narrow, branch and loop outward from a continuous ejecta near the rim (Shoemaker 1962; Melosh 1989).

Direction: rays from the ejecta curtain: clumps of ejecta thrown along paths whose ranges and azimuths scatter, each landing as a chain of secondaries with its own bright streaks, the whole fading with maturity.

When: the earlier-phases revisit.

Closed in the thirtieth milestone (F91). A crater's ejecta curtain gathers its mass into some thirteen jets; beyond the blanket a ray is the share of the ground the jets' secondaries excavate from under the weathered skin, falling with distance by Elliott et al. (2018)'s secondary-crater model and scaled to their lunar ray lengths, so it narrows to its jet's midline, breaks into segments and fades. Noise bends the jets' paths and bunches their fragments. The blanket's lobes and the loop craters' secondaries follow the same jets. An icy surface's lag sets the skin the secondaries dig through where it is thicker than the Moon's, so Callisto's rays are short or absent under its thick low-latitude lag. Still open: rays bright by composition, where a crater throws one unit's material over another (Lichtenberg's highland rays over mare, Hawke et al. 2004), which needs the primary ejecta's share along the ray; the secondaries' herringbone chevrons, below the grid at the map scales; the skin depth of rock bodies other than the Moon and every body's fragment population, which set their rays' reach (the lunar ray lengths are the calibration, below); the rays' fading, still the provisional ray lifetime rather than the maturing of the ground the secondaries dug, which waits on Mercury's weathering timescales; and the deep material a crater digs up, which lies in a soft ramp to the blanket's edge (Mercury's dark halos).

### T5. The dichotomy as a hyperbolic tangent

The primordial crust's hemispheric dichotomy is half its amplitude times tanh(3 cos θ) about an axis: a smooth slope over thousands of kilometres. On Mars (20 km of crust, some 2.4 km of relief) it does not show under the basins, with no boundary scarp, no fretted terrain and no lowland plains younger than the highlands.

It is also a rule where it should be an outcome. `tectonics_dichotomy_amplitude` gives every stagnant-lid rock body a contrast of 20 km in crust thickness, an episodic lid 5 km, an icy body 2 km and a mobile lid none, by the regime's name (S1), so Mercury, the Moon and Mars carry the same one. Mars's own parameters do not guarantee it (the user asked, 2026-10-16). Both leading explanations are contingent: a giant impact (the Borealis basin, 10,600 by 8,500 km, Andrews-Hanna et al. 2008; an impactor of 1,600 to 2,700 km at 6 to 10 km/s and 30 to 60 degrees, Marinova et al. 2008), which late accretion makes a matter of chance, since it was dominated by a few massive bodies and the median largest impactor on Mars in Bottke et al. (2010)'s Monte Carlo was 1,500 to 1,800 km; or degree-one convection, which needs a weak asthenosphere or a viscosity jump in the mid-mantle (Zhong and Zuber 2001; Roberts and Zhong 2006).

Direction: the dichotomy as an outcome. The giant-impact tail of the late-accretion population in the impact stage, a mega-impact taking its own form (the crust removed inside an oblique ellipse, its ejecta thickening the rest, its melt ponding in the basin), so that some bodies have a hemispheric basin, some none, some an off-centre or elliptical one and some a giant basin short of a hemisphere, as the Moon's South Pole-Aitken is; the boundary scarp worn by the processes that wear a scarp, and the lowlands resurfaced by the lavas and sediments that bury their craters to ghosts; an endogenic route once the interior model resolves the planform of the mantle's convection; and for a locked moon of a hot planet, the planet's heating of its near side while its crust forms (Roy et al. 2014). The amplitude by name goes. A preset for a real body may carry a seed whose draw includes its event: one draw of the model, not a pin.

When: the earlier-phases revisit, with Mars.

### T6. Grooves as regular troughs

The user's note (2026-10-19): Ganymede's and Enceladus's grooves are troughs too regular to read as faulted ground, and should be more broken up, undulating and organic. The same process should show on Miranda, as a check on the model. Miranda's coronae are belts of subparallel ridges and troughs that Pappalardo et al. (1997) read as tilt blocks from extension and normal faulting above an upwelling (Arden Corona, local extension near 70 percent; warm buoyant ice in Hammond and Barr's convection): the faulting that makes Ganymede's grooves, concentric about the upwelling where Ganymede's follow its lanes. So one process, and Miranda the second body GROOVE_RIFT_LID_SHARE needs.

Direction: the troughs from the faulting itself, tilt blocks and normal faults whose throw and spacing vary along their strike, segments that end, step over and link, so that a groove undulates and breaks; a Miranda preset, its coronae from its upwellings, its grooves' spacing checked against the measured.

When: with the icy moons' tectonics.

Closed in part by F132, the sixty-third milestone: a groove is a chain of fault segments that come and go, step across where one gives way to the next and taper to their tips, dying out before the lane's bounding faults. Still open: Miranda, whose coronae come from upwellings (the user asked whether they should be driven by Ganymede's lanes: they should not; the outline and stress field are a corona's, as on Venus, the fabric this faulting), and the grooves' relief on Enceladus, tens of metres where its sulci's ridges stand up to 1 to 1.5 km, a check on GROOVE_RIFT_LID_SHARE and GROOVE_RELIEF_SHARE on a second body.

Miranda closed in part by F133, the sixty-fourth milestone: a Miranda preset, and its coronae over the upwellings of the convection its resonance's tidal heating drove in its ice, their belts of the same fault segments concentric with outlines the convection's cells give. The spacing passes its check on a second body: Miranda's coronae's faults come out 4.9 to 5.3 km apart on the lid of its episode, against Arden's 5 km and the 5 to 10 km Hammond and Barr (2014) tabulate for its coronae. The relief does not, under the lanes' law: a share of the spacing gives the belts a few hundred metres where Arden's ridges and troughs stand up to 2 km, the 70 percent extension Pappalardo et al. (1997) find in its tilt blocks. A corona's blocks tilt instead until their faults lock, from the 50 degrees Arden's broke at to the 30 at which friction locks a normal fault, which gives Miranda's belts 1.7 km. The lanes keep their share, calibrated on Ganymede. Still open: one law for both, the relief from the strain the faults take up, a corona's to lock-up and a lane's what its rift opens, with Enceladus's grooves' relief.

### T7. Tiger stripes as a set per region

The user's note (2026-10-19): the height map draws Enceladus's polar terrain as prominent, fine fingerprint patterns, and each occurrence of the terrain carries its own set of four to nine stripes, where the real terrain is a network of interconnected, webbing canyons across all of it. Losing the signature four or five parallel fractures at the south pole is a trade the user accepts for natural terrain. The stripes are also coloured light blue, a false colour.

Direction: the stripes as one network of fractures across the terrain the plume's heat renews, linking, branching and ending where the stress does; their relief at the canyons' scale rather than a fingerprint's; their colour their ice's.

When: with T6 and the user's further notes.

## Calibrations on one body (principle 15)

Each constant here was fitted so that one body comes out right. Each should become, or be shown to be, a constant of its process, with a check on what it gives a body between the references.

| Constant | Fitted to | Note |
|---|---|---|
| ICE_RESURFACING_FLUX, ICE_RESURFACING_FLUX_FLOOR | Ganymede: two thirds renewed, the shell left near 40 km | a constant of renewal per excess flux, or the strain an episode's extension opens |
| FOCUSING_PERIOD_FACTOR | Enceladus: episodes 0.52 Gyr apart (refitted from 10 to 40 when F68 thinned its shell to 22 km) | O'Neill and Nimmo 2010 give 0.1 to 1 Gyr from the shell's convection; the direction is a recurrence from the heating cycle of the Dione resonance, not a factor on the conductive time |
| POLAR_FOCUS_FLUX, POLAR_THINNING_MAX | Enceladus | how strongly tidal heating focuses at a pole |
| PLUME_FALLOUT_RATE | Enceladus: about 0.1 mm/yr near the stripes | scale by the plume's mass flux, which the vent sets |
| CHAOS_FRACTION_MAX, CHAOS_FLUX | Europa's share of chaos (CHAOS_FLUX refitted from 0.05 to 0.012 W/m² when F75 brought Europa's tidal flux from 76 to 14 mW/m²) | replaced under S4 |
| LINEA_WIDTH_EXPONENT, LINEA_MARGIN_ONSET_GYR, LINEA_STAIN, LINEA_STRIPE_DIM, LINEA_STAIN_RESIDUAL | stated choices for Europa's lineae (F78): the heavy tail of widths, the time margins take to show, the stain's strength, the medial stripe's dimming and the ancient band's residue | the time scales of the margins' deposit and fading (Fagents et al. 2000's cryoclastic margins; sputtering and frost) |
| CRACK_DENSITY_FACTOR | Europa: one crack per 31 km of transect (refitted from 3 to 12 when F75 brought its shell from 7.7 to 30 km) | crack spacing from the brittle layer's mechanics, the spacing at which a crack's stress shadow lets the next open |
| IO_PATERAE_PER_TW | Io's count over its tidal power | per heat-pipe power once S3 is done |
| ICE_LAG_SHED_START_DEG, ICE_LAG_SHED_FULL_DEG | Callisto: crater walls of 10 to 20 degrees bright over dark plains (Galileo PIA01297, PIA01649; Moore et al. 1999) | the stability of a sublimation lag on a slope, from its grains' cohesion and the rate the ice under it sublimes |
| LAVA_BATCH_SHARE | the Moon: neighbouring mare units hard to tell apart by albedo alone (Iqbal et al. 2019) inside the factor of 1.4 its major units span (Denevi et al. 2009) | each eruption's departure from its source, from the magma's own fractionation |
| IO_PELE_CLASS_FRACTION | about six large plumes at a time on Io (Geissler et al. 2004) | the share of centres whose magma comes up gas-rich, from the depth and temperature of their sources |
| IO_SHIELD_COLUMNS | Io's red-brown plains beginning near 30 degrees (Geissler et al. 1999; Williams et al. 2011) | the column that stops the plasma, from its ions' energies and the gas's stopping power |
| GROOVE_RIFT_LID_SHARE | Ganymede's 7 km grooves | a process constant, calibrated once; checked on a second body (F133): Miranda's coronae's faults 4.9 to 5.3 km apart on the lid of its episode against Arden's 5 km and 5 to 10 km for its coronae (Hammond and Barr 2014) |
| GROOVE_RELIEF_SHARE | Ganymede's grooves, 250 m at 7 km (Squyres 1981) | fails its check on Miranda (F133): a few hundred metres where Arden's tilt blocks stand up to 2 km at 70 percent extension; the coronae take their relief from their blocks' tilt at lock-up instead, and a lane's wants the strain its rift opens |
| CENTRAL_PEAK_HEIGHT | Copernicus | check how it scales with gravity |
| RAY_LENGTH_COEFF_KM, RAY_LENGTH_EXPONENT, RAY_SKIN_DEPTH_M | Elliott et al. (2018)'s median lunar ray lengths, 6.59 R^1.27 km, for their 10 cm skin (F91) | the skin and fragment population per body: the length goes as the skin depth to about -1.16, which an icy surface's lag already sets |
| CURTAIN_SPREAD, CURTAIN_SLOPE, CURTAIN_ORDERS | the Moon's rays: 10.7 a crater, 22 degrees wide at two radii (Elliott et al. 2018) | the target's relief at the rim, which makes the rays (Sabuwala et al. 2018: the count goes as the crater over the relief's wavelength) |
| TERRACE_FULL_COLLAPSE | a stated choice: a lunar crater's wall terraced in full from about 60 km (F92) | check the share of the wall terraced along a sweep of lunar and Mercurian craters |
| POROUS_DISSIPATION | Choblet et al. 2017's model of Enceladus's core, 10 to 30 GW: 28.5 GW, which with its shell's own 0.6 GW (F75) holds its shell at 22.7 km (refitted from 15 GW when the shell's share went from the class k₂/Q's 15 GW to F75's) | a constant of the dissipation per unit of porosity; fails its check on Mimas (F135): its eccentricity of 0.0196 gives it 228 GW, a shell 1.7 km thick over an ocean and a surface renewed every 0.1 Gyr, where Lainey et al. 2024 find a shell 20 to 30 km thick over a young ocean under an old surface |
| RING_STRAIN | the strain of Tyre's and Callanish's inner graben (Singer et al. 2023) | check against the graben depths of Valhalla and Ganymede's furrows once measured |
| PIT_DIAMETER_SHARE, PIT_GROWTH, PIT_DEPTH_SHARE, DOME_PIT_SHARE, DOME_RISE | Ganymede's pits and domes (Alzate and Barlow 2011; Bray et al. 2012; Schenk 1993) | check against Callisto's dome craters and Ceres's pits (Occator, Dantu) |
| MANTLE_SPREAD | a stated choice: the length over which a mantle evens the ground, twice its thickness | check against the rim profiles of Enceladus's muted craters (Kirchoff and Schenk 2009) |
| FLOW_NARROWING_SLOPE | Tharsis: a 5 degree flank narrows a flow to about a third of its width on the plains (Peters et al. 2021: the narrowest flows on Olympus Mons's flanks) | a flow's width from its lava's yield strength and its effusion rate on the slope (Hulme 1974) |
| FLOW_STEEP_SLOPE | a stated choice: a degree, above the plains' slopes, so that the plains' floods and Io's flows keep the steering they had | the slope at which gravity's push across a front outweighs the bumps it rides over; check against flow paths on a shield's flanks |
| ICY_FLUX, ICY_FLUX_RISE | Nesvorný et al. 2023's comet rates on the moons of Jupiter, Saturn and Uranus, each over its planet's focusing (F135) | fitted on three planets; within a factor of 1.3 on six moons; check against Triton's and Pluto's craters, beyond the planets it was fitted on |
| ICY_EARLY_TOTAL | Wong et al. 2023's comets of 1 km over the history of Mimas and Enceladus, over Nesvorný et al.'s rates (F135) | passes on Mimas, its counts at 10 and 20 km within 7 percent; check Iapetus's oldest plains and Callisto's |
| ROCKY_BELT_FLUX, ROCKY_BELT_RAMP, ROCKY_OUTER_FADE | Mars about twice the Moon (Ivanov 2001); the giant planets' moons struck by comets, not asteroids (Zahnle et al. 2003) | checked on Ceres (F135): 13 craters over 100 km against the 16 counted; Mars, its basins from the Mars production function continued past its range, takes two impacts wider than itself after its crust solidifies and starts its record at 3.95 Ga |
| ICE_HOLD_RADIUS_FRACTION, ROCK_HOLD_RADIUS_FRACTION | the largest craters icy moons keep (Odysseus 0.84 R) and Rheasilvia on Vesta (1.93 R) | the hold from the energy that disrupts a body against its gravity and strength (Movshovitz et al. 2015) |
| EARLY_TIME_EXPONENTS | stated choices: the orbital period at the snow line, the disc's solids over the area within it (F135) | no body to check; around a red dwarf the snow line where its planets formed is set by its far brighter youth, not by its luminosity today |

Constants marked provisional for Phases 5 and 7 (atmosphere, albedo, imagery) are left to those phases.

## Checks against measurement

### M1. Enceladus's ice shell (closed by F68)

The model derives a mean shell of 44 km for Enceladus: a conductive shell under 18.9 mW/m², from 14.9 GW of tidal heating, inside validate's range of 5 to 40 GW. Gravity, shape and libration give a mean of about 18 to 25 km, about 5 km at the south pole (Čadek et al. 2016; Beuthe et al. 2016; Hemingway and Mittal 2019), and conducting through a shell that thin loses some 25 to 40 GW.

The consequence shows on the maps. Under a 44 km shell no crater relaxes, so a 140 km crater stays crisp, where the real cratered plains have nothing crisp above about 40 km and their craters are relaxed in step with the heat flow (Bland et al. 2012; the reference map in research/reference/maps).

Direction: the heat budget, never a pinned shell. Tidal heating in a porous, unconsolidated core can supply 10 to 30 GW and keep the ocean warm from below (Choblet et al. 2017); with it the shell follows from the heat it must carry. Add the shell thickness to validate so that the result is checked.

When: with the ring systems of large impacts, whose palimpsests and relaxed rings depend on the same heat flow history.

Closed by F68, the seventeenth milestone: tides dissipated in a porous rocky interior, its porosity from the pressure in it, bring Enceladus's tides to 30 GW and its shell to 22 km, and validate checks the shell against 15 to 30 km. Its craters relax more beside its terrains and across its plains; the rims of the largest stand on its 12 km lid, which holds wavelengths of a few kilometres up.

### M2. Palimpsests are dark

A palimpsest is a basin whose relief has relaxed away and whose centre stays bright: Valhalla's is 650 km across, Ganymede's Memphis Facula is a bright disk (Jones et al. 2003). The model relaxes the old basins flat (F39, F70) and breaks their lithosphere into rings (F69), but their melt sheets and deep ejecta darken back to the lag with the rest of the surface, so a palimpsest four billion years old is as dark as the ground about it. The physics: a basin's melt sheet is clean ice laid from below, hundreds of metres to kilometres thick, far thicker than the layer gardening mixes, so it stays clean as the lanes do and darkens only by its own impurities' lag; an ordinary crater's melt and ejecta are thin and gardened into the regolith.

When: Phase 5, with the lag of icy surfaces and Callisto's bright craters.

### M3. Enceladus's relaxed craters

Enceladus's cratered plains hold craters relaxed by 40 to 90 percent, which needs a heat flux of 150 mW/m² or more where and when they relaxed (Bland et al. 2012), and muted by its plumes' fallout (Kirchoff and Schenk 2009). The model's heat flux is near 38 mW/m² over its whole history, and its elastic lid of 12 km holds every load shorter than about 45 km, so its craters relax only beside its terrains while each episode's thinned shell lasts (F66). Its fresh complex craters are deep, as low gravity makes them (F71, Herschel on Mimas), and the largest of them keep their terraced walls, softened only by the fallout (F74), which is thin away from the plumes' poles.

Direction: the heat history that relaxed the real craters. An eccentricity cycle in the resonance with Dione is not it: Enceladus does not oscillate about its tidal equilibrium in the model of Ojakangas and Stevenson (Meyer and Wisdom 2008), and its shell's own tides are a small share of its heat (F75). Episodes of focused heating whose thinning reaches the cratered plains between the terrains, or a past resonance that raised its eccentricity for a time (Meyer and Wisdom 2008 find a temporary capture into the 4:2 resonance with Dione), remain.

When: with the earlier-phases revisit's thermal histories, or sooner if the Europa milestone revisits the orbital histories.

In the sixty-sixth milestone (F135) the outer system's own comets strike the icy moons, and Mimas, the least altered of the small moons, keeps its counts at 10 and 20 km within 7 percent; Enceladus's oldest plains, from its record's start at 4.25 Ga, keep 1.5 to 3.7 times its counts at 5 to 20 km (Kirchoff and Schenk's counts as Wong et al. 2023 give them), and their craters stand sharp. Its plains are younger than its record, 4.10 Ga by Wong et al.'s isochrons, and relaxed: the heat history this entry wants would give both. Miranda's cratered terrain, mantled and subdued on Voyager's images, likewise keeps more and sharper craters than its map.

### M4. Enceladus's old fractures

Enceladus's cratered plains are cut by fractures. The class rule's Love number gave its shell an h₂ ten times Beuthe's (2018) and a diurnal stress of 194 kPa, which cracked the whole shell; with the shell's own Love number (F75) its uniform shell is stressed by 20 kPa (Beuthe: about 15 kPa in a uniform 23 km shell), under the 40 kPa at which ice cracks, and only its thinned south polar shell and the zones its past episodes thinned are stressed past it (F76). Its old plains therefore lose the cracks the diurnal tide had drawn there.

Direction: the stresses that crack a shell beyond the diurnal tide: the tension of a thickening shell as its ocean freezes and pressurizes (Manga and Wang 2007; Rudolph et al. 2022), which each episode's end and each fall of the heat flux bring; true polar wander as the body reorients to each new episode's pole (Matsuyama and Nimmo 2008); and non-synchronous rotation. Kihoulou et al. (2025) find extension in the upper shell as Europa's shell thickens after a pulse of eccentricity, and compression as it thins.

Since F82 (2026-10-15, the twenty-fourth milestone) a shell that thickens stands at the cracking stress in tension. Enceladus's mean shell, held by its porous core's steady heat (F68), thickens by only 0.6 km over its life, some 650 Pa of tension, so its plains still do not crack; the thickening that would crack them is local, the zone refreezing at each focusing episode's end, which F82's global history of the shell does not see.

When: with a local history of the shell under the focusing episodes.

### M5. Europa's old bands are bright

Galileo's young bands are dark grey wedges and its old ones brighter than the plains around them; the reference shows bright bands in the south. A band (F83) is grey when young and fades with the salts to the ancient band's faint stain, never brighter than the plains.

Direction: the brightening of old band ice from what it is made of and how it weathers, against the plains' radiolytic darkening.

When: Phase 5's materials work.

### M6. Europa's cracks follow the tide closely

Europa's cracks, under the tide, the shell's freezing and its heterogeneity (F82), still trace the tidal field's great arcs more plainly than the real moon's, whose lineae cross more and undulate more at every scale. The user found the look acceptable (2026-10-15).

Direction: the stresses that free cracks from the tide's pattern beyond those modelled: the interaction of cracks with the open cracks about them for the narrow fabric (F81), whose tracing answers no other crack; non-synchronous rotation's rate against the margins' life; true polar wander.

When: the earlier-phases revisit, or sooner if the user asks.

### M7. Io's white plains

Io's white and grey materials, moderate-to-coarse-grained SO2, cover about 27 percent of the surface, mostly at low latitudes (Geissler et al. 1999), and its white diffuse deposits often surround the margins of lava flows (Williams et al. 2011). F84's frost comes only from plumes' rings and young flows' halos and sublimes at Ingersoll's rate wherever it lies, so frost shows over about 5 percent of the surface; its white is also brighter and more opaque than the reference's.

Direction: the frost cycle of the sublimation atmosphere, which takes SO2 from warm ground and lays it on cold, and the albedo feedback that keeps a bright frost cold, so that frost seeded by deposits persists in fields; coarse-grained, translucent frost apart from a fresh deposit's fine white.

When: with Phase 5's materials, or sooner with the flows.

Since F89 the plains' fallout varies with the volcanic centres and its clumps, which mottles the plains where the plasma alters them, beyond about 30 degrees. At low latitudes the plains are one sulfur under a fall of sulfur, unaltered, and stay nearly smooth: their white, yellow and grey mottling waits on this frost. The plume deposits' own mottling (F84, IO_DEPOSIT_MOTTLING) reaches its stated 0.35 only at the noise's extremes and is about 5 percent typically; F89's clumping of the background is drawn at a unit spread instead.

Tried after F89 and not kept (2026-10-16; the experiment is not on main): the albedo feedback, each frost column subliming at the vapour pressure of its own temperature, thin frost the reference frost over the sulfur plains and darker over darker ground, thick frost fresh frost's; and a frost cycle, what sublimes condensing again over the region about it. SO2's vapour pressure falls some 200-fold for 15 K, so the feedback decides everything. With fresh frost bright enough to be the coldest ground, as McEwen et al. (1988) find Io's SO2-rich ground (temperatures varying by 40 K with albedo under the same sun), every plume's ring became a cold trap and stayed white through the history, half the surface white. Calibrating fresh frost's bolometric albedo to 0.78 brought the frost to Io's 27 percent, but the condensation, spread evenly over 200 to 500 km, built its fields at 30 to 50 degrees of latitude, where bare sulfur is colder than dim frost at the equator and where Io is poor in SO2 (McEwen et al. 1988), and left the equator with rings. Sintering, fresh frost coarsening to the reference in ten years, cancelled the feedback. What is missing: SO2 laid broader than a plume's ring (its gas condensing across the deposit, the atmosphere collapsing each night), the frost's temperature through the day with its thermal inertia, and condensation that goes to the coldest ground the atmosphere reaches, so that bright frost at the equator traps it before the mid-latitudes do. These belong with Phase 5's volatile transport.

### M8. Io's flow coverage

Io's flow fields cover 28 percent of its surface, its plains 66 (Williams et al. 2011). F85's fields cover 20 percent and the plains 70: the fields' areas are stated shares of their reach, and the number of fields three generations a patera. Since F87 and F88, with the lava that ponds in a patera counted as its floor, the fields cover 16 percent and the plains 79. Since F117, the flows floods from the paterae's rims that run down and pond in the closed lows they reach, they cover 19.9 percent and the plains 74.5, the patera floors 2.3.

Direction: the fields' areas from the lava erupted over their lives, the eruptive share of the magma rate (Spencer et al. 2020: most of Io's magma intrudes) over a field's thickness, so that the coverage follows from the heat the lava carries out.

When: with S3, the regime's switch.

### M9. A rise stays where its plume rose

Mars's rise came out at 83 degrees south, its shields round the pole. A load the size of Tharsis turns the body until it lies on the equator, which is why Tharsis does (true polar wander: Bouley et al. 2016). ExoGen turns a body for a light terrain, Enceladus's thinned pole going to the pole (F53), but not for a heavy one.

Direction: the body's reorientation from its loads, positive and negative, the rise carried to the equator as it grows, and everything older turned with it.

When: the earlier-phases revisit, with Mars's rise.

Closed in part by F104, the thirty-seventh milestone: the body turns until its axis of greatest moment of inertia is its spin axis, from its rises' mass against the rotational bulge its elastic lithosphere holds frozen, so Mars's rise is carried from 83 degrees south to the equator and Venus's three rises to within 12 degrees of it. Still open: the negative loads stay with F53's own turn (Enceladus's thinned pole), and the rises' shields and floods add no mass of their own; the latitudes things had before the turn will matter when the climate's processes come in Phase 5.

### M10. Venus's deformation belts and rifts

Venus's plains are cut by belts of ridges and fractures and by rift zones (chasmata) that join its highlands and coronae into a network, the bright web of the radar maps; coronae string along the rifts. ExoGen's Venus has none: 506 coronae stand alone as thin rings, the tesserae are plateaus of noise, and the plains between are smooth under their wrinkle ridges.

Direction: the lid's stress from the mantle's convection, rifting where it pulls apart and folding where it pushes, coronae where plume heads meet it, as annuli of fractures and ridges round a dome or sag, and tesserae as the oldest crust deformed in the belts.

When: the earlier-phases revisit, with Venus.

### M11. The rocky bodies' colours

Mars's albedo is one grey, Venus's another, and Mercury's flooded floors are flat tan patches in dark halos, where its fresh craters' ejecta should be bright. Mars's red is the oxidation of its dust and Mercury's tones its materials' and their maturing; the user left Mercury's colour for a later pass (Phase 3).

Direction: the weathering and transport processes of Phase 5 (oxidation, dust, maturity), with the composition of each unit.

When: Phase 5, or a colour pass the user asks for.

### M12. Rayed craters on the icy moons

With F91's rays, Ganymede has 111 craters young enough to keep rays (under the 1.5 Gyr ray lifetime), twelve of them over 86 km, and its rays brighten 28 percent of its surface by more than a tenth of the fresh contrast, densest at middle and high latitudes where its lag is thin. Callisto's thick low-latitude lag hides its rays there, but poleward of about 55 degrees its lag thinned and rays from craters far away appeared across the cold ground, switching on across a line. Since F93 the lag is summed over unresolved facets whose spread follows the surface's roughness, and the rays come in over 30 to 50 degrees, patchy at first. Since F96 a young crater's discontinuous ejecta run out along its jets, so its spikes go on into its rays with no grey disc between; the high latitudes' rays still lie long and bright over the dark lag, where the real ones would lie on frost (the direction's last point).

Direction: check the counts of rayed craters on Ganymede and Callisto, and their rays' reach, against the Voyager and Galileo maps; the ray lifetime from the maturing of the dug-up ice (gardening, lag, sputtering) rather than RAY_FADING_MAX_GYR; and the skin at high latitudes, where frost and the gardened regolith, not a lag, mature the surface.

When: with the icy moons' crater population, or a colour pass.

### M13. Callisto's troughs and rims stand undegraded

Callisto's ring troughs and crater rims keep the shapes they were drawn with, smoothed by diffusion, where Galileo shows them broken into knobs and short scarps with the dark lag gathered in the lows (PIA01297, PIA01649; Moore et al. 1999). The process is sublimation: ice sublimes from the walls that have shed their lag (F98), so they retreat, and the crests between retreating walls break into knobs (Howard and Moore 2008). A pass breaking crests by the ice their walls sublimed against their height was tried in the thirty-fourth milestone and set aside: at the lag model's rate from the latitude's mean temperature it barely acts outside Callisto's equatorial band, while sublimation follows the daytime peak, some twenty to thirty times faster; at that rate it would strip Ceres, which the model treats as icy and warm, because the model does not place a body's near-surface ice at depth.

Direction: an ice table, ice stable at the surface only where the sublimation it would lose over the surface's life is small and elsewhere below a dry layer as deep as the sublimation and the vapour's diffusion through the layer leave it; sublimation from the diurnal temperatures, the lag (F55, F93) and the shedding (F98) from the same; then the walls' retreat and the crests' breaking into knobs from the ice the walls lose, on every icy body by its own temperature and ice table.

When: Phase 5, with the surface temperatures through the day and the year and the stability of the volatiles.

### M14. A rise's late floods fill its canyons (closed by F108; reopened in part by F116; closed by F119)

On Mars the floods keep coming over the rise to the present, at the pace of the body's magma production, so a canyon opened in the first half of the body's life is mostly filled by the lava that follows it: the deepest, 7.4 km, keeps 2.8 km of relief and a floor of flows, and its landslides (F106) lie under the lava. Valles Marineris is younger than the plains about it and stands open, because Mars's later volcanism gathered at its great shields as its lithosphere thickened, rather than flooding the rise (Hauber et al. 2011 on Tharsis's evolution from plains to central volcanism).

Direction: the share of a body's magma erupted in floods falling as its lithosphere thickens and the share erupted at central vents, the shields, rising, so that late volcanism builds the shields and leaves the canyons open.

When: the earlier-phases revisit, with the Mars dichotomy's lowland burial, which also sets where the floods go.

Closed in the forty-first milestone (F108): a rise's shields take its eruptions as they start, each by its share of their volume, and what they take leaves the floods. Mars's canyons carved after its first shields started keep their relief, its largest a median 2,980 m of its 7,400 m along its trace, against 790 m before; those carved during the plains volcanism stay buried. Still open: when a shield starts is a time drawn between 1.0 and 2.5 Gyr (SHIELD_START_GYR), not a threshold in the lithosphere's thickness, which is what the direction asks; a shield's size comes from its own rule, not from the magma it takes, so Mars's edifices hold 8.3 million km³ of the 61 million km³ the shields take, and the rest is not laid; and the shields' own aprons in the loop wait on M15.

Reopened in part in the forty-eighth milestone (F116): the rise's faults now form as its dome rises, so its canyons open in the Hesperian, during the plains volcanism, and the rise's eruptions before its shields take them fill them in part, as this entry's direction foresaw: at the preview quality the deepest standing canyon stands 1.7 km below its rims, against 2.7 km when the faults ran on to the present. The direction stands: the share of the magma erupted in floods falling as the lithosphere thickens, the shields starting at that threshold and taking the rest.

Closed in the fiftieth milestone (F119), by another control than the direction's. A threshold in the lithosphere's thickness does not hold across the bodies: the Moon's maria erupted as floods through a lithosphere the model makes 210 to 320 km thick, Io erupts floods through a thick one (O'Reilly and Davies 1981), and a threshold at Mars's turn would never start Venus's shields. A volcano's vents turn central where its load's pressure outweighs the regional differential stress and gathers the dikes (Muller, Ito and Martel 2001): a rise's dikes run along the stretching of its growing dome as fissure floods, and its shields gather them as the growth slows, the share by which the growth's rate has fallen from its peak, each shield starting as that share reaches what the larger ones hold. At the standard quality Mars's canyons keep a median half their carved depth below their rims, against 38 percent, and its shields start in order of size from 0.90 to 1.48 Gyr. Still open: the canyons opened early in the dome's growth fill with the floods of that growth; the shields' edifices keep their own size rule, holding less than the magma they take.

### M15. Floods run off the slopes (closed by F109)

The loop's lava runs from its vent, and from a full lake's pass, down the steepest way to the next low, and lays nothing on the way (F107). Lava flows leave their sheets on the slopes they run down: Tharsis's flanks and the plains about its shields are covered by flows hundreds of kilometres long from their vents. In the model, a rise's flanks keep their older ground, graben and craters, except in the lows on them, and the lava gathers in its canyons and at its foot. Lava whose cascade would run past the flows' reach of 1,500 km is not laid at all: 1 percent of the Moon's erupted volume, 9 percent of Mercury's, 18 percent of Mars's and 24 percent of Venus's; since F108, which gives the shields a rise's late eruptions, 12 percent of what Mars's floods erupt and 9 percent of Venus's. A shield's eruptions could not be laid as floods round it for the same reason: on its slopes there are no closed lows to hold them.

Direction: a flow's sheet laid on its way down, a flow's thickness over the cells it crosses and as wide as its supply spreads it, until its volume is spent or it reaches a low, which is the loop's spreading following the flows (W1's open part). The finish would draw it as a sheet over the ground, not a fill to a level, which a sheet on a slope cannot be.

When: the earlier-phases revisit, with M14, which sets where a rise's late lava erupts.

Closed in the forty-second milestone (F109): the lava leaves its sheet, a flow's thickness, on each cell it runs over, as wide as its volume spreads it (a swath 2.2 times as long as wide, from the Roza flow field), until its volume is spent or it reaches a low; the finish draws a sheet over the ground. Less than 1 percent of any body's lava is lost past the reach, and the Moon's lava covers 16.3 percent of its surface against the maria's 17. Still open: the swath's aspect is one flow field's; a flow's width from its supply and the slope (Hulme 1974), its thickness from the lava's yield strength, wait on the same revisit as W1's.

### M16. The finish's floods keep four fifths of the loop's lava

A flood's pass at the finishing resolution takes its coarse smear off where its flows do not reach, and that lava goes nowhere: on the Moon the passes keep about four fifths of the lava the loop laid. A visible flood's flows cover about a third of its runs, where F94 meant seven tenths when its lakes were deep, and a buried flood loses some at its outline and where the fine ground stands above its lake. Where floods stack, the fine ground drifts below the loop's surface. The lakes refill the drift as finer hollows (F95, F110), so it shows as gaps between flows rather than steps, and the Moon's lava covers 15.5 percent of its surface against the 17 percent the maria cover.

Direction: flows that cover their runs, their lobes sized by the lava's supply and the slope (Hulme 1974, W1's open part), with the lava of the gaps on the lobes; then each pass keeps its flood's volume. Spreading the lost lava evenly over the cells laid was tried in the forty-third milestone and made the plains grainy.

When: the earlier-phases revisit, with W1.

### M17. Tharsis's long radial graben come from dikes the model lacks (closed in part by F113)

The longest graben systems on Mars, Memnonia, Sirenum, Icaria, Thaumasia and Claritas Fossae, radiate 3,000 to 4,000 km from Tharsis (Wilson and Head 2002), and Venus carries 163 radial systems of graben, fissures and fractures, radii averaging some 325 km and some beyond 2,000 km, at least 72 percent of them made by dike swarms (Grosfils and Head 1994). Mars bears up to 16 giant dike swarms, most about Tharsis (Ernst, Grosfils and Mège 2001); dikes lie beneath the graben south of Arsia Mons (Goudy and Schultz 2005), though whether one swarm or several, and how much of the radial faulting is the regional stress instead, is debated. The model's graben are tectonic only: their lengths follow fault populations' power law, averaging 85 km on Mars, and the systems thousands of kilometres long do not form.

Direction: radiating dike swarms from a rise's plume and its volcanoes, each dike propagating laterally along the stress trajectories as far as its magma carries it, a graben above it as wide as twice the depth to its top (Mastin and Pollard 1988) and as long as the dike, with pit chains where the ground collapses into it (Mège and Masson 1996); the swarms' magma from the plume's intrusive budget.

When: the next milestone.

Closed in part in the forty-sixth milestone (F113): a rise's plume head and its volcanoes send out radiating dike swarms, a graben above each dike as long as the dike, so Mars carries 68 graben longer than 1,000 km and Venus fans of long graben about its volcanoes. Still open: the swarms' counts, reaches and fans are stated choices rather than the plume's intrusive magma budget and the dikes' overpressure; the pit chains that form above dikes; and dike swarms about the volcanoes and coronae off the rises, which on Venus are most of its radial systems.

### M18. A rise's axis is drawn, not computed

F114 gives each rise an axis, its volcanoes lining up on it and its stress a pure shear about it, tension across the axis, fading beyond the foot as the load's own stress does, so its faults gather along the axis and spare the ground across it, as Tharsis's dike swarms run out to the south-west and north-east along the Tharsis Montes' trend (Orlov et al. 2023). The axis's azimuth is drawn per rise, as the first dikes across its crest break its symmetry, and the pure shear's share is a stated choice. A regional stress that does not fade with the rise's own turns dikes into its strike far from their source (Muller and Pollard 1977); the model's, tied to the rise, cannot, and one fading as the inverse of the distance, tried, kept the ground along the axis faulting out to 2.3 rise radii, with canyons alone on the plains. The regional stress has sources the model could compute: the body's reorientation under the rise (F104), which stresses the whole lithosphere in a pattern set by the pole's path (Melosh 1980), late on Mars driven by the Tharsis Montes' own lavas (Kite et al. 2009); the crustal dichotomy, its thick southern crust against the thin northern; and the other rises' loads. The centres of Tharsis's faulting also moved with time, from Claritas in the Noachian to Pavonis Mons in the Amazonian (Anderson et al. 2001, from 24,452 structures in five stages), where the model's rise has one centre for its whole life. And the sectors the graben survive in are partly burial: the young lavas west of Tharsis cover its older graben there; the model's floods and shields bury older graben, but which sectors they spare follows their flows' paths, not the axis.

Direction: the regional stress from those sources, and the axis from them where they set it, drawn only where nothing does; the dikes turning into the regional strike far from their source; the rise's centre of activity moving over its life.

When: with the rise's flexure and the dichotomy, the revisit plan's fifth item.

### M19. A volcano's load in the stress field, and its lava at its foot

A large volcano is a load on the lithosphere: beneath and about it the ground is compressed, and in the flexural arch beyond, stretched, so regional graben end at or curve round the great shields and concentric graben ring the loads (Comer, Solomon and Head 1985; McGovern and Solomon 1993; at Alba Patera the volcano's own stress shaped the graben, Cailleau et al. 2003). The model's stress field (F105, F114) is the rise's alone: a fault younger than a shield draws over it as it would over the plains, its floor the dropped flank (F116 has made such faults few on Mars, the faults forming as the dome rises, before the shields). And a shield's outline is its own rule: its lava buries what it grows over (F115), but lava reaching a depression at the edifice's foot, a crater the flank has half buried, would run in and pond across it, which the outline does not do.

Direction: each shield's load in the stress field by the lithosphere's flexure at its time, with the rise's; the edifice's foot from where its lava goes, ponding in the depressions it reaches.

When: with the rise's flexure, the revisit plan's fifth item.

### M20. Io's flow fields are not the floods' process (closed in part by F117)

Io's flow fields share the flow tracer the floods and the shields' flows were built from (F59, F94, F101), but not the floods' later physics: a flow ends where it would climb out of the hollow it has ponded in, rather than filling the hollow and spilling on (F107), and it lays no sheet on the slopes it runs down (F109) nor lobes from each run (F110); there is no flood process on Io at all, its lava being the paterae's flow fields alone. One process on every body (principle 11) asks that Io's eruptions run through the floods' process, with Io's own rates and volumes.

Direction: Io's eruptions as floods from its paterae, by its heat budget's magma, filling and spilling over its plains and paterae as the floods do elsewhere, the flow fields' look (F59 to F62) carried over.

When: its own milestone, after this one.

Closed in part in the forty-ninth milestone (F117): Io's eruptions are floods from its paterae's rims, filling and spilling over its lows and the paterae they reach, laying their sheets and lobes as every body's floods do, the paterae cut through them after; and the form of every flood's flows follows how fast its lava poured out, compound fields where slow, sheet lobes where fast, so that one process gives Io's compound fields and the maria's and Tharsis's sheets alike. Still open: the eruptions' volumes and number are the old fields' shares of their reach, not Io's heat budget's magma (the direction's "by its heat budget's magma", with audit M8); the flood model's sheet is 50 m thick on every body, where a compound field's flow units are a few metres; and one eruption's lobes keep its one age.

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
