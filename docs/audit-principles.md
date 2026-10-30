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

Partly closed in the fifty-second milestone (F121): the episodic lid is no longer a cut on proxies. A dry rocky body's mantle stores its heat under a lid that breaks where the convecting sublayer's stress along a cell as deep as the mantle exceeds the lid's friction over its thickness (after Moresi and Solomatov 1998): Venus's breaks 21 times and last 0.61 Ga; Mars's, the Moon's and Mercury's never. A lid broken within the last half of its body's life is episodic, its retention and resurfacing ages the time since. Still open: the plates' rule and the ice shells' names; the Io regime (S3); the stagnant bodies' histories, which keep the equilibrium heat flux until the composition model gives the Moon its rock (S2); the melt efficiency by name, so that Venus's melt clock still runs on the episodic lid's constant. Since the ninety-eighth milestone (F185) the thrusts and the plains' ridges no longer switch on the stagnant lid's name: an episodic lid has them between its overturns.

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

### W9. Renewed ground carries a noise unlike its surroundings (closed by F196)

The user's note (2026-10-26): on Enceladus's height map the polar and fossil terrains show a finer noise than the ground about them, which looks smoother; they should take the noise that would be in their place, or have it worked as the cratered ground's is. The renewals laid the base relief as the body formed with it (F141), where the ground about them held it relaxed and levelled by its craters' ejecta. Only the ice shell's renewals do this: lanes, coronae, polar and fossil terrains, on every icy body that has them; no rocky process lays the base relief again.

Closed by F196, the hundred-and-tenth milestone: the loop keeps the base relief as it stands, relaxed and levelled by each crater as the ground is, and the renewals take it so.

### W10. Renewed terrains outlined by their margins

The user's note (2026-10-26): each of Enceladus's terrains has folds or cracks running round its perimeter that trace its exact shape, the terrains are near circles and ovals, the cross-folds stand out against the other folds, and the terrains do not meet one another to make a network as Enceladus's do. The fissures turned round the margin about the pole (F43), the outline was a taper about the pole, the margin belt's folds are a set of their own along the edge (F194), and each terrain's fissures stay within its own thinned zone.

Direction: the user chose three steps: the outlines and the fissures from the field; the margin's folds from the fold fabric, turned and raised by the spreading; the fissures running on between the terrains, the younger cutting the older.

Closed in part by F198, the hundred-and-twelfth milestone: the fissures run on until the thick shell stops them and the outline is warped at every scale, so that no ring of fissures traces a terrain. Still open: the folds and the belt along the margin; the network.

Closed further by F199, the hundred-and-thirteenth milestone: the margin's folds are the fabric's own, raised and widened where the edge faces the terrain's spreading across its fissures, so that no set of folds crosses the fabric or runs round the edge. Still open: the network.

Closed further by F200, the hundred-and-fourteenth milestone: each reorientation fractures the plains, and the terrains' fissures running out of their zones end against the fractures, a network that links the terrains. Still open: the fractures are narrow and shallow, faint at 4,096 pixels; the belts that link Enceladus's terrains on its images, kilometres in relief, want a shell thinner where they formed (audit M3).

### W3. Tessera outlines

`TESSERA_SMOOTHING_PASSES` (20 passes on the loop grid) smooths the old crust's height before the tesserae are chosen as its highest share. The smoothing is about two cells, so on Venus it is about 80 km at standard quality and 170 km at preview: the tesserae's outlines change with the quality.

Direction: smooth over a physical length, such as the lithosphere's flexural wavelength or the scale of the mantle's convection. The tessera texture is the lid's strain since the hundred-and-first milestone (F187, T1).

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
| Coronae | `rock_coronae`, `rock_corona` (stages/corona.rs) | closed by F145 and F146, the seventy-first milestone: a diapir's melt pocket under a lid it has thinned, the corona where its loads break the lid in an annulus; its profile the lid's flexure under its cumulates, its rim, its lava and its head's buoyancy, with a drip at its edge while it lasts |
| Tessera texture | `TesseraFabric` (stages/tessera.rs) | closed by F187, the hundred-and-first milestone: the strain the mantle's flow put on the young lid where its drag broke it, folds and graben across the principal strains at the brittle layer's and the elastic lid's buckling and necking wavelengths, as high as the layer's friction and strength let them stand |
| Rises | `dome` (stages/plume.rs) | a plume's buoyancy under the lid's flexure |
| Shields | `shield_height`, `shield_flows` (stages/plume.rs) | partly by F101 and F102, the thirty-fifth milestone: the surface is the eruptions since the edifice last grew by three sheets, compound flow fields from fissures round the caldera, on the flanks and on the rift zone, whose eruptions run the farthest and spread as aprons, the flows narrowing on the steep flanks; the summit sinks in collapses that wander, shrink and deepen; the edifice stretches along its rift and its foot is lobed. Partly by F115, the forty-eighth: its lava buries the ground it grows over, craters and all, as a deposit buries what it falls on. Still drawn: the edifice's profile, one power law at five degrees on every body, where its flows' lengths and the lithosphere's support would set it (Venus's large volcanoes are broad and gentle); the foot's lobing; the vents' shares and the run-outs, where the supply and the effusion rates would set them; the edifice's outline its own rule, where its lava would pond in the depressions at its foot (audit M19) |
| Canyons | `trough_depth_at`, `stress_trace`, `trough_throw` (stages/plume.rs) | partly by F105, the thirty-eighth milestone: the faults start where the rise's stress is tensile, strike across its most tensile direction and die out where the tension fades, so the flanks carry radial graben, the crest a maze and the canyons run out from it; they form as the rise's load grows. Partly by F111, the forty-fourth: a fault grows in straight segments that link, turning toward the stress field and stepping en-echelon where they meet; a graben is as wide as its faulted layer sets it, its throw peaking in each segment and dipping where they link; a fault whose throw would exceed the layer is a rift, its segments half-graben whose master side flips. Partly by F112, the forty-fifth: a graben keeps its throw along its length, a confined fault's flat-topped profile; its walls are rough fault traces, so it pinches and swells, and every fault's trace wanders, self-affine; parallel graben stand apart by a few thicknesses of their faulted layer. Partly by F113, the forty-sixth: dike swarms radiate from the rise and its volcanoes, a graben above each dike. Partly by F114, the forty-seventh: a rise grows along an axis, its volcanoes on it and its stress a pure shear about it added to the load's, so the faults beyond its foot gather along the axis and spare the ground across it. Partly by F116, the forty-eighth: the faults form as the dome rises, the dike swarms with the magma that feeds them. Partly by F147 and F148, the seventy-second: the faults a lid's own strain opens, from its change of radius, are drawn by the same generator, and a wide trough keeps its floor across the turns of its trace. Still drawn: the profile across a trough (a flat floor between tanh walls, no step faults, tilted blocks or raised flanks, audit M23); the faulted layer drawn per fault rather than as a field over the rise; the rise's axis drawn (audit M18); a rift's depth and width from its length, not the strain where it runs, so a long fault far out along a rise's axis can open a canyon alone on the plains; and the stress field's shape, a dome's bending and a load's hoop tension rather than a thin shell solved for the rise's load and support |
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

### T5. The dichotomy as a hyperbolic tangent (closed in part by F168)

The primordial crust's hemispheric dichotomy is half its amplitude times tanh(3 cos θ) about an axis: a smooth slope over thousands of kilometres. On Mars (20 km of crust, some 2.4 km of relief) it does not show under the basins, with no boundary scarp, no fretted terrain and no lowland plains younger than the highlands.

It is also a rule where it should be an outcome. `tectonics_dichotomy_amplitude` gives every stagnant-lid rock body a contrast of 20 km in crust thickness, an episodic lid 5 km, an icy body 2 km and a mobile lid none, by the regime's name (S1), so Mercury, the Moon and Mars carry the same one. Mars's own parameters do not guarantee it (the user asked, 2026-10-16). Both leading explanations are contingent: a giant impact (the Borealis basin, 10,600 by 8,500 km, Andrews-Hanna et al. 2008; an impactor of 1,600 to 2,700 km at 6 to 10 km/s and 30 to 60 degrees, Marinova et al. 2008), which late accretion makes a matter of chance, since it was dominated by a few massive bodies and the median largest impactor on Mars in Bottke et al. (2010)'s Monte Carlo was 1,500 to 1,800 km; or degree-one convection, which needs a weak asthenosphere or a viscosity jump in the mid-mantle (Zhong and Zuber 2001; Roberts and Zhong 2006).

Direction: the dichotomy as an outcome. The giant-impact tail of the late-accretion population in the impact stage, a mega-impact taking its own form (the crust removed inside an oblique ellipse, its ejecta thickening the rest, its melt ponding in the basin), so that some bodies have a hemispheric basin, some none, some an off-centre or elliptical one and some a giant basin short of a hemisphere, as the Moon's South Pole-Aitken is; the boundary scarp worn by the processes that wear a scarp, and the lowlands resurfaced by the lavas and sediments that bury their craters to ghosts; an endogenic route once the interior model resolves the planform of the mantle's convection; and for a locked moon of a hot planet, the planet's heating of its near side while its crust forms (Roy et al. 2014). The amplitude by name goes. A preset for a real body may carry a seed whose draw includes its event: one draw of the model, not a pin.

When: the earlier-phases revisit, with Mars.

Closed in part in the eighty-fourth milestone (F168): the crust starts without a dichotomy and the amplitude by name is gone. A rocky body holds a basin as wide as a hemisphere, as Mars kept Borealis, and the impacts wider than the stamp draws take a form of their own: the crust excavated inside an elliptical rim down to what the melt sheet makes, the excavated crust thickening the crust beyond the rim, the ground where the crust's buoyancy holds it, the floor level, the older relief under the ejecta buried. Some bodies take one, some none, some several, some a basin short of a hemisphere. Still open: the boundary scarp worn by the surface processes (Phase 5), the lowlands' sediments; the endogenic route, degree-one convection, once the interior model resolves the mantle's planform; a locked moon's near side heated by its hot planet while its crust forms (Roy et al. 2014), which the Moon's near and far sides need (deferred by the user); the giant basins' outlines broken up, held back by the user until the processes that wear them are in.

### T6. Grooves as regular troughs

The user's note (2026-10-19): Ganymede's and Enceladus's grooves are troughs too regular to read as faulted ground, and should be more broken up, undulating and organic. The same process should show on Miranda, as a check on the model. Miranda's coronae are belts of subparallel ridges and troughs that Pappalardo et al. (1997) read as tilt blocks from extension and normal faulting above an upwelling (Arden Corona, local extension near 70 percent; warm buoyant ice in Hammond and Barr's convection): the faulting that makes Ganymede's grooves, concentric about the upwelling where Ganymede's follow its lanes. So one process, and Miranda the second body GROOVE_RIFT_LID_SHARE needs.

Direction: the troughs from the faulting itself, tilt blocks and normal faults whose throw and spacing vary along their strike, segments that end, step over and link, so that a groove undulates and breaks; a Miranda preset, its coronae from its upwellings, its grooves' spacing checked against the measured.

When: with the icy moons' tectonics.

Closed in part by F132, the sixty-third milestone: a groove is a chain of fault segments that come and go, step across where one gives way to the next and taper to their tips, dying out before the lane's bounding faults. Still open: Miranda, whose coronae come from upwellings (the user asked whether they should be driven by Ganymede's lanes: they should not; the outline and stress field are a corona's, as on Venus, the fabric this faulting), and the grooves' relief on Enceladus, tens of metres where its sulci's ridges stand up to 1 to 1.5 km, a check on GROOVE_RIFT_LID_SHARE and GROOVE_RELIEF_SHARE on a second body.

Closed in part by F194, the hundred-and-eighth milestone: the check holds, the lanes' grooves standing as Enceladus's ridged plains do (Bland et al. 2007); the sulci of the reference are the margins of its tectonized terrains, and the margin belt is now the margin's brittle layer buckling under the terrain's spreading, as wide as the margin's lid's flexural parameter, at F187's relief (about a kilometre crest to trough), drawn in the lanes' segments rather than rings. Still open: Miranda.

Miranda closed in part by F133, the sixty-fourth milestone: a Miranda preset, and its coronae over the upwellings of the convection its resonance's tidal heating drove in its ice, their belts of the same fault segments concentric with outlines the convection's cells give. The spacing passes its check on a second body: Miranda's coronae's faults come out 4.9 to 5.3 km apart on the lid of its episode, against Arden's 5 km and the 5 to 10 km Hammond and Barr (2014) tabulate for its coronae. The relief does not, under the lanes' law: a share of the spacing gives the belts a few hundred metres where Arden's ridges and troughs stand up to 2 km, the 70 percent extension Pappalardo et al. (1997) find in its tilt blocks. A corona's blocks tilt instead until their faults lock, from the 50 degrees Arden's broke at to the 30 at which friction locks a normal fault, which gives Miranda's belts 1.7 km. The lanes keep their share, calibrated on Ganymede. Still open: one law for both, the relief from the strain the faults take up, a corona's to lock-up and a lane's what its rift opens, with Enceladus's grooves' relief.

### T7. Tiger stripes as a set per region

The user's note (2026-10-19): the height map draws Enceladus's polar terrain as prominent, fine fingerprint patterns, and each occurrence of the terrain carries its own set of four to nine stripes, where the real terrain is a network of interconnected, webbing canyons across all of it. Losing the signature four or five parallel fractures at the south pole is a trade the user accepts for natural terrain. The stripes are also coloured light blue, a false colour.

Direction: the stripes as one network of fractures across the terrain the plume's heat renews, linking, branching and ending where the stress does; their relief at the canyons' scale rather than a fingerprint's; their colour their ice's.

Closed in part by F193, the hundred-and-seventh milestone: each terrain's fissures are one network, the stripes its greatest, saturated at the shell's thickness and linked by branches that leave fissures as the stripes' forks do, each past episode's meeting one another as the present's do; the folds take their strike from the nearest fissures and stand at the funiscular terrain's 100 m, their depth rising and falling along them; the false colour is gone. Still open: the network reads faintly at a global view, where the reference's dominant belts are Enceladus's sulci, whose relief is T6's open part.

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
| GRABEN_WALL_SEGMENT_WIDTHS, GRABEN_WALL_JOG, GRABEN_WALL_RAMP_WIDTHS, GRABEN_WALL_ROUGH | a wall's straight pieces, their jogs and roughness (F149): Pluto's and Charon's chasmata run straight for a few widths between jogs | the bounding fault's own segments and relays, from the fault's growth by linking (F111), rather than a second set drawn along the wall |
| RAY_LENGTH_COEFF_KM, RAY_LENGTH_EXPONENT, RAY_SKIN_DEPTH_M | Elliott et al. (2018)'s median lunar ray lengths, 6.59 R^1.27 km, for their 10 cm skin (F91) | the skin and fragment population per body: the length goes as the skin depth to about -1.16, which an icy surface's lag already sets |
| CURTAIN_SPREAD, CURTAIN_SLOPE, CURTAIN_ORDERS | the Moon's rays: 10.7 a crater, 22 degrees wide at two radii (Elliott et al. 2018) | the target's relief at the rim, which makes the rays (Sabuwala et al. 2018: the count goes as the crater over the relief's wavelength) |
| TERRACE_FULL_COLLAPSE | a stated choice: a lunar crater's wall terraced in full from about 60 km (F92) | check the share of the wall terraced along a sweep of lunar and Mercurian craters |
| POROUS_DISSIPATION | Choblet et al. 2017's model of Enceladus's core, 10 to 30 GW: 28.5 GW, which with its shell's own 0.6 GW (F75) holds its shell at 22.7 km (refitted from 15 GW when the shell's share went from the class k₂/Q's 15 GW to F75's) | a constant of the dissipation per unit of porosity; passes its check on Mimas since F138: the core dissipates only while its pores hold water, and Mimas's froze within 0.3 Gyr of its forming and have not thawed under its ocean some 20 Myr old, so that its shell is 24 km thick over the young ocean Lainey et al. 2024 find (F135 had given it 228 GW and a shell 1.7 km thick) |
| RING_STRAIN | the strain of Tyre's and Callanish's inner graben (Singer et al. 2023) | check against the graben depths of Valhalla and Ganymede's furrows once measured |
| PIT_DIAMETER_SHARE, PIT_GROWTH, PIT_DEPTH_SHARE, DOME_PIT_SHARE, DOME_RISE | Ganymede's pits and domes (Alzate and Barlow 2011; Bray et al. 2012; Schenk 1993) | check against Callisto's dome craters and Ceres's pits (Occator, Dantu) |
| MANTLE_SPREAD | a stated choice: the length over which a mantle evens the ground, twice its thickness | check against the rim profiles of Enceladus's muted craters (Kirchoff and Schenk 2009) |
| FLOW_NARROWING_SLOPE | Tharsis: a 5 degree flank narrows a flow to about a third of its width on the plains (Peters et al. 2021: the narrowest flows on Olympus Mons's flanks) | a flow's width from its lava's yield strength and its effusion rate on the slope (Hulme 1974) |
| FLOW_STEEP_SLOPE | a stated choice: a degree, above the plains' slopes, so that the plains' floods and Io's flows keep the steering they had | the slope at which gravity's push across a front outweighs the bumps it rides over; check against flow paths on a shield's flanks |
| ICY_FLUX, ICY_FLUX_RISE | Nesvorný et al. 2023's comet rates on the moons of Jupiter, Saturn and Uranus, each over its planet's focusing (F135) | fitted on three planets; within a factor of 1.3 on six moons; check against Triton's and Pluto's craters, beyond the planets it was fitted on |
| ICY_EARLY_TOTAL | Wong et al. 2023's comets of 1 km over the history of Mimas and Enceladus, over Nesvorný et al.'s rates (F135) | passes on Mimas, Tethys and Dione (F136): with Bottke et al. 2024's sizes and nothing fitted, their counts from 20 km up within 1.4 times; Rhea's at 0.6, its record starting with its crust; check Iapetus's oldest plains and Callisto's |
| ICY_IMPACTOR_SFD | Bottke et al. 2024's size distribution, read off their Fig. 1: what their craters on Iapetus and Phoebe require, to 32 km, and their collisional model beyond | a model fitted on two bodies; checked on Saturn's mid-sized moons' large craters (F136), within 1.4 times |
| CRATER_PI_SCALING, TARGET_STRENGTH_ICE, TARGET_STRENGTH_ROCK | cold ice's (Holsapple 2022, as Bottke et al. 2024 adopt it for the large moons); rock's from Holsapple 1993's soft rock, a stated choice | one set for every icy body, where Bottke et al. use a weaker one on the small moons; check against Ganymede's and Callisto's counts, which come out twice their ancient terrains' at 30 to 200 km |
| ROUGHNESS_RMS_SHARE at the largest scale | a saturated cratered surface's rms over its craters' depth, set for the sub-grid roughness (F54) and carried to the primordial relief's largest octave (F137) | the largest scale is not saturated with craters as wide as the largest a body holds; checked by the spectra's shape on Mimas, the Moon and Mars, where Mercury carries too much at its largest scale (0.56 of its relief against its map's 0.39) |
| WARM_START_SHELL | a stated choice: an icy moon's H2O layer liquid after accretion under a crust of 1 km where a warm moon's heat would hold no floating shell (F138) | the early thermal state from the accretion's heat, the short-lived isotopes and the moon's time of forming; Enceladus keeps its ocean only from a warm start, its eccentricity of 0.0047 never melting it frozen through |
| CORE_THERMAL_DIFFUSIVITY | rock's k / ρ c_p, 2.5 W/m/K over 2.5 MJ/m³/K (F138) | the porous interior's own thermal history, its pore water's latent heat and the hydrothermal circulation that carries Enceladus's heat (Choblet et al. 2017), for the time a thawing interior takes |
| PUMPED_ECCENTRICITY_MAX | a safety cap: 0.2, near where Nakajima et al. (2019) find Enceladus's orbit crossing Tethys's (F138) | each resonance's own largest eccentricity, from the dynamics of its pair |
| LANE_RENEWAL_GYR | Ganymede's two thirds renewed by its resonance episode (F140), as ICE_RESURFACING_FLUX was calibrated | a renewal by the heat the lid carries above the floor, linear in time; fails its check on Mimas: its young ocean's 13 mW/m² over 33 Myr renew a tenth of it in lanes, where Lainey et al. 2024 find no surface expression; the lanes want the strain of the lid, and a melting ocean contracts the body |
| Flexural support (F139) | a thin plate's flexural wavelength 2π (D / ρ g)^¼ (Turcotte and Schubert 2002), the share settling from two thirds to four thirds of it | the viscoelastic relaxation of a lid whose viscosity falls with depth, which a sharp elastic thickness stands for; Enceladus's relaxed craters (Bland et al. 2012) want lids of a kilometre or two |
| Resonance supplies | the heat stated for the episode, or what its eccentricity gives a warm moon (F138): Mimas's 16 GW from Meyer and Wisdom 2008's equilibrium heating scaled to Lainey et al. 2017's Saturn, Miranda's 42 GW and Ganymede's 0.67 TW from their eccentricities | the planet's tidal torque on the moons of the resonance chain, the supply each resonance can deliver; Uranus's torque at a Q of 17,500 delivers Ariel and Umbriel some 6.5 GW, a sixth of Miranda's |
| FINE_TEXTURE_SHARE | a stated choice: a renewed terrain keeps the base relief below the length it levels over, and below its fabric's spacing 0.4 of its amplitude, the fabric overprinting what is finer (F141); chosen by eye with the user for the smoother look of new ice between its grooves | the roughness of new ice at scales below the fabric, from the fabric's own process; checked on Mimas, Enceladus and Miranda against the Gaussian look it replaces |
| ELASTIC_ISOTHERM_ROCK | 650 °C (F144): McKenzie, Jackson and Priestley 2005's 600 °C for the oceanic lithosphere, Davaille et al. 2017's 650 °C for Venus, dry diabase's creep at 1e-17 /s under 100 MPa (928 K); checked on Venus's 29 ± 6 km | the isotherm from each body's rheology and its loads' time, wet or dry (the open rheology item) |
| CORONA_HEAD_RADIUS, CORONA_FEED_GYR | Venus's coronae (F145): some 750 of them, about 200 km across, with the heat share and flux of Smrekar and Stofan 1997; the feeding within their cases' 50 to 430 Myr | the diapirs' size and life from the boundary layer they rise from, as the plumes' should be |
| CORONA_UPWELLING_DRAW | Venus's coronae gathered in its upwellings' region, at its hot spots and along its rifts (Stofan et al. 1992: 21, 68 and 11 percent; F145) | the share of the boundary layer's instabilities the large upwellings sweep in, from the flow's speed against the diapirs' rise |
| CORONA_HEAT_SHARE, CORONA_BUOYANCY_FLUX | Smrekar and Stofan 1997: upwelling at coronae some 15 percent of Venus's heat loss, 0.30 Mg/s for an average corona (F145) | the diapirs' heat from the mantle's own convection |
| CORONA_HEAD_ASPECT, CORONA_HEAD_EXCESS, PLUME_HEAD_SHARE | stated choices (F145): a head's buoyant column as deep as it is wide, 100 to 300 K hot (Dombard et al. 2007), a plume's hot head half as wide as its dome | the head's shape and heat from its ascent |
| CORONA_FAULT_STRENGTH | 20 MPa within Dombard et al. 2007's 10 to 80 MPa (F145) | the fracture of the lid's top by depth, its friction against the bending stress |
| CORONA_INTRUSIVE_SHARE, CORONA_VENT_SHARE, CORONA_STACK_DENSITY, CORONA_RING_SHARE | intrusion four fifths of the melt (Dombard et al.; White, Crisp and Spera 2006), half the lava at the vents, cumulates 200 kg/m³ denser than the crust (Dombard et al.), rings 0.15 of the radius (F146) | the lava's run over the ground it floods; the lower crust's flow under the load |
| CORONA_RIM_VARIATION, CORONA_HEAD_WANDER | stated choices (F146): the rim's thickening varying by 0.7 round the edge, the drip in arcs, the outline wandering by 0.15 of the radius | the head's own shape and the crust's flow round the pocket, the arcs from the lid's tearing (Davaille et al. 2017) |
| FAULT_LAYER_PRESSURE, LID_LAYER_SPREAD | Charon's Serenity Chasma, 40 to 50 km wide (Beyer et al. 2017), over ice at 0.29 m/s² (F147); checked on Tethys's Ithaca Chasma, 100 km wide, given 87; the spread a stated choice | the faulted layer's base from the lid's strength envelope, friction against creep at the strain rate the lid stretches at |
| LID_FAULT_MIN_WIDTHS, LID_FAULT_MAX_LENGTH, LID_LAYER_RADIUS_SHARE | restricted faults run two to ten times as long as they are high (Nicol et al. 1996; Soliva and Benedicto 2005), half a great circle, a fifth of the radius against Ithaca's 0.19 (F147) | the faults' lengths from their linkage as the strain grows, the layer's bound from a thick shell's stresses |
| lid compliance (density over stiffness against rock's) | Charon's chasmata, at most 0.016 of their length deep (F147), between the stiffness alone and the frictional stress at a depth | a fault's slip per length from its stress drop and the rock's stiffness, checked on Pluto's and the Moon's graben |
| ICE_HP_FREEZING_CHANGE | between ice III's 2.6 percent at 350 MPa and ice II's 6 to 8 (F147) | each high-pressure ice's density against the water's at the ocean floor's pressure and temperature |
| LID_LAYER_MEDIAN_SHARE, LID_LAYER_SPREAD | a lid's faults cut to half the deepest at the median, by 0.5 in the logarithm (F149): Charon's chasmata 10 to 50 km wide, Serenity the widest; Pluto's fossae a few km to 20 | the depths faults nucleate at and grow to within the brittle layer, from the layer's strength with depth |
| RISE_MEMBRANE_SHARE | a rise's crest tension the membrane stress of its uplift on the sphere, E h / R (F150; Turcotte et al. 1981) | the thin shell solved for the rise's load and support (the direction of the canyons' audit row) |
| UPWELLING_LIMB_REACH, UPWELLING_LIMB_TENSION, UPWELLING_LIMB_WIDTH | Venus's three rises joined by rifts along their upwelling limbs (F150): plumes within four radii together on one upwelling, pulled apart by 0.8 of the lesser crest's tension over a rise radius | the mantle's convection pattern, its upwelling sheets and the traction they put on the lid |
| the GPE margin's plane strain over the flexural wavelength | a stated choice (F150): a plateau's margin spreads down its slope, the stress along the contour Poisson's share of the stress across it | the membrane stress of the relief solved on the shell (Banerdt 1986) |
| FAULT_FLANK_SHARE, GRABEN_HALF_SHARE, GRABEN_THROW_ROUGHNESS | stated choices (F149): a footwall a tenth of its throw up, two graben in five half-graben, the throw rough by 0.2 along the fault | the footwall's flexure from the elastic lid; the floor's tilt from the master fault's share |
| THRUST_STRIKE_SCATTER, THRUST_FIELD_FOLLOWING | stated choices (F154, F155, F157): the field decides a thrust's strike as its differential stress passes a tenth of the lid's strength in compression or of the mean squeeze, whichever is greater; a trace turns toward the strike by 0.35 a step | the lid's strength's own scatter from its fabric, its craters' rims and old faults, which the thrusts would follow in a field alike every way |
| SCARP_REFERENCE_STRAIN | refitted (F155): Mercury's thousand-kilometre scarps at the 0.4 percent its thrusts take up beyond its lid's elastic strain | the longest thrust from the depth of the layer it cuts and the linking of its segments as the strain grows |
| the plains at their strength (F157) | a stated choice with nothing fitted, replacing RIDGE_DRIVE_SHARE (F156): the plains stand at their strength at least, brought there by drives the field does not hold, Venus's climate warming its plains after their renewal (Solomon et al. 1999) and their own cooling; Venus's ridge midpoints 84 percent below its mean height against 93 measured | the drives themselves: the climate's history once the atmosphere has one, the plains' cooling as a thermal stress |
| PLATE_CELLS_PER_PARAMETER, PLATE_FACE_RANGE, LOAD_EPOCH_GYR, LOAD_MIN_INCREMENT | numerical choices (F157): two cells across the plate's flexural parameter, faces 16 to 128, the field read every 0.3 Gyr, a load once 10 m of lava lies anywhere | checks: the shell against the Hankel transform's solutions within a few percent; Venus's thin lids under-resolved at 128 |
| the loads' graben within the faulted layer (F158) | a stated choice: the loads' stretch is small and none of its faults rift | the throw's partition between graben and rifts from the strain the loads open |
| SYSTEM_SPACING, SYSTEM_ARREST_SHARE, SYSTEM_INDEX_WIDTH | stated choices (F160): a neighbour 0.6 to 1.6 of its host's reach aside (the Tharsis plateau's ridges some 30 km apart, Watters 1991), an arrest within a quarter of the reaches, an index of 120 km cells | the spacing from the stress shadow a thrust casts over its faulted layer, as the cracks' spacing from theirs (CRACK_DENSITY_FACTOR) |
| the spread toe (F172) | a stated choice: an edifice past what its crust holds pushes as much of its toe out as the lava it would have stacked above, which with the lava makes its aureole | the flank's spreading on its base by its stress and the base's strength |
| SCARP_WARP, SCARP_WALL_SLOPE_DEG, AUREOLE_LOBES, AUREOLE_RIDGE_SPACING, AUREOLE_RIDGE_SHARE (F172) | stated choices: the scarp's variation round the edifice, its wall's slope, the aureole's lobes and their ridges, after Olympus Mons's | the lobes from the ground's fall, the ridges from the slide's flow |
| CONVECTIVE_TOPOGRAPHY (F173) | the mean of Guimond, Rudge and Shorttle 2022's eleven stagnant-lid simulations, 1.5 to 3.2, applied beyond their Rayleigh numbers in the boundary-layer form that held across them | the support from the flow's own buoyancy through the layer |
| CONVECTION_CELL_WAVELENGTH, CONVECTION_SPECTRUM_WIDTH, CONVECTION_ONSET (F173) | stated choices: the cells as wide as convection's onset wavelength between free-slip boundaries, their power spread about their degree, a layer convecting once several sublayers deep | the planform from the mantle's own heating, rheology and geometry, rolls in a thin shell and plumes in a thick one |
| the convection's pattern (F173) | a draw from the body's seed: the flow is chaotic, its scale, strength and renewal the mantle's | the pattern from the mantle's evolution itself |
| AUREOLE_MOBILITY (F172) | Scheidegger 1973's trend of a rock avalanche's run against its volume as it is quoted, unverified against the paper | a run-out from the slide's flow over its base |
| MAGMA_PROVINCE_FLOOR (F171) | Mars's Early Hesperian plains, some 800 to 1,000 m over most of its northern lowlands (Head, Kreslavsky and Pratt 2002): a province's centre four times the rest | the province as its cumulates' heat in the thermal history (M33) |
| CRUST_FLOW_RELIEF, CRUST_FLOW_MIN_VISCOSITY (F170) | stated choices: the lower crust's viscosity read at the stress a kilometre of relief drives, and never below 10¹⁸ Pa s, a lower crust near its solidus | the stress each band's relief drives in the power-law creep, and the crust's melt |
| CRUST_MAX_THICKNESS (F184) | a stated choice, 65 km, the most a basaltic crust stands before its base turns to eclogite and founders | the transition's depth on the crust's own geotherm, from the gabbro-eclogite phase boundary and its kinetics (Semprich et al. 2025; Namiki and Solomon 1993) |
| the relaxing crust's thickness (F170) | a stated choice: the crust's tenth thickest sets the body's relaxation, band by band, not each cell's own | a relaxation that varies with the crust beneath each cell |
| the rocky relaxation's scope (F170) | a stated scope by regime (S1): a mobile lid's crust is its plates', and Io's heat pipe is not the conductive lid its thermal history describes | the plates' keels (Phase 6) and Io's heat-pipe lid in its thermal history |
| the levelling width (F169) | a stated reading: an ejecta deposit levels the ground over the width of the fresh crater as deep as it is thick, no farther than LANDSLIDE_REACH's longer end times its thickness, the constants of the crater scaling and of the landslides' run-out | the ejecta's flow on landing, ballistic sedimentation's reach |
| MEGA_FLOOR_CRUST_SHARE (F168) | a stated choice: the crust a mega-basin's melt sheet makes, half the starting crust, between Mars's northern lowlands (some 32 km against 58, Neumann et al. 2004) and South Pole–Aitken's thinner floor | the crust the melt sheet's volume and differentiation make |
| MEGA_WALL_SHARE (F168) | a stated choice: a mega-basin's wall, the outer twentieth of its radius, a boundary scarp a few loop cells wide | the wall's collapse at the crust's strength |
| a mega-impact's ejecta (F168) | stated choices: all the excavated crust lands beyond the rim, as a basin's ejecta thins with distance, with no convergence at the antipode; the relief it buries measured against the ground low-passed over the basin's radius | the ejecta's ballistic ranges on a sphere |
| ELLIPTICAL_THRESHOLD_EXPONENT (F167) | the line of Collins et al. (2011)'s Figure 1 as digitised: the threshold angle 90° times the cratering efficiency to the power of −0.73 | the paper's own fit, if it states one |
| ELLIPTICITY_AT_THRESHOLD, ELLIPTICITY_GROWTH (F167) | the elongation's growth below the threshold angle in one suite of simulations on a flat aluminium target (Collins, Miljković and Davison 2013) | a growth law across efficiencies and materials |
| one impact speed for the efficiency (F167) | a stated choice: every population's craters read their efficiency at the body's mean impact speed, though comets strike faster and a moon's debris slower | the efficiency by population, once a crater records the population it came from |
| ROOT_MIN_CONTRAST (F166) | a stated floor: the least density contrast the mantle pushes back with beneath an edifice filling its sag, 100 kg/m³, below any lava a rocky crust buoys up | none needed while lavas stay lighter than the rock beneath |
| the intruded surplus as underplating (F166) | a stated reading: the magma a full edifice cannot lift spreads beneath the rise, lighter than the mantle it replaces, which the dome holds, rather than loading the lid beneath the edifice | the rise's dome from the underplating it holds |
| the plumes' heads under their rises' provinces (F165) | a stated choice: the plumes melt over the area their rises' Gaussian provinces cover, 3.4 percent of Venus and 16 percent of Mars | the plumes' heads from their buoyancy flux |
| the overturn's resurfacing while the ambient mantle melts after it (F165) | a stated choice: the melt that comes before the ambient mantle stops melting after an overturn is the event's, laid by its fill | the event's fill from that melt |
| the regional stress a volcano's weight turns dikes against (F165) | a stated choice: the lid's yield in tension at the rise's peak, falling as its growth wanes, F119's own reading | the rise's stress as the lid's field reads it |
| the shields' supply read against the ground the loop begins on (F164) | a stated choice: the share each rise's shields take of a step's budget is read from the vent weights over the first ground, so that the schedule is known before the loop | the share read from the ground as each step finds it, the shields growing in the loop |
| a full edifice's surplus intrudes, its flanks renewed at the pace of its supply (F164) | stated choices: the magma a full edifice cannot lift intrudes beneath it, and eruptions too small to count keep its flanks young | the share of a full volcano's magma that erupts on its flanks and past its foot, and its intrusions' load (M19) |
| a flow's lava spread evenly along its run (F163) | a stated choice: the share of a shield's lava past its foot, 0.41 for a round edifice, from F101's vents and run-outs, each flow's area spread evenly along its run as its breakouts are drawn by the lengths of the lobes laid | the lava along a flow from its cooling and the slope it runs down |
| RIDGE_REFERENCE_LENGTH, RIDGE_MIN_SIZE | stated choices (F159): a ridge 100 km long stands RIDGE_HEIGHT, its arch and crest no narrower than a quarter of theirs | the ridge's fold from the thrust beneath it and the layered plains' thickness |
| GRABEN_BLOCK_LAYERS, GRABEN_INTERIOR_SHARE | a wide trough's floor in blocks a quarter of its faulted layer wide, taking 0.3 of its throw (F152): Mandjet's tilted blocks 5 to 10 km on a layer of 20 to 40 km, Beyer et al. 2017; the bounding faults carry most of a graben's throw | the blocks' width and rotation from the domino kinematics of the extension over a listric master (Gibbs 1983) |
| GROOVE_END_TAPER_WIDTHS, GROOVE_LANE_TAPER_POWER | a lane's end in the old ground closes linearly over three widths (F142), as a fault's displacement falls to its tip; Ganymede's lanes pinch out in long points | the taper from the rift's own mechanics, the strain falling along the lane's last reach |
| DEBRIS_FRAGMENTS, DEBRIS_BREAK_KM, DEBRIS_SLOPES, DEBRIS_ECCENTRICITY | Mimas's craters beyond the comets' (F136) | checked on Tethys and Dione, within 1.4 times; the scaling with the square of the impact speed, which spares Miranda, is a stated choice with no derivation: standard ejecta scaling makes the mass an impact of a given crater throws out of the moon independent of its speed |
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

In the sixty-ninth milestone (F138) the moons' eccentricities and interiors evolve together, and the episodes of Mimas's and Miranda's heating emerge from their resonances and their ice; Enceladus's do not. Held in its resonance with Dione, its heating stays at the equilibrium the process finds for it, 29 GW, from its warm start to the present: the published dynamics give it no excursion above equilibrium in its past resonances (Meyer and Wisdom 2008 find its capture into the 2:1 e-Enceladus resonance close to equilibrium after a brief capture into the 4:2 mixed resonance, and exclude its having occupied the 3:2 with Mimas), and Goldreich et al. (2025)'s limit cycle, which heats it to some 290 mW/m² every 13 Myr, would relax its plains throughout their history. A past resonance stated for it would heat its plains as Miranda's heats its terrain, but none is known to the dynamics. Since F139 its relaxation also wants lids of a kilometre or two, as Bland et al.'s heat flows give, and the ground beside its terrains, under an episode's thinned shell, relaxes less than F66 made it: its plains keep sharper craters than its map, until a heat history that reaches those flows is found.

Closed in part by F195, the hundred-and-ninth milestone: the Phase 4F gate check found Enceladus's relief at four times the measured, its primordial relief held whole by a lid that held any load shorter than its support, whatever its size; a lid now holds no more than its strength bears, and the warm ice under it flows from under the rest, which halves the relief (-3.2 to +2.8 km) and takes the broad swells. Still open: its craters and its bands under 80 km, in cold ice, which want the heat history; Goldreich et al. (2025)'s limit cycle is the candidate the dynamics give.

Closed further by F197, the hundred-and-eleventh milestone: the lid that holds a load is as old as the load, so that the relaxation law lets Ganymede's craters relax at the heat flows Bland et al. (2017) find for them, and Enceladus's relief eases again (-2.6 to +2.2 km). The limit cycle cannot arise in the model: its heat is made in its porous core (F68, Choblet et al. 2017), whose dissipation does not fall as the shell thickens, and which takes up any supply steadily (116 GW, Goldreich et al.'s resonance-locking supply, holds a steady 5.7 km shell at an eccentricity of 0.0091, against the 20 to 25 km measured). Still open: the plains' craters, which want heat flows above 150 mW/m2 where its lid stays near 7 km; a heat history that brings them under the porous core (its hot spots moving, or its heat varying) is the direction.

### M3b. The coronae's dark belts

Voyager's coronae on Miranda are darker outside and bright within: Inverness's banded outer material and Arden's dark bands are the darkest ground on the moon, and the USGS map splits the coronae into an outer and an inner material. Since F143 the belt keeps the old crust's tone with clean ice on its scarps and the core is clean, so the belt is as dark as the cratered plains and no darker.

Direction: the dark component of the Uranian moons' ice, which darkens the plains to their measured 0.20 (Miranda's show 0.58, F126) and which the belt's faulting would expose or concentrate; with the materials of Phase 5.

### M4. Enceladus's old fractures

Enceladus's cratered plains are cut by fractures. The class rule's Love number gave its shell an h₂ ten times Beuthe's (2018) and a diurnal stress of 194 kPa, which cracked the whole shell; with the shell's own Love number (F75) its uniform shell is stressed by 20 kPa (Beuthe: about 15 kPa in a uniform 23 km shell), under the 40 kPa at which ice cracks, and only its thinned south polar shell and the zones its past episodes thinned are stressed past it (F76). Its old plains therefore lose the cracks the diurnal tide had drawn there.

Direction: the stresses that crack a shell beyond the diurnal tide: the tension of a thickening shell as its ocean freezes and pressurizes (Manga and Wang 2007; Rudolph et al. 2022), which each episode's end and each fall of the heat flux bring; true polar wander as the body reorients to each new episode's pole (Matsuyama and Nimmo 2008); and non-synchronous rotation. Kihoulou et al. (2025) find extension in the upper shell as Europa's shell thickens after a pulse of eccentricity, and compression as it thins.

Since F82 (2026-10-15, the twenty-fourth milestone) a shell that thickens stands at the cracking stress in tension. Enceladus's mean shell, held by its porous core's steady heat (F68), thickens by only 0.6 km over its life, some 650 Pa of tension, so its plains still do not crack; the thickening that would crack them is local, the zone refreezing at each focusing episode's end, which F82's global history of the shell does not see.

When: with a local history of the shell under the focusing episodes.

Closed in part by F200, the hundred-and-fourteenth milestone: each reorientation to a new episode's pole stresses the shell by its figure's change, some hundreds of kilopascals, and its first fractures cut the plains, a few long ones to each reorientation, narrow and shallow for the small strain they take up. Still open: the local refreezing's tension at each episode's end.

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

Since F202 the young flows' SO2 rises in a plume at their fronts and falls in the Prometheus class's ring about them, in place of the halo that hugged their margins; the rings are white wherever their fall outpaces the sublimation, crisp and opaque, the frost's look this entry names.

Since F203 the frost coarsens as it ages and an old frost lets the ground through; since F204 the plains are Io's measured pale yellow rather than a mustard, so the frost stands out from them by what Io's does. Still open: the fed rings are crisp, white wherever their fall outpaces the sublimation, some 35 of them, and the plains carry no white of their own: the frost cycle this entry's direction names.

Closed in part by F205, the hundred-and-eighteenth milestone: the first of the missing pieces above, SO2 laid broader than a plume's ring, its gas settling across its region as the atmosphere collapses each night. Frost shows over 40 percent of the surface, in broad soft fields where the fronts gather and most rings merge with them, and beyond 30 degrees under the plasma's alteration; the body's mean stays Io's measured colour. Tried first and not kept (2026-10-27): each frost subliming against the atmosphere its region's frost holds up, without the albedo feedback; frost rose from 4.4 to 4.6 percent, the gas the model then laid being far under the sublimation away from the rings. Still open: the frost's temperature through the day and the condensation that goes to the coldest ground, and the fields' fine structure.

Since F206 the frost does not hold on the dark lava still showing, which the sun warms far above it, so the flows stand out within the frost fields; the frost's temperature over the bright ground is still its latitude's and slope's alone.

### M8. Io's flow coverage

Io's flow fields cover 28 percent of its surface, its plains 66 (Williams et al. 2011). F85's fields cover 20 percent and the plains 70: the fields' areas are stated shares of their reach, and the number of fields three generations a patera. Since F87 and F88, with the lava that ponds in a patera counted as its floor, the fields cover 16 percent and the plains 79. Since F117, the flows floods from the paterae's rims that run down and pond in the closed lows they reach, they cover 19.9 percent and the plains 74.5, the patera floors 2.3.

Direction: the fields' areas from the lava erupted over their lives, the eruptive share of the magma rate (Spencer et al. 2020: most of Io's magma intrudes) over a field's thickness, so that the coverage follows from the heat the lava carries out.

When: with S3, the regime's switch.

Closed in part by F201, the hundred-and-fifteenth milestone: the fields' lava is the heat budget's eruptive share, shared among the paterae by their floors' areas, each field growing over the years its patera's share takes to erupt it; they cover 27.4 percent of the surface, lava showing over 15 percent of it. Still open: a field's area for its volume is the floods' 50 m sheet (M20), and the eruptive share the floods' one in ten (Crisp 1984), not Io's own (Spencer et al. 2020).

### M9. A rise stays where its plume rose

Mars's rise came out at 83 degrees south, its shields round the pole. A load the size of Tharsis turns the body until it lies on the equator, which is why Tharsis does (true polar wander: Bouley et al. 2016). ExoGen turns a body for a light terrain, Enceladus's thinned pole going to the pole (F53), but not for a heavy one.

Direction: the body's reorientation from its loads, positive and negative, the rise carried to the equator as it grows, and everything older turned with it.

When: the earlier-phases revisit, with Mars's rise.

Closed in part by F104, the thirty-seventh milestone: the body turns until its axis of greatest moment of inertia is its spin axis, from its rises' mass against the rotational bulge its elastic lithosphere holds frozen, so Mars's rise is carried from 83 degrees south to the equator and Venus's three rises to within 12 degrees of it. Still open: the negative loads stay with F53's own turn (Enceladus's thinned pole), and the rises' shields and floods add no mass of their own; the latitudes things had before the turn will matter when the climate's processes come in Phase 5.

### M10. Venus's deformation belts and rifts

Venus's plains are cut by belts of ridges and fractures and by rift zones (chasmata) that join its highlands and coronae into a network, the bright web of the radar maps; coronae string along the rifts. ExoGen's Venus has none: 506 coronae stand alone as thin rings, the tesserae are plateaus of noise, and the plains between are smooth under their wrinkle ridges.

Closed in part by F145 and F146, the seventy-first milestone: the coronae come of the mantle's diapirs, melt pockets under a lid their heads thin, gathered about the rises' upwellings, their profiles the lid's flexure under their loads. Closed further by F150 and F151, the seventy-third: the rifts join the rises along the upwelling limbs between them, 23,000 km of them, and the coronae are drawn along the limbs and the rift zones, two in five on the limbs. Still open: the share on the rifts themselves (M25); the deformation belts of ridges and the wrinkle ridges in the lows, which the same field's compression will place (M72c); the novae still a random share of the coronae; Mercury's three faint domes.

Closed further in the ninety-eighth milestone (F185): the flow drags the lid, squeezing it over its downwellings by up to 140 MPa on Venus; an episodic lid has thrusts, which form over the downwellings of their time, nine tenths of them in the first hundred million years after the overturn, on the plateaus the flow piled; the plains' ridges lean to the downwellings; the novae are the young coronae whose rising domes stretch past the lid's strength, 17 percent of them against some 12 percent catalogued. Still open: the lowland ridge belts (Atalanta, Lavinia), which the drag through a thickened lid no longer drives past the strength; the thrusts on the plateaus as the tesserae's own fabric (W3); the coronae on the rifts (M25); Mercury's three faint domes.

Direction: the lid's stress from the mantle's convection, rifting where it pulls apart and folding where it pushes, coronae where plume heads meet it, as annuli of fractures and ridges round a dome or sag, and tesserae as the oldest crust deformed in the belts.

When: the earlier-phases revisit, with Venus.

### M11. The rocky bodies' colours

Mars's albedo is one grey, Venus's another, and Mercury's flooded floors are flat tan patches in dark halos, where its fresh craters' ejecta should be bright. Mars's red is the oxidation of its dust and Mercury's tones its materials' and their maturing; the user left Mercury's colour for a later pass (Phase 3).

Direction: the weathering and transport processes of Phase 5 (oxidation, dust, maturity), with the composition of each unit.

When: Phase 5, or a colour pass the user asks for.

The silicates' maturing done (F213, 2026-10-28): a silicate's mature colour is its fresh one weathered by its iron, as intensely as its body's weathering runs against its gardening, which darkens and reddens Mercury's iron-poor crust as far as the Moon's highlands, by the measured amounts. Mars's oxidation and dust, Titan's haze fallout, the tholins of the outer bodies and Miranda's darkening wait for Phase 5.

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

### M19. A volcano's load in the stress field, and its lava at its foot (closed by F161, F162, F163, F166)

A large volcano is a load on the lithosphere: beneath and about it the ground is compressed, and in the flexural arch beyond, stretched, so regional graben end at or curve round the great shields and concentric graben ring the loads (Comer, Solomon and Head 1985; McGovern and Solomon 1993; at Alba Patera the volcano's own stress shaped the graben, Cailleau et al. 2003). The model's stress field (F105, F114) is the rise's alone: a fault younger than a shield draws over it as it would over the plains, its floor the dropped flank (F116 has made such faults few on Mars, the faults forming as the dome rises, before the shields). And a shield's outline is its own rule: its lava buries what it grows over (F115), but lava reaching a depression at the edifice's foot, a crater the flank has half buried, would run in and pond across it, which the outline does not do.

The ground does not bend under its loads either: the floods' loads stress the lid (F157) but leave the ground as it lies, and a shield stands on level plains, where the lid would sag into a moat about it and rise in a bulge beyond (Olympus Mons's moat, partly filled by its aureole; McGovern and Solomon 1993), and a filled basin would sink further under its fill, which the lava's subsidence, a twentieth of its thickness (LAVA_SUBSIDENCE), stands in for.

Direction: each shield's load in the stress field by the lithosphere's flexure at its time, with the rise's; the edifice's foot from where its lava goes, ponding in the depressions it reaches; and the plate's deflection under every load (F157's) laid into the ground, a moat about a shield and the maria sunk under their fill.

Closed in part in the seventy-eighth milestone (F161, F162): the shields' loads join the floods' in the field, and every load's deflection is laid into the ground, the edifices filling the sag beneath them so that each stands in a moat with a bulge beyond, and the maria sink under their fill. Still open: the edifice's foot, whose lava would run into its moat and the depressions about it and pond there, as Olympus Mons's aureole and lava plains fill its moat; and the root, the lava that fills a shield's sag, left out of its load.

Closed further in the seventy-ninth milestone (F163): the lava a shield's flows carry past its foot, some two fifths of its lava, starts from along their runs beyond it and runs on down as a flood's does, ponding in the moat, which it mostly fills on Mars, and spilling past it. Still open: the root, the lava that fills a shield's sag, left out of its load.

Closed in the eighty-second milestone (F166): beneath an edifice the lid sinks into the edifice's own lava, which fills the sag as far as the edifice is thick, so that the mantle pushes back there by the difference between the rock's density and the lava's. Mars's shields' moats deepen to some 3 to 4 km, which their lava only partly fills (M30); the magma a full edifice cannot lift spreads beneath the rise as underplating rather than loading the lid beneath it.

When: with the rise's flexure, the revisit plan's fifth item.

### M20. Io's flow fields are not the floods' process (closed in part by F117)

Io's flow fields share the flow tracer the floods and the shields' flows were built from (F59, F94, F101), but not the floods' later physics: a flow ends where it would climb out of the hollow it has ponded in, rather than filling the hollow and spilling on (F107), and it lays no sheet on the slopes it runs down (F109) nor lobes from each run (F110); there is no flood process on Io at all, its lava being the paterae's flow fields alone. One process on every body (principle 11) asks that Io's eruptions run through the floods' process, with Io's own rates and volumes.

Direction: Io's eruptions as floods from its paterae, by its heat budget's magma, filling and spilling over its plains and paterae as the floods do elsewhere, the flow fields' look (F59 to F62) carried over.

When: its own milestone, after this one.

Closed in part in the forty-ninth milestone (F117): Io's eruptions are floods from its paterae's rims, filling and spilling over its lows and the paterae they reach, laying their sheets and lobes as every body's floods do, the paterae cut through them after; and the form of every flood's flows follows how fast its lava poured out, compound fields where slow, sheet lobes where fast, so that one process gives Io's compound fields and the maria's and Tharsis's sheets alike. Still open: the eruptions' volumes and number are the old fields' shares of their reach, not Io's heat budget's magma (the direction's "by its heat budget's magma", with audit M8); the flood model's sheet is 50 m thick on every body, where a compound field's flow units are a few metres; and one eruption's lobes keep its one age.

Closed further by F201, the hundred-and-fifteenth milestone: the fields' number and times follow the heat budget's lava (M8), and each field's lobes keep the times they were laid through its growth. Still open: the 50 m sheet against Io's flow units of metres.

### M21. A lid's extension is not one budget

The lid's change of radius stretches it (F147), and its faults take that strain up, but Ganymede's lanes (F44, F61) and Europa's bands (F83) and cracks open on drivers of their own, the heat excess and the tide, so the extension they take up and the lid's faults' are counted apart: Ganymede's lanes cover most of its surface while its lid's 0.27 percent of extension opens some 800 graben on top, faint at its scale; Europa's opens three in its last 100 Myr beside its bands. One process asks that a lid's extension be spent once, on the structures its state favours: graben in a cold thick lid, lanes where it is warm and weak (Buck 1991's wide rifts), bands where a thin shell parts over its ocean.

Direction: the lanes and bands from the lid's strain with the tide's and the heat's, the graben taking what they do not.

When: with the lanes' revisit (the M72 design's fifth item).

### M22. Pluto's and Mimas's oceans in their histories

Pluto's ocean freezes through in its thermal history, so its lid stretches by 2.4 percent of its area and opens some 1,400 graben across all of it, and its deepest water freezes as ice II from 2.8 Gyr; Pluto's faults are regional, and its lack of compressional features argues that no ice II has formed and its ocean survives (Hammond, Barr and Parmentier 2016), insulated, perhaps, by clathrates the model lacks (Kamata et al. 2019). Mimas's young ocean melts its frozen-through interior from 4.46 Gyr (F138, F140), shortening its lid by 4.8 percent, and refreezes in its last 10 Myr, opening a dozen graben; Mimas shows no tectonics, which is how its ocean is known to be young. Both are the thermal histories', which the faults read as they stand: no body's faults are tuned.

Direction: the ice shell's insulation and its ocean's survival in the thermal history; the share of a freezing or melting ocean's change of volume its pressure takes up rather than the lid (Manga and Wang 2007), which F82 reads for the cracks.

When: with the thermal model's revisit.

### M23. Charon's chasmata lack their fabric (closed in part by F149)

Charon's chasmata came out as long troughs 30 to 66 km wide and some 4 km deep in belts, as Serenity's and Mandjet's are, but each a flat floor between two walls with landslide alcoves. Serenity's walls step down in fault terraces, its floor is broken into blocks, its footwalls are tilted back with raised flanks, and the belt about it is a run of tilted blocks 5 to 6 degrees with 4.5 to 5.5 km of relief (Beyer et al. 2017). The flexure of a plate as thick as the lid the faults cut, or its elastic part, gives a swell of about a kilometre over some 300 km beside each chasma, not Serenity's narrow raised rims, whose 2.5 km plate is thinner than Charon's lid was in its thermal history. And the strain is the same every way, so the belts gather where faults happened to fall, where Charon's runs round the boundary of its southern plain.

Closed in part by F149 (the seventy-second milestone, revised): a graben's walls are scarps in straight pieces with jogs, its floor the old ground let down, its footwalls raised, two in five tilted; the chasmata 24 km wide at the median. Closed further by F152 (the seventy-fourth): the floor of a wide chasma breaks into tilted blocks on a population of interior faults from the throw's partition, as dominoes above a listric master, which step the walls down in terraces and cut Valles Marineris's floors as well. Still open: the belt's own tilted-block fabric between the chasmata, and the belts' direction.

Direction: the fault fabric from the throw's partition on a listric master fault over the lid's ductile base, with rollover and tilted blocks, as Mars's canyons and Noctis Labyrinthus need; the relief's own stress, stretching the highs and compressing the lows (Bilotti and Suppe 1999), to give the strain a direction.

When: the next milestone, with Venus's rifts and Mars's canyons.

### M25. Venus's coronae on its rifts

Two thirds of Venus's coronae lie along its chasmata (Stofan et al. 1992). With the limbs drawing the diapirs (F151), two in five lie on the limbs and a fourteenth within 150 km of a rift. Two things hold the share down. The rifts are a few lines in a belt a rise wide, so a corona on the limb is seldom at a rift; and a diapir makes no corona where the mantle about it melts (F145), while the real coronae crowd the hot upwelling region, the BAT, where the model's mantle is hottest. Thinning the lid over the rifts, the measured state (Smrekar et al. 2023), made it worse under the melting rule, and drawing the diapirs toward the rift zones moved them into the hot ground about the rises, where the loop's count fell by three quarters.

Direction: the melting rule itself, which should give a corona where a rift focuses the upwelling's melt rather than none where the region melts; and the limbs' own heat in the mantle, which would place the diapirs' births along them without a draw.

When: with the corona process's next revision.

Closed in part in the ninety-sixth milestone (F183): where the region melts, a diapir on a limb, where the rifts run, makes a pocket as far as its head at least doubles the region's melt, so that 14 percent of Venus's 841 coronae lie within 150 km of a rift and 38 on a limb, where 7 and 30 did. Drawing the births by the convection's upwellings spread them over the pattern's broad highs and was set aside. Still open: the share on the rifts themselves, which the rifts the loop opens, drawn after the coronae, cannot yet focus; the hot upwelling region the real coronae crowd.

### M26. The field reads the ground before the plains and the lid of the present (closed in part by F157)

The lid's stress field (F150) reads its relief as the ground before the rises, low-passed over the flexural wavelength of the lid's present elastic thickness, and holds no load of the floods. The plains' ridges (F156) stand where the field squeezes them at their time, but the field does not know the plains themselves: the floods' load, which bends the lid under a filled basin and squeezes its interior while stretching its margin (the lunar maria's mascons, Freed et al. 2001), stands in as a drive past the strength, a stated share alike everywhere; and a thick present lid spreads the relief's force thinner and over longer lengths than the thinner lid the plains formed on, so the field's lows are the broad ones of the ground before the floods. Venus's ridges lie 87 percent below its mean height against the 93 percent of its ridged plains measured, the rest on plains the field's lows of the old ground do not hold.

Direction: the relief and the floods' load read through the loop at their times, on the lid's elastic thickness then (F144), the floods' load bending the lid as a load on the plate (F146's flexure, solved on the sphere), so that the maria's ridges and the arcuate graben at their margins come of their own load and the drive's stated share goes.

When: with the rise's flexure, the revisit plan's fifth item.

Closed in part in the seventy-sixth milestone (F157, F158): the loop reads the relief and the lid at epochs, and each load of lava bends the lid of its time as a shallow spherical shell, whose membrane turns a filled basin's thrusts round it and pulls the lid beyond its margin both ways, where its graben open; the ridges read the field alone and the stated share is gone. Still open: a floor remains, the plains standing at their strength, for the drives the field does not hold, Venus's climate above all, whose history the atmosphere held at the present does not carry; the relief's potential energy is still carried by the elastic lid alone; and the rises' own faults read the field of the ground before the loop and the lid of the present, as they are planned before it runs.

### M27. Mercury's thrusts have no direction (closed in part by F160)

Mercury's lid is squeezed alike every way by its shrinking, 313 MPa of elastic stress against a few of its relief's, so its thrusts lie and run at random (F155). Its large thrusts run north and south at low latitudes and somewhat east and west at high ones (Watters et al. 2015), a pattern of its global contraction with its tidal despinning, or with its mantle's convection; neither is in the model, its rotation held at the present (the inputs held constant) and its mantle raising no downwellings. Its thrusts also gather in systems and belts, as its graben do (F112), which the thrusts' population lacks: no thrust runs beside another or stops at an older one.

The plains' ridges lack their systems likewise (F159): the maria's few prominent ridges run in en-echelon systems and rings broken by gaps, where the model's, turned alike by the field, lie as a whorl of arcs filling a basin.

Direction: the despinning's stress from the rotation's history once that has one, and the mantle's downwellings in the field, squeezing the lid above them; the thrusts and the ridges drawn by the faults' shared generator (F147), with their neighbours and their arrest.

When: the despinning with the rotation's history; the generator with the thrusts' next revision.

Closed in part in the seventy-seventh milestone (F160): the thrusts and the plains' ridges gather in systems, half beside one drawn before near where they start and the rest stopping at one at least as long. Since the seventy-sixth (F157) Mercury's field no longer leaves its thrusts at random: the loads of its plains and its relief, carried by the thinner lid of their times, decide most of their strikes (the field's differential stress 9 MPa at the median against a scatter of 4.6), so that they run regionally alike, by 0.3 to 0.75 in sixty-degree bands of longitude where they ran at random before, along its plains' margins and their loads' fields. Still open: the north and south its large thrusts run at low latitudes, which wants the despinning or the mantle's convection.

Closed further in the ninety-sixth milestone (F182): the mantle's flow bends the lid, 7 to 15 MPa of differential stress on Mercury, and loosens its thrusts' trend. Their strikes are its plains' loads', some 70 to 100 MPa through their membrane at its equator, east and west at low latitudes, where its large thrusts run north and south: still open, the despinning, which wants the rotation's history (an input held at the present).

### M24. The lid's faults are not in the loop's history (closed by F153)

A rise's wide troughs are carried by the coarse loop as they grow (F105, F116), so its floods pond in them and its relaxation works on them; the lid's faults (F147) were drawn by the finish alone, among the craters at their times, as the narrow graben are. The loop's floods and relaxation did not see them: on Charon a later flood's lobes piled up at a chasma's end in jagged fronts instead of ponding along it, as Vulcan Planitia's plains embay its chasmata, and old chasmata on the warmer lids kept walls the loop would have relaxed.

Closed in the seventy-fourth milestone (F153): the strain process carries the lid's faults the loop's grid resolves, carving each within the step it forms, and the finish lays its own version over the loop's, as for a rise's canyons.

### M28. The shields' melt is not one budget (closed by F164)

A rise's shields take a share of the floods' eruptions as the dome's growth wanes (F108, F119), and that share leaves the lava's budget whatever the shields do with it. Their edifices grow by their own clock (F120), fastest as they start, when their rise's capture is still nothing, and the lava past their feet follows the edifices (F163). The two do not meet: Mars's shields take 6.6 × 10⁷ km³ of its floods' budget, while their edifices hold 8.3 × 10⁶ km³ and the lava they send past their feet 5.4 × 10⁶ km³; Venus's take 2.7 × 10⁵ km³, while their edifices hold 1.9 × 10⁶ km³ and their feet's lava 1.0 × 10⁶ km³, more than twice its floods' whole budget after the overturn, 1.2 × 10⁶ km³.

Direction: a shield's growth from the supply it takes, its edifice by what it holds and the rest past its foot, so that the melt is one budget and the edifices grow as the supply comes.

Closed in the eightieth milestone (F164): a shield grows from the lava its rise's eruptions bring it, scheduled before the loop, and the floods lose what the shields take; an edifice holds what its flows lay on it until it is full, its flows carrying their share past its foot, and the magma a full edifice cannot lift intrudes beneath it. Mars's edifices fill by 2.9 to 2.6 Ga and 5.7 × 10⁷ km³ intrudes; Venus's shields stand 6 to 15 percent built on its small budget (M29).

### M29. Venus's shields are as large as its supply lets them (closed by F165)

Since F164 a shield grows from the lava its rise brings it, and Venus's shields take 3.4 × 10⁵ km³ after the overturn: 6 to 15 percent of their edifices, 0.4 to 1.6 km high, where Venus's large volcanoes stand one to several kilometres over their plains. Its budget after the overturn is a stated share, 2 percent of the extrusive budget (EPISODIC_POST_EVENT_FRACTION, "the rest having gone into the event"), 0.0024 km³/yr, and its shields take their share only once its rises' growth turns (F119), the floods taking all of it until then. The crater record bounds the plains' share: only 2.5 percent of Venus's craters are embayed by lava (Strom, Schaber and Dawson 1994), and the model's Venus has 5.5 percent of its surface formed after the overturn, its craters with it; a larger budget spent as floods would pass it.

Direction: Venus's eruptions after the overturn from its mantle's history rather than a share, gathered at its rises' volcanoes as the Beta-Atla-Themis region's are, so that the plains stay within the crater record while the shields take most of the lava.

Closed in the eighty-first milestone (F165): the ambient mantle makes the share of the melt its own gate keeps and the plumes the rest under their heads, so that once Venus's ambient mantle stops melting after its overturn its lava comes out at its rises alone; the melt before that is the overturn's own resurfacing. A rise's shields share its lava by the dikes their weights gather, and each rise's leading volcano stands 2.8 to 5.1 km over a few hundred kilometres, as Maat Mons does. Still open: the event lays its resurfacing by its own fill (EPISODIC_FILL_FRACTION), not from the melt its surge brings.

Examined in the ninety-ninth milestone, nothing changed. The surge after Venus's last overturn makes 2.3 km of melt per unit area and erupts 46 m of it (FloodProcess::surge_lava): the eruptive share's contraction shutoff (VOLCANISM_SHUTOFF_STRAIN) reads the body's whole contraction since it formed, 1.7 × 10⁻³, leaving 12 percent of the extrusive share. Read since the lid last formed, as the thrusts' squeeze is (F185), the surge erupts some 300 m, and Venus's floods after its overturn some seven times what they are, 7.7 million km³ against 1.0, its surface younger than the overturn 4.0 percent against 3.3; but its leading volcanoes then fill their edifices, 10.5 km high, and shed their flanks, where Maat Mons stands some 5 km over its plains: the correction overshoots Venus's young volcanism, which points at another limit on it, and was set aside. Laid as the event's lava, 300 m floods Venus to some 200 m below its mean ground and leaves 54 percent of it dry, where its tesserae are 8 percent; three times as much would leave 22 percent dry and ten times 4, and the plains' thickness wants some 2 km. The geological record has intrusion dominating eruption in every setting, the mid-ocean ridges' thin lithosphere included (4:1 to 9:1), so no thin-lid rule lets the surge's melt out. If the overturn renews the surface by the old lid's foundering, the tesserae are the old lid that did not founder: crust buoyant against the lid's thermal excess (ρ_m − ρ_c against ρ α ΔT d / 2, 107 km and 936 K) passes at 9.9 km, 99 percent of the model's Venus, whose crust stands 13.5 km at the median, 24 at the ninetieth percentile and 44 at the ninety-ninth; a crust whose base is hot enough (some 1,220 K for diabase) to decouple from a delaminating mantle lid wants 56 km, none of it. Basalt sinking turns to eclogite and drags itself down, so buoyancy resists only the start of foundering. A crust thicker than 30 km covers 7 percent, near the tesserae's 8, but on no founded criterion.

Direction: the event from the surge's melt and the old lid's foundering together, the tesserae what the foundering leaves, on a rule for the crust's part in foundering that the sources found: its strength coupling it to the sinking lid, the eclogite transition's depth and kinetics in the sinking crust; and the limit on Venus's eruptions after its overturn, which the contraction since the formation of the body stands in for.

### M30. A shield's flanks do not spread into its moat (closed by F172)

Olympus Mons is ringed by a basal scarp up to some 8 km high and by its aureole, lobes of ridged deposits that reach hundreds of kilometres out, which its flanks shed by gravity spreading or landsliding (Lycus Sulci among them), and its moat is thought mostly filled with lava and those deposits (Isherwood et al. 2013). Since F166 the model's shields stand in moats some 3 to 4 km deep on Mars, which the lava their flows send past their feet (F163) only partly fills; their flanks stand as they built them, without a scarp at their foot or the debris that would fill the moat.

Direction: an edifice's flanks fail when the stress in them passes what they hold, by spreading on a weak base or by sliding, as the canyon's walls shed their alcoves (F106's landslides), the debris running out into the moat and the lobes it lays ridged; a basal scarp where the flank's toe is cut back.

Closed in the eighty-eighth milestone (F172): an edifice past the relief its crust holds spreads its flanks with the lava it would have stacked above, its toe breaking away in a scarp and running out with that lava as its aureole's ridged lobes, as far as rock avalanches of their volumes run. Mars's largest shield ends in a scarp 4.5 to 9.9 km high and lays lobes 120 to 350 km long into its moat and beyond. The failure comes once an edifice has filled to the cap, so Venus's shields, which its supply builds only part of the way, keep their flanks. Still open: the lobes' directions from the ground they run down, and the weak base that lets a flank spread, which the model reads only as the crust's strength.

### M31. Venus's highlands stand on the primordial relief (closed in part by F173)

Venus's tessera plateaus stand one to four kilometres above its plains on crust their tectonics thickened, which holds them up as a compensated load: their gravity shows crust as thick as 40 km beneath them against a global mean under 30 km (Simons, Solomon and Hager 1997; Hansen, Banks and Ghent 1999; Maia and Wieczorek 2022). The model's tesserae are the old surface its overturn leaves, standing on the primordial relief (F137) and the craters on it. A hot early crust relaxes that relief away by its lower crust's flow (F170), and Venus came out flat but for its volcanoes, so an episodic lid does not relax.

Direction: the tesserae from the tectonics that thickened their crust, the folds and the shortening of the event or of its downwellings, compensated by their roots, so that the lower crust's flow relaxes them on the time its temperature sets, and Venus's relief relaxes as a stagnant lid's does. The rest of Venus's relief, some nine tenths of its surface, is held up by the mantle's flow (Simons, Solomon and Hager 1997): its rises over the upwellings, its lowlands, Atalanta and Lavinia, over the downwellings; the primordial relief stands in for that support (the user noted it is about right against Venus's map), and the lower crust's flow would flatten it. Its dynamic topography from the mantle's convective planform, which the thermal history does not yet hold, the rises' domes its start (F170 relaxes toward them), before Venus relaxes.

Closed in part in the eighty-ninth milestone (F173): a rocky mantle convects beneath its lid in cells from its depth, and the ground its flow holds up rises over the upwellings and sinks over the downwellings, some 300 m at its root mean square on Venus, from −1,350 to 890 m; the plumes rise at its strongest upwellings, the lower crust relaxes toward it and the lid bears none of it. Still open: Venus relaxing as a stagnant lid does, with the mantle's support as its equilibrium; the tesserae from the tectonics that thickened their crust, which the convection's downwellings would squeeze (the downwelling plateaus of Bindschadler, Schubert and Kaula 1992); and the stresses the flow puts on the lid, for Mercury's thrusts (M27) and Venus's ridges, and its upwellings' heat for the coronae (M25). The upwelling limbs between plumes (F150) follow the pattern's ridges.

Closed further in the ninety-seventh milestone (F182, F184): the flow's stresses bend the lid (F182, the ninety-sixth); the mantle's flow carries the crust beneath the stiff lid toward its downwellings after each overturn, piling it to the depth its base turns to eclogite, and the lower crust's flow spreads it, each cell under its own crust; the tesserae are the plateaus the last epoch piled, some 3.7 km above the plains, and an episodic lid relaxes as a stagnant one does, its primordial relief gone in its first gigayear; the crust that founders past the cap leaves Venus's 16 km thick on average, near the 20 measured. Still open: the tesserae's own fabric from the shortening that piled them (W3, T1); the share of the surface they take, TESSERA_FRACTION, which the surge's lava against the plateaus' height would set (M29); the plateaus' margins, which the hot phase smooths where Ovda's stand steep.

Closed further in the hundredth milestone (F186): the plains onlap the tesserae, the fabric sinking beneath them beyond the edge in an onlap wedge and the lava ponded within it at the edge's level, so that the contact interfingers, ridges standing as inliers and valleys embayed, each cell by the share the lava left standing; the tesserae keep their share. Still open: the tesserae's fabric from the shortening that piled them and the extension after (W3, T1).

Closed further in the hundred-and-first milestone (F187): the tesserae's fabric is the strain the mantle's flow put on their lid where its drag broke it after the overturn before the last, folds and graben at the brittle layer's and the lid's wavelengths. Still open: the strain after the event, which both the tesserae and the plains take and the thrusts and ridges stand for (F185); the lid's own thickening where it is shortened, which the crust's thickness does not count; the shortening the folds cannot take, which thrusts and the lid's thickening would.

Closed further in the hundred-and-second milestone (F188): the plateaus' collapse, the strain the lower crust's flow put on the lid as the relief relaxed over the lid's life, is a later set of graben and folds across the drag's. Still open: the strain after the event; the lid's thickening where it is shortened; the plains' strain from the same collapse, which the event's lava covers but the later relaxation would deform.

Closed further in the hundred-and-third milestone (F189): each layer's fabric, its sets together, stands no higher than the layer's strength holds. Still open: the plateaus' broad internal relief, the basins and swirls 50 to 300 km across in Ovda's interior, which neither the crust the flow piles nor the fabric makes.

Closed further in the hundred-and-fourth milestone (F190): each fold crest and graben is a line along the strain's strike, evenly spaced streamlines that end where the strain's sense changes or their neighbours crowd them, so that crossing sets bound cells. Still open: the collapse's graben stand at most 0.27 km, the lid's tension at the fault layer's depth, against the folds' 0.88, so the tiles over the plateaus' crests show faintly; the tension integrated through the stretched layer would deepen them.

Closed further in the hundred-and-fifth milestone (F191): the tesserae's graben follow the model's graben and the lid's necking its rifts, deepening them from 0.27 km to the graben's and rifts' own bounds, so that the collapse's crossing rifts tile the plateaus' crests; the yield envelope stays set aside (F178). Still open: the plateaus' broad internal relief; the strain after the event.

Closed further in the hundred-and-sixth milestone (F192): the tesserae's graben take the model's graben form, segmented, stepping, scattered in throw and spacing and wandering, a later fault ending at an earlier, the lid rifting across its greatest stretching alone, so that the plateaus' crests no longer carry a grid. Still open: the plateaus' broad internal relief; the strain after the event.

### M32. Valles Marineris is a stretched trough, not a subsided one

Valles Marineris stands some 7 to 10 km below the plateau about it. Its faults dip at least 85°, taking up less than 4.6 km of extension, so that it is a trough the ground subsided into rather than one stretching opened (Andrews-Hanna 2012, part 1); the tension that opened it lay in a narrow belt just south of and along the dichotomy's boundary buried beneath Tharsis, where the dikes it drew cut the whole lithosphere (part 2); and the blocks those dikes cut free sank to their isostatic level, some 0.5 km, and then under the sediments that filled them and by the crust's slow flow some 8.6 ± 5.6 km in all (part 3), the fill since partly eroded. The model's canyons are a rise's normal faults, as deep as a thousandth and a half of their length, which open radially from the rise wherever its tension passes the lid's strength: on Mars's preset (seed 1048) its one carried canyon runs 1,185 km from the rise into the lowland, across the boundary, some 2.8 km deep. The dichotomy's own tension, the spreading of its highlands down their slope (F150), lies along the boundary as part 2's belt does but is some 2 MPa against the rise's 100 to 250 MPa.

Direction: the rise's stresses solved through a lid whose stiffness follows the crust's thickness, so that the boundary focuses them (part 2); the troughs where the dikes cut the lid, subsiding as their blocks lose its support and sink under their fill (part 3), with the sediments and their erosion of Phase 5; the fabric of their walls with M23. Deferred by the user to be done with Phase 5.

### M33. The Moon's maria do not gather in its province (a known limitation)

The Procellarum KREEP Terrane covers some 16 percent of the Moon and holds more than 60 percent of its maria by area (Jolliff et al. 2000). The model's Moon (seed 10) lays lava over 30 percent of its surface (its maria cover some 17) with 20 to 28 percent of it within its province's radius, which covers 12 percent, at any weight the province is given (F171): its basins' thin crust decides where its lava lies.

Tried, not kept (2026-10-28, measured on the Moon's preset at 2,048 pixels):

- **The province as a mantle of its own.** The magma ocean's last liquid, a share of the crust's heat-producing elements, heated a second mantle history beneath a Gaussian province. Its mantle ran only 50 to 60 K hotter, convection being a thermostat, and melted 1.5 to 3 times what the rest did; 16 percent of the lava lay within the province's radius.
- **The heat taken from the cumulates, the liquid gathered.** The last liquid took its heat from the mantle's share rather than the crust's, so that the mantle elsewhere, depleted, stopped melting at about 1.5 Gyr rather than 3.4, as Laneuville et al. (2013) find for the Moon. After the first step the liquid gathered beneath the crust the mantle's flow had thinned (M35's degree one), filling the crust's base as an inverted lake up to its share of the surface: an irregular province like the terrane, and every eruption after 1.5 Gyr within it. But the province's melt, seven parts in eight of it intruding (F174), laid 25 to 35 km of rock beneath the province and raised it 11 km above the rest, even with its loads laid on its own thinner lid. Its lava ran off the plateau, and 17 to 19 percent of the lava lay within it. With the vents' paths read over a reach about the melt's depth instead of the whole body, every region erupted its own melt into its own lows and the maria covered the far side.
- **The liquid's heat in the lid** (Wieczorek and Phillips 2000). The layer's heat bent the lid's conductive profile and thinned it, the mantle beneath no hotter. The volumes were the Moon's (8.6 × 10⁶ km³, lava over 18.6 percent of the surface) and the province melted a billion years longer than the rest, but no more, so 13 percent of the lava lay within it, and the province still stood 3.7 km high. The general form had large side effects. Every crust's heat thinning its lid, and a lid no thicker than its crust not foundering (the province's thin lid otherwise foundered 219 times and froze its mantle), doubled Venus's and Earth's overturns and stopped Io's.
- **The gathering alone,** the province still a weight on the vents: 7.5 percent.

What blocks it:

1. **The hot start.** For its first 1.3 to 1.5 Gyr the mantle melts everywhere, making 40 to 90 percent of all the lava, which floods the far side and its megabasin. The magma ocean whose freezing would set the early mantle is not modelled (M35).
2. **The intrusive share.** A province hot enough to make Procellarum's lava intrudes tens of kilometres of rock beneath itself, as a plateau its lid holds up; the terrane is low ground.
3. **The megabasin's crust.** The megabasin thins the crust to some 5 km, against some 20 km beneath the South Pole–Aitken basin (Wieczorek et al. 2013), and draws every stream's lava.

Direction: the magma ocean's freezing first (M35), setting the early mantle cool and stratified, with the fusible late cumulates gathered with the last liquid; then the province's heat shared between its lid and its mantle as the layer's profile shares it, the intrusions' share and weight where the magma is denser than the crust, and the megabasin's crust against the measured. The tries' code is kept outside the repository; its pieces worth keeping are the gathering beneath the thinned crust, the depleted mantle elsewhere, and the loads laid on the lid they lie on.

A known limitation (the user, 2026-10-29): every physical route tried here and under M35 has failed, so the Moon's maria lying on both hemispheres is recorded as a limitation of the model and not pursued further; it is not a major concern for the Moon's look. The nearside highlands' darkening toward the maria (M36) goes with it.

### M35. The starting state: a finished crust on a mantle with no lid

A rocky body starts with its crust whole, tectonics.initial_crust_thickness (40 km on every rocky body, S2), on a mantle at its starting temperature with no lid, which then grows as a stagnant lid does (F121). For its first 30 Myr the stiff lid is thinner than the crust, so the mantle's flow carries the crust (F184) and piles it over the downwellings by ±3 km of relief at the cells' wavelengths, which a cold crust holds for the age of the solar system: on Mars at degree 2.4, on the Moon at degree 1, on Mercury at degree 9. A magma ocean's crust and its lid form together as it solidifies, over a time the model does not hold, and the Moon's and Mars's degree-one crusts have been laid to such a flow; the amount is uncalibrated.

Direction: the magma ocean's solidification, which sets the crust's thickness by class (S2) and the first lid together; until then the starting crust and lid are inputs. M33's tries (2026-10-28) found the hot start the main thing keeping the Moon's maria from gathering in its province: the mantle melts everywhere for its first 1.3 to 1.5 Gyr, making 40 to 90 percent of the lava, where a frozen magma ocean leaves cool, stratified cumulates, its fusible late ones with its last liquid.

When: with the crust by class (S2), the earlier-phases revisit.

Tried on the Moon's preset (2026-10-29), none kept; the code is kept outside the repository. Together with M33's tries, these are the measurements the next attempt starts from.

- **The ocean's freezing.** The magma ocean that floats plagioclase freezes through its crust as a Stefan front, its residual liquid's latent heat leaving by conduction: some 30 Myr for the Moon, by the share of the ocean that floats plagioclase, and nothing melts beneath it meanwhile. The Moon's crust is already lopsided without it: the mantle's early flow (F184) leaves its hemispheres' mean crusts 14.8 km apart (the Moon's some 20 km). Letting the crust float free of any lid while the ocean froze doubled that to 32 km and sank 8 km of crust into the mantle.
- **The last liquid gathered.** At the ocean's end the last liquid gathers where the crust's base stands highest, averaged over the province's breadth, and lands on the thin side. But the lava within the province's radius fell from 15 to 8 percent of the whole: the vents still follow the thinnest crust anywhere, the basins', and the province adds no heat.
- **A cool start.** The ocean that floats a crust leaves its cumulates at the solidus at the crust's base, some 1,420 K on the Moon, where one that floats nothing leaves a hot mush (1,750 K); by the share that floats plagioclase, so only the Moon changes. The mantle beyond the province then never melts, and the lava falls to 9 percent of the surface.
- **The province's own mantle.** The last liquid takes nine tenths of the mantle's heat-producing elements beneath the province (some 11 percent of the surface): its mantle melts from 0.2 Gyr to the present, peaking at 0.9, and the depleted mantle elsewhere not at all, their melts taken into the vents' one weight. The lava covers 14.5 percent of the surface, 94 percent of it on the province's hemisphere (the nearside holds some 94 percent of the maria; before, 64). But some 16 km of the melt stays beneath the province's centre, and the lid, which deflects some 0.4 km under a 15 MPa load, holds it up as a plateau 15.7 km high; the lava fills a ring of basins 40 to 80 degrees out.
- **The liquid's heat in its layer at the crust's base.** In the lid once the lid grows past it, the heat thins the lid and leaves through the surface. The province melts from 0.2 to 2.1 Gyr, 99.7 percent of the lava on its hemisphere, but the lava covers 4.6 percent: in the stagnant lid's scaling the mantle's flux into the lid does not depend on the lid, so heat within the lid never warms the mantle beneath it. The province's own history must founder only with the body's lid: its lid, thinner than its crust, broke again and again on its own and cooled its mantle to 1,000 K, and a general rule that a lid of crust cannot founder would stop Io's overturns.
- **The giant basin.** On the preset's seed the 2,576 km basin strikes at 0.415 Gyr, 71 degrees from the province gathered at 0.05 Gyr, and its ejecta raises the province's centre 8 km whichever way the heat lies; the South Pole–Aitken basin lies some 120 degrees from the Procellarum KREEP Terrane.

Direction: the province's loads laid on its own hotter, thinner lid, so that the intrusions beneath it sink as the terrane's low ground has them, or the lid's heat warming the mantle beneath it (a stagnant lid whose flux reads the heat within it). The cool start, the depleted mantle elsewhere, the regional histories foundering with the body's, and the melts taken into one vent weight are the pieces worth keeping.

The Moon's part closed with M33 as a known limitation (2026-10-29); the starting state stays an input until the crust by class (S2).

## Validation (principle 15)

### M34. A rise is a dome, not a pile (closed by F174)

Tharsis's long-wavelength gravity and Mars's shape are explained by its load of some 3 × 10⁸ km³ of magma (Phillips et al. 2001), most of it in place before about 3.7 Ga (Hynek et al. 2011), which buried the ancient terrain beneath it. The model's rise is a dome its plume holds up, 0.28 of the relief the crust holds and no volume of its own; its lava comes only through the floods, spread over the body by their vents' weights (F171), and the shields, and the magma the full edifices cannot take, some 5.7 × 10⁷ km³ on Mars (F164), is lost rather than thickening the crust beneath them. On Mars's preset (seed 1048, F173) the rise stands over an old basin some 2,400 km across, whose floor it lifts some 2 km while the basin keeps its depth and its shields' moats sink in it, where Tharsis's lavas would have filled it first.

Direction: the plume's melt under its head builds the rise: the intrusions it makes, the full edifices' lost magma among them, thicken the crust beneath it and raise it as Airy's isostasy does, and its lavas pond in the rise's lows, the floods' budget kept as one, so that a rise is a pile as well as a dome.

When: the next milestone (the user's choice, 2026-10-22).

Closed in the ninetieth milestone (F174): the melt is one budget, what stays below intruding the crust beneath where it rose from and the full edifices' surplus beneath their rises, the ground rising as the crust thickens and the intrusions' weight bending the lid. 3.2 × 10⁸ km³ intrudes beneath Mars's rise, its flanks some 2 km higher and the old basin raised some 2.5 km. The lavas ponding in a rise's own lows followed in the ninety-second milestone (F176): the vents read the ground less what the plumes and the mantle's flow hold up, and the rise draws some two fifths more of the floods, the basin beneath it taking a kilometre more of lava, which the lid's sagging under a load so wide leaves raising its floor some 0.2 km, as Mars's filled basins still stand deep. Still open: the dome's height, still the share Tharsis's whole height was fitted to (RISE_HEIGHT_FRACTION), which its plume head's buoyancy should give now that the pile stands beside it.

### M36. The Moon's highlands stand too dark against its rays

The user's note at the Phase 4F gate check (2026-10-26): the contrast between the Moon's highlands and its rayed ejecta is too strong; the reference map shows the highlands a lighter shade. The bright rims of its fresh craters, which the reference shows too, stay.

Direction: the highlands' mature albedo and the rays' fresh one against the measured, from the maturing of the regolith rather than a tone.

When: with the Moon's gate items.

The rays' half done (F212, 2026-10-28): what a ray's secondaries throw out far from its crater is mostly the skin's mature soil, so the rays fade from their crater to the ground, the brightest half percent of the Moon at 0.228 against the reference's 0.195 (0.254 before) and its 95th percentile at 0.153 against 0.163. The highlands' own albedo, the median 0.105 against the reference's 0.116, waits for the silicates' space-weathered colours. The rays still read as smooth spokes from a bright disc near their craters, where Tycho's are fans of fine streaks: where the jets excavate more than RAY_CLUSTER_CAP the ray is drawn smooth (F128), and the clusters that would streak it lie farther out, now faint.

The highlands' colours done (F213, 2026-10-28): the mature colours are weathered by each material's iron, the farside highlands' pure mature ground at 0.149 against the reference's 0.13 to 0.14, South Pole-Aitken's noritic floor at 0.080 against 0.080 and the maria at 0.052 against 0.053. The brightest half percent, the rays' near parts, stands at 0.197 against the reference's 0.195 (0.228 before). The highlands now read darker as a whole for another reason, the median 0.089 against the reference's 0.109: the model carries more of the intermediate shades, 36 percent of its surface between 0.065 and 0.10 against the reference's 24 and 17 percent between 0.10 and 0.13 against 36, the noritic deep material its basins excavate from below 25 km (materials.deep_bedrock_depth, provisional) lying in wide dark halos about them and its megabasin, and the lava mixed sideways about its maria. On the real Moon the basins' ejecta are mostly feldspathic and South Pole-Aitken's floor is the noritic province.

Direction: the deep material's depth from the crust's structure (the upper crust's feldspathic thickness over a noritic lower crust) and the basins' ejecta from the depths their excavation reaches, with M33's megabasin.

The halos done (F215, 2026-10-28): a crater brings up the deep bedrock by the share of its excavation cavity below it, at a depth that follows the crust under it, so the basins' ejecta are mostly upper crust, Orientale's size a quarter norite, and South Pole-Aitken's province a partial share. The farside highlands' median is 0.133 against the reference's 0.125 (0.114 before). The nearside highlands now stand too bright, 0.150 against 0.114, and the median at 0.130 against 0.109: the model's maria are not gathered on one hemisphere (M33), and its highlands take none of the basalt that impacts in the maria throw onto them, which darkens the real nearside's highlands and its intermediate shades (24 percent of the reference's surface at 0.065 to 0.10, 6 percent of the model's).

Direction: the maria's ejecta spread over the highlands about them, with M33's maria; the upper crust's own variation in iron (the outer feldspathic highlands' FeO above the core's).

The rings and floors done (F216, 2026-10-28): the deep share lands nearest the rim and floors the crater, so the basins wear collars of norite and the larger craters noritic floors; the median is 0.114 against the reference's 0.109 and the farside highlands' 0.127 against 0.125. The user's note (2026-10-28): the darker material about the maria, grading into the highlands, gave a natural transition that the real Moon's map shows and that F215 took away; the nearside highlands still stand at 0.147 against 0.114, the shades at 0.10 to 0.13 over 16 percent of the surface against 36. That transition is the maria's own basalt thrown onto the highlands by the impacts in them, the direction above.

Tried, not kept (2026-10-29): the basalt sprayed by the craters smaller than the finish draws. The impacts in the maria throw basalt onto the highlands by the production of craters too small to draw, each spreading its ejecta as the cube of its range: some 0.003 on the highlands near the contacts, the production function run past its measured sizes to reach even that.

The buried maria done (F217, 2026-10-29): lava a crater's ejecta covers is buried, and later craters dig it up by the share of their cavities below it, the drawn ones as dark-haloed craters over 1.55 percent of the Moon's surface about its maria and the smaller ones stirring up 0.07 percent; the burials lie deep, half below 1.1 km. The highlands' profile toward the maria moves by 0.003 (0.133, 0.138 and 0.138 at 25 to 200 km against the reference's 0.102, 0.112 and 0.120). Neither the spray nor the buried maria makes the gradient: cryptomaria are hidden by definition, and the darkening most likely is the nearside's crust itself, richer in iron (Jolliff et al. 2000) about the maria gathered there, so that the distance from a mare largely measures the nearside against the farside. That waits on M33's province and M35's start.

The gate's re-check (the user, 2026-10-29): the rays still stand brighter against the ground than the reference's; the maria's interiors show flows of different ages with mottled contacts, as they should. The highlands' even tone lacks the reference's mottle, which does not follow their height on the real maps: it is to come from noise modulated by a physical field (principle 16; M41), later. The nearside's darkening toward the maria is closed with M33.

The rays' contrast done in part (F220, F221, 2026-10-30): the regression the user saw came with F212, which freshened nothing where a jet excavated under half the ground and so took away the clusters' fine streaks; a secondary's depth now goes with its distance and its jet's speed, so the sparse clusters freshen their ground again. Measured in crater radii against the ground 10.6 to 13 radii out, ten young model craters matched eight young reference craters beyond 2.8 radii but stood at 2.08 from 1.4 to 2.8 radii against 1.10 to 1.66: the rays now freshen only the ground the continuous deposit leaves, and the deposit by its share of the crater's own ejecta, so the ring stands at 1.70 and the brightest half percent of the surface at 0.202 against 0.192. Still open: the near ring against the reference's median of 1.37 (Tycho's 1.22 is in part its dark halo of impact melt, which the model lacks); and most of the maria's bright excess, 0.84 percent of the dark ground 50 percent brighter than the ground about it against 0.23, which is the rims of craters standing in and about the lava, not the rays.

Tried, not kept (2026-10-30; F222, the hundred-and-thirty-third milestone, reverted at the user's verdict): a crater in lava kept the lava units but for the share of its cavity in the deep bedrock, where a crater that reached the deep bedrock by any share takes the deep units and loses its lava. The maria's ground standing more than 1.5 times the ground about it fell from 1.29 to 1.02 percent against the reference's 0.30, and the craters younger than the lava with bright rims from 10 to 4 percent; but the coarse loop's history shifted with it, moving the maria's outlines in places, and on Mercury the lava, brighter than its crust, showed in the craters on its plains, its mean albedo up 2 percent. The user found the terrain on both bodies more cohesive without it.

### M37. Io's layered plains

Some 4 percent of Io's surface is layered plains, layers of the plains isolated from those beneath them by bounding scarps, irregular and fretted (Williams et al. 2011; McCauley et al. 1979), which sapping by liquid SO2 erodes back, with lobed deposits at their feet (Moore et al. 2001); Io's colour mosaics show them as plateaus with bright rims. The model's plains are one surface, buried smooth (F86), its plateaus only the mountains' layered form (F87).

Direction: the plains' layers stripped back from their scarps by sapping, the SO2 within them escaping where a scarp exposes it and the debris spreading in lobes at its foot, from the scarps that the paterae's walls, the mountains' feet and the flows' fronts start, with frost at the scarps' feet; the scarps' heights and their rates of retreat from Io's imagery.

When: with Phase 5's erosion, its mass wasting and sublimation textures (decided 2026-10-27, against building it with every number stated).

### M38. Mercury's rays

Mercury's reference map shows its young craters' rays across the planet, Hokusai's and Debussy's the longest, and bright halos about many craters. The model's Mercury has all but none: its rays fade in 0.15 Gyr (RAY_FADING_GYR_1AU over its solar flux), and after F212 0.1 percent of its surface is brightened past a tenth of the fresh contrast. Braden and Robinson (2013) compared the maturing of Mercury's and the Moon's regolith; its darker, iron-poor surface may also lower the contrast a fresh exposure has.

Direction: the rays' fading on Mercury from its regolith's maturing (its micrometeorite and solar-wind fluxes against the gardening that brings fresh soil up) and its fresh contrast, with its albedo (the colours step: the model's 0.083 against 0.142).

When: with the colours step of the Moon's gate items.

The fresh contrast done (F213, 2026-10-28): Mercury's regolith weathers 4.4 times as intensely as the Moon's, so its crust matures from 0.17 to 0.077 and a fresh exposure stands 2.2 times as bright as the ground, as on the Moon; its mean albedo is 0.073, the measured comparison on one photometric standard (Braden et al. 2011) putting its mature ground at 0.55 of the lunar highlands', where the geometric albedo of 0.142 holds its opposition surge. Its rays still fade in 0.15 Gyr, 0.1 percent of its surface brightened past a tenth of the fresh contrast, where Braden and Robinson (2013) find its optical maturation up to four times as fast as the Moon's, not the 6.7 the solar flux gives.

The fading done (F214, 2026-10-28): a fresh surface matures by its weathering, the solar wind and micrometeorite melt, and by its gardening, Mercury's 2.4 times as fast as the Moon's where the sunlight alone gave 6.7; its rays last 0.42 Gyr and its ground matures in 0.12 Gyr, so it carries 0.58 of the Moon's young craters per area against the measured half (Braden et al. 2012; Chabot et al. 2009), 20 rayed craters at a 2048-pixel map where it had 5, its colours unchanged. Two things still keep its map from the reference's look. Its bright young craters are mostly smaller than the cells of a map at that width (7.5 km), so they show only at higher resolution. And the dark halos its craters dig from its graphite-bearing deep crust stand out as round blots with ragged edges across its map, where on the reference the low-reflectance material is diffuse, about some basins and craters, and the bright craters and rays catch the eye.

Direction: the deep crust's halos on Mercury: its depth and share against the measured low-reflectance material's extent (Denevi et al. 2009), and their mixing into the ground by later ejecta and gardening, as the Moon's noritic halos (M36).

The blots done (F215, 2026-10-28): a crater brings up the deep crust by the share of its cavity below it, so Mercury's halos are partial and fade into its ground, the surface at least half deep material 14 percent against Denevi et al.'s (2009) 15 (51 before), and its rays and bright craters catch the eye. Its bright young craters smaller than the map's cells stay open.

The user's note (2026-10-29): Mercury's rays look about right; what it lacks are the small bright-walled craters the Moon's map shows. Three things keep them from it. Its complex craters, shallower than the Moon's (Pike 1988), have walls near 19 degrees at 30 km, under the regolith's stable 20, where the Moon's stand near 25. Its finish draws craters from 7.5 km at the standard resolution, past its transition at 7.2 km, so it draws no simple craters, whose walls shed in full, where the Moon draws them from 5.3 to 16.5 km; the craters smaller than the finish draws shed nothing on any body. And its measured transition is 10 to 12 km (Barnouin et al. 2012), some 1.8 times Mars's at the same gravity, a difference neither the targets' properties nor the projectiles' speeds account for, where the model's transition, by gravity alone, gives both 7.2 km. Both reference maps also show their craters' walls partly by their shading. The terraces done (F218, 2026-10-29): a complex crater's wall sheds as read along its terraces, its scarps some 3.6 times as steep as its mean slope, and Mercury's steep share rises from 0.74 to 1.11 percent, its complex craters faintly bright-walled. Still open: Mercury's transition against the measured, and the walls of the craters smaller than the finish draws.

Mercury's transition tried (2026-10-29), not adopted. Susorney et al. (2016) put it at 11.7 ± 1.2 km from MESSENGER's altimetry and imaging, where Talpe (2012) found 8 km from the depth-to-diameter break alone. The one mechanism with numbers is the impact speed. Silber et al. (2017) model lunar craters whose walls collapse while the rock is fluidized, which lasts as long as blocks of the impactor's size take to ring down. A small, fast impactor's collapse stops sooner, so the onset of complex craters grows as the speed to the 0.47: 19.1 km at 10 km/s and 26.4 km at 20 km/s, their 6 and 15 km/s runs bracketing the same power. Read from the Moon at its 15.2 km/s, with the ice factor refitted to 0.14 to keep Ganymede and Callisto at their 3 km (Schenk et al. 2004), it gives:

- Mercury 11.2 km, its 40 km craters deepened from 1.8 to 2.4 km against Pike's (1988) 2.2;
- Mars 6.2 km, against about 6 to 8 measured;
- Ceres 9.3 km (14.1 now), against 7.5 to 12 (Hiesinger et al. 2016);
- Europa 3.6 km (3.1 now), against about 5;
- Pluto and Charon, struck at about 2 km/s, 2.5 and 5.0 km (6.5 and 13.9 now), against 11 to 12.5 and 13.5 to 16 (Robbins et al. 2021).

The icy bodies' transitions show no sign of the speed: Charon's, scaled by gravity, is Ganymede's at a tenth of the speed. A second law to save the rule, colder ice standing stronger, has no measurement to set it by, and the targets' strength (Pike 1980, 1988; Barnouin et al. 2012) has none for Mercury against Mars. The user chose to record it and keep the transition by gravity alone. Still open: the walls of the craters smaller than the finish draws.

The fine walls done (F219, 2026-10-29): the craters too small to draw shed on 0.56 percent of Mercury's surface and the Moon's, brightening the mean albedo by about 0.1 percent, too little to see, and Mercury's walls, drawn and not, cover 1.64 percent at the standard resolution and 1.55 at the high. What Mercury's map lacks against the Moon's is drawn simple craters, whose walls shed in full: its transition by gravity alone, 7.2 km, keeps them under the 7.5 km the finish draws at the standard resolution, where the Moon's run from 5.3 to 16.5 km. That stays with the transition, tried above and not adopted.

The detail back (F216, 2026-10-28): the user saw Mercury lose its contrast and detail under F215, its map 1.26 from its darkest tenth to its brightest against 1.46 on the user's true-colour map, which shows it brown-grey as the Moon and of low contrast. With the deep share landing nearest the rim and flooring the larger craters its contrast is 1.37, 1.23 at some 90 km against 1.33, its larger craters dark-floored and dark-rimmed.

### M39. Enceladus's fractures do not form a global network (a known limitation)

The gate's re-check (2026-10-29), the user's verdict: Enceladus is near passing. Its terrains blend into one seamless look (their ridged margins were taken out because they drew the terrains' outlines), and each terrain's fractures and grooves look right on their own. But they stay within their own terrain: on the reference map they cross from one terrain into the next and join into a network over the whole body. That alone would make Enceladus pass.

Direction: the fractures and grooves read from one global stress and strain history, so that a set that forms in one terrain carries on into the next where the stress allows, thinning out over ground that resists it, and crosses or meets the sets of other episodes.

When: before Phase 4F closes.

Tried (2026-10-29), not kept; a known limitation (the user's choice beforehand: try it, and mark it a known limitation if it does not work). A zone's fissure stopped where its zone's taper fell to 0.3, and every crack ran only while the stress across it passed the 40 kPa that starts one, where Enceladus's tide stresses its thick shell by some 20 kPa. A crack that has started runs on while the stress intensity at its tip passes the ice's fracture toughness, some 0.1 MPa m^0.5 (Schulson and Duval 2009): an edge crack opening down against the overburden feels at its tip K = sqrt(pi d) (1.1215 s - 0.683 rho g d) (Tada, Paris and Irwin 2000), at most 0.98 s^1.5 / sqrt(rho g), so a running crack needs some 10 kPa on Enceladus and 23 on Europa. With that, the zones' fissures ran on beyond them as fractures of the tide, each relieving the membrane stress for half its length as F200's do: 153 of them, 20,300 km in all, about 130 km each, linking the terrains in the data. But it changed nothing visible:

- 143 of the 153 are dated to past episodes, and their troughs, 100 to 350 m deep, lie under the plume's fallout laid since (334 m at the median, F74), as the order of events has it;
- the present episode's 10 do not stand out among the terrains' own lineations;
- the reorientations' fractures, run on to half the body round, relieved the shell in three (some 20 to a reorientation before);
- Enceladus's finish took some 400 s more (931 s against 520).

The network the reference shows wants troughs whose relief outlasts the fallout, kilometres wide and hundreds of metres deep, and so a stretch of the shell of about a percent, where the model's stresses stretch it by some 0.04 percent. Kept outside the repository with the measurements, for when the physics is final and principle 16 decides what noise may carry.

### M40. Io's flow fields do not merge

The gate's re-check (2026-10-29), the user's verdict: Io passes. The flows' shapes are among the model's strongest results and are to stay. But each lobe stands apart from the next, where on the reference the flows of a field join into broad regions of graded tone. Pele's red rings are a little fainter than the reference's, and clear features of the albedo; strengthening them is optional.

Looked into (2026-10-30), not changed; measured at 2,048 pixels within 45 degrees of the equator. The fields' placement is not what differs: the fresh-flow unit covers 27.4 percent of the surface (Williams et al. 2011: 27.8) in broad, connected sheets with lobate margins, as the reference's fields are, and the ground darker than 0.9 of the median within 300 km of it covers 16.7 percent against the reference's 16.4, in regions as large or larger. What differs is how the sheets show. A third of the flow cells show under 5 percent of their lava, buried to the plains' colour: not by the far fallout, which would leave 99.6 percent of them showing more than 60 percent, but by the plumes at the young flows' fronts, whose fall over the history's last 500 years, 0.2 to 0.9 optical depths a year near the flows, coats the older lava over its whole age, 1,000 to 100,000 years, against the rough lava's 10,000 optical depths (and stands where today's fronts are, where the fronts of the older lava's own time stood elsewhere on its patera's fields). Each lobe, laid at one time, is dark or gone, with a narrow band of grey between, and a field shows as a shrub of its younger lobes. Tripling the rough lava's depth, an experiment, left a ninth of the flow cells buried, but the fields stayed shrubs of separate lobes, each one tone with a sharp edge. The reference's mosaic lifts its blacks (its darkest 0.05 percent of the tropics at 0.36 of the median, where Io's dark paterae reflect about a tenth of its plains, as the model's do), so its tones are compared by their structure only. Left for the noise pass at the end of the surface phases (principle 16; the user, 2026-10-30): the tone within a field, a compound field's surface being a mosaic of breakouts of every age below the drawn lobes (principle 16); the metre-scale relief of inflated flows, which would keep old flows grey under metres of fallout; and the fronts' fall over the older lava's age.

Direction: the lobes of a field connecting where they meet, as later flows spread over and between earlier ones, keeping each lobe's shape.

### M41. The regional variation the physics does not make

Principle 16 (the user, 2026-10-29): where the modelled physics has gone as far as it can, the rest of the variation a real surface shows comes from noise modulated by a physical field, at every scale up to the global look. Two cases at the gate: the lunar highlands' mottle of shades, which on the real maps does not follow their height, and Io's grey-green regions. The user's suggestion for the Moon: a field of the crust's buried composition beside a global noise field, so that the pattern convinces without being uniform; for Io, a field such as the density of its flows.

Direction: one mechanism for every body, a noise field whose amplitude and tone a named physical field modulates, applied to the materials' colours.

When: one pass over every body at the end of the surface phases, after Phase 6, calibrated against many bodies at once, where other bodies' own regional variation joins Io's (the user, 2026-10-30; first set after Phase 4F's last checks, 2026-10-29).

### M42. Compression's mountain belts (moved to Phases 5 and 6)

PLAN.md's Phase 4F asks compression to run as one process from wrinkle ridges through lobate scarps to thrust-bounded mountain belts where the shortening is largest: a range like Phlegra Montes or the Thaumasia highlands on Mars, and Venus's belts round Lakshmi Planum with the episodic regime.

Measured (2026-10-30, at 1,024 pixels): the shortening the lid's thrusts take up, twice its overstress over its stiffness (F155), reaches 3.2 percent at most, on Mars at its rise's foot, past 1 percent over 0.3 percent of its area; Mercury 1.2 percent, all of it about 0.3 to 1.2 from its contraction; Venus and the Moon half a percent. That thickens a crust 50 km thick by a few hundred metres, scarps and ridges, not ranges. Venus's highest ground, 10.2 km, is a volcano; its plateaus, piled over the downwellings by the flow (F184), stand 2 to 5 km. A belt needs a push that goes on after the lid breaks, so that the shortening piles into the broken zone.

Tried, not kept (2026-10-30): the lid goes with the flow where the drag breaks it (F185, F187), its crust with it, first at the flow's speed and then as fast as keeps its stress at its strength, the drag being what the lid bears standing still. Both swept Venus's crust into its downwellings, where it piled past the eclogite cap and foundered: its mean crust fell from 35 km to under one within 0.4 Gyr, most of it in the first 50 Myr while the young lid was thin. Over the upwellings the broken lid stretches and gives up its crust, and nothing makes new crust there from the melt the stretching brings up: a lid that moves with the flow wants crust made at its rifts and lost at its convergence, the machinery of Phase 6. Mars's Tharsis is held up by its plume and does not spread; the Thaumasia highlands as the toe of a spreading plateau want the weak layer beneath it, the salt and ice of Montgomery et al. (2009), Phase 5's volatiles; Phlegra Montes's origin is debated. The plateaus' own collapse after each overturn (F188) presses their feet smoothly; concentrating it into belts would be a new process of the strain's localisation, and is not the leading account of Lakshmi's belts.

Moved by the user (2026-10-30): Venus's belts with Phase 6, Mars's with Phase 5 (PLAN.md).

### V1. Interpolated bodies in validation

`exogen validate` checks the reference bodies only. Principle 15 asks for bodies between them as well. The first sweeps to carry:

- an Enceladus scaled in radius (in the test suite now, S4 and S8);
- a rock body from the Moon's tidal flux to Io's (S3);
- an ice body from Ganymede's shell to Europa's (S4 to S6);
- a body from Mars to Venus in radius and heat flux (S1);
- a surface from dry to wet at Venus's and Earth's size (S1).

Each sweep prints the count or share of each feature class along the path and fails on a step where no named threshold is crossed.

When: each milestone adds the sweep for what it touches.
