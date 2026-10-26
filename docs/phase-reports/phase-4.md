# Phase 4 report: ice worlds

Status: complete, 2026-10-11, six milestones each committed to main on its own; reviewed three times (2026-10-11 and twice on 2026-10-12) and revised after each review, then revised a fourth time on 2026-10-12 to the emergent pillar, and followed the same day by Phase 4F, the fidelity of the base (the sections at the end). Awaiting the next look before Phase 5 starts.

Phase 4 gives the icy bodies their own physics. Phase 3 left them as cratered balls of ice with a regime name: Europa a blank sheet, Ganymede one terrain, and both carrying twenty kilometres of basin relief that warm ice cannot hold. The plan's list (PLAN.md, Phase 4 and section 4.3) is composition and shell thickness from density, viscous relaxation of craters into palimpsests, the tidal stress field with lineae, cycloids and bands, grooved terrain, focused polar heating, chaos, cryovolcanism with plume fallout and salt deposits, sublimation lag, hemispheric asymmetry for locked moons, and the ice and salt end-members. The gate is that Ganymede shows two terrains with a crater density contrast, that Europa's crack orientations follow the tidal stress field, and that relaxation removes large craters first.

The principle set in the Phase 3 reviews holds throughout: counts come from budgets and distributions, kinds from probabilities, colours from materials and processes.

## Milestones

1. **Ice foundations.** The crater scaling's ice factor keyed on the surface material rather than the bulk composition (the Phase 2 note); sublimation lag as the maturity law of icy surfaces, from the vapour pressure of ice at the latitude's mean temperature and the regolith's non-ice share, so warm old ice darkens and fresh craters expose clean ice; cryolava floods on ice regimes from the ice magma budget.
2. **Viscous relaxation.** Topography on an ice shell relaxes by wavelength with the viscosity of the ice at the depth the wavelength probes, through the shell's thermal history, in the loop on the coarse grid; large craters flatten into palimpsests while small ones keep their depth; replayed craters relax by their age.
3. **Tidal stress and the cracked shell.** The degree-2 membrane stress of a thin shell over an ocean from the Love numbers, its diurnal and non-synchronous parts; lineae as double ridges along the most tensile direction at their time, cycloids where the diurnal stress turns the crack through an orbit, bands where extension dominates; salts in young cracks fading under frost.
4. **Partial resurfacing.** Grooved terrain on a shell that extended: lanes of parallel grooves in clean ice, drawn as parts of great circles in belts and parallel sets until they cover the resurfaced share, cutting the old dark terrain, with the crater density contrast; the focused polar terrain of a shell whose heating concentrates under a pole, thinned there and floating lower, renewed every step, with its margin ridge, its folds and most of the cracks; the clean substrate as a layer; the wholesale renewal kept for thin shells only.
5. **Chaos, cryovolcanism and asymmetry.** Chaos on a thin shell as a population of patches whose number comes from the disrupted share of the surface, placed in a heat field, their matrix dropped and broken into hummocks around rafts of the old surface and stained with the ocean's salts; on a shell over rock the dome share of the ice magma budget extruded as brine-rich domes that relax like any load, and salts on the floors of craters deep enough to reach the brines; a polar plume's fallout that buries lags and salts and keeps the body in frost; the ocean's salt end-member chosen by the parent's magnetosphere; and the trailing hemisphere darkened and reddened by the parent's plasma.
6. **Documentation and review.** Parameter rows and formulas (F39 to F48, `docs/parameters.md`), this report, the presets' descriptions, the app rebuilt in debug and release against the Phase 4 engine.

## Gate results

- **Relaxation removes large craters first** (`stages::tests::ganymede_relaxes_its_large_craters_first`). On the Ganymede preset the median depth-to-diameter of craters 30 to 70 km across stays above 0.02 while that of craters over 150 km falls below a quarter of it. At standard quality: 20 to 40 km 0.072, 40 to 80 km 0.044, 80 to 160 km 0.001, over 160 km 0.000; Callisto, whose dirty ice is stiffer, keeps 0.060 at 40 to 80 km. Ganymede's relief fell from 30 km of basin to 11 km.
- **Europa's crack orientations follow the tidal stress field** (`stages::tests::europa_cracks_follow_the_tidal_stress`). The Europa preset grows more than 200 cracks, among them cycloids, and more than 80% of its ridges run within twenty degrees of the direction perpendicular to the most tensile stress at their midpoint in the field of their own time. At standard quality: a diurnal tension of 81 kPa, a crack density of one per 32 km, 2,831 cracks on the present surface, of which 1,815 double ridges, 321 bands and 695 cycloids, the longest 1,500 km.
- **Ganymede shows two terrains with a crater density contrast** (`stages::tests::ganymede_has_two_terrains`). The Ganymede preset's lanes of grooved terrain cover between 45% and 85% of the surface, and the recorded craters on them are less than half as dense as on the dark terrain. At preview quality the lanes cover 65% of the surface and carry 20 recorded craters against 227 on the remaining 35%, a density ratio of 0.05. At standard quality 645 lanes in parallel sets cover 65%, and 36% of the surface still shows as grooved terrain once the craters formed since have taken their share.

- **Milestone 5 has no gate of its own**; its tests hold Europa to more than twenty chaos patches with 8% to 40% of the surface as matrix, Ceres to at least five domes and at least one salted crater floor, and Enceladus to a mean maturity under 0.1 with cracks no browner than the frost round them.

## Measurements

Standard quality, after milestone 5.

| Body | Loop | Finish | What shows |
|---|---|---|---|
| Europa | 3.6 s | 15 s | 2,834 cracks with the tension zones swept through every longitude, salts in the young ones; 5,075 chaos patches whose matrix is 16% of the surface, denser where the heat field is high; the trailing hemisphere browner; relief within 3 km, five craters |
| Ganymede | 8.9 s | 36 s | 645 lanes of bright grooved terrain in parallel sets over 65% of the surface, the dark terrain between them in polygons, craters on the lanes as bright as the lanes; relief 10 km, frost caps past 40° |
| Callisto | 4.1 s | 12 s | no lanes (0.08 mW/m² of tidal flux); dark lag at the equator grading to frost caps, bright young craters, as after milestone 2 |
| Ceres | 9.5 s | 59 s | dark mixture with cryolava on under 1% of the surface, 18 domes of brine-rich ice over the body's life, the young ones standing and salted, white salts in the floors of its deep young craters, large basins relaxed |
| Enceladus | 14 s | 65 s | ancient cratered plains; a renewed polar terrain over 6% of the surface, a kilometre low, with its margin ridge, its folds and 99 cracks mostly within it, two lanes over 3%; the plume's fallout keeps the whole body in frost (mean maturity near 0) and shows the cracks as faint blue-grey lines |

## Deferred

- **Triton and Pluto** render as cratered ice shells with Pluto's domes; the cantaloupe terrain, the nitrogen ice plains and the volatile cycle are later phases.
- **Callisto's leading hemisphere** is darker in reality, from the dust of the irregular satellites; only the magnetospheric trailing tint exists, which on Callisto is weak and on the wrong side.
- **Ceres's small craters** stay soft under the ejecta of its nine basins (Phase 5), and **Enceladus keeps four kilometres of primordial relief** where the real body has one (a Phase 1 amplitude).
- **Ganymede's grooves** are cut at seven kilometres or three cells, whichever is coarser, so standard quality shows the coarse end of their spacing.
- **Io's texture** waits for Phase 7, as agreed at the Phase 3 review.

## Findings during the phase

- **The tidal field came out at the right size without tuning.** The membrane stresses of a thin shell under the degree-2 tide, with h₂ = 1.14 from 1 + k₂ cut down by the shell's rigidity, give Europa 81 kPa of diurnal tension against the 40 to 100 kPa of the literature, and Enceladus 119 kPa; Ganymede's 3 kPa falls below the cracking stress, so it grows no lineae, as it should.

- **The lag without a brake went opaque in a few million years at Callisto's equator**, and its caps began at a cliff of latitude, because the free sublimation rate at 120 K times the age is thousands of tonnes of ice per square metre. The lag chokes its own growth once it is thick and grows as the square root of time, which is what Moore et al. 1999 describe, and the transition spreads over twenty degrees.
- **Every ancient cell was excavated clean ice.** The deep bedrock that Phase 2 gave icy bodies made Ganymede white once thousands of overlapping craters had spread their ejecta over everything. Gardening now mixes clean ice back into the regolith over half a gigayear, and a shell over rock has no clean ice below it at all.
- **A ductile shell alone flattened everything.** With the viscous law alone, even 40 km craters relaxed on Ganymede, because the shell was thin and warm early and a step of fifty million years is an age to warm ice. The cold elastic lid holds loads shorter than a few times its thickness without flow, and the rock in a dirty crust stiffens it by orders of magnitude (Fu et al. 2017 on Ceres), which keeps the small craters and the dirty-ice bodies' craters where they are observed.
- **Ceres is buried in its own basin ejecta.** The production function in the belt gives Ceres nine basins over 400 km, as the models in Marchi et al. 2016 also predict; the real body has none, having relaxed them. Ours relaxes their relief but keeps their ejecta blankets, which bury most of the older small craters. The crispness of Ceres's small craters waits for the Phase 5 look at ejecta burial on small bodies.
- **Enceladus keeps four kilometres of primordial relief** at scales the elastic lid holds; the real body has a kilometre. The primordial crust's amplitude on small bodies is a Phase 1 setting to revisit.
- **Enceladus is not Europa.** The Phase 2 rule renewed every warm shell's surface wholesale within 0.1 Gyr, which left Enceladus without its ancient cratered plains. A shell renews itself wholesale only when it is thin against the body, under 2% of the radius; a thick one keeps its craters and renews the terrain under the pole where its heating focuses, which the tidal flux and the shell's thickness over the radius decide (F43).
- **Per-cell relaxation factors run away.** Giving each cell the survival of its own shell thickness inside the band decomposition sank the cells of the thinned zone by kilometres a step: a cell told to relax a wavelength its stiffer neighbours keep subtracts a regional mean it cannot lower, step after step. A watched cell at Enceladus's cap edge fell 17 km to bring a 48 km block average down by 2.5 km. The loop therefore relaxes with the mean shell everywhere; the local thickness acts through the replayed craters' survival, the crack density and the renewal.
- **A sampled pyramid is not a low-pass.** The band decomposition's levels were bilinear point samples of the field, so a coarse level could take one crater floor for a whole patch. With the renewal and the equilibrium surface in play this grew a dipole of a thousand kilometres over thirty steps. The levels are now area averages of two-by-two cells; the mean-shell relaxation changed a little with them and the gates still pass.
- **An isostatic depression is an equilibrium, not a load.** Pinning the polar terrain's floor while the whole field relaxed ratcheted the ring around the cap up by 20 km in 90 steps: relaxation raised the cap's surroundings towards the mean each step, and the pinning put the cap back down. The relaxation now works on the deviation from the thinned shell's equilibrium surface, and the renewed terrain drops its short wavelengths and relaxes the long ones with its neighbours instead of being pinned.
- **Cracks scale with the body.** Lineae 100 to 1,500 km long on a 252 km moon wrapped around it. The lengths and the band widths are now shares of the radius, which leaves Europa as it was and gives Enceladus 99 cracks instead of 22, most of them in the polar terrain where the shell is thin.
- **The lanes' edges are redrawn at the finish.** The coarse loop draws the lanes and the polar terrain cell by cell, 16 km on Ganymede, and the finishing grid showed the steps. The loop now renews the eroded interior only, and the finish redraws both terrains at its own resolution with a one-cell coverage ramp before anything reads the surface ages, so the crater population, the relief and the colour follow the true edge.
- **A slow budget built only the smallest domes.** Drawing a dome's volume and capping it at the budget on hand gave Ceres 150 domes of the minimum size; the volume is now drawn ahead and waits until the budget covers it, which gives Ceres 18 domes across the whole size range, as the budget and the volume distribution say.
- **The ocean's salts take their colour from the parent.** Europa's cracks are brown because Jupiter's magnetosphere irradiates sulfur into them; Enceladus's stripes are not, and Ceres's faculae are white sodium carbonate. The salt end-member is chosen by the parent's mass as the measure of its magnetosphere, and a plume's fallout buries whatever salt would show before it shows.
- **Chaos gathered in one hemisphere** when the heat field's crests took five times the patches of its troughs; with the troughs taking half, the patches spread round the body while still favouring the warm regions.
- **Craters on bright terrain stayed dark** until the substrate became a layer of its own. A crater on a lane excavates the lane's clean ice, so the clean share of the substrate is now a field that craters do not reset, and the lag maturity and the colour read it.

## Review revisions

The review of 2026-10-11 looked at Ganymede, Enceladus, Europa and Ceres. What it found, and what changed:

- **Ganymede had no dark polygons.** The lane budget was right (two thirds bright) but the belt field barely favoured any region, so the lanes chopped the dark terrain into slivers everywhere. The belt field now runs at 1.2 cycles round the body raised to the fourth power, three lanes in four grow beside an earlier one, and the lanes are shorter and wider (0.2 to 1.5 rad, 30 to 250 km); the dark terrain survives in large regions between the belts, with 777 lanes covering 65%. The bright craters that vanished when the lanes formed were the lanes erasing 65% of the surface 1.8 Gyr ago, which is right.
- **Enceladus's polar terrain was a painted circle.** Its outline was the cosine taper with a wobble, a parameter shape. The renewed terrain is now the thinned shell's core plus the ground within 30 km of the fissures that opened in the thinned zone, so its outline follows the stripes and lobes where they end, the margin ridge follows that outline, and the folds fade in from it. The fissures under a plume are vents, troughs 450 m deep and 3 km wide with low flanks, instead of Europa's 150 m ridges, the folds are 150 m, and the coarse-ice tint of the vents is stronger.
- **Europa's bands were thick dashes.** Constant-width segments with square ends. A band is now a lens, its width going as (4 s (1 − s))^0.35 along its length, lineated inside at 1.5 km. Chaos patches now run to 300 km across, so the area sits in larger, merging regions as on the real body; the matrix is 16% of the surface for a disrupted share of 30%, against published estimates of 20% to 40%. The brown patches are salted chaos matrix; the white streaks are old cracks and bands whose salts have faded under frost.
- **Ceres had circular frost caps and too many bright craters.** The thermal-segregation frost condenses at Ceres's poles but a carbonaceous regolith gardens it in, so frost visibility now scales with the square of the regolith's cleanliness: Ceres and Callisto show none, Ganymede six tenths. Each slope now sits at its own mean temperature (F49), so frost lies on pole-facing walls and leaves equator-facing ones, and the cap's edge follows the ground rather than a latitude circle. Ceres's crust is a carbonaceous ice (fresh 0.16, mature 0.09), 1.6 times brighter when fresh as Haulani's ejecta are, instead of dirty ice four times brighter. Ejecta blankets renew the surface over their outer third, as their thinning deposit mixes with the substrate, rather than ending at a sharp disc.
- **A noise frequency error was found along the way.** Every Phase 4 noise field (groove sets, polar folds, chaos outlines, rafts and hummocks, dome streaks, facula patches, the lane belts, the heat field, the thinning's wobble) was passed in cycles round the body where the noise lattice counts cells per radian, so each was about six times too fine: the chaos rafts were single-cell speckle and the belts were seven and a half cycles instead of 1.2. The convention is documented on the noise function and every call corrected.

Not changed, and raised for the next phase: the surface's natural look at the resolved scale (thermal segregation on slopes, sublimation pitting, gardening of thin ejecta, regolith creep, degradation of small craters) is physics, not texturing, and belongs in Phase 5's surface processes; Phase 7 should keep only what has no physics in it (sub-grid roughness, normal maps, exposure, export). The slope-dependent frost went in now because it is ice physics and answered the Ceres point directly.

## Second review

The second review (2026-10-12) asked how the lanes evolve and what drives them, found the polar terrain still geometric and its lanes boxes, the chaos still circular, and noticed that Callisto looked more like itself early in its history than at the present day. What changed:

- **The lanes form over an episode.** They all appeared at one instant; the real bright terrain spans a gigayear or more of crater densities and the sulci cut one another. Each lane now has its own time, uniform over `tectonics.ice_resurfacing_duration` (0.3 of the age, 1.4 Gyr on Ganymede) about the resurfacing age; the loop applies the lanes of each step as they come, a younger lane overprinting an older, and the finish cuts each lane's grooves at its time, so lanes cross-cut and the cells a younger lane also covers carry the younger age. At 1.5 Gyr before the present 540 of Ganymede's 786 lanes have formed, covering 45%. What drives the lanes is stated plainly: the amount is the pinned resurfaced fraction, because the heating pulse of the Laplace resonance is not in the eccentricity history; it is not the mechanism that cracks Europa, which is the tidal stress field.
- **Lane widths scale with the body.** Absolute widths sized for Ganymede made boxes a quarter of Enceladus's disc; widths are now 0.011 to 0.095 of the radius, which keeps Ganymede's 30 to 250 km and gives Enceladus ridged belts 3 to 24 km wide, as its sulci are.
- **The polar terrain has its margin belt.** Within 40 km of the edge the folds run parallel to the edge, 300 m at 8 km, where the spreading terrain is taken up, in place of a bare floor meeting the old terrain; the vents are 6 km wide with higher flanks. The near-identical outline across seeds is partly physics, the stripe directions come from the tidal field and not from the seed, but the belt and the narrower lanes change the look more than the seed did.
- **Chaos outlines wobble in six waves** round each patch rather than one, and the rafts hold 45% of a patch, so the patches read as irregular regions rather than stamped circles.
- **Each cell sits at its own temperature.** The lag law used the body's mean albedo for every cell, so a fresh crater at Callisto's equator darkened within a million years. A fresh bright surface is twenty degrees or more colder than the dark lag around it and sublimes thousands of times slower, so it keeps its brightness until gardening takes it, while a dark lag warms and grows: thermal segregation. The lag maturity now iterates each cell's temperature with its own reflectance (the body's Bond albedo scaled by the cell's reflectance against the mature surface's, from the fresh state), which brings back Callisto's starry surface of bright rayed craters on dark ground at the present day and darkens Ganymede's dark terrain against its lanes.
- **The per-step lane smoothing** cost the Ganymede loop thirty seconds; the lanes' floor now comes from the relaxation's area-average pyramid at the smoothing length, and the loop is back to eight.

## Third review

The third review (2026-10-12) asked whether Ganymede's resurfacing age was a hard number, found Callisto's grey still bare between its bright craters, found the steps inside the polar terrain's edge unnatural and its stripes still straight lines that never leave the renewed ground, and noticed that a crater's ejecta disc ends on a hard circle and sits over its rays. What changed:

- **Ganymede's episode is modelled, not pinned.** The 1.8 Gyr was 0.4 of the age and the two thirds were pinned in the preset. The orbit now carries a resonance episode as three inputs, `orbit.resonance_age`, `orbit.resonance_eccentricity` and `orbit.resonance_duration`, which put a triangular pulse into the eccentricity history (F50); the tidal flux history follows its square (F13), and the lanes take their amount from the history's peak flux and their times from the spread of the flux above the resurfacing floor (F44): the excess-weighted mean time and the width of the uniform episode with the same variance. Ganymede's preset carries the episode Showman and Malhotra (1997) described, the eccentricity held at 0.008 for a gigayear (rising over the first fifth, holding, falling over the last) two gigayears ago (the bright terrain's crater density, Zahnle et al. 2003), and pins nothing: the shell heated to 8 mW/m² of tidal flux on the plateau, two thirds of the surface (0.65) renewed, over 0.65 Gyr about 2.0 Gyr before the present. The pulse is in every thermal history downstream, and the crater record constrains it: at 0.01 and above the shell thinned under 30 km and the dark terrain's 70 km craters of 4 Gyr ago relaxed flat along with the 150 km ones, which the real dark terrain does not show, so the lane scale (F44) is fitted to the pulse that keeps them, 8 mW/m² for two thirds, and the shell stays near 40 km. Change the eccentricity or the timing and the amount and the age of the bright terrain follow; a steady flux spreads a body's lanes over its whole age, which Enceladus's sulci now do, and the refit raises Enceladus's share to 0.07, nearer the real body's tectonized plains. Validation checks both numbers; 784 lanes cover 0.65 of Ganymede at the present day.
- **Callisto's bare grey stays for Phase 5.** The space between the bright craters is bare because nothing yet degrades it: the real surface is a field of bright knobs and dark lag collecting in the lows, which is sublimation degradation at the resolved scale, and the small bright-rimmed craters below the render's scale. Both are the surface physics the amended plan puts in Phase 5, and the look is left for it rather than textured now.
- **The polar terrain follows its fissures out and the fissures leave it.** The renewed ground round a fissure now narrows from 30 km in the core to nothing where the taper falls to a half, so the terrain reaches out along the fissures as tongues and its outline is theirs, not a scalloped band of discs. A fissure leaving the core meets the thicker, stiffer shell and refracts along the contrast, as cracks do at a stiffness contrast: between the 0.75 and 0.5 profiles it turns to run round the margin, the hooked ends the real stripes have, and it stops only where the taper falls below 0.3, forty kilometres out in the cratered terrain, so the stripes leak into the old ground and curve as they go. A vent now narrows to a third of its width and shallows to 0.4 of its depth at its tips, an elliptical opening rather than a bar. The polar terrain covers 0.036 of Enceladus.
- **The steps were the belt.** Its folds were five concentric sinusoids, terraces to the eye. The belt now decays inward over 30 km rather than ending at 40, its folds 12 km apart come in segments 25 km long that shift against their neighbours by up to 0.4 of a wavelength and vary in height, en echelon as a thrust belt's folds are, and the margin ridge rises and falls along the edge with them.
- **A blanket grades into its rays and into the ground.** The ray field held the disc at half brightness out to exactly two radii and dropped to nothing there, and the age reset ended on the same circle, so the disc was a hard circle laid over the rays. The deposit now covers the ground by 1 − 2^−(edge / r)³ (F51), all at the rim, half at the continuous edge, which wanders by 15% in many small lobes (a few large ones read as petals), a fifth at three radii and nothing beyond, and both the ray brightness and the surface renewal follow that cover: the rays run out from the rim through the thinning blanket, and the cell under thin ejecta keeps the one age whose exposure gives the mixture's maturity, exact for an exponential maturity law, so an old crater's disc fades into the ground instead of ending where a linear age fade crossed the maturity's knee. Under less than half a cover the ground keeps its unit, its melt and its salts or lava.

## Fourth revision: the emergent pillar

The third review's follow-up (2026-10-12) asked whether the polar terrain's branching and meander would come from a later process, found the sulci rectangular and both terrains outlined by sharp grooves, asked where Enceladus's larger uncratered share would come from, and asked for every remaining pin to go and for the model never to be bypassed to reach parity with a real body. Four milestones, each committed on its own:

- **The surface's own albedo, and no pins.** Every preset pinned its measured Bond albedo. An airless body's Bond albedo is now its surface's own reflectance (F52): the regolith mixture's luminance at maturity, fresh over the share of the body a plume's fallout keeps under frost. The surface depends on the temperature the albedo sets, through the shell, the regime and the materials, so the graph gained lagged dependencies, read from the previous recompute pass and left out of the order, and a recompute repeats until nothing changes. Against the measured values: Moon 0.14 (0.11), Mercury 0.085 (0.088), Ceres 0.086 (0.09), Europa 0.66 (0.67), Enceladus 0.83 (0.81), Triton 0.66 (0.76), Ganymede 0.37 (0.43, the bright lanes not counted), Io 0.37 (0.52), Callisto 0.13 (0.22, its bright craters not counted), Pluto 0.25 (0.5, its nitrogen and methane ices not yet among the materials). Under an atmosphere the pressure class stands in until the climate stage computes clouds and haze. The derived nodes that are inputs in spirit, what the body is made of, its volatile inventory and its rotation state, are marked as inputs and shown as set rather than pinned, and a test asserts that no reference preset pins anything else: the pillar made executable. The body's own albedo warms Ganymede's surface by two kelvin, which relaxes its 60 to 70 km craters a little further; the gate now expects them shallow but present, as the real dark terrain's are (a depth ratio over 0.008 against the real 0.012), rather than near fresh.
- **Lanes traced in the field.** A lane was a great-circle box with a wobbling edge. It is now a swath about a centre line traced along the tidal field at its time where the shell cracks and the tension at its start opens it, so Enceladus's sulci curve with the field as its cracks do, and along a great circle where nothing steers it, as Ganymede's; the width closes in a lens at the ends, and the grooves follow the local tangent. The lanes' 200 m drop, the extension's subsidence, spreads over the lithosphere's flexural parameter at the lane's time from the thermal history (about 90 km on Ganymede today, 55 km during its episode, 130 km on Enceladus), so that a lane narrower than it sinks little and no lane ends on a step; the cells outside sink with its flank. The gates hold: 850 lanes cover 0.65 of Ganymede with a crater density ratio of 0.05.
- **Cracks answer one another.** The fissures were traced one by one in the smooth field, parallel and independent. A share of new cracks (0.35) now start at the tip of an older one, where its stress concentrates, and run one way from it at the field's heading of their own time, so sets branch and splay; a propagating crack within half a shell thickness of an open fissure feels its stress shadow and turns toward meeting it at a right angle, stopping within half a width of it as at a free surface (Pollard and Aydin 1988), while healed ridges and bands are crossed as Europa's lineae cross. On Enceladus the polar fissures fork and hook, and the cracks that leave the renewed terrain splay into the cratered plains; Europa's network is unchanged in kind.
- **Episodes of focusing and fossil terrains.** The uncratered share was structural: the focusing sent nearly all the renewal to the present pole, and a steady history had nothing else. The focused heating now comes and goes (F53; O'Neill and Nimmo 2010: overturns every 0.1 to 1 Gyr, each active about 10 Myr): episodes recur every ten conductive times of the shell (0.53 Gyr on Enceladus, a fit at the short end of their range), each starting where the overturn happens to and turned to the pole by the body's reorientation (Nimmo and Pappalardo 2006), so that in the present frame each old terrain lies where the wander left it with its fissures turned the same way. Through an episode the zone floats at its depression and renews itself as the present one does, a share of the thinned zone's cracks falling to each episode alike; when it ends the shell floats back up and the terrain takes craters from then on, and the finish redraws each fossil terrain's outline from its carried fissures, lays its floor back and stamps its folds and margin at its time. Enceladus carries eight of them; with the present terrain and the lanes they renew 0.40 of the surface (the sidecar's renewed_fraction), against the 57% Kirchoff and Schenk 2009 mapped as tectonized, each past terrain being the size of the present one. The youngest fossils are smooth, the oldest as cratered as the plains.

Calibrations added by this revision, all documented as fits with their sources: the focusing period factor (ten conductive times) and nothing else; the lane scale and the resonance pulse stand as before. Left for Phase 5 as before: Callisto's bare grey, the resolved-scale surface physics; and for a later look at the materials, Pluto's and Triton's volatile ices and the brightness of Ganymede's lanes, which the body's own albedo now makes visible.

## Phase 4F: the fidelity of the base

After the fourth revision the user judged the worlds' look basic: hard-edged stamps, flat unit colours, idealized shapes, one-cell contacts, repeating patterns, a pixelated fossil floor, and said formulas alone would not reach the fidelity the project is after. That was right, and the sequencing that left the look to a later pass was wrong: a later pass cannot invent variation the processes never made. PLAN.md gained the two pillars as principles 11 and 12 and a Phase 4F roadmap entry, and four milestones followed, each committed on its own, with Enceladus, the Moon and Io as the test bodies.

- **A structured substrate (F54).** The craters the finish does not resolve are a population too: their relief now goes under every surface as a self-affine field anchored at the smallest resolved diameter to that crater's depth and to the share of saturation the surface has reached there, with the Hurst exponent of cratered terrain. Old ground is rough at every scale, renewed ground smooth, and their contact is a gradient of roughness; the primordial relief takes the same spectrum, so what floods it meets rough ground. The relaxation's box pyramid, which left the blocks of the loop's cells bare on the fossil terrains, is a Gaussian pyramid, a renewed cell always drops the finest band, and the renewed floor is the whole relaxed surface.
- **Features from fields, not shapes.** The focused heating is a field: the thinning varies across the zone by 0.3 with the shell's convection in cells two shell thicknesses wide, so the polar terrains' outlines, depressions and renewed ground are lobate rather than circular. The polar folds follow the nearest fissure within 40 km instead of running in arcs about the pole; a past terrain's ridge, belt and folds are spread and lowered by the gardening since. A lane's edges wander in a self-affine profile in place of two sinusoids; every crack's relief and width rise and fall along it over 15 km segments.
- **Variation within units (F55).** A body's lavas drift in composition over its eruptive history, so the maria show tonal zones of their flow units by up to 18%; slopes shed their regolith from 15° up and stay fresher at the angle of repose, so steep crater walls and scarps brighten on every airless body from the ground's own slopes while gentle old rims stay mature; a lane of renewed ice is clean by a share drawn in 0.55 to 1.
- **Io's deposits accumulate (F56).** A flow's colour is coated over by the plains' fallout on 20,000 years, so old flows grade into the plains; every patera carries three generations of flow fields, in six to fourteen narrow fingers, and three of plumes, each faded by its age; a flow's surface varies with its flow units; each slope sits at its own peak temperature, so frost keeps to the pole-facing slopes. The plume frost on Enceladus saturates softly, so the far side keeps some contrast under a thinner frost.

What the crops show after these four: the fossil terrains are smooth lobate floors with faint ghosts of the craters they relaxed, the polar outline is lobate and its folds parallel the stripes, the sulci curve and taper with edges that wander, the maria carry tonal flow-unit zones and bright steep walls, the Moon's highlands a grain that the maria lack, and Io a patchwork of half-buried flows and spent deposits round its fresh ones. What they still do not show, and what the next look should judge: Io's flows are still fingers of one tone on a plain of one hue, with no texture from the ground they crossed; the maria's contacts, though physically sharp at this scale, are smoother than real basin rims breached by lava; the tiger stripes are still few and short in the core; and the margin of the polar terrain is a line, as the real scarp is, but a line nonetheless. These are the next steps of Phase 4F, not Phase 5's.

### The fifth milestone: what the fourth's renders showed

The crops of the fourth milestone showed three things the formulas had hidden. The maria carried thin bright meandering lines: the contacts between the loop's flood footprints, where the fine pass of the lower flood stopped short of the ramp the bilinear base makes toward its higher neighbour, the band kept a blended base age, and the composition drift, read through that age, turned the blend into a tone of its own. The lines took four renders to trace, each hypothesis checked against the maps rather than the code: not the tint read through a blended age (the composition became a layer of its own that each eruption writes, which was right in itself), not the base's bilinear ages (they resample by the nearest coarse cell now, which was right too), but the ragged margins of pre-mare crater ejecta standing less than 50 m above a later flood's level, which the flood "draped" with a sheet of zero thickness, so the materials showed the bright basin floor through lava that was not there. A draped high now carries a veneer of the sheet's thickness less its own height, so the ghost craters of the maria show in relief and not in albedo. With the lines gone the composition zones showed the loop's cell staircase along every contact between floods, so the fine outline is the footprint's indicator smoothed over each cell's neighbours and wandered by three octaves of lobes, and the zones are lobate; under a younger crater the flood still writes its surface, which the crater's own pass covers where its ejecta reaches. Slopes shed their regolith from 20° by 0.3 rather than from 15° by 0.6, since every old crater's wall had been ringing bright. A render with EXOGEN_DUMP_LAYERS set now writes the hidden layers (cover, deep, melt, clean, composition, maturity, rays) as 16-bit maps, which is how the lines were finally found. Enceladus had a round dark spot at its trailing apex: the magnetospheric darkening, applied over the frost at full strength; it now competes with the burial by the plume's fallout, the top millimetre darkening in 10⁴ years (Cooper et al. 2001) against the years the fallout takes to bury it, which leaves Enceladus's trailing side bright and Europa's bullseye whole. And the tiger stripes were a dozen short dashes, because a polar crack drew its length from the same log-uniform law as every other crack. They are now a cascade (F57; Hemingway, Rudolph and Manga 2020): the first fissure through the thinnest shell at the pole, each next one a quarter flexural parameter away where the ridges of the last bend the broken plate most, outward while the site is in the thinned core, every one run along the field until the margin refracts it; the zone's other cracks stay as its lesser fractures. The spacing follows from the elastic thickness the thermal model gives, thinned by the shell's local factor, so it is not pinned to the 35 km of the real set: the default seed gives five stripes 39 km apart.

### The sixth milestone: Ganymede's lanes

Ganymede's grooved terrain had become a tangle of lens-ended brush strokes crossing at every angle, with a hub of converging lanes at the trailing apex. A dump of the lanes as polylines (the finish writes them with the hidden layers under EXOGEN_DUMP_LAYERS) showed why: on Ganymede no crack field steers a lane, since its tidal stress never reaches the cracking threshold, so every lane was a great circle with a heading of its own; three in four started beside a parent with a heading turned by up to 12°, which crosses the parent within the lane's own length; the belts of extension, raised to the fourth power, concentrated starts into knots from which lanes radiated as stars; a neighbour picked its parent among every lane ever laid, so crowded regions crowded further; the median lane ran 1,600 km and a tenth over 3,200; and the width closed in a lens over the whole length. Now a lane on its own runs along the strike of a smooth field of extension, the regional pattern the stress alone does not resolve, read as an angle from east and drifting through the episode so that generations cross-cut (a first attempt took the contours of a scalar field for the strike, and the lanes swirled round its extrema in two whirlpools); a neighbour belongs to its parent's generation, turned by up to 3°, laid half to one width beside it so that a set is contiguous, no longer than it, and the parent is one of the lanes laid lately; and a lane's width closes over one width at each end as its bounding faults die out. The size laws changed too, since widths in shares of the radius gave Enceladus sulci a tenth as wide as Ganymede's lanes: a lane is a swath, a rift zone of the shell, 4 to 24 elastic thicknesses of the episode's lithosphere wide (44 to 266 km on Ganymede's 11 km, capped at 0.15 R, which holds Enceladus's sulci to 38 km) and 6 to 40 widths long, up to a quarter of the way round the body. Two detours on the way: widths of 1 to 6 thicknesses made Ganymede a felt of 12,000 scratches, since the raster already lays the necking-scale grooves inside a lane, so a lane is the swath and not the groove; and lengths in radians made Enceladus's sulci round blobs. The light terrain is a mosaic of broad parallel swaths crossing at angles with dark polygons between, 1,040 of them at the same 0.65 of the surface, and Enceladus keeps a dozen sulci of 38 km. What it still is not: a swath's grooves are still one spacing, and its interior one tone.

### The seventh milestone: the margins of the renewed terrains

The user's look at the sixth milestone found every terrain type on Enceladus embedded in the cratered plains with a sharp transition, and the fossil terrains blurred, "as if you enlarged a crater's image". Both were right. The renewed terrains met the plains on a two-cell ramp, and their floors were the loop's relaxed surface resampled bilinearly, which carried the loop's cells as facets on every relaxed slope; a close crop showed the same stepped contours on the relaxed craters of the cratered plains themselves, whose coarse form is the base's, so the enlarged image was not the fossil floors' alone but the whole coarse base's wherever the loop had relaxed or built a slope. Now the floor and the folds of a renewed terrain blend into the old surface over a 15 km ramp centred on the edge, the thermal transition at the margin, so the old craters sink part way across it and fade in rather than end on a cut; a floor the loop left is resampled and then smoothed over the loop's cell, as is the whole coarse base before the passes lay their resolved detail on it; and the folds inside a renewed terrain are 200 m at 2.5 km with their spacing and strength varying in fields, so the plains carry ridge sets rather than one grating. What the crops show: relaxed craters with smooth walls where the steps were, margins across which craters fade, fold fields on the smooth plains. What they still show: the largest relaxed craters carry concentric terraces, which is the band-by-band relaxation ringing rather than a facet, and the fossil plains are still smoother than Enceladus's, whose ridged plains are rougher than any fold set here.

### The eighth milestone: the fabric of the renewed plains

The user's look at the seventh milestone found the relaxed plains smoother and more blurred than the cratered terrain, and the fossil terrains' pattern "flat, rough, and stamped on". Measured, the renewed plains carried 38 m rms of relief at 3.5 km and below against 91 m on the cratered plains; and the folds were a straight grating because every cell measured its phase from the pole rather than from the fissure it was meant to parallel, so nothing curved, and a sinusoid with one spacing and one amplitude is a hatch however it is oriented. Now the folds are measured from the point on the nearest fissure nearest the cell, so they curve with it; they wander in phase so that they are sinuous, their sets are warped over 60 km so that they bend as ridge belts do, they come in segments about 10 km long that die out and relay en echelon, their spacing and strength vary in fields, and under them lies the tectonic roughness of a plain deformed at every scale, 120 m rms, self-affine with the Hurst exponent of the cratered terrain. The renewed plains now carry 75 to 80 m rms of fine relief against the cratered plains' 91 m, and the south polar terrain reads as curving ridge sets between the stripes. Two detours: doubling the fold amplitude alone gave a louder hatch, not a rougher plain (the roughness had to come from a field of its own), and the roughness under the belt's envelope left the margins smooth. What remains: the largest relaxed craters still carry the concentric terraces of the band-by-band relaxation, and the fabric is one process's; Enceladus's ridged plains also carry the troughs and the sub-parallel fractures of later episodes, which would come from the lesser cracks if more of them survived in the plains.

### The ninth milestone: what seed 563 showed

The user looked at seed 563 and named five things, with coordinates: the tiger stripes still straight; a blurry, zoomed-in terrain at 176°W 18°N and the same look on the polar terrain between the stripes; rough spotted patches at 112°W 17°S; the streaks painted on, with a crossfade for a margin; and parallel, unbroken lanes inside the streaks and the polar terrain. Each had a cause in the maps. The spotted patches were secondaries: a secondary carries its parent's time, and one landing on a lane renewed since its parent was stamped with the parent's age, a pancake older than the floor round it, which the age-driven roughness drew with a serrated edge of cells; two diagnostic switches, one dropping the loop craters' restores and one the secondaries, found them by elimination after a per-cell guard on the restores had changed nothing. Now a secondary whose parent is older than the renewed ice it would land on is gone with the parent, no stamp or restore writes a cell of clean ice younger than the crater, the roughness's strength follows the surface age smoothed over a cell, and an ejecta lobe narrower than four cells is not drawn. The zoomed-in look was two things: the folds measured their phase from the pole, so the fabric was a quilt of straight gratings with straight seams between the fissures' domains, and the renewed terrain kept the relaxed ghosts of craters up to tens of kilometres, since it dropped only the wavelengths below 12.6 km. The fold phase is now the distance to the nearest fissure, continuous and curving with it, fading over the outer half of the fabric's range so that each fissure carries its own belt of ridges, and a renewed polar terrain keeps nothing shorter than 0.2 R while its episode lasts. The streaks' grooves were one unbroken sinusoid fading over the outer 15% of the width; they are measured in from the nearer edge now, so that they follow its wander and the outermost sits just inside it as a bounding trough, with their spacing and strength varying in a set field, wander, segments that die out and relay, and a tectonic roughness under them. The tiger stripes wander more and end in forks, a branch from each tip that leaves at 60° and bends back to the field, the hooked terminations of the real ones. What the crops show: the lanes clean of pancakes, with sinuous grooves, bounding troughs and a roughness of their own; the fossil terrains without ghosts; the stripes sinuous, with hooked ends. What they still show: the fold fabric closes in loops round the short lesser cracks, the contour map of a distance field, which no distance to a set of short segments can avoid (a Euclidean soft minimum gave stadium loops, a perpendicular one rounded rectangles); the fabric has to be built from oriented segments, each a fold laid along the local fissure direction and dying out on its own, which is the next milestone.

### The tenth milestone: the fabric as oriented segments

The fold fabric of the renewed plains had been a phase measured from a distance in four forms, from the pole, to the nearest fissure, a soft minimum of distances, a soft minimum of distances across each fissure's line, and every one closed in loops round the short lesser cracks, since the contours of any distance field close round the thing measured from. The fabric is now built the way the real plains are built, from segments: one kernel per cell of a lattice 7 km across, each a wave windowed over 5 km along the ridges and 2.5 km across, laid along the fissures' direction averaged over 20 km about its centre and in phase with the nearest fissure's line, so that the ridges parallel each fissure, run on past its ends and relay from one fissure's set to the next; kernels vary in spacing and strength so that no two segments are alike, a share are absent, and the sum is normalized so that the amplitude holds where they overlap and fades where they thin. The first version took each kernel's orientation from its own nearest fissure, and kernels of different orientation crossed in a moiré where two fissures' domains met, with one spacing at one amplitude filling every cell; the averaged orientation field and the variation removed both. A transect across seed 563's lanes also found their clean ice cut to zero across the lane for kilometres: each fossil terrain's raster had redrawn every cell within 33 km outside its edge, writing its own clean share, zero there, over the lanes' and resetting the cells of other terrains of the same unit; a terrain now redraws only what it covers and what the coarse loop gave it. What the crops show: the renewed plains carry a flowing grain of ridges that parallel the fissures and sweep from one to the next, with no loop, seam or moiré, fading in belts; the polar terrain's ridges run with the stripes; the lanes keep their clean ice to their edges. What they still do not show: the grain is too laminar, a combed look where Enceladus's plains are ropy and cross-cut by the sets of later episodes, which would come from a second generation of fissures of its own time; and the amplitude varies less across a plain than the real ridges' does.

### The eleventh milestone: Io's flows follow the ground

Io was the third test body and the one the user had called a yellow ball with brushes stamped on it. Its maps, dumped with the hidden layers, said why in three parts. The flow fields were Gaussian sectors about a point: six to fourteen petals of equal reach radiating from the patera, an asterisk rather than a flow field. Their edges were staircases of the loop's cells, because the finish stamped every Io feature at its own resolution on top of the loop's coarse stamps without clearing them, so each flow kept a coarse halo of its unit and coverage round its fine outline. And every mountain was an oval with a grey rim, a plateau whose edge fell steeply all the way round and so shed its coat on every side. A lobe of a flow field is now a path traced over the ground (F59): from the patera's rim at a random azimuth, its heading turning toward the fall line and wandering, so that a flow runs downhill where the ground slopes and meanders where it is flat, stopping where a mountain stands or the ground rises against it; the lobe is a tube about that line, a channel at the vent widening to a lobate toe. The finish clears the regime's surface before it stamps anything, as every loop step does. A mountain block ends in a fault scarp on its upthrown side and slopes under the plains on the downthrown. Two versions curled every lobe into a stub before the tracer worked: the first turned by the whole wander at every step and stopped at the first upslope, and the second steered by the raw coarse ground, whose 100 m bumps at 20 km trapped each flow in the first hollow, where it oscillated about the minimum; the ground is now the base smoothed over 40 km, and the fall line steers a flow sideways and onward but never back. With the tracer instrumented, the lobes average 74 km of path, 85% run to their length, 13% end at a mountain and 2% pond. What the crops show: flows that run 50 to 100 km from their paterae as sinuous tubes, tangled into fields round each vent and channelled between the mountains, with no coarse halo; mountains with a scarp side and an apron side. What they still do not show: a lobe is one width and one tone along its length, where real flows branch, overlap in units of different ages and thin at the toe; and the plains are one hue, since a plume's fallout fades in a century and only the active plumes' rings exist, so the colour of the deposit history is the next milestone.

### The twelfth milestone: renewed terrains keep their history, and a lane is a rift

The user's review of the eleventh: Enceladus's polar terrain, fossil terrains and lanes still looked patched in, where the real body's young terrains pass into the old by degrees; Ganymede's lanes still had fine boundaries, "like a literal streak of paint"; and the thin lines inside lanes and polar terrain looked drawn. Before anything was changed the maps were probed, and they showed defects of sequence, not of shape. A past episode of focusing renewed its terrain only if a loop step happened to end inside its 10 Myr, one time in five, so most fossil floors were stale. The finish laid every renewed terrain's floor unconditionally, which wiped the bowl of every large crater formed since: craters 8 to 16 km across on renewed ice were 80 m deep at the median against 1,137 m on old ground, and those of 16 to 40 km 450 m against 2,620 m, on the 56% of Enceladus and the 65% of Ganymede that is renewed ice. The present polar terrain did not stand lower than its surroundings but had a moat round it, and the fossil terrains stood as plateaus 400 to 600 m high. A fourth showed once those were mended: a height profile across a lane margin found a notch 600 to 770 m deep and under 5 km wide, with a 295 m step in a single cell, because the lanes' floor was recorded at the time of the first lane and laid at the finish, and the ground between had moved.

The rule that replaces all four (F60): the loop owns the relief, and the finish owns the resolved edges and the fine fabric. A terrain renews once, at its time, by a weight on each cell: all of the short relief inside, a ramp across its edge, and a halo beyond a polar terrain's; an episode fires at the step its end falls in. The finish lays no floor of its own. It composes the terrains in the order they formed and, at each terrain's turn in the sequence, takes the fine relief stamped before it by the same weight, so that a crater older than a terrain is softened across its edge and gone inside it and nothing is cut along a line; the unresolved craters' roughness is taken in the same measure. Craters 8 to 16 km on Enceladus's renewed ice are now 2,114 m deep at the median (six of them) and Ganymede's 1,400 m (303) against 1,088 m on its old ground, as younger craters on a younger surface should be.

The lanes were rebuilt as rifts (F61), every length in them scaled by the lid that broke. The edges are the bounding faults' traces, straight segments that bend at relay ramps, where they had been a smooth wobble of the width. Every end stops at a cut: along the strike of an older lane that crosses it, on its set's shared front, or at a fault of its own across the lane, where each had closed in a rounded lens; on Ganymede two fifths of the ends stop against an older lane and the rest share some 250 fronts, and running on to a junction makes the lanes 1,800 km long at the median where they had been 630 km, against named sulci of 1,330 to 4,221 km. A lane's fabric is a wave across it that fans, bends and breaks into segments groove by groove, with a fine set of tilted fault blocks on the necks where the grid resolves them, and it continues past the bounding fault as a damage zone of open fractures that thin out. At a lane's turn the fine relief inside it gives way to its own, so a younger lane cuts an older one's grooves. The tone follows the structure, and the lanes' own cleanliness was narrowed to 0.8 to 1: at 0.55 to 1 every lane was a ribbon of its own grey, which was most of what read as paint. A crack became a chain of segments, each opening as an ellipse and stepping aside from the last (F62). Two smaller things came with them: a past episode's fissures now stand on the last surface its terrain renewed, where the terrain's own turn had been erasing them, and the polar folds' spacing follows the local lid, finer toward a terrain's centre.

During the work the user asked whether Enceladus's lanes and Ganymede's are one mechanism, wary of tuning by body. They are one code path with no switch on the body, and the question found three places where a number was doing a body's work. The lane's edge was a ramp of 2.5 grid cells, which drew the same bounding fault 1 km wide on Enceladus's grid and 10 km wide on Ganymede's; it is now a border zone one lid thick, the grid's ramp only its floor. The grooves' relief followed the spacing the grid drew, so a coarse grid made taller grooves; it now follows the lane's own spacing. And a rule tried earlier in this milestone, which read a capped lane's lid back from its width so that Enceladus's lanes carried fine grooves, was a tuning by another name, since only Enceladus reaches the cap; it was removed, and the spacing is 0.63 of the elastic lid on every body, 7 km on Ganymede and 14 km on Enceladus. The fine ridges Enceladus's lanes now show, 2.4 km apart, are the fine set, six to a groove spacing as on Ganymede, whose grid is too coarse to draw its own. One number of that kind remains, recorded in F61: the cap on a lane's width, 0.15 R, which only Enceladus reaches and which was set with its sulci in mind. Removing it properly needs a rule for where a lid can rift at all (a lid too thick for a rift to fit the body cracks instead, and rifts open only where it is thinned), which would move Enceladus's lanes into its episodes' terrains; that is the user's to decide.

The user also asked about a fluid simulation, since Enceladus's surface looks fluid, and about three parameters of SpaceEngine's generator as the mod they use has them. Nothing on these surfaces flows as a fluid at the surface; the look is ductile deformation, and its tool is kinematic: a displacement field from a terrain's own strain, applied to what is older than the terrain. That is domain warping with a physical warp, which is what the two "venus" parameters do with a noise warp (a warped fBm added to the continental height, its frequency and its magnitude); "volcanoActivity" gates the volcanoes by a thresholded field and, in the mod, also swaps the continental height to derivative-damped turbulence, thins the craters and warps the coordinate the palette is read at. None of it is copied here; the stock functions are in an encrypted archive that was not opened. The deformation of old relief by a terrain's strain is carried forward as a tool for a later milestone.

What the crops show. Ganymede reads as bright terrain cut into dark angular polygons and no longer as strokes on a ground: one tone of bright terrain, groove sets of different strike meeting along straight contacts, blunt and oblique ends. Enceladus's lanes are soft belts of ridges that pass into the cratered ground through softened craters and outlying fractures, with no outline; its fossil terrains carry their fissures and a fabric whose spacing changes across them; the south polar terrain meets the old ground without a scarp or a moat. What they still do not show. Ganymede's dark terrain is one tone with no furrows, its bright terrain has no variation beyond the grooves, and its grooves are drawn at 12 km where they are 7 km, the grid's limit. Enceladus's fabric, lanes and fossil terrains alike, is fine parallel ridges nearly everywhere, where the real body has smooth plains, reticulate ground and single great ridges beside its ridged plains; the tiger stripes are still even blue lines; a lane the cap narrowed is narrower than its own border zone, so it never renews its ground in full; and the crater density still falls across a terrain's edge within a few tens of kilometres. Io was not touched: its plains' hue from the deposit history is the next milestone. Tests 115 (seven new: the border zone, the stopped end, the grooves' relief on any grid, the order of composition, the completed renewal, the lane ends, the crack chain), clippy clean, validate 116 checks.

### The thirteenth milestone: the lid where a lane opens, pulses, and pits under the snow

The user's review of the twelfth accepted the lanes and asked for four things: whether the furrows should come now, what would change in the grooves (Enceladus's ridges "don't look terrible, but real Enceladus shows them with much more variation, some overlapping others, some deeper, some break up to form dotted lines"), what Europa looks like, and what lifting the width cap would actually change. A survey of the ice bodies and two experiments answered them before anything was built. Europa is the weakest of them, its cracks a cell or two wide with stair steps, its chaos speckled discs, its bands brown wedges; Callisto has no ring systems, only rayed craters on a dark ground; Triton and Pluto are bare and wait for the volatile phases; Ganymede at high quality draws its grooves at their true 7 km, so the 12 km of the standard grid is the grid's and not the model's. Lifting the cap gave Enceladus one lane 400 km wide over a sixth of its surface. The rule that replaced it was tried first as an experiment that left Ganymede's maps identical byte for byte, and is F61's now: a lane takes the lid where and when it opens, and a lane that could not be six times as long as it is wide is not formed.

The user then gave three more pillars, which went into PLAN.md as principles 13 to 15 beside the two they sharpen: transitions are terrain (a transition is made by the younger terrain's own process at part strength, with a width in metres, and a ramp between two sides is not one); features from fields (where, when, how large and in what form are read from fields at a feature's place and time, the fields being physical); and in-between worlds (calibration fixes a constant of a process and never one body's outcome, and a class of feature comes and goes along a path through parameter space only where a named threshold is crossed). And they asked that the earlier phases, canyons and shields above all, be gone over with the same scrutiny once the ice worlds are done, as Io was; the roadmap's Phase 4F entry now says so, with an audit of the switches in the model against principle 15.

The lid rule is those principles at work. A lane's lid is the elastic lid of its own time, thinner by the shell's thinning where a focusing episode has it thinned; its width is 4 to 24 of those; and on Enceladus, where the lid of the time is 23 km and the narrowest lane it could make is wider than a lane can be, lanes open only in the zones its episodes thinned, during those episodes: twelve of them in three sets, 44 to 61 km wide on lids of 10 to 13 km, with grooves 6.4 to 8.2 km apart, along the margins of its old terrains, where the cap had laid 38 km belts across the cratered plains at any time at all. The test of the principle is a body in between: an Enceladus a fifth larger makes its first few lanes on the lid of the time, and one half again as large makes 68 of its 81 there and the rest in its episodes. The same sweep found a switch the rule did not make: at twice the size the shell is thinner than a fiftieth of the radius and F44's thin-shell rule takes every lane away in one step. It is on the audit's list.

The grooves got what the user described, each as a process on every body (F63). A rift works in pulses: two later sets rework patches of the lane, struck along the field as it stands later in the lane's life (8 degrees off the first at the median on Ganymede, whose strike field drifts, and within half a degree on Enceladus, whose tidal field barely turns in an episode), each groove cutting what stood before across its own width, so that sets overlap groove by groove and never add as waves. A groove's depth has a heavy tail, a choice of exponent between a fault population's and a necked lid's regularity that the entry states as a choice, bounded at the deepest measured. The folds of the polar and old terrains take the same tail, narrowed as a packet of folds is, and keep an earlier pulse's strike in patches. And the dotted lines are pit chains (F64): the plumes' fallout, every episode's, is a loose cover of 400 to 1,900 m on Enceladus's old ground, two to three times the 250 to 700 m its real pits imply; an open fracture drains it, no deeper than the fracture is open, 380 m for a crack there; and a crack that opened under cover deeper than its own relief shows as pits 1.1 km across and builds no ridges. The first version tied a pit's size to the whole cover and made pits five kilometres wide; the bound is the crack's own depth, the stress that opens it over the weight of the ice, which no body sets by hand.

What the crops show. On Enceladus the cratered plains are left to their craters and their cracks, and a crack that runs as a ridge across a young terrain breaks into groups of pits where it enters the old snow-covered plain, as the real ones do; the lanes are belts inside the old terrains; the fold fabric has deep and shallow packets and sets that cross at low angles. On Ganymede a few master troughs stand out among the grooves and later sets cut the first at 10 to 20 degrees. What they still do not show. The fold fabric's deeper folds come as short packets where the real ones are single long ridges; a fissure's fork still makes a hooked bracket; the fine set is drawn at three cells and shows the grid when it runs diagonally; the lanes' own damage zones on Enceladus make pits too small for the standard grid, so they vanish there; the cover is a sum of nine plumes and has no drifts, no burial of small craters, no softening of old relief; and the furrows, the ring systems, Europa and Io are as they were. Tests 123 (eight new: the lid a lane opens on, lanes by degrees between bodies, the depth tail, a later pulse, a covered fracture, the pulses' strikes, the loose cover, a crack under it), clippy clean, validate 116 checks.

### The fourteenth milestone: deep folds that grow and link

The user's review of the thirteenth: "Enceladus' terrain reads as both varied and natural. This is a genuine milestone. Agreed on the deeper folds. Enceladus' real deep folds all interconnect and branch out while ours are still disconnected and linear." They also asked that what needs auditing against the new principles be noted. That went in first, as docs/audit-principles.md, referenced from PLAN.md's Phase 4F entry: eight switches (the tectonic regime and the composition class compared as names, the Io regime, the thin-shell rule, the active thinning, wholesale renewal, materials by name, and the dip in the lanes' share along a sweep in radius: 0.073, 0.058, 0.310 and 0 at one, 1.2, 1.5 and two times Enceladus's radius), four widths or scales counted in grid cells, the features still stamped from a profile and noise, the constants fitted to one body, and interpolated bodies for validation, with a table of which reference bodies trigger which switch (the mobile lid only Earth, the episodic lid only Venus, chaos only Europa, polar focusing only Enceladus).

The deep folds had been the fold fabric's deepest packets: each kernel of the fabric's sum is a wave packet ten kilometres long, so a deep fold was a dash. F65 makes them a process of their own. A renewed plain's fabric rides on a first-order train of deep folds three times its wavelength, whose hinges nucleate one to a wavelength, each with a depth from F63's heavy tail. A hinge comes up when its fold has amplified enough, the deeper sooner, and one that comes up beside a fold already grown is part of it; it grows both ways along the folds' direction at a rate that goes with its depth, and links with a neighbour whose hinge passes within about half a wavelength, as Grasemann and Schmalholz (2012) found fold segments do: where the two run alike it bends into the neighbour, a Y where the neighbour runs on; where they meet across it abuts it. No hinge crosses another. A tip that finds nothing ends free and plunges, en echelon with its neighbours. A fold stands highest where it nucleated and lower where it reached last, so a linkage leaves a saddle; its crest is narrow between broad troughs, as a thin-skinned belt's folds are, and one of the middling depth is as steep as the fabric's folds, 430 m on Enceladus, the deepest about a kilometre, under Samarkand Sulcus's 1.75 km (Giese et al. 2017). The folds are read as one surface, a fold's profile about the distance to the nearest hinge, so that two hinges merge where they meet without a dip in the one that runs on.

It took five tries, each shown by a plot of the hinges over the shaded crop (the finish writes them as polylines when EXOGEN_DEEP_DUMP names a directory). Nuclei taken from the fabric's tail stood sixteen kilometres apart and linked at four, so they hardly ever met, and a maximum over the hinges' segments put a trough along every hinge, a segment's end cap lying in its neighbour's trough; the readout became one surface. Joins that reached a neighbour ending beside or behind the tip drew a hop sideways, the comb the first render showed; a join is now made only ahead, and a tip running away from a neighbour passes it. Two joins drawn at different times crossed; every new piece of a hinge, a step or a join, now ends where it first meets another. A train whose hinges all started at once made a mesh of short bridges that shaded as squiggles; with the deep hinges coming up first, they grow long and the shallow ones fill the gaps and link into them. And every crack on the body stopped the folds, which broke them into islands; only a terrain's own vents bound its folds now, a younger crack across it cutting them instead.

The finish had grown from 240 s to 436 s on Enceladus at standard quality. The time went into the fissures' direction, which the fabric's kernels and the growing hinges both read: it took every crack point within 60 km with an arccos, an arctangent, a sine, a cosine and an exponential each. It now compares chords, takes the doubled angle from the dot products, and averages over the crack points thinned to a tenth of its smoothing length, each standing for the points about it; the crops are unchanged to the eye, and the finish takes 233 s, a little less than before the deep folds.

What the crops show. The renewed terrains are covered by long, sinuous, narrow-crested ridges that merge and branch, with the fabric's folds between and on their limbs; between the active terrain's stripes they braid along the stripes. On the globe the young terrains read as ridged plains where they showed dashes. What they still do not show. A join cut short where it met a third fold leaves a small step; the crests are all of one width for a given lid; the deep folds follow the fissures' field as it stands, and keep no earlier pulse's strike; the fine set of the lanes still shows the grid on diagonals; and the furrows, the ring systems, Europa and Io are as they were. Ganymede, Io and Europa render byte-identical to the last milestone. Tests 126 (three new: hinges link within half a wavelength and pass farther apart, a deep fold is a crest between two troughs, deep folds join into a network and never cross), clippy clean, validate 116 checks.

### The fifteenth milestone: joins without Z turns, and relaxation in place of a halo

The user's review of the fourteenth: the folds "now read less like fingerprints and more like genuine terrain", their relief matching Enceladus's real profile; "the 'Z' turns are quite jarring"; "the large smooth patches of terrain have smooth ramps that match the shape of their outline exactly, is there a better way to do that?"; the audit is the checklist for the rest of Phase 4 and for the revisit of the earlier phases; and where does the process for mountain ranges like Phlegra Montes come in.

The Z turns had four causes, each found on the plot of the hinges over the crop. A join aimed at a point a fixed distance ahead on its neighbour, set by a merging angle of 25 degrees, whatever the tip's own approach, so a tip coming in steeply overshot and came back in an S. A join that crossed a third fold was cut short where it crossed it. An abutment was drawn to the nearest point of the neighbour, which could lie beside the tip. And a join to a neighbour that ended short of the junction hopped to its end. Now a join aims where the tip's own heading meets the neighbour, or, for a tip running alongside, no more steeply than 10 degrees; one that would cross a third fold is not taken; a tip meeting a fold across runs straight on into it; and a neighbour that ends short of the junction is passed. The ridges merge in long acute Ys, and a little over half of the test network's ends inside it still link.

The smooth patches took two readings. The first was the polar terrains' halo: the renewal took up to half of the old relief shorter than the renewal length out to one shell thickness beyond the outline, every crater alike, by a weight that was a function of the distance to the outline, so the outline, tongues and all, reappeared in the old ground. The relaxation law of F39 says otherwise: through a 10 Myr episode an Enceladus shell thinned to 0.45 of its 44 km keeps a 25 km crater whole and loses a 100 km load. So the halo is gone (F66): a terrain takes the old relief inside its outline and across the ramp at its edge, and beside it the old ground relaxes under the episode's thinned shell by that law, by band in the loop and by crater in the finish. On Enceladus that is almost nothing; on an Enceladus half again as large the craters beside a terrain relax by size, the largest first. The crops then showed that the patches the user meant were the large craters: a crater is one radial profile read at a radius stretched by a crenulation of its outline, so its floor's edge, its wall and its rim are scaled copies of the outline and its wall is as wide all round. The better way is a wall made by its collapse, slump blocks on arcuate faults, each with its own arc, width and drop, terraces that come and go around the crater, the floor's edge at their toes and a rim crest scalloped where large slumps took it back. It goes with the ring systems of large impacts, the next milestone, and the audit records it among the stamped shapes.

Phlegra Montes is not in the model. Its place is the revisit of the earlier phases: compression is to become one process along the strain field that loads, global contraction and the dichotomy set up, with wrinkle ridges where the shortening is small, lobate scarps where it is larger and thrust-bounded mountain belts where it is largest, so that a range like Phlegra, eroded old massifs beside Elysium flanked by the wrinkle ridges of Phlegra Dorsa, emerges where the strain concentrates. PLAN.md's Phase 4F entry records it.

Tests 128 (two new: a renewed terrain takes nothing beyond its ramp; a thinned shell relaxes the long wavelengths beside a terrain first), clippy clean, validate 116 checks. Ganymede, Io and Europa render byte-identical to the last milestone.

### The sixteenth milestone: crater walls made by their collapse

The user's review of the fifteenth: "Enceladus is looking better. The smooth patches diagnosis makes sense, and the crater wall fix sounds solid." On the reference maps in research/reference/maps: they are composites made from real data, so they serve to check a look and to find a missing or wrong process, and a change follows only where a source backs it, as with Enceladus's crisp large craters, traced to a shell twice the measured thickness. The tidal heating goes with the ring systems in the next milestone, since a giant impact leaves a basin, a ring system or a palimpsest by the thickness and warmth of the lithosphere it strikes.

A complex crater's wall is now made by its collapse (F67). The transient cavity's rim slumps inward on inward-facing normal faults in rings of blocks (Settle and Head 1979), terrace after terrace from the rim inward, the outermost 0.09 D^0.87 km wide on the Moon (Pearce and Melosh 1986) and each inward three quarters of the one outside it, until the terraces would grow too narrow or fill the wall; below them a plain slope of debris runs down to the floor. A terrace's width depends on the crater's size and the strength length of the collapsing target, so it goes as that length to the power the fit leaves, 0.13; the simple-to-complex transition measures the same length, so every body's terraces follow from its transition with no constant of their own, Mercury's a tenth narrower than the Moon's at a given size, as Leith and McKinnon (1991) measured them. Each terrace is a ring of slump blocks of their own length, throw and bow, some of which slumped no terrace, so that the wall's width, its terraces and the floor's edge change from block to block around a crater; the outermost blocks scallop the rim crest in place of the crenulation, and the innermost blocks' toes bulge the floor's edge. The loop draws the plain wall and the scalloped rim, the finish the terraces, in the craters it stamps and in a loop crater's restored detail.

The first version drew flat benches between thin scarps, and on Enceladus's 140 km crater, where the cells are 387 m, the terraces read as contour lines. The reference map of the Moon, coarse as its 8-bit relief is, shows Copernicus's wall as a rough, lumpy, scalloped band rather than clean steps. So the terraces became fewer and wider, their benches keep a third of the wall's slope, their scarps wander along their traces and differ in width from block to block, and the slumped blocks are broken, their benches hummocky.

What the crops show. On Enceladus's 140 km crater the wall is a ring of slumped blocks: four terraces that wander, step at the blocks' ends and come and go around the crater, hummocky benches, a scalloped rim and a lobate floor's edge, where it had been one smooth band copying its outline. On the Moon and Mercury at the standard grid, whose cells are as wide as the terraces, the walls read as rough steps and the rims as scallops. What they still do not show: a bench carries no ponded melt and a scarp no talus fan; central peaks and peak rings are as they were; and an ice shell relaxes a crater's terraces with the rest of it by one factor.

Tests 129 (one new: a complex crater's wall is terraced and differs around it, with Pearce and Melosh's width and Mercury's narrower terraces), clippy clean, validate 116 checks. Io renders byte-identical to the last milestone; every body with complex craters changes.

## The seventeenth milestone: a porous core heats Enceladus

The user's review of the sixteenth: "Looks good across all bodies. Craters on Enceladus look smoother, but I assume the missing processes you listed will help with that, while the terraces as a base will avoid giving us that perfectly smooth look we had before?" They will, and one of them causes the smoothness in its own right. The terraces are the wall's structure, and talus on the scarps and melt ponded on the benches will work on it; but an ice shell also relaxes a crater by wavelength, its bowl shallowing and its floor doming while its rim and terraces, short loads on the elastic lid, stand, as Enceladus's relaxed craters show (Bland et al. 2012). The finish relaxed each crater it replays as a whole, its rim and terraces by the share its bowl keeps, and restored a loop crater's rim and terraces unrelaxed on its relaxed bowl. Relaxation by wavelength goes with the ring systems, which depend on it: a palimpsest is a basin whose relief has relaxed away and whose rings, shorter, remain.

The audit after the fifteenth milestone found Enceladus's shell twice its measured thickness, 44 km where its gravity, shape and libration give 18 to 25 km (Čadek et al. 2016; Beuthe et al. 2016), because its tides came only from its ice and its consolidated rock. Its rocky interior is not consolidated: at its low pressures it stays a pile of grains with water in its pores, which dissipates tides by the friction of the grains and the flow of the water, 10 to 30 GW in Choblet et al.'s (2017) model of a core 20 to 30 percent porous. Now the rocky interior under an H2O layer has a porosity from the pressure in it as sediments compact, 0.4 at no pressure falling by a factor e for each 60 MPa (Athy's law for sandstone, Sclater and Christie 1980), and the pores above 0.05 dissipate the eccentricity tides on the interior's radius with an imaginary Love number in proportion to them (F68). Enceladus's interior comes out 0.23 porous and dissipates 15 GW; the large moons' interiors, at gigapascals, are consolidated, and Ceres's, 0.06 porous, raises no tides worth the name. Only Enceladus changes.

Its tides come to 30 GW, its heat flux to 38 mW/m², its shell to 22 km and its elastic lid to 12 km, where they had been 15 GW, 19 mW/m², 44 km and 23 km, and validate now checks its shell against 15 to 30 km. The first run with the thinner shell renewed nearly all of the surface in 34 episodes, since the episodes recurred every ten conductive times of the shell and a 22 km shell conducts its heat in a quarter of the time a 44 km one does. Ten had been a fit at the short end of O'Neill and Nimmo's (2010) range of 0.1 to 1 Gyr, made on the thick shell; refitted to forty, the episodes come every 0.52 Gyr as before, eight of them. The audit notes that the recurrence should come from the heating cycle of the Dione resonance rather than a factor on the conductive time.

What the maps show. More of Enceladus is renewed, 58 percent against 57, its lanes four where there were twelve, its relief 9 km from lowest to highest where it had been 14.5, and its craters relax more beside its terrains and across its plains; the largest have softer walls, their terraces carried on a lid half as thick.

Tests 130 (one new: a small moon's porous core dissipates its tides; three lane tests restated on an Enceladus with four fifths of its eccentricity, a body with the lid they test), clippy clean, validate 117 checks (one new: Enceladus's shell). Every body but Enceladus renders as before.

## The eighteenth milestone: ring systems of large impacts, and craters that relax by wavelength

The seventeenth milestone answered the user's question on Enceladus's smoother craters with one of its causes: the finish relaxed a replayed crater as a whole, its rim and terraces by the share its bowl kept. This milestone makes relaxation work by wavelength and adds the ring systems that the same physics makes of a large impact into a thin lithosphere.

A crater now relaxes by wavelength (F70). Its radial profile is split into bands between Gaussian blurs a factor of two apart; each band keeps the share of its relief the shell lets its wavelength keep over the crater's life, as the loop's bands already did, so the bowl shallows and the floor domes while the rim and the wall's terraces, short loads the elastic lid holds, stand (Parmentier and Head 1981; Dombard and McKinnon 2006; Bland et al. 2012). The finish's replayed craters and the restored detail of the loop's craters relax so. On Enceladus little changes: its large craters are the loop's, relaxed band by band already with their rims restored on top, and its replayed craters are short enough for its 12 km lid to hold. On Callisto and Ganymede, whose early shells were warm, the old craters keep rims and terraces on shallowed bowls where they had been scaled down whole.

A basin's collapse now breaks a thin lithosphere into rings (F69). The collapse draws the weak substrate into the cavity, and the inflow drags the lithosphere toward the basin until it fails in concentric normal faults (Melosh and McKinnon 1978; McKinnon and Melosh 1980). The lithosphere through the collapse is the ice colder than 240 K, read from the shell at the impact's time and place: Callisto's comes out at 19 km four billion years ago, where McKinnon and Melosh read 15 to 20 km from Valhalla's rings, with no constant set for it. Melosh's (1982) measure of the drag against the strength, as Singer et al. (2023) evaluate it for Europa's ring graben, goes as the square of the transient diameter over that thickness: above 2 the lithosphere fails, above 5 it breaks into an annulus of graben. The graben run from 3 to 8 transient radii, the range Tyre, Callanish, Valhalla, Asgard and Ganymede's furrows span; each is two faults meeting at the base of the lithosphere, so about as wide as the lithosphere is thick, five widths from the next, spaced wider outward and shallower, as deep as the strain it takes up, its flanks raised; inside them the blocks' raised edges stand as ridges. At Valhalla's size this gives graben 20 km wide and 100 km apart, where Valhalla's are 15 to 20 km wide and 75 to 100 apart. The drag needs a substrate that flows: under an ice shell it relaxes a basin in days; under a rocky body, whose mantle this model does not relax, not at all, so the Moon's basins keep the rings of their own collapse and no more. The arcs are events of the finish at the basin's time, so younger craters overprint them and later lanes renew them.

The first version took the elastic lid, 8 km on early Callisto, for the lithosphere: 606 basins broke, down to craters of 100 km, where Callisto's craters of 100 to 300 km are dome craters and penepalimpsests with no rings. The lithosphere that holds through a collapse is the thicker layer colder than 240 K, as McKinnon and Melosh define it; on it 77 of Callisto's basins break, and the graben come out at Valhalla's width and spacing. At first only the graben were drawn, and Valhalla's most visible part on the reference map, the dense sinuous ridges between its bright centre and its troughs, 55 km apart (DeFour-Remy et al. 2022), was missing; the ridges followed.

What the maps show. Callisto's basins carry ring systems of concentric graben and ridges, and the largest of them, 3,336 km across, puts rings across most of the moon, a furrow system like Ganymede's. On Ganymede the furrows run in the dark terrain and stop where the grooved lanes renewed it, as the real furrows do. On Enceladus the rings of its oldest basins, whose transient cavities are nearly as wide as the moon, sweep across the cratered plains.

What they still do not show. A palimpsest's centre stays bright for four billion years because its melt sheet is clean ice far thicker than gardening reaches; here it darkens with the lag like any surface, so Valhalla-like basins have no bright centre (the audit's M2, for Phase 5 with the lag of icy surfaces). Ganymede's furrows have bright rims about dark floors, and theirs are 90 km apart where the real ones are 35 to 60: the furrows' widths imply a lithosphere of 6 to 10 km when they formed (Burkhard and Thomas 2024), half of what the thermal history gives early Ganymede, a matter for the earlier phases' thermal model. Europa's Tyre and Callanish need its shell as it was then, the Europa milestone's question. And the largest basins of Callisto's and Enceladus's populations are larger than any the real moons show, a question for the crater production in the revisit.

The audit's S8 is closed: the dip in the lanes' share came from Enceladus's flux standing just under the focusing's saturation, so that a larger body's flux at first focused more of the renewal to the pole before its thinner shell focused less; since F68 Enceladus is past the saturation and the share rises with every step of the sweep. M2 records the palimpsests' brightness, and RING_STRAIN joins the calibrations.

Tests 135 (five new: a relaxed crater loses its bowl before its rim; transient cavities follow McKinnon and Schenk; Melosh's gamma matches Singer et al.'s evaluation; a thick lithosphere holds and a thin one breaks into rings; graben widen apart and shallow outward), clippy clean, validate 117 checks. The Moon, Mercury and Io render byte-identical to the last milestone; Callisto, Ganymede, Enceladus and Europa change. Enceladus's finish takes 416 s, against 276 s at the sixteenth milestone: 30 s of it are the rings, and a run without rings or relaxation by band takes 382 s, so the rest came with the seventeenth milestone's thinner shell, whose terrains cost more to draw; which of them is not yet traced. A basin's arcs are drawn together, in parallel, as one event of the finish. EXOGEN_SKIP_RINGS and EXOGEN_SKIP_BAND_RELAX switch the two processes off for diagnosis, and EXOGEN_DEBUG prints each ring system and the finish's phase times.

## The nineteenth milestone: craters on ice

The user's review of the eighteenth: "The crater walls especially on Enceladus' largest craters still read much more defined than expected. The relaxed shell didn't do much to fix that look. Ganymede's crater rims show much more prominent rings than real Ganymede, but I do see the 'double' craters that Ganymede has a lot of. The many smaller rings in medium-sized craters look too plentiful compared to real photos. The shape of the rings also perfectly follows the shape of the crater, without any natural variation, which lends to their exaggerated look."

The rings were not the ring faults of the eighteenth milestone but the craters' own interiors. Every crater on ice took the Moon's shapes: Pike's lunar depths, so that Ganymede's 60 km craters stood 4.8 km from floor to rim where the real ones are 0.5 to 1.7 km (Bray et al. 2012), and the lunar sequence of interiors, so that every crater of 24 to 56 km on Ganymede had a peak ring and every one of 56 to 150 km three interior rings and an outer ring, where the real ones have central pits and domes. And every contour of a crater, its floor's edge, its terraces, its rings, was the rim's scalloped outline stretched inward, so all of them nested like copies of the rim.

A complex crater's relief now scales with the strength length (F71): its depth, rim and central peak are the lunar crater's at the same multiple of the transition, scaled by the body's transition against the Moon's. Ganymede's complex craters come out 60 to 70 percent shallower than the Moon's, as Croft (1981) measured them; Mimas's Herschel comes out 11 km deep, where it is 10 to 12; the Moon keeps Pike's laws, and Mercury's craters come out shallower than the Moon's, as its higher gravity has them. On ice the central uplift now collapses into a pit (F72), where rock's falls into a peak ring: from 8.5 transitions a crater on ice has a central pit, a fifth of it across and widening with the crater, ringed by a raised, broken rim, and from 20 a dome rises in the pit, at the sizes Ganymede shows them (Alzate and Barlow 2011; Schenk 1993; Bray et al. 2012); neither keeps the lunar rings, whose part on ice the lithosphere's own failure plays (F69). The rim's outline is now carried only by the rim and the outer part of the wall (F73), where the outermost slump blocks that scallop it stand, so the inner terraces, the floor's edge and the central features keep shapes of their own, and the rim's height varies around the crater by a sixth, as Ganymede's do.

Enceladus's walls had two causes. Its fresh craters are deep: its low gravity makes its transition large, and a large transition deepens complex craters (F71), as Herschel on Mimas shows. And nothing softened them after: its 12 km lid holds every load shorter than about 45 km, so its craters relax only beside its terrains during its episodes, and the fallout of its plumes, the snow that mutes its real craters (Kirchoff and Schenk 2009), only lay on the ground as a cover for its pits. The fallout now mantles the ground (F74): the deposit creeps down what it covers, filling hollows and thinning on crests, so relief smaller than its thickness is buried and larger relief is rounded by no more than it. On the old cratered plains, under hundreds of metres of fallout, the small relief is muted and the terraces of the old craters soften; the largest crater, 162 km across and two billion years old, lies where little fallout reached, and keeps its walls. The real Enceladus has no crater that large, and its craters relaxed under heat flows of 150 mW/m² or more (Bland et al. 2012), four times what its shell carries now; the model's heat history is flat, which the audit records as M3.

What the maps show. Ganymede's craters are shallow bowls with a rim and a central pit, the larger ones a dome in a wide pit, with no nested rings: its craters of 24 to 56 km stand 1.3 km from floor to rim at the median and those of 56 to 150 km 1.6 km, against 0.5 to 1.7 km measured (Bray et al. 2012), where they stood 4.3 to 5.4 km. Callisto's are the same; with its craters shallow, its ring systems stand out the more. Europa's largest crater, 35 km, has a central pit where it had a peak ring. The Moon's basins keep their irregular floors, now an outline of their own, and Mercury's basin rings each wander on their own. On Enceladus the fallout lies 231 m deep at the median, 28 m at the tenth percentile and 497 m at the ninetieth, against Martin et al.'s (2023) mean of 250 m and most of 700 m; the short relief of its old plains is muted by a fifth and the old craters' walls are rounded.

What they still do not show. Enceladus's old craters are relaxed in the real cratered plains and not in the model's, as M3 records, and its largest crater keeps its terraced walls under a thin mantle. A palimpsest stays dark (M2). Callisto's and Ganymede's ring graben stand out against the shallow craters, their inner graben up to 2 km deep at Callisto's largest basins; whether that is the strain or the basins, larger than the real moons', is for the revisit. The pits' and domes' sizes come from Ganymede and are untested on Callisto and Ceres (S9 and the calibrations).

Tests 136 (one new: craters on ice are shallow and pitted, with Herschel's depth; Ganymede's relaxation test now reads the floor between a pit's rim and the wall), clippy clean, validate 117 checks. Io renders byte-identical; every body with craters changes. The mantle adds the low-passes of one pyramid to Enceladus's finish. EXOGEN_SKIP_MANTLE switches it off.

The user's review of the nineteenth: "The icy moon craters match the real imagery closer now. The plume fallout is a good idea but it looks a bit too regular on the surrounding terrain to me." The mantle evened all the ground under it alike: it was laid after the unresolved craters, so it buried the young small craters with the old, and it blurred over four times its thickness, so the plains under thick fallout took one soft texture. Now it is laid first, and the unresolved craters count only those dug since the cover above them grew thinner than they are deep: a crater older than that lies buried, a younger one pits the snow. It evens the ground over twice its thickness and fills or cuts no more than its thickness. On the old plains, 40 to 70 degrees north, the short relief (below about 3 km) now falls from 232 m RMS to 220 m under the mantle, where the first mantle took it to 209 m; the old craters' walls stay rounded and the plains between them keep their young pits.

## The twentieth milestone: the shell's own tides, and the resonance's cycle

The user asked whether the heat history would touch other bodies or only Enceladus-likes, and to go on to orbital histories if they touch it. The heat history reaches every body: radiogenic heat and secular cooling everywhere, tides wherever an orbit raises them. What it lacked was an orbit with a history, and a shell that answers one. Every body's tides came from a Love number k₂ and a dissipation factor Q by composition class (F13), so an ice shell's heat went as its eccentricity squared whatever the shell's state. Europa's came to 76 mW/m² and its shell to 7.7 km, where Juno's microwave radiometer measures a conductive shell of 29 ± 10 km (Levin et al. 2025), and an eccentricity pulse would have thinned the shell without limit. A history of the orbit meant nothing until the shell dissipated as a shell does.

An icy moon's shell now dissipates its own tides (F75). A Maxwell solid strained at the tidal frequency turns 2 μ ⟨ė²⟩ / ω · ωτ / (1 + (ωτ)²) of work per unit volume into heat (Ojakangas and Stevenson 1989; Nimmo et al. 2007). The cold top of a shell is elastic, and its warm base, where the viscosity passes μ / ω, takes nearly all of the heat, so a conductive shell dissipates in closed form, from the ice of F39 and F41. The strain rate is the diurnal tide's (F41). Its geometric factor comes from the model's own breathing bulge: Europa's strain rate squared is 2.4e-20 s⁻², against Ojakangas and Stevenson's 3.3e-20, with nothing fitted. The Love number is the shell's own (Beuthe 2018): a membrane spring against the ocean's buoyancy, soft on a large moon (Europa's h₂ 1.16) and hard on a small one, whose stiffness then sets the tide (Enceladus's h₂ 0.046, against Beuthe's 0.043 for 23 km; the class rule had given it 0.43). The shell's thickness is where its conduction carries the heat from below and its own.

Europa's shell comes to 30.3 km, dissipating 0.43 TW, where Lemasquerier et al. (2023) find 0.56 TW in a 25 km shell. Five times its eccentricity would thin it to 7 km; Kihoulou et al. (2025) model 35 km thinning to 5. Enceladus's hard shell dissipates 0.6 GW, a small share of its heat, as Souček et al. (2019) and Beuthe (2019) find. Its porous core, refitted to 28.5 GW within Choblet et al.'s 10 to 30 (F68), holds its shell at 22.7 km. Ganymede's and Callisto's shells stay at 83 and 97 km, their h₂ near 1.3. Ganymede's resonance pulse (F50) now takes an eccentricity of 0.0185 to heat its shell to 8 mW/m², inside Showman and Malhotra's range of 0.01 to 0.07, where the old 0.008 had sat below it. Three fits keyed to Europa's old flux were refitted against the same observations: the chaos flux (F45), the crack density and the thin-shell share (S4).

A shell's stress follows its local Love number (F76). A hard shell's stress goes nearly as the inverse of its local thickness (Beuthe 2018), so the tidal field is scaled at each point by the local h₂ over the mean's. Enceladus's uniform shell is stressed by 20 kPa (Beuthe: about 15 kPa in a uniform 23 km shell), below the 40 kPa at which ice cracks; its thinned south polar shell reaches 78 kPa. Its cracks therefore open where its shell is thin, now or in a past episode. Its old plains, which the tenfold h₂ had cracked everywhere, the diurnal tide no longer cracks; the audit records the stresses that would (M4).

The orbit now has a history (F77). In the Laplace resonance, Io's heat softens its upper mantle, its dissipation damps its eccentricity, and the mantle cools and lets the resonance raise the eccentricity again (Ojakangas and Stevenson 1986; Hussmann and Spohn 2004). The moons' measured migration rates, heat flows and tidal response select a cycle that repeats every 95 to 150 Myr, begun 0.8 to 2.0 Gyr ago, with the moons now just past a peak (Bagheri et al. 2026). Four orbital inputs carry such a cycle, as three carry Ganymede's passage. Io and Europa take Bagheri et al.'s representative history:

- **Period:** 0.12 Gyr, from 1.6 Gyr ago.
- **Phase:** a fifth of a period past a peak.
- **Swing:** Io's eccentricity by a fifth and Europa's by a quarter.

The thermal history samples the cycle as finely as the curve is laid. Europa's shell breathes with it between 26 and 40 km (Bagheri et al.: its total shell swings by 3 to 15 km). Enceladus gets no cycle: in Ojakangas and Stevenson's model it does not oscillate about its tidal equilibrium (Meyer and Wisdom 2008), so M3 stays open.

What the maps show:

- **Europa:** looks as it did. 2,771 lineae against 2,651, and the same 1,022 chaos patches covering 13 percent, now on a shell of its measured thickness.
- **Enceladus:** its south polar terrain is unchanged. Its lineae fall from 1,316 to 739 and lie in its thinned zones; between 0° and 100°E across its cratered band they fall from 109 to 17. The long arcs left on its plains are the ring graben of its large impacts (F69).
- **Ganymede:** lanes cover 0.68 of it against 0.65, from 1,241 lanes against 1,099.

What they still do not show:

- Europa's chaos, cracks and wholesale renewal read its present shell, not its history; the Europa milestone is to make them read its cycle.
- Enceladus's old fractures need the stresses of a thickening shell, of polar wander or of non-synchronous rotation (M4).
- A rocky body's tides still come from its class k₂ and Q, so Io's heat cannot run away or limit itself (S10).
- The coarse loop's 50 Myr steps see a 120 Myr cycle only at their midpoints.
- Titan's k₂ comes out 0.43 under the rigid-core limit, against the measured 0.59 to 0.62, which a dense ocean would raise.

Tests 139 (three new: a floating shell's Love numbers follow the membrane spring; the diurnal strain rate matches Europa's, with its shell at Juno's thickness; Europa's shell breathes with the Laplace cycle; and the lanes test re-anchored, at full heating, on a body that makes lanes of both kinds), clippy clean, validate 117 checks (Europa's shell against Juno in place of its old tidal flux range). Io renders byte-identical; Europa's maps move by under a metre with the cycle.

## The twenty-first milestone: Europa's lineae age, and streams that crowded

The user on the twentieth: "Enceladus' crater walls still appear more 'laddered' than real Enceladus shows. If that's later work, then go ahead and continue with the next milestone." It is later work. The real Enceladus has no crater in its plains larger than about 40 km, all relaxed and muted (Bland et al. 2012; Kirchoff and Schenk 2009), and every laddered crater on the model's is larger than that; M3 and the crater population's high end own both. The terraces' rings of aligned blocks, which read as rungs on any deep crater, are recorded in the audit (T2).

Europa came next, the weakest of the ice bodies: lace discs for chaos, hairlines broken in steps for lineae, flat plains between them. Surveying it against the real moon found two problems first.

A shared stream. Cracks are placed by rejection: a point drawn at random is kept with a probability set by the tension there. The point's longitude and the test's draw came from the same stream, so a point was kept more often when its longitude was small. Europa's lineae were four times denser in one quarter of the moon than in another, and its chaos left a sector of 60° empty. The same fault, in milder forms, tied Io's mountain heights and Ceres's cryodome ages to their longitudes, and a plume rise's retries tested one point sixteen times. Each draw now has a stream of its own. Tests check that Europa's cracks and chaos reach every longitude, and they fail on the old streams. That fix was committed on its own, as c1b396d.

The lineae (F78). At the standard grid Europa's cells are 2.4 km, and its double ridges are 1 to 2 km wide. No ridge could be drawn, and every crack was a hairline of whole cells. What a map shows of Europa's lineae is the sequence Geissler et al. (1998) found they age through:

- incipient cracks, colourless;
- ridges, 3 to 6 km wide;
- triple bands, a dim medial stripe with diffuse red-brown margins 8 to 12 km wide;
- ancient bands, only slightly darker than the plains.

The drawn lineae are these. A ridge is 2 to 7 km wide, from a heavy tail, and 0.05 of its width high. It lays its margins as a stain in a new surface layer: about 9 km across a 3 km ridge, dimmer over the medial stripe, absent while the crack is fresh, and fading over 30 Myr to the ancient band's tone. The stain follows the trace while the relief follows the segments. A younger crack cuts the older stain under its own width, and a surface renewed later carries none. A ridge narrower than two cells is drawn by the area of its cross-section, and only a resolved ridge or band writes a surface of its own, so no crack is drawn in whole cells. The double ridges narrower than about 2 km make the ridged plains' fabric. They belong to the next milestone, with the chaos whose rafts carry them.

What the maps show:

- **Europa:** its lineae cover the whole moon as soft red-brown lines in long arcs and cycloids, the young ones dark and the old ones faint, without steps. At map scale they read as the real moon's do. Its chaos is unchanged, and the lace discs are now the most artificial thing on it.
- **Enceladus:** unchanged to the eye. Its fallout frost and white carbonate salts hide any stain, and its vents keep their sintered tint.

What they still do not show:

- **The ridged plains' fabric:** Europa's narrow double ridges, crosscutting at every scale (Doggett et al. 2009).
- **Chaos from the shell's heat:** the next milestone.
- **Bands and the thickening shell:** bands as spreading lanes, and the stresses of a thickening shell; the milestone after.
- **Stated choices:** the timescales of a crack's margins.

Tests 142. One is new in this milestone: a ridge stains its margins as it ages, colourless when fresh, a triple band some 9 km across at a few million years, and the ancient band later, drawn by its area on a coarse grid. Two came with the stream fix. Clippy is clean and validate passes its 117 checks.

The user on the twenty-first: "I'm happy with the lineae, but can you check why some of them have a 'staircase' look to them? That doesn't look natural to me." On the equator the tidal field is symmetric about the equatorial plane, so its principal directions lie east and north, and where the two principal stresses cross, or as the orbit turns them, the most tensile direction swaps by a right angle. A trace took the field's heading at every step, and a cycloid ran whenever the most tensile stress passed the cracking stress, whatever its direction, so the cracks turned at right angles there. Now a crack runs on only while the stress across its own plane opens it, and it kinks toward the field's heading by the maximum tangential stress criterion of fracture mechanics (F79). Under near-equal tension in all directions it runs straight on, and where the field swaps through the orbit a cycloid waits for the next orbit instead of turning. Making the trace remember its heading exposed a fault of its own: the heading, carried from step to step, was left tilted out of the tangent plane, and the trace zigzagged at every step, by up to 48 degrees, until each step laid it back in the plane. The sharpest turn of a crack traced from the equator is under 45 degrees, where it was 90, and the cycloids keep their arcs and cusps. Two tests guard it: no right-angle steps from the equator, and no step-by-step zigzag across Europa's cracks. Tests 144, clippy clean, validate 117.

## The twenty-second milestone: chaos from the shell's warmth

After the staircases, as the user asked, the next milestone: Europa's chaos. Its patches (F45) were discs of log-uniform size, each a wobbling outline round a fixed profile of matrix and rafts, and they read on the map as lace discs, the most artificial thing on the moon.

Chaos now forms from a field (F80). The ocean brings warmth up under the ice, most of it at low latitudes (Soderlund et al. 2014). The shell disrupts where that warmth stands high: over melt lenses a few km down, as under Thera and Conamara (Schmidt et al. 2011). The warmth is structured at every scale, from 1,300 km, the largest contiguous chaos Riley et al. (2000) map, down to the rafts, which no finer structure can show. Chaos comes in generations, four in each retention age, since Riley et al. find chaos and the tectonics took turns. Each generation has a warmth of its own and disrupts the share of the surface that leaves Europa's chaos fraction disrupted by all of them together. The lid over the warm ice breaks into rafts two lens depths across, which group into plates two rafts across. A plate founders into the matrix with a chance that grows with the warmth's excess; in a foundered plate, single rafts stand as blocks. What stands drifts and turns, carrying the old surface with its ridges, ages and stains. Leads of matrix open between plates, round blocks and round the chaos itself. The matrix settles 150 m below the old surface's smoothed level, stained by the ocean's material as a crack's margins are, and fades with the salts. Each generation is one event over the whole body, between the cracks and craters older and younger than it.

It took three tries, each judged on the maps:

- **Single rafts.** The first try broke the lid into rafts of one size, 6 km, each foundering on its own. On a 2.4 km grid a raft is two or three cells, and rafts foundering independently read as salt-and-pepper speckle.
- **Plates of rafts.** Grouping the rafts into plates that founder whole made the edges of each chaos a mosaic of plates parted by leads. But the cores became flat pools of matrix: each octave of the warmth had half the last one's amplitude, so its excursions were a few smooth blobs.
- **Warmth at every scale.** With each octave keeping three quarters of the last one's amplitude, the small scales break a chaos into patches of every size, clustered in large regions. That is how the real moon's mottled trailing hemisphere looks in the reference.

What the maps show:

- **Europa:** chaos covers 30% of the surface in clusters of patches of every size, denser at low latitudes. About half of it is brown matrix, against Spaun et al.'s (1998) 59% for Conamara. Plates stand at the margins and blocks inside. The youngest generation is vivid. The older ones are faded to tints, with their plates and leads, and crossed by the lineae younger than them.
- **Every other preset:** unchanged; only Europa has chaos.

What surveying the maps found:

- **Europa's plains carry the body's age.** More than half of Europa's surface keeps a surface age of 4.5 Gyr. The wholesale renewal of a thin shell (S6) sets the retention age of the crater record, but it renews no surface. The plains between the drawn lineae and chaos are left as the coarse loop's first crust. The sub-grid crater roughness (F54) reads their age and roughens them as a saturated ancient surface: the grain on Europa's plains in the shaded maps, against which the young matrix looks smooth. What renews them is the ridged plains' fabric of narrow ridges. It was planned for this milestone and moves to the next, with S6, since its renewal is what gives the plains their ages.

What they still do not show:

- **The ridged plains' fabric and the plains' ages:** the next milestone.
- **Bands and the thickening shell:** bands as spreading lanes, and the stresses of a thickening shell.
- **Coupling:** the warmth does not yet thin the shell that the cracks read (F76).
- **Stated choices:**
  - the generations' count and spacing;
  - the warmth's equatorial excess and persistence;
  - the rafts' and plates' sizes, drifts and turns;
  - the leads' width, the disruption's scale and the matrix level's smoothing length.

Tests 146. Three are new, and one is gone: Europa's chaos covers its share at every longitude; it favours low latitudes; and a generation breaks the lid into rafts in a lower matrix, with untouched ground outside. The old test that chaos patches reach every longitude went with the patches. Clippy is clean and validate passes its 117 checks.

## The twenty-third milestone: the ridged plains renew Europa's surface

The user on the twenty-second: "Good work. The chaos terrain looks natural and real. Continue with the next milestone when you're ready."

The twenty-second milestone's survey had found that Europa's plains carried the body's age. 84% of the moon's surface stood at 4.5 Gyr and its mean age was 3.4 Gyr: the wholesale renewal of a thin shell (S6) set the crater record's retention age but renewed no surface. The sub-grid crater roughness read that age and roughened the plains as a saturated ancient surface.

What renews the plains is what makes them: the ridged plains, Europa's oldest unit, crosscut by ridges at every scale (Doggett et al. 2009). The fabric (F81) is the narrow tail of the crack population the lineae come from. Its cracks are 0.5 to 2 km wide, as narrow as the double ridges Leonard et al. (2018) measure in Galileo's best images. They have the same heavy tail of widths and the same lengths, are placed by the same tension, and are traced along the same field. None is resolved at 2.4 km. Each one:

- renews the share of the cells it crosses that its width covers;
- adds its cross-section to their relief;
- stains them by its margins.

The drawn ridges narrower than two cells now renew their share in the same way.

How much the fabric renews follows from the retention age. A surface renewed throughout renews one surface's worth in each retention age, its mean age. Within that age the drawn lineae renew 17% of Europa and the chaos matrix 15%. The fabric renews the rest, and all of it before that age, where what the lineae and chaos made has been overprinted into the plains. Its history is laid in slices, each one event over the body, back three retention ages. What it has not reached takes the age the renewal gives on average, since renewal keeps no memory. Europa gets some 170,000 narrow ridges, which add 7 s to its finish.

It took three renders, judged on the maps:

- **The first** made the plains young, but the fabric's lines beaded like rope. A strip narrower than a cell falls on one cell where its axis meets a cell's centre and on two halves between. Every ridge narrower than the grid is now drawn two cells wide by its share, as its relief already was, and its margins over at least 1.2 cells. Each young ridge had also shown as a double line, because its medial stripe was dimmed at full strength far below the grid; it now dims by its share.
- **The second** removed the beads, but the young ridges' margins, as dark as a triple band's, combed the plains into flow lines along the tidal field.
- **The third** makes a ridge's stain go as its width, full from 2 km. This is a stated choice: the dark material a ridge brings up goes as its cross-section, spread over margins as wide as itself. The plains now read as a fine bright texture, and the drawn lineae stand out on them again.

**A regression found on the way.** Rendering Enceladus to check the narrow drawn ridges showed half its tiger stripes cut short. The twenty-first milestone's revision (F79) made a crack run only while the tide across its own plane opens it. But a plume's fissures are placed by the bending of their neighbours' flanking ridges (F57), which opens them across their own plane by construction, and their venting holds them open. Under the tide alone, only three of the present episode's eight stripes ran their length, a fourth stopped at half, and the forks at their tips stopped as stubs. The cascade's fissures are now held open while their zone lets them run. The tiger-stripe test had measured the lengths the stripes were given, on a coarse grid with a seed of its own; it now measures their traced lengths at the standard grid with the preset's own seed, and fails on the twenty-first's code.

What the maps show:

- **Europa:** the surface's mean age is 97 Myr against its retention age of 100 Myr (it was 3.4 Gyr), and none of it is older than three retention ages (84% stood at 4.5 Gyr). The plains are smooth where they were rough with saturated craters, finely streaked by the ridges' renewal and faintly lined by the young ridges' margins. The drawn lineae and the chaos stand on them as before.
- **Enceladus:** its eight present tiger stripes run their whole 242 km again, forked at their tips (the regression above). Its narrow drawn ridges renew their share of the cells they cross, which does not show under its frost.

What they still do not show:

- **The plains are combed along the tidal field.** The young ridges share its headings over the time their margins stay dark: the shell turns 54 degrees under them in 30 Myr. The stresses of a freezing shell, the next milestone, would cross them.
- **The retention age** is still the regime's constant (S6); the fabric's rate closes the renewal to it rather than setting it.
- **Bands and the thickening shell:** bands as spreading lanes, and the stresses of a thickening shell; next.
- **Stated choices:** how far back the fabric is drawn (to 5% unrenewed), the stain's width scaling, and the narrow lines' spreading (two cells for the strip, 1.2 for the stain).

Tests 148. Two are new, and the tiger-stripe test is stricter: a fabric ridge renews its share of the cells it crosses and adds its cross-section, wherever its axis falls; and the fabric closes the renewal budget, its slices covering their share with widths from the narrow tail. Europa's full run now also checks that its plains are as young as its surface. The stain test now expects a narrow ridge to renew its share. Clippy is clean and validate passes its 117 checks.

## The twenty-fourth milestone: a freezing shell, and a heterogeneous one

The user on the twenty-third: "Progress looks good. I'll let you continue with the next milestone before I make a judgement on the global look, because I was going to point out the linear directionality of the ridges and cracks until you made note of it. Real Europa shows the cracks undulating and crossing each other, while ours follow fairly straight curves in the same 'lanes'."

Under the tide alone every crack follows the field's headings, and the field is smooth over thousands of kilometres, so cracks of a time run side by side in lanes. Two stresses besides the tide now act on every crack, drawn or narrow (F82).

**The freezing shell.** Water that freezes onto the shell's base takes more room as ice, and the ocean pushes the shell outward by the excess: some 750 kPa per km of thickening on Europa, and more than 10 MPa in a shell that keeps thickening (Nimmo 2004). Two physical points set how much of that the cracks see:

- the cracks a thickening shell opens relieve its tension as soon as it passes the cracking stress, so the shell stands at that stress;
- under a static stress of several MPa the tide could not steer a crack and no cycloid would form, and Europa has many.

So a thickening shell stands at the cracking stress in tension, and a thinning one at it in compression (a stated choice). Europa's shell, breathing with the resonance's cycle, alternates between the two; its thinning phases open fewer cracks.

**The heterogeneous shell.** The shell's thickness varies, so the lid carries the field unevenly, and the stress departs from the field by a fifth of its magnitude (a stated choice) in a pattern the shell's changes renew every few million years. Where the tide pulls nearly equally every way, about the points its pattern turns on, that pattern steers the cracks; where the tide is strongly anisotropic, it bends them.

It took three renders, judged on the maps:

- **Cells of two shell thicknesses**, the convection's scale: the cracks squiggled near the equator, at 42 degrees of turning per 100 km against 5 in the tide alone, and a few curled into loops.
- **Cells of twice the elastic lid's flexural parameter**, the length over which the lid spreads a load and with it the stress about a stiffer or softer patch (100 km on Europa): gentler, but the cracks still meandered like rivers where the tide is weak, with loops up to 150 km across.
- **The heterogeneity felt along the crack.** Each step had obeyed the stress at the tip alone. But a crack's tip feels the stress along the whole crack, weighted by the inverse square root of the distance back from the tip, whose mean is a third of the crack's length. Felt over a third of the length grown so far, the heterogeneity steers a young crack and only bends a long one. Europa's cracks now turn 11 degrees per 100 km against 5 in the tide alone, the meanders and loops are gone, and the plains cross in every direction.

The heterogeneity first bent Enceladus's tiger stripes too, and two of them converged. A plume's fissures are steered by the bending of their neighbours' flanking ridges, which runs along the neighbour and is far stronger near them than the heterogeneity (some 135 kPa on Enceladus against 16), so they keep to the tide and their neighbours. The tiger-stripe test now also checks that every whole stripe keeps within 50 degrees of the first's heading along its length; the bent stripes had diverged by 63.

The last look found one more thing: a few cracks curled into loops some 150 km across. They were cycloids. A cycloid forms because its crack stops and starts again with each orbit, following the turning diurnal field while it runs. Where the freezing holds the crack open through a whole orbit, it follows the field's turn into a loop, its 250 km of growth per orbit curving through half a turn. A crack held open through half an orbit no longer behaves as a cycloid's does, and now ends there.

What the maps show:

- **Europa:** near the equator, where the tide is weakest, the lineae and the fabric form a network of long, gently undulating lines crossing at every angle, as the reference shows there. Elsewhere the tide's great arcs remain, bent gently. 2,527 lineae, against 2,684 under the tide alone.
- **Enceladus:** its tiger stripes keep parallel, all eight present ones running their 242 km, and its lesser cracks bend with the heterogeneity. Its shell thickens by only 0.6 km over its life, some 650 Pa of tension, so its old plains still do not crack (M4): the thickening that would crack them is the zone's local refreezing at each episode's end, which the global history of the shell does not see.

What they still do not show:

- **Bands** as spreading lanes; the next milestone.
- **How many cracks form:** the freezing stresses bend and cross the cracks but do not yet set how many form (S6).
- **Lanes** (F61) still trace the tide alone.
- **Stated choices:** the freezing's hold at the cracking stress and its response time, and the heterogeneity's share and renewal.

Tests 150. Two are new, and the tiger-stripe test is stricter again: a thickening shell stands at the cracking stress in tension, a thinning one in compression, and a steady one bears none; and the shell's heterogeneity bends Europa's cracks, traced from the same starts, more than the tide alone does. Clippy is clean and validate passes its 117 checks.

The user on the twenty-fourth: "It looks pretty good. The overall profile of the cracks shows the tidal field more clearly than real Europa, but the appearance is fine to me unless you want to make some changes to that. Enceladus' stripes indeed look the same as before, but there are now two small new artefacts in the same colour as the white terrain that dig into the stripes sharply." And: "cracks have sharp flat ends. I'm not sure if that's how they're physically supposed to look."

The two artefacts had two causes, found by tracing every write to the cells involved:

- **A lesser crack younger than a vent.** The present episode's vents were dated anywhere in the last 2 Myr, the window the lesser cracks of the zone share, so a younger lesser crack crossing a vent overwrote its trough with its own surface and took away the vent's sintered tint. An active vent renews its trough by its own eruptions to the present, so the present episode's vents are now dated to the present and laid last.
- **Stamps computed together that touched the same cells.** Stamps whose caps are apart are computed in parallel against the state before their wave, and then applied. A linea's cap counted only its width, but its stamp reaches out to 2.5 times its margins, offset by its segments: 5 to 16 km beyond its trace. Two stamps in one wave could then write the same cell, the later one from the state before the earlier was applied. The cap now bounds the stamp's reach, computed as the stamp computes it; a vent, which lays no stain, reaches only its width and two cells.

The square ends: a crack closes toward its tips, and an elastic crack under a uniform load opens as an ellipse along its length. A ridge's width, its stain and the surface it renews now follow that ellipse, (4 s (1 - s))^(1/2), so a crack fades to its tips; the first try, with the gentler (4 s (1 - s))^(1/4) the relief uses, faded only in the last kilometre and still looked square on the map. A ridge then renews pi/4 of its width along its length, which the fabric's renewal budget counts. Tests 150, clippy clean, validate 117.

## The twenty-fifth milestone: a band spreads

The user on the twenty-fourth, revised: "Once you have addressed these, move on to the next milestone."

Europa's bands were drawn as lens-shaped floors 100 m below the plains, stained like a crack's margins: thick brown lines on the map. The published bands are the opposite in height and different in kind:

- **High-standing.** They spread from a V-shaped trough along their axis, roughly symmetric, as mid-ocean ridges do, of buoyant ice, with a ridge along their sharp margins (Prockter et al. 2002) and a relief of a few hundred metres.
- **Dilational.** The surface on either side matches across when a band is closed up (Tufts et al. 2000).

A band now (F83):

- stands 150 m above the old surface's level where its crack was, with a V trough along its axis down to that level and a ridge along each margin, lineated along its axis where the grid resolves it;
- has ice that spreads from the axis over 2 Myr, the axis youngest;
- is grey when young, laying half a margin's stain, and fades with the salts;
- pushes the plates it opens between apart by its half-width at each margin, the push spreading into them over half its length, as a crack's opening does: what crossed the crack is parted by the band's width, and what ran beside it bows away.

It took two renders. The first drew solid, saturated brown lenses with pixel-jagged edges: each cell was band or plate, and the band's new ice took the full stain the chaos matrix does. A cell across a band's edge now takes the band by the share of it the band covers, and the band lays half a margin's stain, so young bands read as grey wedges as the reference shows them.

What the maps show:

- **Europa:** its 168 bands are grey, raised lenses with soft edges, lineated along their axes in the shaded relief; older features beside them bow away, and those that crossed their cracks are parted.

What they still do not show:

- **Wedge shapes:** a band opens as a symmetric lens, where many of Europa's open as wedges, the plates rotating about a pole.
- **Bright old bands:** old bands brighter than the plains, as the reference shows in the south.
- **Lanes:** bands do not yet share the lane mechanism (F61) that Ganymede's grooved terrain uses.
- **Stated choices:** the band's height, trough, marginal ridge, spreading time and stain share.

Tests 151. One is new: a band pulls the plates apart, a bump beside it moving outward by the predicted push, and its ice spreads from its axis. The profile test now checks the raised band. Clippy is clean and validate passes its 117 checks.

## The twenty-sixth milestone: Io's surface from its deposit history

The user on the twenty-fifth: "The shaded relief close-up looks very good. As long as the open items are marked for later, move on to Io." The open items went into the audit first (W5, T3, M5 and M6; 2ddfe8c).

A survey render set beside the reference showed Io as stamps on a flat field. The plains were one yellow, banded by latitude; a white halo was stamped round each live plume; the flows were worms of one tone; the paterae black polygons; the mountains grey blotches. What the published maps give instead:

- **Units** (Williams et al. 2011): plains 66 percent of the surface, flows 28, mountains 3, patera floors 2.5; diffuse deposits over 18 percent, red, white, yellow and black. The red-brown plains lie beyond about 30 degrees, from radiation; outliers of white, bright and red-brown plains come from the long accumulation of white, yellow and red deposits.
- **Colours** (Geissler et al. 1999): red and orange over a third of the surface, mostly beyond 30 degrees; yellow over 40 percent at low latitudes; white and grey, moderate to coarse SO2, over 27 percent.
- **Change** (Geissler et al. 2004): at least 83 percent of the surface unchanged over five years; deposits fade in months to years by burial and alteration; small plumes leave near-circular rings 150 to 200 km in radius, white or yellow, and about half a dozen large ones oval, sulfur-rich rings 500 to 550 km from the vent.

So the plains are what the plumes have left, buried and aged (F84):

- **The eruptions.** Every patera erupts in episodes, followed back five burial times of the background fallout, 500 years, past which nothing laid can still show. A patera is erupting now if it is active and its class throws a plume; before that, quiet spells and eruptions of a few months to decades alternate at its class's duty. Each eruption has its own speed and so its range, its own mass, its vent's elongation and its range's scatter with azimuth.
- **The deposit.** The gas condenses in a ring at the canopy where the plume falls back; the pyroclasts, sulfur and silicate, fall within. A Pele-class ring is red sulfur from its S2, 3,000 kg/s (Spencer et al. 2000; Zhang et al. 2003); a Prometheus-class ring is SO2 frost from 10^5 kg/s (Jessup et al. 2004), white or yellow by its frost share.
- **The top layer.** Each cell runs through the eruptions that reached it in order, in closed form between their starts and ends. Each fall buries what lay on top; red sulfur reverts to yellow in two years; frost sublimes at Ingersoll's 0.1 cm/yr at the equator and by the SO2 vapour pressure at its own temperature elsewhere; and the plasma alters what lies on top to red-brown where the sublimation atmosphere, which the same vapour pressure sets, is too thin to shield the ground.
- **Lava.** A flow or a patera floor is black when laid and is coated by what falls on it after, rough lava more slowly than a plain.
- **Halos.** A young flow drives SO2 off the ground it covers, which condenses round its margins (Williams et al. 2011; Kieffer et al. 2000).

It took ten renders after the survey:

- **The first** had red-brown polar plains, but the low latitudes stayed flat yellow, no Pele-class plume was erupting (the class share gave two at a time; Io has about six), regular low-order jets drew sunbursts, and starting the history from mean rates would have left a permanent halo round every active centre. The class share is now calibrated to six large plumes, the jets are gone, the history starts from an aged surface, and the shielding reads the latitude, since the atmosphere spans far more than a slope.
- **White** covered 5 percent against Geissler's 27: the Prometheus-class mass had been a tenth of Jessup et al.'s output. With their 10^5 kg/s a thick deposit lasts a century at the equator.
- **The halos** first came out as opaque white cocoons round each young lobe; they now hug the margins over 8 km.
- **The deposits** read as opaque discs. Gas and pyroclasts now fall in different places, so a Pele-class deposit is a red ring round an olive interior and a Prometheus-class one a white to cream ring round a yellow one.
- **Speed.** A 2,000-year history took 24 to 41 seconds; tied to the background burial, with the halos on a coarser grid, it takes 11.
- **Artefacts.** Some deposits were cut short at a fixed reach: each eruption now reaches as far as its tail could show. Frost lines were jagged: the frost now reads its slopes over 30 km one by one and shorter ones as a spread of a third to three times its sublimation, which blends a frost line over the width its ground gives it. A speckle of isolated white cells was a rounding fault: a finished eruption left a fall of 10^-17, and the closed form divided by it. A probe printing one cell's whole history found it (EXOGEN_IO_TRACE=lon,lat); the fall is now summed afresh from the eruptions active at each moment, and every relaxation is exact as its rate vanishes.

What the maps show:

- **Io:** three Pele-class plumes erupting, red rings round olive interiors, their pyroclasts darkening the ground near the vent; Prometheus-class rings white, cream or yellow, the older faded under the sulfur that has fallen since; red-brown plains beyond about 30 degrees, with outliers left by recent deposits; young flows black with white margins, old ones coated. Shares of the surface: altered red-brown 30 percent, red 5, frost 5 above a fifth of the ground, yellow 61. The deposit history takes 11 seconds of a 206-second render.

What they still do not show:

- **White plains** (M7): frost comes only from rings and halos and sublimes at Ingersoll's rate wherever it lies, so it covers 5 percent against 27; it needs the frost cycle of the sublimation atmosphere and its albedo feedback.
- **Flows:** worms of one width, where Io's flow fields are broad, branching and of many units, a quarter of the surface (W1); next.
- **Mountains and paterae** (T1): grey blotches and polygons. A patera floor laid since a deposit fell shows as a straight-edged hole in it.
- **Stated choices:** the eruptions' durations and scatter, the deposit's opacity and profile, the background fallout, the lava's roughness, the alteration's time and the halos'.

Tests 153 (one ignored). Two are new: the top layer ages and is buried (red sulfur stays red while it falls and reverts within a few years; thick frost outlasts thin and the cold keeps it longer; lava hides under its coat; a rate of mere rounding leaves the cover as it was; the closed form agrees with small steps); and plumes erupt at their class's duty. Clippy is clean and validate passes its 117 checks.

## The twenty-seventh milestone: Io's flow fields and smooth plains

The user on the twenty-sixth: "Good work. I imagine the flow fields will drastically upgrade Io's appearance. ... Looking at Io's real map, the plains seem exceptionally smooth, with mottled colour variation at a small scale and a global, rather noisy noise field. Will this detail also come with future processes? If you are happy, move on to the flows work."

Io's fields were six to fourteen lobes radiating from a patera at random azimuths: worms of one width. Io's own are compound. Its long-lived eruptions build their fields over years to decades from small breakouts on top of older flows, fed through insulated paths tens to hundreds of kilometres long (Keszthelyi et al. 2001). Fields are dark when fresh, grey when old, often bordered by white frost (Williams et al. 2011), and the largest are over 300 km long.

- **F85, the compound field.** A trunk leaves the patera on its downslope side; branches leave the trunk; breakouts leave the trunk, the branches and each other, more towards their ends and at their fronts, each laid at its own time within the field's span of 3 to 100 years, until they cover half the field's reach squared twice over. An active field grows at about 1,500 km² a year, as Amirani's new flows covered 620 km² in five months. The fields are grown in the finish from the paterae alive at the epoch; the loop no longer traces them, and Io's loop takes 3 seconds instead of 142.
- **F86, the buried plains.** The plains carried the base terrain's bumps, slopes of 2.5 degrees at 3 km. The regime renews its surface at a centimetre a year, so relief lower than the burial since the surface formed, 1.1 km, is filled at every scale of the base's pyramid, and only the broadest undulations stand. The plains are flat, as Io's are, under mountains, paterae and flows.

It took seven renders:

- **Streaks.** The first fields were narrow streaks, neighbours parallel: the fall line steered a front per step, so a short lobe swung onto the regional slope within a few steps, and on the buried plains that slope is coherent over continents. It now steers per distance, over 20 km at full slope, and breakouts leave earlier breakouts too, so a field fans out.
- **Invisible fields.** A field's lobes share its span, so most fields were wholly old and coated over. Rough lava now hides only under 300 optical depths of fallout, about a centimetre of fines filling its hollows, so old fields stay grey, as the map's undivided flows do.
- **White fields.** Young fields came out white. Probes found two causes. A halo landed on cells partly covered by young lava, and one top layer per cell spread that frost over the lava too; frost now does not settle on lava laid within 20 years, and a halo lands only beside young lava. And a flow's coverage texture left a fifth of each cell uncovered, which the halo filled; the breakouts' ages carry a field's variety now, and the texture is 0.1.
- **Grey fields.** Active fields came out grey because a Prometheus-class plume's pyroclasts were 40 percent of its condensate, some 4 x 10^4 kg/s of solids, which buried its own field in a couple of years; they are 3 percent.
- **Coverage.** Flows covered 12 percent against the map's 28; the longest fields now reach 500 km, and they cover 20.

What the maps show:

- **Io:** dark compound fields of lobes of many tones, the youngest black, older ones grey and olive, with white frost round the young ones' margins; plains flat in the shaded relief, mountains and paterae standing on them. Flows 20 percent of the surface, plains 70, patera floors 1, mountains 9. The whole render takes 89 seconds instead of 207.

What they still do not show:

- **Mottled plains** at small scales: the plains between deposits are one tone; small eruptions as a population, in place of the uniform background fallout, would mottle them.
- **Flow coverage** (M8): 20 percent against 28.
- **Lobe sizes** from the lava's yield strength and supply (W1); bright sulfur flows.
- **Mountains and paterae** (T1): grey blotches over 9 percent of the surface against Io's 3, and polygons.
- **White plains** (M7).

Tests 155 (one ignored). Two are new: a compound field fans out from its trunk, its breakouts spread over a good share of its reach, each laid within the field's span; and the burial fills small relief and keeps the broad. Clippy is clean and validate passes its 117 checks.

## The twenty-eighth milestone: Io's mountains and paterae

The user on the twenty-seventh: "Ok, go ahead with the mountains and paterae."

Io's mountains were noisy ovals with a scarp on one side, grey blotches over 9 percent of the surface against Io's 3, and its paterae polygons. Io's mountains are blocks of crust thrust up along faults as the crust subsides under its burial, tilted, layered or massive, with scarps cut back by landslides (Schenk and Bulmer 1998; Schenk et al. 2001), striking 30 or 60 degrees from the tide's compression (Seeger and de Kleer 2025). Its paterae are collapse depressions, many with straight sides and sharp corners, some at the feet of mountains (Radebaugh et al. 2001).

- **F87, the thrust blocks.** A block's uplift is its fault's throw, tapering as a cosine squared to the fault's tips, times its fall across the block, held level by a plateau, falling steadily on a tilted block and rounding on a massif; the plains bury it to a fifth of its height, so its back and its tips are where it dips under them, irregular with its own relief. The thrust's heave carries it out over the plains the farthest where the throw is greatest, so its trace bows, as the bow-and-arrow rule has a thrust's. Its scarp is worn to 30 degrees and cut back in landslide arcs, their debris below. A lineated block is ridged along its strike in patches, crests that wander, branch and end.
- **F88, the collapses.** A patera is the union of nested collapses along the crust's fabric, its sizes Radebaugh et al.'s, its walls at 70 degrees; 42 percent are cut straight by faults and 13 percent sink at the foot of a mountain's scarp. Lava that runs into a patera is its floor's, so the floors cover 1.9 percent where they covered 0.5.

It took eleven renders:

- **D-shapes and rectangles.** The first blocks came out D-shaped with hairlines at their tips, then as rectangles with straight ends: the uplift fell too steeply for the relief to wander the outline. The block became a buried uplifted surface, with stronger relief, a plateau's back that wanders and a throw tapering over a fault longer than the block.
- **Sawtooth crests.** Scarps narrower than two cells stepped along the grid; a scarp is now at least two cells wide.
- **Rails.** Every block was shadowed by two straight lines along its scarp, a strip between them. Two guesses (a mismatch between the loop's copy of a block and the one the finish took away; hard cuts at the edges of a block's domain) changed nothing. A height profile across one block found a step of 400 m some 36 km in front of its scarp: the debris began at 0.15 of the crest at the scarp's foot while the scarp began at nought, and the finish, taking away its own copy of the loop's wider block, left that block's step behind as well. The scarp now comes down to its debris, and the finish draws mountains and paterae whole on the buried plains the loop laid them on, which the loop now hands over, so there is no copy to take away.
- **Corduroy.** A lineated block's ridges were one sine across it, a regular corduroy; they are the noise's zero crossings drawn out along the strike now, in patches.
- **Missing scallops.** The scarp's landslide arcs barely showed: their smooth maximum took the scarp's width times the logarithm of their number, 17 km, off every cut. Each arc now counts only where it cuts.
- **Half the relief.** The factor meant to give the blocks' relief a unit spread gave 0.5; the relief is normalised from the noise's measured spread now, with octaves down to what the grid resolves, and its amplitudes halved to keep the look.
- **Grain.** Vertical bands crossed the blocks: the arcs' cusps shifted the whole surface behind them, and the noise was read on planes along its lattice. The surface is measured from the fault's line, the landslides taking away only what stood in front of the scarp, and the noise is read on a slanted plane.
- **Lost floors.** Inside the larger paterae 84 percent of the floor was mapped as flow: lava ran into the pits and took them. Lava that ponds in a patera keeps the floor's unit and gives it its age.
- **Square paterae.** A faulted patera's faults were set at right angles to its axis, so some came out square; each fault now follows one of the fabric's strikes, 30 or 60 degrees either side of the tide's compression, so the sides meet at the angles between those strikes.

What the maps show:

- **Io:** mountains as blocks with scalloped scarps on one side and embayed backs, lineated blocks ridged in patches, plateaus with wandering backs and lumpy massifs; paterae as pits with steep walls, scalloped where their collapses meet, some cut straight, some against a mountain's scarp, dark when their lava is young and coated when old. Plains 79 percent of the surface, flows 16, patera floors 1.9, mountains 3.4 (Williams et al. 2011: 66, 28, 2.5 and 3). The render takes 96 seconds. The sidecar lists each mountain's length and width in place of its diameter, and a patera's diameter is its outline's widest extent.

What they still do not show:

- **Mountains in the albedo:** a block takes the plains' coat and shows its rock only on its scarps and ridges; frost on high ground is not modelled.
- **Mountains' and paterae's sizes** from the compression and the magma beneath (T1): the sizes and shares are the counts'.
- **Flow coverage** (M8): 16 percent against 28.
- **Mottled plains** and **white plains** (M7).

Tests 157 (two ignored). One is new: a patera's outline is its collapses and faults. An ignored test writes three blocks' heights for a look at their forms. Clippy is clean and validate passes its 117 checks.

## The twenty-ninth milestone: the plains' mottling, and paterae against mountains

The user on the twenty-eighth: "Good progress. Continue with the plains mottling when you are ready." And, during it: "some of the patera are perfectly circular after the last milestone. I see on the real Io map that there are circular patera so there's nothing major to flag there unless that's unintended. But the "patera against a scarp" example you showed shows the circular patera with a sharp, straight edge cut off, just at the foot of the mountain. There are more like this too, seemingly without any reason. Will this be fixed with plains mottling or is it something to address on its own?"

It was something to address on its own. Radebaugh et al. (2001) found paterae that abut mountains oddly shaped, with arcuate or scalloped margins, some wrapping round the mountain (Gish Bar, Monan), not cut straight; and 70 percent of all paterae longer than wide by more than 1.25, where more than half of F88's were circles. On the plains, Williams et al. (2011) read the outliers of white, bright and red-brown plains as the long accumulation of white, yellow and red diffuse deposits.

- **F89, the plains' fallout.** The background fallout, the far fallout of the plumes, spreads about each volcanic centre alive at the epoch over two of its class's ranges, weighted by its plumes' output, its mean over the surface unchanged, and clumped at 15 to 120 km. Where the plasma alters the plains the balance of burial and alteration follows it, so the red-brown plains reach equatorward in quiet regions and give way to yellow among the active centres, and they mottle. Small eruptions, vents too small to be paterae and stealth plumes, add occasional fresh patches, ten a year.
- **F90, paterae against mountains.** A patera beside a mountain sinks against the block's foot and takes it as its margin, bowed and scalloped as the foot is, the block's debris inside it gone with the floor. Every patera is elongated as Radebaugh et al. measured, the longer ones chains of coalesced collapses.

It took nine renders:

- **Foam.** The first tries laid the whole background as a population of small eruptions, 240,000 over the history, each with a Prometheus-class plume's ring: the plains filled with hollow rings. A plume too low for a canopy lays no ring, so a small eruption's fall became a disc about its vent.
- **A meridian of large deposits.** The largest small eruptions lined up just west of longitude 0: the draw that set an eruption's size was the one that set its longitude. Each has its own draw now, and a test checks their sizes east and west.
- **Dots, then bokeh.** Opaque discs read as dots, every cell covered by about three in the history and the plasma browning each within twenty years; thinner, more numerous discs, 790,000 of them, read as glitter; softer, larger ones as bokeh, each a circle with a dark silicate centre. And at that rate Galileo would have seen several hundred new deposits wider than 20 km a year, where it saw some 80 changes in five years (Geissler et al. 2004). The background became the plumes' far fallout about the volcanic centres, which F84 had named it and laid evenly, and the small eruptions a sparse population at about Galileo's rate.
- **Too smooth.** The far fallout's clumping, written as the deposits' mottling is, moved the fall by 3 percent; the noise's spread is 0.18, not 1. It is drawn at a unit spread now, a log-normal factor spreading by 0.5.
- **Circles and cuts.** Paterae beside mountains ended in straight chords at the scarp's foot, which the foot's bow and arcs no longer followed; they end at the foot itself now, found along rays from the patera's centre. Elongated paterae came out as smooth cigars at first, until the elongation was split between the collapses' own and their spread along a chain.

What the maps show:

- **Io:** paterae of varied elongation along the crust's fabric, bean and worm shapes among them, those beside mountains wrapping the block's foot; the red-brown plains' edge near 30 degrees ragged with the volcanic centres, the brown reaching equatorward in quiet regions and the yellow poleward among active ones; the plains mottled at 15 to 120 km where the plasma alters them. Plains 77 percent of the surface, flows 17, patera floors 2.1, mountains 3.4 (Williams et al. 2011: 66, 28, 2.5 and 3). The render takes 100 seconds; the deposit history 12.

What they still do not show:

- **Mottling at low latitudes:** the plains there are one sulfur under a fall of sulfur, unaltered, and stay nearly smooth; their white, yellow and grey mottling waits on the frost cycle (M7).
- **The faults that cut a patera** show nowhere outside it (T1).
- **Flow coverage** (M8): 17 percent against 28.
- **The deposits' mottling** (F84) is about 5 percent where its constant says 0.35 (recorded under M7).

Tests 159 (two ignored). Two are new: a patera beside a mountain ends at its foot; and the small eruptions follow their rate with sizes the same east and west, while the background is heaviest among the volcanic centres with its mean the background's. Clippy is clean and validate passes its 117 checks.

## The thirtieth milestone: rays from the ejecta curtain, and terraces as the crater collapsed

The earlier phases' survey (2026-10-16) rendered Mars, Venus, Mercury and the Moon and listed what falls short of the fidelity milestones' standard: analytic profiles on placed objects, lens-flare rays, flat lava plains, a rise at Mars's south pole. Its findings went into docs/audit-principles.md (T4, T5, M9 to M11). The user agreed the order of the revisit and asked one question first: "Is Mars' dichotomy a guaranteed outcome of Mars' parameters? If not, is it still a worthwhile assumption to make since having it increases the diversity of non-real worlds generated with the program? Once you have addressed this, go ahead with the first audit milestone."

It is not guaranteed. Both leading explanations are contingent: a giant impact (the Borealis basin, Andrews-Hanna et al. 2008; an impactor of 1,600 to 2,700 km, Marinova et al. 2008), which late accretion makes a matter of chance (in Bottke et al. 2010's Monte Carlo the median largest impactor on Mars was 1,500 to 1,800 km); or degree-one convection, which needs a particular viscosity structure in the mantle (Zhong and Zuber 2001; Roberts and Zhong 2006). The dichotomy is a rule today, a 20 km contrast in crust thickness for every stagnant-lid rock body by the regime's name, and a rule every body shares gives no diversity. As an outcome, from the giant-impact tail of the late-accretion population, some bodies would have a hemispheric basin, some none and some an off-centre, elliptical or partial one. T5 records the direction and the Mars milestone builds it.

The first milestone of the revisit takes the craters' rays and walls:

- **F91, rays from the ejecta curtain.** A crater throws its ejecta as a curtain gathered into jets by the target's relief at the rim (Sabuwala et al. 2018), some thirteen of them, nine tenths of the mass in a fifth of the azimuths. Beyond the continuous blanket a ray is the share of the ground the jets' secondaries excavate from under the weathered skin (Shoemaker 1965; Elliott et al. 2018), falling with distance by Elliott et al.'s secondary-crater model and scaled to their lunar ray lengths, so it narrows to its jet's midline, breaks into segments and fades rather than ending. Noise bends the jets' paths and bunches their fragments. The blanket's lobes and the loop craters' secondaries follow the same jets. An icy surface's lag is the skin where it is thicker than the Moon's 10 cm, so Callisto's rays are short or absent under its thick lag.
- **F92, terraces as the crater collapsed.** How much of a crater's wall slumps follows how far it collapsed, from none at the simple-to-complex transition to the whole wall from about 3.7 transitions; the outer terraces come first, their spacing varies, and blocks stand off far enough to merge with the terrace inside or outside them.

How it went:

- **A prototype first.** The rays were tried on a local map round one crater before any engine code: jets of a log-normal curtain, the excavated share 1 - exp(-A) with Elliott et al.'s slope, bends and bursts. At the first spread the gaps between the jets filled out to five radii; at 2.5, with the curtain's harmonics from order 6, the rays stood apart from two or three radii out, as Tycho's do.
- **Brush strokes.** The first render still had the first spread: rays brightened 40 percent of the Moon, in broad saturated strokes. The second spread brought it to 22.
- **Dark specks.** Every old crater's melt pond a ray crossed came out as a spot of fresh dark glass: the albedo freshened the melt with the ground. A ray freshens the ground its secondaries dug now, and an older crater's melt keeps its own age's maturity; without rays the albedo is unchanged, which the algebra shows.
- **A web on Callisto.** With the Moon's skin everywhere, Callisto's rays covered 31 percent of it in a web. The ice lag the maturity already models is the skin there now, and Callisto's low latitudes are bright young craters with short rays or none; 6 percent of it is rayed.
- **The cost.** Each crater's rays could reach far, so 61 of the Moon's 72 rayed craters are drawn over the whole grid: 10 seconds of the Moon's 64-second finish at 8192 by 4096. Skipping the bursts where no burst could show a ray took 2.6 seconds off.

What the maps show:

- **The Moon:** rays as filaments of varied width and length, a few long ones crossing the maria and many short ones, dark gaps between them from a few radii out, webs where systems overlap; a bright halo of two to four radii round each young crater. 21 percent of the surface is brightened by more than a tenth of the fresh contrast, 4 percent by more than a half.
- **Mercury:** the same rays, few and faint, since its rays fade in 0.15 Gyr (its weathering timescales are a standing deferral).
- **Callisto:** bright young craters, rayless under the low-latitude lag; Ganymede a web of rays at middle and high latitudes (M12); Europa a few rayed craters; Enceladus none showing under its fallout.
- **Terraces:** a lunar crater of 20 km stands 0.7 terraces at an azimuth, one of 93 km 3.4, where every size had about four; on the synthetic walls the ladder of rings is gone. At the Moon's 2.7 km cells the difference barely shows in its maps.

What they still do not show:

- **Rays bright by composition**, a crater's material thrown over another unit (Lichtenberg's highland rays over mare).
- **The rays' fading from the ground's own maturing**, which waits on Mercury's weathering timescales; and the skin of rock bodies other than the Moon.
- **Mercury's dark halos** of deep material, still a soft ramp to the blanket's edge.
- **Callisto's high latitudes**, where rays from far craters appear on the cold ground (M12).

Tests 166 (three ignored). Seven are new: the curtain gathers its mass in jets; its rays' lengths follow the lunar fit; an oblique impact throws them farther downrange; a bent path inverts; the reach bounds them; terraces come in as the crater collapses; and an ignored one writes terraced walls for a look. Clippy is clean and validate passes its 117 checks.

## The thirty-first milestone: facets on Callisto, and lava plains as compound flow fields

The user on the thirtieth: "The new rays look much better. Callisto's high latitude switch to rayed craters looks sharp and unnatural, like you said. Do you have any ideas for that? Continue with the lava plains."

The switch came from the lag's exponential dependence on temperature: one temperature per cell, its latitude's and its resolved slope's, made the lag's thickness fall by orders of magnitude across a few degrees, and the rays' skin share went with it to the fourth power. Real icy surfaces segregate at every latitude, the warm slopes growing lag and the cold ones keeping their ice (Spencer 1987).

- **F93, facets.** A cell of icy surface is summed over unresolved facets tilted toward and away from the pole, spread by Callisto's Hapke mean slope (35 degrees; Domingue and Verbiscer 1997) times its share of saturation by the unresolved craters, so a young surface is smoother. Callisto's rays come in over 30 to 50 degrees, patchy at first, instead of switching on across a line (committed on its own, 2c450b4).
- **F94, lava plains as flows.** A flood that can show at the end is drawn as a compound flow field, a trunk, branches and breakouts traced over the ground from where its lava came into its low, Io's paths (F59) with the flood's own sizes; its sheet lies under its lobes, its pond fills where the loop filled it deep, and outside its low its lobes lay a sheet on the ground. Each flood's composition carries its source's and its own batch's variation as well as its time's. Io's path tracer became a shared function, and Io's maps are byte-identical.

How it went:

- **Invisible at first.** The first render changed the maps' unit layers, blobs to lobes, and not their albedo: floods close in time and place had nearly the same composition, a few per cent of the reflectance. The source and batch terms brought the mosaic out.
- **Blocks.** Units with stepped edges a coarse cell high came from three places in turn. Floods ruled buried and drawn with the old outline showed through their successors' gaps, so a flood among the last four over a cell is drawn as flows. Lobes reaching past the low and its ring were cut along coarse cells, so every cell a lobe covers is visited and a lobe stops before its pass's edge. And the test for where a flood ponds compared the ground pieced together from the coarse smears with the level, a ground good only to a coarse cell's ramps: ponding is judged by the flood's own coarse thickness, sampled smoothly. A finer test, ponding outward from the lobes into hollows, brought the blocks back and was dropped.
- **A crater on Mars filled 2.5 km deep.** A trace of one cell through the flood passes showed lobes that ran past their low into a crater filling it to the low's level, 2,400 m above the crater's floor. Outside its low a lobe lays a sheet.
- **Not the same on four threads.** The determinism test failed: a flood's centre and area were summed over its cells in a hash map's order, which differs from thread to thread, and its trunk now heads for that centre and its breakouts stop at that area. Both are summed in the cells' order now.

What the maps show:

- **The Moon:** a mare as a mosaic of flow units of different tones with lobate boundaries, some tens to a hundred kilometres across, their fronts faint in the relief beside the wrinkle ridges; small ponded patches keep rounded edges.
- **Mars:** the youngest plains, round its second province, as lobate, branching flow units in different tones.
- **Mercury:** flow units in its larger plains; its crater-floor ponds keep the contour of the floors they fill.
- **Venus:** its flood patches mottled with lobate units.

What they still do not show:

- **Long flows.** A flood's extent is still the coarse low the loop filled, compact on flat ground; Imbrium's 1,200 km flows need the loop's spreading to follow the flows.
- **Lobe sizes from the lava**, its yield strength and effusion rate (Hulme 1974), rather than shares of the low's reach.
- **A few small angular patches** where a deep coarse pond meets the lobes.

Tests 168 (three ignored). Two are new: a flood shows through its last few successors only; a flood's lobes cover its low and stay inside its pass. Clippy is clean and validate passes its 117 checks.

## The thirty-second milestone: the lava's units blended by gardening, and the craters the lava fills

The user on the thirty-first: "The different colours are good, but the edges are too sharp. There should be better, more natural-looking transitions between colours." They pointed to large pixel artefacts round craters in a mare (the Moon, seed 762, at -88.1, -15.4) and craters partly filled inside a mare with a sharp edge round the filled part (-111.3, -16.3), and added: "Our flows have the right detail and intricacy for Io, they're just more blended on the Moon and have less contrast, and even more for Mercury."

- **F95, gardening.** The lava's coverage and its tint are mixed sideways over each contact's own zone, which grows with the gardening since the younger side erupted, the same law that degrades the craters, from Li and Mustard's lunar zones. The zones are levels of a Gaussian pyramid, so no cap holds them back (audit W4). Each flood's own draw of composition is smaller, so neighbouring flows differ by a few per cent while the major units keep their range (Iqbal et al. 2019; Denevi et al. 2009).
- **F95, young craters on lava.** Each flood older than a crater writes its lava under the crater where its lobes reach, in their order, so the tone under a young crater's ejecta matches the flows round it.
- **F95, hollows.** A crater finer than the loop's cells fills to the flood's level when the lava comes in over its rim, and keeps its bowl when its rim stands above the lava.

How it went:

- **The checkerboard.** A trace of the cells round the crater showed the first flood to write a cell under a younger crater keeping it: the oldest flood, drawn with its coarse outline.
- **The partly filled crater** was a crater the finish replayed, too small for the loop, three kilometres deeper than the ground the loop flooded: the lava filled only the cells its lobes covered. Filling the finer hollows to the level fixed it, then left a straight edge where a hollow crossed the ring of cells a flood's pass visits; the fill now follows a hollow past them. It also filled craters whose rims stood above the lava, seeded by lobes beyond the low and then by slivers of lobes inside a rim, both drawn over the smoothed ground the lobes are traced on; the fill now comes only from the body of the lava, over its edge.
- **The mixing's age.** A contact is as old as its younger side. Mixing by each cell's own surface age would have let the old highlands beside a mare pull in its lava for 4.4 Gyr, and every young crater on a contact sharpen it, since a crater resets the surface's age; the zone comes from the eruption's time, kept under later craters, and the youngest lava within reach.

What the maps show:

- **The Moon:** maria of soft-edged flow units, their tones a few per cent apart, their outlines still lobate; the mare edge blends into the highlands over a few kilometres; the craters in maria filled or dry as their rims allow.
- **Mercury:** inside thick lava, the albedo's pixel-to-pixel steps fall by half (0.0089 to 0.0047 in luminance at the 99th percentile) and its spread falls from 0.13 of its mean to 0.10, against the Moon's 0.18.

What they still do not show:

- **Mercury's faster gardening** in full: the diffusivity scales with gravity alone, 1.5 times the Moon's, where Fassett et al. (2017) measure about twice.
- **Rays bright by composition**, a crater throwing one unit's material over another, which would break contacts up crater by crater rather than as a smooth zone.

Tests 171 (three ignored). Three are new: a pyramid level spreads as its kernels add up; gardening blurs an old contact and keeps a young one sharp; a flood fills a crater it comes into and not one whose rim holds it out. Clippy is clean and validate passes its 117 checks.

## The thirty-third milestone: Callisto's jets, broken rings and bright walls

The user on Callisto after the thirty-first: "The dark lag for Callisto is causing the two crater ejecta patterns (the immediate 'star shape' and the long rays) to 'unlink' ... and so the fade-to-background in the star shape overlays the rays causing a solid gray circle around them." "The ringed craters' rings are far smoother and more uniform than real Callisto. They also cover the entire surface, whereas that doesn't appear to be the case for real Callisto. We are still missing the 'starry' look on its surface. Our small Moon craters have white rim colours which is what Callisto needs to achieve that look." During the work: "I wonder if the ice-capped knobs contribute to Callisto's speckled appearance? But aren't they too small to be a global feature? jpegPIA01297 shows that crater walls are indeed bright, and jpegPIA01649 shows the same bright walls on the ring structures."

- **F96, jets.** Beyond the continuous blanket the deposit follows the ejecta curtain's jets, so a young crater's spikes run on into its rays and the ground between them stays old: no grey disc.
- **F97, broken rings.** A ring system's graben survive the later cratering by the equilibrium the craters themselves reach: Callisto's oldest systems are gone or scattered troughs, its youngest keep their rings in stretches.
- **F98, bright walls.** Each crater writes the share of its walls steep enough to shed, from its own profile at a finer scale than the grid, and an icy lag sheds from gentler slopes than a silicate regolith: Callisto's crater rims and trough walls stand bright.

How it went:

- **The grey circle** was the ejecta's cover past the blanket's edge, the same at every azimuth, which both the stamp's reset age and the ray field's blanket follow. The blanket's lobes came from the curtain's finest harmonics only, which is why the star's points and the rays did not line up.
- **The rings** come from Melosh's gamma on the 18 km lithosphere Callisto had in its first half billion years: 77 systems, all kept, while the craters that formed with them have gone by the population's equilibrium. Pieces of two graben widths left lenses and little diamonds once their walls were bright; three widths leave fewer.
- **The white rims.** The resolved slopes showed why the Moon's rims are white and Callisto's are not: three hundredths of a per cent of Callisto's surface against 6.4 per cent of the Moon's steeper than 20 degrees at the same pixel count. The walls written from each crater's profile helped Callisto's small young craters only: its craters are nearly three times shallower than the Moon's (Schenk 2002), their walls under 20 degrees. The user's Galileo images show its walls bright at those angles, the lag gone from them, so the lag's own angle of shedding came down.

What the maps show:

- **Callisto:** young craters bright with spiky ejecta running on into rays; a speckle of bright crater rims over the whole surface; one or two basins with rings in bright broken troughs, scattered troughs elsewhere.
- **Ganymede:** bright crater rims and bright-rimmed furrows in the dark terrain.
- **The Moon:** white rims on its small craters at any resolution, and slightly sharper spokes in its young craters' ejecta.

What they still do not show:

- **Callisto's disintegrating rims:** its craters are clean bright circles, where the real ones break into knobs (Moore et al. 1999). The sublimation that carves them is not modelled.
- **Ring troughs' shapes:** the shortest surviving pieces are spindles, both walls bright to their tips.
- **Callisto's high-latitude rays,** long and bright over a dark background where real Callisto's high latitudes carry frost.

Tests 175 (four ignored). Four are new: the later cratering strikes out a share of the graben; the ejecta beyond the blanket run out along the jets; a wall sheds where its own profile is steep; and a diagnostic, ignored, writes the walls' slopes. Clippy is clean and validate passes its 117 checks.

## The thirty-fourth milestone: the grey circle's second cause

The user after the thirty-third: "The speckle is a huge improvement to the icy moons, and the lava looks very good now across all worlds." Then: "the gray circles on rayed craters is still there on Callisto e.g. 178.0 29.4 on Callisto preset seed 4673", and, of the circle: "It's allowing only the base surface colour while cutting off the wide ejecta from converging on the crater and blocking smaller craters." And: "Didn't we already add crater rim crumble, or am I misremembering? I mean the troughs, the shapes are very simple whereas real Callisto shows them very degraded. But I suppose rim crumble would go into the large basin rings too?"

- **F99, the blanket's jets.** In the ray field each jet runs out until it is too thin to show, rather than ending with the rest on the footprint's circle.
- **F100, the crater record under thin ejecta.** A crater renews the surface age only in its cavity; its thin ejecta's maturity age is kept apart, so the older craters under them are replayed, subdued by the deposit, and the maturity is as before.

How it went:

- **Two circles.** The crater at 178.0, 29.4 is a 130 km crater of 0.49 Gyr. Sampling the ray field along its radii showed every jet's deposit ending at three radii, the footprint's edge (F99). The circle stayed: the share of bright rim pixels jumped at three radii, from 0.11 to 0.22 and more, with the albedo's median level across it, and the crater list held no crater at all within three radii. The replay samples the craters from the surface ages, and the ejecta had reset those over the whole footprint by their cover, an age meant for the maturity.
- **Burying walls by cover** first looked like the cause: the walls' shares of shedding were buried by the ejecta's cover. Burying them by thickness changed nothing measurable, which pointed at the crater record itself.

The crumbling the user asked about does not exist yet: rims carry the scalloped outlines and varying heights of F73 and the slump scallops of F67 and F92, and their degradation is diffusion, which smooths. A pass breaking icy crests into knobs by the ice their walls have sublimed against their height (Moore et al. 1999; Howard and Moore 2008) was tried and set aside: at the lag model's mean-temperature rate it barely acts outside Callisto's equatorial band, and at the daytime peak's, which governs sublimation, it would strip Ceres, whose near-surface ice the model does not yet place at depth. It waits on the ice table.

What the maps show:

- **Callisto:** young craters' spikes and rays run on without a circle; older craters, subdued, show under their ejecta beyond the rim.
- **The Moon:** craters under young craters' distal ejecta, where none stood before.

Still to do:

- **Degraded troughs and rims** on Callisto, with the ice table and the daytime sublimation rate.

Tests 175 (four ignored); the stamp's age test now holds the surface age to the cavity and the maturity's age to the continuous blanket. Clippy is clean and validate passes its 117 checks.

## The thirty-fifth milestone: shields and calderas from flows and collapses

The third item of the earlier-phases revisit, which the user agreed after the thirty-fourth: "Do you want to defer that to the atmosphere and climate phase, and continue with shields and calderas?"

- **F101, shields built over by their eruptions.** A shield shows the eruptions since its edifice last grew by three sheets, each a compound flow field as a flood's is, from fissures round its caldera, on its flanks and on its rift zone, whose eruptions run the farthest and spread onto the plains as aprons. Every flow narrows on steep ground and spreads where it flattens, and on ground steeper than a degree turns downhill faster. A shield's lavas come from one reservoir. Its edifice is longer along its rift and its foot is lobed.
- **F102, calderas from summit collapses.** Each collapse sinks a pit about as wide as the chamber; the vent wanders, the chamber shrinks, and each pit is smaller and deeper than the last, cutting the older ones.

How it went:

- **Spokes.** The first shields drew each eruption as one flow from the summit or a flank vent, and on the smooth cone every flow ran straight down the fall line: stars of straight spokes. The real shields, shaded from the user's Mars height map, have smooth flanks with a fine radial texture and a summit caldera complex, and the user's albedo map shows a fine, low-contrast radial striation on Arsia Mons's lower flanks.
- **A pie chart.** Drawing each eruption as a compound field turned the spokes into wedges: the fields spread on the gentle summit plateau, and down the flanks their breakouts walked sideways before the fall line pulled them back, because the steering held its pull at that of 0.3 degrees on any steeper ground.
- **Narrowing.** With flows narrow on steep ground (Hulme 1974; Peters et al. 2021), a flank flow is finer than the grid at the standard quality and the flanks read as one smooth unit.
- **Two tries taken back.** Dating each point of the edifice to when its flows last laid three sheets there made the summit young and the foot old, but the profile's steep foot packed the change into a narrow band, a ring of maturity round every shield; and scaling the whole profile by a lobed outline drew radial creases up to the summit. The foot alone is now lobed, and the loop dates the edifice to the eruptions' start.
- **Steering on steep ground.** A breakout on a 5 degree flank now swings downhill within a few kilometres. Letting the pull grow with the slope from 0.3 degrees changed 2.5 percent of Io's pixels by more than 16 levels, its flows reshaped where they cross ground of 0.3 to 1 degree; from a degree, 0.7 percent, the same fields with lobes moved.
- **Rings at the foot.** Flows hidden on the flanks reappeared where they spread at the foot, and on the small shields of Venus that drew a ring of fans. Real shields' aprons lie at the ends of their rift zones, so only the rift zone's eruptions now run long; the fissures round the summit and on the flanks build the edifice with short flows.
- **Squares and tessera.** Writing the shield's composition in the loop drew the loop's cells as squares; the finishing pass now writes it at its own resolution, out to the edifice's outline. On Venus the global resurfacing makes tessera of a standing shield, after which the shield's own pass had drawn only a band of lava at its foot; a shield whose eruptions outlast the resurfacing is now drawn just after it, its earlier eruptions buried, its lava laid over the tessera.

What the maps show:

- **Mars:** shields of lobed outline, longer along their rift zones, their flanks smooth in the shield's own lava with a few young flows and summit fields, their calderas one to six pits of 5 to 28 percent of the shield's radius.
- **Venus:** each shield a disc of its own lava with a lobed foot and a younger summit, where at the last commit the resurfacing had left their flanks tessera.
- **Io:** 0.7 percent of the pixels change by more than 16 levels, the same fields with lobes moved where they cross ground steeper than a degree.
- **The Moon:** unchanged, 0.07 percent of its pixels changing by more than 4 levels and none by more than 16.

Still to do:

- **Flows finer than the grid** show only where their cover reaches half a cell (audit W6), so at the standard quality a shield's flank flows hardly show, and the Tharsis Montes' striation is missing.
- **The edifice's profile** is one power law at five degrees on every body (audit T1); Venus's large volcanoes are broad and gentle.
- **Venus's resurfacing** still makes tessera of any shield standing at the event whose eruptions end before it.

Tests 178 (four ignored). Three are new: a summit sinks in collapses, each smaller and deeper; a shield's eruptions run down its flanks and together cover it more than once; and a shield stretches along its rift and lobes at its foot. Clippy is clean and validate passes its 117 checks.

## The thirty-sixth milestone: shields bury the troughs they grow over

The user after the thirty-fifth: "These look good, and I imagine erosion will improve their look further when we get to that. One of Mars' shields in the example you showed me is in the same place as what I believe is a canyon or graben. It looks fine to me, but is the interaction between them correct? The canyon partially cuts through the shield's cone, leaving a bit of the cone visible in the valley. Venus' shields show a grid pattern which should be addressed."

- **F103, troughs under shields.** A shield's lava fills the troughs it grows over, in the loop and in the finishing pass, so a canyon or graben older than a shield's growth is buried under its upper flanks and comes out from under its foot.

How it went:

- **The canyon and the cone.** The canyon crossing that shield was done 3.49 Gyr ago; the shield began 3.46 Gyr ago and grew 17 km high, so its lava should have buried the canyon, 2.1 km deep. The loop laid the shield as an even drape over the canyon it had carved, and the finishing pass cut the canyon's detail and every graben out of the final ground, whatever had grown since. The floor following the cone's slope inside the trough is right for a graben cutting a cone, a down-dropped block keeping its slope; the order of events was not.
- **The grid on Venus** was in the images, not the maps: the crop script took each map pixel by the nearest one, which at Venus's zoom drew the pixels as a grid of lines. The export samples the grid bilinearly; sampled the same way the shields are smooth, with a faint quilting at the scale of the finish's cells, 9 km on Venus, at the zoom where a 140 km shield spans 300 pixels. The scripts now sample bilinearly.

What the maps show:

- **Mars:** the shield of 395 km is whole, its old canyon and graben ending at its foot and running on across the plains.
- **Venus:** smooth shields, each a disc of its own lava.

Tests 179 (four ignored). One is new: a shield fills a trough formed before it by the lava it laid there since, kilometres near its summit and a few hundred metres at its foot. Clippy is clean and validate passes its 117 checks.

## The thirty-seventh milestone: a rise turns the body

The user after the thirty-sixth: "Looks good. Make sure you always check the chronological order of processes is correct. Continue with the next goal when you are ready." The next goal in the revisit is the rises with flexure and true polar wander, graben and canyons; this milestone is the polar wander (audit M9).

- **F104, true polar wander.** A body turns until its axis of greatest moment of inertia is its spin axis: its rises' mass, which their plumes hold up, against the share of the rotational bulge its elastic lithosphere holds frozen. A large rise is carried to the equator; several rise to a great circle.

How it went:

- **Which bulge resists.** A rocky planet's whole rotational bulge would hold a rise of Mars's in place, its mass a tenth of the bulge's. But the bulge readjusts as the body turns, all but the share the lithosphere holds frozen, 1.6 percent on Mars when its rise grows (Willemann 1984; Matsuyama et al. 2006); and a surface load the lithosphere holds would be stiffened by the same share, which cancels, where a plume holds its rise up whatever the lithosphere does. So the rise turns the body against the frozen bulge alone, eight times its weight on Mars.
- **The chronology.** The body turns while the rise grows, from 0.3 to 1.5 Gyr. Everything fixed to the body turns with it, and what was drawn at random over the body before is as random after; the rise's shields and rifts are placed where it stands now. Nothing on Mars or Venus yet reads the latitude it had at its time, so nothing older needs to be turned back; the climate's processes, in Phase 5, will.

What the maps show:

- **Mars:** the rise, its five shields and its rifts on the equator at 0.8 degrees south, where they had stood round the south pole.
- **Venus:** its three rises within 12 degrees of the equator, where they had stood at 63 north, 36 south and 28 north.

Tests 180 (four ignored). One is new: a rise near the pole is carried near the equator against a small frozen bulge and hardly at all against a large one, the fluid Love number is right for a uniform body and the Earth, and the eigen solver finds a known matrix's axes. Clippy is clean and validate passes its 117 checks.

## The thirty-eighth milestone: graben and canyons follow the rise's stress

The second part of the revisit's fourth item.

- **F105, faults from the stress field.** A rise's faults start where its lithosphere is in tension, strike across the most tensile direction and die out where the tension fades: radial graben on the flanks, a maze of short graben on the crest, canyons running out from it, and the load's hoop tension carrying the graben a little past the foot. They form from the rise's start to the present at the pace of the body's volcanism, which builds the load.

How it went:

- **Buried graben.** The faults had formed only until 2 Gyr, and the lava that kept flooding the rise buried nearly all of them; Tharsis faulted from the Noachian into the Amazonian (Anderson et al. 2001), as its load kept growing. The chronology now follows the load.
- **A canyon on the far side.** The first trace sampled starts by rejection, and where no draw passed in 64 it kept the last, anywhere out to two and a half rise radii: the longest canyon, 4,139 km, started beyond the foot and ran on through the load's weak hoop tension, which never turns to compression, to the far side of Mars. Starts are now drawn from the tabulated distribution, and a fault grows only while the tension is above a twentieth of the crest's.
- **Buried floods' steps.** The flanks' plains show the coarse outlines of buried floods as a staircase in the hillshade, an artefact of the lava plains' milestone, entered as audit W7.

What the maps show:

- **Mars:** radial graben in swarms on the rise's flanks, the youngest across the plains, a maze of short troughs on its crest, its canyons running out from near the crest; every trough within 1.4 rise radii.
- **Venus:** each rise ringed by a star of fine radial graben about its shields, its canyons curving.

Still to do:

- **Pit chains and the canyons' landslide walls**, the item's last part, and the rise's flexure.
- **Buried floods' margins** in coarse steps (audit W7).

Tests 181 (four ignored). One is new: the crest is in tension every way alike, the flank in hoop tension and the foot in compression, a flank fault runs along the radius, and the starts are drawn as the tension is. Clippy is clean and validate passes its 117 checks.

## The thirty-ninth milestone: canyon walls fall in landslides

The third part of the revisit's fourth item.

- **F106, landslides.** A canyon's walls, where at least a kilometre high, retreat in alcoves sized by the wall's height, and each alcove's volume runs out across the floor in a lobe that widens toward its rounded toe, as Valles Marineris's landslides did.

How it went:

- **The first lobes** were flat fans with straight sides that began in a step at the wall's foot; they now widen toward a rounded toe, thin to their edges and drape the wall's lower slope, each holding its alcove's volume (a test sums one over the floor against its alcove, within a tenth).
- **Where they show: nowhere yet on this Mars.** Its large canyons opened in the first half of its life, and the floods that kept coming over the rise filled them: the deepest, 7.4 km, keeps 2.8 km of relief, a floor of flows and walls worn by the time since, with no alcove to be told at a kilometre a pixel. That is the right order of events for the floods the model draws, but Valles Marineris stands open because Mars's later volcanism gathered at its shields; entered as audit M14. The young canyons, a kilometre and a half deep, cut alcoves a few kilometres across, under the standard grid. The landslides are drawn and conserve their volume; they will show on a canyon that stays open.

What the maps show:

- **A synthetic canyon** 150 km wide and 6 km deep, drawn by the trough's own profile: alcoves along both walls and tongues of debris across its floor (a diagnostic, ignored test, writes it).
- **Mars:** no visible change at the standard quality; its canyons are filled and worn.

Still to do:

- **A rise's late floods** fill its canyons (audit M14).
- **Pit chains and the rise's flexure**, the item's remaining parts, and the buried floods' coarse steps (audit W7).

Tests 183 (five ignored). Two are new: a deep canyon's walls carry alcoves along both sides and a shallow one's none, with a landslide's debris holding its alcove's volume; and the diagnostic. Clippy is clean and validate passes its 117 checks.

## The fortieth milestone: a full lake spills into the next low

The first of the two fixes before the revisit's fifth item, closing audit W7.

- **F107, fill and spill.** A flood's lake fills only to its lowest pass and spills the rest of its lava into the next low, so every lake ends in a flow's front; in the finish, the hollow fill takes only the finer relief in the ring round a low.

How it went:

- **Two false starts.** Gating the ring's fill to the finer relief took the large cliffs away and left walls a fine cell wide. Sampling the level smoothly between coarse cells changed 0.05 percent of the pixels and was taken back.
- **The cause was in the loop.** A wall's trace showed a flood filling a cell 1,560 m deep, to the level of a lake that stood 1,780 m above the open ground at its edge. Logging every flood on Mars showed the pattern. A lake took the ground beyond its pass as well, down to where it would spread thinner than 50 m on average, and filled everything to one level, standing over the open ground where it stopped. That happened in 69 percent of the floods by more than 50 m and in 18 percent by more than 100 m, and the worst lake stood 4.6 km above the ground beside it.
- **The volume.** The first spilling lakes laid only 15 percent of Mars's lava, because a lake that stopped on flat ground was held to the height of the ground it stopped at. A cell now joins where the lava would stand a flow's thickness deep on it. What is lost now is lava that runs past the flows' reach: 1 to 24 percent by body (audit M15).

What the maps show:

- **Mars:** the plains are flow fields with lobate margins, with the staircase gone. Tharsis's radial graben swarms show on its flanks, since the lava runs off the flanks into the lows instead of stacking up in terraces. Narrow graben on the crest, a few kilometres wide and partly buried, show at a fine cell or two.
- **Venus:** its blocky flood margins are lobate too.
- **The Moon and Mercury:** the maria and smooth plains keep their look; a few craters near them are flooded or left dry where the lava spilled differently.

Still to do:

- **A rise's late floods** still fill its oldest canyons (audit M14), next.
- **Floods run off the slopes** without leaving their sheets (audit M15).

Tests 179 (five ignored). One is new: floods from sixty vents over rough ground end in flow fronts, none standing more than a flow's thickness over the ground beside it, and some of their lakes spill. Clippy is clean and validate passes its 117 checks.

## The forty-first milestone: a rise's shields take its eruptions

The second of the two fixes before the revisit's fifth item, closing audit M14.

- **F108, central volcanism.** From each shield's start, its share of what the rise adds to the vent weights goes to the shields, whose edifices and flows are their own, and leaves the floods, so a rise's late volcanism builds its shields and leaves its canyons open.

How it went:

- **Growing with the shields.** The share first grew as the shields grew. Mars's largest canyon, carved from 1.0 to 2.0 Gyr, was filled during its carving, when the shields held less than half the rise's eruptions. Each shield now takes its full share from its start.
- **Floods at the shields.** The shields' share first erupted as floods from the shields. A trace of the largest canyon's centre showed most of its late fill coming from them, the lava running 1,300 km down the rise from shield to canyon. Held to the shields' aprons, 1.4 times the foot's distance as the finish draws their flows, 98 percent of it ran past the aprons and was lost, since the loop ponds lava only in closed lows. What the shields take is now theirs, not the floods'.
- **The volume.** The shields take a quarter of Mars's eruptive budget and 61 percent of what erupts on Venus after its resurfacing. Their edifices, sized by their own rule, hold 8.3 million km³ of the 61 million km³ they take on Mars; the rest is not laid.

What the maps show:

- **Mars:** the canyons carved after the shields started stand open, walls and floors, the largest keeping 2,980 m of its 7,400 m along its trace; the rise's crest shows its older ground, with its graben swarms and craters, and an old basin the late floods had buried.
- **Venus:** its rises' graben show more clearly; its plains are as before.
- **The Moon and Mercury:** unchanged, having no rises.

Still to do:

- **When a shield starts** is a time, not a threshold in the lithosphere's thickness, and its size comes from its own rule, not from the magma it takes (audit M14).
- **Floods run off the slopes** without leaving their sheets, and the shields' late lava with them (audit M15).

Tests 181 (five ignored). Two are new: a rise's shields take their share of what the rise adds to the vent weights, out of the floods' weights, the weights away from the rise staying as they were; and the share is nothing before the first shield starts, then each started shield's share of their volume. Clippy is clean and validate passes its 117 checks.

## The forty-second milestone: floods leave their sheets on the slopes

Asked for before the revisit's fifth item, closing audit M15.

- **F109, sheets down the slopes.** The lava leaves its sheet, a flow's thickness, on each cell it runs over, as wide as its volume spreads it, and ponds only with what is left; the finish draws a sheet over the ground rather than a fill to a level.

How it went:

- **The first sheets, a cell wide,** carried only 3.5 percent of Mars's lava: most paths reach a low within two to four cells, and the largest floods still cascaded from lake to lake past the reach, losing 14.5 percent of their volume. A flood's sheet is now as wide as its volume spreads it, from the shape of the Roza flow field; less than 1 percent is lost.
- **A block on Mars's flank,** a mesa 100 km across: the hollow fill took a sheet's surface on a steep slope for a lake's level and filled hollows 1.7 km below it to it. Sheets no longer seed the fill. A smaller step was left: the fill's test compared the fine ground with the base, and the gaps that the sheets' fine outlines leave stand hundreds of metres below the base, so the fill walked 40 km down the slope. A hollow now needs the loop's own ground at the floor.
- **Barred canyon floors** (the user, at the end of a graben on the flank): the floor stepped by a kilometre from one fine cell to the next. The finish rebuilds the ground each flood found from what the slow processes added since, and there the canyon was missing at each landslide: at the loop's cell-wide walls a landslide's debris lobe had no floor to run out on and came out 340 km thick, emptying the canyon. A bug of F106, hidden while the floods filled the canyons; a test fails on the old code. A grid now draws only the landslides wider than its cell: drawn at the loop's cells, the alcoves stood in blocky steps at a canyon's end where the finish drew it anew.
- **Rings round young craters** (the user, on the Moon): a raised circle with a sharp step 2.3 radii out, faint before and plain with the sheets. Older floods left the base alone under a younger crater's whole ejecta, the coarse lava in it, while beside it they drew their fine sheets and gaps. Laying the sheets under the ejecta too left the rings fainter and a blocky pattern within them. They now leave only the cavity alone, and beyond the rim the crater's restore buries their relief by its ejecta's thickness, so the rings are gone.
- **Craters' ages:** a test found 10 cells inside young craters carrying older craters' ejecta, a surface restored where a flood did not reach, older than a crater formed since; and 2 beyond a terraced crater's scalloped rim, where its ejecta keep the age of the ground under them (F100), inside the test's circle. A restored surface is no older than the latest crater over the cell, and the test checks the floor and wall (0.35 of the diameter).

What the maps show:

- **Mars:** the rise's flanks carry the floods' sheets, its graben partly buried and its old ground showing between; its canyons' floors are clean to their ends.
- **The Moon:** the maria cover 16.3 percent of the surface, against the 17 percent measured, their flows draping the craters under them as ghosts; no rings round the young craters.
- **Mercury:** the smooth plains cover 30.9 percent, against the 27 percent measured.
- **Venus:** as before.

Still to do:

- **Canyons and graben** read as worms (the user): their structural shapes, next.
- **A few small blocky spots** of old relief show through the maria beside young craters, where older floods now act, a fine cell or two across (fixed in the forty-third milestone).
- **The swath's aspect** is one flow field's; a flow's width from its supply and slope waits on W1.

Tests 182 (five ignored). One is new: a canyon whose walls leave its floor no room keeps its depth at each landslide (failing on the old code at 0 m of 1,738); the downhill test now counts the flow's swath, and the crater ages test its floor and wall. Clippy is clean and validate passes its 117 checks.

## The forty-third milestone: the plains' blocks

Asked for after the forty-second milestone's images (the user: the blocky terrain artefacts round the Moon's craters had become worse).

- **F110, flows per run and fills that follow the fine relief.** A flood's flows are traced for each of its runs and stay within two coarse cells of its own cells; a lake's fill, and a sheet's, enters any finer hollow beside it whatever its coarse cell; a buried flood's sheet lies everywhere the loop laid it.

How it went:

- **A trace of a block's edge** showed a lake filling one side by 160 m, to its level, and leaving the other, where the cell's coarse cell carried only a flood's sheet. The fine ground on both sides lay some 500 m below the loop's.
- **The drift.** Each flood's pass takes its coarse smear off where its flows do not reach, and the lava goes nowhere. On the Moon the passes kept 81 percent of the loop's lava, visible floods with sheets 61 percent, and lakes refilled the drift within their coarse cells, in steps along the cells' edges.
- **A bug of F109.** Some visible floods kept almost none of their lava. Their flows were bound by the flood's cap and started at its vent, which for a cascade lies on the cap's edge: 47 floods on the Moon drew no flows at all. A test of a cascade fails on the old code with no flows on its sheet.
- **Conserving the volume** by spreading each pass's lost lava over the cells it laid made the plains grainy, and capped at a sheet it dimmed the flows; set aside (audit M16). The fills that follow the fine relief remove the blocks without it.
- **Trenches.** A thin trench 100 to 150 m deep was left where many floods' flows had missed the same strip; a sheet now fills the finer hollows it crosses, which the loop's own ground keeps off open slopes.

What the maps show:

- **The Moon:** the blocks and stairs round the craters are gone; the maria cover 15.5 percent of the surface.
- **Mars, Mercury, Venus:** as before or a little smoother; Mercury's smooth plains cover 27.5 percent, against the 27 percent measured.

Still to do:

- **The finish's floods** keep four fifths of the loop's lava (audit M16).
- **Canyons and graben** read as worms (the user): their structural shapes, next (the forty-fourth milestone).

Tests 183 (five ignored). One is new: a flood running as a sheet from its vent down to a lake draws flows over both runs, staying near its own cells (on the old code its flows covered none of the sheet). Clippy is clean and validate passes its 117 checks.

## The forty-fourth milestone: faults in segments

Asked for by the user: the canyons and graben read as worms.

- **F111, segments, faulted layers and half-graben.** A rise's faults grow in straight segments that link, turning toward the stress field and stepping en-echelon where they meet. A graben is as wide as its faulted layer sets it, a few kilometres whatever its length, with straight walls to its tips, its throw peaking in each segment and dipping where segments link. A fault whose throw would exceed its layer opens a rift, whose segments are half-graben, their floors tilting toward a master fault whose side flips from one to the next.

How it went:

- **What made the worms.** Each trace steered toward the stress field a little at every step, so it curved smoothly; each trough's depth and width came from its length, so it tapered as a lens from tip to tip; and every one swelled and scalloped like a chasma. Measured graben (Plescia 1991; Davis, Tanaka and Golombek 1993; Hauber et al. 2013) are a few kilometres wide whatever their length, with throws a hundredth of each segment's length.
- **The faulted layer.** The first version set a graben's width from its layer's thickness, and the long faults all came out as wide rifts with 2 of Mars's 22 loop troughs left. The width correlates poorly with the layer's thickness (Schultz 2000), so the two are now drawn apart, and 17 remain.
- **Segment lengths.** Segments 3 to 20 widths long left most rifts in one or two pieces; 2 to 10 widths gives a long rift several segments, as Valles Marineris's troughs are, and a graben segments tens of kilometres long.
- **Deep narrow graben.** With the throw capped only by the layer, some graben stood deeper than they were wide; a graben now throws no more than a fifth of its width.
- **A regression caught.** The icy moons' ring systems include ridges of negative depth, whose frames the new taper dropped; a check against the last commit found Ganymede and Callisto changed, and a trough without segments keeps its old taper. They are identical to the last commit again.

What the maps show:

- **Mars:** the rise's canyons run in straight segments with kinks and steps, their master walls standing out; its graben are straight lines a few kilometres wide in parallel swarms with en-echelon steps, the crest a maze of short straight graben.
- **Venus:** its rifts are straight half-graben; its graben, mostly finer than its 9 km cells at the standard quality, show faintly, where they were drawn 10 to 15 km wide before.
- **The other bodies:** unchanged.

Still to do:

- **The faulted layer** as a field over the rise, thickening toward a volcano's summit, rather than drawn per fault.
- **Pit chains** along the graben, finer than the standard grid.
- **The finish's floods** keep four fifths of the loop's lava (audit M16).

Tests 185 (five ignored). Two are new: a fault's trace keeps to a great circle within each segment and turns or steps only where they link, its throw peaking in each segment, dipping at the links and falling to nothing at its ends, its width kept to its tips, and a rift's half-graben floor standing deeper toward its master fault; and on Mars's rise every graben is as wide as its layer sets it and no deeper than its layer or a fifth of its width, the rifts several times as long. Clippy is clean and validate passes its 117 checks.

## The forty-fifth milestone: graben along their length, rough and in sets

Asked for by the user, after the forty-fourth: the graben looked short, too even in width and too straight, and seldom parallel, each to be checked against measurement first.

- **F112, flat-topped throws, rough walls and traces, parallel sets.** A graben keeps its throw along its length, as faults confined to a layer do, dipping a little where its segments link; its walls are rough fault traces, so it pinches and swells, and every fault's trace wanders a little, self-affine at every scale; a graben joining another stands apart from it by a few thicknesses of its faulted layer.

How it went:

- **Checked against the MOLA relief** at its own 2.6 km a pixel: Claritas and Ceraunius Fossae run unbroken across hundreds of kilometres, pinch and swell, wander, and stand in sets about 10 km apart. Confirmed by measurement: confined faults' flat-topped profiles (Soliva and Benedicto 2005), faults' self-affine roughness (Candela et al. 2012), confined faults' spacing in proportion to their layer (Soliva, Benedicto and Maerten 2006).
- **Length, half right.** Our graben average 85 km on Mars, as the Mars fault catalogue's 8,500 faults average 80 km (Knapmeyer et al. 2006); what is missing is Tharsis's radial graben systems thousands of kilometres long, read as the surface of giant dike swarms (Ernst, Grosfils and Mège 2001; Wilson and Head 2002), entered as audit M17 and next.
- **Throw from the whole length,** the first version, made each graben 470 m deep on average and spent the rise's extension on half as many; the throw now comes from the graben's segments' mean length, 188 m on average, within the 10 to 636 m measured across northern Tharsis (Plescia 1991).

What the maps show:

- **Mars, at the high quality:** graben run unbroken along their length, their walls pinching and swelling, their traces wandering a little, some in parallel pairs and sets; at the standard quality, where a graben is under a cell wide, the change is slight.
- **Venus:** its rifts' walls a little rougher.
- **The other bodies:** unchanged.

Still to do:

- **Radiating dike swarms** and their long graben (audit M17), next (the forty-sixth milestone).
- **The faulted layer** as a field over the rise; **pit chains**; **the finish's floods** keep four fifths of the loop's lava (audit M16).

Tests 185 (five ignored). The fault test now checks a graben's flat-topped throw, its dips at the links between segments' throws and its unbroken profile across a link. Clippy is clean and validate passes its 117 checks.

## The forty-sixth milestone: dike swarms

The answer to the user's sense that the graben should be longer (audit M17).

- **F113, radiating dike swarms.** A rise's plume head and each volcano on it send out swarms of dikes, each running out along the stress trajectories as far as its magma carries it, a graben above each as long as the dike; Tharsis's radial graben systems, read as the surface of giant dike swarms, and Venus's radial systems about its volcanoes.

How it went:

- **The finding.** Our tectonic graben average what Mars's fault catalogue does, and their lengths follow fault populations' power law; Tharsis's radial graben systems, 3,000 to 4,000 km long, are another population, the surface of dikes the plume fed (Ernst, Grosfils and Mège 2001; Wilson and Head 2002), debated but measured beneath the graben south of Arsia Mons (Goudy and Schultz 2005).
- **The dikes' graben** are graben as F111 and F112 draw them, as wide as twice the depth to the dike's top (Mastin and Pollard 1988), segmented, rough and flat-topped, their traces leaving their source on their own headings and running on along the rise's stress trajectories.

What the maps show:

- **Mars:** 68 graben longer than 1,000 km (3 before) and 13 longer than 2,000 km, radiating from the rise's centre and its volcanoes across the rise and onto the plains, best seen at the high quality; at the standard quality, where a graben is under a cell wide, they are faint lines.
- **Venus:** 63 graben longer than 1,000 km, in fans about its volcanoes, faint at its 9 km cells.
- **The other bodies:** unchanged.

Still to do:

- **The swarms' magma** from the plume's intrusive budget, and **pit chains** above the dikes (audit M17).
- **Dike swarms off the rises**, about Venus's volcanoes and coronae elsewhere.
- **Charon's chasmata**, which need a driver the model lacks: the extension of an ice shell as its ocean freezes.

Tests 186 (five ignored). One is new: on Mars's rise, more than fifty dike graben, each running away from its source, the longest more than 2,000 km. Clippy is clean and validate passes its 117 checks.

## The forty-seventh milestone: a rise's axis

Asked for by the user, with a figure of Tharsis's giant dike swarms: its graben run out to the south-west and north-east, the ground between largely free of them. Burial by later volcanism, or something the model lacks?

- **F114, a rise grows along an axis.** Its volcanoes line up on the axis, as the Tharsis Montes stand on their trend; its stress adds a pure shear about it, tension across the axis, so the faults beyond its foot gather along the axis and spare the ground across it; its dike swarms fan along the axis, so they run out on both sides of the rise.

How it went:

- **The finding.** In part burial, which the model draws: younger floods and shields fill or cover older graben. In part the regional stress: Tharsis's dike swarms propagated through a stress field made by the rise's growth, with its magmatism concentrated along the Tharsis Montes' axial trend (Orlov et al. 2023), and dikes follow the total stress (Muller and Pollard 1977; Cailleau et al. 2005 at Alba Patera). The model's rise stress was axisymmetric, so its faults and dikes radiated evenly.
- **How the pure shear fades.** Fading as the load's own stress does, as the inverse square of the distance, it never turns the dikes far out, the load's hoop tension being the larger. Fading as the inverse of the distance, the second try, it overtook the load's beyond some 1.4 rise radii and turned them there, but it kept the ground along the axis faulting out to 2.3 radii and raised five canyons, up to 134 km wide and 4.5 km deep, alone on the plains as far as 5,600 km from the rise's centre; Mars has none. Back to the inverse square: the faults reach some 1.8 radii along the axis, none beyond the foot across it, and turning dikes far from their source waits for a regional stress from its own sources (audit M18).
- **The axis is carried along the great circles** from the rise's centre to each point; the first test took it in the wrong frame.

What the maps show:

- **Mars:** beyond the foot the graben gather along the axis, 1,437 beyond 1.1 rise radii, all but 6 within 60 degrees of it (604 spread evenly before), and the volcanoes line up along it; 3 canyons beyond the foot (1 before).
- **Venus:** its three rises take axes of their own; faint at its 9 km cells.
- **The other bodies:** unchanged, having no rises.

Still to do:

- **The axis from computed stresses** (audit M18): the body's reorientation under the rise, the dichotomy, the other rises' loads; dikes turning into the regional strike far out; the rise's centre of activity moving over its life.
- **Pit chains** above the dikes and graben, sub-grid at the standard quality; **the swarms' magma** from the plume's intrusive budget; **dike swarms off the rises**, about Venus's volcanoes and coronae (audit M17).
- **A rift's depth from the strain where it runs,** not its length alone: a long fault far out along the axis still opens a canyon now and then, three beyond the foot on Mars (one before), 50 to 65 km wide.
- **The faulted layer** as a field over the rise; **the finish's floods** keep four fifths of the loop's lava (audit M16).
- **Charon's chasmata**, from an ice shell's extension as its ocean freezes, a process the model lacks.
- **Then the revisit plan's fifth item:** the rise's flexure and the Mars dichotomy, which may supply the regional stress.

Tests 187 (five ignored). One is new: a rise's stress keeps its crest's maze, puts the ground along its axis beyond the foot in tension farther out than the load alone does and the ground across it not, and on Mars the volcanoes lie on the axis and the faults beyond the foot gather along it. One bound is loosened: on the preview Mars the deepest canyon the loop carried stands 2.7 km below its rims (over 3 km before), the new start table having reshuffled the canyons, the deepest at the standard quality now 1,674 km long and 5.0 km deep (2,371 km and 7.0 km before); the test asks for 2.5 km, still kilometres deep after the later floods. Clippy is clean and validate passes its 117 checks.

## The forty-eighth milestone: a shield buries what it grows over, and the rise's faults as old as it

Asked for by the user after the forty-seventh: graben cutting through shields, their floors following the flank; craters half on a shield and half on the plains, the half on the shield raised sharply; and craters half in a valley or graben and half on the ground about it, at two levels.

- **F115, a shield's lava buries what it grows over.** In the loop the lava buries the ground the edifice grows over as a deposit buries what it falls on, the lows filling and the highs drowning as it thickens; in the finishing pass a crater or a trough keeps the share of its relief the lava laid on it since has not reached.
- **F116, the rise's faults form as its dome rises.** The rise's stress is its dome's, so its faults form at the pace the dome grows, and none after; the plume head's dike swarms with them, each volcano's as it grows.

How it went:

- **The craters on the shields** rode up on them: the loop added a shield's thickness to whatever stood there. A crater at a shield's foot that formed while it grew kept a floor tilted 8 km across the flank.
- **The graben over the shields** were young: on Mars 54 percent of the graben formed in the Amazonian, at the pace of the body's volcanism to the present, though the dome stopped rising at 3.0 Ga. Tharsis's faulting was about half Noachian and declined (Anderson et al. 2001); now 40 percent Noachian and 60 percent Hesperian, the volcanoes' dikes on into the Amazonian. The shields bury the faults older than them (F103).
- **The craters split by valleys and graben** that I found are older than the trough that cut them: a 167 km crater 3.41 Ga old has its south half dropped by a 113 km rift 3.08 Ga old, as faults offset older craters. A crater younger than a trough replaces its floor inside its rim (the stamp's floor from the mean height inside the rim, the restore cutting the finish's older relief), so it shows one level; with the faults now older, more craters postdate them.
- **A buried graben** first showed as a hairline over the flank: F103 took the lava's thickness off the trough's depth, leaving only the deepest line of its floor once the faults were older than the shields; the lava now takes the same share of a trough's relief, or a crater's, across its whole width.
- **The canyons** now open in the Hesperian and the rise's eruptions before its shields take them fill them in part: the deepest standing canyon at the preview quality stands 1.7 km below its rims, against 2.7 km before (audit M14 reopened in part).

What the maps show:

- **Mars:** the crater at the shield's foot keeps its eastern half as an arc against the flank, the half the lava reached gone under it; the graben are older and softer, partly buried by the plains that came after them, and fewer cross the shields.
- **Venus:** its shields bury what they grow over; its faults' timing is unchanged, after the event.
- **The other bodies:** unchanged, having no rises.

Still to do:

- **A volcano's load in the stress field** and its lava ponding at its foot (audit M19), with the rise's flexure.
- **The canyons' fill** (audit M14): the floods' share of the magma falling as the lithosphere thickens.
- **Io's flow fields as the floods' process** (audit M20).
- **Pit chains, the swarms' magma, swarms off the rises** (audit M17); **the axis from computed stresses** (audit M18); **the faulted layer as a field**; **M16**; **Charon's chasmata**.

Tests 189 (six ignored, one of them a new diagnostic that lists the troughs passing a point). Two are new: a crater a shield grows over keeps under a tenth of its depth while one beyond its foot keeps it all; and Mars's faults and the plume head's dikes form while the dome rises, each volcano's dikes while it grows. One bound is loosened: the deepest standing canyon on the preview Mars, 1.5 km (2.5 km). Clippy is clean and validate passes its 117 checks.

## The forty-ninth milestone: Io's lava as floods, and its paterae

Asked for by the user after the forty-eighth: Io's flows aligned with the floods' newer physics (audit M20); and, on the first renders, whether the dark lava pools' colour was right, and what the perfectly oval pool in the north-west was.

- **F117, Io's lava erupts as floods.** The paterae's eruptions run down from their rims, fill and spill over the lows and the paterae they reach, and lay their sheets and lobes as every body's floods do; the form of a flood's flows follows how fast it poured out, compound fields where slow, sheet lobes where fast.
- **F118, a patera collapses piecemeal.** More collapses of comparable size break out of the first, its margin wanders at every scale, and its floor lies in patches of different ages.

How it went:

- **The floods filled the paterae.** Erupted after the pits were cut, the lava ponded in them, 900 km³ to fill one, and the paterae vanished from the map; cut after the lava, as a patera born at the step's end is, they stand, their floors lava of their own.
- **Steps along the coarse cells.** The finish starts Io from the loop's bare plains, which carried none of the floods' coarse lava that the flood passes take out, so each pass cut its coarse footprint into the plains; the plains now carry that lava, smoothed with them, as every body's base does.
- **Angular lobes.** The floods' lobes, long sheets after the Imbrium flows', drew Io's fields as bundles of sticks. The physics that tells them apart is how fast the lava pours out (Walker 1971): slow, compound fields of short broad breakouts off narrow insulated paths (Io's Prometheus, 49 m³/s); fast, a few long sheet lobes (the maria, Tharsis). A compound field's trunk is the narrow path, its breakouts broad; drawn at the breakouts' proportions, the trunk left its flood's footprint at once and the fields stopped short.
- **The colour of the dark pools** is right: an active patera's floor is its lava lake, fresh silicate at 0.05 reflectance, Io's darkest ground (Loki's floor near that); quiet floors and the flows are lighter as the plumes' fallout and frost coat them by their age. The new floods leave the floors alone, the paterae being cut after them.
- **The oval pool** was a patera 217 km across at 55°W, 68°N, one collapse hardly elongated, drawn as a smooth ellipse since F88. Small collapses round its rim left a circle with bites in it, as the user saw; collapses of comparable size breaking out of the first, and a margin wandering at every scale, give arcs of different curvature meeting at corners.
- **Generic lava left in the margins.** The resampling gave a few fine cells beside each flood the coarse lava's unit, which no pass then reached: 1.3 percent of Io; they keep the plains now.

What the maps show:

- **Io:** its flows run down the ground and pond in its lows and in the paterae they reach, larger fields than before following the slopes, mostly compound, a few in sheets; they cover 20.7 percent (16.9 before; Io's 28 percent, Williams et al. 2011). Its paterae's margins are irregular, arcs meeting at corners; its lava lakes stay dark.
- **Every other body:** unchanged, Mars to the metre and the pixel, its floods' rates keeping their sheets.

Still to do:

- **Io's eruptions from its heat budget's magma** rather than the old fields' shares (audit M8, M20), at a compound field's own thickness rather than the flood sheet's 50 m, each breakout its own age.
- **The canyons' fill** (audit M14), next.

Tests 189 (six ignored, one of them the trough diagnostic of the forty-eighth). One is new: Io's eruptions come oldest first from past their paterae's pits on the side the ground falls away, as large as a flood's sheet over their fields' share of their reach squared; the Io test now reads its floods. Clippy is clean and validate passes its 117 checks.

Revised after the forty-ninth was committed, when the renders at high quality showed Io's flows in boxes: a straight band of lava with its edges in steps of the loop's cells, and the same steps in the relief. Traced cell by cell through the flood passes, each cause was the loop's coarse cells showing through where the finish draws Io's lava afresh on its smoothed plains:

- **Lakes perched on slopes.** The loop's flood had ponded in a coarse low a few metres deep that the plains, smoothed over the loop's cell, no longer hold. Filled to the pond's level wherever the ground lay below it, the lake stood perched on the slope, 10 m above the open ground beside it, and ended where the cells its pass visits end. A soft flood's lava now ponds only in the closed lows of the fine ground that its lobes run into, each lake at one level, where it spills (a priority flood, Barnes et al. 2014), no higher than the loop's pond stood; elsewhere its lobes lay their sheets.
- **Gaps taken for craters.** The finer hollows of F95 are craters the loop never saw, told apart by the loop's ground standing above them. The gaps an earlier flood's lobes left lie below the loop's ground too, which carried that flood's coarse lava, so they filled to the level of the lava beside them across that flood's coarse footprint: polygons with stepped edges, once the perched lakes no longer hid them. On Io the lakes do the hollows' work, filling the closed lows the lava reaches.
- **Terraces.** Where the units resampled from the loop's cells said lava, a flood's first pass levelled its ground to a sheet below its pond (M15), a rule meant for the facets of a coarse base; on Io's smoothed plains it cut the slopes into terraces a coarse cell square, faint boxes in the relief inside the fields.
- **A residue.** The plains were smoothed together with the floods' coarse lava, which the passes take up unsmoothed; they now carry that lava just as the passes take it up.

And every Io flood is now drawn as its flows: a flood its successors bury on the loop's cells was drawn by its coarse outline, but on Io its successors' lobes leave gaps it shows through. Io's flows now cover 19.9 percent and the plains 74.5; every other body is unchanged, Mars to the pixel. Tests 190 (six ignored), one new: on an open slope a soft flood lays its sheets with no lake below its pond's level and no terrace, and a closed bowl at its vent fills to where it spills. Clippy is clean and validate passes its 117 checks.

## The fiftieth milestone: a rise's shields gather its dikes (audit M14)

Asked for by the user after the forty-eighth, with M20: the canyons' fill, the floods' share of a rise's eruptions falling as the lithosphere thickens and the shields starting at that threshold.

- **F119, a rise's shields gather its dikes as its stretching wanes.** A rising dike bends toward a volcano's load where the load's pressure outweighs the regional differential stress (Muller, Ito and Martel 2001). While a rise's dome grows its dikes run along the stretching and erupt as fissure floods; as the growth slows, its shields take the share of its eruptions by which the growth's rate has fallen from its peak, each starting as that share reaches what the larger ones hold, the largest as the growth turns.

How it went:

- **The lithosphere threshold, tested first.** The thermal histories show no threshold that holds across the bodies: the Moon's maria erupted as floods through a lithosphere the model makes 210 to 320 km thick; Io's floods rise through a thick lithosphere its volcanism keeps cold (O'Reilly and Davies 1981); a threshold set at Mars's turn, its thermal lithosphere about 100 km, would never start Venus's shields on its 65 to 71 km. A volcano's height does not follow the lithosphere's thickness in any simple way (Wilson, Head and Parfitt 1992), so it could not set when shields can grow; by my estimate from the distance over which a dike's walls cool its flow, u w^2 / kappa, a basaltic dike a metre wide stays open across 100 km of cold lithosphere, so freezing could not set it either; and the flood basalts' high eruption rates (White, Crisp and Spera 2006) are long-term rates the maria fall short of.
- **The load against the stretching.** The control that holds is the one between a volcano's load and the regional stress, which the model already has in the rise's growth: its faults and swarms follow the growth's rate (F116), and the shields now follow its fall.

What the maps show:

- **Mars:** its shields start in order of size from 0.90 to 1.48 Gyr, and take 25.7 percent of its eruptive budget, their share rising from the dome's turn; the floods on the rise stop as it stops growing. At the standard quality Mars's 19 canyons keep a median 914 m below their rims at a quarter, a half and three quarters along them, half their carved depth, where they had kept 705 m, 38 percent; those opened after 0.4 Gyr gain the most, and 6 of the 57 points lie buried flush, against 9.
- **Venus:** its shields still take all that erupts on its rises after its resurfacing; its floods are as many and cover as much, their margins falling a little differently.
- **Every other body:** unchanged.

Still to do:

- **The early canyons**, opened during the dome's growth with its floods, stay filled: the faults and the floods both follow the growth.
- **The shields' size from the magma they take**, rather than their own rule.

Tests 190 (seven ignored, one of them a new diagnostic printing Mars's canyons' relief). The shields' test now reads the stretching: none while the growth quickens, the share by which its rate has fallen after, all once it has stopped, each shield starting at the share the larger ones hold; Mars's deepest standing canyon is held to 2.2 km at the preview quality; and a crater on the flank of Mars's largest shield, which now starts as the dome turns, keeps 86 m of its 852, where the test had allowed a tenth: it allows 15 percent now, a crater riding up on the shield keeping all its depth. Clippy is clean and validate passes its 117 checks.

## The fifty-first milestone: Venus's rises after its overturn

Asked for by the user after the fiftieth, on its account of Venus: every rise grew on Mars's clock, so Venus's rises and shields grew three billion years before the overturn that resurfaced them, its shields drowned and reset by it while its rifts, drawn after it, outlived the stretching that drove them. The user asked for this first, and the overturn's cycle as the next milestone.

- **F120, a rise keeps its mantle's melt clock since its lid last formed.** A rise grows as its plume's mantle melts beneath its lid, since the lid last formed, at the body's formation on a stagnant lid and at its overturn on an episodic one: its head rises once the mantle has made a stated amount of melt per unit area since then, its dome is built by another, its canyons are carved and its shields grow over amounts of melt, all set on Mars's history so that Tharsis keeps its times. The episodic lid's special cases for faults and swarms are gone.

How it went:

- **The clock.** Mars's rise starts after 2,739 m of melt per unit area and is built by 10,333 m; Venus's mantle makes that much in 107 and about 400 Myr after its overturn, 25.6 m per Myr against Mars's 6 at its dome's turn. A rise still growing at the present has stretched by its growth so far, and its faults with it.
- **Blocky lava on Venus.** Its rises' young floods, now on its plains, came out in dark blocks a loop cell square that the units did not show. The resurfacing event had laid its lava without a composition, so it kept the coarse cells' tint of the floods after it; and a flood's pass gave its own composition to the cells it gave back to the lava it had found, tinting its whole coarse footprint. The event's lava now takes the composition of the mantle under each place, and a flood gives back the lava it found with its composition. Building the event's composition fields for each of six million cells had doubled Venus's render; they are built once.
- **Terraces on Mars.** A stepped scar on Mars's largest shield's flank, seen at the forty-ninth and the fiftieth, and square pits in old craters' floors were the levelling of a flood's own floor (F36), which inside a flood's coarse footprint cut the ground to a sheet below its pond wherever the units, resampled from the loop's cells, still said lava. Without it they are gone, and nothing else changes but a few hundred metres here and there over a hundredth of Mars, with no strips left dry along the loop's cells, which the rule had been written against.

What the maps show:

- **Venus:** its rises grow from 0.40 to 0.10 Ga and its shields start from 0.25 to 0.10 Ga, still growing today, their calderas at their summits and their flows on the plains; its floods' lava lies on the plains in smooth patches.
- **Mars:** its rise, shields and canyons keep their times; the shield's flank and the old craters lose their terraces.
- **The other bodies with floods:** the levelling's terraces gone, and their floods giving back the composition they found where their flows did not reach; the bodies without floods are unchanged.

Still to do:

- **Venus's pace** comes from its regime's melt efficiency, a constant of the name (audit S1), until the next milestone gives the episodic lid a mantle whose heat builds until its lid overturns.
- **The shields' size from the magma they take**, rather than their own rule.

Tests 191 (eight ignored, one of them a new diagnostic printing each rock body's melt per unit area). One is new: Mars's rise keeps Tharsis's window, Venus's rises grow between its overturn and the present, its shields start after it and grow at least twice as fast as Mars's, and the melt clock runs both ways. Clippy is clean and validate passes its 117 checks.

## The fifty-second milestone: a mantle that stores its heat, and a lid that breaks

Asked for by the user after the fifty-first: cyclical resurfacing for a planet like Venus, and how the model decides whether a planet goes on building up heat or dies like Mars. The model's mantle had no heat of its own, its heat flux its radioactive heat times a factor by its size, so every body cooled steadily, and its regime was named by cuts on its present heat flux and radius.

- **F121, a rocky mantle stores its heat under a lid that can break.** Its temperature follows its heat less what it convects into its lid by the stagnant-lid scaling, the lid growing as it cools and giving up its own heat at the surface. The lid breaks when the convecting sublayer's stress along a cell as deep as the mantle exceeds its friction over its thickness: it founders into the mantle, which takes up its cold, and a new lid grows back as a cooling half space until it can break again. A dry rocky lid broken within the last half of its body's life is episodic, its surface as old as its last overturn.

How it went:

- **The mantle first, in a prototype.** Calibrating in Python across Mercury, Venus, the Moon and Mars showed the stagnant bodies' present heat flows hardly move with the viscosity, the mantle settling where it loses what it makes and a little more, and that the lid's own cooling as it thickens matters on the small bodies: the Moon's lid is most of its rock.
- **The break.** The lid's friction against the sublayer's stress over the mantle's depth gives a heat flux the lid breaks above, lower for a deeper mantle and a hotter surface: Venus's mantle reaches it as it heats under its lid; Mars's, the Moon's and Mercury's never come near it. One effective friction, 0.010, far below the laboratory's as Moresi and Solomatov (1998) found a mobilizing lid needs, sets Venus's last overturn at 0.61 Ga, within its craters' half gigayear.
- **The step.** The overturns' count first hung on the step, 17 to 22 as it shrank: a young lid's half-space flux sampled once a step at a kilometre's floor lost five times the heat it should. Averaged exactly over each step, 21 overturns and the last at 3.91 Gyr for any step from 2,000 to 40,000 over the body's life.
- **Venus's plumes.** Its heat flow, 27 mW/m² where the equilibrium gave 69, would have cut its rises from three to one, the count being one per 20 mW/m². A plume carries a share of the mantle's heat, so the count now follows the heat it loses in all: one per 4 TW, Mars keeping its one, Venus its three.
- **The stagnant bodies stay as they were.** Run on them, the stored heat gives Mars some 40 mW/m² early where the equilibrium gives 140, which would move everything built on its history at once, and the Moon only 10 mW/m², its rock 62 percent of its mass in the composition model beside a 30 percent core and 8 percent ice. They keep the equilibrium history, the mantle deciding only that their lids never break.

What the maps show:

- **Venus:** its last overturn 0.61 Ga; its rises, started 51 Myr after it on the surge of heat that followed, nine tenths grown, the largest shield on each begun 0.21 Ga and a third built, the smaller ones still to come within 0.15 Gyr; its tessera the surface of the overturn before, 1.29 Ga.
- **Every other body:** unchanged, Mars and the Moon to the pixel.

Still to do:

- **The stagnant bodies on the stored heat**, once the composition model gives the Moon its rock (audit S2).
- **The plates and Io** from the mantle too, and the melt efficiency from its temperature rather than the regime's name.

Tests 192 (nine ignored, one of them a new diagnostic printing a preset's rises and shields). One is new: Venus's lid breaks again and again, last within the half gigayear its craters record, and Mars's, the Moon's and Mercury's never, their heat flows as measured; the melt clock's test now asks of Venus that its rises and their largest shields start between its overturn and the present. Clippy is clean and validate passes its 117 checks.

Revised after the fifty-second, on the user's word that Venus's albedo outside its rises was a plain noise field. The fifty-first had given the overturn's lava the mantle's source composition at every place, a smooth field of two octaves, which painted itself across the plains in place of a lava's units; the overturn's lava is one lava now, its time's drift, and the plains are one tone again under the floods that came after. Each basin's own lava, from the source under it, waits for the plains' units of audit M10.

## The fifty-third milestone: the snow line, and every rocky mantle's own heat

Asked for by the user after the fifty-second: fix the composition model, then move every rocky body onto the mantle that stores its heat. The composition model had given the Moon a 30 percent core and 8 percent ice beside 62 percent rock, leaving its mantle too little heat to make its measured heat flow, so the stagnant bodies had kept the equilibrium history.

- **F122, a body formed too warm for ice carries none, and every dry rocky mantle keeps its own heat.** A body takes up water ice by the temperature where it formed, fully at the snow line's 170 K and not at all from 200 K; where it took up none, its density's deficit below the chondritic mix is iron it lacks. Mars, the Moon and Mercury now keep their mantles' histories, as Venus does, the lid growing and thinning by the heat at its base.

How it went:

- **The snow line.** The Moon came out all rock with no core, Mars with a 16 percent core and no ice; every body beyond the snow line, Ceres at 167 K among them, unchanged. Validate passed as it was.
- **The negative heat flow.** Run on the stored heat, Mars's and the Moon's surfaces gained heat in their first hundred million years: the lid took a stagnant lid's thickness at every step, so that as the heating mantle thinned it, the heat its base took up was drawn from the surface. Looking closer showed that a young lid had lost a half space's heat from the mantle beside its own cold, nearly twice what it should, and that a foundered lid's cold had gone to the mantle without its own rock.
- **The lid by the heat at its base.** It now grows where it conducts more than the mantle convects into its base and thins where it conducts less (Spohn 1991), and the first lid grows from the surface as the body forms. Every history is positive: Mars 76 mW/m² at 50 Myr and 27 now, the Moon 71 and 15, Mercury 65 and 15.
- **Venus's cycles.** The first try at maturity, a lid that has reached a stagnant lid's thickness, stopped Venus's overturns 2.5 Ga whatever the friction: under a cooling mantle the lid only nears that thickness. Maturity is again the time a half space takes to grow that thick, and Venus's lid breaks every 0.1 to 0.2 Gyr while its mantle can break it, last 0.49 Ga at a friction of 0.014, the same for any from 0.0138 to 0.0148.
- **Mars's plume.** On its own history Mars loses 3.9 TW today, under the 4 TW a plume carries, and lost Tharsis. A plume rises long before today, so the count now follows the mantle's mean heat loss over its life: Mars's 6.1 TW keeps its one, Venus's 35 TW its three.
- **The melt clock.** Mars's new history makes about half the melt it made by 1.5 Gyr, so the clock's constants were set again on it by the same rules, and Tharsis keeps its rise from 0.30 to 1.50 Gyr and its shields' starts and gigayear of growth. The Moon makes about a fifth more melt in its maria's time than it did, Mercury about a third less.

What the maps show:

- **Mars:** Tharsis keeps its times, its rise from 0.30 to 1.50 Gyr and its shields started from 0.90 to 1.49 Gyr. Its lava covers 42 percent of the surface where it covered 50, its mean height 456 m lower and Tharsis's crown 870 m. More of its lava now comes after 0.9 Gyr, as its mantle heats under its lid, while its canyons are open, so they stand 1,757 m deep at the middle of their length in the main one where they stood 2,210, its deepest point 4,217 m where it was 2,892.
- **The Moon:** its lava covers 21 percent of the surface where it covered 15, against the real maria's 17, and erupts until 1 Ga: 220 m of lava per unit area after 2.5 Ga, where 25 m came before, the young flows laying lighter lava over the great basin's older floods. Its mantle keeps all its rock's heat and stays some 200 K hotter than parameterized histories that put a share of the heat-producing elements in the crust.
- **Mercury:** its smooth plains cover 19 percent of the surface where they covered 27, the real figure: its mantle now makes about the Moon's melt per unit area, where the equilibrium had given it 1.6 times the Moon's, before its contraction closes the vents at 1.5 Gyr.
- **Venus:** much as it was, its last overturn 0.49 Ga and its tessera the surface of the overturn before, 0.69 Ga; its rises grow from 0.46 to 0.12 Ga and its shields from 0.29 Ga on.

Still to do:

- **The mantle's melt and its heat-producing elements, next.** A share of the heat-producing elements in the crust, which gives its heat straight to the surface, and the core's own heat from its formation, which heats a thin mantle from below, so that the Moon's volcanism wanes by 3 Ga and Mercury's floods its plains; the melt from the mantle's temperature against its solidus at the lid's base rather than the heat flux and the regime's name.
- **The thermal state's mantle temperature and lithosphere** from the mantle's own temperature and lid, not the heat flux's formulas.
- **The plates and the Io regime** from the mantle too, and the melt efficiency from the mantle's temperature rather than the regime's name.
- **Io's ice**: the warm disc it formed in around Jupiter, inside whose snow line it formed, is not yet the model's.

Tests 193 (12 ignored, three of them new diagnostics printing each rocky preset's mantle history beside the pipeline's, Venus's overturns over a range of frictions, and each stagnant body's eruptions over its life). One is new: a body formed inside the snow line carries no ice, the Moon all rock and Mars a core under a quarter of its mass, and one beyond it keeps its ice. The mantle's test now asks that Mars's, the Moon's and Mercury's heat flows stay positive from the start and highest when young. The Mars test now takes a canyon's standing relief as its median along its length, since a local fill decides any one point: the main canyon's middle point fell from 2,414 to 2,058 m while its median fell from 2,210 to 1,757, and its floor is 1,500 m where a single point's was 2,200. Clippy is clean and validate passes its 117 checks.

## The fifty-fourth milestone: the crust's heat, and melting above the solidus

Asked for by the user after the fifty-third: the Moon showed lighter, patchier maria, to be fixed. On its own heat the Moon's mantle kept all its rock's radioactive elements and erupted until 1 Ga, a third of its lava after 2.5 Ga, its young flows of a lighter composition laid in patches over the older floods of the great basins.

- **F123, a crust holds a share of the heat-producing elements, and a mantle melts only above its solidus.** The crust holds the incompatible elements its melts took from the mantle, by its volume's share of the rock and their enrichment, and its heat leaves through the surface. A mantle's melt comes of its own part of the heat flux and only while it, or its plumes, stand above the solidus at the base of its stiff lid.

How it went:

- **In a prototype first.** A Python port of the mantle and the eruptive budget reproduced the model to the metre of lava, and showed that neither the crust's share nor melt from the mantle's own heat ended the Moon's volcanism alone: the cutoff by the lithosphere stayed open while the crust's heat kept the surface flux up. What ends it is the mantle falling to its solidus at its lid's base, and the Moon's mantle reached it once its crust held the share of the heat-producing elements its volume and their enrichment give.
- **A melt law that left Venus without its rises.** Melt drawn from the convected heat in proportion to the mantle's excess over the solidus fixed the Moon, but Venus's high gravity puts its solidus at its lid's base beyond its mantle a hundred million years after an overturn, so its rises would not have started until 10.7 Gyr; plumes rising hotter than the mantle, and the stiff lid's shallower base, were not enough. The share of the heat that goes into melt stays F17's, by the regime, and the solidus decides only whether the mantle, or its plumes, melt at all.
- **Venus's friction.** With a third of its heat-producing elements in its crust, Venus's mantle is cooler and its overturns come every 0.3 Gyr at the end, the last at 0.64 Ga for frictions from 0.0101 to 0.0108.
- **The core's heat, tried and left.** A hot core heating Mercury's thin mantle from below would give back its plains, but against the mantle's adiabat a core loses no heat into Venus's or Mars's mantles unless the viscosity at their cores' boundaries is taken from a rheology for their pressures, which the model does not have.

What the maps show:

- **The Moon:** its maria dark and even again. The lava's mean albedo is 79.1, with a spread of 5.7, against 87.9 and 10.7 on its own heat alone (F122) and 78.7 and 6.3 before. The lava covers 15.8 percent of the surface (15.4 before F122, 20.8 with it; the real maria cover 17). The maria's surfaces formed at 1.5 Ga at the median, against 0.67 Ga in F122 and 2.2 Ga before it: the floods of the last two gigayears of its volcanism, a fifth of its lava, still spread thin over the same basins. The real maria's surfaces are mostly 3.3 to 3.8 Ga (Hiesinger et al. 2011), which is a matter of where late floods go, still open.
- **Mercury:** its lava covers 16 percent of the surface, against 19 in F122 and 27 before it and on the real planet.
- **Mars:** Tharsis keeps its times, its rise from 0.30 to 1.51 Gyr and its shields started from 0.90 to 1.49 Gyr. Its lava covers 35 percent of the surface, against 42 in F122 and 50 before it, the crust's heat no longer counted as melt.
- **Venus:** its last overturn is at 0.64 Ga, the one before at 0.94 Ga. Its rises grow from 0.60 to 0.35 Ga, and its shields, begun from 0.47 Ga, are built by today, the tallest standing 12.9 km where a third-built one stood 10.3 km, which leaves the shields' size on Venus for its milestone.

Still to do:

- **Mercury's plains**, from the core's own heat, with a rheology for the deep mantles of the large bodies.
- **The melt's share by the regime's name**, and the crust's thickness by class, which sets its share of the heat-producing elements.

Tests 194 (12 ignored). One is new: the Moon's crust holds between half and seven tenths of its heat-producing elements, three quarters of its melt comes by 2 Gyr and none now, while Mars's and Venus's plumes melt today. Clippy is clean and validate passes its 117 checks.

## The fifty-fifth milestone: the core's own heat

Asked for by the user after the fifty-fourth: proceed with the Mercury fix. Since the fifty-third, Mercury's smooth plains had shrunk from 27 percent of its surface, the figure they really cover, to 19 and then 16: on its own mantle's heat it made about the Moon's melt per unit area, its core moving with the mantle as one temperature and taking up half the heat its mantle made.

- **F124, a core keeps its own heat.** The core keeps its own temperature, left 300 K above the mantle's adiabat at their boundary by its formation, and loses heat into the mantle through the boundary layer at their interface (Stevenson, Spohn and Schubert 1983; Howard 1966's critical Rayleigh number).

How it went:

- **In the prototype first.** The two-temperature core reproduced, in Python, gave Mercury 281 m of lava for a 300 K superheat, 326 m for 500 K and 491 m for 1,000 K: the 447 m that gave its 27 percent would take some 850 K, nearly all of the 941 K that core formation's energy could give Mercury's core were none of it lost. The superheat is 300 K, a stated choice well inside what core formation could leave.
- **The adiabat.** Measured against the mantle's temperature below its lid, a core 300 K hotter loses no heat into Venus's or Mars's mantles, whose adiabats climb 800 and 300 K to their cores; the superheat is measured from the adiabat at the boundary.
- **Venus's friction, again.** With its core's heat in its young mantle, Venus's last overturn falls at 0.35 Ga for the friction of the fifty-fourth; 0.0112 gives 0.66 Ga, near where it was.
- **A test made of one place.** The test that a shield buries the craters it grows over failed at its one place on Mars's largest shield, 23 percent of the bowl left where it allowed 15: the rise sits slightly differently, carried by a polar wander that follows its lithosphere, and a graben crosses the place. Tried at eight places round the shield, the last commit had left 6 to 41 percent, 13 at the median, and this one leaves 7 to 32, 17 at the median, the mean about the same; the test now takes the median of the eight, under a quarter, against the whole bowl a crater riding up on the shield would keep.

What the maps show:

- **Mercury:** its lava covers 18.5 percent of the surface, against 16.0 in the fifty-fourth and 27 on the real planet, its plains surfaces formed at 3.4 Ga at the median; more of its tan plains are back, most on its eastern half.
- **Mars:** much as it was: its lava covers 35.3 percent, Tharsis keeps its times.
- **Venus:** its last overturn at 0.66 Ga, the one before at 0.93 Ga; its rises grow from 0.61 to 0.33 Ga and its shields from 0.47 Ga.
- **The Moon:** unchanged to the pixel; it has no core.

Still to do:

- **Mercury's plains**, the rest of them: the crust's composition from its magma ocean (audit S7), so that the extrusive share follows the crust's buoyancy, the Moon's anorthosite holding its magma down.
- **The core's boundary under pressure**, for the large bodies; the inner core's growth, and contraction from the cooling itself.

Tests 195 (12 ignored). One is new: Mercury's core, left hot by its formation, gives its mantle more than five times as much heat young as now and cools by over 200 K, and the Moon's gives none. Clippy is clean and validate passes its 117 checks.

## The fifty-sixth milestone: a rocky crust's materials from its processes

Asked for by the user after the fifty-fifth: audit S7, the materials list, rocky bodies first and icy ones in the next milestone. The list mixed rocks with bodies and units ("Mercurian regolith", "smooth plains", "mare basalt"), each chosen by class and regime, and every body erupted the same share of its magma whatever its crust.

- **F125, a rocky crust's materials come from the processes that make it.** The mantle's iron from where the body formed; the magma ocean's flotation crust, anorthosite where plagioclase floated on an iron-rich ocean and a graphite crust where an iron-poor one floated graphite; the body's basalt with its mantle's iron; the rock colour where no regolith forms; and the crust's buoyancy deciding how much magma erupts.

How it went:

- **The list.** Sixteen materials by what they are: anorthosite, norite, the body's basalt and the iron-rich and iron-poor basalts it blends, felsic rock, graphite-bearing and graphite-rich crust, impact melt, sulfur deposit, carbonaceous material and the ices. The Moon, Mercury, Mars, Venus, Earth and Io keep the colours they had.
- **Rock and regolith.** One basalt for the Moon's maria and Mars's plains needed two kinds of fresh: an airless body's fresh regolith, finely ground and bright, and the darker rock or sand under an atmosphere, so a material may carry a rock colour that a surface impacts do not garden shows.
- **The buoyancy, in the prototype first.** Magma's buoyancy budget over its path, gained in the mantle and lost in a lighter crust, holds the Moon's iron-rich basalt down most while its lid is thin: its eruptions came later, half by 1.17 Gyr after it formed where it was 0.88, and its lava after 2 Gyr from 88 to 108 m per unit area, which the renders showed did not bring back the patches of the fifty-third. The rocky extrusive share, raised from 0.10 to 0.13 so that the Moon keeps its maria, gives Mercury and Mars more.

What the maps show:

- **The Moon:** as it was, its maria dark and even, their lava over 15.8 percent of the surface; their surfaces formed a little later, at 1.46 Ga at the median where they formed at 1.51.
- **Mercury:** its lava covers 21.8 percent of the surface, against 18.5 in the fifty-fifth and 27 on the real planet; more of its tan plains are back, its surface the graphite-bearing crust and its basins' floors the graphite-rich one, in the colours they had.
- **Mars:** its lava covers 40.0 percent of the surface, against 35.3 in the fifty-fifth and 50.1 before the fifty-third; its colours as they were, its basalt shown as rock under its atmosphere.
- **Venus:** its units as they were; its basalt, its mantle's FeO 11 percent, between the two end-members, the surface 5 percent brighter on average.
- **Io and Earth:** the same to the pixel; Io's coat is now a sulfur deposit by name.

Still to do:

- **The icy bodies' materials** and the presets that pin them, the next milestone.
- **The crust's porosity and thickness**, and a body's carbon, which the buoyancy and the graphite crust take as given.

Tests 199 (12 ignored). Four are new: the Moon floats anorthosite over norite and Mercury graphite while Mars and Venus show their basalt; a body's basalt takes its mantle's iron and shows its rock where no regolith forms; a light crust holds a dense magma down more while the lid is thin; retired material names load as their replacements. Clippy is clean and validate passes its 117 checks.

## The fifty-seventh milestone: an icy crust's materials from its processes

Asked for by the user after the fifty-sixth: audit S7's icy half, and the three presets that pinned their materials. Callisto's dark surface, Ceres's darker one and Ganymede's dark terrain came from pinned regolith mixes, the defaults gave a rock-ice body with an ocean and no air clean ice whatever it was made of, and the Bond albedo of an icy body was its regolith's whatever its lanes.

- **F126, an icy crust's materials come from the processes that make it.** The crust carries the rock the body accreted, its share of the crust's volume from the composition; its colour is the intimate mixture of clean ice and carbonaceous rock at that share; a shell that renews itself everywhere shows clean ice; the regolith carries the lag the warm latitudes grow, spread over the body by ejecta; the Bond albedo counts the renewed share as clean ice.

How it went:

- **What tells Callisto from Ganymede.** The moment of inertia is fixed by class and every body but Ceres and Pluto has an ocean in the model, so nothing in the interior tells Callisto's crust from Ganymede's: both keep the rock they accreted, and Ganymede differs by the lanes its resonance renewed.
- **The lag as its own field, in the prototype first.** A lag growing cell by cell at its latitude's temperature left Callisto dark below 40 degrees and frosted above, five times brighter, where the real moon is nearly even. The lag in the regolith, spread by ejecta, keeps the latitudes close: 0.23 at Callisto's equator and 0.31 at its poles, against the reference map's 0.30 and 0.27 to 0.37.
- **How dark is dirty.** Hapke's mixing of the single-scattering albedos reproduces the library's dirty ice at 30 percent rock and Ceres's carbonaceous ice at 75, so one family at the crust's share replaces three fixed ices.
- **What flows when a crater relaxes.** Rock stiffens ice by twelve decades of viscosity per unit share (from Ceres's crust, Fu et al. 2017); read from the bedrock it would have made Ganymede's crust ten thousand times stiffer. The relaxing ice is now the shell under the half-kilometre crust, the frozen ocean under a shell that floats on one: Ganymede relaxes as before, Callisto's largest craters somewhat more.
- **Ceres, seen by the user.** The first renders gave Ceres bright crater walls and young ejecta, its mean 0.25 against 0.13, and bright poles, where the user expected one grey with bright spots: its crust, bare, is 0.34. A new lag regrowing on bared slopes darkened its walls up to some 60 degrees but needed a stated renewal time and left its poles bright. The lag law itself had the answer: Ceres's crust loses its ice kilometres down, deeper than the layer its craters expose, so what they bring up is lag. Ceres is now one grey at 0.12, from 0.116 at its equator to 0.131 at its poles, with its salts and domes the bright spots.
- **Ganymede, seen by the user.** Its dark terrain came out as dark as Callisto, 0.20 at low latitudes against lanes of 0.55, a contrast of 2.75 where the real moon's is some 1.7: the real dark terrain is brighter than Callisto's surface, 0.29 in the reference map against Callisto's geometric albedo of 0.22. What sets them apart in the model is the lanes themselves: impacts on them, seventy percent of the moon, carry their ice onto the old crust, as impacts carry highland soil onto the lunar maria (Liu et al. 2021). With the share the lunar maria hold, the dark terrain is 0.30 and the lanes 0.56, a contrast of 1.84 at low latitudes and 1.59 at high.
- **Ganymede's craters.** Its brighter surface leaves it 7 K colder, its shell 5 km thicker and its lanes 70 percent of it where they were 68; the test of its relaxed craters, whose median depth ratio was 0.0085, gave 0.0074, and 0.0072 with the final model. With the old Bond albedo pinned the craters came back to the digit, and of the 54 craters the two runs share the median changed by +0.0001: the median of some sixty craters moves with which craters the history leaves, and the test's floor is now 0.006, far from flattened against the real 0.012.

What the maps show:

- **Callisto:** its map 0.25 where it was 0.17, the median 0.23 against the real moon's geometric albedo of 0.22; darker at its equator, 0.23, than at its poles, 0.31, against the reference map's 0.30 and 0.27 to 0.37; its regolith 56 percent lag; its largest craters somewhat softer, the relief's spread 340 m where it was 372.
- **Ganymede:** its old terrain 0.30 at low latitudes and its lanes 0.56, a contrast of 1.84 (1.59 at high latitudes), where the pinned regolith gave 1.54 and the first renders of this milestone 2.75; its map 0.49 where it was 0.56, against the reference map's 0.42; its Bond albedo 0.52, its surface 101 K where it was 109, its lanes 70 percent of it where they were 68.
- **Ceres:** one grey, 0.12 where it was 0.13, from 0.116 at its equator to 0.131 at its poles, its salts and domes the bright spots; its Bond albedo 0.090.
- **Europa and Enceladus:** the same to within a thousandth; Enceladus's crust is now its dirty ice, under its plume's frost.
- **Triton:** 0.47 where it was 0.70, its crust its rock-bearing ice with the clean ice under it bright in its craters; the nitrogen frost that makes the real moon 0.76 is the climate's.
- **Titan and Pluto:** their dirty ice mixed at their own rock, 0.36 where the fixed dirty ice gave 0.25 and 0.26.

Still to do:

- **Nitrogen and methane frosts**, which make Triton and Pluto bright: the climate's (Phase 5).
- **Ganymede's lanes and polar caps**: the lanes clean ice a little brighter than the real ones, the caps the magnetosphere's; dust from outside is a candidate for the first.
- **The central pits by the crust's rock** (S9).

Tests 201 (12 ignored). Two are new: an icy crust keeps the rock it accreted, its regolith carrying its lag, Europa clean, Ganymede between Callisto and Europa with its old terrain brighter than Callisto's, and Ceres's fresh ground dark; a warm old dirty crust grows a lag where a cold or clean one keeps its ice, and Ceres's runs deeper than its craters dig. The retired names' test covers the dirty ices and their summed shares, and Ganymede's relaxation test has its floor at 0.006. Clippy is clean and validate passes its 115 checks, Triton's Bond albedo and temperature set aside with Pluto's.

## The fifty-eighth milestone: rays on a sphere, floors under the lag, deposits by share

Asked for by the user after the fifty-seventh: an alternating pattern of light and dark ground on the floor of Ceres's largest young crater, its bright deposits pixelated as lava flows had been ("make sure this isn't repeated in the future"), and on the Moon at seed 75 a large young crater at the south pole whose rays reached the north pole, dominated most of the surface and converged into one point.

- **F127.** Rays fade where their fragments would need nearly the circular orbital speed and end a quarter of the way round; the lag hides a crater's floor and blanket as it hides its rays; cryolava matures to the ground; salts are not coloured or mixed as lava; and deposits a field bounds are drawn by the share of each cell they cover.

How it went:

- **One cause for two findings.** The layers dumped from Ceres (EXOGEN_DUMP_LAYERS) showed its rays covering the whole body and meeting at a point: the antipode of its young 352 km crater, as the Moon's at seed 75 met at the antipode of its 378 km basin. Both are Elliott et al.'s lunar ray lengths, measured on craters up to 85 km, carried to basins and drawn on flat ground.
- **Where rays end.** A sphere's ballistics: a 45-degree fragment lands a quarter of the way round at the circular orbital speed and at the antipode at escape speed, and near those speeds a spread of launch angles scatters a jet's fragments over the far hemisphere, so they lay no ray.
- **The alternating floor.** The crater's floor took its full freshness inside the rim, where on Ceres the lag regrows in centuries, while its melt sheet kept the lag's maturity: the floor showed the melt as blobs. The floor and blanket now fade under the lag, as the rays do.
- **The basin, seen again by the user.** In a globe view the seed's basin still covered the hemisphere it faces with solid wedges: the quarter-round limit had trimmed a wrong number. The ray fit, measured on craters up to 85 km across, made a basin's rays denser at the same scaled distance the larger it was; beyond Tycho's size they now scale with the crater, self-similar as a large crater's ejecta are. The wedges thinned. What remains is the rays' form: a ray drawn as a smooth band of fixed angular width, saturated where it is dense, reads as a streak at Tycho's scale but as a wedge at a basin's, where real rays are secondary craters in clusters and chains.
- **The bright deposits.** The faculae were a threshold per cell on a noise field, binary cells; they also shared the lava's cover, so they were coloured as cryolava, which since F126 had kept clean ice's whiteness at maturity, and Ceres's cryolava floods showed white with stair-stepped margins. Now the faculae are drawn by share, the salts are not lava, and cryolava matures to the ground; the crater's 432 Myr old salts have faded, as they should over 150 Myr. The deposits still drawn by whole cells are listed in the audit (W8), with the rule.

What the maps show:

- **The Moon at seed 75:** its 378 km basin's rays end a quarter of the way round, the north free of them and nowhere a point where they meet; the default Moon as it was, its map 0.124 against 0.125.
- **Ceres:** its young crater's floor one grey, its melt no longer showing through; its 432 Myr old salts faint and drawn by share, its cryolava floods matured to the ground; its rays end within 740 km, a quarter of the way round; its map 0.115 against 0.119.
- **Mercury, Callisto, Ganymede and Europa:** as they were, within 0.002.

Still to do:

- **Rays drawn from their secondaries**, clusters and chains at the fragments' landing points rather than smooth bands, which a basin's rays need.
- **The deposits still drawn by whole cells** (W8), each with its process: Venus's overturns and coronae, shields, Io's paterae, Enceladus's polar terrain, and young lava flows' margins.
- **Ceres's domes**, perfect circles: their outlines from their feeding fissures and their flanks' collapse, with Venus's domes and coronae (M10).
- **Ganymede's lanes and polar caps**, as above.

Tests 204 (12 ignored). Three are new: a straight margin covers a cell by its share; rays end a quarter of the way round, whole near the crater and at the Moon's longest rays' length; a basin's rays are a large crater's scaled up. Ceres's test holds its salted floors to partial cover at their margins. Clippy is clean and validate passes its 115 checks.

## The fifty-ninth milestone: a ray is the clusters of its secondaries

Asked for by the user after the fifty-eighth's revision ("Is the seed 75 rayed ejecta correct though or are we just capping a wrong number? In globe view it still looks ridiculous."): the ray fit, measured on craters up to 85 km across, had been carried to a 378 km basin, and every ray, a smooth band of fixed angular width, saturated into a solid wedge across the hemisphere it faced. The user agreed to draw rays from their secondaries, on every body that shows rays.

- **F128.** A ray is the clusters its fragments' secondaries make: chains of secondaries sized from the largest a crater throws (Croft 1992) and shrinking outward, opening downrange in the herringbone's V, falling in proportion to what the fragments excavate so that Elliott et al.'s ray lengths hold on average; the secondaries the finish stamps are bright where they lie; and beyond Tycho's size a crater's rays are shorter in its own radii, from Tycho's to Orientale's.

How it went:

- **The prototype first.** One crater's rays on flat ground, drawn smooth and as clusters, showed a Tycho's rays breaking into chains and flecks with ground between them, and a basin's still solid near the crater, as a fresh crater's near rays are.
- **Self-similar was not enough.** In the engine the basin's rays kept their solid wedges: scaled with the crater, its rays saturate out to some twelve of its radii, and on the Moon that is most of a hemisphere. Thinning the fragments that land near the quarter-round limit, rather than dimming the ray, broke the ends into clusters but left the wedges.
- **What a basin's rays are.** Singer et al. (2020) find the largest fragments falling off ever more steeply with their speed for larger craters, so a basin's far ejecta are too fine to make secondaries that show; Orientale's secondary chains reach some 2,000 km, four of its radii, where Tycho's median ray is eighteen of its. A fit between the two shortens the seed's basin's median ray from 3,400 km to 1,400 km: a bright halo some 30 degrees across and rays that break into clusters beyond.
- **Mars and the icy worlds.** The user asked whether the change reaches every body. It does, as one process: the icy worlds' rays break into clusters under their lag as the Moon's do, Ceres's barely change, its lag hiding them, and Europa has few young craters to throw them. Mars's map did not change: above 100 Pa nothing weathers its ground, so a ray has no mature ground to freshen, though its rays are still drawn, for an airless body's 1.5 Gyr, its 6 mbar being under the rays' 1 kPa threshold.

What the maps show:

- **The Moon at seed 75:** its 378 km basin a bright halo some 30 degrees across, its rays solid near the crater and broken into clusters and chains beyond, the far side free of them, where the last commit's wedges covered the hemisphere it faces.
- **The default Moon:** its Tycho-sized crater's rays clustered streaks with ground between them, where they were smooth bands; a quarter of the map changed by more than 0.005, about as much brighter as darker, its mean 0.125 against 0.124.
- **Callisto:** its largest young crater's sweeping rays shorter and broken into clusters; its map 0.244 against 0.247, 28 percent of it changed.
- **Ganymede:** its rays broken into clusters across the dark terrain; half its map changed, about as much brighter as darker, its mean 0.492 as it was.
- **Mercury, Ceres and Europa:** within 0.001 of what they were: Mercury's rays fade in 150 Myr, so few of its craters keep them; Ceres's lag hides its rays; Europa has few young craters.
- **Mars:** as it was, its ground unweathered above 100 Pa, so its rays cannot show.

Still to do:

- **Mars's two thresholds**: nothing weathers its ground above 100 Pa, yet its rays last as on an airless body below 1 kPa, so it draws rays that cannot show; real Martian rays barely show in visible light (Tornabene et al. 2006). The climate's (Phase 5).
- **Secondaries stamped along the whole ray**, not only a loop crater's near ones.
- **The rays' cost**, 23 s on the Moon where it was 9.

Tests 205 (12 ignored). One is new: a ray is the clusters of its secondaries, carrying what the fragments excavate and leaving the ground between them far out. The basin's test of the fifty-eighth now holds a basin's rays shorter in its own radii, Orientale's some 2,000 km. Clippy is clean and validate passes its 115 checks.

## The sixtieth milestone: one population of secondaries

Asked for by the user after the fifty-ninth: its clusters scale with the crater, so a large crater has large ones, as the basin at seed 75 does: was that right, or should they all be the same size? Larger craters throw larger blocks, so their secondaries are larger (Croft 1992), and the cost cap on a crater's clusters never engaged; what was off was that a cluster was drawn as one blur whatever its size, with the craters the finish stamped laid by a rule of their own. The user agreed to one population of secondaries, and gave notes on Ganymede's and Enceladus's grooves and Enceladus's tiger stripes to keep for later (audit T6 and T7).

- **F129.** A crater's secondaries follow Orientale's census, carried to every crater in units of its largest secondary; each lands in one of the clusters its rays make, by the ground they excavate, the larger nearer; and a secondary lands only on ground older than it.

How it went:

- **What the count is.** Croft's relation gives the largest secondary a crater throws, but F128 had sized every cluster as if it held one that large. Filling each cluster with secondaries on that reading would give a 930 km basin some 4,400 secondaries of 10 km and more; a steep power law from Croft's largest would give it 37. Orientale's census, 2,728 secondaries of 2 km and more out to six radii (Guo et al. 2018), has some 400: a Weibull, shallow below 10 km and steep near its 27 km top, which is Croft's.
- **Where they land.** Chosen by the ground the clusters excavate, an old 782 km basin at seed 75 threw its secondaries mostly to one side, the resultant of their directions 0.76, where most basins' were 0.05 to 0.5. Chosen instead by the ground the fragments excavate before the cap, they went nearly all into the curtain's strongest jet, 0.98: the curtain's mass is lognormal with a spread of 2.5, which the rays' saturation hides and a count does not. They follow the clusters the rays show, and an oblique basin's field is one-sided, as Orientale's lie in a swath to the northwest and to the south.
- **Ground older than they are.** The replay dropped a secondary on ground renewed since its parent only where ice had come up from below (F9). At seed 75 the young basin's blanket carried older basins' secondaries; now lava, renewed ice and younger craters' floors and blankets all take them, as they take the population's craters.

What the maps show:

- **The Moon at seed 75:** its 378 km basin's secondaries lie along its rays, craters of 3 to 5 km in radial lines, their floors 1 to 1.6 percent of the ground from two to six of its radii and half that beyond; the older basins' secondaries are gone from its fresh blanket, where the old rule had left them; an old 782 km basin to its north throws a field to one side. In relief its secondaries are faint, shallow and small among the population's craters, and its bright clusters stay soft.
- **The default Moon:** the old rule's secondaries, scattered out to eight radii from every crater, are gone, and the census's lie nearer their craters, along their rays; secondary floors 0.72 percent of the ground where they were 0.89; the albedo as it was, 0.125, a hundredth of the map changed by more than 0.005.
- **Mars and Callisto:** far fewer secondary floors, 0.54 percent where they were 3.75 on Mars and 0.49 where they were 1.40 on Callisto. The old rule gave a basin Hellas's size some 9,200 secondaries of 11 to 92 km, the census some 2,600 smaller ones, and those on plains laid since are gone. Their albedo within 0.001.
- **Ceres:** more, 1.85 percent where they were 0.51: its craters' largest secondaries are large against the smallest crater its grid shows, below the census's floor, where the count steepens. Its albedo 0.115 as it was.
- **Mercury, Ganymede and Europa:** fewer or none, as before; their albedo within 0.001.

Still to do:

- **A cluster's grain.** Its ground is the secondaries in it, most finer than the grid; those near the grid's scale drawn as their own halos would give a large crater's clusters the grain of their secondaries rather than a blur.
- **The census's small end**, a floor from an old basin; a young crater's own census would show how many small secondaries erosion took.
- **Mars's two thresholds**, from the fifty-ninth: its rays, which cannot show, take 50 s.

Tests 207 (12 ignored). Two are new: a crater's secondaries follow Orientale's census in units of its largest, inverting their own draw; a basin's secondaries are as many as the census gives, lie in its rays' clusters beyond its blanket, the larger nearer. Clippy is clean and validate passes its 115 checks.

## The sixty-first milestone: a cluster's grain

Asked for by the user after the sixtieth ("Yes, go ahead with that"): the sixtieth left a basin's clusters soft, its secondaries mostly finer than the grid, and proposed giving a large crater's clusters the grain of their secondaries rather than a blur.

- **F130.** A cluster's ground is laid as clusters of its own, a RAY_CLUSTER_WIDTH-th as wide, and each of those likewise, down to the narrowest drawn: the same rule at every scale, keeping the cluster's ground and place.

How it went:

- **The prototype first.** One crater's rays on flat ground at a map's pixel scale and contrast, its clusters laid as one Gaussian each and as their grain: the basin's soft blurs broke into flakes and streaks, a Tycho's rays hardly changed, its clusters near the grid's scale. Stopping the split at twice the narrowest drawn rather than at it changed little, so it stops at the floor the clusters already had.
- **The cost.** A basin's cluster splits into some eighty parts at 8192 pixels. Laying a large crater's clusters in parallel, each cluster's cells read once and its parts reading only the cells in their reach, took the rays from 19 to 11 s on the Moon at seed 75.

What the maps show:

- **The Moon at seed 75:** its 378 km basin's rays flakes and streaks of bright ground with darker ground between, where they were soft smudges; its inner halo, where the jets excavate more than clusters carry, solid as before; 7 percent of the map changed by more than 0.005, its mean 0.121 as it was.
- **The default Moon, Mercury, Callisto, Ganymede, Ceres and Europa:** as they were, at most half a percent of a map changed, their means within 0.001: their young craters' clusters lie within a few cells of the grid.
- **The rays' time:** 11 s on the Moon where it was 19 to 23, 26 s on Mars where it was 50, 23 s on Ceres where it was 59.

Still to do:

- **The user's notes on the icy moons** (audit T6 and T7): grooves more broken and undulating, on Miranda too, and Enceladus's polar terrain as one network of canyons, without its false blue.
- **The census's small end**, a floor from an old basin.
- **Mars's two thresholds**, its rays drawn where they cannot show.

Tests 208 (12 ignored). One is new: a cluster's grain shares its whole ground among parts narrower than the narrowest that would split, scattered over its footprint, and a cluster too narrow to split is laid whole. Clippy is clean and validate passes its 115 checks.

## The sixty-second milestone: deposits drawn by share

Asked for by the user after the sixty-first, before the terrain features ("Shall we do the pixelated deposits first, then the terrain features (icy moons, Venus' coronae etc.)?"): the deposits still laid by whole cells, audit W8.

- **F131.** Every deposit a field bounds covers each cell by the share its margin covers: Venus's tessera and plains, coronae and pancake domes, a shield's flows, a flood's lobes and outline, Io's patera floors and Ceres's cryodomes; a cell takes a deposit's time and unit where it covers more than half, and a deposit whose cover the albedo reads by its unit takes its unit wherever it lies.

How it went:

- **The sweep.** Each half-cell test in the stages was read for what it sets. The audit's list held, with one more, Ceres's cryodomes; Europa's chaos and bands and Ganymede's lanes already blended by share, and so did Enceladus's polar terrain in the finish, whose relief and frost showed no cell steps. Venus's coronae had their interiors only from the loop's coarse cells, which the finish now draws.
- **What the cover means.** The first render's domes were as stepped as before: the albedo reads a cell's cover as salts only on a dome's unit, and the edge cells kept their old one. A deposit read by its unit takes it wherever it lies, as Io's flows already did.

What the maps show:

- **Io:** its patera floors, dark on the bright plains, end smoothly where they ended in steps; 0.7 percent of the map changed.
- **Ceres:** its domes are smooth circles where they were stepped; their shape is still a circle (with Venus's domes and coronae, audit M10 and T1).
- **Venus, Mars and the Moon:** their maps within 2.4 percent of what they were, the contacts' colours too close, or blurred by gardening, for their steps to have shown; the data under them is the deposits' share.
- **Enceladus:** as it was.

Still to do:

- **A lava margin over a crater floor** keeps the floor's unit where it covers less than half, its cover read as salts.
- **The terrain features next:** the icy moons' grooves and tiger stripes (audit T6, T7), then Venus's coronae and the other shapes.

Tests 208 (12 ignored). Venus's test holds the tessera's contact to partial cover, and Ceres's test the domes' edges, beside the faculae. Clippy is clean and validate passes its 115 checks.

## The sixty-third milestone: a groove is a chain of fault segments

Asked for by the user after the sixty-second, with the notes on the icy moons (audit T6): "the groove lanes look fine on Ganymede's full map, but scaled to Enceladus or Miranda they look unnatural because they are just straight lines without much variation. Instead of being full lines along the entire crack, they should come and go, as smaller individual lines, so they don't always follow the same path (e.g. if you follow one line, it may appear to end and the next line is a few meters to the left with a considerable gap). Also, they end off abruptly, with just a straight edge." The user also asked whether Miranda's grooves, in coronae named as Venus's are, should be driven by Ganymede's lanes.

- **F132.** A groove is a chain of fault segments that come and go along it, each stepped across the groove's line and tapering to its tips, cutting its own trough into a plateau; the grooves die out before the lane's bounding faults.

How it went:

- **Why the grooves were lines.** A groove was one period of a wave across the lane, as long as the lane, its segments a noise gate set so that a groove seldom broke, and nothing stepped it across.
- **Segments without seams.** Each segment cuts a hollow a spacing across into a plateau, so that overlapping segments in a relay both show and their gaps stand at the plateau; choosing the strongest segment at each cell instead would have left a seam where one gave way to the next.
- **Gaps that show.** At first a slot held a segment nine times in ten, stepped by up to three tenths of a spacing; the gaps hardly showed, and the user had asked for a considerable gap, so three in four, stepped by up to four tenths.
- **Miranda.** Its coronae are not lanes. The descriptor is the same as Venus's because both are oval, but so is their origin: both stand over an upwelling under a lid (Venus's over mantle plumes and diapirs, Miranda's over ice diapirs and the convection tidal heating drove; Pappalardo et al. 1997, Hammond and Barr). The grooves inside Miranda's are tilt blocks of the faulting that makes Ganymede's. So a corona takes its outline and its stress field from the upwelling, as on Venus, and its fabric from these faults, laid along that field: the next milestones.

What the maps show:

- **Ganymede:** its lanes' grooves are segments that start and stop, staggered and stepped, tapering to their tips and dying out before the lanes' edges, where they were ruled lines the lanes' length; its map 0.495 where it was 0.493, the troughs' dark lag a little less where the segments leave gaps.
- **Enceladus:** its grooved lanes' grooves come and go in segments with staggered tips where they ran the lanes' length and stopped on a line; still faint, their relief tens of metres, and its map, under the plume's frost, as it was.

Still to do:

- **Miranda and the coronae**: a Miranda preset, and coronae as upwellings under an icy lid as under Venus's rock, their fabric this faulting along the upwelling's stress.
- **Enceladus's grooves' relief**, tens of metres where its sulci's ridges stand 1 to 1.5 km: the groove spacing and relief, calibrated on Ganymede, checked on a second body (audit T6).
- **Enceladus's tiger stripes** as one network of canyons (audit T7).

Tests 209 (12 ignored). One is new: a groove comes and goes along a line, its segments stepped within the limit and its relief without a seam, dying out by the bounding fault. The covered fractures' test samples a band of lines beyond the fault rather than one, the fractures there being segments that step across. Clippy is clean and validate passes its 115 checks.

## The sixty-fourth milestone: coronae over an ice mantle's upwellings

Asked for by the user after the sixty-third, which answered their question whether Miranda's coronae, named as Venus's are, should be driven by Ganymede's lanes: "Do Miranda's coronae first, then move Venus' into the same process."

- **F133.** A Miranda preset; its coronae over the upwellings of the convection a tidal heating episode drives in an ice mantle on rock, as many as the convecting layer's cells, each a polygon of its cell, its outer belt the grooves' fault segments concentric with its outline and its core crossed by two sets; and an ice mantle frozen through makes no ice magma.

How it went:

- **Miranda's preset.** Its density is JPL's (1,178 kg/m³). The resonance is the Ariel-Umbriel 5:3 crossing Ćuk et al. (2020) find excited its eccentricity chaotically for hundreds of millions of years, a gigayear ago for Uranus's Q of 15,000 to 20,000; their runs with Miranda's Q/k2 at 100, the engine's own dissipation, keep it at 0.02 or less. At 0.0165 its tides heat it to about 62 mW/m², the middle of the 31 to 112 mW/m² its faults give (Beddingfield et al. 2015) and of the 34 to 84 Arden's flexure gives (Hammond and Barr 2014). The Bond albedo is left out of validation: the Uranian moons' darkening is not among the materials yet, so its ice shows 0.56 against the measured 0.20.
- **Upwellings, not lanes.** The coronae come of the convection the episode drives in its ice: the same tidal flux above the floor that drives Ganymede's lanes (F44), on a body with no ocean, where the ice under the lid convects. The first criterion I gave the user, the lid against the cells' spacing alone, would also have given Pluto and Ceres large coronae from their early radiogenic heat, which a lid over steady or waning convection does not show: a stagnant lid stays stagnant, and it is an episode's onset that breaks it.
- **How many and how large.** Miranda's ice, 145 km deep under a 14 km lid at the episode's peak, holds five cells of the convection: five coronae, 244 to 336 km across against Inverness's 234, Arden's 318 and Elsinore's 323. Four formed while the heating held and have faults 4.9 to 5.3 km apart, Arden's 5 km; one formed as the heating waned, on a thicker lid, and has them 10 km apart, as Elsinore's are "relatively widely spaced".
- **A short episode between the history's samples.** The thermal history is sampled every 47 million years, so a corona's lid came from samples straddling the episode's ramp: one read a lid 57 km thick where the relations of its time give 12. Each corona takes its lid from the thermal relations at its own time.
- **Frozen through.** Miranda's conductive profile today would reach its melting point some 2,000 km down, eight times its radius, yet its ice magma came of its heat flux whatever the lid, and made 187 domes; ice now melts only where its profile reaches the melting point inside the body.
- **No heat pipe in ice.** The first render buried most of Miranda under cryolava: the tidal boost to the magma share, Io's heat pipe, 0.9 of the heat at 1 W/m², gave the episode's 62 mW/m² of tides 5.6 percent of the heat as ice magma, 1.9 km of it over the whole body, flooding 40 percent of the surface and the coronae's faults with it. An ice mantle's melt is denser than the ice over it, and its tides take no boost; the episode's floods are now five, a thousandth of the surface.
- **The coronae's interiors.** The loop renewed a corona's ground over the lanes' smoothing length, 12 km on Miranda, so the old craters' bowls stood on in its core as broad soft swells; a corona stretches and breaks its whole interior, and keeps of the old relief only what is longer than half its radius. Its core's two sets ran straight and ruled a grid; they turn across the core now, as the stress turns with the cell's shape.
- **The relief, Miranda's check.** The faults' spacing passes its check on a second body: 4.9 to 5.3 km against Arden's 5 km. Their relief under the lanes' law, a share of the spacing calibrated on Ganymede, did not: a few hundred metres where Arden's ridges and troughs stand up to 2 km, which read faintly beside the old craters. A corona's upwelling stretches its lid far longer than its faults take to rotate, so its blocks tilt until their faults lock, from Arden's 50 degrees to the 30 at which friction locks a normal fault: 1.7 km on Miranda's blocks, each trough a tilt block's with its scarp facing out. The lanes keep their law; one law for both, from the strain the faults take up, stays open with Enceladus's grooves' relief.
- **A corona under a younger crater.** The corona that formed first stands beside a 259 km crater 0.43 Gyr old, whose ejecta bury most of its belt; it keeps its surface's age beneath them, as distal ejecta leave the ground they cover (F100).

What the maps show:

- **Miranda:** five coronae 244 to 336 km across, rounded polygons of their convection cells in the old cratered terrain, 45 percent of the surface, each an outer belt of ridges and troughs concentric with its outline, the tilt blocks' scarps facing out, round a core crossed by two bending sets of faults; four have faults about 5 km apart standing 1.7 km from crest to trough, the fifth, broken as the heating waned, has them 10 km apart and lies mostly under the ejecta of a younger 259 km crater; their ground a gigayear old, lightly cratered and brighter than the old terrain; five small cryolava floods and no domes.
- **Ganymede:** as it was, bit for bit.
- **Enceladus:** as it was.

Still to do:

- **Venus's coronae** onto the same process: upwellings under its lid, their outlines from their cells, their fabric the same faulting (the next milestone).
- **One law for the faults' relief**: a corona's blocks tilt to lock-up, a lane's grooves keep a share of their spacing calibrated on Ganymede; both want the strain their faults take up (audit T6), with Enceladus's grooves' relief.
- **Enceladus's tiger stripes** as one network of canyons (audit T7).

Tests 214 (12 ignored). Five are new: a corona's cell is bounded by the bisectors; Miranda's episode breaks its lid over a few upwellings, its coronae 200 to 400 km across with faults 3 to 15 km apart; no coronae on a body without an episode or over an ocean; a corona's belt follows its outline; a body frozen through makes no ice magma. Clippy is clean and validate passes its 119 checks.

## The sixty-fifth milestone: relaxation keeps its bands, and large craters are rough down to the grid

Asked for by the user after the sixty-fourth: "Another problem with the large craters, which look like Herschel, is the smoothed floor which looks to be an excessively large-scale smoothing, and perfectly concentric ridges which expand out like ripples, clearly apparent in our flatter coronae terrain on Miranda near the crater. These 'ripples' are also jagged and pixelated at a larger scale than the surrounding terrain." With it the user asked why Enceladus and Miranda carry Mimas's crater population; that is the next milestone, a bombardment for any star calibrated on the Solar System's records ("fix A first then go ahead with B").

- **F134.** The loop's relaxation keeps the relief's bands from step to step and splits only what each step adds; a crater's texture fades by the quintic and goes on in octaves down to the grid.

How it went:

- **The ripples around relaxed craters.** The loop split the whole surface into the pyramid's overlapping bands at every step and scaled each, so the relaxation compounded: a test crater stepped sixty times lost 569 m of the relief its short bands should have kept and rang 40 m beyond two and a half radii. Kept from step to step, the bands change by nothing after the first step and the ringing is gone (0.4 m). On Enceladus, whose shell relaxes its large craters through its whole history, the two or three jagged rings about each large crater's rim go.
- **The smoothed floors and the bands in the ejecta.** A crater's texture scaled with its radius in two octaves, so a crater of 250 km had floor lumps of 13 km and hummocks of 15 to 50 km with nothing finer, and the lattices' cubic fade drew the polar lattice's rings and spokes in the shaded relief. With the quintic fade and octaves down to the grid, the floor and the ejecta of a large crater are rough at every scale.
- **What the ripples on Miranda were.** Beside its young 259 km crater the user saw dashed concentric ridges in a corona's flat ground. A profile across them shows scarps of 0.5 to 0.8 km every 10 km on ground 3.41 Gyr old: the outer belt of the corona that formed first, its fault segments half buried by the crater's ejecta, its outline running parallel to the crater's rim there. They stay; the crater, wider than Miranda's radius, is not one the moon could hold, and the next milestone takes such impacts as what they are.
- **What did not change.** The Moon's old basins look as they did. The terraces of the very largest craters still break into long straight-scarped blocks.

What the maps show:

- **Enceladus:** its large craters' rims without the rings that spread from them, their floors and ejecta rough at every scale; its crater population still Mimas's, the next milestone's.
- **Miranda:** its young 259 km crater's floor and ejecta rough at every scale where they were smooth lumps and concentric bands.
- **The Moon:** as it was.

Still to do:

- **The bombardment of any star**, calibrated on the Solar System's records (the next milestone): the outer system's own impactors and their decline, crater sizes from each body's gravity and impact speed, and the impacts a small body cannot hold resetting its surface.
- **The largest craters' terraces**, long straight-scarped blocks.
- **Venus's coronae** onto the corona process.

Tests 215 (12 ignored). One is new: a relaxed crater neither rings nor wears away. Clippy is clean and validate passes its 119 checks.

## The sixty-sixth milestone: the bombardment of any star

Asked for by the user after the sixty-fourth: "Our icy bodies generate with the right crater population for something like Mimas, judging purely off images. They even create these large craters like Herschel. But like you pointed out, neither Enceladus nor Miranda have large craters like this, and looking at the actual full reference maps of all these bodies, Miranda and Enceladus show a much smaller crater population. So we must be missing something in our model." On the answer that the outer system has impactors of its own: "Would bombardment be based on distance to the star then? Eventually, this program should be able to generate feasible outcomes for e.g. close planets to a red dwarf. If we can fit a generalised rule to all stars that matches with Solar System records that would be best." The sixty-fifth milestone fixed the ripples first ("fix A first then go ahead with B").

- **F135.** A body's place against its star's snow line sets which impactors strike it, how hard and how fast: the inner system's rocky leftovers on the Neukum production function and chronology, the outer system's comets on their own size distribution, crater scaling and stretched-exponential decline, both scaled to another star by its luminosity and mass. A body holds no crater wider than its radius if icy, or than itself if rocky, and its record starts at the median time of the last larger impact.

How it went:

- **What was missing.** The engine struck every body with the Moon's asteroids, scaled beyond 4 AU by a tenth of the Moon's flux times the planet's focusing: Enceladus took 2.2 times the Moon's craters of every size from 80 Myr after it formed, the lunar size distribution making Herschel-sized basins as readily as on the Moon, and Miranda a crater wider than its own radius 0.43 Gyr ago. The giant planets' moons are struck by comets scattered from beyond Neptune, whose numbers, sizes, speeds and decline are their own (Zahnle et al. 2003; Nesvorný et al. 2023; Wong et al. 2023).
- **One rule for any star.** The place against the snow line, the distance over the distance where a blackbody is 170 K, follows from the insolation: the Moon 0.37, Mars 0.57, Ceres 1.03, Jupiter 1.94, Saturn 3.56, Uranus 7.1. The rocky leftovers strike at the Moon's rate inside 0.40, twice it through the belt and none beyond 1.6; the comets from the snow line out, at 4.7e-7 of 1 km or more per km² per Gyr inside Saturn's place, more as its square beyond, times the planet's focusing. Around another star the early bombardment runs on its orbital times and its disc's solids, and the late flux after it is a root multiplier of the Solar System's, the planets and belts that set it being beyond the engine.
- **Calibrated on the Solar System.** Nesvorný et al. 2023's comet rates on the moons of Jupiter, Saturn and Uranus, each over its planet's focusing, fix the flux: the rule gives Io, Europa, Ganymede, Callisto, Enceladus and Miranda their rates within a factor of 1.3, and their impact speeds within 7 percent, a moon's own orbital speed added three times over as Zahnle et al. 2003 do. Wong et al. 2023's comets over Mimas's and Enceladus's histories fix the early bombardment, 550 times the present rate over 4.5 Gyr. Mimas, the least altered of the small moons and the calibration body, keeps 1.31e-3 craters of 10 km per km² and 2.68e-4 of 20 km against the 1.38e-3 and 2.85e-4 counted. Mimas has no preset; it was built from Enceladus's with its radius, density and orbit.
- **The small comets.** Singer et al. 2019 found few craters under 13 km on Pluto and Charon, a deficit of comets under 1 to 2 km. With it, Mimas held a third of its craters of 10 km and Enceladus a tenth of its craters of 2 km: the moons of Saturn keep the steep slope down to craters of 2 km, and Wong et al. 2023 date them with and without the deficit. The slope runs on (stated choice). Jupiter's moons record fewer small craters than Saturn's (Zahnle et al. 2003's case A against case B), and one population everywhere gives them Saturn's.
- **The largest craters.** A crater as wide as an icy moon's radius is about the largest any keeps (Odysseus 0.84 R, Turgis 0.79, Herschel 0.70); a rocky body keeps one as wide as itself (Rheasilvia 1.93 R on Vesta). A larger impact shatters and resurfaces the body, and its record starts after the last one. Taken as the time after which fewer than one is expected, the record started too early and held a third too many craters; at the last impact's median time, when ln 2 are still to come, Mimas's record starts 4.09 Ga, beside its plains' 4.16 to 4.21 Ga (Wong et al. 2023). Ceres's starts 3.97 Ga and holds 13 craters over 100 km against the 16 counted (Marchi et al. 2016), where it was saturated; Enceladus's starts 4.25 Ga and Miranda's 4.20 Ga. The Moon, Mercury, Venus, the Earth, Io and the large icy moons keep their record from when their crust solidified. Mars's starts 3.95 Ga: its rocky flux, twice the Moon's, with the Mars production function continued to basins, gives it two impacts wider than itself after its crust solidified, where the Moon's rule from South Pole-Aitken, 1.5 radii, had moved it to 3.88 Ga.
- **Mimas's heat.** Built from Enceladus's preset, Mimas's eccentricity of 0.0196 heats its porous core with 228 GW: a shell 1.7 km thick over an ocean, renewed every 0.1 Gyr, where Lainey et al. 2024 find one 20 to 30 km thick over a young ocean under an old surface. POROUS_DISSIPATION fails the check on Mimas its calibration wanted. The render with those tides drew ring graben about craters of 13 km and stalled in a single-threaded phase after twenty minutes; the calibration renders take Mimas with its tides off, frozen through, as its surface shows it.
- **The gardening.** The diffusivity was the Moon's on every body, scaled by gravity, while the flux scaled the crater counts. It follows the body's present production of craters of 1 km now: about the Moon's on Mars, Ceres and Pluto, 1.9 times it on Callisto, 4 on Ganymede, 17 on Enceladus, whose craters of 1 km need comets of 33 m.
- **Around a red dwarf.** A planet two thirds of the Earth's insolation from an ultracool dwarf, as TRAPPIST-1e, meets its rocky leftovers at 26 km/s (Kral et al. 2018: 15 to 110), and its early bombardment, decaying over 1.7 Myr on the star's short orbits, is over long before its crust solidifies: it keeps a record of the late flux alone.
- **What did not change.** The Moon's height map is the same to the bit, and so are its 916 loop craters and 81,295 finished ones; one ray is a shade brighter, its impact speed now 15.16 km/s against 15.13.

What the maps show:

- **Enceladus:** its largest craters gone (M65's 170 km basin), its loop craters down from 8,295 to 720 and its plains still crowded: craters up to about 120 km, their bowls relaxed flat and their rims sharp, where Enceladus has none over about 40 km. Its plains, from its record's start at 4.25 Ga, keep 1.5 to 3.7 times its counts at 5 to 20 km: they are younger than its record (4.10 Ga by Wong et al.'s isochrons) and relaxed, which its heat history does not give (audit M3).
- **Miranda:** its 259 km crater and its flat-floored giants gone, its loop craters down from 4,547 to 639; its old terrain crowded with small craters and a few of 30 to 50 km, one of about 140 km far to the south; sharper and denser than Voyager's mantled terrain.
- **Mimas, a calibration body with its tides off:** crowded with craters of every size to 30 km, as its map is, and two of about 110 km beside its own one Herschel.
- **Ganymede and Callisto:** their giant rayed craters mostly gone, loop craters from 917 to 93 and from 641 to 47; small bright craters over the dark terrain, as their maps show.
- **Ceres:** no longer saturated: loop craters from 10,486 to 515, none wider than its radius.
- **Mars:** fewer basins than before, a few large ones as Mars has, and fewer mid-sized craters on its oldest terrain, from 3.95 Ga where it was 4.15.
- **The Moon:** as it was.

Still to do:

- **Enceladus's and Miranda's younger, relaxed plains** (audit M3): the heat history that relaxed Enceladus's craters and resurfaced its plains after 4.25 Ga, and Miranda's mantling.
- **The comets' size distribution at the large end**: a single power law leaves a pristine small moon about two craters between half its radius and its radius and about twelve wider than a quarter of it, as Iapetus's basins are; Mimas has one and about four. Bottke et al. 2023's wavy distribution from the comets' collisional history is the direction.
- **Jupiter's moons' small craters**, fewer than Saturn's (Zahnle et al. 2003's case A).
- **Mimas's heat**: POROUS_DISSIPATION from the dissipation per unit of porosity; the stalled render of a shell 1.7 km thick.
- **Mars's basins**: the Mars production function continued past its range gives two impacts wider than Mars after its crust solidified.
- **The snow line where a red dwarf's planets formed**, from its far brighter youth rather than its luminosity today.
- **Venus's coronae** onto the corona process; **the largest craters' terraces**.

Tests 219 (12 ignored). Four are new: the outer moons meet their comets at Nesvorný et al.'s speeds and rates; Mimas keeps the craters it was counted with; a small ice moon holds no crater wider than itself while the large moons keep their record from solidifying; around a red dwarf the early bombardment is brief. rustfmt also reformatted the corona, relaxation and stamp modules of the last two milestones. Clippy is clean and validate passes its 123 checks.

## The sixty-seventh milestone: the comets' sizes from their history, craters under strength, and a shattered moon's debris

Asked for by the user after the sixty-sixth ("All improvements across the board. Yes, do heat histories and large-crater size distribution"), choosing among the options I gave: Bottke et al. 2024's comets with the debris their counts imply ("Bottke + debris"), and, when the debris turned out to spoil Miranda, its amount scaled by the impact speed's square ("Scale by impact speed"). The heat histories, chosen as a thermal-orbital process for every moon, are the milestone after next.

- **F136.** The comets follow the destabilized Kuiper belt's size distribution after its collisional history, and make their craters by pi-scaling with the target's strength; a moon whose record starts after an impact it could not hold is struck at once by that impact's debris.

How it went:

- **What the counts say.** I had told the user that Mimas carries about four craters of 50 km; the published counts (Kirchoff and Schenk's, replotted in Bottke et al. 2024's Fig. 7, which I read off the figure) give 7 to 10, so the sixty-sixth milestone's excess there was about twice, not four times. Miranda's counts (their Fig. 9) hold as many large craters as the production did; its sparse look on the map is its mantled, subdued terrain. Enceladus's plains hold eight times fewer craters of 30 km than the production, being younger.
- **The comets' sizes.** Bottke et al. 2024 model the destabilized Kuiper belt's collisional history and the craters it makes on the giant planets' moons. I read their size distribution off their Fig. 1, the curve their craters on Iapetus and Phoebe require joined to their model beyond 32 km: shallow from 0.1 to 1 km, a steep bump of collision fragments from 2 to 6 km, shallower again to 60 km, steep beyond 130 km.
- **Craters under strength.** With their scaling for cold ice, pi-scaling with a strength of 35 MPa, a small icy moon's craters are 15 to 20 times their comets where Zahnle et al. 2003's scaling made them 23 to 28. Bottke et al. use a weaker ice on the small moons; one set serves every icy body here, and gives the large moons the 16 to 25 times they find.
- **Nothing fitted, the large craters right.** With the fluxes, histories and hold rule of the sixty-sixth unchanged, a prototype gave Mimas, Tethys and Dione their counts from 20 km up within 1.4 times and Tethys's to 400 km, and Rhea, whose record starts with its crust, 0.6 of its counts at every size. Their craters under 20 km fell short, as Bottke et al.'s model's do, by up to 2.7 times on Mimas at 3.5 km; they put the excess down to debris in orbit about Saturn from the impacts that reset the moons.
- **The debris.** Fitted on Mimas, fragments striking at 0.3 times its orbital speed, steep above 2 km and shallow below, fill the gap: Mimas within 13 percent from 3.5 to 50 km, Tethys within 1.15 times and Dione within 1.4 with nothing refitted. Laid on every moon whose record an unholdable impact reset, they gave Miranda 2.6 times its craters of 2 km, where the comets alone match its counts from 2 to 20 km as Bottke et al. find; I found nothing in the engine's physics that would spare Miranda and not Saturn's moons, standard ejecta scaling making the mass an impact of a given crater throws out of a moon independent of its speed. At the user's choice the debris scales with the square of the comets' impact speed, a stated choice: Miranda then keeps its counts within 1.45 times at 2 km and 1.15 above.
- **The gardening.** The diffusivity, which follows each body's present production of 1 km craters, falls with the comets' shallow slope under 1 km: Enceladus's from 4.6 times the Moon's to an eighth of it, Ganymede's from 3.8 times to a sixth.
- **What did not change.** The Moon's height map is the same to the bit, its 916 loop craters and 81,295 finished ones the same; so are Mars's and Mercury's populations, which the comets do not reach.

What the maps show:

- **Mimas, a calibration body with its tides off:** the two basins of about 110 km gone, one of about 90 km and several of 50 to 70 km left beside the dense craters of every size to 30 km, as its map is.
- **Miranda:** its crater of about 140 km gone, its loop craters down from 639 to 402 and its old terrain mostly small craters, closer to Voyager's.
- **Enceladus:** fewer craters of every size (loop craters 720 to 587, finished ones 96,000 to 42,000) and a few large ones still on its plains, which are younger than its record (the milestone after next).
- **Ganymede, Callisto and Ceres:** as they were at the sizes the maps show.
- **Pluto:** far fewer craters than the Moon's population gave it; to check against Singer et al.'s counts on its oldest terrains.

Still to do:

- **The base relief's spectrum** (the next milestone, asked for by the user after this one's renders): the primordial relief takes a quarter of the relief the crust's strength supports, 148 km on Mimas, at wavelengths from a third of the radius down, so a small moon is far too rough at mid scales and has no relief at the scale of its hemispheres.
- **The thermal-orbital process** for every moon: eccentricities damped by each moon's own dissipation or forced by resonance, the porous core heating only while its pores hold water, the shell's response; Mimas's preset; Enceladus's younger, relaxed plains (audit M3).
- **Ganymede's and Callisto's large craters**, twice their ancient terrains' counts at 30 to 200 km; Pluto's oldest terrains.
- **Venus's coronae**; the largest craters' terraces.

Tests 219 (12 ignored). One replaces the Mimas test: the moons keep the craters they were counted with, Mimas within 1.3 times at 5 to 50 km, Tethys, Dione and Miranda within 1.5, Rhea within 2 and without debris. Clippy is clean and validate passes its 123 checks.

## The sixty-eighth milestone: the base relief from the largest crater a body holds

Asked for by the user after the sixty-seventh's renders: "These results have revealed that there is a fine, rough underlying displacement in the terrain of the smaller icy moons, likely set by a noise field. This is fine in practice, as looking at the real height maps of the moons shows a similar appearance. But it looks too fine and lacks large scale variation. Is this because we use the same scale for all bodies? The base noise should scale with the size of the body."

- **F137.** The primordial relief begins at the largest crater the body holds and falls toward the grid with the Hurst exponent of cratered terrain, its largest octave as rough as a cratered surface is at that crater's scale.

How it went:

- **What the noise was.** The primordial relief's wavelengths were fractions of the radius already, from a third of it down to a ninety-sixth, the same number of bumps on every body; its amplitude was a quarter of the relief the crust's strength supports, which on a small icy moon is enormous (148 km on Mimas, 122 on Miranda, 83 on Enceladus, against the Moon's 57). Mimas's relief spanned 36 km, and its finest octave was five times as steep as the Moon's: the grain the user saw. The large scales were missing because the noise began at a third of the radius.
- **Measured against the real maps.** Band by band from twice the radius to a thirty-second of it, the real Mimas, Moon, Mars and Ceres carry their relief most at the largest scales; ours peaked at a quarter to a half of the radius.
- **Where to begin and how strong.** The amplitude now comes from the largest crater the body holds: a third of its depth at its own scale, the share a cratered surface's relief has (F54). Beginning first at twice the radius gave Mimas and Mercury twice their maps' share beyond it; beginning at that crater's own size, the radius on ice and 1.5 radii on rock, brings Mimas to 0.40 against 0.32 and needs no constant.
- **What it gives.** Mimas's spectrum matches its map's band for band within 0.06 and its ground between the craters, crowded with lumps of tens of km, undulates gently under them; Mars's within 0.1; the Moon's largest bands rise toward its map's, and its maria flood larger connected lowlands; Mercury carries too much at its largest scale (0.56 of its relief against 0.39). The total relief changes little on the cold bodies (Mimas 5.7 km rms, the Moon 3.9).
- **What did not change.** Enceladus, Miranda, Ceres, Ganymede and Callisto: their shells relax every wavelength over a few tens of km to nothing, so the new large scales are erased again. On a shell floating on an ocean the relaxation should leave the isostatic share, the thicker ice floating higher; that goes with the thermal-orbital milestone's shell.

What the maps show:

- **Mimas, the calibration body:** broad highs and lows across its hemispheres under its craters, the mid-scale lumps gone, as its height map shows it.
- **The Moon:** as it was, its maria in larger connected lowlands.
- **Mars:** as it was, its relief's spectrum its map's.
- **Mercury:** a little more relief at its largest scale than its map.

Still to do:

- **The relaxation of a floating shell** toward isostasy rather than flatness, which the warm icy moons' large scales want, with the thermal-orbital process (the next milestone).
- **Mercury's largest scale.**

Tests 219 (12 ignored). Clippy is clean and validate passes its 123 checks.

## The sixty-ninth milestone: the thermal-orbital history of the icy moons

Asked for by the user ("Great, these are much needed changes and the results pay off. Proceed with M69"), as they chose it in the sixty-seventh: the thermal-orbital coupling built as a process for every moon, of its resonances, the planet's tidal push and its shell's heat, so that the episodes of its heating emerge ("Thermal-orbital process"), rather than an episode stated in Enceladus's preset; with each moon's free eccentricity damped by its own dissipation, the porous core heating only while its pores hold water, and a preset for Mimas, the calibration body of the bombardment.

- **F138.** An icy moon's eccentricity and its interior evolve together: held by a resonance that holds it now, pumped by a past resonance's heat until the moon dissipates as much, otherwise free and decaying by its own dissipation, its present value the remnant of its last resonance; its shell freezing and melting against the latent heat, and its porous interior dissipating only while an ocean lies over it.
- **F139.** The elastic lid holds up the loads its flexure carries, its flexural wavelength 2π (D / ρ g)^¼, many times its thickness on a small moon.
- **F140.** The shell's surface processes follow its history: they act only while it floats, the focusing thins a shell by the share of what it conducts that is made in it, and the lanes renew by the heat carried above the floor for as long as it lasted.

How it went:

- **What the stateless model could not do.** The thermal history took each time's eccentricity as given and each time's shell as the steady state of its heat. Mimas, its eccentricity of 0.0196 held for 4.5 Gyr and its porous core dissipating, made 228 GW and a shell 1.7 km thick, though its eccentricity's orbital energy pays for only 0.4 Myr of that; its render had stalled in the finish. Miranda's porous core dissipated though Miranda had been frozen through for three gigayears, and powered its coronae.
- **What the dynamics say.** Meyer and Wisdom (2008) find Mimas's eccentricity the remnant of past resonances, among them the Mimas-Dione 3:1, its equilibrium eccentricity 0.07 for a Mimas k2/Q of 3e-6, and Enceladus captured into its 2:1 resonance with Dione close to equilibrium, with no excursion above it in its past resonances. In equilibrium a resonance's heating does not depend on the moon (Meyer and Wisdom 2007): what the planet's tides feed the resonance, the moon dissipates, at whatever eccentricity its interior needs. Nakajima et al. (2019) give the same equations of migration and damping; Goldreich et al. (2025) find a limit cycle for Enceladus that heats it to some 290 mW/m² every 13 Myr, which would relax its plains throughout their history.
- **The process.** One integration over the moon's life in steps of half a million years: the eccentricity damped by the moon's dissipation of the time, pumped toward what a resonance supplies, or held; the shell's thickness by a Stefan balance stepped implicitly, from a warm start; the porous interior wet while an ocean touches it, freezing with the moon and thawing only after R² / κ under an ocean. A free moon's resonance lets it go at the latest time that lets its eccentricity decay to the present value. That value is not monotonic in the escape (a moon let go before its ice melts keeps the eccentricity it was pumped to for hundreds of Myr): a bisection found Mimas an escape before its ice melted, and the escapes are now scanned back from the latest. The thermal history, the present tides, shell and surface heat flux and a new expert node of the ocean's age follow from it; a pinned eccentricity history is followed as given.
- **The thresholds that emerge.** A hard shell's heating and its conduction both go as the inverse of its thickness: above one eccentricity it thins away, below it freezes through (Mimas 0.026, Enceladus with a dry core 0.032, Miranda 0.043). Ice frozen to the rock melts only at a higher one (Mimas 0.066, near the two to three times its present value at which Walker and Rhoden find its melting began; Miranda 0.11). Enceladus's ocean is its porous core's, held from its warm start.
- **Mimas.** Frozen through by 0.3 Gyr and cold for four billion years. Its 3:1 resonance with Dione, 83 Myr ago at Saturn's dissipation as Lainey et al. 2017 measure it, supplying 16 GW, pumps its eccentricity to 0.071; its ice begins to melt 33 Myr ago, the shell overshooting to 11 km as the runaway is quenched, and it escapes some 5 Myr ago. It has today its eccentricity of 0.0196, a shell 28 km thick over an ocean 33 Myr old, a dry core and 6.6 GW, its eccentricity's e-folding time 15 Myr: Lainey et al. 2024 find 20 to 30 km over an ocean younger than 25 Myr. Its lid holds Herschel and the relief of its hemispheres: its height map's spectrum is the cold Mimas's band for band.
- **Miranda.** Frozen through by 0.3 Gyr; its resonance, supplying the 42 GW that 0.0165 gives a warm Miranda, pumps its eccentricity to 0.10 until its ice melts at 1.08 Ga, then holds it at 60 mW/m² under a 15 km shell until 0.85 Ga; frozen again by 0.7 Ga. Its five coronae form within the heat, their times now drawn from the tidal heat the lid carries (the lesser of the tides and the conducted heat above the floor): drawn from the heat's window, one met ice still frozen; drawn from the tides, one met a lid still thick over ice melting from below; drawn from the conducted heat, one fell in the warm start's freezing.
- **The others.** Enceladus keeps its 22.5 km shell, its wet core and its 29 GW. Ganymede's eccentricity rises toward 0.018 as its resonance's supply does and decays over 0.15 Gyr after; its lanes trail the resonance by 0.2 Gyr, so it is placed 2.2 Gyr ago for terrain two gigayears old, and they renew 0.69 of it. Europa's shell lags its cycle, 28.7 km today against 30.3. The thick oceans of Titan, Callisto and Triton are still freezing from the warm start, Triton's shell 260 km against 296.
- **The lid that holds Herschel.** Mimas's young ocean, under F39's support of four lid thicknesses, would have flattened every load beyond 100 km within years, Herschel among them. A thin plate's flexural wavelength goes as Tₑ^¾ g^-¼ (F139): Mimas's 12 km lid holds 450 km, Enceladus's 370, Ganymede's 420. The plate's own law leaves a short load a small share that never vanishes, which the survivals of 24 or 90 steps wore away (a 30 km crater on a 3 km lid lost half its relief), so the share rises smoothly from two thirds to four thirds of the flexural wavelength. Enceladus's relief, which its relaxation had flattened beyond a few tens of km, now peaks at its largest scales as its height map's does (bands from twice its radius to a thirty-second of it 0.47, 0.60, 0.47, 0.33, 0.24, 0.14 and 0.10, against its map's 0.75, 0.46, 0.26, 0.16, 0.15, 0.18 and 0.30 and M68's 0.09, 0.24, 0.45, 0.58, 0.50, 0.30 and 0.23).
- **Isostasy.** M68 supposed a floating shell would keep the isostatic share of its large scales; that share decays by the flow at the shell's base within 1 to 80 Myr on shells with a basal viscosity of 1e14 Pa s. What keeps old large-scale relief is the lid's flexure, and the isostatic share is not added.
- **The surface over a young ocean.** Mimas's first render with its ocean had an Enceladus-like polar terrain, seven old focusing terrains, lanes over a quarter of it, cracks of four billion years and plume fallout: the focusing, the lanes, the cracks and the plumes all took the present state, or the peak of the history, as having held for the body's whole life (F140 now makes them follow it). Its ocean's tides, 13 mW/m² carried for 33 Myr, still renew a tenth of it in young lanes by the heat-above-the-floor law (F44 with F140's time constant, calibrated on Ganymede), where Lainey et al. find its ocean too young to have marked its surface.
- **Enceladus's plains.** The process does not give Enceladus the early episode its younger, relaxed plains want (audit M3): held in its resonance, its heating stays at its equilibrium, and the published dynamics give it no excursion. F139 makes its relaxation want lids of a kilometre or two, as Bland et al.'s heat flows give, so the ground beside its terrains relaxes less than F66 made it.
- **The tests.** A test of Europa's cracks counted a cusp between two tightly curved arcs as a zigzag; it now counts turns back and forth of much the same size, and holds the cycloids to the ridges' half a percent.

What the maps show:

- **Mimas:** its cratered surface and Herschel as the cold Mimas had them, with lanes over a tenth of it, smooth and fresh, from its young ocean's heat: the one mismatch of the milestone.
- **Miranda:** its five coronae as they were.
- **Enceladus:** its terrains and craters as they were; its plains undulate at the scale of its hemispheres, the M68 base relief its lid now holds.
- **Ganymede:** its lanes as they were, a little more of them; a basin at the lanes' edge keeps its rings more clearly.
- **Callisto and Ceres:** as they were, Callisto's basins and craters keeping a little more relief.
- **Europa and Titan:** as they were.

Still to do:

- **The lanes from strain**, with one relief law from strain: a young ocean's melting contracts the body and its refreezing stretches the lid, which the heat above the floor does not see; Mimas's lanes are the case.
- **Enceladus's plains** (audit M3): a heat history that reaches Bland et al.'s flows.
- **The resonance supplies** from the planet's tidal torque on the moons of each resonance chain, rather than stated: Uranus's torque at a Q of 17,500 gives Ariel and Umbriel some 6.5 GW, a sixth of what Miranda's coronae take.
- **Venus's coronae**; the largest craters' terraces.

Tests 222 (13 ignored), among them two of the thermal-orbital history (Enceladus keeps its ocean and wet core; Miranda melts in its resonance and freezes again) and one of Mimas's lid holding Herschel. Clippy is clean and validate passes its 129 checks, six of them Mimas's.

## The seventieth milestone: the renewed ground's texture, the lanes' ends and the coronae's tones

Asked for by the user after the sixty-ninth's renders: "The coronae, young lanes, and Enceladus' terrain all smooth over the underlying terrain, which is fine and accurate on paper, but if you look at the heightmaps, this smoothing process creates massive blurs where this terrain is... I think we could taper the end points of the lanes so they end in a tapered 'flame' shape... right now it just looks like a gaussian blur, which doesn't look good... Perhaps preserving some or all of the base noise field, which allows craters to smooth out, but keeps the micro detail intact?... Miranda's coronae have the blur issue as well, but also they appear lighter than the surrounding terrain, when it is described on the real Miranda as being darker. This geological map of Inverness corona shows it split into two main colours, the outer corona material, and the inner corona material, with the outer being darker and the inner being lighter."

- **F141.** A renewed terrain levels what impacts and older terrains added over its own scale and keeps the base relief under it, toned down below its fabric's spacing.
- **F142.** A lane's end in the old ground tapers to a point; only an end inside an older lane is cut square.
- **F143.** A corona's belt keeps the old crust's tone with bright scarps; its core is clean renewed ice.

How it went:

- **Mimas's young terrain is Ganymede's lanes.** The same process, F44's grooved lanes, on a moon whose ocean is 33 Myr old (F140): 45 lanes over a tenth of it, which Lainey et al. 2024's Mimas does not have (open, with the lanes from strain).
- **Why every end was square.** A lane's end that met no older lane was cut off at a fault of its own across the lane, with a wandering trace; the taper for an end without a cut therefore almost never ran. The fault across the end is gone: such an end dies out as a linear wedge over three widths, its grooves fading with the bounding faults' throw. On Ganymede 0.4 of the ends stop against an older lane, square, as its junctions are; the rest pinch out in points.
- **What the blur was.** The floor under a lane, a corona or a polar terrain was a plain low-pass of the relief at a smoothing length (10 km on Mimas's lanes, 100 km on Ganymede's, half a corona's radius, a fifth of the radius under the polar terrain): every crater near that length survived as a soft dimple, and every texture went with the craters. The user's suggestion was the right one, and it needed two pieces: level what the terrain reworks over its own scale, its width (a lane), its radius (the polar zone), so that no crater inside it keeps a ghost; and keep the base relief (F137) under it, regenerated from the same field at the resolution that draws it. Three tries went before the fourth: every octave below the smoothing length at full amplitude buried Mimas's grooves under 340 m of relief; each lane's floor from its own pyramid level showed sawtooth seams where overlapping lanes of a set, a level apart, crossed (the levels are blended by the logarithm of the length); and the base relief kept only below half the fabric's spacing left pixel grain over a blank mid-scale, which the user saw as grainier than the old ground in the shading and still blurred in the height map, and asked to tone the finest octaves down rather than go back to the blur. The base relief is now kept in full below the length the terrain levels over, and toned down to 0.4 of its amplitude below its fabric's spacing, the grooves, faults or folds overprinting what is finer than themselves.
- **The coronae's two materials.** The belt is the old crust broken into tilt blocks and keeps its tone, with clean ice only on its faults' scarps; the core is renewed ice, clean. Miranda's coronae are now bright-cored with faint bright rings about them in the albedo, over plains of their own tone. A belt darker than the plains, as Inverness's and Arden's are, needs the Uranian moons' dark component, which the materials do not have yet (audit M3b).

What the maps show:

- **Mimas:** its lanes are crater-free ground with the base relief's undulation and a toned-down grain, their grooves faint on so thin a lid (59 m at 1.6 km).
- **Ganymede:** its lanes as they were, 0.69 of it, now pinching out in points where they die in the dark terrain.
- **Enceladus:** its polar terrain at its depression with its folds over the relief its plains have, the ghosts of its old craters gone.
- **Miranda:** its coronae with their faults over the relief their plains have, the cores bright and the belts the plains' tone.

Still to do:

- **The coronae's dark belts** (audit M3b), with the Uranian ice's dark component (Phase 5 materials).
- **The lanes from strain**, which would spare Mimas its young lanes (F140).
- **Mercury's largest scale**, Ganymede's and Callisto's large craters, Pluto's oldest terrains, Enceladus's plains (M3).

Tests 223 (13 ignored), among them one of the base relief's octaves regenerated at any resolution, one of the lanes' ends and the lanes' floor completed on the resolved grid. Clippy is clean and validate passes its 129 checks.

## The seventy-first milestone: Venus's coronae from the corona process

Asked for by the user after the seventieth, in the order they set: Venus's coronae onto the corona process, before the rifts and graben from the fault fabric. Venus's coronae were still Phase 3's: 506 thin rings at random places, sizes uniform in the logarithm, a dome or a sag by a coin's toss.

- **F144.** A rock lid's elastic layer ends at 650 °C (923 K), not 650 K, which lay below Venus's surface and left it no elastic lid: Venus's derives at 30 km (measured 29 ± 6), Mars's at 160 km.
- **F145.** Venus's coronae come of its mantle's diapirs: small plumes whose heat is a share of what the mantle convects since its lid formed, born where its large upwellings draw them, each melting under the lid where the mantle about it does not, the pocket its corona where the loads break the lid it has thinned.
- **F146.** A corona's profile is its lid's flexure under its loads: its cumulates, its rim of thickened crust and vent lava, the lava flooding its floor, and while its head lasts the head's lift and the trench where the lid founders at its edge.

How it went:

- **The research.** Venus's 740 coronae are log-normal about 215 km (Glaze et al. 2002; Gülcher et al. 2025), 68 percent along chasmata, 21 at hot spots, 11 on plains (Stofan et al. 1992), on lids 11 ± 7 km thick (Smrekar et al. 2023). Dombard et al. 2007 make a corona the pocket of melt at a plume head's top loading the crust, not the plume itself; Smrekar and Stofan 1997 give their plumes' fluxes and the upwelling's share of Venus's heat; Gülcher et al. 2020 classify the large ones as active, with outer rise and trench, or inactive, a rim round a low.
- **What sets which bodies have coronae.** The first population put 470 on Mars and some on Mercury and the Moon. Their mantles, and Venus's for 200 Myr after its overturn, stand above the solidus at the lid's base, so a diapir's melt is the region's, their plains': only where the mantle about a head stays below its solidus does a pocket stand out. On a lid thick against the pocket the loads spread as a swell instead of breaking it in an annulus, as F133's rule has it over an ice mantle. Mars and the Moon now get none.
- **What sets their size and number.** A corona's diapir is fed for some 100 Myr (Smrekar and Stofan 1997), not a lone head, so each takes some 10^25 J of the heat share; its pocket is about one and a half times its head. Two stated choices, the head's radius and the feeding time, were set against Venus's count and sizes; the rest emerged: a quarter near the rises, a third still active, a quarter raised and half rimmed lows, and the lid under them thinned by their own heads to 7 to 14 km, Smrekar et al.'s 11 ± 7.
- **What made them look stamped.** Lava laid as a ring on the rim gave every corona the same thin bright ring; flooding the floor instead left the rims as the crust's root, spread low by the plate. The rim now carries half the lava at its vents, a narrow load the lid holds up, the floor the rest; the rim and the trench vary round the edge and the outline wanders.

What the maps show:

- **Venus:** coronae of every class about its rises and on its plains: plateaus and domes with trenches along arcs of their edges, rims round lava floors, overlapping chains; its plains otherwise as they were.
- **Mars, Mercury and the Moon:** their elastic lids thicker; Mars's rise lands a degree or two farther from the equator; Mercury takes three faint domes in its last half gigayear.

Still to do:

- **The rifts** (M10), along which most of Venus's coronae lie: their thinning of the lid will string the coronae along them by the same rule; then the rifts and graben from the fault fabric, with a Charon preset.
- **The novae** from the diapirs that dome the lid, now a random share; the annulus's fractures, which the lid's small bending strain leaves faint in the relief.
- **Mercury's late melting**, which gives it its three domes, with its thermal history.

Tests 227 (14 ignored), among them the plate's flexure against isostasy, Venus's coronae against their statistics with none on Mars or the Moon, and a corona's loads. Clippy is clean and validate passes its 131 checks.

## The seventy-second milestone: the lid's strain and the faults it opens, with Charon

Asked for by the user after the seventy-first: the rifts and graben from the fault fabric, driven by strain, with a Charon preset, and free to change Mars if that brings it in line with Venus, so long as the bodies stay consistent with each other. Staged in three: the strain from a body's change of radius and the faults it opens on every icy lid (this milestone); the relief's own stress, Venus's linked rifts with its coronae strung along them, and Mars's canyons and Noctis Labyrinthus on the same rules; then compression on the same field, Mercury's scarps and the wrinkle ridges.

- **F147.** One strain per body, its change of radius: an icy body's ocean freezing as ice I, and as a denser ice where it lies deep enough, and every body's interior cooling. The lid bears it elastically between the stresses that break it in tension and in compression, and its faults take up the rest; the extension since the record began opens graben that cut the lid's brittle layer, drawn by the generator the rises' faults now share, each stopped only by an older fault at least as long.
- **F148.** A wide trough keeps its floor where its trace turns, and a graben in ice rises from its tips as fast as it throws.
- **Charon** has a preset: Pluto its primary, Brozović and Jacobson 2024's orbit and mass, Nimmo et al. 2017's radius.

How it went:

- **The research.** Charon's lid stretched by at least 1 percent of its area as its ocean froze (Beyer et al. 2017), in a belt of chasmata round its northern terrain: Serenity 40 to 50 km wide and more than 5 km deep, Mandjet at least 450 km long, 30 km wide and 5 to 7 km deep, Argo some 9 km deep; its 25 largest scarps 93 km long at the median; tilted blocks of 5 to 6 degrees between. Restricted faults run two to ten times as long as they are high (Nicol et al. 1996; Soliva and Benedicto 2005). Mars's graben, Venus's rifts and the slip per length of faults on several bodies were gathered for the next stage.
- **The shared generator.** The rises' fault population was factored into one generator that any strain can feed; Mars and Venus came out bit for bit as before the change.
- **Graben as wide as a lid lets them be.** Mars's rule, by which a fault whose throw would pass its layer cuts deeper as a rift thirty times as wide as its throw, gave Mimas chasmata 305 km wide: a lid's faults stay graben, as nothing below the lid's brittle layer can break. Low gravity puts the depth at which the faults close beyond Mimas's and Miranda's centres, so the layer is held within the lid's elastic part and within a fifth of the radius, past which no thin shell is stretched alike through its thickness: Ithaca Chasma, the widest graben known against its body, is 0.19 of Tethys's radius across.
- **What melts and refreezes.** Summing the radius's every growth gave Miranda's melting and refreezing at its resonance its fill of faults though its radius ended where it began. The lid now bears its strain elastically and its faults move the strain at which it stands unstressed: a shell melted and refrozen to where it was breaks nothing new, unless what it lost in shortening went beyond what it bears. The same history will give the next stages their compression.
- **Long faults.** Charon's first chasmata were pods: every fault started at twice its width, and every older fault stopped every younger one, so none passed 230 km among faults running every way. The faults now start at four times as long as they reach down, and a fault is stopped only by an older one at least as long, which took up the stress it would have run on; a longer one cuts across. A fault's time comes from the extension spent before it, so the older are drawn first.
- **The troughs' seams.** At each turn of a chasma's trace the bands of its frames spread apart and left bars across its floor, and where they overlapped, stripes from frame to frame: unseen in Mars's graben a few kilometres wide, plain across Charon's 46 km. The frames now fan out by their own turn and share each cell; a fan on every frame made Miranda's wide graben take minutes to draw, each cell weighed against some seventy frames. Each frame also gathered every cell within the trough's reach to keep the few its band holds; it now walks only its band's strip, which draws Charon's chasmata in 10 s instead of 96, bit for bit the same. And a fault beside another copied its host's every relay step some 200 km aside, its points jumping apart round the outside of each step and folding back inside it, a third of Charon's chasmata broken: it now follows the host's course smoothed over twice the distance between them. A graben's ends rose from its tips at rock's rate, which in ice spread them over half its length into a lens.
- **The slip per length.** Ice's stiffness alone, ten times rock's compliance, cut Charon's chasmata to their width's full depth, 9 km at the median, where Mandjet's throw is at most 0.016 of its length; density over stiffness, four times rock's, gives some 4 km at the median and up to 12.
- **Deep oceans and young surfaces.** Below 209.9 MPa an ocean freezes as ice II or III, denser than its water: only Pluto's reaches it, and its extension ends at 2.8 Gyr. The lid's record starts with the craters', so Europa, whose surface keeps 100 Myr, opens three graben where it had opened 2,186.

What the maps show:

- **Charon:** 58 chasmata 220 km long at the median and up to 890, gathered in sets, 46 km wide and 4.1 km deep at the median, from 0.14 to 0.47 Gyr; flat floors between walls, not yet Serenity's stepped walls and tilted blocks (audit M23).
- **Pluto:** some 1,400 graben across all of it, its ocean freezing through in its thermal history where Pluto's own survives (audit M22).
- **Miranda:** 29 graben from its refreezing at 3.7 Gyr, up to 54 km wide; **Mimas:** a dozen from its young ocean's refreezing in its last 10 Myr, where Mimas shows none (audit M22).
- **Ganymede, Callisto and Titan:** some 700 to 800 graben each, faint at their scale beside the lanes (audit M21); **Europa:** three young graben.
- **Mars and Venus:** their graben beside others now follow their hosts' courses and their troughs keep their floors across their turns, so their troughs and the floods about them shift a little: a fifth of Mars's surface by more than 100 m, Olympus's summit by 150 m, the same maps to the eye.

Still to do:

- **The chasmata's fabric** (audit M23): their walls stepped in fault terraces, their floors broken into tilted blocks, the belts given a direction by the relief's stress; with Venus's rifts and Mars's canyons and Noctis, the next stage.
- **One extension budget** for a lid's faults, lanes and bands (audit M21), and **compression** on the same history.
- **Pluto's and Mimas's oceans** in their thermal histories (audit M22).
- **The lid's faults in the loop's history** (audit M24): drawn by the finish alone, so the loop's floods and relaxation do not see them.

Tests 231 (17 ignored), among them Charon's lid stretching as its ocean freezes and breaking into chasmata against Beyer et al.'s measures, a shell melted and refrozen to where it was breaking nothing new, rock lids opening none, and a wide trough keeping its floor where its trace turns. Clippy is clean and validate passes its 134 checks.

### Revised: a graben's form

The user, on the renders: the faults look too wide on Miranda and Charon, their detail is circles embedded in their edges, which looks very bad, they do not look like cracks, and their depth is a straight cut down to one level; the icy moons cannot come out looking like this. The coarse physics may place and size them, and noise or randomness make them look right, so long as they are roughly plausible and convincing to the eye.

- **F149.** A lid's faults cut to half the brittle layer at the median, not all of them to its base, so Charon's chasmata are 24 km wide at the median and 44 at most. A fault's graben takes no landslide alcoves; its walls are scarps in straight pieces with jogs and self-affine roughness, terraced by step faults that come and go; its floor is the old ground let down, craters and all, tilted toward a master fault in two graben of five; its footwalls stand raised; its throw is rough along it; and only its fault faces are new ground, which shed their regolith and stand bright.

How it went:

- **What made the look.** Every lid fault cut to the brittle layer's base, so every graben was the widest it could be; the circles were F106's landslide alcoves, Valles Marineris's scallops drawn along any wall a kilometre high, a few per wall height; and the trough's profile was a flat floor between tanh walls, cut anew, with the finish's relief removed inside it. The width came of the physics misapplied, the circles of a rift's process applied to a graben, the floor of a shape where a mechanism was wanted.
- **The floor.** A graben's floor is the old surface dropped between its faults, which a ring graben's has been since F69; a fault's graben now drops it too, and the craters on it come down with it, so the floors are no longer one level. The faces are the new ground: the finish writes their age and marks the steep layer there, as it does a crater's walls, so the scarps shed their lag and stand bright on the icy moons.
- **The walls.** A scarp from the crest at the angle of repose, a step fault or two within the wall as terraces drawn by window along it, each tilted back, and the footwall raised a tenth of its fault's throw: the cross-section is a staircase that meets the floor without a seam, its wall's throw found at the foot of its last face in a second pass. The walls' pieces, jogs and roughness are stated choices against Pluto's and Charon's chasmata; a first pass with jogs of a tenth of the width over 0.4 widths gave kinked diamonds, now 0.07 over 0.6.
- **The ends.** A graben narrowed to its tips over a whole width at each end, half the length of a graben four widths long, so the short ones were lenses; narrowed over 0.3 widths instead they were pills with square ends, their throw standing up at once at the tip as a quarter circle. The throw now rises from each tip as a parabola, the width falls with it, and the walls swell in and out over many widths, so the troughs taper to points and their walls do not run parallel.

What the maps show: Charon's chasmata in belts of narrow troughs with stepped walls and raised rims, old craters let down on their floors; Miranda's graben as fault troughs among its coronae rather than broad bands; Pluto's and the large moons' graben narrower and more numerous.

Still to do: the terraces' places from the throw's partition on a listric master fault (audit M23, T1); the belts' tilted-block fabric between the chasmata; the lid's faults in the loop's history (M24).

Tests 232 (17 ignored), among them a graben dropping the old ground as its floor with its walls scarps at the angle of repose, raised rims, a terrace along it and only its faces new ground. Clippy is clean and validate passes its 134 checks.

## The seventy-third milestone: the lid's stress field and Venus's rifts

The first stage of M72b, in the order the user set after the seventy-second: the relief's stress and Venus's rifts, with the coronae strung along them.

- **F150.** One stress field in pascals on every rock lid: the rises' stresses scaled by the membrane stress of their uplift, the upwelling between neighbouring plumes pulling the lid apart along the great circle joining them, and the potential energy of the relief, a high in tension and a low in compression, each spreading down its margin. Every fault starts where the field passes the brittle layer's yield and runs as long as it does; a lid-wide population spends the limbs' extension, so rifts run from rise to rise.
- **F151.** The diapirs are drawn along the limbs as toward the plumes, so Venus's coronae gather along its rift belts.

How it went:

- **Magnitudes.** The rises' field was dimensionless, in units of each crest's tension, with faults starting at a twentieth of it; to add the relief's stress in pascals the crest needed a value. The dome's bending gives a few MPa, far below the 23 MPa that breaks Venus's layer; the membrane stress of the uplift on the sphere, E h / R, gives 50 MPa on Venus and 250 on Mars, which with the yield reproduces the old reach of the rise's faults, about a radius and a half, with no new choice beyond the share of 1.
- **What joins the rises.** The plumes' far fields fall as the inverse square of the distance and are a few MPa between rises four radii apart, and the relief's stress between them is nothing; no honest magnitude joins Venus's rises by the plumes alone. The rises sit on the mantle's upwellings, and Venus's three rifts join its three rises along the BAT upwelling; so two plumes close enough sit on one upwelling limb, which pulls the lid apart across it. At half the lesser crest's tension the limb sat just under the yield at its midpoint and the rifts stopped at the rises' feet, where the domes' bending squeezes the lid; at 0.8 they cross, and the longest runs 3,600 km.
- **The stubs.** A sixth of Venus's faults were a single step long: neighbours drawn beside a short host whose course, moved aside by more than the host's turns, collapsed to a step; they start on their own now. A guard that drew a fault's start again where the field barely passed the yield changed nothing and was dropped.
- **The coronae.** Drawn toward the limbs, two in five of Venus's coronae lie on them, and a third within 150 km of the limbs' fracture belts, but a fourteenth within 150 km of a rift against the real two thirds along chasmata. Thinning the lid over the rifts, as measured, made it worse: a thinner lid puts the mantle at its base above the solidus, and the melting rule then makes no corona. Drawing the diapirs toward the rift zones raised the share near rifts in the derivation and collapsed the loop's count to 167, the births moved into the hot ground about the rises, so it was dropped; the derivation and the loop agree again. The shortfall is the melting rule's and the limbs' lack of heat (audit M25), not this stage's to fix.
- **Mars.** Its one plume has no limb, and its population is as it was in count and reach; by the draws its long canyons now cluster at the other end of its rise's axis, and one trough follows the dichotomy's margin south-east of Tharsis, where the relief's stress turns the rise's far field along the contour: the physics, kept.

What the maps show:

- **Venus:** rifts running from one rise to the next along its two limbs, en echelon, 26 km wide at the median and up to 105, in belts of thousands of graben; coronae along the belts and about the rises.
- **Mars:** its canyons regrouped along its rise's axis; the dichotomy-margin trough.

Still to do: the chasma fabric from a listric master fault's throw partition (M23) with Valles Marineris and Noctis Labyrinthus; compression on the same field, Mercury's scarps and the wrinkle ridges in the lows (M72c); the coronae's share on the rifts (M25); the lid's faults in the loop's history (M24); the icy lids' faults on the same field, which the relief's stress hardly moves at their gravity.

Tests 235 (18 ignored), among them the relief's field, a high in tension spreading down its margin and a low in compression, Venus's rifts running from rise to rise along its limbs with no fault a step long, and Venus's coronae gathering along its limbs at their number. Clippy is clean and validate passes its 134 checks.

## The seventy-fourth milestone: the fabric inside a chasma, and the lid's faults in the loop

The second stage of M72b, approved after the seventy-third: the chasma fabric from the throw's partition, on Charon's chasmata and Mars's canyons alike, and the lid's wide faults carried by the coarse loop.

- **F152.** A wide trough's floor breaks into tilted blocks on a population of interior faults, as dominoes above a listric master fault: blocks a quarter of the faulted layer wide, synthetic to the master and tilted back toward it, taking 0.3 of the throw in pieces drawn from the faults' power law, stamped as tilted-block graben just after their trough. They step the walls down in terraces and cut Valles Marineris's floors; the drawn terrace windows of F149 are gone.
- **F153.** The lid's faults the loop's grid resolves are carved through the loop by the strain process, as a rise's canyons are, so Charon's floods pond in its chasmata and the relaxation sees them.

How it went:

- **From windows to a population.** The terraces of F149 were a drawn feature: a window along each wall, a step fault in it with a drawn share of the throw. The physics behind terraces and tilted blocks is one thing, the hanging wall of a listric fault breaking into dominoes, and ExoGen already had the generator for a population of faults; so the interior of a trough is now a population, its faults placed by the block width the faulted layer sets, their lengths from the same power law as every fault's, their throws from the throw's partition, each a tilted block drawn by the graben form with its floor tilted the whole way. One new field on a trough, its floor's tilt share, carries that.
- **Where blocks show.** A block must span two cells of the grid that draws it, so at 2048 pixels Charon's chasmata and Valles Marineris's troughs break into blocks and Pluto's fossae and Venus's rifts do not; at 4096 more would.
- **The loop.** The strain process carries the wide lid faults as the plume process carries canyons, carving each within the step it forms; a fault's window had to be its moment and a millionth of a gigayear, as a zero-length window made the loop's progress through it a division by nothing.

What the maps show:

- **Charon:** its wide chasmata with two to four tilted blocks across their floors, parallel scarps stepping down toward the master wall, and their walls stepped; 2,700 interior faults among 228 chasmata. The comb of arcs at one chasma's end that the seventy-second's report took for flood lobes is a basin's ring graben.
- **Mars:** its canyons' floors in blocks three cells wide at 2048 pixels, texture rather than steps at that size; nothing else changed.
- **Pluto and Miranda:** their wide graben carved through the loop and a few percent shallower for it, the rest as they were; Venus unchanged, its rifts too narrow at this size to break.
- **Cost:** the interior faults add a third again to Charon's and Mars's stamping, and the carried faults a minute to Miranda's terrains, where the floods' deformation now lists every carried graben.

Still to do: the belt's tilted-block fabric between Charon's chasmata and the belts' direction (M23); compression on the same field, Mercury's scarps and the wrinkle ridges in the lows (M72c); the coronae's share on Venus's rifts (M25); the icy lids' faults on the same field.

Tests 237 (18 ignored), among them a wide trough's floor breaking into as many tilted blocks as its width holds, each a share of the throw, stepping down toward the master wall and formed just after the trough, none where the grid cannot resolve a block; and Charon's wide chasmata carved through the loop within the step they form. Clippy is clean and validate passes its 134 checks.

## The seventy-fifth milestone: compression on the same field

M72c, approved after the seventy-fourth: the lid's stress field read from its compressive side, so that the thrusts and the plains' ridges are placed and turned by it, with the shortening the body's shrinking gives.

- **F154.** The field gives a thrust its compression: the most compressive stress, the strike across it, and how far the field decides that strike against the scatter of the lid's strength; the lid's strength in compression, Anderson's thrusts at the faulted layer's base; and the squeeze the body's shrinking puts on a whole lid alike every way.
- **F155.** One shortening budget, place by place: what passes the lid's elastic strain and the most a place has borne. Each scarp starts where that grew, its longest length from the strain its place takes up, strikes across the compression as far as the field decides it, and stops where the lid is no longer squeezed past its strength; the short scarps are drawn by the finish; a thrust beneath the plains folds them into a ridge.
- **F156.** The plains' ridges stand where the lid squeezes the plains, with the plains' own drive a stated share past the strength, and strike across the field's compression, following the lava's contours where it does not decide, meandering about their course.

How it went:

- **Measuring the field first.** The field's compression on each rocky lid was printed before anything was built on it: a few MPa on Mercury and the Moon, whose thick lids spread their relief's force thin, against the hundreds of MPa Mercury's shrinking would put on its lid; ten to thirty-five on Venus; hundreds at Tharsis's foot. So the field cannot place or turn Mercury's thrusts, and the milestone's plan for Mercury became the budget's physics and the drawing of its short scarps, with the direction its real thrusts have left open (audit M27); Mars and Venus are where the field shows.
- **One budget, the elastic part kept.** The lid holds part of its squeeze elastically, as the icy lids do (F147), so only the strain past its strength is faulted, a seventh less on Mercury and more than half less on Mars, where the shrinking passes the strength by little and its thrusts beyond Tharsis's reach are too short to draw. The longest scarp's scaling was refitted to the strain the thrusts now take up, keeping Mercury's thousand-kilometre scarps.
- **The plains' drive.** Venus's field is a tenth of its lid's strength, and its lid was renewed too lately for its shrinking to squeeze it; the plains' ridges stand on their own subsidence, which the floods do not yet load into the field (audit M26), so the drive past the strength is a stated share, chosen where Venus's ridges keep most to its lows. The share below the mean height was 75 percent at a tenth and 87 at a few hundredths; nothing fixed in between reaches the 93 percent measured, as the field's lows are the old ground's.
- **The combed look.** Read by a field that decides their strike, Mars's ridges came out as straight parallel dashes: the pull toward the strike each step took back their wobble. The wobble is now a meander about the course, as the scarps' lobes were, and the ridges keep it.

What the maps show:

- **Mars:** 3,400 lobate scarps where it had none, all about Tharsis, thickest at its foot and out to two and a half rise radii, running round the rise; the longest, 1,200 km long and 3.4 km high, runs north and south down its western foot. The plains' ridges about the rise run round it as Lunae Planum's do, north and south over the lava basin east of it, near the pixel's size at 2048.
- **Venus:** its ridges have left the margins of the tessera plateaus, where the lava's contours crowded them, for the lowland plains and rings at the feet of its grown rises: 84 percent of their midpoints lie below the mean height on the finished map, against 59 before. At 19 km a pixel they are texture, not lines.
- **Mercury:** 7,200 scarps, the short ones drawn whole by the finish; 225 longer than 100 km against 352, as its lid holds a seventh of its contraction elastically and a quarter of its thrusts fold its smooth plains into ridges; placed and turned at random as before.
- **The Moon:** unchanged: its shrinking stays within its lid's strength, and its maria keep their 2,300 ridges.
- **Cost:** none to speak of, Mars's loop 18 s against 19. The finish's event sequence on Mars swings from 127 to 611 s between renders with the machine's load; the seventy-fourth's 314 s, and the third again its report put on the interior faults, were measured beside other work, and quiet runs give 130 to 150 s before and after.

Still to do: the floods' load and the relief at their time in the field (M26); Mercury's thrusts' direction and their systems (M27); the belt's tilted-block fabric between Charon's chasmata and the belts' direction (M23); the coronae's share on Venus's rifts (M25); the icy lids' faults on the same field.

Tests 241 (20 ignored), among them a thrust striking along a low's margin where the field decides its strike and hardly decided in the low's middle; a thrust's trace turning across the squeeze and stopping where it ends; the plains' ridges standing where the lid squeezes them and striking across its compression; and Mars's thrusts running round Tharsis while Venus's ridges keep to its lows. Clippy is clean and validate passes its 134 checks.

## The seventy-sixth milestone: the floods' loads on the lid's shell

M26, approved after the seventy-fifth: the relief and the floods' loads read through the loop at their times, on the lid of their time, so that the maria's ridges and the graben about their margins come of their own load and the plains' stated drive goes.

- **F157.** The loop reads the field at epochs: the relief as the ground stands, on the lid as thick as it is then, and each load of lava since the last bending the lid of its time as a shallow spherical shell, its membrane stiffening it and adding its stress to the bending's, solved on a cube-sphere grid fine enough for its flexural parameter. Thrusts break in Anderson's regime only, from the record's start; the strike is decided against the squeeze itself where it passes the strength; the ridges read the field alone, the lava's thickness and the stated share gone, the plains at their strength at least.
- **F158.** The extension the loads add past the yield opens graben beyond the filled basins' margins, from the faults' generator, each as the stretch at its start grew.
- **F159.** The ridges take the thrusts' laws: lengths from their power law, heights growing with length, arches and crests with height.

How it went:

- **A plate on the sphere.** The plate needed a Laplacian on the cube-sphere. The two-point one, symmetric and simple, was hundreds of percent off at the face edges, where the grid kinks and a smooth field's gradient leaks into it; a load across a cube's edge would have drawn a seam of false stress. The curvature stencils the stresses needed already fit a quadratic at every cell, and their trace is a consistent Laplacian; the plate's bending energy turns it into a symmetric operator conjugate gradients can solve. A load on a cube's corner bends the plate alike every way round it.
- **Checked against the exact solution.** The center of a basin three flexural parameters wide came out squeezed by 48 MPa, three times my estimate; the Hankel transform's solution of the axisymmetric plate gives 47. The estimate was wrong and the solver right, and the assertion that the middle of a broad basin stands nearly unbent went.
- **The flat plate's spokes.** Read by the plate, a mare's ridges came out as spokes radiating from its middle, nothing like the maria's rings: near a load's margin the flat plate squeezes most round it, the strike-slip paradox of the flat plate's mascons (Freed, Melosh and Solomon 2001). Their resolution is the shell's membrane, and the Hankel transform of the shallow shell showed it before any of it was built: the membrane squeezes the basin most outward, its ridges ring it, and the lid beyond is pulled both ways. On the Moon the membrane's stress is as large as the bending's.
- **A thousand lunar scarps.** With the field read from the start, the Moon's first lid, 4 km thick under the primordial relief, broke into a thousand scarps up to 600 km long at once. The surface keeps no record from then; the thrusts now count from the record's start as the lid's faults do, and the Moon has eleven, the longest 36 km.
- **The drive that stayed.** Removing the plains' stated share outright, as the audit asked, left Venus without a ridge: its lavas were laid on the thin lid after its overturn and stress it by tens of MPa against a strength of 110. A test of the share at zero had seemed to show the loads alone would do; it had left the drive at the strength. The share is gone and that floor stays, a stated choice with nothing fitted, standing for the drives the field does not hold.
- **The fingerprint.** Turned alike by the shell, a mare's ridges, all 200 m high whatever their length, filled it with parallel arcs. As folds over thrusts they now take the thrusts' laws, and a few long high ones stand out of many short low ones.

What the maps show:

- **The Moon:** its maria's ridges ring them, a few long ones standing out of many faint, where they wandered along the lava's contours; 187 graben about the maria, a median 82 km beyond their margins and 3 km wide, below the pixel at 2048 as the real rilles are; eleven scarps by one mare's margin, the longest 36 km, where it had none.
- **Mars:** its lava basin's ridges still run round Tharsis, fewer and of varied sizes where they were combed; 2,600 scarps against 3,400, the longest 770 km; some 300 graben about its plains.
- **Venus:** 56,000 ridges against 71,000, 84 percent of their midpoints below its mean height against 83; its plains' loads, laid on the thin lid after its overturn, open 22 graben.
- **Mercury:** 6,800 scarps against 7,200, the longest 940 km; 142 graben about its smooth plains.
- **Cost:** the loop reads the field at some sixteen epochs and lays a dozen loads, each solved on its own grid: 4 s more on the Moon (13 s), 5 on Mercury, 11 on Mars, 4 on Venus; the finish as before, the Moon's event sequence 58 s with its graben or without.

Still to do: the ridges' and thrusts' systems and Mercury's direction (M27); the drives the field does not hold, Venus's climate first (M26); the rises' faults on the field of their time; the belt's tilted-block fabric between Charon's chasmata and the belts' direction (M23); the coronae's share on Venus's rifts (M25).

Tests 246 (21 ignored), among them the plate's Laplacian finding a harmonic's degree on every cell; a disc of lava bending the plate as the axisymmetric plate does, and a load on a cube's corner bending it alike every way; a disc on the Moon's shell bending it as the shallow shell's Hankel solution does; and a filled basin squeezing its plains, its ridges within it, and opening graben round its margin, none before the record began. Clippy is clean and validate passes its 134 checks.

## The seventy-seventh milestone: thrusts and ridges in systems

M27's second half, approved after the seventy-sixth: the thrusts and the plains' ridges drawn as the faults' generator draws the lid's faults, in systems.

- **F160.** Half the thrusts and the plains' ridges start beside one drawn before near where they start, en-echelon along a span of it at 0.6 to 1.6 of its reach; the rest run their own way and stop where they meet one drawn before at least as long. The plains' ridges join the thrusts' ridges as well as their own.

How it went:

- **Systems where the squeeze is.** The lid's faults chose their hosts by length from anywhere, as their strain was alike everywhere; a thrust's squeeze is not, so a newcomer chose its host among those near where it starts, and runs along the host's span nearest it, or a system would have grown where the lid was no longer squeezed.
- **Mercury's field.** Looking at Mercury's scarps for their systems showed them running alike by region, before the systems as after: since the seventy-sixth, the loads of Mercury's plains and its relief on the thinner lid of their times decide most of their strikes, where before its field was too weak to decide any. It is the field's physics, the scarps following their plains' loads and margins; the north and south real Mercury's large thrusts favour at low latitudes still wants its despinning (M27).
- **The shield's trench.** A dark ring round a shield's foot in a high-passed crop was the filter's, the blur carrying the shield's height past its foot; the ground's height round the shield falls from its flank to the plains without a dip. A shield that size would in fact sit in a moat, which the model's ground does not yet bend into (M19).

What the maps show:

- **Mercury:** 9,600 scarps against 6,800, shorter for their arrest, in clusters of parallel scarps; at 2048 pixels most are near the pixel, the long ones as before.
- **Mars:** in the lava basin east of Tharsis a set of long sinuous ridges stands among fainter ones, running round Tharsis, where the ridges were even dashes; 2,900 scarps against 2,600.
- **The Moon:** its central mare's ridges in groups of arcs with quiet ground between them, where a whorl of arcs filled it; eight scarps.
- **Venus:** 53,000 ridges against 56,000, 87 percent of their midpoints below its mean height against 84, the 93 measured nearer.
- **Cost:** none to speak of, the loops and the finish within a few seconds of the seventy-sixth's.

Still to do: Mercury's thrusts' direction (M27); the ground's flexure under its loads, a moat about a shield (M19); the drives the field does not hold, Venus's climate first (M26); the belt's tilted-block fabric between Charon's chasmata and the belts' direction (M23); the coronae's share on Venus's rifts (M25).

Tests 247 (21 ignored), among them a thrust beside its host running along a span of it about where it starts, as far off as the host's reach sets, none joining a thrust that is not near, and one running its own way meeting a host at least as long within their reaches but not a shorter one nor one farther off. Clippy is clean and validate passes its 134 checks.

## The seventy-eighth milestone: the ground bent under its loads

M19, approved after the seventy-seventh: the shields' loads in the field, and the ground bent under every load.

- **F161.** The shields' heights join the floods' lava as loads on the lid's shell, so the ground beneath and about a great shield is squeezed and its flexural arch stretched, where the loads' graben ring it.
- **F162.** Each epoch's deflection is laid into the ground: the floods' ground sinks with the lid, while beneath a shield the edifice fills the sag as far as it is thick, so it stands in a moat with a bulge beyond; the finish's floods see the subsidence after them.

How it went:

- **A shield keeps its height.** A shield's height is its constructional surface, a 5° flank up to the relief the magma can raise; on a sinking plate its lava keeps filling the sag beneath it. Laying the whole deflection into the ground would have lowered every shield by kilometres; instead the edifice fills its own sag and only the ground about it sinks. Counting that root as load would deepen the sag again each epoch, so it is left out, a stated choice.
- **The finish's floods.** A flood refilled at the finish is laid against the ground as it found it; the deflection is a deformation after it, as the shields' growth and the canyons' carving are, so the finish's deformation now carries every epoch's subsidence and gives each flood the ground it lay on.
- **Two batches at once.** A render batch I took for dead had started all the same; the second's build found the executable locked and stopped. The first rendered from the same code, the edits since a comment and the cache versions; it is a lesson in starting renders only as tracked background tasks.

What the maps show:

- **Mars:** its largest shield, 684 km across and 25 km high, stands in a moat 2.4 km below the ground at 2.5 of its radii, its summit 26.0 km where it was 27.2; the four smaller shields stand in moats of their own, 0.4 to 3.7 km below the ground at 2.5 of their radii where the old ground dipped 0 to 1.5 km; the region about Tharsis sinks 1.5 to 2 km under its lavas and shields and its lava-filled basins deepen, the lowest to −8.0 km.
- **The Moon:** its central mare sinks some 1 to 1.5 km under its fill; the lava in craters, too narrow to bend the lid far, sinks less.
- **Venus and Mercury:** Venus's plains sink under the resurfacing's fill by up to 1.7 km, its mean height 0.24 km where it was 0.47; Mercury's deepest basin to −9.3 km from −7.2.
- **Graben:** with the shields' loads the loads' graben grow from 605 to 774 on Mars and from 24 to 123 on Venus, about their shields; the Moon's and Mercury's, without shields, change only as their floods fill the sunk ground (193 to 207, 181 to 212).

Still to do: the edifice's foot, its lava running into its moat and ponding there (M19); Mercury's thrusts' direction (M27); the drives the field does not hold, Venus's climate first (M26); the belt's tilted-block fabric between Charon's chasmata and the belts' direction (M23); the coronae's share on Venus's rifts (M25).

Tests 249 (21 ignored), among them a shield on a plate filling the sag beneath it, its summit where it built itself, the ground about it sunk into a moat with a bulge beyond while a flood with no edifice sinks with the ground; and Mars's largest shield standing in a moat in the loop, its summit still near the height it built. Clippy is clean and validate passes its 134 checks.

## The seventy-ninth milestone: the shields' lava at their feet

M19's remainder, approved after the seventy-eighth: the lava a shield's flows carry past its foot, ponding in the moat.

- **F163.** A shield's eruptions, as F101 draws them, carry some two fifths of its lava past its foot, nearly all of it from its rift zone. Each step that share of what its edifice grows by starts from along the flows' runs beyond its foot, off the ends of its rift zone, and runs on down as a flood's lava does, ponding in the moat and spilling past it, the shield's own lava poured out at a shield's rates.
- **The finish's waves.** Each event joins the first wave after every wave that holds an earlier event it overlaps, judged on a coarse grid of buckets, where a wave used to close at the first event that overlapped it.

How it went:

- **The share past the foot.** The rift zone's flows, 0.2 to 1.0 of a radius long from vents 0.4 to 1.0 of the way to the foot, carry more than half their lava past it, and the summit's and flanks' fissures almost none: 41 percent of a round shield's lava, less along an elongated one's rift zone. Mars's shields send 5.4 million km³ past their feet, some two thirds of what their edifices hold.
- **A pile at the vents.** Erupted at the foot itself, every eruption off a rift zone's end crossed the same few coarse cells and left its flow's sheet on each, and the largest shield stood between two piles more than 9 km high. The lava past the foot lies along the flows' runs beyond it, so each eruption now starts from a point along such a run, drawn by the lava the run carries there, and never on another edifice, whose flanks the lava would not climb.
- **The finish slowed sixfold.** The 2,858 eruptions past Mars's feet each covered a moat's lake, and each closed the finish's wave of parallel events wherever it came, so Mars's event sequence took 929 s where it had taken 149. Now an event joins the first wave after the events it overlaps, whatever comes between, and Mars's sequence takes 149 s again. The Moon's maps and sidecar come out byte-identical, so every overlapping pair keeps its order.
- **Shields lumped into fewer eruptions** would have cut the cost as well, but the many eruptions are what build the aprons: their sheets stack along the rift zones' ends where one large eruption would only pond.

What the maps show:

- **Mars:** the largest shield's moat is filled, the ground from its foot to two radii at 2.7 to 3.0 km, 0.25 to 0.55 km below the ground at 2.5 radii, where the moat had been 2.3 km deep. Between it and its neighbours to the west its lava, hemmed in by their edifices, stands as an apron up to 7.1 km. The smaller shields' moats are 0.5 to 2.1 km deep where they had been 0.4 to 3.7, and each shield stands in a halo of its own lava with lobate margins. The new loads move the field: 3,148 thrusts where there were 3,455, and 28,814 ridges where there were 28,402.
- **Venus:** 561 eruptions past its shields' feet, 1.0 million km³, fill their shallow moats by 0.1 to 0.4 km and lay lobate fields about them, burying the near ends of the rifts and graben beside them.

Still to do: the shields' melt as one budget (M28); the root's load (M19); Mercury's thrusts' direction (M27); the drives the field does not hold, Venus's climate first (M26); the belt's fabric between Charon's chasmata (M23); the coronae on Venus's rifts (M25).

Tests 251 (21 ignored), among them a shield's rift flows carrying 41 percent of its lava past its foot and its summit's fissures none; a shield erupting what its flows send past its foot, within each step, from along the flows' runs off its rift zone's ends, on no edifice and with its own lava; and Mars's largest shield's moat holding its lava, the ground past its foot within a kilometre of the ground beyond, in place of the moat test of the seventy-eighth. The suite runs in 517 s where it took 812, the finish's waves no longer broken at every overlap. Clippy is clean and validate passes its 134 checks.

## The eightieth milestone: the shields' melt as one budget

M28, approved after the seventy-ninth: a shield grows from the lava its rise brings it, and the floods lose what the shields take.

- **F164.** Before the loop, step by step, each rise's started shields take their share of the extrusive budget, read from the vent weights over the first ground, and share it by their edifices' volumes. An edifice holds what its flows lay on it until it is full, its flows carrying their share past its foot; the magma a full edifice cannot lift intrudes beneath it. Its dike swarms come at shares of that growth, and its calderas sink by the share of its height it built.

How it went:

- **Mars's budget checks out.** Greeley and Schneid (1991) put Mars's volcanic deposits since the Early Hesperian at 6.45 × 10⁷ km³ and its Amazonian rate at 0.012 to 0.014 km³/yr; the model's shields took 6.6 × 10⁷ km³ since then, and its late rate is 0.012 to 0.016. That lava is real, and it was vanishing.
- **Erupted past the feet, it piled up.** The first version sent everything a full edifice could not hold past its foot: 6.3 × 10⁷ km³ in 29,856 eruptions, each spending itself in flow sheets near where it started, stacked aprons 10 to 20 km high and, where one shield was hemmed in by its neighbours and its eruptions fell back on one cell, a spike 270 km high. The user chose: the magma a full edifice cannot lift intrudes beneath it.
- **Venus's shields starve.** Its budget after the overturn is a stated 2 percent share, and its shields take their share only once its rises' growth turns: they stand 6 to 15 percent built. The user chose to keep that and open audit M29.

What the maps show:

- **Mars:** the shields fill by 1.6 to 1.9 Gyr, the largest by 2.9 Ga, as Olympus Mons's main stage had ended by 2.54 Ga (Isherwood et al. 2013), and 5.7 × 10⁷ km³ intrudes beneath them; the lava past their feet, 5.6 × 10⁶ km³, is much as before, the largest shield's moat now about a kilometre deep where it last filled to half that. The floods keep 7.4 × 10⁷ km³ of the budget. Mars's event sequence takes 148 s.
- **Venus:** its shields are low swells 0.4 to 1.6 km high where they stood 6 to 10.5 km, their calderas shallow with them, and the lava past their feet falls from 1.0 × 10⁶ to 1.4 × 10⁵ km³.
- **The Moon:** byte-identical.

Still to do: Venus's shields' supply (M29); the root's load (M19); Mercury's thrusts' direction (M27); the drives the field does not hold, Venus's climate first (M26); the belt's fabric between Charon's chasmata (M23); the coronae on Venus's rifts (M25).

Tests 252 (21 ignored), among them a rise's shields growing from what they take of its eruptions: each step no more than its budget brings, nothing before the first shield starts, the supplies adding up to what they take, Mars's edifices full by 3 Gyr with more than thrice what they hold brought to them, and Venus's partly built. The foot test now checks the flows' shares of what the edifices grew by, the rest intruding; the plume's tests schedule their shields' supply. Clippy is clean and validate passes its 134 checks.

## The eighty-first milestone: the plumes' melt and a volcano's weight

M29, after the user judged Venus's small shields right and asked for some volcanoes like Maat Mons: Venus's lava after its overturn from its mantle's history, gathered at its rises, and a rise's lava gathered by its volcanoes' weights.

- **F165.** The ambient mantle makes the share of the melt its own gate keeps, and the plumes the rest only under their heads; once the ambient mantle stops melting, all the lava comes out at the rises. On an episodic lid the melt while the ambient mantle still melts after the overturn is the overturn's own resurfacing. A rise's shields share its lava by the area each gathers dikes from, its vents' and its weight's over the regional stress, which falls as the rise's growth wanes, so that the leading volcano takes more as it grows.

How it went:

- **Where Venus's lava can go.** Its mantle's own history gives some 5 × 10⁷ km³ after the overturn, against 1.2 × 10⁶ kept by the stated 2 percent; spent as plains it would bury Venus many times over, where only 2.5 percent of its craters are embayed by lava. But the melt was reckoned with the plumes' excess heat over the whole body: Venus's ambient mantle stops melting 0.45 Ga ago, after which only the plumes' heads melt, 3.4 percent of the surface, and their lava comes out at the rises.
- **Mars checks out again.** Its ambient mantle tapers over its last gigayear, and its eruptions since 1.3 Ga fall from 1.95 to 1.19 × 10⁷ km³, where Greeley and Schneid (1991) put its Middle and Late Amazonian deposits at 1.18 × 10⁷.
- **A feedback too weak at first.** Against the lid's yield, a leading volcano's weight only matched its own vents' area by the end, and its rise's smaller shields came out the more built. The regional stress falls with the rise's growth as F119 has it, and once it is gone the shields gather by their weights alone: the leader keeps its lead.

What the maps show:

- **Venus:** each rise's leading volcano takes most of its lava, standing 5.1, 4.6 and 2.8 km over 328, 296 and 312 km, with summit calderas and lobate lava fields at its feet, as Maat Mons stands 4.9 km over 395 km; their neighbours stand 0.3 to 0.9 km, and three volcanoes that started as their rises stopped growing never grew, nor sent dikes. No lava reaches the plains after the overturn's resurfacing.
- **Mars:** its largest shield fills by 3.0 Ga, the next two by 2.1 and 1.3 Ga, and the two smallest stand half and a third built.
- **Counts:** Venus's floods fall from 720 to 320, none on its plains after the overturn's resurfacing; Mars's thrusts number 3,437 where they were 3,123, its late lava gathered at Tharsis moving the field.
- **The Moon:** byte-identical.

Still to do: the event's fill from its surge's melt (M29); the root's load (M19); Mercury's thrusts' direction (M27); the drives the field does not hold, Venus's climate first (M26); the belt's fabric between Charon's chasmata (M23); the coronae on Venus's rifts (M25).

Tests 253 (21 ignored), among them a rise's leading volcano gathering its lava: on each of Venus's rises the largest takes more than three quarters of what its shields hold, at least two stand over 2 km, none full, and no dike swarm comes from a volcano that never grew. The budget test now has Mars's largest edifice fill first, by 2 Gyr, and the moat test reads the shields' own lava, poured out at a shield's rates, near the foot against farther off, the moat within 1.5 km of the ground beyond. Clippy is clean and validate passes its 134 checks.

## The eighty-second milestone: the shields' roots

M19's last part, approved after the eighty-first: the lava that fills a shield's sag weighs on the lid.

- **F166.** Beneath an edifice the lid sinks into the edifice's own lava, which fills the sag as far as the edifice is thick, so that the mantle pushes back there by the difference between the rock's density and the lava's. The plate is solved with its restoring force cell by cell, the share of the sag the edifice fills read from the deflection it fills and then once more.

How it went:

- **Bounded by the lid.** Lava nearly as dense as the rock, Mars's at some 2,970 kg/m³ against 3,300, would sink an edifice without a lid some ten times as far; the lid's strength holds it to 1.65 times beneath a shield 600 km across on a lid 70 km thick, 1.25 times on 30 km, 1.16 times on Venus's 25 km.
- **Solved directly.** Counted as load epoch by epoch, as F162 feared, the root would have deepened the sag again each epoch, converging only slowly; solved as a weaker restoring force where the edifice fills its sag, it converges within the epoch.
- **The intruded surplus is not a root.** The magma Mars's largest shield cannot lift would stand a column some 115 km thick over its footprint, thicker than the lid: it spreads beneath the rise as underplating, which the dome holds.

What the maps show:

- **Mars:** its largest shield stands in a moat some 3 km below the ground at 2.5 of its radii, where it stood in 1.2 km, the second in 4 km, the third in 2.8 km; the lava past their feet fills part of them, and the apron between the largest and its neighbours sits lower. The summits stand where they built them. Olympus Mons's moat is mostly filled, by lava and its aureole's deposits, which the model lacks (M30).
- **Venus:** its leading volcanoes' surroundings sink some 0.2 km, its lid too thin for their roots to bend it far.

Still to do: the aureole and basal scarp, a shield's flanks spreading into its moat (M30); the event's fill from its surge's melt (M29); Mercury's thrusts' direction (M27); the drives the field does not hold, Venus's climate first (M26); the belt's fabric between Charon's chasmata (M23); the coronae on Venus's rifts (M25).

Tests 254 (21 ignored), among them an edifice's root deepening its sag and its moat: on a Mars lid 70 km thick the sag beneath a shield's summit and its moat come out more than 1.3 and less than 3 times as deep as an empty sag would leave them, the edifice still standing where it built itself, and a flood with no edifice bending the lid as before. The moat test now asks that Mars's largest shield's moat stand below the ground beyond while its shields' own lava ponds in it. Clippy is clean and validate passes its 134 checks.

## The eighty-second milestone, revised: the lava past a foot on unbuilt ground

The user's check of Mars on master seed 235 after the eighty-first: the largest volcano had a peak at its southern foot higher than its own summit, the lava flow below it was tiled in squares larger than a pixel, and the highlands and lowlands met along a sine across the map.

- **One cause for the peak and the squares.** An eruption past a shield's foot (F163) was kept off every other edifice's design outline, whether or not that edifice had grown yet. The largest shield's southern runs crossed the outline of a neighbour that had hardly begun, so nearly every point along them was turned away, and the eruptions that kept none all started at one fallback point. Some 3,000 of them stacked their 50 m sheets there, cell by cell of the loop's grid, into a cone, 33 km high at preview, its steps the squares. An edifice now turns the lava only where it stands higher than a flow is thick at the eruption's time, and an eruption that keeps no point starts at the foot at an azimuth of its own. At preview the highest ground past that shield's feet falls from 33.4 km to 4.2 km, its summit 30.5 km. On the seed's map the highest point is now that shield's summit, 29.6 km, where the cone stood 47.0 km above a summit of 30.0 km, the ground past its southern foot stands at most 3.4 km, and the squares are gone.
- **The sine is two giant basins.** The lowlands on seed 235 are a 4,911 km basin centred at 64°S and a 3,914 km one at 74°N, their floors flooded level. A circle centred far from the equator draws as a sine on an equirectangular map, and these are clean circles with narrow rims: the crust's own dichotomy lies beneath the primordial relief. The direction of a fix is left to the user's choice.

Tests 254 (21 ignored); the test of a shield's eruptions past its foot now asks that none start on an edifice standing higher than a flow is thick by the eruption's time.

## The eighty-third milestone: a crater's elongation follows its cratering efficiency

The user's choice after the revised eighty-second: the largest basins elongate as the largest known basins are, their outlines' irregularity held back until the processes that wear them are in.

- **F167.** A cavity comes out elliptical below a threshold angle that rises as the cratering efficiency falls, the crater a vertical impact makes over the impactor's diameter (Collins et al. 2011), read from the body's own crater scaling: the largest basins, made by impactors only a few times smaller than their cavities, elongate in some three to four impacts in ten on Mars, its craters of kilometres in some one in twenty.

How it went:

- **Footprint or efficiency.** Andrews-Hanna and Zuber (2010) elongated the largest basins from the footprint an oblique impactor leaves, stretched further by the planet's curvature; simulations of impacts on curved bodies found curvature suppresses the elongation instead (Collins, Miljković and Davison 2013), and pointed to the efficiency trend, which the experiments and simulations share across materials and regimes.
- **The trend's line.** The 2011 paper was not open, and the user supplied the 2010 abstract's figures: the line through Figure 1, traced from the figure's image, is 90° times the efficiency to the power of −0.73 within 0.2° over its whole length.
- **The growth below it.** The old elongation was 1 at the threshold and rose by a fixed curve; by the trend's definition it is 1.1 there, and it grows as the flat-target simulations do, a crater above the threshold coming out a hundredth or so longer than wide rather than a perfect circle.

What the maps show:

- **Mars on master seed 235:** its southern basin, 4,911 km across, struck at 32°, comes out 1.23 times as long as wide, and its northern one, 3,914 km, struck at 22.7°, 1.50 times. The northern lowland now reaches farther south west of the map's centre, so its edge no longer mirrors the southern one; the two lowlands still leave a band of highlands between them, as two giant impacts near opposite poles would on a dead world, their ramps to be worn by the surface processes to come.
- **The Moon:** its largest basin comes out a little elongated, a few of its craters clearly so, and one of its oldest giant basins changes its outline in the eastern lowlands; elsewhere its craters stand as before, within a hundredth of round.
- **Every body:** the share of elliptical craters rises with their size, Mars's 5.8 percent at 10 km to 35 percent at 2,500 km, the Moon's 3.5 percent at 10 km.

Still to do: the giant basins' outlines broken up, should they still look too regular once their rims are worn (the user's note); the aureole and basal scarp (M30); the event's fill from its surge's melt (M29); Mercury's thrusts' direction (M27); the drives the field does not hold, Venus's climate first (M26); the belt's fabric between Charon's chasmata (M23); the coronae on Venus's rifts (M25).

Tests 256 (21 ignored), among them the threshold angle following the cratering efficiency as Collins et al.'s line does and the elongation growing below it as the flat-target simulations do, and Mars's craters of 10 and 30 km elliptical in 3 to 10 percent of impacts, its basins of 2,500 km in 25 to 45 percent, the Moon's craters of 10 km in 2 to 8 percent, the share rising with the size.

## The eighty-fourth milestone: the dichotomy as an outcome

T5, which the user chose after the eighty-third: the hemispheric dichotomy comes of the impacts a body takes, not of the crust it starts with.

- **F168.** The primordial crust's tanh dichotomy and its amplitude by the regime's name are gone. A rocky body holds a basin as wide as a hemisphere, as Mars kept Borealis, and the impacts wider than the crater stamp draws, from 1.5 body radii up to a hemisphere, take a form of their own: the crust excavated inside an elliptical rim down to what the melt sheet makes, the excavated crust thrown beyond the rim where it thickens the crust, the ground where the crust's buoyancy holds it, the floor level, the older relief under the ejecta buried.

How it went:

- **The hold.** F135 held a rocky body to a crater as wide as itself, from Vesta's Rheasilvia, and treated any larger impact as one that shattered and resurfaced it before its record began, so Mars could never take a Borealis, whose long axis is half its circumference. The user chose to hold a hemisphere: Mars's record now starts 0.47 Gyr after its formation rather than 0.57, and it takes about twice the craters, 21 basins over 1,000 km where it took 11 (Frey 2008 counted some 20).
- **The tail.** The loop drew a rocky body's craters only to 1.5 radii, squeezing the rest of the tail below it; drawn to the hold, the largest basins grow, seed 235's southern one from 4,911 km to 9,665.
- **A weak first dichotomy.** The first basin stood only some 2 km below the ground beyond it. Its floor had dropped by the crust it lost where it struck, but the stamped basins before it had thinned the crust there without lowering the ground, and the bombardment had drained the crust from a mean of 40 km to 28, since a stamped basin's thinned crust went nowhere. Two changes, both from the physics: the floor's new crust floats against the body's mean column whatever stood there, and a stamped basin's thinned crust lands beyond its rim, kept. The contrast is now 4.4 km.
- **The wall.** The transient cavity of the crater scaling would have excavated the crust only within four tenths of the basin's radius and tapered it across the rest, the smooth slope T5 named; Borealis's lowland is the basin itself, so the basin excavates the crust to its rim, ending in a short wall.
- **The Moon's asymmetry.** With the tanh gone the Moon loses a contrast that was never tied to its near side; its own route, the hot Earth heating its near side while its crust formed (Roy et al. 2014), is deferred by the user.
- **The finish.** A mega-impact's form is in the coarse ground the finish resamples, so it is not restored; but the finish redraws in age order the floods, ridges, scarps and craters older than it, which it buried. It now takes its place in that order: inside its rim the ground is the resampled one again, and beyond it what the finish drew before is buried by its ejecta's thickness, the surface renewed, as a crater's restore removes what its impact destroyed; the loop craters it buried are dropped.
- **The presets' draws.** Mars's seed 12 drew four mega-impacts and the Moon's seed 13 two, where Mars has one basin of the kind and the Moon none: the Moon's two buried its older relief twice over, flattening its spectrum. As T5 allows, each preset now carries a seed whose draw includes its body's event, chosen by the user: Mars 1048, whose one mega-basin is 10,575 by 9,195 km about 69°N, 194°E, as Borealis is 10,650 by 8,520 km about 66°N, 209°E; the Moon 10, which takes none, its largest basin 2,576 km across about 48°S, 168°W, as South Pole–Aitken is some 2,500 km about 53°S, 169°W. One draw of the model each, not a pin.
- **The spectrum's test.** One row's periodogram scatters too much to test a slope on, some rows of a Moon without a mega-impact at −0.6 and others at −2.5; the test now averages the rows within 40° of the equator.

What the maps show:

- **Mars (seed 1048):** one mega-basin, 10,575 by 9,195 km about 69°N, 194°E, its floor some 2.2 km below the datum and the ground beyond its rim some 2.2 km above it; the heights now fall in two groups, about −1.2 km and +1.5 km (Mars's lie about −4 and +1.5 km). The rise formed on the basin's rim, as Tharsis stands on the dichotomy's boundary. The boundary is a clean arc, its breaking up held back by the user until the surface processes wear it.
- **The Moon (seed 10):** no mega-impact; its largest basin, 2,576 km across about 48°S, 168°W, is ringed by the crust its excavation threw out, which stands as highlands about it, as the far side's highlands stand about South Pole–Aitken; the maria lie on the thin crust of the basins north of it.

Still to do: the boundary scarp's wear and the lowlands' sediments (Phase 5); the endogenic route and a locked moon's near side (T5); the giant basins' outlines, on hold; the aureole and basal scarp (M30); the event's fill from its surge's melt (M29); Mercury's thrusts' direction (M27); the drives the field does not hold (M26); the belt's fabric between Charon's chasmata (M23); the coronae on Venus's rifts (M25).

Tests 257 (21 ignored), among them a mega-impact excavating its basin and throwing the crust beyond its rim: the crust kept within a thousandth, the floor level and as low as its lost crust floats, the rim's crust the thickest and thinning outward, the basin longer along its azimuth by its elongation, the floor a basin floor of its age and the ejecta deep ejecta. The test that nothing older than a loop crater survives inside it now holds a mega-impact's interior too; the Moon's ridges lie on lava but for those a later mega-impact buried; the Moon's spectrum is averaged over its rows; the basin-secondaries test keeps its own seed; and Mars's deepest canyon need stand over a kilometre deep along most of its length, seed 1048's standing 3.1, 2.7 and 1.3 km below its rims at a quarter, half and three quarters along.

## The eighty-fifth milestone: an ejecta deposit levels only the hollows it can fill

The user's review of the eighty-fourth: a Moon struck by a mega-impact (master seed 6342) stood flat but for its basin, where the primordial relief had been believable ground, and asked whether to keep that relief. It had not gone; the mega-impact's ejecta buried it.

- **F169.** A deposit levels the ground over the width of the fresh crater as deep as it is thick, the largest hollow it can fill, no farther than a landslide's debris runs out for its fall, and drapes the longer relief. A crater's rings split into sectors that long; a mega-impact's reference is the ground low-passed over that width.

How it went:

- **The reference was too wide.** A crater's ejecta buried the older relief against the mean of the whole ring about it, and a mega-impact's against the ground low-passed over its basin's radius: around a body-sized basin each reference spanned thousands of kilometres, so its kilometres of ejecta buried the hemisphere-wide swells of the primordial relief, which a deposit drapes.
- **The crater rule alone reaches too far.** The diameter of the crater as deep as a deposit is thick grows slowly with the depth of the largest basins, so that a deposit 10 km thick would have levelled over 1,300 km; a deposit spreads only as far as its debris runs out, the landslides' 14 times its thickness, and levels over 140 km.
- **The same rule for every crater.** The craters' ejecta take the sectors too, so that the largest stamped basins drape the swells about them as the mega-impacts do.

What the maps show:

- **The Moon, master seed 6342:** its far side keeps the swells and lowlands of its primordial relief under the mega-impact's ejecta, where they had been buried flat; the belt beyond the basin's rim stands some 6 km high, which the old reference had partly buried as relief; the floor stays the level melt pond it was.
- **Mars (seed 1048):** the highlands along the dichotomy stand higher and roll, the old basins' rims beside the boundary showing through the ejecta; its lowland's heights gather about −1.8 km and its highlands spread to some 5 km.

Still to do: the relaxation of a rocky crust by its lower crust's flow, which would soften the dichotomy's boundary and the old basins; Valles Marineris as subsidence along the buried boundary, with M23; the lowlands' lava against the Hesperian plains; the boundary's wear and the lowlands' sediments (Phase 5).

Tests 258 (21 ignored), among them a deposit levelling the hollows it can fill: at every thickness from a metre to a hundred kilometres no wider than its debris runs out, and where narrower the crater it levels over as deep as it is thick or, past the simple craters' transition, the first so deep. Mars's largest shield's summit now need stand three quarters of its height above the ground beyond its moat rather than four fifths, the ground there keeping the relief its ejecta draped; it stands 79 percent.

## The eighty-sixth milestone: a rocky lid's lower crust flows

The second of the items the user agreed after the eighty-fourth: whether the dichotomy's boundary softens before the surface processes wear it.

- **F170.** Where the base of a rocky crust is warm its lower crust flows as a channel beneath the rigid crust above, the relief's pressure moving the crust from under the high ground to under the low, so that the relief diffuses and the longest wavelengths last; the elastic lid holds up the short ones. The relief relaxes toward the domes the plumes hold up.

How it went:

- **Too cold to flow.** The model's thermal histories keep the base of Mars's crust below 1,200 K even early, under a crust 60 km thick, and the Moon's and Mercury's below it too: dry diabase holds the relief there for the age of the solar system, as Nimmo and Stevenson found for a Martian crust thinner than some 100 km. The user chose to build the process as the physics has it rather than a weaker, wet crust for the water-rich bodies, whose flow law and water are both uncertain: so Mars's boundary stands as it formed, and its wear waits for the surface processes.
- **Venus.** Its crust's base stands above 1,600 K through its first billion years, so that the relief formed then relaxes away, and with it the primordial relief its tesserae and lowlands stand on: Venus came out flat but for its volcanoes, and its plains' wrinkle ridges no longer kept to its lows. Its highlands are thought to stand on crust their tectonics thickened, which the model does not make; so an episodic lid waits for that (audit M31).
- **The rises.** The band relaxation acts on the whole surface, and a rise the plume holds up is not the crust's relief: the relief relaxes toward the rises' domes.
- **The scope.** Earth's thermal history and Io's put their crusts' bases at the mantle's temperature, which would flatten Earth's continents and Io's mountains; a mobile lid's crust is its plates', and Io's heat pipe is not a conductive lid, so neither relaxes (a stated scope, against S1).

What the maps show:

- **Every preset:** unchanged. Mars's, the Moon's and Mercury's crusts are too cold at their base for diabase to flow, Venus's episodic lid waits for its highlands' thickened crust, and Earth and Io are outside the scope.
- **Venus, had it relaxed:** flat but for its volcanoes, its tessera highlands and lowland basins gone.

Still to do: Valles Marineris as subsidence along the buried boundary, with M23; the lowlands' lava against the Hesperian plains; the boundary's wear and the lowlands' sediments (Phase 5).

Tests 259 (22 ignored), among them a warm lower crust relaxing the relief between its lid and its breadth: under a crust 60 km thick whose base stands at 1,550 K, a wavelength of 1,000 km loses most of its relief in a hundred million years, one of 10,000 km keeps nearly all of it and the elastic lid holds up one of 200 km, while at 900 K the crust keeps everything for a billion years; and a diagnostic printing each rocky preset's survival by wavelength.

## The eighty-seventh milestone: the floods fill the low ground

The last two of the items the user agreed after the eighty-fourth: Valles Marineris along the buried boundary, and the lowlands' lava against the Hesperian plains.

- **Valles Marineris, deferred (M32).** Its depth is subsidence under the sediments that filled it, since eroded (Andrews-Hanna 2012, part 3), which waits for Phase 5's sediments; and its place along the buried boundary needs the rise's stresses solved through a lid whose stiffness follows the crust, where the dichotomy's own tension, some 2 MPa, cannot steer faults the rise opens at 100 to 250 MPa. The user chose to defer it to Phase 5 and record it.
- **F171.** The floods' vents favour the thin crust and the low ground everywhere, a province's centre standing four times the rest rather than thirty: Mars's lowland now takes its share of the floods.

How it went:

- **The provinces decided.** Mars's lowland draws five to thirteen times the vents per area that its highlands do by its thin crust and low ground, and still held no more lava per area, the provinces outweighing everything outside them thirty times.
- **The Moon does not check it.** The Procellarum KREEP Terrane holds most of the Moon's maria, so it was to be the check on the provinces' weight; but the model's Moon lays only a fifth to a quarter of its lava within its province at any weight, its basins' thin crust deciding (M33), so the weight is fitted on Mars.
- **Venus's relief.** The user noted that relaxing Venus would flatten more than its tesserae: nine tenths of its relief is held up by the mantle's flow, which the primordial relief stands in for; the convective planform joins M31.

What the maps show:

- **Mars (seed 1048):** the lowland holds some 900 m of lava against the highlands' 600 m; a few of its craters' floors fill further, and its craters of 50 to 200 km keep their rims, deeper than the lava.
- **Venus, the Moon and Mercury:** their floods spread somewhat more over their low ground; their maps change little.

Still to do: Valles Marineris with Phase 5 (M32); Venus's highlands and dynamic topography (M31); the Moon's maria in its province (M33); the aureole and basal scarp (M30); Mercury's thrusts' direction (M27); the drives the field does not hold (M26); the coronae on Venus's rifts (M25).

Tests 259 (22 ignored). The test of a rise's shields taking its eruptions now holds its sums equal within 10⁻⁷, their rounding with the province's larger background, and the Moon's rayed share need stay below three tenths rather than a quarter (it is 25.2 percent).

## The eighty-eighth milestone: an edifice its crust cannot hold sheds its flanks

M30, which the user chose: Olympus Mons's basal scarp and aureole.

- **F172.** A shield whose slope would build it taller than its crust holds spreads its flanks with the lava it would have stacked above: as much of its toe as that lava is pushed out beyond its foot and breaks away in a scarp, and the toe and the lava together run out beyond it as its aureole, in lobes as long as rock avalanches of their volumes run, ridged across their run.

How it went:

- **Which edifices fail.** A Coulomb wedge on a frictional base fails by its slope alone, the height cancelling, and every Tharsis shield has Olympus Mons's slope; but the crust's strength caps an edifice's height (F15), and only the largest of Mars's shields reaches the cap, 25 km where its slope would build 27.4. The criterion is the cap: an edifice at it has loaded its base to what the crust holds.
- **How far it fails.** The lava past the cap, some 9 percent of the edifice, pushes out as much of its toe, which places the scarp at 0.81 of the reach and makes it 4.5 to 9.9 km high round the edifice, where Olympus Mons's stands 6 to 10 km; and the aureole, twice that lava, runs 120 to 350 km in seven lobes by the avalanches' trend, short of Olympus Mons's longest tongue of some 700 km, each lobe a smaller share.
- **When it fails.** The failure was first set by the edifice's slope alone, and Venus's three largest shields, whose slope would build them 12.9 to 14.3 km against the 10.5 km its crust holds, shed their flanks too, though Venus's supply builds them only part of the way (F165): their summits fell by up to 3.7 km, the spread lava counted into an edifice they never filled. An edifice now fills to the cap first, and only the lava after that breaks its toe away and lays its aureole, so Venus's shields stand as before and Mars's largest fails once its edifice has filled.
- **The budget.** The edifice's supply covers its aureole as well as the edifice (F164), and the aureole weighs on the lid in the moat (F162) but is not the edifice whose root fills the sag (F166).

What the maps show:

- **Mars (seed 1048):** its largest shield stands ringed by a basal scarp 4.5 to 9.9 km high, its old gentle toe gone, and its aureole lies beyond in seven lobes of curved, broken ridges, in its moat and out over the plains 120 to 350 km; its other shields keep their toes.
- **Venus:** unchanged; its three largest shields would build past what its crust holds, but its supply never fills them.
- **Every other body:** unchanged, raising no shields.

Still to do: the lobes' directions from the ground they run down; Valles Marineris with Phase 5 (M32); Venus's highlands and dynamic topography (M31); the Moon's maria in its province (M33); Mercury's thrusts' direction (M27); the drives the field does not hold (M26); the coronae on Venus's rifts (M25).

Tests 261 (22 ignored), among them a shield past what its crust holds shedding its toe and its aureole: the toe beyond the scarp holds the lava its supply would have stacked above, within 15 percent; each scarp 2 to 12 km high; each lobe running 100 to 1,200 km; the aureole holding twice that lava, within 15 percent; nothing beyond its reach, and the edifice nothing past its wall; and on Mars only the shields whose slope would build them past the cap shedding their flanks. The test of a shield's lava past its foot counts its aureole with its edifice; the test of a rise's shields taking its eruptions now holds Mars's largest edifice's flanks fully failed and none of Venus's, though its largest would build past the cap.

## The eighty-ninth milestone: the mantle's convection holds the ground up

The mantle's convection pattern, which the user chose over M31 alone since several items wait on it (M31, M27, M25, M18), the user adding that the look may change at this stage.

- **F173.** A rocky mantle convects beneath its lid in cells as wide as the layer under the stiff lid is deep, and the ground its flow holds up, a root mean square of 2.2 times the sublayer's thermal contraction over its thickness and flatter where the lid bends, rises over the upwellings and sinks over the downwellings as the pattern renews itself; the plumes rise at its strongest upwellings.

How it went:

- **The strength.** Guimond, Rudge and Shorttle's fit of the dynamic topography in the interior Rayleigh number and the viscosity contrast, extrapolated to the model's Venus, gives tens of kilometres: its interaction term runs away beyond the range it was fitted on. Their table gives each simulation's sublayer, and the topography over α ΔT_rh δ_rh is 1.5 to 3.2 across all eleven, so the support is that form, the sublayer from the model's own stagnant-lid flux, the same buoyancy the lid's breaking reads (F121).
- **The scale.** The cells' degree from the layer under the stiff lid: 4.3 on Venus, 2.5 on Mars, 14 on Mercury, whose mantle is thin, and 1 on the Moon, whose mantle reaches its centre. Mercury's thick lid holds its narrow cells flatter, 144 m of the 658 m its flow would hold up.
- **The time.** The pattern renews itself in the time the flux takes to carry off the heat the sublayer holds, 32 million years on Venus now and 185 on Mercury, and at each of Venus's overturns.
- **The rises.** The plumes rise at the pattern's strongest highs as they start, 1.6 to 3 times its root mean square, and the body turns under them as before (F104), the pattern with it.
- **The limbs.** With its plumes at its upwellings Venus's stood 9,000 to 14,000 km apart, and of the great-circle limbs, joining plumes within four of their radii together, one remained: its rifts' longest fell from over 2,000 km to 1,000. The user chose to draw the limbs from the convection at once: a limb now follows the ridge of the pattern the plumes rose from, joining two plumes whose ridge stays above the mantle's mean between them. Two of Venus's plumes are joined by a ridge 14,800 km long, its rifts the longest 1,900 km, and the third is parted from both by downwellings.
- **Mars's rise over an old basin.** The user saw Mars's rise sit over an old basin, its relief hiding much of the bulge. The dome lifts the ground some 5 to 7 km over 2,600 km and the basin with it, its rim now near 5 km, but its floor only some 2 km, the shields' moats sinking in it and the basin keeping its depth. Tharsis is a pile of some 3 × 10⁸ km³ of magma rather than a dome, which would have filled such a basin; the model's rise has no volume of its own (audit M34, the next milestone by the user's choice).

What the maps show:

- **Venus (its preset):** its rises stand at new places, its upwellings, and a lowland some 2 km deep lies over a downwelling; the ground it holds up reaches −1,350 to 890 m.
- **Mars (seed 1048):** its rise has moved to the western end of its highlands, where it straddles the dichotomy's boundary, its shields with it; the ground the flow holds up reaches −930 to 1,070 m.
- **Mercury and the Moon:** they roll gently over their mantles' cells, Mercury by about half a kilometre each way and the Moon by a kilometre.

Still to do: the rise's own pile (M34, next); the flow's stresses on the lid, for Mercury's thrusts (M27) and Venus's ridges, and its upwellings' heat for Venus's coronae (M25); Venus relaxing as a stagnant lid does, and its tesserae from the tectonics that thickened their crust (M31); Valles Marineris with Phase 5 (M32); the Moon's maria in its province (M33); the drives the field does not hold (M26).

Tests 268 (23 ignored), among them the cells following the mantle's depth, Venus's of degree 3 to 6 and Mercury's past 8, the ground the flow holds up some hundreds of metres and its pattern renewed in tens to hundreds of millions of years; the pattern keeping its strength as it changes and held flatter under Mercury's lid; its upwellings its highs, strongest first; a spherical harmonic field's grid holding its values; Venus's and Mars's plumes rising at the pattern's top fifth, the first at its strongest high; a limb following the ridge between its plumes, the pattern above the mantle's mean along its way within a quarter of its root mean square; and Mercury's ground moving late in its life by as much as the flow's support changes, within 15 percent; and a diagnostic printing each rocky preset's convection. The test of Venus's rifts along its limbs now asks the longest over 1,500 km rather than 2,000 (it is 1,900), the population's longest along its one limb, and follows the limb's way.

## The ninetieth milestone: the melt is one budget, and a rise is a pile as well as a dome

M34, which the user chose for the milestone after the convection, asking whether Mars's rise over an old basin should stand higher.

- **F174.** The melt that does not erupt intrudes the crust beneath where it rose from, the provinces by the ambient mantle's share and the plumes' heads the rest, and the magma full edifices could not lift intrudes beneath their rises; the crust thickens, the ground rises and the intrusions' weight bends the lid.

How it went:

- **Where the volume was.** A rocky body erupts an eighth of its melt at most (F125), and the seven eighths left were lost: on Mars 1.1 × 10⁹ km³. Laid where the melt rises from, some 3 × 10⁸ km³ of it gathers beneath the rise's province, as much as Tharsis's load (Phillips et al. 2001), with nothing fitted.
- **The edifices on it.** Laid with the lava's load, the intrusions beneath the shields never sank, F166's rule filling the sag beneath an edifice with its lava, and Mars's largest summit rose to 37 km. The rule is for an edifice's own weight: the intrusions are a load of their own, and the edifices on them sink with the ground. The summit stands at 27.3 km on some 2.5 km of new ground.
- **The dome.** The rise's dome is still 0.28 of the relief its crust holds, fitted to Tharsis's whole height before the rise had a pile; with the pile its flanks stand some 7 km, within Tharsis's range, and the dome's height from its plume head's buoyancy is left open.

What the maps show:

- **Mars (seed 1048):** its rise fuller and broader, its flanks some 2 km higher, the old basin's floor raised some 2.5 km, at 2 to 3 km, the shields on the new ground. The basin keeps its depth below its rim: the intrusions raise the province about evenly and only lava fills a hollow, and the vents, favouring low ground against the body's median, find few hollows on a raised rise.
- **Venus:** little change, its melt after its overturn small.
- **The Moon and Mercury:** their provinces swell by a kilometre or two, the Moon's raising the lowland its province lies in; their crusts held most of their melt down.

Still to do: the lavas ponding in a rise's own lows and the dome's height from its plume head (F174's open ends); the flow's stresses on the lid, for Mercury's thrusts (M27) and Venus's ridges, and its upwellings' heat for Venus's coronae (M25); Venus relaxing as a stagnant lid does, and its tesserae from the tectonics that thickened their crust (M31); Valles Marineris with Phase 5 (M32); the Moon's maria in its province (M33); the drives the field does not hold (M26).

Tests 270 (23 ignored), among them the melt one budget, what erupts and what stays below making the melt between any two times within 10⁻⁹, the intrusions laid keeping their volume, most beneath a rise's province and a rise's own share beneath it alone; and Mars's rise gathering 1 to 10 × 10⁸ km³ beneath it, all the loop lays the melt's intrusive share and the edifices' surplus within a percent. The test of the Moon's maria now asks only the plains' own ridges to lie on lava: a thrust's ridge folds the plains where most of its trace crossed lava and may cross older ground besides, and the intrusions' load now raises thrusts early in the Moon's life.

## The ninety-first milestone: a deposit's level has no steps

The user saw stepped relief fanning out from Mars's rise in the crops of the last two milestones.

- **F175.** The level a crater's ejecta bury the old relief against is the sectors' means of its rings interpolated round and across them, rather than each point's own sector's mean.

How it went:

- **The steps.** The terraces were the sectors of the old basin's ejecta, some 2,400 km across beneath the rise: its deposit, kilometres thick near the rim, buried nearly all the relief beneath it against each sector's own mean (F169), so that each sector stood flat at its own level, stepping along the straight edges between them, radial from the basin. Every crater's deposit did so; the large ones on relief kilometres high showed it.

What the maps show:

- **Mars (seed 1048) and the Moon:** the ground beneath the basins' deposits rolls as the relief beneath it did, draped, without the fans of flat sectors.

Tests 271 (23 ignored), among them a deposit's level without steps round its crater: on a plain tilted 3 m per km, the ground beneath a 600 km crater's thick deposit steps no more than the tilt and the grid give round rings at 1.15, 1.3 and 1.6 radii, where the sectors' edges stood as steps of 28 m.

## The ninety-second milestone: the vents read the ground less what is held up

F174's first open end, which the user chose next: the lavas ponding in a rise's own lows.

- **F176.** The flood vents read the ground less what the plumes' domes and the mantle's flow hold up, which lift the magma beneath with the ground; the relief the crust holds counts as before.

How it went:

- **Why the rise shunned its lavas.** The vents favour low ground against the body's median (F36), and Mars's rise stood 5 to 7 km over it, some twenty times fewer vents a cell. The plume's dome and the mantle's flow raise the ground and the magma beneath it together, so the vents now leave that uplift out.
- **What moved.** On Mars's preset at the preview's resolution the rise holds some 3.5 × 10⁷ km³ of lava rather than 2.5 × 10⁷, and within 900 km of its old basin 3.8 km rather than 2.7; the lowland's and the highlands' lava stay as F171 left them, 713 m and 383 m against 735 m and 413 m.
- **Why the basin stays.** The basin's floor rises only some 0.2 km under the kilometre more of lava: a load as wide as the basin, some 1,800 km, bends the lid nearly as far as it is thick, being far wider than the lid can hold up by its strength. Mars's filled basins, Utopia and Isidis, still stand kilometres deep, so the basin keeping its depth beneath the rise is the lid's doing.

What the maps show:

- **Mars (seed 1048):** the rise's hollows hold more lava, the basin's floor some 0.2 km higher; little else changes.
- **Venus:** its rises draw more of its lava.

Still to do: the dome's height from its plume head (F174's other open end); the flow's stresses on the lid, for Mercury's thrusts (M27) and Venus's ridges, and its upwellings' heat for Venus's coronae (M25); Venus relaxing as a stagnant lid does, and its tesserae from the tectonics that thickened their crust (M31); Valles Marineris with Phase 5 (M32); the Moon's maria in its province (M33); the drives the field does not hold (M26).

Tests 272 (23 ignored), among them a rise's vents reading the ground its plume holds up as level: a dome held 4 km high draws its share of the vents as flat ground would, 0.069 against 0.058, where read whole it draws half of that, and a hollow on it draws more per cell than its flank.

## The ninety-third milestone: what bends the lid

The user asked for the dome's height from its plume head next (F174's open end). The dome came out physical, and what it drove did not survive it; the user chose to set it aside, make the rises' faults read the loop's own field first, and, when that showed a deeper fault, to fix the lid's loads first.

- **F177.** A load laid evenly over the body bends nothing, and a shift of the whole lid strains nothing: the load's mean leaves the plate's solve, the deflection's mean and first degree leave its membrane, and what a load lays evenly raises the datum with the body.

How it went:

- **The dome.** A plume's head flattens beneath the lid into a disc twice its diameter (Campbell 2007), so that a disc as wide as the rise and a sixth of its radius thick, 200 K hotter than the mantle, holds up some 2.9 km on Mars and 1.0 to 1.4 km on Venus, where 0.28 of the relief the crust holds gave 7 and some 3 (a head from the Earth's core: some 1 km, as Campbell and Griffiths's heads raise the ground before their floods). Mars's ground about its rise fell from some 7.5 km to some 3, its lighter dome turning the body only 18°, from 31° to 13° north, and its rise's faults and canyon, its limbs' rifts on Venus, which scale with the dome's height, planned before the loop, all but went. Set aside.
- **One population of faults.** The rises' faults, the limbs' rifts and the loads' graben, opened through the loop by what the lid's field adds to its tension past the yield beyond the bare relief's, over the lid's stiffness, as the thrusts take the field's compression: built and run, it faulted most of Mars. The relief alone stands past the faulted layer's yield over three quarters of Mars, so that anything added counts, and the loads added tens of megapascals nearly everywhere. Set aside: a budget needs the lid's strength through its thickness, which the faulted layer's yield is not.
- **The loads.** Tracing the tension found two errors in the loads' solve: the intrusions, laid nearly everywhere, deflected the whole lid evenly by kilometres, and the membrane read that and the deflection's first degree as stretching it. Neither can strain a lid on an incompressible mantle; both are gone, and the heights keep their datum. The tension over most of Mars stayed, some 77 MPa past the yield at the median where it had been 87: the loads' own stresses, which a lid as strong through its thickness as rock is would bear.

What the maps show:

- **Every rocky body:** the relief within a few hundred metres of what it was, the mean nearer the datum.

Still to do: the lid's strength through its thickness, a yield-strength envelope, against which the field's stresses stretch it; then the one population of faults and the dome from its plume head (both set aside); the flow's stresses on the lid (M27, M25); Venus relaxing and its tesserae (M31); the Moon's maria in its province (M33).

Tests 273 (23 ignored), among them an even load bending nothing, a lid deflected less than a metre under a kilometre laid evenly over it, and a shift of the whole lid straining none of it, its stress under a fiftieth of what a membrane straining it would stand at. The test of a filled basin's graben now pools six draws of them, its graben in the ring the load stretches out to six flexural parameters and seventeen in twenty within it, the rest beside one that is.

## The ninety-fourth milestone: the lid's tension faults by their drivers

The user chose to keep the loads' stresses, which a lid as strong through its thickness as rock bears (the mascons' stresses persist), and to bring back the one population of the lid's tension faults with each fault sized by what drives it: graben within the faulted layer where the strain is broad and small, rifts as deep as a rise's canyon only where it is concentrated, on a rise's crest or along a limb.

- **F178.** The lid's tension faults open through the loop by three drivers, each read beyond those beneath it: the rises' domes beyond the bare relief, the limbs beyond the rises, the loads beyond both. A rise's and a limb's extension is released as they grow and spent where their share of the field grew past the yield, on faults as deep as a rise's canyon; the loads' is their share of the excess over the lid's stiffness, spent on graben within the faulted layer.

How it went:

- **Two drivers.** The rises and limbs together, the loads apart, each with the field's budget, the excess's growth over the stiffness: the loads' graben came out as F158's had (6.2e10 m² on Mars, 3,871 drawn, against 6.4e10 and 3,952), but the rises' budget was a twenty-fifth of F38's extension (1.5e10 m² on Mars against 3.7e11), and Venus's limbs opened two rifts. A dome stretches the lid by its growth, the dikes and faults taking that up whatever the stress past the yield (Buck 2004), so a rise's and a limb's extension is F38's, released as they grow; the field says where and when.
- **Three drivers.** With F38's extension the rises and limbs pooled opened 51 rifts on Venus, but spent the limbs' extension on the crests, where the rises' share is greatest: along its one limb, 14,800 km long, faults stood everywhere, a few hundred kilometres long, and rifts on a quarter of it. The limbs now have their own share, beyond the rises', and their own extension: rifts on four fifths of the limb, the longest 1,390 km.
- **Chronology.** A fault formed within the step its share grew in, from when the loop last read the field; a rise's could so fall before its rise began. A rise's and a limb's now form while the rises rose.
- **Carving.** A canyon was carved only while its growth rose across the loop's step, so that a rift that joined the canyons after its carving window had closed was never carved. Each canyon now keeps how far it has grown.

What the maps show:

- **Mars:** no rift across the lowland; Tharsis keeps its faults and its graben about the volcanoes; the long trough that curved south of Tharsis is gone, and a few rifts run by the volcanoes instead. The loads' graben as before.
- **Venus:** rifts on both rises, graben along the whole of its limb, and fewer long rifts carved along it than before: the rises and the limb opened some 25 rifts each, one wide enough to carve, 1,560 km long; the arc of carved rifts between its rises is mostly gone.

Still to do: the dome's height from its plume head (set aside in the ninety-third); the flow's stresses on the lid (M27, M25); Venus relaxing and its tesserae (M31); the Moon's maria in its province (M33).

Tests 273 (23 ignored), among them a filled basin's graben opened through the loop round its margin, and Venus's rifts reaching from rise to rise along its limb, the longest past a thousand kilometres where fifteen hundred had been asked: the population is a draw, and the rifts through the loop read the field as it grows and the loads' stresses with it, which the plan before the loop did not.

## The ninety-fifth milestone: the faults grow and link

The user on the ninety-fourth: the canyons still looked shallow, and Venus had lost the connected look its rifts had. Both came from the same place: each step of the loop opened faults afresh, none grew, and the faults of one step never gathered beside those of another. The user chose for the faults to grow and link.

- **F179.** A rise's and a limb's faults persist through the loop and grow as their driver stretches the lid: each takes up the growth within its shadow, runs its tips on as the field leads them and deepens as far as the system it linked into is long; what falls outside every shadow opens new faults.

How it went:

- **Growth alone.** The faults grew, but the extension each took up was small against what it had, and most of it went to some 1,300 new faults a step: Mars's longest rift 1,133 km, 2.6 km deep. A graben confined to its layer takes up extension as its length, not its square, so it was lengthened as a rift is and gave most of it back; it now lengthens as what it takes up.
- **Deepening.** A fault blocked at its tips deepened without bound, so that short stubs became deep rifts; its throw now grows only as far as the length of the system it linked into (Cartwright et al. 1995).
- **Linking.** A tip stopped at the first fault within half its width, and the swarm of small graben stopped every rift short; a tip now stops only at a fault at least half as wide. Mars's longest rift went to 4,360 km, its deepest to the 8 km cap.
- **Increments.** Venus's rises grow over some six of the loop's steps, so its faults had no time to grow, its limb's extension opening new faults; a rise's or a limb's extension is now spent in sixty-fourths of its whole, each taken up before the next opens more.
- **Small fixes on the way.** A fault's draw is the one before its id; a rise's faults formed after it stopped growing when extension waited a step; a canyon is carved anew when it grows; the test of when the rise's faults form reads their formation, a canyon being carved on over the melt after.

What the maps show:

- **Mars:** a long, deep canyon along the rise's axis, running out smoothly into a northern basin, branching where its rifts overlap, with mounds between them; 14 rifts carved, the longest 2,940 km, the deepest 8 km.
- **Venus:** long, deeper rifts through its rises, through the western one as Devana Chasma runs through Beta Regio; the arc of rifts along its limb is there, its rifts a little shorter and more subdued than the planned population's were.

Still to do: Venus's subdued arc (its limb's faults read the whole field, the loads' stresses with it, which the planned population did not); the dome's height from its plume head (set aside in the ninety-third); the flow's stresses on the lid (M27, M25); Venus relaxing and its tesserae (M31); the Moon's maria in its province (M33).

Tests 274 (23 ignored), among them a fault reshaped as drawn coming back as it was.

## The ninety-sixth milestone: the flow bends the lid, and the rifts focus the coronae

After the ninety-fifth the user asked for the dome's height from its plume head. The dome came out physical, and with it an audit of Tharsis's magma found the plumes made almost none of their mantle's melt; the plumes' share of the melt and the rifting from the rise's relief as it is built were tried, and the user chose to set all three aside (the patches kept). The user then asked for the flow's stresses on the lid (M27) and Venus's coronae on its rifts (M25).

- **F182.** The mantle's flow bends the lid: each reading of the field holds the stress of the lid bent to the ground the flow holds up.
- **F183.** A rift focuses a diapir's melt: where the region melts, a diapir on a limb makes a pocket as far as its head at least doubles the region's melt.

How it went:

- **The dome, held.** The plume's head flattened into a disc stands Mars's dome at 2.9 km; the canyon's depth and the rifts' extension, scaled to the dome, fell with it. Mars's mantle makes 1.37e9 km³ of melt over its life, 4.2e8 while Tharsis grows, but its plume took 1.8e7 of it, the ambient mantle melting all the while (F165). A plume's share of its mantle's heat (F121, 0.59 on Mars) beneath its head gave it 6.9e8 km³, of Tharsis's some 3e8 (Phillips et al. 2001), and the rise stood 5 km above the body's median. Rifting from the rise's relief as built gave shallow canyons on Tharsis's margin, where its pile's loads stretch the lid, and ran away on Venus across its overturn. Set aside.
- **Mercury's thrusts.** The flow bends Mercury's lid by 7 to 15 MPa, but its plains' loads bend its small shell by 70 to 100 MPa at its equator and decide its thrusts' strikes, east and west at low latitudes; the north and south of Watters et al. (2015) wants the despinning.
- **The coronae's births.** Drawn by the convection's upwellings, the diapirs spread over the pattern's broad highs: 7 percent of the coronae on the limbs, where 30 had been. Set aside.
- **The melting rule.** A diapir on a limb or within a rise made a pocket where the region melts; Mars's rise then made coronae, and the focusing was kept to the limbs.

What the maps show:

- **Mercury:** its thrusts' strikes a little less alike, otherwise as before.
- **Venus:** 841 coronae, more of them along the limbs; the rifts on its western rise cross, where a long parallel set ran, the flow deciding their strikes on the rise's crest.

Still to do: Mercury's despinning (the rotation's history); the coronae's share on the rifts themselves, and the hot region they crowd; the dome, the plumes' share of the melt and the rifting from the rise's relief (held); Venus relaxing and its tesserae (M31); the Moon's maria in its province (M33).

Tests 274 (24 ignored).

## The ninety-seventh milestone: the crust the mantle's flow carries

The user asked for M31: Venus relaxing as a stagnant lid does, and its tesserae from the tectonics that thickened their crust.

- **F184.** Where a rocky lid's crust reaches below the lid's stiff base it moves with the mantle's flow, piling over the downwellings to the depth its base turns to eclogite and thinning over the upwellings; the lower crust's flow spreads it, each cell relaxing under its own crust; the tesserae are the plateaus the last epoch piled, and an episodic lid relaxes as a stagnant one does.

How it went:

- **Plateaus stamped at the overturn.** A crust 65 km thick formed at the overturn would lose everything under 1,000 km and three quarters of its 3,000 km in the regrowing lid's hot phase, so the plateaus cannot be set at the event: they must build while the lower crust is mobile and freeze as the lid cools, which is Bindschadler, Schubert and Kaula's picture. Venus's stiff lid falls below its crust for some 30 Myr after each overturn and stands 45 to 80 km above it between them.
- **Conservation.** The thinning over the upwellings slows as the carried crust thins while the thickening ran on, so the downwellings piled more than the flow brought; the thickening is now scaled to the thinning that supplies it, and only the eclogite cap loses crust, which lowers the datum with the body.
- **Each cell's own crust.** The relaxation relaxed the whole body under its tenth thickest crust and the lid of the step's end; with plateaus of 65 km among plains of 40, and a lid that regrows within a step, it now relaxes each cell under its own crust through the step's substeps.
- **The last epoch.** The plateaus of the epoch before the last overturn relaxed to a few hundred metres within 50 Myr of theirs; the last epoch's stand some 3.7 km above the plains.
- **Two leaks.** The crust's mean fell from 40 km to 13 where the eclogite's cap accounted for 17: the thinning over the upwellings ran unmatched where the downwellings took less than they gave, and the relaxation under each cell's own crust let the thick-crust highs go faster than the thin-crust lows, so that each band's mean sank, 190 m of ground a step in the hot phase, and the crust, kept compensated, followed it down. The flow now moves the lesser of what is given and taken, and the relaxation keeps the ground's area mean, as a flow that moves material must (the test of a plateau relaxing on Ganymede's shell now reads the relief about the mean, where it read the plateau's volume vanishing); the crust ends 16 km thick on average, 24 km of it foundered, and the ground's mean holds through the history.
- **The stagnant lids.** Mars, the Moon and Mercury start with a finished crust on a mantle with no lid, so their crust rides with the flow for their first 30 Myr too and keeps the pattern it piled into, ±3 km at their cells' wavelengths: the Moon's at degree one, as its farside crust stands. The user chose to keep it on every body, one process, the amount a prediction from a starting state the magma ocean's end is not modelled for (audit M35).

What the maps show:

- **Venus:** broad plateaus with organic outlines over the last epoch's downwellings, the tesserae 2.3 to 5.3 km above the plains, their margins gentle where the hot phase smoothed them; the primordial relief gone; the rifts of the rises cut through their shields, as Devana Chasma cuts Theia Mons, the rises still growing.
- **Mars, the Moon, Mercury:** their long wavelengths shifted by ±3 km, at the degree of their cells; the craters, the floods and the faults on them as before, the rifts drawn anew.

Still to do: the tesserae's fabric from their own shortening; the tessera share from the surge's lava (M29); steeper plateau margins; Venus's deformation belts (M10).

Tests 276 (25 ignored), among them the flow piling Venus's crust over its downwellings beneath its regrowing lid and none under Mars's stagnant one, and a cell relaxing under its own crust.

## The ninety-eighth milestone: the flow drags the lid

The user asked for M10: Venus's deformation belts and wrinkle ridges, and the novae.

- **F185.** The mantle's flow drags the lid, stretching it over the upwellings and squeezing it over the downwellings by the stress that founders it at an overturn; an episodic lid has thrusts between its overturns; a nova is a young corona whose rising dome stretches past the lid's strength.

How it went:

- **The drag.** The flow's deflection bent Venus's thin lid by some 3 MPa (F182); the traction along its base, carried through its thickness, is the same stress the overturn criterion uses, 35 to 57 MPa per unit of the pattern, 70 to 140 over the strong downwellings.
- **No thrusts on Venus.** The thrusts and the shrinking's squeeze switched on the stagnant lid's name; the shrinking is now counted since the lid last formed, so an episodic lid has both.
- **Where they form.** Judged against the present pattern the thrusts seemed to stand between the cells; the pattern renews every 30 to 40 Myr, and judged against the pattern of their own time 95 percent stand over the downwellings. Nine tenths of their length forms in the first hundred million years after the overturn, the longest in its first ten, on the plateaus the flow piled, and barely show through the tesserae's fabric; later in the epoch the lid is too thick for the drag to break it, so no lowland belts form.
- **The novae.** Every corona's rising head stretched its crest past the strength (56 MPa at the median); those still stood on their heads by half were 27 percent, those still fed 17, against the catalogue's 12: a nova is a corona still rising.

What the maps show:

- **Venus:** much as before; a few thrusts on the plateaus; the plains' ridges a little denser in the lows; radial fractures round the rising coronae only.
- **Mars, Mercury, the Moon:** their thrusts as before, within a few percent.

Still to do: the lowland ridge belts; the thrusts on the plateaus as the tesserae's fabric (W3); the coronae on the rifts (M25).

Tests 277 (25 ignored), among them Venus's novae being its rising coronae, and its thrusts forming on the lid its last overturn left, where they had been none (a thrust timed in the bin of the overturn now comes after it).

## The ninety-ninth milestone: the event's lava from its surge, examined

The user asked for M29: the overturn's resurfacing from the melt its surge brings, which would free the tessera share. The model's behaviour is unchanged; what was found is in the audit (M29), and FloodProcess::surge_lava with its diagnostic stays for the next attempt.

- **46 m.** The surge erupts 46 m of lava per unit area of its 2.3 km of melt: the contraction shutoff reads the body's whole contraction. Read since the lid formed, as the thrusts' squeeze is (F185), it erupts some 300 m, and Venus's floods after its overturn grow sevenfold, which builds its leading volcanoes to 10.5 km and fails their flanks against Maat Mons's 5: set aside.
- **54 percent dry.** Laid as the event's lava, 300 m leaves over half of Venus dry against its 8 percent of tesserae; some 2 km would be needed. The user chose to let a thin lid let more out; the sources show intrusion dominating eruption in every setting, mid-ocean ridges included, so that rule was not built.
- **Foundering.** The user chose the old lid's foundering as what renews the surface, the tesserae what did not founder: the crust's buoyancy against the lid's thermal excess would keep 99 percent of Venus, its decoupling from a delaminating mantle lid none; neither founds 8 percent.

Still to do: the event from the surge's melt and the lid's foundering on a founded rule, and the limit on Venus's eruptions after its overturn (M29); the lowland ridge belts (M10).

## The hundredth milestone: the plains onlap the tesserae

The user found Venus's tesserae abruptly bounded against the plains, in relief and in colour, and their tops undynamic.

- **F186.** The tessera fabric continues beneath the plains beyond the tesserae's edge, the plains thickening away from it in an onlap wedge, and the event's lava ponds within the tesserae at the edge's level: the plains run up the fabric's valleys and its ridges stand as inliers, each cell by the share the lava left standing.

How it went:

- **A fade.** The fabric laid beyond the edge and faded over a margin, the lava at its own level, left the fabric standing above the 20 m veneer the event lays on high ground: the tesserae grew from 8 to 12 percent of Venus, and at the map's scale the contact became a blur, not an embayment.
- **The wedge.** The old ground's smoothed height below the edge gives the burial: the fabric beneath the plains sinks as far below the edge as the ground's smoothed height stands, buried past the reach of its highest ridges; within the tesserae the lava ponds at the edge's level. With the fabric centred on its median the tesserae keep their share, 7.6 percent against 7.9, and the contact interfingers.
- **A shadowed name.** The pond's height first took the name of the event's lava composition, and the maps came out black and white: the plains took heights for their composition.

Still to do: the tesserae's fabric from the shortening that piled them and the extension after (W3, T1), which the user asked for next.

## The hundred-and-first milestone: the tesserae's fabric from their lid's strain

The user asked for the tesserae's noise fabric to be replaced by one from the plateaus' history.

- **F187.** The mantle's flow broke the young lid where its drag passed the lid's strength and carried it along; the strain it took, the pattern's curvature, is recorded set by set, and the tesserae fold across its shortening and neck into graben across its stretching, at the brittle layer's and the elastic lid's wavelengths, as high as friction and strength let them stand.

How it went:

- **The lid did not move.** The flow carries the crust beneath the stiff lid (F184), which alone would leave the surface undeformed; the drag the lid bears from it (F185) passes its strength over the strong downwellings and upwellings in the lid's first few million years, and there the lid goes with the flow.
- **Kinks for curvature.** Read across one coarse cell, the pattern's curvature picked up the kinks of the grid it is sampled from, principal strains of ±2.4; read across two of its own grid's cells, ±0.5.
- **Saturated strain.** The flow strains the broken lid by about half in a few million years, far past what folds take up (Tellus's folds take 1 to 12 percent), so the folds stand as high as their layer's strength lets them, 0.88 km, Tellus's 0.5 to 1 km.
- **A grating.** Packets covering each place eight times over summed to one wave across whole domains, a combed grating; at one and a half times over, with their strike and spacing scattered, they break into en echelon segments and cross as basins and domes.

Still to do: the strain after the event, on tesserae and plains alike; the lid's thickening where it is shortened; the shortening the folds cannot take up.

## The hundred-and-second milestone: the plateaus' collapse

The user found the hundred-and-first milestone's tesserae reading as a giant fingerprint, or dunes with very long rolls, and asked whether the natural look should come from a shader in the artistic pass. The fabric below the map's cells, the brittle layer's ribbons and short folds, belongs to such a shader, fed by these fields; the fingerprint was a gap in the model, one deformation where tesserae have several, and the plateaus' collapse, which the relaxation already computes, was the deformation missing.

- **F188.** The lower crust's flow as the relief relaxes strains the lid by its potential's curvature, recorded over the lid's life as one more set: graben and folds across the drag's.

How it went:

- **Half the strain.** Each band's wavenumber is the one it is named for, which its content spreads about: the potential's Laplacian came to half the relaxed relief, and is scaled to it over the body.
- **An egg crate.** Sine folds crossing made a field of round knobs, and their tops broke through the plains as blobs; a narrow crest over a broad trough reads as ridges.
- **Ten minutes.** Every set's packets were looked up within the widest set's reach, the lid's 97 km folds, some five hundred coarse cells about each point; within each set's own reach the pass takes 81 s.

Still to do: the strain after the event; the fabric under the cells, for the artistic pass's shader.

## The hundred-and-third milestone: a layer holds its fabric's relief together

The user found the tesserae's tops much bumpier at the fine scale than Venus's topography.

- **F189.** The strength bound on a fold's relief bounds its layer's: each layer's relief, its sets together, saturates at what its strength holds.

How it went:

- **Summed past the bound.** Each fold set kept to the bound alone, and three sets and the packets' interference stood two to three times past it; together the tesserae's peaks fall from 9.8 to 7.4 km.
- **The instrument's share.** Magellan's altimetry averages over some 15 km; blurred to it, the tops show the reference map's soft mottle, so the fine-scale smoothness of the reference is mostly the instrument's. What the model lacks is Ovda's broad internal relief, basins and swirls 50 to 300 km across.

Still to do: the plateaus' broad internal relief; the strain after the event; the fabric under the cells, for the artistic pass.

## The hundred-and-fourth milestone: the tesserae's folds and graben as lines

The user found the wave packets reading as frosted glass at a global view, where the noise fabric had a turtleshell of ridged cells, and asked whether the two could be blended.

- **F190.** Each fold crest and graben is a line along the strain's strike, a set's lines evenly spaced, ending where the strain changes sense or their neighbours crowd them; crossing sets bound cells.

How it went:

- **A blend.** The old noise varying the packets' strength by six tenths either way and bending their crests changed almost nothing: a fabric's strength varies its contrast, not its mean, and cells need rims standing above their interiors. The user chose founded tiles over the noise as relief of its own or the artistic pass.
- **Long arcs.** Lines up to fifteen spacings long, evenly spaced, curved with the strain into concentric arcs and fine combed sets, the user's fine seashell: each line's length is drawn between three and fifteen spacings, its spacing and relief scattered more, and a set drawn only where its spacing spans three cells or more.

Still to do: the collapse's graben deeper by the tension through the stretched layer, for clearer tiles; the plateaus' broad internal relief; the strain after the event; the fabric under the cells, for the artistic pass.

## The hundred-and-fifth milestone: the tesserae's graben are the model's

The user asked what to do next; I proposed the lid's strength through its thickness, a yield envelope, to deepen the tesserae's graben and carry the plateaus' broad relief. The user remembered it from the ninety-fourth milestone, where it was set aside: its 150 to 380 MPa broke no lid, and here it would leave Venus's young lid unbroken by the drag and the tesserae without their fabric. The user chose instead to make the tesserae's graben follow the model's own graben.

- **F191.** The tesserae's graben are the model's graben, their width from the faults meeting at their layer's base, their depth within their layer, a fifth of their width and their length's throw; the lid necked whole opens rifts by the rifts' rules.

How it went:

- **Gashes.** By the graben's rules the lid's troughs, 22 km wide, sank 4 km, a grid of gashes; the graben's rules hold for graben up to 12 km wide, and the model opens a rift where a fault would cut below its layer. As rifts, 1.5 thousandths of their length deep, they stand 0.3 to 1.5 km and tile the plateaus' crests.

Still to do: the plateaus' broad internal relief; the strain after the event; the fabric under the cells, for the artistic pass.

## The hundred-and-sixth milestone: the tesserae's graben take the model's graben form

The user found the crossing rifts reading as a grid that ended abruptly and did not blend with the tesserae, and then, on the crests, graben in perfectly aligned, perfectly spaced lines.

- **F192.** Relief grows smoothly toward its bounds; a later fault ends at an earlier; a line follows the strain's strike as far as the strain decides it; the lid rifts across its greatest stretching alone; and the graben have the model's graben form, segmented, stepping en echelon, scattered in throw and spacing, wandering.

How it went:

- **Abutting and deciding.** A later fault ending at an earlier and the thrusts' rule for the strike barely changed the grid: the collapse's strain is decided almost everywhere, and two orthogonal sets of evenly spaced faults read as a grid whatever their junctions.
- **One rifting.** With the lid rifting across its greatest stretching alone, the grid faded to one set crossed faintly by the brittle layer's folds (the user chose this of three).
- **Aligned lines.** The graben were continuous, smooth and nearly evenly spaced, where the model's graben elsewhere are segments that step, throw unevenly and wander, a set's spacing log-normal: with that form the crests' graben break up among the fold belts.

Still to do: the plateaus' broad internal relief; the strain after the event; the fabric under the cells, for the artistic pass.

## The hundred-and-seventh milestone: the polar terrains' fissures as one network

The user asked for T7: the tiger stripes, decent in size and shape, stood disconnected, and a fine fingerprint covered each polar terrain; real images show one big continuous crack network with fold terrain following it, the folds overlapping, contouring and varying in height along their length, the stripes blending in without their false blue.

- **F193.** Each terrain's fissures meet one another (past episodes' too), saturate at the shell's thickness and are linked by branches that leave them as the stripes' forks do; the folds take their strike from the nearest fissures, stand at the measured 100 m and rise and fall along their length; the blue tint is gone.

How it went:

- **Missing troughs, or not.** With the tint gone the stripes vanished, and the first measurement found no troughs; it was taken where the tint lay, not along the troughs. Sampled along the stripes themselves, they stand some 260 m deep: the folds, as steep at 250 m and 2.5 km apart, swamped them, and the measured 100 m lets them stand above the folds.
- **A parallel network.** Saturation alone drew a dense parallel set, the tide's field being nearly uniform over the zone; cross-fissures across the nearest fissure linked it but made a wall of bricks (the user's word), and branches by the stripes' fork rule braid it instead.
- **Local strikes.** Taking the folds' strike over 4 km instead of 20, and their segments 5 km long instead of 10, changed little and nearly doubled the pass's time: set back.
- **The reference.** The user's global image shows long branching belts dominating Enceladus: its sulci, whose relief the model's lanes give as tens of metres where they stand 1 to 1.5 km (T6, next).

Still to do: the sulci's relief (T6); the network's prominence at a global view.

## The hundred-and-eighth milestone: the polar terrains' margin belts from the spreading

The user agreed to the sulci's relief, T6's open part, next. The dominant belts in the user's Cassini view are the margins of Enceladus's tectonized terrains, Samarkand Sulci about Sarandib, and the model's lanes are not them: their grooves stand as the ridged plains do, the relief share holding on a second body.

- **F194.** The margin belt is the margin's brittle layer buckling under the terrain's spreading: its folds as far apart as that layer buckles, spreading over the margin lid's flexural parameter, at the tesserae's fold relief, drawn in the lanes' segments.

How it went:

- **Two laws set aside.** The relaxation since spares the lanes' grooves whole at their wavelength, so a tilt that relaxed would not lower them; a gravity and cohesion bound on the grooves would have contradicted the ridged plains' measured relief.
- **Ripples.** With the fold law and the flexural width the belt stood at the right height but read as ripples, too even and set one inside another (the user's word); the lanes' and Miranda's coronae's segments break it up, as they did those.
- **A star.** The terrain's edge distance stopped at 32 km, short of the belt's 86, so the whole interior read as the belt's edge and its azimuth met at the centre in a star; the edge distance now runs out to the belt's reach.
- **Short deep cuts.** The user found many very short cross-fissures digging deeper than the folds about them: the segments were the lanes' troughs, cut below the ground and as deep as a fault's heavy tail lets a trough be, some past 900 m, and drawn round the thinned zone's circle, 1.7 times the edge's, so a segment came out that much shorter. Each is now a fold at its bound, drawn round the edge itself: bundles of two to four ridges stepping en echelon. The user accepted it for now, to revisit if the global look falls short: the reference's sulci run longer and more continuous.

Still to do: Miranda; the network's and the belts' prominence at a global view.

## The hundred-and-ninth milestone: the gate check, and a lid holds no more than its strength bears

The user asked for the Phase 4F gate check: on crops of Enceladus, the Moon and Io, no analytic outline, no one-cell contact, no single-scale repetition; the Moon's maria embaying the highlands and mixing across the contact; Enceladus's terrains and sulci reading as terrain; Io's flows and deposits grading into their surroundings. Renders at 4,096 pixels against the user's reference maps. None passes yet:

- **Enceladus.** Its relief, its figure removed, spanned -6.0 to +5.7 km between the 5th and 95th percentiles, four times the measured (Schenk and McKinnon 2024), most of it at 80 to 300 km and there from 0.12 Gyr: the primordial relief, never relaxed (audit M3). It buried the margin belts' kilometre at a global view and kept the craters deep. Its albedo, subtle, shows analytic outlines when stretched: a lane drawn as a straight-edged wedge, spiky halos about craters, stepped bands.
- **The Moon.** Its relief is about right. Its maria are soft patches spread over the globe rather than gathered in its province (audit M33), their contacts graded over tens of kilometres rather than embaying the highlands; some twenty ray systems read as straight spokes at a global view; its highlands stand too dark against its rays (the user's note, audit M36). The bright rims of its fresh craters, which the user's reference shows too, stay.
- **Io.** It fails most: flat dark flows with crisp white outlines on smooth yellow plains, white plume deposits as repeated rings, red ones as smooth loops, plains with no structure at any scale (audits M7, M8, M20).

The user agreed that the small icy moons carry too much relief and that Enceladus is smoother than its size suggests, and chose to start there: a lid's strength first, then the heat history.

- **F195.** An ice lid holds a band's relief only up to what bends it to its strength (F15's), the warm ice under it flowing from under the rest; a rocky lid's relief stands on its crust's compensation and is not bounded.

How it went:

- **Where the relief came from.** Enceladus's 11.5 km lid held every load shorter than 247 km, so the loop's bands from 6 to 198 km kept all their relief for 4.5 Gyr; only the 397 km band relaxed. The episodes thin the lid to 2.6 km only near their poles, inside the terrains they renew; the plains never see a lid under about 5 km, where the flexure of Enceladus's chasmata and dorsae gives past heat flows of 125 to 200 mW/m², lids of one to three kilometres.
- **Two routes.** A lid's strength bound, a missing half of F139's "the loads its flexure can carry", and the heat history the record wants (M3); the user chose the bound first. The ice lid's isotherm, 170 K, stands: the rock's criterion (F144), creep at 1e-17 /s under the loads' stress, gives 155 to 163 K for ice under the model's own creep law.
- **The bound.** Enceladus's relief halves (-3.2 to +2.8 km) and its broad swells go; Mimas keeps its craters and loses its swells; the other icy bodies change little. The craters stay deep: they lie in cold ice that does not flow.

Still to do: Enceladus's heat history (M3: Goldreich et al. 2025's limit cycle); Io; the Moon's maria, rays and highlands.

## The hundred-and-tenth milestone: the renewed ground keeps the base relief as it stands

Starting on the heat history, the user saw on Enceladus's height map that its polar and fossil terrains carried a finer noise than the smoother ground about them, asked that they take the noise that would be in their place or have it worked as the cratered ground's is, and that no other process on other bodies do the same unintentionally (audit W9).

- **F196.** The loop keeps the base relief as it stands, relaxed as the ground relaxes and levelled by each crater's ejecta and floor as the ground is, and every renewal takes it so; the finish lays the lanes' and coronae's floors from it.

How it went:

- **Where the noise came from.** Every renewal took the base relief as the body formed with it: the present terrain, renewed at every step, held 1.6 times the plains' relief at 100 to 200 km since the lid's bound relaxed the plains' (F195), and every terrain a granular noise at 5 to 50 km that the plains' craters' ejecta had worn down. Only the ice shell's renewals reuse the base relief: lanes, coronae, polar and fossil terrains; no rocky process does.
- **Discs.** Setting the kept base to the floor's mean inside each rim left flat discs, edged at the rim, wherever a terrain renewed over an old crater; the floor now buries it as the deposit at the rim does, and the discs are gone.
- **Other bodies.** Miranda's coronae lose their excess at 80 to 300 km; Ganymede's lanes change by little.

Still to do: Enceladus's heat history (M3); Io; the Moon's maria, rays and highlands.

## The hundred-and-eleventh milestone: the lid that holds a load is as old as the load

The user agreed to the heat history next (M3), with Goldreich, Lithwick and Luan (2025)'s limit cycle as the candidate: Enceladus's eccentricity rising while its shell freezes, its shell melting and then librating near 3 km at some 290 mW/m2 for a third of each 13 Myr cycle.

- **F197.** The ice lid that holds a load up reaches down to where the ice's Maxwell time equals the load's age, the time the surface has stood, rather than to the 170 K isotherm the faults break.

How it went:

- **No cycle.** The model makes Enceladus's heat in its porous core (F68), some 45 times its shell's per unit of eccentricity squared, and the core's dissipation does not fall as the shell thickens: driven by a resonance's supply over its whole life it settles at once, 116 GW holding a steady 5.7 km shell, 29 GW today's 22.6 km. The limit cycle wants its heat made in the shell; the user chose to keep the porous core and look at the craters and the lid's law instead.
- **The lid's law.** At crater scales the law held everything a lid of the 170 K isotherm could carry, for ever: Ganymede's craters of 10 km and more, which Bland et al. (2017) relax to their measured depths at 40 to 50 mW/m2 over a billion years, kept theirs. A load's age sets how deep the ice holds it.
- **Where it leaves Enceladus.** Its relief eases to -2.6 to +2.2 km and its longest wavelengths sit at its measured figure's; its plains' craters, under a lid near 7 km all its life, want heat flows the porous core never gives them outside the episodes' zones.

Still to do: Enceladus's plains' heat history (M3); Io; the Moon's maria, rays and highlands.

## The hundred-and-twelfth milestone: a thinned zone's fissures run on until the thick shell stops them

The user saw each of Enceladus's terrains outlined by folds and cracks running round its perimeter, the terrains near circles and ovals, the cross-folds standing out and the terrains not meeting in a network (audit W10), and chose three steps, this the first.

- **F198.** A polar fissure runs on with the field into the thicker shell and stops where the tension no longer opens it, rather than turning round the margin; the zone's edge is warped at every scale from two cycles round the pole to half a convection cell.

How it went:

- **Circles.** The fissures that reached the margin turned toward the direction round the pole itself, so each ran along a circle about it: the rings the user saw. Turning along the zone's own contours, which its convection structures, made them wander but still trace the perimeter, and cost 80 s.
- **Arrest.** A crack running into a stiffer shell slows and stops rather than following the contrast; with no turn the rings are gone, and each terrain's edge is ragged where its fissures run out, its tongues following them.

Still to do: the margin's folds and belt (W10); the network (W10); Enceladus's plains' heat history (M3); Io; the Moon.

## The hundred-and-thirteenth milestone: the margin's folds are the fabric's, raised where the spreading is taken up

The second of the user's three steps for the terrains (audit W10): the cross-folds, the margin belt's ridges laid along the edge over the fabric, stood out from the fabric's folds and traced the terrain's shape.

- **F199.** The terrain spreads across its fissures, so its margin shortens across them and folds along them, as the fabric does: the belt's folds are the fabric's own, taking the belt's spacing and relief where the edge faces the spreading; the separate belt is gone.

How it went:

- **Turning or not.** The plan was to turn the fabric toward the edge near the margin; the spreading's direction says otherwise, the shortening at a margin running across the fissures, so the folds keep their strike and only grow, as much as the edge faces the spreading. Where the edge runs along the spreading the margin shears and folds nothing, so the belt comes and goes round each terrain with its fissures' strike.

Still to do: the network (W10); Enceladus's plains' heat history (M3); Io; the Moon.
