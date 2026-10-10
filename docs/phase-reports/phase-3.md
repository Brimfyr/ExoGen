# Phase 3 report: interior, volcanism, stagnant and episodic lid

Status: complete, 2026-10-11. Awaiting review before Phase 4 starts.

Phase 3 gives every airless rocky world its interior history and lets that history shape the surface: cooling shrinks the body into thrust scarps, magma floods the basins and builds rises, shields and canyons, a foundering lid renews a whole surface, and a tidally heated moon resurfaces itself faster than impacts can mark it. Eight milestones went in, each committed to main on its own: the parallel finishing pass, the thermal history, contraction scarps, flood volcanism, lava subsidence with wrinkle ridges, plume rises, the episodic lid, and the Io regime. The finishing pass that Phase 2 built, coarse evolution plus resolved detail, carries all of it: every new process runs inside the loop at the step it happens and is sharpened against the resolved terrain afterwards, in time order with the craters.

## Gate results

The plan's gate was that the Mars preset produces a rise, a canyon system and giant shields consistent with the relief limit, that the Venus preset's crater population is sparse and uniform, and that the Io preset has no craters and hundreds of paterae. All three hold: Mars's tallest shield reaches 22.9 km under a 25 km relief limit; Venus carries 25 craters above the recording size, fresh, spread evenly and none older than the event; Io has 403 paterae and no crater. Each gate below is a test in `crates/exogen-core/src/stages`, run at preview quality on the reference presets; the full suite is 87 tests and `exogen validate` passes its 104 reference checks.

- **Mercury takes up its contraction in scarps.** The scarps stamped in the loop carry more than 30% of the shortening budget 2 ΔR / R × area (the rest is in scarps shorter than the loop grid), the longest exceeds 300 km and the highest 1.5 km, and the Moon has far less shortening per area. At standard quality Mercury carries 2,272 scarps, the longest 843 km, from 5.7 km of radial contraction (Byrne et al. 2014 give 5 to 7 km).
- **The Moon floods its thin-crust basins.** Lava covers between 5 and 35% of the surface, most of it on the thin-crust half, and the finishing pass keeps the loop's extent within 5%. At standard quality the maria cover 14.8% of the Moon (17% in reality) from 7,731 eruptions totalling 1.6e7 km³ (about 1e7 km³ in reality), with 1,561 wrinkle ridges on them.
- **Mars raises a rise with shields and a canyon.** One rise stands more than 0.3 of its height above the ground beyond it, the tallest shield more than half its height above its foot, the canyon floor more than 3 km below its rims, and the Moon raises none. At standard quality: a 4,400 km rise 5.9 km high, five shields (the largest 523 km across and 22.9 km high, with three calderas), a canyon 3,484 km long, 200 km wide and 7 km deep, 68 graben.
- **Venus renews its surface and keeps tessera.** The event lies at the retention age, the tessera make up 4 to 14% of the surface and keep an older age, more than 70% of the plains carry the event's age, there are more than 50 coronae and some domes, and no crater is older than the event. At standard quality: tessera 7.9%, plains 91.8%, 506 coronae, 152 domes, three rises with nine shields and one chasma, 78,000 wrinkle ridges, 765 flows since the event.
- **Io is craterless with paterae and mountains.** More than 100 paterae and 20 mountains, flow fields and plume rings, not one crater, fresh flows between 0.5 and 30% of the surface, and a mean colour that is yellow. At standard quality: 403 paterae (85 active), 149 mountains up to 15 km, 140 flow fields, 25 plume rings.
- **Unit gates.** The thermal history reproduces the present-day interior nodes within 2% and the young Moon is hotter with a thinner lithosphere; a scarp stamped on a flat grid steps up behind its face and tapers at its ends; a lava lake runs downhill, fills the pit and keeps its volume, and overflows a pass into the next pit; ridges laid out on a filled basin run along its contours at the planned length per area; depression filling levels pits and keeps slopes; the shield, trough, corona, dome, patera and block profiles have their shapes; the loop population still follows the production function with the steady-flux cap.

## Deliverables

- **Thermal history** (`thermal.rs`): heat flux, tidal flux from the eccentricity history, mantle temperature, lithosphere and elastic thickness, a relief limit that follows the elastic thickness, cumulative contraction and magma production, sampled through the loop and exposed as curve nodes.
- **Contraction scarps** (`stages/tectonics.rs`): a power-law population of lobate thrust scarps whose shortening equals the strain, each stamped at the step it forms so it cuts what is older and is cut by what is younger, with a 12° face, a long back slope and three harmonics of lobing along the trace; the finish restores the face. The longest scarp scales with the square of the strain, which gives Mercury thousand-kilometre ruptures and the Moon kilometre-scale ones below map resolution. The contraction law gained a 0.7 Gyr onset and a mantle term that scales down for bodies that accreted cool inside, so the Moon's record is a few hundred metres (Watters et al. 2010).
- **Flood volcanism** (`stages/volcanism.rs`): the extrusive share of the magma budget, cut back once the lithosphere passes 400 km or contraction closes the pathways, erupts from vents that favour thin crust and low ground; large impacts now thin the crust beneath them, so basins draw the eruptions. Each eruption fills one lake from the nearest low by priority fill, crossing passes into deeper ground, spreading no thinner than a flow and no further than 1,500 km, and subsiding by 5% of its thickness, which shows buried rims as ghost rings. Floods reset the surface age and unit; craters on lava take lava units. The finishing pass keeps a ledger of every flood's coarse smear and refills each flood to its level against the resolved ground inside a lobate outline of its footprint, so contacts are sharp and flow fronts lobate. Mare basalt and Mercury's brighter smooth plains join the materials.
- **Wrinkle ridges** (`stages/ridges.rs`): laid out at the end of the loop where the lava fill thins fastest and round a rise's load, following the thickness contours with inertia and wobble, a broad arch with a narrow wandering crest, standing 0.2 Gyr after the last flow they cross and cut by younger craters.
- **Plume rises** (`stages/plume.rs`): a dome grown over the first 1.5 Gyr that lifts what is on it, radial graben, a canyon system of overlapping scalloped troughs on a wide rise, and shields with 5° slopes capped by the strength-limited relief, with nested calderas and radial flow texture cut at the finish. Eruptions favour the rise.
- **Episodic lid** (`stages/episodic.rs`): the event fills every depression 70% of the way to its spill point and keeps the highest-standing 8% as tessera with a fabric of crossing ridges, refined at the finish against the resolved ground; afterwards a small share of the budget erupts, and coronae (annular ridges with concentric fractures round a dome or sag) and pancake domes dot the plains. Rises form on episodic lids too, narrower than on a thick stagnant lid.
- **Io regime** (`stages/io.rs`): on a rocky body above 0.5 W/m² of tidal flux the surface is laid down in the last step: paterae concentrated at low latitudes, a fifth active with black floors; dark flow fields with soft lobate margins; plume rings of red sulfur and white frost; elongated tilted mountain blocks with ridged texture; sulfur plains reddening towards the poles with patches of frost. Deposits and flows blend over the plains by a coverage carried in the melt field. No crater is kept.
- **Finishing pass**: parallel waves of non-overlapping stamps, byte-identical to the sequential pass; events of every kind interleaved in time order; floods sampled through a sparse padded field, which brought Mars from 147 to 50 seconds.
- **Outputs**: the sidecar lists scarps (with strike), the largest four thousand floods and their count, ridges, rises, shields, canyons, graben, coronae, domes, paterae, mountains, the lava area fraction, the tessera fraction and the resurfacing age; the exporter adds `units.png` and `surface_age.png`; the shaded preview is exposed to a common mean under a 30° sun.
- **Parameter graph**: 20 new nodes, documented in `docs/parameters.md` with formulas F35 to F37 and a revised F16.

## Measurements

Standard quality, 8192 × 4096 output, the presets' own seeds.

| Body | Loop | Finish | Signature features |
|---|---|---|---|
| Moon | 3.0 s | 18.9 s | maria 14.8% of the surface, 1,561 ridges, no visible scarps |
| Mercury | 3.3 s | 20.4 s | 2,272 scarps to 843 km, smooth plains 16.6%, 3,958 ridges |
| Mars | 7.0 s | 50.2 s | rise, 5 shields, canyon 3,484 × 200 × 7 km, 68 graben, plains 36.8%, 12,777 ridges |
| Venus | 1.5 s | 5.0 s | tessera 7.9%, 506 coronae, 152 domes, 3 rises, 1 chasma, 78,183 ridges |
| Io | 0.7 s | 2.2 s | 403 paterae, 149 mountains, 140 flow fields, 25 plume rings, 0 craters |

Mercury's smooth plains come to 17% against 27% on the real planet; the budget is right within the uncertainty of the extrusive share and the shut-off by contraction, which ends Mercury's volcanism 3.5 Gyr ago as observed. Mars's 160 craters above 300 km are what the lunar production function scaled by area and flux predicts; the real planet shows some 20 of them because its early erosion and the burial of the lowlands hid the rest, which waits for Phase 5.

## Findings during the phase

- **Lakes, not ponds.** The first flood model let lava pond in the nearest pit, spill and run on; on cratered ground it hopped through tiny pits until the hop limit and dropped most of its volume, and the chains of filled pits read as staircases. A single hydrostatic lake per eruption, grown from the nearest low in height order, fills a basin fed many times as one lake and keeps every cubic metre. A reach-limited version in between filled sub-regions of a basin to different levels and left the hummocks between them dry.
- **The ledger.** Refilling floods at the finish against the final coarse state fails when floods overlap, because the base carries every later flood's smear. Keeping a running sum of the smears already taken out, and leaving cells under younger craters to the base by the craters' fine footprints rather than the base's blurred ages, made overlapping floods compose exactly.
- **Restores and later surfaces.** A loop crater's fine restore used to rewrite its floor and ejecta units over cells that lava had covered later; the flood pass then took those cells for younger craters and left holes below the lava. The restore now leaves the surface of cells resurfaced after the crater alone.
- **Late giants.** The lunar production function extrapolated to the steady era put a 1,400 km basin on Mars 2.5 Gyr ago and 550 km craters on its shields. After the heavy bombardment the asteroid belt has no bodies left for such craters: the steady flux now falls two slopes steeper above 150 km and makes nothing over 600 km (F37).
- **Blocky margins.** Every coarse footprint refined at the finish needs a smooth outline, or its edge shows the loop's cells. Floods use the smeared indicator of their footprint at a threshold below a half with lobes of noise; troughs and calderas are carried at the loop's wall width and refined with the same width in the subtraction, or the residual ramp stays.
- **Contraction timing.** With contraction starting at formation the scarps formed under the heavy bombardment and were chopped to pieces; the onset and the cool-interior scaling moved the record to where the planets keep it.

## Deferred from the Phase 3 list

- Lava channels and rilles: a few kilometres wide, below the standard finishing cell on Venus and Mars; Phase 7 detail synthesis.
- Stratovolcanoes belong to arcs, so to the mobile lid in Phase 6.
- Plains flow-field texture, lake levels across a pass inside one lake (the hydrostatic approximation of F36), and a more varied mare colour by titanium.
- The relaxation of early basins under the hot lithosphere, noted in Phase 2, now has the lithosphere history to draw on but is not yet applied.

## Notes for Phase 4

- Ice worlds get their own regimes; the flood process and the Io process are gated by composition and tidal flux and will not run on them. Cryovolcanism should reuse the lake fill with the cryolava material.
- The thermal history's relief limit follows the elastic thickness; for ice shells the elastic thickness node already uses the 170 K isotherm.
- Mercury's weathering timescales (0.15 Gyr from the inverse-square scaling) still want calibration when Mercury's colour is reviewed; space weathering now saturates in 0.3 Gyr at 1 AU.
