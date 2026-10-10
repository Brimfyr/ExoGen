# Phase 2 report: airless cratered worlds

Status: complete, 2026-10-10. Awaiting review before Phase 3 starts.

Gate from PLAN.md: statistical tests of section 6 for crater size-frequency distribution, depth-diameter and power spectral slope pass on the Moon preset; a 1,500 km and a 6,000 km body produced from the same inputs differ in crater morphology as the scaling laws predict.

## Gate results

| Gate item | Evidence |
|---|---|
| Crater size-frequency distribution | `crater::production::tests::population_matches_the_isochron_on_a_young_surface_and_saturates_on_an_old_one`: counts above 40, 80 and 160 km on a 3.3 Gyr surface match the isochron, and on a 4.5 Gyr surface they match the Hartmann equilibrium at 5% of geometric saturation. `stages::impacts::tests::loop_population_follows_the_production_function`: over 60 seeds the loop population's share above twice the threshold, its count in the last 1 Gyr and the rarity of young giant basins follow the production function, and the whole-age basin count is Moon-like. `diameter_sampling_follows_the_cumulative_distribution` checks the sampler itself. `stages::tests::the_moon_is_heavily_cratered_and_earth_is_not` checks the regime contrast, the rayed fraction and the albedo range |
| Depth-diameter | `crater::morphology::tests::classes_and_depths_follow_pike`: a 5 km simple crater has a depth of one fifth its diameter; a Copernicus-sized crater has the Pike depth, rim height and central peak; Schrödinger-sized is peak-ring; Orientale-sized is multi-ring at the basin depth cap. `profile_is_continuous_at_the_rim_and_measures_its_own_depth`, `diffusion_shallows_a_small_crater` and `crater::stamp::tests::a_stamped_fresh_crater_has_the_right_depth_and_resets_ages_under_its_ejecta` cover the profile and the stamp |
| Power spectral slope | `stages::tests::moon_height_spectrum_has_a_terrain_like_slope`: the equatorial row of the Moon preset's height map has a log-log periodogram slope between -1.3 and -3.8 over wavenumbers 4 to 128 |
| Two bodies, same inputs | `crater::morphology::tests::gravity_moves_the_transition`: a 12 km crater is simple at 1.38 m/s² and complex with a central peak at 9.2 m/s². The CLI comparison below |

Engine test suite: 59 tests, all passing. `exogen validate`: 98 checks over 14 presets, 0 failed. Clippy with warnings denied and rustfmt are clean on all three crates.

### Two bodies from the same inputs

Both project files are the Moon preset with the radius changed, the density left at 3,344 kg/m³, the composition class pinned to rock and the regime pinned to stagnant lid, rendered at preview quality with the same seed.

| Body | Gravity | Simple to complex | Peak-ring onset | Multi-ring onset | Loop threshold | Classes in the sidecar list |
|---|---|---|---|---|---|---|
| 1,500 km | 1.40 m/s² | 19.1 km | 162 km | 381 km | 73.6 km | 32 to 128 km: all complex. 128 to 256 km: 45 complex, 88 peak-ring. 256 to 512 km: 27 peak-ring, 9 multi-ring. Above 512 km: all multi-ring |
| 6,000 km | 5.61 m/s² | 4.77 km | 40.5 km | 95.3 km | 294 km | All 853 listed craters (128 km and up) multi-ring |

The transitions scale with the inverse of gravity, as the scaling law states, and the emplaced populations follow them. With the composition class left derived, the 6,000 km body decompresses its lunar density to 2.79 g/cm³ and is classed rock-ice, which applies the ice factor and puts the simple-to-complex transition at 0.71 km. That is the bulk composition model doing what the spec says, but see the Phase 4 note below.

## Deliverables

| Item | Where | Notes |
|---|---|---|
| Production and chronology | `crates/exogen-core/src/crater/production.rs` | Neukum production function and chronology, flux scale per body, record start 4.2 Gyr before present, Hartmann saturation with the surviving fraction, per-cell Poisson sampling per log-diameter bin with truncated-exponential formation times |
| Morphology | `crater/morphology.rs` | Pike dimensions (depth, rim height, central peak, floor), class transitions scaling with the inverse of gravity with an ice factor, peak-ring and multi-ring multiples, profile tables with diffusion pre-degradation, measured depth |
| Stamping | `crater/stamp.rs` | Footprint search on the cubesphere across face edges and poles, floor replacement and ejecta addition, surface age and unit reset under the ejecta, ellipticity and azimuth for oblique impacts |
| Loop emplacement | `stages/impacts.rs` | Craters above the loop threshold sampled each step and emplaced oldest first at the step they form; retention cutoff for regimes that resurface wholesale |
| Finishing pass | `stages/mod.rs` | Resample, restore the loop craters at finishing resolution, sample and replay the smaller population oldest first with pre-degradation, secondaries, rays, optical maturity, materials. Epoch rendering works through the same pass |
| Primordial crust | `stages/primordial_crust.rs` | Low-degree harmonic crustal thickness with a dichotomy term, Airy isostasy, two-level crust for mobile-lid bodies; replaces the Phase 1 placeholder terrain |
| Materials | `stages/materials.rs` | Ten end-members with fresh and mature colours and a roughness value; bedrock and regolith-mix controls; albedo blends the bedrock's fresh colour to the regolith mix's mature colour by maturity reduced by rays; continental and oceanic units; liquid colour by depth |
| Parameters | `params/registry.rs`, `params/derive.rs` | New nodes: cratering loop threshold, minimum diameter, retention age, empirical saturation, secondary fraction, oblique fraction; surface space-weathering timescale; materials bedrock, regolith mix and colour space. Every preset has its own master seed |
| Units | `units.rs` | Primordial, continental, oceanic, crater floor, ejecta, basin floor, secondary field; legend written to the sidecar |
| Export, CLI, app | `export.rs`, `crates/exogen-cli`, `crates/exogen-app` | Sidecar lists the loop craters and the larger replayed ones with position, diameter, age and class, plus the replayed count; `render --seed`; app preview channels for maturity and rays |
| Validation | `validate.rs` | Retention age and bedrock rows for ten bodies |

## Measurements

Development machine, 6-core laptop CPU, release build, standard quality (loop face 256 by 90 steps, finishing face 1024, export 2048 px):

| Preset | Loop craters | Replayed craters | Loop | Finishing | Export | Total |
|---|---|---|---|---|---|---|
| Moon | 950 | 66,795 | 0.91 s | 9.16 s | 0.37 s | 10.44 s |
| Mercury | 724 | 92,228 | 0.96 s | 9.68 s | 0.35 s | 10.99 s |
| Callisto | 480 | 52,912 | 0.77 s | 7.33 s | 0.40 s | 8.50 s |

Preview quality renders take 2.3 to 3.7 s, so the app scrubber stays responsive. Finishing time is almost all crater stamping and scales with the replayed count.

Moon population at standard quality: 34 basins above 300 km and 6 above 1,000 km from the loop, against about 43 and a handful in the lunar catalogue; preview quality with the same seed gives 44 and 6, the spread being Poisson. The youngest crater above 150 km is 186 km at 0.52 Gyr, which the production function allows at the few-tens-of-percent level. Rendering the Moon at 3 Gyr before present replays 9,868 craters instead of 24,455 at preview quality, with every crater younger than 1.4 Gyr at that epoch.

## Findings during the phase

- Extrapolating the chronology to the body's formation gave three times the Moon's basin count and frequent young giant basins. The crater record now starts 4.2 Gyr before present, the age of the oldest dated lunar surfaces; the cumulative count is flat and the rate zero before then. A 60-seed statistical test guards the population.
- Earth was flooded with secondaries from ancient basins the regime would have erased. Loop craters older than the retention age are no longer emplaced on regimes that resurface wholesale (mobile lid, episodic lid, tidally heated ice shells); stagnant lids keep theirs, since erosion removes small craters long before large ones.
- The volcanic burial rule had set the Moon's retention age to 14 Myr. Burial now applies only above a tidal flux of 0.5 W/m², the Io class.
- Rays have a constant linear width per ray and a length proportional to the square root of 85 km over the diameter, so large craters have relatively shorter rays. The ray lifetime is capped at 1.5 Gyr: Callisto's derived fading time otherwise made the ray pass take two minutes. Rays need at least six finishing cells.
- The space-weathering timescale is capped at 4 Gyr; Callisto came out uniform otherwise. The material rule changed from a single mixture to the bedrock's fresh colour blended to the regolith mix's mature colour, so fresh craters show the excavated rock.
- All presets shared the master seed 1 and so shared their crater positions; each now has its own.
- Mercury's derived space-weathering and ray-fading timescales are both 0.15 Gyr from the inverse-square distance scaling, so nearly every crater is mature and the albedo map is close to uniform. Real Mercury keeps bright rayed craters for a few hundred million years. Calibrate both when Mercury is completed in Phase 3.
- The determinism test holds with the crater stages: identical loop and finishing hashes on 1 and 4 threads.

## Deferred from the Phase 2 list

- Regolith roughness: the end-members carry a roughness value; the map is exported with the studio set in Phase 7.
- Colour integration against colour matching functions and custom end-members with spectra: Phase 7 colour management. The library holds integrated colours now.
- Impactor-to-crater scaling: the engine samples crater diameters from the crater production function directly. The pi-scaling spike in `research/crater_scaling.py` is the reference for when impactor populations are specified.
- Infill of crater floors is implicit in the oldest-first replay; volcanic infill is Phase 3.
- Disk spill of the stage cache and keyframe compaction: still in-memory only.

## Notes for Phase 3

- Flood basalts by topographic fill will darken basin floors; the basin floor unit exists to receive them. A fill must reset the surface age in the loop the way a crater stamp does, so that the age-ordered replay puts ghost craters under the fill and younger craters on top.
- Mercury's timescales and fresh-to-mature contrast, as above. Scarps are Phase 3.
- Finishing is dominated by stamping. Large bodies at standard quality will be the first place a parallel or tiled stamping pass is needed.

## Notes for Phase 4

- The ice factor on the crater scaling keys on the bulk composition class, so a dense body with an ice shell, or a rock-ice body whose surface is pinned to rock, gets the wrong transition. Key it on the surface material (the bedrock node) when the ice-world material model lands.
- The compressible composition model noted in the Phase 1 report affects the class of large bodies, as the 6,000 km test body showed.
