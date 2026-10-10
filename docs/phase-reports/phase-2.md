# Phase 2 report: airless cratered worlds

Status: complete, 2026-10-10, revised the same day after review. Awaiting review before Phase 3 starts.

Gate from PLAN.md: statistical tests of section 6 for crater size-frequency distribution, depth-diameter and power spectral slope pass on the Moon preset; a 1,500 km and a 6,000 km body produced from the same inputs differ in crater morphology as the scaling laws predict.

Two checks were added during the review: craters erase each other, measured as the share of loop craters that lose their relief and the count of recognisable basins; and bodies of any age get the same early bombardment at formation and the steady flux after it.

## Gate results

| Gate item | Evidence |
|---|---|
| Crater size-frequency distribution | `crater::production::tests::population_matches_the_isochron_on_a_young_surface_and_saturates_on_an_old_one`: counts above 40, 80 and 160 km on a 3.3 Gyr surface match the isochron, and on a 4.5 Gyr surface they match the Hartmann equilibrium at 5% of geometric saturation. `stages::impacts::tests::loop_population_follows_the_production_function`: over 60 seeds the loop population's share above twice the threshold, its count in the last 1 Gyr and the rarity of young giant basins follow the production function, and the whole-age basin count is Moon-like. `diameter_sampling_follows_the_cumulative_distribution` checks the sampler. `stages::tests::the_moon_is_heavily_cratered_and_earth_is_not` checks the regime contrast, the rayed fraction and the albedo range |
| Flux history | `crater::production::tests::flux_history_matches_the_neukum_chronology`: the lunar record integrated from solidification equals the chronology's value 4.2 Gyr before present, the last 1 Gyr and 3.3 Gyr match the spike values, and the formation-time inverse round-trips. `bodies_of_any_age_get_the_early_bombardment_at_formation`: a 10 Gyr record is a 6 Gyr record plus 4 Gyr of steady flux |
| Depth-diameter | `crater::morphology::tests::classes_and_depths_follow_pike`: a 5 km simple crater has a depth of one fifth its diameter; a Copernicus-sized crater has the Pike depth, rim height and central peak; Schrödinger-sized is peak-ring with one ring; Orientale-sized is multi-ring at the basin depth cap with three rings at the Outer Rook, Inner Rook and inner-ring diameters. `profile_is_continuous_at_the_rim_and_measures_its_own_depth`, `diffusion_shallows_a_small_crater` and `crater::stamp::tests::a_stamped_fresh_crater_has_the_right_depth_and_resets_ages_under_its_ejecta` cover the profile and the stamp |
| Craters erode each other | `crater::stamp::tests::ejecta_buries_older_relief`: a small crater under a large crater's ejecta loses a good part of its relief, and without burial it only gets draped. `stages::tests::later_impacts_erase_some_loop_craters_and_most_basins_stay_visible`: between 2% and 50% of the loop craters lose their relief, while at least half of the basins stay recognisable with kilometres of relief. `stages::tests::nothing_older_than_a_loop_crater_survives_inside_it`: no cell inside a loop crater is older than the crater |
| Oblique impacts | `crater::morphology::tests::impact_angles_follow_sin_two_theta_and_few_cavities_are_elliptical`: mean angle 45°, the share below 12° is sin² 12°, and 5% of cavities are elliptical. `crater::stamp::tests::oblique_ejecta_is_thicker_downrange_and_absent_uprange`: a vertical stamp is symmetric to 5%; a 10° stamp is 1.3 times thicker downrange, below 0.4 uprange and raised to the sides |
| Power spectral slope | `stages::tests::moon_height_spectrum_has_a_terrain_like_slope`: the equatorial row of the Moon preset's height map has a log-log periodogram slope between -1.3 and -3.8 over wavenumbers 4 to 128 |
| Two bodies, same inputs | `crater::morphology::tests::gravity_moves_the_transition`: a 12 km crater is simple at 1.38 m/s² and complex with a central peak at 9.2 m/s². The CLI comparison below |

Engine test suite: 65 tests, all passing. `exogen validate`: 101 checks over 14 presets, 0 failed. Clippy with warnings denied and rustfmt are clean on all three crates.

### Two bodies from the same inputs

Both project files are the Moon preset with the radius changed, the density left at 3,344 kg/m³, the composition class pinned to rock and the regime pinned to stagnant lid, rendered at preview quality with the same seed.

| Body | Gravity | Simple to complex | Peak-ring onset | Multi-ring onset | Loop threshold | Classes in the sidecar list |
|---|---|---|---|---|---|---|
| 1,500 km | 1.40 m/s² | 19.1 km | 162 km | 381 km | 73.6 km | 32 to 128 km: all complex. 128 to 256 km: 45 complex, 88 peak-ring. 256 to 512 km: 27 peak-ring, 9 multi-ring. Above 512 km: all multi-ring |
| 6,000 km | 5.61 m/s² | 4.77 km | 40.5 km | 95.3 km | 294 km | All 853 listed craters (128 km and up) multi-ring |

The transitions scale with the inverse of gravity, as the scaling law states, and the emplaced populations follow them. With the composition class left derived, the 6,000 km body decompresses its lunar density to 2.79 g/cm³ and is classed rock-ice, which applies the ice factor and puts the simple-to-complex transition at 0.71 km. That is the bulk composition model doing what the spec says, but see the Phase 4 note below.

### Bodies of any age

The same Callisto and Mercury project files with the body age pinned, preview quality, one seed. The record begins at the solidification time after formation, so an older body has the same early bombardment and more steady flux; the basin counts differ only by the Poisson draw.

| Body | Age | Basins above 300 km formed | Visible | Oldest crater |
|---|---|---|---|---|
| Callisto | 4.52 Gyr | 71 | 65 | 4.27 Gyr |
| Callisto | 6.0 Gyr | 65 | 60 | 5.75 Gyr |
| Mercury | 4.52 Gyr | 55 | 53 | 4.16 Gyr |
| Mercury | 10.0 Gyr | 45 | 44 | 9.64 Gyr |

## Deliverables

| Item | Where | Notes |
|---|---|---|
| Production and flux history | `crates/exogen-core/src/crater/production.rs` | Neukum production function; the flux history in time since formation, a steady rate plus an early component decaying with a 144 Myr e-fold, which is the Neukum chronology rewritten with the present 4.5 Gyr after the Moon formed; the record starts at the solidification time; Hartmann saturation with the surviving fraction; per-cell Poisson sampling per log-diameter bin with formation times that follow the rate |
| Morphology | `crater/morphology.rs` | Pike dimensions, class transitions scaling with the inverse of gravity with an ice factor, a peak ring at half the radius, multi-ring basins with a ring for each sqrt(2) of diameter above the onset at successive 1/sqrt(2) of the rim radius, profile tables with diffusion pre-degradation, impact angles from sin 2θ and cavity elongation below the elliptical threshold |
| Stamping | `crater/stamp.rs` | Footprint search across face edges and poles; floor replacement and ejecta addition; burial of the relief under ejecta by its thickness; azimuthal ejecta asymmetry for oblique impacts; scalloped rims; a per-crater noise lattice for hummocky ejecta, basin sculpture and rough floors; impact melt placement; surface age and unit reset under the ejecta; a loop-crater restore that returns the cavity to the coarse base and buries the finer relief |
| Loop emplacement | `stages/impacts.rs` | Craters above the loop threshold sampled each step and emplaced oldest first; retention cutoff for regimes that resurface wholesale; residual relief measured at the end of the loop and recognisable craters marked |
| Finishing pass | `stages/mod.rs` | Resample; one age-ordered sequence of loop-crater restores and replayed craters with each crater pre-degraded by the gardening that followed it; secondaries in radial chains and clusters; rays with the oblique asymmetry, a ragged continuous-ejecta edge and brightness mottling; optical maturity; materials. Epoch rendering works through the same pass |
| Primordial crust | `stages/primordial_crust.rs` | Low-degree harmonic crustal thickness with a dichotomy term, Airy isostasy, two-level crust for mobile-lid bodies |
| Materials | `stages/materials.rs` | Thirteen end-members with fresh and mature colours and a roughness value; bedrock and regolith-mix controls; the deep bedrock that large craters excavate, with ice brought up from under a lag darkening back to the lag; impact melt that converges on the regolith with age; continental and oceanic units; liquid colour by depth |
| Parameters | `params/registry.rs`, `params/derive.rs` | New nodes: cratering loop threshold, minimum diameter, retention age, empirical saturation, secondary fraction, oblique fraction, early flux ratio and decay; surface solidification time and space-weathering timescale; materials bedrock, regolith mix, deep bedrock, deep bedrock depth and colour space. The inert late-spike node was removed. Every preset has its own master seed |
| Units | `units.rs` | Primordial, continental, oceanic, crater floor, ejecta, basin floor, secondary field, deep-excavating floor, deep-sourced ejecta, impact melt; legend in the sidecar |
| Export, CLI, app | `export.rs`, `crates/exogen-cli`, `crates/exogen-app` | Sidecar lists loop craters and the larger replayed ones with position, diameter, age, class, impact angle, residual relief and visibility; writes retry briefly when a viewer holds a file; `render --seed`; app preview channels for maturity and rays |
| Validation | `validate.rs` | Retention age, bedrock, deep bedrock and solidification rows |

## Measurements

Development machine, 6-core laptop CPU, release build, standard quality (loop face 256 by 90 steps, finishing face 1024, export 2048 px):

| Preset | Loop craters | Replayed craters | Basins above 300 km formed, visible | Loop craters still visible | Loop | Finishing | Total |
|---|---|---|---|---|---|---|---|
| Moon | 953 | 61,335 | 31, 28 | 799 of 953 | 1.1 s | 17.0 s | 18.5 s |
| Callisto | 690 | 61,619 | 67, 63 | 585 of 690 | 1.1 s | 17.1 s | 18.6 s |
| Mercury | 491 | 76,198 | 59, 46 | 380 of 491 | 1.1 s | 18.0 s | 19.6 s |

The finishing pass took 7 to 10 s before the review; burial, the age-ordered restores, texture and the wider pre-degradation blurs added the rest. It is single-threaded. Preview quality takes 2 to 4 s.

The Moon's basin count sits below the lunar catalogue's roughly 43 above 300 km at this seed; the 60-seed test puts the expectation at 43 with Poisson spread. Six to twenty percent of the loop craters lose their relief to later cavities and ejecta, which is the overprinting and burial part of crater-on-crater erasure; basins are too deep for burial alone and wait for relaxation under the hot early lithosphere, which the Phase 3 heat budget provides.

## Findings during the phase

- Extrapolating the chronology to the body's formation gave three times the Moon's basin count. The first fix was a hard start of the record 4.2 Gyr before present; the review replaced it with the flux history in time since formation and a derived, pinnable solidification time, so a 6 Gyr Callisto or a 10 Gyr Mercury gets its heavy bombardment at formation and the steady flux after.
- The finishing pass restored all loop craters first and then stamped the replayed population in its own order, so secondaries of older craters wrote old surface ages into fresh floors and read as dark specks. Loop craters now take their place in one age-ordered sequence, and a test guards it.
- Ejecta was draped over older relief; it now buries it by its thickness against the azimuthal mean, lows filling faster than highs. Pre-degradation follows the gardening that came after each crater, so craters from the early bombardment are heavily softened and recent ones stay crisp, with the lunar diffusivity at the present flux as the anchor.
- Every crater had symmetric ejecta and rays. Impact angles now follow sin 2θ; only the 5% below about 12° elongate, but half of all impacts thicken their ejecta downrange, and the lowest open an uprange zone of avoidance and side lobes.
- Basins were smooth bowls with one ring and every world showed two colours. Multi-ring basins now carry their rings, rims are scalloped, ejecta and floors are textured, secondaries form chains, and a deep-bedrock material gives the Moon noritic basin floors and ejecta, Mercury low-reflectance patches and Callisto bright clean-ice craters; impact melt mottles fresh floors.
- Earth was flooded with secondaries from ancient basins the regime would have erased. Loop craters older than the retention age are not emplaced on regimes that resurface wholesale; stagnant lids keep theirs.
- The volcanic burial rule had set the Moon's retention age to 14 Myr; burial now applies only above a tidal flux of 0.5 W/m², the Io class.
- Rays have a constant linear width and a length proportional to the square root of 85 km over the diameter; the ray lifetime is capped at 1.5 Gyr and the space-weathering timescale at 4 Gyr.
- Mercury's derived space-weathering and ray-fading timescales are both 0.15 Gyr from the inverse-square distance scaling, so nearly every crater is mature and the albedo map stays low in contrast. Calibrate both when Mercury is completed in Phase 3.
- Callisto's 4 Gyr weathering cap leaves an age contrast between its oldest terrain and younger ejecta that the real lag-covered surface does not show; the Phase 4 lag model replaces it.
- The determinism test holds through all of this: identical loop and finishing hashes on 1 and 4 threads.

## Second revision

A second look at the renders after the review found four more faults, all fixed:

- The map sets in a folder could disagree: the exporter wrote the height maps first, so a run that failed on a file held open by a viewer left a new height map beside old albedo and surface-age maps. The exporter now checks that every output can be opened before it writes any.
- Wedge patterns radiating from crater centres, read as rays from a point, in both the floor texture and the melt. The texture lattice began at the floor radius, so every floor cell clamped to one lattice row and the pattern depended on azimuth alone. The lattice now covers the whole interior, and a test checks that the floor texture varies with radius.
- Hard-edged clumps of melt and a hemisphere-scale two-tone from the deep-ejecta unit of the largest basin. Impact melt and deep-sourced material are now continuous layers: melt is a soft fraction, and the deep fraction follows ballistic mixing, all primary ejecta near the rim and all local material at the continuous edge, with later craters inheriting the material they excavate. Basin sculpture is finer and weaker and applies to multi-ring basins only.
- Rings that were scaled copies of the rim. Each ring now has its own scalloped shape and a few gaps around the azimuth, and degradation widens and lowers it. A first attempt made them too low, wide and broken, and let four rings converge toward the centre; they are now narrower, taller, mostly continuous and at most three, the innermost at about a third of the rim radius as at Orientale. Fully accurate ring prominence waits for relaxation in Phase 3.
- The floor texture still converged to a point at crater centres after the lattice fix, because a polar lattice shrinks its cells without limit toward the centre. The floor now uses a planar lattice on the crater's tangent plane, continuous through the centre, for both roughness and melt.
- Central peaks were a single bump about one finishing cell wide, erased by degradation and burial. They are now a cluster of one to three massifs with a basal radius of 0.22 crater radii, the measured peak diameters, placed by the stamp and degraded with everything else.

Standard-quality renders now take 22 to 28 s; the finishing pass remains the cost centre and the parallel stamping noted below is the answer.

## Deferred from the Phase 2 list

- Regolith roughness: the end-members carry a roughness value; the map is exported with the studio set in Phase 7.
- Colour integration against colour matching functions and custom end-members with spectra: Phase 7 colour management. The library holds integrated colours now.
- Impactor-to-crater scaling: the engine samples crater diameters from the crater production function directly. The pi-scaling spike in `research/crater_scaling.py` is the reference for when impactor populations are specified.
- Explicit neighbour degradation after Minton and colleagues 2019, each crater degrading its surroundings in proportion to its radius squared: the flux-scaled gardening and burial cover the dominant effects; this is the refinement if equilibrium depth-diameter statistics need it.
- Disk spill of the stage cache and keyframe compaction: still in-memory only.

## Notes for Phase 3

- Flood basalts by topographic fill darken basin floors; the basin floor and deep floor units and the crater records mark where. A fill must reset the surface age in the loop the way a crater stamp does, so that the age-ordered replay puts ghost craters under the fill and younger craters on top.
- The heat budget gives the lithosphere through time: derive the solidification time from it and relax the early basins on the hot lithosphere, which is what should bring the visible basin count down from the formed count.
- Mercury's timescales and fresh-to-mature contrast, as above. Scarps are Phase 3.
- The finishing pass is the cost centre: stamping is sequential because stamps overlap. Partition the sphere into tiles with halos or batch non-overlapping craters to parallelise it.

## Notes for Phase 4

- The ice factor on the crater scaling keys on the bulk composition class, so a dense body with an ice shell, or a rock-ice body whose surface is pinned to rock, gets the wrong transition. Key it on the surface material when the ice-world material model lands.
- The compressible composition model noted in the Phase 1 report affects the class of large bodies, as the 6,000 km test body showed.
- Relaxation of icy basins into palimpsests and the lag timescale that replaces the weathering cap.
