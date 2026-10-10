# Phase 1 report: foundation

Status: complete, 2026-10-10. Awaiting review before Phase 2 starts.

Gate from PLAN.md: a project file round-trips; the same seed produces identical output twice; derived values for the reference bodies match the Phase 0 tables; scrubbing the timeline shows the placeholder terrain changing step by step.

## Gate results

| Gate item | Evidence |
|---|---|
| Project file round-trips | `project::tests::project_round_trips_through_toml_and_rebuilds_the_same_graph`: every node value identical after save and load, pinned state preserved, newer-version files refused |
| Same seed, identical output | `stages::tests::the_same_seed_gives_identical_output_on_different_thread_counts`: loop and finishing hashes equal on 1 and 4 threads |
| Reference bodies match | `exogen validate`: 81 checks over 14 presets, 0 failed; also a unit test |
| Timeline shows change | `stages::tests::keyframes_change_over_time_and_the_cache_is_used`; the app scrubber shows each stored keyframe |

Engine test suite: 45 tests, all passing. Clippy with warnings denied and rustfmt are clean on all three crates.

## Deliverables

| Item | Where | Notes |
|---|---|---|
| Parameter graph | `crates/exogen-core/src/params/` | About 130 nodes across 13 groups. Derived by default, pin to override, reset to derive; multipliers; stage outputs; plausibility rules; canonical hashing for cache keys. A test asserts every derivation reads only the nodes it declares |
| Derivations | `params/derive.rs`, `constants.rs` | Formulas F1 to F33 of the spec in the groups implemented; every constant with a source; provisional constants marked with the phase that calibrates them |
| Presets | `crates/exogen-core/presets/*.toml` | Fourteen reference bodies as project files |
| Project files | `project.rs` | TOML with version, name, overrides; unknown ids reported, newer versions refused |
| Fields and layers | `fields.rs` | Per-cell fields, the nine-layer surface state, seamless resampling, deterministic reductions, hypsometric fill |
| Random numbers | `rng.rs` | Counter-based draws keyed by seed, stream, index and draw |
| Noise and harmonics | `noise.rs`, `spharm.rs` | Seeded gradient noise on the sphere; real spherical harmonics to degree 8 |
| Stage runner | `stages/mod.rs` | Coarse loop with keyframes, finishing pass, in-memory cache keyed by parameter hashes, progress and cancellation |
| Placeholder terrain | `stages/placeholder_terrain.rs` | Primordial crust with isostasy and a two-level crust for mobile-lid bodies, drifting perturbation, patchy resurfacing that stamps surface age, diffusion |
| Hydrosphere fill | `stages/mod.rs` | Working liquid filled to its inventory volume; sea level and covered fraction written back as stage outputs |
| Export | `export.rs` | Height EXR float32 and 16-bit PNG, albedo EXR float16 and sRGB PNG, liquid mask, surface age EXR, sidecar JSON with every parameter and its state |
| CLI | `crates/exogen-cli` | `render`, `derive`, `validate`, `preset list/show/new`, and the Phase 0 grid spike |
| Desktop app | `crates/exogen-app` | Parameter panel with groups, search, expert toggle, tooltips with basis and source, pin and reset; pipeline panel with run, cancel, progress, stage state and plausibility notes; map and software-rendered globe with channel selection, cursor readout, sun direction, drag and zoom; time scrubber with on-demand keyframe previews; project load and save; export; log. All computation on a worker thread |

## Measurements

Development machine, 6-core laptop CPU, debug-opt build:

| Case | Loop | Finishing | Export at 2048 px | Total |
|---|---|---|---|---|
| Earth, preview quality (loop face 128 x 30 steps, finishing face 512) | 0.22 s | 0.38 s | 0.35 s | 0.94 s |
| Mars, standard quality (loop face 256 x 90 steps, finishing face 1024) | 1.9 s | 0.83 s | 0.39 s | 3.1 s |

Earth at preview quality gets a liquid fraction of 0.74 with the sea level 2.5 km above the mean solid surface, which is the right order for a 2.7 km global water layer over a bimodal crust.

## Findings during the phase

- Stellar relations were switched from the older fits to Eker et al. 2018 after the source check; the Gough luminosity law, the Gladman locking formula and the radiogenic shares were confirmed against sources.
- The Jeans retention threshold evaluated at the surface temperature was far too lenient; it now uses an exobase temperature model and the six-times-thermal-speed rule (λ above 36). Provisional, Phase 5.
- The water latent heat is an effective mean over the liquid range so the boiling point lands within 0.3 K.
- The determinism test caught a parallel floating-point sum in the hypsometric fill. All reductions now use fixed-size chunks summed in order.
- A blank mobile-lid world needs a two-level crust to have continents; without it, Earth's water volume drowns almost everything. Plate tectonics replaces this in Phase 6.
- The bulk composition model (uncompressed-density mixing) classifies all fourteen reference bodies correctly but gives Europa an H2O layer of 220 km against 80 to 170 km measured. A compressible three-layer model is Phase 4 work.

## Deferred from the Phase 1 list

- Disk spill of the stage cache and keyframe compaction: in-memory only, four cache entries, at most 24 keyframes. Needed when loops get heavier in Phase 2.
- Tiling: not needed until the detail pass in Phase 7.
- Curve editor in the app: history curves display but are not editable; Phase 3 when the eccentricity history matters.
- File dialogs: paths are typed. Phase 9.

## Notes for Phase 2

- The placeholder process shows coarse-cell blockiness in its resurfacing pattern after resampling. It is a placeholder; the cratering stage replaces it and small features are emplaced at the finishing resolution.
- The reference-check table in `validate.rs` is the place to add crater statistics when the cratering stage lands.
- The app's preview finishing at face size 256 takes well under a second per keyframe, so the scrubber feels immediate once a keyframe is cached.
