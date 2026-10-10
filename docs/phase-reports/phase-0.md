# Phase 0 report: research spikes and decisions

Status: complete, 2026-10-10. Awaiting review before Phase 1 starts.

Gate from PLAN.md: spikes reproduce the published curves; decisions recorded.

## Deliverables

| Item | Where | Result |
|---|---|---|
| Stack and grid spike | `crates/exogen-core/src/grid/`, `crates/exogen-cli`, `crates/exogen-app` | Equiangular cubesphere with exact cell areas, geometric edge neighbours, seamless bilinear resampling to equirectangular; CLI spike and egui-on-wgpu window both run |
| Parameter graph specification | `docs/parameters.md` | Every user-visible node with unit, kind, derivation and source; 34 numbered formulas; plausibility rules |
| Architecture decision records | `docs/adr/0001` to `0005` | Grid, time model, stage caching, determinism, stack |
| Crater scaling spike | `research/crater_scaling.py` | See below |
| Climate spike | `research/ebm_climate.py` | See below |
| Constants re-check | deferred by decision | Each stage's constants are verified against sources when that stage is implemented |
| CI | `.github/workflows/ci.yml` | fmt, clippy with warnings denied, tests on Windows and Linux |

## Grid spike results

Face size 256 (393,216 cells), resampled to 2048 x 1024:

| Measure | Value |
|---|---|
| Total solid angle | 12.566371 sr, equal to 4π to the printed precision (exact closed form) |
| Max error, face interior | 1.8e-4 |
| Max error, cells on a face edge | 1.9e-4 |
| Max error, cells at a cube corner | 7.3e-4 |
| Fill time | 54 ms |
| Resample time | 114 ms |

Edge error equals interior error, so the ghost-ring interpolation is seamless. Corner cells carry a first-order term confined to the eight cube corners; tests pin it below 0.05 and the resampling error is confirmed second order (falls by about four per halving of cell size). Nine unit tests cover index round trip, right-handed face frames, locate round trip, exact area, area ratio, interior and edge neighbours, convergence and ghost accuracy.

## Crater scaling spike results

- Neukum chronology: N(1 km) of 3.2e-3 per km² at 3.3 Gyr and 6.3e-2 at 4.0 Gyr, matching mare and highland densities.
- Transition diameter 1/g scaling within 15% for the Moon, Mars, Earth and Venus. Mercury under-predicted (7.2 against 10.3 km), a known velocity effect; the engine will use per-class anchors with 1/g interpolation.
- Atmospheric cut-off power law gives Venus 2.4 km and Mars 13 m against observed 1.5 to 3 km and metres. Titan 0.6 km is low and is re-fit in Phase 5.
- Pi-scaling reproduces Meteor Crater within the impactor-size uncertainty and gives 111 km for a 10 km stony impactor at 20 km/s.
- Pike's complex depth law matches Copernicus, Theophilus and Langrenus within 10%; the central-peak coefficient over-predicts and is fitted in Phase 2.
- Linear diffusion alone degrades a 1 km crater more slowly than observed populations, as expected; overprinting and burial are separate processes in the replay and the diffusivity is calibrated in Phase 2.

## Climate spike results

- Earth 287.7 K global mean, Hadley edge 30.5°; Mars 212 K with an 84° cell; Venus 737 K uniform; Titan 94 K uniform. All within a kelvin of the targets once the grey optical depth is set from the observed surface temperature, which is the calibration path F21 specifies.
- The ice-albedo feedback has a snowball branch; the climate stage must start warm or ramp insolation. Recorded as a requirement.
- Locked Earth twin at 1 bar: the Williams-Kasting transport scaling over-predicts transport uncapped. A cap of 10 x Earth's coefficient gives a 63 K day-night contrast, inside the range published GCMs give; 100 x gives 8 K. The cap is a Phase 5 calibration parameter.

## Notes for Phase 1

- egui 0.36 changed its API from what most published examples show: the application trait's entry point is `ui(&mut self, ui, frame)` and panels are `Panel::top(id).show(ui, ...)`. The app crate now follows that.
- The desktop app builds on its own target directory during development to avoid lock contention with engine builds; CI builds the workspace in one go.
- The grid module is spike quality: stencil operators with proper ghost handling, cell edge lengths and the tiled detail layout are Phase 1 work.
- Carry-forward calibration items: Mercury transition diameter, Titan cut-off, central-peak height, crater diffusivity, EBM warm start, locked-world transport cap, Earth ice line (needs the seasonal cycle and threshold tuning), Mars polar temperature (needs CO2 condensation from the volatile stage).
