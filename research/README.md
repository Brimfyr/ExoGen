# Research scripts

Python scripts that check equations against published anchor values before they are
implemented in the engine. They are not part of the build. Plots and reports go to
`research/output/`, which is not tracked.

Setup (once):

```
py -3.10 -m venv research/.venv
research/.venv/Scripts/python -m pip install numpy scipy matplotlib
```

Run:

```
research/.venv/Scripts/python research/crater_scaling.py
research/.venv/Scripts/python research/ebm_climate.py
```

## crater_scaling.py (Phase 0 spike, 2026-10-10)

Checks the Neukum production function and chronology, pi-group scaling, Pike morphometry,
McGetchin ejecta, the gravity scaling of the simple-to-complex transition, the atmospheric
crater cut-off and diffusional degradation.

Findings:

- Chronology gives N(1 km) of 3.2e-3 per km² at 3.3 Gyr and 6.3e-2 at 4.0 Gyr, matching
  mare and highland densities. The polynomial and the chronology agree at 1 Gyr to 2.5%.
- The 1/g scaling of the transition diameter reproduces the Moon, Mars, Earth and Venus
  within 15%. Mercury is under-predicted (7.2 km against 10.3 km observed), a known
  effect of its high impact velocity; the engine uses per-class anchors with 1/g
  interpolation between them. The ice factor 0.15 gives 2.8 to 3.2 km for the Galilean
  satellites against 2.5 km observed.
- The column-mass power law anchored at Earth gives a Venus cut-off of 2.4 km (observed
  1.5 to 3 km) and Mars 13 m. Titan comes out at 0.6 km, low against the observed dearth
  of craters below several km, which erosion and burial also cause; re-fit in Phase 5.
- Pi-scaling with Collins et al. 2005 gives Meteor Crater 1.55 km for a 50 m iron
  impactor at 12.8 km/s (observed 1.19 km; a 40 m impactor fits) and 111 km for a 10 km
  stony impactor at 20 km/s. Impactor sizes for a 100 km crater range from 5.4 km on
  Mercury to 11 km on Mars, driven by impact velocity.
- Pike's complex depth law matches Copernicus, Theophilus and Langrenus within 10%. Tycho
  is deeper (young) and Petavius shallower (old, floor-fractured), which is the degradation
  spread the engine must reproduce. The central-peak height law over-predicts Copernicus
  (1.9 km against 1.2 km): fit the coefficient in Phase 2.
- Linear diffusion with the lunar diffusivity of 5.5 m² per Myr reduces a 1 km crater
  from d/D 0.24 to 0.14 over 3.6 Gyr. Observed small-crater populations degrade faster
  than that, because overlapping impacts and ejecta burial add to diffusion; the replay
  handles those separately, and the diffusivity is calibrated in Phase 2 against published
  degradation sequences.

Feeds: Phase 2 constants file, cratering stage, validation tests.

## ebm_climate.py (Phase 0 spike, 2026-10-10)

One-dimensional energy balance on the sphere in the angle the forcing depends on (latitude
for rotators with annual-mean insolation, angle from the substellar point for locked
bodies), grey greenhouse, Williams-Kasting transport scaling, Held-Hou Hadley width.

Findings:

- With the grey optical depth set from the observed surface temperature, the model
  returns Earth 287.7 K (equator 296, 60° 277, pole 270), Mars 212 K (equator 224, pole
  183), Venus 737 K uniform, Titan 94 K uniform. Earth's polar temperatures are too warm
  because the warm start keeps the poles ice-free; the ice line needs the seasonal cycle
  and a lower ice threshold, which is Phase 5 calibration. Mars's poles are too warm
  because CO2 condensation is not in this spike; it belongs to the volatile stage.
- Held-Hou gives an Earth Hadley edge of 30.5°, Mars 84° (a hemisphere-wide cell, as
  observed), Venus and Titan a single cell to the pole.
- The ice-albedo feedback has two branches. A cold start lands on the snowball branch;
  the engine must start warm or ramp insolation. Recorded as a requirement for the
  climate stage.
- For a locked Earth twin on a 10 day orbit with 1 bar, the Williams-Kasting scaling
  over-predicts transport when uncapped: at 100 x Earth's coefficient the day-night
  contrast is 8 K; at 10 x it is 63 K (substellar 325 K, antistellar 261 K), which is in
  the range published GCMs give. The cap is a calibration parameter for Phase 5, fitted
  against Joshi et al. 1997, Yang et al. 2013 and Haqq-Misra et al. 2018.
- The airless locked case gives a substellar 381 K, a floor-limited night side and water
  ice beyond 75° from the substellar point.

Feeds: Phase 5 climate stage, parameter spec F21, F22, F32, F33.
