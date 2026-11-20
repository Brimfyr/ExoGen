# ExoGen

Process-based planetary map generator: full-planet texture sets for any body gravity has rounded, made by simulating the processes that shape real worlds.

ExoGen produces texture sets (albedo, height, materials, masks, climate and ice layers) by simulating the processes that shape real planetary surfaces: impact cratering, tectonics, volcanism, climate, erosion, volatile transport and weathering, scaled to the body's size, composition, thermal state, atmosphere and orbital environment. The maps are for physically based renderers, game engines and studio pipelines.

## Status

In development, not yet released. The design and the full roadmap are in [PLAN.md](PLAN.md); progress milestone by milestone is in the phase reports, [docs/phase-reports/](docs/phase-reports/).

**Now:** 5C, water and ice at work.

| Phase | What it builds | State |
|---|---|---|
| 0 | Research spikes and decisions | Done |
| 1 | Foundation: the parameter graph, the grid and its fields, the stage runner, export, the command line and the desktop app | Done |
| 2 | Airless cratered worlds | Done |
| 3 | Interiors, volcanism, stagnant and episodic lids | Done |
| 4 | Ice worlds | Done |
| 4F | The earlier phases' processes revised for fidelity on every body | Done |
| 5 | Atmosphere, climate, volatiles and erosion | Under way |
| 5A | Climate on every body: sunlight, the energy balance, days and seasons, the spin's history | Done |
| 5B | Volatiles: frost and caps, sea ice, the climate through the star's life, the atmosphere's loss | Done |
| 5C | Water and ice at work: wind, moisture and clouds, snow, the seas' currents and ice, rivers and lakes, glaciers, slope wasting, Titan's methane cycle | Under way |
| 5D | Dust, sand and chemistry | Planned |
| 5E | The icy bodies' first shell from the melt at their birth, dust from outside, and surface physics on every body | Planned |
| 5F | Io, Thaumasia and molten surfaces | Planned |
| 5G | The weather model, cloud types and their maps | Planned |
| 6 | Plate tectonics and the rocky first crust | Planned |
| 6N | The noise pass: field-modulated noise over every body, calibrated on all of them together | Planned |
| 6R | Equilibrium figures, rings and multiple light sources | Planned |
| 6G | Deep atmospheres, from hycean worlds to giants | Planned |
| 7 | Detail synthesis to 16k and 32k, studio export | Planned |
| 8 | Vegetation and biomes | Planned |
| 9 | Polish and release | Planned |

This section changes when a phase's or a sub-phase's gate passes, not with each milestone; a test holds its states to PLAN.md's.

A command-line tool renders the map set for a project or preset. A desktop app shows the parameter graph, a preview globe and the pipeline's state.

```
cargo build --release -p exogen-cli    # the binary is target/release/exogen
exogen preset list
exogen derive --preset Titan
exogen render --preset Mars --width 2048 --out out/Mars
exogen render --preset Earth --width 2048 --draft --out out/Earth   # quicker climate for iterating
exogen validate
```

## Outputs

Today: equirectangular maps (float EXR and 16-bit PNG) of albedo, height, a shaded preview, a snapshot at the chosen season, the liquid's depth and mask, the year's mean, warmest and coldest temperature, precipitation, aridity, wind, clouds, frost, ice, sea ice, the cryosphere's season and thickness, rivers, lakes, surface age and geological units. A JSON sidecar records every parameter, derived value, seed, legend and unit, and the atmosphere's profile by height (its mean and its warmest and coldest tenths: pressure, temperature, density, vapour, Rayleigh scattering and haze).

Planned: cube maps beside the equirectangular layout for every map (standard, and equi-angular to match the simulation's grid); roughness, specular and normal maps; cloud layers, and the clouds' climatology through the year for engines that raise their own weather; an emissive map for surfaces hot enough to glow; ring textures; and in the sidecar the hazes by height and place, the gases' visible absorption, and the air's glow and aurorae.

## Intended scope

- Every body gravity has rounded, from small icy moons through Earth-likes and magma worlds to hycean worlds, sub-Neptunes and giants up to the deuterium-burning limit, with their rings. One model covers them all: no class of body is a separate branch, and features change by degrees as the inputs change.
- Every physical parameter is derived from the star, orbit, body and atmosphere inputs, and any parameter can be unlocked and set by hand.
- Terrain and albedo come from a time-stepped simulation on a sphere, not from layered noise. Physics decides what is where and why. Noise carries that to a realistic surface at every scale but never decides alone where a feature is or what a surface is made of.
- Reference presets reproduce the character of Mercury, Venus, Earth, the Moon, Mars, Io, Europa, Ganymede, Callisto, Titan, Enceladus, Mimas, Miranda, Triton, Pluto, Charon and Ceres, validated against published scaling laws, crater statistics and measured climates. A tidally locked ocean world (TRAPPIST-1e) is checked against published climate models.
- Out of scope: irregular bodies such as Proteus or Hyperion, whose shapes cannot be derived; 3D terrain meshes; browser deployment.
- Desktop application and command-line tool. Runs locally.
- Stack: Rust, CPU compute first (deterministic), GPU compute later.

## Licence

GPL-3.0-or-later. See [LICENSE](LICENSE).

ExoGen is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. It is distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose. See the GNU General Public License for more details.
