# ExoGen

Process-based planetary map generator: full-planet texture sets for any body gravity has rounded, made by simulating the processes that shape real worlds.

ExoGen produces texture sets (albedo, height, materials, masks, climate and ice layers) by simulating the processes that shape real planetary surfaces: impact cratering, tectonics, volcanism, climate, erosion, volatile transport and weathering, scaled to the body's size, composition, thermal state, atmosphere and orbital environment. The maps are for physically based renderers, game engines and studio pipelines.

## Status

In development. The design and roadmap are in [PLAN.md](PLAN.md).

- **Built:** airless cratered worlds; interiors, volcanism and stagnant or episodic lids; ice worlds; the climate on every body (energy balance, seasons, days and nights, the star's history); frost, caps and the atmosphere's loss; the wind, moisture, snow, the seas' currents and sea ice; rivers and lakes.
- **Under way:** water and ice at work (glaciers and slope wasting still to come), then dust and sand, surface physics, Io, and weather and clouds.
- **Planned:** plate tectonics; equilibrium figures, rings and multiple light sources; magma worlds; deep atmospheres from hycean worlds to gas giants; detail synthesis to 16k and 32k; vegetation.

A command-line tool renders the map set for a project or preset. A desktop app shows the parameter graph, a preview globe and the pipeline's state.

```
cargo build --release -p exogen-cli    # the binary is target/release/exogen
exogen preset list
exogen derive --preset Titan
exogen render --preset Mars --width 2048 --out out/Mars
exogen validate
```

## Outputs

Today: equirectangular maps (float EXR and 16-bit PNG) of albedo, height, a shaded preview, a snapshot at the chosen season, the liquid's depth and mask, the year's mean, warmest and coldest temperature, precipitation, aridity, wind, frost, ice, sea ice, the cryosphere's season and thickness, rivers, lakes, surface age and geological units. A JSON sidecar records every parameter, derived value, seed, legend and unit.

Planned: cube maps beside the equirectangular layout for every map (standard, and equi-angular to match the simulation's grid); roughness, specular and normal maps; cloud layers; an emissive map for surfaces hot enough to glow; ring textures; and an atmosphere profile by height in the sidecar.

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
