# ExoGen

Process-based planetary map generator for solid-surface worlds.

ExoGen produces full-planet texture sets (albedo, height, roughness, specular, masks and cloud layers) by simulating the processes that shape real planetary surfaces: impact cratering, tectonics, volcanism, climate, erosion, volatile transport and weathering, scaled to the body's size, composition, thermal state and orbital environment. Outputs are equirectangular maps at up to 16k for use in physically based renderers.

Status: planning. There is no code yet. The design is in [PLAN.md](PLAN.md).

## Intended scope

- Desktop application and command-line tool. Runs locally. No browser.
- Every physical parameter is derived from star, orbit, body and atmosphere inputs, and any parameter can be unlocked and set by hand.
- Terrain and albedo come from a time-stepped simulation on a sphere, not from layered noise. Noise is used only for initial conditions and for fine detail conditioned on the simulated structure.
- Reference presets reproduce the character of Mercury, Venus, Earth, the Moon, Mars, Io, Europa, Ganymede, Callisto, Titan, Enceladus, Triton, Pluto and Ceres, validated against published scaling laws and crater statistics.
- Planned stack: Rust, with CPU compute first and GPU compute later.

## Licence

GPL-3.0-or-later. See [LICENSE](LICENSE).

ExoGen is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. It is distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose. See the GNU General Public License for more details.
