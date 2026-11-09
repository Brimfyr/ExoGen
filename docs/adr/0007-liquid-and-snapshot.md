# ADR 0007: The liquid apart from the albedo, and the snapshot

Status: accepted, 2026-11-09.

## Context

The user looked at the full map sets of Earth, Venus, Mars and TRAPPIST-1e (2026-11-09) and found:

- the shaded preview showed a uniform noise texture over the seas, where a liquid's surface is flat: the preview lit the ground's height, and under a sea that is the seabed's relief;
- the albedo drew a cell under a liquid in the liquid's colour alone, a tint from deep to shallow by its depth: the ground's own colour under it was never computed, so that the ground under a beach turned the shallow tint's blue, the seabed's fine relief showed in the coastal water's colour, and an engine could not set its own sea level over the ground;
- nothing drew the ices that come and go over the albedo, so the albedo maps looked bare of the season's ice (ADR 0006 leaves that composite to an engine, and baked frames to Phase 7).

The user chose to separate the liquid from the albedo, as the ices already are, so that a designer can use the dry albedo with a sea level set in the engine, or composite the two.

## Decision

1. `albedo` is the ground's own colour everywhere, under a liquid too. No process under a liquid is modelled yet: a seabed shows the colour its unit's ground has, the oceanic crust its basalt. The perennial caps stay in it as ground (F241).
2. `liquid`, a float32 EXR with the channel `depth` and a 16-bit grey PNG from 0 to the deepest: the liquid's depth over the ground, m, seas and lakes alike, its surface flat; 0 where none lies. The sidecar's `liquid` block holds the working liquid, the sea's surface above the datum (none where only lakes lie), its deep reflectance (linear RGB), its light's diffuse attenuation in each channel, 1/m, the channels' wavelengths, and how to compose it: where the depth H is above 0, ground of albedo A shows through the liquid as R = R_deep + (A - R_deep) exp(-2 K H) in each channel (Maritorena, Morel and Gentili 1994). Water's attenuation is pure water's absorption at the sRGB primaries' dominant wavelengths, 612, 549 and 464 nm (Pope and Fry 1997); the deep reflectances are those the albedo drew the liquids with; the other liquids' attenuation is stated, grey, until measured.
3. `snapshot`, a float16 EXR RGB and an 8-bit sRGB PNG: the surface as it looks at the snapshot's season. The albedo under the liquid by that model, then the ice that lies on each cell at the snapshot's phase as the cryosphere's maps draw it (`climate::Window::at_phase`, ADR 0006), hiding 1 - exp(-m / cover mass) of what lies under it. It is composed on the model's grid and read between its cells as the albedo is: evaluated, then filtered. It is one baked frame; Phase 7's bake gives any phase.
4. `shaded` is the snapshot, exposed to a common mean, lit over the surface: the ground's height where it is dry, the liquid's flat surface where a liquid lies. Its exaggeration is still read from the ground's slopes, so that a sea's flat surface does not steepen the land's.
5. `mask_liquid` stays: where the depth is above 0.

## Consequences

- A pipeline that wants the old look composes `albedo` and `liquid` by the sidecar's model, at ExoGen's sea level or its own; `snapshot` is that composite at ExoGen's level, with the snapshot's ices over it.
- Under a liquid the ground is the model's ground, not a sediment: seabeds stay basalt until processes under the liquid (sediment, ooze) are modelled.
- The desktop app's views follow: the snapshot (shown first), the albedo (the ground's own), and the shaded relief (the snapshot lit over its surface).
- Liquids other than water carry stated optics until measured ones replace them.
