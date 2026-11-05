# ADR 0006: The cryosphere's export

Status: accepted, 2026-11-01.

## Context

The ices that come and go through the year (F246) were exported as one RGBA map, `cryosphere.exr` and a 16-bit `cryosphere.png`: R the phase an ice freezes at, G the phase it melts out at, B its thickness at its peak, A the phase of its peak, the ices named in prose in the sidecar. The user's design was an example, open to change while the format is young. A live test of the map in an engine, and a review written during it (2026-11-01, checked against the code and a fresh render), found:

- a loader that premultiplies by alpha damages the phases: on Mars 28.5 percent of frosted pixels have an 8-bit alpha of 26 or less;
- which ice a pixel holds was implied by the liquid mask, sampled bilinearly where the phases are read at the nearest cell, so that the two disagree along coasts;
- the ices' density, the mass that hides the ground and the shape of their year existed only in the source;
- the sidecar's parameters were labelled in display units but written in SI: Mars's year read as 59 million days;
- the frost's windows at its margins copied the last frosted latitude's whole window, so that a band of the cap appeared and vanished at once (F252).

## Decision

1. Three maps of one meaning each, none with an alpha channel:
   - `cryosphere_window`, a float32 EXR with channels R, G and B and a 16-bit RGB PNG: R the phase an ice freezes at, G the phase it melts out at, B the phase of its peak. Phases are shares of the cycle's time after the northern spring equinox, the year or a resonant body's whole repeat. R is 0 and G 1 where an ice lies all year; all three are 0 where none lies. Read at the nearest cell per pixel; never filtered.
   - `cryosphere_thickness`, a float32 EXR with the channel `thickness` and a 16-bit grey PNG from 0 to the sidecar's `thickest_m`: the ice's thickness at its peak, m. It is read bilinearly, but only where the nearest cell's window lies.
   - `cryosphere_ice`, an 8-bit grey PNG: which ice lies, its index among the sidecar's `cryosphere.ices` from 1, 0 none. Read at the nearest cell; never filtered.
2. The sidecar's `cryosphere` block holds:
   - the phase origin and how to read the maps (`evaluation`);
   - the snapshot's phase, the moment its maps (`frost`, `sea_ice`) were taken, the end of the step that holds its season (F284), and `thickest_m`;
   - for each ice its index, name, species, where it lies (land or sea), colour (linear reflectance), density, its cover mass (the mass that hides 1 − 1/e of what lies under it), its lying mass (where half is hidden, which opens and closes its window) and its year curve's exponents, `rise` and `fall`;
   - the clock: the cycle's length, its orbits and steps, the solar day, the obliquity, the eccentricity, the star's longitude at perihelion and the phases of the equinoxes and solstices;
   - the year step by step: each step's phase, the star's longitude, declination and flux, the surface pressure, and each ice's mean over each hemisphere.
3. An engine reads the maps per texel at their resolution and filters the result. An ice lies at phase p where p falls in its window, from R round to G. Its mass is its thickness at its peak times its density, times (t / rise)^a from R to its peak's phase B, t the time since R and rise the time from R to B, then times 1 − (g / fall)^b from B to G, g the time since B and fall the time from B to G, (a, b) its curve; all year its peak where R is 0 and G 1. It hides 1 − exp(−m / cover mass) of what lies under it and stands m / density thick. `climate::Window::at_phase` is the reference; `docs/cryosphere-vectors.json` gives readings to check a port against.
4. The phases keep one origin, the northern spring equinox, for every ice and hemisphere. Filtering the window is the mistake; an origin per hemisphere would soften its damage but not cure it, and baked frames (Phase 7's `exogen bake`) serve the pipelines that must filter.
5. Every map's record in the sidecar gives its `filter`, `linear` or `nearest`, and where it has several channels each one's meaning, unit and range. The parameters are written in their display units with the `scale` to SI. The feature lists go to `features.json`.

## Consequences

- An engine needs the three maps and the sidecar, and no constant from the source.
- `cryosphere.exr` and `cryosphere.png` are gone; tools written for them move to the three maps.
- The year curve is a fit, one shape per ice over its places; the year's series gives the hemispheres' exact means.
- Perennial ice is drawn at its peak all year: its seasonal thinning is not in the maps.
- A pipeline that filters textures needs baked frames, which wait for Phase 7.
