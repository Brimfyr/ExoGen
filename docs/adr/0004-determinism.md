# ADR 0004: Determinism and random numbers

Status: accepted, 2026-10-10.

## Context

A studio that renders a planet must be able to regenerate the same maps months later, possibly on another machine, at a higher resolution, or for a different epoch of the same world. Seeds alone do not give that: thread scheduling, hash-map iteration order, platform math libraries and GPU reduction order all change results. The product also promises that changing the crater seed does not change the continents.

## Decision

1. **The CPU path is the reference and is bit-exact** for a given software version across runs, thread counts and supported platforms. The GPU path is visually identical but not guaranteed bit-exact; the sidecar records which path produced each map.
2. **Random numbers are counter-based, not sequential.** Every random draw is a pure function of (master seed, stage id, feature or cell index, draw index) through a counter-based generator (Philox or ChaCha-based). No stage carries a stream whose position depends on execution order. This makes results independent of threading and allows per-tile detail synthesis to be reproduced in isolation.
3. **Per-stage sub-seeds** are derived from the master seed by hashing the stage id. Any sub-seed can be overridden by the user, matching the parameter lock model.
4. **Resolution independence.** Statistical populations (craters, vents, dunes) are sampled in physical units on the sphere from the seed, then rasterised at whatever resolution is in use. Raising the resolution adds smaller members to a population; it never changes the existing ones.
5. **No order-dependent containers** in simulation code: ordered maps or sorted vectors instead of hash maps, and reductions in a fixed order. Parallel reductions use a deterministic tree.
6. **Platform math.** Transcendental functions in simulation kernels come from a pure-Rust math library (the `libm` crate) rather than the platform C library, so Windows, Linux and macOS agree. Fused multiply-add is used only where written explicitly with `mul_add`.
7. **Project files store everything needed to reproduce**: parameters with lock states, multipliers, seeds, history curves, output settings and the software version. The CLI refuses to render a project from a newer version than itself.

## Consequences

- Determinism is tested: the suite renders reference presets twice with different thread counts and compares byte for byte, and a cross-platform CI job compares hashes.
- Counter-based generation costs a little more per draw than a sequential generator; this is negligible against the physics.
- GPU acceleration is an optimisation that can be turned off, never a requirement for reproducibility.
- Any change to a stage's behaviour changes its code version (ADR 0003), so a given project file plus version string identifies an output exactly.
