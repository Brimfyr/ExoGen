# ADR 0003: Stage graph, caching and invalidation

Status: accepted, 2026-10-10.

## Context

A full run takes minutes. Users will change a material colour, a cloud parameter or the output epoch far more often than they change the star or the plate count. Without caching, every change re-runs tectonics. The UI also needs to show which stages are current, which are stale and what re-running will cost.

## Decision

1. The pipeline is an explicit **directed acyclic graph of stages**. Each stage declares the parameter nodes, seeds and upstream stage outputs it reads, and the fields and feature lists it writes. Stages cannot read anything they did not declare; the runner enforces this.
2. Each stage output is keyed by a **content hash** of: the declared parameter values after derivation, the declared seeds, the upstream output hashes, the stage's code version string, and the resolution. Equal key means the cached output is reused.
3. Cache entries live in memory up to a budget and spill to a per-project cache directory on disk. Keyframes of the coarse loop are cache entries like any other output.
4. Changing a parameter marks exactly the stages whose keys change as stale. The UI colours stale stages, shows the last run time of each, and offers "run to here" and "invalidate from here".
5. The runner executes stale stages in topological order, in a worker thread, reporting progress per stage. A run can be cancelled between stages, and between loop steps inside the coarse loop.
6. Stage code versions are bumped by hand when a stage's behaviour changes, so old caches are not reused across behaviour changes. The test suite includes a check that every stage has a version string and that it changed when the stage's source changed in a release.

## Consequences

- Late-stage iteration (materials, clouds, detail, export) is interactive after the first run.
- Declared inputs make stages unit-testable in isolation with synthetic upstream outputs.
- The hash must be stable across runs and machines: parameters are hashed in a canonical serialisation, floats as their exact bit pattern.
- Disk cache can grow large; the project settings hold a size limit with least-recently-used eviction.
- Headless CLI runs use the same runner and cache, so a render farm job that changes one parameter between frames pays only for what changed.
