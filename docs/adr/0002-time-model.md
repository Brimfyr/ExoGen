# ADR 0002: Time model, coarse forward loop and age-ordered replay

Status: accepted, 2026-10-10.

## Context

The product's defining property is that surfaces are the result of processes acting over time. A naive implementation would step every process at full resolution over the body's whole age, which is too slow at 4k and above. A purely staged implementation (tectonics, then craters, then erosion, each once) loses the interactions that make old surfaces look old: lava flooding earlier craters, scarps cutting craters, rivers breaching rims. The design must keep those interactions while bounding cost, and it must support previewing and exporting the body at any epoch.

## Decision

1. **Coarse forward loop.** Simulated time runs from the body's formation to the present in N steps (default 90 steps of 50 Myr) on a coarse cubesphere (F=256 default, 512 high). Each step samples the history curves, then runs interior, tectonics, volcanism, large impacts, climate (every k steps) and coarse erosion, then writes bookkeeping and a keyframe.
2. **Split rule.** A feature is simulated inside the loop when it spans at least about eight loop cells; smaller features are recorded or sampled for the finishing pass. The loop resolution is therefore a quality dial that moves features between the two populations without changing their statistics.
3. **Age-ordered replay.** The finishing pass, at the full simulation resolution, resamples the coarse state at the output epoch and emplaces every sub-threshold feature oldest first. Each feature is degraded by the process history of its cell since its own age (diffusion, infill, breaching, deformation, relaxation) using the rate fields the loop produced, and warped by the deformation accumulated since its age. Younger features overprint older ones because of the ordering.
4. **Keyframes and epochs.** The loop stores a keyframe per step. An output epoch parameter selects the keyframe the finishing pass starts from. The preview scrubber applies a cheap preview-resolution finishing pass to any keyframe on demand.
5. **History curves.** Stellar luminosity, eccentricity, atmospheric pressure, volatile inventory and impact flux may vary over simulated time as piecewise-linear curves with constant or derived defaults.

## Consequences

- Interactions between large features are real; interactions involving small features are approximated by integrated rates and corrected by a bounded number of full-resolution erosion iterations.
- Memory for keyframes: about 5 MB each at F=256 with key layers at 16-bit precision; 450 MB for a full history in memory; disk-backed above that.
- The loop must record per-step plate assignments and a coarse non-rigid displacement field so the replay can compose deformation cheaply.
- Every feature type that can be replayed needs an analytical degradation model; diffusion is exact (pre-blur with σ² = 2∫κ dt), the others are approximations documented per process.
- Any epoch can be exported at full quality for the cost of one normal render, and sequences of epochs are a batch job.
