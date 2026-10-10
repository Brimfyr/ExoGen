# ADR 0001: Simulation grid

Status: accepted, 2026-10-10.

## Context

Every process in ExoGen acts on fields defined over a sphere. The choice of grid decides whether erosion, flow routing and cratering have artefacts at the poles, whether the finite-volume weather model is practical, how fields map to GPU textures, and how the final equirectangular export is produced. Candidates were the equirectangular grid, the cubesphere, the icosahedral grid and HEALPix.

Equirectangular has a singularity and a strong cell-area gradient toward the poles, which breaks flow routing and biases crater statistics. The icosahedral grid is the most uniform but has no natural 2D texture layout and needs index tables for every neighbour lookup. HEALPix is equal-area but awkward for stencil operations. The cubesphere has mild area variation, direct 2D neighbour access within a face, six 2D textures on the GPU, and a long history in finite-volume atmospheric dynamics.

## Decision

1. All gridded fields live on an **equiangular cubesphere**: six faces of F x F cells, with the equiangular mapping so that the cell-area ratio between face centre and corner is about 1.3. The same grid code serves every resolution: the coarse forward loop (F=256 to 512), the climate and weather models (F=64 to 256), the finishing pass (F=1024 to 2048) and the detail pass (tiles of an F=4096 to 8192 virtual grid).
2. Neighbour access across face edges and corners goes through precomputed tables. Stencil operators (gradient, Laplacian, flow direction) are written once against those tables.
3. **Plate tectonics alone uses a point representation**: a Fibonacci sphere of 10k to 100k points with its spherical Voronoi tessellation, because plates are regions whose material is advected, and a point set handles that without remeshing. Plate state is rasterised to the cubesphere by spherical interpolation at each loop step.
4. **Equirectangular is an export projection only**. Resampling to it happens at the end with area-aware filtering. Gradients for normals and slopes are computed on the cubesphere before projection, so polar pixels are correct.

## Consequences

- No pole handling anywhere in the process code.
- One grid implementation to test, with resolution as a parameter.
- Cell size varies by about 15% across a face. Processes whose rates depend on cell size (diffusion, flow accumulation) use the true cell area and edge lengths, never a constant.
- Face boundaries are a correctness risk for stencils and flow routing. The grid module ships with tests that run every operator on a field with no visible face seams (for example a spherical harmonic) and assert seam error below a tolerance.
- The export step is a resampling, so export resolution is independent of simulation resolution.
