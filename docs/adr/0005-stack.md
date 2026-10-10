# ADR 0005: Language, libraries and workspace layout

Status: accepted, 2026-10-10.

## Context

The project needs a single shippable binary per platform, multi-threaded CPU numerics now and portable GPU compute later, a dense technical desktop UI, and a research-heavy development process in which models are tuned and validated against published data. The candidates were Rust, C++ and Python. Python is the fastest to iterate on but the wrong production choice for packaging and 8k performance. C++ is viable with no capability advantage here and more cost in build system, memory safety and packaging.

## Decision

1. **Rust**, stable toolchain, with the project licensed GPL-3.0-or-later. Every dependency must be GPL-3.0-compatible; MIT and Apache-2.0 crates qualify.
2. **Workspace layout**:
   - `exogen-core`: the engine as a library. Parameter graph, grids, fields, stages, materials, export, validation. No UI dependency.
   - `exogen-cli`: headless rendering, derivation reports, validation runs, preset management, sequences.
   - `exogen-app`: desktop UI on `eframe`/`egui` with `wgpu`, with `egui_dock` for panels and `egui_plot` for diagnostics.
   - `exogen-py` (optional, later): thin PyO3 bindings for analysis notebooks. No production code in Python.
3. **Compute**: CPU first through `rayon`. GPU kernels through `wgpu` compute shaders (WGSL) for the hot stages once the CPU reference exists. CPU fallback is mandatory for every GPU kernel.
4. **Libraries**: `image` for PNG and TIFF, `exr` for OpenEXR, `serde` with `toml` for project files and `serde_json` for sidecars, `libm` for portable math in kernels, `glam` or `nalgebra` for small linear algebra, `rand_chacha` or a Philox implementation for counter-based randomness.
5. **Research scripts** live in `research/` and use Python 3.10 with numpy, scipy and matplotlib in a local virtual environment. They validate equations and produce comparison plots; they are not part of the build.
6. **Quality gates**: `cargo fmt`, `cargo clippy` with warnings denied in CI, unit tests per module, the validation suite behind `exogen validate`, and GitHub Actions on Windows, Linux and macOS once there is a build to run.
7. **Platform order**: Windows first (development machine), then macOS and Linux, which the stack supports without architectural change.

## Consequences

- One engine, two front ends, plus optional bindings. No logic lives in the UI.
- Rust compile times slow iteration compared with Python; the research scripts absorb the equation-tuning iteration, and the engine is built once a model is settled.
- `egui` is less polished than Qt for complex layout. The UI is deliberately dense and simple, which `egui` handles well.
- WGSL kernels must be written to match the CPU reference within visual tolerance and are tested against it.
