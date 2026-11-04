// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! ExoGen engine.
//!
//! The engine is a library with no UI dependency. The desktop application and the
//! command-line tool are thin clients over it. See `PLAN.md` and `docs/adr/` in the
//! repository for the design.

pub mod climate;
pub mod constants;
pub mod convection;
pub mod crater;
pub mod currents;
pub mod drift;
pub mod export;
pub mod fields;
pub mod flexure;
pub mod fluvial;
pub mod geom;
pub mod grid;
pub mod insolation;
pub mod mantle;
pub mod moisture;
pub mod noise;
pub mod orbital;
pub mod params;
pub mod presets;
pub mod project;
pub mod rivers;
pub mod rng;
pub mod shade;
pub mod spharm;
pub mod spin;
pub mod stages;
pub mod thermal;
pub mod units;
pub mod validate;
pub mod wind;
mod sweeps;

/// Version of the engine, written into project files and sidecars.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
