// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! ExoGen engine.
//!
//! The engine is a library with no UI dependency. The desktop application and the
//! command-line tool are thin clients over it. See `PLAN.md` and `docs/adr/` in the
//! repository for the design.

pub mod constants;
pub mod crater;
pub mod export;
pub mod fields;
pub mod geom;
pub mod grid;
pub mod mantle;
pub mod noise;
pub mod orbital;
pub mod params;
pub mod presets;
pub mod project;
pub mod rng;
pub mod shade;
pub mod spharm;
pub mod stages;
pub mod thermal;
pub mod units;
pub mod validate;

/// Version of the engine, written into project files and sidecars.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
