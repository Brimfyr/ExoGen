// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! Impact cratering: production and chronology, morphology, and stamping (PLAN.md 4.5).

pub mod morphology;
pub mod production;
pub mod stamp;

pub use morphology::{CraterClass, Dimensions, ProfileTable, Scaling, dimensions};
pub use production::{CraterSample, PopulationConfig, ProductionFunction, sample_population};
pub use stamp::{StampParams, Stamper, stamp};

/// A crater that was emplaced, kept for the sidecar, for rays and secondaries, and for
/// re-stamping resolved detail at the finishing resolution.
#[derive(Clone, Debug)]
pub struct CraterRecord {
    pub center: [f64; 3],
    pub diameter: f64,
    pub time_gyr: f64,
    pub azimuth: f64,
    pub ellipticity: f64,
    pub class: CraterClass,
    /// True for craters emplaced inside the coarse loop.
    pub from_loop: bool,
    pub id: u64,
}
