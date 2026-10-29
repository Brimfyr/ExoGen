// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! Impact cratering: production and chronology, morphology, and stamping (PLAN.md 4.5).

pub mod curtain;
pub mod morphology;
pub mod production;
pub mod stamp;

pub use curtain::{Curtain, Rays};
pub use morphology::{
    CraterClass, Dimensions, Levelling, ProfileTable, Scaling, dimensions, ellipticity_for,
    impact_angle,
};
pub use production::{
    Bombardment, BombardmentParams, CraterSample, Early, FluxHistory, ImpactorCraters,
    ImpactorSizes, PopulationConfig, ProductionFunction, Source, Target, hold_diameter,
    sample_population,
};
pub use stamp::{
    CellWrite, CraterTexture, DeepShare, EdgeShape, RestoreParams, StampOutput, StampParams,
    Stamper, burial_factor, crenulation, ejecta_coverage, ejecta_modulation, excavated_below,
    restore, restore_compute, shed_area, stamp, stamp_compute, with_stamper,
};

/// A crater that was emplaced, kept for the sidecar, for rays and secondaries, and for
/// re-stamping resolved detail at the finishing resolution.
#[derive(Clone, Debug)]
pub struct CraterRecord {
    pub center: [f64; 3],
    pub diameter: f64,
    pub time_gyr: f64,
    /// Downrange direction, radians in the local east-north frame.
    pub azimuth: f64,
    /// Impact angle from horizontal, radians.
    pub angle: f64,
    pub ellipticity: f64,
    pub class: CraterClass,
    /// True for craters emplaced inside the coarse loop.
    pub from_loop: bool,
    /// Rim-to-floor relief left on the coarse surface at the end of the loop, metres;
    /// the fresh relief for replayed craters.
    pub relief: f32,
    /// Whether the crater is still recognisable: residual relief above a fraction of the
    /// fresh relief.
    pub visible: bool,
    /// Formed in lava plains, so its floor and ejecta are made of them.
    pub on_lava: bool,
    pub id: u64,
}
