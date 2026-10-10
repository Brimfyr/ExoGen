// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! Geological unit labels written into the `unit_id` layer and exported with a legend.

pub const UNIT_NONE: u16 = 0;
pub const UNIT_PRIMORDIAL: u16 = 1;
pub const UNIT_CONTINENTAL: u16 = 2;
pub const UNIT_OCEANIC: u16 = 3;
pub const UNIT_CRATER_FLOOR: u16 = 4;
pub const UNIT_EJECTA: u16 = 5;
pub const UNIT_BASIN_FLOOR: u16 = 6;
pub const UNIT_SECONDARY_FIELD: u16 = 7;

pub fn unit_name(id: u16) -> &'static str {
    match id {
        UNIT_PRIMORDIAL => "primordial crust",
        UNIT_CONTINENTAL => "continental crust",
        UNIT_OCEANIC => "oceanic crust",
        UNIT_CRATER_FLOOR => "crater floor",
        UNIT_EJECTA => "continuous ejecta",
        UNIT_BASIN_FLOOR => "basin floor",
        UNIT_SECONDARY_FIELD => "secondary crater field",
        _ => "none",
    }
}

/// All known units with their names, for the sidecar legend.
pub fn legend() -> Vec<(u16, &'static str)> {
    (0..=UNIT_SECONDARY_FIELD)
        .map(|id| (id, unit_name(id)))
        .collect()
}
