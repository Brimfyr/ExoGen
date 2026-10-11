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
/// Floors and ejecta of craters that excavated below the surface layer into the deep
/// bedrock.
pub const UNIT_DEEP_FLOOR: u16 = 8;
pub const UNIT_DEEP_EJECTA: u16 = 9;
/// Lava plains, and the floors and ejecta of craters formed in them.
pub const UNIT_LAVA: u16 = 10;
pub const UNIT_LAVA_CRATER_FLOOR: u16 = 11;
pub const UNIT_LAVA_EJECTA: u16 = 12;
/// Deformed old crust that survived an episodic lid's resurfacing.
pub const UNIT_TESSERA: u16 = 13;
/// The Io regime: patera floors, fresh silicate flows, plume deposits, mountains.
pub const UNIT_PATERA_FLOOR: u16 = 14;
pub const UNIT_IO_FLOW: u16 = 15;
pub const UNIT_PLUME_RED: u16 = 16;
pub const UNIT_PLUME_WHITE: u16 = 17;
pub const UNIT_MOUNTAIN: u16 = 18;
/// An ice shell's cracks: double ridges and cycloids, and dilational bands.
pub const UNIT_LINEA: u16 = 19;
pub const UNIT_BAND: u16 = 20;
/// Lanes of grooved terrain, and the renewed terrain under a pole where the shell is thin.
pub const UNIT_GROOVED: u16 = 21;
pub const UNIT_POLAR_TERRAIN: u16 = 22;
/// The matrix of a chaos patch, and a cryovolcanic dome.
pub const UNIT_CHAOS: u16 = 23;
pub const UNIT_CRYODOME: u16 = 24;

pub fn unit_name(id: u16) -> &'static str {
    match id {
        UNIT_PRIMORDIAL => "primordial crust",
        UNIT_CONTINENTAL => "continental crust",
        UNIT_OCEANIC => "oceanic crust",
        UNIT_CRATER_FLOOR => "crater floor",
        UNIT_EJECTA => "continuous ejecta",
        UNIT_BASIN_FLOOR => "basin floor",
        UNIT_SECONDARY_FIELD => "secondary crater field",
        UNIT_DEEP_FLOOR => "deep-excavating crater floor",
        UNIT_DEEP_EJECTA => "deep-sourced ejecta",
        UNIT_LAVA => "lava plains",
        UNIT_LAVA_CRATER_FLOOR => "crater floor in lava",
        UNIT_LAVA_EJECTA => "lava ejecta",
        UNIT_TESSERA => "tessera",
        UNIT_PATERA_FLOOR => "patera floor",
        UNIT_IO_FLOW => "fresh silicate flow",
        UNIT_PLUME_RED => "plume deposit, sulfur",
        UNIT_PLUME_WHITE => "plume deposit, frost",
        UNIT_MOUNTAIN => "mountain block",
        UNIT_LINEA => "linea",
        UNIT_BAND => "band",
        UNIT_GROOVED => "grooved terrain",
        UNIT_POLAR_TERRAIN => "polar terrain",
        UNIT_CHAOS => "chaos",
        UNIT_CRYODOME => "cryovolcanic dome",
        _ => "none",
    }
}

/// All known units with their names, for the sidecar legend.
pub fn legend() -> Vec<(u16, &'static str)> {
    (0..=UNIT_CRYODOME).map(|id| (id, unit_name(id))).collect()
}
