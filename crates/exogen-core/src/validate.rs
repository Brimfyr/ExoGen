// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! Validation against the reference bodies (PLAN.md section 6).
//!
//! Each check derives a quantity for a preset and compares it with the measured value and
//! a tolerance. The CLI prints the table; the test suite requires every check to pass, so a
//! change to a derivation that breaks a known body is caught immediately.

use crate::params::ParamGraph;
use crate::presets;

#[derive(Clone, Debug)]
pub enum Expected {
    /// Numeric value within [lo, hi] in stored units.
    Range(f64, f64),
    /// Exact string (enum) value.
    Str(&'static str),
    Bool(bool),
}

#[derive(Clone, Debug)]
pub struct Check {
    pub preset: &'static str,
    pub node: &'static str,
    pub expected: Expected,
    /// What the measured value is and where it comes from.
    pub note: &'static str,
}

#[derive(Clone, Debug)]
pub struct CheckResult {
    pub check: Check,
    pub actual: String,
    pub pass: bool,
}

fn range(preset: &'static str, node: &'static str, lo: f64, hi: f64, note: &'static str) -> Check {
    Check {
        preset,
        node,
        expected: Expected::Range(lo, hi),
        note,
    }
}

fn text(preset: &'static str, node: &'static str, s: &'static str, note: &'static str) -> Check {
    Check {
        preset,
        node,
        expected: Expected::Str(s),
        note,
    }
}

fn flag(preset: &'static str, node: &'static str, b: bool, note: &'static str) -> Check {
    Check {
        preset,
        node,
        expected: Expected::Bool(b),
        note,
    }
}

/// The reference table. Tolerances are stated in the note where they are not obvious.
pub fn reference_checks() -> Vec<Check> {
    vec![
        // Earth
        range(
            "Earth",
            "body.gravity",
            9.72,
            9.92,
            "9.82 m/s² from mass and mean radius, 1%",
        ),
        range(
            "Earth",
            "body.escape_velocity",
            11_080.0,
            11_300.0,
            "11.19 km/s, 1%",
        ),
        text("Earth", "body.composition_class", "rock", ""),
        text("Earth", "rotation.mode", "free", ""),
        range(
            "Earth",
            "atmosphere.equilibrium_temperature",
            252.0,
            257.0,
            "254.6 K at Bond albedo 0.30",
        ),
        range(
            "Earth",
            "atmosphere.surface_temperature_mean",
            284.0,
            292.0,
            "288 K, 4 K",
        ),
        range(
            "Earth",
            "interior.surface_heat_flux",
            0.065,
            0.110,
            "87 mW/m², 25%",
        ),
        text("Earth", "interior.tectonic_regime", "mobile_lid", ""),
        text("Earth", "volatiles.working_liquid", "water", ""),
        range(
            "Earth",
            "climate.hadley_edge_latitude",
            26.0,
            38.0,
            "about 30°",
        ),
        range(
            "Earth",
            "interior.max_relief",
            8_000.0,
            11_000.0,
            "Everest 8.8 km above sea level, 9.5 km design value",
        ),
        range(
            "Earth",
            "cratering.impact_velocity",
            15_000.0,
            24_000.0,
            "mean about 19 to 20 km/s",
        ),
        range(
            "Earth",
            "cratering.retention_age",
            0.05,
            0.2,
            "mobile lid keeps craters for about 100 Myr; Phase 2 rule",
        ),
        text(
            "Earth",
            "materials.bedrock",
            "felsic",
            "continental default; oceanic crust uses basalt by unit",
        ),
        range(
            "Moon",
            "cratering.retention_age",
            4.0,
            5.0,
            "airless stagnant lid keeps craters for its whole age",
        ),
        text("Moon", "materials.bedrock", "anorthosite", ""),
        text("Moon", "materials.deep_bedrock", "norite", ""),
        range(
            "Moon",
            "surface.solidification_time",
            0.25,
            0.35,
            "record from 4.2 Gyr",
        ),
        range("Mercury", "cratering.retention_age", 4.0, 5.0, ""),
        range(
            "Mercury",
            "interior.contraction_total",
            5_000.0,
            7_000.0,
            "Byrne et al. 2014",
        ),
        range(
            "Moon",
            "interior.contraction_total",
            100.0,
            600.0,
            "Watters et al. 2010",
        ),
        range(
            "Mars",
            "interior.contraction_total",
            1_000.0,
            4_000.0,
            "Nahm and Schultz 2011",
        ),
        text("Mercury", "materials.bedrock", "mercurian_regolith", ""),
        text(
            "Mercury",
            "materials.deep_bedrock",
            "low_reflectance_material",
            "",
        ),
        range(
            "Venus",
            "cratering.retention_age",
            0.3,
            0.7,
            "global resurfacing 300 to 700 Myr ago",
        ),
        text("Venus", "materials.bedrock", "basalt", ""),
        range(
            "Mars",
            "cratering.retention_age",
            3.0,
            4.0,
            "thin atmosphere; highlands retain craters to about 3.5 Gyr",
        ),
        text("Mars", "materials.bedrock", "basalt", ""),
        range(
            "Io",
            "cratering.retention_age",
            0.0,
            0.01,
            "resurfaced faster than craters form",
        ),
        text("Io", "materials.bedrock", "sulfur_coated", ""),
        range(
            "Europa",
            "cratering.retention_age",
            0.05,
            0.2,
            "surface age 40 to 90 Myr",
        ),
        text("Europa", "materials.bedrock", "water_ice", ""),
        range(
            "Ganymede",
            "cratering.retention_age",
            4.0,
            5.0,
            "dark terrain is ancient; grooved terrain is Phase 4",
        ),
        range(
            "Callisto",
            "cratering.retention_age",
            4.0,
            5.0,
            "the most heavily cratered surface known",
        ),
        range(
            "Enceladus",
            "cratering.retention_age",
            3.0,
            5.0,
            "the cratered plains are ancient; the south polar terrain renews itself on its own",
        ),
        // Mars
        range("Mars", "body.gravity", 3.64, 3.78, "3.71 m/s²"),
        range(
            "Mars",
            "body.escape_velocity",
            4_950.0,
            5_100.0,
            "5.03 km/s",
        ),
        text("Mars", "body.composition_class", "rock", ""),
        text("Mars", "rotation.mode", "free", ""),
        range(
            "Mars",
            "atmosphere.surface_temperature_mean",
            207.0,
            219.0,
            "213 K, 6 K",
        ),
        range(
            "Mars",
            "interior.surface_heat_flux",
            0.015,
            0.045,
            "20 to 30 mW/m² (InSight)",
        ),
        text("Mars", "interior.tectonic_regime", "stagnant_lid", ""),
        text(
            "Mars",
            "volatiles.working_liquid",
            "none",
            "water is ice at 213 K and 6 mbar",
        ),
        range(
            "Mars",
            "interior.max_relief",
            20_000.0,
            30_000.0,
            "Olympus Mons 22 to 26 km",
        ),
        range(
            "Mars",
            "cratering.simple_complex_transition",
            5_000.0,
            9_000.0,
            "6 to 7 km observed, 7.2 km from 1/g",
        ),
        // Moon
        range("Moon", "body.gravity", 1.59, 1.65, "1.62 m/s²"),
        text("Moon", "rotation.mode", "locked", ""),
        text("Moon", "body.composition_class", "rock", ""),
        range(
            "Moon",
            "interior.surface_heat_flux",
            0.009,
            0.021,
            "12 to 16 mW/m² (Apollo)",
        ),
        range(
            "Moon",
            "atmosphere.surface_temperature_mean",
            262.0,
            278.0,
            "equilibrium about 270 K at Bond albedo 0.11",
        ),
        range(
            "Moon",
            "cratering.simple_complex_transition",
            15_000.0,
            18_000.0,
            "15 to 18 km",
        ),
        // Mercury
        range("Mercury", "body.gravity", 3.63, 3.77, "3.70 m/s²"),
        text(
            "Mercury",
            "body.composition_class",
            "iron_rich",
            "core about 70% of the mass",
        ),
        range(
            "Mercury",
            "atmosphere.surface_temperature_mean",
            425.0,
            450.0,
            "equilibrium about 437 K at Bond albedo 0.088",
        ),
        text("Mercury", "interior.tectonic_regime", "stagnant_lid", ""),
        range(
            "Mercury",
            "interior.contraction_total",
            4_000.0,
            8_000.0,
            "5 to 7 km of radial contraction",
        ),
        // Venus
        range("Venus", "body.gravity", 8.78, 8.96, "8.87 m/s²"),
        range(
            "Venus",
            "atmosphere.equilibrium_temperature",
            225.0,
            234.0,
            "229 K at Bond albedo 0.76",
        ),
        range(
            "Venus",
            "atmosphere.surface_temperature_mean",
            715.0,
            760.0,
            "737 K, 3%",
        ),
        text("Venus", "body.composition_class", "rock", ""),
        text("Venus", "interior.tectonic_regime", "episodic_lid", ""),
        range(
            "Venus",
            "interior.max_relief",
            9_000.0,
            12_500.0,
            "Maxwell Montes 11 km",
        ),
        range(
            "Venus",
            "atmosphere.crater_cutoff_diameter",
            1_500.0,
            3_200.0,
            "no craters below 1.5 to 3 km",
        ),
        // Io
        range("Io", "body.gravity", 1.76, 1.84, "1.80 m/s²"),
        text("Io", "rotation.mode", "locked", ""),
        text("Io", "body.composition_class", "rock", ""),
        range(
            "Io",
            "interior.tidal_flux",
            1.0,
            4.0,
            "about 2 W/m², factor of two",
        ),
        range(
            "Io",
            "interior.magma_production_rate",
            5_000.0,
            50_000.0,
            "500 to 1500 km³/yr is 16,000 to 48,000 m³/s",
        ),
        // Europa
        range("Europa", "body.gravity", 1.29, 1.34, "1.31 m/s²"),
        text("Europa", "rotation.mode", "locked", ""),
        text("Europa", "body.composition_class", "rock_ice", ""),
        range(
            "Europa",
            "interior.tidal_flux",
            0.03,
            0.30,
            "0.05 to 0.3 W/m² estimates",
        ),
        flag("Europa", "interior.has_subsurface_ocean", true, ""),
        text("Europa", "interior.tectonic_regime", "ice_shell_ocean", ""),
        range(
            "Europa",
            "body.ice_layer_thickness",
            60_000.0,
            250_000.0,
            "80 to 170 km H2O layer; the provisional mixing model gives 220 km",
        ),
        // Ganymede
        range("Ganymede", "body.gravity", 1.40, 1.46, "1.43 m/s²"),
        text("Ganymede", "rotation.mode", "locked", ""),
        text("Ganymede", "body.composition_class", "rock_ice", ""),
        flag("Ganymede", "interior.has_subsurface_ocean", true, ""),
        range(
            "Ganymede",
            "cratering.simple_complex_transition",
            2_000.0,
            3_500.0,
            "2 to 3 km on ice",
        ),
        // Callisto
        range("Callisto", "body.gravity", 1.21, 1.26, "1.24 m/s²"),
        text("Callisto", "rotation.mode", "locked", ""),
        text("Callisto", "body.composition_class", "rock_ice", ""),
        // Titan
        range("Titan", "body.gravity", 1.33, 1.38, "1.35 m/s²"),
        text("Titan", "rotation.mode", "locked", ""),
        text("Titan", "body.composition_class", "rock_ice", ""),
        range(
            "Titan",
            "atmosphere.surface_temperature_mean",
            90.0,
            98.0,
            "94 K, 4 K",
        ),
        text("Titan", "volatiles.working_liquid", "methane_ethane", ""),
        range(
            "Titan",
            "climate.hadley_edge_latitude",
            80.0,
            90.0,
            "single cell to the pole",
        ),
        // Enceladus
        range("Enceladus", "body.gravity", 0.108, 0.118, "0.113 m/s²"),
        text("Enceladus", "rotation.mode", "locked", ""),
        text("Enceladus", "body.composition_class", "ice_dominated", ""),
        range(
            "Enceladus",
            "interior.tidal_power",
            5.0e9,
            4.0e10,
            "16 GW observed, factor of two",
        ),
        range(
            "Enceladus",
            "tectonics.polar_thinning",
            0.5,
            0.8,
            "a shell a few kilometres thick under the south pole against twenty elsewhere",
        ),
        // Triton
        range("Triton", "body.gravity", 0.76, 0.80, "0.78 m/s²"),
        text("Triton", "rotation.mode", "locked", ""),
        range(
            "Triton",
            "atmosphere.surface_temperature_mean",
            33.0,
            42.0,
            "38 K",
        ),
        // Pluto
        range("Pluto", "body.gravity", 0.60, 0.64, "0.62 m/s²"),
        text(
            "Pluto",
            "rotation.mode",
            "free",
            "locked to Charon in reality, which is not modelled",
        ),
        text("Pluto", "body.composition_class", "rock_ice", ""),
        range(
            "Pluto",
            "atmosphere.surface_temperature_mean",
            34.0,
            48.0,
            "about 40 K at Bond albedo 0.5",
        ),
        // Ceres
        range("Ceres", "body.gravity", 0.27, 0.29, "0.28 m/s²"),
        text("Ceres", "rotation.mode", "free", ""),
        text("Ceres", "body.composition_class", "rock_ice", ""),
        range(
            "Ceres",
            "atmosphere.surface_temperature_mean",
            155.0,
            175.0,
            "about 168 K",
        ),
    ]
}

/// Run every reference check.
pub fn run_reference_checks() -> Vec<CheckResult> {
    let mut graphs: std::collections::BTreeMap<&'static str, ParamGraph> = Default::default();
    let mut out = Vec::new();
    let mut checks = reference_checks();
    checks.sort_by_key(|c| {
        presets::PRESET_NAMES
            .iter()
            .position(|n| *n == c.preset)
            .unwrap_or(usize::MAX)
    });
    for check in checks {
        let g = graphs.entry(check.preset).or_insert_with(|| {
            let p = presets::preset(check.preset)
                .unwrap_or_else(|| panic!("no preset {}", check.preset));
            p.to_graph().0
        });
        let (actual, pass) = match &check.expected {
            Expected::Range(lo, hi) => {
                let v = g.f64(check.node);
                (format_number(v), v >= *lo && v <= *hi)
            }
            Expected::Str(s) => {
                let v = g.str(check.node).to_string();
                let pass = v == *s;
                (v, pass)
            }
            Expected::Bool(b) => {
                let v = g.bool(check.node);
                (v.to_string(), v == *b)
            }
        };
        out.push(CheckResult {
            check,
            actual,
            pass,
        });
    }
    out
}

/// Compact number formatting for reports.
pub fn format_number(v: f64) -> String {
    if v == 0.0 {
        return "0".to_string();
    }
    let a = v.abs();
    if !(1e-3..1e6).contains(&a) {
        format!("{v:.3e}")
    } else if a >= 100.0 {
        format!("{v:.1}")
    } else if a >= 1.0 {
        format!("{v:.3}")
    } else {
        format!("{v:.4}")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_reference_check_passes() {
        let results = run_reference_checks();
        let failures: Vec<String> = results
            .iter()
            .filter(|r| !r.pass)
            .map(|r| {
                format!(
                    "{} {} = {} (expected {:?}; {})",
                    r.check.preset, r.check.node, r.actual, r.check.expected, r.check.note
                )
            })
            .collect();
        assert!(
            failures.is_empty(),
            "{} reference checks failed:\n{}",
            failures.len(),
            failures.join("\n")
        );
    }
}
