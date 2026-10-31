// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! Plausibility rules (docs/parameters.md). Rules warn and explain; they never block.

use super::ParamGraph;
use super::derive::physics::volatile;
use super::derive::{PHASE_LIQUID, PHASE_LIQUID_AND_SOLID};
use crate::constants::{AU, JEANS_RETENTION_THRESHOLD, THERMAL_TIDE_PRESSURE};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Severity {
    Note,
    Warning,
}

#[derive(Clone, Debug)]
pub struct Warning {
    pub severity: Severity,
    /// The parameter the message refers to.
    pub node: &'static str,
    pub text: String,
}

pub fn evaluate(g: &ParamGraph) -> Vec<Warning> {
    let mut out = Vec::new();
    let note = |out: &mut Vec<Warning>, node, text| {
        out.push(Warning {
            severity: Severity::Note,
            node,
            text,
        })
    };
    let warn = |out: &mut Vec<Warning>, node, text| {
        out.push(Warning {
            severity: Severity::Warning,
            node,
            text,
        })
    };

    let rho = g.f64("body.density");
    if !(900.0..=8000.0).contains(&rho) {
        warn(
            &mut out,
            "body.density",
            format!("Density {:.0} kg/m³ is outside 900 to 8000 kg/m³.", rho),
        );
    }

    let r = g.f64("body.radius");
    if r < 4.0e5 {
        note(
            &mut out,
            "body.radius",
            format!(
                "Radius {:.0} km is below 400 km; the body is unlikely to be spheroidal.",
                r / 1e3
            ),
        );
    }

    let star_age = g.f64("star.age");
    let t_ms = g.f64("star.main_sequence_lifetime");
    if star_age > t_ms {
        warn(
            &mut out,
            "star.age",
            format!("Star age {star_age:.2} Gyr exceeds the main-sequence lifetime {t_ms:.2} Gyr."),
        );
    }

    let escaping = g.str("atmosphere.retention_warning");
    if !escaping.is_empty() {
        warn(
            &mut out,
            "atmosphere.composition",
            format!(
                "Jeans parameter below {JEANS_RETENTION_THRESHOLD:.0} at the exobase temperature {:.0} K for {escaping}. These species are unlikely to survive {:.1} Gyr.",
                g.f64("atmosphere.exobase_temperature"),
                g.f64("body.age")
            ),
        );
    }

    let p = g.f64("atmosphere.surface_pressure");
    if let Some(inv) = g.get("volatiles.inventory").as_map() {
        for (species, amount) in inv {
            if *amount <= 0.0 {
                continue;
            }
            let Some(v) = volatile(species) else { continue };
            if p < v.p_triple {
                note(
                    &mut out,
                    "atmosphere.surface_pressure",
                    format!(
                        "Surface pressure {:.4} bar is below the {species} triple point ({:.4} bar); no liquid {species} can exist.",
                        p / 1e5,
                        v.p_triple / 1e5
                    ),
                );
            }
        }
    }

    // F238: the present spin against the history its tides allow.
    let age = g.f64("body.age");
    let mode = g.str("rotation.mode");
    if mode == "free" {
        let t_lock = g.f64("rotation.locking_timescale");
        if t_lock < age {
            if g.f64("atmosphere.surface_pressure") >= THERMAL_TIDE_PRESSURE {
                note(
                    &mut out,
                    "rotation.period",
                    format!(
                        "Its primary's tide would have locked its spin within {t_lock:.2} Gyr; a thick atmosphere's thermal tides can hold a spin asynchronous against it (Leconte et al. 2015), which the model does not compute."
                    ),
                );
            } else {
                warn(
                    &mut out,
                    "rotation.period",
                    format!(
                        "Its primary's tide would have locked its spin within {t_lock:.2} Gyr, less than its age of {age:.2} Gyr: a free spin needs a reason the model does not have."
                    ),
                );
            }
        }
    } else if g.f64("rotation.capture_time") < 0.0 {
        let born = g.f64("rotation.primordial_period") / 3600.0;
        warn(
            &mut out,
            "rotation.period",
            format!(
                "Its tides would not bring the {born:.1} h spin it was born with to its present {mode} state within its age of {age:.2} Gyr."
            ),
        );
    }
    if g.f64("satellite.mass") > 0.0 {
        let t_sat = g.f64("satellite.locking_timescale");
        let mass = g.f64("body.mass") + g.f64("satellite.mass");
        let distance = g.f64("satellite.semi_major_axis");
        let orbit =
            std::f64::consts::TAU * (distance.powi(3) / (crate::constants::G * mass)).sqrt();
        let synchronous = (g.f64("rotation.period") / orbit - 1.0).abs() < 1e-3;
        if t_sat < age && !synchronous {
            warn(
                &mut out,
                "satellite.mass",
                format!(
                    "Its satellite's tide would have locked its spin to the satellite's orbit within {t_sat:.3} Gyr, less than its age of {age:.2} Gyr."
                ),
            );
        }
    }

    let column = g.f64("atmosphere.column_mass");
    let gravity = g.f64("body.gravity");
    if column > 1.0e5 && gravity < 3.0 {
        note(
            &mut out,
            "atmosphere.surface_pressure",
            format!(
                "Column mass {column:.3e} kg/m² on gravity {gravity:.2} m/s² has no solar-system analogue apart from Titan."
            ),
        );
    }

    if let Some(states) = g.get("volatiles.phase_state").as_map() {
        let any_liquid = states
            .values()
            .any(|c| *c == PHASE_LIQUID || *c == PHASE_LIQUID_AND_SOLID);
        let has_inventory = g
            .get("volatiles.inventory")
            .as_map()
            .is_some_and(|m| m.values().any(|a| *a > 0.0));
        if has_inventory && !any_liquid && p > 0.0 {
            let t = g.f64("atmosphere.surface_temperature_mean");
            note(
                &mut out,
                "volatiles.inventory",
                format!(
                    "No inventory species is liquid at the mean surface temperature {t:.0} K and pressure {:.3} bar.",
                    p / 1e5
                ),
            );
        }
    }

    let d_au = g.f64("orbit.star_distance") / AU;
    if g.str("orbit.primary") == "planet" && d_au < 0.1 {
        note(
            &mut out,
            "parent.semi_major_axis",
            format!("Parent planet at {d_au:.3} AU; tides from the star will also matter."),
        );
    }

    out
}
