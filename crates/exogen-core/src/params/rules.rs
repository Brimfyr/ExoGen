// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! Plausibility rules (docs/parameters.md). Rules warn and explain; they never block.

use super::ParamGraph;
use super::derive::physics::volatile;
use super::derive::{PHASE_LIQUID, PHASE_LIQUID_AND_SOLID};
use crate::constants::{AU, JEANS_RETENTION_THRESHOLD};

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

    if g.is_pinned("rotation.mode") && g.str("rotation.mode") == "free" {
        let t_lock = g.f64("rotation.locking_timescale");
        let age = g.f64("body.age");
        if t_lock < age {
            warn(
                &mut out,
                "rotation.mode",
                format!(
                    "Locking timescale {t_lock:.2} Gyr is shorter than the body age {age:.2} Gyr, but rotation is pinned free."
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
