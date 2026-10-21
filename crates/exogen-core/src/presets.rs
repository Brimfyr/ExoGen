// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! Reference-body presets. Each is a project file embedded in the binary. They pin the
//! observed inputs of a solar-system body so that the derived values can be checked
//! against the measured ones (PLAN.md section 6) and so that users have known starting
//! points.

use crate::project::Project;

macro_rules! presets {
    ($( $name:literal => $file:literal ),* $(,)?) => {
        /// Preset names in display order.
        pub const PRESET_NAMES: &[&str] = &[$($name),*];

        fn preset_text(name: &str) -> Option<&'static str> {
            match name {
                $( $name => Some(include_str!(concat!("../presets/", $file))), )*
                _ => None,
            }
        }
    };
}

presets! {
    "Earth" => "earth.toml",
    "Mars" => "mars.toml",
    "Moon" => "moon.toml",
    "Mercury" => "mercury.toml",
    "Venus" => "venus.toml",
    "Io" => "io.toml",
    "Europa" => "europa.toml",
    "Ganymede" => "ganymede.toml",
    "Callisto" => "callisto.toml",
    "Titan" => "titan.toml",
    "Enceladus" => "enceladus.toml",
    "Mimas" => "mimas.toml",
    "Miranda" => "miranda.toml",
    "Triton" => "triton.toml",
    "Pluto" => "pluto.toml",
    "Charon" => "charon.toml",
    "Ceres" => "ceres.toml",
}

/// Load a preset by name.
pub fn preset(name: &str) -> Option<Project> {
    let text = preset_text(name)?;
    Some(
        Project::from_toml(text)
            .unwrap_or_else(|e| panic!("embedded preset {name} is invalid: {e}")),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_preset_parses_and_only_sets_known_nodes() {
        for name in PRESET_NAMES {
            let p = preset(name).unwrap();
            let (g, unknown) = p.to_graph();
            assert!(unknown.is_empty(), "{name}: unknown ids {unknown:?}");
            assert!(g.f64("body.gravity") > 0.0);
        }
    }
}
