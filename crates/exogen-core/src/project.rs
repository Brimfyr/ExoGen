// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! Project files (ADR 0003, ADR 0004): everything needed to reproduce a world.
//!
//! A project stores the engine version, a name, and the parameter overrides: modified
//! roots and multipliers, and pinned derived nodes with their values. Everything else is
//! derived on load. A project written by a newer engine than the one reading it is refused.

use std::collections::BTreeMap;
use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::params::{ParamGraph, Value};

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct Project {
    pub version: String,
    pub name: String,
    #[serde(default)]
    pub description: String,
    #[serde(default)]
    pub values: BTreeMap<String, Value>,
}

#[derive(Debug, thiserror::Error)]
pub enum ProjectError {
    #[error("could not read {path}: {source}")]
    Read {
        path: String,
        source: std::io::Error,
    },
    #[error("could not write {path}: {source}")]
    Write {
        path: String,
        source: std::io::Error,
    },
    #[error("project file is not valid TOML: {0}")]
    Parse(#[from] toml::de::Error),
    #[error("could not serialise project: {0}")]
    Serialise(#[from] toml::ser::Error),
    #[error("project was written by engine version {0}, newer than this engine ({1})")]
    NewerVersion(String, String),
}

impl Project {
    /// Capture the current overrides of a graph.
    pub fn from_graph(name: impl Into<String>, graph: &ParamGraph) -> Self {
        Self {
            version: crate::VERSION.to_string(),
            name: name.into(),
            description: String::new(),
            values: graph.overrides(),
        }
    }

    /// Build a graph from this project. Returns the graph and the ids it did not recognise.
    pub fn to_graph(&self) -> (ParamGraph, Vec<String>) {
        let mut g = ParamGraph::new();
        let unknown = g.apply_overrides(&self.values);
        (g, unknown)
    }

    pub fn from_toml(text: &str) -> Result<Self, ProjectError> {
        let p: Project = toml::from_str(text)?;
        if version_tuple(&p.version) > version_tuple(crate::VERSION) {
            return Err(ProjectError::NewerVersion(
                p.version.clone(),
                crate::VERSION.to_string(),
            ));
        }
        Ok(p)
    }

    pub fn to_toml(&self) -> Result<String, ProjectError> {
        Ok(toml::to_string_pretty(self)?)
    }

    pub fn load(path: &Path) -> Result<Self, ProjectError> {
        let text = std::fs::read_to_string(path).map_err(|source| ProjectError::Read {
            path: path.display().to_string(),
            source,
        })?;
        Self::from_toml(&text)
    }

    pub fn save(&self, path: &Path) -> Result<(), ProjectError> {
        let text = self.to_toml()?;
        std::fs::write(path, text).map_err(|source| ProjectError::Write {
            path: path.display().to_string(),
            source,
        })
    }
}

fn version_tuple(v: &str) -> (u64, u64, u64) {
    let mut parts = v.split('.').map(|p| p.trim().parse::<u64>().unwrap_or(0));
    (
        parts.next().unwrap_or(0),
        parts.next().unwrap_or(0),
        parts.next().unwrap_or(0),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn project_round_trips_through_toml_and_rebuilds_the_same_graph() {
        let mut g = ParamGraph::new();
        g.set("body.radius", 3.3895e6).unwrap();
        g.set("body.density", 3934.0).unwrap();
        g.set("atmosphere.surface_pressure", 610.0).unwrap();
        g.set("rotation.mode", "free").unwrap();
        let mut inv = BTreeMap::new();
        inv.insert("H2O".to_string(), 25.0);
        g.set("volatiles.inventory", Value::Map(inv)).unwrap();

        let p = Project::from_graph("Mars-like", &g);
        let text = p.to_toml().unwrap();
        let back = Project::from_toml(&text).unwrap();
        assert_eq!(back, p);
        let (g2, unknown) = back.to_graph();
        assert!(unknown.is_empty());
        for def in g.defs() {
            assert_eq!(
                g.get(def.id),
                g2.get(def.id),
                "{} differs after round trip",
                def.id
            );
        }
        assert!(g2.is_pinned("rotation.mode"));
    }

    #[test]
    fn newer_projects_are_refused_and_unknown_ids_are_reported() {
        let text = "version = \"99.0.0\"\nname = \"x\"\n";
        assert!(matches!(
            Project::from_toml(text),
            Err(ProjectError::NewerVersion(..))
        ));
        let text = "version = \"0.0.1\"\nname = \"x\"\n[values]\n\"no.such_node\" = 1.0\n\"body.radius\" = 1.0e6\n";
        let p = Project::from_toml(text).unwrap();
        let (g, unknown) = p.to_graph();
        assert_eq!(unknown, vec!["no.such_node".to_string()]);
        assert_eq!(g.f64("body.radius"), 1.0e6);
    }
}
