// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! The parameter graph (PLAN.md 3.3, docs/parameters.md).
//!
//! Every parameter is a node. Root nodes hold inputs with defaults. Derived nodes compute
//! their value from upstream nodes and are recomputed whenever anything changes; setting
//! one pins it, and downstream nodes keep deriving from the pinned value. Multipliers are
//! artistic scale factors that are never derived. Stage outputs are filled in by the
//! pipeline and are read-only in the UI.
//!
//! Recomputation is a full pass over the topological order. With a few hundred nodes this
//! costs microseconds, which keeps the logic simple and makes every value consistent after
//! every change. A node may also declare lagged dependencies: the far side of a feedback
//! (the surface's reflectance sets the temperature that sets the surface), read from the
//! previous pass and left out of the order so the graph stays acyclic; recomputation
//! repeats its pass until nothing changes, a few passes at most.

pub mod derive;
pub mod registry;
pub mod rules;
pub mod value;

use std::collections::{BTreeMap, HashMap};

pub use value::{Curve, Value, ValueType};

/// What a node is.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Kind {
    /// An input with a default.
    Root,
    /// Computed from upstream nodes; may be pinned.
    Derived,
    /// Artistic scale factor, default 1, never derived.
    Multiplier,
    /// Written by the pipeline; read-only.
    StageOutput,
}

/// Whether a node is shown in the normal parameter panel or only in expert mode.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Visibility {
    Normal,
    Expert,
}

/// A default for a root or multiplier node, expressible as a constant.
#[derive(Clone, Copy, Debug)]
pub enum Default {
    F64(f64),
    Bool(bool),
    Str(&'static str),
    Map(&'static [(&'static str, f64)]),
    Curve(f64),
    None,
}

impl Default {
    pub fn value(&self) -> Option<Value> {
        match self {
            Default::F64(v) => Some(Value::F64(*v)),
            Default::Bool(b) => Some(Value::Bool(*b)),
            Default::Str(s) => Some(Value::Str((*s).to_string())),
            Default::Map(entries) => {
                let mut m = BTreeMap::new();
                for (k, v) in *entries {
                    m.insert((*k).to_string(), *v);
                }
                Some(Value::Map(m))
            }
            Default::Curve(v) => Some(Value::Curve(Curve::constant(*v))),
            Default::None => None,
        }
    }
}

/// Context handed to derivation functions: read-only access to upstream values. In tests
/// the context can record which nodes a derivation read, to check declared dependencies.
pub struct Ctx<'a> {
    graph: &'a ParamGraph,
    accessed: Option<&'a std::cell::RefCell<std::collections::BTreeSet<String>>>,
}

impl<'a> Ctx<'a> {
    pub fn new(graph: &'a ParamGraph) -> Self {
        Self {
            graph,
            accessed: None,
        }
    }

    pub fn tracking(
        graph: &'a ParamGraph,
        accessed: &'a std::cell::RefCell<std::collections::BTreeSet<String>>,
    ) -> Self {
        Self {
            graph,
            accessed: Some(accessed),
        }
    }

    pub fn value(&self, id: &str) -> &Value {
        if let Some(set) = self.accessed {
            set.borrow_mut().insert(id.to_string());
        }
        self.graph.get(id)
    }

    pub fn f64(&self, id: &str) -> f64 {
        self.value(id)
            .as_f64()
            .unwrap_or_else(|| panic!("node {id} is not numeric"))
    }

    pub fn bool(&self, id: &str) -> bool {
        self.value(id)
            .as_bool()
            .unwrap_or_else(|| panic!("node {id} is not boolean"))
    }

    pub fn str(&self, id: &str) -> &str {
        self.value(id)
            .as_str()
            .unwrap_or_else(|| panic!("node {id} is not a string"))
    }

    pub fn map(&self, id: &str) -> &BTreeMap<String, f64> {
        self.value(id)
            .as_map()
            .unwrap_or_else(|| panic!("node {id} is not a map"))
    }

    pub fn curve(&self, id: &str) -> &Curve {
        self.value(id)
            .as_curve()
            .unwrap_or_else(|| panic!("node {id} is not a curve"))
    }

    pub fn map_get(&self, id: &str, key: &str) -> f64 {
        self.map(id).get(key).copied().unwrap_or(0.0)
    }
}

/// Static definition of a node.
pub struct NodeDef {
    pub id: &'static str,
    pub group: &'static str,
    pub label: &'static str,
    /// Display unit suffix, for example "km" or "bar".
    pub unit: &'static str,
    /// Display value = stored value / scale.
    pub scale: f64,
    pub kind: Kind,
    pub visibility: Visibility,
    pub ty: ValueType,
    pub default: Default,
    pub deps: &'static [&'static str],
    /// Upstream nodes read from the previous recompute pass: the far side of a feedback,
    /// left out of the order so the graph stays acyclic and settled by repeated passes.
    pub lagged: &'static [&'static str],
    pub derive: Option<fn(&Ctx) -> Value>,
    pub doc: &'static str,
    pub source: &'static str,
    /// A derived node whose value is an input in spirit, a default the preset is expected
    /// to set, such as what the body is made of: setting it is giving an input, not
    /// pinning an outcome.
    pub input: bool,
}

impl NodeDef {
    pub fn display(&self, v: f64) -> f64 {
        v / self.scale
    }

    pub fn from_display(&self, v: f64) -> f64 {
        v * self.scale
    }
}

/// Lock state of a derived node.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LockState {
    Derived,
    Set,
}

#[derive(Debug, thiserror::Error)]
pub enum SetError {
    #[error("unknown parameter {0}")]
    Unknown(String),
    #[error("parameter {0} is written by the pipeline and cannot be set")]
    ReadOnly(String),
    #[error("value has the wrong type for parameter {0}")]
    WrongType(String),
}

/// Passes a recompute makes at most while lagged feedbacks settle.
const RECOMPUTE_PASSES: usize = 4;

/// Equality that lets a NaN settle (a NaN is never equal to itself).
fn same_value(a: &Value, b: &Value) -> bool {
    match (a, b) {
        (Value::F64(x), Value::F64(y)) => x == y || (x.is_nan() && y.is_nan()),
        _ => a == b,
    }
}

/// The live graph: definitions, values and lock states.
pub struct ParamGraph {
    defs: &'static [NodeDef],
    index: HashMap<&'static str, usize>,
    order: Vec<usize>,
    values: Vec<Value>,
    state: Vec<LockState>,
}

impl Clone for ParamGraph {
    fn clone(&self) -> Self {
        Self {
            defs: self.defs,
            index: self.index.clone(),
            order: self.order.clone(),
            values: self.values.clone(),
            state: self.state.clone(),
        }
    }
}

impl std::default::Default for ParamGraph {
    fn default() -> Self {
        Self::new()
    }
}

impl ParamGraph {
    /// A graph with every root at its default and every derived node computed.
    pub fn new() -> Self {
        let defs = registry::nodes();
        let mut index = HashMap::with_capacity(defs.len());
        for (i, d) in defs.iter().enumerate() {
            assert!(
                index.insert(d.id, i).is_none(),
                "duplicate node id {}",
                d.id
            );
        }
        let order = topological_order(defs, &index);
        let mut values = Vec::with_capacity(defs.len());
        let mut state = Vec::with_capacity(defs.len());
        for d in defs {
            let v = match d.kind {
                Kind::Root | Kind::Multiplier => d
                    .default
                    .value()
                    .unwrap_or_else(|| panic!("root node {} has no default", d.id)),
                Kind::Derived | Kind::StageOutput => {
                    d.default.value().unwrap_or(placeholder(&d.ty))
                }
            };
            values.push(d.ty.normalise(v));
            state.push(LockState::Derived);
        }
        let mut g = Self {
            defs,
            index,
            order,
            values,
            state,
        };
        g.recompute();
        g
    }

    pub fn defs(&self) -> &'static [NodeDef] {
        self.defs
    }

    pub fn def(&self, id: &str) -> &NodeDef {
        &self.defs[self.idx(id)]
    }

    pub fn has(&self, id: &str) -> bool {
        self.index.contains_key(id)
    }

    fn idx(&self, id: &str) -> usize {
        *self
            .index
            .get(id)
            .unwrap_or_else(|| panic!("unknown parameter {id}"))
    }

    pub fn get(&self, id: &str) -> &Value {
        &self.values[self.idx(id)]
    }

    pub fn f64(&self, id: &str) -> f64 {
        self.get(id).as_f64().unwrap_or(f64::NAN)
    }

    pub fn str(&self, id: &str) -> &str {
        self.get(id).as_str().unwrap_or("")
    }

    pub fn bool(&self, id: &str) -> bool {
        self.get(id).as_bool().unwrap_or(false)
    }

    pub fn curve(&self, id: &str) -> Option<&Curve> {
        self.get(id).as_curve()
    }

    pub fn is_pinned(&self, id: &str) -> bool {
        let i = self.idx(id);
        self.defs[i].kind == Kind::Derived && self.state[i] == LockState::Set
    }

    /// Whether a root or multiplier differs from its default.
    pub fn is_modified(&self, id: &str) -> bool {
        let i = self.idx(id);
        match self.defs[i].kind {
            Kind::Root | Kind::Multiplier => {
                self.defs[i]
                    .default
                    .value()
                    .map(|d| self.defs[i].ty.normalise(d))
                    != Some(self.values[i].clone())
            }
            Kind::Derived => self.state[i] == LockState::Set,
            Kind::StageOutput => false,
        }
    }

    /// Set a root, multiplier or derived node. Setting a derived node pins it.
    pub fn set(&mut self, id: &str, v: impl Into<Value>) -> Result<(), SetError> {
        let i = *self
            .index
            .get(id)
            .ok_or_else(|| SetError::Unknown(id.to_string()))?;
        let d = &self.defs[i];
        if d.kind == Kind::StageOutput {
            return Err(SetError::ReadOnly(id.to_string()));
        }
        let v = d.ty.normalise(v.into());
        if !d.ty.accepts(&v) {
            return Err(SetError::WrongType(id.to_string()));
        }
        self.values[i] = v;
        if d.kind == Kind::Derived {
            self.state[i] = LockState::Set;
        }
        self.recompute();
        Ok(())
    }

    /// Write a stage output. Only the pipeline calls this.
    pub fn set_stage_output(&mut self, id: &str, v: impl Into<Value>) {
        let i = self.idx(id);
        debug_assert_eq!(self.defs[i].kind, Kind::StageOutput);
        self.values[i] = self.defs[i].ty.normalise(v.into());
    }

    /// Unpin a derived node, or restore a root or multiplier to its default.
    pub fn reset(&mut self, id: &str) {
        let i = self.idx(id);
        match self.defs[i].kind {
            Kind::Derived => self.state[i] = LockState::Derived,
            Kind::Root | Kind::Multiplier => {
                if let Some(d) = self.defs[i].default.value() {
                    self.values[i] = self.defs[i].ty.normalise(d);
                }
            }
            Kind::StageOutput => {}
        }
        self.recompute();
    }

    /// Recompute every unpinned derived node in topological order, and again while any
    /// value changed, up to a few passes, so that the feedbacks the lagged dependencies
    /// express settle. Without lagged edges the second pass changes nothing.
    pub fn recompute(&mut self) {
        for _ in 0..RECOMPUTE_PASSES {
            if !self.recompute_pass() {
                break;
            }
        }
    }

    /// One pass over the order; whether any derived value changed.
    fn recompute_pass(&mut self) -> bool {
        let mut changed = false;
        for k in 0..self.order.len() {
            let i = self.order[k];
            let d = &self.defs[i];
            if d.kind != Kind::Derived || self.state[i] == LockState::Set {
                continue;
            }
            let Some(f) = d.derive else { continue };
            let v = {
                let ctx = Ctx::new(self);
                f(&ctx)
            };
            let v = d.ty.normalise(v);
            debug_assert!(d.ty.accepts(&v), "derivation of {} produced {:?}", d.id, v);
            if !same_value(&self.values[i], &v) {
                self.values[i] = v;
                changed = true;
            }
        }
        changed
    }

    /// Everything a project file needs to recreate this graph: modified roots and
    /// multipliers, and pinned derived nodes with their values.
    pub fn overrides(&self) -> BTreeMap<String, Value> {
        let mut out = BTreeMap::new();
        for (i, d) in self.defs.iter().enumerate() {
            if self.is_modified(d.id) {
                out.insert(d.id.to_string(), self.values[i].clone());
            }
        }
        out
    }

    /// Apply a set of overrides. Unknown ids are returned rather than failing, so a project
    /// written by a newer version degrades gracefully.
    pub fn apply_overrides(&mut self, overrides: &BTreeMap<String, Value>) -> Vec<String> {
        let mut unknown = Vec::new();
        for (id, v) in overrides {
            let Some(&i) = self.index.get(id.as_str()) else {
                unknown.push(id.clone());
                continue;
            };
            let d = &self.defs[i];
            if d.kind == Kind::StageOutput {
                continue;
            }
            let v = d.ty.normalise(v.clone());
            if d.ty.accepts(&v) {
                self.values[i] = v;
                if d.kind == Kind::Derived {
                    self.state[i] = LockState::Set;
                }
            } else {
                unknown.push(id.clone());
            }
        }
        self.recompute();
        unknown
    }

    /// Stable 64-bit hash of the named nodes' values, for stage cache keys.
    pub fn hash_of(&self, ids: &[&str]) -> u64 {
        let mut bytes = Vec::new();
        for id in ids {
            bytes.extend_from_slice(id.as_bytes());
            bytes.push(0);
            self.get(id).write_canonical(&mut bytes);
        }
        fnv1a64(&bytes)
    }

    /// Hash of every node in a group.
    pub fn hash_of_group(&self, group: &str) -> u64 {
        let ids: Vec<&str> = self
            .defs
            .iter()
            .filter(|d| d.group == group)
            .map(|d| d.id)
            .collect();
        self.hash_of(&ids)
    }

    /// Groups in definition order, without duplicates.
    pub fn groups(&self) -> Vec<&'static str> {
        let mut out: Vec<&'static str> = Vec::new();
        for d in self.defs {
            if !out.contains(&d.group) {
                out.push(d.group);
            }
        }
        out
    }

    /// Nodes of a group in definition order.
    pub fn nodes_in(&self, group: &str) -> impl Iterator<Item = &'static NodeDef> + '_ {
        let group = group.to_string();
        self.defs.iter().filter(move |d| d.group == group)
    }

    /// Plausibility warnings for the current state.
    pub fn warnings(&self) -> Vec<rules::Warning> {
        rules::evaluate(self)
    }
}

fn placeholder(ty: &ValueType) -> Value {
    match ty {
        ValueType::F64 { .. } => Value::F64(0.0),
        ValueType::Bool => Value::Bool(false),
        ValueType::Enum(variants) => Value::Str(variants[0].to_string()),
        ValueType::Text => Value::Str(String::new()),
        ValueType::Curve => Value::Curve(Curve::constant(0.0)),
        ValueType::Map(_) => Value::Map(BTreeMap::new()),
        ValueType::List => Value::List(Vec::new()),
    }
}

/// Kahn's algorithm over the declared dependencies. Panics on unknown dependencies or
/// cycles, which are programming errors caught by the registry test.
fn topological_order(defs: &[NodeDef], index: &HashMap<&'static str, usize>) -> Vec<usize> {
    let n = defs.len();
    let mut indegree = vec![0usize; n];
    let mut dependents: Vec<Vec<usize>> = vec![Vec::new(); n];
    for (i, d) in defs.iter().enumerate() {
        for dep in d.deps {
            let j = *index
                .get(dep)
                .unwrap_or_else(|| panic!("node {} depends on unknown node {}", d.id, dep));
            dependents[j].push(i);
            indegree[i] += 1;
        }
    }
    let mut ready: Vec<usize> = (0..n).filter(|&i| indegree[i] == 0).collect();
    ready.reverse();
    let mut order = Vec::with_capacity(n);
    while let Some(i) = ready.pop() {
        order.push(i);
        for &k in &dependents[i] {
            indegree[k] -= 1;
            if indegree[k] == 0 {
                ready.push(k);
            }
        }
    }
    assert_eq!(order.len(), n, "parameter graph has a dependency cycle");
    order
}

/// FNV-1a over bytes, 64 bits. Stable across platforms and versions.
pub fn fnv1a64(bytes: &[u8]) -> u64 {
    let mut h = 0xcbf2_9ce4_8422_2325u64;
    for &b in bytes {
        h ^= b as u64;
        h = h.wrapping_mul(0x0000_0100_0000_01b3);
    }
    h
}
