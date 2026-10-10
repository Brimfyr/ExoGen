// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! Parameter values and their types.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

/// A piecewise-linear history over simulated time. Times are in Gyr since the body formed
/// and must be ascending. Evaluation clamps outside the defined range.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Curve {
    pub points: Vec<(f64, f64)>,
}

impl Curve {
    pub fn constant(v: f64) -> Self {
        Self {
            points: vec![(0.0, v)],
        }
    }

    pub fn new(points: Vec<(f64, f64)>) -> Self {
        debug_assert!(
            points.windows(2).all(|w| w[0].0 <= w[1].0),
            "curve times must ascend"
        );
        Self { points }
    }

    pub fn eval(&self, t: f64) -> f64 {
        match self.points.len() {
            0 => f64::NAN,
            1 => self.points[0].1,
            _ => {
                if t <= self.points[0].0 {
                    return self.points[0].1;
                }
                let last = self.points[self.points.len() - 1];
                if t >= last.0 {
                    return last.1;
                }
                let k = self.points.partition_point(|p| p.0 <= t);
                let (t0, v0) = self.points[k - 1];
                let (t1, v1) = self.points[k];
                if t1 > t0 {
                    v0 + (v1 - v0) * (t - t0) / (t1 - t0)
                } else {
                    v1
                }
            }
        }
    }

    pub fn is_constant(&self) -> bool {
        self.points.len() <= 1 || self.points.iter().all(|p| p.1 == self.points[0].1)
    }
}

/// A parameter value. Enumerations and free text are both strings; the node type says which.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(untagged)]
pub enum Value {
    F64(f64),
    Bool(bool),
    Str(String),
    Curve(Curve),
    Map(BTreeMap<String, f64>),
    List(Vec<String>),
}

impl Value {
    pub fn as_f64(&self) -> Option<f64> {
        match self {
            Value::F64(v) => Some(*v),
            Value::Bool(b) => Some(if *b { 1.0 } else { 0.0 }),
            _ => None,
        }
    }

    pub fn as_bool(&self) -> Option<bool> {
        match self {
            Value::Bool(b) => Some(*b),
            Value::F64(v) => Some(*v != 0.0),
            _ => None,
        }
    }

    pub fn as_str(&self) -> Option<&str> {
        match self {
            Value::Str(s) => Some(s),
            _ => None,
        }
    }

    pub fn as_curve(&self) -> Option<&Curve> {
        match self {
            Value::Curve(c) => Some(c),
            _ => None,
        }
    }

    pub fn as_map(&self) -> Option<&BTreeMap<String, f64>> {
        match self {
            Value::Map(m) => Some(m),
            _ => None,
        }
    }

    pub fn as_list(&self) -> Option<&[String]> {
        match self {
            Value::List(l) => Some(l),
            _ => None,
        }
    }

    /// Canonical bytes for hashing: a type tag followed by the exact contents. Floats hash
    /// by bit pattern so that equal values hash equal on every platform.
    pub fn write_canonical(&self, out: &mut Vec<u8>) {
        match self {
            Value::F64(v) => {
                out.push(1);
                out.extend_from_slice(&v.to_bits().to_le_bytes());
            }
            Value::Bool(b) => {
                out.push(2);
                out.push(*b as u8);
            }
            Value::Str(s) => {
                out.push(3);
                out.extend_from_slice(&(s.len() as u64).to_le_bytes());
                out.extend_from_slice(s.as_bytes());
            }
            Value::Curve(c) => {
                out.push(4);
                out.extend_from_slice(&(c.points.len() as u64).to_le_bytes());
                for (t, v) in &c.points {
                    out.extend_from_slice(&t.to_bits().to_le_bytes());
                    out.extend_from_slice(&v.to_bits().to_le_bytes());
                }
            }
            Value::Map(m) => {
                out.push(5);
                out.extend_from_slice(&(m.len() as u64).to_le_bytes());
                for (k, v) in m {
                    out.extend_from_slice(&(k.len() as u64).to_le_bytes());
                    out.extend_from_slice(k.as_bytes());
                    out.extend_from_slice(&v.to_bits().to_le_bytes());
                }
            }
            Value::List(l) => {
                out.push(6);
                out.extend_from_slice(&(l.len() as u64).to_le_bytes());
                for s in l {
                    out.extend_from_slice(&(s.len() as u64).to_le_bytes());
                    out.extend_from_slice(s.as_bytes());
                }
            }
        }
    }
}

impl From<f64> for Value {
    fn from(v: f64) -> Self {
        Value::F64(v)
    }
}

impl From<bool> for Value {
    fn from(v: bool) -> Self {
        Value::Bool(v)
    }
}

impl From<&str> for Value {
    fn from(v: &str) -> Self {
        Value::Str(v.to_string())
    }
}

impl From<String> for Value {
    fn from(v: String) -> Self {
        Value::Str(v)
    }
}

impl From<Curve> for Value {
    fn from(v: Curve) -> Self {
        Value::Curve(v)
    }
}

/// The type a node's value must have.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ValueType {
    /// A real number with a soft range for the UI and plausibility checks.
    F64 {
        min: f64,
        max: f64,
    },
    Bool,
    /// One of a fixed set of variants.
    Enum(&'static [&'static str]),
    /// Free text, display only.
    Text,
    /// A history curve over simulated time in Gyr.
    Curve,
    /// A map from species or component name to a number; keys limited to the given set.
    Map(&'static [&'static str]),
    /// A list of strings (species names, feature ids).
    List,
}

impl ValueType {
    /// Whether a value is acceptable for this type (variant and enum membership).
    pub fn accepts(&self, v: &Value) -> bool {
        match (self, v) {
            (ValueType::F64 { .. }, Value::F64(x)) => x.is_finite(),
            (ValueType::Bool, Value::Bool(_)) => true,
            (ValueType::Enum(variants), Value::Str(s)) => variants.contains(&s.as_str()),
            (ValueType::Text, Value::Str(_)) => true,
            (ValueType::Curve, Value::Curve(c)) => !c.points.is_empty(),
            (ValueType::Curve, Value::F64(x)) => x.is_finite(),
            (ValueType::Map(keys), Value::Map(m)) => m.keys().all(|k| keys.contains(&k.as_str())),
            (ValueType::List, Value::List(_)) => true,
            _ => false,
        }
    }

    /// Coerce a loosely typed value (for example a number given for a curve) into the exact
    /// representation this type stores.
    pub fn normalise(&self, v: Value) -> Value {
        match (self, v) {
            (ValueType::Curve, Value::F64(x)) => Value::Curve(Curve::constant(x)),
            (ValueType::Bool, Value::F64(x)) => Value::Bool(x != 0.0),
            (_, v) => v,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn curve_interpolates_and_clamps() {
        let c = Curve::new(vec![(0.0, 1.0), (2.0, 3.0), (4.0, 3.0)]);
        assert_eq!(c.eval(-1.0), 1.0);
        assert_eq!(c.eval(1.0), 2.0);
        assert_eq!(c.eval(3.0), 3.0);
        assert_eq!(c.eval(9.0), 3.0);
        assert!(Curve::constant(2.0).is_constant());
    }

    #[test]
    fn canonical_bytes_distinguish_values() {
        let mut a = Vec::new();
        let mut b = Vec::new();
        Value::F64(1.0).write_canonical(&mut a);
        Value::F64(-1.0).write_canonical(&mut b);
        assert_ne!(a, b);
        let mut c = Vec::new();
        Value::F64(1.0).write_canonical(&mut c);
        assert_eq!(a, c);
    }

    #[test]
    fn values_round_trip_through_toml() {
        let mut m = BTreeMap::new();
        m.insert("N2".to_string(), 0.78);
        let values = vec![
            Value::F64(2.5),
            Value::Bool(true),
            Value::Str("locked".into()),
            Value::Curve(Curve::new(vec![(0.0, 1.0), (4.5, 2.0)])),
            Value::Map(m),
        ];
        for v in values {
            let mut table = BTreeMap::new();
            table.insert("v", v.clone());
            let text = toml::to_string(&table).unwrap();
            let back: BTreeMap<String, Value> = toml::from_str(&text).unwrap();
            assert_eq!(back["v"], v, "{text}");
        }
    }
}
