// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! The parameter panel: groups of rows, each with label, value, unit and lock state.
//! Derived values are dimmed; clicking "pin" turns one into an editable override; "reset"
//! returns it to derivation. Hovering a label shows the basis and the source.

use std::collections::BTreeMap;

use eframe::egui;
use exogen_core::params::{Kind, NodeDef, ParamGraph, Value, ValueType, Visibility};

/// Draw the panel. Returns true if any value changed.
pub fn show(ui: &mut egui::Ui, graph: &mut ParamGraph, expert: bool, search: &str) -> bool {
    let mut changed = false;
    let search = search.trim().to_lowercase();
    let groups = graph.groups();
    for group in groups {
        let defs: Vec<&'static NodeDef> = graph
            .nodes_in(group)
            .filter(|d| expert || d.visibility != Visibility::Expert)
            .filter(|d| {
                search.is_empty()
                    || d.id.contains(&search)
                    || d.label.to_lowercase().contains(&search)
            })
            .collect();
        if defs.is_empty() {
            continue;
        }
        let default_open = matches!(group, "body" | "atmosphere") || !search.is_empty();
        egui::CollapsingHeader::new(group)
            .default_open(default_open)
            .show(ui, |ui| {
                egui::Grid::new(format!("grid-{group}"))
                    .num_columns(4)
                    .spacing([8.0, 2.0])
                    .striped(true)
                    .show(ui, |ui| {
                        for def in defs {
                            if row(ui, graph, def) {
                                changed = true;
                            }
                            ui.end_row();
                        }
                    });
            });
    }
    changed
}

fn state_of(graph: &ParamGraph, def: &NodeDef) -> &'static str {
    match def.kind {
        Kind::Root => {
            if graph.is_modified(def.id) {
                "set"
            } else {
                "default"
            }
        }
        Kind::Derived => {
            if graph.is_pinned(def.id) {
                "pinned"
            } else {
                "derived"
            }
        }
        Kind::Multiplier => "multiplier",
        Kind::StageOutput => "stage output",
    }
}

fn tooltip(ui: &mut egui::Ui, graph: &ParamGraph, def: &NodeDef) {
    ui.set_max_width(360.0);
    ui.label(egui::RichText::new(def.id).monospace());
    if !def.doc.is_empty() {
        ui.label(def.doc);
    }
    if !def.source.is_empty() {
        ui.label(format!("Source: {}", def.source));
    }
    if !def.deps.is_empty() {
        ui.label(format!("Derived from: {}", def.deps.join(", ")));
    }
    ui.label(format!("State: {}", state_of(graph, def)));
}

/// One row. Returns true when the graph changed.
fn row(ui: &mut egui::Ui, graph: &mut ParamGraph, def: &'static NodeDef) -> bool {
    let editable = match def.kind {
        Kind::Root | Kind::Multiplier => true,
        Kind::Derived => graph.is_pinned(def.id),
        Kind::StageOutput => false,
    };
    let dim = !editable;
    let label_text = if dim {
        egui::RichText::new(def.label).weak()
    } else {
        egui::RichText::new(def.label)
    };
    ui.label(label_text)
        .on_hover_ui(|ui| tooltip(ui, graph, def));

    let mut changed = false;
    let value = graph.get(def.id).clone();
    match (&def.ty, value) {
        (ValueType::F64 { min, max }, Value::F64(v)) => {
            if editable {
                let mut shown = def.display(v);
                let speed = (shown.abs() * 0.01)
                    .max(def.display(*max - *min) * 1e-4)
                    .max(1e-6);
                let resp = ui.add(
                    egui::DragValue::new(&mut shown)
                        .speed(speed)
                        .max_decimals(6),
                );
                if resp.changed() {
                    let _ = graph.set(def.id, def.from_display(shown));
                    changed = true;
                }
            } else {
                ui.label(
                    egui::RichText::new(format_num(def.display(v)))
                        .weak()
                        .monospace(),
                );
            }
        }
        (ValueType::Bool, Value::Bool(b)) => {
            if editable {
                let mut bb = b;
                if ui.checkbox(&mut bb, "").changed() {
                    let _ = graph.set(def.id, bb);
                    changed = true;
                }
            } else {
                ui.label(egui::RichText::new(if b { "yes" } else { "no" }).weak());
            }
        }
        (ValueType::Enum(variants), Value::Str(s)) => {
            if editable {
                let mut sel = s.clone();
                egui::ComboBox::from_id_salt(def.id)
                    .selected_text(sel.clone())
                    .show_ui(ui, |ui| {
                        for v in variants.iter() {
                            ui.selectable_value(&mut sel, (*v).to_string(), *v);
                        }
                    });
                if sel != s {
                    let _ = graph.set(def.id, sel.as_str());
                    changed = true;
                }
            } else {
                ui.label(egui::RichText::new(s).weak());
            }
        }
        (ValueType::Text, Value::Str(s)) => {
            ui.label(egui::RichText::new(if s.is_empty() { "-".to_string() } else { s }).weak());
        }
        (ValueType::Curve, Value::Curve(c)) => {
            let text = if c.is_constant() {
                format!("constant {}", format_num(c.points[0].1))
            } else {
                format!(
                    "curve, {} points, {} to {}",
                    c.points.len(),
                    format_num(c.points[0].1),
                    format_num(c.points[c.points.len() - 1].1)
                )
            };
            ui.label(egui::RichText::new(text).weak());
        }
        (ValueType::Map(keys), Value::Map(m)) => {
            if editable {
                if map_editor(ui, def, keys, &m)
                    .map(|new| {
                        let _ = graph.set(def.id, Value::Map(new));
                    })
                    .is_some()
                {
                    changed = true;
                }
            } else {
                let text = if m.is_empty() {
                    "-".to_string()
                } else {
                    m.iter()
                        .map(|(k, v)| format!("{k} {}", format_num(*v)))
                        .collect::<Vec<_>>()
                        .join(", ")
                };
                ui.label(egui::RichText::new(text).weak());
            }
        }
        (ValueType::List, Value::List(l)) => {
            ui.label(
                egui::RichText::new(if l.is_empty() {
                    "-".to_string()
                } else {
                    l.join(", ")
                })
                .weak(),
            );
        }
        (_, other) => {
            ui.label(egui::RichText::new(format!("{other:?}")).weak());
        }
    }

    ui.label(egui::RichText::new(def.unit).weak().small());

    match def.kind {
        Kind::Derived => {
            if graph.is_pinned(def.id) {
                if ui
                    .small_button("reset")
                    .on_hover_text("Return to the derived value")
                    .clicked()
                {
                    graph.reset(def.id);
                    changed = true;
                }
            } else if ui
                .small_button("pin")
                .on_hover_text("Set this value by hand; downstream nodes keep deriving from it")
                .clicked()
            {
                let v = graph.get(def.id).clone();
                let _ = graph.set(def.id, v);
                changed = true;
            }
        }
        Kind::Root | Kind::Multiplier => {
            if graph.is_modified(def.id) {
                if ui
                    .small_button("default")
                    .on_hover_text("Restore the default")
                    .clicked()
                {
                    graph.reset(def.id);
                    changed = true;
                }
            } else {
                ui.label("");
            }
        }
        Kind::StageOutput => {
            ui.label(egui::RichText::new("pipeline").weak().small());
        }
    }
    changed
}

/// Inline editor for a species map. Returns a new map when something changed.
fn map_editor(
    ui: &mut egui::Ui,
    def: &NodeDef,
    keys: &[&str],
    current: &BTreeMap<String, f64>,
) -> Option<BTreeMap<String, f64>> {
    let mut result = None;
    ui.vertical(|ui| {
        let summary = if current.is_empty() {
            "empty".to_string()
        } else {
            current
                .iter()
                .map(|(k, v)| format!("{k} {}", format_num(*v)))
                .collect::<Vec<_>>()
                .join(", ")
        };
        egui::CollapsingHeader::new(summary)
            .id_salt(format!("map-{}", def.id))
            .show(ui, |ui| {
                for key in keys {
                    let mut v = current.get(*key).copied().unwrap_or(0.0);
                    ui.horizontal(|ui| {
                        ui.label(egui::RichText::new(*key).monospace());
                        let speed = (v.abs() * 0.01).max(1e-4);
                        if ui
                            .add(egui::DragValue::new(&mut v).speed(speed).max_decimals(6))
                            .changed()
                        {
                            let mut m = current.clone();
                            if v > 0.0 {
                                m.insert((*key).to_string(), v);
                            } else {
                                m.remove(*key);
                            }
                            result = Some(m);
                        }
                    });
                }
            });
    });
    result
}

pub fn format_num(v: f64) -> String {
    exogen_core::validate::format_number(v)
}
