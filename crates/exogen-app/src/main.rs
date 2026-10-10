// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! ExoGen desktop application (PLAN.md section 7): parameter panel on the left, preview in
//! the centre with a time scrubber, pipeline and plausibility on the right, log at the
//! bottom. All computation runs on a worker thread.

mod params_panel;
mod preview;
mod worker;

use std::path::PathBuf;
use std::sync::Arc;

use eframe::egui;
use exogen_core::params::ParamGraph;
use exogen_core::project::Project;
use exogen_core::stages::{FinishOutput, PipelineResult, Progress, Runner};

use preview::{Channel, EquirectRgb, GlobeView};
use worker::{Job, Msg, Worker};

const PREVIEW_WIDTH: usize = 1024;
const PREVIEW_FACE: u32 = 256;
const GLOBE_SIZE: usize = 640;

fn main() -> eframe::Result {
    let options = eframe::NativeOptions {
        renderer: eframe::Renderer::Wgpu,
        viewport: egui::ViewportBuilder::default()
            .with_title("ExoGen")
            .with_inner_size([1500.0, 900.0]),
        ..Default::default()
    };
    eframe::run_native(
        "ExoGen",
        options,
        Box::new(|cc| Ok(Box::new(App::new(&cc.egui_ctx)))),
    )
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum View {
    Map,
    Globe,
}

struct App {
    graph: ParamGraph,
    project_name: String,
    project_path: String,
    preset_choice: String,
    expert: bool,
    search: String,

    worker: Worker,
    progress: Option<Progress>,
    log: Vec<String>,
    result: Option<Arc<PipelineResult>>,
    result_loop_key: u64,
    result_finish_key: u64,
    previews: std::collections::HashMap<usize, Arc<FinishOutput>>,
    pending_preview: Option<usize>,

    channel: Channel,
    view: View,
    keyframe: Option<usize>,
    globe: GlobeView,
    image: Option<Arc<EquirectRgb>>,
    image_key: (Channel, Option<usize>, u64),
    map_texture: Option<egui::TextureHandle>,
    globe_texture: Option<egui::TextureHandle>,
    globe_dirty: bool,

    export_dir: String,
    export_width: usize,
}

impl App {
    fn new(ctx: &egui::Context) -> Self {
        let ctx2 = ctx.clone();
        let worker = Worker::start(move || ctx2.request_repaint());
        let (graph, _) = exogen_core::presets::preset("Earth").unwrap().to_graph();
        let mut app = Self {
            graph,
            project_name: "Earth".to_string(),
            project_path: "project.toml".to_string(),
            preset_choice: "Earth".to_string(),
            expert: false,
            search: String::new(),
            worker,
            progress: None,
            log: Vec::new(),
            result: None,
            result_loop_key: 0,
            result_finish_key: 0,
            previews: Default::default(),
            pending_preview: None,
            channel: Channel::Albedo,
            view: View::Map,
            keyframe: None,
            globe: GlobeView::default(),
            image: None,
            image_key: (Channel::Albedo, None, 0),
            map_texture: None,
            globe_texture: None,
            globe_dirty: true,
            export_dir: "out".to_string(),
            export_width: 4096,
        };
        app.log("ExoGen started. Earth preset loaded. Press Run.");
        app
    }

    fn log(&mut self, s: impl Into<String>) {
        self.log.push(s.into());
        if self.log.len() > 500 {
            self.log.remove(0);
        }
    }

    fn run(&mut self) {
        if self.worker.busy {
            return;
        }
        self.progress = None;
        self.log(format!(
            "run: loop face {} x {} steps, finishing face {}",
            self.graph.f64("detail.loop_resolution"),
            exogen_core::stages::LOOP_STEPS[match self.graph.str("detail.quality") {
                "preview" => 0,
                "high" => 2,
                _ => 1,
            }],
            self.graph.f64("detail.finishing_resolution")
        ));
        self.worker.submit(Job::Run(self.graph.clone()));
    }

    fn stale(&self) -> bool {
        self.result.is_some()
            && (Runner::loop_key(&self.graph) != self.result_loop_key
                || Runner::finish_key(&self.graph, self.result.as_ref().unwrap().loop_output.hash)
                    != self.result_finish_key)
    }

    fn handle_messages(&mut self) {
        for m in self.worker.poll() {
            match m {
                Msg::Progress(p) => self.progress = Some(p),
                Msg::RunDone(r) => {
                    for (stage, d) in &r.timings {
                        self.log(format!("{stage}: {d:.2?}"));
                    }
                    if r.loop_cached {
                        self.log("coarse loop: cached");
                    }
                    if r.finish_cached {
                        self.log("finishing: cached");
                    }
                    self.log(format!(
                        "done: sea level {} m, liquid fraction {:.3}, {} keyframes",
                        params_panel::format_num(r.finish.sea_level),
                        r.finish.ocean_fraction,
                        r.loop_output.keyframes.len()
                    ));
                    self.result_loop_key = Runner::loop_key(&self.graph);
                    self.result_finish_key = Runner::finish_key(&self.graph, r.loop_output.hash);
                    self.graph.set_stage_output(
                        "volatiles.sea_level",
                        if r.finish.sea_level.is_finite() {
                            r.finish.sea_level
                        } else {
                            0.0
                        },
                    );
                    self.graph
                        .set_stage_output("volatiles.ocean_fraction", r.finish.ocean_fraction);
                    self.previews.clear();
                    self.keyframe = None;
                    self.result = Some(r);
                    self.progress = None;
                    self.image = None;
                }
                Msg::PreviewDone(index, f) => {
                    self.previews.insert(index, f);
                    if self.pending_preview == Some(index) {
                        self.pending_preview = None;
                    }
                    self.image = None;
                }
                Msg::ExportDone(s) => self.log(s),
                Msg::Error(e) => {
                    self.log(format!("error: {e}"));
                    self.progress = None;
                    self.pending_preview = None;
                }
            }
        }
    }

    /// The finished state to show: the final output, or a scrubbed keyframe preview.
    fn current_output(&mut self) -> Option<Arc<FinishOutput>> {
        let result = self.result.clone()?;
        match self.keyframe {
            None => Some(result.finish.clone()),
            Some(i) => {
                if let Some(f) = self.previews.get(&i) {
                    Some(f.clone())
                } else {
                    if self.pending_preview.is_none() && !self.worker.busy {
                        self.pending_preview = Some(i);
                        self.worker.submit(Job::Preview {
                            graph: self.graph.clone(),
                            loop_output: result.loop_output.clone(),
                            index: i,
                            face: PREVIEW_FACE,
                        });
                    }
                    None
                }
            }
        }
    }

    fn ensure_image(&mut self, ctx: &egui::Context) {
        let Some(out) = self.current_output() else {
            return;
        };
        let key = (self.channel, self.keyframe, out.hash);
        if self.image.is_some() && self.image_key == key {
            return;
        }
        let img = Arc::new(preview::build_equirect(&out, self.channel, PREVIEW_WIDTH));
        let color = preview::to_color_image(&img, self.globe.exposure);
        self.map_texture =
            Some(ctx.load_texture("preview-map", color, egui::TextureOptions::LINEAR));
        self.image = Some(img);
        self.image_key = key;
        self.globe_dirty = true;
    }

    fn ensure_globe(&mut self, ctx: &egui::Context) {
        if !self.globe_dirty {
            return;
        }
        let Some(img) = &self.image else { return };
        let color = preview::render_globe(img, &self.globe, GLOBE_SIZE);
        self.globe_texture =
            Some(ctx.load_texture("preview-globe", color, egui::TextureOptions::LINEAR));
        self.globe_dirty = false;
    }

    fn load_preset(&mut self, name: &str) {
        if let Some(p) = exogen_core::presets::preset(name) {
            let (g, _) = p.to_graph();
            self.graph = g;
            self.project_name = p.name.clone();
            self.log(format!("loaded preset {name}"));
        }
    }

    fn load_project(&mut self) {
        match Project::load(&PathBuf::from(&self.project_path)) {
            Ok(p) => {
                let (g, unknown) = p.to_graph();
                for id in unknown {
                    self.log(format!("unknown parameter {id} ignored"));
                }
                self.graph = g;
                self.project_name = p.name.clone();
                self.log(format!("loaded {}", self.project_path));
            }
            Err(e) => self.log(format!("error: {e}")),
        }
    }

    fn save_project(&mut self) {
        let p = Project::from_graph(self.project_name.clone(), &self.graph);
        match p.save(&PathBuf::from(&self.project_path)) {
            Ok(()) => self.log(format!("saved {}", self.project_path)),
            Err(e) => self.log(format!("error: {e}")),
        }
    }

    fn export(&mut self) {
        let Some(finish) = self.result.as_ref().map(|r| r.finish.clone()) else {
            self.log("nothing to export: run first");
            return;
        };
        if self.worker.busy {
            return;
        }
        self.log(format!(
            "export: {}x{} to {}",
            self.export_width,
            self.export_width / 2,
            self.export_dir
        ));
        self.worker.submit(Job::Export {
            graph: self.graph.clone(),
            finish,
            dir: PathBuf::from(&self.export_dir),
            width: self.export_width,
            name: self.project_name.clone(),
        });
    }

    // ------------------------------------------------------------------ panels

    fn left_panel(&mut self, ui: &mut egui::Ui) {
        ui.horizontal(|ui| {
            ui.label("Preset");
            let mut choice = self.preset_choice.clone();
            egui::ComboBox::from_id_salt("preset")
                .selected_text(choice.clone())
                .show_ui(ui, |ui| {
                    for name in exogen_core::presets::PRESET_NAMES {
                        ui.selectable_value(&mut choice, (*name).to_string(), *name);
                    }
                });
            if choice != self.preset_choice {
                self.preset_choice = choice.clone();
                self.load_preset(&choice);
            }
            ui.checkbox(&mut self.expert, "expert");
        });
        ui.horizontal(|ui| {
            ui.label("Find");
            ui.add(egui::TextEdit::singleline(&mut self.search).desired_width(160.0));
            if ui.small_button("clear").clicked() {
                self.search.clear();
            }
        });
        ui.separator();
        egui::ScrollArea::vertical()
            .auto_shrink([false, false])
            .show(ui, |ui| {
                params_panel::show(ui, &mut self.graph, self.expert, &self.search);
            });
    }

    fn right_panel(&mut self, ui: &mut egui::Ui) {
        ui.heading("Project");
        ui.horizontal(|ui| {
            ui.label("Name");
            ui.add(egui::TextEdit::singleline(&mut self.project_name).desired_width(150.0));
        });
        ui.horizontal(|ui| {
            ui.label("File");
            ui.add(egui::TextEdit::singleline(&mut self.project_path).desired_width(150.0));
        });
        ui.horizontal(|ui| {
            if ui.button("Load").clicked() {
                self.load_project();
            }
            if ui.button("Save").clicked() {
                self.save_project();
            }
        });
        ui.separator();

        ui.heading("Pipeline");
        let stale = self.stale();
        ui.horizontal(|ui| {
            let run_label = if self.result.is_none() {
                "Run"
            } else if stale {
                "Run (parameters changed)"
            } else {
                "Run again"
            };
            if ui
                .add_enabled(!self.worker.busy, egui::Button::new(run_label))
                .clicked()
            {
                self.run();
            }
            if ui
                .add_enabled(self.worker.busy, egui::Button::new("Cancel"))
                .clicked()
            {
                self.worker.cancel();
            }
        });
        if let Some(p) = &self.progress {
            ui.add(egui::ProgressBar::new(p.fraction).text(format!("{}: {}", p.stage, p.message)));
        }
        egui::Grid::new("stages")
            .num_columns(2)
            .spacing([12.0, 2.0])
            .show(ui, |ui| {
                let loop_state = match &self.result {
                    None => "not run".to_string(),
                    Some(_) => {
                        if Runner::loop_key(&self.graph) != self.result_loop_key {
                            "stale".to_string()
                        } else {
                            "current".to_string()
                        }
                    }
                };
                ui.label("coarse loop");
                ui.label(loop_state);
                ui.end_row();
                let finish_state = match &self.result {
                    None => "not run".to_string(),
                    Some(r) => {
                        if Runner::finish_key(&self.graph, r.loop_output.hash)
                            != self.result_finish_key
                        {
                            "stale".to_string()
                        } else {
                            "current".to_string()
                        }
                    }
                };
                ui.label("finishing");
                ui.label(finish_state);
                ui.end_row();
                ui.label("export");
                ui.label("on request");
                ui.end_row();
            });
        ui.separator();

        ui.heading("Export");
        ui.horizontal(|ui| {
            ui.label("Directory");
            ui.add(egui::TextEdit::singleline(&mut self.export_dir).desired_width(150.0));
        });
        ui.horizontal(|ui| {
            ui.label("Width");
            egui::ComboBox::from_id_salt("export-width")
                .selected_text(self.export_width.to_string())
                .show_ui(ui, |ui| {
                    for w in [1024usize, 2048, 4096, 8192, 16384] {
                        ui.selectable_value(&mut self.export_width, w, w.to_string());
                    }
                });
            if ui
                .add_enabled(
                    self.result.is_some() && !self.worker.busy,
                    egui::Button::new("Export"),
                )
                .clicked()
            {
                self.export();
            }
        });
        ui.separator();

        ui.heading("Plausibility");
        let warnings = self.graph.warnings();
        if warnings.is_empty() {
            ui.label(egui::RichText::new("No notes.").weak());
        } else {
            egui::ScrollArea::vertical()
                .max_height(220.0)
                .show(ui, |ui| {
                    for w in warnings {
                        let prefix = match w.severity {
                            exogen_core::params::rules::Severity::Warning => "Warning",
                            exogen_core::params::rules::Severity::Note => "Note",
                        };
                        ui.label(format!("{prefix}: {}", w.text))
                            .on_hover_text(w.node);
                    }
                });
        }
    }

    fn center_panel(&mut self, ui: &mut egui::Ui) {
        let ctx = ui.ctx().clone();
        ui.horizontal(|ui| {
            ui.selectable_value(&mut self.view, View::Map, "Map");
            ui.selectable_value(&mut self.view, View::Globe, "Globe");
            ui.separator();
            ui.label("Channel");
            let mut ch = self.channel;
            egui::ComboBox::from_id_salt("channel")
                .selected_text(ch.label())
                .show_ui(ui, |ui| {
                    for c in Channel::ALL {
                        ui.selectable_value(&mut ch, c, c.label());
                    }
                });
            self.channel = ch;
            if self.view == View::Globe {
                ui.separator();
                ui.label("Sun az");
                if ui
                    .add(
                        egui::Slider::new(&mut self.globe.sun_azimuth, -3.14..=3.14)
                            .show_value(false),
                    )
                    .changed()
                {
                    self.globe_dirty = true;
                }
                ui.label("el");
                if ui
                    .add(
                        egui::Slider::new(&mut self.globe.sun_elevation, -1.5..=1.5)
                            .show_value(false),
                    )
                    .changed()
                {
                    self.globe_dirty = true;
                }
            }
            if let Some(img) = &self.image {
                if self.channel != Channel::Albedo && self.channel != Channel::Unit {
                    ui.separator();
                    ui.label(
                        egui::RichText::new(format!(
                            "range {} to {}",
                            params_panel::format_num(img.range.0 as f64),
                            params_panel::format_num(img.range.1 as f64)
                        ))
                        .weak(),
                    );
                }
            }
        });

        // Time scrubber.
        if let Some(r) = self.result.clone() {
            let n = r.loop_output.keyframes.len();
            ui.horizontal(|ui| {
                ui.label("Time");
                let mut idx = self.keyframe.unwrap_or(n - 1);
                let resp = ui.add(egui::Slider::new(&mut idx, 0..=n - 1).show_value(false));
                let t = r.loop_output.keyframes[idx].time_gyr;
                ui.label(egui::RichText::new(format!("{t:.2} Gyr after formation")).monospace());
                if resp.changed() {
                    self.keyframe = if idx == n - 1 { None } else { Some(idx) };
                }
                if self.keyframe.is_some() {
                    ui.label(
                        egui::RichText::new(format!("preview at face size {PREVIEW_FACE}")).weak(),
                    );
                    if ui.small_button("present").clicked() {
                        self.keyframe = None;
                    }
                } else {
                    ui.label(
                        egui::RichText::new(format!(
                            "final output, face size {}",
                            r.finish.grid.face_size()
                        ))
                        .weak(),
                    );
                }
            });
        }

        self.ensure_image(&ctx);
        if self.view == View::Globe {
            self.ensure_globe(&ctx);
        }

        let avail = ui.available_size();
        match self.view {
            View::Map => {
                if let Some(tex) = &self.map_texture {
                    let w = avail.x.min(avail.y * 2.0).max(64.0);
                    let sized = egui::load::SizedTexture::new(tex.id(), egui::vec2(w, w / 2.0));
                    let resp = ui.add(egui::Image::from_texture(sized).sense(egui::Sense::hover()));
                    if let (Some(pos), Some(out)) = (resp.hover_pos(), self.current_output()) {
                        let rect = resp.rect;
                        let fx = ((pos.x - rect.left()) / rect.width()).clamp(0.0, 1.0) as f64;
                        let fy = ((pos.y - rect.top()) / rect.height()).clamp(0.0, 1.0) as f64;
                        let lon = -180.0 + 360.0 * fx;
                        let lat = 90.0 - 180.0 * fy;
                        let d = exogen_core::grid::direction_from_lon_lat(
                            lon.to_radians(),
                            lat.to_radians(),
                        );
                        let cell = out.grid.locate(d).0;
                        let h = out.height.get(cell);
                        let liquid = out.layers.liquid.get(cell);
                        let age = out.layers.surface_age.get(cell);
                        let unit = out.layers.unit_id.get(cell);
                        resp.on_hover_text(format!(
                            "lon {lon:.1}  lat {lat:.1}\nheight {:.0} m  liquid {:.0} m\nsurface age {age:.2} Gyr  unit {unit}",
                            h, liquid
                        ));
                    }
                } else {
                    self.empty_hint(ui);
                }
            }
            View::Globe => {
                if let Some(tex) = &self.globe_texture {
                    let s = avail.x.min(avail.y).max(64.0);
                    let sized = egui::load::SizedTexture::new(tex.id(), egui::vec2(s, s));
                    let resp = ui
                        .add(egui::Image::from_texture(sized).sense(egui::Sense::click_and_drag()));
                    if resp.dragged() {
                        let d = resp.drag_delta();
                        self.globe.yaw -= d.x * 0.01;
                        self.globe.pitch = (self.globe.pitch + d.y * 0.01).clamp(-1.5, 1.5);
                        self.globe_dirty = true;
                    }
                    if resp.hovered() {
                        let scroll = ui.input(|i| i.smooth_scroll_delta.y);
                        if scroll != 0.0 {
                            self.globe.zoom =
                                (self.globe.zoom * (1.0 + scroll * 0.002)).clamp(0.5, 8.0);
                            self.globe_dirty = true;
                        }
                    }
                } else {
                    self.empty_hint(ui);
                }
            }
        }
    }

    fn empty_hint(&self, ui: &mut egui::Ui) {
        ui.centered_and_justified(|ui| {
            let text = if self.worker.busy {
                "Running."
            } else if self.result.is_none() {
                "No output yet. Press Run."
            } else {
                "Preparing preview."
            };
            ui.label(egui::RichText::new(text).weak());
        });
    }

    fn bottom_panel(&mut self, ui: &mut egui::Ui) {
        egui::ScrollArea::vertical()
            .stick_to_bottom(true)
            .auto_shrink([false, false])
            .show(ui, |ui| {
                for line in &self.log {
                    ui.label(egui::RichText::new(line).monospace().small());
                }
            });
    }
}

impl eframe::App for App {
    fn ui(&mut self, ui: &mut egui::Ui, _frame: &mut eframe::Frame) {
        self.handle_messages();
        if self.worker.busy {
            ui.ctx()
                .request_repaint_after(std::time::Duration::from_millis(100));
        }
        egui::Panel::left("params")
            .default_size(420.0)
            .show(ui, |ui| self.left_panel(ui));
        egui::Panel::right("pipeline")
            .default_size(300.0)
            .show(ui, |ui| self.right_panel(ui));
        egui::Panel::bottom("log")
            .default_size(110.0)
            .show(ui, |ui| self.bottom_panel(ui));
        egui::CentralPanel::default().show(ui, |ui| self.center_panel(ui));
    }
}
