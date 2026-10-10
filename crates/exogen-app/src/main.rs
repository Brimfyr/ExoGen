// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! ExoGen desktop application. During Phase 0 this is the grid spike: one cubesphere
//! field resampled to an equirectangular image and shown as a texture through egui on wgpu.

use eframe::egui;
use exogen_core::grid::CubeSphere;
use exogen_core::grid::projection::{test_field, to_equirectangular};

fn main() -> eframe::Result {
    let options = eframe::NativeOptions {
        renderer: eframe::Renderer::Wgpu,
        viewport: egui::ViewportBuilder::default()
            .with_title("ExoGen grid spike")
            .with_inner_size([1100.0, 620.0]),
        ..Default::default()
    };
    eframe::run_native(
        "ExoGen grid spike",
        options,
        Box::new(|cc| Ok(Box::new(App::new(&cc.egui_ctx)))),
    )
}

struct App {
    face_size: u32,
    width: usize,
    texture: Option<egui::TextureHandle>,
    last_build_ms: f64,
}

impl App {
    fn new(ctx: &egui::Context) -> Self {
        let mut app = Self {
            face_size: 256,
            width: 1024,
            texture: None,
            last_build_ms: 0.0,
        };
        app.rebuild(ctx);
        app
    }

    fn rebuild(&mut self, ctx: &egui::Context) {
        let t0 = std::time::Instant::now();
        let grid = CubeSphere::new(self.face_size);
        let values: Vec<f32> = grid
            .cells()
            .map(|c| test_field(grid.center(c)) as f32)
            .collect();
        let height = self.width / 2;
        let img = to_equirectangular(&grid, &values, self.width, height);
        let (mut lo, mut hi) = (f32::INFINITY, f32::NEG_INFINITY);
        for &v in &img {
            lo = lo.min(v);
            hi = hi.max(v);
        }
        let scale = if hi > lo { 255.0 / (hi - lo) } else { 0.0 };
        let gray: Vec<u8> = img
            .iter()
            .map(|&v| ((v - lo) * scale).round() as u8)
            .collect();
        let color = egui::ColorImage::from_gray([self.width, height], &gray);
        self.texture = Some(ctx.load_texture("grid-spike", color, egui::TextureOptions::LINEAR));
        self.last_build_ms = t0.elapsed().as_secs_f64() * 1000.0;
    }
}

impl eframe::App for App {
    fn ui(&mut self, ui: &mut egui::Ui, _frame: &mut eframe::Frame) {
        let ctx = ui.ctx().clone();
        egui::Panel::top("controls").show(ui, |ui| {
            ui.horizontal(|ui| {
                ui.label("Face size (cells)");
                let mut f = self.face_size;
                egui::ComboBox::from_id_salt("face-size")
                    .selected_text(f.to_string())
                    .show_ui(ui, |ui| {
                        for v in [64u32, 128, 256, 512, 1024] {
                            ui.selectable_value(&mut f, v, v.to_string());
                        }
                    });
                if f != self.face_size {
                    self.face_size = f;
                    self.rebuild(&ctx);
                }
                ui.separator();
                ui.label(format!("Build {:.1} ms", self.last_build_ms));
                ui.separator();
                ui.label(format!(
                    "Cells {}",
                    6 * self.face_size as usize * self.face_size as usize
                ));
            });
        });
        egui::CentralPanel::default().show(ui, |ui| {
            if let Some(tex) = &self.texture {
                let avail = ui.available_size();
                let w = avail.x.min(avail.y * 2.0);
                let sized = egui::load::SizedTexture::new(tex.id(), egui::vec2(w, w / 2.0));
                ui.add(egui::Image::from_texture(sized));
            }
        });
    }
}
