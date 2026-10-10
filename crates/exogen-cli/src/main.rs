// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! ExoGen command-line tool: derivation reports, validation, presets, rendering and,
//! during the early phases, the spikes.

use std::path::PathBuf;

use clap::{Parser, Subcommand};
use exogen_core::params::{Kind, ParamGraph, Value, ValueType, Visibility};
use exogen_core::project::Project;
use exogen_core::validate::format_number;

#[derive(Parser)]
#[command(
    name = "exogen",
    version,
    about = "Process-based planetary map generator"
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Print every parameter of a project with its derived value, unit and state.
    Derive {
        /// Project file. Omit to use a preset or the defaults.
        project: Option<PathBuf>,
        /// Start from a named preset instead of a file.
        #[arg(long)]
        preset: Option<String>,
        /// Only this group (star, orbit, body, rotation, interior, atmosphere, volatiles, ...).
        #[arg(long)]
        group: Option<String>,
        /// Include expert nodes.
        #[arg(long)]
        expert: bool,
    },
    /// Run the pipeline for a project and write the map set.
    Render {
        /// Project file. Omit to use a preset or the defaults.
        project: Option<PathBuf>,
        /// Start from a named preset instead of a file.
        #[arg(long)]
        preset: Option<String>,
        /// Output directory.
        #[arg(long, default_value = "out")]
        out: PathBuf,
        /// Equirectangular width in pixels; overrides output.width.
        #[arg(long)]
        width: Option<usize>,
        /// preview, standard or high; overrides detail.quality.
        #[arg(long)]
        quality: Option<String>,
        /// Gyr before present to render; overrides output.epoch.
        #[arg(long)]
        epoch: Option<f64>,
    },
    /// Check the derived values of the reference bodies against measured values.
    Validate,
    /// List, show or copy presets.
    Preset {
        #[command(subcommand)]
        action: PresetAction,
    },
    /// Phase 0 grid spike: fill a cubesphere with an analytic field, resample it to an
    /// equirectangular image, report the error and write PNGs.
    SpikeGrid {
        #[arg(long, default_value_t = 256)]
        face_size: u32,
        #[arg(long, default_value_t = 2048)]
        width: usize,
        #[arg(long, default_value = "out")]
        out: PathBuf,
    },
}

#[derive(Subcommand)]
enum PresetAction {
    /// List preset names.
    List,
    /// Print a preset's project file.
    Show { name: String },
    /// Write a preset to a project file to start from.
    New { name: String, out: PathBuf },
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    match Cli::parse().command {
        Command::Derive {
            project,
            preset,
            group,
            expert,
        } => derive(project, preset, group, expert),
        Command::Render {
            project,
            preset,
            out,
            width,
            quality,
            epoch,
        } => render(project, preset, out, width, quality, epoch),
        Command::Validate => validate(),
        Command::Preset { action } => preset_cmd(action),
        Command::SpikeGrid {
            face_size,
            width,
            out,
        } => spike_grid(face_size, width, &out),
    }
}

fn load_graph(
    project: Option<PathBuf>,
    preset: Option<String>,
) -> Result<(ParamGraph, String), Box<dyn std::error::Error>> {
    if let Some(path) = project {
        let p = Project::load(&path)?;
        let (g, unknown) = p.to_graph();
        for id in unknown {
            eprintln!("warning: unknown parameter {id} ignored");
        }
        return Ok((g, p.name));
    }
    if let Some(name) = preset {
        let p = exogen_core::presets::preset(&name)
            .ok_or_else(|| format!("no preset named {name}; run `exogen preset list`"))?;
        let (g, _) = p.to_graph();
        return Ok((g, p.name));
    }
    Ok((ParamGraph::new(), "defaults".to_string()))
}

fn format_value(g: &ParamGraph, id: &str) -> String {
    let def = g.def(id);
    match g.get(id) {
        Value::F64(v) => {
            let shown = def.display(*v);
            if def.unit.is_empty() {
                format_number(shown)
            } else {
                format!("{} {}", format_number(shown), def.unit)
            }
        }
        Value::Bool(b) => if *b { "yes" } else { "no" }.to_string(),
        Value::Str(s) => {
            if s.is_empty() {
                "-".to_string()
            } else {
                s.clone()
            }
        }
        Value::Curve(c) => {
            if c.is_constant() {
                format!("constant {}", format_number(c.points[0].1))
            } else {
                let first = c.points[0];
                let last = c.points[c.points.len() - 1];
                format!(
                    "curve {} at {:.2} Gyr to {} at {:.2} Gyr",
                    format_number(first.1),
                    first.0,
                    format_number(last.1),
                    last.0
                )
            }
        }
        Value::Map(m) => {
            if m.is_empty() {
                "-".to_string()
            } else {
                m.iter()
                    .map(|(k, v)| format!("{k} {}", format_number(*v)))
                    .collect::<Vec<_>>()
                    .join(", ")
            }
        }
        Value::List(l) => {
            if l.is_empty() {
                "-".to_string()
            } else {
                l.join(", ")
            }
        }
    }
}

fn derive(
    project: Option<PathBuf>,
    preset: Option<String>,
    group: Option<String>,
    expert: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    let (g, name) = load_graph(project, preset)?;
    println!("{name}");
    for grp in g.groups() {
        if group.as_deref().is_some_and(|w| w != grp) {
            continue;
        }
        println!("\n[{grp}]");
        for def in g.nodes_in(grp) {
            if def.visibility == Visibility::Expert && !expert {
                continue;
            }
            let state = match def.kind {
                Kind::Root => {
                    if g.is_modified(def.id) {
                        "set"
                    } else {
                        "default"
                    }
                }
                Kind::Derived => {
                    if g.is_pinned(def.id) {
                        "pinned"
                    } else {
                        "derived"
                    }
                }
                Kind::Multiplier => "multiplier",
                Kind::StageOutput => "stage output",
            };
            let value = format_value(&g, def.id);
            let kind_note = match def.ty {
                ValueType::Enum(_) | ValueType::Text | ValueType::Bool => String::new(),
                _ => String::new(),
            };
            println!("  {:<38} {:<42} {:<12}{}", def.id, value, state, kind_note);
        }
    }
    let warnings = g.warnings();
    if !warnings.is_empty() {
        println!("\n[plausibility]");
        for w in warnings {
            println!("  {:?}: {} ({})", w.severity, w.text, w.node);
        }
    }
    Ok(())
}

struct ConsoleProgress {
    last: std::sync::Mutex<(String, f32)>,
}

impl exogen_core::stages::ProgressSink for ConsoleProgress {
    fn report(&self, p: exogen_core::stages::Progress) {
        let mut last = self.last.lock().unwrap();
        if last.0 != p.stage || p.fraction - last.1 >= 0.1 || p.fraction >= 1.0 {
            println!(
                "  {:<12} {:3.0}%  {}",
                p.stage,
                p.fraction * 100.0,
                p.message
            );
            *last = (p.stage.to_string(), p.fraction);
        }
    }
}

fn render(
    project: Option<PathBuf>,
    preset: Option<String>,
    out: PathBuf,
    width: Option<usize>,
    quality: Option<String>,
    epoch: Option<f64>,
) -> Result<(), Box<dyn std::error::Error>> {
    let (mut g, name) = load_graph(project, preset)?;
    if let Some(q) = quality {
        g.set("detail.quality", q.as_str())?;
    }
    if let Some(e) = epoch {
        g.set("output.epoch", e)?;
    }
    if let Some(w) = width {
        g.set("output.width", w as f64)?;
    }
    for w in g.warnings() {
        println!("note: {}", w.text);
    }
    println!(
        "{name}: loop face {} x {} steps, finishing face {}",
        g.f64("detail.loop_resolution"),
        exogen_core::stages::LOOP_STEPS[match g.str("detail.quality") {
            "preview" => 0,
            "high" => 2,
            _ => 1,
        }],
        g.f64("detail.finishing_resolution")
    );
    let sink = ConsoleProgress {
        last: std::sync::Mutex::new((String::new(), -1.0)),
    };
    let mut runner = exogen_core::stages::Runner::new();
    let t0 = std::time::Instant::now();
    let result = runner.run(&g, &sink)?;
    for (stage, d) in &result.timings {
        println!("  {stage}: {d:.2?}");
    }
    let settings = exogen_core::export::ExportSettings {
        width: g.f64("output.width") as usize,
        dir: out.clone(),
        project_name: name,
    };
    let t1 = std::time::Instant::now();
    let sidecar = exogen_core::export::export(&g, &result.finish, &settings)?;
    println!("  export: {:.2?}", t1.elapsed());
    println!(
        "wrote {} maps to {} ({}x{}), sea level {} m, liquid fraction {:.3}, total {:.2?}",
        sidecar.maps.len(),
        out.display(),
        sidecar.width,
        sidecar.height,
        format_number(result.finish.sea_level),
        result.finish.ocean_fraction,
        t0.elapsed()
    );
    Ok(())
}

fn validate() -> Result<(), Box<dyn std::error::Error>> {
    let results = exogen_core::validate::run_reference_checks();
    let mut failed = 0usize;
    let mut current = "";
    for r in &results {
        if r.check.preset != current {
            current = r.check.preset;
            println!("\n{current}");
        }
        let mark = if r.pass { "ok  " } else { "FAIL" };
        let expected = match &r.check.expected {
            exogen_core::validate::Expected::Range(lo, hi) => {
                format!("{} to {}", format_number(*lo), format_number(*hi))
            }
            exogen_core::validate::Expected::Str(s) => s.to_string(),
            exogen_core::validate::Expected::Bool(b) => b.to_string(),
        };
        println!(
            "  {mark} {:<40} {:<14} expected {:<22} {}",
            r.check.node, r.actual, expected, r.check.note
        );
        if !r.pass {
            failed += 1;
        }
    }
    println!("\n{} checks, {} failed", results.len(), failed);
    if failed > 0 {
        std::process::exit(1);
    }
    Ok(())
}

fn preset_cmd(action: PresetAction) -> Result<(), Box<dyn std::error::Error>> {
    match action {
        PresetAction::List => {
            for name in exogen_core::presets::PRESET_NAMES {
                let p = exogen_core::presets::preset(name).unwrap();
                println!("{name:<10} {}", p.description);
            }
        }
        PresetAction::Show { name } => {
            let p = exogen_core::presets::preset(&name)
                .ok_or_else(|| format!("no preset named {name}"))?;
            print!("{}", p.to_toml()?);
        }
        PresetAction::New { name, out } => {
            let mut p = exogen_core::presets::preset(&name)
                .ok_or_else(|| format!("no preset named {name}"))?;
            p.name = format!("{name} copy");
            p.save(&out)?;
            println!("wrote {}", out.display());
        }
    }
    Ok(())
}

fn spike_grid(
    face_size: u32,
    width: usize,
    out: &PathBuf,
) -> Result<(), Box<dyn std::error::Error>> {
    use exogen_core::grid::projection::{
        angle_to_nearest_cube_corner, test_field, to_equirectangular,
    };
    use exogen_core::grid::{CubeSphere, direction_from_lon_lat};
    use std::f64::consts::{FRAC_PI_2, PI};

    let t0 = std::time::Instant::now();
    let grid = CubeSphere::new(face_size);
    let values: Vec<f32> = grid
        .cells()
        .map(|c| test_field(grid.center(c)) as f32)
        .collect();
    let total_area: f64 = grid.cells().map(|c| grid.cell_area(c)).sum();
    println!(
        "cubesphere f={face_size}: {} cells, total solid angle {:.6} (4 pi = {:.6}), fill {:?}",
        grid.cell_count(),
        total_area,
        4.0 * PI,
        t0.elapsed()
    );

    let height = width / 2;
    let t1 = std::time::Instant::now();
    let img = to_equirectangular(&grid, &values, width, height);
    println!("resampled to {width}x{height} in {:?}", t1.elapsed());

    let bands = 9usize;
    let mut band_max = vec![0.0f64; bands];
    let mut band_sq = vec![0.0f64; bands];
    let mut band_n = vec![0usize; bands];
    let mut class_max = [0.0f64; 3];
    let mut errors = vec![0.0f32; width * height];
    let delta = grid.delta();
    for py in 0..height {
        let lat = FRAC_PI_2 - PI * (py as f64 + 0.5) / height as f64;
        let band = (((lat.abs() / FRAC_PI_2) * bands as f64) as usize).min(bands - 1);
        for px in 0..width {
            let lon = -PI + 2.0 * PI * (px as f64 + 0.5) / width as f64;
            let d = direction_from_lon_lat(lon, lat);
            let err = (img[py * width + px] as f64 - test_field(d)).abs();
            errors[py * width + px] = err as f32;
            band_max[band] = band_max[band].max(err);
            band_sq[band] += err * err;
            band_n[band] += 1;
            let (c, _, _) = grid.locate(d);
            let last = face_size - 1;
            let class = if angle_to_nearest_cube_corner(d) < 1.5 * delta {
                2
            } else if c.i == 0 || c.j == 0 || c.i == last || c.j == last {
                1
            } else {
                0
            };
            class_max[class] = class_max[class].max(err);
        }
    }
    println!(
        "max error: face interior {:.2e}, edge cells {:.2e}, corner cells {:.2e}",
        class_max[0], class_max[1], class_max[2]
    );
    println!("error by |latitude| band (degrees): max / rms");
    for b in 0..bands {
        let lo = 90.0 * b as f64 / bands as f64;
        let hi = 90.0 * (b + 1) as f64 / bands as f64;
        println!(
            "  {lo:4.0}-{hi:3.0}: {:.2e} / {:.2e}",
            band_max[b],
            (band_sq[b] / band_n[b].max(1) as f64).sqrt()
        );
    }

    std::fs::create_dir_all(out)?;
    write_gray_png(&out.join("spike_grid_field.png"), &img, width, height)?;
    write_gray_png(&out.join("spike_grid_error.png"), &errors, width, height)?;
    println!("wrote {}", out.display());
    Ok(())
}

fn write_gray_png(
    path: &std::path::Path,
    data: &[f32],
    width: usize,
    height: usize,
) -> Result<(), Box<dyn std::error::Error>> {
    let (mut lo, mut hi) = (f32::INFINITY, f32::NEG_INFINITY);
    for &v in data {
        lo = lo.min(v);
        hi = hi.max(v);
    }
    let scale = if hi > lo { 65535.0 / (hi - lo) } else { 0.0 };
    let pixels: Vec<u16> = data
        .iter()
        .map(|&v| ((v - lo) * scale).round() as u16)
        .collect();
    let img = image::ImageBuffer::<image::Luma<u16>, Vec<u16>>::from_raw(
        width as u32,
        height as u32,
        pixels,
    )
    .ok_or("image buffer size mismatch")?;
    img.save(path)?;
    Ok(())
}
