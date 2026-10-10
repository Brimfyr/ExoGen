// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! ExoGen command-line tool. Headless rendering, derivation reports, validation and,
//! during Phase 0, the spikes.

use std::path::PathBuf;

use clap::{Parser, Subcommand};

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
    /// Phase 0 grid spike: fill a cubesphere with an analytic field, resample it to an
    /// equirectangular image, report the error and write PNGs.
    SpikeGrid {
        /// Cells along each cube face edge.
        #[arg(long, default_value_t = 256)]
        face_size: u32,
        /// Output image width in pixels; height is half the width.
        #[arg(long, default_value_t = 2048)]
        width: usize,
        /// Output directory.
        #[arg(long, default_value = "out")]
        out: PathBuf,
    },
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    match Cli::parse().command {
        Command::SpikeGrid {
            face_size,
            width,
            out,
        } => spike_grid(face_size, width, &out),
    }
}

fn spike_grid(
    face_size: u32,
    width: usize,
    out: &PathBuf,
) -> Result<(), Box<dyn std::error::Error>> {
    use exogen_core::grid::projection::{test_field, to_equirectangular};
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

    // Error statistics by latitude band, and by position relative to face edges and corners.
    use exogen_core::grid::projection::angle_to_nearest_cube_corner;
    let bands = 9usize;
    let mut band_max = vec![0.0f64; bands];
    let mut band_sq = vec![0.0f64; bands];
    let mut band_n = vec![0usize; bands];
    let mut class_max = [0.0f64; 3]; // interior, edge, corner
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
