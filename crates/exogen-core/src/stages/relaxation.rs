// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! Viscous relaxation of an ice shell's topography (PLAN.md 4.3 and Phase 4).
//!
//! Ice flows. A load of wavelength λ on a shell relaxes on a timescale τ = 4π η / (ρ g λ)
//! with the viscosity of the ice at the depth the load reaches into, about a quarter of
//! the wavelength, and that ice is the warmer the deeper, along the shell's conductive
//! profile, down to the melting point at its base. So the longest wavelengths go first
//! and fastest: a basin a few hundred kilometres across flattens within a million years
//! on a thin warm shell and leaves a palimpsest, a bright circular ghost with its rings,
//! while a crater of twenty kilometres keeps its depth for the age of the solar system.
//! The loop applies the decay by wavelength band on the coarse grid at every step, with
//! the shell as the thermal history has it then, so the ancient craters of a shell that
//! was once thin are the flat ones; the finish restores the resolved rims on the flattened
//! forms, and relaxes the replayed craters by the same law over their own lives.

use crate::constants::*;
use crate::fields::Field;
use crate::grid::CubeSphere;

/// Viscosity of ice at a temperature, Pa s: an Arrhenius law about the melting point
/// (Goldsby and Kohlstedt 2001, with the basal viscosity of a convecting shell).
pub fn ice_viscosity(t: f64) -> f64 {
    let t = t.max(ICE_VISCOSITY_FLOOR_K);
    ICE_VISCOSITY_MELT
        * (ICE_CREEP_ACTIVATION / GAS_CONSTANT * (1.0 / t - 1.0 / ICE_MELTING_POINT)).exp()
}

/// The relaxation timescale of a wavelength on a shell, Gyr: the ice at the depth the
/// load reaches, no deeper than the shell, along a conductive profile from the surface
/// temperature to the melting point.
pub fn relaxation_time_gyr(
    wavelength: f64,
    shell_thickness: f64,
    surface_temperature: f64,
    gravity: f64,
) -> f64 {
    let depth = (RELAXATION_DEPTH_FRACTION * wavelength).min(shell_thickness.max(1.0));
    let frac = (depth / shell_thickness.max(1.0)).clamp(0.0, 1.0);
    // The conductivity of ice falls as 1/T, so the conductive profile is exponential in
    // depth rather than linear: T(z) = T_s (T_m / T_s)^(z/d).
    let t = surface_temperature * (ICE_MELTING_POINT / surface_temperature.max(1.0)).powf(frac);
    let eta = ice_viscosity(t);
    let tau_s = 4.0 * std::f64::consts::PI * eta / (DENSITY_CRUST_ICE * gravity * wavelength.max(1.0));
    tau_s / SECONDS_PER_GYR
}

/// The fraction of a wavelength's relief that survives a time on a shell.
pub fn survival(wavelength: f64, dt_gyr: f64, shell: f64, surface_temperature: f64, gravity: f64) -> f64 {
    let tau = relaxation_time_gyr(wavelength, shell, surface_temperature, gravity);
    if tau <= 0.0 {
        0.0
    } else {
        (-dt_gyr / tau).exp()
    }
}

/// A pyramid of low-passed copies of a field, each level smoothing over twice the
/// scale of the last, by resampling down and up again.
fn pyramid(field: &Field<f32>, levels: usize) -> Vec<Field<f32>> {
    let grid = field.grid;
    let mut out = Vec::with_capacity(levels);
    let mut face = grid.face_size();
    for _ in 0..levels {
        face = (face / 2).max(4);
        let coarse = CubeSphere::new(face);
        out.push(field.resample(coarse).resample(grid));
    }
    out
}

/// Relax a height field over a step: each wavelength band, from twice the cell up to
/// the whole body, keeps the share of its relief the shell's viscosity lets it keep.
pub fn relax(
    height: &mut Field<f32>,
    body_radius: f64,
    dt_gyr: f64,
    shell: f64,
    surface_temperature: f64,
    gravity: f64,
) -> bool {
    let grid = height.grid;
    let cell = grid.delta() * body_radius;
    let levels = ((grid.face_size() as f64).log2().floor() as usize).saturating_sub(2).max(1);
    // Survival per band, from the finest the pyramid resolves to the longest. Nothing to
    // do when even the longest wavelength keeps everything.
    let factors: Vec<f32> = (0..=levels)
        .map(|l| {
            let wavelength = 2.0 * cell * (2u32.pow(l as u32 + 1) as f64);
            survival(wavelength, dt_gyr, shell, surface_temperature, gravity) as f32
        })
        .collect();
    if factors.iter().all(|&f| f > 0.999) {
        return false;
    }
    let lows = pyramid(height, levels);
    let h = height.as_mut_slice();
    let n = h.len();
    // Bands: h - low_0, low_0 - low_1, ..., low_{L-1}; the last is the global shape.
    for k in 0..n {
        let mut v = 0.0f32;
        let mut prev = h[k];
        for (l, low) in lows.iter().enumerate() {
            let band = prev - low.as_slice()[k];
            v += band * factors[l];
            prev = low.as_slice()[k];
        }
        v += prev * factors[levels];
        h[k] = v;
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn warm_shells_relax_long_wavelengths_first() {
        // Ganymede today: a 100 km shell at 106 K. A 300 km load goes within a million
        // years, a 20 km crater outlasts the solar system.
        let (shell, t_s, g) = (83e3, 106.0, 1.43);
        assert!(relaxation_time_gyr(300e3, shell, t_s, g) < 1e-3);
        assert!(relaxation_time_gyr(20e3, shell, t_s, g) > 100.0);
        // Europa's thin shell relaxes even 30 km craters within its surface age.
        assert!(relaxation_time_gyr(30e3, 8e3, 92.0, 1.31) < 0.05);
        // Relaxation keeps a small bump and flattens a broad one.
        let grid = CubeSphere::new(64);
        let r = 2.6e6;
        let mut h = Field::from_fn(grid, |_, d| {
            let broad = (d[0] * 1.0 - 1.0).abs() < 0.2;
            let fine = ((d[1] * 200.0).sin() * (d[2] * 200.0).sin()) as f32;
            if broad { 1000.0 + 10.0 * fine } else { 10.0 * fine }
        });
        let before: f32 = h.as_slice().iter().map(|v| v.abs()).sum();
        let changed = relax(&mut h, r, 0.05, 83e3, 106.0, 1.43);
        assert!(changed);
        let after: f32 = h.as_slice().iter().map(|v| v.abs()).sum();
        assert!(after < 0.5 * before, "{after} vs {before}");
    }
}
