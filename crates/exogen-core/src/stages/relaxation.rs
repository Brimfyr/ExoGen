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
pub fn ice_viscosity(t: f64, melting_point: f64) -> f64 {
    let t = t.max(ICE_VISCOSITY_FLOOR_K);
    ICE_VISCOSITY_MELT
        * (ICE_CREEP_ACTIVATION / GAS_CONSTANT * (1.0 / t - 1.0 / melting_point.max(t))).exp()
}

/// The shell as relaxation sees it at a time: its thickness, the temperatures at its
/// top and base, and the gravity.
#[derive(Clone, Copy, Debug)]
pub struct Shell {
    pub thickness: f64,
    /// The cold elastic lid, which holds short loads up without flow.
    pub elastic: f64,
    pub surface_temperature: f64,
    pub basal_temperature: f64,
    pub gravity: f64,
    /// Stiffening of the ice by the rock in it, a factor on the viscosity.
    pub viscosity_scale: f64,
}

/// The relaxation timescale of a wavelength on a shell, Gyr: the ice at the depth the
/// load reaches, no deeper than the shell, along a conductive profile from the surface
/// temperature to the melting point.
pub fn relaxation_time_gyr(wavelength: f64, shell: &Shell) -> f64 {
    let d = shell.thickness.max(1.0);
    let depth = (RELAXATION_DEPTH_FRACTION * wavelength).min(d);
    let frac = (depth / d).clamp(0.0, 1.0);
    // The conductivity of ice falls as 1/T, so the conductive profile is exponential in
    // depth rather than linear: T(z) = T_s (T_m / T_s)^(z/d).
    let t_s = shell.surface_temperature.max(1.0);
    let t = t_s * (shell.basal_temperature.max(t_s) / t_s).powf(frac);
    let eta = ice_viscosity(t, shell.basal_temperature) * shell.viscosity_scale.max(1.0);
    let tau_s = 4.0 * std::f64::consts::PI * eta
        / (DENSITY_CRUST_ICE * shell.gravity * wavelength.max(1.0));
    tau_s / SECONDS_PER_GYR
}

/// The share of a load that rides on the ductile ice rather than the elastic lid: none
/// below the lid's flexural wavelength, all well beyond it.
pub fn ductile_share(wavelength: f64, shell: &Shell) -> f64 {
    let support = ELASTIC_SUPPORT_FACTOR * shell.elastic.max(0.0);
    if support <= 0.0 {
        return 1.0;
    }
    let x = wavelength / support;
    // Smooth from 0 at the flexural wavelength to 1 at twice it.
    ((x - 1.0).clamp(0.0, 1.0)).powi(2) * (3.0 - 2.0 * (x - 1.0).clamp(0.0, 1.0))
}

/// The fraction of a wavelength's relief that survives a time on a shell.
pub fn survival(wavelength: f64, dt_gyr: f64, shell: &Shell) -> f64 {
    let ductile = ductile_share(wavelength, shell);
    if ductile <= 0.0 {
        return 1.0;
    }
    let tau = relaxation_time_gyr(wavelength, shell);
    let flowed = if tau <= 0.0 {
        0.0
    } else {
        (-dt_gyr / tau).exp()
    };
    1.0 - ductile * (1.0 - flowed)
}

/// The fraction of a crater's relief that survives from its formation to the present,
/// the shell changing as the thermal history has it, sampled in steps.
pub fn survival_since(
    wavelength: f64,
    formed_gyr: f64,
    present: f64,
    shell_at: &dyn Fn(f64) -> Shell,
) -> f64 {
    if present <= formed_gyr {
        return 1.0;
    }
    let steps = 24;
    let dt = (present - formed_gyr) / steps as f64;
    let mut kept = 1.0;
    for i in 0..steps {
        let t = formed_gyr + (i as f64 + 0.5) * dt;
        kept *= survival(wavelength, dt, &shell_at(t));
        if kept < 1e-6 {
            return 0.0;
        }
    }
    kept
}

/// A pyramid of low-passed copies of a field, each level smoothing over twice the
/// scale of the last, by resampling down and up again.
fn pyramid(field: &Field<f32>, levels: usize) -> Vec<Field<f32>> {
    let grid = field.grid();
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
pub fn relax(height: &mut Field<f32>, body_radius: f64, dt_gyr: f64, shell: &Shell) -> bool {
    let grid = height.grid();
    let cell = grid.delta() * body_radius;
    let levels = ((grid.face_size() as f64).log2().floor() as usize)
        .saturating_sub(2)
        .max(1);
    // Survival per band, from the finest the pyramid resolves to the longest. Nothing to
    // do when even the longest wavelength keeps everything.
    let factors: Vec<f32> = (0..=levels)
        .map(|l| {
            let wavelength = 2.0 * cell * (2u32.pow(l as u32 + 1) as f64);
            survival(wavelength, dt_gyr, shell) as f32
        })
        .collect();
    if factors.iter().all(|&f| f > 0.999) {
        return false;
    }
    if std::env::var_os("EXOGEN_DEBUG").is_some() {
        eprintln!(
            "relax shell {:.0} m elastic {:.0} m survival by band {:?}",
            shell.thickness,
            shell.elastic,
            factors
                .iter()
                .map(|f| (f * 100.0).round() / 100.0)
                .collect::<Vec<_>>()
        );
    }
    let lows = pyramid(height, levels);
    let h = height.as_mut_slice();
    // Bands: h - low_0, low_0 - low_1, ..., low_{L-1}; the last is the global shape.
    for (k, cell) in h.iter_mut().enumerate() {
        let mut v = 0.0f32;
        let mut prev = *cell;
        for (l, low) in lows.iter().enumerate() {
            let band = prev - low.as_slice()[k];
            v += band * factors[l];
            prev = low.as_slice()[k];
        }
        v += prev * factors[levels];
        *cell = v;
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
        let ganymede = Shell {
            thickness: 83e3,
            elastic: 32e3,
            surface_temperature: 106.0,
            basal_temperature: 273.0,
            gravity: 1.43,
            viscosity_scale: 1.0,
        };
        assert!(relaxation_time_gyr(300e3, &ganymede) < 1e-3);
        assert!(relaxation_time_gyr(20e3, &ganymede) > 100.0);
        // Europa's thin shell relaxes even 30 km craters within its surface age.
        let europa = Shell {
            thickness: 8e3,
            elastic: 3e3,
            surface_temperature: 92.0,
            basal_temperature: 273.0,
            gravity: 1.31,
            viscosity_scale: 1.0,
        };
        assert!(relaxation_time_gyr(30e3, &europa) < 0.05);
        assert!(survival_since(30e3, 4.0, 4.5, &|_| europa) < 1e-3);
        assert!(survival_since(20e3, 0.0, 4.5, &|_| ganymede) > 0.95);
        // The elastic lid holds a 60 km crater on today's Ganymede however warm the base.
        assert!(survival_since(60e3, 0.0, 4.5, &|_| ganymede) > 0.95);
        assert!(survival_since(400e3, 4.0, 4.5, &|_| ganymede) < 0.05);
        // Relaxation keeps a small bump and flattens a broad one.
        let grid = CubeSphere::new(64);
        let r = 2.6e6;
        let mut h = Field::from_fn(grid, |_, d| {
            let broad = (d[0] * 1.0 - 1.0).abs() < 0.2;
            let fine = ((d[1] * 200.0).sin() * (d[2] * 200.0).sin()) as f32;
            if broad {
                1000.0 + 10.0 * fine
            } else {
                10.0 * fine
            }
        });
        let before: f32 = h.as_slice().iter().map(|v| v.abs()).sum();
        let changed = relax(&mut h, r, 0.05, &ganymede);
        assert!(changed);
        let after: f32 = h.as_slice().iter().map(|v| v.abs()).sum();
        assert!(after < 0.5 * before, "{after} vs {before}");
    }
}
