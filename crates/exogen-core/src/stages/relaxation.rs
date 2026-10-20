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
//!
//! The bands are kept from step to step (F134). The pyramid's bands overlap in wavelength,
//! so splitting the whole surface afresh at every step and scaling each band compounds:
//! a wavelength shared between a band that survives and one that relaxes loses its share
//! of the second every step, so that after many steps relief well below the shortest
//! relaxing wavelength is gone too, and the repeated filter rings outward from every edge
//! in concentric ripples. Each step instead splits only what changed since the last, adds
//! it to the bands the relief already holds, and scales those, so that every part of the
//! relief decays in the bands it had when it came, by the product of their survivals.

use crate::constants::*;

use crate::fields::Field;
use crate::grid::{Cell, CubeSphere};
use crate::stages::episodic::smooth;

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

/// The share of a wavelength's relief that a thinned shell, `factor` of the shell's
/// thickness and lid where it stands, lets survive for a time beyond what the shell
/// itself would (F66): the ground beside a renewed terrain, under the episode's thinned
/// shell while the episode lasts.
pub fn thinned_survival(wavelength: f64, dt_gyr: f64, shell: &Shell, factor: f64) -> f64 {
    let thinned = Shell {
        thickness: shell.thickness * factor,
        elastic: shell.elastic * factor,
        ..*shell
    };
    let base = survival(wavelength, dt_gyr, shell);
    if base <= 0.0 {
        return 1.0;
    }
    (survival(wavelength, dt_gyr, &thinned) / base).clamp(0.0, 1.0)
}

/// The fraction of a crater's relief that survives from its formation to the present,
/// the shell changing as the thermal history has it, sampled in steps.
pub fn survival_since(
    wavelength: f64,
    formed_gyr: f64,
    present: f64,
    shell_at: &dyn Fn(f64) -> Shell,
    local_factor: f64,
) -> f64 {
    if present <= formed_gyr {
        return 1.0;
    }
    let steps = 24;
    let dt = (present - formed_gyr) / steps as f64;
    let mut kept = 1.0;
    for i in 0..steps {
        let t = formed_gyr + (i as f64 + 0.5) * dt;
        let s = shell_at(t);
        let local = Shell {
            thickness: s.thickness * local_factor,
            elastic: s.elastic * local_factor,
            ..s
        };
        kept *= survival(wavelength, dt, &local);
        if kept < 1e-6 {
            return 0.0;
        }
    }
    kept
}

/// The field low-passed over a length, m: the Gaussian pyramid's first level at least
/// that long, expanded back to the grid.
pub fn low_pass(field: &Field<f32>, cell: f64, length: f64) -> Vec<f32> {
    let grid = field.grid();
    let levels = ((grid.face_size() as f64).log2().floor() as usize)
        .saturating_sub(2)
        .max(1);
    let level = (0..levels)
        .find(|&l| 2.0 * cell * (2u32.pow(l as u32 + 1) as f64) >= length)
        .unwrap_or(levels - 1);
    let lows = pyramid(field, level + 1);
    match lows.last() {
        Some(low) => low.as_slice().to_vec(),
        None => field.as_slice().to_vec(),
    }
}

/// Low-passed copies of a field from one Gaussian pyramid, each with its length, m, as
/// `low_pass` reads a level's, up to the first level at least `max_length` long.
pub fn low_passes(field: &Field<f32>, cell: f64, max_length: f64) -> Vec<(f64, Vec<f32>)> {
    let grid = field.grid();
    let levels = ((grid.face_size() as f64).log2().floor() as usize)
        .saturating_sub(2)
        .max(1);
    let length = |l: usize| 2.0 * cell * (2u32.pow(l as u32 + 1) as f64);
    let level = (0..levels)
        .find(|&l| length(l) >= max_length)
        .unwrap_or(levels - 1);
    pyramid(field, level + 1)
        .into_iter()
        .enumerate()
        .map(|(l, f)| (length(l), f.as_slice().to_vec()))
        .collect()
}

/// How many levels a pyramid `depth` deep holds on a grid: the halving stops once a face
/// is odd or down to four cells.
fn pyramid_depth(grid: &CubeSphere, depth: usize) -> usize {
    let mut f = grid.face_size();
    let mut n = 0;
    while n < depth && f.is_multiple_of(2) && f >= 8 {
        f /= 2;
        n += 1;
    }
    n
}

/// Passes of the smoothing kernel before each halving and after each expansion.
const PYRAMID_BLUR_PASSES: usize = 2;

/// A Gaussian pyramid of low-passed copies of a field, each level smoothing over twice
/// the scale of the last: blurred before each halving, so the coarse level holds no
/// aliased detail, and expanded back level by level with a blur after each step, so that
/// no block or tent of a coarse grid shows in the result (a box pyramid left the blocks
/// of the loop's cells in every renewed floor).
fn pyramid(field: &Field<f32>, levels: usize) -> Vec<Field<f32>> {
    let (grids, coarse_levels) = pyramid_coarse(field, levels);
    coarse_levels
        .iter()
        .enumerate()
        .map(|(l, c)| pyramid_expand(&grids, c, l))
        .collect()
}

/// The coarse levels of a Gaussian pyramid, each blurred and halved from the last, with
/// the grids from the field's own down: for expanding one level at a time where every
/// level at full resolution at once would not fit.
pub(crate) fn pyramid_coarse(
    field: &Field<f32>,
    levels: usize,
) -> (Vec<CubeSphere>, Vec<Field<f32>>) {
    let mut grids = vec![field.grid()];
    let mut coarse_levels: Vec<Field<f32>> = Vec::with_capacity(levels);
    let mut current = field.clone();
    for _ in 0..levels {
        let Some(coarse) = halve(&blurred(&current, PYRAMID_BLUR_PASSES)) else {
            break;
        };
        grids.push(coarse.grid());
        coarse_levels.push(coarse.clone());
        current = coarse;
    }
    (grids, coarse_levels)
}

/// Level `l` of a pyramid expanded back to the field's grid, blurred after each step.
pub(crate) fn pyramid_expand(grids: &[CubeSphere], level: &Field<f32>, l: usize) -> Field<f32> {
    let mut f = level.clone();
    for g in grids[..=l].iter().rev() {
        f = blurred(&f.resample(*g), PYRAMID_BLUR_PASSES);
    }
    f
}

/// The Gaussian spread of a pyramid's level `l`, in cells of the field's grid: each pass
/// of the kernel adds a quarter of a cell's variance at its level, each halving a two-cell
/// box's quarter, and each expansion the bilinear tent's, a sixth of the coarser spacing
/// squared.
pub(crate) fn pyramid_spread(l: usize) -> f64 {
    let per_level = 0.5 * PYRAMID_BLUR_PASSES as f64 + 0.25 + 4.0 / 6.0;
    (per_level * (4f64.powi(l as i32 + 1) - 1.0) / 3.0).sqrt()
}

/// The field smoothed by passes of the 1/8, 1/2, 1/8 kernel on its own grid.
pub(crate) fn blurred(field: &Field<f32>, passes: usize) -> Field<f32> {
    let grid = field.grid();
    let v = smooth(&grid, field.as_slice(), passes);
    let mut out = Field::new(grid);
    out.as_mut_slice().copy_from_slice(&v);
    out
}

/// The field averaged onto the grid with half the cells per face edge, each coarse cell
/// the mean of its four fine cells: an area average, so that no single cell stands for a
/// patch. None once a face is odd or down to four cells.
pub(crate) fn halve(field: &Field<f32>) -> Option<Field<f32>> {
    let grid = field.grid();
    let f = grid.face_size();
    if !f.is_multiple_of(2) || f < 8 {
        return None;
    }
    let coarse = CubeSphere::new(f / 2);
    let src = field.as_slice();
    let mut out = Field::new(coarse);
    for (k, v) in out.as_mut_slice().iter_mut().enumerate() {
        let c = coarse.cell(k);
        let mut sum = 0.0f32;
        for (di, dj) in [(0, 0), (1, 0), (0, 1), (1, 1)] {
            sum += src[grid.index(Cell {
                face: c.face,
                i: 2 * c.i + di,
                j: 2 * c.j + dj,
            })];
        }
        *v = 0.25 * sum;
    }
    Some(out)
}

/// The relief a relaxing surface holds, band by band (F134): the bands of everything the
/// surface has taken in, each scaled by the survivals of the steps since, and the surface
/// as the last step left it, against which the next step finds what changed.
#[derive(Default, Clone)]
pub struct Bands {
    bands: Vec<Vec<f32>>,
    last: Option<Vec<f32>>,
}

impl Bands {
    /// Take in what changed in `height` since the last step, split into the pyramid's
    /// `levels` bands and the long remainder, and added to the bands held.
    fn absorb(&mut self, height: &Field<f32>, levels: usize) {
        let grid = height.grid();
        let n = grid.cell_count();
        let mut delta = Field::<f32>::new(grid);
        match &self.last {
            Some(last) if last.len() == n => {
                for ((d, h), l) in delta
                    .as_mut_slice()
                    .iter_mut()
                    .zip(height.as_slice())
                    .zip(last)
                {
                    *d = h - l;
                }
            }
            _ => delta.as_mut_slice().copy_from_slice(height.as_slice()),
        }
        let lows = pyramid(&delta, levels);
        let count = lows.len() + 1;
        if self.bands.len() != count || self.bands.iter().any(|b| b.len() != n) {
            self.bands = vec![vec![0.0; n]; count];
        }
        for (k, &d) in delta.as_slice().iter().enumerate() {
            let mut prev = d;
            for (l, low) in lows.iter().enumerate() {
                let v = low.as_slice()[k];
                self.bands[l][k] += prev - v;
                prev = v;
            }
            self.bands[count - 1][k] += prev;
        }
    }

    /// Scale the bands held by each band's factor, times a cell's own for the bands that
    /// have one (`own(l, k)`), and write their sum back to `height`, which the next step
    /// measures its change from.
    fn decay(
        &mut self,
        height: &mut Field<f32>,
        factors: &[f32],
        own: impl Fn(usize, usize) -> f32,
    ) {
        let count = self.bands.len();
        let h = height.as_mut_slice();
        for (k, cell) in h.iter_mut().enumerate() {
            let mut v = 0.0f32;
            for (l, band) in self.bands.iter_mut().enumerate() {
                let f = factors[l.min(factors.len() - 1)] * if l + 1 < count { own(l, k) } else { 1.0 };
                band[k] *= f;
                v += band[k];
            }
            *cell = v;
        }
        self.last = Some(h.to_vec());
    }
}

/// Relax a height field over a step: each wavelength band, from twice the cell up to
/// the whole body, keeps the share of its relief the shell's viscosity lets it keep.
/// What varies across the shell for one relaxation step: the equilibrium surface the
/// topography relaxes towards, since an isostatic depression is not a load; and the
/// cells of a renewed terrain, which give up the wavelengths shorter than its renewal
/// length, m, by the weight each carries (1 inside, less across its edge and its halo),
/// and relax the rest with their neighbours.
///
/// The survival factors themselves are those of the mean shell everywhere. Giving a soft
/// zone its own factors beside a stiff one does not work in a band decomposition: a soft
/// cell told to relax a wavelength its neighbours keep subtracts a regional mean it cannot
/// lower, step after step, and sinks many times its offset (the Phase 4 report). The local
/// thickness acts through the replayed craters' survival and the crack density instead.
#[derive(Default, Clone, Copy)]
pub struct Local<'a> {
    pub equilibrium: Option<&'a [f32]>,
    pub renew: Option<(&'a [f32], f64)>,
    /// While a terrain renews, the factor its episode's thinning puts on the shell under
    /// each cell, and for how long, Gyr: the ground beside the terrain relaxes the
    /// wavelengths shorter than the renewal length under that thinner shell (F66), each
    /// by its own survival.
    pub halo: Option<(&'a [f32], f64)>,
}

/// The step, on the bands the surface holds (F134). Returns the relaxed surface, which is
/// the renewed terrain's floor where it is renewed, when there is such a terrain.
pub fn relax(
    bands: &mut Bands,
    height: &mut Field<f32>,
    body_radius: f64,
    dt_gyr: f64,
    shell: &Shell,
    local: &Local,
) -> Option<Vec<f32>> {
    let grid = height.grid();
    let cell = grid.delta() * body_radius;
    if let Some(eq) = local.equilibrium {
        for (h, e) in height.as_mut_slice().iter_mut().zip(eq) {
            *h -= e;
        }
    }
    let result = relax_deviation(bands, height, cell, dt_gyr, shell, local);
    if let Some(eq) = local.equilibrium {
        for (h, e) in height.as_mut_slice().iter_mut().zip(eq) {
            *h += e;
        }
    }
    result.map(|mut floor| {
        if let Some(eq) = local.equilibrium {
            for (f, e) in floor.iter_mut().zip(eq) {
                *f += e;
            }
        }
        floor
    })
}

fn relax_deviation(
    held: &mut Bands,
    height: &mut Field<f32>,
    cell: f64,
    dt_gyr: f64,
    shell: &Shell,
    local: &Local,
) -> Option<Vec<f32>> {
    let grid = height.grid();
    let depth = ((grid.face_size() as f64).log2().floor() as usize)
        .saturating_sub(2)
        .max(1);
    let levels = pyramid_depth(&grid, depth).max(1);
    // Survival per band, from the finest the pyramid resolves to the longest. Nothing to
    // do when even the longest wavelength keeps everything and nothing is renewed: what
    // changed waits in the surface for the next step that relaxes.
    let wavelengths: Vec<f64> = (0..=levels)
        .map(|l| 2.0 * cell * (2u32.pow(l as u32 + 1) as f64))
        .collect();
    let factors: Vec<f32> = wavelengths
        .iter()
        .map(|&w| survival(w, dt_gyr, shell) as f32)
        .collect();
    if local.renew.is_none() && factors.iter().all(|&f| f > 0.999) {
        return None;
    }
    // A renewed terrain keeps nothing shorter than the first band at least its renewal
    // length long, and never the finest band: new ice has no cell-scale relief, whatever
    // the grid's cell.
    let renew_level = local.renew.map(|(_, length)| {
        wavelengths
            .iter()
            .position(|&w| w >= length)
            .unwrap_or(levels)
            .max(1)
            .min(levels - 1)
    });
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
    // What a cell keeps of the bands shorter than the renewal length: none inside a
    // renewed terrain, by its weight across the edge, and beside it what the episode's
    // thinner shell lets each band keep (F66), from a table over the thinning factor.
    let kept = |k: usize| -> (usize, f32) {
        match (local.renew, renew_level) {
            (Some((weights, _)), Some(level)) => (level, 1.0 - weights[k].clamp(0.0, 1.0)),
            _ => (0, 1.0),
        }
    };
    const STEPS: usize = 100;
    let thinned: Vec<Vec<f32>> = match (local.halo, renew_level) {
        (Some((_, dt)), Some(level)) if dt > 0.0 => (0..level)
            .map(|l| {
                (0..=STEPS)
                    .map(|i| {
                        thinned_survival(wavelengths[l], dt, shell, i as f64 / STEPS as f64) as f32
                    })
                    .collect()
            })
            .collect(),
        _ => Vec::new(),
    };
    let beside = |l: usize, k: usize| -> f32 {
        match (local.halo, thinned.get(l)) {
            (Some((factors, _)), Some(table)) => {
                let x = factors[k].clamp(0.0, 1.0) * STEPS as f32;
                let i = (x.floor() as usize).min(STEPS - 1);
                let f = x - i as f32;
                table[i] + (table[i + 1] - table[i]) * f
            }
            _ => 1.0,
        }
    };
    // Bands: h - low_0, low_0 - low_1, ..., low_{L-1}; the last is the global shape, held
    // from step to step. A renewed cell gives up the bands shorter than its renewal length
    // by its weight and relaxes the rest like its neighbours, so that it never stands apart
    // from them.
    held.absorb(height, depth);
    held.decay(height, &factors, |l, k| {
        let (level, keep) = kept(k);
        if l < level { keep * beside(l, k) } else { 1.0 }
    });
    // The whole relaxed surface, for the renewed terrain's floor: defined everywhere, so
    // that a floor laid at a finer resolution, where the terrain reaches past the loop's
    // cells, has something continuous to lay.
    renew_level.map(|_| height.as_slice().to_vec())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_pyramid_level_spreads_as_its_kernels_add_up() {
        // A point spread through each level of the pyramid has the spread pyramid_spread
        // gives it, which the gardening's mixing reads to pick its levels (F95).
        let grid = CubeSphere::new(128);
        for l in 0..4 {
            let mut worst: f64 = 0.0;
            for (i0, j0) in [(64, 64), (65, 64), (65, 65)] {
                let mut field = Field::<f32>::new(grid);
                field.as_mut_slice()[grid.index(Cell {
                    face: 0,
                    i: i0,
                    j: j0,
                })] = 1.0;
                let (grids, coarse) = pyramid_coarse(&field, l + 1);
                let f = pyramid_expand(&grids, &coarse[l], l);
                let (mut w, mut m1, mut m2) = (0.0, 0.0, 0.0);
                for i in 0..128 {
                    for j in 0..128 {
                        let v = f.as_slice()[grid.index(Cell { face: 0, i, j })] as f64;
                        w += v;
                        m1 += v * i as f64;
                        m2 += v * (i * i) as f64;
                    }
                }
                let spread = (m2 / w - (m1 / w).powi(2)).sqrt();
                worst = worst.max((spread / pyramid_spread(l) - 1.0).abs());
            }
            assert!(worst < 0.1, "level {l}: off by {worst:.3}");
        }
    }

    #[test]
    fn a_thinned_shell_relaxes_the_long_wavelengths_beside_a_terrain_first() {
        // F66: under an Enceladus shell thinned to 0.45 of its 44 km for an episode of 10
        // Myr, a 100 km load relaxes and a 12 km crater keeps its relief; the shell as it
        // is relaxes neither.
        let shell = Shell {
            thickness: 44e3,
            elastic: 23e3,
            surface_temperature: 75.0,
            basal_temperature: 273.0,
            gravity: 0.113,
            viscosity_scale: 1.0,
        };
        assert!(thinned_survival(100e3, 0.01, &shell, 0.45) < 0.05);
        assert!(thinned_survival(12e3, 0.01, &shell, 0.45) > 0.99);
        assert!((thinned_survival(100e3, 0.01, &shell, 1.0) - 1.0).abs() < 1e-9);
        // A thinner shell, a larger body's, relaxes a crater of 25 km beside the terrain.
        let thin = Shell {
            thickness: 13.3e3,
            elastic: 6.9e3,
            ..shell
        };
        assert!(thinned_survival(25e3, 0.01, &thin, 0.3) < 0.05);
    }

    /// F134: a crater whose long bands relax and short bands hold, stepped many times with
    /// the same survivals, is relaxed once and then left alone: the bands it keeps do not
    /// wear away step by step, and no ripples spread out from its rim. Split afresh at every
    /// step, as the bands were, the same relief wears down and rings.
    #[test]
    fn a_relaxed_crater_neither_rings_nor_wears_away() {
        let grid = CubeSphere::new(128);
        let center = [1.0f64, 0.0, 0.0];
        let radius = 0.12;
        let crater = Field::from_fn(grid, |_, d| {
            let r = crate::geom::angle_between(center, d) / radius;
            let bowl = if r < 1.0 { -1000.0 * (1.0 - r * r) } else { 0.0 };
            let rim = 300.0 * (-((r - 1.0) / 0.12).powi(2)).exp();
            (bowl + rim) as f32
        });
        let depth = ((grid.face_size() as f64).log2().floor() as usize).saturating_sub(2);
        let count = pyramid_depth(&grid, depth) + 1;
        // The three finest bands hold, the rest go at once.
        let factors: Vec<f32> = (0..count).map(|l| if l < 3 { 1.0 } else { 0.0 }).collect();
        let run = |steps: usize, fresh: bool| -> Vec<f32> {
            let mut h = crater.clone();
            let mut held = Bands::default();
            for _ in 0..steps {
                if fresh {
                    held = Bands::default();
                }
                held.absorb(&h, depth);
                held.decay(&mut h, &factors, |_, _| 1.0);
            }
            h.as_slice().to_vec()
        };
        let once = run(1, false);
        let many = run(60, false);
        let worst = once
            .iter()
            .zip(&many)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0f32, f32::max);
        assert!(worst < 1e-2, "kept bands changed by {worst} m over 60 steps");
        // Far outside the crater, beyond two radii, nothing.
        let far = |h: &[f32]| -> f32 {
            (0..grid.cell_count())
                .filter(|&k| {
                    crate::geom::angle_between(center, grid.center(grid.cell(k))) > 2.5 * radius
                })
                .map(|k| h[k].abs())
                .fold(0.0f32, f32::max)
        };
        assert!(far(&many) < 30.0, "ripples of {} m beyond the crater", far(&many));
        // The old way, split afresh at every step: the kept relief wears down.
        let fresh = run(60, true);
        let lost = once
            .iter()
            .zip(&fresh)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0f32, f32::max);
        eprintln!(
            "kept bands changed {worst} m; beyond the crater {} m held, {} m split afresh; lost {lost} m",
            far(&many),
            far(&fresh)
        );
        assert!(lost > 10.0 * worst.max(1.0), "split afresh: {lost} m");
    }

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
        assert!(survival_since(30e3, 4.0, 4.5, &|_| europa, 1.0) < 1e-3);
        assert!(survival_since(20e3, 0.0, 4.5, &|_| ganymede, 1.0) > 0.95);
        // The elastic lid holds a 60 km crater on today's Ganymede however warm the base.
        assert!(survival_since(60e3, 0.0, 4.5, &|_| ganymede, 1.0) > 0.95);
        assert!(survival_since(400e3, 4.0, 4.5, &|_| ganymede, 1.0) < 0.05);
        // Where the shell is a fifth as thick the same crater goes.
        assert!(survival_since(60e3, 0.0, 4.5, &|_| ganymede, 0.2) < 0.5);
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
        relax(&mut Bands::default(), &mut h, r, 0.05, &ganymede, &Local::default());
        let after: f32 = h.as_slice().iter().map(|v| v.abs()).sum();
        assert!(after < 0.5 * before, "{after} vs {before}");
    }
}
