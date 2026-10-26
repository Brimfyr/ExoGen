// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! Viscous relaxation of an ice shell's topography (PLAN.md 4.3 and Phase 4), and of a
//! rocky lid's by the flow of its lower crust (F170).
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

/// Bumped when the behaviour changes, so caches are not reused across versions.
pub const VERSION: u64 = 3;

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
    /// A rocky lid's crust (F170): `thickness` is the crust's, `basal_temperature` its
    /// base's, and the relief relaxes by the flow of its lower crust.
    pub rock: bool,
    /// The stress at which the lid breaks in bending, Pa, the crust's strength (F15): it
    /// holds up no more of a load than bends it that far (F195).
    pub strength: f64,
}

/// The viscosity of a rocky lid's lower crust at a temperature, Pa s, under the stress a
/// kilometre of relief drives on a body of this gravity (F170).
pub fn crust_viscosity(t: f64, gravity: f64) -> f64 {
    let (a, n, q) = DIABASE_CREEP;
    let stress = DENSITY_CRUST_ROCK * gravity * CRUST_FLOW_RELIEF;
    let rate = a * (stress / 1e6).powf(n) * (-q / (GAS_CONSTANT * t.max(1.0))).exp();
    if rate <= 0.0 {
        return f64::INFINITY;
    }
    (stress / (2.0 * rate)).max(CRUST_FLOW_MIN_VISCOSITY)
}

/// The relaxation timescale of a wavelength on a rocky lid, Gyr (F170). The lower crust
/// flows where it is warm, a channel at the crust's base as thick as the depth over which
/// its viscosity falls by e, δ = R T² / (Q dT/dz) along the crust's conductive profile,
/// beneath crust it holds rigid; the pressure of the relief above drives it, moving the
/// crust from under the high ground to under the low, so that, the crust kept compensated,
/// the relief diffuses with κ = (Δρ/ρₘ) ρ g δ³ / η and a wavelength decays as (λ/2π)² / κ:
/// the longest wavelengths last (McKenzie, Nimmo and others 2000; Nimmo and Stevenson 2001).
pub fn rock_relaxation_time_gyr(wavelength: f64, shell: &Shell) -> f64 {
    let (_, _, q) = DIABASE_CREEP;
    let t_base = shell.basal_temperature.max(shell.surface_temperature + 1.0);
    let gradient = (t_base - shell.surface_temperature) / shell.thickness.max(1.0);
    let channel = (GAS_CONSTANT * t_base * t_base / (q * gradient)).min(shell.thickness);
    let eta = crust_viscosity(t_base, shell.gravity);
    let airy = (DENSITY_ROCK - DENSITY_CRUST_ROCK) / DENSITY_ROCK;
    let kappa = airy * DENSITY_CRUST_ROCK * shell.gravity * channel.powi(3) / eta;
    if kappa <= 0.0 {
        return f64::INFINITY;
    }
    let k = wavelength / std::f64::consts::TAU;
    k * k / kappa / SECONDS_PER_GYR
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

/// The flexural wavelength of the elastic lid, m (F139): 2π (D / ρ g)^¼, its flexural
/// rigidity D = E Tₑ³ / 12 (1 − ν²) (Turcotte and Schubert 2002, section 3.17). It goes
/// as Tₑ^¾ g^-¼, so that a small moon's lid holds up loads many times its thickness:
/// Mimas's 12 km lid 450 km, Enceladus's 12 km 370 km, Ganymede's 32 km 420 km.
pub fn flexural_wavelength(shell: &Shell) -> f64 {
    let te = shell.elastic.max(0.0);
    let (modulus, poisson, density) = if shell.rock {
        (ROCK_YOUNG_MODULUS, ROCK_POISSON, DENSITY_CRUST_ROCK)
    } else {
        (ICE_YOUNG_MODULUS, ICE_POISSON, DENSITY_CRUST_ICE)
    };
    let rigidity = modulus * te.powi(3) / (12.0 * (1.0 - poisson * poisson));
    std::f64::consts::TAU * (rigidity / (density * shell.gravity.max(1e-9))).powf(0.25)
}

/// The share of a load that rides on the ductile ice rather than the elastic lid (F139).
/// A thin plate's degree of compensation is 1 / (1 + (λ_f / λ)⁴), half at the flexural
/// wavelength; the share is taken as rising smoothly from none at two thirds of it to all
/// at four thirds, so that the share the lid holds of a short load, small but never nil in
/// the plate's law, is not worn away by the survivals of step after step.
pub fn ductile_share(wavelength: f64, shell: &Shell) -> f64 {
    let support = 2.0 / 3.0 * flexural_wavelength(shell);
    if support <= 0.0 {
        return 1.0;
    }
    let x = (wavelength / support - 1.0).clamp(0.0, 1.0);
    x * x * (3.0 - 2.0 * x)
}

/// The share of a wavelength's relief riding on the ductile ice or lower crust that is
/// still there after a time: what the flow has not yet carried away.
pub fn unflowed(wavelength: f64, dt_gyr: f64, shell: &Shell) -> f64 {
    let tau = if shell.rock {
        rock_relaxation_time_gyr(wavelength, shell)
    } else {
        relaxation_time_gyr(wavelength, shell)
    };
    if tau <= 0.0 {
        0.0
    } else {
        (-dt_gyr / tau).exp()
    }
}

/// The fraction of a wavelength's relief that survives a time on a shell.
pub fn survival(wavelength: f64, dt_gyr: f64, shell: &Shell) -> f64 {
    let ductile = ductile_share(wavelength, shell);
    if ductile <= 0.0 {
        return 1.0;
    }
    1.0 - ductile * (1.0 - unflowed(wavelength, dt_gyr, shell))
}

/// The most relief of a wavelength an ice lid holds up, m (F195): the weight of a
/// sinusoidal load h bends the lid by a moment ρ g h / k², which brings its surface to a
/// stress 6 M / Tₑ²; past the lid's strength it breaks, and what it cannot hold rides on
/// the ice below and flows away as that ice does. A rocky lid's relief stands on its
/// crust's compensation (F170), not on the lid's bending, so none is bounded there.
pub fn lid_relief(wavelength: f64, shell: &Shell) -> f64 {
    if shell.rock {
        return f64::INFINITY;
    }
    let k = std::f64::consts::TAU / wavelength.max(1.0);
    shell.strength * (shell.elastic * k).powi(2)
        / (6.0 * DENSITY_CRUST_ICE * shell.gravity.max(1e-9))
}

/// A band's relief at a cell after a time on a shell (F195): the lid holds its share up
/// (`1 - ductile`), no more than its strength bears, approached softly as a layer's relief
/// is (F189); the rest rides on the ductile ice and keeps what the flow leaves
/// (`unflowed`). With a strength past any load this is `survival` times the relief.
pub fn kept_relief(relief: f32, held_share: f64, unflowed: f64, most: f64) -> f32 {
    let held = held_share * relief as f64;
    let carried = if most.is_finite() && most > 0.0 {
        most * (held / most).tanh()
    } else if most.is_finite() {
        0.0
    } else {
        held
    };
    (carried + (relief as f64 - carried) * unflowed) as f32
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

/// A cell's value low-passed over `length`, from `low_passes`' levels: between the two
/// levels that bracket the length, by its logarithm, so that two lengths near each other
/// give near the same floor and no seam runs between terrains of different widths (F141).
pub fn blend_levels(levels: &[(f64, Vec<f32>)], length: f64, k: usize) -> f32 {
    let i = levels.iter().position(|(l, _)| *l >= length);
    match i {
        None => levels.last().map_or(0.0, |(_, v)| v[k]),
        Some(0) => levels[0].1[k],
        Some(i) => {
            let (l0, v0) = (&levels[i - 1].0, &levels[i - 1].1);
            let (l1, v1) = (&levels[i].0, &levels[i].1);
            let t = ((length / l0).ln() / (l1 / l0).ln()).clamp(0.0, 1.0) as f32;
            v0[k] + (v1[k] - v0[k]) * t
        }
    }
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
    /// When kept, the potential of the flow that relaxed the relief, m³ (F188): each band's
    /// relaxed relief over its wavenumber squared, summed over the steps since it was
    /// cleared, whose curvature is the flow's strain times the crust's root share and
    /// thickness.
    pub flow: Option<Vec<f32>>,
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

    /// Relax the bands held, each cell's relief in a band to what `keep(l, shape, k,
    /// relief)` leaves of it, `shape` for the last band, the global shape, and write their
    /// sum back to `height`, which the next step measures its change from.
    fn decay(
        &mut self,
        height: &mut Field<f32>,
        lengths: &[f64],
        keep: impl Fn(usize, bool, usize, f32) -> f32,
    ) {
        let count = self.bands.len();
        let h = height.as_mut_slice();
        let mut flow = self.flow.take().filter(|f| f.len() == h.len());
        let reach: Vec<f32> = (0..count)
            .map(|l| {
                let w = lengths.get(l).or(lengths.last()).copied().unwrap_or(0.0);
                (w / std::f64::consts::TAU).powi(2) as f32
            })
            .collect();
        for (k, cell) in h.iter_mut().enumerate() {
            let mut v = 0.0f32;
            let mut potential = 0.0f32;
            for (l, band) in self.bands.iter_mut().enumerate() {
                let kept = keep(l, l + 1 == count, k, band[k]);
                potential += (kept - band[k]) * reach[l];
                band[k] = kept;
                v += band[k];
            }
            *cell = v;
            if let Some(flow) = flow.as_mut() {
                flow[k] += potential;
            }
        }
        self.flow = flow;
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
    /// A rocky lid's crust under each cell (F184): each cell relaxes under its own crust,
    /// the lid as it stood through the step's substeps.
    pub rock: Option<RockCrust<'a>>,
    pub renew: Option<(&'a [f32], f64)>,
    /// While a terrain renews, the factor its episode's thinning puts on the shell under
    /// each cell, and for how long, Gyr: the ground beside the terrain relaxes the
    /// wavelengths shorter than the renewal length under that thinner shell (F66), each
    /// by its own survival.
    pub halo: Option<(&'a [f32], f64)>,
    /// The base relief's bands shorter than the renewal length (F141), which a renewed
    /// cell takes in the measure it gave up what stood there: the texture under the new
    /// terrain, where the low-pass alone left a blur of the old craters.
    pub fine: Option<&'a [f32]>,
}

/// A rocky lid's crust for relaxation (F184): its thickness under each cell, m, the shell
/// of the lid at a time over a crust of a thickness, and the step's substeps, (time,
/// duration) Gyr, over which the lid changes.
#[derive(Clone, Copy)]
pub struct RockCrust<'a> {
    pub crust: &'a [f32],
    pub shell_at: &'a (dyn Fn(f64, f64) -> Shell + Sync),
    pub steps: &'a [(f64, f64)],
}

/// The share of a wavelength's relief a rocky lid lets survive over a step's substeps under
/// a crust of `thickness`, m (F184): the survivals of the lid as it stood at each.
pub fn rock_survival(wavelength: f64, rock: &RockCrust, thickness: f64) -> f64 {
    rock.steps
        .iter()
        .map(|&(t, dt)| survival(wavelength, dt, &(rock.shell_at)(t, thickness)))
        .product()
}

/// Whether a shell leaves a wavelength's relief as it is over a time, whatever the relief:
/// it keeps all of it, and the ice under the lid would not flow from under what the lid
/// could not hold (F195).
fn band_still(wavelength: f64, dt_gyr: f64, shell: &Shell) -> bool {
    survival(wavelength, dt_gyr, shell) as f32 > 0.999
        && unflowed(wavelength, dt_gyr, shell) as f32 > 0.999
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
    // Nothing to do when no band relaxes and nothing is renewed: leave the surface as it is,
    // not even shifted to its equilibrium and back.
    if local.renew.is_none() {
        let depth = ((grid.face_size() as f64).log2().floor() as usize)
            .saturating_sub(2)
            .max(1);
        let levels = pyramid_depth(&grid, depth).max(1);
        // Under each cell's own crust, the thickest relaxes most.
        let thickest = local
            .rock
            .as_ref()
            .map(|r| r.crust.iter().fold(0.0f32, |m, &c| m.max(c)) as f64);
        let still = (0..=levels).all(|l| {
            let w = 2.0 * cell * (2u32.pow(l as u32 + 1) as f64);
            match (&local.rock, thickest) {
                (Some(r), Some(c)) => rock_survival(w, r, c) as f32 > 0.999,
                // An ice lid holds a load only up to its strength (F195): still only where
                // the ice below would not flow from under what the lid cannot hold.
                _ => band_still(w, dt_gyr, shell),
            }
        });
        if still {
            return None;
        }
    }
    if let Some(eq) = local.equilibrium {
        for (h, e) in height.as_mut_slice().iter_mut().zip(eq) {
            *h -= e;
        }
    }
    // The ground's mean is its height above the body's mean radius, no load: what the
    // relaxation takes from the mean it gives back, so that a datum's move (F177, F184) is
    // not relaxed away.
    let mean_of = |h: &[f32]| -> f64 {
        let (mut sum, mut area) = (0.0f64, 0.0f64);
        for (k, &v) in h.iter().enumerate() {
            let a = grid.cell_area(grid.cell(k));
            sum += v as f64 * a;
            area += a;
        }
        sum / area.max(1e-30)
    };
    let mean = mean_of(height.as_slice());
    let result = relax_deviation(bands, height, cell, dt_gyr, shell, local);
    let shifted = (mean - mean_of(height.as_slice())) as f32;
    for h in height.as_mut_slice().iter_mut() {
        *h += shifted;
    }
    if let Some(eq) = local.equilibrium {
        for (h, e) in height.as_mut_slice().iter_mut().zip(eq) {
            *h += e;
        }
    }
    result.map(|mut floor| {
        for f in floor.iter_mut() {
            *f += shifted;
        }
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
    // Under a rocky lid's crust each cell keeps each band by its own crust's survival
    // (F184), from a table over the crusts the lid holds, the global shape under the
    // median crust.
    let tables: Option<(Vec<Vec<f32>>, f32, f32)> = local.rock.as_ref().map(|r| {
        let (lo, hi) = r
            .crust
            .iter()
            .fold((f32::MAX, f32::MIN), |(a, b), &c| (a.min(c), b.max(c)));
        let tables = wavelengths
            .iter()
            .map(|&w| {
                (0..=ROCK_CRUST_TABLE)
                    .map(|j| {
                        let c = lo + (hi - lo) * j as f32 / ROCK_CRUST_TABLE as f32;
                        rock_survival(w, r, c as f64) as f32
                    })
                    .collect()
            })
            .collect();
        (tables, lo, hi)
    });
    let under = |l: usize, thickness: f32| -> f32 {
        match &tables {
            Some((tables, lo, hi)) => {
                let x = ((thickness - lo) / (hi - lo).max(1e-3)).clamp(0.0, 1.0)
                    * ROCK_CRUST_TABLE as f32;
                let i = (x.floor() as usize).min(ROCK_CRUST_TABLE - 1);
                let f = x - i as f32;
                tables[l][i] + (tables[l][i + 1] - tables[l][i]) * f
            }
            None => 1.0,
        }
    };
    let factors: Vec<f32> = match (&local.rock, &tables) {
        (Some(r), Some(_)) => {
            let mut sorted = r.crust.to_vec();
            sorted.sort_by(f32::total_cmp);
            let median = sorted[sorted.len() / 2];
            (0..wavelengths.len())
                .map(|l| if l + 1 < wavelengths.len() { 1.0 } else { under(l, median) })
                .collect()
        }
        _ => wavelengths
            .iter()
            .map(|&w| survival(w, dt_gyr, shell) as f32)
            .collect(),
    };
    // An ice lid holds each band's share up to its strength and lets the rest ride on the
    // ice below (F195): per band, the share the lid holds, what the flow leaves of the rest
    // and the most relief the lid holds.
    let lid: Option<Vec<(f64, f64, f64)>> = local.rock.is_none().then(|| {
        wavelengths
            .iter()
            .map(|&w| {
                (
                    1.0 - ductile_share(w, shell),
                    unflowed(w, dt_gyr, shell),
                    lid_relief(w, shell),
                )
            })
            .collect()
    });
    let still = match &tables {
        Some((tables, _, _)) => tables.iter().all(|t| t.iter().all(|&f| f > 0.999)),
        None => wavelengths.iter().all(|&w| band_still(w, dt_gyr, shell)),
    };
    if local.renew.is_none() && still {
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
            "relax shell {:.0} m elastic {:.0} m survival by band {:?}, the lid holding at most {:?} m",
            shell.thickness,
            shell.elastic,
            factors
                .iter()
                .map(|f| (f * 100.0).round() / 100.0)
                .collect::<Vec<_>>(),
            lid.as_ref().map(|lid| {
                lid.iter()
                    .map(|&(_, _, most)| most.min(1e6).round())
                    .collect::<Vec<_>>()
            })
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
    // Under an ice lid each band keeps what the lid holds of it and what the flow leaves of
    // the rest (F195), beside a renewing terrain under the episode's thinner lid, which
    // holds as the square of its thickness.
    held.absorb(height, depth);
    held.decay(height, &wavelengths, |l, shape, k, relief| {
        let b = l.min(factors.len() - 1);
        let (level, keep) = kept(k);
        let thin = l < level && !shape;
        let relaxed = match &lid {
            Some(lid) => {
                let (share, unflowed, most) = lid[b];
                let lid_factor = match local.halo {
                    Some((thinning, _)) if thin => thinning[k].clamp(0.0, 1.0) as f64,
                    _ => 1.0,
                };
                kept_relief(relief, share, unflowed, most * lid_factor * lid_factor)
            }
            None => relief * factors[b],
        };
        if shape {
            return relaxed;
        }
        let own = if thin { keep * beside(l, k) } else { 1.0 };
        match &local.rock {
            Some(r) => relaxed * own * under(l, r.crust[k]),
            None => relaxed * own,
        }
    });
    // What a renewed cell gave up of the short bands it takes back from the base relief
    // (F141): laid now, it relaxes from now on like anything new.
    if let (Some(fine), Some((weights, _))) = (local.fine, local.renew) {
        for (k, h) in height.as_mut_slice().iter_mut().enumerate() {
            let w = weights[k].clamp(0.0, 1.0);
            if w > 0.0 {
                *h += fine[k] * w;
            }
        }
    }
    // The whole relaxed surface, for the renewed terrain's floor: defined everywhere, so
    // that a floor laid at a finer resolution, where the terrain reaches past the loop's
    // cells, has something continuous to lay.
    renew_level.map(|_| height.as_slice().to_vec())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// F188: the flow's potential the bands keep has minus the relaxed relief for its
    /// Laplacian, to within the bands' spread of wavelengths, which the strain's reading
    /// scales out (crust::StrainRecord::add_collapse).
    #[test]
    fn the_relaxing_flow_s_potential_has_the_relaxed_relief_for_its_laplacian() {
        let grid = CubeSphere::new(64);
        let center = grid.center(crate::grid::Cell { face: 0, i: 32, j: 32 });
        let sigma = 6.0 * grid.delta();
        let mut h = Field::<f32>::new(grid);
        for (k, v) in h.as_mut_slice().iter_mut().enumerate() {
            let x = crate::geom::angle_between(grid.center(grid.cell(k)), center) / sigma;
            *v = (1000.0 * (-0.5 * x * x).exp()) as f32;
        }
        let before = h.as_slice().to_vec();
        let depth = ((grid.face_size() as f64).log2().floor() as usize).saturating_sub(2);
        let levels = pyramid_depth(&grid, depth).max(1);
        let lengths: Vec<f64> = (0..=levels)
            .map(|l| 2.0 * grid.delta() * (2u32.pow(l as u32 + 1) as f64))
            .collect();
        let mut held = Bands {
            flow: Some(vec![0.0; grid.cell_count()]),
            ..Default::default()
        };
        held.absorb(&h, depth);
        held.decay(&mut h, &lengths, |_, _, _, v| 0.5 * v);
        let flow = held.flow.unwrap();
        for (i, j) in [(32, 32), (38, 32), (32, 41)] {
            let c = grid.cell(grid.index(crate::grid::Cell { face: 0, i, j }));
            let k = grid.index(c);
            let side = |di: i32, dj: i32| flow[grid.index(grid.neighbor(c, di, dj))] as f64;
            let laplacian = (side(1, 0) + side(-1, 0) + side(0, 1) + side(0, -1)
                - 4.0 * flow[k] as f64)
                / (grid.delta() * grid.delta());
            let relaxed = (h.as_slice()[k] - before[k]) as f64;
            let ratio = -laplacian / relaxed;
            assert!(
                ratio > 0.3 && ratio < 3.0,
                "({i}, {j}): laplacian {laplacian}, relaxed {relaxed}"
            );
        }
    }

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

    /// F170: a rocky lid's lower crust flows where it is warm. Under a crust 60 km thick on
    /// Mars whose base stands at 1,550 K, a wavelength of 1,000 km loses most of its relief
    /// in a hundred million years, one of 10,000 km keeps nearly all of it, and the elastic
    /// lid 30 km thick holds up one of 200 km; at 900 K the crust keeps everything for a
    /// billion years.
    #[test]
    fn a_warm_lower_crust_relaxes_the_relief_between_its_lid_and_its_breadth() {
        let hot = Shell {
            thickness: 60e3,
            elastic: 30e3,
            surface_temperature: 214.0,
            basal_temperature: 1_550.0,
            gravity: 3.71,
            viscosity_scale: 1.0,
            rock: true,
            strength: STRENGTH_ROCK,
        };
        assert!(survival(1_000e3, 0.1, &hot) < 0.5, "{}", survival(1_000e3, 0.1, &hot));
        assert!(survival(10_000e3, 0.1, &hot) > 0.95, "{}", survival(10_000e3, 0.1, &hot));
        assert_eq!(survival(200e3, 0.1, &hot), 1.0);
        let cold = Shell {
            basal_temperature: 900.0,
            ..hot
        };
        assert!(survival(1_000e3, 1.0, &cold) > 0.999, "{}", survival(1_000e3, 1.0, &cold));
        // The viscosity falls with the temperature, to the floor near the solidus.
        assert!(crust_viscosity(1_100.0, 3.71) > 1e3 * crust_viscosity(1_400.0, 3.71));
        assert_eq!(crust_viscosity(2_000.0, 3.71), CRUST_FLOW_MIN_VISCOSITY);
    }

    /// F184: under a rocky lid each cell relaxes by its own crust's survival over the step's
    /// substeps: a crust 65 km thick on Venus just after its overturn loses its thousand
    /// kilometres where one 40 km thick keeps half of them, and later neither relaxes.
    #[test]
    fn a_cell_relaxes_under_its_own_crust() {
        let (g, _) = crate::presets::preset("Venus").unwrap().to_graph();
        let thermal = crate::thermal::ThermalHistory::from_graph(&g);
        let surface = g.f64("atmosphere.surface_temperature_mean");
        let gravity = g.f64("body.gravity");
        let shell_at = move |t: f64, crust: f64| {
            let state = thermal.at(t);
            let share = (crust / state.lithosphere.max(1.0)).min(1.0);
            Shell {
                thickness: crust,
                elastic: state.elastic,
                surface_temperature: surface,
                basal_temperature: surface + (state.mantle_temperature - surface) * share,
                gravity,
                viscosity_scale: 1.0,
                rock: true,
                strength: STRENGTH_ROCK,
            }
        };
        let overturn = crate::thermal::ThermalHistory::from_graph(&g)
            .mantle
            .unwrap()
            .last_overturn()
            .unwrap();
        let steps: Vec<(f64, f64)> = (0..10)
            .map(|i| (overturn + 0.005 * (i as f64 + 0.5), 0.005))
            .collect();
        let crust = [40e3f32, 65e3];
        let rock = RockCrust {
            crust: &crust,
            shell_at: &shell_at,
            steps: &steps,
        };
        let thin = rock_survival(1_000e3, &rock, 40e3);
        let thick = rock_survival(1_000e3, &rock, 65e3);
        assert!(thick < 0.1 && thin > 0.3, "{thick} {thin}");
        let later: Vec<(f64, f64)> = (0..10)
            .map(|i| (overturn + 0.3 + 0.005 * (i as f64 + 0.5), 0.005))
            .collect();
        let rock = RockCrust {
            steps: &later,
            ..rock
        };
        assert!(rock_survival(1_000e3, &rock, 65e3) > 0.99);
    }

    #[test]
    fn a_thinned_shell_relaxes_the_long_wavelengths_beside_a_terrain_first() {
        // F66 with F139: under Enceladus's 22 km shell thinned to a fifth for an episode of
        // 10 Myr, its lid's flexural wavelength falls from 370 to 110 km: a 300 km load
        // relaxes and a 12 km crater keeps its relief; the shell as it is holds up most of
        // the 300 km load and all of the crater.
        let shell = Shell {
            thickness: 22e3,
            elastic: 11.6e3,
            surface_temperature: 61.0,
            basal_temperature: 273.0,
            gravity: 0.113,
            viscosity_scale: 1.0,
            rock: false,
            strength: STRENGTH_ICE,
        };
        let flexural = flexural_wavelength(&shell);
        assert!((330e3..410e3).contains(&flexural), "{flexural}");
        assert!(thinned_survival(300e3, 0.01, &shell, 0.2) < 0.05);
        assert!(thinned_survival(12e3, 0.01, &shell, 0.2) > 0.99);
        assert!((thinned_survival(300e3, 0.01, &shell, 1.0) - 1.0).abs() < 1e-9);
        assert!(survival(300e3, 0.01, &shell) > 0.6);
        assert!(survival(12e3, 4.5, &shell) > 0.999);
        // Thinned to a tenth, its lid of 1.2 km lets a crater of 60 km relax in part and one
        // of 100 km wholly, as Enceladus's relaxed craters want lids of a kilometre or two
        // (Bland et al. 2012).
        assert!(thinned_survival(60e3, 0.01, &shell, 0.1) < 0.8);
        assert!(thinned_survival(100e3, 0.01, &shell, 0.1) < 0.05);
    }

    /// F139: Mimas's lid over its young ocean, 12 km thick, holds Herschel up: its flexural
    /// wavelength on a body of so little gravity is some 450 km.
    #[test]
    fn a_small_moon_s_lid_holds_its_largest_crater() {
        let mimas = Shell {
            thickness: 24e3,
            elastic: 12.3e3,
            surface_temperature: 69.0,
            basal_temperature: 273.0,
            gravity: 0.064,
            viscosity_scale: 1.0,
            rock: false,
            strength: STRENGTH_ICE,
        };
        let flexural = flexural_wavelength(&mimas);
        assert!((400e3..500e3).contains(&flexural), "{flexural}");
        assert!(survival(130e3, 0.02, &mimas) > 0.95);
        assert!(survival(260e3, 0.02, &mimas) > 0.85);
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
            let bowl = if r < 1.0 {
                -1000.0 * (1.0 - r * r)
            } else {
                0.0
            };
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
                held.decay(&mut h, &[], |l, _, _, v| v * factors[l.min(factors.len() - 1)]);
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
        assert!(
            worst < 1e-2,
            "kept bands changed by {worst} m over 60 steps"
        );
        // Far outside the crater, beyond two radii, nothing.
        let far = |h: &[f32]| -> f32 {
            (0..grid.cell_count())
                .filter(|&k| {
                    crate::geom::angle_between(center, grid.center(grid.cell(k))) > 2.5 * radius
                })
                .map(|k| h[k].abs())
                .fold(0.0f32, f32::max)
        };
        assert!(
            far(&many) < 30.0,
            "ripples of {} m beyond the crater",
            far(&many)
        );
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
        // Ganymede today: an 83 km shell at 106 K. A 300 km load would go within a million
        // years were its lid not holding it, a 20 km crater outlasts the solar system.
        let ganymede = Shell {
            thickness: 83e3,
            elastic: 32e3,
            surface_temperature: 106.0,
            basal_temperature: 273.0,
            gravity: 1.43,
            viscosity_scale: 1.0,
            rock: false,
            strength: STRENGTH_ICE,
        };
        assert!(relaxation_time_gyr(300e3, &ganymede) < 1e-3);
        assert!(relaxation_time_gyr(20e3, &ganymede) > 100.0);
        // A thin shell, an 8 km one on Europa, flows fast under even a 30 km crater, but its
        // 3 km lid holds the crater up (F139): only loads beyond its 74 km flexural
        // wavelength settle within the surface's age.
        let europa = Shell {
            thickness: 8e3,
            elastic: 3e3,
            surface_temperature: 92.0,
            basal_temperature: 273.0,
            gravity: 1.31,
            viscosity_scale: 1.0,
            rock: false,
            strength: STRENGTH_ICE,
        };
        assert!(relaxation_time_gyr(30e3, &europa) < 0.05);
        assert!(survival_since(30e3, 4.0, 4.5, &|_| europa, 1.0) > 0.95);
        assert!(survival_since(300e3, 4.0, 4.5, &|_| europa, 1.0) < 0.05);
        assert!(survival_since(20e3, 0.0, 4.5, &|_| ganymede, 1.0) > 0.95);
        // The elastic lid holds a 60 km crater on today's Ganymede however warm the base,
        // and lets a basin of 1000 km settle.
        assert!(survival_since(60e3, 0.0, 4.5, &|_| ganymede, 1.0) > 0.95);
        assert!(survival_since(1000e3, 4.0, 4.5, &|_| ganymede, 1.0) < 0.05);
        // Where the shell is a fifth as thick a 200 km load goes.
        assert!(survival_since(200e3, 0.0, 4.5, &|_| ganymede, 0.2) < 0.5);
        assert!(survival_since(200e3, 0.0, 4.5, &|_| ganymede, 1.0) > 0.9);
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
        // The relief about the mean, by area, and the mean: the flow spreads the plateau over
        // the sphere and keeps what it is made of (F184).
        let measure = |h: &Field<f32>| -> (f32, f32) {
            let (mut sum, mut area) = (0.0f64, 0.0f64);
            for (k, &v) in h.as_slice().iter().enumerate() {
                let a = grid.cell_area(grid.cell(k));
                sum += v as f64 * a;
                area += a;
            }
            let mean = sum / area;
            let spread: f64 = h
                .as_slice()
                .iter()
                .enumerate()
                .map(|(k, &v)| (v as f64 - mean).abs() * grid.cell_area(grid.cell(k)))
                .sum::<f64>()
                / area;
            (spread as f32, mean as f32)
        };
        let (before, mean) = measure(&h);
        relax(
            &mut Bands::default(),
            &mut h,
            r,
            0.05,
            &ganymede,
            &Local::default(),
        );
        let (after, kept) = measure(&h);
        assert!(after < 0.5 * before, "{after} vs {before}");
        assert!((kept - mean).abs() < 1.0, "{kept} vs {mean}");
    }

    /// F195: an ice lid holds a band's relief only up to what bends it to its strength, and
    /// the warm ice under it flows from under the rest; a relief well within the bound, or
    /// a lid of no bound, keeps its survival.
    #[test]
    fn an_ice_lid_holds_no_more_relief_than_its_strength_bears() {
        let enceladus = Shell {
            thickness: 22e3,
            elastic: 11.5e3,
            surface_temperature: 75.0,
            basal_temperature: 273.0,
            gravity: 0.113,
            viscosity_scale: 1.0,
            rock: false,
            strength: STRENGTH_ICE,
        };
        // About a kilometre at 200 km on Enceladus, four times that at half the wavelength.
        let most = lid_relief(200e3, &enceladus);
        assert!(most > 0.8e3 && most < 1.3e3, "{most}");
        assert!((lid_relief(100e3, &enceladus) / most - 4.0).abs() < 1e-9);
        let share = 1.0 - ductile_share(200e3, &enceladus);
        let flow = unflowed(200e3, 0.05, &enceladus);
        let kept = kept_relief(4e3, share, flow, most) as f64;
        assert!(kept > 0.0 && kept <= most, "{kept}");
        let survives = survival(200e3, 0.05, &enceladus);
        let small = kept_relief(10.0, share, flow, most) as f64;
        assert!((small / 10.0 - survives).abs() < 0.01, "{small}");
        let free = kept_relief(4e3, share, flow, f64::INFINITY) as f64;
        assert!((free / 4e3 - survives).abs() < 1e-6, "{free}");
        // A rocky lid's relief stands on its crust's compensation: no bound.
        let rock = Shell {
            rock: true,
            strength: STRENGTH_ROCK,
            ..enceladus
        };
        assert!(lid_relief(200e3, &rock).is_infinite());
    }
}
