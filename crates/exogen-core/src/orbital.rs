// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! The thermal-orbital history of an icy moon locked to its planet (F138).
//!
//! A moon's eccentricity and its interior change together. The heat its tides make is the
//! eccentricity's orbital energy, G M m e² / 2a, spent: a free eccentricity decays by the
//! moon's own dissipation, its e-folding time 2E / P, which for Mimas with the porous core
//! and thin shell a steady state would give it is 0.4 Myr, against Meyer and Wisdom's
//! (2008) 325 Myr for a cold Mimas. A resonance that holds a moon now holds its forced
//! eccentricity. One that captured it in the past supplies the energy the planet's tides
//! feed it, and the eccentricity rises until the moon dissipates as much: in equilibrium
//! its heating is the resonance's, whatever the moon is made of (Meyer and Wisdom 2007),
//! and a moon that dissipates little is pumped to a higher eccentricity than one that
//! dissipates much. After the resonance lets go the eccentricity is free again.
//!
//! The moon answers with its ice. Its shell thickens and thins with the heat that reaches
//! it, against the latent heat of the water it freezes and melts, from a warm start: the
//! H2O layer liquid under the shell its heat would hold, or a thin crust where it would
//! hold none. Its porous rocky interior dissipates only while its pores hold water (F68),
//! which they do while an ocean touches it: they freeze when the moon freezes through, and
//! thaw again only once an ocean has lain over it for the time heat takes to diffuse
//! through it, R² / κ. The shell's own tides give it a threshold: a hard shell's heating
//! and its conduction both go as the inverse of its thickness, so that above one
//! eccentricity it thins away and below it freezes through, and ice frozen to the rock,
//! whose tides are those of an elastic sphere, melts only at a higher one still. So a
//! moon frozen through is held frozen by any eccentricity a resonance does not pump past
//! that, and one whose porous interior holds water keeps an ocean at a far lower one: the
//! state a moon is in depends on its history, and the episodes of its heating come of its
//! resonances and its ice together.

use crate::constants::*;
use crate::params::derive::physics::{radiogenic_per_kg, tidal_power};
use crate::params::value::Curve;
use crate::params::{Ctx, ParamGraph};
use crate::stages::tidal::ShellTides;

/// A resonance's thermal-orbital cycle (F77), as the factor it puts on the eccentricity at
/// a time since formation: one before its onset, then swinging by its amplitude about one,
/// peaking a whole number of periods from its last peak and growing to its full swing over
/// its first period.
#[derive(Clone, Copy, Debug)]
pub struct ResonanceCycle {
    pub period: f64,
    pub onset: f64,
    pub amplitude: f64,
    pub last_peak: f64,
}

impl ResonanceCycle {
    /// From the cycle's period, Gyr, how long before the present it began, its amplitude
    /// and the share of a period since its last peak; none when any of these is nil.
    pub fn new(
        period: f64,
        onset_before: f64,
        amplitude: f64,
        phase: f64,
        age: f64,
    ) -> Option<Self> {
        let amplitude = amplitude.clamp(0.0, 0.9);
        let before = onset_before.min(age);
        if period <= 0.0 || amplitude <= 0.0 || before <= 0.0 {
            return None;
        }
        Some(Self {
            period,
            onset: age - before,
            amplitude,
            last_peak: age - phase.clamp(0.0, 1.0) * period,
        })
    }

    pub fn factor(&self, t: f64) -> f64 {
        if t < self.onset {
            return 1.0;
        }
        let grown = ((t - self.onset) / self.period).min(1.0);
        let angle = std::f64::consts::TAU * (t - self.last_peak) / self.period;
        1.0 + self.amplitude * grown * angle.cos()
    }
}

/// A past resonance (F50, F138): it captures the moon at `start` and lets it go at `end`,
/// Gyr since formation, supplying its eccentricity `supply` W, which rises from nothing over
/// `rise` Gyr.
#[derive(Clone, Copy, Debug)]
pub struct Episode {
    pub start: f64,
    pub end: f64,
    pub rise: f64,
    pub supply: f64,
}

impl Episode {
    fn supply_at(&self, t: f64) -> f64 {
        if t < self.start || t >= self.end {
            return 0.0;
        }
        self.supply * ((t - self.start) / self.rise.max(1e-12)).min(1.0)
    }
}

/// What drives the eccentricity: its present value, whether a resonance holds it now
/// (forced, with the resonance's cycle), and a past resonance.
#[derive(Clone, Debug)]
pub struct Forcing {
    pub eccentricity_now: f64,
    pub forced: bool,
    pub cycle: Option<ResonanceCycle>,
    pub episode: Option<Episode>,
}

/// The interior an icy moon's thermal-orbital history evolves.
#[derive(Clone, Debug)]
pub struct OrbitalInputs {
    pub age_gyr: f64,
    pub tides: ShellTides,
    pub area: f64,
    /// The tides the porous interior dissipates while its pores hold water, per e², W (F68).
    pub core_per_e2: f64,
    /// How long an ocean must lie over the porous interior before its pores thaw, Gyr.
    pub core_thaw: f64,
    pub rock_mass: f64,
    pub secular_fraction: f64,
    /// G M m / (2a): the orbital energy of the eccentricity per e², J.
    pub energy_per_e2: f64,
}

/// The shell and the porous interior of an icy moon through its life, sampled every
/// `step` Gyr from formation.
#[derive(Clone, Debug)]
pub struct ShellHistory {
    pub step: f64,
    pub eccentricity: Vec<f64>,
    pub thickness: Vec<f64>,
    pub wet: Vec<bool>,
}

/// One sample of the history.
#[derive(Clone, Copy, Debug)]
pub struct ShellState {
    pub eccentricity: f64,
    pub thickness: f64,
    pub wet: bool,
}

impl ShellHistory {
    /// The state at a time since formation: the shell and eccentricity interpolated, the
    /// pores' water as at the sample before.
    pub fn at(&self, t: f64) -> ShellState {
        let n = self.thickness.len();
        let x = (t / self.step).clamp(0.0, (n - 1) as f64);
        let i = (x.floor() as usize).min(n.saturating_sub(2));
        let f = (x - i as f64).clamp(0.0, 1.0);
        let j = (i + 1).min(n - 1);
        ShellState {
            eccentricity: self.eccentricity[i] + (self.eccentricity[j] - self.eccentricity[i]) * f,
            thickness: self.thickness[i] + (self.thickness[j] - self.thickness[i]) * f,
            wet: self.wet[i],
        }
    }

    /// How long the shell has floated on an ocean without freezing through, up to the
    /// present, Gyr, for an H2O layer `layer` m thick.
    pub fn ocean_age(&self, layer: f64) -> f64 {
        let n = self.thickness.len();
        let frozen = self.thickness.iter().rposition(|&d| d >= layer - 1.0);
        match frozen {
            Some(k) => (n - 1 - k) as f64 * self.step,
            None => (n - 1) as f64 * self.step,
        }
    }

    pub fn now(&self) -> ShellState {
        let n = self.thickness.len();
        ShellState {
            eccentricity: self.eccentricity[n - 1],
            thickness: self.thickness[n - 1],
            wet: self.wet[n - 1],
        }
    }
}

impl OrbitalInputs {
    /// The interior of an icy moon locked to its planet, from the graph; none for any
    /// other body (`ShellTides::from_graph`).
    pub fn from_graph(graph: &ParamGraph) -> Option<Self> {
        let tides = ShellTides::from_graph(graph)?;
        Some(Self::build(tides, &|id| graph.f64(id)))
    }

    /// The same from the graph context, mirroring `from_graph`.
    pub fn from_ctx(c: &Ctx) -> Option<Self> {
        let tides = crate::params::derive::shell_tides(c)?;
        Some(Self::build(tides, &|id| c.f64(id)))
    }

    fn build(tides: ShellTides, f: &dyn Fn(&str) -> f64) -> Self {
        let radius = f("body.radius");
        let area = 4.0 * std::f64::consts::PI * radius * radius;
        let core_radius = (radius - f("body.ice_layer_thickness")).max(0.0);
        let excess = (f("interior.core_porosity") - POROSITY_CONSOLIDATED).max(0.0);
        let core_per_e2 = if excess > 0.0 {
            tidal_power(
                POROUS_DISSIPATION * excess,
                1.0,
                f("orbit.primary_mass"),
                core_radius,
                f("orbit.mean_motion"),
                1.0,
                f("orbit.semi_major_axis"),
            )
        } else {
            0.0
        };
        let mass = f("body.mass");
        Self {
            age_gyr: f("body.age"),
            tides,
            area,
            core_per_e2,
            core_thaw: core_radius * core_radius / CORE_THERMAL_DIFFUSIVITY / GYR,
            rock_mass: mass * f("body.rock_mass_fraction"),
            secular_fraction: f("interior.secular_cooling_fraction"),
            energy_per_e2: G * f("orbit.primary_mass") * mass
                / (2.0 * f("orbit.semi_major_axis").max(1.0)),
        }
    }

    /// The rock's radiogenic heat and secular cooling at a time since formation, W/m², as
    /// the thermal history has them.
    pub fn basal_flux(&self, t: f64) -> f64 {
        radiogenic_per_kg(t) * self.rock_mass / self.area * (1.0 + self.secular_fraction)
    }

    /// The tides the moon dissipates in a state, per e², W.
    fn per_e2(&self, thickness: f64, wet: bool, floating: bool) -> f64 {
        self.tides.heating(thickness, 1.0, floating) * self.area
            + if wet { self.core_per_e2 } else { 0.0 }
    }

    /// The tides the moon dissipates at an eccentricity in the steady state of a warm moon,
    /// its porous interior holding water and its shell where its conduction carries the
    /// heat, W: the stateless model of F68 and F75.
    pub fn warm_power(&self, eccentricity: f64, t: f64) -> f64 {
        let core = self.core_per_e2 * eccentricity * eccentricity;
        let shell = self
            .tides
            .equilibrium(self.basal_flux(t) + core / self.area, eccentricity);
        core + shell.heating * self.area
    }

    pub fn layer(&self) -> f64 {
        self.tides.layer.min(self.tides.radius)
    }

    /// The tides of a state at an eccentricity, W/m².
    pub fn tidal_flux(&self, s: &ShellState, eccentricity: f64) -> f64 {
        let floating = s.thickness < self.layer() - 1.0;
        let core = if s.wet {
            self.core_per_e2 / self.area
        } else {
            0.0
        };
        (core + self.tides.heating(s.thickness, 1.0, floating)) * eccentricity * eccentricity
    }

    /// The heat leaving the surface in a state, W/m²: what the shell conducts while it
    /// floats, and the heat made below where the moon is frozen through.
    pub fn surface_flux(&self, s: &ShellState, eccentricity: f64, t: f64) -> f64 {
        if s.thickness < self.layer() - 1.0 {
            self.tides.conduction / s.thickness.max(1.0)
        } else {
            self.basal_flux(t) + self.tidal_flux(s, eccentricity)
        }
    }

    /// The eccentricity over the moon's life with its interior, from the forcing: a free
    /// moon's last resonance lets it go at the latest time that lets its eccentricity decay
    /// to its present value by the present. The present value does not fall steadily as
    /// the escape comes earlier: a moon let go before its ice melts keeps the eccentricity
    /// it was pumped to for a long time, one let go after decays it fast, so the escapes are
    /// scanned back from the latest before the last is narrowed down.
    pub fn eccentricity_history(&self, forcing: &Forcing) -> Curve {
        let mut run = self.run(forcing, None);
        let held = forcing.forced || forcing.cycle.is_some();
        let target = forcing.eccentricity_now * (1.0 + ESCAPE_TOLERANCE);
        let present = |r: &ShellHistory| r.eccentricity.last().copied().unwrap_or(0.0);
        if let (false, Some(ep)) = (held, forcing.episode)
            && present(&run) > target
        {
            let ended = |end: f64| Forcing {
                episode: Some(Episode { end, ..ep }),
                ..forcing.clone()
            };
            let span = ep.end - ep.start;
            let (mut lo, mut hi) = (ep.start, ep.end);
            for k in 1..=ESCAPE_SCAN {
                let end = ep.end - span * k as f64 / ESCAPE_SCAN as f64;
                if present(&self.run(&ended(end), None)) <= target {
                    lo = end;
                    break;
                }
                hi = end;
            }
            for _ in 0..ESCAPE_ITERATIONS {
                let mid = 0.5 * (lo + hi);
                if present(&self.run(&ended(mid), None)) > target {
                    hi = mid;
                } else {
                    lo = mid;
                }
            }
            run = self.run(&ended(lo), None);
        }
        compress(run.step, &run.eccentricity)
    }

    /// The shell and the porous interior through the moon's life under an eccentricity
    /// history.
    pub fn shell_history(&self, eccentricity: &Curve) -> ShellHistory {
        let forcing = Forcing {
            eccentricity_now: eccentricity.eval(self.age_gyr),
            forced: true,
            cycle: None,
            episode: None,
        };
        self.run(&forcing, Some(eccentricity))
    }

    /// The forward integration: the eccentricity from the forcing, or following a history
    /// when one is given, and the shell and the pores' water.
    fn run(&self, forcing: &Forcing, follow: Option<&Curve>) -> ShellHistory {
        let age = self.age_gyr.max(1e-6);
        let n = ((age / ORBITAL_STEP_GYR).ceil() as usize).max(1);
        let dt = age / n as f64;
        let seconds = dt * GYR;
        let layer = self.layer();
        let latent = DENSITY_CRUST_ICE * ICE_LATENT_HEAT;
        let held = forcing.forced || forcing.cycle.is_some();
        let e_now = forcing.eccentricity_now;
        let baseline = |t: f64| match forcing.cycle {
            Some(c) => e_now * c.factor(t) / c.factor(age).max(1e-6),
            None => e_now,
        };
        // A free moon keeps its present eccentricity from its last escape on: the remnant of
        // an excitation the forcing does not state.
        let floor_from = forcing.episode.map_or(0.0, |ep| ep.end);
        // A free moon with a resonance stated has, before it, no eccentricity: whatever it
        // formed with, its warm interior, its pores full of water, damps within a few Myr.
        // One with none stated keeps its present eccentricity from the start.
        let e0 = match follow {
            Some(c) => c.eval(0.0).max(0.0),
            None if held => baseline(0.0),
            None if forcing.episode.is_some() => 0.0,
            None => e_now,
        };
        // The warm start: the H2O layer liquid under the shell the heat of a warm moon would
        // hold, or a thin crust where it would hold none.
        let start = self.tides.equilibrium(
            self.basal_flux(0.0) + self.core_per_e2 * e0 * e0 / self.area,
            e0,
        );
        let mut d = if start.floating {
            start.thickness
        } else {
            WARM_START_SHELL.min(layer)
        };
        let mut wet = true;
        let mut floating_since = Some(0.0);
        let mut x = e0 * e0;
        let mut out = ShellHistory {
            step: dt,
            eccentricity: Vec::with_capacity(n + 1),
            thickness: Vec::with_capacity(n + 1),
            wet: Vec::with_capacity(n + 1),
        };
        out.eccentricity.push(e0);
        out.thickness.push(d);
        out.wet.push(wet);
        for k in 0..n {
            let t = k as f64 * dt;
            let tm = t + 0.5 * dt;
            let floating = d < layer - 1.0;
            if floating {
                let since = *floating_since.get_or_insert(t);
                if !wet && t - since >= self.core_thaw {
                    wet = true;
                }
            } else {
                wet = false;
                floating_since = None;
            }
            let e = match follow {
                Some(c) => c.eval(t + dt).max(0.0),
                None => {
                    let p_hat = self.per_e2(d, wet, floating).max(1e-30);
                    let decay = (-p_hat / self.energy_per_e2 * seconds).exp();
                    let supply = forcing.episode.map_or(0.0, |ep| ep.supply_at(tm));
                    if supply > 0.0 {
                        let target = (supply / p_hat).min(PUMPED_ECCENTRICITY_MAX.powi(2));
                        x = target + (x - target) * decay;
                    } else {
                        x *= decay;
                    }
                    let floor = if held {
                        baseline(t + dt).powi(2)
                    } else if tm >= floor_from {
                        e_now * e_now
                    } else {
                        0.0
                    };
                    x = x.max(floor);
                    x.sqrt()
                }
            };
            // The shell: what reaches its base and what it makes itself, against what it
            // conducts, the difference freezing or melting it.
            let below = self.basal_flux(tm)
                + if wet {
                    self.core_per_e2 * e * e / self.area
                } else {
                    0.0
                };
            if floating {
                d = self.stefan_step(d, e, below, latent, seconds, layer);
            } else {
                let excess =
                    below + self.tides.heating(layer, e, false) - self.tides.conduction / layer;
                if excess > 0.0 {
                    d = (layer - excess * seconds / latent).max(ORBITAL_SHELL_MIN);
                }
            }
            out.eccentricity.push(e);
            out.thickness.push(d);
            out.wet.push(wet);
        }
        out
    }

    /// One implicit step of a floating shell's thickness: the root of
    /// ρ L (d₁ − d) / Δt = k ΔT / d₁ − q − H(d₁) nearest the shell as it is, or frozen
    /// through when even the whole layer would freeze within the step.
    fn stefan_step(
        &self,
        d: f64,
        e: f64,
        below: f64,
        latent: f64,
        seconds: f64,
        layer: f64,
    ) -> f64 {
        let k = self.tides.conduction;
        let g = |d1: f64| {
            latent * (d1 - d) / seconds - (k / d1 - below - self.tides.heating(d1, e, true))
        };
        let (mut lo, mut hi);
        if g(d) > 0.0 {
            // Thinning: down until the residual changes sign.
            hi = d;
            lo = d / 1.1;
            while g(lo) > 0.0 {
                if lo <= ORBITAL_SHELL_MIN {
                    return ORBITAL_SHELL_MIN;
                }
                hi = lo;
                lo = (lo / 1.1).max(ORBITAL_SHELL_MIN);
            }
        } else {
            // Thickening: up until it changes sign, or frozen through.
            lo = d;
            hi = (d * 1.1).min(layer);
            while g(hi) < 0.0 {
                if hi >= layer {
                    return layer;
                }
                lo = hi;
                hi = (hi * 1.1).min(layer);
            }
        }
        for _ in 0..40 {
            let mid = 0.5 * (lo + hi);
            if g(mid) > 0.0 {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        0.5 * (lo + hi)
    }
}

impl Forcing {
    /// The forcing from the graph context, for a moon whose interior is `inputs`.
    pub fn from_ctx(c: &Ctx, inputs: &OrbitalInputs) -> Self {
        Self::build(inputs, &|id| c.f64(id), c.bool("orbit.resonance_forced"))
    }

    pub fn from_graph(graph: &ParamGraph, inputs: &OrbitalInputs) -> Self {
        Self::build(
            inputs,
            &|id| graph.f64(id),
            graph.bool("orbit.resonance_forced"),
        )
    }

    fn build(inputs: &OrbitalInputs, f: &dyn Fn(&str) -> f64, forced: bool) -> Self {
        let age = inputs.age_gyr;
        let cycle = ResonanceCycle::new(
            f("orbit.cycle_period"),
            f("orbit.cycle_onset"),
            f("orbit.cycle_amplitude"),
            f("orbit.cycle_phase"),
            age,
        );
        let before = f("orbit.resonance_age");
        let width = f("orbit.resonance_duration");
        let eccentricity = f("orbit.resonance_eccentricity");
        let heating = f("orbit.resonance_heating");
        let episode =
            (width > 0.0 && before > 0.0 && before <= age && (eccentricity > 0.0 || heating > 0.0))
                .then(|| {
                    let centre = age - before;
                    let half = 0.5 * width;
                    // Stated by its heat, the resonance supplies that; stated by its eccentricity,
                    // what that eccentricity gives a warm moon, as the dynamics that found it
                    // assumed of the moon.
                    let supply = if heating > 0.0 {
                        heating * inputs.area
                    } else {
                        inputs.warm_power(eccentricity, centre)
                    };
                    Episode {
                        start: (centre - half).max(0.0),
                        end: (centre + half).min(age),
                        rise: RESONANCE_RAMP_SHARE * half,
                        supply,
                    }
                })
                .filter(|ep| ep.end > ep.start && ep.supply > 0.0);
        Self {
            eccentricity_now: f("orbit.eccentricity"),
            forced,
            cycle,
            episode,
        }
    }
}

/// The samples as a curve, a run of equal values kept by its ends.
fn compress(step: f64, values: &[f64]) -> Curve {
    let n = values.len();
    let mut points = Vec::with_capacity(n);
    for (k, &v) in values.iter().enumerate() {
        let same_before = k > 0 && values[k - 1] == v;
        let same_after = k + 1 < n && values[k + 1] == v;
        if !(same_before && same_after) {
            points.push((k as f64 * step, v));
        }
    }
    Curve::new(points)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::presets::preset;

    fn history(name: &str) -> (crate::params::ParamGraph, OrbitalInputs, ShellHistory) {
        let (g, _) = preset(name).unwrap().to_graph();
        let o = OrbitalInputs::from_graph(&g).unwrap();
        let h = o.shell_history(g.curve("orbit.eccentricity_history").unwrap());
        (g, o, h)
    }

    /// Diagnostic: each icy moon's eccentricity, shell, pores and heat through its life;
    /// run with `--ignored --nocapture`.
    #[test]
    #[ignore]
    fn print_histories() {
        for name in crate::presets::PRESET_NAMES {
            let (g, _) = preset(name).unwrap().to_graph();
            let Some(o) = OrbitalInputs::from_graph(&g) else {
                continue;
            };
            let e = g.curve("orbit.eccentricity_history").unwrap();
            let h = o.shell_history(e);
            let age = o.age_gyr;
            println!("{name}: {} points", e.points.len());
            for before in [
                4.5, 4.0, 3.0, 2.5, 2.4, 2.3, 2.2, 2.1, 2.0, 1.8, 1.6, 1.5, 1.4, 1.2, 1.15, 1.1,
                1.05, 1.0, 0.95, 0.9, 0.85, 0.8, 0.7, 0.5, 0.2, 0.1, 0.08, 0.06, 0.05, 0.04, 0.03,
                0.02, 0.01, 0.005, 0.0,
            ] {
                let t = (age - before).max(0.0);
                let s = h.at(t);
                println!(
                    "  {before:5.3} Ga: e {:.4} shell {:6.1} km wet {:5} tides {:7.2} heat {:7.2} mW/m2",
                    s.eccentricity,
                    s.thickness / 1e3,
                    s.wet,
                    o.tidal_flux(&s, s.eccentricity) * 1e3,
                    o.surface_flux(&s, s.eccentricity, t) * 1e3
                );
            }
        }
    }

    /// Enceladus, held in its resonance with Dione from its warm start, keeps its porous
    /// interior wet and its shell near 23 km through its life (F68, F138).
    #[test]
    fn enceladus_keeps_its_ocean_and_wet_core() {
        let (g, _, h) = history("Enceladus");
        let age = g.f64("body.age");
        for t in [0.5, 2.0, age] {
            let s = h.at(t);
            assert!(s.wet, "{t}");
            assert!((18e3..28e3).contains(&s.thickness), "{t}: {}", s.thickness);
        }
    }

    /// Miranda freezes through within half a gigayear of forming; its resonance pumps its
    /// eccentricity until its ice melts, holds it at the heat the coronae's faults record
    /// (31 to 112 mW/m², Beddingfield et al. 2015), and lets it freeze again (F138).
    #[test]
    fn miranda_melts_in_its_resonance_and_freezes_again() {
        let (g, o, h) = history("Miranda");
        let age = g.f64("body.age");
        let layer = o.layer();
        assert!(h.at(0.5).thickness >= layer - 1.0);
        let heated = (0..60)
            .map(|i| age - 1.15 + 0.3 * i as f64 / 60.0)
            .filter(|&t| {
                let s = h.at(t);
                o.tidal_flux(&s, s.eccentricity) > 0.031
            })
            .count();
        assert!(heated >= 30, "{heated} of 60 samples heated");
        let now = h.now();
        assert!(now.thickness >= layer - 1.0 && !now.wet);
    }
}
