// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! Thermal evolution through the body's life (PLAN.md 4.2): heat flux, mantle temperature,
//! lithosphere and elastic thickness, relief limit, cumulative contraction and magma
//! production, sampled from formation to the present with the relations of the present-day
//! interior nodes applied at each time. The coarse loop reads it at every step, and the
//! scarp, volcanism and relaxation processes draw their rates from it.

use crate::constants::*;
use crate::params::ParamGraph;
use crate::params::derive::physics::{
    conductive_thickness, contraction, ice_conductivity, magma_rate, radiogenic_per_kg,
    relief_limit,
};
use crate::params::value::Curve;

/// Samples over the body's life.
pub const THERMAL_SAMPLES: usize = 96;

/// Everything the thermal evolution depends on, gathered from the graph.
#[derive(Clone, Debug)]
pub struct ThermalInputs {
    pub age_gyr: f64,
    pub radius: f64,
    pub gravity: f64,
    pub uncompressed_density: f64,
    pub rock_mass: f64,
    pub secular_fraction: f64,
    /// Present tidal flux and the eccentricity it was computed for.
    pub tidal_flux_now: f64,
    pub eccentricity_now: f64,
    pub eccentricity_history: Curve,
    /// An icy moon's shell, which dissipates the tides of each time at the thickness of
    /// that time (F75), and the present tides of its porous interior, W/m² (F68); the
    /// shell's share of the present tidal flux is then recomputed, not scaled.
    pub shell: Option<crate::stages::tidal::ShellTides>,
    pub core_tidal_flux_now: f64,
    pub icy: bool,
    pub surface_temperature: f64,
    pub melting_point: f64,
    pub crust_strength: f64,
    pub magma_efficiency_base: f64,
}

/// The interior at one time.
#[derive(Clone, Copy, Debug, Default)]
pub struct ThermalState {
    pub time_gyr: f64,
    /// Surface heat flux, W/m².
    pub heat_flux: f64,
    pub tidal_flux: f64,
    pub mantle_temperature: f64,
    pub lithosphere: f64,
    pub elastic: f64,
    /// Relief the lithosphere supports, metres.
    pub relief_limit: f64,
    /// Cumulative radial contraction since formation, metres.
    pub contraction: f64,
    /// Magma production, m³/s.
    pub magma_rate: f64,
}

#[derive(Clone, Debug)]
pub struct ThermalHistory {
    pub samples: Vec<ThermalState>,
}

/// The melt a body's mantle has made per unit area since it formed, m, through its life
/// (F120): the clock a plume's rise and its volcanoes keep, growing as fast as their
/// mantle melts. Past the present it runs on at the present's rate.
#[derive(Clone, Debug)]
pub struct MeltClock {
    times: Vec<f64>,
    melt: Vec<f64>,
}

impl MeltClock {
    /// The rate at the end, m per Gyr.
    fn last_rate(&self) -> f64 {
        let n = self.times.len();
        if n < 2 {
            return 0.0;
        }
        (self.melt[n - 1] - self.melt[n - 2]) / (self.times[n - 1] - self.times[n - 2]).max(1e-12)
    }

    /// The melt per unit area made by a time, m.
    pub fn at(&self, t: f64) -> f64 {
        let n = self.times.len();
        if t <= self.times[0] {
            return 0.0;
        }
        if t >= self.times[n - 1] {
            return self.melt[n - 1] + self.last_rate() * (t - self.times[n - 1]);
        }
        let k = self.times.partition_point(|&x| x <= t).clamp(1, n - 1);
        let f = (t - self.times[k - 1]) / (self.times[k] - self.times[k - 1]).max(1e-12);
        self.melt[k - 1] + f * (self.melt[k] - self.melt[k - 1])
    }

    /// The time by which the melt per unit area reached an amount, Gyr; never, past the
    /// present, if the mantle no longer melts.
    pub fn time_of(&self, m: f64) -> f64 {
        let n = self.times.len();
        if m <= 0.0 {
            return self.times[0];
        }
        if m >= self.melt[n - 1] {
            let rate = self.last_rate();
            return if rate > 0.0 {
                self.times[n - 1] + (m - self.melt[n - 1]) / rate
            } else {
                f64::INFINITY
            };
        }
        let k = self.melt.partition_point(|&x| x < m).clamp(1, n - 1);
        let span = self.melt[k] - self.melt[k - 1];
        if span <= 0.0 {
            return self.times[k];
        }
        self.times[k - 1] + (m - self.melt[k - 1]) / span * (self.times[k] - self.times[k - 1])
    }

    /// The rate at a time, m per Gyr.
    pub fn rate(&self, t: f64) -> f64 {
        let h = 0.01;
        (self.at(t + h) - self.at((t - h).max(0.0))) / (t + h - (t - h).max(0.0))
    }
}

impl ThermalInputs {
    pub fn from_graph(graph: &ParamGraph) -> Self {
        let icy = matches!(
            graph.str("body.composition_class"),
            "rock_ice" | "ice_dominated"
        );
        let regime = graph.str("interior.tectonic_regime");
        let base = match regime {
            "mobile_lid" => MAGMA_EFFICIENCY_MOBILE,
            "episodic_lid" => MAGMA_EFFICIENCY_EPISODIC,
            "ice_shell_ocean" | "ice_shell_rock" => MAGMA_EFFICIENCY_ICE,
            _ => MAGMA_EFFICIENCY_STAGNANT,
        };
        let melting_point = if graph
            .get("volatiles.inventory")
            .as_map()
            .and_then(|m| m.get("NH3").copied())
            .unwrap_or(0.0)
            > 0.0
        {
            T_EUTECTIC_AMMONIA_WATER
        } else {
            T_MELT_WATER
        };
        Self {
            age_gyr: graph.f64("body.age"),
            radius: graph.f64("body.radius"),
            gravity: graph.f64("body.gravity"),
            uncompressed_density: graph.f64("body.uncompressed_density"),
            rock_mass: graph.f64("body.mass") * graph.f64("body.rock_mass_fraction"),
            secular_fraction: graph.f64("interior.secular_cooling_fraction"),
            tidal_flux_now: graph.f64("interior.tidal_flux"),
            eccentricity_now: graph.f64("orbit.eccentricity"),
            eccentricity_history: graph
                .curve("orbit.eccentricity_history")
                .cloned()
                .unwrap_or_else(|| Curve::constant(graph.f64("orbit.eccentricity"))),
            shell: crate::stages::tidal::ShellTides::from_graph(graph),
            core_tidal_flux_now: graph.f64("interior.core_tidal_power")
                / (4.0 * std::f64::consts::PI * graph.f64("body.radius").powi(2)),
            icy,
            surface_temperature: graph.f64("atmosphere.surface_temperature_mean"),
            melting_point,
            crust_strength: graph.f64("interior.crust_strength"),
            magma_efficiency_base: base,
        }
    }
}

impl ThermalHistory {
    pub fn from_graph(graph: &ParamGraph) -> Self {
        Self::new(&ThermalInputs::from_graph(graph))
    }

    pub fn new(inp: &ThermalInputs) -> Self {
        // As finely as the eccentricity history changes, so that a resonance's cycle
        // (F77) is resolved, and no less finely than the base sampling.
        let finest = inp
            .eccentricity_history
            .points
            .windows(2)
            .map(|w| w[1].0 - w[0].0)
            .filter(|dt| *dt > 1e-9)
            .fold(f64::INFINITY, f64::min);
        let n = if finest.is_finite() {
            ((inp.age_gyr / finest).ceil() as usize).clamp(THERMAL_SAMPLES, THERMAL_SAMPLES_MAX)
        } else {
            THERMAL_SAMPLES
        };
        let mut samples: Vec<ThermalState> = (0..=n)
            .map(|k| {
                let t = inp.age_gyr * k as f64 / n as f64;
                state_at(inp, t)
            })
            .collect();
        // A thin early elastic lithosphere carries less relief than the present one: the
        // strength limit is scaled by the elastic thickness relative to the present value,
        // so the present matches the max_relief node. Provisional until a flexural model.
        let elastic_now = samples.last().map_or(1.0, |s| s.elastic).max(1.0);
        for s in &mut samples {
            s.relief_limit *= (s.elastic / elastic_now).clamp(RELIEF_SUPPORT_MIN, 1.0);
        }
        Self { samples }
    }

    /// The state at a time since formation, interpolated between samples.
    pub fn at(&self, time_gyr: f64) -> ThermalState {
        let last = self.samples.len() - 1;
        if last == 0 {
            return self.samples[0];
        }
        let span = self.samples[last].time_gyr.max(1e-12);
        let x = (time_gyr / span * last as f64).clamp(0.0, last as f64);
        let i = (x.floor() as usize).min(last - 1);
        let f = x - i as f64;
        let (a, b) = (&self.samples[i], &self.samples[i + 1]);
        let lerp = |p: f64, q: f64| p + (q - p) * f;
        ThermalState {
            time_gyr,
            heat_flux: lerp(a.heat_flux, b.heat_flux),
            tidal_flux: lerp(a.tidal_flux, b.tidal_flux),
            mantle_temperature: lerp(a.mantle_temperature, b.mantle_temperature),
            lithosphere: lerp(a.lithosphere, b.lithosphere),
            elastic: lerp(a.elastic, b.elastic),
            relief_limit: lerp(a.relief_limit, b.relief_limit),
            contraction: lerp(a.contraction, b.contraction),
            magma_rate: lerp(a.magma_rate, b.magma_rate),
        }
    }

    /// The melt the mantle has made per unit area since formation (F120), for a body of
    /// this surface area, m².
    pub fn melt_clock(&self, area: f64) -> MeltClock {
        let mut times = Vec::with_capacity(self.samples.len());
        let mut melt = Vec::with_capacity(self.samples.len());
        let mut m = 0.0;
        for (k, s) in self.samples.iter().enumerate() {
            if k > 0 {
                let a = &self.samples[k - 1];
                m += 0.5 * (a.magma_rate + s.magma_rate) * (s.time_gyr - a.time_gyr) * GYR
                    / area.max(1.0);
            }
            times.push(s.time_gyr);
            melt.push(m);
        }
        MeltClock { times, melt }
    }

    /// A curve of one quantity over the body's life.
    pub fn curve(&self, f: impl Fn(&ThermalState) -> f64) -> Curve {
        Curve::new(self.samples.iter().map(|s| (s.time_gyr, f(s))).collect())
    }
}

/// The interior at `t` Gyr since formation from the present-day relations.
pub fn state_at(inp: &ThermalInputs, t: f64) -> ThermalState {
    let area = 4.0 * std::f64::consts::PI * inp.radius * inp.radius;
    let radiogenic = radiogenic_per_kg(t) * inp.rock_mass / area * (1.0 + inp.secular_fraction);
    // Tidal power scales with the square of the eccentricity (F13), at a fixed interior;
    // an icy shell's own share is that of the shell the heat of the time leaves (F75).
    let e = inp.eccentricity_history.eval(t).max(0.0);
    let scale = if inp.eccentricity_now > 0.0 {
        (e / inp.eccentricity_now).powi(2)
    } else {
        1.0
    };
    let tidal = match &inp.shell {
        Some(shell) => {
            let core = inp.core_tidal_flux_now * scale;
            core + shell.equilibrium(radiogenic + core, e).heating
        }
        None => inp.tidal_flux_now * scale,
    };
    let heat_flux = radiogenic + tidal;
    let mantle_temperature = if inp.icy {
        inp.melting_point
    } else {
        MANTLE_T_BASE + MANTLE_T_SCALE * (heat_flux / EARTH_HEAT_FLUX).max(0.0).sqrt()
    };
    let dt = (mantle_temperature - inp.surface_temperature).max(1.0);
    let lithosphere = if inp.icy {
        conductive_thickness(
            ice_conductivity(0.5 * (mantle_temperature + inp.surface_temperature)),
            dt,
            heat_flux,
        )
    } else {
        conductive_thickness(K_ROCK, dt, MANTLE_FLUX_FRACTION * heat_flux)
    }
    .min(inp.radius);
    let iso = if inp.icy {
        ELASTIC_ISOTHERM_ICE
    } else {
        ELASTIC_ISOTHERM_ROCK
    };
    let elastic = lithosphere * ((iso - inp.surface_temperature) / dt).clamp(0.0, 1.0);
    let rho = if inp.icy {
        DENSITY_CRUST_ICE
    } else {
        DENSITY_CRUST_ROCK
    };
    let relief_limit = relief_limit(inp.crust_strength, rho, inp.gravity);
    let tidal_efficiency = MAGMA_EFFICIENCY_TIDAL_MAX * (tidal / MAGMA_TIDAL_FLUX_SCALE).min(1.0);
    let efficiency = inp.magma_efficiency_base.max(tidal_efficiency.max(0.0));
    let heat_per_m3 = if inp.icy {
        MAGMA_HEAT_PER_M3_ICE
    } else {
        MAGMA_HEAT_PER_M3_ROCK
    };
    ThermalState {
        time_gyr: t,
        heat_flux,
        tidal_flux: tidal,
        mantle_temperature,
        lithosphere,
        elastic,
        relief_limit,
        contraction: contraction(inp.radius, inp.uncompressed_density, t),
        magma_rate: magma_rate(heat_flux, area, efficiency, heat_per_m3),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::presets::preset;

    /// Diagnostic: each rock body's melt per unit area since formation, m, at a few times
    /// (the calibration of F120); run with `--ignored --nocapture`.
    #[test]
    #[ignore]
    fn print_melt_per_area() {
        for name in ["Mars", "Venus", "Moon", "Mercury", "Earth"] {
            let (g, _) = preset(name).unwrap().to_graph();
            let h = ThermalHistory::from_graph(&g);
            let r = g.f64("body.radius");
            let area = 4.0 * std::f64::consts::PI * r * r;
            let age = g.f64("body.age");
            let n = 4000;
            let mut m = 0.0;
            let mut marks = vec![];
            let times = [0.3, 0.4, 0.9, 1.0, 1.2, 1.3, 1.4, 1.5, 4.0, 4.5];
            let mut next = 0;
            for k in 0..n {
                let t0 = age * k as f64 / n as f64;
                let t1 = age * (k + 1) as f64 / n as f64;
                let rate = h.at(0.5 * (t0 + t1)).magma_rate;
                m += rate * (t1 - t0) * 3.156e16 / area;
                while next < times.len() && t1 >= times[next] - 1e-9 {
                    marks.push((times[next], m));
                    next += 1;
                }
            }
            println!(
                "{name} ({}): {}",
                g.str("interior.tectonic_regime"),
                marks
                    .iter()
                    .map(|(t, m)| format!("{t}: {m:.0}"))
                    .collect::<Vec<_>>()
                    .join(", ")
            );
        }
    }

    /// At the present the history reproduces the present-day nodes it was built from.
    #[test]
    fn history_matches_the_present_day_nodes() {
        for name in ["Earth", "Moon", "Europa"] {
            let (g, _) = preset(name).unwrap().to_graph();
            let h = ThermalHistory::from_graph(&g);
            let now = h.at(g.f64("body.age"));
            let close = |a: f64, b: f64| (a - b).abs() <= 0.02 * b.abs().max(1e-9);
            assert!(
                close(now.heat_flux, g.f64("interior.surface_heat_flux")),
                "{name} flux {} vs {}",
                now.heat_flux,
                g.f64("interior.surface_heat_flux")
            );
            assert!(close(
                now.mantle_temperature,
                g.f64("interior.mantle_temperature")
            ));
            assert!(close(
                now.lithosphere,
                g.f64("interior.lithosphere_thickness")
            ));
            assert!(close(now.contraction, g.f64("interior.contraction_total")));
            assert!(
                close(now.magma_rate, g.f64("interior.magma_production_rate")),
                "{name} magma {} vs {}",
                now.magma_rate,
                g.f64("interior.magma_production_rate")
            );
        }
    }

    /// Heat flux was higher early, the lithosphere thinner, and contraction grows.
    #[test]
    fn the_young_moon_was_hotter_and_contraction_accumulates() {
        let (g, _) = preset("Moon").unwrap().to_graph();
        let h = ThermalHistory::from_graph(&g);
        let early = h.at(0.5);
        let now = h.at(g.f64("body.age"));
        assert!(early.heat_flux > 2.0 * now.heat_flux);
        assert!(early.lithosphere < 0.6 * now.lithosphere);
        assert!(early.relief_limit <= now.relief_limit);
        let mut last = -1.0;
        for s in &h.samples {
            assert!(s.contraction >= last);
            last = s.contraction;
        }
        // The Moon's total contraction is a few hundred metres; Mercury's kilometres.
        assert!(now.contraction > 100.0 && now.contraction < 1_000.0);
    }

    /// A scripted eccentricity history drives the tidal flux with its square.
    #[test]
    fn eccentricity_history_scales_tidal_heating() {
        let (mut g, _) = preset("Io").unwrap().to_graph();
        let e = g.f64("orbit.eccentricity");
        let age = g.f64("body.age");
        g.set(
            "orbit.eccentricity_history",
            Curve::new(vec![(0.0, 2.0 * e), (age, e)]),
        )
        .unwrap();
        let h = ThermalHistory::from_graph(&g);
        let early = h.at(0.0);
        let now = h.at(age);
        assert!((early.tidal_flux / now.tidal_flux - 4.0).abs() < 0.1);
    }

    /// The Laplace resonance's cycle (F77) swings Europa's eccentricity by a quarter about
    /// its mean every 0.12 Gyr since 1.6 Gyr ago, the present just past a peak and
    /// falling, and its shell thins and thickens with it (F75).
    #[test]
    fn europa_s_shell_breathes_with_the_laplace_cycle() {
        let (g, _) = preset("Europa").unwrap().to_graph();
        let age = g.f64("body.age");
        let e = g.curve("orbit.eccentricity_history").unwrap();
        assert!((e.eval(age) - 0.0094).abs() < 1e-6);
        let mean = e.eval(1.0);
        let expected = 0.0094 / (1.0 + 0.25 * (std::f64::consts::TAU * 0.2).cos());
        assert!((mean - expected).abs() < 1e-5, "{mean}");
        let (lo, hi) = (0..=240)
            .map(|i| e.eval(age - 0.24 * i as f64 / 240.0))
            .fold((1.0f64, 0.0f64), |(a, b), x| (a.min(x), b.max(x)));
        assert!(
            (hi / mean - 1.25).abs() < 0.02 && (lo / mean - 0.75).abs() < 0.02,
            "{lo} {hi}"
        );
        assert!(e.eval(age - 0.005) > e.eval(age));
        let h = ThermalHistory::from_graph(&g);
        let (thin, thick) = (0..=96)
            .map(|i| h.at(age - 0.24 * i as f64 / 96.0).lithosphere)
            .fold((f64::INFINITY, 0.0f64), |(a, b), x| (a.min(x), b.max(x)));
        eprintln!("Europa's shell through the cycle: {thin:.0} to {thick:.0} m");
        assert!(thin < 27e3 && thick > 36e3, "{thin} {thick}");
    }

    /// Ganymede's resonance episode (F50) pulses its eccentricity and its tidal flux two
    /// gigayears ago, and the lanes take their amount and their time from the pulse.
    #[test]
    fn ganymede_resonance_episode_heats_the_shell_two_gigayears_ago() {
        let (g, _) = preset("Ganymede").unwrap().to_graph();
        let age = g.f64("body.age");
        let e = g.curve("orbit.eccentricity_history").unwrap();
        assert!((e.eval(0.0) - 0.0013).abs() < 1e-6);
        assert!((e.eval(age) - 0.0013).abs() < 1e-6);
        assert!((e.eval(age - 2.0) - 0.0185).abs() < 1e-6);
        // The plateau holds the forced value over the middle three fifths.
        assert!((e.eval(age - 2.25) - 0.0185).abs() < 1e-6);
        assert!(e.eval(age - 2.45) < 0.012 && e.eval(age - 2.45) > 0.0013);
        let h = ThermalHistory::from_graph(&g);
        let peak = h.samples.iter().map(|s| s.tidal_flux).fold(0.0, f64::max);
        assert!(peak > 0.007 && peak < 0.009, "peak tidal flux {peak}");
        assert!(h.at(age).tidal_flux < 0.001);
        let fraction = g.f64("tectonics.ice_resurfaced_fraction");
        assert!((0.55..=0.75).contains(&fraction), "fraction {fraction}");
        let t = g.f64("tectonics.ice_resurfacing_age");
        assert!((1.9..=2.1).contains(&t), "age {t}");
        let d = g.f64("tectonics.ice_resurfacing_duration");
        assert!((0.3..=0.8).contains(&d), "duration {d}");
    }
}
