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
        let n = THERMAL_SAMPLES;
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

    /// A curve of one quantity over the body's life.
    pub fn curve(&self, f: impl Fn(&ThermalState) -> f64) -> Curve {
        Curve::new(self.samples.iter().map(|s| (s.time_gyr, f(s))).collect())
    }
}

/// The interior at `t` Gyr since formation from the present-day relations.
pub fn state_at(inp: &ThermalInputs, t: f64) -> ThermalState {
    let area = 4.0 * std::f64::consts::PI * inp.radius * inp.radius;
    let radiogenic = radiogenic_per_kg(t) * inp.rock_mass / area * (1.0 + inp.secular_fraction);
    // Tidal power scales with the square of the eccentricity (F13).
    let e = inp.eccentricity_history.eval(t).max(0.0);
    let tidal = if inp.eccentricity_now > 0.0 {
        inp.tidal_flux_now * (e / inp.eccentricity_now).powi(2)
    } else {
        inp.tidal_flux_now
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
}
