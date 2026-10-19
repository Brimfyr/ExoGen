// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! A rocky mantle that stores its heat under a lid that can break (F121).
//!
//! The mantle below the lid has one temperature, which its radioactive and tidal heat raise
//! and the heat it convects into the lid's base lowers; its heat capacity is the rock's and
//! the core's, less the lid's, which grows as the mantle cools. Under a stagnant lid the
//! convected flux follows the boundary-layer scaling of stagnant-lid convection, the
//! viscosity Arrhenian in the mantle's temperature; the surface loses that flux and the heat
//! the lid gives up as it thickens. The lid breaks when the stress the convecting sublayer
//! puts on it, along a cell as long as the mantle is deep, exceeds the lid's frictional
//! strength over its thickness (after Moresi and Solomatov 1998's comparison of convective
//! and yield stresses): it founders into the mantle, which takes up its cold, and a new lid
//! grows as a cooling half space until it is as thick as a stagnant lid would be, when it
//! can break again (Turcotte 1993: the interior reheats under a thickening lid until the lid
//! fails). Whether a body's lid ever breaks, and how often, comes out of its size, its heat,
//! its surface temperature and its lid's friction.

use crate::constants::*;
use crate::params::derive::physics::radiogenic_per_kg;
use crate::params::value::Curve;

/// What the mantle's history depends on.
#[derive(Clone, Debug)]
pub struct MantleInputs {
    pub age_gyr: f64,
    pub radius: f64,
    pub gravity: f64,
    pub mass: f64,
    pub rock_fraction: f64,
    pub core_fraction: f64,
    pub surface_temperature: f64,
    /// Present tidal flux, W/m², the eccentricity it was computed for, and the
    /// eccentricity's history; the tides scale with the eccentricity squared (F13).
    pub tidal_flux_now: f64,
    pub eccentricity_now: f64,
    pub eccentricity_history: Curve,
}

impl MantleInputs {
    pub fn from_graph(graph: &crate::params::ParamGraph) -> Self {
        Self {
            age_gyr: graph.f64("body.age"),
            radius: graph.f64("body.radius"),
            gravity: graph.f64("body.gravity"),
            mass: graph.f64("body.mass"),
            rock_fraction: graph.f64("body.rock_mass_fraction"),
            core_fraction: graph.f64("body.core_mass_fraction"),
            surface_temperature: graph.f64("atmosphere.surface_temperature_mean"),
            tidal_flux_now: graph.f64("interior.tidal_flux"),
            eccentricity_now: graph.f64("orbit.eccentricity"),
            eccentricity_history: graph
                .curve("orbit.eccentricity_history")
                .cloned()
                .unwrap_or_else(|| Curve::constant(graph.f64("orbit.eccentricity"))),
        }
    }
}

/// The mantle through the body's life, one sample per step: each step's mean surface heat
/// flux, W/m², and the mantle's temperature and the lid's thickness at its end.
#[derive(Clone, Debug, Default)]
pub struct MantleHistory {
    pub times: Vec<f64>,
    pub heat_flux: Vec<f64>,
    pub temperature: Vec<f64>,
    pub lid: Vec<f64>,
    /// When the lid broke and foundered, Gyr after formation, oldest first.
    pub overturns: Vec<f64>,
}

impl MantleHistory {
    /// The surface heat flux at a time, W/m², from the steps about it.
    pub fn heat_flux_at(&self, t: f64) -> f64 {
        let n = self.times.len();
        if n == 0 {
            return 0.0;
        }
        let k = self.times.partition_point(|&x| x < t).min(n - 1);
        self.heat_flux[k]
    }

    /// The mantle's temperature at a time, K.
    pub fn temperature_at(&self, t: f64) -> f64 {
        let n = self.times.len();
        if n == 0 {
            return 0.0;
        }
        let k = self.times.partition_point(|&x| x < t).min(n - 1);
        self.temperature[k]
    }

    /// The last time the lid foundered, if it ever did.
    pub fn last_overturn(&self) -> Option<f64> {
        self.overturns.last().copied()
    }
}

/// The stagnant-lid convective flux at a mantle temperature, W/m² (Solomatov and Moresi
/// 2000's form): a k (rho g alpha / (kappa eta))^(1/3) (R T² / E)^(4/3), the last factor the
/// rheological temperature scale over which the viscosity changes by a factor e.
fn stagnant_flux(t: f64, gravity: f64) -> f64 {
    let eta = MANTLE_VISCOSITY_REF
        * (MANTLE_ACTIVATION / GAS_CONSTANT * (1.0 / t - 1.0 / MANTLE_VISCOSITY_REF_K)).exp();
    let rheological = GAS_CONSTANT * t * t / MANTLE_ACTIVATION;
    STAGNANT_LID_PREFACTOR
        * K_ROCK
        * (DENSITY_ROCK * gravity * MANTLE_EXPANSIVITY / (MANTLE_DIFFUSIVITY * eta)).cbrt()
        * rheological.powf(4.0 / 3.0)
}

/// The flux above which the lid breaks, W/m²: where the sublayer's buoyancy stress
/// rho g alpha dT_rh delta_rh, over a cell as long as the mantle is deep, matches the lid's
/// friction mu rho g delta² / 2, with delta_rh = k dT_rh / q and delta = k dT / q.
fn breaking_flux(t: f64, surface: f64, depth: f64) -> f64 {
    let rheological = GAS_CONSTANT * t * t / MANTLE_ACTIVATION;
    LID_FRICTION * K_ROCK * (t - surface).powi(2)
        / (2.0 * MANTLE_EXPANSIVITY * rheological * rheological * depth.max(1.0))
}

/// The mantle's history from formation to the present (F121).
pub fn evolve(inp: &MantleInputs) -> MantleHistory {
    let r = inp.radius;
    let area = 4.0 * std::f64::consts::PI * r * r;
    let mantle = inp.mass * inp.rock_fraction;
    let core = inp.mass * inp.core_fraction;
    let core_radius = (3.0 * core / (4.0 * std::f64::consts::PI * DENSITY_IRON)).cbrt();
    let depth = (r - core_radius).max(1.0);
    let capacity = (mantle * MANTLE_HEAT_CAPACITY + core * CORE_HEAT_CAPACITY) / area;
    let ts = inp.surface_temperature;
    // The lid's volume per unit area, m, and the heat it holds below the mantle's, J/m².
    let shell = |d: f64| (r.powi(3) - (r - d.min(depth)).powi(3)) / (3.0 * r * r);
    let deficit =
        |t: f64, d: f64| DENSITY_ROCK * MANTLE_HEAT_CAPACITY * shell(d) * (t - ts).max(0.0) / 2.0;
    let tides = |t: f64| {
        let e = inp.eccentricity_history.eval(t).max(0.0);
        let scale = if inp.eccentricity_now > 0.0 {
            (e / inp.eccentricity_now).powi(2)
        } else {
            1.0
        };
        inp.tidal_flux_now * scale
    };
    let n = MANTLE_STEPS;
    let dt = inp.age_gyr / n as f64;
    let mut out = MantleHistory {
        times: Vec::with_capacity(n),
        heat_flux: Vec::with_capacity(n),
        temperature: Vec::with_capacity(n),
        lid: Vec::with_capacity(n),
        overturns: Vec::new(),
    };
    let mut temp = MANTLE_START_K;
    let mut lid = (K_ROCK * (temp - ts) / stagnant_flux(temp, inp.gravity)).min(0.99 * depth);
    let mut broke: Option<f64> = None;
    let mut t = 0.0;
    let mut k = 0;
    while k < n {
        let heat = radiogenic_per_kg(t + 0.5 * dt) * mantle / area + tides(t + 0.5 * dt);
        let q_stagnant = stagnant_flux(temp, inp.gravity);
        let steady = (K_ROCK * (temp - ts) / q_stagnant).min(0.99 * depth);
        // A lid that broke grows back as a cooling half space, its thickness the
        // conductive one its flux gives, sqrt(pi kappa tau), and the flux it loses over the
        // step the half space's, averaged exactly (a flux sampled once a step had made the
        // overturns' count hang on the step).
        let (grown, young, mature) = match broke {
            Some(t0) => {
                let ta = (t - t0).max(0.0) * GYR;
                let tb = ta + dt * GYR;
                let pk = std::f64::consts::PI * MANTLE_DIFFUSIVITY;
                let young = K_ROCK * (temp - ts) * 2.0 * (tb.sqrt() - ta.sqrt())
                    / (pk.sqrt() * (tb - ta));
                ((pk * tb).sqrt().min(steady), young, (pk * ta).sqrt() >= steady)
            }
            None => (steady, 0.0, true),
        };
        if mature && q_stagnant > breaking_flux(temp, ts, depth) {
            // The lid founders, and the mantle takes up its cold.
            let convecting = (capacity - DENSITY_ROCK * MANTLE_HEAT_CAPACITY * shell(lid))
                .max(0.05 * capacity);
            temp -= deficit(temp, lid) / convecting;
            broke = Some(t);
            lid = 0.0;
            out.overturns.push(t);
            continue;
        }
        let base = if grown < steady {
            young.max(q_stagnant)
        } else {
            q_stagnant
        };
        let convecting =
            (capacity - DENSITY_ROCK * MANTLE_HEAT_CAPACITY * shell(grown)).max(0.05 * capacity);
        let next = temp + (heat - base) * dt * GYR / convecting;
        // The surface loses what reaches the lid's base and what the lid gives up as its
        // cold grows.
        let surface = base + (deficit(next, grown) - deficit(temp, lid)) / (dt * GYR);
        temp = next;
        lid = grown;
        t += dt;
        k += 1;
        out.times.push(t);
        out.heat_flux.push(surface);
        out.temperature.push(temp);
        out.lid.push(lid);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::presets::preset;

    fn history(name: &str) -> (MantleHistory, f64) {
        let (g, _) = preset(name).unwrap().to_graph();
        (evolve(&MantleInputs::from_graph(&g)), g.f64("body.age"))
    }

    #[test]
    fn venus_s_lid_breaks_again_and_again_and_the_small_bodies_never() {
        // F121: Venus's mantle heats under its lid until the lid founders, again and again,
        // the last time within the half gigayear its craters record; the smaller bodies'
        // mantles cool under lids that never break, their present heat flows as measured.
        let (venus, age) = history("Venus");
        assert!(venus.overturns.len() >= 3, "{:?}", venus.overturns);
        let last = age - venus.last_overturn().unwrap();
        assert!(last > 0.3 && last < 0.7, "Venus last overturned {last} Gyr ago");
        for (name, lo, hi) in [("Mars", 0.015, 0.045), ("Moon", 0.009, 0.021), ("Mercury", 0.008, 0.03)] {
            let (h, _) = history(name);
            assert!(h.overturns.is_empty(), "{name}: {:?}", h.overturns);
            let now = *h.heat_flux.last().unwrap();
            assert!(now > lo && now < hi, "{name}: {now} W/m²");
        }
    }
}
