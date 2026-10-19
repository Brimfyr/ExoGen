// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! A rocky mantle that stores its heat under a lid that can break (F121).
//!
//! The mantle below the lid has one temperature, which its radioactive and tidal heat raise
//! and the heat it convects into the lid's base lowers; its heat capacity is the rock's and
//! the core's, less the share of the lid's that moves with it. The convected flux follows
//! the boundary-layer scaling of stagnant-lid convection, the viscosity Arrhenian in the
//! mantle's temperature. The lid conducts from the mantle's temperature at its base to the
//! surface's; it grows where it conducts more than the mantle convects into its base and
//! thins where it conducts less, by the heat its base takes up or gives (a Stefan
//! condition, as in parameterized histories since Schubert, Cassen and Young 1979), so that
//! a lid that has just formed grows as a cooling half space and an old one follows the
//! thickness a stagnant lid would have; the surface loses what the body's energy says it
//! does. The first lid grows from the surface as the body forms. The lid breaks when the
//! stress the convecting sublayer puts on it, along a cell as long as the mantle is deep,
//! exceeds the lid's frictional strength over its thickness (after Moresi and Solomatov
//! 1998's comparison of convective and yield stresses), but only once it has grown as thick
//! as a stagnant lid since it formed: until then the cold boundary layer is the lid itself,
//! growing, with no sublayer under it. A lid that breaks founders into the mantle, which
//! takes up its cold, and a new lid grows until it can break again (Turcotte 1993: the
//! interior reheats under a thickening lid until the lid fails). Whether a body's lid ever
//! breaks, and how often, comes out of its size, its heat, its surface temperature and its
//! lid's friction.
//!
//! The crust holds a share of the heat-producing elements, which the melts that made it took
//! from the mantle, and its heat leaves through the surface without heating the mantle
//! (F123). The mantle melts only while it stands above the rock's solidus at the base of its
//! stiff lid, where the mantle rising under the lid stops, so that as its lid thickens and it
//! cools its melting ends; a body's plumes rise hotter than the mantle about them and go on
//! melting where it no longer does.

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
    /// The crust's thickness, m, whose volume against the mantle's sets its share of the
    /// heat-producing elements (F123).
    pub crust_thickness: f64,
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
            crust_thickness: graph.f64("tectonics.initial_crust_thickness"),
        }
    }

    /// The share of the heat-producing elements the crust holds (F123): its volume's share
    /// of the rock, phi, weighted by the enrichment of incompatible elements in the melts
    /// that made it, Lambda phi / (Lambda phi + 1 - phi).
    pub fn crust_share(&self) -> f64 {
        let r = self.radius;
        let core = self.mass * self.core_fraction;
        let core_radius = (3.0 * core / (4.0 * std::f64::consts::PI * DENSITY_IRON)).cbrt();
        let rock = (r.powi(3) - core_radius.powi(3)).max(1.0);
        let crust = r.powi(3) - (r - self.crust_thickness.clamp(0.0, r - core_radius)).powi(3);
        let phi = (crust / rock).clamp(0.0, 1.0);
        CRUST_ENRICHMENT * phi / (CRUST_ENRICHMENT * phi + 1.0 - phi)
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
    /// The crust's share of the heat-producing elements.
    pub crust_share: f64,
    /// The body's surface gravity, m/s², and temperature, K, for the stiff lid's base.
    pub gravity: f64,
    pub surface_temperature: f64,
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

    /// Whether the mantle melts at a time, 0 to 1 (F123): how far the mantle rising under the
    /// lid, or on a body with plumes its plumes, stand above the solidus at the stiff lid's
    /// base, over one rheological temperature scale. The stiff lid ends where the conductive
    /// profile reaches LID_BASE_SCALES rheological scales below the mantle's temperature,
    /// where the mantle stops flowing (Davaille and Jaupart 1993).
    pub fn melt_gate_at(&self, t: f64, plumes: bool) -> f64 {
        let n = self.times.len();
        if n == 0 {
            return 0.0;
        }
        let k = self.times.partition_point(|&x| x < t).min(n - 1);
        let temp = self.temperature[k];
        let ts = self.surface_temperature;
        let rheological = GAS_CONSTANT * temp * temp / MANTLE_ACTIVATION;
        let stiff = self.lid[k] * (temp - LID_BASE_SCALES * rheological - ts).max(0.0)
            / (temp - ts).max(1.0);
        let solidus = solidus(DENSITY_ROCK * self.gravity * stiff);
        let excess = if plumes {
            temp + PLUME_EXCESS_K - solidus
        } else {
            temp - solidus
        };
        (excess / rheological).clamp(0.0, 1.0)
    }

    /// The surface heat flux averaged over the body's life, W/m².
    pub fn mean_heat_flux(&self) -> f64 {
        if self.heat_flux.is_empty() {
            return 0.0;
        }
        self.heat_flux.iter().sum::<f64>() / self.heat_flux.len() as f64
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

/// The dry peridotite solidus at a pressure, K (Hirschmann 2000's fit, P in GPa).
pub fn solidus(pressure: f64) -> f64 {
    let p = pressure / 1e9;
    273.15 + SOLIDUS_C0 + SOLIDUS_C1 * p + SOLIDUS_C2 * p * p
}

/// Whether the lid breaks: where the sublayer's buoyancy stress rho g alpha dT_rh delta_rh,
/// delta_rh = k dT_rh / q, over a cell as long as the mantle is deep, exceeds the lid's
/// friction mu rho g d² / 2 over its thickness d. At a stagnant lid's thickness,
/// d = k dT / q, this is a flux above mu k dT² / (2 alpha dT_rh² D).
fn breaks(t: f64, depth: f64, lid: f64, q: f64, friction: f64) -> bool {
    let rheological = GAS_CONSTANT * t * t / MANTLE_ACTIVATION;
    2.0 * MANTLE_EXPANSIVITY * rheological * rheological * K_ROCK * depth.max(1.0)
        > friction * lid * lid * q
}

/// The mantle's history from formation to the present (F121).
pub fn evolve(inp: &MantleInputs) -> MantleHistory {
    evolve_with(inp, LID_FRICTION)
}

/// The mantle's history under a lid of this effective friction.
fn evolve_with(inp: &MantleInputs, friction: f64) -> MantleHistory {
    let r = inp.radius;
    let area = 4.0 * std::f64::consts::PI * r * r;
    let mantle = inp.mass * inp.rock_fraction;
    let core = inp.mass * inp.core_fraction;
    let core_radius = (3.0 * core / (4.0 * std::f64::consts::PI * DENSITY_IRON)).cbrt();
    let depth = (r - core_radius).max(1.0);
    let capacity = (mantle * MANTLE_HEAT_CAPACITY + core * CORE_HEAT_CAPACITY) / area;
    let ts = inp.surface_temperature;
    let rho_c = DENSITY_ROCK * MANTLE_HEAT_CAPACITY;
    // The lid's diffusivity, its conductivity over the heat capacity its deficit counts, so
    // that its growth and its heat agree.
    let kappa = K_ROCK / rho_c;
    // The lid's volume per unit area, m, and the heat it holds below the mantle's, J/m², its
    // temperature rising evenly from the surface's to the mantle's.
    let shell = |d: f64| (r.powi(3) - (r - d.min(depth)).powi(3)) / (3.0 * r * r);
    let deficit = |t: f64, d: f64| rho_c * shell(d) * (t - ts).max(0.0) / 2.0;
    // The body's heat per unit area, J/m²: all of it at the mantle's temperature, less the
    // lid's deficit.
    let energy = |t: f64, d: f64| capacity * t - deficit(t, d);
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
    let crust_share = inp.crust_share();
    let mut out = MantleHistory {
        times: Vec::with_capacity(n),
        heat_flux: Vec::with_capacity(n),
        temperature: Vec::with_capacity(n),
        lid: Vec::with_capacity(n),
        overturns: Vec::new(),
        crust_share,
        gravity: inp.gravity,
        surface_temperature: ts,
    };
    let mut temp = MANTLE_START_K;
    // The first lid grows from the surface as the body's magma ocean freezes.
    let mut lid = 0.0_f64;
    let mut formed = 0.0;
    let mut t = 0.0;
    let mut k = 0;
    let step = dt * GYR;
    while k < n {
        // The crust's share of the radioactive heat leaves through the surface; the rest and
        // the tides heat the mantle.
        let radioactive = radiogenic_per_kg(t + 0.5 * dt) * mantle / area;
        let heat = (1.0 - crust_share) * radioactive + tides(t + 0.5 * dt);
        let q = stagnant_flux(temp, inp.gravity);
        let steady = (K_ROCK * (temp - ts) / q).min(0.99 * depth);
        // Grown as thick as a stagnant lid: by the time a half space would have, since the
        // lid under a cooling mantle only nears that thickness.
        let mature = std::f64::consts::PI * kappa * (t - formed) * GYR >= steady * steady;
        if mature && breaks(temp, depth, lid, q, friction) {
            // The lid founders, and the whole mantle takes up its cold.
            temp -= deficit(temp, lid) / capacity;
            lid = 0.0;
            formed = t;
            out.overturns.push(t);
            continue;
        }
        // The mantle keeps its heat less what it convects into the lid's base; the lid's
        // deficit moves with its base's temperature by half its volume.
        let convecting = (capacity - 0.5 * rho_c * shell(lid)).max(0.05 * capacity);
        let next = temp + (heat - q) * step / convecting;
        // The lid's base takes up the heat it conducts beyond what reaches it, so its square
        // grows at 4 kappa (1 - d q / (k dT)) per its base's share of the surface: a half
        // space's from nothing, none at a stagnant lid's thickness, shrinking past it.
        let dtemp = (temp - ts).max(1.0);
        let base_share = ((r - lid.min(depth)) / r).powi(2).max(1e-3);
        let square =
            lid * lid + step * 4.0 * kappa * (1.0 - lid * q / (K_ROCK * dtemp)) / base_share;
        let grown = square.max(0.0).sqrt().min(0.99 * depth);
        // The surface loses the heat made less the heat the body keeps, and the crust's.
        let surface =
            heat - (energy(next, grown) - energy(temp, lid)) / step + crust_share * radioactive;
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

    /// Diagnostic: each rocky preset's mantle history beside the thermal history the
    /// pipeline reads, mW/m² and K, at a few times; run with `--ignored --nocapture`.
    #[test]
    #[ignore]
    fn print_mantle_histories() {
        let mut csv = String::from("body,time_gyr,heat_flux,temperature,lid\n");
        for name in ["Mercury", "Venus", "Earth", "Moon", "Mars", "Io"] {
            let (g, _) = preset(name).unwrap().to_graph();
            let h = evolve(&MantleInputs::from_graph(&g));
            let th = crate::thermal::ThermalHistory::from_graph(&g);
            let age = g.f64("body.age");
            let r = g.f64("body.radius");
            let mean = h.heat_flux.iter().sum::<f64>() / h.heat_flux.len() as f64;
            println!(
                "{name}: rock {:.3} core {:.3}, overturns {}, regime {}, mean loss {:.2} TW",
                g.f64("body.rock_mass_fraction"),
                g.f64("body.core_mass_fraction"),
                h.overturns.len(),
                g.str("interior.tectonic_regime"),
                mean * 4.0 * std::f64::consts::PI * r * r / 1e12
            );
            for k in 0..h.times.len() {
                csv += &format!(
                    "{name},{:.5},{:.6},{:.1},{:.0}\n",
                    h.times[k], h.heat_flux[k], h.temperature[k], h.lid[k]
                );
            }
            for t in [0.05, 0.1, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, age] {
                let k = h.times.partition_point(|&x| x < t).min(h.times.len() - 1);
                println!(
                    "  {t:>5.2} Gyr: mantle {:>6.1} mW/m² {:>6.0} K lid {:>5.0} km | pipeline {:>6.1} mW/m² {:>6.0} K",
                    1e3 * h.heat_flux_at(t),
                    h.temperature_at(t),
                    h.lid[k] / 1e3,
                    1e3 * th.at(t).heat_flux,
                    th.at(t).mantle_temperature
                );
            }
        }
        if let Ok(path) = std::env::var("EXOGEN_MANTLE_CSV") {
            std::fs::write(path, csv).unwrap();
        }
    }

    /// Diagnostic: Venus's overturns under a range of frictions; run with `--ignored
    /// --nocapture`.
    #[test]
    #[ignore]
    fn print_friction_scan() {
        let (g, _) = preset("Venus").unwrap().to_graph();
        let inp = MantleInputs::from_graph(&g);
        let age = g.f64("body.age");
        // EXOGEN_MU_START and EXOGEN_MU_STEP narrow the scan.
        let start = std::env::var("EXOGEN_MU_START")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(0.004);
        let stride = std::env::var("EXOGEN_MU_STEP")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(0.0005);
        for i in 0..=40 {
            let mu: f64 = start + stride * i as f64;
            let h = evolve_with(&inp, mu);
            let ago: Vec<String> = h
                .overturns
                .iter()
                .rev()
                .take(4)
                .map(|t| format!("{:.2}", age - t))
                .collect();
            println!(
                "mu {mu:.4}: {:>3} overturns, last Ga {:?}, now {:.1} mW/m²",
                h.overturns.len(),
                ago,
                1e3 * h.heat_flux.last().unwrap()
            );
        }
    }

    #[test]
    fn venus_s_lid_breaks_again_and_again_and_the_small_bodies_never() {
        // F121: Venus's mantle heats under its lid until the lid founders, again and again,
        // the last time within the half gigayear its craters record; the smaller bodies'
        // mantles cool under lids that never break, their present heat flows as measured.
        let (venus, age) = history("Venus");
        assert!(venus.overturns.len() >= 3, "{:?}", venus.overturns);
        let last = age - venus.last_overturn().unwrap();
        assert!(
            last > 0.3 && last < 0.7,
            "Venus last overturned {last} Gyr ago"
        );
        for (name, lo, hi) in [
            ("Mars", 0.015, 0.045),
            ("Moon", 0.009, 0.021),
            ("Mercury", 0.008, 0.03),
        ] {
            let (h, _) = history(name);
            assert!(h.overturns.is_empty(), "{name}: {:?}", h.overturns);
            let now = *h.heat_flux.last().unwrap();
            assert!(now > lo && now < hi, "{name}: {now} W/m²");
            // F122: the surface loses heat from the start, a lid thinning under a heating
            // mantle only as fast as its base takes up the heat, and lost most when young.
            let least = h.heat_flux.iter().copied().fold(f64::INFINITY, f64::min);
            assert!(least > 0.0, "{name}: {least} W/m²");
            assert!(
                h.heat_flux[0] > 2.0 * now,
                "{name}: young {} now {now}",
                h.heat_flux[0]
            );
        }
    }
}
