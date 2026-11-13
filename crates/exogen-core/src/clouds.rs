//! The clouds (M66).
//!
//! The working liquid's clouds through the cycle's year at its places, from what its air holds
//! and how it moves, as a coarse climate model's diagnostic scheme reads them (Slingo 1987, as
//! ECMWF ran it): cloud from the convection its rain shows, a tower over a share of the place
//! and an anvil spread from it once the convection is strong; cloud in the layers of air that
//! rise, by how near its column stands to saturation; and low cloud under the inversion a
//! cold surface makes beneath a warmer free air, by its strength (Wood and Bretherton 2006),
//! where the air beneath it is moist. The free air over a place stands at the warmest of the
//! latitudes its overturning cell spans, as weak temperature gradients hold it (Sobel et al.
//! 2001), and over its own latitude beyond the cell.
//!
//! Each kind reflects the star's light by its optical depth, its water's path over its drops'
//! size (Stephens 1978), as a conservative layer of its drops' asymmetry, and sends its
//! infrared to space from its top and down to the ground from its base at its emissivity: the
//! kinds overlap at random, the higher hiding the lower from space and the lower the higher
//! from the ground.
//!
//! They are a diagnostic for now: the climate keeps the clouds' stated law (F329) until the
//! planet's albedo and greenhouse come from the clouds themselves (M67; the user, 2026-11-13).

use crate::constants::*;
use crate::params::ParamGraph;
use crate::params::derive::physics::{
    mean_heat_capacity, mean_molecular_mass, vapour_pressure, volatile,
};
use crate::radiation::{Column, Heights, Outgoing, Sunlight};

/// The clouds' kinds at a place and a step, each by the share of the place it covers.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Kinds {
    /// The convection's towers, and the anvils they spread at the tropopause.
    pub tower: f64,
    pub anvil: f64,
    /// The layers of rising air near saturation, their tops in the middle air.
    pub layer: f64,
    /// The low cloud under an inversion.
    pub low: f64,
}

impl Kinds {
    /// The share of the place any of them covers, overlapping at random.
    pub fn cover(&self) -> f64 {
        1.0 - (1.0 - self.tower) * (1.0 - self.anvil) * (1.0 - self.layer) * (1.0 - self.low)
    }
}

/// The clouds' covers at a place (M66): its convection's rain and its large-scale rain, mm a
/// day, its column's saturation, its inversion's estimated strength, K, and the relative
/// humidity of the air under the inversion.
///
/// The convection's cloud C = a + b ln P, at most CLOUD_CONVECTIVE's most (Slingo 1987); a
/// share CLOUD_TOWER of it the towers, and the anvils 2 (C - 0.3) of the place where it passes
/// 0.3 (ECMWF 1991, which began them at 0.4, a step this begins nowhere). The large-scale rain's
/// systems cover a share of the place by the same law, their cloud the layers' (M66, the law
/// read as the share of a place its raining systems cover; Slingo's own layers came from each
/// model layer's humidity, which a column's mean saturation, some 0.7 to 0.8 under clouds that
/// stand saturated in their own layers, cannot show); and the layers' own cloud ((r - r_c) /
/// (1 - r_c))², r the column's saturation outside the towers and r_c CLOUD_LAYER (ECMWF 1991,
/// eq. 2.55, without its rise, which a year's mean wind does not hold: the storms lift the air,
/// not the mean flow). The low cloud 0.14 + 0.060 EIS (Wood and Bretherton 2006), none where the
/// air under the inversion stands under half saturated, all of it from seven tenths (ECMWF
/// 1991's dry inversions over deserts and the winter pole).
pub fn kinds(convective: f64, large: f64, saturation: f64, eis: f64, humidity: f64) -> Kinds {
    let (a, b, most) = CLOUD_CONVECTIVE;
    let law = |rain: f64| {
        if rain > 0.0 {
            (a + b * rain.ln()).clamp(0.0, most)
        } else {
            0.0
        }
    };
    let deep = law(convective);
    let tower = CLOUD_TOWER * deep;
    let (from, slope) = CLOUD_ANVIL;
    let anvil = (slope * (deep - from)).clamp(0.0, 1.0);
    let outside = (saturation - tower).clamp(0.0, 1.0);
    let critical = CLOUD_LAYER;
    let humid = ((outside - critical) / (1.0 - critical)).max(0.0);
    let layer = 1.0 - (1.0 - law(large)) * (1.0 - (humid * humid).min(1.0));
    let (base, per_kelvin) = CLOUD_LOW;
    let (dry, moist) = CLOUD_LOW_HUMIDITY;
    let wet = ((humidity - dry) / (moist - dry)).clamp(0.0, 1.0);
    let low = (base + per_kelvin * eis).clamp(0.0, 1.0) * wet;
    Kinds {
        tower,
        anvil,
        layer,
        low,
    }
}

/// What the clouds read of the body and its air (M66): its clear air's streams over its
/// ground's temperature and over grounds at each height (F293, M55), without the clouds'
/// stated law; how far its ice brightens toward ICE_ALBEDO under its star's light (F250); its
/// lapse rate, K/m, its scale height, m, its tropopause's temperature, its
/// skin's, K; its surface pressure, Pa, gravity, m/s², the air's gas constant and heat
/// capacity, J/kg/K; and its water.
#[derive(Clone, Debug)]
pub struct Sky {
    pub clear: Outgoing,
    pub brightening: f64,
    pub heights: Heights,
    pub lapse: f64,
    pub scale_height: f64,
    pub tropopause: f64,
    pub pressure: f64,
    pub gravity: f64,
    pub gas_constant: f64,
    pub heat_capacity: f64,
    pub water: &'static Volatile,
}

impl Sky {
    /// A body's sky for its water's clouds under its energy balance `balance`: none without air
    /// or with another working liquid (its other condensables' clouds wait on their own
    /// milestone, the user, 2026-11-12).
    pub fn from_graph(graph: &ParamGraph, balance: &crate::climate::EnergyBalance) -> Option<Self> {
        if graph.str("volatiles.working_liquid") != "water" {
            return None;
        }
        let pressure = graph.f64("atmosphere.surface_pressure");
        if pressure <= 0.0 {
            return None;
        }
        let empty = std::collections::BTreeMap::new();
        let composition = graph
            .get("atmosphere.composition")
            .as_map()
            .unwrap_or(&empty);
        let gravity = graph.f64("body.gravity");
        let t_eq = graph.f64("atmosphere.equilibrium_temperature");
        let skin = t_eq * SKIN_SHARE.powf(0.25);
        let water = balance.water;
        let column = Column::new(composition, pressure, gravity, water, skin)
            .without_clouds()
            .with_sunlight(Sunlight::of(
                composition,
                (pressure, gravity),
                graph.f64("orbit.insolation") / 4.0,
                graph.f64("atmosphere.bond_albedo"),
                graph.f64("star.xuv_relative"),
            ));
        let (low, high) = ((t_eq / 3.0).max(5.0), 4.0 * t_eq + 100.0);
        Some(Self {
            clear: Outgoing::of(&column, low, high),
            brightening: balance.brightening,
            heights: Heights::of(&column, low, high),
            lapse: graph.f64("climate.lapse_rate"),
            scale_height: graph.f64("atmosphere.scale_height"),
            tropopause: skin,
            pressure,
            gravity,
            gas_constant: GAS_CONSTANT / mean_molecular_mass(composition),
            heat_capacity: mean_heat_capacity(composition),
            water: volatile("H2O")?,
        })
    }

    /// The air's temperature at `ratio` of the surface pressure over air at `air` K at the
    /// ground: down its lapse rate to its tropopause, and no colder.
    fn at(&self, air: f64, ratio: f64) -> f64 {
        let z = -self.scale_height * ratio.max(1e-9).ln();
        (air - self.lapse * z).max(self.tropopause.min(air))
    }

    /// The share of the surface pressure the tropopause stands at over air at `air` K.
    fn tropopause_ratio(&self, air: f64) -> f64 {
        let z = ((air - self.tropopause) / self.lapse.max(1e-9)).max(0.0);
        (-z / self.scale_height.max(1.0)).exp()
    }

    /// The saturation mixing ratio of the water over `t` K at `p` Pa.
    fn mixing(&self, t: f64, p: f64) -> f64 {
        let e = vapour_pressure(self.water, t).min(0.5 * p);
        self.gas_constant / self.water.r_specific * e / (p - e)
    }

    /// The estimated inversion strength over a surface at `surface` K under free air that
    /// stands at `free` K at its own ground (M66): Wood and Bretherton (2006)'s EIS, LTS -
    /// Gamma_m(850) (z700 - LCL), its LTS the potential temperature of the free air at
    /// CLOUD_LEVELS's 700 hPa share of the surface pressure over the surface's, the free air
    /// down its lapse rate; its moist adiabat's potential temperature gradient Gamma_m at the
    /// mean of the two temperatures at the 850 hPa share; its lifting condensation level 125 m
    /// a kelvin of the surface air's dew point depression at CLOUD_SURFACE_HUMIDITY (Lawrence
    /// 2005).
    pub fn eis(&self, surface: f64, free: f64) -> f64 {
        let (x700, x850) = CLOUD_LEVELS;
        let (g, r, cp) = (self.gravity, self.gas_constant, self.heat_capacity);
        let z700 = r * surface / g * (1.0 / x700).ln();
        let t700 = free - self.lapse * r * free / g * (1.0 / x700).ln();
        let theta700 = t700 * (1.0 / x700).powf(r / cp);
        let lts = theta700 - surface;
        let t850 = 0.5 * (surface + t700);
        let p850 = x850 * self.pressure;
        let q = self.mixing(t850, p850);
        let (l, rv) = (self.water.l_vap, self.water.r_specific);
        let gamma = g / cp
            * (1.0 - (1.0 + l * q / (r * t850)) / (1.0 + l * l * q / (cp * rv * t850 * t850)));
        let dew = 1.0 / (1.0 / surface - rv / l * CLOUD_SURFACE_HUMIDITY.ln());
        let lcl = CLOUD_LCL_PER_KELVIN * (surface - dew).max(0.0);
        lts - gamma * (z700 - lcl)
    }

    /// The share of the surface pressure a place's lifting condensation level stands at, over
    /// air at `air` K at CLOUD_SURFACE_HUMIDITY.
    fn base_ratio(&self, air: f64) -> f64 {
        let (l, rv) = (self.water.l_vap, self.water.r_specific);
        let dew = 1.0 / (1.0 / air - rv / l * CLOUD_SURFACE_HUMIDITY.ln());
        let lcl = CLOUD_LCL_PER_KELVIN * (air - dew).max(0.0);
        (-lcl / self.scale_height.max(1.0)).exp()
    }

    /// What clouds of `kinds` change at a place (M66), the air over it at `air` K, its ground
    /// at `ground` K of albedo `albedo`, `sea` the share of it the sea covers: the planet's
    /// albedo there, what it lets out to space, W/m², and what its air sends down to its ground,
    /// W/m², each with the clouds less without.
    pub fn effects(&self, kinds: &Kinds, air: f64, ground: f64, albedo: f64, sea: f64) -> Effects {
        let sigma = |t: f64| SIGMA_SB * t.powi(4);
        // The kinds' water and optics.
        let (liquid, ice) = CLOUD_RADIUS;
        let drop = liquid.1 + sea.clamp(0.0, 1.0) * (liquid.0 - liquid.1);
        let top = self.tropopause_ratio(air);
        let base = self.base_ratio(air);
        let tower_path = CLOUD_TOWER_WATER * self.pressure * (base - top).max(0.0) / self.gravity;
        let tau_liquid = |path: f64| 3.0 * path / (2.0 * self.water.liquid_density * drop);
        let (g_liquid, g_ice) = (
            CLOUD_ASYMMETRY.0 + CLOUD_ASYMMETRY.1 * drop * 1e6,
            CLOUD_ASYMMETRY.2 + CLOUD_ASYMMETRY.3 * ice * 1e6,
        );
        let ice_path =
            CLOUD_ANVIL_DEPTH / (CLOUD_ICE_DEPTH.0 + CLOUD_ICE_DEPTH.1 / (ice * 1e6)) * 1e-3;
        let ice_tau = CLOUD_ANVIL_DEPTH;
        let reflect = |tau: f64, g: f64| {
            let x = 0.75 * (1.0 - g) * tau;
            x / (1.0 + x)
        };
        let r_tower = reflect(tau_liquid(tower_path), g_liquid);
        let r_anvil = reflect(ice_tau, g_ice);
        let r_layer = reflect(tau_liquid(CLOUD_LIQUID_PATH), g_liquid);
        let r_low = r_layer;
        let e_liquid = |path: f64| 1.0 - (-CLOUD_LIQUID_EMISSION * path * 1e3).exp();
        let e_ice = 1.0
            - (-DIFFUSIVITY
                * (CLOUD_ICE_EMISSION.0 + CLOUD_ICE_EMISSION.1 / (ice * 1e6))
                * ice_path
                * 1e3)
                .exp();
        let (e_tower, e_layer) = (e_liquid(tower_path), e_liquid(CLOUD_LIQUID_PATH));
        // The star's light: the clouds over the ground as one conservative layer reflecting
        // what any of them reflects at random overlap, the light reflected between it and the
        // ground.
        let k = kinds;
        let reflected = 1.0
            - (1.0 - k.tower * r_tower)
                * (1.0 - k.anvil * r_anvil)
                * (1.0 - k.layer * r_layer)
                * (1.0 - k.low * r_low);
        let a = albedo.clamp(0.0, 1.0);
        let planet = reflected + (1.0 - reflected).powi(2) * a / (1.0 - reflected * a).max(1e-9);
        // To space: from the tops, the higher hiding the lower.
        let out_from = |ratio: f64| {
            let t = self.at(air, ratio);
            sigma(t) * self.heights.share(t, ratio)
        };
        let clear = self.clear.air_up(air) + self.clear.window(air) * sigma(ground);
        let (mid_top, low_top) = CLOUD_TOPS;
        let high = 1.0 - (1.0 - k.tower * e_tower) * (1.0 - k.anvil * e_ice);
        let (mid, low) = (k.layer * e_layer, k.low * e_layer);
        let olr = high * out_from(top)
            + (1.0 - high)
                * (mid * out_from(mid_top)
                    + (1.0 - mid) * (low * out_from(low_top) + (1.0 - low) * clear));
        // To the ground: from the bases through the air's window, the lower hiding the higher.
        let window = self.clear.window(air);
        let (mid_base, anvil_base) = CLOUD_BASES;
        let near = 1.0 - (1.0 - k.low * e_layer) * (1.0 - k.tower * e_tower);
        let down = window
            * (near * sigma(self.at(air, base))
                + (1.0 - near)
                    * (mid * sigma(self.at(air, mid_base))
                        + (1.0 - mid) * k.anvil * e_ice * sigma(self.at(air, anvil_base))));
        Effects {
            albedo: planet - a,
            olr: olr - clear,
            down,
        }
    }
}

/// What clouds change at a place (M66): the planet's albedo there, what it lets out to space,
/// W/m², and what its air sends down to its ground, W/m², each with the clouds less without.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Effects {
    pub albedo: f64,
    pub olr: f64,
    pub down: f64,
}

#[cfg(test)]
mod tests {
    use super::*;

    /// M66: every kind covers between none and all of its place, whatever its causes, and a
    /// place with no rain, a column far from saturation and dry air under its inversion has
    /// none.
    #[test]
    fn the_kinds_cover_no_more_than_their_place() {
        for convective in [0.0, 0.1, 1.0, 5.0, 30.0, 300.0] {
            for large in [0.0, 0.5, 5.0, 100.0] {
                for saturation in [0.0, 0.5, 0.8, 1.0, 1.3] {
                    for eis in [-15.0, 0.0, 8.0, 30.0] {
                        for humidity in [0.0, 0.6, 1.0, 1.4] {
                            let k = kinds(convective, large, saturation, eis, humidity);
                            for v in [k.tower, k.anvil, k.layer, k.low, k.cover()] {
                                assert!(
                                    (0.0..=1.0).contains(&v),
                                    "{v} at {convective} {large} {saturation} {eis} {humidity}"
                                );
                            }
                        }
                    }
                }
            }
        }
        assert_eq!(kinds(0.0, 0.0, 0.5, 10.0, 0.3).cover(), 0.0);
    }

    /// M66: each kind grows with its cause: the towers and their anvils with the convection's
    /// rain, the anvils none under weak convection; the layers with the large-scale rain and
    /// the column's saturation; the low cloud with its inversion's strength and the air's
    /// humidity under it, none in dry air whatever the inversion.
    #[test]
    fn each_kind_grows_with_its_cause() {
        let rains = [0.0, 0.5, 1.0, 2.0, 5.0, 10.0, 30.0, 100.0];
        for w in rains.windows(2) {
            let (a, b) = (
                kinds(w[0], 0.0, 0.5, 0.0, 0.0),
                kinds(w[1], 0.0, 0.5, 0.0, 0.0),
            );
            assert!(b.tower >= a.tower && b.anvil >= a.anvil, "{w:?}");
            let (a, b) = (
                kinds(0.0, w[0], 0.5, 0.0, 0.0),
                kinds(0.0, w[1], 0.5, 0.0, 0.0),
            );
            assert!(b.layer >= a.layer, "{w:?}");
        }
        assert!(kinds(100.0, 0.0, 0.5, 0.0, 0.0).tower > kinds(2.0, 0.0, 0.5, 0.0, 0.0).tower);
        assert_eq!(kinds(1.0, 0.0, 0.5, 0.0, 0.0).anvil, 0.0);
        assert!(kinds(30.0, 0.0, 0.5, 0.0, 0.0).anvil > 0.5);
        let saturations = [0.0, 0.5, 0.7, 0.8, 0.9, 1.0];
        for w in saturations.windows(2) {
            let (a, b) = (
                kinds(0.0, 0.0, w[0], 0.0, 0.0),
                kinds(0.0, 0.0, w[1], 0.0, 0.0),
            );
            assert!(b.layer >= a.layer, "{w:?}");
        }
        assert!(kinds(0.0, 0.0, 1.0, 0.0, 0.0).layer > 0.9);
        for w in [-5.0, 0.0, 4.0, 8.0, 12.0].windows(2) {
            let (a, b) = (
                kinds(0.0, 0.0, 0.5, w[0], 0.9),
                kinds(0.0, 0.0, 0.5, w[1], 0.9),
            );
            assert!(b.low > a.low, "{w:?}");
        }
        for w in [0.3, 0.5, 0.6, 0.7, 0.9].windows(2) {
            let (a, b) = (
                kinds(0.0, 0.0, 0.5, 6.0, w[0]),
                kinds(0.0, 0.0, 0.5, 6.0, w[1]),
            );
            assert!(b.low >= a.low, "{w:?}");
        }
        assert_eq!(kinds(0.0, 0.0, 0.5, 20.0, 0.4).low, 0.0);
    }

    /// M66: the Earth-like's clouds reflect its light, keep its heat from space and send it
    /// down: none of them change nothing; a low deck reflects much and keeps little, a high
    /// anvil keeps far more from space and sends less down; a deck over ice adds less to its
    /// albedo than over the sea; an airless body has no sky for clouds.
    #[test]
    #[cfg_attr(not(feature = "whole-body"), ignore = "whole body")]
    fn the_earth_like_s_clouds_reflect_its_light_and_keep_its_heat() {
        let graph = crate::presets::preset("Earth").unwrap().to_graph().0;
        let balance = crate::climate::EnergyBalance::from_graph(&graph);
        let sky = Sky::from_graph(&graph, &balance).expect("a sky");
        let none = sky.effects(&Kinds::default(), 288.0, 288.0, 0.1, 1.0);
        assert!(none.albedo.abs() < 1e-12 && none.olr.abs() < 1e-9 && none.down.abs() < 1e-9);
        let low = Kinds {
            low: 1.0,
            ..Kinds::default()
        };
        let anvil = Kinds {
            anvil: 1.0,
            ..Kinds::default()
        };
        let (deck, high) = (
            sky.effects(&low, 290.0, 292.0, SEA_ALBEDO, 1.0),
            sky.effects(&anvil, 300.0, 301.0, SEA_ALBEDO, 1.0),
        );
        println!("Earth-like: a low deck {deck:?}, an anvil {high:?}");
        assert!(
            deck.albedo > 0.2 && deck.olr < 0.0 && deck.down > 0.0,
            "{deck:?}"
        );
        assert!(high.olr < deck.olr - 50.0, "{high:?} {deck:?}");
        assert!(deck.down > high.down, "{high:?} {deck:?}");
        let over_ice = sky.effects(&low, 260.0, 255.0, ICE_ALBEDO, 0.0);
        assert!(over_ice.albedo < deck.albedo - 0.1, "{over_ice:?}");
        let moon = crate::presets::preset("Moon").unwrap().to_graph().0;
        let airless = crate::climate::EnergyBalance::from_graph(&moon);
        assert!(Sky::from_graph(&moon, &airless).is_none());
    }
}
