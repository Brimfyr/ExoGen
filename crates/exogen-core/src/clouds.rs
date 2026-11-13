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
    /// The convection's towers, and the anvils they spread at the tropopause, `medium` of
    /// them medium-thick and the rest thin cirrus (M67).
    pub tower: f64,
    pub anvil: f64,
    pub medium: f64,
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
/// The convection's high cloud a + b ln P by its climatological rain, its towers c P of the
/// place, no more than all of it, and its anvils the rest of it, overlapping the towers at
/// random (CLOUD_HIGH, CLOUD_THICK; M67). The large-scale rain's systems cover a share of the
/// place by Slingo (1987)'s law for convective cloud, C = a + b ln P at most CLOUD_CONVECTIVE's
/// most, their cloud the layers' (M66, the law
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
    let (h0, h1) = CLOUD_HIGH;
    let high = if convective > 0.0 {
        (h0 + h1 * convective.ln()).clamp(0.0, 1.0)
    } else {
        0.0
    };
    let tower = (CLOUD_THICK * convective).clamp(0.0, high);
    let anvil = if tower < 1.0 {
        (1.0 - (1.0 - high) / (1.0 - tower)).max(0.0)
    } else {
        0.0
    };
    let outside = (saturation - tower).clamp(0.0, 1.0);
    let critical = CLOUD_LAYER;
    let humid = ((outside - critical) / (1.0 - critical)).max(0.0);
    let layer = 1.0 - (1.0 - law(large)) * (1.0 - (humid * humid).min(1.0));
    let (base, per_kelvin) = CLOUD_LOW;
    let (dry, moist) = CLOUD_LOW_HUMIDITY;
    let wet = ((humidity - dry) / (moist - dry)).clamp(0.0, 1.0);
    let low = (base + per_kelvin * eis).clamp(0.0, 1.0) * wet;
    let (m0, m1, _) = CLOUD_ANVIL_MEDIUM;
    let medium = if anvil > 0.0 {
        (m0 + m1 * convective).clamp(0.0, 1.0)
    } else {
        0.0
    };
    Kinds {
        tower,
        anvil,
        medium,
        layer,
        low,
    }
}

/// The clouds over a place at a step as the year reads them (M67), each kind in the order
/// tower, anvil, layer, low cloud: the share of the place it covers; its optical depth in the
/// star's light scaled by its drops' asymmetry, (1 - g) tau, which sets its reflection; its
/// optical depth delta-scaled, (1 - g²) tau (Joseph et al. 1976), which sets how its
/// reflection rises as the sun sinks; and its emissivity in the infrared. And the share of
/// the light the clouds absorb (M67), and the share of the air's absorbers of the light above
/// the clouds that reflect it, the kinds' tops by what each reflects.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Deck {
    pub cover: [f32; 4],
    pub scaled: [f32; 4],
    pub depth: [f32; 4],
    pub emissivity: [f32; 4],
    pub absorbs: f32,
    pub above: f32,
}

impl Deck {
    /// The share of the star's light the clouds reflect with the sun at height cosine `mu`
    /// (M67): each kind a conservative layer by the delta-Eddington approximation, (s + (2/3 -
    /// mu)(1 - exp(-tau' / mu))) / (4/3 + s), s its scaled depth and tau' its delta-scaled
    /// depth (Joseph et al. 1976; at mu = 2/3, the diffuse light's s / (4/3 + s), M66's); the
    /// kinds overlapping at random.
    pub fn reflectance(&self, mu: f64) -> f64 {
        let mu = mu.clamp(1e-3, 1.0);
        let mut through = 1.0;
        for i in 0..4 {
            let (c, s, t) = (
                self.cover[i] as f64,
                self.scaled[i] as f64,
                self.depth[i] as f64,
            );
            if c <= 0.0 || s <= 0.0 {
                continue;
            }
            let r = ((s + (2.0 / 3.0 - mu) * (1.0 - (-t / mu).exp())) / (4.0 / 3.0 + s))
                .clamp(0.0, 1.0);
            through *= 1.0 - c * r;
        }
        1.0 - through
    }

    /// The share of the place any kind covers, overlapping at random.
    pub fn cover(&self) -> f64 {
        1.0 - self.cover.iter().map(|&c| 1.0 - c as f64).product::<f64>()
    }
}

/// What the clouds read of the body and its air (M66): its clear air's streams over its
/// ground's temperature and over grounds at each height (F293, M55), without the clouds'
/// stated law, its vapour held as its year's is (M59); how far its ice brightens toward
/// ICE_ALBEDO under its star's light (F250); its lapse rate, K/m, its scale height, m, its
/// tropopause's temperature, its skin's, K; its surface pressure, Pa, gravity, m/s², the air's
/// gas constant and heat capacity, J/kg/K; its water; and the ground under its air (M67), its
/// land's own albedo, its bedrock's mature colour's, and the share of the light its air's
/// molecules reflect (F293).
#[derive(Clone, Debug)]
pub struct Sky {
    pub clear: std::sync::Arc<Outgoing>,
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
    pub ground: f64,
    pub rayleigh: f64,
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
            .with_vapour_most(balance.vapour_most)
            .with_sunlight(Sunlight::of(
                composition,
                (pressure, gravity),
                graph.f64("orbit.insolation") / 4.0,
                graph.f64("atmosphere.bond_albedo"),
                graph.f64("star.xuv_relative"),
            ));
        let (low, high) = ((t_eq / 3.0).max(5.0), 4.0 * t_eq + 100.0);
        let library = crate::stages::materials::Library::from_graph(graph);
        let luminance = |c: [f32; 3]| (c[0] + c[1] + c[2]) as f64 / 3.0;
        let ground = library
            .member(graph.str("materials.bedrock"))
            .or_else(|| library.member("basalt"))
            .map_or(ALBEDO_BENEATH, |m| luminance(m.mature));
        let rayleigh = crate::radiation::rayleigh_reflectance(
            crate::radiation::rayleigh_depth(composition, graph.f64("atmosphere.column_mass")),
            graph.f64("star.teff"),
        );
        Some(Self {
            clear: std::sync::Arc::new(Outgoing::of(&column, low, high)),
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
            ground,
            rayleigh,
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

    /// The deck clouds of `kinds` make at a place (M66, M67), the air over it at `air` K,
    /// `sea` the share of it the sea covers: each kind's water held to what its air gives it,
    /// its share of it liquid by its level's temperature (CLOUD_PHASE), its drops' size; its
    /// optical depth by its water's path over its drops' size (Stephens 1978), its drops'
    /// asymmetry, its emissivity, and the light it absorbs (CLOUD_ABSORPTION).
    pub fn deck(&self, kinds: &Kinds, air: f64, sea: f64) -> Deck {
        let (liquid, ice) = CLOUD_RADIUS;
        let drop = liquid.1 + sea.clamp(0.0, 1.0) * (liquid.0 - liquid.1);
        let top = self.tropopause_ratio(air);
        let base = self.base_ratio(air);
        // Each kind's water, kg/m²: the towers' a share of what their ascent condenses from
        // their base to the tropopause, the decks' the adiabatic water of their depth at their
        // level (CLOUD_DEPTH), the anvils' their stated depth's.
        let mass = self.pressure * (base - top).max(0.0) / self.gravity;
        let at_base = self.mixing(self.at(air, base), base * self.pressure);
        let tower_path = CLOUD_TOWER_SHARE * at_base * mass;
        let (mid_top, low_top) = CLOUD_TOPS;
        let (mid_base, _) = CLOUD_BASES;
        let layer_level = 0.5 * (mid_top + mid_base);
        // The low cloud's base its air's lifting condensation level, the layers' their stated
        // base.
        let low_base = base.max(low_top);
        let deck_path = |ratio: f64, depth: f64| {
            let t = self.at(air, ratio);
            0.5 * CLOUD_ADIABATIC * self.condensing(t, ratio * self.pressure) * depth * depth
        };
        let (low_depth, layer_depth) = CLOUD_DEPTH;
        // The optics per kg/m² of liquid and of ice.
        let tau_liquid = 3.0 / (2.0 * self.water.liquid_density * drop);
        let tau_ice = (CLOUD_ICE_DEPTH.0 + CLOUD_ICE_DEPTH.1 / (ice * 1e6)) * 1e3;
        let (g_liquid, g_ice) = (
            CLOUD_ASYMMETRY.0 + CLOUD_ASYMMETRY.1 * drop * 1e6,
            CLOUD_ASYMMETRY.2 + CLOUD_ASYMMETRY.3 * ice * 1e6,
        );
        let (k_liquid, k_ice) = (
            CLOUD_LIQUID_EMISSION * 1e3,
            DIFFUSIVITY * (CLOUD_ICE_EMISSION.0 + CLOUD_ICE_EMISSION.1 / (ice * 1e6)) * 1e3,
        );
        let anvil_path = CLOUD_ANVIL_DEPTH / tau_ice;
        let (cold, warm) = CLOUD_PHASE;
        let liquid_at = |ratio: f64| ((self.at(air, ratio) - cold) / (warm - cold)).clamp(0.0, 1.0);
        // Each kind's optical depth, asymmetry, emissivity and absorption from its water.
        let optics = |path: f64, liquid: f64| {
            let (tl, ti) = (liquid * path * tau_liquid, (1.0 - liquid) * path * tau_ice);
            let tau = tl + ti;
            let g = if tau > 0.0 {
                (tl * g_liquid + ti * g_ice) / tau
            } else {
                g_liquid
            };
            let e = 1.0 - (-(liquid * k_liquid + (1.0 - liquid) * k_ice) * path).exp();
            let (a, b) = CLOUD_ABSORPTION;
            (tau, g, e, a * path / (path + b))
        };
        // The anvils' thin cirrus and their medium-thick cloud as one kind: its reflection
        // the two's by their shares, its depth the one that reflects it in the diffuse light,
        // its emissivity and its absorption the two's (M67).
        let anvil = {
            let medium = kinds.medium.clamp(0.0, 1.0);
            let (thin, thick) = (
                optics(anvil_path, 0.0),
                optics(CLOUD_ANVIL_MEDIUM.2 / tau_ice, 0.0),
            );
            let diffuse = |(tau, g, _, _): (f64, f64, f64, f64)| {
                let s = (1.0 - g) * tau;
                s / (4.0 / 3.0 + s)
            };
            let r = (1.0 - medium) * diffuse(thin) + medium * diffuse(thick);
            let tau = 4.0 / 3.0 * r / (1.0 - r).max(1e-9) / (1.0 - g_ice);
            let mix = |a: f64, b: f64| (1.0 - medium) * a + medium * b;
            (tau, g_ice, mix(thin.2, thick.2), mix(thin.3, thick.3))
        };
        let parts = [
            optics(tower_path, 1.0),
            anvil,
            optics(deck_path(mid_base, layer_depth), liquid_at(layer_level)),
            optics(deck_path(low_base, low_depth), liquid_at(low_top)),
        ];
        let cover = [kinds.tower, kinds.anvil, kinds.layer, kinds.low].map(|c| c.clamp(0.0, 1.0));
        let through: f64 = (0..4).map(|i| 1.0 - cover[i] * parts[i].3).product();
        // The share of the air's absorbers above the clouds that reflect the light, each kind's
        // top by what it reflects in the diffuse light: the towers' at the tropopause, the
        // anvils' within their layer, the decks' at theirs (VAPOUR_SCALE).
        let tops = [top, (top * CLOUD_BASES.1.max(top)).sqrt(), mid_top, low_top];
        let (mut weight, mut over_them) = (0.0, 0.0);
        for i in 0..4 {
            let s = (1.0 - parts[i].1) * parts[i].0;
            let r = cover[i] * s / (4.0 / 3.0 + s);
            weight += r;
            over_them += r * tops[i].clamp(0.0, 1.0).powf(VAPOUR_SCALE);
        }
        Deck {
            cover: cover.map(|c| c as f32),
            scaled: std::array::from_fn(|i| ((1.0 - parts[i].1) * parts[i].0) as f32),
            depth: std::array::from_fn(|i| ((1.0 - parts[i].1 * parts[i].1) * parts[i].0) as f32),
            emissivity: std::array::from_fn(|i| parts[i].2 as f32),
            absorbs: (1.0 - through) as f32,
            above: if weight > 0.0 {
                (over_them / weight) as f32
            } else {
                0.0
            },
        }
    }

    /// The water a cloud's air condenses as it rises at `t` K and `p` Pa, kg/m³ a metre (M67),
    /// by `adiabatic_water`.
    fn condensing(&self, t: f64, p: f64) -> f64 {
        let air = (self.gravity, self.gas_constant, self.heat_capacity);
        adiabatic_water(self.water, air, t, p)
    }

    /// What a deck changes of a place's infrared (M66, M67), the air over it at `air` K and
    /// its ground sending up `emitted` W/m²: what the place lets out to space, from the
    /// clouds' tops, the higher hiding the lower, and what its air sends down to its ground,
    /// from their bases through the air's window, the lower hiding the higher, each with the
    /// clouds less without, W/m².
    pub fn infrared(&self, deck: &Deck, air: f64, emitted: f64) -> (f64, f64) {
        let sigma = |t: f64| SIGMA_SB * t.powi(4);
        let [c_tower, c_anvil, c_layer, c_low] = deck.cover.map(|c| c as f64);
        let [e_tower, e_anvil, e_layer, e_low] = deck.emissivity.map(|e| e as f64);
        if c_tower + c_anvil + c_layer + c_low <= 0.0 {
            return (0.0, 0.0);
        }
        let top = self.tropopause_ratio(air);
        let base = self.base_ratio(air);
        let out_from = |ratio: f64| {
            let t = self.at(air, ratio);
            sigma(t) * self.heights.share(t, ratio)
        };
        let clear = self.clear.air_up(air) + self.clear.window(air) * emitted;
        let (mid_top, low_top) = CLOUD_TOPS;
        let (mid_base, anvil_base) = CLOUD_BASES;
        // The towers let out from the tropopause, the anvils from the middle of their layer in
        // the pressure's logarithm (M67): their tops stand under the cold tropopause, near
        // 220 K over Earth's tropics (Hartmann and Larson 2002), and a thin cloud lets out from
        // within it.
        let anvil_level = (top * anvil_base.max(top)).sqrt();
        let (tower, anvil) = (c_tower * e_tower, c_anvil * e_anvil);
        let (mid, low) = (c_layer * e_layer, c_low * e_low);
        let olr = tower * out_from(top)
            + (1.0 - tower)
                * (anvil * out_from(anvil_level)
                    + (1.0 - anvil)
                        * (mid * out_from(mid_top)
                            + (1.0 - mid) * (low * out_from(low_top) + (1.0 - low) * clear)));
        let window = self.clear.window(air);
        let near = 1.0 - (1.0 - low) * (1.0 - c_tower * e_tower);
        let down = window
            * (near * sigma(self.at(air, base))
                + (1.0 - near)
                    * (mid * sigma(self.at(air, mid_base))
                        + (1.0 - mid) * c_anvil * e_anvil * sigma(self.at(air, anvil_base))));
        (olr - clear, down)
    }

    /// What clouds of `kinds` change at a place (M66), the air over it at `air` K, its ground
    /// at `ground` K of albedo `albedo`, `sea` the share of it the sea covers, the sun at the
    /// diffuse light's mean height: the planet's albedo there, its air's molecules scattering
    /// over its clouds and ground (M67), what it lets out to space, W/m², and what its air
    /// sends down to its ground, W/m², each with the clouds less without.
    pub fn effects(&self, kinds: &Kinds, air: f64, ground: f64, albedo: f64, sea: f64) -> Effects {
        let deck = self.deck(kinds, air, sea);
        let reflected = deck.reflectance(2.0 / 3.0);
        let a = albedo.clamp(0.0, 1.0);
        let over = crate::radiation::over;
        let planet = over(self.rayleigh, over(reflected, a));
        let (olr, down) = self.infrared(&deck, air, SIGMA_SB * ground.powi(4));
        Effects {
            albedo: planet - over(self.rayleigh, a),
            olr,
            down,
        }
    }
}

/// The water a cloud's air condenses as it rises at `t` K and `p` Pa, kg/m³ a metre (M67), in
/// an air of `gravity` m/s², gas constant and heat capacity J/kg/K: its adiabatic water's rise
/// (Albrecht et al. 1990), the saturation vapour's density times (p / (p - e)) (L Gamma_m / (R_v
/// T²) - g / (R T)), Gamma_m the moist adiabat's lapse rate, L the latent heat of the vapour's
/// condensing, over ice below its triple point.
pub fn adiabatic_water(
    water: &Volatile,
    (gravity, gas_constant, heat_capacity): (f64, f64, f64),
    t: f64,
    p: f64,
) -> f64 {
    let e = vapour_pressure(water, t).min(0.5 * p);
    let l = if t >= water.t_triple {
        water.l_vap
    } else {
        water.l_sub
    };
    let (g, r, cp, rv) = (gravity, gas_constant, heat_capacity, water.r_specific);
    let q = r / rv * e / (p - e);
    let lapse = g / cp * (1.0 + l * q / (r * t)) / (1.0 + l * l * q / (cp * rv * t * t));
    let vapour = e / (rv * t);
    (vapour * p / (p - e) * (l * lapse / (rv * t * t) - g / (r * t))).max(0.0)
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

    /// M66, M67: each kind grows with its cause: the towers and the high cloud they and their
    /// anvils make with the convection's rain, the anvils none under weak convection; the layers with the large-scale rain and
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
            let high = |k: &Kinds| 1.0 - (1.0 - k.tower) * (1.0 - k.anvil);
            assert!(b.tower >= a.tower && high(&b) >= high(&a), "{w:?}");
            let (a, b) = (
                kinds(0.0, w[0], 0.5, 0.0, 0.0),
                kinds(0.0, w[1], 0.5, 0.0, 0.0),
            );
            assert!(b.layer >= a.layer, "{w:?}");
        }
        assert!(kinds(100.0, 0.0, 0.5, 0.0, 0.0).tower > kinds(2.0, 0.0, 0.5, 0.0, 0.0).tower);
        assert_eq!(kinds(0.5, 0.0, 0.5, 0.0, 0.0).anvil, 0.0);
        assert!(kinds(18.0, 0.0, 0.5, 0.0, 0.0).anvil > 0.4);
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

    /// M67: the water a cloud's air condenses as it rises, as Albrecht et al. (1990)'s
    /// adiabatic rate gives it in Earth's air: within a tenth of 2.64, 2.18 and 1.85 g/m³ a km
    /// over liquid at 300 K and 900 hPa, 285 K and 900 hPa and 285 K and 700 hPa, and of 0.92 and
    /// 0.082 over ice at 260 and 230 K and 900 hPa (Murphy and Koop 2005's saturation; the
    /// research notes' table), Albrecht's measured 2.5 at some 289 K and 1000 hPa within a
    /// tenth; falling with the cold; and a low deck of CLOUD_DEPTH at its adiabaticity holding
    /// the 70 g/m² the in-cloud paths most often take at 285 K and 850 hPa, a 230 K one under
    /// 5 g/m².
    #[test]
    fn a_cloud_holds_the_water_its_air_condenses() {
        let water = volatile("H2O").unwrap();
        let air = (9.81, 287.05, 1004.0);
        let rate = |t: f64, p: f64| adiabatic_water(water, air, t, p) * 1e6;
        for (t, p, expected) in [
            (300.0, 9.0e4, 2.64),
            (285.0, 9.0e4, 2.18),
            (285.0, 7.0e4, 1.85),
            (260.0, 9.0e4, 0.92),
            (230.0, 9.0e4, 0.082),
            (289.0, 1.0e5, 2.5),
        ] {
            let found = rate(t, p);
            println!("{t} K, {p} Pa: {found:.3} g/m³ a km against {expected}");
            assert!(
                (found / expected - 1.0).abs() < 0.1,
                "{t} K {p} Pa: {found} against {expected}"
            );
        }
        for w in [220.0, 240.0, 260.0, 280.0, 300.0].windows(2) {
            assert!(rate(w[1], 8.0e4) > rate(w[0], 8.0e4), "{w:?}");
        }
        let (low, _) = CLOUD_DEPTH;
        let deck =
            |t: f64| 0.5 * CLOUD_ADIABATIC * adiabatic_water(water, air, t, 8.5e4) * low * low;
        let (warm, cold) = (deck(285.0), deck(230.0));
        println!(
            "a low deck holds {:.1} g/m² at 285 K, {:.2} at 230 K",
            warm * 1e3,
            cold * 1e3
        );
        assert!((0.055..0.085).contains(&warm), "{warm}");
        assert!(cold < 0.005, "{cold}");
    }

    /// M67: a deck reflects none of the light where it covers none, the more the thicker at
    /// every height of the sun, the diffuse light's as M66's two-stream layer did, a thick
    /// deck the more as the sun sinks, and a deck over half the place half what it reflects
    /// over all of it; never more than all of it.
    #[test]
    fn a_deck_reflects_more_as_it_thickens_and_the_sun_sinks() {
        let g = 0.85;
        let deck = |cover: f64, tau: f64| Deck {
            cover: [0.0, 0.0, cover as f32, 0.0],
            scaled: [0.0, 0.0, ((1.0 - g) * tau) as f32, 0.0],
            depth: [0.0, 0.0, ((1.0 - g * g) * tau) as f32, 0.0],
            emissivity: [0.0; 4],
            absorbs: 0.0,
            above: 0.0,
        };
        assert_eq!(deck(0.0, 10.0).reflectance(0.5), 0.0);
        for mu in [0.02, 0.1, 0.3, 0.5, 2.0 / 3.0, 1.0] {
            let mut last = 0.0;
            for tau in [0.05, 0.1, 0.5, 1.0, 3.0, 10.0, 30.0, 100.0, 1e3] {
                let r = deck(1.0, tau).reflectance(mu);
                assert!((0.0..=1.0).contains(&r), "{mu} {tau}: {r}");
                assert!(r >= last, "{mu} {tau}: {r} under {last}");
                last = r;
            }
        }
        let x = 0.75 * (1.0 - g) * 10.0;
        let diffuse = deck(1.0, 10.0).reflectance(2.0 / 3.0);
        assert!((diffuse - x / (1.0 + x)).abs() < 1e-6, "{diffuse}");
        let mut last = 0.0;
        for mu in [1.0, 0.7, 0.5, 0.3, 0.1, 0.02] {
            let r = deck(1.0, 20.0).reflectance(mu);
            assert!(r > last, "{mu}: {r} under {last}");
            last = r;
        }
        let half = deck(0.5, 10.0).reflectance(0.5);
        assert!((half - 0.5 * deck(1.0, 10.0).reflectance(0.5)).abs() < 1e-6);
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
