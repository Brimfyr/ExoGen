// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! Icy surfaces (PLAN.md Phase 4): what the surface of an ice world does on its own.
//!
//! The maturity of an icy surface is its sublimation lag. Ice sublimes at the latitude's
//! mean temperature and leaves its non-ice share behind, a dust that darkens the surface
//! once it is a few centimetres thick; a warm, old, dirty surface goes dark (Callisto), a
//! cold or young one stays bright (Enceladus, Europa), and a fresh crater that brings
//! clean ice up from under the lag is bright until its own lag forms. The poles, which the
//! mean insolation hardly warms, keep their frost.

use crate::constants::*;
use crate::fields::{Field, LayerSet};
use crate::noise::GradientNoise;
use crate::params::ParamGraph;
use crate::rng::Rng;

/// Whether the body's surface is ice, whatever its bulk is made of.
pub fn is_ice_surface(graph: &ParamGraph) -> bool {
    crate::stages::materials::is_ice(graph.str("materials.bedrock"))
}

/// Vapour pressure of water ice, Pa (Murphy and Koop 2005, valid above 110 K and close
/// enough below for a rate that is vanishing there anyway).
pub fn ice_vapour_pressure(t: f64) -> f64 {
    if t < 20.0 {
        return 0.0;
    }
    (9.550_426 - 5_723.265 / t + 3.530_68 * t.ln() - 0.007_283_32 * t).exp()
}

/// Free sublimation flux of water ice at a temperature, kg/m²/s: the Hertz-Knudsen rate
/// at the ice's vapour pressure.
pub fn ice_sublimation_flux(t: f64) -> f64 {
    if t < 20.0 {
        return 0.0;
    }
    ice_vapour_pressure(t) * (WATER_MOLECULE_MASS / (2.0 * std::f64::consts::PI * K_B * t)).sqrt()
}

/// Diurnal-mean temperature of an airless surface at a latitude: the mean insolation,
/// S cos(lat) / π, radiated away, no lower than the polar floor.
pub fn mean_surface_temperature(insolation: f64, bond_albedo: f64, cos_lat: f64) -> f64 {
    let absorbed = insolation * (1.0 - bond_albedo) * cos_lat.max(0.0) / std::f64::consts::PI;
    (absorbed / SIGMA_SB)
        .powf(0.25)
        .max(ICE_POLAR_TEMPERATURE_FLOOR_K)
}

/// The non-ice share of the regolith mix, which is what the lag is made of.
pub fn regolith_impurity(graph: &ParamGraph) -> f64 {
    match graph.get("materials.regolith_mix").as_map() {
        Some(mix) if !mix.is_empty() => {
            let total: f64 = mix.values().sum();
            if total <= 0.0 {
                return 1.0;
            }
            mix.iter()
                .map(|(id, share)| share * crate::stages::materials::non_ice_fraction(id))
                .sum::<f64>()
                / total
        }
        _ => crate::stages::materials::non_ice_fraction(graph.str("materials.bedrock")),
    }
}

/// Maturity of a cell of icy surface at an exposure: the lag its sublimation left,
/// which chokes its own growth once it is thick (the ice under it sublimes through it
/// ever more slowly, so the lag grows as the square root of time), opaque past the
/// reference mass; and the gardening that mixes clean ice back into the regolith.
pub fn lag_maturity_at(temperature: f64, impurity: f64, exposure_gyr: f64) -> f64 {
    let exposure = exposure_gyr.max(0.0);
    let supply = ice_sublimation_flux(temperature) * exposure * SECONDS_PER_GYR * impurity
        / LAG_OPAQUE_KG_M2;
    let thickness = (1.0 + 2.0 * supply).sqrt() - 1.0;
    let lag = 1.0 - (-thickness).exp();
    let gardened = 1.0 - (-exposure / ICE_EJECTA_FADE_GYR).exp();
    1.0 - (1.0 - lag) * (1.0 - gardened)
}

/// Thermal segregation: the share of a latitude covered by frost that condensed there
/// out of what the warm latitudes sublimed over the body's life.
pub struct FrostCover {
    insolation: f64,
    albedo: f64,
    supply: f64,
}

impl FrostCover {
    pub fn new(graph: &ParamGraph) -> Self {
        let insolation = graph.f64("orbit.insolation");
        let albedo = graph.f64("atmosphere.bond_albedo");
        let t_eq = mean_surface_temperature(insolation, albedo, 1.0);
        let sublimed = ice_sublimation_flux(t_eq) * graph.f64("body.age") * SECONDS_PER_GYR;
        Self {
            insolation,
            albedo,
            supply: (sublimed / FROST_SUPPLY_KG_M2).min(1.0),
        }
    }

    /// Frost share at a direction on the sphere.
    pub fn at(&self, d: [f64; 3]) -> f64 {
        let cos_lat = (1.0 - d[2] * d[2]).max(0.0).sqrt();
        let t = mean_surface_temperature(self.insolation, self.albedo, cos_lat);
        let persist = 1.0 / (1.0 + ((t - FROST_HALF_K) / FROST_WIDTH_K).exp());
        FROST_COVER_MAX * persist * self.supply
    }
}

/// The maturity field of an icy surface from the surface ages. Ice renewed from below,
/// the grooved lanes and the polar terrain, is clean whatever the regolith around it.
pub fn lag_maturity(graph: &ParamGraph, layers: &LayerSet, present: f64) -> Field<f32> {
    let grid = layers.grid;
    let insolation = graph.f64("orbit.insolation");
    let albedo = graph.f64("atmosphere.bond_albedo");
    let impurity = regolith_impurity(graph);
    let clean = crate::stages::materials::non_ice_fraction("water_ice");
    let clean_f = layers.clean.as_slice();
    Field::from_fn(grid, |c, d| {
        let cos_lat = (1.0 - d[2] * d[2]).max(0.0).sqrt();
        let t = mean_surface_temperature(insolation, albedo, cos_lat);
        let exposure = present - layers.surface_age.get(c) as f64;
        let share = clean_f[grid.index(c)].clamp(0.0, 1.0) as f64;
        let impurity = impurity + (clean - impurity) * share;
        lag_maturity_at(t, impurity, exposure) as f32
    })
}

/// The shell's thickness across the body (F43): thinned under one pole where the tidal
/// heating focuses, as a share of the mean, tapering to the mean away from it; the
/// terrain over the thinnest part is renewed.
#[derive(Clone, Debug)]
pub struct ShellMap {
    pub pole: Option<[f64; 3]>,
    pub thinning: f64,
    pub radius_rad: f64,
    /// The field that wobbles the thinning's edge round the pole.
    noise: GradientNoise,
}

impl ShellMap {
    pub fn new(graph: &ParamGraph) -> Self {
        let thinning = graph.f64("tectonics.polar_thinning").clamp(0.0, 0.95);
        let seed = graph.f64("seed.tectonics") as u64;
        let noise = GradientNoise::new(seed ^ 0x5027_B0B1);
        if thinning <= 0.0 {
            return Self {
                pole: None,
                thinning: 0.0,
                radius_rad: 0.0,
                noise,
            };
        }
        // Which pole is a coin the body tossed long ago.
        let rng = Rng::new(seed ^ 0x5027_0000);
        let south = rng.bits(0, 0) & 1 == 0;
        Self {
            pole: Some([0.0, 0.0, if south { -1.0 } else { 1.0 }]),
            thinning,
            radius_rad: POLAR_TERRAIN_RADIUS_DEG.to_radians(),
            noise,
        }
    }

    /// The taper's radius in a direction, rad: the mean, wobbled round the pole.
    pub fn radius_at(&self, d: [f64; 3]) -> f64 {
        self.radius_rad * (1.0 + POLAR_MARGIN_WOBBLE * self.noise.fbm(d, 2.0, 2, 2.0, 0.5))
    }

    /// Angular radius of the active terrain in a direction, rad.
    pub fn active_radius_at(&self, d: [f64; 3]) -> f64 {
        if self.radius_rad <= 0.0 {
            return 0.0;
        }
        self.active_radius_rad() * self.radius_at(d) / self.radius_rad
    }

    /// The focusing profile at a direction: 1 at the pole, 0 beyond the taper.
    pub fn profile(&self, d: [f64; 3]) -> f64 {
        let Some(pole) = self.pole else {
            return 0.0;
        };
        let x = crate::geom::angle_between(pole, d) / self.radius_at(d).max(1e-9);
        if x >= 1.0 {
            0.0
        } else {
            0.5 * (1.0 + (std::f64::consts::PI * x).cos())
        }
    }

    /// The local shell thickness over the mean.
    pub fn factor(&self, d: [f64; 3]) -> f64 {
        1.0 - self.thinning * self.profile(d)
    }

    /// Whether the thinning is deep enough for the terrain over it to renew itself.
    pub fn has_active_terrain(&self) -> bool {
        self.pole.is_some() && self.thinning > POLAR_ACTIVE_THINNING
    }

    /// Angular radius of the active terrain, rad: where the profile crosses the active
    /// share.
    pub fn active_radius_rad(&self) -> f64 {
        (2.0 * POLAR_ACTIVE_PROFILE - 1.0).clamp(-1.0, 1.0).acos() / std::f64::consts::PI
            * self.radius_rad
    }

    /// Whether the terrain at a direction is the renewed polar terrain.
    pub fn active(&self, d: [f64; 3]) -> bool {
        self.has_active_terrain() && self.profile(d) > POLAR_ACTIVE_PROFILE
    }

    /// How much of a cell at a direction the active terrain covers, 0 to 1: a one-cell
    /// ramp across its edge, for an edge drawn without stair steps.
    pub fn active_coverage(&self, d: [f64; 3], body_radius: f64, cell: f64) -> f64 {
        let Some(pole) = self.pole else {
            return 0.0;
        };
        if !self.has_active_terrain() {
            return 0.0;
        }
        let angle = crate::geom::angle_between(pole, d);
        ((self.active_radius_at(d) - angle) * body_radius / cell.max(1.0) + 0.5).clamp(0.0, 1.0)
    }

    /// How far below the mean surface the thinner shell floats, m: the isostatic
    /// depression of the thinning, less what the lid holds up.
    pub fn depression(&self, d: [f64; 3], mean_thickness: f64) -> f64 {
        let thinned = mean_thickness * self.thinning * self.profile(d);
        POLAR_DEPRESSION_FRACTION * (DENSITY_ICE - DENSITY_CRUST_ICE) / DENSITY_ICE * thinned
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ice_sublimes_fast_at_callisto_and_not_at_europa() {
        // Murphy and Koop: 2.8e-11 Pa at 115 K, about 1e-6 Pa at 140 K.
        let p115 = ice_vapour_pressure(115.0);
        assert!(p115 > 1e-11 && p115 < 1e-10, "{p115}");
        let j = ice_sublimation_flux(115.0);
        assert!(j > 1e-14 && j < 2e-13, "{j}");
        assert!(ice_sublimation_flux(92.0) < 1e-17);
        // A dirty surface at Callisto's equatorial mean goes opaque over the age; a
        // clean one at Europa's stays bright; the poles stay bright anywhere.
        let callisto = mean_surface_temperature(50.25, 0.22, 1.0);
        let europa = mean_surface_temperature(50.25, 0.67, 1.0);
        assert!(callisto > 110.0 && callisto < 125.0, "{callisto}");
        assert!(europa < 100.0, "{europa}");
        assert!(lag_maturity_at(callisto, 0.7, 4.5) > 0.95);
        assert!(lag_maturity_at(callisto, 0.7, 0.05) < 0.6);
        assert!(lag_maturity_at(europa, 0.02, 0.06) < 0.15);
        // Old polar ice is gardened but carries no lag; the equator carries both.
        let pole = lag_maturity_at(mean_surface_temperature(50.25, 0.22, 0.05), 0.7, 4.5);
        assert!(
            pole > 0.9 && pole < lag_maturity_at(callisto, 0.7, 4.5),
            "{pole}"
        );
        let mid = lag_maturity_at(mean_surface_temperature(50.25, 0.22, 0.5), 0.7, 0.2);
        let eq = lag_maturity_at(callisto, 0.7, 0.2);
        assert!(mid < eq, "{mid} vs {eq}");
        // Frost: a Ganymede-like body wears caps past 40° and nothing at the equator.
        let frost = FrostCover {
            insolation: 50.25,
            albedo: 0.43,
            supply: 1.0,
        };
        let lat = |deg: f64| [deg.to_radians().cos(), 0.0, deg.to_radians().sin()];
        assert!(frost.at(lat(0.0)) < 0.05 && frost.at(lat(70.0)) > 0.3);
        assert!(frost.at(lat(70.0)) > frost.at(lat(45.0)));
    }
}
