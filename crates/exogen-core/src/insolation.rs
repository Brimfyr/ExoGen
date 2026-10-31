// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! The star's light on a body (PLAN.md 4.6, Phase 5A, F234): the insolation of any rotation
//! state from one integral, the star's flux on the ground averaged over its path across the
//! sky. A body that spins many times an orbit sees the star cross every longitude alike, so
//! its sunlight is the day's mean at each latitude (Berger 1978) taken through the seasons
//! its obliquity and eccentricity make; a body whose spin is in an exact resonance with its
//! orbit sees the same sky every few orbits, so its sunlight keeps the longitudes the
//! resonance fixes: the substellar hemisphere of a body locked to its star, Mercury's hot
//! longitudes under its 3:2 resonance. A moon's sky is its planet's: its year and its
//! seasons are its planet's, its day its own orbit.
//!
//! By Kepler's second law the year's mean is a mean over the star's longitude divided by
//! (1 - e^2)^(1/2), whatever the season of perihelion (Ward 1974), so the seasons' timing
//! enters only a resonant body's sunlight, through where its hot longitudes face at
//! perihelion.

use std::f64::consts::{PI, TAU};

use rayon::prelude::*;

use crate::constants::*;
use crate::fields::Field;
use crate::grid::{CubeSphere, lon_lat_from_direction};
use crate::params::ParamGraph;

/// The year's sunlight over a grid (F235): its mean flux and the sun's mean height.
#[derive(Clone, Debug)]
pub struct Sunlight {
    /// The year's mean flux on the ground, W/m².
    pub flux: Field<f32>,
    /// The cosine of the star's zenith angle averaged over the year, each moment weighted by
    /// the flux it brings: 1 with the star always overhead, 0 where it never shines.
    pub height: Field<f32>,
}

/// The star's path across a body's sky: everything its insolation depends on.
#[derive(Clone, Copy, Debug)]
pub struct Sky {
    /// The star's flux at the semi-major axis of the orbit about it, W/m².
    pub flux: f64,
    /// The eccentricity of that orbit (a moon's planet's).
    pub eccentricity: f64,
    /// The spin axis's mean tilt from that orbit's normal, rad.
    pub obliquity: f64,
    /// The half-angle of the cone the axis precesses on about its mean tilt, faster than the
    /// loop's steps, rad: a moon's orbit circling its planet's equator.
    pub swing: f64,
    /// The star's longitude seen from the body at perihelion, rad, 0 at the northern spring
    /// equinox.
    pub perihelion_season: f64,
    /// Sidereal spins per orbit about the star.
    pub spins: f64,
    /// An exact spin-orbit resonance with the star: spins and orbits to repeat.
    pub resonance: Option<(u32, u32)>,
    /// The longitude under the star at perihelion, rad, for a resonant body.
    pub subsolar_longitude: f64,
}

impl Sky {
    pub fn from_graph(graph: &ParamGraph) -> Self {
        let spins = graph.f64("rotation.spins_per_year");
        Self {
            flux: graph.f64("orbit.insolation"),
            eccentricity: graph.f64("orbit.star_eccentricity"),
            obliquity: graph.f64("rotation.effective_obliquity").to_radians(),
            swing: graph.f64("rotation.obliquity_swing").to_radians(),
            perihelion_season: graph.f64("orbit.perihelion_season").to_radians(),
            spins,
            resonance: if graph.f64("rotation.resonance") > 0.0 {
                spin_orbit_resonance(spins)
            } else {
                None
            },
            subsolar_longitude: graph.f64("climate.subsolar_longitude").to_radians(),
        }
    }

    /// The star's flux and declination at a mean anomaly from perihelion, for a tilt (F236).
    pub fn at_mean_anomaly(&self, m: f64, eps: f64) -> (f64, f64) {
        let e = self.eccentricity;
        let big_e = eccentric_anomaly(m, e);
        let r = 1.0 - e * big_e.cos();
        let nu = 2.0
            * ((1.0 + e).sqrt() * (big_e / 2.0).sin())
                .atan2((1.0 - e).sqrt() * (big_e / 2.0).cos());
        let ls = nu + self.perihelion_season;
        let decl = (eps.sin() * ls.sin()).clamp(-1.0, 1.0).asin();
        (self.flux / (r * r), decl)
    }

    /// The obliquities the axis takes as it precesses on its cone, sampled evenly through
    /// the precession: one when it does not swing.
    pub fn obliquities(&self) -> Vec<f64> {
        if self.swing <= 0.0 {
            return vec![self.obliquity];
        }
        let (sa, ca) = self.obliquity.sin_cos();
        let (sb, cb) = self.swing.sin_cos();
        (0..INSOLATION_SWING_SAMPLES)
            .map(|j| {
                let node = (j as f64 + 0.5) / INSOLATION_SWING_SAMPLES as f64 * TAU;
                (ca * cb - sa * sb * node.cos()).clamp(-1.0, 1.0).acos()
            })
            .collect()
    }

    /// The year's mean flux on the ground at a latitude and longitude, W/m².
    pub fn annual_mean_at(&self, lat: f64, lon: f64) -> f64 {
        self.flux * self.shares_at(lat, lon).0
            / (1.0 - self.eccentricity * self.eccentricity).sqrt()
    }

    /// The sun's year-mean height at a latitude and longitude, weighted by its flux (F235).
    pub fn annual_height_at(&self, lat: f64, lon: f64) -> f64 {
        let (one, two) = self.shares_at(lat, lon);
        if one > 0.0 { two / one } else { 0.0 }
    }

    /// The year's mean share of the overhead flux at a point, and of that share times the
    /// zenith angle's cosine, over the axis's precession.
    fn shares_at(&self, lat: f64, lon: f64) -> (f64, f64) {
        let tilts = self.obliquities();
        let each: Vec<(f64, f64)> = match self.resonance {
            None => tilts
                .iter()
                .map(|&eps| zonal_annual_shares(lat, eps))
                .collect(),
            Some((_, orbits)) => tilts
                .iter()
                .map(|&eps| path_shares(&self.resonant_path(eps, orbits), lat, lon))
                .collect(),
        };
        average_pairs(&each)
    }

    /// The year's mean flux on the ground over a grid, W/m².
    pub fn annual_mean(&self, grid: CubeSphere) -> Field<f32> {
        self.annual_sunlight(grid).flux
    }

    /// The year's sunlight over a grid, its flux and the sun's mean height, from a table in
    /// latitude (a fast rotator) or in latitude and longitude (a resonant body) read
    /// bilinearly, so that the fields are the same at every resolution.
    pub fn annual_sunlight(&self, grid: CubeSphere) -> Sunlight {
        let scale = self.flux / (1.0 - self.eccentricity * self.eccentricity).sqrt();
        let tilts = self.obliquities();
        let n_lat = INSOLATION_TABLE_LATITUDES;
        let lat_of = |i: usize| -PI / 2.0 + PI * i as f64 / (n_lat - 1) as f64;
        let height = |(one, two): (f64, f64)| if one > 0.0 { two / one } else { 0.0 };
        match self.resonance {
            None => {
                let table: Vec<(f64, f64)> = (0..n_lat)
                    .into_par_iter()
                    .map(|i| {
                        let each: Vec<(f64, f64)> = tilts
                            .iter()
                            .map(|&eps| zonal_annual_shares(lat_of(i), eps))
                            .collect();
                        average_pairs(&each)
                    })
                    .collect();
                let shares: Vec<f64> = table.iter().map(|p| p.0).collect();
                let heights: Vec<f64> = table.iter().map(|&p| height(p)).collect();
                Sunlight {
                    flux: Field::from_fn(grid, |_, d| {
                        let (_, lat) = lon_lat_from_direction(d);
                        (scale * read_row(&shares, lat)) as f32
                    }),
                    height: Field::from_fn(grid, |_, d| {
                        let (_, lat) = lon_lat_from_direction(d);
                        read_row(&heights, lat) as f32
                    }),
                }
            }
            Some((_, orbits)) => {
                let n_lon = INSOLATION_TABLE_LONGITUDES;
                let paths: Vec<Vec<PathStep>> = tilts
                    .iter()
                    .map(|&eps| self.resonant_path(eps, orbits))
                    .collect();
                let table: Vec<(f64, f64)> = (0..n_lat * n_lon)
                    .into_par_iter()
                    .map(|k| {
                        let (i, j) = (k / n_lon, k % n_lon);
                        let lon = -PI + TAU * j as f64 / n_lon as f64;
                        let each: Vec<(f64, f64)> = paths
                            .iter()
                            .map(|p| path_shares(p, lat_of(i), lon))
                            .collect();
                        average_pairs(&each)
                    })
                    .collect();
                let shares: Vec<f64> = table.iter().map(|p| p.0).collect();
                let heights: Vec<f64> = table.iter().map(|&p| height(p)).collect();
                Sunlight {
                    flux: Field::from_fn(grid, |_, d| {
                        let (lon, lat) = lon_lat_from_direction(d);
                        (scale * read_grid(&shares, n_lon, lat, lon)) as f32
                    }),
                    height: Field::from_fn(grid, |_, d| {
                        let (lon, lat) = lon_lat_from_direction(d);
                        read_grid(&heights, n_lon, lat, lon) as f32
                    }),
                }
            }
        }
    }

    /// The star's path through a resonant body's sky over the orbits it takes to repeat,
    /// sampled evenly in time: at each step its flux over the mean's, its declination, and
    /// the hour angle at longitude zero.
    pub(crate) fn resonant_path(&self, eps: f64, orbits: u32) -> Vec<PathStep> {
        let e = self.eccentricity;
        let n = INSOLATION_ORBIT_SAMPLES * orbits as usize;
        let (se, ce) = eps.sin_cos();
        let ls_p = self.perihelion_season;
        // The rotation's phase puts the chosen longitude under the star at perihelion.
        let ra_p = (ls_p.sin() * ce).atan2(ls_p.cos());
        let theta0 = ra_p - self.subsolar_longitude;
        (0..n)
            .map(|j| {
                let m = (j as f64 + 0.5) / INSOLATION_ORBIT_SAMPLES as f64 * TAU;
                let big_e = eccentric_anomaly(m, e);
                let r = 1.0 - e * big_e.cos();
                let nu = 2.0
                    * ((1.0 + e).sqrt() * (big_e / 2.0).sin())
                        .atan2((1.0 - e).sqrt() * (big_e / 2.0).cos());
                // Without the mean distance's factor: the caller divides by (1 - e^2)^(1/2)
                // for every body alike, and the time mean of 1/r^2 is that factor.
                let flux = (1.0 - e * e).sqrt() / (r * r);
                let ls = nu + ls_p;
                let decl = (se * ls.sin()).clamp(-1.0, 1.0).asin();
                let ra = (ls.sin() * ce).atan2(ls.cos());
                let theta = theta0 + self.spins * m;
                PathStep {
                    flux,
                    sin_decl: decl.sin(),
                    cos_decl: decl.cos(),
                    hour: theta - ra,
                }
            })
            .collect()
    }
}

/// One step of the star's path: its flux, declination and the hour angle at longitude zero.
#[derive(Clone, Copy, Debug)]
pub(crate) struct PathStep {
    /// The flux over its mean, times (1 - e^2)^(1/2).
    pub(crate) flux: f64,
    pub(crate) sin_decl: f64,
    pub(crate) cos_decl: f64,
    pub(crate) hour: f64,
}

/// The mean of the flux on the ground at a point over a path, as a share of the flux
/// overhead at the mean distance, and of that share times the zenith angle's cosine.
fn path_shares(path: &[PathStep], lat: f64, lon: f64) -> (f64, f64) {
    let (sp, cp) = lat.sin_cos();
    let (one, two) = path.iter().fold((0.0, 0.0), |(one, two), s| {
        let cz = (sp * s.sin_decl + cp * s.cos_decl * (s.hour + lon).cos()).max(0.0);
        (one + s.flux * cz, two + s.flux * cz * cz)
    });
    (one / path.len() as f64, two / path.len() as f64)
}

/// The mean of pairs.
fn average_pairs(each: &[(f64, f64)]) -> (f64, f64) {
    let n = each.len().max(1) as f64;
    let (a, b) = each.iter().fold((0.0, 0.0), |(a, b), p| (a + p.0, b + p.1));
    (a / n, b / n)
}

/// The day's mean flux on the ground at a latitude, as a share of the flux overhead, with
/// the star at a declination: the integral of the sun's height from sunrise to sunset
/// (Berger 1978), whole days of light or night inside the polar circles.
pub fn diurnal_share(lat: f64, decl: f64) -> f64 {
    diurnal_shares(lat, decl).0
}

/// The day's mean share of the overhead flux, and of that share times the zenith angle's
/// cosine (F235), from sunrise to sunset.
pub fn diurnal_shares(lat: f64, decl: f64) -> (f64, f64) {
    let (sp, cp) = lat.sin_cos();
    let (sd, cd) = decl.sin_cos();
    let (a, b) = (sp * sd, cp * cd);
    let h0 = if a >= b {
        PI
    } else if a <= -b {
        0.0
    } else {
        (-a / b).acos()
    };
    let (sh, ch) = h0.sin_cos();
    let one = (h0 * a + b * sh) / PI;
    let two = (a * a * h0 + 2.0 * a * b * sh + b * b * (h0 + sh * ch) / 2.0) / PI;
    (one.max(0.0), two.max(0.0))
}

/// The year's mean of the day's share at a latitude on a circular orbit, the star's
/// longitude sampled evenly; an eccentric orbit's is this over (1 - e^2)^(1/2).
pub fn zonal_annual_share(lat: f64, obliquity: f64) -> f64 {
    zonal_annual_shares(lat, obliquity).0
}

/// The year's means of the day's two shares at a latitude on a circular orbit (F235).
pub fn zonal_annual_shares(lat: f64, obliquity: f64) -> (f64, f64) {
    let n = INSOLATION_SEASON_SAMPLES;
    let s = obliquity.sin();
    let each: Vec<(f64, f64)> = (0..n)
        .map(|j| {
            let ls = (j as f64 + 0.5) / n as f64 * TAU;
            diurnal_shares(lat, (s * ls.sin()).clamp(-1.0, 1.0).asin())
        })
        .collect();
    average_pairs(&each)
}

/// The eccentric anomaly of a mean anomaly, by Newton's method from E = M + e sin M.
pub fn eccentric_anomaly(m: f64, e: f64) -> f64 {
    let mut big_e = m + e * m.sin();
    for _ in 0..50 {
        let step = (big_e - e * big_e.sin() - m) / (1.0 - e * big_e.cos());
        big_e -= step;
        if step.abs() < 1e-14 {
            break;
        }
    }
    big_e
}

/// An exact spin-orbit resonance: the fewest orbits after which a spin of `spins` sidereal
/// turns an orbit shows the star the same face again, with the turns it makes in them, for
/// the low orders tides capture a spin into (Goldreich and Peale 1966); none otherwise.
pub fn spin_orbit_resonance(spins: f64) -> Option<(u32, u32)> {
    if !spins.is_finite() || spins <= 0.0 {
        return None;
    }
    (1..=RESONANCE_MAX_ORBITS).find_map(|q| {
        let turns = spins * q as f64;
        let p = turns.round();
        (p >= 1.0 && p <= RESONANCE_MAX_SPINS as f64 && (turns - p).abs() < RESONANCE_DRIFT_TURNS)
            .then_some((p as u32, q))
    })
}

/// A table over latitude from pole to pole, read linearly.
pub(crate) fn read_row(table: &[f64], lat: f64) -> f64 {
    let n = table.len();
    let x = ((lat + PI / 2.0) / PI * (n - 1) as f64).clamp(0.0, (n - 1) as f64);
    let i = (x.floor() as usize).min(n - 2);
    let f = x - i as f64;
    table[i] * (1.0 - f) + table[i + 1] * f
}

/// A table over latitude (rows, pole to pole) and longitude (columns, from -180°, wrapping),
/// read bilinearly.
pub(crate) fn read_grid(table: &[f64], n_lon: usize, lat: f64, lon: f64) -> f64 {
    let n_lat = table.len() / n_lon;
    let y = ((lat + PI / 2.0) / PI * (n_lat - 1) as f64).clamp(0.0, (n_lat - 1) as f64);
    let i = (y.floor() as usize).min(n_lat - 2);
    let fy = y - i as f64;
    let x = (lon + PI).rem_euclid(TAU) / TAU * n_lon as f64;
    let j = (x.floor() as usize) % n_lon;
    let fx = x - x.floor();
    let j1 = (j + 1) % n_lon;
    let at = |r: usize, c: usize| table[r * n_lon + c];
    let top = at(i, j) * (1.0 - fx) + at(i, j1) * fx;
    let bottom = at(i + 1, j) * (1.0 - fx) + at(i + 1, j1) * fx;
    top * (1.0 - fy) + bottom * fy
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sky(e: f64, eps_deg: f64, spins: f64) -> Sky {
        Sky {
            flux: 1361.0,
            eccentricity: e,
            obliquity: eps_deg.to_radians(),
            swing: 0.0,
            perihelion_season: 0.0,
            spins,
            resonance: spin_orbit_resonance(spins),
            subsolar_longitude: 0.0,
        }
    }

    /// The area-weighted mean over a grid.
    fn global_mean(f: &Field<f32>) -> f64 {
        let g = f.grid();
        let (mut sum, mut area) = (0.0, 0.0);
        for c in g.cells() {
            let a = g.cell_area(c);
            sum += f.get(c) as f64 * a;
            area += a;
        }
        sum / area
    }

    /// F234: the year's sunlight at the Earth's equator and poles, against a brute-force
    /// integration of the star's path with the day resolved (scratchpad insol_ref.py:
    /// 415.89 and 172.35 W/m² at 1361), and the pole's closed form S sin(eps) / pi.
    #[test]
    fn the_earth_s_equator_and_poles_take_their_year_s_sunlight() {
        let earth = sky(0.0167, 23.44, 366.256);
        assert!(earth.resonance.is_none());
        let equator = earth.annual_mean_at(0.0, 0.0);
        let pole = earth.annual_mean_at(PI / 2.0, 0.0);
        assert!((equator / 415.89 - 1.0).abs() < 1e-3, "equator {equator}");
        assert!((pole / 172.35 - 1.0).abs() < 1e-3, "pole {pole}");
        let closed = 1361.0 * 23.44f64.to_radians().sin() / PI / (1.0 - 0.0167f64.powi(2)).sqrt();
        assert!(
            (pole / closed - 1.0).abs() < 1e-4,
            "pole {pole} closed form {closed}"
        );
    }

    /// F234: every body takes a quarter of its year's mean flux on average, whatever its
    /// tilt, spin or resonance.
    #[test]
    fn every_body_takes_a_quarter_of_its_flux_on_average() {
        let grid = CubeSphere::new(32);
        let mut triton = sky(0.0086, 151.68, 3_000.0);
        triton.swing = 23.0f64.to_radians();
        for (name, s) in [
            ("Earth", sky(0.0167, 23.44, 366.256)),
            ("Mars", sky(0.0934, 25.19, 669.6)),
            ("Mercury", sky(0.2056, 0.03, 1.5)),
            ("locked", sky(0.05, 0.0, 1.0)),
            ("Pluto", sky(0.2488, 119.6, 14_178.0)),
            ("Triton", triton),
        ] {
            let f = s.annual_mean(grid);
            let expect = s.flux / 4.0 / (1.0 - s.eccentricity * s.eccentricity).sqrt();
            let mean = global_mean(&f);
            assert!(
                (mean / expect - 1.0).abs() < 2e-3,
                "{name}: {mean} against {expect}"
            );
        }
    }

    /// F234: past about 54° of obliquity the poles take more of the year's sunlight than the
    /// equator (Ward 1974; insol_ref.py: 54° gives the two within 0.2 percent).
    #[test]
    fn a_steep_axis_warms_the_poles_past_the_equator() {
        let at = |eps: f64| {
            let s = sky(0.0, eps, 1000.37);
            (s.annual_mean_at(0.0, 0.0), s.annual_mean_at(PI / 2.0, 0.0))
        };
        let (eq, pole) = at(50.0);
        assert!(eq > pole * 1.05, "50°: equator {eq} pole {pole}");
        let (eq, pole) = at(60.0);
        assert!(pole > eq * 1.05, "60°: equator {eq} pole {pole}");
        let (eq, pole) = at(119.6);
        assert!(pole > eq * 1.05, "Pluto's 119.6°: equator {eq} pole {pole}");
    }

    /// F234: Mercury's 3:2 resonance keeps its hot longitudes, under the star at perihelion,
    /// at 2.53 times the year's sunlight of the warm longitudes between them (insol_ref.py:
    /// 0.42641 and 0.16868 of the flux at the mean distance).
    #[test]
    fn mercury_s_resonance_keeps_its_hot_longitudes() {
        let mercury = sky(0.2056, 0.0, 1.5);
        assert_eq!(mercury.resonance, Some((3, 2)));
        let hot = mercury.annual_mean_at(0.0, 0.0) / 1361.0;
        let warm = mercury.annual_mean_at(0.0, PI / 2.0) / 1361.0;
        let other_hot = mercury.annual_mean_at(0.0, PI) / 1361.0;
        assert!((hot / 0.42641 - 1.0).abs() < 2e-3, "hot {hot}");
        assert!((warm / 0.16868 - 1.0).abs() < 2e-3, "warm {warm}");
        assert!((other_hot / hot - 1.0).abs() < 1e-3);
        // The field read from its table agrees with the direct integral.
        let f = mercury.annual_mean(CubeSphere::new(64));
        let g = f.grid();
        let (c, _, _) = g.locate([1.0, 0.0, 0.0]);
        assert!((f.get(c) as f64 / 1361.0 / hot - 1.0).abs() < 5e-3);
    }

    /// F234: a body locked to its star faces it: the substellar point takes the full flux,
    /// 60° from it half, the far side none (insol_ref.py; an eccentric orbit's libration
    /// lights the terminator a little, 0.064 of the flux at e = 0.1).
    #[test]
    fn a_body_locked_to_its_star_faces_it() {
        let locked = sky(0.0, 0.0, 1.0);
        assert_eq!(locked.resonance, Some((1, 1)));
        let at = |s: &Sky, lon: f64| s.annual_mean_at(0.0, lon.to_radians()) / 1361.0;
        assert!((at(&locked, 0.0) - 1.0).abs() < 1e-3);
        assert!((at(&locked, 60.0) - 0.5).abs() < 1e-3);
        assert!(at(&locked, 180.0) < 1e-9);
        let eccentric = sky(0.1, 0.0, 1.0);
        assert!((at(&eccentric, 90.0) / 0.06373 - 1.0).abs() < 1e-2);
    }

    /// Prints each preset's year of sunlight: the mean flux by latitude on the prime meridian
    /// and, for a resonant body, around its equator.
    #[test]
    #[ignore = "diagnostic: the presets' insolation by latitude"]
    fn print_insolation() {
        for name in crate::presets::PRESET_NAMES {
            let sky = Sky::from_graph(&crate::presets::preset(name).unwrap().to_graph().0);
            let row: Vec<String> = [-90.0, -60.0, -30.0, 0.0, 30.0, 60.0, 90.0]
                .iter()
                .map(|lat: &f64| format!("{:.4}", sky.annual_mean_at(lat.to_radians(), 0.0)))
                .collect();
            println!(
                "{name:10} tilt {:6.2}° swing {:5.2}° spins {:10.3} resonance {:?}: {}",
                sky.obliquity.to_degrees(),
                sky.swing.to_degrees(),
                sky.spins,
                sky.resonance,
                row.join(" ")
            );
            if sky.resonance.is_some() {
                let around: Vec<String> = [0.0, 45.0, 90.0, 135.0, 180.0]
                    .iter()
                    .map(|lon: &f64| format!("{:.1}", sky.annual_mean_at(0.0, lon.to_radians())))
                    .collect();
                println!(
                    "{:10} around the equator from 0° by 45°: {}",
                    "",
                    around.join(" ")
                );
            }
        }
    }

    /// F234: the presets' skies: Mercury in its 3:2, a moon on its planet's year with its
    /// planet's tilt, the Moon at its own 1.54°, Triton's axis swinging.
    #[test]
    fn the_presets_see_their_star_s_path() {
        let graph = |name: &str| crate::presets::preset(name).unwrap().to_graph().0;
        let mercury = Sky::from_graph(&graph("Mercury"));
        assert_eq!(mercury.resonance, Some((3, 2)));
        let earth = Sky::from_graph(&graph("Earth"));
        assert!(earth.resonance.is_none() && (earth.spins - 366.26).abs() < 0.05);
        let titan_graph = graph("Titan");
        let year_days = titan_graph.f64("orbit.year") / 86_400.0;
        assert!(
            (year_days / 10_759.0 - 1.0).abs() < 0.01,
            "Titan's year {year_days} d"
        );
        let titan = Sky::from_graph(&titan_graph);
        assert!(titan.resonance.is_none());
        assert!((titan.obliquity.to_degrees() - 26.73).abs() < 1e-9);
        assert!((titan.eccentricity - 0.0565).abs() < 1e-9);
        let moon = Sky::from_graph(&graph("Moon"));
        assert!((moon.obliquity.to_degrees() - 1.54).abs() < 1e-9);
        let triton = Sky::from_graph(&graph("Triton"));
        let tilts: Vec<f64> = triton
            .obliquities()
            .iter()
            .map(|e| e.to_degrees())
            .collect();
        let (lo, hi) = tilts
            .iter()
            .fold((180.0f64, 0.0f64), |(a, b), &t| (a.min(t), b.max(t)));
        assert!(lo > 128.0 && hi < 175.0, "Triton's tilts {lo} to {hi}");
    }

    /// F235: the sun stands higher over the tropics than the poles: on the Earth its
    /// year-mean height, weighted by the flux, is 0.754 at the equator and 0.313 at the poles
    /// (ebm_proto3.py); a body locked to its star has it overhead at the substellar point.
    #[test]
    fn the_sun_stands_higher_over_the_tropics() {
        let earth = sky(0.0167, 23.44, 366.256);
        let equator = earth.annual_height_at(0.0, 0.0);
        let pole = earth.annual_height_at(PI / 2.0 - 0.004, 0.0);
        assert!((equator - 0.754).abs() < 0.005, "equator {equator}");
        assert!((pole - 0.313).abs() < 0.01, "pole {pole}");
        let locked = sky(0.0, 0.0, 1.0);
        assert!((locked.annual_height_at(0.0, 0.0) - 1.0).abs() < 1e-6);
        assert!((locked.annual_height_at(0.0, 60f64.to_radians()) - 0.5).abs() < 1e-6);
        assert_eq!(locked.annual_height_at(0.0, PI), 0.0);
        // The field's table agrees with the direct integral.
        let f = earth.annual_sunlight(CubeSphere::new(32));
        let g = f.height.grid();
        let (c, _, _) = g.locate([1.0, 0.0, 0.0]);
        assert!((f.height.get(c) as f64 - equator).abs() < 0.01);
    }

    /// F234: no resonance for a fast spin or a spin that drifts against its orbit; the low
    /// orders tides capture.
    #[test]
    fn resonances_are_exact_and_of_low_order() {
        assert_eq!(spin_orbit_resonance(1.0), Some((1, 1)));
        assert_eq!(spin_orbit_resonance(1.499_997), Some((3, 2)));
        assert_eq!(spin_orbit_resonance(2.5), Some((5, 2)));
        assert_eq!(spin_orbit_resonance(366.256), None);
        assert_eq!(spin_orbit_resonance(0.9246), None);
        assert_eq!(spin_orbit_resonance(1.51), None);
    }
}
