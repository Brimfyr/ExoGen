// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! Crater production functions, chronology, saturation and population sampling
//! (PLAN.md 4.5, docs/parameters.md F31).
//!
//! Counts are cumulative per km² with diameters in km, following the Neukum convention.
//! The crater record begins [`CHRONOLOGY_CLIP_GYR`] before present: the chronology is flat
//! earlier, because its exponential term is not calibrated there and crusts were still
//! forming. With this start the formed basin counts of a Moon-like body match the Moon's.
//!
//! Saturation follows the standard picture: new craters obliterate old ones, so a surface
//! approaches an equilibrium density that is a fixed fraction of geometric saturation. Per
//! size bin the surviving density is `n_eq (1 - exp(-n_prod / n_eq))`, and the formation
//! times of survivors follow a truncated exponential in accumulated production, which is
//! sampled directly so that old surfaces carry mostly young survivors at small sizes and
//! old craters at large sizes, as observed.

use rayon::prelude::*;

use crate::constants::*;
use crate::fields::Field;
use crate::grid::CubeSphere;
use crate::rng::Rng;

/// A production function: polynomial size distribution plus chronology coefficients.
#[derive(Clone, Copy, Debug)]
pub struct ProductionFunction {
    pub coeffs: &'static [f64; 12],
    /// N(>= 1 km) = a (exp(b T) - 1) + c T per km², T in Gyr.
    pub chronology: (f64, f64, f64),
}

impl ProductionFunction {
    pub fn lunar() -> Self {
        Self {
            coeffs: &NPF_LUNAR,
            chronology: (NEUKUM_A, NEUKUM_B, NEUKUM_C),
        }
    }

    pub fn mars() -> Self {
        Self {
            coeffs: &NPF_MARS,
            chronology: MARS_CHRONOLOGY,
        }
    }

    pub fn by_name(name: &str) -> Self {
        match name {
            "neukum_mars" => Self::mars(),
            _ => Self::lunar(),
        }
    }

    /// Cumulative density N(>= D) per km² for a 1 Gyr surface, D in km. Outside the fitted
    /// range the polynomial is continued with the slope at the range edge.
    pub fn cumulative_1gyr(&self, d_km: f64) -> f64 {
        let (lo, hi) = (NPF_RANGE_KM.0, NPF_RANGE_KM.1);
        let x = d_km.clamp(lo, hi).log10();
        let mut acc = 0.0;
        for (j, a) in self.coeffs.iter().enumerate() {
            acc += a * x.powi(j as i32);
        }
        if d_km < lo || d_km > hi {
            let slope = self.log_slope(d_km.clamp(lo, hi));
            acc += slope * (d_km.log10() - x);
        }
        10f64.powf(acc)
    }

    /// d log N / d log D at a diameter.
    pub fn log_slope(&self, d_km: f64) -> f64 {
        let x = d_km.clamp(NPF_RANGE_KM.0, NPF_RANGE_KM.1).log10();
        let mut s = 0.0;
        for (j, a) in self.coeffs.iter().enumerate().skip(1) {
            s += a * j as f64 * x.powi(j as i32 - 1);
        }
        s
    }

    /// Chronology: N(>= 1 km) per km² for a surface of age `t_gyr`. The crater record
    /// begins [`CHRONOLOGY_CLIP_GYR`] before present: the function is flat beyond it, so
    /// nothing formed earlier is counted. The chronology is not calibrated there and crusts
    /// were still forming; holding the early rate constant instead triples the basin count
    /// of a Moon-like body.
    pub fn n1(&self, t_gyr: f64) -> f64 {
        let (a, b, c) = self.chronology;
        let t = t_gyr.clamp(0.0, CHRONOLOGY_CLIP_GYR);
        a * ((b * t).exp() - 1.0) + c * t
    }

    /// Production rate of craters >= 1 km per km² per Gyr at `t_gyr` before present; zero
    /// before the start of the record.
    pub fn rate1(&self, t_gyr: f64) -> f64 {
        if !(0.0..=CHRONOLOGY_CLIP_GYR).contains(&t_gyr) {
            return 0.0;
        }
        let (a, b, c) = self.chronology;
        a * b * (b * t_gyr).exp() + c
    }

    /// Scale factor turning the 1 Gyr density into the density for age `t_gyr`.
    pub fn age_factor(&self, t_gyr: f64) -> f64 {
        self.n1(t_gyr) / self.n1(1.0)
    }

    /// Inverse of the chronology: the age at which N(>= 1 km) equals `n`.
    pub fn age_for_n1(&self, n: f64) -> f64 {
        if n <= 0.0 {
            return 0.0;
        }
        let (mut lo, mut hi) = (0.0f64, 20.0f64);
        if self.n1(hi) < n {
            return hi;
        }
        for _ in 0..60 {
            let mid = 0.5 * (lo + hi);
            if self.n1(mid) < n {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    }

    /// Cumulative density N(>= D) per km² for a surface of age `t_gyr`.
    pub fn cumulative(&self, d_km: f64, t_gyr: f64) -> f64 {
        self.cumulative_1gyr(d_km) * self.age_factor(t_gyr)
    }

    /// Sample a diameter (km) above `d_min_km` from the cumulative distribution, capped at
    /// `d_max_km`, by bisection on the cumulative density.
    pub fn sample_diameter_above(&self, u: f64, d_min_km: f64, d_max_km: f64) -> f64 {
        let n_min = self.cumulative_1gyr(d_min_km);
        let n_max = self.cumulative_1gyr(d_max_km);
        // Target cumulative value between n_max and n_min; larger u gives a smaller crater.
        let target = n_max + (n_min - n_max) * u.clamp(0.0, 1.0);
        let (mut lo, mut hi) = (d_min_km.ln(), d_max_km.ln());
        for _ in 0..50 {
            let mid = 0.5 * (lo + hi);
            if self.cumulative_1gyr(mid.exp()) > target {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        (0.5 * (lo + hi)).exp()
    }
}

/// Cumulative geometric saturation N(>= D) per km², D in km.
pub fn geometric_saturation(d_km: f64) -> f64 {
    GEOMETRIC_SATURATION / (d_km * d_km)
}

/// Equilibrium (empirical saturation) density of a size bin per km².
pub fn bin_equilibrium(d_lo_km: f64, d_hi_km: f64, fraction: f64) -> f64 {
    fraction * (geometric_saturation(d_lo_km) - geometric_saturation(d_hi_km))
}

/// Surviving density given production density and equilibrium density.
pub fn survivors(n_prod: f64, n_eq: f64) -> f64 {
    if n_eq <= 0.0 {
        return 0.0;
    }
    n_eq * (1.0 - (-n_prod / n_eq).exp())
}

/// Knuth's Poisson sampler for small means; a normal approximation above 40.
pub fn poisson(rng: &Rng, index: u64, draw_base: u64, lambda: f64) -> (u32, u64) {
    if lambda <= 0.0 {
        return (0, draw_base);
    }
    if lambda > 40.0 {
        let z = rng.normal(index, draw_base);
        return (
            ((lambda + z * lambda.sqrt()).round().max(0.0)) as u32,
            draw_base + 2,
        );
    }
    let l = (-lambda).exp();
    let mut k = 0u32;
    let mut p = 1.0;
    let mut draw = draw_base;
    loop {
        p *= rng.uniform(index, draw);
        draw += 1;
        if p <= l {
            break;
        }
        k += 1;
        if k > 1000 {
            break;
        }
    }
    (k, draw)
}

/// A sampled crater, before stamping.
#[derive(Clone, Debug)]
pub struct CraterSample {
    pub center: [f64; 3],
    pub diameter: f64,
    /// Formation time, Gyr since the body formed.
    pub time_gyr: f64,
    pub azimuth: f64,
    pub ellipticity: f64,
    pub key: u64,
}

/// Inputs to the statistical population sampler.
#[derive(Clone, Debug)]
pub struct PopulationConfig {
    pub pf: ProductionFunction,
    /// Impact flux relative to the lunar production function, including the density
    /// multiplier.
    pub flux_scale: f64,
    pub body_radius: f64,
    /// Smallest and largest diameters sampled here, metres.
    pub d_min: f64,
    pub d_max: f64,
    pub bins_per_octave: u32,
    /// Fraction of geometric saturation at equilibrium.
    pub saturation_fraction: f64,
    /// Longest effective surface age for crater retention, Gyr.
    pub retention_age: f64,
    pub oblique_fraction: f64,
    pub seed: u64,
}

/// Sample the surviving sub-threshold crater population at `present` (Gyr since formation)
/// given the time each cell was last resurfaced. Deterministic per cell and bin.
pub fn sample_population(
    grid: &CubeSphere,
    surface_age: &Field<f32>,
    present: f64,
    cfg: &PopulationConfig,
) -> Vec<CraterSample> {
    if cfg.d_max <= cfg.d_min {
        return Vec::new();
    }
    let octaves = (cfg.d_max / cfg.d_min).log2();
    let nbins = ((octaves * cfg.bins_per_octave as f64).ceil() as usize).max(1);
    let ratio = (cfg.d_max / cfg.d_min).powf(1.0 / nbins as f64);
    let bins: Vec<(f64, f64, f64)> = (0..nbins)
        .map(|b| {
            let lo = cfg.d_min * ratio.powi(b as i32);
            (lo, lo * ratio, (lo * lo * ratio).sqrt())
        })
        .collect();
    // Per-bin production density per km² for a 1 Gyr surface, and equilibrium density.
    let bin_prod_1gyr: Vec<f64> = bins
        .iter()
        .map(|(lo, hi, _)| {
            cfg.flux_scale
                * (cfg.pf.cumulative_1gyr(lo / 1e3) - cfg.pf.cumulative_1gyr(hi / 1e3)).max(0.0)
        })
        .collect();
    let bin_eq: Vec<f64> = bins
        .iter()
        .map(|(lo, hi, _)| bin_equilibrium(lo / 1e3, hi / 1e3, cfg.saturation_fraction))
        .collect();
    let n1_1gyr = cfg.pf.n1(1.0);
    let rng = Rng::new(cfg.seed);
    let delta = grid.delta();
    let r2_km2 = cfg.body_radius * cfg.body_radius / 1e6;

    let per_cell: Vec<Vec<CraterSample>> = (0..grid.cell_count())
        .into_par_iter()
        .map(|k| {
            let cell = grid.cell(k);
            let area_km2 = grid.cell_area(cell) * r2_km2;
            let age = (present - surface_age.as_slice()[k] as f64).clamp(0.0, cfg.retention_age);
            if age <= 0.0 {
                return Vec::new();
            }
            let age_factor = cfg.pf.age_factor(age);
            let (a0, b0) = grid.angles(cell);
            let mut out = Vec::new();
            for (b, (lo, hi, _mid)) in bins.iter().enumerate() {
                let n_prod = bin_prod_1gyr[b] * age_factor;
                let n_eq = bin_eq[b];
                let lambda = survivors(n_prod, n_eq) * area_km2;
                if lambda < 1e-9 {
                    continue;
                }
                let draw_base = (b as u64) << 20;
                let (count, mut draw) = poisson(&rng, k as u64, draw_base, lambda);
                for _ in 0..count {
                    let u_pos_a = rng.uniform(k as u64, draw);
                    let u_pos_b = rng.uniform(k as u64, draw + 1);
                    let u_d = rng.uniform(k as u64, draw + 2);
                    let u_t = rng.uniform(k as u64, draw + 3);
                    let u_az = rng.uniform(k as u64, draw + 4);
                    let u_ob = rng.uniform(k as u64, draw + 5);
                    draw += 6;
                    let center = CubeSphere::direction(
                        cell.face,
                        a0 - 0.5 * delta + u_pos_a * delta,
                        b0 - 0.5 * delta + u_pos_b * delta,
                    );
                    // Diameter within the bin from the local power-law slope.
                    let slope = -cfg.pf.log_slope((lo * hi).sqrt() / 1e3).max(0.5);
                    let q = (hi / lo).powf(-slope);
                    let diameter = lo * (1.0 - u_d * (1.0 - q)).powf(-1.0 / slope);
                    // Accumulated production after formation, truncated exponential.
                    let f = if n_eq > 0.0 {
                        -n_eq * (1.0 - u_t * (1.0 - (-n_prod / n_eq).exp())).ln()
                    } else {
                        u_t * n_prod
                    };
                    let after_factor = if bin_prod_1gyr[b] > 0.0 {
                        f / bin_prod_1gyr[b]
                    } else {
                        0.0
                    };
                    let age_before_present = cfg.pf.age_for_n1(after_factor * n1_1gyr).min(age);
                    let time_gyr = present - age_before_present;
                    let (azimuth, ellipticity) = if u_ob < cfg.oblique_fraction {
                        (
                            u_az * 2.0 * std::f64::consts::PI,
                            1.15 + 0.35 * rng.uniform(k as u64, draw + 6),
                        )
                    } else {
                        (0.0, 1.0)
                    };
                    out.push(CraterSample {
                        center,
                        diameter,
                        time_gyr,
                        azimuth,
                        ellipticity,
                        key: (k as u64) << 24 | ((b as u64) << 12) | (out.len() as u64),
                    });
                }
            }
            out
        })
        .collect();
    let mut all: Vec<CraterSample> = per_cell.into_iter().flatten().collect();
    // Oldest first; ties broken by key so the order is deterministic.
    all.sort_by(|a, b| a.time_gyr.total_cmp(&b.time_gyr).then(a.key.cmp(&b.key)));
    all
}

/// Expected number of craters at or above `d_min` (metres) over the whole body formed
/// between `t0` and `t1` Gyr since formation, for a body of age `age` and the given flux.
pub fn expected_count_above(
    pf: &ProductionFunction,
    flux_scale: f64,
    body_radius: f64,
    d_min: f64,
    age: f64,
    t0: f64,
    t1: f64,
) -> f64 {
    let area_km2 = 4.0 * std::f64::consts::PI * body_radius * body_radius / 1e6;
    let tau0 = (age - t0).max(0.0);
    let tau1 = (age - t1).max(0.0);
    let per_km2 = pf.cumulative_1gyr(d_min / 1e3) * (pf.n1(tau0) - pf.n1(tau1)) / pf.n1(1.0);
    (area_km2 * flux_scale * per_km2).max(0.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn chronology_matches_the_spike_values_and_is_monotonic() {
        let pf = ProductionFunction::lunar();
        assert!((pf.n1(1.0) - 8.38e-4).abs() < 1e-6);
        assert!((pf.n1(3.3) - 3.23e-3).abs() < 1e-4);
        assert!((pf.n1(4.0) - 6.28e-2).abs() < 2e-3);
        let mut last = 0.0;
        for i in 0..=84 {
            let t = i as f64 * 0.05;
            let n = pf.n1(t);
            assert!(n >= last);
            last = n;
            assert!((pf.age_for_n1(n) - t).abs() < 1e-6, "inverse at {t}");
        }
        // Flat before the start of the record.
        assert_eq!(pf.n1(4.5), pf.n1(CHRONOLOGY_CLIP_GYR));
        assert_eq!(pf.rate1(4.5), 0.0);
        assert!((pf.cumulative_1gyr(1.0) - 8.17e-4).abs() < 1e-5);
    }

    #[test]
    fn diameter_sampling_follows_the_cumulative_distribution() {
        let pf = ProductionFunction::lunar();
        let rng = Rng::new(3);
        let n = 20_000;
        let d_min = 20.0;
        let mut above_50 = 0;
        for i in 0..n {
            let d = pf.sample_diameter_above(rng.uniform(i, 0), d_min, 2000.0);
            assert!(d >= d_min - 1e-9);
            if d >= 50.0 {
                above_50 += 1;
            }
        }
        let expected = pf.cumulative_1gyr(50.0) / pf.cumulative_1gyr(d_min);
        let got = above_50 as f64 / n as f64;
        assert!(
            (got - expected).abs() < 0.02,
            "got {got} expected {expected}"
        );
    }

    #[test]
    fn population_matches_the_isochron_on_a_young_surface_and_saturates_on_an_old_one() {
        let grid = CubeSphere::new(64);
        let radius = 1.7374e6;
        let cfg = PopulationConfig {
            pf: ProductionFunction::lunar(),
            flux_scale: 1.0,
            body_radius: radius,
            d_min: 40_000.0,
            d_max: 320_000.0,
            bins_per_octave: 4,
            saturation_fraction: EMPIRICAL_SATURATION_FRACTION,
            retention_age: 10.0,
            oblique_fraction: 0.0,
            seed: 11,
        };
        let area_km2 = 4.0 * std::f64::consts::PI * radius * radius / 1e6;
        for (surface_time, present, check_saturation) in [(1.2f64, 4.5f64, false), (0.0, 4.5, true)]
        {
            let ages = Field::filled(grid, surface_time as f32);
            let pop = sample_population(&grid, &ages, present, &cfg);
            let age = present - surface_time;
            for d in [40.0, 80.0, 160.0] {
                let count = pop.iter().filter(|c| c.diameter >= d * 1e3).count() as f64;
                let prod =
                    cfg.pf.cumulative(d, age) * area_km2 - cfg.pf.cumulative(320.0, age) * area_km2;
                let sat = EMPIRICAL_SATURATION_FRACTION
                    * (geometric_saturation(d) - geometric_saturation(320.0))
                    * area_km2;
                let expected = if check_saturation {
                    sat.min(prod)
                } else {
                    prod
                };
                let tol = 0.3 * expected + 3.0 * expected.sqrt();
                assert!(
                    (count - expected).abs() < tol,
                    "D>={d} km: count {count} expected {expected} (prod {prod}, sat {sat})"
                );
            }
            // Ages must lie within the surface age and be oldest-first.
            assert!(pop.windows(2).all(|w| w[0].time_gyr <= w[1].time_gyr));
            assert!(
                pop.iter()
                    .all(|c| c.time_gyr >= surface_time - 1e-9 && c.time_gyr <= present + 1e-9)
            );
        }
    }
}
