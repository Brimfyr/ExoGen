// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! Crater production functions, flux history, saturation and population sampling
//! (PLAN.md 4.5, docs/parameters.md F31).
//!
//! Counts are cumulative per km² with diameters in km, following the Neukum convention.
//! The flux history is written in time since the body formed: a steady rate plus an early
//! component decaying from formation, which is the Neukum chronology rewritten so that a
//! body of any age gets its heavy bombardment at formation and the steady flux ever after.
//! Nothing is counted before the surface solidifies.
//!
//! Saturation follows the standard picture: new craters obliterate old ones, so a surface
//! approaches an equilibrium density that is a fixed fraction of geometric saturation. Per
//! size bin the surviving density is `n_eq (1 - exp(-n_prod / n_eq))`, and the formation
//! times of survivors follow a truncated exponential in accumulated production, which is
//! sampled directly so that old surfaces carry mostly young survivors at small sizes and
//! old craters at large sizes, as observed.

use rayon::prelude::*;

use super::morphology::{ellipticity_for, impact_angle};
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

    /// The flux history of this production function's chronology in time since formation,
    /// for a surface that starts retaining craters `record_start` Gyr after formation.
    pub fn flux_history(&self, record_start: f64) -> FluxHistory {
        FluxHistory::from_chronology(self.chronology, record_start)
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

/// Impact flux history of a body in time since its formation (docs/parameters.md F31).
///
/// The rate of craters at or above 1 km per km² per Gyr is a steady rate plus an early
/// component that decays exponentially from formation. With the Neukum coefficients it is
/// the chronology's derivative written in time since formation, the present being
/// [`NEUKUM_REFERENCE_AGE_GYR`] after formation. Nothing is counted before the surface
/// solidifies at `record_start`.
#[derive(Clone, Copy, Debug)]
pub struct FluxHistory {
    /// Steady rate of craters >= 1 km per km² per Gyr at the reference flux.
    pub steady: f64,
    /// Early rate at formation relative to the steady rate.
    pub early_ratio: f64,
    /// E-folding time of the early component, Gyr.
    pub early_decay: f64,
    /// Time since formation when the surface begins retaining craters, Gyr.
    pub record_start: f64,
}

impl FluxHistory {
    /// From chronology coefficients (a, b, c) of N1(T) = a (exp(b T) - 1) + c T, with T in
    /// Gyr before present.
    pub fn from_chronology(chronology: (f64, f64, f64), record_start: f64) -> Self {
        let (_, b, c) = chronology;
        Self {
            steady: c,
            early_ratio: early_flux_ratio(chronology),
            early_decay: 1.0 / b,
            record_start,
        }
    }

    /// Rate of craters >= 1 km per km² per Gyr at `tau` Gyr since formation.
    pub fn rate1(&self, tau: f64) -> f64 {
        if tau < self.record_start {
            return 0.0;
        }
        self.steady * (1.0 + self.early_ratio * (-tau / self.early_decay).exp())
    }

    /// N(>= 1 km) per km² formed between `tau0` and `tau1` Gyr since formation.
    pub fn n1_between(&self, tau0: f64, tau1: f64) -> f64 {
        let a = tau0.max(self.record_start);
        let b = tau1.max(a);
        let early = self.early_ratio
            * self.early_decay
            * ((-a / self.early_decay).exp() - (-b / self.early_decay).exp());
        self.steady * ((b - a) + early)
    }

    /// Gyr of steady flux equivalent to the production between `tau0` and `tau1`.
    pub fn equivalent_gyr(&self, tau0: f64, tau1: f64) -> f64 {
        if self.steady <= 0.0 {
            return 0.0;
        }
        self.n1_between(tau0, tau1) / self.steady
    }

    /// The time within `tau0` to `tau1` by which the fraction `u` of that interval's
    /// production has arrived, so draws follow the rate within a step.
    pub fn time_at_fraction(&self, tau0: f64, tau1: f64, u: f64) -> f64 {
        let a = tau0.max(self.record_start).min(tau1);
        let total = self.equivalent_gyr(a, tau1);
        if total <= 0.0 {
            return tau1;
        }
        let target = u.clamp(0.0, 1.0) * total;
        let (mut lo, mut hi) = (a, tau1);
        for _ in 0..50 {
            let mid = 0.5 * (lo + hi);
            if self.equivalent_gyr(a, mid) < target {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    }

    /// The time since formation, before `present`, after which `remaining` Gyr of
    /// equivalent flux still arrive; the record start when that is the whole record.
    pub fn time_with_remaining(&self, present: f64, remaining: f64) -> f64 {
        let (mut lo, mut hi) = (self.record_start.min(present), present);
        if remaining >= self.equivalent_gyr(lo, present) {
            return lo;
        }
        if remaining <= 0.0 {
            return hi;
        }
        for _ in 0..60 {
            let mid = 0.5 * (lo + hi);
            if self.equivalent_gyr(mid, present) > remaining {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    }
}

/// Early rate at formation relative to the steady rate, for chronology coefficients.
pub fn early_flux_ratio(chronology: (f64, f64, f64)) -> f64 {
    let (a, b, c) = chronology;
    a * b * (b * NEUKUM_REFERENCE_AGE_GYR).exp() / c
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
    /// Downrange direction, radians in the local east-north frame.
    pub azimuth: f64,
    /// Impact angle from horizontal, radians.
    pub angle: f64,
    pub ellipticity: f64,
    pub key: u64,
}

/// Inputs to the statistical population sampler.
#[derive(Clone, Debug)]
pub struct PopulationConfig {
    pub pf: ProductionFunction,
    pub flux: FluxHistory,
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
    let rng = Rng::new(cfg.seed);
    let delta = grid.delta();
    let r2_km2 = cfg.body_radius * cfg.body_radius / 1e6;

    let per_cell: Vec<Vec<CraterSample>> = (0..grid.cell_count())
        .into_par_iter()
        .map(|k| {
            let cell = grid.cell(k);
            let area_km2 = grid.cell_area(cell) * r2_km2;
            // The surface records craters since it formed, or since the retention age.
            let start = (surface_age.as_slice()[k] as f64).max(present - cfg.retention_age);
            let exposure = cfg.flux.equivalent_gyr(start, present);
            if exposure <= 0.0 {
                return Vec::new();
            }
            let (a0, b0) = grid.angles(cell);
            let mut out = Vec::new();
            for (b, (lo, hi, _mid)) in bins.iter().enumerate() {
                let n_prod = bin_prod_1gyr[b] * exposure;
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
                    // Gyr of steady-flux equivalent still to come after formation.
                    let remaining = if bin_prod_1gyr[b] > 0.0 {
                        f / bin_prod_1gyr[b]
                    } else {
                        0.0
                    };
                    let time_gyr = cfg.flux.time_with_remaining(present, remaining).max(start);
                    let angle = impact_angle(u_ob);
                    let azimuth = u_az * 2.0 * std::f64::consts::PI;
                    let ellipticity = ellipticity_for(angle, cfg.oblique_fraction);
                    out.push(CraterSample {
                        center,
                        diameter,
                        time_gyr,
                        azimuth,
                        angle,
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
/// between `t0` and `t1` Gyr since formation.
pub fn expected_count_above(
    pf: &ProductionFunction,
    flux: &FluxHistory,
    flux_scale: f64,
    body_radius: f64,
    d_min: f64,
    t0: f64,
    t1: f64,
) -> f64 {
    let area_km2 = 4.0 * std::f64::consts::PI * body_radius * body_radius / 1e6;
    let per_km2 = pf.cumulative_1gyr(d_min / 1e3) * flux.equivalent_gyr(t0, t1);
    (area_km2 * flux_scale * per_km2).max(0.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn flux_history_matches_the_neukum_chronology() {
        let pf = ProductionFunction::lunar();
        let flux = pf.flux_history(0.3);
        let present = NEUKUM_REFERENCE_AGE_GYR;
        // The whole record is the chronology's value 4.2 Gyr before present.
        let n_42 = NEUKUM_A * ((NEUKUM_B * 4.2).exp() - 1.0) + NEUKUM_C * 4.2;
        let whole = flux.n1_between(0.0, present);
        assert!((whole / n_42 - 1.0).abs() < 1e-6, "{whole} vs {n_42}");
        // The last 1 Gyr and the last 3.3 Gyr match the spike values.
        assert!((flux.n1_between(present - 1.0, present) - 8.38e-4).abs() < 1e-6);
        assert!((flux.n1_between(present - 3.3, present) - 3.23e-3).abs() < 1e-4);
        // Nothing before the record starts; a heavy flux right after it.
        assert_eq!(flux.rate1(0.2), 0.0);
        assert!(flux.rate1(0.3) > 1_000.0 * flux.rate1(present));
        // The inverse of the remaining production round-trips.
        for t in [0.3, 0.5, 1.0, 2.5, 4.0, 4.4] {
            let remaining = flux.equivalent_gyr(t, present);
            let back = flux.time_with_remaining(present, remaining);
            assert!((back - t).abs() < 1e-6, "inverse at {t}: {back}");
        }
        assert!((pf.cumulative_1gyr(1.0) - 8.17e-4).abs() < 1e-5);
    }

    /// A body of any age gets the early bombardment at formation and the steady flux after
    /// it: the record of a 10 Gyr body is that of a 6 Gyr body plus 4 Gyr of steady flux.
    #[test]
    fn bodies_of_any_age_get_the_early_bombardment_at_formation() {
        let flux = ProductionFunction::lunar().flux_history(0.3);
        let n6 = flux.n1_between(0.0, 6.0);
        let n10 = flux.n1_between(0.0, 10.0);
        assert!(((n10 - n6) / (4.0 * flux.steady) - 1.0).abs() < 1e-9);
        // The early phase holds most of a 4.5 Gyr record and a smaller share of a 10 Gyr one.
        let early = flux.n1_between(0.0, 1.5);
        assert!(early / flux.n1_between(0.0, 4.5) > 0.98);
        assert!(early / n10 < early / flux.n1_between(0.0, 4.5));
        // The early rate at formation is the chronology's exponential term.
        let r = early_flux_ratio((NEUKUM_A, NEUKUM_B, NEUKUM_C));
        assert!((r - 1.57e4).abs() < 2e2, "{r}");
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
            flux: ProductionFunction::lunar().flux_history(0.3),
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
            let exposure = cfg.flux.equivalent_gyr(surface_time, present);
            for d in [40.0, 80.0, 160.0] {
                let count = pop.iter().filter(|c| c.diameter >= d * 1e3).count() as f64;
                let prod = (cfg.pf.cumulative_1gyr(d) - cfg.pf.cumulative_1gyr(320.0))
                    * exposure
                    * area_km2;
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
