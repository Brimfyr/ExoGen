// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! Real spherical harmonics up to a modest degree, for the low-degree primordial crust
//! field (PLAN.md 4.1) and for tests of grid operators.

use std::f64::consts::PI;

use crate::rng::Rng;

/// Associated Legendre functions P_l^m(x) for all l <= lmax and 0 <= m <= l, without the
/// Condon–Shortley phase, by the standard recurrences. Returned as a flat table indexed by
/// `index(l, m)`.
fn legendre_table(lmax: usize, x: f64) -> Vec<f64> {
    let n = (lmax + 1) * (lmax + 2) / 2;
    let mut p = vec![0.0; n];
    let idx = |l: usize, m: usize| l * (l + 1) / 2 + m;
    let somx2 = (1.0 - x * x).max(0.0).sqrt();
    p[idx(0, 0)] = 1.0;
    for m in 1..=lmax {
        p[idx(m, m)] = p[idx(m - 1, m - 1)] * (2 * m - 1) as f64 * somx2;
    }
    for m in 0..lmax {
        p[idx(m + 1, m)] = x * (2 * m + 1) as f64 * p[idx(m, m)];
    }
    for m in 0..=lmax {
        for l in (m + 2)..=lmax {
            p[idx(l, m)] = (x * (2 * l - 1) as f64 * p[idx(l - 1, m)]
                - (l + m - 1) as f64 * p[idx(l - 2, m)])
                / (l - m) as f64;
        }
    }
    p
}

fn factorial_ratio(l: usize, m: usize) -> f64 {
    // (l - m)! / (l + m)!
    let mut r = 1.0;
    for k in (l - m + 1)..=(l + m) {
        r /= k as f64;
    }
    r
}

/// A real spherical-harmonic expansion: coefficients for (l, m) with m in -l..=l, where
/// negative m denotes the sine term.
#[derive(Clone, Debug)]
pub struct Harmonics {
    pub lmax: usize,
    /// Indexed by `l² + (m + l)` for m in -l..=l: degrees below l occupy the first l² slots.
    pub coeffs: Vec<f64>,
}

impl Harmonics {
    pub fn zeros(lmax: usize) -> Self {
        Self {
            lmax,
            coeffs: vec![0.0; (lmax + 1) * (lmax + 1)],
        }
    }

    #[inline]
    pub fn index(l: usize, m: i64) -> usize {
        l * l + (m + l as i64) as usize
    }

    pub fn set(&mut self, l: usize, m: i64, v: f64) {
        self.coeffs[Self::index(l, m)] = v;
    }

    /// Random coefficients with a power spectrum proportional to `(l (l + 1))^(-slope / 2)`
    /// for l from `lmin` to `lmax`, unit variance per degree before weighting.
    pub fn random(rng: &Rng, lmin: usize, lmax: usize, slope: f64) -> Self {
        let mut h = Self::zeros(lmax);
        let mut k = 0u64;
        for l in lmin.max(1)..=lmax {
            let weight = ((l * (l + 1)) as f64).powf(-0.5 * slope) / (2 * l + 1) as f64;
            for m in -(l as i64)..=(l as i64) {
                h.set(l, m, rng.normal(k, 0) * weight.sqrt());
                k += 1;
            }
        }
        h
    }

    /// Evaluate at a unit direction.
    pub fn eval(&self, d: [f64; 3]) -> f64 {
        let cos_theta = d[2].clamp(-1.0, 1.0);
        let phi = d[1].atan2(d[0]);
        let p = legendre_table(self.lmax, cos_theta);
        let idx = |l: usize, m: usize| l * (l + 1) / 2 + m;
        let mut sum = 0.0;
        for l in 0..=self.lmax {
            for m in 0..=l {
                let norm = ((2 * l + 1) as f64 / (4.0 * PI) * factorial_ratio(l, m)).sqrt();
                let plm = p[idx(l, m)] * norm;
                if m == 0 {
                    sum += self.coeffs[Self::index(l, 0)] * plm;
                } else {
                    let s2 = 2f64.sqrt();
                    sum +=
                        self.coeffs[Self::index(l, m as i64)] * s2 * plm * (m as f64 * phi).cos();
                    sum += self.coeffs[Self::index(l, -(m as i64))]
                        * s2
                        * plm
                        * (m as f64 * phi).sin();
                }
            }
        }
        sum
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn low_degree_harmonics_match_closed_forms() {
        // Y_1^0 = sqrt(3 / 4 pi) z ; Y_2^0 = sqrt(5 / 16 pi) (3 z² - 1)
        let mut h = Harmonics::zeros(2);
        h.set(1, 0, 1.0);
        let d = [0.3, 0.4, (1.0f64 - 0.25).sqrt()];
        assert!((h.eval(d) - (3.0 / (4.0 * PI)).sqrt() * d[2]).abs() < 1e-12);
        let mut h2 = Harmonics::zeros(2);
        h2.set(2, 0, 1.0);
        assert!(
            (h2.eval(d) - (5.0 / (16.0 * PI)).sqrt() * (3.0 * d[2] * d[2] - 1.0)).abs() < 1e-12
        );
        // Real Y_1^1 (cos) = sqrt(3 / 4 pi) x
        let mut h3 = Harmonics::zeros(1);
        h3.set(1, 1, 1.0);
        assert!((h3.eval(d) - (3.0 / (4.0 * PI)).sqrt() * d[0]).abs() < 1e-12);
    }

    #[test]
    fn random_field_has_the_expected_spread() {
        let h = Harmonics::random(&Rng::new(5), 1, 8, 2.0);
        let rng = Rng::new(9);
        let (mut sum, mut sq) = (0.0, 0.0);
        let n = 5000u64;
        for i in 0..n {
            let v = h.eval(rng.unit_direction(i, 0));
            sum += v;
            sq += v * v;
        }
        let mean = sum / n as f64;
        let rms = (sq / n as f64).sqrt();
        assert!(mean.abs() < 0.1, "{mean}");
        assert!(rms > 0.05 && rms < 2.0, "{rms}");
    }
}
