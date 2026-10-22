// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! Real spherical harmonics up to a modest degree, for the low-degree primordial crust
//! field (PLAN.md 4.1), the mantle's convection (F173) and for tests of grid operators.

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

    /// Random coefficients whose power per degree is `power(l)` for l from 1 to `lmax`, the
    /// coefficients of a degree sharing it, drawn from the stream `draw` (each draw takes
    /// two of the generator's draws).
    pub fn random_spectrum(
        rng: &Rng,
        draw: u64,
        lmax: usize,
        power: impl Fn(usize) -> f64,
    ) -> Self {
        let mut h = Self::zeros(lmax);
        let mut k = 0u64;
        for l in 1..=lmax {
            let weight = power(l).max(0.0) / (2 * l + 1) as f64;
            for m in -(l as i64)..=(l as i64) {
                h.set(l, m, rng.normal(k, 2 * draw) * weight.sqrt());
                k += 1;
            }
        }
        h
    }

    /// The field's root mean square over the sphere: the coefficients' sum of squares over
    /// 4π, the harmonics being orthonormal.
    pub fn rms(&self) -> f64 {
        (self.coeffs.iter().map(|c| c * c).sum::<f64>() / (4.0 * PI)).sqrt()
    }

    /// Each degree's coefficients times `factor(l)`.
    pub fn scale_degrees(&mut self, factor: impl Fn(usize) -> f64) {
        for l in 0..=self.lmax {
            let f = factor(l);
            for m in -(l as i64)..=(l as i64) {
                self.coeffs[Self::index(l, m)] *= f;
            }
        }
    }

    /// The field on a grid of `nlat` latitude rows and `nlon` longitude columns, a row's
    /// Legendre functions summed for each order and the orders then over the row, so that
    /// the cost is the degrees squared per row rather than per point.
    pub fn grid(&self, nlat: usize, nlon: usize) -> LatLon {
        let lmax = self.lmax;
        let idx = |l: usize, m: usize| l * (l + 1) / 2 + m;
        let mut norm = vec![0.0; (lmax + 1) * (lmax + 2) / 2];
        for l in 0..=lmax {
            for m in 0..=l {
                let s = if m == 0 { 1.0 } else { 2f64.sqrt() };
                norm[idx(l, m)] =
                    s * ((2 * l + 1) as f64 / (4.0 * PI) * factorial_ratio(l, m)).sqrt();
            }
        }
        let mut cs = vec![0.0; (lmax + 1) * nlon];
        let mut sn = vec![0.0; (lmax + 1) * nlon];
        for m in 0..=lmax {
            for j in 0..nlon {
                let phi = LatLon::longitude(j, nlon);
                cs[m * nlon + j] = (m as f64 * phi).cos();
                sn[m * nlon + j] = (m as f64 * phi).sin();
            }
        }
        let mut values = vec![0.0f32; nlat * nlon];
        let (mut a, mut b) = (vec![0.0; lmax + 1], vec![0.0; lmax + 1]);
        for i in 0..nlat {
            let p = legendre_table(lmax, LatLon::colatitude(i, nlat).cos());
            for m in 0..=lmax {
                let (mut sa, mut sb) = (0.0, 0.0);
                for l in m..=lmax {
                    let plm = p[idx(l, m)] * norm[idx(l, m)];
                    sa += self.coeffs[Self::index(l, m as i64)] * plm;
                    if m > 0 {
                        sb += self.coeffs[Self::index(l, -(m as i64))] * plm;
                    }
                }
                a[m] = sa;
                b[m] = sb;
            }
            for j in 0..nlon {
                let mut v = 0.0;
                for m in 0..=lmax {
                    v += a[m] * cs[m * nlon + j] + b[m] * sn[m * nlon + j];
                }
                values[i * nlon + j] = v as f32;
            }
        }
        LatLon { nlat, nlon, values }
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

/// A field on a grid of latitude rows and longitude columns, the nodes at colatitude
/// (i + 0.5) π / nlat and longitude −π + (j + 0.5) 2π / nlon.
#[derive(Clone, Debug)]
pub struct LatLon {
    pub nlat: usize,
    pub nlon: usize,
    pub values: Vec<f32>,
}

impl LatLon {
    pub fn colatitude(i: usize, nlat: usize) -> f64 {
        (i as f64 + 0.5) * PI / nlat as f64
    }

    pub fn longitude(j: usize, nlon: usize) -> f64 {
        -PI + (j as f64 + 0.5) * 2.0 * PI / nlon as f64
    }

    /// The unit direction of a node.
    pub fn direction(&self, i: usize, j: usize) -> [f64; 3] {
        let (st, ct) = Self::colatitude(i, self.nlat).sin_cos();
        let (sp, cp) = Self::longitude(j, self.nlon).sin_cos();
        [st * cp, st * sp, ct]
    }

    pub fn at(&self, i: usize, j: usize) -> f64 {
        self.values[i * self.nlon + j] as f64
    }

    /// The value at a unit direction, bilinear between the nodes about it, the rows beyond
    /// the first and last held at them.
    pub fn sample(&self, d: [f64; 3]) -> f64 {
        let theta = d[2].clamp(-1.0, 1.0).acos();
        let phi = d[1].atan2(d[0]);
        let last = (self.nlat - 1) as f64;
        let fi = (theta / PI * self.nlat as f64 - 0.5).clamp(0.0, last);
        let fj = (phi + PI) / (2.0 * PI) * self.nlon as f64 - 0.5;
        let i0 = (fi.floor() as usize).min(self.nlat.saturating_sub(2));
        let ti = (fi - i0 as f64).clamp(0.0, 1.0);
        let i1 = (i0 + 1).min(self.nlat - 1);
        let j0f = fj.floor();
        let tj = fj - j0f;
        let j0 = (j0f as i64).rem_euclid(self.nlon as i64) as usize;
        let j1 = (j0 + 1) % self.nlon;
        let row = |i: usize| (1.0 - tj) * self.at(i, j0) + tj * self.at(i, j1);
        (1.0 - ti) * row(i0) + ti * row(i1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_grid_holds_the_harmonics_values_at_its_nodes_and_between_them() {
        let h = Harmonics::random_spectrum(&Rng::new(3), 1, 12, |l| 1.0 / l as f64);
        let g = h.grid(48, 96);
        let mut worst = 0.0f64;
        for (i, j) in [(0, 0), (5, 17), (23, 40), (47, 95), (30, 3)] {
            worst = worst.max((g.at(i, j) - h.eval(g.direction(i, j))).abs());
        }
        assert!(worst < 1e-5 * h.rms().max(1.0), "{worst}");
        let rng = Rng::new(11);
        let mut far = 0.0f64;
        for n in 0..200u64 {
            let d = rng.unit_direction(n, 0);
            far = far.max((g.sample(d) - h.eval(d)).abs());
        }
        assert!(far < 0.1 * h.rms(), "{far} of {}", h.rms());
    }

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
