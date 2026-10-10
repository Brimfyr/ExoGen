// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! Crater classes, dimensions and radial profiles (PLAN.md 4.5).
//!
//! Dimensions follow Pike's lunar fits in kilometres, converted to metres, with the
//! transition diameters supplied by the parameter graph so that gravity and target scaling
//! enter in one place. Profiles are radial height functions relative to the pre-impact
//! reference surface, continuous through the rim into the ejecta.

use crate::constants::*;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CraterClass {
    Simple,
    Complex,
    PeakRing,
    MultiRing,
}

impl CraterClass {
    pub fn name(self) -> &'static str {
        match self {
            CraterClass::Simple => "simple",
            CraterClass::Complex => "complex",
            CraterClass::PeakRing => "peak_ring",
            CraterClass::MultiRing => "multi_ring",
        }
    }
}

/// Body-dependent scaling inputs, all in SI.
#[derive(Clone, Copy, Debug)]
pub struct Scaling {
    pub simple_complex: f64,
    pub peak_ring: f64,
    pub multi_ring: f64,
    pub gravity: f64,
}

impl Scaling {
    pub fn classify(&self, diameter: f64) -> CraterClass {
        if diameter < self.simple_complex {
            CraterClass::Simple
        } else if diameter < self.peak_ring {
            CraterClass::Complex
        } else if diameter < self.multi_ring {
            CraterClass::PeakRing
        } else {
            CraterClass::MultiRing
        }
    }
}

/// Crater dimensions in metres.
#[derive(Clone, Copy, Debug)]
pub struct Dimensions {
    pub class: CraterClass,
    pub radius: f64,
    pub depth: f64,
    pub rim_height: f64,
    pub floor_radius: f64,
    pub peak_height: f64,
    pub peak_radius: f64,
    /// Ridge rings inside the rim as (radius, height), outermost first.
    pub rings: [(f64, f64); MAX_RINGS],
    pub ring_count: usize,
    pub ejecta_radius: f64,
    pub ejecta_exponent: f64,
}

fn pow_law(law: (f64, f64), d_km: f64) -> f64 {
    law.0 * d_km.powf(law.1)
}

pub fn dimensions(diameter: f64, s: &Scaling) -> Dimensions {
    let class = s.classify(diameter);
    let d_km = diameter / 1e3;
    let radius = 0.5 * diameter;
    let g_scale = (GRAVITY_MOON / s.gravity.max(0.05)).powf(0.25);
    match class {
        CraterClass::Simple => Dimensions {
            class,
            radius,
            depth: pow_law(PIKE_SIMPLE_DEPTH, d_km) * 1e3,
            rim_height: pow_law(PIKE_SIMPLE_RIM, d_km) * 1e3,
            floor_radius: 0.0,
            peak_height: 0.0,
            peak_radius: 0.0,
            rings: [(0.0, 0.0); MAX_RINGS],
            ring_count: 0,
            ejecta_radius: EJECTA_RADIUS_MULTIPLE * radius,
            ejecta_exponent: EJECTA_EXPONENT,
        },
        _ => {
            let depth =
                (pow_law(PIKE_COMPLEX_DEPTH, d_km) * 1e3).min(MAX_BASIN_DEPTH_LUNAR * g_scale);
            let rim_height =
                (pow_law(PIKE_COMPLEX_RIM, d_km) * 1e3).min(MAX_RIM_HEIGHT_LUNAR * g_scale);
            let floor_radius = (0.5 * pow_law(FLOOR_DIAMETER, d_km) * 1e3)
                .min(0.75 * radius)
                .max(0.3 * radius);
            let (peak_height, peak_radius) = if class == CraterClass::Complex {
                (pow_law(CENTRAL_PEAK_HEIGHT, d_km) * 1e3, 0.15 * radius)
            } else {
                (0.0, 0.0)
            };
            // A peak ring at half the radius; multi-ring basins add a ring for each
            // sqrt(2) of diameter above the onset, each at 1/sqrt(2) of the last.
            let mut rings = [(0.0, 0.0); MAX_RINGS];
            let ring_count = match class {
                CraterClass::PeakRing => {
                    rings[0] = (0.5 * radius, RING_HEIGHT_FRACTION * depth);
                    1
                }
                CraterClass::MultiRing => {
                    let extra = ((diameter / s.multi_ring).ln() / std::f64::consts::SQRT_2.ln())
                        .floor()
                        .max(0.0) as usize;
                    let n = (1 + extra).min(MAX_RINGS);
                    for (k, ring) in rings.iter_mut().enumerate().take(n) {
                        *ring = (
                            radius / std::f64::consts::SQRT_2.powi(k as i32 + 1),
                            RING_HEIGHT_FRACTION * depth * RING_HEIGHT_DECAY.powi(k as i32),
                        );
                    }
                    n
                }
                _ => 0,
            };
            let (ejecta_radius, ejecta_exponent) = if class == CraterClass::MultiRing {
                (BASIN_EJECTA_RADIUS_MULTIPLE * radius, BASIN_EJECTA_EXPONENT)
            } else {
                (EJECTA_RADIUS_MULTIPLE * radius, EJECTA_EXPONENT)
            };
            Dimensions {
                class,
                radius,
                depth,
                rim_height,
                floor_radius,
                peak_height,
                peak_radius,
                rings,
                ring_count,
                ejecta_radius,
                ejecta_exponent,
            }
        }
    }
}

/// Impact angle from horizontal for a uniform draw `u`: isotropic impactors strike a
/// sphere with probability density sin 2θ (Gilbert 1893; Shoemaker 1962), whose
/// cumulative distribution is sin²θ. The most probable angle is 45°.
pub fn impact_angle(u: f64) -> f64 {
    u.clamp(0.0, 1.0).sqrt().asin()
}

/// Elongation of the cavity for an impact angle. Cavities stay circular down to the angle
/// below which `elliptical_fraction` of all craters form (about 12° for 5%), then elongate
/// toward grazing incidence.
pub fn ellipticity_for(angle: f64, elliptical_fraction: f64) -> f64 {
    let threshold = elliptical_fraction.clamp(0.0, 1.0).sqrt().asin();
    if threshold <= 0.0 || angle >= threshold {
        return 1.0;
    }
    1.0 + ELLIPTICITY_MAX_EXCESS * ((threshold - angle) / threshold).powf(1.5)
}

fn smoothstep(t: f64) -> f64 {
    let t = t.clamp(0.0, 1.0);
    t * t * (3.0 - 2.0 * t)
}

/// Fresh radial profile height at distance `r` from the centre, relative to the reference
/// surface. Inside the rim the value is the surface itself; outside it is ejecta thickness.
pub fn fresh_height(dim: &Dimensions, r: f64) -> f64 {
    let rr = dim.radius;
    if r >= rr {
        if r >= dim.ejecta_radius {
            return 0.0;
        }
        // Rings, including a multi-ring basin's outer ring, are added by the stamp with
        // their own shapes.
        return dim.rim_height * (r / rr).powf(dim.ejecta_exponent);
    }
    match dim.class {
        CraterClass::Simple => {
            let x = r / rr;
            -dim.depth + (dim.depth + dim.rim_height) * x * x
        }
        _ => {
            let mut h = if r <= dim.floor_radius {
                -dim.depth
            } else {
                let s = smoothstep((r - dim.floor_radius) / (rr - dim.floor_radius).max(1.0));
                -dim.depth + (dim.depth + dim.rim_height) * s
            };
            if dim.peak_height > 0.0 {
                let sigma = 0.5 * dim.peak_radius;
                h += dim.peak_height * (-(r / sigma).powi(2)).exp();
            }
            h
        }
    }
}

/// A radial profile sampled into a table and optionally diffused.
#[derive(Clone, Debug)]
pub struct ProfileTable {
    pub dim: Dimensions,
    /// Height at `k * step` for k in 0..n, covering 0 to the ejecta radius.
    pub values: Vec<f64>,
    pub step: f64,
    /// Radial blur width applied for degradation, metres.
    pub sigma: f64,
}

impl ProfileTable {
    pub const SAMPLES: usize = 128;

    /// Build the fresh profile and apply topographic diffusion for an elapsed `kappa_t`
    /// (m²) by a radial Gaussian blur in the table, which is a close approximation for the
    /// small blur widths that matter at map scale.
    pub fn new(dim: Dimensions, kappa_t: f64) -> Self {
        let n = Self::SAMPLES;
        let step = dim.ejecta_radius / (n - 1) as f64;
        let fresh: Vec<f64> = (0..n)
            .map(|k| fresh_height(&dim, k as f64 * step))
            .collect();
        // Two-dimensional diffusion for kappa t has sigma² = 2 kappa t per axis; a blur along
        // the radial table alone needs twice that variance to lower a bowl's centre by the
        // same amount, which is the quantity that matters for the degradation state.
        let sigma = (4.0 * kappa_t.max(0.0)).sqrt();
        let values = if sigma < 0.25 * step {
            fresh
        } else {
            let w = ((3.0 * sigma / step).ceil() as i64).min(n as i64);
            let weights: Vec<f64> = (0..=w)
                .map(|j| (-(j as f64 * step / sigma).powi(2) * 0.5).exp())
                .collect();
            (0..n as i64)
                .map(|k| {
                    let mut num = 0.0;
                    let mut den = 0.0;
                    for j in -w..=w {
                        let idx = (k + j).abs().min(2 * (n as i64 - 1) - (k + j).abs());
                        let idx = idx.clamp(0, n as i64 - 1) as usize;
                        let weight = weights[j.unsigned_abs() as usize];
                        num += weight * fresh[idx];
                        den += weight;
                    }
                    num / den
                })
                .collect()
        };
        Self {
            dim,
            values,
            step,
            sigma,
        }
    }

    /// True when the blur is so wide that no relief of the crater remains.
    pub fn is_erased(&self) -> bool {
        self.sigma > 2.0 * self.dim.radius
    }

    /// Height at radius `r` (metres) by linear interpolation; zero beyond the ejecta.
    #[inline]
    pub fn height(&self, r: f64) -> f64 {
        if r >= self.dim.ejecta_radius || r < 0.0 {
            return 0.0;
        }
        let x = r / self.step;
        let k = x.floor() as usize;
        if k + 1 >= self.values.len() {
            return self.values[self.values.len() - 1];
        }
        let f = x - k as f64;
        self.values[k] * (1.0 - f) + self.values[k + 1] * f
    }

    /// Rim-to-floor depth of the (possibly degraded) profile.
    pub fn measured_depth(&self) -> f64 {
        let lo = self.values.iter().cloned().fold(f64::INFINITY, f64::min);
        let hi = self
            .values
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        hi - lo
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lunar() -> Scaling {
        Scaling {
            simple_complex: 16_500.0,
            peak_ring: 140_000.0,
            multi_ring: 330_000.0,
            gravity: 1.62,
        }
    }

    #[test]
    fn classes_and_depths_follow_pike() {
        let s = lunar();
        let simple = dimensions(5_000.0, &s);
        assert_eq!(simple.class, CraterClass::Simple);
        assert!((simple.depth / 5_000.0 - 0.2).abs() < 0.02);
        let copernicus = dimensions(93_000.0, &s);
        assert_eq!(copernicus.class, CraterClass::Complex);
        assert!(
            (copernicus.depth - 4_090.0).abs() < 150.0,
            "{}",
            copernicus.depth
        );
        assert!((copernicus.rim_height - 1_440.0).abs() < 100.0);
        assert!(copernicus.peak_height > 900.0 && copernicus.peak_height < 1_500.0);
        let schrodinger = dimensions(320_000.0, &s);
        assert_eq!(schrodinger.class, CraterClass::PeakRing);
        let orientale = dimensions(930_000.0, &s);
        assert_eq!(orientale.class, CraterClass::MultiRing);
        assert!(orientale.depth <= MAX_BASIN_DEPTH_LUNAR + 1.0);
        // Three rings inside the Cordillera rim, at the Outer Rook, Inner Rook and inner
        // ring radii (620, 480 and 320 km across).
        assert_eq!(orientale.ring_count, 3);
        let diameters: Vec<f64> = orientale.rings[..3]
            .iter()
            .map(|r| 2.0 * r.0 / 1e3)
            .collect();
        for (got, want) in diameters.iter().zip([620.0, 480.0, 320.0]) {
            assert!((got / want - 1.0).abs() < 0.08, "{diameters:?}");
        }
        assert_eq!(schrodinger.ring_count, 1);
    }

    #[test]
    fn profile_is_continuous_at_the_rim_and_measures_its_own_depth() {
        let s = lunar();
        for d in [3_000.0, 30_000.0, 200_000.0, 600_000.0] {
            let dim = dimensions(d, &s);
            let inside = fresh_height(&dim, dim.radius * 0.999);
            let outside = fresh_height(&dim, dim.radius * 1.001);
            assert!(
                (inside - outside).abs() < 0.02 * dim.rim_height + 1.0,
                "D={d}: {inside} vs {outside}"
            );
            let table = ProfileTable::new(dim, 0.0);
            let measured = table.measured_depth();
            let expected = dim.depth + dim.rim_height;
            assert!(
                (measured - expected).abs() < 0.03 * expected + 1.0,
                "D={d}: {measured} vs {expected}"
            );
        }
    }

    #[test]
    fn diffusion_shallows_a_small_crater() {
        let s = lunar();
        let dim = dimensions(1_000.0, &s);
        let fresh = ProfileTable::new(dim, 0.0).measured_depth();
        let old = ProfileTable::new(dim, 20_000.0).measured_depth();
        assert!(old < 0.7 * fresh, "fresh {fresh} old {old}");
    }

    /// Impact angles follow sin 2θ: a mean of 45°, with the share below 12° being sin² 12°,
    /// and only the elliptical fraction of cavities elongated.
    #[test]
    fn impact_angles_follow_sin_two_theta_and_few_cavities_are_elliptical() {
        let n = 20_000;
        let (mut sum, mut below_12, mut elliptical) = (0.0, 0, 0);
        for i in 0..n {
            let a = impact_angle((i as f64 + 0.5) / n as f64);
            sum += a;
            if a < 12f64.to_radians() {
                below_12 += 1;
            }
            if ellipticity_for(a, 0.05) > 1.0 {
                elliptical += 1;
            }
        }
        let mean_deg = (sum / n as f64).to_degrees();
        assert!((mean_deg - 45.0).abs() < 0.5, "{mean_deg}");
        let f12 = below_12 as f64 / n as f64;
        assert!(
            (f12 - 12f64.to_radians().sin().powi(2)).abs() < 0.005,
            "{f12}"
        );
        let fe = elliptical as f64 / n as f64;
        assert!((fe - 0.05).abs() < 0.005, "{fe}");
        assert_eq!(ellipticity_for(0.8, 0.05), 1.0);
        assert!(ellipticity_for(0.02, 0.05) > 1.5);
    }

    #[test]
    fn gravity_moves_the_transition() {
        let big = Scaling {
            simple_complex: 2_900.0,
            peak_ring: 24_650.0,
            multi_ring: 58_000.0,
            gravity: 9.2,
        };
        let small = Scaling {
            simple_complex: 19_400.0,
            peak_ring: 164_900.0,
            multi_ring: 388_000.0,
            gravity: 1.38,
        };
        let d = 12_000.0;
        assert_eq!(dimensions(d, &small).class, CraterClass::Simple);
        assert_eq!(dimensions(d, &big).class, CraterClass::Complex);
        assert!(dimensions(d, &big).peak_height > 0.0);
        assert_eq!(dimensions(d, &small).peak_height, 0.0);
    }
}
