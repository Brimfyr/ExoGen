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
    /// On ice, where rock makes a peak ring: the uplift collapsed into a pit (F72).
    CentralPit,
    /// On ice, where rock makes a multi-ring basin: a dome risen in the pit (F72).
    CentralDome,
}

impl CraterClass {
    pub fn name(self) -> &'static str {
        match self {
            CraterClass::Simple => "simple",
            CraterClass::Complex => "complex",
            CraterClass::PeakRing => "peak_ring",
            CraterClass::MultiRing => "multi_ring",
            CraterClass::CentralPit => "central_pit",
            CraterClass::CentralDome => "central_dome",
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
    /// The target is ice: its larger complex craters have central pits and domes (F72).
    pub icy: bool,
}

impl Scaling {
    pub fn classify(&self, diameter: f64) -> CraterClass {
        if diameter < self.simple_complex {
            CraterClass::Simple
        } else if diameter < self.peak_ring {
            CraterClass::Complex
        } else if diameter < self.multi_ring {
            if self.icy {
                CraterClass::CentralPit
            } else {
                CraterClass::PeakRing
            }
        } else if self.icy {
            CraterClass::CentralDome
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
    /// The width of the outermost wall terrace, m (F67); none for a simple crater.
    pub terrace_width: f64,
    pub peak_height: f64,
    pub peak_radius: f64,
    /// A central pit's radius, its raised rim's outer radius, m, its depth below the floor
    /// and its rim's height, m (F72); none outside the pit classes.
    pub pit_radius: f64,
    pub pit_rim_radius: f64,
    pub pit_depth: f64,
    pub pit_rim_height: f64,
    /// A dome's radius and its height above the pit's floor, m (F72).
    pub dome_radius: f64,
    pub dome_height: f64,
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
    // A complex crater's relief scales with the strength length its transition measures
    // (F71): the lunar crater of the same multiple of the transition, scaled by the
    // body's transition against the Moon's.
    let strength = (s.simple_complex / TERRACE_REFERENCE_TRANSITION).max(1e-3);
    let lunar = |law: (f64, f64)| pow_law(law, d_km / strength) * 1e3 * strength;
    match class {
        CraterClass::Simple => Dimensions {
            class,
            radius,
            depth: pow_law(PIKE_SIMPLE_DEPTH, d_km) * 1e3,
            rim_height: pow_law(PIKE_SIMPLE_RIM, d_km) * 1e3,
            floor_radius: 0.0,
            terrace_width: 0.0,
            peak_height: 0.0,
            peak_radius: 0.0,
            pit_radius: 0.0,
            pit_rim_radius: 0.0,
            pit_depth: 0.0,
            pit_rim_height: 0.0,
            dome_radius: 0.0,
            dome_height: 0.0,
            rings: [(0.0, 0.0); MAX_RINGS],
            ring_count: 0,
            ejecta_radius: EJECTA_RADIUS_MULTIPLE * radius,
            ejecta_exponent: EJECTA_EXPONENT,
        },
        _ => {
            let depth = lunar(PIKE_COMPLEX_DEPTH).min(MAX_BASIN_DEPTH_LUNAR * strength);
            let rim_height = lunar(PIKE_COMPLEX_RIM).min(MAX_RIM_HEIGHT_LUNAR * strength);
            let floor_radius = (0.5 * pow_law(FLOOR_DIAMETER, d_km) * 1e3)
                .min(0.75 * radius)
                .max(0.3 * radius);
            // The outermost terrace's width, scaled from the Moon's by the strength
            // length the transition diameter measures (F67).
            let terrace_width = pow_law(TERRACE_WIDTH, d_km)
                * 1e3
                * (s.simple_complex / TERRACE_REFERENCE_TRANSITION)
                    .max(1e-3)
                    .powf(1.0 - TERRACE_WIDTH.1);
            let uplift = lunar(CENTRAL_PEAK_HEIGHT);
            let (peak_height, peak_radius) = if class == CraterClass::Complex {
                (uplift, PEAK_RADIUS_FRACTION * radius)
            } else {
                (0.0, 0.0)
            };
            // On ice the uplift collapsed into a pit, and in the largest craters a dome
            // rose in it (F72).
            let pitted = matches!(class, CraterClass::CentralPit | CraterClass::CentralDome);
            let (pit_radius, pit_rim_radius, pit_depth, pit_rim_height) = if pitted {
                let share = (PIT_DIAMETER_SHARE
                    * (diameter / s.peak_ring.max(1.0)).max(1.0).powf(PIT_GROWTH))
                .min(PIT_SHARE_MAX);
                (
                    0.5 * share * diameter,
                    0.5 * (PIT_RIM_DIAMETER_SHARE * diameter).max(1.5 * share * diameter),
                    PIT_DEPTH_SHARE * depth,
                    uplift,
                )
            } else {
                (0.0, 0.0, 0.0, 0.0)
            };
            let (dome_radius, dome_height) = if class == CraterClass::CentralDome {
                (DOME_PIT_SHARE * pit_radius, pit_depth + DOME_RISE * depth)
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
                terrace_width,
                peak_height,
                peak_radius,
                pit_radius,
                pit_rim_radius,
                pit_depth,
                pit_rim_height,
                dome_radius,
                dome_height,
                rings,
                ring_count,
                ejecta_radius,
                ejecta_exponent,
            }
        }
    }
}

/// The transient cavity's diameter for a final crater: a complex crater's from McKinnon
/// and Schenk's (1995) scaling; below the transition the rim's collapse widens the cavity
/// by the share the scaling gives at the transition, so the two meet.
pub fn transient_diameter(diameter: f64, s: &Scaling) -> f64 {
    let (k, p) = TRANSIENT_COMPLEX;
    let d_sc = s.simple_complex.max(1.0) / 1e3;
    let d = diameter.max(0.0) / 1e3;
    let tr = if d >= d_sc {
        (d * d_sc.powf(p - 1.0) / k).powf(1.0 / p)
    } else {
        d * k.powf(-1.0 / p)
    };
    tr * 1e3
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
            let h = if r <= dim.floor_radius {
                -dim.depth
            } else {
                let s = smoothstep((r - dim.floor_radius) / (rr - dim.floor_radius).max(1.0));
                -dim.depth + (dim.depth + dim.rim_height) * s
            };
            // The central peak is added by the stamp as a cluster of massifs.
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

    /// The degraded profile relaxed by wavelength (F70): each band between Gaussian blurs
    /// a factor of two apart, from the table's step up to RELAX_LONGEST_BLUR of the radius,
    /// keeps the share `kept` gives its wavelength, and the smooth remainder the share of
    /// the longest band. What is finer than the step stands.
    pub fn relaxed(dim: Dimensions, kappa_t: f64, kept: &dyn Fn(f64) -> f64) -> Self {
        let base = Self::new(dim, kappa_t);
        let mut values = base.values.clone();
        let longest = RELAX_LONGEST_BLUR * dim.radius;
        let mut blur = base.step;
        let mut finer = Self::new(dim, kappa_t + 0.25 * blur * blur);
        loop {
            let lost = 1.0 - kept(RELAX_BAND_WAVELENGTH * blur).clamp(0.0, 1.0);
            let next = 2.0 * blur;
            if next > longest {
                for (v, f) in values.iter_mut().zip(&finer.values) {
                    *v -= lost * f;
                }
                break;
            }
            let coarser = Self::new(dim, kappa_t + 0.25 * next * next);
            for ((v, f), c) in values.iter_mut().zip(&finer.values).zip(&coarser.values) {
                *v -= lost * (f - c);
            }
            finer = coarser;
            blur = next;
        }
        Self { values, ..base }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Diagnostic: the steepest wall of each profile and the share of the radius steeper
    /// than 10 degrees, on Callisto and the Moon, fresh and degraded, for the slopes from
    /// which the ground sheds (F98).
    #[test]
    #[ignore]
    fn write_wall_slopes() {
        let callisto = Scaling {
            simple_complex: 3_244.0,
            peak_ring: 27_577.0,
            multi_ring: 64_887.0,
            gravity: 1.235,
            icy: true,
        };
        for (name, s, kappa) in [("callisto", callisto, 4.8), ("moon", lunar(), 5.5)] {
            for d in [5e3, 10e3, 20e3, 40e3, 80e3] {
                let dim = dimensions(d, &s);
                for eq_gyr in [0.0, 1.0, 5.7, 30.0, 300.0] {
                    let t = ProfileTable::new(dim, kappa * eq_gyr * 1e3);
                    let mut max = 0.0f64;
                    let mut wide = 0.0;
                    for k in 1..t.values.len() {
                        let a = ((t.values[k] - t.values[k - 1]).abs() / t.step)
                            .atan()
                            .to_degrees();
                        max = max.max(a);
                        if a > 10.0 {
                            wide += t.step;
                        }
                    }
                    eprintln!(
                        "{name} D {:>3.0} km E {eq_gyr:>5}: depth {:>5.0} m, steepest {max:>4.1} deg, steeper than 10 deg over {:.2} R",
                        d / 1e3,
                        t.measured_depth(),
                        wide / dim.radius
                    );
                }
            }
        }
    }

    #[test]
    fn transient_cavities_follow_mckinnon_and_schenk() {
        // Tyre, 38 km across on Europa (transition 3 km): Singer et al. 2023 give 23 km.
        let europa = Scaling {
            simple_complex: 3_050.0,
            peak_ring: 25_922.0,
            multi_ring: 60_993.0,
            gravity: 1.315,
            icy: true,
        };
        let tyre = transient_diameter(38e3, &europa);
        assert!((21e3..26e3).contains(&tyre), "{tyre}");
        // The simple and complex branches meet at the transition.
        let below = transient_diameter(3_049.0, &europa);
        let above = transient_diameter(3_051.0, &europa);
        assert!((below - above).abs() < 5.0, "{below} {above}");
    }

    #[test]
    fn a_relaxed_crater_loses_its_bowl_before_its_rim() {
        // F70: the long wavelengths go first, so the floor rises toward the surface while
        // the rim, a short wavelength, stands; a shell that keeps every wavelength leaves
        // the profile as it was.
        let dim = dimensions(40_000.0, &lunar());
        let fresh = ProfileTable::new(dim, 0.0);
        let relaxed =
            ProfileTable::relaxed(dim, 0.0, &|l| if l > 2.0 * dim.radius { 0.1 } else { 1.0 });
        let rim = |t: &ProfileTable| t.values.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        assert!(
            relaxed.height(0.0) > 0.5 * fresh.height(0.0),
            "floor {} against {}",
            relaxed.height(0.0),
            fresh.height(0.0)
        );
        assert!(
            rim(&relaxed) > 0.8 * rim(&fresh),
            "rim {} against {}",
            rim(&relaxed),
            rim(&fresh)
        );
        let kept = ProfileTable::relaxed(dim, 0.0, &|_| 1.0);
        assert!(
            kept.values
                .iter()
                .zip(&fresh.values)
                .all(|(a, b)| (a - b).abs() < 1e-9)
        );
    }

    fn lunar() -> Scaling {
        Scaling {
            simple_complex: 16_500.0,
            peak_ring: 140_000.0,
            multi_ring: 330_000.0,
            gravity: 1.62,
            icy: false,
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
    fn craters_on_ice_are_shallow_and_pitted() {
        // F71: Ganymede's complex craters are 60 to 70 percent shallower than the Moon's
        // (Croft 1981); its central pit craters are 0.5 to 1.7 km deep (Bray et al. 2012).
        let ganymede = Scaling {
            simple_complex: 2_812.0,
            peak_ring: 23_904.0,
            multi_ring: 56_246.0,
            gravity: 1.428,
            icy: true,
        };
        let moon = Scaling {
            simple_complex: TERRACE_REFERENCE_TRANSITION,
            peak_ring: 139_900.0,
            multi_ring: 329_100.0,
            gravity: 1.62,
            icy: false,
        };
        let g40 = dimensions(40e3, &ganymede);
        let m40 = dimensions(40e3, &moon);
        assert!((0.5e3..1.7e3).contains(&g40.depth), "{}", g40.depth);
        assert!(g40.depth < 0.4 * m40.depth, "{} {}", g40.depth, m40.depth);
        // The Moon keeps Pike's laws.
        assert!((m40.depth - 1.044e3 * 40f64.powf(0.301)).abs() < 1.0);
        // F72: a 40 km crater on Ganymede has a central pit a fifth of it across and no
        // peak or rings; a 70 km one a dome in its pit; the Moon's do not.
        assert_eq!(g40.class, CraterClass::CentralPit);
        assert_eq!(g40.peak_height, 0.0);
        assert_eq!(g40.ring_count, 0);
        let share40 = g40.pit_radius / g40.radius;
        assert!((0.15..0.3).contains(&share40), "{share40}");
        assert!(g40.pit_depth > 0.0 && g40.pit_rim_height > 0.0 && g40.dome_radius == 0.0);
        let g70 = dimensions(70e3, &ganymede);
        assert_eq!(g70.class, CraterClass::CentralDome);
        assert!(g70.dome_radius > 0.0 && g70.dome_height > g70.pit_depth);
        assert_eq!(g70.ring_count, 0);
        // The pit widens with the crater: a third of a 100 km crater (Bray et al. 2012).
        let g100 = dimensions(100e3, &ganymede);
        let share100 = g100.pit_radius / g100.radius;
        assert!(
            share100 > share40 && (0.28..0.4).contains(&share100),
            "{share100}"
        );
        assert_eq!(m40.pit_radius, 0.0);
        // Mimas's Herschel, 139 km across on a body whose transition is near 63 km, is 10
        // to 12 km deep.
        let mimas = Scaling {
            simple_complex: 63e3,
            peak_ring: 535e3,
            multi_ring: 1_260e3,
            gravity: 0.064,
            icy: true,
        };
        let herschel = dimensions(139e3, &mimas);
        assert!((9e3..13e3).contains(&herschel.depth), "{}", herschel.depth);
    }

    #[test]
    fn gravity_moves_the_transition() {
        let big = Scaling {
            simple_complex: 2_900.0,
            peak_ring: 24_650.0,
            multi_ring: 58_000.0,
            gravity: 9.2,
            icy: false,
        };
        let small = Scaling {
            simple_complex: 19_400.0,
            peak_ring: 164_900.0,
            multi_ring: 388_000.0,
            gravity: 1.38,
            icy: false,
        };
        let d = 12_000.0;
        assert_eq!(dimensions(d, &small).class, CraterClass::Simple);
        assert_eq!(dimensions(d, &big).class, CraterClass::Complex);
        assert!(dimensions(d, &big).peak_height > 0.0);
        assert_eq!(dimensions(d, &small).peak_height, 0.0);
    }
}
