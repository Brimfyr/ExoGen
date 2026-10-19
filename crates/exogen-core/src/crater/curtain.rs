// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! The ejecta curtain of one crater and the rays it leaves (F91).
//!
//! An impact throws its ejecta out as a curtain whose mass is not even round the crater:
//! the target's relief and heterogeneity at the rim gather it into jets (Sabuwala et al.
//! 2018 made rays in a granular bed whose surface undulated, and a uniform blanket where
//! it was smooth). The curtain's mass by azimuth is the exponential of a spread of
//! harmonics of that unresolved heterogeneity, so that most of it leaves in some ten jets
//! (Elliott et al. 2018 counted 10.7 rays a crater on the Moon) with finer jets inside
//! them. Its finest harmonics also lobe the edge of the continuous blanket, so the blanket
//! and its rays share one outline.
//!
//! Beyond the blanket a ray is the share of the ground the jets' fragments excavate from
//! under the space-weathered skin as secondaries (Shoemaker 1965; Elliott et al. 2018):
//! one minus the exponential of the area they excavate, which falls with distance as the
//! fragments thin out and their secondaries shrink. A ray therefore narrows to its jet's
//! midline, breaks into segments and fades rather than ending. Elliott et al.'s
//! secondary-crater model sets the fall and their lunar ray lengths its scale. Noise bends
//! each jet's path and bunches its fragments into bursts, standing for the azimuthal
//! momentum and the clumping the curtain's heterogeneity gives them (PLAN principle 12).

use crate::constants::{
    CONTINUOUS_EJECTA_MULTIPLE, CURTAIN_EDGE_ORDER, CURTAIN_ORDERS, CURTAIN_SLOPE, CURTAIN_SPREAD,
    RAY_CLUMPING, RAY_DENSITY_SLOPE, RAY_DOWNRANGE_SPEED, RAY_LENGTH_COEFF_KM, RAY_LENGTH_EXPONENT,
    RAY_LENGTH_MAX_RADIUS_KM, RAY_VISIBLE, RAY_WARP,
};
use crate::crater::stamp::ejecta_modulation;
use crate::rng::Rng;

/// Bins of the mass table round the crater.
const TABLE: usize = 1024;
/// The bends of a jet's path: harmonic orders round the crater.
const WARP_ORDERS: [f64; 4] = [3.0, 5.0, 8.0, 13.0];
/// The bursts along the jets: lattice cells per turn and per e-fold of distance, and the
/// weight of each octave.
const CLUMP_OCTAVES: [(u32, f64, f64); 3] = [(256, 8.0, 1.0), (512, 16.0, 0.6), (1024, 32.0, 0.35)];

/// One harmonic of the curtain's heterogeneity.
#[derive(Clone, Copy)]
struct Harmonic {
    order: f64,
    amplitude: f64,
    phase: f64,
}

/// The harmonics of a crater's curtain: every order in CURTAIN_ORDERS, amplitudes falling
/// as the order to the CURTAIN_SLOPE, each scattered by a factor of 0.5 to 1.5.
fn harmonics(seed: u64) -> Vec<Harmonic> {
    let rng = Rng::new(seed ^ 0xC0A7_A1A5);
    let (lo, hi) = CURTAIN_ORDERS;
    (lo..=hi)
        .map(|m| {
            let k = m as u64;
            Harmonic {
                order: m as f64,
                amplitude: (m as f64).powf(-CURTAIN_SLOPE) * (0.5 + rng.uniform(k, 0)),
                phase: 2.0 * std::f64::consts::PI * rng.uniform(k, 1),
            }
        })
        .collect()
}

/// The sum of the harmonics `hs` over `bins` azimuths, by rotating each one's phasor.
fn tabulate(hs: &[Harmonic], bins: usize) -> Vec<f64> {
    let mut out = vec![0.0; bins];
    let step = 2.0 * std::f64::consts::PI / bins as f64;
    for h in hs {
        let (ds, dc) = (h.order * step).sin_cos();
        let (mut s, mut c) = h.phase.sin_cos();
        for v in out.iter_mut() {
            *v += h.amplitude * c;
            let c2 = c * dc - s * ds;
            s = s * dc + c * ds;
            c = c2;
        }
    }
    out
}

#[inline]
fn lookup(table: &[f64], phi: f64) -> f64 {
    let n = table.len();
    let x = phi.rem_euclid(2.0 * std::f64::consts::PI) / (2.0 * std::f64::consts::PI) * n as f64;
    let i = (x.floor() as usize) % n;
    let f = x - x.floor();
    table[i] * (1.0 - f) + table[(i + 1) % n] * f
}

/// The lobes of a continuous blanket's edge, at unit spread: the curtain's harmonics from
/// CURTAIN_EDGE_ORDER up to `max_order`, so the edge is ragged in the same fine jets that
/// start the rays, in many small lobes rather than a few petals.
pub fn edge_lobes(seed: u64, max_order: f64, bins: usize) -> Vec<f64> {
    let hs: Vec<Harmonic> = harmonics(seed)
        .into_iter()
        .filter(|h| h.order >= CURTAIN_EDGE_ORDER && h.order <= max_order)
        .collect();
    if hs.is_empty() {
        return vec![0.0; bins];
    }
    let power: f64 = hs.iter().map(|h| 0.5 * h.amplitude * h.amplitude).sum();
    tabulate(&hs, bins)
        .into_iter()
        .map(|v| v / power.sqrt())
        .collect()
}

/// The curtain's mass by azimuth, mean 1, from its harmonics up to `max_order`: the jets
/// a grid resolves, the finer ones averaged out (F96). The spread is the whole curtain's,
/// so a resolved jet carries the share of the mass it carries in the rays.
pub fn jet_mass(seed: u64, max_order: f64, bins: usize) -> Vec<f64> {
    let all = harmonics(seed);
    let power: f64 = all.iter().map(|h| 0.5 * h.amplitude * h.amplitude).sum();
    let hs: Vec<Harmonic> = all.into_iter().filter(|h| h.order <= max_order).collect();
    if hs.is_empty() {
        return vec![1.0; bins];
    }
    let mass: Vec<f64> = tabulate(&hs, bins)
        .into_iter()
        .map(|v| (CURTAIN_SPREAD * v / power.sqrt()).exp())
        .collect();
    let mean = mass.iter().sum::<f64>() / bins as f64;
    mass.into_iter().map(|m| m / mean).collect()
}

/// The ejecta curtain of one crater.
pub struct Curtain {
    /// The mass thrown per unit azimuth, mean 1 before the oblique modulation.
    mass: Vec<f64>,
    /// The median mass of the jets, the maxima of the mass above its mean.
    median_jet: f64,
    /// The largest mass times the reach its speed gives, for bounding the rays.
    peak_reach: f64,
    azimuth: f64,
    /// Amplitude of the downrange speed-up.
    speed: f64,
    /// Per bend: order, amplitude per e-fold, phase, rate of turn with ln r.
    warp: [(f64, f64, f64, f64); 4],
    seed: u64,
}

impl Curtain {
    /// The curtain of the crater `seed` struck from `azimuth` downrange at `angle` from
    /// the horizontal.
    pub fn new(seed: u64, azimuth: f64, angle: f64) -> Self {
        let hs = harmonics(seed);
        let power: f64 = hs.iter().map(|h| 0.5 * h.amplitude * h.amplitude).sum();
        let n = tabulate(&hs, TABLE);
        let mut mass: Vec<f64> = n
            .iter()
            .map(|v| (CURTAIN_SPREAD * v / power.sqrt()).exp())
            .collect();
        let mean = mass.iter().sum::<f64>() / TABLE as f64;
        let speed = RAY_DOWNRANGE_SPEED * angle.cos().max(0.0);
        let step = 2.0 * std::f64::consts::PI / TABLE as f64;
        for (b, m) in mass.iter_mut().enumerate() {
            *m *= ejecta_modulation(b as f64 * step, azimuth, angle) / mean;
        }
        let mut jets: Vec<f64> = (0..TABLE)
            .filter(|&b| {
                let m = mass[b];
                m > 1.0 && m > mass[(b + TABLE - 1) % TABLE] && m >= mass[(b + 1) % TABLE]
            })
            .map(|b| mass[b])
            .collect();
        jets.sort_by(|a, b| a.total_cmp(b));
        let median_jet = if jets.is_empty() {
            1.0
        } else {
            jets[jets.len() / 2]
        };
        let peak_reach = (0..TABLE)
            .map(|b| {
                let s = 1.0 + speed * (b as f64 * step - azimuth).cos();
                mass[b] * s.powf(2.0 * RAY_DENSITY_SLOPE)
            })
            .fold(0.0, f64::max);
        let rng = Rng::new(seed ^ 0xBE4D_7A75);
        let mut warp = [(0.0, 0.0, 0.0, 0.0); 4];
        for (i, (w, &k)) in warp.iter_mut().zip(WARP_ORDERS.iter()).enumerate() {
            let i = i as u64;
            *w = (
                k,
                RAY_WARP * (3.0 / k) * (0.5 + rng.uniform(i, 0)),
                2.0 * std::f64::consts::PI * rng.uniform(i, 1),
                0.5 + rng.uniform(i, 2),
            );
        }
        Self {
            mass,
            median_jet,
            peak_reach,
            azimuth,
            speed,
            warp,
            seed,
        }
    }

    /// The mass thrown toward the launch azimuth `phi`, mean about 1.
    #[inline]
    pub fn mass(&self, phi: f64) -> f64 {
        lookup(&self.mass, phi)
    }

    /// The largest mass round the crater.
    pub fn max_mass(&self) -> f64 {
        self.mass.iter().cloned().fold(0.0, f64::max)
    }

    /// The median mass of the jets.
    pub fn median_jet(&self) -> f64 {
        self.median_jet
    }

    /// How much faster than the mean the ejecta toward `phi` leave: downrange faster,
    /// uprange slower, by RAY_DOWNRANGE_SPEED at a grazing angle and not at all straight
    /// down (Anderson et al. 2003).
    #[inline]
    pub fn speed(&self, phi: f64) -> f64 {
        1.0 + self.speed * (phi - self.azimuth).cos()
    }

    /// The launch azimuth of the fragments that land at azimuth `phi` and log distance
    /// `x` = ln(r / R): the jets' paths bend past the continuous blanket.
    #[inline]
    pub fn launch_azimuth(&self, phi: f64, x: f64) -> f64 {
        let grow = (x - CONTINUOUS_EJECTA_MULTIPLE.ln()).max(0.0);
        if grow <= 0.0 {
            return phi;
        }
        let bend: f64 = self
            .warp
            .iter()
            .map(|&(k, a, p, w)| a * (k * phi + p + w * x).sin())
            .sum();
        phi + bend * grow
    }

    /// The azimuth at log distance `x` of the path launched toward `launch`: the inverse
    /// of `launch_azimuth`, by Newton's method where the bends do not fold the paths over.
    pub fn path_azimuth(&self, launch: f64, x: f64) -> f64 {
        let grow = (x - CONTINUOUS_EJECTA_MULTIPLE.ln()).max(0.0);
        let mut phi = launch;
        for _ in 0..12 {
            let f = self.launch_azimuth(phi, x) - launch;
            if f.abs() < 1e-12 {
                break;
            }
            let slope = 1.0
                + grow
                    * self
                        .warp
                        .iter()
                        .map(|&(k, a, p, w)| a * k * (k * phi + p + w * x).cos())
                        .sum::<f64>();
            phi -= if slope > 0.2 { f / slope } else { 0.5 * f };
        }
        phi
    }

    /// The bunching of the fragments thrown toward `launch` that land at log distance
    /// `x`, mean about 1: log-normal, in bursts stretched along the path.
    pub fn clumping(&self, launch: f64, x: f64) -> f64 {
        let turns = launch.rem_euclid(2.0 * std::f64::consts::PI) / (2.0 * std::f64::consts::PI);
        let rng = Rng::new(self.seed ^ 0xC1E4_B0B5);
        let mut g = 0.0;
        let mut power = 0.0;
        for (o, &(period, per_efold, weight)) in CLUMP_OCTAVES.iter().enumerate() {
            let u = turns * period as f64;
            let v = x.max(0.0) * per_efold;
            let (iu, iv) = (u.floor(), v.floor());
            let smooth = |t: f64| t * t * (3.0 - 2.0 * t);
            let (su, sv) = (smooth(u - iu), smooth(v - iv));
            let cu0 = (iu as i64).rem_euclid(period as i64) as u64;
            let cu1 = (cu0 + 1) % period as u64;
            let cv = iv as u64;
            let at = |cu: u64, cv: u64| {
                2.0 * rng.uniform(((o as u64) << 56) | (cu << 28) | (cv & 0x0FFF_FFFF), 0) - 1.0
            };
            let a = at(cu0, cv) + (at(cu1, cv) - at(cu0, cv)) * su;
            let b = at(cu0, cv + 1) + (at(cu1, cv + 1) - at(cu0, cv + 1)) * su;
            g += weight * (a + (b - a) * sv);
            power += weight * weight / 3.0;
        }
        let g = g / power.sqrt();
        (RAY_CLUMPING * g - 0.5 * RAY_CLUMPING * RAY_CLUMPING).exp()
    }
}

/// The strongest burst the clumping can make: every octave at its lattice's extreme.
fn strongest_burst() -> f64 {
    let weights: f64 = CLUMP_OCTAVES.iter().map(|o| o.2).sum();
    let power: f64 = CLUMP_OCTAVES.iter().map(|o| o.2 * o.2 / 3.0).sum();
    (RAY_CLUMPING * weights / power.sqrt() - 0.5 * RAY_CLUMPING * RAY_CLUMPING).exp()
}

/// The rays of one crater: the share of the ground its curtain's fragments excavated.
pub struct Rays {
    curtain: Curtain,
    radius: f64,
    /// The excavated area at one crater radius from a jet of unit mass.
    scale: f64,
    /// The area below which a jet leaves less than RAY_VISIBLE of the ground excavated
    /// even in its strongest burst.
    least: f64,
}

/// The median ray length, m, of a crater of radius `radius`, m: Elliott et al. (2018)'s
/// lunar fit, RAY_LENGTH_COEFF_KM (R / km)^RAY_LENGTH_EXPONENT, up to the largest crater
/// it was measured on, and beyond, self-similar, as many crater radii as that crater's
/// (F127).
pub fn ray_length(radius: f64) -> f64 {
    let fit = |r_km: f64| RAY_LENGTH_COEFF_KM * r_km.powf(RAY_LENGTH_EXPONENT);
    let r_km = radius / 1e3;
    let km = if r_km <= RAY_LENGTH_MAX_RADIUS_KM {
        fit(r_km)
    } else {
        fit(RAY_LENGTH_MAX_RADIUS_KM) * r_km / RAY_LENGTH_MAX_RADIUS_KM
    };
    km * 1e3
}

impl Rays {
    /// The rays of a crater of radius `radius`, m: scaled so that its median jet's
    /// fragments excavate half the ground at its median ray length (`ray_length`), the
    /// share Elliott et al. (2018) took as visible.
    pub fn new(curtain: Curtain, radius: f64) -> Self {
        let length = ray_length(radius);
        let scale = std::f64::consts::LN_2 / curtain.median_jet()
            * (length / radius).powf(RAY_DENSITY_SLOPE);
        Self {
            curtain,
            radius,
            scale,
            least: -(1.0 - RAY_VISIBLE).ln() / strongest_burst(),
        }
    }

    /// The curtain the rays came from.
    pub fn curtain(&self) -> &Curtain {
        &self.curtain
    }

    /// The share of the ground excavated to fresh material at distance `r`, m, and
    /// azimuth `phi` from the centre, where the skin the secondaries dig through leaves
    /// the share `skin` of the area they would excavate under the Moon's
    /// (RAY_SKIN_DEPTH_M); nothing inside the rim.
    pub fn excavated(&self, r: f64, phi: f64, skin: f64) -> f64 {
        if r <= self.radius {
            return 0.0;
        }
        let x = (r / self.radius).ln();
        let launch = self.curtain.launch_azimuth(phi, x);
        let reach = self.curtain.speed(launch).powi(2);
        let area = self.scale
            * skin
            * self.curtain.mass(launch)
            * (r / (self.radius * reach)).powf(-RAY_DENSITY_SLOPE);
        // Where no burst could make the jet show, its bursts are not worked out.
        if area < self.least {
            return 0.0;
        }
        1.0 - (-area * self.curtain.clumping(launch, x)).exp()
    }

    /// The distance, m, beyond which no jet excavates RAY_VISIBLE of the ground even in
    /// its strongest burst.
    pub fn reach(&self) -> f64 {
        self.radius
            * (self.scale * self.curtain.peak_reach / self.least).powf(1.0 / RAY_DENSITY_SLOPE)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The local maxima of a curtain's mass above its mean, the jets.
    fn jets(c: &Curtain) -> usize {
        (0..TABLE)
            .filter(|&b| {
                let m = c.mass[b];
                m > 1.0 && m > c.mass[(b + TABLE - 1) % TABLE] && m >= c.mass[(b + 1) % TABLE]
            })
            .count()
    }

    #[test]
    fn the_curtain_gathers_its_mass_in_jets() {
        let mut counts = Vec::new();
        let mut top_shares = Vec::new();
        for seed in 0..200u64 {
            let c = Curtain::new(seed, 0.0, std::f64::consts::FRAC_PI_2);
            let mean = c.mass.iter().sum::<f64>() / TABLE as f64;
            assert!((mean - 1.0).abs() < 1e-9, "mean {mean}");
            counts.push(jets(&c) as f64);
            // The mass the strongest fifth of the azimuths carries.
            let mut m = c.mass.clone();
            m.sort_by(|a, b| b.total_cmp(a));
            top_shares.push(m[..TABLE / 5].iter().sum::<f64>() / TABLE as f64);
        }
        let mean_count = counts.iter().sum::<f64>() / counts.len() as f64;
        let mean_share = top_shares.iter().sum::<f64>() / top_shares.len() as f64;
        eprintln!(
            "jets {mean_count:.1}, the strongest fifth of the azimuths carries {mean_share:.2}"
        );
        // Some ten to twenty jets, most of the mass in a fifth of the azimuths.
        assert!(mean_count > 8.0 && mean_count < 30.0, "jets {mean_count}");
        assert!(mean_share > 0.5, "top fifth carries {mean_share}");
    }

    /// The farthest distance along the path launched toward `launch` at which the
    /// excavated share is above a half: the ray's length, its gaps passed over as Elliott
    /// et al. passed over them where the ray went on.
    fn visible_length(c: &Curtain, rays: &Rays, launch: f64, radius: f64) -> f64 {
        let mut r = 2.0 * radius;
        let mut last = r;
        while r < 500.0 * radius {
            let phi = c.path_azimuth(launch, (r / radius).ln());
            if rays.excavated(r, phi, 1.0) > 0.5 {
                last = r;
            }
            r *= 1.02;
        }
        last
    }

    #[test]
    fn ray_lengths_follow_the_lunar_fit() {
        // The median over craters of the median length of their jets, against Elliott et
        // al.'s 6.59 R^1.27 km, for vertical impacts (no clumping bias: lengths are read
        // where the jets' midlines cross a half).
        for radius_km in [2.0f64, 10.0, 42.0] {
            let radius = radius_km * 1e3;
            let mut medians = Vec::new();
            for seed in 0..40u64 {
                let rays = Rays::new(Curtain::new(seed, 0.0, std::f64::consts::FRAC_PI_2), radius);
                let c = rays.curtain();
                let mut lengths: Vec<f64> = (0..TABLE)
                    .filter(|&b| {
                        let m = c.mass[b];
                        m > 1.0
                            && m > c.mass[(b + TABLE - 1) % TABLE]
                            && m >= c.mass[(b + 1) % TABLE]
                    })
                    .map(|b| {
                        let launch = b as f64 / TABLE as f64 * 2.0 * std::f64::consts::PI;
                        visible_length(c, &rays, launch, radius)
                    })
                    .collect();
                lengths.sort_by(|a, b| a.total_cmp(b));
                medians.push(lengths[lengths.len() / 2]);
            }
            medians.sort_by(|a, b| a.total_cmp(b));
            let got = medians[medians.len() / 2] / 1e3;
            let want = RAY_LENGTH_COEFF_KM * radius_km.powf(RAY_LENGTH_EXPONENT);
            eprintln!("R {radius_km} km: median ray {got:.0} km against {want:.0} km");
            assert!(
                (got / want - 1.0).abs() < 0.35,
                "R {radius_km} km: median ray {got:.0} km against {want:.0} km"
            );
        }
    }

    /// F127: beyond the largest crater measured, a crater's rays are its scaled up: as
    /// many crater radii long, and as dense at the same distance in crater radii.
    #[test]
    fn a_basin_s_rays_are_a_large_crater_s_scaled_up() {
        let edge = RAY_LENGTH_MAX_RADIUS_KM * 1e3;
        let ratio = |r: f64| ray_length(r) / r;
        assert!((ratio(4.0 * edge) - ratio(edge)).abs() < 1e-9);
        assert!(ratio(0.5 * edge) < ratio(edge));
        let curtain = || Curtain::new(7, 0.0, std::f64::consts::FRAC_PI_2);
        let (crater, basin) = (Rays::new(curtain(), edge), Rays::new(curtain(), 4.0 * edge));
        for x in [3.0, 10.0, 30.0] {
            let a = crater.excavated(x * edge, 0.3, 1.0);
            let b = basin.excavated(x * 4.0 * edge, 0.3, 1.0);
            assert!((a - b).abs() < 1e-9, "{x}: {a} {b}");
        }
    }

    #[test]
    fn an_oblique_impact_throws_its_rays_farther_downrange() {
        let radius = 10e3;
        let mut down = 0.0;
        let mut up = 0.0;
        for seed in 0..40u64 {
            let rays = Rays::new(Curtain::new(seed, 0.0, 40f64.to_radians()), radius);
            // The excavated share summed over a ring at ten radii, by half.
            for b in 0..360 {
                let phi = (b as f64).to_radians();
                let e = rays.excavated(10.0 * radius, phi, 1.0);
                if phi.cos() > 0.0 {
                    down += e;
                } else {
                    up += e;
                }
            }
        }
        eprintln!("excavated at ten radii: downrange {down:.1}, uprange {up:.1}");
        assert!(down > 1.5 * up, "downrange {down:.1}, uprange {up:.1}");
    }

    #[test]
    fn a_bent_path_inverts() {
        let c = Curtain::new(7, 0.3, 1.0);
        for &x in &[0.5f64, 1.5, 3.0] {
            for b in 0..16 {
                let launch = b as f64 * 0.4;
                let phi = c.path_azimuth(launch, x);
                let back = c.launch_azimuth(phi, x);
                assert!(
                    (back - launch).abs() < 1e-6,
                    "x {x}: {launch} -> {phi} -> {back}"
                );
            }
        }
    }

    #[test]
    fn the_reach_bounds_the_rays() {
        let radius = 20e3;
        for seed in 0..10u64 {
            let rays = Rays::new(Curtain::new(seed, 1.0, 0.6), radius);
            let reach = rays.reach();
            for b in 0..720 {
                let phi = b as f64 / 720.0 * 2.0 * std::f64::consts::PI;
                let e = rays.excavated(reach * 1.01, phi, 1.0);
                assert!(e < RAY_VISIBLE, "seed {seed}: {e} beyond the reach");
            }
        }
    }
}
