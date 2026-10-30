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
    RAY_CLUMPING, RAY_CLUSTER_ASPECT, RAY_CLUSTER_CAP, RAY_CLUSTER_LEAST, RAY_CLUSTER_MAX,
    RAY_CLUSTER_SPLAY, RAY_CLUSTER_WIDTH, RAY_DENSITY_SLOPE, RAY_DOWNRANGE_SPEED,
    RAY_LENGTH_COEFF_KM, RAY_LENGTH_EXPONENT, RAY_LENGTH_LARGE_EXPONENT, RAY_LENGTH_MAX_RADIUS_KM,
    RAY_VISIBLE, RAY_WARP, SECONDARY_CENSUS_COUNT, SECONDARY_CENSUS_FLOOR_KM,
    SECONDARY_CENSUS_LARGEST_KM, SECONDARY_CENSUS_SHAPE, SECONDARY_MAX_COEFF_KM,
    SECONDARY_MAX_EXPONENT, SECONDARY_SFD_SLOPE, SECONDARY_SIZE_FALL,
};
use crate::crater::production::poisson;
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
    /// The median ray length, m (`ray_length`).
    length: f64,
}

/// The median ray length, m, of a crater of radius `radius`, m: Elliott et al. (2018)'s
/// lunar fit, RAY_LENGTH_COEFF_KM (R / km)^RAY_LENGTH_EXPONENT, up to the largest crater
/// it was measured on, and beyond, growing more slowly than the crater, as the radius to
/// RAY_LENGTH_LARGE_EXPONENT (F128).
pub fn ray_length(radius: f64) -> f64 {
    let fit = |r_km: f64| RAY_LENGTH_COEFF_KM * r_km.powf(RAY_LENGTH_EXPONENT);
    let r_km = radius / 1e3;
    let km = if r_km <= RAY_LENGTH_MAX_RADIUS_KM {
        fit(r_km)
    } else {
        fit(RAY_LENGTH_MAX_RADIUS_KM)
            * (r_km / RAY_LENGTH_MAX_RADIUS_KM).powf(RAY_LENGTH_LARGE_EXPONENT)
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
            length,
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

    /// The distance, m, beyond which no jet's fragments dig below the Moon's skin
    /// (`fresh_share`), so that no ray shows: the fastest jet's end.
    pub fn fresh_reach(&self) -> f64 {
        let fastest = (0..TABLE)
            .map(|b| {
                self.curtain
                    .speed(b as f64 / TABLE as f64 * std::f64::consts::TAU)
            })
            .fold(1.0, f64::max);
        self.length * fastest * fastest
    }

    /// The distance, m, beyond which no jet excavates RAY_VISIBLE of the ground even in
    /// its strongest burst.
    pub fn reach(&self) -> f64 {
        self.radius
            * (self.scale * self.curtain.peak_reach / self.least).powf(1.0 / RAY_DENSITY_SLOPE)
    }

    /// The area the fragments thrown toward `launch` excavate per unit area at distance
    /// `r`, m, before the skin and the bursts: what `excavated` saturates.
    fn launch_area(&self, r: f64, launch: f64) -> f64 {
        let reach = self.curtain.speed(launch).powi(2);
        self.scale
            * self.curtain.mass(launch)
            * (r / (self.radius * reach)).powf(-RAY_DENSITY_SLOPE)
    }

    /// The area the fragments excavate per unit area at distance `r`, m, and azimuth `phi`
    /// from the centre, before the skin and the bursts (F128); none inside the rim.
    pub fn area(&self, r: f64, phi: f64) -> f64 {
        if r <= self.radius {
            return 0.0;
        }
        let launch = self.curtain.launch_azimuth(phi, (r / self.radius).ln());
        self.launch_area(r, launch)
    }

    /// The share of what the secondaries excavate at distance `r`, m, and azimuth `phi`
    /// that comes from beneath the weathered skin (F212 rays): a secondary digging to a depth d
    /// through a skin h throws out its column, h / d of it the skin's own mature soil, so
    /// that the fresh share is 1 - h / d, "the fresh material diluted by the weathered
    /// material" (Elliott et al. 2018). The depth goes as the secondaries' size, and that
    /// as the speed their fragments left at (F220): it falls with distance as the largest
    /// secondary does (SECONDARY_SIZE_FALL), a faster jet carrying the same fragments
    /// farther by its speed squared, however many fragments the jet throws, so that a sparse
    /// cluster digs as deep as a dense one at the same distance. A jet of the mean speed digs
    /// only the skin at the median ray length (`ray_length`, where Elliott et al.'s median
    /// jet excavates half the ground), so that a ray is bright near its crater and fades to
    /// the ground at its end. `skin` is the share of the Moon's area a thicker skin leaves
    /// (F91), which ends the rays sooner as it lowers the area they excavate.
    pub fn fresh_share(&self, r: f64, phi: f64, skin: f64) -> f64 {
        if r <= self.radius || skin <= 0.0 {
            return 0.0;
        }
        let launch = self.curtain.launch_azimuth(phi, (r / self.radius).ln());
        let end =
            self.length * self.curtain.speed(launch).powi(2) * skin.powf(1.0 / RAY_DENSITY_SLOPE);
        (1.0 - (r / end).powf(SECONDARY_SIZE_FALL)).max(0.0)
    }

    /// The footprint, m across and along the radial, of a cluster of secondaries at
    /// distance `r`, m: RAY_CLUSTER_WIDTH of the largest there wide, RAY_CLUSTER_ASPECT
    /// times as long, and no narrower than `least`.
    pub fn cluster_size(&self, r: f64, least: f64) -> (f64, f64) {
        let w =
            (RAY_CLUSTER_WIDTH * largest_secondary(2.0 * self.radius, r / self.radius)).max(least);
        (w, RAY_CLUSTER_ASPECT * w)
    }

    /// The clusters of secondaries the jets' fragments make between the rim and `reach`, m
    /// (F128): each carries the area its footprint excavates, and they fall at random in
    /// proportion to the area the fragments excavate there, up to RAY_CLUSTER_CAP per
    /// area, each on its jet's bent path, so that on average they excavate what the
    /// fragments do, and a ray far out breaks into its clusters with ground between them.
    /// `least` is the narrowest cluster drawn; past RAY_CLUSTER_MAX the clusters are
    /// larger and fewer; `thin` is the share of the fragments that still land by distance,
    /// m, the rest scattered (F127). The same `seed` gives the same clusters.
    pub fn clusters(
        &self,
        seed: u64,
        reach: f64,
        least: f64,
        thin: &dyn Fn(f64) -> f64,
    ) -> Vec<Cluster> {
        const AZIMUTHS: usize = 360;
        let x_max = (reach / self.radius).max(1.0).ln();
        let bins = ((x_max / 0.04).ceil() as usize).max(1);
        let dx = x_max / bins as f64;
        let dphi = 2.0 * std::f64::consts::PI / AZIMUTHS as f64;
        let mut expected = vec![0.0f64; bins * AZIMUTHS];
        for (i, row) in expected.chunks_mut(AZIMUTHS).enumerate() {
            let r = self.radius * ((i as f64 + 0.5) * dx).exp();
            let (w, l) = self.cluster_size(r, least);
            let bin = r * r * dx * dphi;
            let landing = thin(r);
            for (j, n) in row.iter_mut().enumerate() {
                let a =
                    (landing * self.launch_area(r, (j as f64 + 0.5) * dphi)).min(RAY_CLUSTER_CAP);
                if a >= RAY_CLUSTER_LEAST {
                    *n = a * bin / (w * l);
                }
            }
        }
        let total: f64 = expected.iter().sum();
        let grow = (total / RAY_CLUSTER_MAX as f64).max(1.0).sqrt();
        let rng = Rng::new(seed ^ 0xC1A5_7E25);
        let mut out = Vec::new();
        for (k, &n) in expected.iter().enumerate() {
            if n <= 0.0 {
                continue;
            }
            let key = k as u64;
            let (count, mut draw) = poisson(&rng, key, 0, n / (grow * grow));
            let (i, j) = (k / AZIMUTHS, k % AZIMUTHS);
            for _ in 0..count {
                let x = (i as f64 + rng.uniform(key, draw)) * dx;
                let launch = (j as f64 + rng.uniform(key, draw + 1)) * dphi;
                let size = grow * (0.7 + 0.6 * rng.uniform(key, draw + 2));
                draw += 3;
                let r = self.radius * x.exp();
                let (w, l) = self.cluster_size(r, least);
                out.push(Cluster {
                    r,
                    phi: self.curtain.path_azimuth(launch, x),
                    width: w * size,
                    length: l * size,
                    amount: w * l * size * size,
                });
            }
        }
        out
    }
}

/// The largest secondary, m, a crater `diameter` m across makes at `r_over_radius` crater
/// radii (F128): Croft (1992)'s 0.07 (D / km)^0.87 km at the continuous blanket's edge,
/// falling past it as the distance to -SECONDARY_SIZE_FALL.
pub fn largest_secondary(diameter: f64, r_over_radius: f64) -> f64 {
    let at_edge = SECONDARY_MAX_COEFF_KM * (diameter / 1e3).powf(SECONDARY_MAX_EXPONENT) * 1e3;
    at_edge
        * (r_over_radius / CONTINUOUS_EJECTA_MULTIPLE)
            .max(1.0)
            .powf(-SECONDARY_SIZE_FALL)
}

/// The expected number of secondaries at least `d`, m, across that a crater throws whose
/// largest secondary is `largest`, m (F129): Orientale's census in units of its largest, a
/// Weibull above the census's floor, steepening below it to the slope of Tycho's
/// secondaries.
pub fn secondaries_above(d: f64, largest: f64) -> f64 {
    let (n0, scale, floor) = census();
    let x = d / largest;
    if x >= floor {
        n0 * (-(x / scale).powf(SECONDARY_CENSUS_SHAPE)).exp()
    } else {
        SECONDARY_CENSUS_COUNT * (x / floor).powf(-SECONDARY_SFD_SLOPE)
    }
}

/// The size, m, with `count` secondaries at least as large among those a crater throws
/// whose largest secondary is `largest`, m: the inverse of `secondaries_above`, by which a
/// uniform share of the count draws a secondary's size.
pub fn secondary_size_at(count: f64, largest: f64) -> f64 {
    let (n0, scale, floor) = census();
    let count = count.max(1e-12);
    let x = if count <= SECONDARY_CENSUS_COUNT {
        scale
            * (n0 / count)
                .ln()
                .max(0.0)
                .powf(1.0 / SECONDARY_CENSUS_SHAPE)
    } else {
        floor * (count / SECONDARY_CENSUS_COUNT).powf(-1.0 / SECONDARY_SFD_SLOPE)
    };
    x * largest
}

/// Orientale's census as a Weibull in units of its largest secondary: the count it tends to
/// below, its scale, both from one secondary at the largest and the census's count at its
/// floor, and that floor.
fn census() -> (f64, f64, f64) {
    let k = SECONDARY_CENSUS_SHAPE;
    let floor = SECONDARY_CENSUS_FLOOR_KM / SECONDARY_CENSUS_LARGEST_KM;
    let scale = ((1.0 - floor.powf(k)) / SECONDARY_CENSUS_COUNT.ln()).powf(1.0 / k);
    ((1.0 / scale).powf(k).exp(), scale, floor)
}

/// One cluster of secondaries a ray's fragments make (F128).
#[derive(Clone, Copy, Debug)]
pub struct Cluster {
    /// Distance from the crater's centre, m, and the azimuth it lies at.
    pub r: f64,
    pub phi: f64,
    /// Across and along the radial, m.
    pub width: f64,
    pub length: f64,
    /// The ground it excavates, m².
    pub amount: f64,
}

/// A part of a cluster's ground (F130): where it lies from the cluster's centre, m, along
/// the radial outward and across it; its spread along and across, m; and the ground it
/// excavates, m².
#[derive(Clone, Copy, Debug)]
pub struct Grain {
    pub u: f64,
    pub v: f64,
    pub su: f64,
    pub sv: f64,
    pub amount: f64,
}

impl Cluster {
    /// The cluster's ground laid as its grain (F130). A cluster is RAY_CLUSTER_WIDTH of its
    /// largest secondaries across, and each lies among smaller ones by the same rule at its
    /// own scale: the ground is laid as clusters a RAY_CLUSTER_WIDTH-th as wide, scattered
    /// over the footprint as the cluster lays its ground and sharing it by their own
    /// footprints, and each of those likewise, down to the narrowest drawn, `least`, m.
    /// One too narrow to split is laid whole. The same `seed` gives the same grain.
    pub fn grain(&self, least: f64, seed: u64) -> Vec<Grain> {
        let rng = Rng::new(seed ^ 0x6A41_7C3B);
        let mut out = Vec::new();
        let mut stack = vec![(0.0, 0.0, self.width, self.length, self.amount)];
        let mut draw = 0u64;
        while let Some((u, v, w, l, amount)) = stack.pop() {
            let (su, sv) = (l / 2.5, w / 2.5);
            let ws = w / RAY_CLUSTER_WIDTH;
            if ws < least {
                out.push(Grain {
                    u,
                    v,
                    su,
                    sv,
                    amount,
                });
                continue;
            }
            let ls = RAY_CLUSTER_ASPECT * ws;
            let (n, next) = poisson(&rng, 0, draw, w * l / (ws * ls));
            draw = next;
            let mut parts = Vec::with_capacity(n.max(1) as usize);
            let mut total = 0.0;
            for _ in 0..n.max(1) {
                let du = su * rng.normal(0, draw);
                let dv =
                    sv * (1.0 + RAY_CLUSTER_SPLAY * (du / su).max(0.0)) * rng.normal(0, draw + 2);
                let size = 0.7 + 0.6 * rng.uniform(0, draw + 4);
                draw += 5;
                total += ws * ls * size * size;
                parts.push((u + du, v + dv, ws * size, ls * size));
            }
            for (pu, pv, pw, pl) in parts {
                stack.push((pu, pv, pw, pl, amount * pw * pl / total));
            }
        }
        out
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
    fn a_sparse_cluster_digs_as_deep_as_a_dense_one() {
        // F220: what the secondaries throw out is fresh by their depth, which goes with the
        // distance and the jet's speed, not with how many land there: at the same distance
        // the lightest and the heaviest azimuth of a vertical impact's curtain freshen
        // alike, less with distance, and not at all past the median ray length.
        let radius = 10e3;
        let rays = Rays::new(Curtain::new(3, 0.0, std::f64::consts::FRAC_PI_2), radius);
        let length = ray_length(radius);
        let phis: Vec<f64> = (0..360).map(|k| (k as f64).to_radians()).collect();
        let by_mass = |better: fn(f64, f64) -> bool| {
            phis.iter()
                .copied()
                .reduce(|a, b| {
                    let (ma, mb) = (rays.curtain().mass(a), rays.curtain().mass(b));
                    if better(mb, ma) { b } else { a }
                })
                .unwrap()
        };
        let (heavy, light) = (by_mass(|b, a| b > a), by_mass(|b, a| b < a));
        let r = 0.5 * length;
        let (dense, sparse) = (
            rays.fresh_share(r, heavy, 1.0),
            rays.fresh_share(r, light, 1.0),
        );
        assert!(
            rays.area(r, light) < std::f64::consts::LN_2,
            "the light azimuth excavates under half the ground there"
        );
        assert!(
            dense > 0.2 && (dense - sparse).abs() < 1e-9,
            "{dense} {sparse}"
        );
        assert!(rays.fresh_share(2.0 * radius, heavy, 1.0) > dense);
        assert_eq!(rays.fresh_share(1.01 * length, heavy, 1.0), 0.0);
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

    /// F128: beyond the largest crater measured a crater's rays are shorter in its own
    /// radii the larger it is: an Orientale's reach some 2,000 km, about four of its
    /// radii, where Tycho's median ray is eighteen of its.
    #[test]
    fn a_basin_s_rays_are_shorter_in_its_radii() {
        let edge = RAY_LENGTH_MAX_RADIUS_KM * 1e3;
        let ratio = |r: f64| ray_length(r) / r;
        assert!((ratio(edge) - 18.2).abs() < 0.2, "{}", ratio(edge));
        assert!(ratio(0.5 * edge) < ratio(edge) && ratio(4.0 * edge) < ratio(edge));
        let orientale = ray_length(465e3) / 1e3;
        assert!((1_800.0..2_200.0).contains(&orientale), "{orientale}");
        // At its own scaled distance a basin's rays are thinner than a large crater's.
        let curtain = || Curtain::new(7, 0.0, std::f64::consts::FRAC_PI_2);
        let (crater, basin) = (Rays::new(curtain(), edge), Rays::new(curtain(), 4.0 * edge));
        assert!(basin.area(10.0 * 4.0 * edge, 0.3) < crater.area(10.0 * edge, 0.3));
    }

    /// F128: a ray is the clusters of its secondaries. They carry what the fragments
    /// excavate, up to RAY_CLUSTER_CAP per area; far out their footprints leave most of
    /// the ground between them; the same seed gives the same clusters.
    #[test]
    fn a_ray_is_the_clusters_of_its_secondaries() {
        let radius = 40e3;
        let rays = Rays::new(Curtain::new(11, 0.0, std::f64::consts::FRAC_PI_2), radius);
        let reach = 60.0 * radius;
        let clusters = rays.clusters(11, reach, 1e3, &|_| 1.0);
        let again = rays.clusters(11, reach, 1e3, &|_| 1.0);
        assert_eq!(clusters.len(), again.len());
        assert!(
            clusters
                .iter()
                .zip(&again)
                .all(|(a, b)| a.r == b.r && a.phi == b.phi)
        );
        let (lo, hi) = (10.0 * radius, 20.0 * radius);
        let carried: f64 = clusters
            .iter()
            .filter(|c| c.r >= lo && c.r < hi)
            .map(|c| c.amount)
            .sum();
        let (n_r, n_phi) = (200, 720);
        let (dr, dphi) = (
            (hi - lo) / n_r as f64,
            2.0 * std::f64::consts::PI / n_phi as f64,
        );
        let mut expected = 0.0;
        for i in 0..n_r {
            let r = lo + (i as f64 + 0.5) * dr;
            for j in 0..n_phi {
                let phi = (j as f64 + 0.5) * dphi;
                expected += rays.area(r, phi).min(RAY_CLUSTER_CAP) * r * dr * dphi;
            }
        }
        assert!(
            (carried / expected - 1.0).abs() < 0.15,
            "carried {carried:.3e} against {expected:.3e}"
        );
        let far: f64 = clusters
            .iter()
            .filter(|c| c.r > 30.0 * radius)
            .map(|c| c.width * c.length)
            .sum();
        let ring = std::f64::consts::PI * ((60.0 * radius).powi(2) - (30.0 * radius).powi(2));
        assert!(
            far < 0.2 * ring,
            "far clusters cover {:.3} of their ring",
            far / ring
        );
    }

    /// F130: a cluster's grain shares its ground among parts narrower than the narrowest
    /// that would split, scattered over its footprint; the same seed gives the same grain;
    /// one too narrow to split is laid whole.
    #[test]
    fn a_cluster_s_grain_shares_its_ground() {
        let c = Cluster {
            r: 1e6,
            phi: 0.3,
            width: 30e3,
            length: 90e3,
            amount: 27e8,
        };
        let least = 1.6e3;
        let g = c.grain(least, 7);
        let again = c.grain(least, 7);
        assert_eq!(g.len(), again.len());
        assert!(
            g.iter()
                .zip(&again)
                .all(|(a, b)| a.u == b.u && a.amount == b.amount)
        );
        let total: f64 = g.iter().map(|p| p.amount).sum();
        assert!(
            (total / c.amount - 1.0).abs() < 1e-9,
            "{total} of {}",
            c.amount
        );
        assert!(g.len() > 20, "{} parts", g.len());
        assert!(g.iter().all(|p| 2.5 * p.sv < RAY_CLUSTER_WIDTH * least));
        let (su, sv) = (c.length / 2.5, c.width / 2.5);
        let inside: f64 = g
            .iter()
            .filter(|p| {
                p.u.abs() < 3.0 * su && p.v.abs() < 3.0 * sv * (1.0 + 3.0 * RAY_CLUSTER_SPLAY)
            })
            .map(|p| p.amount)
            .sum();
        assert!(inside > 0.95 * c.amount, "{:.3} inside", inside / c.amount);
        let narrow = Cluster {
            width: 4e3,
            length: 12e3,
            amount: 48e6,
            ..c
        };
        let whole = narrow.grain(least, 7);
        assert_eq!(whole.len(), 1);
        assert_eq!(whole[0].amount, narrow.amount);
    }

    /// F129: a crater's secondaries follow Orientale's census in units of its largest: 2,728
    /// of 2 km and more, one of 27 km and some 400 of 10 km and more at Orientale's own; a
    /// handful over 2.7 km for Tycho; and the draw inverts the count.
    #[test]
    fn the_secondaries_follow_orientales_census() {
        let largest = SECONDARY_CENSUS_LARGEST_KM * 1e3;
        let at = |km: f64| secondaries_above(km * 1e3, largest);
        assert!(
            (at(2.0) / SECONDARY_CENSUS_COUNT - 1.0).abs() < 1e-9,
            "{}",
            at(2.0)
        );
        assert!((at(27.0) - 1.0).abs() < 1e-9, "{}", at(27.0));
        assert!(at(10.0) > 300.0 && at(10.0) < 500.0, "{}", at(10.0));
        assert!((at(1.999_999) / at(2.000_001) - 1.0).abs() < 1e-4);
        assert!(at(1.0) > 3.0 * at(2.0) && at(5.0) < at(4.0));
        // Croft's largest is the census's at Orientale's 930 km.
        let orientale = largest_secondary(930e3, CONTINUOUS_EJECTA_MULTIPLE);
        assert!((orientale / largest - 1.0).abs() < 0.02, "{orientale}");
        let tycho = secondaries_above(2.7e3, largest_secondary(85e3, CONTINUOUS_EJECTA_MULTIPLE));
        assert!(tycho > 3.0 && tycho < 15.0, "{tycho}");
        for count in [0.5, 1.0, 37.0, 900.0, 2728.0, 5000.0, 1e5] {
            let d = secondary_size_at(count, largest);
            let back = secondaries_above(d, largest);
            assert!(
                (back / count - 1.0).abs() < 1e-6,
                "{count}: {d} gives {back}"
            );
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
