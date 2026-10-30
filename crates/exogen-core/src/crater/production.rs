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
//! A body is struck by two populations (F135): the rocky leftovers of the inner system,
//! whose craters follow the Neukum production function and decline with the Neukum
//! chronology, and the icy planetesimals scattered from beyond the giant planets, whose
//! numbers by size follow their collisional history, whose craters follow each body's own
//! scaling, and whose flux declines as a stretched exponential from formation and leaks on
//! at a steady rate (F136). A moon shattered by an impact it could not hold is struck
//! besides by that impact's debris. A `Bombardment` holds a body's sources, each with its
//! own production and history.
//!
//! Saturation follows the standard picture: new craters obliterate old ones, so a surface
//! approaches an equilibrium density that is a fixed fraction of geometric saturation. Per
//! size bin the surviving density is `n_eq (1 - exp(-n_prod / n_eq))`, and the formation
//! times of survivors follow a truncated exponential in accumulated production, which is
//! sampled directly so that old surfaces carry mostly young survivors at small sizes and
//! old craters at large sizes, as observed.

use rayon::prelude::*;

use super::morphology::{elliptical_threshold, ellipticity_for, impact_angle};
use crate::constants::*;
use crate::fields::Field;
use crate::grid::CubeSphere;
use crate::rng::Rng;

/// A production function: polynomial size distribution plus chronology coefficients, or
/// the craters a population of impactors makes on the body (F135, F136).
#[derive(Clone, Copy, Debug)]
pub struct ProductionFunction {
    pub coeffs: &'static [f64; 12],
    /// N(>= 1 km) = a (exp(b T) - 1) + c T per km², T in Gyr.
    pub chronology: (f64, f64, f64),
    /// A population of impactors and the craters they make, in place of the polynomial.
    pub impactors: Option<ImpactorCraters>,
}

/// How a population's members are spread in size (F136).
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ImpactorSizes {
    /// The destabilized Kuiper belt's comets (ICY_IMPACTOR_SFD), counted against those at
    /// or above 1 km.
    Comets,
    /// A shattered moon's fragments (DEBRIS_SLOPES), counted against those at or above
    /// DEBRIS_BREAK_KM.
    Debris,
}

impl ImpactorSizes {
    /// Members at or above a size, km, per one at or above the reference size, and the
    /// cumulative slope there.
    pub fn share(&self, d_km: f64) -> (f64, f64) {
        let x = d_km.max(1e-12).log10();
        match self {
            ImpactorSizes::Comets => {
                let t = &ICY_IMPACTOR_SFD;
                let (head, tail) = ICY_IMPACTOR_END_SLOPES;
                let n = t.len();
                if x <= t[0].0 {
                    return (10f64.powf(t[0].1 - head * (x - t[0].0)), head);
                }
                if x >= t[n - 1].0 {
                    return (10f64.powf(t[n - 1].1 - tail * (x - t[n - 1].0)), tail);
                }
                let k = t.partition_point(|p| p.0 <= x).clamp(1, n - 1);
                let (a, b) = (t[k - 1], t[k]);
                let slope = -(b.1 - a.1) / (b.0 - a.0);
                (10f64.powf(a.1 - slope * (x - a.0)), slope)
            }
            ImpactorSizes::Debris => {
                let (small, large) = DEBRIS_SLOPES;
                let r = d_km.max(1e-12) / DEBRIS_BREAK_KM;
                if r < 1.0 {
                    (r.powf(-small), small)
                } else {
                    (r.powf(-large), large)
                }
            }
        }
    }
}

/// What a crater forms in, and how fast it is struck (F136).
#[derive(Clone, Copy, Debug)]
pub struct Target {
    /// m/s².
    pub gravity: f64,
    /// m/s.
    pub speed: f64,
    /// Pa.
    pub strength: f64,
    /// kg/m³.
    pub density: f64,
    /// The simple-to-complex transition, m.
    pub simple_complex: f64,
}

impl Target {
    /// The transient crater an impactor of diameter `d` m makes, m (CRATER_PI_SCALING).
    pub fn transient(&self, d: f64) -> f64 {
        let (k, mu) = CRATER_PI_SCALING;
        let v2 = self.speed.max(1.0).powi(2);
        let gravity = self.gravity * d / (2.0 * v2);
        let strength = (self.strength / (self.density.max(1.0) * v2)).powf(0.5 * (2.0 + mu));
        k * d * (gravity + strength).powf(-mu / (2.0 + mu))
    }

    /// The final crater an impactor of diameter `d` m makes, m: the transient cavity's
    /// collapse, the inverse of `morphology::transient_diameter`.
    pub fn crater(&self, d: f64) -> f64 {
        let (k, p) = TRANSIENT_COMPLEX;
        let tr = self.transient(d) / 1e3;
        let d_sc = self.simple_complex.max(1.0) / 1e3;
        let tr_sc = d_sc * k.powf(-1.0 / p);
        let km = if tr < tr_sc {
            tr * k.powf(1.0 / p)
        } else {
            k * tr.powf(p) * d_sc.powf(1.0 - p)
        };
        km * 1e3
    }
}

/// The impact angle below which a body's craters come out elliptical, by their diameter
/// (F167): from the cratering efficiency, the transient cavity a vertical impact makes over
/// the impactor's diameter, read from the body's crater scaling at its mean impact speed
/// for every population that strikes it (a stated choice), tabulated by the crater's
/// diameter.
#[derive(Clone, Copy, Debug)]
pub struct Oblique {
    /// The threshold angle, radians, for a crater of log10 diameter
    /// `MAP_START + k * MAP_STEP`, km.
    map: [f64; MAP_POINTS],
}

impl Oblique {
    pub fn new(target: &Target) -> Self {
        let mut map = [0.0; MAP_POINTS];
        for (k, m) in map.iter_mut().enumerate() {
            let crater = 10f64.powf(MAP_START + MAP_STEP * k as f64) * 1e3;
            let (mut lo, mut hi) = (1e-4f64.ln(), 1e8f64.ln());
            for _ in 0..80 {
                let mid = 0.5 * (lo + hi);
                if target.crater(mid.exp()) < crater {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            let impactor = (0.5 * (lo + hi)).exp();
            *m = elliptical_threshold(target.transient(impactor) / impactor);
        }
        Self { map }
    }

    /// The same threshold for every crater: `threshold` radians.
    pub fn fixed(threshold: f64) -> Self {
        Self {
            map: [threshold; MAP_POINTS],
        }
    }

    /// The threshold angle, radians, for a crater of `diameter` m.
    pub fn threshold(&self, diameter: f64) -> f64 {
        let x = ((diameter.max(1e-9) / 1e3).log10() - MAP_START) / MAP_STEP;
        let x = x.clamp(0.0, (MAP_POINTS - 1) as f64);
        let i = (x.floor() as usize).min(MAP_POINTS - 2);
        let f = x - i as f64;
        self.map[i] * (1.0 - f) + self.map[i + 1] * f
    }

    /// The elongation of a crater of `diameter` m struck at `angle` radians from horizontal.
    pub fn ellipticity(&self, diameter: f64, angle: f64) -> f64 {
        ellipticity_for(angle, self.threshold(diameter))
    }
}

/// The craters a population of impactors makes on a body (F135, F136): how many of its
/// members strike, by size, and the crater each makes there.
#[derive(Clone, Copy, Debug)]
pub struct ImpactorCraters {
    /// Comets: those at or above 1 km per km² per Gyr of the steady flux. Debris: the
    /// fragments at or above the break per km² over the whole rain.
    pub per_gyr: f64,
    pub sizes: ImpactorSizes,
    /// The crater an impactor of 1 km makes, km.
    pub crater_1km: f64,
    /// log10 of the impactor, km, that makes a crater of log10 diameter
    /// `MAP_START + k * MAP_STEP`, km.
    map: [f64; MAP_POINTS],
}

const MAP_START: f64 = -3.0;
const MAP_STEP: f64 = 0.125;
const MAP_POINTS: usize = 65;

impl ImpactorCraters {
    pub fn new(per_gyr: f64, sizes: ImpactorSizes, target: &Target) -> Self {
        let mut map = [0.0; MAP_POINTS];
        for (k, m) in map.iter_mut().enumerate() {
            let crater = 10f64.powf(MAP_START + MAP_STEP * k as f64) * 1e3;
            let (mut lo, mut hi) = (1e-4f64.ln(), 1e8f64.ln());
            for _ in 0..80 {
                let mid = 0.5 * (lo + hi);
                if target.crater(mid.exp()) < crater {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            *m = (0.5 * (lo + hi)).exp().log10() - 3.0;
        }
        Self {
            per_gyr,
            sizes,
            crater_1km: target.crater(1e3) / 1e3,
            map,
        }
    }

    /// The map's segment for a crater and the position along it.
    fn segment(d_km: f64) -> (usize, f64) {
        let x = (d_km.max(1e-12).log10() - MAP_START) / MAP_STEP;
        let i = x.floor().clamp(0.0, (MAP_POINTS - 2) as f64) as usize;
        (i, x - i as f64)
    }

    /// The impactor that makes a crater, km.
    pub fn impactor(&self, d_km: f64) -> f64 {
        let (i, f) = Self::segment(d_km);
        10f64.powf(self.map[i] + f * (self.map[i + 1] - self.map[i]))
    }

    /// Craters at or above a diameter, km, per km² per Gyr of the steady flux (comets) or
    /// over the whole rain (debris).
    pub fn cumulative(&self, d_km: f64) -> f64 {
        self.per_gyr * self.sizes.share(self.impactor(d_km)).0
    }

    /// d log N / d log D at a diameter.
    pub fn log_slope(&self, d_km: f64) -> f64 {
        let (i, _) = Self::segment(d_km);
        let grows = (self.map[i + 1] - self.map[i]) / MAP_STEP;
        -self.sizes.share(self.impactor(d_km)).1 * grows
    }
}

impl ProductionFunction {
    pub fn lunar() -> Self {
        Self {
            coeffs: &NPF_LUNAR,
            chronology: (NEUKUM_A, NEUKUM_B, NEUKUM_C),
            impactors: None,
        }
    }

    pub fn mars() -> Self {
        Self {
            coeffs: &NPF_MARS,
            chronology: MARS_CHRONOLOGY,
            impactors: None,
        }
    }

    /// A population of impactors on a body (F135, F136); the chronology is unused, the
    /// source's history being its own.
    pub fn from_impactors(impactors: ImpactorCraters) -> Self {
        Self {
            coeffs: &NPF_LUNAR,
            chronology: (0.0, 1.0, impactors.per_gyr),
            impactors: Some(impactors),
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
        if let Some(imp) = &self.impactors {
            return imp.cumulative(d_km);
        }
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
        if let Some(imp) = &self.impactors {
            return imp.log_slope(d_km);
        }
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
        self.sample_with(|d| self.cumulative_1gyr(d), u, d_min_km, d_max_km)
    }

    /// The steady flux's production: the lunar function, falling steeper above the knee
    /// where the asteroid belt runs out of large bodies.
    pub fn cumulative_steady_1gyr(&self, d_km: f64) -> f64 {
        let n = self.cumulative_1gyr(d_km);
        // The outer system's steady leak keeps the scattered disc's own distribution.
        if self.impactors.is_none() && d_km > STEADY_FLUX_KNEE_KM {
            n * (d_km / STEADY_FLUX_KNEE_KM).powf(-STEADY_FLUX_EXTRA_SLOPE)
        } else {
            n
        }
    }

    /// A diameter from the steady flux's production, km, capped at the largest crater a
    /// stray body makes.
    pub fn sample_steady_diameter_above(&self, u: f64, d_min_km: f64, d_max_km: f64) -> f64 {
        let d_max_km = if self.impactors.is_some() {
            d_max_km
        } else {
            d_max_km.min(STEADY_FLUX_MAX_KM).max(d_min_km)
        };
        self.sample_with(|d| self.cumulative_steady_1gyr(d), u, d_min_km, d_max_km)
    }

    fn sample_with(
        &self,
        cumulative: impl Fn(f64) -> f64,
        u: f64,
        d_min_km: f64,
        d_max_km: f64,
    ) -> f64 {
        let n_min = cumulative(d_min_km);
        let n_max = cumulative(d_max_km);
        // Target cumulative value between n_max and n_min; larger u gives a smaller crater.
        let target = n_max + (n_min - n_max) * u.clamp(0.0, 1.0);
        let (mut lo, mut hi) = (d_min_km.ln(), d_max_km.ln());
        for _ in 0..50 {
            let mid = 0.5 * (lo + hi);
            if cumulative(mid.exp()) > target {
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
    /// The weight of the steady flux: 1 for a population that leaks on after its early
    /// decline, 0 for a moon's debris, which comes all at once (F136).
    pub late: f64,
    /// The early component, declining from formation.
    pub early: Early,
    /// Time since formation when the surface begins retaining craters, Gyr.
    pub record_start: f64,
}

/// How a population's early flux declines from formation (F31, F135).
#[derive(Clone, Copy, Debug)]
pub enum Early {
    /// The rocky leftovers of the inner system: an early rate `ratio` times the steady one
    /// at formation, e-folding over `decay` Gyr (the Neukum chronology).
    Exponential { ratio: f64, decay: f64 },
    /// The outer system's scattered planetesimals: `total` Gyr of steady flux in all, of
    /// which the share exp(-(t / tau)^beta) is still to come t Gyr after formation (a
    /// stretched exponential, Wong et al. 2023).
    Stretched { total: f64, beta: f64, tau: f64 },
    /// A shattered moon's debris: `total` units from `start`, falling off over `tau` Gyr
    /// (F136).
    Burst { total: f64, start: f64, tau: f64 },
}

impl Early {
    /// The early rate at a time in units of the steady rate.
    fn rate(&self, tau: f64) -> f64 {
        match *self {
            Early::Exponential { ratio, decay } => ratio * (-tau / decay.max(1e-12)).exp(),
            Early::Stretched {
                total,
                beta,
                tau: scale,
            } => {
                let t = tau.max(1e-9);
                let x = (t / scale.max(1e-12)).powf(beta);
                total * beta * x * (-x).exp() / t
            }
            Early::Burst {
                total,
                start,
                tau: scale,
            } => {
                if tau < start {
                    0.0
                } else {
                    total * (-(tau - start) / scale.max(1e-12)).exp() / scale.max(1e-12)
                }
            }
        }
    }

    /// The early production between two times in Gyr of steady flux.
    fn between(&self, a: f64, b: f64) -> f64 {
        match *self {
            Early::Exponential { ratio, decay } => {
                ratio * decay * ((-a / decay.max(1e-12)).exp() - (-b / decay.max(1e-12)).exp())
            }
            Early::Stretched { total, beta, tau } => {
                let left = |t: f64| (-(t.max(0.0) / tau.max(1e-12)).powf(beta)).exp();
                total * (left(a) - left(b))
            }
            Early::Burst { total, start, tau } => {
                let left = |t: f64| (-(t - start).max(0.0) / tau.max(1e-12)).exp();
                total * (left(a) - left(b))
            }
        }
    }
}

impl FluxHistory {
    /// From chronology coefficients (a, b, c) of N1(T) = a (exp(b T) - 1) + c T, with T in
    /// Gyr before present.
    pub fn from_chronology(chronology: (f64, f64, f64), record_start: f64) -> Self {
        let (_, b, c) = chronology;
        Self {
            steady: c,
            late: 1.0,
            early: Early::Exponential {
                ratio: early_flux_ratio(chronology),
                decay: 1.0 / b,
            },
            record_start,
        }
    }

    /// Rate of craters >= 1 km per km² per Gyr at `tau` Gyr since formation.
    pub fn rate1(&self, tau: f64) -> f64 {
        if tau < self.record_start {
            return 0.0;
        }
        self.steady * (self.late + self.early.rate(tau))
    }

    /// Share of the flux at a time that is the early, basin-forming component.
    pub fn early_share(&self, tau: f64) -> f64 {
        let e = self.early.rate(tau);
        if self.late + e <= 0.0 {
            return 1.0;
        }
        e / (self.late + e)
    }

    /// N(>= 1 km) per km² formed between `tau0` and `tau1` Gyr since formation.
    pub fn n1_between(&self, tau0: f64, tau1: f64) -> f64 {
        let a = tau0.max(self.record_start);
        let b = tau1.max(a);
        self.steady * (self.late * (b - a) + self.early.between(a, b))
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

/// One population striking a body (F135): its craters' production, its history, and its
/// flux relative to the production function, the density multiplier included.
#[derive(Clone, Copy, Debug)]
pub struct Source {
    pub pf: ProductionFunction,
    pub flux: FluxHistory,
    pub scale: f64,
}

/// The populations striking a body (F135).
#[derive(Clone, Debug, Default)]
pub struct Bombardment {
    pub sources: Vec<Source>,
}

impl Bombardment {
    /// One source alone.
    pub fn single(pf: ProductionFunction, flux: FluxHistory, scale: f64) -> Self {
        Self {
            sources: vec![Source { pf, flux, scale }],
        }
    }

    /// Craters at or above a diameter, km, per km² per Gyr of every source's steady flux,
    /// or over a moon's whole debris rain: the shape of the production by size.
    pub fn cumulative_1gyr(&self, d_km: f64) -> f64 {
        self.sources
            .iter()
            .map(|s| s.scale * s.pf.cumulative_1gyr(d_km))
            .sum()
    }

    /// Craters at or above a diameter, km, per km² per Gyr of the steady flux the sources
    /// leak on with: the production today, a moon's debris long gone (F136).
    pub fn steady_1gyr(&self, d_km: f64) -> f64 {
        self.sources
            .iter()
            .map(|s| s.scale * s.flux.late * s.pf.cumulative_1gyr(d_km))
            .sum()
    }

    /// Craters at or above a diameter, km, per km² formed between two times.
    pub fn density_above(&self, d_km: f64, tau0: f64, tau1: f64) -> f64 {
        self.sources
            .iter()
            .map(|s| s.scale * s.pf.cumulative_1gyr(d_km) * s.flux.equivalent_gyr(tau0, tau1))
            .sum()
    }

    /// Gyr of steady flux equivalent to the production between two times, each source
    /// weighted by its craters of 1 km: the clock of the gardening the small craters do.
    pub fn equivalent_gyr(&self, tau0: f64, tau1: f64) -> f64 {
        let steady = self.steady_1gyr(1.0);
        if steady <= 0.0 {
            return 0.0;
        }
        self.density_above(1.0, tau0, tau1) / steady
    }

    /// Expected number of craters at or above `d_min` (metres) over the whole body formed
    /// between `t0` and `t1` Gyr since formation, from every source.
    pub fn expected_count_above(&self, body_radius: f64, d_min: f64, t0: f64, t1: f64) -> f64 {
        self.sources
            .iter()
            .map(|s| expected_count_above(&s.pf, &s.flux, s.scale, body_radius, d_min, t0, t1))
            .sum()
    }

    /// The latest record start of the sources, Gyr since formation.
    pub fn record_start(&self) -> f64 {
        self.sources
            .iter()
            .map(|s| s.flux.record_start)
            .fold(0.0, f64::max)
    }

    /// The time since formation, before `present`, after which `remaining` craters per km²
    /// of a bin still form, given each source's craters per km² per Gyr of its steady flux
    /// in that bin; the record's start when that is the whole record.
    pub fn time_with_remaining(&self, present: f64, remaining: f64, bin: &[f64]) -> f64 {
        let after = |t: f64| -> f64 {
            self.sources
                .iter()
                .zip(bin)
                .map(|(s, &b)| b * s.flux.equivalent_gyr(t, present))
                .sum()
        };
        let (mut lo, mut hi) = (self.record_start().min(present), present);
        if remaining >= after(lo) {
            return lo;
        }
        if remaining <= 0.0 {
            return hi;
        }
        for _ in 0..60 {
            let mid = 0.5 * (lo + hi);
            if after(mid) > remaining {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    }
}

/// The scalars a body's bombardment is built from (F135, F136), as the parameter graph
/// derives them: the rocky leftovers' production, flux relative to the Moon's and early
/// decline; the outer system's comets, their early total and timescale; a shattered moon's
/// debris; when the record starts; the density multiplier.
#[derive(Clone, Copy, Debug)]
pub struct BombardmentParams {
    pub rocky: ProductionFunction,
    pub rocky_scale: f64,
    pub rocky_ratio: f64,
    pub rocky_decay: f64,
    pub comets: Option<ImpactorCraters>,
    pub comet_total: f64,
    pub comet_tau: f64,
    pub debris: Option<ImpactorCraters>,
    pub record_start: f64,
    pub density: f64,
}

impl BombardmentParams {
    /// From the graph's nodes, read by `f64_of`: the rocky leftovers' production function
    /// by name, how far its target is ice (audit S9), when the record starts, and whether a moon's
    /// debris is counted (not when the record's start is being found, which it follows).
    pub fn from_nodes(
        f64_of: impl Fn(&str) -> f64,
        production: &str,
        target_ice: f64,
        record_start: f64,
        with_debris: bool,
    ) -> Self {
        let target = |speed: f64| Target {
            gravity: f64_of("body.gravity"),
            speed,
            strength: super::by_ice(TARGET_STRENGTH_ROCK, TARGET_STRENGTH_ICE, target_ice),
            density: f64_of("cratering.target_density"),
            simple_complex: f64_of("cratering.simple_complex_transition"),
        };
        let comets = f64_of("cratering.icy_flux");
        let fragments = if with_debris {
            f64_of("cratering.debris_fragments")
        } else {
            0.0
        };
        Self {
            rocky: ProductionFunction::by_name(production),
            rocky_scale: f64_of("cratering.flux_relative_to_moon"),
            rocky_ratio: f64_of("cratering.early_flux_ratio"),
            rocky_decay: f64_of("cratering.early_flux_decay"),
            comets: (comets > 0.0).then(|| {
                ImpactorCraters::new(
                    comets,
                    ImpactorSizes::Comets,
                    &target(f64_of("cratering.impact_velocity")),
                )
            }),
            comet_total: f64_of("cratering.icy_early_total"),
            comet_tau: f64_of("cratering.icy_early_timescale"),
            debris: (fragments > 0.0).then(|| {
                ImpactorCraters::new(
                    fragments,
                    ImpactorSizes::Debris,
                    &target(f64_of("cratering.debris_velocity")),
                )
            }),
            record_start,
            density: f64_of("multiplier.crater_density"),
        }
    }

    pub fn build(&self) -> Bombardment {
        let mut sources = Vec::new();
        if self.rocky_scale > 0.0 {
            sources.push(Source {
                pf: self.rocky,
                flux: FluxHistory {
                    steady: self.rocky.chronology.2,
                    late: 1.0,
                    early: Early::Exponential {
                        ratio: self.rocky_ratio,
                        decay: self.rocky_decay,
                    },
                    record_start: self.record_start,
                },
                scale: self.rocky_scale * self.density,
            });
        }
        if let Some(comets) = self.comets.filter(|c| c.per_gyr > 0.0) {
            let pf = ProductionFunction::from_impactors(comets);
            sources.push(Source {
                pf,
                flux: FluxHistory {
                    steady: pf.cumulative_1gyr(1.0),
                    late: 1.0,
                    early: Early::Stretched {
                        total: self.comet_total,
                        beta: ICY_EARLY_BETA,
                        tau: self.comet_tau,
                    },
                    record_start: self.record_start,
                },
                scale: self.density,
            });
        }
        if let Some(debris) = self.debris.filter(|d| d.per_gyr > 0.0) {
            sources.push(Source {
                pf: ProductionFunction::from_impactors(debris),
                flux: FluxHistory {
                    steady: 1.0,
                    late: 0.0,
                    early: Early::Burst {
                        total: 1.0,
                        start: self.record_start,
                        tau: DEBRIS_DECAY_GYR,
                    },
                    record_start: self.record_start,
                },
                scale: self.density,
            });
        }
        Bombardment { sources }
    }
}

/// The diameter of the largest crater a body holds, m (F135): an icy body's radius, a rocky
/// body's diameter, a body of rock and ice between them by its target's ice degree (audit
/// S9). A larger impact shatters and resurfaces it.
pub fn hold_diameter(radius: f64, ice: f64) -> f64 {
    2.0 * radius * super::by_ice(ROCK_HOLD_RADIUS_FRACTION, ICE_HOLD_RADIUS_FRACTION, ice)
}

/// Inputs to the statistical population sampler.
#[derive(Clone, Debug)]
pub struct PopulationConfig {
    /// The populations striking the body, their fluxes including the density multiplier.
    pub bombardment: Bombardment,
    pub body_radius: f64,
    /// Smallest and largest diameters sampled here, metres.
    pub d_min: f64,
    pub d_max: f64,
    pub bins_per_octave: u32,
    /// Fraction of geometric saturation at equilibrium.
    pub saturation_fraction: f64,
    /// Longest effective surface age for crater retention, Gyr.
    pub retention_age: f64,
    /// The threshold angles of elliptical craters (F167).
    pub oblique: Oblique,
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
    // Per-bin production density per km² for a 1 Gyr surface of each source's steady flux,
    // and equilibrium density.
    let sources = &cfg.bombardment.sources;
    let bin_prod_1gyr: Vec<Vec<f64>> = bins
        .iter()
        .map(|(lo, hi, _)| {
            sources
                .iter()
                .map(|s| {
                    s.scale
                        * (s.pf.cumulative_1gyr(lo / 1e3) - s.pf.cumulative_1gyr(hi / 1e3)).max(0.0)
                })
                .collect()
        })
        .collect();
    // The bins' local slopes, from the sum of the sources' productions.
    let bin_slope: Vec<f64> = bins
        .iter()
        .map(|(lo, hi, _)| {
            let n_lo = cfg.bombardment.cumulative_1gyr(lo / 1e3);
            let n_hi = cfg.bombardment.cumulative_1gyr(hi / 1e3);
            if n_lo > 0.0 && n_hi > 0.0 && n_lo > n_hi {
                (n_lo / n_hi).ln() / (hi / lo).ln()
            } else {
                2.0
            }
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
            let exposure: Vec<f64> = sources
                .iter()
                .map(|s| s.flux.equivalent_gyr(start, present))
                .collect();
            if exposure.iter().all(|&e| e <= 0.0) {
                return Vec::new();
            }
            let (a0, b0) = grid.angles(cell);
            let mut out = Vec::new();
            for (b, (lo, hi, _mid)) in bins.iter().enumerate() {
                let n_prod: f64 = bin_prod_1gyr[b]
                    .iter()
                    .zip(&exposure)
                    .map(|(p, e)| p * e)
                    .sum();
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
                    let slope = bin_slope[b].max(0.5);
                    let q = (hi / lo).powf(-slope);
                    let diameter = lo * (1.0 - u_d * (1.0 - q)).powf(-1.0 / slope);
                    // Accumulated production after formation, truncated exponential.
                    let f = if n_eq > 0.0 {
                        -n_eq * (1.0 - u_t * (1.0 - (-n_prod / n_eq).exp())).ln()
                    } else {
                        u_t * n_prod
                    };
                    // The craters of the bin still to come after formation, f per km², and
                    // the time after which they come, from every source's history.
                    let time_gyr = cfg
                        .bombardment
                        .time_with_remaining(present, f, &bin_prod_1gyr[b])
                        .max(start);
                    let angle = impact_angle(u_ob);
                    let azimuth = u_az * 2.0 * std::f64::consts::PI;
                    let ellipticity = cfg.oblique.ellipticity(diameter, angle);
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
            bombardment: Bombardment::single(
                ProductionFunction::lunar(),
                ProductionFunction::lunar().flux_history(0.3),
                1.0,
            ),
            body_radius: radius,
            d_min: 40_000.0,
            d_max: 320_000.0,
            bins_per_octave: 4,
            saturation_fraction: EMPIRICAL_SATURATION_FRACTION,
            retention_age: 10.0,
            oblique: Oblique::fixed(0.0),
            seed: 11,
        };
        let area_km2 = 4.0 * std::f64::consts::PI * radius * radius / 1e6;
        for (surface_time, present, check_saturation) in [(1.2f64, 4.5f64, false), (0.0, 4.5, true)]
        {
            let ages = Field::filled(grid, surface_time as f32);
            let pop = sample_population(&grid, &ages, present, &cfg);
            let source = &cfg.bombardment.sources[0];
            let exposure = source.flux.equivalent_gyr(surface_time, present);
            for d in [40.0, 80.0, 160.0] {
                let count = pop.iter().filter(|c| c.diameter >= d * 1e3).count() as f64;
                let prod = (source.pf.cumulative_1gyr(d) - source.pf.cumulative_1gyr(320.0))
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
