//! The surface layer's stability (M63).
//!
//! The turbulence that carries heat and vapour between a ground and the air over it grows
//! where the ground warms its air from below and dies away where the air stands the warmer,
//! an inversion. The bulk transfer coefficient over its neutral is read from the bulk
//! Richardson number between the ground and the air at the wind's reference height, Ri_b =
//! g z (theta_air - theta_ground) / (T U²), by Monin and Obukhov's similarity with the
//! profile functions measured over flat ground. Where the ground is the warmer, the Dyer and
//! Hicks forms, phi_M = (1 - 16 zeta)^(-1/4) and phi_H = (1 - 16 zeta)^(-1/2), fitted to
//! the Kansas measurements (Paulson 1970, Dyer 1974, Högström 1988). Where the air is the
//! warmer, Holtslag and De Bruin's (1988) and Beljaars and Holtslag's (1991), fitted to the
//! Cabauw mast's measurements up to a zeta of 10, linear near neutral as Webb's (1970) and
//! with no critical Richardson number past which the turbulence stops: the exchange falls by
//! degrees however strong the inversion. The linear forms alone stop it dead at a bulk
//! Richardson number of 0.2 (Grachev and Fairall 1997).
//!
//! Heat and momentum share one roughness, that of the stated neutral coefficient's log
//! profile, C_N = (kappa / ln(z / z0))²: with Lambda = ln(z / z0) and zeta0 = zeta z0 / z,
//!
//! C_H / C_N = Lambda² / ((Lambda - psi_M(zeta) + psi_M(zeta0)) (Lambda - psi_H(zeta) +
//! psi_H(zeta0))), Ri_b = zeta (Lambda - psi_H(zeta) + psi_H(zeta0)) / (Lambda - psi_M(zeta)
//! + psi_M(zeta0))²,
//!
//! read the other way, from the Richardson number, through a table over zeta, along which
//! the Richardson number rises.

use crate::constants::*;
use std::sync::OnceLock;

/// What a body's surface layer reads of it (M63): its gravity, m/s², and its air's heat
/// capacity, J/kg/K.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SurfaceLayer {
    pub gravity: f64,
    pub heat_capacity: f64,
}

impl SurfaceLayer {
    /// The bulk Richardson number between a ground at `ground` K and its air at `air` K at the
    /// wind's reference height, under a wind of `wind` m/s: between their potential
    /// temperatures, the air's the warmer by the dry adiabat's g z / c_p, over the air's
    /// temperature.
    pub fn richardson(&self, ground: f64, air: f64, wind: f64) -> f64 {
        let z = WIND_REFERENCE_HEIGHT;
        let lift = air - ground + self.gravity / self.heat_capacity * z;
        self.gravity * z * lift / (air.max(1.0) * (wind * wind).max(1e-300))
    }

    /// The exchange of heat and vapour between a ground at `ground` K and its air at `air` K
    /// under a wind of `wind` m/s, over its neutral; the neutral's where the body's gravity
    /// or its air's heat capacity is not known.
    pub fn exchange(&self, ground: f64, air: f64, wind: f64) -> f64 {
        if self.gravity <= 0.0 || self.heat_capacity <= 0.0 {
            return 1.0;
        }
        over_neutral(self.richardson(ground, air, wind))
    }
}

/// The integrated profile functions psi_M and psi_H at zeta = z / L (Paulson 1970): the Dyer
/// and Hicks forms below zero, Beljaars and Holtslag's (1991) above.
pub fn profiles(zeta: f64) -> (f64, f64) {
    if zeta < 0.0 {
        let x = (1.0 - PROFILE_UNSTABLE * zeta).powf(0.25);
        let half = |v: f64| ((1.0 + v) / 2.0).ln();
        let m = 2.0 * half(x) + half(x * x) - 2.0 * x.atan() + std::f64::consts::FRAC_PI_2;
        (m, 2.0 * half(x * x))
    } else {
        let [a, b, c, d] = PROFILE_STABLE;
        let tail = b * (zeta - c / d) * (-d * zeta).exp() + b * c / d;
        let m = -(a * zeta + tail);
        let h = -((1.0 + 2.0 * a * zeta / 3.0).powf(1.5) + tail - 1.0);
        (m, h)
    }
}

/// The bulk Richardson number and the exchange over its neutral at zeta = z / L.
pub fn at_zeta(zeta: f64) -> (f64, f64) {
    let log = VON_KARMAN / AIR_EXCHANGE_COEFFICIENT.sqrt();
    let (pm, ph) = profiles(zeta);
    let (qm, qh) = profiles(zeta * (-log).exp());
    let (m, h) = (log - pm + qm, log - ph + qh);
    (zeta * h / (m * m), log * log / (m * h))
}

/// The table of `at_zeta` over zeta, from the most unstable to the most stable,
/// SURFACE_LAYER_STEPS to a decade of |zeta| from SURFACE_LAYER_ZETA's least to its most
/// either side of neutral.
fn table() -> &'static [(f64, f64)] {
    static TABLE: OnceLock<Vec<(f64, f64)>> = OnceLock::new();
    TABLE.get_or_init(|| {
        let (least, unstable, stable) = SURFACE_LAYER_ZETA;
        let side = |most: f64| {
            let decades = (most / least).log10();
            let n = (decades * SURFACE_LAYER_STEPS as f64).round() as usize;
            (0..=n).map(move |k| least * 10f64.powf(decades * k as f64 / n as f64))
        };
        let mut zetas: Vec<f64> = side(unstable).map(|z| -z).collect();
        zetas.reverse();
        zetas.push(0.0);
        zetas.extend(side(stable));
        zetas.into_iter().map(at_zeta).collect()
    })
}

/// The exchange of heat and vapour over its neutral at the bulk Richardson number `ri`, read
/// between the table's entries; past its ends its last.
pub fn over_neutral(ri: f64) -> f64 {
    let t = table();
    let i = t.partition_point(|&(r, _)| r < ri).clamp(1, t.len() - 1);
    let ((r0, f0), (r1, f1)) = (t[i - 1], t[i]);
    let s = ((ri - r0) / (r1 - r0)).clamp(0.0, 1.0);
    f0 + s * (f1 - f0)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// M63: along the table the Richardson number rises with zeta, so that it is read the
    /// other way, and the exchange falls.
    #[test]
    fn the_exchange_falls_as_the_air_stands_the_more_stable() {
        let t = table();
        assert!(t.windows(2).all(|w| w[1].0 > w[0].0));
        assert!(t.windows(2).all(|w| w[1].1 < w[0].1));
        let (lo, hi) = (t[0].0, t[t.len() - 1].0);
        println!("the table from Ri {lo:.1} to {hi:.1}");
        assert!(lo < -1000.0 && hi > 1000.0, "{lo} {hi}");
    }

    /// M63: near neutral the exchange falls with the Richardson number as the measured
    /// profiles' slopes give it: over a stable layer by phi_M's and phi_H's 5 zeta each, some
    /// ten times the number; over an unstable one it grows by phi_M's 4 and phi_H's 8, some
    /// twelve times; neutral at none.
    #[test]
    fn near_neutral_the_exchange_follows_the_measured_slopes() {
        let ri = 1e-3;
        let stable = (1.0 - over_neutral(ri)) / ri;
        let unstable = (over_neutral(-ri) - 1.0) / ri;
        println!("the exchange's slopes near neutral: {stable:.2} stable, {unstable:.2} unstable");
        assert!((over_neutral(0.0) - 1.0).abs() < 1e-9);
        assert!((stable - 10.0).abs() < 0.5, "{stable}");
        assert!((unstable - 12.0).abs() < 0.6, "{unstable}");
    }

    /// M63: no critical Richardson number: past the 0.2 at which the linear profiles stop the
    /// turbulence the exchange keeps on, falling by degrees: a fifth of its neutral at 0.2, a
    /// thirtieth at 1, under a thousandth only past 5.
    #[test]
    fn the_exchange_dies_away_by_degrees_under_an_inversion() {
        let at = over_neutral;
        println!(
            "the exchange over its neutral at Ri 0.2: {:.4}, 0.5: {:.4}, 1: {:.5}, 10: {:.2e}, 100: {:.2e}",
            at(0.2),
            at(0.5),
            at(1.0),
            at(10.0),
            at(100.0)
        );
        assert!(at(0.2) > 0.15 && at(0.5) > 0.05 && at(1.0) > 0.02);
        assert!(at(10.0) < 1e-3 && at(100.0) > 0.0);
        let mut last = 1.0;
        for k in 1..=400 {
            let now = at(0.05 * k as f64);
            assert!(now > 0.0 && now < last, "{k}");
            last = now;
        }
    }

    /// M63: a ground under its air: the Richardson number by the potential temperatures, the
    /// air of Earth 0.1 K the colder at neutral, and the exchange over a cold ground under a
    /// light wind the less than under a strong one, over a warm one the more than neutral.
    #[test]
    fn a_ground_under_its_air_reads_its_stability() {
        let earth = SurfaceLayer {
            gravity: 9.81,
            heat_capacity: 1004.0,
        };
        let neutral = 288.0 - 9.81 / 1004.0 * WIND_REFERENCE_HEIGHT;
        assert!(earth.richardson(288.0, neutral, 5.0).abs() < 1e-12);
        let (calm, windy) = (
            earth.exchange(285.0, 288.0, 1.0),
            earth.exchange(285.0, 288.0, 8.0),
        );
        println!("Earth's ground 3 K under its air: {calm:.4} at 1 m/s, {windy:.4} at 8 m/s");
        assert!(calm < windy && windy < 1.0);
        assert!(earth.exchange(293.0, 288.0, 2.0) > 1.0);
        assert_eq!(SurfaceLayer::default().exchange(250.0, 300.0, 1.0), 1.0);
    }
}
