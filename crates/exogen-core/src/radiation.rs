// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! The air's radiation (F293): what a column of air over the ground lets out to space at the
//! ground's temperature, and what its molecules scatter back of the star's light.
//!
//! The column stands on the adiabat of its air from the ground up, moist where it holds water,
//! its vapour's latent heat slowing its cooling as the vapour condenses (the pseudo-adiabat,
//! Holton 2004), dry where the air is too cold to hold any, never colder than the frost point
//! of its bulk gas, where that gas condenses and its latent heat holds the air on the frost
//! curve (Kasting 1991), nor than the skin over it, the grey skin's 2^(-1/4) of the body's
//! equilibrium temperature. Its water vapour lies where its temperatures hold it, each level
//! saturated alike (F235's relative humidity). Its radiation leaves by the two-stream equations over
//! wavenumber, each layer emitting by its temperature into the transmission above it at the
//! diffusivity 1.66 (Elsasser 1942):
//!
//! ```text
//! OLR = pi sum_nu [ B(nu, T_g) e^(-D tau_g) + sum_layers B(nu, T_l) (e^(-D tau_top) - e^(-D tau_bottom)) ] d nu
//! ```
//!
//! Carbon dioxide's 15 um band, the bending mode's, absorbs by its lines' wings falling off
//! exponentially from its centre, kappa_0 e^(-|nu - nu_0| / l) per kg/m² at a reference
//! pressure and broadened in proportion to the pressure (Williams 2025, from HITRAN2016), so
//! that its band widens with the logarithm of its column and saturates; carbon dioxide's
//! molecules colliding absorb in the far infrared as the square of their density (Ho et al.
//! 1971), which a thick cold atmosphere's window below 250 cm^-1 closes with. Water vapour,
//! methane and nitrogen's collisions stay grey, as the zero-dimensional climate had them
//! (F21, F235), and the sulfuric clouds an atmosphere's sulfur dioxide makes, grey and well
//! mixed: water set by Earth's mean temperature and the slope of its outgoing radiation,
//! methane and nitrogen by Titan's measured mean, the clouds by Venus's.

use std::collections::BTreeMap;
use std::sync::OnceLock;

use crate::constants::*;
use crate::params::derive::physics::{
    mean_heat_capacity, mean_molecular_mass, molar_mass, saturation_temperature, vapour_pressure,
    volatile, water_column,
};

/// Planck's radiance at a wavenumber, cm^-1, and a temperature, K: W/m²/sr per cm^-1.
pub fn planck(nu: f64, t: f64) -> f64 {
    if t <= 0.0 {
        return 0.0;
    }
    RADIATION_C1 * nu.powi(3) / (RADIATION_C2_CM_K * nu / t).exp_m1()
}

/// The wavenumbers the column's radiation is summed over, cm^-1, with their widths: finely
/// over carbon dioxide's band, where its absorption changes by a factor e every 12 cm^-1,
/// coarsely elsewhere; with carbon dioxide's band's absorption per kg/m² at the reference
/// pressure and its molecules' collisions' per amagat² cm at each.
struct Spectrum {
    nu: Vec<f64>,
    width: Vec<f64>,
    band: Vec<f64>,
    collisions: Vec<f64>,
}

fn spectrum() -> &'static Spectrum {
    static SPECTRUM: OnceLock<Spectrum> = OnceLock::new();
    SPECTRUM.get_or_init(|| {
        let (mut nu, mut width) = (Vec::new(), Vec::new());
        for &(from, to, step) in &RADIATION_BINS {
            let mut x = from;
            while x < to - 1e-9 {
                let w = step.min(to - x);
                nu.push(x + 0.5 * w);
                width.push(w);
                x += w;
            }
        }
        let (kappa, centre, wing) = CO2_BAND;
        let band = nu
            .iter()
            .map(|&n: &f64| kappa * (-(n - centre).abs() / wing).exp())
            .collect();
        let (integral, from, to, peak) = CO2_COLLISIONS;
        // A band peaking at `peak`, (nu / peak) e^(1 - nu / peak), its integral over the
        // measured span the measured one.
        let shape = |x: f64| (x + 1.0) * (-x).exp();
        let over = peak * std::f64::consts::E * (shape(from / peak) - shape(to / peak));
        let collisions = nu
            .iter()
            .map(|&n| integral / over * (n / peak) * (1.0 - n / peak).exp())
            .collect();
        Spectrum {
            nu,
            width,
            band,
            collisions,
        }
    })
}

/// A column of air over the ground, as its outgoing radiation reads it (F293).
#[derive(Clone, Debug)]
pub struct Column {
    /// The surface pressure, Pa, and the gravity, m/s².
    pressure: f64,
    gravity: f64,
    /// The air's gas constant, J/kg/K, its heat capacity, J/kg/K, and its dry adiabat's
    /// exponent, R / c_p; water's molecules' mass against the air's.
    gas_constant: f64,
    heat: f64,
    adiabat: f64,
    epsilon: f64,
    /// Carbon dioxide's share of the air's mass and of its molecules.
    co2_mass: f64,
    co2_volume: f64,
    /// The bulk gas, whose frost holds the air, and its share of the molecules.
    bulk: Option<(&'static Volatile, f64)>,
    /// The skin's temperature over the air, K.
    skin: f64,
    /// The grey absorbers' depths at the surface, and how they thin with the pressure:
    /// methane's and nitrogen's collisions', and the sulfuric clouds'.
    grey: Vec<(f64, f64)>,
    /// The water a warm air can hold vapour from, metres of global layer, and its vapour's
    /// grey depth's law, a ln(1 + u / u0) of its column u, kg/m² (F235).
    water: f64,
    water_tau: (f64, f64),
}

impl Column {
    /// The column of an air of `composition` by volume at `pressure`, Pa, under `gravity`,
    /// over a body with `water` metres of global layer of it, its skin at `skin` K.
    pub fn new(
        composition: &BTreeMap<String, f64>,
        pressure: f64,
        gravity: f64,
        water: f64,
        skin: f64,
    ) -> Self {
        let total: f64 = composition.values().sum::<f64>().max(1e-300);
        let mu = mean_molecular_mass(composition);
        let gas_constant = GAS_CONSTANT / mu;
        let heat = mean_heat_capacity(composition);
        let adiabat = gas_constant / heat;
        let share = |s: &str| composition.get(s).copied().unwrap_or(0.0) / total;
        let mass = |s: &str| share(s) * molar_mass(s) / mu;
        let column = pressure / gravity.max(1e-300);
        let bulk = bulk_of(composition);
        let u = |s: &str| mass(s) * column;
        let grey = vec![
            (TAU_CH4.0 * u("CH4").powf(TAU_CH4.1), TAU_CH4.1),
            (TAU_N2_CIA.0 * u("N2").powf(TAU_N2_CIA.1), TAU_N2_CIA.1),
            (SULFURIC_CLOUD_DEPTH * u("SO2"), 1.0),
        ];
        Self {
            pressure,
            gravity,
            gas_constant,
            heat,
            adiabat,
            epsilon: molar_mass("H2O") / mu,
            co2_mass: mass("CO2"),
            co2_volume: share("CO2"),
            bulk,
            skin,
            grey,
            water,
            water_tau: TAU_H2O,
        }
    }

    /// The adiabat's exponent, d ln T / d ln p, at `t` K and `p` Pa: the dry air's, slowed where
    /// the column holds water by its saturated vapour's latent heat,
    /// (1 + L r / (R T)) / (1 + L² r eps / (c_p R T²)), r its mixing ratio.
    fn lapse(&self, t: f64, p: f64) -> f64 {
        let Some(water) = volatile("H2O").filter(|_| self.water > 0.0) else {
            return self.adiabat;
        };
        let e = vapour_pressure(water, t).min(0.5 * p);
        let r = self.epsilon * e / (p - e);
        let l = if t >= water.t_triple {
            water.l_vap
        } else {
            water.l_sub
        };
        let rt = self.gas_constant * t;
        self.adiabat * (1.0 + l * r / rt) / (1.0 + l * l * r * self.epsilon / (self.heat * rt * t))
    }

    /// The air's temperatures at `levels` over ground at `t` K, the first level the ground's:
    /// up the adiabat by midpoint steps in the pressure's logarithm, never colder than the bulk
    /// gas's frost nor the skin; and the water vapour above each level, by the saturated
    /// vapour each holds, in the units of its ground's (none without water).
    fn profile(&self, t: f64, levels: &[f64]) -> (Vec<f64>, Vec<f64>) {
        let n = levels.len();
        let mut adiabat = vec![t; n];
        for k in 0..n - 1 {
            let step = (levels[k + 1] / levels[k]).ln();
            let mid = adiabat[k] * (0.5 * step * self.lapse(adiabat[k], levels[k])).exp();
            let p_mid = (levels[k] * levels[k + 1]).sqrt();
            adiabat[k + 1] = adiabat[k] * (step * self.lapse(mid, p_mid)).exp();
        }
        let temps: Vec<f64> = adiabat
            .iter()
            .zip(levels)
            .map(|(&a, &p)| {
                let frost = self
                    .bulk
                    .map_or(0.0, |(v, x)| saturation_temperature(v, x * p));
                a.max(frost).max(self.skin)
            })
            .collect();
        let mut vapour = vec![0.0; n];
        if let Some(water) = volatile("H2O").filter(|_| self.water > 0.0) {
            let q = |k: usize| vapour_pressure(water, temps[k]).min(0.5 * levels[k]) / levels[k];
            for k in (0..n - 1).rev() {
                vapour[k] = vapour[k + 1] + 0.5 * (q(k) + q(k + 1)) * (levels[k] - levels[k + 1]);
            }
        }
        (temps, vapour)
    }

    /// The outgoing radiation, W/m², of the column over ground at `t` K.
    pub fn olr(&self, t: f64) -> f64 {
        self.olr_at(t, 1.0)
    }

    /// The outgoing radiation, W/m², of the column over ground at `t` K standing at `ratio` of
    /// the surface pressure: the air above it, and of the water vapour the air over the datum
    /// at `t` holds (F235), the share its profile holds above it.
    pub fn olr_at(&self, t: f64, ratio: f64) -> f64 {
        let ground = self.pressure * ratio;
        if ground <= 0.0 || t <= 0.0 {
            return SIGMA_SB * t.max(0.0).powi(4);
        }
        let n = RADIATION_LEVELS;
        // Crowded toward the ground as the square of their rank, where the water vapour's
        // opacity lies within a fifth of the air's scale height.
        let levels: Vec<f64> = (0..n)
            .map(|k| {
                let u = k as f64 / (n - 1) as f64;
                ground * (-RADIATION_DECADES * std::f64::consts::LN_10 * u * u).exp()
            })
            .collect();
        let (temps, vapour) = self.profile(t, &levels);
        // Carbon dioxide's band's path above each level, kg/m² at the reference pressure, and
        // its molecules' collisions', amagat² cm.
        let band: Vec<f64> = levels
            .iter()
            .map(|&p| self.co2_mass * p * p / (2.0 * self.gravity * CO2_BAND_PRESSURE))
            .collect();
        let mut collisions = vec![0.0; n];
        let density = self.co2_volume * AMAGAT_TEMPERATURE / AMAGAT_PRESSURE;
        for k in (0..n - 1).rev() {
            let (p0, p1) = (levels[k + 1], levels[k]);
            let mean = 0.5 * (p0 / temps[k + 1] + p1 / temps[k]);
            collisions[k] = collisions[k + 1]
                + density * density * 100.0 * self.gas_constant / self.gravity * mean * (p1 - p0);
        }
        // The water vapour above each level, kg/m²: the datum's column (F235) by the share of
        // its vapour above the level; each depth by the saturating law of the water above it.
        let datum = if ratio == 1.0 || vapour[0] <= 0.0 {
            vapour[0]
        } else {
            let datum_levels: Vec<f64> = levels.iter().map(|&p| p / ratio).collect();
            self.profile(t, &datum_levels).1[0]
        };
        let column = if datum > 0.0 {
            water_column(t, self.water) / datum
        } else {
            0.0
        };
        let (a, u0) = self.water_tau;
        let grey: Vec<f64> = levels
            .iter()
            .zip(&vapour)
            .map(|(&p, &v)| {
                let x = p / self.pressure;
                let dry: f64 = self.grey.iter().map(|&(tau, b)| tau * x.powf(b)).sum();
                dry + a * (column * v / u0).ln_1p()
            })
            .collect();
        let spectrum = spectrum();
        let mut total = 0.0;
        for (i, &nu) in spectrum.nu.iter().enumerate() {
            let depth = |k: usize| {
                band[k] * spectrum.band[i] + collisions[k] * spectrum.collisions[i] + grey[k]
            };
            let mut above = (-DIFFUSIVITY * depth(n - 1)).exp();
            let mut out = 0.0;
            for k in (0..n - 1).rev() {
                let below = (-DIFFUSIVITY * depth(k)).exp();
                out += planck(nu, 0.5 * (temps[k] + temps[k + 1])) * (above - below);
                above = below;
            }
            out += planck(nu, t) * above;
            total += out * spectrum.width[i];
        }
        std::f64::consts::PI * total
    }

    /// The ground's temperature at which the column lets out `absorbed` W/m², K: never below
    /// the frost point of its bulk gas at the surface, where that gas's frost holds it.
    pub fn surface_temperature(&self, absorbed: f64) -> f64 {
        let floor = self
            .bulk
            .map_or(0.0, |(v, x)| saturation_temperature(v, x * self.pressure));
        balance(|t| self.olr(t), absorbed).max(floor)
    }
}

/// The air's bulk gas, its largest but its water vapour, where it has a frost, and its share of
/// the molecules.
fn bulk_of(composition: &BTreeMap<String, f64>) -> Option<(&'static Volatile, f64)> {
    let total: f64 = composition.values().sum::<f64>().max(1e-300);
    composition
        .iter()
        .filter(|(s, _)| s.as_str() != "H2O")
        .max_by(|a, b| a.1.total_cmp(b.1))
        .and_then(|(s, f)| volatile(s).map(|v| (v, f / total)))
}

/// The temperature at which `olr`, rising with it, reaches `absorbed`, K: bracketed from the
/// black body's, then by the secant held within the bracket (Illinois), in the logarithm.
fn balance(olr: impl Fn(f64) -> f64, absorbed: f64) -> f64 {
    let black = (absorbed.max(1e-30) / SIGMA_SB).powf(0.25);
    let f = |x: f64| olr(x.exp()) - absorbed;
    let (mut a, mut b) = (black.ln(), black.ln());
    let (mut fa, mut fb) = (f(a), f(a));
    for _ in 0..RADIATION_SOLVE_STEPS {
        if fa <= 0.0 && fb >= 0.0 && b > a {
            break;
        }
        if fa > 0.0 {
            (b, fb) = (a, fa);
            a -= 0.4;
            fa = f(a);
        } else {
            (a, fa) = (b, fb);
            b += 0.4;
            fb = f(b);
        }
    }
    if fa > 0.0 || fb < 0.0 {
        return if fa > 0.0 { a.exp() } else { b.exp() };
    }
    let mut side = 0;
    for _ in 0..RADIATION_SOLVE_STEPS {
        let c = (a * fb - b * fa) / (fb - fa);
        let fc = f(c);
        if fc.abs() < 1e-6 * absorbed || (b - a) < 1e-7 {
            return c.exp();
        }
        if fc < 0.0 {
            (a, fa) = (c, fc);
            if side == -1 {
                fb *= 0.5;
            }
            side = -1;
        } else {
            (b, fb) = (c, fc);
            if side == 1 {
                fa *= 0.5;
            }
            side = 1;
        }
    }
    ((a * fb - b * fa) / (fb - fa)).exp()
}

/// The share of a black body's radiation a column lets out at each ground temperature (F293),
/// tabulated over the logarithm of the temperature and read between the entries.
#[derive(Clone, Debug)]
pub struct Outgoing {
    from: f64,
    step: f64,
    share: Vec<f64>,
}

impl Outgoing {
    /// None held back: an airless body's.
    pub fn none() -> Self {
        Self {
            from: 0.0,
            step: 1.0,
            share: vec![1.0],
        }
    }

    /// The column's share from `low` to `high` K.
    pub fn of(column: &Column, low: f64, high: f64) -> Self {
        use rayon::prelude::*;
        let (from, to) = (low.max(1.0).ln(), high.max(low + 1.0).ln());
        let n = ((to - from) / OUTGOING_STEP).ceil() as usize + 1;
        let step = (to - from) / (n - 1) as f64;
        let share = (0..n)
            .into_par_iter()
            .map(|k| {
                let t = (from + step * k as f64).exp();
                column.olr(t) / (SIGMA_SB * t.powi(4))
            })
            .collect();
        Self { from, step, share }
    }

    /// The table from its entries, (temperature, share).
    pub fn from_points(points: &[(f64, f64)]) -> Self {
        match points {
            [] => Self::none(),
            [only] => Self {
                from: only.0.max(1.0).ln(),
                step: 1.0,
                share: vec![only.1],
            },
            _ => {
                let from = points[0].0.max(1.0).ln();
                let step =
                    (points[points.len() - 1].0.max(1.0).ln() - from) / (points.len() - 1) as f64;
                Self {
                    from,
                    step,
                    share: points.iter().map(|p| p.1).collect(),
                }
            }
        }
    }

    /// The entries, (temperature, share).
    pub fn points(&self) -> Vec<(f64, f64)> {
        self.share
            .iter()
            .enumerate()
            .map(|(k, &s)| ((self.from + self.step * k as f64).exp(), s))
            .collect()
    }

    /// The share at `t` K, held at the table's ends past them.
    pub fn share(&self, t: f64) -> f64 {
        if self.share.len() == 1 {
            return self.share[0];
        }
        let x =
            ((t.max(1e-9).ln() - self.from) / self.step).clamp(0.0, (self.share.len() - 1) as f64);
        let k = (x.floor() as usize).min(self.share.len() - 2);
        let f = x - k as f64;
        self.share[k] + f * (self.share[k + 1] - self.share[k])
    }

    /// The outgoing radiation at `t` K, W/m².
    pub fn olr(&self, t: f64) -> f64 {
        self.share(t) * SIGMA_SB * t.powi(4)
    }

    /// Whether the column holds anything back.
    pub fn holds(&self) -> bool {
        self.share.iter().any(|&s| s < 1.0 - 1e-9)
    }
}

/// The share a column lets out at each ground temperature and each height of the ground
/// (F293), as the downscaling reads it: over the logarithms of the temperature and of the
/// ground's share of the surface pressure, read between the entries.
#[derive(Clone, Debug)]
pub struct Heights {
    from: (f64, f64),
    step: (f64, f64),
    count: (usize, usize),
    share: Vec<f64>,
}

impl Heights {
    /// The column's share from `low` to `high` K, its ground from HEIGHTS_RATIO's span of
    /// the surface pressure's logarithm.
    pub fn of(column: &Column, low: f64, high: f64) -> Self {
        use rayon::prelude::*;
        let (t0, t1) = (low.max(1.0).ln(), high.max(low + 1.0).ln());
        let nt = ((t1 - t0) / HEIGHTS_STEP.0).ceil() as usize + 1;
        let (r0, r1) = HEIGHTS_RATIO;
        let nr = ((r1 - r0) / HEIGHTS_STEP.1).ceil() as usize + 1;
        let step = ((t1 - t0) / (nt - 1) as f64, (r1 - r0) / (nr - 1) as f64);
        let share = (0..nt * nr)
            .into_par_iter()
            .map(|i| {
                let (k, j) = (i / nr, i % nr);
                let t = (t0 + step.0 * k as f64).exp();
                let ratio = (r0 + step.1 * j as f64).exp();
                column.olr_at(t, ratio) / (SIGMA_SB * t.powi(4))
            })
            .collect();
        Self {
            from: (t0, r0),
            step,
            count: (nt, nr),
            share,
        }
    }

    /// The share at `t` K over ground at `ratio` of the surface pressure, held at the
    /// table's edges past them.
    pub fn share(&self, t: f64, ratio: f64) -> f64 {
        let (nt, nr) = self.count;
        let at = |x: f64, from: f64, step: f64, n: usize| {
            let x = ((x - from) / step).clamp(0.0, (n - 1) as f64);
            let k = (x.floor() as usize).min(n.saturating_sub(2));
            (k, x - k as f64)
        };
        let (k, fk) = at(t.max(1e-9).ln(), self.from.0, self.step.0, nt);
        let (j, fj) = at(ratio.max(1e-300).ln(), self.from.1, self.step.1, nr);
        let v = |k: usize, j: usize| self.share[k.min(nt - 1) * nr + j.min(nr - 1)];
        let low = v(k, j) + fj * (v(k, j + 1) - v(k, j));
        let high = v(k + 1, j) + fj * (v(k + 1, j + 1) - v(k + 1, j));
        low + fk * (high - low)
    }
}

/// The share of the star's light, a black body at `teff` K, a layer of molecules scattering
/// with an optical depth `tau` at 550 nm reflects (F293): Rayleigh's, as the wavelength to the
/// minus fourth, each wavelength by the conservative two-stream layer's 3 tau / 4 /
/// (1 + 3 tau / 4), weighted by the star's light.
pub fn rayleigh_reflectance(tau: f64, teff: f64) -> f64 {
    if tau <= 0.0 || teff <= 0.0 {
        return 0.0;
    }
    let n = RAYLEIGH_WAVELENGTHS;
    let (mut light, mut reflected) = (0.0, 0.0);
    for k in 0..n {
        let um =
            RAYLEIGH_SPAN.0 + (RAYLEIGH_SPAN.1 - RAYLEIGH_SPAN.0) * (k as f64 + 0.5) / n as f64;
        let b = um.powi(-5) / (RADIATION_C2_UM_K / (um * teff)).exp_m1();
        let t = tau * (0.55 / um).powi(4);
        light += b;
        reflected += b * 0.75 * t / (1.0 + 0.75 * t);
    }
    reflected / light.max(1e-300)
}

/// The optical depth at 550 nm of the molecules of an air of `composition` by volume, `column`
/// kg/m² of it: each gas's cross-section from its refractive index and its King factor,
/// 24 pi³ / (N_L² lambda⁴) ((n² - 1) / (n² + 2))² F_K (Bates 1984; Sneep and Ubachs 2005).
pub fn rayleigh_depth(composition: &BTreeMap<String, f64>, column: f64) -> f64 {
    let total: f64 = composition.values().sum::<f64>();
    if total <= 0.0 || column <= 0.0 {
        return 0.0;
    }
    let mu = mean_molecular_mass(composition);
    let molecules = column / (mu / AVOGADRO);
    let wavenumber: f64 = 1.0 / 550e-9;
    composition
        .iter()
        .map(|(s, f)| {
            let (refraction, king) = RAYLEIGH_GASES
                .iter()
                .find(|g| g.0 == s.as_str())
                .map_or((0.0, 1.0), |g| (g.1, g.2));
            let lorentz = 2.0 * refraction / 3.0;
            let sigma = 24.0 * std::f64::consts::PI.powi(3) * wavenumber.powi(4)
                / (LOSCHMIDT * LOSCHMIDT)
                * lorentz
                * lorentz
                * king;
            f / total * sigma
        })
        .sum::<f64>()
        * molecules
}

/// The albedo of a layer reflecting `r` of the light, conservative, over ground of albedo
/// `beneath`: what it reflects and what the ground returns through it, again and again.
pub fn over(r: f64, beneath: f64) -> f64 {
    r + (1.0 - r) * (1.0 - r) * beneath / (1.0 - r * beneath).max(1e-300)
}

/// The haze a methane atmosphere makes under its star's ultraviolet, the haze's index (F293):
/// methane's share of the molecules, the star's XUV against the Sun's at Earth, and the
/// pressure in bar up to two.
pub fn tholin_production(methane: f64, xuv: f64, pressure: f64) -> f64 {
    methane * xuv * 50.0 * (pressure / 1e5).min(2.0)
}

/// The Bond albedo of a body under an atmosphere (F293): its air's Rayleigh scattering of its
/// star's light over the albedo beneath, the sulfuric clouds its sulfur dioxide makes over
/// that, a conservative layer reflecting s / (1 + s) of the light, s growing with the sulfur
/// dioxide's column; its haze absorbing on each pass through it.
pub fn albedo(composition: &BTreeMap<String, f64>, column: f64, teff: f64, haze: f64) -> f64 {
    let scattered = over(
        rayleigh_reflectance(rayleigh_depth(composition, column), teff),
        ALBEDO_BENEATH,
    );
    let total: f64 = composition.values().sum::<f64>().max(1e-300);
    let sulfur = composition.get("SO2").copied().unwrap_or(0.0) / total * molar_mass("SO2")
        / mean_molecular_mass(composition)
        * column;
    let s = SULFURIC_CLOUD_REFLECTION * sulfur;
    over(s / (1.0 + s), scattered) * (-2.0 * HAZE_ABSORPTION * haze).exp()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::ParamGraph;

    fn entry(g: &ParamGraph, id: &str, key: &str) -> f64 {
        g.get(id)
            .as_map()
            .and_then(|m| m.get(key).copied())
            .unwrap_or(0.0)
    }

    fn graph(name: &str) -> ParamGraph {
        crate::presets::preset(name).unwrap().to_graph().0
    }

    /// A preset's column at its own pressure under the sunlight it absorbs at `albedo`, and that
    /// sunlight, W/m².
    fn column_of(g: &ParamGraph, albedo: f64) -> (Column, f64) {
        let composition = g.get("atmosphere.composition").as_map().unwrap().clone();
        let water = entry(g, "volatiles.inventory", "H2O")
            + if composition.get("H2O").copied().unwrap_or(0.0) > 0.0 {
                1.0
            } else {
                0.0
            };
        let absorbed = g.f64("orbit.insolation") / 4.0 * (1.0 - albedo);
        let skin = (absorbed / SIGMA_SB).powf(0.25) * 0.5f64.powf(0.25);
        let column = Column::new(
            &composition,
            g.f64("atmosphere.surface_pressure"),
            g.f64("body.gravity"),
            water,
            skin,
        );
        (column, absorbed)
    }

    fn halve(f: impl Fn(f64) -> f64, mut lo: f64, mut hi: f64) -> f64 {
        // f rising through zero between lo and hi, halved in the logarithm
        for _ in 0..40 {
            let mid = (lo * hi).sqrt();
            if f(mid) < 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        (lo * hi).sqrt()
    }

    fn rayleigh_of(g: &ParamGraph) -> f64 {
        let composition = g.get("atmosphere.composition").as_map().unwrap().clone();
        over(
            rayleigh_reflectance(
                rayleigh_depth(&composition, g.f64("atmosphere.column_mass")),
                g.f64("star.teff"),
            ),
            ALBEDO_BENEATH,
        )
    }

    /// F293: carbon dioxide's band widens with the logarithm of its column: over a 288 K
    /// ground under an Earth-like air, each doubling of its carbon dioxide from Earth's takes
    /// within a fifth of the same from the outgoing radiation, some 5 W/m² as a dry air's
    /// line-by-line forcing is (Jeevanjee et al. 2021); the grey law it replaced grew as the
    /// column to the 0.86.
    #[test]
    fn carbon_dioxide_s_band_saturates() {
        let olr = |x: f64| {
            let composition: BTreeMap<String, f64> =
                [("N2".to_string(), 1.0 - x), ("CO2".to_string(), x)]
                    .into_iter()
                    .collect();
            Column::new(&composition, 101_325.0, 9.81, 0.0, 214.0).olr(288.0)
        };
        let doublings: Vec<f64> = (0..4)
            .map(|k| {
                let x = 4.2e-4 * 2f64.powi(k);
                olr(x) - olr(2.0 * x)
            })
            .collect();
        println!("each doubling's forcing from 420 ppm: {doublings:.2?} W/m²");
        for d in &doublings {
            assert!(
                (*d - doublings[0]).abs() < 0.2 * doublings[0],
                "{doublings:?}"
            );
        }
        assert!(doublings[0] > 4.0 && doublings[0] < 6.5, "{}", doublings[0]);
    }

    /// F293: Earth's air's Rayleigh depth at 550 nm, measured 0.0973 (Bucholtz 1995), and its
    /// albedo, 0.30 measured, under the Sun; under a red dwarf at 2,566 K the same air
    /// scatters a sixteenth as much of its light.
    #[test]
    fn earth_s_air_scatters_as_measured() {
        let earth = graph("Earth");
        let composition = earth
            .get("atmosphere.composition")
            .as_map()
            .unwrap()
            .clone();
        let tau = rayleigh_depth(&composition, earth.f64("atmosphere.column_mass"));
        let (sun, red) = (
            rayleigh_reflectance(tau, 5_772.0),
            rayleigh_reflectance(tau, 2_566.0),
        );
        let albedo = earth.f64("atmosphere.bond_albedo");
        println!(
            "Earth: tau {tau:.4}, reflects {sun:.4} of the Sun's light, {red:.4} of a red dwarf's; albedo {albedo:.3}"
        );
        assert!((tau - 0.0973).abs() < 0.005, "{tau}");
        assert!((albedo - 0.30).abs() < 0.01, "{albedo}");
        assert!(red < sun / 8.0, "{red} {sun}");
    }

    /// F293: Mars through its thicker air (F263, F264) stays frozen, its greenhouse saturating
    /// as its scattering grows and its carbon dioxide condensing aloft: under 230 K, the most
    /// any carbon dioxide and water atmosphere gives it (Kasting 1991; Forget et al. 2013;
    /// Wordsworth et al. 2013), at 3.5 and 3.8 Ga, where the grey law gave some 289 K.
    #[test]
    fn early_mars_stays_frozen_under_its_thicker_air() {
        let mars = graph("Mars");
        let age = mars.f64("body.age");
        let (pressures, temperatures) = (
            mars.curve("atmosphere.pressure_history").unwrap(),
            mars.curve("atmosphere.surface_temperature_history")
                .unwrap(),
        );
        for ago in [3.8, 3.5, 3.0, 2.0, 0.0] {
            let t = age - ago;
            println!(
                "Mars {ago} Ga ago: {:.3} bar, {:.1} K",
                pressures.eval(t) / 1e5,
                temperatures.eval(t)
            );
        }
        for ago in [3.8, 3.5] {
            let t = temperatures.eval(age - ago);
            assert!(t < 230.0 && t > 180.0, "{ago} Ga: {t}");
        }
        let now = temperatures.eval(age);
        let mean = mars.f64("atmosphere.surface_temperature_mean");
        assert!((now - mean).abs() < 1.0, "{now} {mean}");
    }

    /// F293: the anchors, found: Titan's methane and nitrogen at its measured 93.7 K under its
    /// measured Bond albedo 0.265 (Fulchignoni et al. 2005); Earth's water at its
    /// zero-dimensional 291.6 K with its outgoing radiation rising 2.1 W/m² per kelvin (F235);
    /// Venus's sulfuric clouds at its 737 K; Titan's haze at its 0.265.
    #[test]
    #[ignore]
    fn calibrate_the_greenhouse_s_anchors() {
        let titan = graph("Titan");
        let (base, absorbed) = column_of(&titan, 0.265);
        let scale = halve(
            |f| {
                let mut c = base.clone();
                c.grey[0].0 *= f;
                c.grey[1].0 *= f;
                c.surface_temperature(absorbed) - 93.7
            },
            1e-3,
            1e3,
        );
        println!(
            "Titan: methane's and nitrogen's depths x {scale:.5}: TAU_CH4.0 {:.5e}, TAU_N2_CIA.0 {:.5e}",
            TAU_CH4.0 * scale,
            TAU_N2_CIA.0 * scale
        );
        let earth = graph("Earth");
        let (mut base, absorbed) = column_of(&earth, rayleigh_of(&earth));
        base.grey[1].0 *= scale;
        let at = |u0: f64| {
            let a = halve(
                |a| {
                    let mut c = base.clone();
                    c.water_tau = (a, u0);
                    c.surface_temperature(absorbed) - 291.6
                },
                1e-3,
                1e2,
            );
            let mut c = base.clone();
            c.water_tau = (a, u0);
            (a, (c.olr(292.6) - c.olr(290.6)) / 2.0)
        };
        let u0 = halve(|u| 2.1 - at(u).1, 1e-3, 1e2);
        let (a, slope) = at(u0);
        println!(
            "Earth: albedo {:.4}, water a {a:.5}, u0 {u0:.5} kg/m², slope {slope:.3}",
            rayleigh_of(&earth)
        );
        let venus = graph("Venus");
        let scattered = rayleigh_of(&venus);
        let clouds = halve(|r| over(r, scattered) - 0.76, 1e-6, 0.999);
        let so2_column = {
            let composition = venus.get("atmosphere.composition").as_map().unwrap();
            let total: f64 = composition.values().sum();
            composition["SO2"] / total * molar_mass("SO2") / mean_molecular_mass(composition)
                * venus.f64("atmosphere.column_mass")
        };
        println!(
            "Venus: scattered {scattered:.4}, its clouds reflect {clouds:.4}: SULFURIC_CLOUD_REFLECTION {:.5e}",
            clouds / (1.0 - clouds) / so2_column
        );
        let a_v = 0.76;
        let (mut base, absorbed) = column_of(&venus, a_v);
        base.grey[1].0 *= scale;
        let so2 = base.grey[2].0 / SULFURIC_CLOUD_DEPTH;
        let k = halve(
            |k| {
                let mut c = base.clone();
                c.grey[2].0 = k * so2;
                c.surface_temperature(absorbed) - 737.0
            },
            1e-6,
            1e4,
        );
        println!("Venus: albedo {a_v:.4}, SO2 {so2:.2} kg/m², clouds {k:.5e} per kg/m²");
        let haze = tholin_production(
            entry(&titan, "atmosphere.composition", "CH4"),
            titan.f64("star.xuv_relative"),
            titan.f64("atmosphere.surface_pressure"),
        );
        let scattered = rayleigh_of(&titan);
        let absorption = -(0.265 / scattered).ln() / (2.0 * haze);
        println!(
            "Titan: scattered {scattered:.4}, haze {haze:.5}, HAZE_ABSORPTION {absorption:.5}"
        );
    }
}
