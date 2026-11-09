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
//! pressure and broadened in proportion to the pressure, and water vapour's rotation and
//! vibration-rotation bands alike, with its continuum by its own pressure and its temperature
//! (Williams 2026, from HITRAN2016 and MT_CKD), so that their bands widen with the logarithm of
//! their columns and saturate; carbon dioxide's molecules colliding absorb in the far infrared as
//! the square of their density (Ho et al. 1971), which a thick cold atmosphere's window below 250
//! cm^-1 closes with, and nitrogen's, methane's and hydrogen's by HITRAN's tables of their pairs
//! (M55). The sulfuric clouds an atmosphere's sulfur dioxide makes stay grey and well mixed, set
//! by Venus's measured mean; the clouds and the other absorbers of a water-bearing air, grey and
//! well mixed, a power of its vapour's column, set by Earth's, until the weather's clouds (5G).
//! Each layer's source is linear in its depth (Clough et al. 1992), and the column gives what it
//! lets out to space, the share of the ground's own radiation that reaches space through it, and
//! what its air sends down to the ground (M55).

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
    /// Water vapour's lines per kg/m² at the reference pressure, and its continuum per kg/m²
    /// at its reference state (M55).
    water: Vec<f64>,
    continuum: Vec<f64>,
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
        let water = nu.iter().map(|&n| water_lines(n)).collect();
        let continuum = nu
            .iter()
            .map(|&n| {
                H2O_CONTINUUM
                    .iter()
                    .find(|&&(to, _)| n < to)
                    .map_or(0.0, |&(_, k)| k)
            })
            .collect();
        Spectrum {
            nu,
            width,
            band,
            collisions,
            water,
            continuum,
        }
    })
}

/// The gases' collisions (M55): each pair's absorption, cm^-1 amagat^-2, on the spectrum's
/// wavenumbers at each of its tables' temperatures, from HITRAN's tables coarse-grained by 10
/// cm^-1 (data/collisions.txt): nitrogen's (Chistikov et al. 2019), nitrogen and methane's
/// (Finenko et al. 2022, fitted to Cassini's spectra of Titan), methane's (Borysow and
/// Frommhold 1987), nitrogen and hydrogen's (Borysow and Frommhold 1986).
struct Collisions {
    pairs: Vec<Pair>,
}

struct Pair {
    gases: (String, String),
    temperatures: Vec<f64>,
    /// The absorption at each temperature on each of the spectrum's wavenumbers.
    absorption: Vec<Vec<f64>>,
}

impl Pair {
    /// The table's temperatures about `t`, K, the lower's index and the share of the upper's,
    /// held at the ends.
    fn bracket(&self, t: f64) -> (usize, f64) {
        let ts = &self.temperatures;
        if ts.len() < 2 || t <= ts[0] {
            return (0, 0.0);
        }
        let j = ts.partition_point(|&x| x <= t).clamp(1, ts.len() - 1) - 1;
        (j, ((t - ts[j]) / (ts[j + 1] - ts[j])).clamp(0.0, 1.0))
    }
}

fn collision_tables() -> &'static Collisions {
    static COLLISIONS: OnceLock<Collisions> = OnceLock::new();
    COLLISIONS.get_or_init(|| {
        let spectrum = spectrum();
        let mut pairs: Vec<Pair> = Vec::new();
        let mut rows: Vec<(f64, Vec<f64>)> = Vec::new();
        let mut gases = (String::new(), String::new());
        let mut temperatures: Vec<f64> = Vec::new();
        let finish = |gases: &(String, String), temperatures: &[f64], rows: &[(f64, Vec<f64>)]| {
            if rows.is_empty() {
                return None;
            }
            // Each wavenumber's absorption between the table's, none past its end.
            let at = |j: usize, nu: f64| {
                let k = rows.partition_point(|r| r.0 <= nu);
                if k == 0 || k >= rows.len() {
                    return if k == 0 { rows[0].1[j] } else { 0.0 };
                }
                let (a, b) = (&rows[k - 1], &rows[k]);
                a.1[j] + (nu - a.0) / (b.0 - a.0) * (b.1[j] - a.1[j])
            };
            Some(Pair {
                gases: gases.clone(),
                temperatures: temperatures.to_vec(),
                absorption: (0..temperatures.len())
                    .map(|j| spectrum.nu.iter().map(|&nu| at(j, nu)).collect())
                    .collect(),
            })
        };
        for line in include_str!("../data/collisions.txt").lines() {
            let mut words = line.split_whitespace();
            match words.next() {
                None | Some("#") => {}
                Some(w) if w.starts_with('#') => {}
                Some("pair") => {
                    pairs.extend(finish(&gases, &temperatures, &rows));
                    rows.clear();
                    let a = words.next().unwrap_or_default().to_string();
                    let b = words.next().unwrap_or_default().to_string();
                    gases = (a, b);
                }
                Some("T") => temperatures = words.filter_map(|w| w.parse().ok()).collect(),
                Some(nu) => {
                    if let Ok(nu) = nu.parse::<f64>() {
                        rows.push((nu, words.filter_map(|w| w.parse().ok()).collect()));
                    }
                }
            }
        }
        pairs.extend(finish(&gases, &temperatures, &rows));
        Collisions { pairs }
    })
}

/// Water vapour's lines at a wavenumber, cm^-1, per kg/m² at the reference pressure (M55): its
/// rotation band's below the split, its vibration-rotation band's above it to the end.
fn water_lines(nu: f64) -> f64 {
    let (rotation, flat, l) = H2O_ROTATION;
    let (vibration, from, to, below, above) = H2O_VIBRATION;
    if nu < flat {
        rotation
    } else if nu < H2O_LINES_SPLIT {
        rotation * (-(nu - flat) / l).exp()
    } else if nu < from {
        vibration * (-(from - nu) / below).exp()
    } else if nu < to {
        vibration
    } else if nu < H2O_LINES_END {
        vibration * (-(nu - to) / above).exp()
    } else {
        0.0
    }
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
    /// The pairs of its gases whose collisions absorb, each its index among the tables' and
    /// the product of the two gases' shares of the molecules (M55).
    pairs: Vec<(usize, f64)>,
    /// The water a warm air can hold vapour from, metres of global layer (F235).
    water: f64,
    /// The clouds' grey depth's law, a u^b of the column's vapour u, kg/m² (M55).
    cloud: (f64, f64),
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
        let grey = vec![(SULFURIC_CLOUD_DEPTH * u("SO2"), 1.0)];
        // The pairs of its gases whose collisions absorb, by the product of their shares of
        // the molecules (M55).
        let pairs = collision_tables()
            .pairs
            .iter()
            .enumerate()
            .map(|(k, p)| (k, share(&p.gases.0) * share(&p.gases.1)))
            .filter(|&(_, x)| x > 0.0)
            .collect();
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
            pairs,
            water,
            cloud: WATER_CLOUD_DEPTH,
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
    /// gas's frost nor the skin; the water vapour above each level, by the saturated vapour
    /// each holds, in the units of its ground's; and each level's saturated vapour's share of
    /// the pressure (none without water).
    fn profile(&self, t: f64, levels: &[f64]) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
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
        let (mut vapour, mut humid) = (vec![0.0; n], vec![0.0; n]);
        if let Some(water) = volatile("H2O").filter(|_| self.water > 0.0) {
            for (k, h) in humid.iter_mut().enumerate() {
                *h = vapour_pressure(water, temps[k]).min(0.5 * levels[k]) / levels[k];
            }
            for k in (0..n - 1).rev() {
                let mean = 0.5 * (humid[k] + humid[k + 1]);
                vapour[k] = vapour[k + 1] + mean * (levels[k] - levels[k + 1]);
            }
        }
        (temps, vapour, humid)
    }

    /// The outgoing radiation, W/m², of the column over ground at `t` K.
    pub fn olr(&self, t: f64) -> f64 {
        self.olr_at(t, 1.0)
    }

    /// The outgoing radiation, W/m², of the column over ground at `t` K standing at `ratio` of
    /// the surface pressure: the air above it, and of the water vapour the air over the datum
    /// at `t` holds (F235), the share its profile holds above it.
    pub fn olr_at(&self, t: f64, ratio: f64) -> f64 {
        self.streams_at(t, ratio).olr
    }

    /// The column's streams over ground at `t` K, W/m² (M55): what it lets out to space, the
    /// share of the ground's own radiation that reaches space through it, and what its air
    /// sends down to the ground.
    pub fn streams(&self, t: f64) -> Streams {
        self.streams_at(t, 1.0)
    }

    /// The column's streams over ground at `t` K standing at `ratio` of the surface pressure,
    /// as `olr_at` reads the column. The air's downward stream sums each layer's emission into
    /// the transmission below it to the ground, as the outgoing stream sums it into the
    /// transmission above it to space.
    fn streams_at(&self, t: f64, ratio: f64) -> Streams {
        let ground = self.pressure * ratio;
        if ground <= 0.0 || t <= 0.0 {
            return Streams {
                olr: SIGMA_SB * t.max(0.0).powi(4),
                window: 1.0,
                down: 0.0,
            };
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
        let (temps, vapour, humid) = self.profile(t, &levels);
        // Carbon dioxide's band's path above each level, kg/m² at the reference pressure, and
        // its molecules' collisions', amagat² cm.
        let band: Vec<f64> = levels
            .iter()
            .map(|&p| self.co2_mass * p * p / (2.0 * self.gravity * CO2_BAND_PRESSURE))
            .collect();
        let mut collisions = vec![0.0; n];
        let density = self.co2_volume * AMAGAT_TEMPERATURE / AMAGAT_PRESSURE;
        // Each layer's path for the molecules' collisions at unit shares, amagat² cm (M55).
        let mut pressed = vec![0.0; n - 1];
        for k in (0..n - 1).rev() {
            let (p0, p1) = (levels[k + 1], levels[k]);
            let mean = 0.5 * (p0 / temps[k + 1] + p1 / temps[k]);
            let unit = (AMAGAT_TEMPERATURE / AMAGAT_PRESSURE).powi(2) * 100.0 * self.gas_constant
                / self.gravity
                * mean
                * (p1 - p0);
            pressed[k] = unit;
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
        // The vapour's share of saturation, its column against the saturated column's
        // (M55): the saturated vapour's mass above a level eps / g of its profile's.
        let saturation = column * self.gravity / self.epsilon;
        // The clouds, well mixed in the air's mass, a u^b of the vapour u the datum's column
        // holds, kg/m² (M55).
        let held = column * datum;
        let clouds = if held > 0.0 {
            self.cloud.0 * held.powf(self.cloud.1)
        } else {
            0.0
        };
        let grey: Vec<f64> = levels
            .iter()
            .map(|&p| {
                let x = p / self.pressure;
                let dry: f64 = self.grey.iter().map(|&(tau, b)| tau * x.powf(b)).sum();
                dry + clouds * x
            })
            .collect();
        // Water vapour's paths through each layer (M55): its lines' broadened by the pressure,
        // kg/m² at the reference pressure, and its continuum's by its own pressure and its
        // temperature, kg/m² at the reference state.
        let (reference, warmth, sigma) = H2O_CONTINUUM_STATE;
        let water: Vec<(f64, f64)> = (0..n - 1)
            .map(|k| {
                let (p0, p1) = (levels[k], levels[k + 1]);
                let mass = column * (vapour[k] - vapour[k + 1]);
                let pressed = 0.5 * (humid[k] * p0 + humid[k + 1] * p1);
                let lines = column * pressed * (p0 - p1) / CO2_BAND_PRESSURE;
                let own = saturation * pressed / reference;
                let cold = (sigma * (warmth - 0.5 * (temps[k] + temps[k + 1]))).exp();
                (lines, mass * own * cold)
            })
            .collect();
        let spectrum = spectrum();
        // The depth above each level at each wavenumber: what lies over the top level, and
        // each layer's own depth summed down from it, each absorber's at the layer's own state.
        let layers: Vec<Layer> = (0..n - 1)
            .map(|k| Layer {
                band: band[k] - band[k + 1],
                collisions: collisions[k] - collisions[k + 1],
                grey: grey[k] - grey[k + 1],
                water: water[k].0,
                continuum: water[k].1,
                pairs: self
                    .pairs
                    .iter()
                    .map(|&(p, x)| {
                        let warmth = 0.5 * (temps[k] + temps[k + 1]);
                        let (j, w) = collision_tables().pairs[p].bracket(warmth);
                        (p, j, w, x * pressed[k])
                    })
                    .collect(),
            })
            .collect();
        let (mut total, mut through, mut ground_light, mut down) = (0.0, 0.0, 0.0, 0.0);
        let (mut depths, mut light) = (vec![0.0; n], vec![0.0; n]);
        for (i, &nu) in spectrum.nu.iter().enumerate() {
            depths[n - 1] = band[n - 1] * spectrum.band[i]
                + collisions[n - 1] * spectrum.collisions[i]
                + grey[n - 1];
            for k in (0..n - 1).rev() {
                depths[k] = depths[k + 1] + layers[k].depth(spectrum, i);
            }
            for (b, &temp) in light.iter_mut().zip(&temps) {
                *b = planck(nu, temp);
            }
            let (mut out, mut sent) = (0.0, 0.0);
            // Each layer emits up at its top and down at its bottom, its source linear in its
            // depth between its levels (Clough et al. 1992): the mean of its levels' where it is
            // thin, the emitting level's where it is thick.
            for k in 0..n - 1 {
                let x = DIFFUSIVITY * (depths[k] - depths[k + 1]);
                let (emits, slope) = linear_source(x);
                let up = light[k + 1] * emits + (light[k] - light[k + 1]) * slope;
                let dn = light[k] * emits + (light[k + 1] - light[k]) * slope;
                out += up * (-DIFFUSIVITY * depths[k + 1]).exp();
                sent += dn * (-DIFFUSIVITY * (depths[0] - depths[k])).exp();
            }
            let above = (-DIFFUSIVITY * depths[0]).exp();
            let own = planck(nu, t);
            out += own * above;
            total += out * spectrum.width[i];
            through += own * above * spectrum.width[i];
            ground_light += own * spectrum.width[i];
            down += sent * spectrum.width[i];
        }
        let pi = std::f64::consts::PI;
        Streams {
            olr: pi * total,
            window: through / ground_light.max(1e-300),
            down: pi * down,
        }
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

/// A layer of slant depth `x`'s emission at its edge per its source there, its emissivity
/// 1 - e^(-x), and per its source's change across it, (1 - e^(-x)) / x - e^(-x): its source
/// linear in its depth (Clough et al. 1992), x / 2 - x² / 3 where it is thin.
fn linear_source(x: f64) -> (f64, f64) {
    if x < 1e-4 {
        return (x - 0.5 * x * x, 0.5 * x - x * x / 3.0);
    }
    let through = (-x).exp();
    let emits = 1.0 - through;
    (emits, emits / x - through)
}

/// A layer of the column between two levels, as its absorbers' paths through it read it: carbon
/// dioxide's band's, kg/m² at the reference pressure, its molecules' collisions', amagat² cm, and
/// the grey depth.
struct Layer {
    band: f64,
    collisions: f64,
    grey: f64,
    /// Water vapour's lines' path, kg/m² at the reference pressure, and its continuum's, kg/m²
    /// at its reference state (M55).
    water: f64,
    continuum: f64,
    /// Each pair's collisions: its table, the table's temperatures about the layer's, the
    /// lower's index and the upper's share, and its path, amagat² cm (M55).
    pairs: Vec<(usize, usize, f64, f64)>,
}

impl Layer {
    /// The layer's optical depth at the spectrum's `i`th wavenumber.
    fn depth(&self, spectrum: &Spectrum, i: usize) -> f64 {
        let tables = &collision_tables().pairs;
        let pairs: f64 = self
            .pairs
            .iter()
            .map(|&(p, j, w, path)| {
                let k = &tables[p].absorption;
                let upper = k.get(j + 1).map_or(k[j][i], |row| row[i]);
                path * (k[j][i] + w * (upper - k[j][i]))
            })
            .sum();
        self.band * spectrum.band[i]
            + self.collisions * spectrum.collisions[i]
            + self.grey
            + self.water * spectrum.water[i]
            + self.continuum * spectrum.continuum[i]
            + pairs
    }
}

/// A column's streams over its ground (M55).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Streams {
    /// What it lets out to space, W/m².
    pub olr: f64,
    /// The share of the ground's own radiation that reaches space through it.
    pub window: f64,
    /// What its air sends down to the ground, W/m².
    pub down: f64,
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
/// tabulated over the logarithm of the temperature and read between the entries; with the share
/// of the ground's own radiation that reaches space through it and the share its air sends down
/// to the ground (M55).
#[derive(Clone, Debug)]
pub struct Outgoing {
    from: f64,
    step: f64,
    share: Vec<f64>,
    window: Vec<f64>,
    down: Vec<f64>,
}

impl Outgoing {
    /// None held back: an airless body's.
    pub fn none() -> Self {
        Self {
            from: 0.0,
            step: 1.0,
            share: vec![1.0],
            window: vec![1.0],
            down: vec![0.0],
        }
    }

    /// The column's shares from `low` to `high` K.
    pub fn of(column: &Column, low: f64, high: f64) -> Self {
        use rayon::prelude::*;
        let (from, to) = (low.max(1.0).ln(), high.max(low + 1.0).ln());
        let n = ((to - from) / OUTGOING_STEP).ceil() as usize + 1;
        let step = (to - from) / (n - 1) as f64;
        let streams: Vec<(f64, f64, f64)> = (0..n)
            .into_par_iter()
            .map(|k| {
                let t = (from + step * k as f64).exp();
                let black = SIGMA_SB * t.powi(4);
                let s = column.streams(t);
                (s.olr / black, s.window, s.down / black)
            })
            .collect();
        // Past water's runaway plateau, where its vapour's growth stalls the outgoing radiation
        // (Simpson 1927; Nakajima et al. 1992), held at its highest (M55).
        let mut highest = 0.0f64;
        let share = streams
            .iter()
            .enumerate()
            .map(|(k, s)| {
                let black = SIGMA_SB * (from + step * k as f64).exp().powi(4);
                highest = highest.max(s.0 * black);
                highest / black
            })
            .collect();
        Self {
            from,
            step,
            share,
            window: streams.iter().map(|s| s.1).collect(),
            down: streams.iter().map(|s| s.2).collect(),
        }
    }

    /// The table from its entries, (temperature, share), its window and downward shares
    /// those of a column that lets the ground's radiation through and sends none down until
    /// they are given (`with_streams`).
    pub fn from_points(points: &[(f64, f64)]) -> Self {
        let share: Vec<f64> = points.iter().map(|p| p.1).collect();
        let (from, step) = match points {
            [] => return Self::none(),
            [only] => (only.0.max(1.0).ln(), 1.0),
            _ => {
                let from = points[0].0.max(1.0).ln();
                let last = points[points.len() - 1].0.max(1.0).ln();
                (from, (last - from) / (points.len() - 1) as f64)
            }
        };
        Self {
            from,
            step,
            window: vec![1.0; share.len()],
            down: vec![0.0; share.len()],
            share,
        }
    }

    /// The table with its window's and downward shares' entries, (temperature, share), on its
    /// own temperatures.
    pub fn with_streams(mut self, window: &[(f64, f64)], down: &[(f64, f64)]) -> Self {
        if window.len() == self.share.len() && down.len() == self.share.len() {
            self.window = window.iter().map(|p| p.1).collect();
            self.down = down.iter().map(|p| p.1).collect();
        }
        self
    }

    /// The entries, (temperature, share).
    pub fn points(&self) -> Vec<(f64, f64)> {
        self.entries(&self.share)
    }

    /// The window's entries, (temperature, share of the ground's radiation reaching space).
    pub fn window_points(&self) -> Vec<(f64, f64)> {
        self.entries(&self.window)
    }

    /// The downward entries, (temperature, share of a black body's radiation sent down).
    pub fn down_points(&self) -> Vec<(f64, f64)> {
        self.entries(&self.down)
    }

    fn entries(&self, values: &[f64]) -> Vec<(f64, f64)> {
        values
            .iter()
            .enumerate()
            .map(|(k, &s)| ((self.from + self.step * k as f64).exp(), s))
            .collect()
    }

    /// The share at `t` K, held at the table's ends past them.
    pub fn share(&self, t: f64) -> f64 {
        self.read(&self.share, t)
    }

    /// The share of the ground's own radiation that reaches space through the column over
    /// ground at `t` K.
    pub fn window(&self, t: f64) -> f64 {
        self.read(&self.window, t)
    }

    /// What the column's air over ground at `t` K sends down to it, W/m².
    pub fn down(&self, t: f64) -> f64 {
        self.read(&self.down, t) * SIGMA_SB * t.powi(4)
    }

    /// What the column's air over ground at `t` K lets out to space of its own, W/m²: the
    /// outgoing radiation less the ground's through the window.
    pub fn air_up(&self, t: f64) -> f64 {
        (self.share(t) - self.window(t)) * SIGMA_SB * t.powi(4)
    }

    fn read(&self, values: &[f64], t: f64) -> f64 {
        if values.len() == 1 {
            return values[0];
        }
        let x = ((t.max(1e-9).ln() - self.from) / self.step).clamp(0.0, (values.len() - 1) as f64);
        let k = (x.floor() as usize).min(values.len() - 2);
        let f = x - k as f64;
        values[k] + f * (values[k + 1] - values[k])
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

/// The share of its star's light over the lit hemisphere a layer of grey depth `tau` absorbs on
/// the light's way down (M55), its mean over the sun's heights, 1 - 2 E3(tau).
pub fn hemisphere_absorbed(tau: f64) -> f64 {
    if tau <= 0.0 {
        return 0.0;
    }
    let n = HEMISPHERE_HEIGHTS;
    let through: f64 = (0..n)
        .map(|k| {
            let mu = (k as f64 + 0.5) / n as f64;
            2.0 * mu * (-tau / mu).exp()
        })
        .sum::<f64>()
        / n as f64;
    1.0 - through
}

/// The share of its star's light a methane haze of index `haze` absorbs high over the air and
/// its stratosphere sends back to space (M55; McKay et al. 1991's anti-greenhouse), over the lit
/// hemisphere.
pub fn haze_returned(haze: f64) -> f64 {
    HAZE_RETURNED * hemisphere_absorbed(HAZE_SUNLIGHT * haze)
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
        let s = base.streams(93.7);
        let black = SIGMA_SB * 93.7f64.powi(4);
        println!(
            "Titan at 93.7 K: olr {:.3}, window {:.4} (3.8%), down {:.3} (3.92 to 3.98), net at the ground {:.3} (0.379); the mean under all {absorbed:.3} W/m2 absorbed {:.2} K",
            s.olr,
            s.window,
            s.down,
            black - s.down,
            base.surface_temperature(absorbed)
        );
        // The haze's depth in the sunlight: its stratosphere's 1.5 of 3.6 W/m2; and the share
        // of it sent back to space that keeps Titan's 93.7 K.
        let index = tholin_production(
            entry(&titan, "atmosphere.composition", "CH4"),
            titan.f64("star.xuv_relative"),
            titan.f64("atmosphere.surface_pressure"),
        );
        let tau = halve(|t| hemisphere_absorbed(t) - 1.5 / 3.6, 1e-6, 1e2);
        let light = absorbed / (1.0 - 0.265);
        let high = hemisphere_absorbed(tau) * light;
        let share = halve(
            |f| 93.7 - base.surface_temperature(absorbed - f * high),
            1e-4,
            1.0,
        );
        println!(
            "Titan: haze index {index:.5}, its depth in the sunlight {tau:.5}: HAZE_SUNLIGHT {:.5}; it absorbs {high:.3} of {light:.3} W/m2, sends {share:.4} of it back (HAZE_RETURNED): {:.3} W/m2; the anti-greenhouse {:.2} K",
            tau / index,
            share * high,
            base.surface_temperature(absorbed) - 93.7
        );
        let earth = graph("Earth");
        let (mut base, absorbed) = column_of(&earth, rayleigh_of(&earth));
        base.cloud = (0.0, 1.0);
        let clear = base.surface_temperature(absorbed);
        let s = base.streams(288.0);
        println!(
            "Earth clear: {clear:.2} K; at 288 K olr {:.2} (CERES 267), down {:.2} (314), window {:.4}, slope {:.3}",
            s.olr,
            s.down,
            s.window,
            (base.olr(289.0) - base.olr(287.0)) / 2.0
        );
        // The clouds' a that keeps Earth's zero-dimensional 291.6 K, and their b that keeps its
        // outgoing radiation's slope, 2.1 W/m² per kelvin.
        let with = |cloud: (f64, f64)| {
            let mut c = base.clone();
            c.cloud = cloud;
            c
        };
        let at = |b: f64| {
            let a = halve(
                |a| with((a, b)).surface_temperature(absorbed) - 291.6,
                1e-6,
                1e3,
            );
            let c = with((a, b));
            (a, (c.olr(292.6) - c.olr(290.6)) / 2.0)
        };
        for b in [0.0, 0.1, 0.2, 0.3, 0.5] {
            let (a, slope) = at(b);
            let c = with((a, b));
            println!(
                "  clouds a u^{b}: a {a:.4}, slope {slope:.3}, olr at 250 K {:.1}, 230 K {:.1}",
                c.olr(250.0),
                c.olr(230.0)
            );
        }
        let b = 1.0 - halve(|x| at(1.0 - x).1 - 2.1, 1e-3, 1.0);
        let (k, _) = at(b);
        let m0 = b;
        let (c, cs) = (with((k, b)), base.streams(291.6));
        let s = c.streams(291.6);
        println!(
            "Earth: albedo {:.4}, clouds {k:.5} u^{m0:.4}: olr {:.2} (clear {:.2}: effect {:.2}, CERES 28), down {:.2} (clear {:.2}; 342, 314), window {:.4}, slope {:.3}",
            rayleigh_of(&earth),
            s.olr,
            cs.olr,
            cs.olr - s.olr,
            s.down,
            cs.down,
            s.window,
            (c.olr(292.6) - c.olr(290.6)) / 2.0
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
        let (base, absorbed) = column_of(&venus, a_v);
        let so2 = base.grey[0].0 / SULFURIC_CLOUD_DEPTH;
        let k = halve(
            |k| {
                let mut c = base.clone();
                c.grey[0].0 = k * so2;
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

    /// M55: the column's streams meet what is measured. Earth's clear air at 288 K lets out
    /// CERES's 267 W/m² within 8 (Wild et al. 2019); with its stated clouds its air over its
    /// 291.6 K ground sends down the measured 342 W/m² within 10 (Kato et al. 2018 give 6 of
    /// uncertainty), its window between the cloudy 0.055 and the clear 0.17 of its ground's
    /// light (Costa and Shine 2012). Titan's air over its 93.7 K ground lets 2 to 6 percent of
    /// the ground's light through to space (3.8 on the Huygens probe's profile), sends down
    /// 3.7 to 4.1 W/m², and leaves its ground losing 0.2 to 0.6 W/m² net (0.379, McKay et al.
    /// 1989), where the grey column's window was 0.30. Its haze sends back to space between
    /// half and all of the 1.5 of 3.6 W/m² its stratosphere takes (McKay et al. 1991); Earth's
    /// air none. A clear air lets all its ground's light through and sends none down.
    #[test]
    fn the_air_s_streams_meet_the_measured() {
        let earth = || {
            let g = graph("Earth");
            let (column, _) = column_of(&g, g.f64("atmosphere.bond_albedo"));
            let mut clear = column.clone();
            clear.cloud = (0.0, 0.0);
            let (s, c) = (column.streams(291.6), clear.streams(288.0));
            assert!((c.olr - 267.0).abs() < 8.0, "clear {:.1}", c.olr);
            assert!((s.down - 342.0).abs() < 10.0, "down {:.1}", s.down);
            assert!((0.055..0.17).contains(&s.window), "window {:.3}", s.window);
            assert_eq!(g.f64("atmosphere.sunlight_returned"), 0.0);
        };
        let titan = || {
            let g = graph("Titan");
            let (column, _) = column_of(&g, g.f64("atmosphere.bond_albedo"));
            let s = column.streams(93.7);
            let net = SIGMA_SB * 93.7f64.powi(4) - s.down;
            assert!((0.02..0.06).contains(&s.window), "window {:.4}", s.window);
            assert!((3.7..4.1).contains(&s.down), "down {:.3}", s.down);
            assert!((0.2..0.6).contains(&net), "net {net:.3}");
            let returned = g.f64("atmosphere.sunlight_returned") * 3.6;
            assert!(
                (0.75..1.5).contains(&returned),
                "returned {returned:.3} of 3.6 W/m²"
            );
        };
        let clear = || {
            let argon: BTreeMap<String, f64> = [("Ar".to_string(), 1.0)].into_iter().collect();
            let s = Column::new(&argon, 1e5, 9.81, 0.0, 200.0).streams(288.0);
            assert!(
                (s.window - 1.0).abs() < 1e-9 && s.down.abs() < 1e-9,
                "{s:?}"
            );
        };
        crate::testkit::checks(&[
            ("Earth's streams", &earth),
            ("Titan's streams", &titan),
            ("a clear air's streams", &clear),
        ]);
    }

    /// M55: Earth's outgoing radiation over its zonal means' temperatures, against the zonal
    /// law Graves, Lee and North 1993 fit, 203.3 + 2.09 (T - 273.15) W/m².
    #[test]
    #[ignore]
    fn print_earth_s_outgoing_radiation() {
        let g = graph("Earth");
        let (column, _) = column_of(&g, g.f64("atmosphere.bond_albedo"));
        let mut clear = column.clone();
        clear.cloud = (0.0, 1.0);
        for t in [
            220.0, 230.0, 240.0, 250.0, 260.0, 270.0, 280.0, 291.6, 300.0, 310.0,
        ] {
            let law = 203.3 + 2.09 * (t - 273.15);
            println!(
                "{t:.1} K: olr {:.1} (clear {:.1}), zonal law {law:.1}; down {:.1}; vapour {:.2} kg/m2",
                column.olr(t),
                clear.olr(t),
                column.streams(t).down,
                water_column(t, column.water)
            );
        }
    }

    /// M55: the streams of each body's column at its mean temperature.
    #[test]
    #[ignore]
    fn print_the_streams() {
        for name in ["Earth", "Titan", "Venus", "Mars", "Triton", "Pluto"] {
            let g = graph(name);
            let (column, absorbed) = column_of(&g, g.f64("atmosphere.bond_albedo"));
            let t = g.f64("atmosphere.surface_temperature_mean");
            for dt in [-10.0, 0.0, 10.0] {
                let s = column.streams(t + dt);
                let b = SIGMA_SB * (t + dt).powi(4);
                println!(
                    "{name} {:.1} K: absorbed {absorbed:.2}, black {b:.2}, olr {:.3}, window {:.4}, up from the air {:.3}, down {:.3} ({:.4} of black), net at the ground {:.3}",
                    t + dt,
                    s.olr,
                    s.window,
                    s.olr - s.window * b,
                    s.down,
                    s.down / b,
                    b - s.down
                );
            }
        }
    }

    /// F295: the column's spectrum holds a black body's light: an air that neither absorbs nor
    /// emits lets out the ground's sigma T⁴, within a percent, from the coldest ground a
    /// body's table reads, 5 K, to 2,000 K, where a hot world's thin air would stand (the
    /// tests' tiers); the spectrum from 1 to 5,000 cm^-1 let out 32% too much at 5 K and lost
    /// a twentieth at 1,000 K, a half at 2,000 K.
    #[test]
    fn a_clear_air_lets_out_a_black_body_s_light() {
        let clear: BTreeMap<String, f64> = [("Ar".to_string(), 1.0)].into_iter().collect();
        let column = Column::new(&clear, 100.0, 9.81, 0.0, 100.0);
        let shares: Vec<(f64, f64)> = (0..=40)
            .map(|k| {
                let t = 5.0 * 400f64.powf(k as f64 / 40.0);
                (t, column.olr(t) / (SIGMA_SB * t.powi(4)))
            })
            .collect();
        for &(t, share) in &shares {
            println!("{t:7.1} K: {share:.5}");
        }
        for &(t, share) in &shares {
            assert!((share - 1.0).abs() < 0.01, "{t:.1} K lets out {share:.4}");
        }
    }

    /// Conservative scattering layers stack as one of their summed depth: a layer reflecting
    /// s / (1 + s) of the light, its two-stream law, over one reflecting q / (1 + q) reflects
    /// (s + q) / (1 + s + q), in either order; over a black ground a layer reflects its own
    /// share, over a white one all (the tests' tiers).
    #[test]
    fn scattering_layers_stack_as_one_of_their_summed_depth() {
        let r = |s: f64| s / (1.0 + s);
        let depths = [0.0, 1e-3, 0.05, 0.3, 1.0, 4.0, 50.0];
        for &s in &depths {
            assert!((over(r(s), 0.0) - r(s)).abs() < 1e-12, "{s}");
            assert!((over(r(s), 1.0) - 1.0).abs() < 1e-12, "{s}");
            for &q in &depths {
                let (one, other) = (over(r(s), r(q)), over(r(q), r(s)));
                assert!((one - r(s + q)).abs() < 1e-12, "{s} {q}: {one}");
                assert!((one - other).abs() < 1e-12, "{s} {q}: {one} {other}");
            }
        }
    }

    /// An air's molecules scatter each on its own: its Rayleigh depth is its gases' depths
    /// each over its own column, and grows as the column; what it reflects of its star's light
    /// grows with the depth and with the star's warmth, bluer, never reaching all of it (the
    /// tests' tiers).
    #[test]
    fn an_air_s_molecules_scatter_each_on_its_own() {
        let mix: BTreeMap<String, f64> = [("N2", 0.6), ("CO2", 0.3), ("CH4", 0.1)]
            .into_iter()
            .map(|(s, f)| (s.to_string(), f))
            .collect();
        let column = 1e4;
        let mu = mean_molecular_mass(&mix);
        let depth = rayleigh_depth(&mix, column);
        let parts: f64 = mix
            .iter()
            .map(|(s, &f)| {
                let pure: BTreeMap<String, f64> = [(s.clone(), 1.0)].into_iter().collect();
                rayleigh_depth(&pure, column * f * molar_mass(s) / mu)
            })
            .sum();
        assert!((parts - depth).abs() < 1e-12 * depth, "{parts} {depth}");
        let twice = rayleigh_depth(&mix, 2.0 * column);
        assert!(
            (twice - 2.0 * depth).abs() < 1e-12 * depth,
            "{twice} {depth}"
        );
        let mut last = 0.0;
        for tau in [1e-4, 1e-2, 0.1, 1.0, 10.0, 100.0] {
            let r = rayleigh_reflectance(tau, 5_772.0);
            assert!(r > last && r < 1.0, "{tau}: {r} after {last}");
            last = r;
        }
        let mut last = 0.0;
        for teff in [2_500.0, 3_500.0, 5_772.0, 10_000.0] {
            let r = rayleigh_reflectance(0.1, teff);
            assert!(r > last && r < 1.0, "{teff} K: {r} after {last}");
            last = r;
        }
    }

    /// The outgoing radiation rises with the ground's warmth and falls as more air stands
    /// over it: dry airs of carbon dioxide and nitrogen, from a thin Martian air to a
    /// Venusian, each from 150 to 1,000 K, never let out more over a colder ground, to
    /// rounding, nor under a thicker air over a ground its air may stand on, at or over its
    /// bulk gas's frost point; under that the air would freeze out on the ground, and a
    /// thicker air's frost point holds its column the warmer (the tests' tiers).
    #[test]
    fn the_outgoing_radiation_rises_with_the_warmth_and_falls_with_the_air() {
        for x in [4.2e-4, 0.05, 0.95] {
            let composition: BTreeMap<String, f64> =
                [("N2".to_string(), 1.0 - x), ("CO2".to_string(), x)]
                    .into_iter()
                    .collect();
            let temperatures: Vec<f64> = (0..=24)
                .map(|k| 150.0 * (1_000f64 / 150.0).powf(k as f64 / 24.0))
                .collect();
            let (olr, floors): (Vec<Vec<f64>>, Vec<f64>) = [600.0, 1e4, 1e5, 1e6, 9.2e6]
                .iter()
                .map(|&p| {
                    let column = Column::new(&composition, p, 9.81, 0.0, 150.0);
                    let floor = column
                        .bulk
                        .map_or(0.0, |(v, share)| saturation_temperature(v, share * p));
                    let row = temperatures.iter().map(|&t| column.olr(t)).collect();
                    (row, floor)
                })
                .unzip();
            for (i, row) in olr.iter().enumerate() {
                for k in 1..row.len() {
                    assert!(
                        row[k] >= row[k - 1] * (1.0 - 1e-12),
                        "CO2 {x}, air {i}: {:.1} K lets out {:.3}, {:.1} K {:.3}",
                        temperatures[k - 1],
                        row[k - 1],
                        temperatures[k],
                        row[k]
                    );
                    if i > 0 && temperatures[k] >= floors[i] {
                        assert!(
                            row[k] <= olr[i - 1][k] * (1.0 + 1e-9),
                            "CO2 {x}, {:.1} K: air {i} lets out {:.3}, the thinner {:.3}",
                            temperatures[k],
                            row[k],
                            olr[i - 1][k]
                        );
                    }
                }
            }
        }
    }
}
