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
    heat_capacity, mean_heat_capacity, mean_molecular_mass, molar_mass, saturation_temperature,
    vapour_pressure, volatile, water_column,
};

/// Planck's radiance at a wavenumber, cm^-1, and a temperature, K: W/m²/sr per cm^-1.
pub fn planck(nu: f64, t: f64) -> f64 {
    if t <= 0.0 {
        return 0.0;
    }
    RADIATION_C1 * nu.powi(3) / (RADIATION_C2_CM_K * nu / t).exp_m1()
}

/// Planck's function's change per kelvin at `nu` cm^-1 and `t` K.
fn planck_slope(nu: f64, t: f64) -> f64 {
    let x = RADIATION_C2_CM_K * nu / t.max(1e-300);
    if t <= 0.0 || x > 700.0 {
        return 0.0;
    }
    let less = x.exp_m1();
    RADIATION_C1 * nu.powi(3) / less * x * (less + 1.0) / less / t
}

/// The wavenumbers the column's radiation is summed over, cm^-1, with their widths: finely
/// over carbon dioxide's band, where its absorption changes by a factor e every 12 cm^-1,
/// coarsely elsewhere; with carbon dioxide's band's absorption per kg/m² at the reference
/// pressure at each.
struct Spectrum {
    nu: Vec<f64>,
    width: Vec<f64>,
    band: Vec<f64>,
    /// Water vapour's lines per kg/m² at the reference pressure, and its continuum per kg/m²
    /// at its reference state (M55).
    water: Vec<f64>,
    continuum: Vec<f64>,
    /// Methane's 7.7 um band per kg/m² at the reference pressure, and the haze's thermal depth
    /// per its depth at HAZE_INFRARED's reference wavenumber (M65).
    methane: Vec<f64>,
    haze: Vec<f64>,
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
        // Methane's band by its medians' table, each 25 cm^-1's (M65).
        let methane = nu
            .iter()
            .map(|&n| {
                CH4_BAND
                    .iter()
                    .find(|&&(to, _)| n < to)
                    .map_or(0.0, |&(_, k)| k)
            })
            .collect();
        // The haze absorbing in the far infrared in proportion to the wavenumber (M65).
        let haze = nu.iter().map(|&n| n / HAZE_INFRARED_REFERENCE).collect();
        Spectrum {
            nu,
            width,
            band,
            water,
            continuum,
            methane,
            haze,
        }
    })
}

/// The gases' collisions (M55): each pair's absorption, cm^-1 amagat^-2, on the spectrum's
/// wavenumbers at each of its tables' temperatures, from HITRAN's tables coarse-grained by 10
/// cm^-1 (data/collisions.txt): nitrogen's (Chistikov et al. 2019), nitrogen and methane's
/// (Finenko et al. 2022, fitted to Cassini's spectra of Titan), methane's (Borysow and
/// Frommhold 1987), nitrogen and hydrogen's (Borysow and Frommhold 1986); carbon dioxide's over
/// its five spans from the far infrared to 4,490 cm^-1 (M65: Gruszka and Borysow 1997, Baranov
/// et al. 2003 and 2004, Terragni et al. 2025, Tran et al. 2024).
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
        bracket(&self.temperatures, t)
    }
}

/// A table's temperatures `ts` about `t`, K: the lower's index and the share of the upper's,
/// held at the ends.
fn bracket(ts: &[f64], t: f64) -> (usize, f64) {
    if ts.len() < 2 || t <= ts[0] {
        return (0, 0.0);
    }
    let j = ts.partition_point(|&x| x <= t).clamp(1, ts.len() - 1) - 1;
    (j, ((t - ts[j]) / (ts[j + 1] - ts[j])).clamp(0.0, 1.0))
}

/// A gas's lines (M65), from HITRAN's (data/*_lines.txt): on each of the spectrum's
/// wavenumbers from the table's first, at each of its temperatures, the lines' median absorption
/// within 25 cm^-1 of their centres at the reference pressure, m²/kg, broadened in proportion to
/// the pressure, and their mean, the most a bin's broadened lines absorb; and their far wings,
/// sub-Lorentzian, also in proportion to the pressure, without bound (carbon dioxide's).
struct Lines {
    temperatures: Vec<f64>,
    median: Table,
    mean: Table,
    wing: Table,
    /// The first of the spectrum's wavenumbers the table holds.
    from: usize,
}

/// A line table's absorption at each of its temperatures on each of the spectrum's
/// wavenumbers, with its logarithms, read between its temperatures.
struct Table {
    values: Vec<Vec<f64>>,
    logs: Vec<Vec<f64>>,
}

impl Table {
    fn new(values: Vec<Vec<f64>>) -> Self {
        let logs = values
            .iter()
            .map(|row| {
                row.iter()
                    .map(|&v| if v > 0.0 { v.ln() } else { 0.0 })
                    .collect()
            })
            .collect();
        Self { values, logs }
    }

    /// The absorption at the spectrum's `i`th wavenumber between the `j`th and next
    /// temperatures, `w` of the way, read in the logarithm.
    fn at(&self, j: usize, w: f64, i: usize) -> f64 {
        let lo = self.values[j][i];
        let hi = self.values.get(j + 1).map_or(lo, |row| row[i]);
        if lo > 0.0 && hi > 0.0 {
            let a = self.logs[j][i];
            let b = self.logs.get(j + 1).map_or(a, |row| row[i]);
            (a + w * (b - a)).exp()
        } else {
            lo + w * (hi - lo)
        }
    }
}

/// The gases' lines' tables (M65): carbon dioxide's from 400 cm^-1, its far wings, and from
/// 950 its lines outside its 15 um band (CO2_BAND); water's past its fitted bands
/// (H2O_LINES_END); sulfur dioxide's.
fn line_tables() -> &'static [Lines; 3] {
    static LINES: OnceLock<[Lines; 3]> = OnceLock::new();
    LINES.get_or_init(|| {
        [
            Lines::read(include_str!("../data/co2_lines.txt")),
            Lines::read(include_str!("../data/h2o_lines.txt")),
            Lines::read(include_str!("../data/so2_lines.txt")),
        ]
    })
}

impl Lines {
    /// The table in `text`: its first wavenumber, its temperatures, and each row its bin's
    /// upper edge, then the median, the mean and the wings at each temperature.
    fn read(text: &str) -> Self {
        let spectrum = spectrum();
        let mut temperatures: Vec<f64> = Vec::new();
        let mut first = f64::INFINITY;
        let mut rows: Vec<(f64, Vec<f64>)> = Vec::new();
        for line in text.lines() {
            let mut words = line.split_whitespace();
            match words.next() {
                None => {}
                Some(w) if w.starts_with('#') => {}
                Some("T") => temperatures = words.filter_map(|w| w.parse().ok()).collect(),
                Some("from") => {
                    first = words.next().and_then(|w| w.parse().ok()).unwrap_or(first);
                }
                Some(edge) => {
                    if let Ok(edge) = edge.parse::<f64>() {
                        rows.push((edge, words.filter_map(|w| w.parse().ok()).collect()));
                    }
                }
            }
        }
        let from = spectrum.nu.partition_point(|&nu| nu < first);
        let pick = |m: usize| -> Vec<Vec<f64>> {
            (0..temperatures.len())
                .map(|j| {
                    spectrum
                        .nu
                        .iter()
                        .map(|&nu| {
                            if nu < first {
                                return 0.0;
                            }
                            rows.iter()
                                .find(|r| nu < r.0)
                                .map_or(0.0, |r| r.1.get(3 * j + m).copied().unwrap_or(0.0))
                        })
                        .collect()
                })
                .collect()
        };
        Lines {
            median: Table::new(pick(0)),
            mean: Table::new(pick(1)),
            wing: Table::new(pick(2)),
            temperatures,
            from,
        }
    }

    /// The depth of the gas's lines at the spectrum's `i`th wavenumber over a layer, by its
    /// temperatures about the layer's, `j` and `w`, and its paths, kg/m², `broadened` at the
    /// reference pressure and its own `mass`.
    fn depth(&self, i: usize, (j, w): (usize, f64), broadened: f64, mass: f64) -> f64 {
        if i < self.from {
            return 0.0;
        }
        (broadened * self.median.at(j, w, i)).min(mass * self.mean.at(j, w, i))
            + broadened * self.wing.at(j, w, i)
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
            // Each wavenumber's absorption between the table's, none outside its span (M65:
            // a span may start above the spectrum's first wavenumber).
            let at = |j: usize, nu: f64| {
                let k = rows.partition_point(|r| r.0 <= nu);
                if k == 0 || k >= rows.len() {
                    return 0.0;
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

/// The sunlight a column absorbs and the depths its air takes it by (M65): the light that
/// reaches its top and what the planet absorbs of it at its Bond albedo, W/m² over its whole
/// surface, and the haze's and methane's grey depths over the column in the light; the clouds'
/// and water vapour's are the column's own, from its vapour.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Sunlight {
    pub incident: f64,
    pub absorbed: f64,
    pub haze: f64,
    pub methane: f64,
}

impl Sunlight {
    /// The sunlight an air of `composition` at `pressure`, Pa, under `gravity` absorbs of
    /// `incident` W/m² over its whole surface at its Bond `albedo`, its star's ultraviolet
    /// `xuv` times the Sun's making its haze (M65): the haze's depth HAZE_SUNLIGHT for each
    /// unit of its index, methane's METHANE_SUNLIGHT per kg/m² of it, as the climate's split
    /// reads them.
    pub fn of(
        composition: &BTreeMap<String, f64>,
        (pressure, gravity): (f64, f64),
        incident: f64,
        albedo: f64,
        xuv: f64,
    ) -> Self {
        let total: f64 = composition.values().sum::<f64>().max(1e-300);
        let methane = composition.get("CH4").copied().unwrap_or(0.0) / total;
        let mass = methane * molar_mass("CH4") / mean_molecular_mass(composition);
        Self {
            incident,
            absorbed: incident * (1.0 - albedo).max(0.0),
            haze: HAZE_SUNLIGHT * tholin_production(methane, xuv, pressure),
            methane: METHANE_SUNLIGHT * mass * pressure / gravity.max(1e-300),
        }
    }
}

/// A column's radiative-convective profiles at its levels over the datum, at the nodes of its
/// ground's temperature found so far (M65).
type Nodes = std::sync::Arc<std::sync::Mutex<Vec<(i64, Vec<f64>)>>>;

/// The nodes found for the column `key` names (M65): kept for every column alike, so that the
/// graph's columns over one air, and its recomputes, solve its profiles once; the last
/// RC_COLUMNS_KEPT columns' kept.
fn node_store(key: Vec<u64>) -> Nodes {
    use std::sync::Mutex;
    static STORE: Mutex<Vec<(Vec<u64>, Nodes)>> = Mutex::new(Vec::new());
    let mut store = STORE.lock().unwrap_or_else(|e| e.into_inner());
    if let Some((_, nodes)) = store.iter().find(|(k, _)| *k == key) {
        return nodes.clone();
    }
    if store.len() >= RC_COLUMNS_KEPT {
        store.remove(0);
    }
    let nodes = Nodes::default();
    store.push((key, nodes.clone()));
    nodes
}

/// A column of air over the ground (F293), its profile radiative-convective (M65).
#[derive(Clone, Debug)]
pub struct Column {
    /// The surface pressure, Pa, and the gravity, m/s².
    pressure: f64,
    gravity: f64,
    /// The air's gas constant, J/kg/K, and its gases' stated heat capacity, J/kg/K, carbon
    /// dioxide's read by its temperature (M65); water's molecules' mass against the air's.
    gas_constant: f64,
    heat: f64,
    epsilon: f64,
    /// Carbon dioxide's share of the air's mass and of its molecules; methane's of its mass
    /// (M65).
    co2_mass: f64,
    co2_volume: f64,
    ch4_mass: f64,
    /// Sulfur dioxide's share of the air's mass (M65).
    so2_mass: f64,
    /// The bulk gas, whose frost holds the air, and its share of the molecules.
    bulk: Option<(&'static Volatile, f64)>,
    /// The grey skin's temperature over the air, K, the profiles' cold start.
    skin: f64,
    /// The pairs of its gases whose collisions absorb, each its index among the tables' and
    /// the product of the two gases' shares of the molecules (M55).
    pairs: Vec<(usize, f64)>,
    /// The water a warm air can hold vapour from, metres of global layer (F235).
    water: f64,
    /// The temperature whose column the air's vapour holds at most, K (M59): a body whose
    /// water has no liquid to evaporate draws its vapour off its ice, however warm the ground
    /// under the column; unbounded where its water is liquid.
    vapour_most: f64,
    /// The clouds' grey depth's law, a u^b of the column's vapour u, kg/m² (M55).
    cloud: (f64, f64),
    /// The haze's thermal depth over the column at HAZE_INFRARED's reference wavenumber (M65).
    haze_ir: f64,
    /// The sulfuric deck its sulfur dioxide makes (M65): its infrared and its sunlight's grey
    /// depths and its acid vapour's share of the molecules.
    deck: (f64, f64, f64),
    /// The sunlight it absorbs (M65).
    sunlight: Sunlight,
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
        let share = |s: &str| composition.get(s).copied().unwrap_or(0.0) / total;
        let mass = |s: &str| share(s) * molar_mass(s) / mu;
        let column = pressure / gravity.max(1e-300);
        let bulk = bulk_of(composition);
        let u = |s: &str| mass(s) * column;
        let deck = (
            SULFURIC_CLOUD_DEPTH * u("SO2"),
            SULFURIC_CLOUD_SUNLIGHT * u("SO2"),
            SULFURIC_ACID_SHARE * share("SO2"),
        );
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
            epsilon: molar_mass("H2O") / mu,
            co2_mass: mass("CO2"),
            co2_volume: share("CO2"),
            ch4_mass: mass("CH4"),
            so2_mass: mass("SO2"),
            bulk,
            skin,
            pairs,
            water,
            vapour_most: f64::INFINITY,
            cloud: WATER_CLOUD_DEPTH,
            haze_ir: 0.0,
            deck,
            sunlight: Sunlight::default(),
        }
    }

    /// The column with its vapour no more than the vapour law's column at `most` K, its
    /// profile's and its adiabat's no more than saturated at it (M59).
    pub fn with_vapour_most(mut self, most: f64) -> Self {
        self.vapour_most = most;
        self
    }

    /// The column absorbing `sunlight` (M65): its air taking it by its absorbers' depths, its
    /// haze's thermal depth HAZE_INFRARED for each unit of its index.
    pub fn with_sunlight(mut self, sunlight: Sunlight) -> Self {
        self.sunlight = sunlight;
        self.haze_ir = HAZE_INFRARED * sunlight.haze / HAZE_SUNLIGHT;
        self
    }

    /// The column without the clouds' stated law (M66): its clear air, whose clouds the
    /// climate lays by place and season.
    pub fn without_clouds(mut self) -> Self {
        self.cloud = (0.0, 0.0);
        self
    }

    /// The sunlight the column absorbs (M65).
    pub fn sunlight(&self) -> Sunlight {
        self.sunlight
    }

    /// The air's heat capacity at `t` K, J/kg/K: its gases' stated ones, carbon dioxide's by
    /// its temperature (M65).
    fn heat_at(&self, t: f64) -> f64 {
        let (a, b, c, d, e) = CO2_SHOMATE;
        let x = t.max(200.0) / 1000.0;
        let molar = a + b * x + c * x * x + d * x * x * x + e / (x * x);
        self.heat + self.co2_mass * (molar / molar_mass("CO2") - heat_capacity("CO2"))
    }

    /// The adiabat's exponent, d ln T / d ln p, at `t` K and `p` Pa: the dry air's, slowed where
    /// the column holds water by its saturated vapour's latent heat,
    /// (1 + L r / (R T)) / (1 + L² r eps / (c_p R T²)), r its mixing ratio; its heat capacity
    /// its temperature's (M65).
    fn lapse(&self, t: f64, p: f64) -> f64 {
        let heat = self.heat_at(t);
        let dry = self.gas_constant / heat;
        let Some(water) = volatile("H2O").filter(|_| self.water > 0.0) else {
            return dry;
        };
        let e = vapour_pressure(water, t.min(self.vapour_most)).min(0.5 * p);
        // Its mixing ratio saturated, no more than the water the air can hold all of makes
        // (M65): Venus's air, far from saturated, rises on its dry adiabat.
        let held = self.water * 1_000.0 * self.gravity / self.pressure.max(1e-300);
        let r = (self.epsilon * e / (p - e)).min(held);
        let l = if t >= water.t_triple {
            water.l_vap
        } else {
            water.l_sub
        };
        let rt = self.gas_constant * t;
        dry * (1.0 + l * r / rt) / (1.0 + l * l * r * self.epsilon / (heat * rt * t))
    }

    /// The bulk gas's frost point at `p` Pa, K; none without one.
    fn frost(&self, p: f64) -> f64 {
        self.bulk
            .map_or(0.0, |(v, x)| saturation_temperature(v, x * p))
    }

    /// The adiabat over ground at `t` K through `levels`, the first level the ground's: by
    /// midpoint steps in the pressure's logarithm.
    fn adiabat_over(&self, t: f64, levels: &[f64]) -> Vec<f64> {
        let n = levels.len();
        let mut adiabat = vec![t; n];
        for k in 0..n - 1 {
            let step = (levels[k + 1] / levels[k]).ln();
            let mid = adiabat[k] * (0.5 * step * self.lapse(adiabat[k], levels[k])).exp();
            let p_mid = (levels[k] * levels[k + 1]).sqrt();
            adiabat[k + 1] = adiabat[k] * (step * self.lapse(mid, p_mid)).exp();
        }
        adiabat
    }

    /// The water vapour above each of `levels` under `temps`, by the saturated vapour each
    /// holds, in the units of its ground's, and each level's saturated vapour's share of the
    /// pressure (none without water): no more of it than the driest level below holds, its
    /// water frozen out at its coldest (M65: the cold trap that keeps a stratosphere dry).
    fn vapour_of(&self, temps: &[f64], levels: &[f64]) -> (Vec<f64>, Vec<f64>) {
        let n = levels.len();
        let (mut vapour, mut humid) = (vec![0.0; n], vec![0.0; n]);
        if let Some(water) = volatile("H2O").filter(|_| self.water > 0.0) {
            let mut driest = f64::INFINITY;
            for (k, h) in humid.iter_mut().enumerate() {
                let most = temps[k].min(self.vapour_most);
                driest = driest.min(vapour_pressure(water, most).min(0.5 * levels[k]) / levels[k]);
                *h = driest;
            }
            for k in (0..n - 1).rev() {
                let mean = 0.5 * (humid[k] + humid[k + 1]);
                vapour[k] = vapour[k + 1] + mean * (levels[k] - levels[k + 1]);
            }
        }
        (vapour, humid)
    }

    /// The sulfuric deck in the air over the datum's ground at `t` K (M65): its base, Pa, where
    /// air rising on the adiabat from the datum's ground first saturates in its acid, read
    /// between the levels in the logarithms, as a cloud's base stands at its condensation
    /// level; its top SULFURIC_DECK_SPAN of the base's pressure; none where the rising air never
    /// saturates in it. Placed along the air's own profile instead, the deck and the profile
    /// moved each other, and the solve wandered between a deck high in a warm air and none. The
    /// deck stands in the air, whatever ground stands under it: a summit's column holds what of
    /// it stands over the summit (a deck placed from the summit's ground rose over every height,
    /// and Venus's 11 km summit stood under more greenhouse than its lowlands).
    fn deck_span(&self, t: f64) -> Option<(f64, f64)> {
        let acid = self.deck.2;
        if acid <= 0.0 {
            return None;
        }
        let levels = &self.levels(1.0);
        let temps = &self.adiabat_over(t, levels);
        let (a, b, t0, tc, c) = SULFURIC_ACID_VAPOUR;
        // The acid's partial pressure over its saturated one, in the logarithm.
        let excess = |p: f64, t: f64| -> f64 {
            let t = t.max(1.0);
            let saturated = a - b / t0
                + b * (-1.0 / t + 1.0 / t0 + c / (tc - t0) * (1.0 + (t0 / t).ln() - t0 / t));
            (acid * p / 1e5).ln() - saturated
        };
        let mut last = excess(levels[0], temps[0]);
        if last >= 0.0 {
            return Some((levels[0], levels[0] / SULFURIC_DECK_SPAN));
        }
        for k in 1..levels.len() {
            let now = excess(levels[k], temps[k]);
            if now >= 0.0 {
                let w = last / (last - now);
                let base = (levels[k - 1].ln() + w * (levels[k] / levels[k - 1]).ln()).exp();
                return Some((base, base / SULFURIC_DECK_SPAN));
            }
            last = now;
        }
        None
    }

    /// The deck's depth `whole` above `p` Pa, spread evenly in the pressure's logarithm between
    /// its base and its top (M65): about evenly by height, as Venus's clouds' measured
    /// extinction stands within a factor of a few from their base at 48 km to their top at 70
    /// (Knollenberg and Hunten 1980, Pioneer Venus's particle sizes), where the air's mass
    /// falls some thirtyfold; spread by the air's mass, the deck held nearly all its depth low
    /// and let out its warm upper levels' light.
    fn deck_above(span: Option<(f64, f64)>, whole: f64, p: f64) -> f64 {
        match span {
            Some((base, top)) if base > top => {
                whole * ((p.min(base) / top).ln().max(0.0) / (base / top).ln())
            }
            _ => 0.0,
        }
    }

    /// The column's cold start over ground at `t` K (F293): up the adiabat, never colder than
    /// the bulk gas's frost nor the skin.
    fn cold_profile(&self, t: f64, levels: &[f64]) -> Vec<f64> {
        self.adiabat_over(t, levels)
            .iter()
            .zip(levels)
            .map(|(&a, &p)| a.max(self.frost(p)).max(self.skin))
            .collect()
    }

    /// The air's temperatures at `levels` over ground at `t` K, the first level the ground's
    /// (M65): the warmer of the adiabat from the ground and the radiative-convective profile
    /// over the datum read between the nodes of its ground's temperature, never colder than
    /// the bulk gas's frost; with the water vapour above each level and each level's saturated
    /// vapour's share of the pressure (none without water).
    fn profile(&self, t: f64, levels: &[f64]) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        let thick = self.thickness();
        if thick >= 1.0 {
            let mut temps = self.cold_profile(t, levels);
            temps[0] = t;
            let (vapour, humid) = self.vapour_of(&temps, levels);
            return (temps, vapour, humid);
        }
        let datum = self.levels(1.0);
        let rc = self.rc_at(t);
        let m = datum.len();
        // The profile at `p` Pa over the datum, read in the pressure's logarithm; nothing
        // under the datum's ground.
        let at = |p: f64| -> f64 {
            if p > datum[0] {
                return f64::NEG_INFINITY;
            }
            if p <= datum[m - 1] {
                return rc[m - 1];
            }
            let k = datum.iter().position(|&q| q < p).unwrap_or(m - 1).max(1);
            let w = (p / datum[k - 1]).ln() / (datum[k] / datum[k - 1]).ln();
            rc[k - 1] + w * (rc[k] - rc[k - 1])
        };
        let mut temps: Vec<f64> = self
            .adiabat_over(t, levels)
            .iter()
            .zip(levels)
            .map(|(&a, &p)| a.max(at(p)).max(self.frost(p)))
            .collect();
        if thick > 0.0 {
            for (x, f) in temps.iter_mut().zip(self.cold_profile(t, levels)) {
                *x += thick * (f - *x);
            }
        }
        temps[0] = t;
        let (vapour, humid) = self.vapour_of(&temps, levels);
        (temps, vapour, humid)
    }

    /// How far the column stands from its radiative-convective profile toward F293's, the
    /// adiabat up to its skin, by its infrared depth (M65): none where its ground's light
    /// reaches space through it, all where so little does that its deep air carries its heat
    /// by convection alone, by degrees between RC_THICK's depths. Its depth -ln of the share of
    /// its ground's light its air lets through over ground at its equilibrium temperature, or
    /// its bulk gas's frost point where that is the warmer, on F293's profile, its deck where
    /// its acid condenses; found once for every column alike.
    pub fn thickness(&self) -> f64 {
        use std::sync::Mutex;
        static KEPT: Mutex<Vec<(Vec<u64>, f64)>> = Mutex::new(Vec::new());
        let key = self.node_key();
        if let Some(&(_, w)) = KEPT
            .lock()
            .unwrap_or_else(|e| e.into_inner())
            .iter()
            .find(|(k, _)| *k == key)
        {
            return w;
        }
        let depth = self.infrared_depth();
        let (from, to) = RC_THICK;
        let x = ((depth - from) / (to - from)).clamp(0.0, 1.0);
        let w = x * x * (3.0 - 2.0 * x);
        let mut kept = KEPT.lock().unwrap_or_else(|e| e.into_inner());
        if kept.len() >= RC_COLUMNS_KEPT {
            kept.remove(0);
        }
        kept.push((key, w));
        w
    }

    /// The column's infrared depth (M65): -ln of the share of its ground's light its air lets
    /// through to space over ground at its equilibrium temperature, or its bulk gas's frost
    /// point where that is the warmer, on F293's profile, its deck where its acid condenses.
    pub fn infrared_depth(&self) -> f64 {
        if self.pressure <= 0.0 {
            return 0.0;
        }
        let t = (self.skin / SKIN_SHARE.powf(0.25))
            .max(self.frost(self.pressure))
            .max(1.0);
        let levels = self.levels(1.0);
        let n = levels.len();
        let temps = self.cold_profile(t, &levels);
        let (vapour, humid) = self.vapour_of(&temps, &levels);
        let column = self.vapour_scale(t, vapour[0]);
        let (layers, top) = self.optics(&levels, &temps, (&vapour, &humid), column, 0.0);
        let spectrum = spectrum();
        let (mut through, mut light) = (0.0, 0.0);
        for (i, &nu) in spectrum.nu.iter().enumerate() {
            let depth: f64 = (0..n - 1)
                .map(|k| layers[k].depth(spectrum, i))
                .sum::<f64>()
                + Self::top_depth(spectrum, i, &top);
            let b = planck(nu, t) * spectrum.width[i];
            light += b;
            through += b * (-DIFFUSIVITY * depth).exp();
        }
        -(through / light.max(1e-300)).max(1e-300).ln()
    }

    /// A thick air's deck's depth above `p` Pa (M65): well mixed in its mass, as F293 held it.
    /// Placed where its acid condenses on F293's profile, the deck rose with its ground and its
    /// top stood at about the same temperature whatever its ground's; held where Venus's
    /// measured clouds stand in its air, it let out the same over every ground.
    fn thick_deck(&self, p: f64) -> f64 {
        self.deck.0 * (p / self.pressure).min(1.0)
    }

    /// The radiative-convective profile over the datum's ground at `t` K, read between the
    /// nodes of its temperature RC_NODE_STEP apart in its logarithm (M65).
    fn rc_at(&self, t: f64) -> Vec<f64> {
        let x = t.max(1.0).ln() / RC_NODE_STEP;
        let k0 = x.floor() as i64;
        let w = x - k0 as f64;
        let store = node_store(self.node_key());
        let (a, b) = (self.rc_node(&store, k0), self.rc_node(&store, k0 + 1));
        a.iter().zip(&b).map(|(p, q)| p + w * (q - p)).collect()
    }

    /// Everything that sets the column's radiative-convective profiles, as bits (M65): its
    /// nodes' key in the store all columns share.
    fn node_key(&self) -> Vec<u64> {
        let mut key: Vec<u64> = [
            self.pressure,
            self.gravity,
            self.gas_constant,
            self.heat,
            self.epsilon,
            self.co2_mass,
            self.co2_volume,
            self.ch4_mass,
            self.so2_mass,
            self.skin,
            self.water,
            self.vapour_most,
            self.cloud.0,
            self.cloud.1,
            self.haze_ir,
            self.deck.0,
            self.deck.1,
            self.deck.2,
            self.sunlight.incident,
            self.sunlight.absorbed,
            self.sunlight.haze,
            self.sunlight.methane,
        ]
        .iter()
        .map(|v| v.to_bits())
        .collect();
        match self.bulk {
            Some((v, x)) => key.extend([v as *const Volatile as usize as u64, x.to_bits()]),
            None => key.push(u64::MAX),
        }
        key.push(self.pairs.len() as u64);
        for &(k, x) in &self.pairs {
            key.extend([k as u64, x.to_bits()]);
        }
        key
    }

    /// The radiative-convective profile over the datum's ground at the `k`th node of its
    /// temperature, e^(k RC_NODE_STEP) K, found once and kept (M65): solved from the cold
    /// start at the node nearest the column's equilibrium temperature, whose fourth power its
    /// skin's is SKIN_SHARE of, or its bulk gas's frost point at its ground where that is the
    /// warmer (a ground its air would freeze out on, Venus's under its 92 bar, cost each
    /// column's cold start its 3,000 steps), and from its neighbour's toward it at every
    /// other, so that the nodes follow each other and stand alike however they are asked for.
    ///
    /// The nodes the column's `store` keeps are found under its lock, one chain at a time: a
    /// column read at many temperatures at once solves each node once.
    fn rc_node(&self, store: &Nodes, k: i64) -> Vec<f64> {
        let mut nodes = store.lock().unwrap_or_else(|e| e.into_inner());
        let find = |nodes: &[(i64, Vec<f64>)], j: i64| {
            nodes
                .iter()
                .find(|(key, _)| *key == j)
                .map(|(_, x)| x.clone())
        };
        if let Some(found) = find(&nodes, k) {
            return found;
        }
        let anchor = ((self.skin / SKIN_SHARE.powf(0.25))
            .max(self.frost(self.pressure))
            .max(1.0)
            .ln()
            / RC_NODE_STEP)
            .round() as i64;
        // Back from `k` toward the anchor to the nearest node found, or the anchor's own cold
        // start; then out again to `k`, each node from the one before it.
        let back = if k > anchor { -1 } else { 1 };
        let (mut j, mut start) = (k, None);
        while j != anchor {
            if let Some(x) = find(&nodes, j + back) {
                start = Some(x);
                break;
            }
            j += back;
        }
        loop {
            let t = (j as f64 * RC_NODE_STEP).exp();
            let (_, temps, _, _) =
                self.radiative_convective(t, start.as_deref(), (RC_SETTLED, RC_ROUNDS));
            nodes.push((j, temps.clone()));
            if j == k {
                return temps;
            }
            start = Some(temps);
            j -= back;
        }
    }

    /// The column's levels over ground standing at `ratio` of the surface pressure, Pa, from
    /// the ground up: crowded toward the ground as the square of their rank, where the water
    /// vapour's opacity lies within a fifth of the air's scale height.
    pub fn levels(&self, ratio: f64) -> Vec<f64> {
        let ground = self.pressure * ratio;
        let n = RADIATION_LEVELS;
        (0..n)
            .map(|k| {
                let u = k as f64 / (n - 1) as f64;
                ground * (-RADIATION_DECADES * std::f64::consts::LN_10 * u * u).exp()
            })
            .collect()
    }

    /// The water vapour's column over the datum's ground at `t` K (F235) over what the
    /// profile's saturated levels hold, `datum` of it: the scale its vapour is read by.
    fn vapour_scale(&self, t: f64, datum: f64) -> f64 {
        if datum > 0.0 {
            water_column(t.min(self.vapour_most), self.water) / datum
        } else {
            0.0
        }
    }

    /// The layers between `levels` under `temps`, their water vapour `vapour` and `humid` read
    /// by `column` (M55, M65): each absorber's path through each, and the paths over the top
    /// level, carbon dioxide's band's, the grey depth, methane's band's and the haze's thermal
    /// depth.
    fn optics(
        &self,
        levels: &[f64],
        temps: &[f64],
        (vapour, humid): (&[f64], &[f64]),
        column: f64,
        thick: f64,
    ) -> (Vec<Layer>, [f64; 4]) {
        let n = levels.len();
        // Carbon dioxide's band's path above each level, kg/m² at the reference pressure, and
        // methane's band's alike.
        let band: Vec<f64> = levels
            .iter()
            .map(|&p| p * p / (2.0 * self.gravity * CO2_BAND_PRESSURE))
            .collect();
        // Each layer's path for the molecules' collisions at unit shares, amagat² cm (M55).
        let mut pressed = vec![0.0; n - 1];
        for k in (0..n - 1).rev() {
            let (p0, p1) = (levels[k + 1], levels[k]);
            let mean = 0.5 * (p0 / temps[k + 1] + p1 / temps[k]);
            pressed[k] = (AMAGAT_TEMPERATURE / AMAGAT_PRESSURE).powi(2) * 100.0 * self.gas_constant
                / self.gravity
                * mean
                * (p1 - p0);
        }
        // The vapour's share of saturation, its column against the saturated column's
        // (M55): the saturated vapour's mass above a level eps / g of its profile's.
        let saturation = column * self.gravity / self.epsilon;
        // The clouds, well mixed in the air's mass, a u^b of the vapour u the datum's column
        // holds, kg/m² (M55).
        let held = column * vapour[0];
        let clouds = if held > 0.0 {
            self.cloud.0 * held.powf(self.cloud.1)
        } else {
            0.0
        };
        let span = self.deck_span(temps[0]);
        let grey: Vec<f64> = levels
            .iter()
            .map(|&p| {
                clouds * p / self.pressure
                    + (1.0 - thick) * Self::deck_above(span, self.deck.0, p)
                    + thick * self.thick_deck(p)
            })
            .collect();
        // The haze's thermal depth above each level, by its measured shape (M65).
        let haze: Vec<f64> = levels
            .iter()
            .map(|&p| self.haze_ir * (p / self.pressure).min(1.0).powf(HAZE_SHAPE))
            .collect();
        let (reference, warmth, sigma) = H2O_CONTINUUM_STATE;
        let tables = line_tables();
        let layers = (0..n - 1)
            .map(|k| {
                let (p0, p1) = (levels[k], levels[k + 1]);
                let mass = column * (vapour[k] - vapour[k + 1]);
                let own_pressure = 0.5 * (humid[k] * p0 + humid[k + 1] * p1);
                let lines = column * own_pressure * (p0 - p1) / CO2_BAND_PRESSURE;
                let own = saturation * own_pressure / reference;
                let warm = 0.5 * (temps[k] + temps[k + 1]);
                let cold = (sigma * (warmth - warm)).exp();
                // Each gas's paths, broadened and its own (M65).
                let air = (band[k] - band[k + 1], (p0 - p1) / self.gravity);
                let paths = [
                    (self.co2_mass * air.0, self.co2_mass * air.1),
                    (lines, mass),
                    (self.so2_mass * air.0, self.so2_mass * air.1),
                ];
                Layer {
                    lines: tables
                        .iter()
                        .zip(paths)
                        .enumerate()
                        .filter(|(_, (_, (b, m)))| *b > 0.0 || *m > 0.0)
                        .map(|(g, (t, (b, m)))| {
                            let (j, w) = bracket(&t.temperatures, warm);
                            (g, j, w, b, m)
                        })
                        .collect(),
                    band: self.co2_mass * (band[k] - band[k + 1]),
                    methane: self.ch4_mass * (band[k] - band[k + 1]),
                    grey: grey[k] - grey[k + 1]
                        + CO2_CONTINUUM * self.co2_volume * self.co2_volume * pressed[k],
                    haze: haze[k] - haze[k + 1],
                    water: lines,
                    continuum: mass * own * cold,
                    pairs: self
                        .pairs
                        .iter()
                        .map(|&(q, x)| {
                            let (j, w) = collision_tables().pairs[q].bracket(warm);
                            (q, j, w, x * pressed[k])
                        })
                        .collect(),
                }
            })
            .collect();
        let top = [
            self.co2_mass * band[n - 1],
            grey[n - 1],
            self.ch4_mass * band[n - 1],
            haze[n - 1],
        ];
        (layers, top)
    }

    /// The depth over the top level at the spectrum's `i`th wavenumber, by its paths `top`.
    fn top_depth(spectrum: &Spectrum, i: usize, top: &[f64; 4]) -> f64 {
        top[0] * spectrum.band[i]
            + top[1]
            + top[2] * spectrum.methane[i]
            + top[3] * spectrum.haze[i]
    }

    /// The upward and downward fluxes, W/m², at each of the column's levels over the datum's
    /// ground at `t` K under the air's temperatures `temps` (M65), with what it lets out to
    /// space and each layer's emissivity weighted by its own radiation over the spectrum.
    pub fn fluxes_over(&self, t: f64, temps: &[f64]) -> (Vec<f64>, Vec<f64>, f64, Vec<f64>) {
        let levels = self.levels(1.0);
        let n = levels.len();
        assert_eq!(temps.len(), n);
        let (vapour, humid) = self.vapour_of(temps, &levels);
        let column = self.vapour_scale(t, vapour[0]);
        let (layers, top) =
            self.optics(&levels, temps, (&vapour, &humid), column, self.thickness());
        let spectrum = spectrum();
        let (mut up, mut down, mut olr) = (vec![0.0; n], vec![0.0; n], 0.0);
        let (mut f_up, mut f_dn, mut tr) = (vec![0.0; n], vec![0.0; n], vec![0.0; n - 1]);
        let (mut emit_up, mut emit_dn) = (vec![0.0; n - 1], vec![0.0; n - 1]);
        let (mut emitted, mut black) = (vec![0.0; n - 1], vec![0.0; n - 1]);
        let mut light = vec![0.0; n];
        for (i, &nu) in spectrum.nu.iter().enumerate() {
            let w = spectrum.width[i];
            for (b, &temp) in light.iter_mut().zip(temps) {
                *b = planck(nu, temp);
            }
            for k in 0..n - 1 {
                let x = DIFFUSIVITY * layers[k].depth(spectrum, i);
                let (emits, slope, through) = linear_source(x);
                let (b0, b1) = (light[k], light[k + 1]);
                tr[k] = through;
                let weight = w * 0.5 * (b0 + b1);
                emitted[k] += emits * weight;
                black[k] += weight;
                emit_up[k] = b1 * emits + (b0 - b1) * slope;
                emit_dn[k] = b0 * emits + (b1 - b0) * slope;
            }
            f_up[0] = planck(nu, t);
            for k in 0..n - 1 {
                f_up[k + 1] = f_up[k] * tr[k] + emit_up[k];
            }
            f_dn[n - 1] = 0.0;
            for k in (0..n - 1).rev() {
                f_dn[k] = f_dn[k + 1] * tr[k] + emit_dn[k];
            }
            for k in 0..n {
                up[k] += f_up[k] * w;
                down[k] += f_dn[k] * w;
            }
            // The air above the top level holds back its share of the top's upward stream.
            let above = Self::top_depth(spectrum, i, &top);
            olr += f_up[n - 1] * (-DIFFUSIVITY * above).exp() * w;
        }
        let pi = std::f64::consts::PI;
        (
            up.iter().map(|u| pi * u).collect(),
            down.iter().map(|d| pi * d).collect(),
            pi * olr,
            emitted
                .iter()
                .zip(&black)
                .map(|(e, b)| e / b.max(1e-300))
                .collect(),
        )
    }

    /// The change of the net upward flux at each of the column's levels over the datum's ground
    /// at `t` K under `temps`, W/m², per kelvin of each level's temperature (M65), the air's
    /// absorbers held as they stand: each level's emission into the streams through the layers
    /// about it, their sources linear in their depths, by its Planck function's slope; the
    /// ground's held.
    fn net_slopes(&self, t: f64, temps: &[f64]) -> Vec<Vec<f64>> {
        let levels = self.levels(1.0);
        let n = levels.len();
        let (vapour, humid) = self.vapour_of(temps, &levels);
        let column = self.vapour_scale(t, vapour[0]);
        let (layers, _) = self.optics(&levels, temps, (&vapour, &humid), column, self.thickness());
        let spectrum = spectrum();
        let mut slopes = vec![vec![0.0; n]; n];
        let mut layer = vec![(0.0, 0.0, 0.0); n - 1];
        let (mut up, mut down) = (vec![0.0; n], vec![0.0; n]);
        let pi = std::f64::consts::PI;
        for (i, &nu) in spectrum.nu.iter().enumerate() {
            for (k, l) in layer.iter_mut().enumerate() {
                *l = linear_source(DIFFUSIVITY * layers[k].depth(spectrum, i));
            }
            for m in 1..n {
                let db = pi * spectrum.width[i] * planck_slope(nu, temps[m]);
                if db == 0.0 {
                    continue;
                }
                // Upward: from the top of the layer below the level, and of the layer above.
                up.iter_mut().for_each(|u| *u = 0.0);
                let (emits, slope, _) = layer[m - 1];
                up[m] = emits - slope;
                for j in m..n - 1 {
                    up[j + 1] = up[j] * layer[j].2 + if j == m { layer[m].1 } else { 0.0 };
                }
                // Downward: from the bottom of the layer above the level, and of the layer below.
                down.iter_mut().for_each(|d| *d = 0.0);
                if m + 1 < n {
                    down[m] = layer[m].0 - layer[m].1;
                }
                for j in (1..=m).rev() {
                    down[j - 1] =
                        down[j] * layer[j - 1].2 + if j == m { layer[m - 1].1 } else { 0.0 };
                }
                for j in 0..n {
                    slopes[j][m] += (up[j] - down[j]) * db;
                }
            }
        }
        slopes
    }

    /// The sunlight each layer over the datum's ground at `t` K under `temps` absorbs, W/m², and
    /// what the ground absorbs (M65): what the planet absorbs laid down through the air over
    /// the lit hemisphere's sun heights, each by the light it brings, as the climate's split
    /// reads it, each layer taking what its absorbing depth takes of it:
    /// the haze's by its measured shape, (p / p_s)^HAZE_ABSORBING_SHAPE of its depth above a
    /// level; methane's and the clouds' by the air's mass; water vapour's band by the vapour
    /// above it, each kilogram by its pressure over WATER_SUNLIGHT_PRESSURE (Lacis and Hansen
    /// 1974's scaled amount), as the cycle's split reads them.
    fn sun_layers(&self, t: f64, temps: &[f64]) -> (Vec<f64>, f64) {
        let levels = self.levels(1.0);
        let n = levels.len();
        let light = self.sunlight;
        if light.absorbed <= 0.0 {
            return (vec![0.0; n - 1], 0.0);
        }
        let (vapour, humid) = self.vapour_of(temps, &levels);
        let column = self.vapour_scale(t, vapour[0]);
        let held = column * vapour[0];
        // The vapour above each level, each kilogram by its pressure, in its profile's units.
        let mut scaled = vec![0.0; n];
        for k in (0..n - 1).rev() {
            let mean = 0.5 * (humid[k] * levels[k] + humid[k + 1] * levels[k + 1]);
            scaled[k] = scaled[k + 1] + mean * (levels[k] - levels[k + 1]);
        }
        let clouds = if held > 0.0 {
            CLOUD_SUNLIGHT * self.cloud.0 * held.powf(self.cloud.1)
        } else {
            0.0
        };
        // The haze absorbs above its base, where its extinction from the top reaches
        // HAZE_ABSORBING_BASE, its extinction HAZE_EXTINCTION for each unit of its index by
        // HAZE_SHAPE; of its absorbing depth, (p / p_base)^HAZE_ABSORBING_SHAPE above p.
        let span = self.deck_span(temps[0]);
        let index = light.haze / HAZE_SUNLIGHT;
        let base = if index > 0.0 {
            self.pressure
                * (HAZE_ABSORBING_BASE / (HAZE_EXTINCTION * index))
                    .min(1.0)
                    .powf(1.0 / HAZE_SHAPE)
        } else {
            self.pressure
        };
        let depth = |k: usize, mu: f64| -> f64 {
            let path = magnification(mu);
            let x = (levels[k] / self.pressure).min(1.0);
            let haze = light.haze * (levels[k] / base).min(1.0).powf(HAZE_ABSORBING_SHAPE);
            let upper = span.map(|(base, top)| (base.min(top * SULFURIC_SUN_SPAN), top));
            let grey = haze
                + (light.methane + clouds) * x
                + Self::deck_above(upper, self.deck.1, levels[k]);
            let wet = if held > 0.0 {
                let u = column * scaled[k] / WATER_SUNLIGHT_PRESSURE / 10.0;
                -(1.0 - water_sunlight(u * path)).max(1e-12).ln()
            } else {
                0.0
            };
            grey / mu + wet
        };
        // Over the lit hemisphere's sun heights, RC_SUN_HEIGHTS of them evenly in their
        // cosine, each by the light it brings, its cosine; the haze's own share apart.
        let haze_only = |k: usize, mu: f64| -> f64 {
            light.haze * (levels[k] / base).min(1.0).powf(HAZE_ABSORBING_SHAPE) / mu
        };
        let (mut layers, mut hazed, mut ground, mut weights) =
            (vec![0.0; n - 1], vec![0.0; n - 1], 0.0, 0.0);
        for j in 0..RC_SUN_HEIGHTS {
            let mu = (j as f64 + 0.5) / RC_SUN_HEIGHTS as f64;
            let depths: Vec<f64> = (0..n).map(|k| depth(k, mu)).collect();
            for k in 0..n - 1 {
                let taken = (-depths[k + 1]).exp() - (-depths[k]).exp();
                // The haze's part of the layer's, by its share of the layer's depth.
                let all = depths[k] - depths[k + 1];
                let own = haze_only(k, mu) - haze_only(k + 1, mu);
                let part = if all > 0.0 {
                    (own / all).clamp(0.0, 1.0)
                } else {
                    0.0
                };
                layers[k] += mu * taken * (1.0 - part);
                hazed[k] += mu * taken * part;
            }
            // What the air over the top level takes, its top layer's, the haze's.
            hazed[n - 2] += mu * (1.0 - (-depths[n - 1]).exp());
            ground += mu * (-depths[0]).exp();
            weights += mu;
        }
        // The haze takes its share of the incident light; what the planet reflects comes off
        // what passes under it, so that the column absorbs what the planet does.
        let incident = light.incident.max(light.absorbed) / weights;
        let haze_total: f64 = hazed.iter().sum::<f64>() * incident;
        let rest: f64 = (layers.iter().sum::<f64>() + ground) * incident;
        let keep = if rest > 0.0 {
            ((light.absorbed - haze_total) / rest).clamp(0.0, 1.0)
        } else {
            0.0
        };
        let scale = if haze_total > light.absorbed && haze_total > 0.0 {
            light.absorbed / haze_total
        } else {
            1.0
        };
        (
            layers
                .iter()
                .zip(&hazed)
                .map(|(l, h)| (l * keep + h * scale) * incident)
                .collect(),
            ground * keep * incident,
        )
    }

    /// The radiative-convective profile over the datum's ground at `t` K (M65): the air's
    /// temperatures at the column's levels in radiative equilibrium, each layer's net infrared
    /// heating balancing the sunlight it absorbs (`sun_layers`), save where that would stand
    /// steeper than the adiabat, which convection holds it to from the ground up (Manabe and
    /// Wetherald 1967), never colder than the bulk gas's frost. Each step relaxes every level
    /// by its heating over its radiation's change per kelvin, the layers about it emitting by
    /// their emissivity, halved; the steps are mixed by Anderson's acceleration over the last
    /// RC_HISTORY of them (Anderson 1965; Walker and Ni 2011), the mixing begun afresh where a
    /// step moves the profile more than the last. Marched from `start`, or the cold start,
    /// until no level moves by `settled` K, at most `rounds` steps; returns the levels, the
    /// profile, the steps taken and the last step's largest change, K.
    pub fn radiative_convective(
        &self,
        t: f64,
        start: Option<&[f64]>,
        (settled, rounds): (f64, usize),
    ) -> (Vec<f64>, Vec<f64>, usize, f64) {
        let levels = self.levels(1.0);
        let x = match start {
            Some(s) if s.len() == levels.len() => s.to_vec(),
            _ => self.cold_profile(t, &levels),
        };
        let (x, taken, change) = self.settle(t, &levels, x, (settled, rounds));
        (levels, x, taken, change)
    }

    /// The profile over the datum's ground at `t` K relaxed from `x` through `levels` toward
    /// its radiative-convective equilibrium (M65): `radiative_convective`'s steps; returns the
    /// profile, the steps taken and the last step's largest change, K.
    fn settle(
        &self,
        t: f64,
        levels: &[f64],
        mut x: Vec<f64>,
        (settled, rounds): (f64, usize),
    ) -> (Vec<f64>, usize, f64) {
        let n = levels.len();
        let frost: Vec<f64> = levels.iter().map(|&p| self.frost(p)).collect();
        // Each level's share of the air's mass, its pressure's span between its neighbours'
        // midpoints.
        let mass: Vec<f64> = (0..n)
            .map(|j| {
                let below = if j == 0 {
                    levels[0]
                } else {
                    0.5 * (levels[j - 1] + levels[j])
                };
                let above = if j + 1 == n {
                    0.0
                } else {
                    0.5 * (levels[j] + levels[j + 1])
                };
                below - above
            })
            .collect();
        // The adiabat's ratio from level `j` to the next up, at `x`'s temperature there.
        let ratio = |x: &[f64], j: usize| -> f64 {
            ((levels[j + 1] / levels[j]).ln() * self.lapse(x[j], levels[j])).exp()
        };
        // Convection (Manabe and Strickler 1964): the layer convecting from the ground raised
        // to its adiabat, the ground warming it; above it, each layer standing steeper than
        // the adiabat mixed to one that holds its heat, sum c_p T dp, and grown down and up
        // while it stands so against its neighbours; no level colder than the bulk gas's frost.
        // Whether level `j` stands no warmer than the adiabat from the level below it, within
        // RC_ADIABAT_SLACK: on it, it convects.
        let colder = |x: &[f64], j: usize| -> bool {
            x[j] < x[j - 1] * ratio(x, j - 1) * (1.0 + RC_ADIABAT_SLACK)
        };
        // Returns what holds each level: the ground's convection (-1), the frost (-2), its own
        // radiation (0), or the convecting layer aloft it belongs to (1, 2, ...).
        let adjust = |x: &mut [f64]| -> Vec<i32> {
            let mut kind = vec![0i32; n];
            x[0] = t;
            kind[0] = -1;
            let mut j = 1;
            while j < n && colder(x, j) {
                x[j] = x[j - 1] * ratio(x, j - 1);
                kind[j] = -1;
                j += 1;
            }
            let mut ground = j;
            let mut layer = 0;
            while j < n {
                if !colder(x, j) {
                    j += 1;
                    continue;
                }
                let (mut a, mut b) = (j - 1, j);
                loop {
                    // The adiabat through a..=b holding the layer's heat.
                    let mut shape = vec![1.0; b - a + 1];
                    for k in a..b {
                        shape[k - a + 1] = shape[k - a] * ratio(x, k);
                    }
                    let heat: f64 = (a..=b).map(|k| mass[k] * x[k]).sum();
                    let weight: f64 = (a..=b).map(|k| mass[k] * shape[k - a]).sum();
                    let base = heat / weight.max(1e-300);
                    for k in a..=b {
                        x[k] = base * shape[k - a];
                    }
                    if a > ground && colder(x, a) {
                        a -= 1;
                    } else if b + 1 < n && colder(x, b + 1) {
                        b += 1;
                    } else {
                        break;
                    }
                }
                // Mixed down to the ground's convection and colder than its adiabat: the
                // ground's, raised to its adiabat, the ground warming it.
                if a == ground && colder(x, a) {
                    for k in a..=b {
                        x[k] = x[k - 1] * ratio(x, k - 1);
                        kind[k] = -1;
                    }
                    let mut k = b + 1;
                    while k < n && colder(x, k) {
                        x[k] = x[k - 1] * ratio(x, k - 1);
                        kind[k] = -1;
                        k += 1;
                    }
                    ground = k;
                    j = k;
                    continue;
                }
                layer += 1;
                kind[a..=b].fill(layer);
                j = b + 1;
            }
            for j in 1..n {
                let floor = frost[j].max(1.0);
                if x[j] <= floor {
                    x[j] = floor;
                    kind[j] = -2;
                }
            }
            kind
        };
        // One relaxation step from `x`.
        let relax = |x: &[f64]| -> Vec<f64> {
            let (up, down, _, emissivity) = self.fluxes_over(t, x);
            let (sun, _) = self.sun_layers(t, x);
            let net: Vec<f64> = up.iter().zip(&down).map(|(u, d)| u - d).collect();
            let mut next = x.to_vec();
            for j in 1..n {
                let below = net[j - 1] - net[j] + sun[j - 1];
                let (above, e_above) = if j + 1 < n {
                    (net[j] - net[j + 1] + sun[j], emissivity[j])
                } else {
                    (0.0, 0.0)
                };
                let e = (emissivity[j - 1] + e_above).max(1e-9);
                let stiffness = 4.0 * SIGMA_SB * x[j].powi(3) * e;
                let step = ((below + above) / stiffness).clamp(-RC_CLAMP, RC_CLAMP);
                next[j] = x[j] + 0.5 * step;
            }
            adjust(&mut next);
            next
        };
        // Each level's heating, W/m²: the net flux into the layers about it and the sunlight
        // they absorb, as the relaxation reads it.
        let heating = |x: &[f64]| -> Vec<f64> {
            let (up, down, _, _) = self.fluxes_over(t, x);
            let (sun, _) = self.sun_layers(t, x);
            let net: Vec<f64> = up.iter().zip(&down).map(|(u, d)| u - d).collect();
            (0..n)
                .map(|j| match j {
                    0 => 0.0,
                    _ if j + 1 < n => net[j - 1] - net[j + 1] + sun[j - 1] + sun[j],
                    _ => net[j - 1] - net[j] + sun[j - 1],
                })
                .collect()
        };
        let mut kind = adjust(&mut x);
        let (mut change, mut taken) = (f64::INFINITY, 0);
        // Newton's steps (M65), the air's absorbers held as they stand at each: the unknowns
        // each level its own radiation holds, and each convecting layer aloft by its base, its
        // levels along its adiabat; their equations each level's heating, and each layer's
        // heating summed, none. Where they stall, the relaxation's steps go on from them.
        let (mut best, mut stalls) = (f64::INFINITY, 0);
        while taken < rounds && stalls < RC_NEWTON_STALLS {
            taken += 1;
            let r = heating(&x);
            let slopes = self.net_slopes(t, &x);
            let slope = |j: usize, m: usize| -> f64 {
                if j + 1 < n {
                    slopes[j - 1][m] - slopes[j + 1][m]
                } else {
                    slopes[j - 1][m] - slopes[j][m]
                }
            };
            let mut unknowns: Vec<Vec<(usize, f64)>> = Vec::new();
            let mut j = 1;
            while j < n {
                match kind[j] {
                    0 => {
                        unknowns.push(vec![(j, 1.0)]);
                        j += 1;
                    }
                    // Held by the ground's convection or the frost, and its radiation would warm
                    // it: free of them.
                    -1 | -2 if r[j] > 0.0 => {
                        unknowns.push(vec![(j, 1.0)]);
                        j += 1;
                    }
                    layer if layer > 0 => {
                        let base = x[j];
                        let mut levels = Vec::new();
                        while j < n && kind[j] == layer {
                            levels.push((j, x[j] / base));
                            j += 1;
                        }
                        unknowns.push(levels);
                    }
                    _ => j += 1,
                }
            }
            let m = unknowns.len();
            if m == 0 {
                change = 0.0;
                break;
            }
            let mut a = vec![vec![0.0; m]; m];
            let mut b = vec![0.0; m];
            for (p, rows) in unknowns.iter().enumerate() {
                b[p] = -rows.iter().map(|&(k, _)| r[k]).sum::<f64>();
                for (q, columns) in unknowns.iter().enumerate() {
                    a[p][q] = rows
                        .iter()
                        .map(|&(k, _)| columns.iter().map(|&(l, w)| slope(k, l) * w).sum::<f64>())
                        .sum();
                }
            }
            let Some(step) = solve_dense(a, b) else {
                break;
            };
            // No level stepped further than RC_CLAMP: where the thin air's own radiation hardly
            // holds it, its step is held, as the relaxation's are, and the rest go on.
            let mut next = x.clone();
            for (levels, d) in unknowns.iter().zip(&step) {
                for &(l, w) in levels {
                    next[l] += (w * d).clamp(-RC_CLAMP, RC_CLAMP);
                }
            }
            kind = adjust(&mut next);
            change = next
                .iter()
                .zip(&x)
                .fold(0.0, |c, (u, v)| f64::max(c, (u - v).abs()));
            x = next;
            if change <= settled {
                return (x, taken, change);
            }
            if change < best {
                (best, stalls) = (change, 0);
            } else {
                stalls += 1;
            }
        }
        let (mut dx, mut df): (Vec<Vec<f64>>, Vec<Vec<f64>>) = (Vec::new(), Vec::new());
        let (mut last_x, mut last_f): (Option<Vec<f64>>, Option<Vec<f64>>) = (None, None);
        // The plain step a mixed one was taken in place of.
        let mut fallback: Option<Vec<f64>> = None;
        // The least change so far and the steps since it fell; past RC_STALE of them the
        // mixing, overshooting convection's whole layers, gives way to plain steps.
        let (mut least, mut stale) = (f64::INFINITY, 0);
        while taken < rounds {
            taken += 1;
            let g = relax(&x);
            let f: Vec<f64> = g.iter().zip(&x).map(|(a, b)| a - b).collect();
            let was = change;
            change = f.iter().fold(0.0, |m, v| f64::max(m, v.abs()));
            // Settled: the profile just stepped from is the equilibrium's.
            if change <= settled {
                break;
            }
            if change < least {
                (least, stale) = (change, 0);
            } else {
                stale += 1;
            }
            if stale > RC_STALE {
                x = g;
                continue;
            }
            // A step that moves the profile more than the last did starts the mixing afresh,
            // and a mixed step that did so is undone for the plain step it was taken in place
            // of.
            if change > was {
                dx.clear();
                df.clear();
                (last_x, last_f) = (None, None);
                if let Some(back) = fallback.take() {
                    x = back;
                    change = was;
                    continue;
                }
            }
            if let (Some(lx), Some(lf)) = (&last_x, &last_f) {
                dx.push(x.iter().zip(lx).map(|(a, b)| a - b).collect());
                df.push(f.iter().zip(lf).map(|(a, b)| a - b).collect());
                if dx.len() > RC_HISTORY {
                    dx.remove(0);
                    df.remove(0);
                }
            }
            last_x = Some(x.clone());
            last_f = Some(f.clone());
            // The mixing: gamma minimising |f - dF gamma|, by its normal equations, a little
            // damped; the next profile g - (dX + dF) gamma.
            let m = df.len();
            let mut next = g.clone();
            if m > 0 {
                let mut a = vec![vec![0.0; m]; m];
                let mut b = vec![0.0; m];
                for p in 0..m {
                    for q in 0..m {
                        a[p][q] = df[p].iter().zip(&df[q]).map(|(u, v)| u * v).sum::<f64>();
                    }
                    b[p] = df[p].iter().zip(&f).map(|(u, v)| u * v).sum::<f64>();
                    a[p][p] = a[p][p] * (1.0 + 1e-8) + 1e-12;
                }
                if let Some(gamma) = solve_dense(a, b) {
                    for k in 0..n {
                        let mix: f64 = (0..m).map(|p| gamma[p] * (dx[p][k] + df[p][k])).sum();
                        // No further from the plain step than a plain step goes.
                        next[k] -= mix.clamp(-RC_CLAMP, RC_CLAMP);
                    }
                }
            }
            adjust(&mut next);
            fallback = (m > 0).then(|| g.clone());
            x = next;
        }
        (x, taken, change)
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
        let levels = self.levels(ratio);
        let n = levels.len();
        let (temps, vapour, humid) = self.profile(t, &levels);
        // The water vapour above each level, kg/m²: the datum's column (F235) by the share of
        // its vapour above the level.
        let datum = if ratio == 1.0 || vapour[0] <= 0.0 {
            vapour[0]
        } else {
            let datum_levels: Vec<f64> = levels.iter().map(|&p| p / ratio).collect();
            self.profile(t, &datum_levels).1[0]
        };
        let column = self.vapour_scale(t, datum);
        let (layers, top) =
            self.optics(&levels, &temps, (&vapour, &humid), column, self.thickness());
        let spectrum = spectrum();
        let (mut total, mut through, mut ground_light, mut down) = (0.0, 0.0, 0.0, 0.0);
        let (mut depths, mut light) = (vec![0.0; n], vec![0.0; n]);
        let mut layer = vec![(0.0, 0.0, 0.0); n - 1];
        for (i, &nu) in spectrum.nu.iter().enumerate() {
            depths[n - 1] = Self::top_depth(spectrum, i, &top);
            for k in (0..n - 1).rev() {
                depths[k] = depths[k + 1] + layers[k].depth(spectrum, i);
            }
            for (b, &temp) in light.iter_mut().zip(&temps) {
                *b = planck(nu, temp);
            }
            let (mut out, mut sent) = (0.0, 0.0);
            // Each layer emits up at its top and down at its bottom, its source linear in its
            // depth between its levels (Clough et al. 1992): the mean of its levels' where it is
            // thin, the emitting level's where it is thick; through the layers above it to
            // space, and below it to the ground.
            let mut to_space = (-DIFFUSIVITY * depths[n - 1]).exp();
            for k in (0..n - 1).rev() {
                let x = DIFFUSIVITY * (depths[k] - depths[k + 1]);
                let (emits, slope, through) = linear_source(x);
                let up = light[k + 1] * emits + (light[k] - light[k + 1]) * slope;
                out += up * to_space;
                to_space *= through;
                layer[k] = (emits, slope, through);
            }
            let mut to_ground = 1.0;
            for (k, &(emits, slope, through)) in layer.iter().enumerate() {
                let dn = light[k] * emits + (light[k + 1] - light[k]) * slope;
                sent += dn * to_ground;
                to_ground *= through;
            }
            let above = to_space;
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
        let floor = self.frost(self.pressure);
        balance(|t| self.olr(t), absorbed).max(floor)
    }
}

/// The solution of the small dense system `a` x = `b` by Gaussian elimination with partial
/// pivoting, none where it is singular.
fn solve_dense(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Option<Vec<f64>> {
    let n = b.len();
    for c in 0..n {
        let pivot = (c..n).max_by(|&i, &j| a[i][c].abs().total_cmp(&a[j][c].abs()))?;
        if a[pivot][c].abs() < 1e-300 {
            return None;
        }
        a.swap(c, pivot);
        b.swap(c, pivot);
        for r in c + 1..n {
            let f = a[r][c] / a[c][c];
            let (above, below) = a.split_at_mut(r);
            for (x, y) in below[0][c..].iter_mut().zip(&above[c][c..]) {
                *x -= f * y;
            }
            b[r] -= f * b[c];
        }
    }
    let mut x = vec![0.0; n];
    for r in (0..n).rev() {
        let s: f64 = (r + 1..n).map(|k| a[r][k] * x[k]).sum();
        x[r] = (b[r] - s) / a[r][r];
    }
    x.iter().all(|v| v.is_finite()).then_some(x)
}

/// A layer of slant depth `x`'s emission at its edge per its source there, its emissivity
/// 1 - e^(-x), and per its source's change across it, (1 - e^(-x)) / x - e^(-x): its source
/// linear in its depth (Clough et al. 1992), x / 2 - x² / 3 where it is thin; and its
/// transmission, e^(-x).
fn linear_source(x: f64) -> (f64, f64, f64) {
    let through = (-x).exp();
    if x < 1e-4 {
        return (x - 0.5 * x * x, 0.5 * x - x * x / 3.0, through);
    }
    let emits = 1.0 - through;
    (emits, emits / x - through, through)
}

/// A layer of the column between two levels, as its absorbers' paths through it read it: carbon
/// dioxide's band's, kg/m² at the reference pressure, and the grey depth.
struct Layer {
    band: f64,
    /// Methane's band's path, kg/m² at the reference pressure, and the haze's thermal depth at
    /// HAZE_INFRARED's reference wavenumber (M65).
    methane: f64,
    haze: f64,
    grey: f64,
    /// Water vapour's lines' path, kg/m² at the reference pressure, and its continuum's, kg/m²
    /// at its reference state (M55).
    water: f64,
    continuum: f64,
    /// Each pair's collisions: its table, the table's temperatures about the layer's, the
    /// lower's index and the upper's share, and its path, amagat² cm (M55).
    pairs: Vec<(usize, usize, f64, f64)>,
    /// Each gas's lines (M65): its table, the table's temperatures about the layer's, the
    /// lower's index and the upper's share, and the gas's paths, kg/m², broadened at the
    /// reference pressure and its own mass.
    lines: Vec<(usize, usize, f64, f64, f64)>,
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
        let tables = line_tables();
        let lines: f64 = self
            .lines
            .iter()
            .map(|&(g, j, w, broadened, mass)| tables[g].depth(i, (j, w), broadened, mass))
            .sum();
        lines
            + self.band * spectrum.band[i]
            + self.methane * spectrum.methane[i]
            + self.haze * spectrum.haze[i]
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

/// The frost point of the bulk gas of an air of `composition` at `pressure` Pa, K: none
/// without one (M65).
pub fn frost_point(composition: &BTreeMap<String, f64>, pressure: f64) -> f64 {
    bulk_of(composition).map_or(0.0, |(v, x)| saturation_temperature(v, x * pressure))
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

/// The temperature at which `olr`, rising with it, reaches `absorbed`, K: the warmest such
/// (M65), where an air lets out about what it absorbs high over a span of its ground's warmths
/// (Venus's deck), its balance there its present air's; bracketed over BALANCE_SCAN steps of
/// BALANCE_STEP in the logarithm up from the black body's temperature, the highest crossing
/// among them, or past them up or down from it; then by the secant held within the bracket
/// (Illinois), in the logarithm.
fn balance(olr: impl Fn(f64) -> f64, absorbed: f64) -> f64 {
    let black = (absorbed.max(1e-30) / SIGMA_SB).powf(0.25);
    let f = |x: f64| olr(x.exp()) - absorbed;
    let scan: Vec<(f64, f64)> = (0..=BALANCE_SCAN)
        .map(|k| {
            let x = black.ln() + BALANCE_STEP * k as f64;
            (x, f(x))
        })
        .collect();
    let crossing = scan
        .windows(2)
        .rposition(|w| w[0].1 <= 0.0 && w[1].1 >= 0.0);
    let top = scan[BALANCE_SCAN];
    let ((mut a, mut fa), (mut b, mut fb)) = match crossing {
        _ if top.1 < 0.0 => (top, top),
        Some(k) => (scan[k], scan[k + 1]),
        None => (scan[0], scan[0]),
    };
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
        // Held to no more than at any warmer temperature (M65): past water's runaway plateau,
        // where its vapour's growth stalls the outgoing radiation (Simpson 1927; Nakajima et
        // al. 1992), at the plateau's; where an air lets out about what it absorbs high over a
        // span of its ground's warmths (Venus's deck), at that span's least, so that its balance
        // is the warmest, as the zero-dimensional balance's.
        let mut lowest = f64::INFINITY;
        let mut share: Vec<f64> = streams
            .iter()
            .enumerate()
            .rev()
            .map(|(k, s)| {
                let black = SIGMA_SB * (from + step * k as f64).exp().powi(4);
                lowest = lowest.min(s.0 * black);
                lowest / black
            })
            .collect();
        share.reverse();
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

/// The share of its star's light ozone absorbs on a path of `x` cm at STP (M67): Lacis and
/// Hansen 1974's fits over the Sun's spectrum, the Chappuis band's 0.02118 x / (1 + 0.042 x +
/// 0.000323 x²) and the Hartley and Huggins bands' 1.082 x / (1 + 138.6 x)^0.805 + 0.0658 x /
/// (1 + (103.6 x)³).
pub fn ozone_sunlight(x: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    let visible = 0.02118 * x / (1.0 + 0.042 * x + 0.000323 * x * x);
    let ultraviolet =
        1.082 * x / (1.0 + 138.6 * x).powf(0.805) + 0.0658 * x / (1.0 + (103.6 * x).powi(3));
    visible + ultraviolet
}

/// The ozone column of an air whose oxygen stands at `oxygen` of it, cm at STP (M67): by
/// OZONE_COLUMN, none without oxygen.
pub fn ozone_column(oxygen: f64) -> f64 {
    let level = oxygen.max(0.0) / OXYGEN_PRESENT;
    let (first, last) = (OZONE_COLUMN[0], OZONE_COLUMN[OZONE_COLUMN.len() - 1]);
    let du = if level <= first.0 {
        first.1 * level / first.0
    } else if level >= last.0 {
        last.1
    } else {
        let k = OZONE_COLUMN
            .windows(2)
            .position(|w| level <= w[1].0)
            .unwrap_or(0);
        let (a, b) = (OZONE_COLUMN[k], OZONE_COLUMN[k + 1]);
        let f = (level / a.0).ln() / (b.0 / a.0).ln();
        a.1 + f * (b.1 - a.1)
    };
    du * 1e-3
}

/// The share of its star's light water vapour absorbs on a path of `u` cm of precipitable
/// water (M58): Lacis and Hansen 1974's fit over the Sun's spectrum to its bands,
/// 2.9 u / ((1 + 141.5 u)^0.635 + 5.925 u), saturating towards a half.
pub fn water_sunlight(u: f64) -> f64 {
    if u <= 0.0 {
        return 0.0;
    }
    2.9 * u / ((1.0 + 141.5 * u).powf(0.635) + 5.925 * u)
}

/// The magnification of a vertical path towards a sun at height cosine `mu` (Lacis and Hansen
/// 1974), 35 / sqrt(1224 mu² + 1): the slant path through a curved air.
pub fn magnification(mu: f64) -> f64 {
    35.0 / (1224.0 * mu * mu + 1.0).sqrt()
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

    /// M65 exploration: the radiative-convective column over each body's ground, its sunlight
    /// laid down through its air, against its measured profile by height: Huygens (HASI, PDS),
    /// the US Standard Atmosphere 1976, Viking 1 (Seiff and Kirk 1977, afternoon), VIRA (Seiff
    /// et al. 1985).
    #[test]
    #[ignore]
    fn explore_radiative_convective_profiles() {
        let measured: [(&str, &[(f64, f64)]); 4] = [
            (
                "Titan",
                &[
                    (0.0, 93.65),
                    (1.0, 92.4),
                    (5.0, 88.2),
                    (10.0, 83.9),
                    (20.0, 77.2),
                    (30.0, 72.3),
                    (40.0, 70.5),
                    (44.0, 70.43),
                    (60.0, 75.8),
                    (80.0, 125.8),
                    (100.0, 142.0),
                    (150.0, 168.5),
                    (200.0, 179.7),
                ],
            ),
            (
                "Earth",
                &[
                    (0.0, 288.15),
                    (5.0, 255.7),
                    (11.0, 216.65),
                    (20.0, 216.65),
                    (32.0, 228.65),
                    (47.0, 270.65),
                ],
            ),
            (
                "Mars",
                &[
                    (0.0, 237.3),
                    (10.0, 204.9),
                    (20.0, 193.4),
                    (30.0, 181.6),
                    (40.0, 160.8),
                ],
            ),
            (
                "Venus",
                &[
                    (0.0, 735.3),
                    (10.0, 658.0),
                    (20.0, 579.0),
                    (30.0, 495.0),
                    (40.0, 416.0),
                    (50.0, 348.0),
                    (60.0, 263.0),
                    (70.0, 230.0),
                    (80.0, 197.0),
                ],
            ),
        ];
        for (name, data) in measured {
            let graph = crate::presets::preset(name).unwrap().to_graph().0;
            let eb = crate::climate::EnergyBalance::from_graph(&graph);
            let (column, absorbed) = column_of(&graph, eb.albedo);
            let column = column.with_sunlight(Sunlight {
                incident: graph.f64("orbit.insolation") / 4.0,
                absorbed,
                haze: eb.haze,
                methane: eb.methane,
            });
            // The outgoing radiation over grounds about the measured one.
            let scan: Vec<String> = [0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2]
                .iter()
                .map(|f| {
                    let t = data[0].1 * f;
                    format!("{t:.0}:{:.1}", column.olr(t))
                })
                .collect();
            println!("   OLR over grounds: {}", scan.join(" "));
            let started = std::time::Instant::now();
            let ground = column.surface_temperature(absorbed);
            let spent = started.elapsed().as_secs_f64();
            let levels = column.levels(1.0);
            let n = levels.len();
            let (temps, _, _) = column.profile(ground, &levels);
            let (_, down, olr, _) = column.fluxes_over(ground, &temps);
            let (sun, ground_sun) = column.sun_layers(ground, &temps);
            let mut z = vec![0.0; n];
            for k in 1..n {
                let t_mean = 0.5 * (temps[k] + temps[k - 1]);
                z[k] = z[k - 1]
                    + column.gas_constant * t_mean / column.gravity
                        * (levels[k - 1] / levels[k]).ln();
            }
            let incident = graph.f64("orbit.insolation") / 4.0;
            let band = |lo: f64, hi: f64| -> f64 {
                (0..n - 1)
                    .filter(|&k| z[k] >= lo * 1e3 && z[k] < hi * 1e3)
                    .map(|k| sun[k])
                    .sum::<f64>()
                    / incident
            };
            println!(
                "== {name}: ground {ground:.2} K where measured {:.2}, found in {spent:.2} s; absorbs {absorbed:.2} W/m2, its air {:.2}, its ground {ground_sun:.2}; OLR {olr:.2}; down {:.2}",
                data[0].1,
                sun.iter().sum::<f64>(),
                down[0]
            );
            println!(
                "   sunlight of the incident: above 150 km {:.3}, 80-150 {:.3}, air below 80 {:.3}, ground {:.3} (Titan measured 0.20, 0.20, 0.29, 0.11)",
                band(150.0, 1e9),
                band(80.0, 150.0),
                band(0.0, 80.0),
                ground_sun / incident
            );
            println!(
                "   below 80 km: 0-20 {:.3}, 20-44 {:.3}, 44-80 {:.3}",
                band(0.0, 20.0),
                band(20.0, 44.0),
                band(44.0, 80.0)
            );
            let at = |h: f64| -> f64 {
                let k = z
                    .iter()
                    .position(|&zz| zz >= h * 1e3)
                    .unwrap_or(n - 1)
                    .max(1);
                let w = ((h * 1e3 - z[k - 1]) / (z[k] - z[k - 1]).max(1e-9)).clamp(0.0, 1.0);
                temps[k - 1] + w * (temps[k] - temps[k - 1])
            };
            // The measured profile on the column's levels, by its heights: what it lets out and
            // sends down to its measured ground.
            {
                let measured_at = |h: f64| -> f64 {
                    let k = data
                        .iter()
                        .position(|&(hh, _)| hh >= h)
                        .unwrap_or(data.len() - 1);
                    if k == 0 {
                        return data[0].1;
                    }
                    let ((h0, t0), (h1, t1)) = (data[k - 1], data[k]);
                    if h >= h1 {
                        return t1;
                    }
                    t0 + (h - h0) / (h1 - h0) * (t1 - t0)
                };
                // Heights of the levels by the measured profile itself.
                let mut zm = vec![0.0; n];
                let mut tm = vec![data[0].1; n];
                for k in 1..n {
                    let guess = measured_at(zm[k - 1] / 1e3);
                    zm[k] = zm[k - 1]
                        + column.gas_constant * guess / column.gravity
                            * (levels[k - 1] / levels[k]).ln();
                    tm[k] = measured_at(zm[k] / 1e3);
                }
                let (up_m, down_m, olr_m, _) = column.fluxes_over(data[0].1, &tm);
                {
                    let (sun_m, ground_m) = column.sun_layers(data[0].1, &tm);
                    let edges = [
                        0.0, 2.0, 5.0, 10.0, 20.0, 30.0, 40.0, 60.0, 80.0, 150.0, 1e9,
                    ];
                    let mut sums = vec![(0.0, 0.0); edges.len() - 1];
                    for j in 0..n - 1 {
                        let h = 0.5 * (zm[j] + zm[j + 1]) / 1e3;
                        let r = edges
                            .windows(2)
                            .position(|e| h >= e[0] && h < e[1])
                            .unwrap();
                        let ir = (up_m[j] - down_m[j]) - (up_m[j + 1] - down_m[j + 1]);
                        sums[r].0 += ir;
                        sums[r].1 += sun_m[j];
                    }
                    let line: String = edges
                        .windows(2)
                        .zip(&sums)
                        .map(|(e, (ir, sun))| {
                            format!(" {}-{} km: ir {ir:+.3} sun {sun:.3};", e[0], e[1])
                        })
                        .collect();
                    let ground_net = SIGMA_SB * data[0].1.powi(4) - down_m[0];
                    println!(
                        "   measured profile's heating, W/m2:{line} ground: sun {ground_m:.3}, ir out {ground_net:.3}"
                    );
                }
                println!(
                    "   the measured profile over its {:.2} K ground lets out {olr_m:.3} W/m2 and sends down {:.3} (absorbed {absorbed:.3})",
                    data[0].1, down_m[0]
                );
            }
            for &(h, t) in data {
                println!("  {h:6.1} km  measured {t:7.2}  column {:7.2}", at(h));
            }
        }
    }

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
            + crate::params::derive::physics::air_vapour_water(
                &composition,
                g.f64("atmosphere.surface_pressure"),
                g.f64("body.gravity"),
            );
        let absorbed = g.f64("orbit.insolation") / 4.0 * (1.0 - albedo);
        let skin = (absorbed / SIGMA_SB).powf(0.25) * 0.5f64.powf(0.25);
        let (pressure, gravity) = (g.f64("atmosphere.surface_pressure"), g.f64("body.gravity"));
        let column =
            Column::new(&composition, pressure, gravity, water, skin).with_sunlight(Sunlight::of(
                &composition,
                (pressure, gravity),
                g.f64("orbit.insolation") / 4.0,
                albedo,
                g.f64("star.xuv_relative"),
            ));
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

    /// M67: ozone takes Lacis and Hansen (1974)'s share of the Sun's light, 0.0221 of it
    /// through 300 DU from overhead, the more the longer its path; its column by the air's
    /// oxygen as WACCM6 gives it (Cooke et al. 2022), 279 DU at the present level, 18 at a
    /// thousandth of it, none without oxygen, the more the more oxygen; Earth's air taking some
    /// 10 W/m² of its light by it over the lit hemisphere (the notes' 10 to 11.5).
    #[test]
    fn ozone_takes_its_share_of_the_light() {
        assert!(
            (ozone_sunlight(0.3) - 0.0221).abs() < 5e-4,
            "{}",
            ozone_sunlight(0.3)
        );
        assert_eq!(ozone_sunlight(0.0), 0.0);
        let mut last = 0.0;
        for x in [1e-4, 1e-3, 0.01, 0.1, 0.3, 1.0, 3.0] {
            let a = ozone_sunlight(x);
            assert!(a > last && a < 0.2, "{x}: {a}");
            last = a;
        }
        assert!((ozone_column(OXYGEN_PRESENT) - 0.279).abs() < 1e-9);
        assert!((ozone_column(0.001 * OXYGEN_PRESENT) - 0.018).abs() < 1e-9);
        assert_eq!(ozone_column(0.0), 0.0);
        let mut last = -1.0;
        for level in [0.0, 1e-4, 1e-3, 3e-3, 0.01, 0.03, 0.1, 0.5, 1.0, 2.0] {
            let c = ozone_column(level * OXYGEN_PRESENT);
            assert!(c >= last, "{level}: {c}");
            last = c;
        }
        // Over the body, its lit hemisphere's area even in the sun's height cosine: S0 / 2
        // times the integral of mu A(mu).
        let n = 2000;
        let taken = (0..n)
            .map(|k| {
                let mu = (k as f64 + 0.5) / n as f64;
                mu * ozone_sunlight(magnification(mu) * ozone_column(OXYGEN_PRESENT))
            })
            .sum::<f64>()
            / n as f64
            * 1361.0
            / 2.0;
        println!("Earth's ozone takes {taken:.1} W/m²");
        assert!((8.5..12.0).contains(&taken), "{taken}");
    }

    /// F293: carbon dioxide's band widens with the logarithm of its column: over a 288 K
    /// ground under an Earth-like air, each doubling of its carbon dioxide from Earth's takes
    /// within a fifth of the same from the outgoing radiation; the grey law it replaced grew as
    /// the column to the 0.86. The first within a quarter of the band's own law, 2 l ln 2 pi
    /// (B(nu0, T_s) - B(nu0, T_strat)), its absorption falling by e every l from its centre nu0
    /// (Jeevanjee et al. 2021, some 5 W/m² over the Earth's own stratosphere): the column's
    /// stratosphere, unwarmed by ozone, its own (M65).
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
        let composition: BTreeMap<String, f64> = [
            ("N2".to_string(), 1.0 - 4.2e-4),
            ("CO2".to_string(), 4.2e-4),
        ]
        .into_iter()
        .collect();
        let column = Column::new(&composition, 101_325.0, 9.81, 0.0, 214.0);
        let (_, profile, _, _) = column.radiative_convective(288.0, None, (RC_SETTLED, RC_ROUNDS));
        let stratosphere = profile[profile.len() - 1];
        let (_, centre, l) = CO2_BAND;
        let law = 2.0
            * l
            * std::f64::consts::LN_2
            * std::f64::consts::PI
            * (planck(centre, 288.0) - planck(centre, stratosphere));
        println!("the band's law over its {stratosphere:.1} K stratosphere: {law:.2} W/m²");
        assert!(
            (doublings[0] / law - 1.0).abs() < 0.25,
            "{} {law}",
            doublings[0]
        );
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
        println!(
            "Titan: haze index {index:.5}, its depth in the sunlight {tau:.5}: HAZE_SUNLIGHT {:.5}; it absorbs {high:.3} of {light:.3} W/m2, its stratosphere's own; the ground {:.2} K (93.65)",
            tau / index,
            base.surface_temperature(absorbed)
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
        // M65: the deck's depth is stated, Venus's ground its column's.
        let (base, absorbed) = column_of(&venus, a_v);
        let ground = base
            .clone()
            .with_sunlight(Sunlight {
                incident: venus.f64("orbit.insolation") / 4.0,
                absorbed,
                ..Sunlight::default()
            })
            .surface_temperature(absorbed);
        println!(
            "Venus: albedo {a_v:.4}, deck {:.2} deep in the infrared: its ground {ground:.1} K (735.3)",
            base.deck.0
        );
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
    /// the ground's light through to space (3.8 on the Huygens probe's profile), where the grey
    /// column's window was 0.30, and sends down some 3.98 W/m², its ground losing 0.379 net
    /// (McKay et al. 1989): its column's sends down 4.30 and leaves 0.06 (M65), its troposphere
    /// 4 to 7 K warm over Huygens's (4.20 under Huygens's own profile, its haze's thermal depth
    /// with it), a known limitation until methane's moist convection and its sunlight's bands. Its haze sends back to space between
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
        };
        let titan = || {
            let g = graph("Titan");
            let (column, _) = column_of(&g, g.f64("atmosphere.bond_albedo"));
            let s = column.streams(93.7);
            let net = SIGMA_SB * 93.7f64.powi(4) - s.down;
            assert!((0.02..0.06).contains(&s.window), "window {:.4}", s.window);
            assert!((3.7..4.4).contains(&s.down), "down {:.3}", s.down);
            assert!((0.0..0.6).contains(&net), "net {net:.3}");
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
    /// Venusian, each from 150 to 1,000 K, over a ground its air may stand on, at or over its
    /// bulk gas's frost point, never let out more over a colder ground, to rounding, nor under
    /// a thicker air; under that the air would freeze out on the ground, its profile held at
    /// its frost, and a thicker air's frost point holds its column the warmer (the tests'
    /// tiers).
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
                    if temperatures[k - 1] < floors[i] {
                        continue;
                    }
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

    /// M65: the net flux's slopes Newton's steps read meet the net flux's own differences: over
    /// a thin dry air of carbon dioxide, whose absorbers change least with the temperature they
    /// are held at, each level's own and its neighbours' within a fifth wherever they matter.
    #[test]
    fn the_net_flux_s_slopes_meet_its_differences() {
        let composition: BTreeMap<String, f64> =
            [("CO2".to_string(), 0.95), ("N2".to_string(), 0.05)]
                .into_iter()
                .collect();
        let column = Column::new(&composition, 1e4, 3.71, 0.0, 180.0);
        let (_, x, _, _) = column.radiative_convective(240.0, None, (RC_SETTLED, RC_ROUNDS));
        let net = |x: &[f64]| -> Vec<f64> {
            let (up, down, _, _) = column.fluxes_over(240.0, x);
            up.iter().zip(&down).map(|(u, d)| u - d).collect()
        };
        let slopes = column.net_slopes(240.0, &x);
        let base = net(&x);
        let largest = slopes.iter().flatten().fold(0.0f64, |m, v| m.max(v.abs()));
        for m in [2, 8, 16, 24, 32, 40] {
            let mut y = x.clone();
            y[m] += 0.01;
            let moved = net(&y);
            for j in m.saturating_sub(1)..=(m + 1).min(x.len() - 1) {
                let difference = (moved[j] - base[j]) / 0.01;
                if difference.abs() < 1e-3 * largest {
                    continue;
                }
                assert!(
                    (slopes[j][m] / difference - 1.0).abs() < 0.2,
                    "level {m}, net at {j}: {} where {difference}",
                    slopes[j][m]
                );
            }
        }
    }

    /// M65: columns alike share their radiative-convective nodes, found once; a column under
    /// other sunlight, or with other clouds, keeps its own.
    #[test]
    fn columns_alike_share_their_nodes_and_unlike_ones_do_not() {
        let composition: BTreeMap<String, f64> =
            [("N2".to_string(), 0.95), ("CH4".to_string(), 0.05)]
                .into_iter()
                .collect();
        let make = |absorbed: f64| {
            Column::new(&composition, 1.5e5, 1.35, 0.0, 70.0).with_sunlight(Sunlight {
                incident: 3.7,
                absorbed,
                haze: 0.5,
                methane: 0.1,
            })
        };
        let (a, b, other) = (make(2.7), make(2.7), make(2.6));
        let mut cloudy = a.clone();
        cloudy.cloud = (0.5, 0.5);
        let store = node_store(a.node_key());
        assert!(std::sync::Arc::ptr_eq(&store, &node_store(b.node_key())));
        assert!(!std::sync::Arc::ptr_eq(
            &store,
            &node_store(other.node_key())
        ));
        assert!(!std::sync::Arc::ptr_eq(
            &store,
            &node_store(cloudy.node_key())
        ));
        let found = a.olr(90.0);
        let kept = store.lock().map_or(0, |n| n.len());
        assert!(kept > 0);
        assert_eq!(b.olr(90.0), found);
        assert_eq!(store.lock().map_or(0, |n| n.len()), kept);
    }
}
