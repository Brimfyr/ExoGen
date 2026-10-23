// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! The mantle's convection under a rocky lid (F173).
//!
//! Beneath a stagnant or episodic lid the mantle convects in the layer between the stiff
//! lid's base and the core, the cold sublayer under the lid dripping into its downwellings
//! and its upwellings rising between them (Solomatov and Moresi 2000). The pattern is
//! chaotic, so it is drawn from the body's seed as the craters' places are; its scale, its
//! strength and how fast it changes are the mantle's own (F121). Its cells are as wide as
//! CONVECTION_CELL_WAVELENGTH times the layer's depth: many small ones in a thin shell such as
//! Mercury's, a few large ones in a thick one, as King (2008) found Mercury's thin mantle
//! convecting in many narrow rolls where Venus's, the Earth's and Mars's rise in a few
//! plumes. The ground its flow holds up, raised over the upwellings and drawn down over the
//! downwellings, has a root mean square of CONVECTIVE_TOPOGRAPHY times α ΔT_rh δ_rh, the
//! sublayer's thermal contraction over its thickness (Guimond, Rudge and Shorttle 2022's
//! simulations), which an elastic lid holds flatter at the wavelengths it can bend over. The
//! pattern renews itself in the time the flux takes to carry off the heat its sublayer holds,
//! and at each overturn of an episodic lid; between renewals one state turns into the next.

use std::f64::consts::{FRAC_PI_2, TAU};
use std::sync::OnceLock;

use crate::constants::*;
use crate::geom::Turn;
use crate::mantle::{MantleHistory, stagnant_flux};
use crate::params::Reads;
use crate::rng::Rng;
use crate::spharm::{Harmonics, LatLon};
use crate::thermal::ThermalHistory;

pub const VERSION: u64 = 1;

/// The mantle's convection under its lid at a time.
#[derive(Clone, Copy, Debug, Default)]
pub struct Cells {
    /// The temperature drop across the sublayer under the stiff lid, K, and its
    /// thickness, m.
    pub drop: f64,
    pub sublayer: f64,
    /// The convecting layer's depth, m, and the radius of its middle.
    pub depth: f64,
    pub middle: f64,
    /// The degree of its cells.
    pub degree: f64,
    /// The root mean square of the ground its flow holds up, m.
    pub support: f64,
    /// The time its pattern takes to renew itself, Gyr.
    pub renewal: f64,
}

/// The convection under the lid at a time (F173): the sublayer from the mantle's
/// temperature and its stagnant-lid flux, the layer from the stiff lid's base to the core.
pub fn cells(mantle: &MantleHistory, t: f64) -> Cells {
    let temp = mantle.temperature_at(t);
    if temp <= 0.0 {
        return Cells::default();
    }
    let stiff = mantle.stiff_lid_at(t);
    let flux = stagnant_flux(temp, mantle.gravity).max(1e-9);
    let drop = LID_BASE_SCALES * GAS_CONSTANT * temp * temp / MANTLE_ACTIVATION;
    let sublayer = K_ROCK * drop / flux;
    let depth = (mantle.depth - stiff).max(1.0);
    let middle = (mantle.radius - stiff - 0.5 * depth).max(1.0);
    let degree = (TAU * middle / (CONVECTION_CELL_WAVELENGTH * depth) - 0.5).max(1.0);
    // The layer convects once it is as deep as several of its sublayers.
    let (lo, hi) = CONVECTION_ONSET;
    let x = ((depth / sublayer - lo) / (hi - lo)).clamp(0.0, 1.0);
    let runs = x * x * (3.0 - 2.0 * x);
    Cells {
        drop,
        sublayer,
        depth,
        middle,
        degree,
        support: CONVECTIVE_TOPOGRAPHY * MANTLE_EXPANSIVITY * drop * sublayer * runs,
        renewal: DENSITY_ROCK * MANTLE_HEAT_CAPACITY * drop * sublayer / flux / SECONDS_PER_GYR,
    }
}

/// One state of the pattern: its coefficients, and its values on a grid once asked for.
#[derive(Clone, Debug)]
struct Draw {
    time: f64,
    harmonics: Harmonics,
    grid: OnceLock<LatLon>,
}

impl Draw {
    fn grid(&self) -> &LatLon {
        self.grid.get_or_init(|| {
            let nlat = (4 * self.harmonics.lmax).clamp(32, 4 * CONVECTION_MAX_DEGREE);
            self.harmonics.grid(nlat, 2 * nlat)
        })
    }
}

/// A state of the pattern (F173): its power spread about the cells' degree, of unit root
/// mean square, and the wavelengths an elastic lid this thick bends over held flatter, a
/// wavelength's deflection under the flow's push 1 / (1 + D k⁴ / ρ g) of the mantle's own.
fn state(rng: &Rng, draw: u64, cells: &Cells, elastic: f64, mantle: &MantleHistory) -> Harmonics {
    let width = CONVECTION_SPECTRUM_WIDTH;
    let degree = cells.degree;
    let lmax = ((degree * (3.0 * width).exp()).ceil() as usize).clamp(2, CONVECTION_MAX_DEGREE);
    let mut h = Harmonics::random_spectrum(rng, draw, lmax, |l| {
        let x = (l as f64 / degree).ln() / width;
        (-0.5 * x * x).exp()
    });
    let rms = h.rms();
    if rms > 0.0 {
        h.scale_degrees(|_| 1.0 / rms);
    }
    let rigidity = ROCK_YOUNG_MODULUS * elastic.max(0.0).powi(3)
        / (12.0 * (1.0 - ROCK_POISSON * ROCK_POISSON));
    let radius = mantle.radius.max(1.0);
    let weight = DENSITY_ROCK * mantle.gravity.max(1e-9);
    h.scale_degrees(|l| {
        let k2 = (l * (l + 1)) as f64 / (radius * radius);
        1.0 / (1.0 + rigidity * k2 * k2 / weight)
    });
    h
}

/// A rocky mantle's convection through the body's life (F173), in the mantle's frame.
#[derive(Clone, Debug)]
pub struct Convection {
    draws: Vec<Draw>,
    /// The mantle's steps, Gyr, and the root mean square of the ground its flow holds up
    /// then, m.
    times: Vec<f64>,
    support: Vec<f64>,
}

impl Convection {
    /// A body's convection, drawn from its seed; none where it keeps no history of its
    /// mantle.
    pub fn from_graph<P: Reads>(graph: &P, thermal: &ThermalHistory) -> Option<Self> {
        let mantle = thermal.mantle.as_ref()?;
        let seed = graph.num("seed.tectonics") as u64 ^ 0xC0_7EC7_1015;
        Some(Self::new(mantle, thermal, seed))
    }

    pub fn new(mantle: &MantleHistory, thermal: &ThermalHistory, seed: u64) -> Self {
        let rng = Rng::new(seed);
        let age = mantle.times.last().copied().unwrap_or(0.0);
        let times = mantle.times.clone();
        let support = times.iter().map(|&t| cells(mantle, t).support).collect();
        let mut draws = Vec::new();
        let mut t = 0.0;
        loop {
            let c = cells(mantle, t);
            let elastic = thermal.at(t.min(age)).elastic;
            draws.push(Draw {
                time: t,
                harmonics: state(&rng, draws.len() as u64, &c, elastic, mantle),
                grid: OnceLock::new(),
            });
            if t > age || draws.len() >= CONVECTION_MAX_DRAWS {
                break;
            }
            let mut next = t + c.renewal.max(1e-4);
            // An overturn renews the pattern.
            if let Some(&o) = mantle.overturns.iter().find(|&&o| o > t + 1e-9 && o < next) {
                next = o;
            }
            t = next;
        }
        Self {
            draws,
            times,
            support,
        }
    }

    /// The number of states the pattern passes through.
    pub fn renewals(&self) -> usize {
        self.draws.len()
    }

    /// The root mean square of the ground the flow holds up at a time, m.
    pub fn support_rms(&self, t: f64) -> f64 {
        let n = self.times.len();
        if n == 0 {
            return 0.0;
        }
        let k = self.times.partition_point(|&x| x < t);
        if k == 0 {
            return self.support[0];
        }
        if k >= n {
            return self.support[n - 1];
        }
        let (t0, t1) = (self.times[k - 1], self.times[k]);
        let u = ((t - t0) / (t1 - t0).max(1e-12)).clamp(0.0, 1.0);
        self.support[k - 1] + u * (self.support[k] - self.support[k - 1])
    }

    /// The pattern at a time and a direction of the mantle's frame, of unit root mean square
    /// but where the lid holds it flatter: the states about the time blended along a quarter
    /// circle, so that it keeps its strength as one turns into the next.
    pub fn pattern(&self, t: f64, d: [f64; 3]) -> f64 {
        let j = self
            .draws
            .partition_point(|w| w.time <= t)
            .saturating_sub(1);
        let a = &self.draws[j];
        match self.draws.get(j + 1) {
            Some(b) => {
                let u = ((t - a.time) / (b.time - a.time).max(1e-12)).clamp(0.0, 1.0);
                let (s, c) = (FRAC_PI_2 * u).sin_cos();
                c * a.grid().sample(d) + s * b.grid().sample(d)
            }
            None => a.grid().sample(d),
        }
    }

    /// The ground the flow holds up at a time and a direction of the mantle's frame, m.
    pub fn support(&self, t: f64, d: [f64; 3]) -> f64 {
        self.support_rms(t) * self.pattern(t, d)
    }

    /// The pattern's highs at a time, strongest first: where the mantle wells up.
    pub fn upwellings(&self, t: f64) -> Vec<([f64; 3], f64)> {
        let (nlat, nlon) = (96usize, 192usize);
        let mut grid = LatLon {
            nlat,
            nlon,
            values: vec![0.0; nlat * nlon],
        };
        for i in 0..nlat {
            for j in 0..nlon {
                grid.values[i * nlon + j] = self.pattern(t, grid.direction(i, j)) as f32;
            }
        }
        let mut highs = Vec::new();
        for i in 0..nlat {
            for j in 0..nlon {
                let v = grid.at(i, j);
                let mut high = true;
                for di in -1i64..=1 {
                    let ii = i as i64 + di;
                    if ii < 0 || ii >= nlat as i64 {
                        continue;
                    }
                    for dj in -1i64..=1 {
                        if di == 0 && dj == 0 {
                            continue;
                        }
                        let jj = (j as i64 + dj).rem_euclid(nlon as i64) as usize;
                        if grid.at(ii as usize, jj) >= v {
                            high = false;
                        }
                    }
                }
                if high {
                    highs.push((grid.direction(i, j), v));
                }
            }
        }
        highs.sort_by(|a, b| b.1.total_cmp(&a.1));
        highs
    }
}

/// The convection under the map: its pattern turned with the body as its rises turned it
/// (F104), the mantle turning with the lid.
#[derive(Clone, Debug)]
pub struct Support {
    pub convection: Convection,
    pub turn: Turn,
}

impl Support {
    /// The ground the flow holds up at a time and a direction of the map, m.
    pub fn at(&self, t: f64, d: [f64; 3]) -> f64 {
        self.convection.support(t, self.turn.undo(d))
    }

    /// The pattern at a time and a direction of the map, of unit root mean square: above
    /// nothing over an upwelling, below it over a downwelling.
    pub fn pattern(&self, t: f64, d: [f64; 3]) -> f64 {
        self.convection.pattern(t, self.turn.undo(d))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::presets::preset;

    fn convection(name: &str) -> (Convection, MantleHistory, f64) {
        let (g, _) = preset(name).unwrap().to_graph();
        let thermal = ThermalHistory::from_graph(&g);
        let mantle = thermal.mantle.clone().unwrap();
        (
            Convection::from_graph(&g, &thermal).unwrap(),
            mantle,
            g.f64("body.age"),
        )
    }

    /// F173: a thick mantle convects in a few large cells and a thin one in many small ones,
    /// Venus's of degree 3 to 6 and Mercury's past 8; the ground the flow holds up some
    /// hundreds of metres at its root mean square, so that its highs and lows reach a
    /// kilometre or so; and the pattern renews itself in tens to hundreds of millions of
    /// years.
    #[test]
    fn the_cells_follow_the_mantle_s_depth() {
        let mut lines = Vec::new();
        for (name, lo, hi) in [
            ("Venus", 3.0, 6.0),
            ("Mars", 1.5, 4.5),
            ("Mercury", 8.0, 25.0),
            ("Moon", 1.0, 2.5),
        ] {
            let (_, mantle, age) = convection(name);
            let c = cells(&mantle, age);
            lines.push(format!(
                "{name}: degree {:.1}, layer {:.0} km, sublayer {:.0} km, support {:.0} m, renewal {:.0} Myr",
                c.degree,
                c.depth / 1e3,
                c.sublayer / 1e3,
                c.support,
                c.renewal * 1e3
            ));
            assert!(
                c.degree >= lo && c.degree <= hi,
                "{}",
                lines.last().unwrap()
            );
            assert!(
                c.support > 100.0 && c.support < 1_500.0,
                "{}",
                lines.last().unwrap()
            );
            assert!(
                c.renewal > 0.005 && c.renewal < 0.6,
                "{}",
                lines.last().unwrap()
            );
        }
        println!("{}", lines.join("\n"));
    }

    /// F173: the pattern keeps its strength as one state turns into the next, changes
    /// smoothly within a renewal, and is held flatter on Mercury, whose cells are as narrow
    /// as its thick lid bends over, than on Venus.
    #[test]
    fn the_pattern_keeps_its_strength_as_it_changes() {
        let rng = Rng::new(5);
        let rms = |c: &Convection, t: f64| {
            let n = 4000u64;
            let sq: f64 = (0..n)
                .map(|i| c.pattern(t, rng.unit_direction(i, 0)).powi(2))
                .sum();
            (sq / n as f64).sqrt()
        };
        let (venus, _, age) = convection("Venus");
        assert!(venus.renewals() > 20, "{}", venus.renewals());
        for t in [0.5, 1.37, 2.9, age] {
            let r = rms(&venus, t);
            assert!(r > 0.6 && r < 1.3, "{t}: {r}");
        }
        let d = rng.unit_direction(99, 0);
        let step = (venus.pattern(2.0, d) - venus.pattern(2.0005, d)).abs();
        assert!(step < 0.1, "{step}");
        let (mercury, _, age) = convection("Mercury");
        let r = rms(&mercury, age);
        assert!(r < rms(&venus, age) && r > 0.05, "{r}");
    }

    /// Diagnostic: each rocky preset's convection now, the ground its flow holds up as the
    /// lid lets it stand; run with `--ignored --nocapture`.
    #[test]
    #[ignore]
    fn print_convection() {
        let rng = Rng::new(17);
        for name in ["Mercury", "Venus", "Moon", "Mars"] {
            let (c, mantle, age) = convection(name);
            let overturn = mantle.last_overturn().map(|o| o + 0.01);
            for t in [0.5, 1.0, 2.0, 3.0]
                .into_iter()
                .chain(overturn)
                .chain([age])
            {
                let cells = cells(&mantle, t);
                let values: Vec<f64> = (0..20_000u64)
                    .map(|i| c.support(t, rng.unit_direction(i, 0)))
                    .collect();
                let pattern: Vec<f64> = (0..20_000u64)
                    .map(|i| c.pattern(t, rng.unit_direction(i, 0)))
                    .collect();
                let rms = (values.iter().map(|v| v * v).sum::<f64>() / values.len() as f64).sqrt();
                let lo = values.iter().cloned().fold(f64::MAX, f64::min);
                let hi = values.iter().cloned().fold(f64::MIN, f64::max);
                let plo = pattern.iter().cloned().fold(f64::MAX, f64::min);
                let phi = pattern.iter().cloned().fold(f64::MIN, f64::max);
                println!(
                    "{name} t {t:.3}: degree {:.1}, depth {:.0} km, sublayer {:.1} km, renewal {:.1} Myr, {} renewals, unbent rms {:.0} m, standing rms {rms:.0} m, {lo:.0} to {hi:.0} m, pattern {plo:.2} to {phi:.2}",
                    cells.degree,
                    cells.depth / 1e3,
                    cells.sublayer / 1e3,
                    cells.renewal * 1e3,
                    c.renewals(),
                    cells.support
                );
            }
        }
    }

    /// F173: the highs are where the pattern stands highest, strongest first.
    #[test]
    fn the_upwellings_are_the_pattern_s_highs() {
        let (venus, _, _) = convection("Venus");
        let highs = venus.upwellings(1.5);
        assert!(highs.len() >= 3, "{}", highs.len());
        for w in highs.windows(2) {
            assert!(w[0].1 >= w[1].1);
        }
        let top = highs[0];
        assert!((venus.pattern(1.5, top.0) - top.1).abs() < 1e-3);
        let rng = Rng::new(8);
        let higher = (0..2000u64)
            .filter(|&i| venus.pattern(1.5, rng.unit_direction(i, 0)) > top.1 + 1e-3)
            .count();
        assert_eq!(higher, 0);
    }
}
