// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! The sea ice's drift (F276).
//!
//! The wind's stress on the ice, rho_a C_a |V| V, against the still water's drag on it,
//! rho_w C_w |u| u (McPhee 1980), and the spin's turn of its mass, m f k × u: its free drift,
//! some 2 percent of the wind and turned to its right in the north, the more so the thicker
//! the ice and the gentler the wind (Nansen's 20 to 40°). Where the wind presses the ice
//! together it resists with its strength, P* h e^{-C (1 - A)} (Hibler 1979), as a fluid that
//! bears a pressure up to its strength and none drawn apart (Flato and Hibler 1992): pressed
//! beyond it, it gives and piles, so that ice blown against a coast thickens toward it until
//! the strength its pile gathers holds the wind; drawn apart, it opens freely. The ice its
//! drift carries moves between the places' seas upwind, their coasts closed to it.

use crate::constants::*;
use std::f64::consts::{PI, TAU};

/// A boundary two places' seas share (F276): the places, west and east of it or south and
/// north, whether it faces north, the length of it the seas share, m, the distance between
/// the places across it, m, and the spin's Coriolis parameter there, 1/s.
#[derive(Clone, Copy, Debug)]
struct Face {
    from: usize,
    to: usize,
    north: bool,
    length: f64,
    distance: f64,
    coriolis: f64,
}

/// The sea ice's drift on the cycle's layout (F276): the boundaries between its places' seas,
/// each place's sea's area, m², the boundaries about each place and whether it stands west or
/// south of each, the sea's water and ice densities, kg/m³, the pressure the ice bore at each
/// place through the last step, N/m, from which the next step's is found, and the drift so
/// far: the ice's mass times its speed across the boundaries and its mass there, the time its
/// steps took, s, their count, and their sub-steps' and sweeps'.
pub struct Drift {
    faces: Vec<Face>,
    area: Vec<f64>,
    about: Vec<Vec<(usize, bool)>>,
    densities: (f64, f64),
    pressure: Vec<f64>,
    carried: f64,
    held: f64,
    spent: f64,
    steps: usize,
    substeps: usize,
    sweeps: usize,
}

impl Drift {
    /// The drift on a layout of `n_lat` latitudes and `n_lon` longitudes whose places' seas
    /// cover the shares `sea` of them, the places the shares `weight` of the sphere, on a body
    /// of `radius` m spinning at `omega` rad/s, its sea's water and ice `densities`, kg/m³.
    pub fn new(
        (sea, weight): (&[f64], &[f64]),
        (n_lat, n_lon): (usize, usize),
        (omega, radius): (f64, f64),
        densities: (f64, f64),
    ) -> Self {
        let n = n_lat * n_lon;
        let (dphi, dlon) = (PI / (n_lat.max(2) - 1) as f64, TAU / n_lon.max(1) as f64);
        let lat = |i: usize| -PI / 2.0 + dphi * i as f64;
        let edges = |i: usize| {
            (
                (lat(i) - dphi / 2.0).max(-PI / 2.0),
                (lat(i) + dphi / 2.0).min(PI / 2.0),
            )
        };
        let a = radius;
        let mut faces = Vec::with_capacity(2 * n);
        for p in 0..n.min(sea.len()) {
            let (i, j) = (p / n_lon, p % n_lon);
            let (lo, hi) = edges(i);
            let east = i * n_lon + (j + 1) % n_lon;
            if east != p && sea[p].min(sea[east]) > 0.0 {
                let mid = 0.5 * (lo + hi);
                faces.push(Face {
                    from: p,
                    to: east,
                    north: false,
                    length: a * (hi - lo) * sea[p].min(sea[east]),
                    distance: a * mid.cos() * dlon,
                    coriolis: 2.0 * omega * mid.sin(),
                });
            }
            let north = p + n_lon;
            if i + 1 < n_lat && sea[p].min(sea[north]) > 0.0 {
                faces.push(Face {
                    from: p,
                    to: north,
                    north: true,
                    length: a * hi.cos() * dlon * sea[p].min(sea[north]),
                    distance: a * dphi,
                    coriolis: 2.0 * omega * hi.sin(),
                });
            }
        }
        let mut about = vec![Vec::new(); n];
        for (f, face) in faces.iter().enumerate() {
            about[face.from].push((f, true));
            about[face.to].push((f, false));
        }
        let share = |v: &[f64], p: usize| v.get(p).copied().unwrap_or(0.0);
        let area = (0..n)
            .map(|p| TAU * a * a * share(weight, p) * share(sea, p))
            .collect();
        Self {
            faces,
            area,
            about,
            densities,
            pressure: vec![0.0; n],
            carried: 0.0,
            held: 0.0,
            spent: 0.0,
            steps: 0,
            substeps: 0,
            sweeps: 0,
        }
    }

    /// The drift through its steps so far: the ice's mean speed across the boundaries, m/s,
    /// its mass's weight; the time the steps took, s; and their mean sub-steps and sweeps.
    pub fn summary(&self) -> (f64, f64, f64, f64) {
        let per = |v: usize| v as f64 / self.steps.max(1) as f64;
        let speed = if self.held > 0.0 {
            self.carried / self.held
        } else {
            0.0
        };
        (speed, self.spent, per(self.substeps), per(self.sweeps))
    }

    /// Drifts the ice `ice`, kg/m² of each place's sea, through a step of `dt` s under the
    /// wind's stress on it `stress`, N/m² east and north at each place: each boundary's free
    /// drift under the stress about it, slowed where the ice beyond it bears a pressure the
    /// ice before it does not, and the ice each boundary passes from its upwind place, in steps
    /// that move no more than OCEAN_CFL of a place's ice.
    pub fn step(&mut self, ice: &mut [f64], stress: &[[f64; 2]], dt: f64) {
        self.step_with(ice, &mut [], stress, dt);
    }

    /// `step`, carrying with the ice what rides on it, `riders`, kg/m² of each place's sea
    /// (M64: the snow on it), by the same shares of each place's: none where `riders` is
    /// not one to a place.
    pub fn step_with(&mut self, ice: &mut [f64], riders: &mut [f64], stress: &[[f64; 2]], dt: f64) {
        let n = self.area.len();
        if ice.len() != n || stress.len() != n || dt <= 0.0 {
            return;
        }
        let started = std::time::Instant::now();
        let (water, density) = self.densities;
        let drag = water * SEA_ICE_WATER_DRAG;
        let cover: Vec<f64> = ice
            .iter()
            .map(|&m| 1.0 - (-m.max(0.0) / SEA_ICE_COVER_MASS).exp())
            .collect();
        // Hibler's strength: the ice's mean thickness over the sea, weakened as its cover
        // opens.
        let strength: Vec<f64> = ice
            .iter()
            .zip(&cover)
            .map(|(&m, &c)| {
                SEA_ICE_STRENGTH * m.max(0.0) / density
                    * (-SEA_ICE_STRENGTH_DECAY * (1.0 - c)).exp()
            })
            .collect();
        // Each boundary's free drift across it, m/s, from the stress about it and the mass of
        // the ice over its own area there; and how much the pressure on its two sides slows
        // it, m/s per N/m, its drag taken at the free drift's speed.
        let mut free = vec![0.0; self.faces.len()];
        let mut give = vec![0.0; self.faces.len()];
        for (f, face) in self.faces.iter().enumerate() {
            let (p, q) = (face.from, face.to);
            let (held, covered) = (ice[p].max(0.0) + ice[q].max(0.0), cover[p] + cover[q]);
            if held <= 0.0 || covered <= 0.0 {
                continue;
            }
            let tau = [
                0.5 * (stress[p][0] + stress[q][0]),
                0.5 * (stress[p][1] + stress[q][1]),
            ];
            let turn = held / covered * face.coriolis;
            let (u, k) = free_drift(tau, drag, turn);
            free[f] = if face.north { u[1] } else { u[0] };
            let slows = k * k + turn * turn;
            if slows > 0.0 {
                give[f] = k / (slows * face.distance * 0.5 * covered);
            }
        }
        // The pressure the ice bears: at each place the one under which the ice its
        // boundaries pass in and out balances, none where it would be drawn apart, its
        // strength at most, beyond which it gives; in over-relaxed sweeps from the last
        // step's, until no place's moves by SEA_ICE_PRESSURE_TOLERANCE of the strongest ice's
        // strength, so that the pressure is the step's own and not the last's carried on.
        let mut p = std::mem::take(&mut self.pressure);
        p.resize(n, 0.0);
        for (v, &s) in p.iter_mut().zip(&strength) {
            *v = v.clamp(0.0, s);
        }
        let pressed: Vec<usize> = (0..n).filter(|&k| strength[k] > 0.0).collect();
        let settled = SEA_ICE_PRESSURE_TOLERANCE * strength.iter().copied().fold(0.0, f64::max);
        for _ in 0..SEA_ICE_PRESSURE_SWEEPS {
            self.sweeps += 1;
            let mut moved = 0.0f64;
            for &k in &pressed {
                let (mut out, mut stiff) = (0.0, 0.0);
                for &(f, west) in &self.about[k] {
                    let face = &self.faces[f];
                    let w = free[f] - give[f] * (p[face.to] - p[face.from]);
                    out += if west { w } else { -w } * face.length;
                    stiff += give[f] * face.length;
                }
                if stiff > 0.0 {
                    let next = p[k] - SEA_ICE_PRESSURE_RELAXATION * out / stiff;
                    let next = next.clamp(0.0, strength[k]);
                    moved = moved.max((next - p[k]).abs());
                    p[k] = next;
                }
            }
            if moved <= settled {
                break;
            }
        }
        let speed: Vec<f64> = self
            .faces
            .iter()
            .enumerate()
            .map(|(f, face)| free[f] - give[f] * (p[face.to] - p[face.from]))
            .collect();
        self.pressure = p;
        // The ice each boundary carries from its upwind place, the share of that place's ice
        // a second, in steps that move no more than OCEAN_CFL of a place's ice, none more than
        // it holds.
        let mut leaving = vec![0.0; n];
        let mut moving: Vec<(usize, usize, f64, f64)> = Vec::new();
        for (face, &w) in self.faces.iter().zip(&speed) {
            let (up, down) = if w > 0.0 {
                (face.from, face.to)
            } else {
                (face.to, face.from)
            };
            if w != 0.0 && self.area[up] > 0.0 {
                let rate = w.abs() * face.length / self.area[up];
                leaving[up] += rate;
                moving.push((up, down, rate, w.abs()));
            }
        }
        let fastest = leaving.iter().copied().fold(0.0, f64::max);
        let subs = ((dt * fastest / OCEAN_CFL).ceil() as usize).clamp(1, OCEAN_MOST_SUBSTEPS);
        let h = dt / subs as f64;
        for m in moving.iter_mut() {
            m.2 *= h / (h * leaving[m.0]).max(1.0);
        }
        let mut mass: Vec<f64> = ice
            .iter()
            .zip(&self.area)
            .map(|(&m, &a)| m.max(0.0) * a)
            .collect();
        let riding = riders.len() == n;
        let mut rode: Vec<f64> = match riding {
            true => riders
                .iter()
                .zip(&self.area)
                .map(|(&m, &a)| m.max(0.0) * a)
                .collect(),
            false => Vec::new(),
        };
        let mut moved = vec![0.0; n];
        let mut moved_riders = vec![0.0; rode.len()];
        for _ in 0..subs {
            moved.fill(0.0);
            moved_riders.fill(0.0);
            for &(up, down, share, speed) in &moving {
                let carried = mass[up] * share;
                moved[up] -= carried;
                moved[down] += carried;
                self.carried += carried * speed;
                self.held += carried;
                if riding {
                    let rides = rode[up] * share;
                    moved_riders[up] -= rides;
                    moved_riders[down] += rides;
                }
            }
            for (m, d) in mass.iter_mut().zip(&moved) {
                *m = (*m + d).max(0.0);
            }
            for (m, d) in rode.iter_mut().zip(&moved_riders) {
                *m = (*m + d).max(0.0);
            }
        }
        for ((m, &held), &a) in ice.iter_mut().zip(&mass).zip(&self.area) {
            if a > 0.0 {
                *m = held / a;
            }
        }
        for ((m, &held), &a) in riders.iter_mut().zip(&rode).zip(&self.area) {
            if a > 0.0 {
                *m = held / a;
            }
        }
        self.spent += started.elapsed().as_secs_f64();
        self.steps += 1;
        self.substeps += subs;
    }
}

/// The free drift under the stress `tau`, N/m² east and north (F276): the drift, m/s east and
/// north, at which the still water's drag on the ice, `drag` |u| u with `drag` its density
/// times its drag coefficient, and the spin's turn of the ice's mass, i `turn` u with `turn`
/// its mass over its own area times the Coriolis parameter, balance it; and the drag's
/// coefficient at that speed, `drag` |u|, kg/m²/s.
fn free_drift(tau: [f64; 2], drag: f64, turn: f64) -> ([f64; 2], f64) {
    let pull = tau[0] * tau[0] + tau[1] * tau[1];
    if pull <= 0.0 || drag <= 0.0 {
        return ([0.0; 2], 0.0);
    }
    // |τ|² = (drag² |u|² + turn²) |u|²: its root for |u|², in the form that keeps where the
    // turn is the larger.
    let b2 = turn * turn;
    let s = 2.0 * pull / (b2 + (b2 * b2 + 4.0 * drag * drag * pull).sqrt());
    let k = drag * s.sqrt();
    let d = k * k + b2;
    (
        [
            (tau[0] * k + tau[1] * turn) / d,
            (tau[1] * k - tau[0] * turn) / d,
        ],
        k,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    /// F276: under a 10 m/s wind's stress the free drift runs at some 2 percent of the wind,
    /// turned to its right in the north and to its left in the south, the more under a
    /// gentler wind.
    #[test]
    fn the_free_drift_runs_at_some_two_percent_of_the_wind_turned_by_the_spin() {
        let drag_air = (VON_KARMAN / (WIND_REFERENCE_HEIGHT / ICE_ROUGHNESS).ln()).powi(2);
        let drag = 1_025.0 * SEA_ICE_WATER_DRAG;
        let f = 2.0 * 7.292e-5 * 75f64.to_radians().sin();
        let mut angles = Vec::new();
        for wind in [3.0, 10.0] {
            let tau = [1.3 * drag_air * wind * wind, 0.0];
            let (north, _) = free_drift(tau, drag, 1_800.0 * f);
            let (south, _) = free_drift(tau, drag, -1_800.0 * f);
            let speed = north[0].hypot(north[1]);
            let angle = north[1].atan2(north[0]).to_degrees();
            println!(
                "{wind} m/s: the ice {:.3} m/s, {:.2} percent, turned {angle:.1}°",
                speed,
                100.0 * speed / wind
            );
            assert!((0.015..0.025).contains(&(speed / wind)), "{speed}");
            assert!(angle < 0.0 && south[1] > 0.0, "{north:?} {south:?}");
            angles.push(-angle);
        }
        assert!(angles[0] > angles[1] && angles[0] > 15.0, "{angles:?}");
    }

    /// M64: what rides on the ice, its snow, drifts with it: a tenth of each place's ice
    /// stays a tenth wherever the ice is blown, and the snow over the sea is kept.
    #[test]
    fn the_snow_rides_with_its_ice() {
        let (n_lat, n_lon) = (91, 90);
        let (dphi, dlon) = (PI / (n_lat - 1) as f64, TAU / n_lon as f64);
        let lat = |i: usize| -PI / 2.0 + dphi * i as f64;
        let lon = |j: usize| -PI + dlon * j as f64;
        let inside = |p: usize| {
            let (la, lo) = (lat(p / n_lon).to_degrees(), lon(p % n_lon).to_degrees());
            (70.0..=86.0).contains(&la) && (0.0..60.0).contains(&lo)
        };
        let n = n_lat * n_lon;
        let sea: Vec<f64> = (0..n).map(|p| if inside(p) { 1.0 } else { 0.0 }).collect();
        let weight: Vec<f64> = (0..n)
            .map(|p| {
                let la = lat(p / n_lon);
                let (lo, hi) = (
                    (la - dphi / 2.0).max(-PI / 2.0),
                    (la + dphi / 2.0).min(PI / 2.0),
                );
                (hi.sin() - lo.sin()) / n_lon as f64
            })
            .collect();
        let mut drift = Drift::new(
            (&sea, &weight),
            (n_lat, n_lon),
            (1.4e-4, 6.371e6),
            (1_025.0, 917.0),
        );
        // Ice thinning eastward, its snow a tenth of it.
        let mut ice: Vec<f64> = (0..n)
            .map(|p| sea[p] * (2_000.0 - 20.0 * lon(p % n_lon).to_degrees()))
            .collect();
        let mut snow: Vec<f64> = ice.iter().map(|m| 0.1 * m).collect();
        let area = drift.area.clone();
        let total = |v: &[f64]| v.iter().zip(&area).map(|(m, a)| m * a).sum::<f64>();
        let before = total(&snow);
        let stress = vec![[0.05, 0.02]; n];
        for _ in 0..60 {
            drift.step_with(&mut ice, &mut snow, &stress, 86_400.0);
        }
        let ratio = (0..n)
            .filter(|&p| ice[p] > 1.0)
            .map(|p| (snow[p] / ice[p] - 0.1).abs())
            .fold(0.0, f64::max);
        println!("the snow's share of its ice strays by {ratio:.2e} at most");
        assert!(ratio < 1e-9, "{ratio}");
        assert!(((total(&snow) - before) / before).abs() < 1e-9);
    }

    /// F276: ice a metre thick over a sea from 0° to 60° east between 70° and 86° north,
    /// blown east at 0.05 N/m² without the spin's turn, drifts against its eastern coast and
    /// piles there, its pressure rising toward the coast as the wind's stress over the sheet:
    /// the sheet holds where that stays under its strength, and beyond it the ice piles to
    /// the wind's stress times its distance from the sheet's edge over P*; the western coast
    /// opens; the ice's mass is kept.
    #[test]
    fn ice_blown_on_a_coast_piles_there_until_its_strength_holds_the_wind() {
        let (n_lat, n_lon) = (91, 90);
        let (dphi, dlon) = (PI / (n_lat - 1) as f64, TAU / n_lon as f64);
        let lat = |i: usize| -PI / 2.0 + dphi * i as f64;
        let lon = |j: usize| -PI + dlon * j as f64;
        let inside = |p: usize| {
            let (la, lo) = (lat(p / n_lon).to_degrees(), lon(p % n_lon).to_degrees());
            (70.0..=86.0).contains(&la) && (0.0..60.0).contains(&lo)
        };
        let n = n_lat * n_lon;
        let sea: Vec<f64> = (0..n).map(|p| if inside(p) { 1.0 } else { 0.0 }).collect();
        let weight: Vec<f64> = (0..n)
            .map(|p| {
                let la = lat(p / n_lon);
                let (lo, hi) = (
                    (la - dphi / 2.0).max(-PI / 2.0),
                    (la + dphi / 2.0).min(PI / 2.0),
                );
                (hi.sin() - lo.sin()) / n_lon as f64
            })
            .collect();
        let (a, density) = (6.371e6, 917.0);
        let mut drift = Drift::new(
            (&sea, &weight),
            (n_lat, n_lon),
            (0.0, a),
            (1_025.0, density),
        );
        let start = 917.0;
        let mut ice: Vec<f64> = sea.iter().map(|&s| s * start).collect();
        let area = drift.area.clone();
        let total = |ice: &[f64]| ice.iter().zip(&area).map(|(m, a)| m * a).sum::<f64>();
        let before = total(&ice);
        let tau = 0.05;
        let stress = vec![[tau, 0.0]; n];
        for _ in 0..400 {
            drift.step(&mut ice, &stress, 86_400.0);
        }
        let after = total(&ice);
        assert!(((after - before) / before).abs() < 1e-9, "{before} {after}");
        // Along the row at 78° north.
        let i = (0..n_lat)
            .min_by_key(|&i| ((lat(i).to_degrees() - 78.0).abs() * 100.0) as i64)
            .unwrap();
        let row: Vec<(f64, f64, f64)> = (0..n_lon)
            .map(|j| i * n_lon + j)
            .filter(|&p| inside(p))
            .map(|p| {
                (
                    lon(p % n_lon).to_degrees(),
                    ice[p] / density,
                    drift.pressure[p],
                )
            })
            .collect();
        for (lo, h, p) in &row {
            println!("{lo:6.1}°E: {h:.2} m, {:.1} kN/m", p / 1e3);
        }
        let (west, coast) = (row.first().unwrap().1, row.last().unwrap().1);
        // The sheet's length from the ice it holds: a metre where the wind's stress over it
        // stays under its strength, h = τ x / P* beyond; and its thickness at the coast.
        let (width, h0) = (60f64.to_radians() * a * lat(i).cos(), start / density);
        let rigid = (SEA_ICE_STRENGTH * h0 / tau).min(width);
        let pile = (rigid * rigid + 2.0 * SEA_ICE_STRENGTH * h0 * (width - rigid) / tau).sqrt();
        let expected = (tau * pile / SEA_ICE_STRENGTH).max(h0);
        println!("the coast's {coast:.2} m against the pile's {expected:.2} m");
        assert!(west < 0.2 * start / density, "{west}");
        assert!(
            coast > 0.7 * expected && coast < 1.3 * expected,
            "{coast} {expected}"
        );
        // The pressure rises toward the coast as the wind's stress over the pile.
        let (x0, x1) = (row.len() / 2, row.len() - 2);
        let rise =
            (row[x1].2 - row[x0].2) / ((row[x1].0 - row[x0].0).to_radians() * a * lat(i).cos());
        println!("its pressure rising {rise:.4} N/m² against the wind's {tau}");
        assert!(rise > 0.6 * tau && rise < 1.4 * tau, "{rise}");
    }

    /// The drift keeps the ice: over seas covering every share of their places, land among
    /// them, under a spin and a wind that turns from place to place and step to step, a season
    /// of steps neither makes nor loses ice and leaves none below none (the tests' tiers).
    #[test]
    fn the_drift_keeps_the_ice() {
        let (n_lat, n_lon) = (19, 36);
        let n = n_lat * n_lon;
        let dphi = PI / (n_lat - 1) as f64;
        let rng = crate::rng::Rng::new(7);
        let weight: Vec<f64> = (0..n)
            .map(|p| {
                let la = -PI / 2.0 + dphi * (p / n_lon) as f64;
                let (lo, hi) = (
                    (la - dphi / 2.0).max(-PI / 2.0),
                    (la + dphi / 2.0).min(PI / 2.0),
                );
                (hi.sin() - lo.sin()) / n_lon as f64
            })
            .collect();
        let sea: Vec<f64> = (0..n as u64)
            .map(|p| match rng.uniform(p, 0) {
                x if x < 0.3 => 0.0,
                x if x < 0.5 => 1.0,
                x => x,
            })
            .collect();
        let mut drift = Drift::new(
            (&sea, &weight),
            (n_lat, n_lon),
            (7.292e-5, 6.371e6),
            (1_025.0, 917.0),
        );
        let mut ice: Vec<f64> = (0..n)
            .map(|p| 3.0 * 917.0 * rng.uniform(p as u64, 1) * sea[p].min(1.0).ceil())
            .collect();
        let area = drift.area.clone();
        let total = |ice: &[f64]| ice.iter().zip(&area).map(|(m, a)| m * a).sum::<f64>();
        let before = total(&ice);
        for step in 0..90u64 {
            let stress: Vec<[f64; 2]> = (0..n as u64)
                .map(|p| {
                    [
                        rng.range(p, 2 + 2 * step, -0.2, 0.2),
                        rng.range(p, 3 + 2 * step, -0.2, 0.2),
                    ]
                })
                .collect();
            drift.step(&mut ice, &stress, 86_400.0);
        }
        let after = total(&ice);
        println!("the ice {before:.6e} kg before, {after:.6e} after");
        assert!(((after - before) / before).abs() < 1e-9, "{before} {after}");
        for (p, &m) in ice.iter().enumerate() {
            assert!(m.is_finite() && m >= 0.0, "place {p}: {m}");
        }
    }
}
