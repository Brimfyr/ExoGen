// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! The liquid in the ground (F435).
//!
//! What a body's land sheds soaks into its ground and flows through it down the slope of its
//! table, the surface below which the ground is full of the liquid, to where the table meets
//! the surface: there it seeps out, into the rivers, the lakes and the sea. The flow is
//! Darcy's, its conductivity K = k rho g / mu, the ground's permeability k by its materials, the
//! liquid's density rho and viscosity mu and the body's gravity g; level through the ground's
//! depth below the table (Dupuit-Forchheimer), its transmissivity T the conductivity summed
//! down from the table: through the loose cover at the cover's own, then through the rock at
//! the rock's near the surface, to the depth below which it falls as the continental crust's
//! does with depth, k = 10^-14 z^-3.2, z in km (Manning and Ingebritsen 1999), the depth read
//! at the pressure the body's own gravity gives it (as Clifford and Parker 2001 read it on
//! Mars). Where the table lies near the surface the ground evaporates from it what the air
//! would take from a wet ground beyond what its soil gives, less as the table lies deeper,
//! e-fold each GROUND_DRYING_DEPTH (Faulk et al.'s ground methane evaporation; F436). At the
//! steady state what soaks in balances what flows out and what the ground gives the air,
//! div(T grad h) + R - E = 0, where the table lies below the surface; where it would rise above
//! it, it stands at the surface and the rest seeps out; the sea and the lakes hold it at their
//! levels and take or give what flows to them (Callaghan, Wickert, Barnes and Austermann 2025's
//! Water Table Model, at its steady state). By finite volumes over each cell's sides, the
//! transmissivity at a side the mean of its two cells', its flow the side's length over the
//! distance between the cells' centres; the cells at the surface found as an active set
//! (Hintermüller, Ito and Kunisch 2002), the rest by conjugate gradients, the transmissivity
//! lagged a step (Picard).

use crate::constants::*;
use crate::fields::det_sum;
use crate::geom::dot;
use crate::grid::CubeSphere;
use rayon::prelude::*;

/// The sides of each cell of a grid and their shape: the cells either side, both ways, each
/// cell's in `start[k]..start[k + 1]`, and each side's length over the distance between the
/// centres of its cells.
#[derive(Clone, Debug)]
pub struct Sides {
    start: Vec<u32>,
    from: Vec<u32>,
    to: Vec<u32>,
    shape: Vec<f64>,
}

impl Sides {
    /// The sides of `grid`'s cells, each cell's four within its face and across its edges; a
    /// side its two cells find differently near a face's corner taken at their mean, one only
    /// one of them finds as that one finds it.
    pub fn new(grid: CubeSphere) -> Self {
        let n = grid.cell_count();
        let delta = grid.delta();
        let found: Vec<[(u32, f64); 4]> = (0..n)
            .into_par_iter()
            .map(|k| {
                let c = grid.cell(k);
                let (a, b) = grid.angles(c);
                let here = grid.center(c);
                let h = 0.5 * delta;
                [(1, 0), (-1, 0), (0, 1), (0, -1)].map(|(di, dj): (i32, i32)| {
                    let m = grid.neighbor(c, di, dj);
                    let there = grid.center(m);
                    // The side: from one corner to the other along the cell's edge that way.
                    let (p, q) = if di != 0 {
                        let x = a + di as f64 * h;
                        (
                            CubeSphere::direction(c.face, x, b - h),
                            CubeSphere::direction(c.face, x, b + h),
                        )
                    } else {
                        let y = b + dj as f64 * h;
                        (
                            CubeSphere::direction(c.face, a - h, y),
                            CubeSphere::direction(c.face, a + h, y),
                        )
                    };
                    (
                        grid.index(m) as u32,
                        arc(p, q) / arc(here, there).max(1e-12),
                    )
                })
            })
            .collect();
        let mut pairs: Vec<(u32, u32, f64)> = Vec::with_capacity(4 * n);
        for (k, sides) in found.iter().enumerate() {
            for &(m, shape) in sides {
                if m as usize != k {
                    pairs.push((k as u32, m, shape));
                }
            }
        }
        // Each side once, its two cells' shapes met at their mean.
        let mut once: Vec<(u32, u32, f64, u32)> = pairs
            .iter()
            .map(|&(k, m, s)| (k.min(m), k.max(m), s, 1))
            .collect();
        once.sort_by_key(|x| (x.0, x.1));
        let mut merged: Vec<(u32, u32, f64)> = Vec::with_capacity(once.len());
        let mut i = 0;
        while i < once.len() {
            let (a, b) = (once[i].0, once[i].1);
            let (mut sum, mut count) = (0.0, 0u32);
            while i < once.len() && once[i].0 == a && once[i].1 == b {
                sum += once[i].2;
                count += once[i].3;
                i += 1;
            }
            merged.push((a, b, sum / count as f64));
        }
        let mut count = vec![0u32; n + 1];
        for &(a, b, _) in &merged {
            count[a as usize + 1] += 1;
            count[b as usize + 1] += 1;
        }
        for k in 0..n {
            count[k + 1] += count[k];
        }
        let start = count.clone();
        let mut fill = count;
        let mut from = vec![0u32; start[n] as usize];
        let mut to = vec![0u32; start[n] as usize];
        let mut shape = vec![0.0; start[n] as usize];
        for &(a, b, s) in &merged {
            for (x, y) in [(a, b), (b, a)] {
                let at = fill[x as usize] as usize;
                (from[at], to[at], shape[at]) = (x, y, s);
                fill[x as usize] += 1;
            }
        }
        Self {
            start,
            from,
            to,
            shape,
        }
    }

    /// The sides of cell `k`: the cell across each and its shape.
    #[cfg(test)]
    fn of(&self, k: usize) -> impl Iterator<Item = (usize, f64)> + '_ {
        let (a, b) = (self.start[k] as usize, self.start[k + 1] as usize);
        (a..b).map(move |s| (self.to[s] as usize, self.shape[s]))
    }
}

/// The angle between two unit directions, rad.
fn arc(p: [f64; 3], q: [f64; 3]) -> f64 {
    let c = dot(p, q);
    let x = [
        p[1] * q[2] - p[2] * q[1],
        p[2] * q[0] - p[0] * q[2],
        p[0] * q[1] - p[1] * q[0],
    ];
    (dot(x, x).sqrt()).atan2(c)
}

/// A cell's ground as the liquid flows through it (F435): its loose cover's thickness, m, and
/// conductivity, m/s; its rock's conductivity near the surface, m/s, and the depth below the
/// surface, m, from which it falls as the crust's does, (z / knee)^-GROUND_CRUST_FALL.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Profile {
    pub cover: f64,
    pub cover_conductivity: f64,
    pub rock_conductivity: f64,
    pub knee: f64,
}

impl Profile {
    /// The transmissivity, m²/s, below a table `d` m under the surface: K_c (b - d) through
    /// the cover above it, and the rock's below the deeper of the table and the cover's base
    /// a, K_r (z_k - a) + K_r z_k / (n - 1) from above its knee z_k, K_r (a / z_k)^-n a / (n - 1)
    /// below it, n GROUND_CRUST_FALL.
    pub fn transmissivity(&self, d: f64) -> f64 {
        let d = d.max(0.0);
        let fall = GROUND_CRUST_FALL - 1.0;
        let rock = |a: f64| {
            let k = self.rock_conductivity;
            if a < self.knee {
                k * (self.knee - a) + k * self.knee / fall
            } else {
                k * (a / self.knee).powf(-GROUND_CRUST_FALL) * a / fall
            }
        };
        if d < self.cover {
            self.cover_conductivity * (self.cover - d) + rock(self.cover)
        } else {
            rock(d)
        }
    }

    /// A cell's ground (F435): its loose cover `cover` m thick, of the rivers' sediments, over
    /// a rock of permeability `rock`, m², near the surface, both frozen through where `frozen`;
    /// the rock's knee where the crust's law, GROUND_CRUST_PERMEABILITY at a kilometre falling as
    /// z^-GROUND_CRUST_FALL, comes down to its own, the depth read at the pressure `pressure`
    /// times the Earth's at the same depth; the liquid flowing through it as `flow` = rho g /
    /// mu, 1/(m s).
    pub fn of(cover: f64, rock: f64, frozen: bool, (flow, pressure): (f64, f64)) -> Self {
        let (loose, rock) = if frozen {
            (GROUND_FROZEN_PERMEABILITY, GROUND_FROZEN_PERMEABILITY)
        } else {
            (GROUND_LOOSE_PERMEABILITY, rock)
        };
        let knee = 1_000.0 * (rock / GROUND_CRUST_PERMEABILITY).powf(-1.0 / GROUND_CRUST_FALL)
            / pressure.max(1e-6);
        Profile {
            cover: cover.max(0.0),
            cover_conductivity: loose * flow,
            rock_conductivity: rock * flow,
            knee,
        }
    }

    /// The most a cell's ground takes in at its top, m/s: its top's conductivity, the
    /// gradient one under gravity alone.
    pub fn intake(&self) -> f64 {
        if self.cover > 0.0 {
            self.cover_conductivity
        } else {
            self.rock_conductivity
        }
    }
}

/// The permeability of a crust of `bedrock`'s kind, m² (F435): a basaltic crust's and
/// sulfur's as volcanic rock's, the rest as crystalline rock's.
pub fn rock_permeability(bedrock: &str) -> f64 {
    match bedrock {
        "basalt" | "sulfur_deposit" => GROUND_VOLCANIC_PERMEABILITY,
        _ => GROUND_CRYSTALLINE_PERMEABILITY,
    }
}

/// Each cell's ground (F435): its loose cover `cover` m thick over a crust of `bedrock`'s
/// kind, frozen through where its ground's mean temperature, `temperature`, K, stands below its
/// `liquid`'s melting point, the liquid's viscosity at that temperature, under `gravity`,
/// m/s², the crust's depth read at the pressure that gravity gives it over the Earth's,
/// GROUND_CRUST_GRAVITY; a cell whose temperature is not known read at the melting point.
pub fn profiles(
    cover: &[f32],
    temperature: &(dyn Fn(usize) -> f64 + Sync),
    bedrock: &str,
    (liquid, gravity): (&Volatile, f64),
) -> Vec<Profile> {
    let rock = rock_permeability(bedrock);
    (0..cover.len())
        .into_par_iter()
        .map(|k| {
            let t = Some(temperature(k))
                .filter(|t| t.is_finite() && *t > 0.0)
                .unwrap_or(liquid.t_triple);
            let flow = liquid.liquid_density * gravity / viscosity(liquid, t);
            let pressure = gravity / GROUND_CRUST_GRAVITY;
            Profile::of(cover[k] as f64, rock, t < liquid.t_triple, (flow, pressure))
        })
        .collect()
}

/// A liquid's viscosity at `t` K, Pa s (F435): its fit, read from its triple point up.
pub fn viscosity(liquid: &Volatile, t: f64) -> f64 {
    let [a, b, c] = liquid.liquid_viscosity;
    a * (b / (t.max(liquid.t_triple) - c)).exp()
}

/// The steady table (F435): its head at each cell, m; what seeps out where it stands at the
/// surface, m³/s, none elsewhere; and what flows through the ground into each cell held at
/// its level, the sea's and the lakes', m³/s, less where the cell gives the ground.
#[derive(Clone, Debug)]
pub struct Table {
    pub head: Vec<f64>,
    pub seep: Vec<f64>,
    pub exchange: Vec<f64>,
    /// Whether the table stands at the surface at each cell.
    pub wet: Vec<bool>,
    /// The level each cell is held at, the sea's or a lake's, m, taking in nothing; not a
    /// number where it is not held.
    pub held: Vec<f64>,
    /// What the ground gives the air at each cell, m³/s.
    pub dried: Vec<f64>,
    /// What the lakes over the table below them give it at each cell, m³/s, where the lakes
    /// stand on it (`drain_through`); none from `table` alone.
    pub leak: Vec<f64>,
    /// The ground's year in brief, where the lakes stand on it (`drain_through`).
    pub brief: Brief,
    /// The active set's rounds through all the passes, and the conjugate gradients' steps
    /// through them.
    pub rounds: usize,
    pub steps: usize,
}

/// The ground's year in brief (F435): what soaks into it, what it gives the air, what seeps out
/// of it at the surface, what flows through it into the lakes, less what they give it, and
/// into the sea, what the lakes over a table below them give it, and what the finish's lakes
/// give it on the rivers' finer grid where its own over a table below them do
/// (`route_through`), m³/s; the share of the land where its table meets the surface; and the
/// rounds its table was found again with the lakes as they stood.
#[derive(Clone, Copy, Debug, Default)]
pub struct Brief {
    pub soaked: f64,
    pub dried: f64,
    pub seeped: f64,
    pub to_lakes: f64,
    pub to_sea: f64,
    pub leaked: f64,
    pub drawn: f64,
    pub wet: f64,
    pub coupled: usize,
}

/// The steady table on `grid` of a body of `radius` m (F435): over a ground whose `surface`
/// is at each cell's height, m, of each cell's `profile`, the cells where `held` is finite held
/// at that head (the sea's and the lakes'), m, under what soaks in at each other, `recharge`,
/// m/s, and what the air would take from its ground wet at the surface beyond what its soil
/// gives it, `demand`, m/s, e^(-d / GROUND_DRYING_DEPTH) of it where the table lies d below the
/// surface (F436); from `start`'s heads where given, else from the table found free.
pub fn table(
    (grid, sides, radius): (CubeSphere, &Sides, f64),
    surface: &[f64],
    held: &[f64],
    (recharge, demand): (&[f64], &[f64]),
    profile: &[Profile],
    start: Option<&[f64]>,
) -> Table {
    let n = grid.cell_count();
    let area: Vec<f64> = (0..n)
        .into_par_iter()
        .map(|k| grid.cell_area(grid.cell(k)) * radius * radius)
        .collect();
    let fixed: Vec<bool> = held.iter().map(|h| h.is_finite()).collect();
    let mut head: Vec<f64> = (0..n)
        .map(|k| {
            if fixed[k] {
                held[k]
            } else {
                start.map_or(surface[k], |s| s[k].min(surface[k]))
            }
        })
        .collect();
    // From a start's cells at the surface; with none, from none: the table found free first,
    // the cells it would raise over the surface then held at it.
    let mut wet: Vec<bool> = (0..n)
        .map(|k| start.is_some() && !fixed[k] && head[k] >= surface[k])
        .collect();
    let source: Vec<f64> = (0..n)
        .map(|k| if fixed[k] { 0.0 } else { recharge[k] * area[k] })
        .collect();
    // What a held cell's bed passes at most, m³/s: what its ground takes in at its top over its
    // area. Where the ground would draw more from it, the table beneath falls from its level,
    // the liquid perched over it giving it that; where the table beneath rises to its level
    // again, it is held again.
    let bed: Vec<f64> = (0..n).map(|k| profile[k].intake() * area[k]).collect();
    let mut perched = vec![false; n];
    let mut fed = source.clone();
    // What the ground gives the air at each cell, m³/s, and its rise with the table, m²/s.
    let drying = |head: &[f64], k: usize| -> (f64, f64) {
        if fixed[k] {
            return (0.0, 0.0);
        }
        let d = (surface[k] - head[k]).max(0.0);
        let dried = demand[k] * area[k] * (-d / GROUND_DRYING_DEPTH).exp();
        (dried, dried / GROUND_DRYING_DEPTH)
    };
    let mut conductance = vec![0.0; sides.to.len()];
    let mut residual = vec![0.0; n];
    let mut rise = vec![0.0; n];
    let (mut rounds, mut steps) = (0, 0);
    let tolerance = 1e-9 * det_sum(n, |k| source[k].abs()).max(1e-30);
    // A cell at the surface leaves it where it would draw more than this from it, m³/s.
    let least = GROUND_DRAW_SHARE * tolerance / 1e-9 / n as f64;
    // Each pass the sides' conductances at their cells' transmissivities as the last pass left
    // the table (Picard); within a pass, the cells at the surface found as an active set.
    for _ in 0..GROUND_PASSES {
        conduct(sides, profile, (surface, &head), &mut conductance);
        let before = head.clone();
        for _ in 0..GROUND_ROUNDS {
            rounds += 1;
            net(sides, &conductance, &head, &fed, &mut residual);
            // What the ground gives the air, its rise with the table taken into the step
            // (Newton's).
            for k in 0..n {
                let (dried, slope) = drying(&head, k);
                residual[k] -= dried;
                rise[k] = slope;
            }
            // The cells at the surface that would draw from it leave it.
            let mut changed = 0usize;
            for k in 0..n {
                if wet[k] && residual[k] < -least {
                    wet[k] = false;
                    changed += 1;
                }
                // The held cells the ground would draw more from than their beds pass perch.
                if fixed[k] && !perched[k] && residual[k] < -bed[k] - least {
                    perched[k] = true;
                    fed[k] = bed[k];
                    residual[k] += bed[k];
                    changed += 1;
                }
            }
            let free: Vec<bool> = (0..n)
                .map(|k| (!fixed[k] || perched[k]) && !wet[k])
                .collect();
            let (delta, taken) = solve(sides, (&conductance, &rise), &free, &residual, tolerance);
            steps += taken;
            for k in 0..n {
                if free[k] {
                    head[k] += delta[k];
                    if perched[k] {
                        if head[k] > held[k] {
                            head[k] = held[k];
                            perched[k] = false;
                            fed[k] = 0.0;
                            changed += 1;
                        }
                    } else if head[k] > surface[k] {
                        head[k] = surface[k];
                        wet[k] = true;
                        changed += 1;
                    }
                }
            }
            if changed == 0 {
                break;
            }
        }
        let most = (0..n)
            .map(|k| (head[k] - before[k]).abs())
            .fold(0.0, f64::max);
        if most < GROUND_SETTLED {
            break;
        }
    }
    // What seeps out at the surface, what the ground gives the air, and what each held cell
    // takes from the ground.
    conduct(sides, profile, (surface, &head), &mut conductance);
    net(sides, &conductance, &head, &fed, &mut residual);
    let dried: Vec<f64> = (0..n).map(|k| drying(&head, k).0).collect();
    for k in 0..n {
        residual[k] -= dried[k];
    }
    let seep: Vec<f64> = (0..n)
        .map(|k| if wet[k] { residual[k].max(0.0) } else { 0.0 })
        .collect();
    let exchange: Vec<f64> = (0..n)
        .map(|k| if fixed[k] { residual[k] - fed[k] } else { 0.0 })
        .collect();
    Table {
        head,
        seep,
        exchange,
        wet,
        held: held.to_vec(),
        dried,
        leak: vec![0.0; n],
        brief: Brief::default(),
        rounds,
        steps,
    }
}

/// A ground the liquid flows through (F435): its cells' sides, each cell's profile, and the
/// sea's level, m, the head the table meets the sea at.
pub struct Ground<'a> {
    pub sides: &'a Sides,
    pub profile: &'a [Profile],
    pub level: f64,
}

/// Where the coupling of the ground and the lakes left off (F435): the table's heads, the
/// levels it was held at, NaN where it was not, and what the lakes gave the ground beneath
/// them, m/s.
#[derive(Clone, Debug)]
pub struct Start {
    pub head: Vec<f64>,
    pub held: Vec<f64>,
    pub focused: Vec<f64>,
}

/// The rivers and the lakes of a year's runoff with the ground between them (F435): what the
/// land sheds, `runoff`, kg/m²/s, soaks into its ground as far as the ground takes it in, the
/// rest running off over it; what soaks in flows through the ground to where its table meets
/// the surface, and seeps out there into the rivers; the lakes stand where the rivers bring
/// their water (Fill-Spill-Merge, F277), the table held at their levels where they cover a cell,
/// taking what flows to them through the ground or giving it, and a lake over ground whose
/// table lies below it gives it what its ground takes in, until the ground beneath it fills;
/// each round the table found with the last round's lakes, until they cover the land as they
/// did and stand where they stood, or GROUND_COUPLING rounds have passed, from `start` where
/// given (the year's last), its end returned. A level held the round before moves halfway to
/// the one found, and a cell joins the held where its lake covers GROUND_HOLD_MOST of it and
/// leaves them where it covers less than GROUND_HOLD_LEAST (numerical: the lake found over the
/// ground's last exchange with another level and extent of it overshoots, a lake through which
/// the ground's water flows filling and draining from round to round). The surface's open water
/// gains `gain`, kg/m²/s; the air would take `demand`, kg/m²/s, from the ground wet at its
/// surface beyond what its soil gives it (`table`); the liquid of `density` kg/m³.
pub fn drain_through(
    grid: CubeSphere,
    height: &[f32],
    sea: &[bool],
    (runoff, gain, demand): (&[f32], &[f32], &[f32]),
    (radius, density): (f64, f64),
    ground: &Ground,
    start: Option<&Start>,
) -> (crate::rivers::Drainage, Table, Start) {
    let n = grid.cell_count();
    let area: Vec<f64> = (0..n)
        .into_par_iter()
        .map(|k| grid.cell_area(grid.cell(k)) * radius * radius)
        .collect();
    let surface: Vec<f64> = height.iter().map(|&h| h as f64).collect();
    let intake: Vec<f64> = ground.profile.iter().map(|p| p.intake()).collect();
    // What soaks in, m/s, and what runs off over the ground, kg/m²/s.
    let soak: Vec<f64> = (0..n)
        .map(|k| {
            if sea[k] {
                0.0
            } else {
                (runoff[k].max(0.0) as f64 / density).min(intake[k])
            }
        })
        .collect();
    let over: Vec<f64> = (0..n)
        .map(|k| {
            if sea[k] {
                0.0
            } else {
                (runoff[k].max(0.0) as f64 - soak[k] * density).max(0.0)
            }
        })
        .collect();
    let wanted: Vec<f64> = (0..n)
        .map(|k| {
            if sea[k] {
                0.0
            } else {
                demand[k].max(0.0) as f64 / density
            }
        })
        .collect();
    let start = start.filter(|s| s.head.len() == n && s.held.len() == n && s.focused.len() == n);
    let mut held: Vec<f64> = (0..n)
        .map(|k| {
            if sea[k] {
                ground.level
            } else {
                start.map_or(f64::NAN, |s| s.held[k])
            }
        })
        .collect();
    let mut focused = start.map_or_else(|| vec![0.0; n], |s| s.focused.clone());
    let mut head: Option<Vec<f64>> = start.map(|s| s.head.clone());
    let mut cover = vec![0.0f32; n];
    let land: f64 = (0..n)
        .filter(|&k| !sea[k])
        .map(|k| area[k])
        .sum::<f64>()
        .max(1e-30);
    let soaked: f64 = (0..n).map(|k| soak[k] * area[k]).sum::<f64>().max(1e-30);
    let mut result = None;
    for round in 0..GROUND_COUPLING {
        let recharge: Vec<f64> = (0..n).map(|k| soak[k] + focused[k]).collect();
        let table = table(
            (grid, ground.sides, radius),
            &surface,
            &held,
            (&recharge, &wanted),
            ground.profile,
            head.as_deref(),
        );
        // The surface's water: what runs off over the ground and what seeps out of it; what
        // the ground gives a cell's surface comes to it whether or not a lake floods it, past
        // what the lake's surface gains (`rivers::drain` holding a lake's cell to its surface's
        // gain where that is the less), and what a lake gives the ground leaves its surface; a
        // lake over ground whose table lies below it gives it what it takes in.
        let lake = |k: usize| !sea[k] && held[k].is_finite();
        let given = |k: usize| {
            if lake(k) {
                table.exchange[k]
            } else {
                table.seep[k]
            }
        };
        let run: Vec<f32> = (0..n)
            .map(|k| (over[k] + given(k).max(0.0) / area[k] * density) as f32)
            .collect();
        let surface_gain: Vec<f32> = (0..n)
            .map(|k| {
                let g = gain[k] as f64 + given(k) / area[k] * density;
                let g = if !lake(k) && !sea[k] && !table.wet[k] {
                    g - intake[k] * density
                } else {
                    g
                };
                g as f32
            })
            .collect();
        let drainage =
            crate::rivers::drain(grid, height, sea, (&run, &surface_gain), (radius, density));
        // The lakes as they stand, the table held at their levels; what a lake gives the
        // ground beneath it where the table lies below it soaks in there.
        let mut next: Vec<f64> = (0..n)
            .map(|k| if sea[k] { ground.level } else { f64::NAN })
            .collect();
        let mut leaked = vec![0.0; n];
        for k in 0..n {
            let l = drainage.lake_of[k];
            if l == crate::rivers::NONE || sea[k] {
                continue;
            }
            let hold = if lake(k) {
                GROUND_HOLD_LEAST
            } else {
                GROUND_HOLD_MOST
            };
            if drainage.cover[k] as f64 >= hold {
                let found = drainage.levels[l as usize] as f64;
                next[k] = if lake(k) {
                    held[k] + GROUND_LEVEL_STEP * (found - held[k])
                } else {
                    found
                };
            }
            if !lake(k) && !table.wet[k] {
                leaked[k] = drainage.cover[k] as f64 * intake[k];
            }
        }
        // Settled when the lakes cover the land as they did and what they give the ground
        // beneath them is what was given, each within GROUND_SETTLED_SHARE, and the levels
        // held both rounds moved less than GROUND_LEVEL_SETTLED.
        let moved: f64 = (0..n)
            .map(|k| (drainage.cover[k] - cover[k]).abs() as f64 * area[k])
            .sum::<f64>()
            / land;
        let given: f64 = (0..n)
            .map(|k| (leaked[k] - focused[k]).abs() * area[k])
            .sum::<f64>()
            / soaked;
        let rose = (0..n)
            .filter(|&k| !sea[k] && next[k].is_finite() && held[k].is_finite())
            .map(|k| (next[k] - held[k]).abs())
            .fold(0.0, f64::max);
        let settled = moved < GROUND_SETTLED_SHARE
            && given < GROUND_SETTLED_SHARE
            && rose < GROUND_LEVEL_SETTLED;
        cover.copy_from_slice(&drainage.cover);
        for k in 0..n {
            focused[k] = 0.5 * (focused[k] + leaked[k]);
        }
        let mut table = table;
        table.leak = (0..n)
            .map(|k| {
                if held[k].is_finite() {
                    0.0
                } else {
                    (recharge[k] - soak[k]) * area[k]
                }
            })
            .collect();
        table.brief = brief(&table, sea, &recharge, &area, round + 1);
        held = next;
        head = Some(table.head.clone());
        result = Some((drainage, table));
        if settled {
            break;
        }
    }
    let (drainage, table) = result.expect("at least one round");
    let end = Start {
        head: table.head.clone(),
        held,
        focused,
    };
    (drainage, table, end)
}

/// The ground's year in brief (F435) from its `table` on cells of `area`, m², the sea where
/// `sea`, under `recharge`, m/s, what soaked in from its surface and what the lakes over a
/// table below them give it (its `leak`); found with the lakes over `coupled` rounds.
fn brief(table: &Table, sea: &[bool], recharge: &[f64], area: &[f64], coupled: usize) -> Brief {
    let n = area.len();
    let land_of = |k: usize| !sea[k] && !table.held[k].is_finite();
    let sum = |f: &dyn Fn(usize) -> f64| (0..n).map(f).sum::<f64>();
    let land = sum(&|k| if land_of(k) { area[k] } else { 0.0 });
    Brief {
        soaked: sum(&|k| {
            if land_of(k) {
                recharge[k] * area[k]
            } else {
                0.0
            }
        }),
        dried: table.dried.iter().sum(),
        seeped: table.seep.iter().sum(),
        to_lakes: sum(&|k| {
            if !sea[k] && table.held[k].is_finite() {
                table.exchange[k]
            } else {
                0.0
            }
        }),
        to_sea: sum(&|k| if sea[k] { table.exchange[k] } else { 0.0 }),
        leaked: table.leak.iter().sum(),
        drawn: 0.0,
        wet: sum(&|k| {
            if land_of(k) && table.wet[k] {
                area[k]
            } else {
                0.0
            }
        }) / land.max(1e-30),
        coupled,
    }
}

/// The rivers and the lakes on a grid with the ground between them (F435), the table found on
/// a coarser grid, `coarse`, the ground's, of each of its cells' `profile`, the sea's `level`,
/// m, where `drain_through` finds it over the ground, the runoff and the gain brought to it by
/// their means over its cells, a coarse cell of land at its land's mean height; the rivers and
/// the lakes on `grid`'s own, the ground giving each cell's surface or taking from it what it
/// gives or takes from its coarse cell's, spread over its cells by their area, whether a lake
/// covers them or not: where the coarse cell is held at the sea's or a lake's level, its share
/// of what its lake gives the ground or the ground gives it, a cell's land giving at most what
/// it sheds; elsewhere what each cell sheds soaking in as far as its coarse cell's ground takes
/// it in, and its share of what seeps out; and the lakes over the cells above the table, where
/// the coarse cells' own lakes stood over it, giving their ground in all what those gave it
/// (F435), each at the same share of what its ground takes in: at first all of it; where they
/// give more, the share at which they would give no more over all the land above the table;
/// then between the two shares at which they gave less and more (a lake the wider for giving
/// less), by regula falsi on their logarithms, until what they give it is the coarse grid's
/// lakes' within GROUND_SETTLED_SHARE of what soaks in, or GROUND_COUPLING passes have run.
/// Its table, on `coarse`.
pub fn route_through(
    grid: CubeSphere,
    height: &[f32],
    sea: &[bool],
    (runoff, gain): (&[f32], &[f32]),
    (radius, density): (f64, f64),
    (coarse, profile, level): (CubeSphere, &[Profile], f64),
) -> (crate::rivers::Rivers, Table) {
    let n = grid.cell_count();
    let to: Vec<usize> = (0..n)
        .into_par_iter()
        .map(|k| coarse.index(coarse.locate(grid.center(grid.cell(k))).0))
        .collect();
    let mean = |values: &dyn Fn(usize) -> f64| -> Vec<f64> {
        let m = coarse.cell_count();
        let (mut sum, mut area) = (vec![0.0; m], vec![0.0; m]);
        for k in 0..n {
            let a = grid.cell_area(grid.cell(k));
            sum[to[k]] += values(k) * a;
            area[to[k]] += a;
        }
        (0..m)
            .map(|c| if area[c] > 0.0 { sum[c] / area[c] } else { 0.0 })
            .collect()
    };
    let coarse_sea: Vec<bool> = mean(&|k| if sea[k] { 1.0 } else { 0.0 })
        .iter()
        .map(|&s| s >= 0.5)
        .collect();
    // A coarse cell of land stands at its land's mean height, not below the sea's level by
    // the sea's floor beside it.
    let (all, land) = (
        mean(&|k| height[k] as f64),
        mean(&|k| if sea[k] { 0.0 } else { height[k] as f64 }),
    );
    let land_share = mean(&|k| if sea[k] { 0.0 } else { 1.0 });
    let coarse_height: Vec<f32> = (0..coarse.cell_count())
        .map(|c| {
            if coarse_sea[c] || land_share[c] <= 0.0 {
                all[c] as f32
            } else {
                (land[c] / land_share[c]) as f32
            }
        })
        .collect();
    let coarse_runoff: Vec<f32> = mean(&|k| if sea[k] { 0.0 } else { runoff[k] as f64 })
        .iter()
        .map(|&r| r as f32)
        .collect();
    let coarse_gain: Vec<f32> = mean(&|k| gain[k] as f64)
        .iter()
        .map(|&g| g as f32)
        .collect();
    let sides = Sides::new(coarse);
    let ground = Ground {
        sides: &sides,
        profile,
        level,
    };
    let m = coarse.cell_count();
    let none = vec![0.0f32; m];
    let (_, first, _) = drain_through(
        coarse,
        &coarse_height,
        &coarse_sea,
        (&coarse_runoff, &coarse_gain, &none),
        (radius, density),
        &ground,
        None,
    );
    let r2 = radius * radius;
    let coarse_area: Vec<f64> = (0..m)
        .map(|c| coarse.cell_area(coarse.cell(c)) * r2)
        .collect();
    // Each coarse cell's land, m².
    let coarse_land: Vec<f64> = (0..m).map(|c| coarse_area[c] * land_share[c]).collect();
    let fine_area: Vec<f64> = (0..n)
        .into_par_iter()
        .map(|k| grid.cell_area(grid.cell(k)) * r2)
        .collect();
    let mut table = first;
    // What the ground gives each coarse cell's land, kg/m²/s, less what its lake gives the
    // ground: none where the sea holds it.
    let given = |c: usize| {
        if coarse_sea[c] || coarse_land[c] <= 0.0 {
            0.0
        } else {
            (table.exchange[c] + table.seep[c]) / coarse_land[c] * density
        }
    };
    // Each cell's surface water over its coarse cell's ground: where the coarse cell is held,
    // its share of what the ground gives the cell or takes from it, at most what it sheds;
    // elsewhere what soaks in as far as the ground takes it in, and its share of what seeps
    // out.
    let run: Vec<f32> = (0..n)
        .map(|k| {
            if sea[k] {
                return 0.0;
            }
            let c = to[k];
            let shed = runoff[k].max(0.0) as f64;
            let left = if table.held[c].is_finite() {
                shed + given(c)
            } else {
                shed - (shed / density).min(profile[c].intake()) * density + given(c)
            };
            left.max(0.0) as f32
        })
        .collect();
    // What each cell's surface gains as a lake: what its surface gains and its share of what
    // the ground gives its coarse cell or takes from it.
    let open: Vec<f64> = (0..n).map(|k| gain[k] as f64 + given(to[k])).collect();
    // The cells above the table where the coarse cell's own lakes stood over it, and what the
    // ground takes in there, kg/m²/s.
    let above: Vec<bool> = (0..n)
        .map(|k| {
            let c = to[k];
            !sea[k] && table.leak[c] > 0.0 && height[k] as f64 > table.head[c]
        })
        .collect();
    let intake: Vec<f64> = (0..m).map(|c| profile[c].intake() * density).collect();
    // What the ground's grid's lakes gave it, m³/s, and what the lakes over all the land above
    // the table would give it at what it takes in: at that share of it they give no more.
    let own = table.brief.leaked;
    let reach: f64 = (0..n)
        .filter(|&k| above[k])
        .map(|k| intake[to[k]] * fine_area[k] / density)
        .sum();
    // The share of what the ground takes in, by its logarithm, and how much more or less the
    // lakes gave it than its grid's did, by its logarithm, where they gave less and more; the
    // side found last.
    let mut share = 1.0f64;
    let mut less: Option<(f64, f64)> = None;
    let mut more: Option<(f64, f64)> = None;
    let mut side = 0i8;
    let mut found = None;
    for _ in 0..GROUND_COUPLING {
        let leak = |k: usize| {
            if above[k] { share * intake[to[k]] } else { 0.0 }
        };
        let surface_gain: Vec<f32> = (0..n).map(|k| (open[k] - leak(k)) as f32).collect();
        let drainage =
            crate::rivers::drain(grid, height, sea, (&run, &surface_gain), (radius, density));
        // What the lakes give their ground, m³/s: what each cell a lake covers gives its lake the
        // less for it (`rivers::drain` holding a lake's cell to its surface's gain where that is
        // the less than what it sheds).
        let drawn: f64 = (0..n)
            .filter(|&k| above[k] && drainage.cover[k] > 0.0)
            .map(|k| {
                let a = fine_area[k];
                let shed = run[k] as f64 * a;
                let without = (open[k] * a - shed).min(0.0);
                let with = ((open[k] - leak(k)) * a - shed).min(0.0);
                drainage.cover[k] as f64 * (without - with) / density
            })
            .sum();
        found = Some((drainage, drawn));
        // Settled where what they give it is its grid's lakes' within GROUND_SETTLED_SHARE of
        // what soaks in, or where they give it less at what it takes in.
        if (drawn - own).abs() < GROUND_SETTLED_SHARE * table.brief.soaked
            || (drawn < own && share >= 1.0)
        {
            break;
        }
        // The next share: where they gave more, the least first, then between the shares at
        // which they gave less and more, where the line between them by their logarithms meets
        // its grid's (regula falsi, the Illinois way).
        let at = (share.ln(), (drawn.max(1e-9 * own) / own).ln());
        if drawn > own {
            more = Some(at);
            if side > 0
                && let Some((x, y)) = less
            {
                less = Some((x, 0.5 * y));
            }
            side = 1;
        } else {
            less = Some(at);
            if side < 0
                && let Some((x, y)) = more
            {
                more = Some((x, 0.5 * y));
            }
            side = -1;
        }
        share = match (less, more) {
            (Some((xl, yl)), Some((xm, ym))) => (xl - yl * (xm - xl) / (ym - yl)).exp(),
            _ => (own / reach).min(1.0),
        };
    }
    let (drainage, drawn) = found.expect("at least one pass");
    table.brief.drawn = drawn;
    let rivers = crate::rivers::from_drainage(grid, sea, &drainage, density);
    (rivers, table)
}

/// Each side's conductance, m²/s: the mean of its two cells' transmissivities below their
/// tables times the side's shape.
fn conduct(sides: &Sides, profile: &[Profile], (surface, head): (&[f64], &[f64]), out: &mut [f64]) {
    let t: Vec<f64> = (0..head.len())
        .into_par_iter()
        .map(|k| profile[k].transmissivity(surface[k] - head[k]))
        .collect();
    out.par_iter_mut().enumerate().for_each(|(s, c)| {
        let (k, m) = (sides.from[s] as usize, sides.to[s] as usize);
        *c = 0.5 * (t[k] + t[m]) * sides.shape[s];
    });
}

/// What flows into each cell, m³/s: what soaks in, and through its sides from its
/// neighbours, sum c (h_m - h_k).
fn net(sides: &Sides, conductance: &[f64], head: &[f64], source: &[f64], out: &mut [f64]) {
    out.par_iter_mut().enumerate().for_each(|(k, o)| {
        let (a, b) = (sides.start[k] as usize, sides.start[k + 1] as usize);
        let mut s = source[k];
        for i in a..b {
            s += conductance[i] * (head[sides.to[i] as usize] - head[k]);
        }
        *o = s;
    });
}

/// The heads' change on the `free` cells that brings their inflow `residual` to none, the
/// others held: sum over sides c (x_k - x_m [m free]) + e_k x_k = r_k, e the outflow's rise
/// with the head, by conjugate gradients preconditioned by the diagonal; and its steps.
fn solve(
    sides: &Sides,
    (conductance, rise): (&[f64], &[f64]),
    free: &[bool],
    residual: &[f64],
    tolerance: f64,
) -> (Vec<f64>, usize) {
    let n = free.len();
    let diagonal: Vec<f64> = (0..n)
        .into_par_iter()
        .map(|k| {
            if free[k] {
                let (a, b) = (sides.start[k] as usize, sides.start[k + 1] as usize);
                (conductance[a..b].iter().sum::<f64>() + rise[k]).max(1e-300)
            } else {
                1.0
            }
        })
        .collect();
    let apply = |x: &[f64], y: &mut [f64]| {
        y.par_iter_mut().enumerate().for_each(|(k, yk)| {
            *yk = if free[k] {
                let (a, b) = (sides.start[k] as usize, sides.start[k + 1] as usize);
                let mut s = rise[k] * x[k];
                for (&m, &c) in sides.to[a..b].iter().zip(&conductance[a..b]) {
                    let m = m as usize;
                    let other = if free[m] { x[m] } else { 0.0 };
                    s += c * (x[k] - other);
                }
                s
            } else {
                0.0
            };
        });
    };
    let mut x = vec![0.0; n];
    let mut r: Vec<f64> = (0..n)
        .map(|k| if free[k] { residual[k] } else { 0.0 })
        .collect();
    let mut z: Vec<f64> = (0..n).map(|k| r[k] / diagonal[k]).collect();
    let mut p = z.clone();
    let mut ap = vec![0.0; n];
    let mut rz = det_sum(n, |k| r[k] * z[k]);
    let mut steps = 0;
    for _ in 0..GROUND_STEPS.max(n) {
        if det_sum(n, |k| r[k].abs()) <= tolerance {
            break;
        }
        steps += 1;
        apply(&p, &mut ap);
        let pap = det_sum(n, |k| p[k] * ap[k]);
        if pap <= 0.0 {
            break;
        }
        let alpha = rz / pap;
        x.par_iter_mut()
            .zip(r.par_iter_mut())
            .enumerate()
            .for_each(|(k, (xk, rk))| {
                *xk += alpha * p[k];
                *rk -= alpha * ap[k];
            });
        z.par_iter_mut()
            .enumerate()
            .for_each(|(k, zk)| *zk = r[k] / diagonal[k]);
        let next = det_sum(n, |k| r[k] * z[k]);
        let beta = next / rz;
        rz = next;
        p.par_iter_mut()
            .enumerate()
            .for_each(|(k, pk)| *pk = z[k] + beta * *pk);
    }
    (x, steps)
}

/// The most rounds of the table found with the last round's lakes in a call, and of the
/// finish's lakes found with what they give the ground spread over them (F435), a numerical
/// bound: the year's calls go on from where the last left off.
const GROUND_COUPLING: usize = 8;
/// The lakes settled when their cover of the land moves less than this share of it in a
/// round, what they give the ground less than this share of what soaks in, and no level held
/// both rounds moves more than GROUND_LEVEL_SETTLED, m; the finish's when what they give the
/// ground is off its grid's by less than this share of what soaks in (F435), numerical
/// tolerances.
const GROUND_SETTLED_SHARE: f64 = 1e-4;
const GROUND_LEVEL_SETTLED: f64 = 0.01;
/// The share of the way a lake's level held the round before moves to the one found (F435),
/// a numerical damping.
const GROUND_LEVEL_STEP: f64 = 0.5;
/// The share of a cell a lake covers above which the table is held at its level, a cell held
/// the round before, and one not (F435), numerical bounds against a lake's edge turning over.
const GROUND_HOLD_LEAST: f64 = 0.25;
const GROUND_HOLD_MOST: f64 = 0.75;
/// The most passes of the lagged transmissivity, and rounds of the active set within each
/// (F435), numerical bounds.
const GROUND_PASSES: usize = 40;
const GROUND_ROUNDS: usize = 50;
/// A cell at the surface leaves it where it would draw from it more than this share of what
/// soaks into a cell in the mean (F435), a numerical tolerance against cells that would draw
/// next to nothing turning over and back.
const GROUND_DRAW_SHARE: f64 = 1e-6;
/// The table settled when no cell moves more than this through a pass, m (F435), a numerical
/// tolerance.
const GROUND_SETTLED: f64 = 1e-3;
/// The conjugate gradients' fewest steps allowed before they stop short of the tolerance,
/// at least the cells' count (F435), a numerical bound.
const GROUND_STEPS: usize = 20_000;

#[cfg(test)]
mod tests {
    use super::*;

    /// A polar continent of colatitude 60 degrees on a body of 1,000 km, a plateau `plateau`
    /// m high falling to its coast over its last 20 degrees, a bowl `bowl` m deep within 12
    /// degrees of its pole; the sea beyond.
    fn continent(grid: CubeSphere, plateau: f64, bowl: f64) -> (Vec<f32>, Vec<bool>) {
        let top = [0.0, 0.0, 1.0];
        let heights: Vec<f32> = (0..grid.cell_count())
            .map(|k| {
                let th = dot(grid.center(grid.cell(k)), top)
                    .clamp(-1.0, 1.0)
                    .acos()
                    .to_degrees();
                let h = if th > 60.0 {
                    -100.0
                } else if th > 40.0 {
                    plateau * (60.0 - th) / 20.0
                } else {
                    plateau
                };
                let dip = if th < 12.0 {
                    bowl * (0.5 + 0.5 * (std::f64::consts::PI * th / 12.0).cos())
                } else {
                    0.0
                };
                (h - dip) as f32
            })
            .collect();
        let sea = heights.iter().map(|&h| h < 0.0).collect();
        (heights, sea)
    }

    /// A hollow 300 m deep within 4 degrees of colatitude 50 on the meridian, its floor below
    /// the sea's level, the coast's slope about it; its centre.
    fn dig_hollow(grid: CubeSphere, heights: &mut [f32], sea: &[bool]) -> [f64; 3] {
        let at = crate::grid::direction_from_lon_lat(0.0, 40f64.to_radians());
        for (k, h) in heights.iter_mut().enumerate() {
            let a = dot(grid.center(grid.cell(k)), at)
                .clamp(-1.0, 1.0)
                .acos()
                .to_degrees();
            if a < 4.0 && !sea[k] {
                *h -= (300.0 * (0.5 + 0.5 * (std::f64::consts::PI * a / 4.0).cos())) as f32;
            }
        }
        at
    }

    /// F435: over a ground that carries what soaks into it to the sea, a kilometre of gravel
    /// under 9.5 mm a year, the table lies below a bowl's floor on the plateau and the bowl
    /// holds no lake, what the land sheds reaching the sea through the ground; over a ground
    /// that carries next to none, the table at the surface, the bowl holds the lake the rivers
    /// alone would give it, within 2 percent; and over a crust that takes in two thirds of what
    /// the land sheds, over rock that carries it away, what runs off its crust into the bowl
    /// standing in a lake perched over the table that gives its ground what the crust takes in
    /// as well as its air what it would take, less than the rivers alone would give it and
    /// within 2 percent of the lake they give it from what runs off the crust, its surface
    /// losing what the crust takes in as well.
    #[test]
    fn a_bowl_above_its_table_holds_no_lake() {
        let grid = CubeSphere::new(48);
        let sides = Sides::new(grid);
        let (radius, density) = (1e6, 1000.0);
        let (heights, sea) = continent(grid, 500.0, 300.0);
        let n = grid.cell_count();
        let runoff = vec![3e-7f32; n];
        let gain = vec![-1e-7f32; n];
        let plain = crate::rivers::drain(grid, &heights, &sea, (&runoff, &gain), (radius, density));
        for (name, k, b) in [
            ("gravel", 1e-3, 1000.0),
            ("tight", 1e-13, 0.0),
            ("crusted", 2e-10, 1.0),
        ] {
            let rock = if name == "crusted" { 1e-3 } else { k };
            let profile = vec![
                Profile {
                    cover: b,
                    cover_conductivity: k,
                    rock_conductivity: rock,
                    knee: 1000.0,
                };
                n
            ];
            let ground = Ground {
                sides: &sides,
                profile: &profile,
                level: 0.0,
            };
            let (d, table, _) = drain_through(
                grid,
                &heights,
                &sea,
                (&runoff, &gain, &vec![0.0f32; n]),
                (radius, density),
                &ground,
                None,
            );
            let area = |k: usize| grid.cell_area(grid.cell(k)) * radius * radius;
            let shed: f64 = (0..n)
                .filter(|&k| !sea[k])
                .map(|k| 3e-7 * area(k) / density)
                .sum();
            let through: f64 = (0..n).filter(|&k| sea[k]).map(|k| table.exchange[k]).sum();
            let pole = (0..n)
                .max_by(|&a, &b| {
                    grid.center(grid.cell(a))[2].total_cmp(&grid.center(grid.cell(b))[2])
                })
                .unwrap();
            println!(
                "{name}: lakes {:.4} of the land (the rivers alone {:.4}); the table {:.1} m at the pole under its floor at {:.1}; {:.4e} m³/s shed, {:.4e} through the ground to the sea, {:.4e} by the rivers",
                d.summary.lakes,
                plain.summary.lakes,
                table.head[pole],
                heights[pole],
                shed,
                through,
                d.summary.to_sea
            );
            if name == "gravel" {
                assert!(d.summary.lakes < 1e-3, "{name}");
                assert!(table.head[pole] < heights[pole] as f64 - 50.0, "{name}");
                assert!(((through - shed) / shed).abs() < 1e-3, "{name}");
            } else if name == "crusted" {
                // The rivers alone over what runs off the crust, the lake's surface losing what
                // the crust takes in as well as what the air takes.
                let perched = crate::rivers::drain(
                    grid,
                    &heights,
                    &sea,
                    (&vec![1e-7f32; n], &vec![-3e-7f32; n]),
                    (radius, density),
                );
                println!(
                    "{name}: the rivers over what runs off the crust, its lake losing what the crust takes in: lakes {:.4}",
                    perched.summary.lakes
                );
                assert!(d.summary.lakes < plain.summary.lakes, "{name}");
                assert!(
                    ((d.summary.lakes - perched.summary.lakes) / perched.summary.lakes).abs()
                        < 0.02,
                    "{name}"
                );
            } else {
                assert!(
                    ((d.summary.lakes - plain.summary.lakes) / plain.summary.lakes).abs() < 0.02,
                    "{name}"
                );
            }
        }
    }

    /// F435: the finish's rivers run over the ground's table, found on a grid half as fine as
    /// theirs: over a kilometre of gravel the bowl on the plateau holds no lake, what the land
    /// sheds reaching the sea through the ground or, from the cells by the coast the sea holds
    /// on the ground's grid, by the rivers, within a thousandth, and its table below all its
    /// land, the sea feeding no seep at the coast; over a ground that carries next to none, the
    /// bowl holds the lake the rivers alone give it, within 2 percent.
    #[test]
    fn the_finish_s_rivers_run_over_the_ground_s_table() {
        let (grid, coarse) = (CubeSphere::new(96), CubeSphere::new(48));
        let (radius, density) = (1e6, 1000.0);
        let (heights, sea) = continent(grid, 500.0, 300.0);
        let n = grid.cell_count();
        let runoff = vec![3e-7f32; n];
        let gain = vec![-1e-7f32; n];
        let plain = crate::rivers::route(grid, &heights, &sea, (&runoff, &gain), (radius, density));
        for (name, k, b) in [("gravel", 1e-3, 1000.0), ("tight", 1e-13, 0.0)] {
            let profile = vec![
                Profile {
                    cover: b,
                    cover_conductivity: k,
                    rock_conductivity: k,
                    knee: 1000.0,
                };
                coarse.cell_count()
            ];
            let (rivers, table) = route_through(
                grid,
                &heights,
                &sea,
                (&runoff, &gain),
                (radius, density),
                (coarse, &profile, 0.0),
            );
            let (lakes, alone) = (rivers.summary.lakes, plain.summary.lakes);
            let given: f64 = (0..n)
                .filter(|&k| !sea[k])
                .map(|k| 3e-7 / density * grid.cell_area(grid.cell(k)) * radius * radius)
                .sum();
            let reached = rivers.summary.to_sea + table.brief.to_sea;
            println!(
                "{name}: lakes {lakes:.4} of the land (the rivers alone {alone:.4}); given {given:.2} m³/s, {:.2} by the rivers and {:.2} through the ground, {:.4} seeping out",
                rivers.summary.to_sea, table.brief.to_sea, table.brief.seeped
            );
            assert!(alone > 0.0);
            if name == "gravel" {
                assert!(lakes < 1e-3, "{name}");
                assert!((reached / given - 1.0).abs() < 1e-3, "{name}");
                assert!(table.brief.seeped < 1e-3 * given, "{name}");
            } else {
                assert!(((lakes - alone) / alone).abs() < 0.02, "{name}");
            }
        }
    }

    /// F435: the finish's rivers and ground keep the land's water. Over a ground that takes in
    /// two thirds of what the land sheds, its table deep under a plateau pitted below its
    /// ground's cells and a bowl: the lakes' surfaces gaining what the land sheds, nothing lost
    /// to the air, all the land sheds reaches the sea through the rivers or the ground within a
    /// hundredth; their surfaces gaining nothing, the lakes give the ground no more than its own
    /// lakes give it, within a hundredth, the pits', which stand on the fine grid alone, giving
    /// it nothing; and a hollow by the coast below the table, over a kilometre of gravel, its
    /// lake fed by the ground, all the land sheds reaching the sea within a hundredth, its lakes
    /// on the rivers' grid giving the ground what its own grid's gave it within a hundredth.
    #[test]
    fn the_finish_keeps_the_land_s_water() {
        let (grid, coarse) = (CubeSphere::new(96), CubeSphere::new(24));
        let (radius, density) = (1e6, 1000.0);
        let (mut heights, sea) = continent(grid, 500.0, 300.0);
        let n = grid.cell_count();
        let mut pits = 0;
        for (k, h) in heights.iter_mut().enumerate() {
            let c = grid.cell(k);
            if *h == 500.0 && (1..3).contains(&(c.i % 4)) && (1..3).contains(&(c.j % 4)) {
                *h -= 30.0;
                pits += 1;
            }
        }
        let runoff = vec![3e-7f32; n];
        let profile = vec![
            Profile {
                cover: 1.0,
                cover_conductivity: 2e-10,
                rock_conductivity: 1e-3,
                knee: 1000.0,
            };
            coarse.cell_count()
        ];
        let given: f64 = (0..n)
            .filter(|&k| !sea[k])
            .map(|k| runoff[k] as f64 / density * grid.cell_area(grid.cell(k)) * radius * radius)
            .sum();
        for (name, g) in [
            ("gaining what the land sheds", 3e-7f32),
            ("gaining nothing", 0.0),
        ] {
            let gain = vec![g; n];
            let (rivers, table) = route_through(
                grid,
                &heights,
                &sea,
                (&runoff, &gain),
                (radius, density),
                (coarse, &profile, 0.0),
            );
            let s = rivers.summary;
            let reached = s.to_sea + table.brief.to_sea;
            println!(
                "lakes {name}: {pits} pits; given {given:.1} m³/s, {:.1} by the rivers and {:.1} through the ground: {:+.4} of it; the lakes' surfaces {:+.1}; lakes {:.4} of the land, the table at the surface over {:.3} of it; {:?}",
                s.to_sea,
                table.brief.to_sea,
                reached / given - 1.0,
                s.lakes_gain,
                s.lakes,
                table.brief.wet,
                table.brief
            );
            assert!(pits > 0 && s.lakes > 0.0 && table.brief.wet < 0.5);
            if g > 0.0 {
                let b = table.brief;
                assert!(
                    ((reached + b.drawn - b.leaked) / given - 1.0).abs() < 0.01,
                    "{name}"
                );
            } else {
                let b = table.brief;
                assert!(-s.lakes_gain <= 1.01 * (b.leaked - b.to_lakes), "{name}");
            }
        }
        // The bowl on the plateau over a ground whose table rises through its floor, found on a
        // grid half as fine: its lake fed by the ground.
        let coarse = CubeSphere::new(48);
        let (heights, sea) = continent(grid, 500.0, 300.0);
        let fed = vec![
            Profile {
                cover: 1000.0,
                cover_conductivity: 1e-4,
                rock_conductivity: 1e-4,
                knee: 100.0,
            };
            coarse.cell_count()
        ];
        let (rivers, table) = route_through(
            grid,
            &heights,
            &sea,
            (&runoff, &runoff),
            (radius, density),
            (coarse, &fed, 0.0),
        );
        let s = rivers.summary;
        let b = table.brief;
        let reached = s.to_sea + b.to_sea;
        println!(
            "the bowl fed by its ground: given {given:.1} m³/s, {:.1} by the rivers and {:.1} through the ground, the lakes giving it {:.1} on the rivers' grid where {:.1} on its own: {:+.4} of it; {:.1} seeps out, the ground gives the lakes {:+.1}; lakes {:.4} of the land",
            s.to_sea,
            b.to_sea,
            b.drawn,
            b.leaked,
            (reached + b.drawn - b.leaked) / given - 1.0,
            b.seeped,
            b.to_lakes,
            s.lakes
        );
        assert!(s.lakes > 0.0 && b.seeped + b.to_lakes > 0.1 * given);
        assert!(b.leaked > 0.0 && (b.drawn / b.leaked - 1.0).abs() < 0.01);
        assert!(((reached + b.drawn - b.leaked) / given - 1.0).abs() < 0.01);
    }

    /// The coupled table's cost on finer grids (a probe, not a check).
    #[test]
    #[ignore]
    fn the_table_s_cost() {
        for face in [64u32, 128] {
            let grid = CubeSphere::new(face);
            let t0 = std::time::Instant::now();
            let sides = Sides::new(grid);
            let built = t0.elapsed().as_secs_f64();
            let (radius, density) = (1e6, 1000.0);
            let (heights, sea) = continent(grid, 500.0, 300.0);
            let n = grid.cell_count();
            let runoff = vec![3e-7f32; n];
            let gain = vec![-1e-7f32; n];
            for (name, k, b) in [("gravel", 1e-3, 1000.0), ("sand", 1e-5, 100.0)] {
                let profile = vec![
                    Profile {
                        cover: b,
                        cover_conductivity: k,
                        rock_conductivity: k,
                        knee: 100.0,
                    };
                    n
                ];
                let ground = Ground {
                    sides: &sides,
                    profile: &profile,
                    level: 0.0,
                };
                let t1 = std::time::Instant::now();
                let (d, table, _) = drain_through(
                    grid,
                    &heights,
                    &sea,
                    (&runoff, &gain, &vec![0.0f32; n]),
                    (radius, density),
                    &ground,
                    None,
                );
                println!(
                    "face {face} {name}: sides {built:.2} s, drained {:.2} s; lakes {:.4}; last table {} rounds, {} steps",
                    t1.elapsed().as_secs_f64(),
                    d.summary.lakes,
                    table.rounds,
                    table.steps
                );
            }
        }
    }

    /// F435: the less the ground carries, the more of the land its table meets the surface
    /// at, from a few cells by the coast over a kilometre of gravel, under 2 percent of the
    /// land, to nine tenths and more over a ground that carries next to none, each
    /// conductivity's no less than the next one's.
    #[test]
    fn the_table_follows_the_ground_where_the_ground_carries_little() {
        let grid = CubeSphere::new(32);
        let sides = Sides::new(grid);
        let (radius, density) = (1e6, 1000.0);
        let (heights, sea) = continent(grid, 500.0, 0.0);
        let n = grid.cell_count();
        let runoff = vec![3e-7f32; n];
        let gain = vec![-1e-7f32; n];
        let mut last = f64::INFINITY;
        for k in [1e-11, 1e-9, 1e-7, 1e-5, 1e-4, 1e-3] {
            let profile = vec![
                Profile {
                    cover: 1000.0,
                    cover_conductivity: k,
                    rock_conductivity: k,
                    knee: 100.0,
                };
                n
            ];
            let ground = Ground {
                sides: &sides,
                profile: &profile,
                level: 0.0,
            };
            let (_, table, _) = drain_through(
                grid,
                &heights,
                &sea,
                (&runoff, &gain, &vec![0.0f32; n]),
                (radius, density),
                &ground,
                None,
            );
            let land = (0..n).filter(|&c| !sea[c]).count() as f64;
            let wet = (0..n).filter(|&c| !sea[c] && table.wet[c]).count() as f64 / land;
            println!("K {k:e} m/s: the table at the surface over {wet:.3} of the land");
            assert!(wet <= last, "{k:e}");
            if k <= 1e-11 {
                assert!(wet >= 0.9);
            }
            if k >= 1e-3 {
                assert!(wet < 0.02);
            }
            last = wet;
        }
    }

    /// F435: over a kilometre of gravel, a hollow by the coast whose floor lies below the
    /// table holds a lake the table feeds, its level no higher than the table about it and over
    /// halfway up to it from its floor, while the bowl on the plateau above the table holds
    /// none.
    #[test]
    fn a_lake_stands_where_the_table_meets_its_floor() {
        let grid = CubeSphere::new(48);
        let sides = Sides::new(grid);
        let (radius, density) = (1e6, 1000.0);
        let (mut heights, sea) = continent(grid, 500.0, 300.0);
        let n = grid.cell_count();
        let at = dig_hollow(grid, &mut heights, &sea);
        let runoff = vec![3e-7f32; n];
        let gain = vec![-1e-7f32; n];
        let profile = vec![
            Profile {
                cover: 1000.0,
                cover_conductivity: 1e-3,
                rock_conductivity: 1e-3,
                knee: 100.0,
            };
            n
        ];
        let ground = Ground {
            sides: &sides,
            profile: &profile,
            level: 0.0,
        };
        let (d, table, _) = drain_through(
            grid,
            &heights,
            &sea,
            (&runoff, &gain, &vec![0.0f32; n]),
            (radius, density),
            &ground,
            None,
        );
        let held: Vec<usize> = (0..n)
            .filter(|&k| !sea[k] && d.lake_of[k] != crate::rivers::NONE && d.cover[k] >= 0.5)
            .collect();
        let near: Vec<usize> = (0..n)
            .filter(|&k| {
                let a = dot(grid.center(grid.cell(k)), at)
                    .clamp(-1.0, 1.0)
                    .acos()
                    .to_degrees();
                held.iter().all(|&h| h != k) && (4.0..7.0).contains(&a) && !sea[k]
            })
            .collect();
        let level = held
            .first()
            .map(|&k| d.levels[d.lake_of[k] as usize] as f64);
        let about = near.iter().map(|&k| table.head[k]).sum::<f64>() / near.len() as f64;
        let pole = (0..n)
            .max_by(|&a, &b| grid.center(grid.cell(a))[2].total_cmp(&grid.center(grid.cell(b))[2]))
            .unwrap();
        println!(
            "the hollow's lake over {} cells at {level:?} m, the table about it {about:.1} m; the plateau's bowl's floor {:.1} m over its table {:.1}, lakes {:.4} of the land",
            held.len(),
            heights[pole],
            table.head[pole],
            d.summary.lakes
        );
        let level = level.expect("a lake in the hollow");
        let floor = held
            .iter()
            .map(|&k| heights[k] as f64)
            .fold(f64::INFINITY, f64::min);
        assert!(
            level <= about + 1.0 && level >= 0.5 * (floor + about),
            "{floor} {level} {about}"
        );
        assert!(d.lake_of[pole] == crate::rivers::NONE || d.cover[pole] < 0.01);
    }

    /// F435: each liquid's viscosity within its fit's stated error of the NIST WebBook's
    /// saturated liquid at the ends of its range.
    #[test]
    fn each_liquid_flows_as_nist_measures_it() {
        let measured = [
            ("H2O", 273.16, 1.7914e-3, 0.005),
            ("H2O", 373.16, 2.8155e-4, 0.005),
            ("CH4", 91.0, 1.9176e-4, 0.003),
            ("CH4", 111.0, 1.1840e-4, 0.003),
            ("C2H6", 91.0, 1.2580e-3, 0.023),
            ("C2H6", 186.0, 1.6390e-4, 0.023),
            ("NH3", 196.0, 5.6403e-4, 0.003),
            ("NH3", 240.0, 2.5104e-4, 0.003),
            ("N2", 64.0, 2.9750e-4, 0.004),
            ("N2", 77.0, 1.6294e-4, 0.004),
            ("SO2", 200.0, 9.6157e-4, 0.011),
            ("SO2", 265.0, 4.0616e-4, 0.011),
            ("CO2", 217.0, 2.5163e-4, 0.034),
        ];
        for (name, t, mu, within) in measured {
            let liquid = VOLATILES.iter().find(|v| v.name == name).unwrap();
            let fit = viscosity(liquid, t);
            println!("{name} at {t} K: {fit:.4e} Pa s, NIST {mu:.4e}");
            assert!((fit / mu - 1.0).abs() <= within, "{name} {t}");
        }
    }

    /// A plain like Titan's on a grid: some 50 m of relief from waves of random directions, its
    /// lowest few percent the liquid, on a body of Titan's radius.
    fn plain(grid: CubeSphere) -> (Vec<f32>, Vec<bool>) {
        let rng = crate::rng::Rng::new(7);
        let waves: Vec<([f64; 3], f64, f64)> = (0..40u64)
            .map(|i| {
                let d = rng.unit_direction(i, 0);
                (
                    d,
                    2.0 + i as f64 * 1.5,
                    rng.range(i, 1, 0.0, std::f64::consts::TAU),
                )
            })
            .collect();
        let heights: Vec<f32> = (0..grid.cell_count())
            .map(|k| {
                let c = grid.center(grid.cell(k));
                let h: f64 = waves
                    .iter()
                    .map(|&(d, f, ph)| 40.0 / f.sqrt() * (f * dot(c, d) + ph).sin())
                    .sum();
                h as f32
            })
            .collect();
        let mut sorted = heights.clone();
        sorted.sort_by(f32::total_cmp);
        let low = sorted[sorted.len() / 40];
        let sea = heights.iter().map(|&h| h < low).collect();
        (heights, sea)
    }

    /// The table on a plain like Titan's (a probe, not a check).
    #[test]
    #[ignore]
    fn the_table_on_a_plain() {
        let grid = CubeSphere::new(128);
        let sides = Sides::new(grid);
        let (radius, density) = (2.575e6, 423.0);
        let (heights, sea) = plain(grid);
        let n = grid.cell_count();
        let runoff = vec![3.8e-8f32; n];
        let gain = vec![-2e-7f32; n];
        let level = (0..n)
            .filter(|&k| sea[k])
            .map(|k| heights[k])
            .fold(f32::MIN, f32::max) as f64;
        for (name, k) in [("ice", 2.4e-8), ("loose", 3.1e-7), ("gravel", 3e-4)] {
            let profile = vec![
                Profile {
                    cover: 0.0,
                    cover_conductivity: k,
                    rock_conductivity: k,
                    knee: 100.0,
                };
                n
            ];
            let ground = Ground {
                sides: &sides,
                profile: &profile,
                level,
            };
            let t0 = std::time::Instant::now();
            let (d, table, _) = drain_through(
                grid,
                &heights,
                &sea,
                (&runoff, &gain, &vec![0.0f32; n]),
                (radius, density),
                &ground,
                None,
            );
            let plain =
                crate::rivers::drain(grid, &heights, &sea, (&runoff, &gain), (radius, density));
            let land = (0..n).filter(|&c| !sea[c]).count() as f64;
            let wet = (0..n).filter(|&c| !sea[c] && table.wet[c]).count() as f64 / land;
            println!(
                "{name}: {:.2} s; lakes {:.4} of the land where the rivers alone {:.4}; the table at the surface over {wet:.3}; last table {} rounds {} steps",
                t0.elapsed().as_secs_f64(),
                d.summary.lakes,
                plain.summary.lakes,
                table.rounds,
                table.steps
            );
        }
    }

    /// F436: the ground gives the air from its table near the surface: over a ground that
    /// carries little, the air asking ten times what soaks in, it gives the air all that soaks
    /// in, its table sinking below the surface until it gives no more, none reaching the sea
    /// through it or seeping out; over a kilometre of gravel, its table lying deep under the
    /// plateau, it gives the air less, the rest reaching the sea; each ground's year balanced
    /// within a thousandth.
    #[test]
    fn a_table_near_the_surface_gives_the_air_what_it_would_take() {
        let grid = CubeSphere::new(32);
        let sides = Sides::new(grid);
        let (radius, density) = (1e6, 1000.0);
        let (heights, sea) = continent(grid, 500.0, 0.0);
        let n = grid.cell_count();
        let runoff = vec![3e-7f32; n];
        let gain = vec![-1e-7f32; n];
        let demand = vec![3e-6f32; n];
        for (name, k, b) in [("tight", 1e-9, 0.0), ("gravel", 1e-3, 1000.0)] {
            let profile = vec![
                Profile {
                    cover: b,
                    cover_conductivity: k,
                    rock_conductivity: k,
                    knee: 100.0,
                };
                n
            ];
            let ground = Ground {
                sides: &sides,
                profile: &profile,
                level: 0.0,
            };
            let (d, table, _) = drain_through(
                grid,
                &heights,
                &sea,
                (&runoff, &gain, &demand),
                (radius, density),
                &ground,
                None,
            );
            let b = table.brief;
            let lakes_out: f64 = 0.0;
            let gone = b.dried + b.seeped + b.to_sea + b.to_lakes + lakes_out;
            println!(
                "{name}: {:.4e} m³/s soaks in, {:.4e} to the air, {:.4e} seeps out, {:.4e} into the sea; the table at the surface over {:.3}; lakes {:.4}",
                b.soaked, b.dried, b.seeped, b.to_sea, b.wet, d.summary.lakes
            );
            assert!(((gone - b.soaked) / b.soaked).abs() < 1e-3, "{name}");
            if name == "tight" {
                assert!(b.dried > 0.999 * b.soaked && b.wet == 0.0, "{name}");
            } else {
                assert!(
                    b.dried < 0.5 * b.soaked && b.to_sea > 0.5 * b.soaked,
                    "{name}"
                );
            }
        }
    }

    /// F435: the rock's transmissivity is its conductivity summed down from the table, within
    /// a part in a thousand of the sum by steps of a metre; its knee, where the crust's law
    /// comes down to crystalline rock's 10^-14.1 m², some 1.07 km down under the Earth's
    /// gravity, 7.3 times deeper under Titan's, as each cell's ground reads it by its body's
    /// gravity; a ground frozen through, its mean below its liquid's melting point, carrying
    /// as continuous permafrost.
    #[test]
    fn the_crust_carries_what_its_depth_law_sums() {
        let earth = Profile::of(0.0, GROUND_CRYSTALLINE_PERMEABILITY, false, (1.0, 1.0));
        let titan = Profile::of(
            0.0,
            GROUND_CRYSTALLINE_PERMEABILITY,
            false,
            (1.0, 1.352 / 9.81),
        );
        println!("knees: {:.0} m, {:.0} m", earth.knee, titan.knee);
        assert!((earth.knee - 1_072.0).abs() < 5.0);
        assert!((titan.knee / earth.knee - 9.81 / 1.352).abs() < 1e-9);
        let water = VOLATILES.iter().find(|v| v.name == "H2O").unwrap();
        let read = |gravity: f64, t: f64| profiles(&[0.0], &|_| t, "felsic", (water, gravity))[0];
        assert!((read(1.352, 290.0).knee / read(9.81, 290.0).knee - 9.81 / 1.352).abs() < 1e-9);
        let (thawed, frozen) = (read(9.81, 274.0), read(9.81, 272.0));
        assert!(thawed.rock_conductivity > 1e-9 && frozen.rock_conductivity < 1e-12);
        for p in [earth, titan] {
            for d in [0.0, 500.0, 3_000.0, 20_000.0] {
                let k = |z: f64| {
                    let k = p.rock_conductivity;
                    if z < p.knee {
                        k
                    } else {
                        k * (z / p.knee).powf(-GROUND_CRUST_FALL)
                    }
                };
                let summed: f64 = (0..2_000_000).map(|i| k(d + i as f64 + 0.5)).sum();
                let tail = k(d + 2e6) * (d + 2e6) / (GROUND_CRUST_FALL - 1.0);
                let law = p.transmissivity(d);
                println!("{d} m: {law:.6e}, summed {:.6e}", summed + tail);
                assert!(((law - summed - tail) / law).abs() < 1e-3, "{d}");
            }
        }
    }

    /// The heads of a polar island of colatitude `cap` under a steady recharge `r`, of
    /// transmissivity `t`, its coast held at 0, on a body of `radius`: the axisymmetric
    /// Poisson's (1 / (a² sin θ)) d/dθ (sin θ dh/dθ) = -R / T, h = (2 R a² / T) ln(cos(θ/2) /
    /// cos(θ0/2)).
    fn island_head(theta: f64, cap: f64, r: f64, t: f64, radius: f64) -> f64 {
        2.0 * r * radius * radius / t * ((0.5 * theta).cos() / (0.5 * cap).cos()).ln()
    }

    /// F435: an island's table under a steady recharge rises as Poisson's mound on a sphere,
    /// its coast where its sea is held, within 2 percent at its top, the recharge all reaching
    /// its sea.
    #[test]
    fn an_island_s_table_rises_as_poisson_s_mound() {
        let grid = CubeSphere::new(48);
        let sides = Sides::new(grid);
        let radius = 1e6;
        let cap = 20f64.to_radians();
        let top = [0.0, 0.0, 1.0];
        let n = grid.cell_count();
        let theta: Vec<f64> = (0..n)
            .map(|k| dot(grid.center(grid.cell(k)), top).clamp(-1.0, 1.0).acos())
            .collect();
        let held: Vec<f64> = theta
            .iter()
            .map(|&th| if th > cap { 0.0 } else { f64::NAN })
            .collect();
        let surface = vec![1e4; n];
        // A cover so thick its transmissivity is near 1 m²/s at any of these depths.
        let profile = vec![
            Profile {
                cover: 1e7,
                cover_conductivity: 1e-7,
                rock_conductivity: 0.0,
                knee: 100.0,
            };
            n
        ];
        let r = 1e-9;
        let recharge = vec![r; n];
        let table = table(
            (grid, &sides, radius),
            &surface,
            &held,
            (&recharge, &vec![0.0; n]),
            &profile,
            None,
        );
        let k_top = (0..n)
            .min_by(|&a, &b| theta[a].total_cmp(&theta[b]))
            .unwrap();
        let t = profile[k_top].transmissivity(surface[k_top] - table.head[k_top]);
        // The coast where the sea is held: the mean colatitude of the held cells beside the
        // land.
        let ring: Vec<f64> = (0..n)
            .filter(|&k| held[k].is_finite() && sides.of(k).any(|(m, _)| held[m].is_nan()))
            .map(|k| theta[k])
            .collect();
        let coast = ring.iter().sum::<f64>() / ring.len() as f64;
        let expected = island_head(theta[k_top], coast, r, t, radius);
        let area: f64 = (0..n)
            .filter(|&k| held[k].is_nan())
            .map(|k| grid.cell_area(grid.cell(k)) * radius * radius)
            .sum();
        let to_sea: f64 = table.exchange.iter().sum();
        println!(
            "top {:.3} m, Poisson's {:.3}; {} rounds, {} steps; recharge {:.4e} m³/s, to the sea {:.4e}",
            table.head[k_top],
            expected,
            table.rounds,
            table.steps,
            r * area,
            to_sea
        );
        assert!(((table.head[k_top] - expected) / expected).abs() < 0.02);
        assert!(((to_sea - r * area) / (r * area)).abs() < 1e-6);
        assert!(table.wet.iter().all(|&w| !w));
    }
}
