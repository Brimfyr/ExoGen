//! The climatological surface wind (F267).
//!
//! The air's flow over the ground through the cycle's year, at the cycle's places: the
//! surface pressure its temperatures and its spin lay out, the ground's drag turning the flow
//! across the isobars, and its strength from the work the atmosphere's heat engine does as it
//! carries its heat from where it is warm to where it is cold.
//!
//! The pressure has two parts. Along each latitude the air's departure from its latitude's
//! mean, a warm place's air lower and a cold place's higher, as a boundary layer's pressure
//! follows its temperature (Lindzen and Nigam 1987): the thermal lows over a summer's land,
//! the highs over a winter's, a locked body's low under its star. And the latitudes' own
//! belts: from the thermal equator, the warmest of the latitudes' means, the pressure rises
//! by the fall of their temperature through the cell that overturns from it, as far as the
//! cell reaches (Held and Hou 1980), its edges at that reach from the equator while its rising
//! branch moves between them (Lindzen and Hou 1988); past them the eddies' belts, westerlies
//! first, alternating at the Rhines scale (Rhines 1975), their strength the one that brings
//! the surface's torque on the air to nothing on each side of the thermal equator, as the
//! atmosphere's angular momentum requires. Where the cell reaches the pole, on a slowly
//! spinning body, its surface branch flows toward the thermal equator all the way, and its
//! own westerlies stand poleward of its trades from a share 1/sqrt(3) of the way to its edge
//! (Held and Hou 1980), where its descending branch brings its high angular momentum down to
//! the ground: a force on the air, not a ridge in its pressure. The two share the torque by
//! the room past the cell's edge for the first of the eddies' belts, the belts coming in by
//! degrees as the room holds them.

use crate::constants::*;
use crate::fields::Field;
use crate::grid::{CubeSphere, lon_lat_from_direction};
use crate::insolation::{read_grid, read_row};
use crate::params::ParamGraph;
use crate::params::derive::physics::{held_hou_reach, mean_molecular_mass};
use rayon::prelude::*;
use std::f64::consts::{FRAC_PI_2, PI};

/// What the wind reads of the body and its air (F267).
#[derive(Clone, Copy, Debug, Default)]
pub struct Air {
    /// The body's spin, rad/s, and its radius, m.
    pub omega: f64,
    pub radius: f64,
    /// The surface pressure, Pa, the air's molecular mass, kg/mol, and its scale height, m.
    pub pressure: f64,
    pub molar_mass: f64,
    pub scale_height: f64,
    /// How far the cell that overturns from the thermal equator reaches from it, rad: the
    /// Held and Hou edge (F32), uncapped.
    pub reach: f64,
}

impl Air {
    pub fn from_graph(graph: &ParamGraph) -> Self {
        let empty = std::collections::BTreeMap::new();
        let composition = graph
            .get("atmosphere.composition")
            .as_map()
            .unwrap_or(&empty);
        let (omega, radius) = (
            graph.f64("rotation.angular_velocity").abs(),
            graph.f64("body.radius"),
        );
        let scale_height = graph.f64("atmosphere.scale_height");
        Self {
            omega,
            radius,
            pressure: graph.f64("atmosphere.surface_pressure").max(0.0),
            molar_mass: mean_molecular_mass(composition),
            scale_height,
            reach: held_hou_reach(
                graph.f64("body.gravity"),
                TROPOPAUSE_SCALE_HEIGHTS * scale_height,
                HELD_HOU_DELTA_H,
                omega,
                radius,
            ),
        }
    }

    /// The depth the ground's drag acts over, m.
    fn surface_layer(&self) -> f64 {
        (SURFACE_LAYER_SHARE * self.scale_height).max(1.0)
    }

    /// The width of the eddies' belt at a latitude, rad of latitude, under a wind of `speed`
    /// m/s: half the Rhines scale's wavelength, pi (2 U / beta)^1/2 over the radius.
    fn belt(&self, lat: f64, speed: f64) -> f64 {
        let beta = 2.0 * self.omega * lat.cos().abs() / self.radius.max(1.0);
        if beta <= 0.0 {
            return f64::INFINITY;
        }
        PI * (2.0 * speed.max(WIND_LEAST_SPEED) / beta).sqrt() / self.radius.max(1.0)
    }

    /// The Rhines phase from latitude `from` to `to`, rad: the integral of the belts'
    /// wavenumber along the meridian, pi per belt.
    fn phase(&self, from: f64, to: f64, speed: f64) -> f64 {
        let scale = (self.radius * self.omega / speed.max(WIND_LEAST_SPEED)).sqrt();
        // Simpson's rule on the square root of the cosine.
        let n = 32;
        let h = (to - from) / n as f64;
        let s: f64 = (0..=n)
            .map(|m| {
                let w = if m == 0 || m == n {
                    1.0
                } else if m % 2 == 1 {
                    4.0
                } else {
                    2.0
                };
                w * (from + h * m as f64).cos().max(0.0).sqrt()
            })
            .sum();
        scale * s * h / 3.0
    }
}

/// A place the wind blows over (F267): its latitude and longitude, rad, its share of the
/// sphere's area, and the shares of it the sea and the caps' ice cover.
#[derive(Clone, Copy, Debug)]
pub struct Site {
    pub lat: f64,
    pub lon: f64,
    pub weight: f64,
    pub sea: f64,
    pub ice: f64,
}

impl Site {
    /// The ground's drag coefficient at the wind's reference height, its sea's, its ice's and
    /// its land's by their shares: (kappa / ln(z / z0))², the law of the wall over each
    /// one's roughness.
    fn drag(&self) -> f64 {
        let law = |z0: f64| (VON_KARMAN / (WIND_REFERENCE_HEIGHT / z0).ln()).powi(2);
        let (sea, ice) = (self.sea.clamp(0.0, 1.0), self.ice.clamp(0.0, 1.0));
        let land = (1.0 - sea - ice).max(0.0);
        sea * law(SEA_ROUGHNESS) + ice * law(ICE_ROUGHNESS) + land * law(LAND_ROUGHNESS)
    }
}

/// The wind through the cycle's year (F267).
#[derive(Clone, Debug)]
pub struct Wind {
    /// The layout's latitudes and longitudes, and the cycle's steps.
    pub layout: (usize, usize),
    pub steps: usize,
    /// The wind at each place at the end of each step, eastward and northward, m/s, step
    /// after step, place after place within each: the tilts' mean.
    pub east: Vec<f32>,
    pub north: Vec<f32>,
    /// At each step: the thermal equator, rad; the trades' ends, where the eddies' or the
    /// cell's own westerlies begin, north and south, rad; the wind's mean speed over the body,
    /// m/s; and the work the atmosphere's heat engine does, W/m² over the body.
    pub equator: Vec<f64>,
    pub ridges: Vec<(f64, f64)>,
    pub speed: Vec<f64>,
    pub work: Vec<f64>,
}

/// The year's wind in brief (F267): its mean speed over the body, m/s, the fastest of its
/// places' year-mean speeds, m/s, the heat engine's mean work, W/m² over the body, the
/// thermal equator's southernmost and northernmost through the year, deg, and the trades'
/// ends north and south in the year's mean, deg.
#[derive(Clone, Copy, Debug)]
pub struct Summary {
    pub speed: f64,
    pub fastest: f64,
    pub work: f64,
    pub equator: (f64, f64),
    pub ridges: (f64, f64),
}

/// One step's wind at one tilt.
struct Blown {
    east: Vec<f64>,
    north: Vec<f64>,
    equator: f64,
    ridges: (f64, f64),
    speed: f64,
    work: f64,
}

impl Wind {
    /// The year's mean at each place: its eastward and northward wind, m/s, and its mean
    /// speed, m/s.
    pub fn year_mean(&self) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        let n = self.layout.0 * self.layout.1;
        let (mut east, mut north, mut speed) = (vec![0.0; n], vec![0.0; n], vec![0.0; n]);
        for k in 0..self.steps {
            for p in 0..n {
                let (u, v) = (self.east[k * n + p] as f64, self.north[k * n + p] as f64);
                east[p] += u;
                north[p] += v;
                speed[p] += u.hypot(v);
            }
        }
        let steps = self.steps.max(1) as f64;
        for x in east.iter_mut().chain(&mut north).chain(&mut speed) {
            *x /= steps;
        }
        (east, north, speed)
    }

    /// The year's mean speed over the body, m/s.
    pub fn mean_speed(&self) -> f64 {
        self.speed.iter().sum::<f64>() / self.speed.len().max(1) as f64
    }

    /// The year in brief.
    pub fn summary(&self) -> Summary {
        let steps = self.steps.max(1) as f64;
        let (.., speed) = self.year_mean();
        let equator = self
            .equator
            .iter()
            .fold((f64::MAX, f64::MIN), |(lo, hi), &e| (lo.min(e), hi.max(e)));
        let ridges = self
            .ridges
            .iter()
            .fold((0.0, 0.0), |(n, s), r| (n + r.0 / steps, s + r.1 / steps));
        Summary {
            speed: self.mean_speed(),
            fastest: speed.iter().fold(0.0, |a: f64, &v| a.max(v)),
            work: self.work.iter().sum::<f64>() / steps,
            equator: (equator.0.to_degrees(), equator.1.to_degrees()),
            ridges: (ridges.0.to_degrees(), ridges.1.to_degrees()),
        }
    }

    /// The year's mean on `grid`, read between the places: its eastward and northward wind
    /// and its mean speed, m/s.
    pub fn fields(&self, grid: CubeSphere) -> [Field<f32>; 3] {
        let (east, north, speed) = self.year_mean();
        let n_lon = self.layout.1;
        let read = |values: &[f64]| {
            Field::from_fn(grid, |_, d| {
                let (lon, lat) = lon_lat_from_direction(d);
                let v = match n_lon {
                    1 => read_row(values, lat),
                    _ => read_grid(values, n_lon, lat, lon),
                };
                v as f32
            })
        };
        [read(&east), read(&north), read(&speed)]
    }
}

/// The wind at each step of the cycle (F267), from the air over each place at the end of each
/// step at each tilt, `air[tilt][step * places + place]`, K, the places, the layout's
/// latitudes and longitudes, and each place's links to its neighbours, the conductance of the
/// air's transport to each, W/m²/K per unit of the place's area. None without air or places.
pub fn blow(
    air: &Air,
    sites: &[Site],
    (n_lat, n_lon): (usize, usize),
    links: &[Vec<(usize, f64)>],
    temperatures: &[&[f64]],
) -> Option<Wind> {
    let n = n_lat * n_lon;
    if air.pressure <= 0.0 || n_lat < 3 || sites.len() != n || temperatures.is_empty() {
        return None;
    }
    let steps = temperatures[0].len() / n;
    if steps == 0 || temperatures.iter().any(|t| t.len() != steps * n) {
        return None;
    }
    let tilts = temperatures.len() as f64;
    let each: Vec<Vec<Blown>> = temperatures
        .iter()
        .map(|series| {
            (0..steps)
                .into_par_iter()
                .map(|k| {
                    step(
                        air,
                        sites,
                        (n_lat, n_lon),
                        links,
                        &series[k * n..(k + 1) * n],
                    )
                })
                .collect()
        })
        .collect();
    let mut wind = Wind {
        layout: (n_lat, n_lon),
        steps,
        east: vec![0.0; steps * n],
        north: vec![0.0; steps * n],
        equator: vec![0.0; steps],
        ridges: vec![(0.0, 0.0); steps],
        speed: vec![0.0; steps],
        work: vec![0.0; steps],
    };
    for tilt in &each {
        for (k, b) in tilt.iter().enumerate() {
            for p in 0..n {
                wind.east[k * n + p] += (b.east[p] / tilts) as f32;
                wind.north[k * n + p] += (b.north[p] / tilts) as f32;
            }
            wind.equator[k] += b.equator / tilts;
            wind.ridges[k].0 += b.ridges.0 / tilts;
            wind.ridges[k].1 += b.ridges.1 / tilts;
            wind.speed[k] += b.speed / tilts;
            wind.work[k] += b.work / tilts;
        }
    }
    Some(wind)
}

/// One step's wind from the air over each place, `t`, K.
fn step(
    air: &Air,
    sites: &[Site],
    (n_lat, n_lon): (usize, usize),
    links: &[Vec<(usize, f64)>],
    t: &[f64],
) -> Blown {
    let n = n_lat * n_lon;
    let a = air.radius.max(1.0);
    let dphi = PI / (n_lat - 1) as f64;
    let dlon = std::f64::consts::TAU / n_lon as f64;
    let lat_of = |i: usize| -FRAC_PI_2 + dphi * i as f64;
    let whole: f64 = sites.iter().map(|s| s.weight).sum::<f64>().max(1e-30);
    // The latitudes' means and the thermal equator, the warmest of them, between the rows.
    let zonal: Vec<f64> = (0..n_lat)
        .map(|i| t[i * n_lon..(i + 1) * n_lon].iter().sum::<f64>() / n_lon as f64)
        .collect();
    let zonal_at = |lat: f64| {
        let x = ((lat + FRAC_PI_2) / dphi).clamp(0.0, (n_lat - 1) as f64);
        let i = (x.floor() as usize).min(n_lat - 2);
        let s = x - i as f64;
        zonal[i] * (1.0 - s) + zonal[i + 1] * s
    };
    let top = (0..n_lat)
        .max_by(|&i, &j| zonal[i].total_cmp(&zonal[j]))
        .unwrap_or(n_lat / 2);
    let mut equator = lat_of(top);
    if top > 0 && top + 1 < n_lat {
        let (l, c, r) = (zonal[top - 1], zonal[top], zonal[top + 1]);
        let curve = l - 2.0 * c + r;
        if curve < 0.0 {
            equator += dphi * (0.5 * (l - r) / curve).clamp(-0.5, 0.5);
        }
    }
    // The cells keep their edges at the Held and Hou reach from the equator while their
    // rising branch moves between them with the thermal equator, as a cell heated off the
    // equator does (Lindzen and Hou 1988), no farther than a row inside them.
    let reach = air.reach.min(FRAC_PI_2);
    let bound = (reach - dphi).clamp(0.0, FRAC_PI_2 - dphi);
    let equator = equator.clamp(-bound, bound);
    // The heat engine's work (Lorenz 2002): the heat the air's transport takes from where it
    // carries heat away and gives where it brings it, at their mean temperatures, its Carnot
    // share; half of it spent against the ground (Peixoto and Oort 1992).
    let (mut taken, mut taken_t, mut given, mut given_t) = (0.0, 0.0, 0.0, 0.0);
    for p in 0..n {
        let q: f64 = links[p].iter().map(|&(o, c)| c * (t[o] - t[p])).sum();
        let w = sites[p].weight;
        if q < 0.0 {
            taken += -q * w;
            taken_t += -q * w * t[p];
        } else {
            given += q * w;
            given_t += q * w * t[p];
        }
    }
    let work = if taken > 0.0 && given > 0.0 {
        let (hot, cold) = (taken_t / taken, given_t / given);
        (taken / whole * (1.0 - cold / hot)).max(0.0)
    } else {
        0.0
    };
    let dissipation = SURFACE_DISSIPATION_SHARE * work;
    let density: Vec<f64> = t
        .iter()
        .map(|&t| air.pressure * air.molar_mass / (GAS_CONSTANT * t.max(1.0)))
        .collect();
    let drag: Vec<f64> = sites.iter().map(Site::drag).collect();
    let rows = |f: &dyn Fn(usize) -> f64| -> Vec<f64> {
        (0..n_lat)
            .map(|i| (0..n_lon).map(|j| f(i * n_lon + j)).sum::<f64>() / n_lon as f64)
            .collect()
    };
    let row_density = rows(&|p| density[p]);
    let row_drag = rows(&|p| drag[p]);
    let row_weight = rows(&|p| sites[p].weight);
    // The air's departure from its latitude's mean, the pressure's along the latitude.
    let departure: Vec<f64> = (0..n).map(|p| t[p] - zonal[p / n_lon]).collect();
    let depth = air.surface_layer();
    let mut speed = WIND_START_SPEED;
    let mut out = Blown {
        east: vec![0.0; n],
        north: vec![0.0; n],
        equator,
        ridges: (equator, equator),
        speed: 0.0,
        work,
    };
    for _ in 0..WIND_ITERATIONS {
        let friction = |drag: f64| drag * speed.max(WIND_LEAST_SPEED) / depth;
        // Each side of the thermal equator: the cell's edge, the trades' end, and the
        // latitudes' pressure gradient, K per metre northward, through the trades by the fall
        // of their temperature, past them the belts that bring the surface's torque to
        // nothing.
        let mut gradient = vec![0.0; n_lat];
        let mut force = vec![0.0; n_lat];
        let mut ridges = (equator, equator);
        for side in [1.0f64, -1.0] {
            // The cell's edge, how far the room past it falls short of the first of the eddies'
            // belts, and where its trades end and its own westerlies begin.
            let edge = side * reach;
            let room = FRAC_PI_2 - side * edge;
            let width = air.belt(edge, speed);
            let held = (1.0 - room / width).clamp(0.0, 1.0);
            let ridge = equator + (edge - equator) * (1.0 - (1.0 - 1.0 / 3f64.sqrt()) * held);
            if side > 0.0 {
                ridges.0 = ridge;
            } else {
                ridges.1 = ridge;
            }
            let span = (edge - equator).abs().max(1e-6);
            let outer = (edge - ridge).abs();
            let rise = (zonal_at(equator) - zonal_at(edge)).max(0.0);
            let trades = 0.5 * PI * rise / (a * span);
            let on = |lat: f64| side * (lat - equator) >= 0.0 && (lat != equator || side > 0.0);
            // The cell's gradient through it, the eddies' belts' past it, and the force of
            // the cell's own westerlies over its outer part, the last two at unit strength.
            let mut trade = vec![0.0; n_lat];
            let mut belt = vec![0.0; n_lat];
            let mut push = vec![0.0; n_lat];
            for i in 1..n_lat - 1 {
                let lat = lat_of(i);
                if !on(lat) {
                    continue;
                }
                let from = side * (lat - equator);
                if from <= span {
                    trade[i] = side * trades * (PI * from / span).sin();
                    let past = side * (lat - ridge);
                    if outer > 1e-6 && past > 0.0 {
                        push[i] = (PI * past / outer).sin() / row_density[i];
                    }
                } else {
                    let (lo, hi) = if side > 0.0 { (edge, lat) } else { (lat, edge) };
                    belt[i] = -side * air.phase(lo, hi, speed).sin() * lat.cos();
                }
            }
            // The surface's torque on the air through each row, rho eps h u cos(lat) over its
            // area, the zonal wind from a gradient and a force by the ground's friction and
            // the spin.
            let torque = |g: &[f64], x: &[f64]| -> f64 {
                (1..n_lat - 1)
                    .map(|i| {
                        let lat = lat_of(i);
                        let f = 2.0 * air.omega * lat.sin();
                        let eps = friction(row_drag[i]);
                        let u = (eps * x[i] - f * g[i] / row_density[i])
                            / (eps * eps + f * f).max(1e-300);
                        row_density[i] * eps * u * lat.cos() * row_weight[i]
                    })
                    .sum()
            };
            let none = vec![0.0; n_lat];
            let pull = torque(&trade, &none);
            let (by_belts, by_cell) = (torque(&belt, &none), torque(&none, &push));
            let (mut belts, mut cell) = (0.0, 0.0);
            if pull < 0.0 {
                let owed = -pull;
                if by_belts > 0.0 && by_cell > 0.0 {
                    (belts, cell) = ((1.0 - held) * owed / by_belts, held * owed / by_cell);
                } else if by_belts > 0.0 {
                    belts = owed / by_belts;
                } else if by_cell > 0.0 {
                    cell = owed / by_cell;
                }
            }
            for i in 0..n_lat {
                gradient[i] += trade[i] + belts * belt[i];
                force[i] += cell * push[i];
            }
        }
        // The wind at each place from the gradient, the latitudes' and its departure's.
        let mut east = vec![0.0; n];
        let mut north = vec![0.0; n];
        for (i, &along) in gradient.iter().enumerate().take(n_lat - 1).skip(1) {
            let lat = lat_of(i);
            let f = 2.0 * air.omega * lat.sin();
            let across = 2.0 * a * lat.cos().max(1e-6) * dlon;
            for j in 0..n_lon {
                let p = i * n_lon + j;
                let gx = if n_lon > 1 {
                    let (w, e) = (
                        i * n_lon + (j + n_lon - 1) % n_lon,
                        i * n_lon + (j + 1) % n_lon,
                    );
                    -(departure[e] - departure[w]) / across
                } else {
                    0.0
                };
                let gy = along - (departure[p + n_lon] - departure[p - n_lon]) / (2.0 * a * dphi);
                let eps = friction(drag[p]);
                let den = (eps * eps + f * f).max(1e-300);
                let x = force[i];
                east[p] = (eps * x - (eps * gx + f * gy) / density[p]) / den;
                north[p] = ((f * gx - eps * gy) / density[p] - f * x) / den;
            }
        }
        // Its strength: the surface spends its share of the heat engine's work, the mean
        // wind the mean of a Rayleigh distribution of speeds, whose cube's mean is 6/pi its
        // mean's cube.
        let spent: f64 = (0..n)
            .map(|p| {
                let v = east[p].hypot(north[p]);
                sites[p].weight * density[p] * drag[p] * RAYLEIGH_CUBE * v * v * v
            })
            .sum::<f64>()
            / whole;
        let scale = if spent > 0.0 && dissipation > 0.0 {
            (dissipation / spent).cbrt()
        } else {
            0.0
        };
        for x in east.iter_mut().chain(&mut north) {
            *x *= scale;
        }
        let mean = (0..n)
            .map(|p| sites[p].weight * east[p].hypot(north[p]))
            .sum::<f64>()
            / whole;
        let settled = (mean - speed).abs() <= WIND_SETTLED * speed.max(WIND_LEAST_SPEED);
        out.east = east;
        out.north = north;
        out.ridges = ridges;
        out.speed = mean;
        speed = mean;
        if settled {
            break;
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A layout's places on a fast spinner's grid, all sea or all land.
    fn sites((n_lat, n_lon): (usize, usize), sea: f64) -> Vec<Site> {
        let dphi = PI / (n_lat - 1) as f64;
        (0..n_lat * n_lon)
            .map(|k| {
                let (i, j) = (k / n_lon, k % n_lon);
                let lat = -FRAC_PI_2 + dphi * i as f64;
                let (lo, hi) = (
                    (lat - dphi / 2.0).max(-FRAC_PI_2),
                    (lat + dphi / 2.0).min(FRAC_PI_2),
                );
                Site {
                    lat,
                    lon: -PI + std::f64::consts::TAU * j as f64 / n_lon as f64,
                    weight: (hi.sin() - lo.sin()) / n_lon as f64,
                    sea,
                    ice: 0.0,
                }
            })
            .collect()
    }

    /// The layout's links as the cycle lays them, at a transport of `d` W/m²/K.
    fn links((n_lat, n_lon): (usize, usize), d: f64) -> Vec<Vec<(usize, f64)>> {
        let (dphi, dlon) = (
            PI / (n_lat - 1) as f64,
            std::f64::consts::TAU / n_lon as f64,
        );
        (0..n_lat * n_lon)
            .map(|k| {
                let (i, j) = (k / n_lon, k % n_lon);
                let lat = -FRAC_PI_2 + dphi * i as f64;
                let (lo, hi) = (
                    (lat - dphi / 2.0).max(-FRAC_PI_2),
                    (lat + dphi / 2.0).min(FRAC_PI_2),
                );
                let band = hi.sin() - lo.sin();
                let mut out = Vec::new();
                if i > 0 {
                    out.push((k - n_lon, d * lo.cos() / (dphi * band)));
                }
                if i + 1 < n_lat {
                    out.push((k + n_lon, d * hi.cos() / (dphi * band)));
                }
                let across = (hi - lo) / (dlon * (0.5 * (lo + hi)).cos() * dlon * band);
                out.push((i * n_lon + (j + n_lon - 1) % n_lon, d * across));
                out.push((i * n_lon + (j + 1) % n_lon, d * across));
                out
            })
            .collect()
    }

    fn earth() -> Air {
        Air {
            omega: 7.292e-5,
            radius: 6.371e6,
            pressure: 101_325.0,
            molar_mass: 0.028_97,
            scale_height: 8_400.0,
            reach: 32f64.to_radians(),
        }
    }

    /// The zonal mean of the eastward and northward wind at a latitude, deg.
    fn zonal(wind: &Wind, lat: f64) -> (f64, f64) {
        let (n_lat, n_lon) = wind.layout;
        let i = ((lat + 90.0) / 180.0 * (n_lat - 1) as f64).round() as usize;
        let (east, north, _) = wind.year_mean();
        let row = i * n_lon..(i + 1) * n_lon;
        (
            east[row.clone()].iter().sum::<f64>() / n_lon as f64,
            north[row].iter().sum::<f64>() / n_lon as f64,
        )
    }

    /// F267: an Earth-like aquaplanet's air, warm at its equator and cold at its poles, blows
    /// in Earth's three belts: easterly trades toward the equator within some 30°, the
    /// westerlies poleward of them to some 60°, easterlies at the poles; some 5 to 9 m/s,
    /// near the 7 m/s of Earth's ocean winds (Wentz et al. 2007), under a heat engine
    /// working at about the 2 to 3 W/m² of Earth's (Peixoto and Oort 1992).
    #[test]
    fn an_earth_like_s_winds_blow_in_three_belts() {
        let layout = (91, 90);
        let sites = sites(layout, 1.0);
        let air: Vec<f64> = sites
            .iter()
            .map(|s| 300.0 - 50.0 * s.lat.sin().powi(2))
            .collect();
        let wind = blow(&earth(), &sites, layout, &links(layout, 0.6), &[&air]).unwrap();
        let (trades, toward) = zonal(&wind, 15.0);
        let (westerlies, away) = zonal(&wind, 45.0);
        let (polar, _) = zonal(&wind, 75.0);
        let (south, _) = zonal(&wind, -45.0);
        assert!(trades < -1.0 && toward < 0.0, "15°N: {trades}, {toward}");
        assert!(westerlies > 1.0 && away > 0.0, "45°N: {westerlies}, {away}");
        assert!(south > 1.0, "45°S: {south}");
        assert!(polar < 0.0, "75°N: {polar}");
        let (mean, work) = (wind.mean_speed(), wind.work[0]);
        assert!((5.0..9.0).contains(&mean), "{mean} m/s");
        assert!((1.5..4.0).contains(&work), "{work} W/m²");
    }

    /// F267: a faster spinner's belts are the narrower, so that it has more of them: at four
    /// times the Earth's spin and an eighth of its cell's reach, the westerlies' belt from
    /// near 8° gives way to easterlies before 40° and to westerlies again before the pole.
    #[test]
    fn a_fast_spinner_has_more_belts() {
        let layout = (91, 90);
        let sites = sites(layout, 1.0);
        let air: Vec<f64> = sites
            .iter()
            .map(|s| 300.0 - 50.0 * s.lat.sin().powi(2))
            .collect();
        let fast = Air {
            omega: 4.0 * 7.292e-5,
            reach: 8f64.to_radians(),
            ..earth()
        };
        let wind = blow(&fast, &sites, layout, &links(layout, 0.6), &[&air]).unwrap();
        let signs: Vec<f64> = (4..44)
            .step_by(2)
            .map(|lat| zonal(&wind, lat as f64 * 2.0).0)
            .collect();
        let turns = signs.windows(2).filter(|w| w[0] * w[1] < 0.0).count();
        assert!(turns >= 3, "{signs:?}");
    }

    /// F267: a body locked to its star draws its air from its night side to the point under
    /// its star: at its equator the wind blows toward the substellar longitude from either
    /// side.
    #[test]
    fn a_locked_body_s_air_blows_toward_its_star() {
        let layout = (46, 90);
        let sites = sites(layout, 1.0);
        let air: Vec<f64> = sites
            .iter()
            .map(|s| 230.0 + 50.0 * (s.lat.cos() * s.lon.cos()).max(0.0))
            .collect();
        let locked = Air {
            omega: TAU_6_DAYS,
            radius: 5.85e6,
            reach: PI,
            ..earth()
        };
        let wind = blow(&locked, &sites, layout, &links(layout, 0.6), &[&air]).unwrap();
        let (n_lat, n_lon) = wind.layout;
        let (east, ..) = wind.year_mean();
        let row = (n_lat / 2) * n_lon;
        let at = |lon: f64| {
            let j = ((lon + 180.0) / 360.0 * n_lon as f64).round() as usize % n_lon;
            east[row + j]
        };
        assert!(
            at(-40.0) > 0.0 && at(40.0) < 0.0,
            "{} {}",
            at(-40.0),
            at(40.0)
        );
        for lat in [-66.0, -30.0, 30.0, 66.0] {
            let (_, north) = zonal(&wind, lat);
            assert!(north * lat < 0.0, "{lat}°: {north}");
        }
    }

    const TAU_6_DAYS: f64 = std::f64::consts::TAU / (6.1 * 86_400.0);

    /// F267: a cell heated far from its equator keeps its edges: where a Mars-like body's air
    /// is warmest near its south pole in its southern summer, its thermal equator stands no
    /// farther south than its cell's reach, 44°, and its winter cell's surface branch crosses
    /// the equator toward the summer, as Mars's solstitial cell spans from the summer's
    /// subtropics to the winter's mid-latitudes in its models (Haberle et al. 1993).
    #[test]
    fn a_cell_heated_off_its_equator_keeps_its_edges() {
        let layout = (91, 90);
        let sites = sites(layout, 0.0);
        let air: Vec<f64> = sites
            .iter()
            .map(|s| 170.0 + 60.0 * (-((s.lat + 70f64.to_radians()) / 0.6).powi(2)).exp())
            .collect();
        let mars = Air {
            omega: 7.088e-5,
            radius: 3.39e6,
            pressure: 610.0,
            molar_mass: 0.0434,
            scale_height: 11_100.0,
            reach: 44f64.to_radians(),
        };
        let wind = blow(&mars, &sites, layout, &links(layout, 0.002), &[&air]).unwrap();
        let equator = wind.equator[0].to_degrees();
        assert!((-44.0..-40.0).contains(&equator), "{equator}°");
        let (_, across) = zonal(&wind, 0.0);
        let (_, north) = zonal(&wind, 20.0);
        assert!(across < 0.0 && north < 0.0, "{across} {north}");
    }

    /// F267: a slowly spinning body's cell reaches its poles: its trades end some 50° from
    /// the equator, where its own westerlies take over, with no easterlies poleward of them,
    /// and its surface branch flows toward the equator all the way.
    #[test]
    fn a_slow_spinner_s_cell_reaches_its_poles() {
        let layout = (91, 90);
        let sites = sites(layout, 0.0);
        let air: Vec<f64> = sites
            .iter()
            .map(|s| 94.0 - 3.0 * s.lat.sin().powi(2))
            .collect();
        let titan = Air {
            omega: 4.56e-6,
            radius: 2.575e6,
            pressure: 146_700.0,
            molar_mass: 0.028,
            scale_height: 20_000.0,
            reach: PI,
        };
        let wind = blow(&titan, &sites, layout, &links(layout, 0.5), &[&air]).unwrap();
        let ridge = wind.ridges[0].0.to_degrees();
        assert!((45.0..60.0).contains(&ridge), "{ridge}°");
        let (low, toward) = zonal(&wind, 20.0);
        let (high, back) = zonal(&wind, 70.0);
        assert!(low < 0.0 && high > 0.0, "{low} {high}");
        assert!(toward < 0.0 && back < 0.0, "{toward} {back}");
    }
}
