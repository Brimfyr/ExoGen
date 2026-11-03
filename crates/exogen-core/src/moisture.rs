//! The moisture (F268).
//!
//! The working liquid's vapour through the cycle's year at its places: a column of vapour over
//! each place, evaporated from its open sea and its wet land at its own wind (F267), carried by
//! the wind and mixed by the eddies that carry the heat, and rained out as the column nears
//! saturation.
//!
//! The column at saturation is the saturated vapour's density at the ground over the vapour's
//! scale height, R_v T² / (L Γ), from Clausius and Clapeyron along the lapse rate: e T / (L Γ),
//! some 30 kg/m² over Earth's 288 K, less over high ground, whose air is the colder. The rain
//! rises with the column's saturation r, the vapour it holds over what it would hold
//! saturated, as Bretherton et al. (2004) measured it over the tropical oceans, exp(15.6 (r -
//! 0.603)) mm a day over their saturated columns of some 65 kg/m², and in proportion to the
//! saturated column elsewhere: no threshold, only a steep rise. The eddies mix the vapour as
//! they mix the heat, at the air's transport over its column's heat capacity. The land holds
//! what rains on it in a bucket of 15 cm of the liquid (Manabe 1969), evaporating at its full
//! rate while three quarters full and in proportion below, and sheds what overflows as
//! runoff. The sea under its ice and the caps' ice give none.

use crate::constants::*;
use crate::fields::Field;
use crate::grid::{CubeSphere, lon_lat_from_direction};
use crate::insolation::{read_grid, read_row};
use crate::params::ParamGraph;
use crate::params::derive::physics::{
    mean_heat_capacity, mean_molecular_mass, vapour_pressure, volatile,
};
use crate::wind::Site;
use std::f64::consts::{PI, TAU};

/// What the moisture reads of the body, its air and its liquid (F268).
#[derive(Clone, Copy, Debug)]
pub struct Air {
    /// The working liquid's vapour.
    pub vapour: &'static Volatile,
    /// The surface pressure at the datum, Pa, the air's molecular mass, kg/mol, its heat
    /// capacity, J/kg/K, and its scale height, m.
    pub pressure: f64,
    pub molar_mass: f64,
    pub heat_capacity: f64,
    pub scale_height: f64,
    /// The body's gravity, m/s², and radius, m, and the air's lapse rate, K/m.
    pub gravity: f64,
    pub radius: f64,
    pub lapse: f64,
    /// The air's transport, W/m²/K, and the share of it a sea's water carries across the
    /// sea's share of a boundary (F245).
    pub transport: f64,
    pub water: f64,
}

impl Air {
    /// A body's air and working liquid, none without either; `transport` the climate's,
    /// W/m²/K, `water` the share of it a sea's water carries.
    pub fn from_graph(graph: &ParamGraph, transport: f64, water: f64) -> Option<Self> {
        let working = graph.str("volatiles.working_liquid");
        let species = crate::stages::caps::liquid_species(working);
        let vapour = volatile(species.first()?)?;
        let pressure = graph.f64("atmosphere.surface_pressure");
        if pressure <= 0.0 || transport <= 0.0 {
            return None;
        }
        let empty = std::collections::BTreeMap::new();
        let composition = graph
            .get("atmosphere.composition")
            .as_map()
            .unwrap_or(&empty);
        Some(Self {
            vapour,
            pressure,
            molar_mass: mean_molecular_mass(composition),
            heat_capacity: mean_heat_capacity(composition),
            scale_height: graph.f64("atmosphere.scale_height"),
            gravity: graph.f64("body.gravity"),
            radius: graph.f64("body.radius"),
            lapse: graph.f64("climate.lapse_rate"),
            transport,
            water,
        })
    }

    /// The latent heat at `t`, J/kg: of vaporisation above the triple point, of sublimation
    /// below it, as the vapour pressure reads it.
    fn latent(&self, t: f64) -> f64 {
        let v = self.vapour;
        if t >= v.t_triple { v.l_vap } else { v.l_sub }
    }

    /// The saturated column over ground whose air stands at `t`, kg/m²: the vapour's
    /// saturated density over its scale height, e T / (L Γ).
    pub fn saturated_column(&self, t: f64) -> f64 {
        let t = t.max(1.0);
        vapour_pressure(self.vapour, t) * t
            / (self.latent(t) * self.lapse.max(MOISTURE_LEAST_LAPSE))
    }

    /// The saturated specific humidity at `t` under `p`, Pa.
    fn saturated_humidity(&self, t: f64, p: f64) -> f64 {
        let e = vapour_pressure(self.vapour, t.max(1.0)).min(MOISTURE_MOST_PARTIAL * p);
        let ratio = GAS_CONSTANT / self.vapour.r_specific / self.molar_mass;
        ratio * e / (p - (1.0 - ratio) * e)
    }

    /// The precipitation from a column of `w` kg/m² whose saturation holds `saturated`,
    /// kg/m²/s (Bretherton et al. 2004).
    pub fn rain(&self, w: f64, saturated: f64) -> f64 {
        if saturated <= 0.0 || w <= 0.0 {
            return 0.0;
        }
        let r = w / saturated;
        let rise = (RAIN_SATURATION_SLOPE * (r - RAIN_SATURATION_ONSET)).min(MOISTURE_MOST_RISE);
        saturated / (RAIN_TIMESCALE_DAYS * SECONDS_PER_DAY) * rise.exp()
    }

    /// The share of the precipitation at `t` that falls as snow: half at a kelvin above the
    /// liquid's freezing point, the rest in a logistic a kelvin wide (Jennings et al. 2018).
    fn snow_share(&self, t: f64) -> f64 {
        let x = (t - self.vapour.t_triple - RAIN_SNOW_OFFSET) / RAIN_SNOW_WIDTH;
        1.0 / (1.0 + x.clamp(-50.0, 50.0).exp())
    }
}

/// The implicit mixing's solver (F268): (diag + links) x = rhs, the links each place's
/// neighbours and its coupling to each, x the start on entry.
pub type Solve<'s> = dyn Fn(&[Vec<(usize, f64)>], &[f64], &[f64], &mut [f64]) + 's;

/// The weather the moisture reads at a step (F268): the air over each place at the datum, K,
/// the wind, m/s, the land's surface and the sea's at the datum, K, the sea's ice, kg/m², and
/// the gusts about the wind.
pub struct Weather<'a> {
    pub air: &'a [f64],
    pub east: &'a [f32],
    pub north: &'a [f32],
    pub land: Vec<f64>,
    pub sea: Vec<f64>,
    pub ice: Vec<f64>,
    /// The body's mean wind at the step, m/s, the weather's gusts about each place's mean
    /// wind.
    pub gust: f64,
}

/// The moisture through the cycle's year (F268).
#[derive(Clone, Debug)]
pub struct Moisture {
    /// The layout's latitudes and longitudes, and the cycle's steps.
    pub layout: (usize, usize),
    pub steps: usize,
    /// The precipitation at each place through each step, kg/m²/s, step after step, place
    /// after place within each, and the share of it that falls as snow.
    pub falling: Vec<f32>,
    pub snow: Vec<f32>,
    /// The year's mean at each place: its precipitation and evaporation over the place, its
    /// potential evaporation, of the open sea over its sea and of wet ground over its land,
    /// and its land's runoff over its land, kg/m²/s; and its column's mean saturation.
    pub precipitation: Vec<f64>,
    pub evaporation: Vec<f64>,
    pub potential: Vec<f64>,
    pub runoff: Vec<f64>,
    pub saturation: Vec<f64>,
    /// The share of each place's land its bucket keeps wet through the year (F269).
    pub wetness: Vec<f64>,
    /// Each place's share of the sphere and of its land.
    pub weights: Vec<f64>,
    pub land: Vec<f64>,
}

/// The year's moisture in brief (F268), kg/m² a year: precipitation and evaporation over the
/// body; over its land, precipitation, evaporation and runoff; and the shares of the land whose
/// precipitation over its potential evaporation stands under 0.2, arid, and under 0.05,
/// hyper-arid (UNEP 1992).
#[derive(Clone, Copy, Debug)]
pub struct Summary {
    pub precipitation: f64,
    pub evaporation: f64,
    pub land_precipitation: f64,
    pub land_evaporation: f64,
    pub runoff: f64,
    pub arid: f64,
    pub hyper_arid: f64,
}

impl Moisture {
    /// The year in brief.
    pub fn summary(&self) -> Summary {
        let year = SECONDS_PER_DAY * 365.25;
        let (mut whole, mut land) = (0.0, 0.0);
        let (mut p, mut e, mut lp, mut le, mut ro) = (0.0, 0.0, 0.0, 0.0, 0.0);
        let (mut arid, mut hyper) = (0.0, 0.0);
        for k in 0..self.precipitation.len() {
            let (w, l) = (self.weights[k], self.weights[k] * self.land[k]);
            whole += w;
            p += w * self.precipitation[k];
            e += w * self.evaporation[k];
            if l > 0.0 {
                land += l;
                lp += l * self.precipitation[k];
                le += l * self.land_evaporation(k);
                ro += l * self.runoff[k];
                let index = self.aridity(k);
                if index < ARID_INDEX {
                    arid += l;
                }
                if index < HYPER_ARID_INDEX {
                    hyper += l;
                }
            }
        }
        let over = |x: f64, a: f64| if a > 0.0 { x / a * year } else { 0.0 };
        Summary {
            precipitation: over(p, whole),
            evaporation: over(e, whole),
            land_precipitation: over(lp, land),
            land_evaporation: over(le, land),
            runoff: over(ro, land),
            arid: if land > 0.0 { arid / land } else { 0.0 },
            hyper_arid: if land > 0.0 { hyper / land } else { 0.0 },
        }
    }

    /// The land's evaporation at place `k` over its land, kg/m²/s: what the place gives that
    /// its land and not its sea gives, in the year's balance its precipitation less its
    /// runoff.
    fn land_evaporation(&self, k: usize) -> f64 {
        (self.precipitation[k] - self.runoff[k]).max(0.0)
    }

    /// The year's precipitation, kg/m² a year, and its aridity index, held under
    /// ARIDITY_MAP_MOST, on `grid`, read between the places.
    pub fn fields(&self, grid: CubeSphere) -> [Field<f32>; 2] {
        let year = SECONDS_PER_DAY * 365.25;
        let n_lon = self.layout.1;
        let rain: Vec<f64> = self.precipitation.iter().map(|&p| p * year).collect();
        let dry: Vec<f64> = (0..rain.len())
            .map(|k| self.aridity(k).min(ARIDITY_MAP_MOST))
            .collect();
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
        [read(&rain), read(&dry)]
    }

    /// The aridity index at place `k`: its precipitation over its potential evaporation.
    pub fn aridity(&self, k: usize) -> f64 {
        if self.potential[k] > 0.0 {
            self.precipitation[k] / self.potential[k]
        } else {
            f64::INFINITY
        }
    }
}

/// The moisture through the cycle (F268): the year run again and again, `steps` steps of a
/// cycle of `year` s, until its land's buckets have come to their year, each step's weather
/// from `weather`; `heights` each place's ground over the sea, m; `solve` the implicit
/// mixing's system, (diag + links) x = rhs, by the cycle's transport solver.
#[allow(clippy::too_many_arguments)]
pub fn rain<'a>(
    air: &Air,
    sites: &[Site],
    (n_lat, n_lon): (usize, usize),
    heights: &[f64],
    (year, steps): (f64, usize),
    weather: &dyn Fn(usize) -> Weather<'a>,
    solve: &Solve<'_>,
) -> Option<Moisture> {
    let n = n_lat * n_lon;
    if n_lat < 3 || sites.len() != n || heights.len() != n || steps == 0 || year <= 0.0 {
        return None;
    }
    let dt = year / steps as f64;
    let a = air.radius.max(1.0);
    let (dphi, dlon) = (PI / (n_lat - 1) as f64, TAU / n_lon as f64);
    let around = if n_lon > 1 { dlon } else { TAU };
    let edges = |i: usize| {
        let lat = -PI / 2.0 + dphi * i as f64;
        (
            (lat - dphi / 2.0).max(-PI / 2.0),
            (lat + dphi / 2.0).min(PI / 2.0),
        )
    };
    let area: Vec<f64> = (0..n)
        .map(|p| {
            let (lo, hi) = edges(p / n_lon);
            a * a * (hi.sin() - lo.sin()) * around
        })
        .collect();
    // The eddies' mixing through each link, per second, the air's transport over its
    // column's heat capacity, the step's share of it for the implicit system.
    let mix: Vec<Vec<(usize, f64)>> = crate::climate::place_links((n_lat, n_lon))
        .iter()
        .enumerate()
        .map(|(p, links)| {
            links
                .iter()
                .map(|&(q, g)| {
                    let wet = sites[p].sea.min(sites[q].sea);
                    let rate = air.transport * g * (1.0 - air.water * wet) * air.gravity
                        / (air.heat_capacity * air.pressure);
                    (q, dt * rate)
                })
                .collect()
        })
        .collect();
    // Each place's ground: its land's share, its height and its air's pressure there.
    let land: Vec<f64> = sites
        .iter()
        .map(|s| (1.0 - s.sea - s.ice).clamp(0.0, 1.0))
        .collect();
    let pressure: Vec<f64> = heights
        .iter()
        .map(|&h| air.pressure * (-h / air.scale_height.max(1.0)).exp())
        .collect();
    let capacity = SOIL_CAPACITY * air.vapour.liquid_density;
    let years = ((MOISTURE_SPIN_UP / year).ceil() as usize)
        .max(MOISTURE_LEAST_YEARS)
        .min((MOISTURE_MOST_STEPS / steps).max(MOISTURE_LEAST_YEARS));
    let first = weather(0);
    let mut w: Vec<f64> = (0..n)
        .map(|p| {
            let t = first.air[p] - air.lapse * heights[p];
            MOISTURE_START_SATURATION * air.saturated_column(t)
        })
        .collect();
    let mut bucket = vec![0.5 * capacity; n];
    let mut out = Moisture {
        layout: (n_lat, n_lon),
        steps,
        falling: vec![0.0; steps * n],
        snow: vec![0.0; steps * n],
        precipitation: vec![0.0; n],
        evaporation: vec![0.0; n],
        potential: vec![0.0; n],
        runoff: vec![0.0; n],
        saturation: vec![0.0; n],
        wetness: vec![0.0; n],
        weights: sites.iter().map(|s| s.weight).collect(),
        land: land.clone(),
    };
    let (mut diag, mut rhs) = (vec![0.0; n], vec![0.0; n]);
    let (mut gives, mut takes) = (vec![0.0; n], vec![0.0; n]);
    let (mut exchange, mut beta) = (vec![0.0; n], vec![0.0; n]);
    let (mut humid, mut sea_q, mut soil_q) = (vec![0.0; n], vec![0.0; n], vec![0.0; n]);
    let (mut east, mut north) = (vec![0.0; n], vec![0.0; n]);
    for round in 0..years {
        let last = round + 1 == years;
        for k in 0..steps {
            let now = weather(k);
            // The ground's air, its saturated column and humidity, and what the sea and the
            // wet land give, gives - takes * w, kg/m²/s, at the wind's speed through the
            // weather about its mean, sqrt(|V|² + U²), U the body's mean wind: the gusts that
            // keep a calm mean wind's ridges evaporating.
            let ground: Vec<f64> = (0..n)
                .map(|p| now.air[p] - air.lapse * heights[p])
                .collect();
            let saturated: Vec<f64> = ground.iter().map(|&t| air.saturated_column(t)).collect();
            for p in 0..n {
                let speed = (now.east[p] as f64)
                    .hypot(now.north[p] as f64)
                    .hypot(now.gust);
                let rho = pressure[p] * air.molar_mass / (GAS_CONSTANT * ground[p].max(1.0));
                exchange[p] = rho * AIR_EXCHANGE_COEFFICIENT * speed;
                beta[p] = (bucket[p] / (BUCKET_FULL_SHARE * capacity)).min(1.0);
                humid[p] = air.saturated_humidity(ground[p], pressure[p]);
                sea_q[p] = air.saturated_humidity(now.sea[p], air.pressure);
                soil_q[p] =
                    air.saturated_humidity(now.land[p] - air.lapse * heights[p], pressure[p]);
                let open = sites[p].sea * (-now.ice[p] / SEA_ICE_COVER_MASS).exp();
                let wet = land[p] * beta[p];
                gives[p] = exchange[p] * (open * sea_q[p] + wet * soil_q[p]);
                takes[p] = exchange[p] * (open + wet) * humid[p] / saturated[p].max(1e-30);
            }
            // The wind carries the vapour across the places' faces, upwind, in steps short
            // enough that no place gives more than it holds.
            for p in 0..n {
                let (i, j) = (p / n_lon, p % n_lon);
                let (lo, hi) = edges(i);
                east[p] = if n_lon > 1 {
                    let e = i * n_lon + (j + 1) % n_lon;
                    0.5 * (now.east[p] + now.east[e]) as f64 * a * (hi - lo)
                } else {
                    0.0
                };
                north[p] = if i + 1 < n_lat {
                    0.5 * (now.north[p] + now.north[p + n_lon]) as f64 * a * hi.cos() * around
                } else {
                    0.0
                };
            }
            let mut leaving = vec![0.0; n];
            for p in 0..n {
                let (i, j) = (p / n_lon, p % n_lon);
                let mut add = |q: usize, f: f64| {
                    if f > 0.0 {
                        leaving[q] += f;
                    }
                };
                if n_lon > 1 {
                    let e = i * n_lon + (j + 1) % n_lon;
                    add(p, east[p]);
                    add(e, -east[p]);
                }
                if i + 1 < n_lat {
                    add(p, north[p]);
                    add(p + n_lon, -north[p]);
                }
            }
            let fastest = (0..n).map(|p| leaving[p] / area[p]).fold(0.0, f64::max);
            let subs = ((dt * fastest / MOISTURE_CFL).ceil() as usize).max(1);
            let h = dt / subs as f64;
            let mut change = vec![0.0; n];
            for _ in 0..subs {
                change.fill(0.0);
                for p in 0..n {
                    let (i, j) = (p / n_lon, p % n_lon);
                    let mut carry = |from: usize, to: usize, f: f64| {
                        let m = if f > 0.0 { f * w[from] } else { f * w[to] };
                        change[from] -= m;
                        change[to] += m;
                    };
                    if n_lon > 1 {
                        carry(p, i * n_lon + (j + 1) % n_lon, east[p]);
                    }
                    if i + 1 < n_lat {
                        carry(p, p + n_lon, north[p]);
                    }
                }
                for p in 0..n {
                    w[p] = (w[p] + h * change[p] / area[p]).max(0.0);
                }
            }
            // The eddies mix it, implicitly.
            diag.fill(1.0);
            rhs.copy_from_slice(&w);
            solve(&mix, &diag, &rhs, &mut w);
            // The ground gives it and it rains, together and implicitly, on the step's end:
            // w' + dt (takes w' + P(w')) = w + dt gives.
            for p in 0..n {
                // The root lies between no rain and all of the column: Newton's steps kept
                // within the bracket, halved where they leave it; the rain the column's balance.
                let start = w[p].max(0.0) + dt * gives[p];
                let balance = |x: f64| x + dt * (takes[p] * x + air.rain(x, saturated[p])) - start;
                let (mut lo, mut hi) = (0.0, start.max(0.0));
                let mut x = w[p].clamp(lo, hi);
                for _ in 0..MOISTURE_NEWTON {
                    let f = balance(x);
                    if f > 0.0 {
                        hi = x;
                    } else {
                        lo = x;
                    }
                    let rate = air.rain(x, saturated[p]);
                    let rising = RAIN_SATURATION_SLOPE
                        * (x / saturated[p].max(1e-30) - RAIN_SATURATION_ONSET)
                        < MOISTURE_MOST_RISE;
                    let slope = 1.0
                        + dt * takes[p]
                        + if rising {
                            dt * rate * RAIN_SATURATION_SLOPE / saturated[p].max(1e-30)
                        } else {
                            0.0
                        };
                    let mut next = x - f / slope;
                    if !(next > lo && next < hi) {
                        next = 0.5 * (lo + hi);
                    }
                    let settled = (next - x).abs() <= 1e-12 * start.max(1e-30);
                    x = next;
                    if settled || hi - lo <= 1e-12 * start.max(1e-30) {
                        break;
                    }
                }
                let fell = ((start - x) / dt - takes[p] * x).max(0.0);
                let given = gives[p] - takes[p] * x;
                w[p] = x;
                let r = x / saturated[p].max(1e-30);
                // The land's bucket: the rain on it, less what it gives the air, over its
                // land; what overflows runs off.
                let mut ran = 0.0;
                if land[p] > 0.0 {
                    let dried = exchange[p] * beta[p] * (soil_q[p] - r * humid[p]);
                    bucket[p] = (bucket[p] + dt * (fell - dried)).max(0.0);
                    if bucket[p] > capacity {
                        ran = (bucket[p] - capacity) / dt;
                        bucket[p] = capacity;
                    }
                }
                if last {
                    out.falling[k * n + p] = fell as f32;
                    out.snow[k * n + p] = air.snow_share(ground[p]) as f32;
                    out.precipitation[p] += fell / steps as f64;
                    out.evaporation[p] += given / steps as f64;
                    out.runoff[p] += ran / steps as f64;
                    out.saturation[p] += r / steps as f64;
                    out.wetness[p] += beta[p] / steps as f64;
                    // The potential evaporation: the open sea's and the land's wet.
                    let potential = exchange[p]
                        * (sites[p].sea * (sea_q[p] - r * humid[p]).max(0.0)
                            + land[p] * (soil_q[p] - r * humid[p]).max(0.0));
                    out.potential[p] += potential / steps as f64;
                }
            }
        }
    }
    Some(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn earth_air() -> Air {
        Air {
            vapour: volatile("H2O").unwrap(),
            pressure: 101_325.0,
            molar_mass: 0.028_97,
            heat_capacity: 1_004.0,
            scale_height: 8_400.0,
            gravity: 9.81,
            radius: 6.371e6,
            lapse: 9.81 / 1_004.0 * MOIST_LAPSE_FACTOR,
            transport: 0.58,
            water: OCEAN_TRANSPORT_SHARE,
        }
    }

    /// F268: the saturated column over Earth's 288 K holds some 30 kg/m², near the 25 kg/m²
    /// its air holds at its mean saturation of some 0.8, and over its tropical oceans' 300 K
    /// some 65; at the onset of Bretherton et al.'s fit a tropical column rains 1 mm a day.
    #[test]
    fn a_column_saturates_and_rains_as_measured() {
        let air = earth_air();
        let mild = air.saturated_column(288.0);
        let warm = air.saturated_column(300.0);
        assert!((25.0..35.0).contains(&mild), "{mild}");
        assert!((55.0..75.0).contains(&warm), "{warm}");
        let onset = air.rain(RAIN_SATURATION_ONSET * 65.0, 65.0) * SECONDS_PER_DAY;
        assert!((onset - 1.0).abs() < 1e-9, "{onset}");
        assert!(air.rain(0.8 * 65.0, 65.0) > 10.0 * air.rain(0.6 * 65.0, 65.0));
        assert!(air.snow_share(260.0) > 0.99 && air.snow_share(280.0) < 0.01);
    }
}
