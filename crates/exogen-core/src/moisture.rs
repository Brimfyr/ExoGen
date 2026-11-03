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
//!
//! Over the relief at a finer grid's scale (F270), air carried over rising ground condenses
//! its vapour at C_w w, w the lift, the wind against the slope, C_w the saturated vapour's
//! density at the ground times the moist adiabat's lapse rate over the air's (Smith and
//! Barstad 2004), while it is saturated: for the share of each step the column rains at its
//! saturation, its rain over the law's at saturation, read at the ground's own height, where
//! the column carried up to it stands the nearer saturation, and all of the step where it
//! would saturate. The water falls once it has turned to rain and fallen, carried downwind
//! meanwhile, and descending air evaporates it again and dries the rain that falls through
//! it. The air carried on is short of what the relief rained beyond the places' rain, and
//! rains the less downwind, by the law's rise, until the rain it holds back has made it good:
//! the rain shadows.

use crate::constants::*;
use crate::fields::{Field, det_sum};
use crate::geom::{dot, local_basis, normalize};
use crate::grid::projection::{PaddedField, sample_bilinear};
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

    /// The moist adiabat's lapse rate at `t` under `p`, K/m: the dry adiabat's, g / c_p, as
    /// much less as the latent heat rising saturated air gives up warms it, (1 + L q / (R_d
    /// T)) / (1 + L² q / (c_p R_v T²)), q the saturated specific humidity.
    fn moist_lapse(&self, t: f64, p: f64) -> f64 {
        let t = t.max(1.0);
        let (q, l) = (self.saturated_humidity(t, p), self.latent(t));
        let dry = GAS_CONSTANT / self.molar_mass;
        let gives = 1.0 + l * q / (dry * t);
        let warms = 1.0 + l * l * q / (self.heat_capacity * self.vapour.r_specific * t * t);
        self.gravity / self.heat_capacity * gives / warms
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
    /// The year's lifting wind at each place (F270) over ground at each of LIFT_LEVELS
    /// heights over the place's, from LIFT_LOWEST of its vapour's scale heights LIFT_STEP
    /// apart, place after place, level after level within each: its wind east and north,
    /// m/s, times C_w, the vapour saturated air condenses for each metre it is lifted, kg/m³,
    /// times the share of the step the air there is saturated; and the same without the wind,
    /// kg/m³, over which the lifting wind is the wind its rain drifts on.
    pub lift: Vec<f64>,
    /// Each place's ground over the sea, m, its vapour's scale height through the year, R_v
    /// T² / (L Γ), m, and its saturated column through the year, kg/m² (F270).
    pub heights: Vec<f64>,
    pub vapour_height: Vec<f64>,
    pub saturated: Vec<f64>,
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
    /// ARIDITY_MAP_MOST, over the relief `above`, the ground's height over the sea, m, of a
    /// body of `radius` m (F270), the potential evaporation read between the places.
    pub fn fields(&self, above: &Field<f32>, radius: f64) -> [Field<f32>; 2] {
        let year = SECONDS_PER_DAY * 365.25;
        let n_lon = self.layout.1;
        let mut rain = self.over_relief(above, radius);
        let dry = Field::from_fn(above.grid(), |c, d| {
            let (lon, lat) = lon_lat_from_direction(d);
            let potential = match n_lon {
                1 => read_row(&self.potential, lat),
                _ => read_grid(&self.potential, n_lon, lat, lon),
            };
            let index = if potential > 0.0 {
                rain.get(c) as f64 / potential
            } else {
                f64::INFINITY
            };
            index.min(ARIDITY_MAP_MOST) as f32
        });
        rain.map_in_place(|v| (v as f64 * year) as f32);
        [rain, dry]
    }

    /// The year's precipitation over the relief `above`, the ground's height over the sea on a
    /// grid finer than the places', m, of a body of `radius` m, kg/m²/s (F270). Air carried
    /// over rising ground condenses at C_w w, w the lift, the wind against the slope; the
    /// condensation is as linear in the wind as in the slope, so that the year's is the
    /// places' lifting wind, read at the ground's height over each place's, against the
    /// ground's slope. Its water falls once it has turned to rain and fallen, carried on the
    /// rain's wind meanwhile, the share fallen by a drift s 1 - (1 + s/d) e^(-s/d), d the
    /// drift through each delay; descending air takes what the rain over it would give, and
    /// where it would take more the rain that falls through it dries by e^(-taken/rain). The
    /// air carried on is short of what the relief rained beyond the places' rain upwind of
    /// it, and rains the less downwind until it has made it good; the body's rain is kept.
    pub fn over_relief(&self, above: &Field<f32>, radius: f64) -> Field<f32> {
        let grid = above.grid();
        let (n_lat, n_lon) = self.layout;
        let cell = grid.delta() * radius;
        let corners = |d: [f64; 3]| {
            let (lon, lat) = lon_lat_from_direction(d);
            crate::climate::reader(n_lat, n_lon, lat, lon)
        };
        let read = |values: &[f64], at: &[(usize, f64); 4]| -> f64 {
            at.iter().map(|&(p, w)| w * values[p]).sum()
        };
        // The lifting wind over ground `h` m over the sea, east and north, kg/m²/s for each
        // unit of slope, and its weight, kg/m³: each place's at the ground's height over its
        // own in its vapour's scale heights, between its levels.
        let lifting = |at: &[(usize, f64); 4], h: f64| -> [f64; 3] {
            let mut out = [0.0; 3];
            for &(p, w) in at.iter().filter(|c| c.1 > 0.0) {
                let u = (h - self.heights[p]) / self.vapour_height[p].max(1.0);
                let x = ((u - LIFT_LOWEST) / LIFT_STEP).clamp(0.0, (LIFT_LEVELS - 1) as f64);
                let j = (x as usize).min(LIFT_LEVELS - 2);
                let f = x - j as f64;
                let base = (p * LIFT_LEVELS + j) * 3;
                for (c, o) in out.iter_mut().enumerate() {
                    *o += w * (self.lift[base + c] * (1.0 - f) + self.lift[base + 3 + c] * f);
                }
            }
            out
        };
        // The lift's condensation at each cell: the lifting wind against the slope, the slope
        // by least squares over the cell's four neighbours.
        let source = Field::from_fn(grid, |c, d| {
            let (east, north) = local_basis(d);
            let h0 = above.get(c) as f64;
            let (mut xx, mut xy, mut yy, mut xh, mut yh) = (0.0, 0.0, 0.0, 0.0, 0.0);
            for (di, dj) in [(1, 0), (-1, 0), (0, 1), (0, -1)] {
                let nb = grid.neighbor(c, di, dj);
                let p = grid.center(nb);
                let dp = [
                    (p[0] - d[0]) * radius,
                    (p[1] - d[1]) * radius,
                    (p[2] - d[2]) * radius,
                ];
                let (x, y) = (dot(dp, east), dot(dp, north));
                let dh = above.get(nb) as f64 - h0;
                xx += x * x;
                xy += x * y;
                yy += y * y;
                xh += x * dh;
                yh += y * dh;
            }
            let det = xx * yy - xy * xy;
            if det <= 0.0 {
                return 0.0;
            }
            let (gx, gy) = ((yy * xh - xy * yh) / det, (xx * yh - xy * xh) / det);
            let [east, north, _] = lifting(&corners(d), h0);
            (east * gx + north * gy) as f32
        });
        // Carried on the rain's wind, the lifting wind over its weight, through its two delays
        // before it falls.
        let rain_wind = |at: &[(usize, f64); 4], h: f64| -> (f64, f64) {
            let [east, north, weight] = lifting(at, h);
            if weight > 0.0 {
                (east / weight, north / weight)
            } else {
                (0.0, 0.0)
            }
        };
        let padded = PaddedField::new(&grid, source.as_slice());
        let fallen = |s: f64, drift: f64| 1.0 - (1.0 + s / drift) * (-s / drift).exp();
        let lifted = Field::from_fn(grid, |c, d| {
            let (ue, un) = rain_wind(&corners(d), above.get(c) as f64);
            let speed = ue.hypot(un);
            let drift = speed * OROGRAPHIC_DELAY_S;
            if drift <= 0.0 || fallen(0.5 * cell, drift) >= OROGRAPHIC_KERNEL_TAIL {
                return source.get(c);
            }
            let (east, north) = local_basis(d);
            let e: [f64; 3] = std::array::from_fn(|i| (ue * east[i] + un * north[i]) / speed);
            let (mut sum, mut done) = (0.0, 0.0);
            for j in 0..OROGRAPHIC_MOST_STEPS {
                let reach = fallen((j as f64 + 0.5) * cell, drift);
                let back = j as f64 * cell / radius;
                let p = normalize(std::array::from_fn(|i| d[i] - back * e[i]));
                sum += (reach - done) * sample_bilinear(&grid, &padded, p) as f64;
                done = reach;
                if done >= OROGRAPHIC_KERNEL_TAIL {
                    break;
                }
            }
            (sum / done) as f32
        });
        // The rain over the relief on the rain `left` it would have: the more by what the
        // lift brings down, and where descent would take more than falls, the less by
        // e^(-taken/rain).
        let over = |left: f64, lifted: f64| {
            if lifted >= 0.0 {
                left + lifted
            } else if left > 0.0 {
                left * (lifted / left).exp()
            } else {
                0.0
            }
        };
        // The air carried on the rain's wind is short of what the relief rained beyond the
        // places' rain upwind of it, and rains the less for it by the law's rise, e^(-a D /
        // W), D the shortfall, W the saturated column, a RAIN_SATURATION_SLOPE, until the rain
        // it holds back has made it good: dD/dt = gained - P (1 - e^(-a D / W)), on a coarser
        // grid of the shortfall's reach, traced upwind along the rain's wind over
        // OROGRAPHIC_SHADOW_STEPS of its cells and carried down it a step at a time,
        // implicitly.
        let shadow = CubeSphere::new(OROGRAPHIC_SHADOW_FACE.min(grid.face_size()));
        let gained = Field::from_fn(grid, |c, d| {
            let background = read(&self.precipitation, &corners(d));
            (over(background, lifted.get(c) as f64) - background) as f32
        });
        let gained = crate::climate::coarsen(&gained, shadow);
        let gained = PaddedField::new(&shadow, gained.as_slice());
        let level = crate::climate::coarsen(above, shadow);
        let level = PaddedField::new(&shadow, level.as_slice());
        let step = shadow.delta() * radius;
        let short = Field::from_fn(shadow, |_, d| {
            let mut path = Vec::with_capacity(OROGRAPHIC_SHADOW_STEPS);
            let mut p = d;
            for _ in 0..OROGRAPHIC_SHADOW_STEPS {
                let (ue, un) = rain_wind(&corners(p), sample_bilinear(&shadow, &level, p) as f64);
                let speed = ue.hypot(un);
                if speed <= 0.0 {
                    break;
                }
                let (east, north) = local_basis(p);
                let e: [f64; 3] = std::array::from_fn(|i| (ue * east[i] + un * north[i]) / speed);
                p = normalize(std::array::from_fn(|i| p[i] - step / radius * e[i]));
                let at = corners(p);
                let (rain, held) = (read(&self.precipitation, &at), read(&self.saturated, &at));
                let gain = sample_bilinear(&shadow, &gained, p) as f64;
                path.push((gain, rain, held.max(1e-30), step / speed));
            }
            let mut lost = 0.0f64;
            for &(gain, rain, held, dt) in path.iter().rev() {
                let rise = RAIN_SATURATION_SLOPE / held;
                let mut x = (lost + dt * gain).max(0.0);
                for _ in 0..OROGRAPHIC_SHADOW_NEWTON {
                    let back = (-rise * x).exp();
                    let f = x - lost - dt * gain + dt * rain * (1.0 - back);
                    x = (x - f / (1.0 + dt * rain * rise * back)).max(0.0);
                }
                lost = x;
            }
            lost as f32
        });
        let short = PaddedField::new(&shadow, short.as_slice());
        let mut rain = Field::from_fn(grid, |c, d| {
            let at = corners(d);
            let background = read(&self.precipitation, &at);
            let held = read(&self.saturated, &at).max(1e-30);
            let lost = (sample_bilinear(&shadow, &short, d) as f64).max(0.0);
            let left = background * (-RAIN_SATURATION_SLOPE * lost / held).exp();
            over(left, lifted.get(c) as f64) as f32
        });
        // The body's rain kept, whatever of the shortfall its reach leaves unmade.
        let n = grid.cell_count();
        let want = det_sum(n, |k| {
            let c = grid.cell(k);
            read(&self.precipitation, &corners(grid.center(c))) * grid.cell_area(c)
        });
        let got = det_sum(n, |k| {
            rain.as_slice()[k] as f64 * grid.cell_area(grid.cell(k))
        });
        let keep = if got > 0.0 { want / got } else { 1.0 };
        rain.map_in_place(|v| (v as f64 * keep) as f32);
        rain
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
    // Each place's vapour scale height through the year, R_v T² / (L Γ), the measure of the
    // lift's levels (F270).
    let lapse = air.lapse.max(MOISTURE_LEAST_LAPSE);
    let scale = |t: f64| air.vapour.r_specific * t * t / (air.latent(t) * lapse);
    let mut vapour_height = vec![0.0; n];
    for k in 0..steps {
        let now = weather(k);
        for (p, v) in vapour_height.iter_mut().enumerate() {
            *v += scale((now.air[p] - air.lapse * heights[p]).max(1.0)) / steps as f64;
        }
    }
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
        lift: vec![0.0; n * LIFT_LEVELS * 3],
        heights: heights.to_vec(),
        vapour_height,
        saturated: vec![0.0; n],
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
                    out.saturated[p] += saturated[p] / steps as f64;
                    out.wetness[p] += beta[p] / steps as f64;
                    // The potential evaporation: the open sea's and the land's wet.
                    let potential = exchange[p]
                        * (sites[p].sea * (sea_q[p] - r * humid[p]).max(0.0)
                            + land[p] * (soil_q[p] - r * humid[p]).max(0.0));
                    out.potential[p] += potential / steps as f64;
                    // The lift's condensation through the step at each level of the ground over
                    // the place's (F270): C_w, the saturated vapour's density at the ground
                    // times the moist adiabat's lapse rate over the air's, the less by e^(-dz/H)
                    // in the colder air over higher ground, while the air there is saturated:
                    // for the share of the step its column, the more saturated by e^(dz/H),
                    // rains at its saturation, all of it once it would saturate.
                    let t = ground[p].max(1.0);
                    let density = vapour_pressure(air.vapour, t) / (air.vapour.r_specific * t);
                    let c_w = density * air.moist_lapse(t, pressure[p]) / lapse;
                    let measure = out.vapour_height[p] / scale(t);
                    for j in 0..LIFT_LEVELS {
                        let dz = (LIFT_LOWEST + LIFT_STEP * j as f64) * measure;
                        let share = (RAIN_SATURATION_SLOPE * (r * dz.exp() - 1.0))
                            .min(0.0)
                            .exp();
                        let lift = share * c_w * (-dz).exp() / steps as f64;
                        let at = (p * LIFT_LEVELS + j) * 3;
                        out.lift[at] += lift * now.east[p] as f64;
                        out.lift[at + 1] += lift * now.north[p] as f64;
                        out.lift[at + 2] += lift;
                    }
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

    /// F270: air carried east over a ridge along a meridian rains on the ridge's west, its
    /// windward side, more than on its east, its lee, and beyond the lee the less for what the
    /// ridge took from it; the rain over the body is kept; and in calm air the rain is the
    /// places' own, read between them.
    #[test]
    fn a_ridge_rains_on_its_windward_side() {
        let (n_lat, n_lon) = (13, 16);
        let n = n_lat * n_lon;
        let day = 1.0 / SECONDS_PER_DAY;
        // Some 1 mm a day everywhere, more toward the equator, under a column that rains at
        // its saturation a twentieth of the time, C_w some 0.007 kg/m³, the wind 10 m/s east.
        let lat = |p: usize| -PI / 2.0 + PI * (p / n_lon) as f64 / (n_lat - 1) as f64;
        let moisture = |wind: f64| Moisture {
            layout: (n_lat, n_lon),
            steps: 1,
            falling: vec![0.0; n],
            snow: vec![0.0; n],
            precipitation: (0..n).map(|p| day * (1.0 + lat(p).cos())).collect(),
            evaporation: vec![day; n],
            potential: vec![2.0 * day; n],
            runoff: vec![0.0; n],
            saturation: vec![0.75; n],
            wetness: vec![0.5; n],
            lift: (0..n * LIFT_LEVELS * 3)
                .map(|i| {
                    let u = LIFT_LOWEST + LIFT_STEP * ((i / 3) % LIFT_LEVELS) as f64;
                    let lift = 0.05 * 0.007 * (-u).exp();
                    [lift * wind, 0.0, lift][i % 3]
                })
                .collect(),
            heights: vec![0.0; n],
            vapour_height: vec![2_400.0; n],
            saturated: vec![25.0; n],
            weights: vec![1.0 / n as f64; n],
            land: vec![1.0; n],
        };
        let grid = CubeSphere::new(64);
        let radius = 6.371e6;
        // A ridge 2 km high and some 300 km wide along the meridian at 30° east.
        let ridge = Field::from_fn(grid, |_, d| {
            let (lon, _) = lon_lat_from_direction(d);
            let x = (lon - 30f64.to_radians()) / 3f64.to_radians();
            (2_000.0 * (-x * x).exp()) as f32
        });
        let total = |f: &Field<f32>| -> f64 {
            (0..grid.cell_count())
                .map(|k| f.as_slice()[k] as f64 * grid.cell_area(grid.cell(k)))
                .sum()
        };
        let calm = moisture(0.0);
        let still = calm.over_relief(&ridge, radius);
        let windy = moisture(10.0).over_relief(&ridge, radius);
        let at = |f: &Field<f32>, lon: f64| {
            let d = crate::grid::direction_from_lon_lat(lon.to_radians(), 0.0);
            f.get(grid.locate(d).0) as f64
        };
        let (windward, lee, beside) = (at(&windy, 27.0), at(&windy, 33.0), at(&windy, -90.0));
        println!("windward {windward:.3e}, lee {lee:.3e}, beside {beside:.3e} kg/m²/s");
        assert!(
            windward > 1.3 * beside && lee < 0.8 * beside,
            "{windward} {lee} {beside}"
        );
        let (a, b) = (total(&windy), total(&still));
        assert!((a / b - 1.0).abs() < 1e-6, "{a} {b}");
        // A step 2 km up at 30° east: its edge rains on the air its descent no longer gives
        // back, and the plateau beyond rains the less for it.
        let step = Field::from_fn(grid, |_, d| {
            let (lon, _) = lon_lat_from_direction(d);
            let x = (lon - 30f64.to_radians()) / 2f64.to_radians();
            (1_000.0 * (1.0 + x.tanh())) as f32
        });
        let stepped = moisture(10.0).over_relief(&step, radius);
        let (edge, plateau) = (at(&stepped, 29.0), at(&stepped, 40.0));
        let beside = at(&stepped, -90.0);
        println!("edge {edge:.3e}, plateau {plateau:.3e}, beside {beside:.3e} kg/m²/s");
        assert!(
            edge > 1.3 * beside && plateau < 0.97 * beside,
            "{edge} {plateau} {beside}"
        );
        // Calm air: the places' rain read between them.
        for k in (0..grid.cell_count()).step_by(97) {
            let d = grid.center(grid.cell(k));
            let (lon, lat) = lon_lat_from_direction(d);
            let read = read_grid(&calm.precipitation, n_lon, lat, lon);
            let got = still.as_slice()[k] as f64;
            assert!((got / read - 1.0).abs() < 1e-5, "{got} {read}");
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
