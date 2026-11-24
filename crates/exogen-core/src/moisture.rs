//! The moisture (F268).
//!
//! The working liquid's vapour through the cycle's year at its places: a column of vapour over
//! each place, evaporated from its open sea and its wet land at its own wind (F267), carried by
//! the wind and mixed by the eddies that carry the heat, and rained out as the column nears
//! saturation.
//!
//! The wind carries it in two parts (M74). The air that flows along the isobars carries its
//! whole column. The drift across them is the low branch of an overturning, whose air aloft
//! returns it, and its vapour is counted with the overturning's heat (M75): the heat an
//! overturning carries between two places, the transport's share that is not the eddies', is
//! the mass it turns over times its gross stability, what its air aloft holds over its air
//! below, and that mass carries the colder place's surface air's vapour to the warmer (Held
//! 2001; Siler et al. 2018): the trades' toward the equator, a sea's toward a land warmer
//! than it, a night's toward its day. Its air aloft holds 1 + CELL_STABILITY of the energy,
//! c_p T + L q, of the air that rises, the air where it rains. The eddies carry
//! the vapour as the latent share of the heat they carry (Frierson et al. 2007): the air they
//! exchange between two places carries its own heat, c_p a kelvin, and its vapour's, L dq/dT,
//! so that the transport D moves D / (c_p + L dq/dT) kg of air a kelvin, and with it the
//! difference of the two places' surface air's vapour; and they carry the share of the
//! transport that is theirs: all of it from the overturning cell's edge on, and less within
//! the cell by degrees to none at the thermal equator, where the cell's own overturning
//! carries the heat aloft and its vapour below (Siler et al. 2018).
//!
//! The open sea evaporates into its surface air, which stands at CLOUD_SURFACE_HUMIDITY of
//! saturation whatever its column holds (M74), as Earth's seas' does within a few hundredths
//! from 60 south to 60 north, the air the sea's clouds and the column's convection already
//! read; the wet land into air at its column's saturation, as before.
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
//! runoff, and of what falls on it the share that falls on ground already full, the square of
//! its fill (M74), which holds a year's balance on Budyko's (1974) curve. The sea under its
//! ice and the caps' ice give none. What a wet ground would give, the land's potential
//! evaporation, is Penman's (1948): the heat its ground has for its air, sensible and
//! latent, shared as a wet surface shares it, and its air's drying power (M75); a dry
//! ground's own heat is no measure of it.
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
    /// How far the overturning cell reaches from the equator, rad (F32), within which the
    /// eddies carry the less of the transport (M74) and the overturning the more (M75).
    pub reach: f64,
    /// The time a column takes to rain its saturated vapour at the law's onset, s: Earth's,
    /// RAIN_TIMESCALE_DAYS, over the share of Earth's measured sunlight its ground absorbs
    /// (F431), set by the year that carries it.
    pub rain_time: f64,
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
        let circulation = crate::wind::Air::from_graph(graph);
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
            reach: circulation.reach,
            rain_time: RAIN_TIMESCALE_DAYS * SECONDS_PER_DAY,
        })
    }

    /// The air's rain time over ground that absorbs `sunlight` W/m² over the body through its
    /// year (F431): Earth's, RAIN_TIMESCALE_DAYS, over its share of Earth's measured
    /// RAIN_SUNLIGHT_REFERENCE. The sunlight its ground absorbs is what its ground gives its
    /// air to drive the convection that rains a column, its evaporation's bound; none, no
    /// rain but near saturation.
    pub fn rain_time_under(sunlight: f64) -> f64 {
        RAIN_TIMESCALE_DAYS * SECONDS_PER_DAY * RAIN_SUNLIGHT_REFERENCE
            / sunlight.max(RAIN_SUNLIGHT_LEAST)
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
    /// kg/m²/s (Bretherton et al. 2004), over the air's rain time (F431).
    pub fn rain(&self, w: f64, saturated: f64) -> f64 {
        if saturated <= 0.0 || w <= 0.0 {
            return 0.0;
        }
        self.rain_over(w, saturated, self.rain_time)
    }

    /// The precipitation from a column of `w` kg/m² whose saturation holds `saturated`,
    /// kg/m²/s, over the rain time `time`, s.
    pub fn rain_over(&self, w: f64, saturated: f64, time: f64) -> f64 {
        if saturated <= 0.0 || w <= 0.0 {
            return 0.0;
        }
        let r = w / saturated;
        let rise = (RAIN_SATURATION_SLOPE * (r - RAIN_SATURATION_ONSET)).min(MOISTURE_MOST_RISE);
        saturated / time * rise.exp()
    }

    /// The share of the precipitation at `t` that falls as snow: half at a kelvin above the
    /// liquid's freezing point, the rest in a logistic a kelvin wide (Jennings et al. 2018).
    fn snow_share(&self, t: f64) -> f64 {
        let x = (t - self.vapour.t_triple - RAIN_SNOW_OFFSET) / RAIN_SNOW_WIDTH;
        1.0 / (1.0 + x.clamp(-50.0, 50.0).exp())
    }

    /// The vapour's scale height over ground whose air stands at `t`, m: R_v T² / (L Γ).
    pub fn vapour_height(&self, t: f64) -> f64 {
        let t = t.max(1.0);
        self.vapour.r_specific * t * t / (self.latent(t) * self.lapse.max(MOISTURE_LEAST_LAPSE))
    }

    /// The energy of a kilogram of surface air at `t` K holding `vapour` kg of vapour, J/kg
    /// (M75): c_p T + L q, its moist static energy at the ground.
    pub fn energy(&self, t: f64, vapour: f64) -> f64 {
        self.heat_capacity * t + self.latent(t) * vapour
    }

    /// An overturning's gross stability, J/kg (M75): what its air aloft holds over its air
    /// below, its air aloft with 1 + CELL_STABILITY of the energy `rising` of the air that
    /// rises, its air below with the energy `low` of the place it leaves; no less than at
    /// the rising end, where a place holds more than the rising air (Siler et al. 2018).
    pub fn cell_stability(&self, rising: f64, low: f64) -> f64 {
        ((1.0 + CELL_STABILITY) * rising - low).max(CELL_STABILITY * rising)
    }

    /// The overturning's share of the air's transport between the latitudes `a` and `b`,
    /// rad, the thermal equator at `equator` (M75): what the eddies leave of it (M74).
    pub fn cell_share(&self, a: f64, b: f64, equator: f64) -> f64 {
        1.0 - 0.5 * (self.eddies_share(a, equator) + self.eddies_share(b, equator))
    }

    /// What a wet ground would give its air, kg/m²/s (M75; Penman 1948): the heat the
    /// ground has for its air, `sensible` and `latent`, both as the vapour their heat would
    /// evaporate, kg/m²/s, shared between the two as a wet surface shares it, by the slope
    /// of the saturated humidity at the air's temperature `t` K and pressure `p` Pa against
    /// c_p / L; and the air's drying power `power`, kg/m²/s, the exchange times its deficit
    /// under saturation. A wet ground's own evaporation, where the ground is wet; none where
    /// the ground takes heat from its air.
    pub fn penman(&self, (t, p): (f64, f64), sensible: f64, latent: f64, power: f64) -> f64 {
        let (heat, drying) = self.penman_shares(t, p);
        (heat * (sensible + latent) + drying * power).max(0.0)
    }

    /// The shares Penman's (1948) wet surface gives the heat its ground has for its air and
    /// its air's drying power at the air's temperature `t` K and pressure `p` Pa: the slope of
    /// the saturated humidity and c_p / L, each over their sum.
    fn penman_shares(&self, t: f64, p: f64) -> (f64, f64) {
        let t = t.max(1.0);
        let l = self.latent(t);
        let slope = self.saturated_humidity(t, p) * l / (self.vapour.r_specific * t * t);
        let psychro = self.heat_capacity / l;
        (slope / (slope + psychro), psychro / (slope + psychro))
    }

    /// The latent heat a kilogram of air at `t` under `p` Pa and saturation `r` carries a
    /// kelvin beside its own (M74), J/kg/K: L dq/dT along the places' air, q its vapour, r
    /// q_s, and dq_s/dT by Clausius and Clapeyron, q_s L / (R_v T²).
    pub fn latent_capacity(&self, t: f64, p: f64, r: f64) -> f64 {
        let t = t.max(1.0);
        let l = self.latent(t);
        r.clamp(0.0, 1.0) * self.saturated_humidity(t, p) * l * l / (self.vapour.r_specific * t * t)
    }

    /// The share of the air's transport its eddies carry at latitude `lat`, rad, the thermal
    /// equator at `equator` (M74): all of it from the overturning cell's edge on, which
    /// stands at the cell's reach from the equator, and within the cell less by degrees to
    /// none at the thermal equator, sin²(pi/2 x) of the way x from the one to the other,
    /// level at both: 0.22, 0.67 and 0.98 at 10, 20 and 30 degrees of Earth's, where Siler et
    /// al. (2018) fit ERA's with 1 - exp(-(sin lat / 0.3)²), 0.29, 0.73 and 0.94. All of it
    /// without a cell.
    pub fn eddies_share(&self, lat: f64, equator: f64) -> f64 {
        if self.reach <= 0.0 {
            return 1.0;
        }
        let side = if lat >= equator { 1.0 } else { -1.0 };
        let span = (side * self.reach - equator).abs().max(1e-6);
        let x = ((lat - equator).abs() / span).clamp(0.0, 1.0);
        (0.5 * PI * x).sin().powi(2)
    }
}

/// The share of its full rate a land's bucket evaporates at, `fill` of its capacity full
/// (F268): all of it from BUCKET_FULL_SHARE, in proportion below (Manabe 1969).
pub fn bucket_gives(fill: f64) -> f64 {
    (fill / BUCKET_FULL_SHARE).clamp(0.0, 1.0)
}

/// The share of the rain and the melt on a land that runs off it, its bucket `fill` of its
/// capacity full (M74): the share of its ground already full, its fill to BUCKET_SHED.
pub fn bucket_sheds(fill: f64) -> f64 {
    fill.clamp(0.0, 1.0).powf(BUCKET_SHED)
}

/// The implicit mixing's solver (F268): (diag + links) x = rhs, the links each place's
/// neighbours and its coupling to each, x the start on entry.
pub type Solve<'s> = dyn Fn(&[Vec<(usize, f64)>], &[f64], &[f64], &mut [f64]) + 's;

/// The weather the moisture reads at a step (F268): the air over each place at the datum, K,
/// the wind, m/s, the land's surface and the sea's at the datum, K, the sea's ice, kg/m², the
/// year's snowpack on the land, its snow at the step's end and what its melt took through the
/// step, kg/m² of the land, NaN where its year's snow is unknown and none without the year's
/// snowpack (F277), and the gusts about the wind.
pub struct Weather<'a> {
    pub air: &'a [f64],
    pub east: &'a [f32],
    pub north: &'a [f32],
    pub land: Vec<f64>,
    pub sea: Vec<f64>,
    pub ice: Vec<f64>,
    pub snow: Option<(Vec<f64>, Vec<f64>)>,
    /// The body's mean wind at the step, m/s, the weather's gusts about each place's mean
    /// wind.
    pub gust: f64,
    /// The part of the wind that crosses the isobars, the drift of the layer the ground's
    /// drag acts over, eastward and northward, m/s (M74), which carries no column: its
    /// vapour is its overturning's (M75); and the thermal equator, rad.
    pub drift: (&'a [f32], &'a [f32]),
    pub equator: f64,
}

/// The moisture through the cycle's year (F268).
#[derive(Clone, Debug)]
pub struct Moisture {
    /// The layout's latitudes and longitudes, and the cycle's steps.
    pub layout: (usize, usize),
    pub steps: usize,
    /// The precipitation at each place through each step, kg/m²/s, step after step, place
    /// after place within each, and the share of it that falls as snow; and its column's
    /// saturation at each (M66), laid out alike.
    pub falling: Vec<f32>,
    pub snow: Vec<f32>,
    pub saturations: Vec<f32>,
    /// The year's mean at each place: its precipitation and evaporation over the place, its
    /// potential evaporation, of the open sea over its sea and of wet ground over its land,
    /// and its land's runoff over its land, kg/m²/s; and its column's mean saturation.
    pub precipitation: Vec<f64>,
    pub evaporation: Vec<f64>,
    pub potential: Vec<f64>,
    pub runoff: Vec<f64>,
    pub saturation: Vec<f64>,
    /// Each place's vapour column through the year, kg/m² (M69).
    pub column: Vec<f64>,
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
    /// The air's rain time its year ran under, s (F431).
    pub rain_time: f64,
    /// The share of each place's land its lakes cover through the year where its runoff
    /// reaches no sea (F433), none where it does; and what spills into them from its
    /// neighbours' lakes less what spills out of them, kg/m²/s over its land.
    pub lakes: Vec<f64>,
    pub spill: Vec<f64>,
}

/// The year's moisture in brief (F268), kg/m² a year: precipitation and evaporation over the
/// body; over its land, precipitation, evaporation and runoff; and the shares of the land whose
/// precipitation over its potential evaporation stands under 0.2, arid, and under 0.05,
/// hyper-arid (UNEP 1992); its air's rain time, days (F431), and the share of its land its
/// lakes cover in the year (F433).
#[derive(Clone, Copy, Debug)]
pub struct Summary {
    pub precipitation: f64,
    pub evaporation: f64,
    pub land_precipitation: f64,
    pub land_evaporation: f64,
    pub runoff: f64,
    pub arid: f64,
    pub hyper_arid: f64,
    pub rain_time_days: f64,
    pub lakes: f64,
}

impl Moisture {
    /// The year in brief.
    pub fn summary(&self) -> Summary {
        let year = SECONDS_PER_DAY * 365.25;
        let (mut whole, mut land) = (0.0, 0.0);
        let (mut p, mut e, mut lp, mut le, mut ro) = (0.0, 0.0, 0.0, 0.0, 0.0);
        let (mut arid, mut hyper, mut lakes) = (0.0, 0.0, 0.0);
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
                lakes += l * self.lakes.get(k).copied().unwrap_or(0.0);
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
            rain_time_days: self.rain_time / SECONDS_PER_DAY,
            lakes: if land > 0.0 { lakes / land } else { 0.0 },
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
    /// body of `radius` m (F270), the potential evaporation read between the places; `above`
    /// the shape the grid resolves, the cells about each cell its samples (F422).
    pub fn fields(&self, above: &Field<f32>, radius: f64) -> [Field<f32>; 2] {
        let year = SECONDS_PER_DAY * 365.25;
        let n_lon = self.layout.1;
        let mut rain = self.over_relief(above, radius);
        let dry = Field::from_fn(above.grid(), |c, d| {
            let (lon, lat) = lon_lat_from_direction(d);
            let potential =
                crate::insolation::smooth_grid(&self.potential, n_lon, lat, lon).max(0.0);
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

    /// The land's water over the relief (F277), kg/m²/s at each cell of the year's precipitation
    /// `rain` there, kg/m² a year (`fields`): what its land sheds, the rain less its places'
    /// land's evaporation read between them, none where that takes all of it, the evaporation
    /// spread from the places with land over the sea's so that a coast reads its land's; and
    /// what an open water's surface gains there, the rain less the potential evaporation read
    /// between the places.
    pub fn land_water(&self, rain: &Field<f32>) -> [Field<f32>; 2] {
        let year = SECONDS_PER_DAY * 365.25;
        let n_lon = self.layout.1;
        let dries: Vec<f64> = (0..self.precipitation.len())
            .map(|k| self.land_evaporation(k))
            .collect();
        let own: Vec<bool> = self.land.iter().map(|&l| l > 0.0).collect();
        let dries = crate::climate::spread_from(&dries, &own, self.layout);
        let read = |table: &[f64], d: [f64; 3]| {
            let (lon, lat) = lon_lat_from_direction(d);
            crate::insolation::smooth_grid(table, n_lon, lat, lon).max(0.0)
        };
        let shed = Field::from_fn(rain.grid(), |c, d| {
            (rain.get(c) as f64 / year - read(&dries, d)).max(0.0) as f32
        });
        let gain = Field::from_fn(rain.grid(), |c, d| {
            (rain.get(c) as f64 / year - read(&self.potential, d)) as f32
        });
        [shed, gain]
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
        // The places' rain read smoothly between them (F273).
        let background = |d: [f64; 3]| {
            let (lon, lat) = lon_lat_from_direction(d);
            crate::insolation::smooth_grid(&self.precipitation, n_lon, lat, lon).max(0.0)
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
            let background = background(d);
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
            let background = background(d);
            let held = read(&self.saturated, &at).max(1e-30);
            let lost = (sample_bilinear(&shadow, &short, d) as f64).max(0.0);
            let left = background * (-RAIN_SATURATION_SLOPE * lost / held).exp();
            over(left, lifted.get(c) as f64) as f32
        });
        // The body's rain kept, whatever of the shortfall its reach leaves unmade.
        let n = grid.cell_count();
        let want = det_sum(n, |k| {
            let c = grid.cell(k);
            background(grid.center(c)) * grid.cell_area(c)
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
    // The air's transport through each link over the step, J/m²/K of the place: the heat
    // the eddies carry between the two a kelvin where all of the transport is theirs.
    let carried: Vec<Vec<(usize, f64)>> = crate::climate::place_links((n_lat, n_lon))
        .iter()
        .enumerate()
        .map(|(p, links)| {
            links
                .iter()
                .map(|&(q, g)| {
                    let wet = sites[p].sea.min(sites[q].sea);
                    (q, dt * air.transport * g * (1.0 - air.water * wet))
                })
                .collect()
        })
        .collect();
    let mut mix = carried.clone();
    let lat_of = |p: usize| -PI / 2.0 + dphi * (p / n_lon) as f64;
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
    // At least the years of the spin-up, then on until the year balances (M63), within the
    // steps the run may take.
    let most = (MOISTURE_MOST_STEPS / steps).max(MOISTURE_LEAST_YEARS);
    let years = ((MOISTURE_SPIN_UP / year).ceil() as usize)
        .max(MOISTURE_LEAST_YEARS)
        .min(most);
    // Each place's vapour scale height through the year, R_v T² / (L Γ), the measure of the
    // lift's levels (F270).
    let lapse = air.lapse.max(MOISTURE_LEAST_LAPSE);
    let marine = CLOUD_SURFACE_HUMIDITY;
    let scale = |t: f64| air.vapour_height(t);
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
        saturations: vec![0.0; steps * n],
        precipitation: vec![0.0; n],
        evaporation: vec![0.0; n],
        potential: vec![0.0; n],
        runoff: vec![0.0; n],
        saturation: vec![0.0; n],
        column: vec![0.0; n],
        wetness: vec![0.0; n],
        lift: vec![0.0; n * LIFT_LEVELS * 3],
        heights: heights.to_vec(),
        vapour_height,
        saturated: vec![0.0; n],
        weights: sites.iter().map(|s| s.weight).collect(),
        land: land.clone(),
        rain_time: air.rain_time,
        lakes: vec![0.0; n],
        spill: vec![0.0; n],
    };
    let (mut diag, mut rhs) = (vec![0.0; n], vec![0.0; n]);
    let (mut gives, mut takes) = (vec![0.0; n], vec![0.0; n]);
    let (mut exchange, mut beta) = (vec![0.0; n], vec![0.0; n]);
    // What each bucket gives through the step, its share of the wet land's, no more than it
    // holds (F433): beta's, the ground's wetness, as its heat for its air reads it.
    let mut draw = vec![0.0; n];
    let mut store = vec![0.0; n];
    // The surface layer over each place's sea and its land (M63), its exchange over its
    // neutral.
    let layer = crate::surface_layer::SurfaceLayer {
        gravity: air.gravity,
        heat_capacity: air.heat_capacity,
    };
    let (mut over_sea, mut over_land) = (vec![1.0; n], vec![1.0; n]);
    let (mut humid, mut sea_q, mut soil_q) = (vec![0.0; n], vec![0.0; n], vec![0.0; n]);
    let (mut east, mut north) = (vec![0.0; n], vec![0.0; n]);
    let (mut heats, mut eddies) = (vec![0.0; n], vec![0.0; n]);
    let mut carries = vec![0.0; n];
    // The rain each place took the step before, where the overturning's air rises, and the
    // air the overturning turns over through each link, kg a step, from its cold place to
    // its warm (M75).
    let mut rained_last = vec![0.0; n];
    let mut flows: Vec<(usize, usize, f64)> = Vec::new();
    let mut last = years <= 1;
    let rain_time = air.rain_time;
    // The lakes (F433): where no place holds a sea its runoff reaches, the share of each
    // place's land its lakes cover, and the year's sums that set it.
    let closed = !sites.iter().any(|s| s.sea > 0.0) && land.iter().any(|&l| l > 0.0);
    let mut lake = vec![0.0; n];
    let mut kept = vec![0.0; n];
    let mut sums = vec![[0.0f64; 3]; n];
    let mut pond = vec![(0.0f64, 0.0f64); n];
    for round in 0.. {
        // The body's rain and evaporation through the year, kg/s.
        let (mut rained, mut gave) = (0.0, 0.0);
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
                beta[p] = bucket_gives(bucket[p] / capacity);
                draw[p] = beta[p];
                humid[p] = air.saturated_humidity(ground[p], pressure[p]);
                sea_q[p] = air.saturated_humidity(now.sea[p], air.pressure);
                let soil = now.land[p] - air.lapse * heights[p];
                soil_q[p] = air.saturated_humidity(soil, pressure[p]);
                // The exchange the surface layer's stability lets over the sea and over the
                // land, each surface against the air over it (M63).
                if sites[p].sea > 0.0 {
                    over_sea[p] = layer.exchange(now.sea[p], now.air[p], speed);
                }
                pond[p] = (0.0, 0.0);
                if land[p] > 0.0 {
                    over_land[p] = layer.exchange(soil, ground[p], speed);
                    // A bucket gives no more through the step than it holds as the step
                    // begins (F433), into air however dry, where its floor gave back what it
                    // lacked.
                    let r = w[p] / saturated[p].max(1e-30);
                    let x = exchange[p] * over_land[p];
                    let most = x * soil_q[p];
                    if most > 0.0 {
                        draw[p] = beta[p].min(bucket[p] / (most * dt));
                    }
                    // Its lakes give their air what a wet ground would (F433), Penman's: the
                    // heat its ground has for its air and its air's drying power shared as a
                    // wet surface shares them, the land's potential (M75), as the finish's
                    // lakes give it; linear in the column's saturation, a - b r, none where
                    // that gives none as the step begins.
                    if closed {
                        let l = air.latent(ground[p].max(1.0));
                        let sensible = air.heat_capacity / l * x * (soil - ground[p]);
                        let (heat, drying) = air.penman_shares(ground[p], pressure[p]);
                        let a = heat * (sensible + x * beta[p] * soil_q[p]) + drying * x * humid[p];
                        let b = (heat * beta[p] + drying) * x * humid[p];
                        if a - b * r > 0.0 {
                            pond[p] = (a, b);
                        }
                    }
                }
                // The open sea into its surface air, at CLOUD_SURFACE_HUMIDITY of saturation
                // whatever its column holds, none into air that holds more than the sea's
                // own; the wet land into air at its column's saturation (M74), its buckets'
                // share of it; its lakes theirs (F433).
                let open = sites[p].sea * (-now.ice[p] / SEA_ICE_COVER_MASS).exp() * over_sea[p];
                let wet = land[p] * (1.0 - lake[p]) * draw[p] * over_land[p];
                let sea = (sea_q[p] - marine * humid[p]).max(0.0);
                let ponds = land[p] * lake[p];
                gives[p] = exchange[p] * (open * sea + wet * soil_q[p]) + ponds * pond[p].0;
                takes[p] =
                    (exchange[p] * wet * humid[p] + ponds * pond[p].1) / saturated[p].max(1e-30);
                // Its surface air's vapour a kilogram of its column's (M74).
                carries[p] = humid[p] / saturated[p].max(1e-30);
            }
            // The wind carries the vapour across the places' faces, upwind, in steps short
            // enough that no place gives more than it holds: the air along the isobars its
            // column (M74); the drift across them none, its vapour its overturning's (M75).
            for p in 0..n {
                let (i, j) = (p / n_lon, p % n_lon);
                let (lo, hi) = edges(i);
                east[p] = if n_lon > 1 {
                    let e = i * n_lon + (j + 1) % n_lon;
                    let across = 0.5 * (now.drift.0[p] + now.drift.0[e]) as f64;
                    let all = 0.5 * (now.east[p] + now.east[e]) as f64;
                    (all - across) * a * (hi - lo)
                } else {
                    0.0
                };
                north[p] = if i + 1 < n_lat {
                    let across = 0.5 * (now.drift.1[p] + now.drift.1[p + n_lon]) as f64;
                    let all = 0.5 * (now.north[p] + now.north[p + n_lon]) as f64;
                    (all - across) * a * hi.cos() * around
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
                    let q = p + n_lon;
                    add(p, north[p]);
                    add(q, -north[p]);
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
                        let up = if f > 0.0 { from } else { to };
                        let m = f * w[up];
                        change[from] -= m;
                        change[to] += m;
                    };
                    if n_lon > 1 {
                        let e = i * n_lon + (j + 1) % n_lon;
                        carry(p, e, east[p]);
                    }
                    if i + 1 < n_lat {
                        carry(p, p + n_lon, north[p]);
                    }
                }
                for p in 0..n {
                    w[p] = (w[p] + h * change[p] / area[p]).max(0.0);
                }
            }
            // The overturning carries it (M75). The heat an overturning carries through a
            // link, the air's transport's share that is not the eddies', is the mass it
            // turns over times its gross stability; that mass carries the cold place's
            // surface air's vapour to the warm place, in steps short enough that no place
            // gives more than it holds. The air that rises is the air where it rains: its
            // energy the rain's own places', the step before; the greatest before any rain.
            {
                let energy = |p: usize, w: &[f64]| air.energy(now.air[p], w[p] * carries[p]);
                let (mut held, mut fell) = (0.0, 0.0);
                for p in 0..n {
                    let r = rained_last[p] * area[p];
                    held += r * energy(p, &w);
                    fell += r;
                }
                let rising = if fell > 0.0 {
                    held / fell
                } else {
                    (0..n).map(|p| energy(p, &w)).fold(f64::MIN, f64::max)
                };
                flows.clear();
                let mut leaving = vec![0.0; n];
                for (p, links) in carried.iter().enumerate() {
                    for &(q, heat) in links {
                        if q <= p {
                            continue;
                        }
                        let (warm, cold) = if now.air[p] >= now.air[q] {
                            (p, q)
                        } else {
                            (q, p)
                        };
                        let share = air.cell_share(lat_of(p), lat_of(q), now.equator);
                        let joules = heat * area[p] * share * (now.air[warm] - now.air[cold]);
                        let mass = joules / air.cell_stability(rising, energy(cold, &w));
                        if mass > 0.0 {
                            leaving[cold] += mass * carries[cold];
                            flows.push((cold, warm, mass));
                        }
                    }
                }
                let fastest = (0..n).map(|p| leaving[p] / area[p]).fold(0.0, f64::max);
                let subs = ((fastest / MOISTURE_CFL).ceil() as usize).max(1);
                let mut change = vec![0.0; n];
                for _ in 0..subs {
                    change.fill(0.0);
                    for &(cold, warm, mass) in &flows {
                        let m = mass / subs as f64 * w[cold] * carries[cold];
                        change[cold] -= m;
                        change[warm] += m;
                    }
                    for p in 0..n {
                        w[p] = (w[p] + change[p] / area[p]).max(0.0);
                    }
                }
            }
            // The eddies mix it, implicitly (M74): the air they exchange between two places,
            // the heat they carry over what a kilogram of it carries a kelvin, its own and
            // its vapour's, carries the difference of the two places' surface air's vapour,
            // x, each column's w its surface air's over `carries`: w / carries + sum k (x -
            // x_q) = w; and theirs is the share of the transport beyond the thermal equator.
            for p in 0..n {
                let r = w[p] / saturated[p].max(1e-30);
                heats[p] = air.heat_capacity + air.latent_capacity(ground[p], pressure[p], r);
                eddies[p] = air.eddies_share(lat_of(p), now.equator);
            }
            for (p, (links, steady)) in mix.iter_mut().zip(&carried).enumerate() {
                for (link, &(q, heat)) in links.iter_mut().zip(steady) {
                    let share = 0.5 * (eddies[p] + eddies[q]);
                    link.1 = heat * share / (0.5 * (heats[p] + heats[q]));
                }
            }
            for p in 0..n {
                diag[p] = 1.0 / carries[p].max(1e-30);
                rhs[p] = w[p];
                w[p] *= carries[p];
            }
            solve(&mix, &diag, &rhs, &mut w);
            for p in 0..n {
                w[p] = (w[p] / carries[p].max(1e-30)).max(0.0);
            }
            // The ground gives it and it rains, together and implicitly, on the step's end:
            // w' + dt (takes w' + P(w')) = w + dt gives.
            for p in 0..n {
                // The root lies between no rain and all of the column: Newton's steps kept
                // within the bracket, halved where they leave it; the rain the column's balance.
                let start = w[p].max(0.0) + dt * gives[p];
                let balance = |x: f64| {
                    x + dt * (takes[p] * x + air.rain_over(x, saturated[p], rain_time)) - start
                };
                let (mut lo, mut hi) = (0.0, start.max(0.0));
                let mut x = w[p].clamp(lo, hi);
                for _ in 0..MOISTURE_NEWTON {
                    let f = balance(x);
                    if f > 0.0 {
                        hi = x;
                    } else {
                        lo = x;
                    }
                    let rate = air.rain_over(x, saturated[p], rain_time);
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
                rained_last[p] = fell;
                let given = gives[p] - takes[p] * x;
                rained += fell * area[p];
                gave += given * area[p];
                w[p] = x;
                let r = x / saturated[p].max(1e-30);
                // The land's bucket: the rain on it, less what it gives the air, over its
                // land; what overflows runs off. Where the year's snowpack is known its snow
                // lies on the land until it melts (F277), in the step by the share of what lay on
                // the snowpack that its melt took, by the air's warmth where none lay, and no
                // deeper than SNOW_MOST_MASS: the bucket takes the rain and the melt.
                let mut ran = 0.0;
                if land[p] > 0.0 {
                    let deficit = soil_q[p] - r * humid[p];
                    let dried = exchange[p] * over_land[p] * draw[p] * deficit;
                    let pack = now
                        .snow
                        .as_ref()
                        .map(|(lying, melted)| (lying[p], melted[p]));
                    let wet = match pack.filter(|&(_, melted)| melted.is_finite()) {
                        Some((lying, melted)) => {
                            let cold = air.snow_share(ground[p]);
                            let share = if lying + melted > 0.0 {
                                melted / (lying + melted)
                            } else {
                                1.0 - cold
                            };
                            store[p] = (store[p] + dt * fell * cold).min(SNOW_MOST_MASS);
                            let melts = store[p] * share.clamp(0.0, 1.0);
                            store[p] -= melts;
                            fell * (1.0 - cold) + melts / dt
                        }
                        None => fell,
                    };
                    // What falls on its ground already full runs off (M74); the bucket takes
                    // the rest, and sheds what overflows it.
                    ran = wet * bucket_sheds(bucket[p] / capacity);
                    bucket[p] = (bucket[p] + dt * (wet - ran - dried)).max(0.0);
                    if bucket[p] > capacity {
                        ran += (bucket[p] - capacity) / dt;
                        bucket[p] = capacity;
                    }
                    // Where its runoff reaches no sea, the lakes' year (F433): the rain on
                    // its land, what its buckets shed and what its lakes would give the air.
                    if closed {
                        let sum = &mut sums[p];
                        sum[0] += wet;
                        sum[1] += ran;
                        sum[2] += (pond[p].0 - pond[p].1 * r).max(0.0);
                    }
                }
                if last {
                    out.falling[k * n + p] = fell as f32;
                    out.snow[k * n + p] = air.snow_share(ground[p]) as f32;
                    out.saturations[k * n + p] = r as f32;
                    out.precipitation[p] += fell / steps as f64;
                    out.evaporation[p] += given / steps as f64;
                    out.runoff[p] += (1.0 - lake[p]) * ran / steps as f64;
                    out.saturation[p] += r / steps as f64;
                    out.column[p] += x / steps as f64;
                    out.saturated[p] += saturated[p] / steps as f64;
                    out.wetness[p] += (lake[p] + (1.0 - lake[p]) * draw[p]) / steps as f64;
                    // The potential evaporation: the open sea's, and what the land would
                    // give wet (M75): Penman's, by the heat its ground has for its air, as
                    // the vapour it would evaporate, and its air's drying power.
                    let wet_land = {
                        let x = exchange[p] * over_land[p];
                        let soil = now.land[p] - air.lapse * heights[p];
                        let l = air.latent(ground[p].max(1.0));
                        let sensible = air.heat_capacity / l * x * (soil - ground[p]);
                        let latent = x * beta[p] * (soil_q[p] - r * humid[p]);
                        let power = x * (1.0 - r).max(0.0) * humid[p];
                        air.penman((ground[p], pressure[p]), sensible, latent, power)
                    };
                    let potential = exchange[p]
                        * sites[p].sea
                        * over_sea[p]
                        * (sea_q[p] - marine * humid[p]).max(0.0)
                        + land[p] * wet_land;
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
        if last {
            out.lakes = lake;
            out.spill = kept;
            break;
        }
        // Each place's lakes cover the share of its land that gives its air what its buckets
        // shed into them and what spills into them, less the rain on them, f A (E - P) =
        // (1 - f) A R + I (F433), A its land; a place whose lakes would gain all lake, its
        // gain spilling to its land's neighbours, on until lakes that can give it take it;
        // what finds none is kept out of the air's cycle, until its air rains no more.
        let mut moved = 0.0f64;
        if closed {
            let year: Vec<[f64; 3]> = sums.iter().map(|s| s.map(|x| x / steps as f64)).collect();
            let reach: Vec<f64> = (0..n).map(|p| land[p] * area[p]).collect();
            // What a place's lakes gain over what they give, kg/s, all lake under its
            // inflow `into`.
            let gain = |p: usize, into: f64| {
                let [fell, _, gives] = year[p];
                into - reach[p] * (gives - fell)
            };
            let neighbours: Vec<Vec<usize>> = carried
                .iter()
                .map(|links| {
                    links
                        .iter()
                        .map(|&(q, _)| q)
                        .filter(|&q| land[q] > 0.0)
                        .collect()
                })
                .collect();
            let whole: f64 = (0..n).map(|p| reach[p] * year[p][0]).sum();
            let mut inflow = vec![0.0; n];
            let mut passed = vec![0.0; n];
            for _ in 0..LAKE_SPILL_PASSES {
                let mut spilt = false;
                for p in (0..n).filter(|&p| reach[p] > 0.0) {
                    let fresh = gain(p, inflow[p]).max(0.0) - passed[p];
                    if fresh > LAKE_SPILL_LEAST * whole && !neighbours[p].is_empty() {
                        let each = fresh / neighbours[p].len() as f64;
                        for &q in &neighbours[p] {
                            inflow[q] += each;
                        }
                        passed[p] += fresh;
                        spilt = true;
                    }
                }
                if !spilt {
                    break;
                }
            }
            for p in 0..n {
                let next = if reach[p] > 0.0 {
                    let [fell, shed, gives] = year[p];
                    let net = gives - fell;
                    let into = (inflow[p] - passed[p]) / reach[p];
                    kept[p] = into;
                    if net > 0.0 {
                        ((shed + into) / (shed + net)).clamp(0.0, 1.0)
                    } else {
                        1.0
                    }
                } else {
                    0.0
                };
                moved = moved.max((next - lake[p]).abs());
                lake[p] = next;
            }
            sums.iter_mut().for_each(|s| *s = [0.0; 3]);
        }
        // The year the moisture keeps once the spin-up's years have run, the column neither
        // fills nor empties over the year and its lakes hold, or the run's steps are spent.
        let settled = (rained - gave).abs() <= MOISTURE_SETTLED * rained.abs().max(gave.abs());
        last = (round + 2 >= years && settled && moved <= MOISTURE_SETTLED) || round + 2 >= most;
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
            reach: 32f64.to_radians(),
            rain_time: RAIN_TIMESCALE_DAYS * SECONDS_PER_DAY,
        }
    }

    /// F270: air carried east over a ridge along a meridian rains on the ridge's west, its
    /// windward side, more than on its east, its lee, and beyond the lee the less for what the
    /// ridge took from it; the rain over the body is kept; and in calm air the rain is the
    /// places' own, read smoothly between them.
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
            saturations: vec![0.75; n],
            precipitation: (0..n).map(|p| day * (1.0 + lat(p).cos())).collect(),
            evaporation: vec![day; n],
            potential: vec![2.0 * day; n],
            runoff: vec![0.0; n],
            saturation: vec![0.75; n],
            column: vec![0.0; n],
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
            rain_time: RAIN_TIMESCALE_DAYS * SECONDS_PER_DAY,
            lakes: vec![0.0; n],
            spill: vec![0.0; n],
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
        // Calm air: the places' rain read smoothly between them (F273).
        for k in (0..grid.cell_count()).step_by(97) {
            let d = grid.center(grid.cell(k));
            let (lon, lat) = lon_lat_from_direction(d);
            let read = crate::insolation::smooth_grid(&calm.precipitation, n_lon, lat, lon);
            let got = still.as_slice()[k] as f64;
            assert!((got / read - 1.0).abs() < 1e-5, "{got} {read}");
        }
    }

    /// M74: a land's bucket through a year of steady rain holds its balance on Budyko's (1974)
    /// curve: what it evaporates of its rain, by its wet ground's evaporation over its rain,
    /// within 0.06 of (x tanh(1 / x) (1 - exp(-x)))^(1/2) from a half to 8, the rest shed;
    /// where the bucket shed only what overflowed it, it evaporated all of its rain from x of
    /// 1 on, 0.31 over the curve.
    #[test]
    fn a_bucket_sheds_as_budyko_s_catchments() {
        let budyko = |x: f64| (x * (1.0 / x).tanh() * (1.0 - (-x).exp())).sqrt();
        assert!((budyko(1.0) - 0.694).abs() < 1e-3 && (budyko(2.0) - 0.894).abs() < 1e-3);
        for x in [0.5f64, 0.75, 1.0, 1.5, 2.0, 4.0, 8.0] {
            // Days of rain, 1 a day, its wet ground's evaporation x a day, the bucket 50
            // days' rain deep.
            let (mut fill, capacity, mut gave) = (0.5f64, 50.0f64, 0.0);
            for day in 0..20_000 {
                let ran = bucket_sheds(fill);
                let dried = x * bucket_gives(fill);
                fill = (fill + (1.0 - ran - dried) / capacity).clamp(0.0, 1.0);
                if day >= 19_000 {
                    gave += dried / 1_000.0;
                }
            }
            let curve = budyko(x);
            println!("x {x}: evaporates {gave:.3} of its rain, Budyko's {curve:.3}");
            assert!((gave - curve).abs() < 0.06, "{x}: {gave} {curve}");
        }
        assert_eq!(bucket_sheds(0.0), 0.0);
        assert_eq!(bucket_sheds(1.0), 1.0);
        assert!((bucket_sheds(0.5) - 0.25).abs() < 1e-12);
        assert!(bucket_gives(0.75) == 1.0 && (bucket_gives(0.375) - 0.5).abs() < 1e-12);
    }

    /// M75: an overturning's gross stability is what its air aloft holds over its air below.
    /// Earth's tropical surface air, 300 K at 0.8 of saturation, holds some 345 kJ/kg, a
    /// seventh of it its vapour's; where it rises its air aloft holds CELL_STABILITY of that
    /// more, some 21 kJ/kg (Siler et al. 2018), and over air 10 K colder as much again as
    /// that air holds less; no less than at the rising end over air that holds more. The
    /// overturning's share of the transport is what the eddies leave of it: all of it at the
    /// thermal equator, some four fifths at 10 degrees, none from the cell's edge on or
    /// without a cell.
    #[test]
    fn an_overturning_s_stability_is_its_air_aloft_s_over_its_air_below() {
        let air = earth_air();
        let humid = |t: f64| CLOUD_SURFACE_HUMIDITY * air.saturated_humidity(t, air.pressure);
        let rising = air.energy(300.0, humid(300.0));
        assert!((338e3..352e3).contains(&rising), "{rising}");
        let vapour = air.latent(300.0) * humid(300.0);
        assert!(
            (0.1..0.17).contains(&(vapour / rising)),
            "{vapour} {rising}"
        );
        let at_rise = air.cell_stability(rising, rising);
        assert!((at_rise - CELL_STABILITY * rising).abs() < 1e-6);
        assert!((19e3..23e3).contains(&at_rise), "{at_rise}");
        let low = air.energy(290.0, humid(290.0));
        let over = air.cell_stability(rising, low);
        assert!((over - (at_rise + rising - low)).abs() < 1e-6);
        assert!(over > 2.0 * at_rise, "{over} {at_rise}");
        // Air that holds more than the rising air: the rising end's stability.
        assert_eq!(air.cell_stability(rising, 1.2 * rising), at_rise);
        // The overturning's share: what the eddies leave.
        let d = |x: f64| x.to_radians();
        assert_eq!(air.cell_share(0.0, 0.0, 0.0), 1.0);
        let ten = air.cell_share(d(10.0), d(10.0), 0.0);
        assert!((0.74..0.82).contains(&ten), "{ten}");
        assert!((ten + air.eddies_share(d(10.0), 0.0) - 1.0).abs() < 1e-12);
        let across = air.cell_share(d(8.0), d(12.0), 0.0);
        let mean = 0.5 * (air.eddies_share(d(8.0), 0.0) + air.eddies_share(d(12.0), 0.0));
        assert!((across + mean - 1.0).abs() < 1e-12);
        assert_eq!(air.cell_share(d(40.0), d(44.0), 0.0), 0.0);
        let none = Air { reach: 0.0, ..air };
        assert_eq!(none.cell_share(0.0, d(4.0), 0.0), 0.0);
    }

    /// M75: what a wet ground would give is Penman's (1948). A wet ground's own evaporation
    /// is its potential: with the heat it gives its air as a ground 2 K over its air gives
    /// it, Penman's returns the evaporation the exchange draws from it, within its
    /// linearisation. A dry ground 15 K over its air gives its air nearly all of its heat as
    /// warmth; wet it would give most of that heat as vapour, far less than the exchange
    /// would draw from a wet surface as hot as the dry ground stands, which no wet ground
    /// would be. None where the ground takes heat from its air and the air is saturated.
    #[test]
    fn a_wet_ground_s_potential_is_penman_s() {
        let air = earth_air();
        let (t, p, r) = (298.0, air.pressure, 0.6);
        let q = |t: f64| air.saturated_humidity(t, p);
        let psychro = air.heat_capacity / air.latent(t);
        // The exchange, kg/m²/s: the air's density, the coefficient and a 5 m/s wind.
        let x = 1.2 * AIR_EXCHANGE_COEFFICIENT * 5.0;
        let power = x * (1.0 - r) * q(t);
        // A wet ground 2 K over its air.
        let wet = x * (q(t + 2.0) - r * q(t));
        let sensible = psychro * x * 2.0;
        let potential = air.penman((t, p), sensible, wet, power);
        assert!((potential / wet - 1.0).abs() < 0.03, "{potential} {wet}");
        // A dry ground 15 K over its air, a twentieth wet.
        let hot = x * (q(t + 15.0) - r * q(t));
        let (sensible, latent) = (psychro * x * 15.0, 0.05 * hot);
        let potential = air.penman((t, p), sensible, latent, power);
        let slope = q(t) * air.latent(t) / (air.vapour.r_specific * t * t);
        let by_hand = (slope * (sensible + latent) + psychro * power) / (slope + psychro);
        assert!(
            (potential / by_hand - 1.0).abs() < 1e-12,
            "{potential} {by_hand}"
        );
        println!(
            "a dry ground's potential {:.2} mm a day, where a wet surface as hot would give {:.2}",
            potential * SECONDS_PER_DAY,
            hot * SECONDS_PER_DAY
        );
        assert!(potential < 0.5 * hot, "{potential} {hot}");
        // No more than the heat it has for its air and its air's drying power together.
        assert!(potential < sensible + latent + power);
        assert!(potential > latent, "{potential} {latent}");
        // A ground colder than a saturated air: none.
        assert_eq!(air.penman((t, p), -psychro * x * 3.0, 0.0, 0.0), 0.0);
    }

    /// M74: the air the eddies exchange carries its vapour's latent heat beside its own:
    /// about as much again at 285 K, a latent share of some half of the heat they carry, as
    /// Earth's atmosphere carries some 0.4 of its heat as vapour near 40 degrees (Trenberth
    /// and Stepaniak 2003); under a sixth at 250 K and some three quarters at 300 K; none in
    /// dry air and in proportion to its saturation. Their share of the transport rises from
    /// none at the thermal equator to all of it at the cell's edge, within 0.08 of Siler et
    /// al.'s (2018) fit to ERA's at 10, 20 and 30 degrees; the cell's edges stand at its reach
    /// from the equator wherever the thermal equator; all of it without a cell.
    #[test]
    fn the_eddies_carry_vapour_as_their_heat_s_latent_share() {
        let air = earth_air();
        let p = air.pressure;
        let share = |t: f64, r: f64| {
            let l = air.latent_capacity(t, p, r);
            l / (air.heat_capacity + l)
        };
        let (cold, mild, warm) = (share(250.0, 0.8), share(285.0, 0.8), share(300.0, 0.8));
        assert!((0.45..0.60).contains(&mild), "{mild}");
        assert!(cold < 1.0 / 6.0 && cold > 0.03, "{cold}");
        assert!((0.65..0.80).contains(&warm), "{warm}");
        assert_eq!(air.latent_capacity(285.0, p, 0.0), 0.0);
        let (half, whole) = (
            air.latent_capacity(285.0, p, 0.4),
            air.latent_capacity(285.0, p, 0.8),
        );
        assert!((whole / half - 2.0).abs() < 1e-12, "{whole} {half}");
        // By hand: r q_s L² / (R_v T²).
        let (t, v) = (285.0f64, air.vapour);
        let by_hand =
            0.8 * air.saturated_humidity(t, p) * v.l_vap * v.l_vap / (v.r_specific * t * t);
        assert!((whole / by_hand - 1.0).abs() < 1e-12, "{whole} {by_hand}");
        // The eddies' share of the transport.
        let reach = air.reach;
        assert_eq!(air.eddies_share(0.0, 0.0), 0.0);
        for (lat, siler) in [(10.0f64, 0.285), (20.0, 0.727), (30.0, 0.938)] {
            let got = air.eddies_share(lat.to_radians(), 0.0);
            let fit = 1.0 - (-(lat.to_radians().sin() / 0.3).powi(2)).exp();
            assert!((fit - siler).abs() < 1e-3, "{fit}");
            assert!((got - siler).abs() < 0.08, "{lat}: {got} {siler}");
            assert!((air.eddies_share(-lat.to_radians(), 0.0) - got).abs() < 1e-12);
        }
        assert!((air.eddies_share(reach, 0.0) - 1.0).abs() < 1e-12);
        assert!((air.eddies_share(1.3, 0.0) - 1.0).abs() < 1e-12);
        let equator = 10f64.to_radians();
        assert_eq!(air.eddies_share(equator, equator), 0.0);
        assert!((air.eddies_share(reach, equator) - 1.0).abs() < 1e-12);
        assert!((air.eddies_share(-reach, equator) - 1.0).abs() < 1e-12);
        // Ten degrees toward the farther edge, less of the way to it.
        let (far, near) = (
            air.eddies_share(0.0, equator),
            air.eddies_share(2.0 * equator, equator),
        );
        assert!(far > 0.0 && far < near && near < 1.0, "{far} {near}");
        let none = Air { reach: 0.0, ..air };
        assert_eq!(none.eddies_share(0.3, 0.0), 1.0);
    }

    /// The year's water over an aquaplanet of 46 latitudes, warm at its equator, its trades
    /// drifting toward it across its isobars at up to `trades` m/s within the cell's reach.
    fn aquaplanet(air: &Air, trades: f64) -> Moisture {
        let reach = 32f64.to_radians();
        strip(
            air,
            |_| (1.0, 0.0),
            |lat| {
                if lat.abs() < reach {
                    -trades * (PI * lat / reach).sin()
                } else {
                    0.0
                }
            },
            true,
        )
    }

    /// A year's moisture over a strip of latitudes under a steady weather (M74): each
    /// latitude's sea's and ice's shares from `ground`, its northward wind, m/s, from `wind`:
    /// a drift across its isobars where `drifts`, along them where not.
    fn strip(
        air: &Air,
        ground: impl Fn(f64) -> (f64, f64),
        wind: impl Fn(f64) -> f64,
        drifts: bool,
    ) -> Moisture {
        strip_over(air, ground, wind, drifts, 1.0)
    }

    /// The strip over a year `years` Earth years long, its 360 steps the longer.
    fn strip_over(
        air: &Air,
        ground: impl Fn(f64) -> (f64, f64),
        wind: impl Fn(f64) -> f64,
        drifts: bool,
        years: f64,
    ) -> Moisture {
        let (n_lat, steps) = (46usize, 360usize);
        let dphi = PI / (n_lat - 1) as f64;
        let sites: Vec<Site> = (0..n_lat)
            .map(|i| {
                let lat = -PI / 2.0 + dphi * i as f64;
                let (lo, hi) = (
                    (lat - dphi / 2.0).max(-PI / 2.0),
                    (lat + dphi / 2.0).min(PI / 2.0),
                );
                let (sea, ice) = ground(lat);
                Site {
                    lat,
                    lon: 0.0,
                    weight: 0.5 * (hi.sin() - lo.sin()),
                    sea,
                    ice,
                }
            })
            .collect();
        let temps: Vec<f64> = sites
            .iter()
            .map(|s| 300.0 - 45.0 * s.lat.sin().powi(2))
            .collect();
        let seas: Vec<f64> = temps.iter().map(|t| t + 1.0).collect();
        let north: Vec<f32> = sites.iter().map(|s| wind(s.lat) as f32).collect();
        let east = vec![0.0f32; n_lat];
        let weather = |_: usize| Weather {
            air: &temps,
            east: &east,
            north: &north,
            land: seas.clone(),
            sea: seas.clone(),
            ice: vec![0.0; n_lat],
            snow: None,
            gust: 6.0,
            drift: (&east, if drifts { &north } else { &east }),
            equator: 0.0,
        };
        // (diag + links) x = rhs along the latitudes, by Thomas's algorithm.
        let solve = |links: &[Vec<(usize, f64)>], diag: &[f64], rhs: &[f64], x: &mut [f64]| {
            let n = diag.len();
            let link = |p: usize, q: usize| links[p].iter().find(|l| l.0 == q).map_or(0.0, |l| l.1);
            let mut b: Vec<f64> = (0..n)
                .map(|p| diag[p] + links[p].iter().map(|l| l.1).sum::<f64>())
                .collect();
            let (mut c, mut r) = (vec![0.0; n], rhs.to_vec());
            for p in 0..n {
                if p + 1 < n {
                    c[p] = -link(p, p + 1);
                }
                if p > 0 {
                    let m = -link(p, p - 1) / b[p - 1];
                    b[p] -= m * c[p - 1];
                    r[p] -= m * r[p - 1];
                }
            }
            x[n - 1] = r[n - 1] / b[n - 1];
            for p in (0..n - 1).rev() {
                x[p] = (r[p] - c[p] * x[p + 1]) / b[p];
            }
        };
        let year = years * 365.25 * SECONDS_PER_DAY;
        let heights = vec![0.0; n_lat];
        rain(
            air,
            &sites,
            (n_lat, 1),
            &heights,
            (year, steps),
            &weather,
            &solve,
        )
        .expect("a year")
    }

    /// F431: a column rains its saturated vapour over Earth's time where its ground absorbs
    /// Earth's measured sunlight, the longer as its ground absorbs the less: Titan's some
    /// 0.35 W/m² some 460 times; and the law's rise its own at any time.
    #[test]
    fn a_column_rains_by_the_sunlight_its_ground_absorbs() {
        let earth = Air::rain_time_under(RAIN_SUNLIGHT_REFERENCE);
        assert!((earth / (RAIN_TIMESCALE_DAYS * SECONDS_PER_DAY) - 1.0).abs() < 1e-12);
        let titan = Air::rain_time_under(0.35);
        assert!(
            (titan / earth - 161.0 / 0.35).abs() < 1e-9,
            "{}",
            titan / earth
        );
        assert!(Air::rain_time_under(0.0).is_finite());
        let air = earth_air();
        for (w, saturated) in [(30.0, 65.0), (5.0, 8.0), (4_000.0, 8_000.0)] {
            let rain = air.rain_over(w, saturated, titan);
            assert!((rain / air.rain_over(w, saturated, earth) - earth / titan).abs() < 1e-12);
            assert!(
                (air.rain(w, saturated) / air.rain_over(w, saturated, earth) - 1.0).abs() < 1e-12
            );
        }
    }

    /// F433: where no sea takes its runoff, a land's runoff stands in its lakes and its lakes
    /// give it back to the air, no more than its rain over the year at any place: a strip of
    /// land under Earth's air, its year balanced, its lakes standing where its buckets shed
    /// and under its rain where they would gain, each a share of its land; where its buckets'
    /// floor gave back the runoff, the land gave its air its rain and its runoff besides.
    /// Over a strip with a sea, none.
    #[test]
    fn a_land_s_runoff_stands_in_its_lakes_where_no_sea_takes_it() {
        let air = earth_air();
        // Its year Earth's, its steps a day; and thirty years long, its steps a month, as
        // coarse as Titan's, where a bucket's step could take more than it holds.
        for years in [1.0, 30.0] {
            let dry = strip_over(&air, |_| (0.0, 0.0), |_| 0.0, false, years);
            let brief = dry.summary();
            println!("a strip of land over {years} years: {brief:?}");
            let (p, e) = (brief.precipitation, brief.evaporation);
            assert!(
                p > 0.0 && (p - e).abs() <= 2.0 * MOISTURE_SETTLED * p.max(e),
                "{p} {e}"
            );
            assert!(dry.lakes.iter().all(|&f| (0.0..=1.0).contains(&f)));
            assert!(dry.lakes.iter().any(|&f| f > 0.0), "{:?}", dry.lakes);
            // Each place gives its air its rain and what spills in, less what spills out,
            // neither more, liquid from nothing, nor less, liquid lost.
            for k in 0..dry.precipitation.len() {
                assert!(
                    dry.wetness[k] >= dry.lakes[k] - 1e-12,
                    "{k}: the year's wetness without its lakes"
                );
                let (fell, gave, into) = (dry.precipitation[k], dry.evaporation[k], dry.spill[k]);
                assert!(
                    (gave - fell - into).abs() <= 0.05 * fell + 1e-12,
                    "{years} years, {k}: gave {gave} of {fell} and {into}"
                );
            }
        }
        let wet = strip(&air, |_| (0.3, 0.0), |_| 0.0, false);
        assert!(wet.lakes.iter().all(|&f| f == 0.0));
    }

    /// M74: the wind carries a place's vapour to the place downwind of it, and a land sheds
    /// the share of its rain that falls on ground already full, in the year's own loop. A
    /// strip of latitudes, its sea open from 10 to 20 degrees north alone, ice south of it and
    /// land north, under a wind blowing north along its isobars from 10 to 50 degrees, with
    /// no transport, neither eddies nor overturning: the land
    /// downwind of the sea takes its rain, tens of times what falls upwind of it, where the
    /// air's first vapour alone rains out; and that land's buckets, never full, shed the
    /// square of their fill of what falls on them, where a bucket that shed only what
    /// overflowed it would shed nothing. And what that land would give wet is Penman's (M75),
    /// by hand from the year's own air, ground and wetness.
    #[test]
    fn a_wind_carries_its_vapour_downwind_to_a_land_that_sheds() {
        let air = Air {
            transport: 0.0,
            ..earth_air()
        };
        let band = |lat: f64, lo: f64, hi: f64| (lo..hi).contains(&lat.to_degrees());
        let year = strip(
            &air,
            |lat| {
                if band(lat, 10.0, 20.0) {
                    (1.0, 0.0)
                } else if lat.to_degrees() < 10.0 {
                    (0.0, 1.0)
                } else {
                    (0.0, 0.0)
                }
            },
            |lat| if band(lat, 10.0, 50.0) { 0.3 } else { 0.0 },
            false,
        );
        let n = year.precipitation.len();
        let dphi = PI / (n - 1) as f64;
        let lat = |i: usize| (-PI / 2.0 + dphi * i as f64).to_degrees();
        let day = SECONDS_PER_DAY;
        // The rain over a band of latitudes, kg/m²/s.
        let rain = |lo: f64, hi: f64| {
            let (mut sum, mut area) = (0.0, 0.0);
            for i in (0..n).filter(|&i| lat(i) >= lo && lat(i) < hi) {
                sum += year.weights[i] * year.precipitation[i];
                area += year.weights[i];
            }
            sum / area
        };
        let (upwind, downwind) = (rain(-20.0, 6.0), rain(22.0, 42.0));
        println!(
            "upwind of the sea {:.3} mm a day, downwind of it {:.3}",
            upwind * day,
            downwind * day
        );
        assert!(downwind * day > 0.2, "{downwind}");
        assert!(downwind > 20.0 * upwind, "{downwind} {upwind}");
        // What the wind carries north across the sea's edge is its upwind place's: all that
        // lies beyond rains what its ground gives it and what the wind brings, the sea's own
        // column at the edge's wind, and not the land's beyond it.
        let i = (0..n).rfind(|&i| lat(i) < 20.0).unwrap();
        let gained: f64 = (i + 1..n)
            .map(|j| year.weights[j] * (year.precipitation[j] - year.evaporation[j]))
            .sum();
        let edge = (-PI / 2.0 + dphi * (i as f64 + 0.5)).cos();
        let carried = |from: usize| 0.3 * edge * year.column[from] / (2.0 * air.radius);
        println!(
            "across the sea's edge the land beyond gains {:.4e}, by hand from upwind {:.4e}, from downwind {:.4e}",
            gained,
            carried(i),
            carried(i + 1)
        );
        assert!((gained / carried(i) - 1.0).abs() < 0.01, "{gained}");
        assert!((gained / carried(i + 1) - 1.0).abs() > 0.1, "{gained}");
        // The land downwind: its buckets under three quarters full, shedding the square of
        // their fill, the fill three quarters of the wetness the year records.
        let mut shed = 0;
        for i in (0..n).filter(|&i| lat(i) >= 22.0 && lat(i) < 42.0) {
            assert_eq!(year.land[i], 1.0);
            let (wet, fell, ran) = (year.wetness[i], year.precipitation[i], year.runoff[i]);
            assert!(wet > 0.0 && wet < 1.0, "{i}: wetness {wet}");
            assert!(ran > 0.0 && ran < fell, "{i}: {ran} of {fell}");
            let by_hand = bucket_sheds(BUCKET_FULL_SHARE * wet);
            println!(
                "{:.0} degrees: wetness {wet:.3}, sheds {:.4} of its rain, by hand {by_hand:.4}",
                lat(i),
                ran / fell
            );
            assert!(
                (ran / fell / by_hand - 1.0).abs() < 0.05,
                "{i}: {ran} {fell} {by_hand}"
            );
            // What it would give wet: the heat its ground, a kelvin over its air, has for
            // its air, and its air's drying power, at the wind's and the gusts' speed.
            let ta = 300.0 - 45.0 * lat(i).to_radians().sin().powi(2);
            let (soil, speed) = (ta + 1.0, 0.3f64.hypot(6.0));
            let rho = air.pressure * air.molar_mass / (GAS_CONSTANT * ta);
            let layer = crate::surface_layer::SurfaceLayer {
                gravity: air.gravity,
                heat_capacity: air.heat_capacity,
            };
            let x = rho * AIR_EXCHANGE_COEFFICIENT * speed * layer.exchange(soil, ta, speed);
            let q = |t: f64| air.saturated_humidity(t, air.pressure);
            let r = year.saturation[i];
            let sensible = air.heat_capacity / air.latent(ta) * x * (soil - ta);
            let latent = x * wet * (q(soil) - r * q(ta));
            let potential = air.penman((ta, air.pressure), sensible, latent, x * (1.0 - r) * q(ta));
            assert!(
                (year.potential[i] / potential - 1.0).abs() < 1e-3,
                "{i}: {} {potential}",
                year.potential[i]
            );
            // Under the draw of a wet surface as warm as its ground, which gives the more.
            assert!(year.potential[i] < x * (q(soil) - r * q(ta)));
            shed += 1;
        }
        assert!(shed >= 4, "{shed}");
    }

    /// M74, M75: the year's vapour over an aquaplanet. Its vapour is kept: what rains is what
    /// its sea gave. Its overturning gathers it at the equator, which rains more than its sea
    /// gives, and its subtropics less: across a latitude within the cell, by hand, the mass
    /// the overturning's heat turns over, its heat over its gross stability, times the colder
    /// latitude's surface air's vapour. The trades' drift across the isobars carries no
    /// column of its own: with the drift or without, the same vapour crosses. Its eddies
    /// carry it out of the cell's reach as the latent share of their heat, some half of it at
    /// 40 degrees; were all of the transport theirs within the cell too, the equator would
    /// lose where it gains. Its sea gives into its surface air, CLOUD_SURFACE_HUMIDITY of
    /// saturation, whatever its column holds.
    #[test]
    fn a_year_s_vapour_goes_by_its_overturning_and_its_eddies() {
        let air = earth_air();
        let day = SECONDS_PER_DAY;
        let year = aquaplanet(&air, 3.0);
        let n = year.precipitation.len();
        let dphi = PI / (n - 1) as f64;
        let lat = |i: usize| (-PI / 2.0 + dphi * i as f64).to_degrees();
        let brief = year.summary();
        assert!(
            (brief.precipitation / brief.evaporation - 1.0).abs() < MOISTURE_SETTLED,
            "{brief:?}"
        );
        // What a band gains from beyond it, kg/m²/s over the band.
        let gain = |m: &Moisture, lo: f64, hi: f64| {
            let (mut sum, mut area) = (0.0, 0.0);
            for i in (0..n).filter(|&i| lat(i) >= lo && lat(i) <= hi) {
                sum += m.weights[i] * (m.precipitation[i] - m.evaporation[i]);
                area += m.weights[i];
            }
            sum / area
        };
        let (equator, subtropics) = (gain(&year, -9.0, 9.0), gain(&year, 14.0, 30.0));
        println!(
            "the equator gains {:.2} mm a day, the subtropics {:.2}; rain {:.2} and {:.2}",
            equator * day,
            subtropics * day,
            year.precipitation[n / 2] * day,
            year.precipitation[n / 2 + 5] * day
        );
        assert!(
            equator * day > 0.3 && subtropics * day < -0.3,
            "{equator} {subtropics}"
        );
        let south_of = |m: &Moisture, edge: f64| -> f64 {
            (0..n)
                .filter(|&i| lat(i) < edge)
                .map(|i| m.weights[i] * (m.evaporation[i] - m.precipitation[i]))
                .sum()
        };
        let links = crate::climate::place_links((n, 1));
        let t = |i: usize| 300.0 - 45.0 * lat(i).to_radians().sin().powi(2);
        // The surface air's vapour and its energy at a latitude, the year's own.
        let vapour =
            |m: &Moisture, i: usize| m.saturation[i] * air.saturated_humidity(t(i), air.pressure);
        // The overturning alone, its cell reaching far past the poles and no wind: across 18
        // degrees north the vapour that crosses toward the equator, by hand.
        let cell = Air { reach: 1e3, ..air };
        let pure = aquaplanet(&cell, 0.0);
        let i = (0..n).rfind(|&i| lat(i) < 18.0).unwrap();
        let g = links[i].iter().find(|l| l.0 == i + 1).unwrap().1;
        let crossed = -south_of(&pure, lat(i) + 1.0);
        let (mut held, mut fell) = (0.0, 0.0);
        for j in 0..n {
            let r = pure.weights[j] * pure.precipitation[j];
            held += r * cell.energy(t(j), vapour(&pure, j));
            fell += r;
        }
        let rising = held / fell;
        let low = cell.energy(t(i + 1), vapour(&pure, i + 1));
        let heat = cell.transport * (1.0 - cell.water) * g * (t(i) - t(i + 1)) * pure.weights[i];
        let by_hand = heat * cell.cell_share(0.0, 0.0, 0.0) / cell.cell_stability(rising, low)
            * vapour(&pure, i + 1);
        println!(
            "the overturning alone: across 18 degrees {:.4e} toward the equator, by hand {:.4e}; its rising air {:.1} kJ/kg, its stability there {:.1}",
            crossed,
            by_hand,
            rising / 1e3,
            cell.cell_stability(rising, low) / 1e3
        );
        assert!(
            crossed > 0.0 && (crossed / by_hand - 1.0).abs() < 0.005,
            "{crossed} {by_hand}"
        );
        // A transport a hundred times as strong would move many times a place's column in a
        // step: carried in shorter steps, its vapour is kept, every column holds some, and
        // its equator gains the more.
        let strong = aquaplanet(
            &Air {
                transport: 100.0 * air.transport,
                ..cell
            },
            0.0,
        );
        let b = strong.summary();
        assert!(
            (b.precipitation / b.evaporation - 1.0).abs() < MOISTURE_SETTLED,
            "{b:?}"
        );
        assert!(strong.column.iter().all(|&w| w.is_finite() && w > 0.0));
        let more = gain(&strong, -9.0, 9.0) / gain(&pure, -9.0, 9.0);
        println!("a hundred times the transport: the equator gains {more:.2} times as much");
        assert!((1.2..100.0).contains(&more), "{more}");
        // The drift across the isobars carries no column: the trades change what the sea
        // gives under them, by their speed, and nothing of what crosses beyond them.
        let calm = aquaplanet(&air, 0.0);
        let beyond = (0..n).rfind(|&i| lat(i) < 40.0).unwrap();
        for m in [&year, &calm] {
            assert!(gain(m, -9.0, 9.0) > 0.0);
        }
        let ratio = gain(&year, -9.0, 9.0) / gain(&calm, -9.0, 9.0);
        println!("with its trades the equator gains {ratio:.3} of what it gains without");
        assert!((0.95..1.08).contains(&ratio), "{ratio}");
        // The eddies' vapour across 40 degrees against the heat they carry: its latent share.
        let i = beyond;
        let g = links[i].iter().find(|l| l.0 == i + 1).unwrap().1;
        let heat = air.transport * (1.0 - air.water) * g * (t(i) - t(i + 1)) * year.weights[i];
        let latent = air.vapour.l_vap * south_of(&year, 40.0);
        // By hand: the air the heat moves, over what a kilogram of it carries a kelvin,
        // carrying the difference of the two latitudes' surface air's vapour.
        let carries = |i: usize| {
            air.heat_capacity + air.latent_capacity(t(i), air.pressure, year.saturation[i])
        };
        let by_hand = air.vapour.l_vap * heat / (t(i) - t(i + 1))
            * (vapour(&year, i) - vapour(&year, i + 1))
            / (0.5 * (carries(i) + carries(i + 1)));
        println!(
            "across 40 degrees the vapour carries {:.2} of the eddies' heat, by hand {:.2}",
            latent / heat,
            by_hand / heat
        );
        assert!((0.15..0.6).contains(&(latent / heat)), "{latent} {heat}");
        assert!((latent / by_hand - 1.0).abs() < 0.1, "{latent} {by_hand}");
        // All of the transport the eddies' within the cell too: no overturning.
        let mixed = aquaplanet(&Air { reach: 0.0, ..air }, 3.0);
        let kept = gain(&mixed, -9.0, 9.0) / equator;
        println!("with the eddies' whole share within the cell the equator gains {kept:.2} of it");
        assert!(kept < 0.0, "{kept}");
        // The sea gives into its surface air whatever its column holds.
        let p = n / 2 + 9;
        let (ta, ts) = (t(p), t(p) + 1.0);
        let rho = air.pressure * air.molar_mass / (GAS_CONSTANT * ta);
        let layer = crate::surface_layer::SurfaceLayer {
            gravity: air.gravity,
            heat_capacity: air.heat_capacity,
        };
        let by_hand = rho
            * AIR_EXCHANGE_COEFFICIENT
            * 6.0
            * layer.exchange(ts, ta, 6.0)
            * (air.saturated_humidity(ts, air.pressure)
                - CLOUD_SURFACE_HUMIDITY * air.saturated_humidity(ta, air.pressure));
        for m in [&year, &pure, &mixed] {
            assert!(
                (m.evaporation[p] / by_hand - 1.0).abs() < 1e-9,
                "{} {by_hand}",
                m.evaporation[p]
            );
        }
        let e = n / 2;
        assert!((pure.saturation[e] - mixed.saturation[e]).abs() > 0.01);
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
