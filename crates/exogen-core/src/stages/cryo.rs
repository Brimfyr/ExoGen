// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! What an ice shell does with its heat and its ocean's salts (PLAN.md 4.3 and Phase
//! 4): chaos where a thin shell melts through, cryovolcanic domes and salt-bright crater
//! floors on a shell over rock, and a polar plume's fallout of fresh ice.
//!
//! Chaos is a population of patches, their number from the share of the surface the
//! shell's heat disrupts and their sizes from a distribution, drawn in a heat field; the
//! matrix inside a patch drops and breaks into hummocks around rafts of the old surface,
//! and the ocean's salts stain it. On a shell over rock the ice magma that does not
//! flood extrudes as domes of brine-rich ice, which relax like any other load, and a
//! crater deep enough to reach the brines gets salts on its floor that darken with age.
//! A plume lays fresh ice down over the body, fastest near its pole: it buries lags and
//! salts before they can show and keeps the surface in frost.

use rayon::prelude::*;

use crate::constants::*;
use crate::crater::{CellWrite, CraterRecord, StampOutput, Stamper, dimensions};
use crate::fields::LayerSet;
use crate::geom::{angle_between, cross, dot, local_basis, normalize, rotate};
use crate::grid::CubeSphere;
use crate::grid::projection::{PaddedField, sample_bilinear};
use crate::noise::GradientNoise;
use crate::params::ParamGraph;
use crate::rng::Rng;
use crate::units::*;

pub const VERSION: u64 = 3;

/// A polar plume's fallout of fresh ice.
pub struct Fallout {
    pole: [f64; 3],
    /// The rate at the pole, m/Gyr.
    rate: f64,
    angle: f64,
}

/// The share of the body a plume's fallout keeps under frost (F47, F52): the frost cover
/// at each distance from the pole, integrated over the sphere.
pub fn mean_frost_share(rate_pole: f64) -> f64 {
    if rate_pole <= 0.0 {
        return 0.0;
    }
    let angle = PLUME_FALLOUT_ANGLE_DEG.to_radians();
    let n = 360;
    let mut sum = 0.0;
    for k in 0..n {
        let theta = (k as f64 + 0.5) / n as f64 * std::f64::consts::PI;
        let rate =
            rate_pole * (PLUME_FALLOUT_FAR + (1.0 - PLUME_FALLOUT_FAR) * (-theta / angle).exp());
        sum += (1.0 - (-rate / PLUME_FROST_RATE).exp()) * theta.sin();
    }
    0.5 * sum * std::f64::consts::PI / n as f64
}

impl Fallout {
    pub fn new(graph: &ParamGraph) -> Option<Self> {
        let rate = graph.f64("volcanism.plume_fallout_rate");
        if rate <= 0.0 {
            return None;
        }
        let pole = crate::stages::ice::ShellMap::new(graph).pole?;
        Some(Self {
            pole,
            rate,
            angle: PLUME_FALLOUT_ANGLE_DEG.to_radians(),
        })
    }

    /// The fallout rate at a direction, m/Gyr.
    pub fn rate_at(&self, d: [f64; 3]) -> f64 {
        let x = angle_between(self.pole, d) / self.angle;
        self.rate * (PLUME_FALLOUT_FAR + (1.0 - PLUME_FALLOUT_FAR) * (-x).exp())
    }

    /// How long the fallout takes to bury a lag or a salt at a direction, Gyr.
    pub fn burial_gyr(&self, d: [f64; 3]) -> f64 {
        PLUME_BURIAL_DEPTH / self.rate_at(d).max(1e-30)
    }

    /// The frost cover the fallout keeps at a direction, 0 to 1.
    pub fn frost_at(&self, d: [f64; 3]) -> f64 {
        // A soft saturation: the fallout wins against the lag and the gardening by
        // degrees, so the far side keeps some of its own contrast under a thinner frost.
        1.0 - (-self.rate_at(d) / PLUME_FROST_RATE).exp()
    }

    /// How much a fissure at a direction is one of the plume's warm vents, 0 to 1.
    pub fn vent_weight(&self, d: [f64; 3]) -> f64 {
        (-angle_between(self.pole, d) / self.angle).exp()
    }
}

/// The plumes' fallout mantles the ground it lies on (F74). Snow laid over billions of
/// years creeps down the slopes of what it covers, filling the hollows and thinning on the
/// crests, so relief smaller than its thickness is buried and what stands larger is
/// rounded by no more than it: Enceladus's old craters are muted (Kirchoff and Schenk
/// 2009). Its loose top, what drains into fractures, is a few hundred metres (Martin et al.
/// 2023); the whole deposit, sintered below that, mantles. The ground is drawn toward its
/// own low-pass over MANTLE_SPREAD times the deposit's thickness `cover`, m, by no more
/// than the thickness either way; the deposit's mean, which only raises the ground, is
/// left out. The mantle buries only what is older than it: the finish lays it before the
/// unresolved craters, which count only those dug since the cover above them grew
/// thinner than they are deep (LooseCover::buried_before).
pub fn mantle(layers: &mut LayerSet, body_radius: f64, cover: &[f32]) {
    let grid = layers.grid;
    let cell = grid.delta() * body_radius;
    let thickest = cover.iter().cloned().fold(0.0f32, f32::max) as f64;
    if thickest <= 0.0 {
        return;
    }
    if std::env::var_os("EXOGEN_DEBUG").is_some() {
        let mut sorted: Vec<f32> = cover.iter().copied().filter(|h| *h > 0.0).collect();
        sorted.sort_by(f32::total_cmp);
        let at = |q: f64| sorted[((sorted.len() - 1) as f64 * q) as usize];
        eprintln!(
            "mantle: fallout {:.0} m at the 10th percentile, {:.0} at the median, {:.0} at the 90th, {:.0} at most",
            at(0.1),
            at(0.5),
            at(0.9),
            thickest
        );
    }
    let lows =
        crate::stages::relaxation::low_passes(&layers.bedrock, cell, MANTLE_SPREAD * thickest);
    if lows.is_empty() {
        return;
    }
    layers
        .bedrock
        .as_mut_slice()
        .par_iter_mut()
        .enumerate()
        .for_each(|(k, z)| {
            let h = cover[k] as f64;
            if h <= 0.0 {
                return;
            }
            let here = *z as f64;
            let want = MANTLE_SPREAD * h;
            // The low-pass at that length, between the pyramid's levels about it.
            let low = match lows.iter().position(|(len, _)| *len >= want) {
                Some(0) => here + (want / lows[0].0) * (lows[0].1[k] as f64 - here),
                Some(i) => {
                    let (l0, l1) = (lows[i - 1].0, lows[i].0);
                    let t = (want / l0).ln() / (l1 / l0).ln();
                    lows[i - 1].1[k] as f64 * (1.0 - t) + lows[i].1[k] as f64 * t
                }
                None => lows[lows.len() - 1].1[k] as f64,
            };
            *z = (here + (low - here).clamp(-h, h)) as f32;
        });
}

/// The loose cover the plumes' fallout has laid on the ground (F64): the snow of the
/// present plume and of each past episode's, each falling off from its own pole as the
/// fallout does (F47), summed over the time a surface has stood.
pub struct LooseCover {
    /// Each plume's pole in the present frame, and when it ran, Gyr.
    plumes: Vec<([f64; 3], f64, f64)>,
    /// The rate at a plume's pole, m/Gyr.
    rate: f64,
    angle: f64,
    /// The depth a crack stays open to, m: what it can drain.
    crack_depth: f64,
}

impl LooseCover {
    /// The plumes that have run by `present`: one at the centre of each episode of
    /// focusing (F53) while it lasted, each venting as the present one does.
    pub fn new(
        graph: &ParamGraph,
        episodes: &[crate::stages::ice::FocusEpisode],
        present: f64,
    ) -> Option<Self> {
        let rate = graph.f64("volcanism.plume_fallout_rate");
        if rate <= 0.0 {
            return None;
        }
        let shell = crate::stages::ice::ShellMap::new(graph);
        let pole = shell.pole?;
        let mut plumes: Vec<([f64; 3], f64, f64)> = episodes
            .iter()
            .filter(|e| e.start < present)
            .map(|e| (e.center, e.start, e.end.min(present)))
            .collect();
        let start = (graph.f64("body.age") - FOCUSING_EPISODE_DURATION_GYR).max(0.0);
        if shell.has_active_terrain() && present > start {
            plumes.push((pole, start, present));
        }
        // A crack stays open down to where the weight of the ice matches the tension
        // that opened it.
        let crack_depth =
            ICE_CRACK_STRESS / (DENSITY_CRUST_ICE * graph.f64("body.gravity").max(1e-6));
        Some(Self::of_plumes(plumes, rate, crack_depth))
    }

    /// The cover of given plumes: each one's pole and the times it ran between, Gyr, all
    /// at one rate at the pole, m/Gyr; and the depth a crack can drain, m.
    pub fn of_plumes(plumes: Vec<([f64; 3], f64, f64)>, rate: f64, crack_depth: f64) -> Self {
        Self {
            plumes,
            rate,
            angle: PLUME_FALLOUT_ANGLE_DEG.to_radians(),
            crack_depth,
        }
    }

    /// How much of a cover a crack drains, m: no more than the crack is deep.
    pub fn drained(&self, cover: f64) -> f64 {
        cover.min(self.crack_depth)
    }

    /// The time before which what was dug at a direction lies buried by `until`, Gyr:
    /// after it no more than `depth` of cover fell. `since` when even all the cover
    /// since then is thinner (F74).
    pub fn buried_before(&self, d: [f64; 3], since: f64, until: f64, depth: f64) -> f64 {
        if self.thickness_between(d, since, until) <= depth {
            return since;
        }
        let (mut lo, mut hi) = (since, until);
        for _ in 0..40 {
            let mid = 0.5 * (lo + hi);
            if self.thickness_between(d, mid, until) > depth {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    }

    /// The cover laid at a direction between two times, m: what lies on a surface that
    /// formed at `since` by the time `until`.
    pub fn thickness_between(&self, d: [f64; 3], since: f64, until: f64) -> f64 {
        self.plumes
            .iter()
            .map(|&(pole, t0, t1)| {
                let span = t1.min(until) - t0.max(since);
                if span <= 0.0 {
                    return 0.0;
                }
                let x = angle_between(pole, d) / self.angle;
                self.rate * (PLUME_FALLOUT_FAR + (1.0 - PLUME_FALLOUT_FAR) * (-x).exp()) * span
            })
            .sum()
    }
}

/// A generation of chaos (F80): at one time the shell disrupts wherever the warmth the
/// ocean brings up under it stands above a threshold. The brittle lid over the warm ice
/// breaks into rafts, grouped into plates, along leads that the matrix fills; a plate
/// founders into the matrix with a chance that grows with the warmth's excess, leaving
/// some of its rafts standing as blocks, and what stands drifts and turns, carrying the
/// old surface on it.
#[derive(Clone, Debug)]
pub struct ChaosGeneration {
    pub time_gyr: f64,
    pub id: u64,
    warmth: Warmth,
    /// The warmth above which the shell disrupts, in spreads of its structured part.
    threshold: f64,
    rafts: Lattice,
    plates: Lattice,
    radius: f64,
}

/// The warmth the ocean brings up under the shell at one time (F80): structured at every
/// scale from the largest chaos down to the rafts, which no finer structure can show,
/// and strongest at low latitudes.
#[derive(Clone, Debug)]
struct Warmth {
    noise: GradientNoise,
    frequency: f64,
    octaves: u32,
    /// The structured part's spread, which the field is measured in.
    spread: f64,
}

/// Points on which a generation's warmth is measured to find its threshold.
const CHAOS_SAMPLES: usize = 1 << 16;

/// The `i`-th of `n` points spread evenly over the sphere (a Fibonacci lattice).
fn even_point(i: usize, n: usize) -> [f64; 3] {
    let z = 1.0 - (2.0 * i as f64 + 1.0) / n as f64;
    let r = (1.0 - z * z).max(0.0).sqrt();
    let phi = i as f64 * std::f64::consts::PI * (3.0 - 5f64.sqrt());
    [r * phi.cos(), r * phi.sin(), z]
}

/// The ocean's heat at a latitude, in spreads of the structured warmth: the square of the
/// cosine of the latitude less its mean over the sphere, two thirds, so the equator
/// stands CHAOS_EQUATOR above the poles.
fn latitude_warmth(d: [f64; 3]) -> f64 {
    CHAOS_EQUATOR * (1.0 - d[2] * d[2] - 2.0 / 3.0)
}

impl Warmth {
    fn structured(&self, d: [f64; 3]) -> f64 {
        self.noise.fbm(
            d,
            self.frequency,
            self.octaves,
            2.0,
            CHAOS_WARMTH_PERSISTENCE,
        )
    }

    fn at(&self, d: [f64; 3]) -> f64 {
        self.structured(d) / self.spread + latitude_warmth(d)
    }
}

/// Jittered points, one in each cell of a cubic lattice in space: each place on the
/// sphere belongs to the nearest, so the points break the sphere into polygons.
#[derive(Clone, Debug)]
struct Lattice {
    seed: u64,
    /// The spacing as a share of the body's radius, and in metres.
    step: f64,
    spacing: f64,
}

/// A lattice point and its key.
type Site = ([f64; 3], u64);

impl Lattice {
    fn new(seed: u64, spacing: f64, body_radius: f64) -> Self {
        Self {
            seed,
            step: spacing / body_radius,
            spacing,
        }
    }

    /// The points of the 27 cells round a place's, and which of them is nearest.
    fn around(&self, p: [f64; 3]) -> ([Site; 27], usize) {
        let h = self.step;
        let base = [
            (p[0] / h).floor() as i64,
            (p[1] / h).floor() as i64,
            (p[2] / h).floor() as i64,
        ];
        let rng = Rng::new(self.seed);
        let mut sites = [([0.0f64; 3], 0u64); 27];
        let mut nearest = (f64::INFINITY, 0usize);
        let mut n = 0;
        for di in -1..=1 {
            for dj in -1..=1 {
                for dk in -1..=1 {
                    let c = [base[0] + di, base[1] + dj, base[2] + dk];
                    let pack = |v: i64| ((v + (1 << 20)) as u64) & 0x1F_FFFF;
                    let key = (pack(c[0]) << 42) | (pack(c[1]) << 21) | pack(c[2]);
                    let bits = rng.bits(key, 0);
                    let jitter = |s: u32| ((bits >> s) & 0x1F_FFFF) as f64 / (1u64 << 21) as f64;
                    let site = [
                        (c[0] as f64 + jitter(0)) * h,
                        (c[1] as f64 + jitter(21)) * h,
                        (c[2] as f64 + jitter(42)) * h,
                    ];
                    let r = [p[0] - site[0], p[1] - site[1], p[2] - site[2]];
                    if dot(r, r) < nearest.0 {
                        nearest = (dot(r, r), n);
                    }
                    sites[n] = (site, key);
                    n += 1;
                }
            }
        }
        (sites, nearest.1)
    }

    fn nearest(&self, p: [f64; 3]) -> Site {
        let (sites, i) = self.around(p);
        sites[i]
    }
}

/// How a raft stands in a generation: on its plate, which stands; on its own, its plate
/// foundered; or foundered with it.
#[derive(Clone, Copy, PartialEq)]
enum Standing {
    Plate,
    Block,
    Foundered,
}

/// The stand of a raft whose disruption is known: its plate founders with a chance of
/// the plate's disruption, and then the raft with a chance of its own.
struct RaftState {
    plate: u64,
    standing: Standing,
}

impl ChaosGeneration {
    fn new(g: u64, time_gyr: f64, seed: u64, body_radius: f64, share: f64) -> ChaosGeneration {
        let spacing = CHAOS_RAFT_LIDS * CHAOS_LENS_DEPTH;
        let rng = Rng::new(seed);
        // Noise frequencies are lattice cells per radian: the radius over the length.
        let octaves = ((CHAOS_LARGEST / spacing).log2().ceil() as u32).max(1);
        let mut warmth = Warmth {
            noise: GradientNoise::new(rng.bits(g, 1)),
            frequency: body_radius / CHAOS_LARGEST,
            octaves,
            spread: 1.0,
        };
        let structured: Vec<f64> = (0..CHAOS_SAMPLES)
            .into_par_iter()
            .map(|i| warmth.structured(even_point(i, CHAOS_SAMPLES)))
            .collect();
        let mean = structured.iter().sum::<f64>() / CHAOS_SAMPLES as f64;
        let variance =
            structured.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / CHAOS_SAMPLES as f64;
        warmth.spread = variance.sqrt().max(1e-12);
        let mut values: Vec<f64> = structured
            .iter()
            .enumerate()
            .map(|(i, v)| v / warmth.spread + latitude_warmth(even_point(i, CHAOS_SAMPLES)))
            .collect();
        values.sort_by(f64::total_cmp);
        let at = (((1.0 - share) * CHAOS_SAMPLES as f64) as usize).min(CHAOS_SAMPLES - 1);
        ChaosGeneration {
            time_gyr,
            id: (0xC4A0 << 32) | g,
            warmth,
            threshold: values[at],
            rafts: Lattice::new(rng.bits(g, 2), spacing, body_radius),
            plates: Lattice::new(rng.bits(g, 3), CHAOS_PLATE_RAFTS * spacing, body_radius),
            radius: body_radius,
        }
    }

    /// How far the shell under a point is disrupted, 0 to 1, or none where its warmth is
    /// below the threshold and it stays whole.
    fn disruption(&self, d: [f64; 3]) -> Option<f64> {
        let excess = self.warmth.at(d) - self.threshold;
        (excess > 0.0).then(|| (excess / CHAOS_DISRUPTION_SPREADS).min(1.0))
    }

    /// How a raft, given by its lattice point, stands, with its own disruption.
    fn raft_state(&self, site: Site, disruption: f64) -> RaftState {
        let (plate, plate_key) = self.plates.nearest(site.0);
        let plate_disruption = self.disruption(normalize(plate)).unwrap_or(0.0);
        let standing = if Rng::new(self.plates.seed).uniform(plate_key, 1) >= plate_disruption {
            Standing::Plate
        } else if Rng::new(self.rafts.seed).uniform(site.1, 1) >= disruption {
            Standing::Block
        } else {
            Standing::Foundered
        };
        RaftState {
            plate: plate_key,
            standing,
        }
    }

    /// The share of the body's surface whose rafts the generation foundered into its
    /// matrix, measured on points spread evenly over it (its leads aside).
    pub fn matrix_share(&self) -> f64 {
        let n = 1 << 14;
        (0..n)
            .into_par_iter()
            .filter(|&i| {
                let (sites, ia) = self.rafts.around(even_point(i, n));
                let a = sites[ia];
                self.disruption(normalize(a.0))
                    .is_some_and(|e| self.raft_state(a, e).standing == Standing::Foundered)
            })
            .count() as f64
            / n as f64
    }

    /// Whether the generation disrupted the shell at a point: whether the raft it lies on
    /// stands over warmth above the threshold.
    pub fn disrupts(&self, d: [f64; 3]) -> bool {
        let site = self.rafts.nearest(d);
        self.disruption(normalize(site.0)).is_some()
    }
}

/// The chaos present at the epoch (F80): the generations within the retention age, one
/// every CHAOS_GENERATIONS-th of it through the body's life, so that the retention age
/// always holds as many, each with a warmth field of its own and a share of the surface
/// that leaves the chaos share disrupted by them all together.
pub fn chaos_population(graph: &ParamGraph, present: f64, seed: u64) -> Vec<ChaosGeneration> {
    let fraction = graph.f64("tectonics.chaos_fraction").clamp(0.0, 0.99);
    if fraction <= 0.0 {
        return Vec::new();
    }
    let body_radius = graph.f64("body.radius");
    let retention = graph.f64("cratering.retention_age").max(1e-6);
    let n = CHAOS_GENERATIONS.max(1) as f64;
    let slot = retention / n;
    let share = 1.0 - (1.0 - fraction).powf(1.0 / n);
    let seed = seed ^ 0xC4A0_5000;
    let first = ((present - retention) / slot).floor().max(0.0) as u64;
    let last = (present / slot).floor().max(0.0) as u64;
    (first..=last)
        .filter_map(|g| {
            let time = (g as f64 + 0.5) * slot;
            (time <= present && time > present - retention)
                .then(|| ChaosGeneration::new(g, time, seed, body_radius, share))
        })
        .collect()
}

fn write(idx: usize, bedrock: f64) -> CellWrite {
    CellWrite {
        idx: idx as u32,
        bedrock: bedrock as f32,
        surface: false,
        age: 0.0,
        unit: 0,
        deep: 0.0,
        melt: 0.0,
        cover: 0.0,
    }
}

/// The stain of a chaos matrix as it stands at an age, Gyr: the ocean's material stains
/// it as it stains a crack's margins (F78), from the start, and fades with the salts to
/// the ancient band's share.
pub fn matrix_stain(age: f64) -> f64 {
    LINEA_STAIN
        * (LINEA_STAIN_RESIDUAL + (1.0 - LINEA_STAIN_RESIDUAL) * (-age / SALT_FADE_GYR).exp())
}

/// A generation of chaos on the resolved grid (F80). Each cell lies on a raft; where the
/// raft's warmth is above the threshold, the matrix fills the share of the cell that a
/// lead covers, or the whole cell where the raft foundered, and the rest shows the old
/// surface from where the raft's plate, or the raft on its own, carried it. Leads open
/// along the edges between plates, round a block, and round the chaos, where the lid
/// breaks from the shell that stays whole; a cell of a plate takes the plate's width of
/// lead and a block's cell the block's. The matrix lies at the old surface's level
/// smoothed over a few rafts, below it by the drop, stained by the ocean's material; what
/// stands keeps the age, the unit and the stain of what it carries.
pub fn chaos_compute(
    layers: &LayerSet,
    grid: &CubeSphere,
    cell: f64,
    present: f64,
    generation: &ChaosGeneration,
) -> StampOutput {
    let bedrock = layers.bedrock.as_slice();
    let ages = layers.surface_age.as_slice();
    let units = layers.unit_id.as_slice();
    let deep = layers.deep.as_slice();
    let melt = layers.melt.as_slice();
    let cover = layers.cover.as_slice();
    let padded_bedrock = PaddedField::new(grid, bedrock);
    let padded_stain = PaddedField::new(grid, layers.stain.as_slice());
    let level = crate::stages::relaxation::low_pass(
        &layers.bedrock,
        cell,
        CHAOS_MATRIX_LEVEL_RAFTS * generation.rafts.spacing,
    );
    let radius = generation.radius;
    // Hummocks finer than the grid average out over a cell, to the share of their
    // relief that one cell's mean keeps.
    let hummocks = GradientNoise::new(generation.id ^ 0x4B0C_0000);
    let hummock_frequency = radius / CHAOS_HUMMOCK_SPACING;
    let hummock = CHAOS_HUMMOCK * (CHAOS_HUMMOCK_SPACING / cell).min(1.0);
    let time = generation.time_gyr as f32;
    let stain_matrix = matrix_stain((present - generation.time_gyr).max(0.0)) as f32;
    let turn_max = CHAOS_RAFT_TURN_DEG.to_radians();
    let writes: Vec<(CellWrite, [f32; 2])> = (0..grid.cell_count())
        .into_par_iter()
        .filter_map(|k| {
            let d = grid.center(grid.cell(k));
            let (sites, ia) = generation.rafts.around(d);
            let a = sites[ia];
            let disruption = generation.disruption(normalize(a.0))?;
            let state = generation.raft_state(a, disruption);
            // What carries the cell's surface: the plate, or the block on its own.
            let (centre, key, lattice, carried) = match state.standing {
                Standing::Plate => {
                    let (p, key) = generation.plates.nearest(a.0);
                    let centre = normalize(p);
                    let e = generation.disruption(centre).unwrap_or(0.0);
                    (centre, key, &generation.plates, e)
                }
                _ => (normalize(a.0), a.1, &generation.rafts, disruption),
            };
            let half_lead = 0.5 * CHAOS_LEAD_SHARE * lattice.spacing;
            // The nearest open edge: the planes halfway to the other rafts' points, the
            // edge open between plates, round a block, and round the chaos.
            let r = [d[0] - a.0[0], d[1] - a.0[1], d[2] - a.0[2]];
            let near = dot(r, r);
            let mut edge = f64::INFINITY;
            for (i, &b) in sites.iter().enumerate() {
                if i == ia {
                    continue;
                }
                let ab = [b.0[0] - a.0[0], b.0[1] - a.0[1], b.0[2] - a.0[2]];
                let len = dot(ab, ab).sqrt();
                if len <= 0.0 {
                    continue;
                }
                let rb = [d[0] - b.0[0], d[1] - b.0[1], d[2] - b.0[2]];
                let distance = (dot(rb, rb) - near) / (2.0 * len) * radius;
                if distance >= edge || distance >= half_lead + 0.5 * cell {
                    continue;
                }
                let open = state.standing != Standing::Plate
                    || match generation.disruption(normalize(b.0)) {
                        None => true,
                        Some(e) => {
                            let other = generation.raft_state(b, e);
                            other.plate != state.plate || other.standing != Standing::Plate
                        }
                    };
                if open {
                    edge = distance;
                }
            }
            // The share of the cell the lead covers: the overlap of the cell's width with
            // the lead's along the line to the edge.
            let lead = ((edge + 0.5 * cell).min(half_lead) - (edge - 0.5 * cell).max(-half_lead))
                .max(0.0)
                / cell;
            let matrix = if state.standing == Standing::Foundered {
                1.0
            } else {
                lead.min(1.0)
            };
            let matrix_bedrock = level[k] as f64 - CHAOS_DROP
                + hummock * hummocks.fbm(d, hummock_frequency, 2, 2.0, 0.5);
            // What stands carries the cell's surface from where its drift and its turn
            // about its centre brought it.
            let rng = Rng::new(lattice.seed);
            let (east, north) = local_basis(centre);
            let heading = 2.0 * std::f64::consts::PI * rng.uniform(key, 2);
            let (s, c) = heading.sin_cos();
            let along = [
                east[0] * c + north[0] * s,
                east[1] * c + north[1] * s,
                east[2] * c + north[2] * s,
            ];
            let drift = CHAOS_RAFT_DRIFT * lattice.spacing * carried * rng.uniform(key, 3);
            let turn = turn_max * carried * (2.0 * rng.uniform(key, 4) - 1.0);
            let back = rotate(d, normalize(cross(centre, along)), -drift / radius);
            let source = rotate(back, centre, -turn);
            let raft_bedrock = sample_bilinear(grid, &padded_bedrock, source) as f64;
            let raft_stain = sample_bilinear(grid, &padded_stain, source);
            let from = grid.index(grid.locate(source).0);
            let mut w = write(k, matrix * matrix_bedrock + (1.0 - matrix) * raft_bedrock);
            w.surface = true;
            if matrix > 0.5 {
                w.age = time;
                w.unit = UNIT_CHAOS;
            } else {
                w.age = ages[from];
                w.unit = units[from];
                w.deep = deep[from];
                w.melt = melt[from];
                w.cover = cover[from];
            }
            let m = matrix as f32;
            let stain = m * stain_matrix + (1.0 - m) * raft_stain;
            Some((w, [stain, 1.0]))
        })
        .collect();
    let (cells, stain): (Vec<CellWrite>, Vec<[f32; 2]>) = writes.into_iter().unzip();
    StampOutput {
        cells,
        composition: None,
        stain: Some(stain),
        steep: None,
        ejecta_age: None,
        compositions: None,
    }
}

/// A cryovolcanic dome: brine-rich ice extruded onto a shell over rock.
#[derive(Clone, Debug)]
pub struct Cryodome {
    pub center: [f64; 3],
    pub radius: f64,
    /// Height as built, m; it relaxes with the shell.
    pub height: f64,
    pub time_gyr: f64,
    pub id: u64,
}

impl Cryodome {
    pub fn cap_radius(&self, body_radius: f64) -> f64 {
        1.3 * self.radius / body_radius
    }
}

/// The domes built through the body's life: per time bin, the dome share of the ice
/// magma budget spent in domes of volumes drawn log-uniform as a share of the body's
/// cube, placed at random.
pub fn cryodome_population(graph: &ParamGraph, present: f64, seed: u64) -> Vec<Cryodome> {
    if graph.str("interior.tectonic_regime") != "ice_shell_rock"
        || !crate::stages::ice::is_ice_surface(graph)
    {
        return Vec::new();
    }
    // The magma rate is kept in m³/s; the budget here is in km³ per Gyr.
    let rate = graph.f64("volcanism.rate") * GYR / 1e9 * CRYODOME_SHARE;
    if rate <= 0.0 {
        return Vec::new();
    }
    let body_radius = graph.f64("body.radius");
    let age = graph.f64("body.age");
    let r_km3 = (body_radius / 1e3).powi(3);
    let (lo, hi) = (
        CRYODOME_VOLUME_SHARE.0 * r_km3,
        CRYODOME_VOLUME_SHARE.1 * r_km3,
    );
    let rng = Rng::new(seed ^ 0xD0CE_5000);
    let mut carry = 0.0;
    // The next dome's volume is drawn ahead and waits until the budget covers it, so
    // that a slow budget builds domes of every size rather than only the smallest.
    let mut pending: Option<f64> = None;
    let mut out = Vec::new();
    for (bin, b0, b1) in crate::stages::time_bins(age, 0.0, present) {
        carry += rate * (b1 - b0);
        let mut j = 0u64;
        loop {
            let draw = 16 + j * 8;
            let volume =
                *pending.get_or_insert_with(|| lo * (hi / lo).powf(rng.uniform(bin, draw)));
            if carry < volume || j >= 1000 {
                break;
            }
            pending = None;
            carry -= volume;
            // V = (2π/3) aspect r³ for a bell of height 2 aspect r.
            let radius =
                (3.0 * volume * 1e9 / (2.0 * std::f64::consts::PI * CRYODOME_ASPECT)).cbrt();
            out.push(Cryodome {
                center: rng.unit_direction(bin, draw + 1),
                radius,
                height: 2.0 * CRYODOME_ASPECT * radius,
                time_gyr: b0 + (b1 - b0) * rng.uniform(bin, draw + 3),
                id: (bin << 32) | j,
            });
            j += 1;
        }
    }
    out
}

/// A dome on the resolved grid at its height now: a bell with streaked flanks, its
/// surface the dome's own, salted.
pub fn cryodome_compute(
    layers: &LayerSet,
    grid: &CubeSphere,
    body_radius: f64,
    dome: &Cryodome,
    height_now: f64,
    stamper: &mut Stamper,
) -> StampOutput {
    let bedrock = layers.bedrock.as_slice();
    let ages = layers.surface_age.as_slice();
    let deep = layers.deep.as_slice();
    let melt = layers.melt.as_slice();
    let noise = GradientNoise::new(dome.id ^ 0xD0CE_0000);
    let streak_frequency = body_radius / (0.2 * dome.radius);
    let time = dome.time_gyr as f32;
    stamper.gather(grid, body_radius, dome.center, dome.radius);
    let mut cells: Vec<CellWrite> = stamper
        .footprint
        .par_iter()
        .map(|&(k, r, _)| {
            let d = grid.center(grid.cell(k));
            let x = (r / dome.radius).min(1.0);
            let bell = (1.0 - x * x).powf(1.5);
            let streaks = CRYODOME_FLANK_ROUGHNESS
                * bell
                * (1.0 - bell)
                * noise.fbm(d, streak_frequency, 2, 2.0, 0.5);
            let mut w = write(k, bedrock[k] as f64 + height_now * (bell + streaks));
            if ages[k] <= time + 1e-6 {
                // The dome's surface by the share of the cell its edifice covers (F131): its
                // salts are the cover the albedo reads on its unit, so it takes the unit and
                // time wherever it lies, as a flow's margin does, its salts the share.
                let share = dome_share(dome, d, body_radius, grid.delta()) as f32;
                if share > 0.0 {
                    w.surface = true;
                    w.age = time;
                    w.unit = UNIT_CRYODOME;
                    w.deep = (1.0 - share) * deep[k];
                    w.melt = (1.0 - share) * melt[k];
                    w.cover = share * CRYODOME_SALT as f32;
                }
            }
            w
        })
        .collect();
    cells.sort_by_key(|c| c.idx);
    StampOutput {
        cells,
        composition: None,
        stain: None,
        steep: None,
        ejecta_age: None,
        compositions: None,
    }
}

/// The share of a cell about `d`, `cell` radians wide, that a cryodome's edifice covers,
/// where its profile stands above a twentieth of its height (F131).
fn dome_share(dome: &Cryodome, d: [f64; 3], body_radius: f64, cell: f64) -> f64 {
    let bell = |e: [f64; 3]| {
        let x = (crate::geom::angle_between(dome.center, e) * body_radius / dome.radius).min(1.0);
        (1.0 - x * x).powf(1.5)
    };
    crate::fields::level_share(&bell, d, 0.05, cell)
}

/// After the stamps: the domes' substrate is partly clean ice, and the floors of
/// craters deep enough to reach the brines carry salts in patches.
pub fn post_pass(
    layers: &mut LayerSet,
    body_radius: f64,
    domes: &[Cryodome],
    craters: &[CraterRecord],
    scaling: &crate::crater::Scaling,
    present: f64,
) {
    let grid = layers.grid;
    if !domes.is_empty() {
        let ages = layers.surface_age.as_slice().to_vec();
        let clean = layers.clean.as_mut_slice();
        crate::crater::with_stamper(&grid, |stamper| {
            for dome in domes {
                stamper.gather(&grid, body_radius, dome.center, dome.radius);
                for &(k, _, _) in &stamper.footprint {
                    // The dome's clean ice by its share of the cell, where nothing younger
                    // has renewed it (F131).
                    if ages[k] > dome.time_gyr as f32 + 1e-3 {
                        continue;
                    }
                    let d = grid.center(grid.cell(k));
                    let share = dome_share(dome, d, body_radius, grid.delta()) as f32;
                    clean[k] = clean[k].max(share * CRYODOME_CLEAN as f32);
                }
            }
        });
    }
    let faculae: Vec<&CraterRecord> = craters
        .iter()
        .filter(|c| {
            c.visible
                && present - c.time_gyr < 4.0 * FACULA_FADE_GYR
                && dimensions(c.diameter, scaling).depth >= FACULA_BRINE_DEPTH
        })
        .collect();
    if faculae.is_empty() {
        return;
    }
    let units = layers.unit_id.as_slice().to_vec();
    let ages = layers.surface_age.as_slice().to_vec();
    let cover = layers.cover.as_mut_slice();
    let noise = GradientNoise::new(0xFAC0_1A00);
    let patch_frequency = body_radius / FACULA_PATCH_SCALE;
    crate::crater::with_stamper(&grid, |stamper| {
        for c in faculae {
            let dim = dimensions(c.diameter, scaling);
            stamper.gather(
                &grid,
                body_radius,
                c.center,
                FACULA_FLOOR_SHARE * dim.radius,
            );
            // The salts lie where the patches rise above their level, drawn by the share
            // of each cell they cover (F127).
            let patch = |p: [f64; 3]| {
                let taper = 1.0
                    - crate::geom::angle_between(c.center, p) * body_radius
                        / (FACULA_FLOOR_SHARE * dim.radius);
                noise.fbm(p, patch_frequency, 2, 2.0, 0.5) + 0.6 * taper
            };
            for &(k, _, _) in &stamper.footprint {
                let floor = units[k] == UNIT_CRATER_FLOOR || units[k] == UNIT_DEEP_FLOOR;
                if !floor || (ages[k] as f64 - c.time_gyr).abs() > 1e-3 {
                    continue;
                }
                let d = grid.center(grid.cell(k));
                let share = crate::fields::level_share(&patch, d, 0.45, grid.delta());
                if share > 0.0 {
                    cover[k] = cover[k].max(share as f32);
                }
            }
        }
    });
}

/// The share of the surface that is chaos matrix.
pub fn chaos_share(layers: &LayerSet) -> f64 {
    let grid = layers.grid;
    let units = layers.unit_id.as_slice();
    let (mut chaos, mut total) = (0.0, 0.0);
    for (k, &u) in units.iter().enumerate() {
        let a = grid.cell_area(grid.cell(k));
        total += a;
        if u == UNIT_CHAOS {
            chaos += a;
        }
    }
    chaos / total.max(1e-300)
}

/// The leading and trailing apexes of a locked moon: the sub-parent point lies on +x and
/// the orbit runs prograde about +z, so the moon moves towards -y and its trailing
/// hemisphere faces +y.
pub const TRAILING_APEX: [f64; 3] = [0.0, 1.0, 0.0];

#[cfg(test)]
mod tests {
    use super::*;

    /// Europa's chaos (F80) covers its share of the surface, at every longitude.
    #[test]
    fn europa_s_chaos_covers_its_share_at_every_longitude() {
        let (g, _) = crate::presets::preset("Europa").unwrap().to_graph();
        let fraction = g.f64("tectonics.chaos_fraction");
        let generations = chaos_population(&g, g.f64("body.age"), 11);
        assert_eq!(generations.len(), CHAOS_GENERATIONS);
        let n = 20_000;
        let mut sectors = [(0usize, 0usize); 6];
        for i in 0..n {
            let d = even_point(i, n);
            let lon = d[1].atan2(d[0]).to_degrees() + 180.0;
            let s = &mut sectors[((lon / 60.0) as usize).min(5)];
            s.1 += 1;
            if generations.iter().any(|c| c.disrupts(d)) {
                s.0 += 1;
            }
        }
        let share = sectors.iter().map(|s| s.0).sum::<usize>() as f64 / n as f64;
        assert!(
            (share / fraction - 1.0).abs() < 0.15,
            "share {share:.3} of {fraction:.3}"
        );
        for s in &sectors {
            assert!(s.0 as f64 / s.1 as f64 > 0.4 * share, "{sectors:?}");
        }
    }

    /// The ocean brings the most heat to the ice at low latitudes, where chaos is most
    /// common (F80).
    #[test]
    fn chaos_favours_low_latitudes() {
        let (g, _) = crate::presets::preset("Europa").unwrap().to_graph();
        let n = 20_000;
        let edge = 40f64.to_radians().sin();
        let (mut low, mut high) = ([0usize; 2], [0usize; 2]);
        for seed in 1..=4 {
            let generations = chaos_population(&g, g.f64("body.age"), seed);
            for i in 0..n {
                let d = even_point(i, n);
                let chaos = generations.iter().any(|c| c.disrupts(d)) as usize;
                let band = if d[2].abs() < edge {
                    &mut low
                } else {
                    &mut high
                };
                band[0] += chaos;
                band[1] += 1;
            }
        }
        let low = low[0] as f64 / low[1] as f64;
        let high = high[0] as f64 / high[1] as f64;
        assert!(
            low > 1.2 * high,
            "{low:.3} within 40 degrees, {high:.3} beyond"
        );
    }

    /// A generation breaks the lid over its warm ice into rafts in a lower matrix (F80):
    /// the rafts keep what they carry, the matrix is renewed at a level the old relief
    /// does not show through and stained by the ocean's material, and the surface
    /// outside is untouched.
    #[test]
    fn a_generation_breaks_the_lid_into_rafts_in_a_matrix() {
        let radius = 100e3;
        let grid = CubeSphere::new(64);
        let cell = grid.delta() * radius;
        let generation = ChaosGeneration::new(0, 0.05, 7, radius, 0.3);
        let mut layers = LayerSet::new(grid);
        for k in 0..grid.cell_count() {
            let d = grid.center(grid.cell(k));
            layers.bedrock.as_mut_slice()[k] = (300.0 * (60.0 * d[0]).sin()) as f32;
            layers.surface_age.as_mut_slice()[k] = 0.01;
            layers.unit_id.as_mut_slice()[k] = UNIT_LINEA;
        }
        let out = chaos_compute(&layers, &grid, cell, 0.06, &generation);
        let stain = out.stain.as_ref().unwrap();
        let region = (0..grid.cell_count())
            .filter(|&k| generation.disrupts(grid.center(grid.cell(k))))
            .count();
        assert_eq!(out.cells.len(), region);
        assert!(
            (0.2..0.45).contains(&(region as f64 / grid.cell_count() as f64)),
            "region {region}"
        );
        let (mut matrix, mut old, mut new) = (0usize, Vec::new(), Vec::new());
        for (w, s) in out.cells.iter().zip(stain) {
            let k = w.idx as usize;
            assert!(generation.disrupts(grid.center(grid.cell(k))));
            if w.unit == UNIT_CHAOS {
                matrix += 1;
                assert_eq!(w.age, 0.05);
                assert!(s[0] > 0.3, "matrix stain {}", s[0]);
                old.push(layers.bedrock.as_slice()[k] as f64);
                new.push(w.bedrock as f64);
            } else {
                assert_eq!(w.unit, UNIT_LINEA);
                assert_eq!(w.age, 0.01);
            }
        }
        let share = matrix as f64 / out.cells.len() as f64;
        assert!((0.2..0.9).contains(&share), "matrix share {share:.3}");
        let spread = |v: &[f64]| {
            let m = v.iter().sum::<f64>() / v.len() as f64;
            (v.iter().map(|x| (x - m).powi(2)).sum::<f64>() / v.len() as f64).sqrt()
        };
        let mean = |v: &[f64]| v.iter().sum::<f64>() / v.len() as f64;
        assert!(
            spread(&new) < 0.5 * spread(&old),
            "matrix relief {:.0} m over {:.0} m",
            spread(&new),
            spread(&old)
        );
        assert!(
            (mean(&new) - mean(&old) + CHAOS_DROP).abs() < 0.3 * CHAOS_DROP,
            "matrix drop {:.0} m",
            mean(&old) - mean(&new)
        );
    }

    #[test]
    fn fallout_falls_off_from_the_pole_and_buries_quickly() {
        let f = Fallout {
            pole: [0.0, 0.0, -1.0],
            rate: 1e5,
            angle: 25f64.to_radians(),
        };
        assert!((f.rate_at([0.0, 0.0, -1.0]) - 1e5).abs() < 1.0);
        let equator = f.rate_at([1.0, 0.0, 0.0]);
        assert!(equator > 1e3 && equator < 1e4, "{equator}");
        assert!(f.burial_gyr([1.0, 0.0, 0.0]) < 1e-5);
        // The frost saturates softly: nearly full at the equator, part of the way on
        // the far side, where the floor of the fallout is the frost threshold.
        assert!(f.frost_at([1.0, 0.0, 0.0]) > 0.9);
        let far = f.frost_at([0.0, 0.0, 1.0]);
        assert!(far > 0.5 && far < 0.75, "{far}");
    }

    /// The loose cover is the fallout of every plume that has run (F64): the present
    /// episode's at the present pole at least, the far share of every episode's far
    /// from the poles, and nothing on a surface that has only just formed.
    #[test]
    fn loose_cover_is_the_sum_of_the_plumes_fallout() {
        let (g, _) = crate::presets::preset("Enceladus").unwrap().to_graph();
        let shell = crate::stages::ice::ShellMap::new(&g);
        let episodes = crate::stages::ice::focusing_episodes(&g, &shell);
        let present = g.f64("body.age");
        let cover = LooseCover::new(&g, &episodes, present).unwrap();
        let pole = shell.pole.unwrap();
        let rate = g.f64("volcanism.plume_fallout_rate");
        let one = rate * FOCUSING_EPISODE_DURATION_GYR;
        let at_pole = cover.thickness_between(pole, 0.0, present);
        assert!(at_pole >= 0.999 * one, "{at_pole} against {one}");
        assert_eq!(cover.thickness_between(pole, present, present), 0.0);
        // The least any old surface carries: the far share of every plume.
        let floor = PLUME_FALLOUT_FAR * one * (episodes.len() + 1) as f64;
        let mut least = f64::INFINITY;
        let mut most: f64 = 0.0;
        for k in 0..2000u64 {
            let d = Rng::new(7).unit_direction(k, 0);
            let h = cover.thickness_between(d, 0.0, present);
            least = least.min(h);
            most = most.max(h);
        }
        eprintln!(
            "enceladus: {} past episodes; loose cover {least:.0} to {most:.0} m on old ground, {at_pole:.0} m at the pole",
            episodes.len()
        );
        assert!(least >= 0.999 * floor, "{least} against {floor}");
        // The same order as the pit chains of its cratered plains have it (Martin et al.
        // 2017: 250 m on average, to 700 m), and two to three times over: nine plumes'
        // fallout, where the real body's history is not known.
        assert!(least > 100.0 && most < 3000.0, "{least} to {most}");
        // A crack there drains a few hundred metres of it.
        assert!(
            (300.0..500.0).contains(&cover.drained(most)),
            "{}",
            cover.drained(most)
        );
    }
}
