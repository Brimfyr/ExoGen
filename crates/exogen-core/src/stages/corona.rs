// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! Coronae (PLAN.md 4.3, F133): the ground a lid breaks over the upwellings beneath it.
//!
//! A heating episode switches convection on under a lid, and the upwellings of the new
//! convection lift and stretch it: Miranda's coronae stand over upwellings in its ice that
//! tidal heating drove (Hammond and Barr 2014), as Venus's stand over plumes and diapirs
//! under its rock. Each upwelling holds a cell of the convection, bounded by the bisectors
//! with its neighbours, so that a cell is a polygon, as Miranda's coronae are, "polygonal
//! to ovoidal" (Inverness a trapezoid with a chevron at one corner). The corona is the part
//! of its cell where the flow under the lid stretches it, out to where the flow turns to
//! shortening it. Its outer belt is a ring of fault segments concentric with its outline,
//! the tilt blocks of Arden Corona (Pappalardo et al. 1997), whose scarps face outward, so
//! that the belt stands above the old ground beyond its bounding faults; its core is
//! crossed by sets of its own, where the lid stretched every way over the upwelling broke
//! in more than one direction. The faults are the grooves' (F132), spaced by the lid of the
//! time and thinner over the upwelling than at its edge, and the corona renews the ground it
//! breaks as a lane does (F60, F61).

use rayon::prelude::*;

use crate::constants::*;
use crate::crater::{CellWrite, StampOutput, Stamper, with_stamper};
use crate::fields::LayerSet;
use crate::geom::{angle_between, azimuth_from, dot, local_basis, normalize, offset};
use crate::grid::CubeSphere;
use crate::noise::GradientNoise;
use crate::params::ParamGraph;
use crate::rng::Rng;
use crate::stages::ice_shell::{LaneFabric, LaneHit, fault_set, flexural_parameter};
use crate::thermal::ThermalHistory;

/// A convecting layer under a lid at a time (F133).
#[derive(Clone, Copy, Debug)]
pub struct Convection {
    /// The conductive lid over the layer and the layer's depth, m.
    pub lid: f64,
    pub depth: f64,
    /// The upwellings' spacing at the lid's base, m, and how many the sphere holds there.
    pub spacing: f64,
    pub upwellings: f64,
}

impl Convection {
    /// An ice layer `ice_layer` m thick on rock under a lid `lid` m thick, on a body of
    /// `radius` m: none where the lid reaches down to the rock. Upwellings stand twice
    /// the layer's depth apart at its mid-depth, which widens to the lid's base as the
    /// radius does, and as many as hexagonal cells of that spacing cover the sphere there.
    pub fn ice(radius: f64, ice_layer: f64, lid: f64) -> Option<Self> {
        let depth = ice_layer.min(radius) - lid;
        if depth <= 0.0 || lid <= 0.0 {
            return None;
        }
        let top = radius - lid;
        let mid = top - 0.5 * depth;
        let spacing = 2.0 * CORONA_CELL_ASPECT * depth * top / mid.max(1.0);
        let upwellings = 8.0 * std::f64::consts::PI / 3f64.sqrt() * (top / spacing).powi(2);
        Some(Self {
            lid,
            depth,
            spacing,
            upwellings,
        })
    }
}

/// The share of an upwelling's stretching a lid takes in faults of its own (F133): all
/// of it where the lid bends over a shorter length than the upwelling is wide (its
/// flexural parameter against `width`, the distance from the upwelling to its cell's
/// edge), falling where the lid's bending spans the upwelling and spreads its lift as a
/// broad swell.
pub fn localized(flexural: f64, width: f64) -> f64 {
    1.0 / (1.0 + (flexural / width.max(1.0)).powf(CORONA_FLEXURE_POWER))
}

/// The heating episode that breaks an ice mantle's lid in coronae (F133): how strongly
/// the tidal heating drives it, above the floor that renews an icy surface (F44), and the
/// convection at the episode's thinnest lid, with the lid's flexural parameter there and
/// the strength of the coronae over its mean upwelling.
#[derive(Clone, Copy, Debug)]
pub struct IceEpisode {
    pub drive: f64,
    pub convection: Convection,
    pub flexural: f64,
    pub strength: f64,
}

impl IceEpisode {
    pub fn new(
        history: &ThermalHistory,
        radius: f64,
        ice_layer: f64,
        gravity: f64,
    ) -> Option<Self> {
        let peak = history
            .samples
            .iter()
            .map(|s| s.tidal_flux)
            .fold(0.0, f64::max);
        let drive = ((peak - ICE_RESURFACING_FLUX_FLOOR) / ICE_RESURFACING_FLUX).clamp(0.0, 1.0);
        if drive <= 0.0 {
            return None;
        }
        let thinnest = history
            .samples
            .iter()
            .filter(|s| s.tidal_flux > ICE_RESURFACING_FLUX_FLOOR)
            .min_by(|a, b| a.lithosphere.total_cmp(&b.lithosphere))?;
        let convection = Convection::ice(radius, ice_layer, thinnest.lithosphere)?;
        let flexural = flexural_parameter(thinnest.elastic, gravity);
        Some(Self {
            drive,
            convection,
            flexural,
            strength: drive * localized(flexural, 0.5 * convection.spacing),
        })
    }
}

/// A corona over an upwelling (F133).
#[derive(Clone, Debug)]
pub struct Corona {
    pub id: u64,
    /// The upwelling's place, and east and north there.
    pub center: [f64; 3],
    frame: ([f64; 3], [f64; 3]),
    /// The bisectors with the neighbouring upwellings that bound its cell: each the unit
    /// normal of the bisecting plane, towards the upwelling, and the angle from the
    /// upwelling to that plane, rad.
    bisectors: Vec<([f64; 3], f64)>,
    /// The angle from the upwelling to the farthest reach of its outline, its wander
    /// included, rad; and the outline's mean radius, m, the length round its belt.
    extent: f64,
    outline_radius: f64,
    /// The upwelling's width: the mean distance from it to its cell's edge, m.
    pub width: f64,
    /// How fully the lid broke over it, 0 to 1: its faults' relief and their renewal of
    /// the ground go with it.
    pub strength: f64,
    /// Its faults' spacing where the lid has its mean thickness, m.
    pub spacing: f64,
    pub time_gyr: f64,
    /// How clean its ice is, 0 to 1.
    pub clean: f64,
    /// The strike of its core's first set of faults, rad from east.
    pub strike: f64,
}

/// A corona's fields of unresolved variation: its faults' sets and wander, the roughness
/// under them, its outline's wander, and the draws of its fault segments.
pub struct CoronaNoise {
    grooves: GradientNoise,
    rough: GradientNoise,
    outline: [GradientNoise; 3],
    draws: Rng,
}

impl CoronaNoise {
    pub fn new(c: &Corona) -> Self {
        let id = c.id ^ 0xC020_0000_0000;
        Self {
            grooves: GradientNoise::new(id ^ 0x6700),
            rough: GradientNoise::new(id ^ 0x0C4A),
            outline: [
                GradientNoise::new(id ^ 0x0071),
                GradientNoise::new(id ^ 0x0072),
                GradientNoise::new(id ^ 0x0073),
            ],
            draws: Rng::new(id ^ 0xD0A5),
        }
    }
}

/// A corona's faults' middling relief as a share of their spacing (F133). The upwelling
/// stretches the lid for as long as the episode lasts, far longer than its faults take to
/// rotate, so its blocks tilt until their faults lock: from the dip they broke at to the
/// dip at which friction locks a normal fault, new faults then taking the strain at the
/// same tilt. A block tilted so stands its width times the sine of the tilt from crest to
/// trough, which is the middling fault's trough below its plateau (one and the mean
/// depth, F132).
pub fn relief_share() -> f64 {
    (CORONA_FAULT_DIP_DEG - CORONA_LOCKUP_DIP_DEG).to_radians().sin() / (1.0 + GROOVE_DEPTH_MEAN)
}

fn smoothstep(x: f64) -> f64 {
    let x = x.clamp(0.0, 1.0);
    x * x * (3.0 - 2.0 * x)
}

/// The share of the way to a bisector a point stands at: 0 at the upwelling and 1 on the
/// plane, linear in the distance along the sphere.
fn way(d: [f64; 3], bisector: ([f64; 3], f64)) -> f64 {
    1.0 - dot(d, bisector.0).clamp(-1.0, 1.0).asin() / bisector.1.max(1e-9)
}

impl Corona {
    /// The corona over the upwelling at `center`, bounded by its neighbours at `others`
    /// (the nearest of them; the antipode for a lone upwelling), with its strength,
    /// spacing, time and the rest to be set.
    pub fn new(id: u64, center: [f64; 3], others: &[[f64; 3]], body_radius: f64) -> Self {
        let mut near: Vec<([f64; 3], f64)> = others
            .iter()
            .map(|&o| (o, angle_between(center, o)))
            .filter(|o| o.1 > 1e-6)
            .collect();
        near.sort_by(|a, b| a.1.total_cmp(&b.1));
        near.truncate(16);
        let bisectors: Vec<([f64; 3], f64)> = near
            .iter()
            .map(|&(o, _)| {
                let n = normalize([center[0] - o[0], center[1] - o[1], center[2] - o[2]]);
                (n, dot(center, n).clamp(-1.0, 1.0).asin())
            })
            .collect();
        let mut corona = Self {
            id,
            center,
            frame: local_basis(center),
            bisectors,
            extent: 0.0,
            outline_radius: 0.0,
            width: 0.0,
            strength: 0.0,
            spacing: 1.0,
            time_gyr: 0.0,
            clean: 1.0,
            strike: 0.0,
        };
        // The cell's edge and the outline round it, ray by ray from the upwelling; a
        // bisector farther than the cell reaches never bounds it.
        let rays = 64;
        let (mut edge_max, mut out_max, mut out_sum, mut edge_sum) = (0.0f64, 0.0f64, 0.0, 0.0);
        for i in 0..rays {
            let phi = 2.0 * std::f64::consts::PI * i as f64 / rays as f64;
            let edge = corona.ray_to(phi, 1.0);
            let outline = corona.ray_to(phi, CORONA_OUTLINE_SHARE);
            edge_max = edge_max.max(edge);
            edge_sum += edge;
            out_max = out_max.max(outline);
            out_sum += outline;
        }
        corona.bisectors.retain(|b| b.1 <= 1.05 * edge_max);
        corona.extent = (out_max * (1.0 + 3.0 * CORONA_OUTLINE_WANDER)).min(std::f64::consts::PI);
        corona.outline_radius = out_sum / rays as f64 * body_radius;
        corona.width = edge_sum / rays as f64 * body_radius;
        corona
    }

    /// The angle from the upwelling, along azimuth `phi`, at which the sharp share of the
    /// way to the cell's edge reaches `share`, rad.
    fn ray_to(&self, phi: f64, share: f64) -> f64 {
        let (east, north) = self.frame;
        let at = |theta: f64| -> f64 {
            let d = offset(self.center, east, north, theta, phi);
            self.bisectors
                .iter()
                .map(|&b| way(d, b))
                .fold(f64::NEG_INFINITY, f64::max)
        };
        let (mut lo, mut hi) = (0.0, std::f64::consts::PI);
        for _ in 0..40 {
            let mid = 0.5 * (lo + hi);
            if at(mid) < share {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    }

    /// The point as the outline's wander moves it.
    fn wandered(&self, d: [f64; 3], body_radius: f64, noise: &CoronaNoise) -> [f64; 3] {
        let scale = (CORONA_OUTLINE_WANDER_SCALE * self.outline_radius).max(1.0);
        let f = body_radius / scale;
        let a = CORONA_OUTLINE_WANDER * self.outline_radius / body_radius;
        let w: Vec<f64> = noise
            .outline
            .iter()
            .map(|n| n.fbm(d, f, 3, 2.0, 0.5))
            .collect();
        normalize([d[0] + a * w[0], d[1] + a * w[1], d[2] + a * w[2]])
    }

    /// Where a point stands in the corona's cell: 0 at the upwelling and 1 on the cell's
    /// edge, the share of the way to the nearest bisector, its corners rounded; and the
    /// length one unit of that share stands for there, m, the bisectors' distances
    /// weighted as the rounding weights them.
    fn position(&self, d: [f64; 3], body_radius: f64, noise: &CoronaNoise) -> (f64, f64) {
        let d = self.wandered(d, body_radius, noise);
        let k = CORONA_CORNER_SHARPNESS;
        let best = self
            .bisectors
            .iter()
            .map(|&b| way(d, b))
            .fold(f64::NEG_INFINITY, f64::max);
        let (mut sum, mut length) = (0.0, 0.0);
        for &b in &self.bisectors {
            let w = (k * (way(d, b) - best)).exp();
            sum += w;
            length += w * b.1;
        }
        (best + sum.ln() / k, length / sum * body_radius)
    }

    /// The outline's mean radius, m.
    pub fn outline_radius(&self) -> f64 {
        self.outline_radius
    }

    /// The length below which the corona takes the old relief inside it, m.
    pub fn smooth_length(&self) -> f64 {
        CORONA_SMOOTH_SHARE * self.outline_radius
    }

    /// The lid the corona broke, m: the layer whose necking spaced its faults (F61).
    pub fn lid(&self) -> f64 {
        self.spacing / (GROOVE_NECKING_WAVELENGTH * GROOVE_RIFT_LID_SHARE)
    }

    /// The width of the border zone across its outline on a grid of a cell, m: its
    /// bounding faults cut the whole lid, as a lane's do.
    pub fn border(&self, cell: f64) -> f64 {
        (GROOVE_BORDER_LIDS * self.lid()).max(GROOVE_EDGE_CELLS * cell.max(1.0))
    }

    /// The spacing of its faults on a grid of a cell, m: no finer than the grid draws.
    fn faults_on(&self, cell: f64) -> f64 {
        self.spacing.max(GROOVE_MIN_CELLS * cell)
    }

    /// How far from the upwelling the corona reaches on a grid of a cell, m: its outline
    /// and the damage zone beyond it.
    pub fn reach(&self, body_radius: f64, cell: f64) -> f64 {
        self.extent * body_radius
            + (GROOVE_DAMAGE_SPACINGS * self.faults_on(cell)).max(0.5 * self.border(cell))
            + 2.0 * cell
    }

    pub fn cap_radius(&self, body_radius: f64, cell: f64) -> f64 {
        self.reach(body_radius, cell) / body_radius
    }

    /// The corona's cover of a cell, 0 to 1: its strength, ramped across the border zone
    /// at its outline.
    fn cover_with(&self, d: [f64; 3], body_radius: f64, cell: f64, noise: &CoronaNoise) -> f64 {
        let (rho, unit) = self.position(d, body_radius, noise);
        let inward = (CORONA_OUTLINE_SHARE - rho) * unit;
        self.strength * smoothstep(inward / self.border(cell) + 0.5)
    }

    /// The cells the corona covers on a grid, with its cover of each.
    pub fn covers(
        &self,
        grid: &CubeSphere,
        body_radius: f64,
        cell: f64,
        stamper: &mut Stamper,
    ) -> Vec<(usize, f64)> {
        let noise = CoronaNoise::new(self);
        stamper.gather(
            grid,
            body_radius,
            self.center,
            self.extent * body_radius + self.border(cell),
        );
        let cells: Vec<usize> = stamper.footprint.iter().map(|f| f.0).collect();
        cells
            .into_par_iter()
            .filter_map(|idx| {
                let c = self.cover_with(grid.center(grid.cell(idx)), body_radius, cell, &noise);
                (c > 0.0).then_some((idx, c))
            })
            .collect()
    }

    /// What the corona makes of a cell (F133), as a lane's fabric: its cover; inside its
    /// outline the relief of its faults, an outer belt of segments concentric with the
    /// outline and a core crossed by two sets, with the footwalls of its outward-facing
    /// bounding faults up; beyond the outline their hanging walls down and the belt's
    /// faults running on as open fractures, fewer and shorter the farther out; the tone of
    /// its ice, darker in the troughs, and the clean ice the fractures show.
    pub fn fabric_at(
        &self,
        d: [f64; 3],
        body_radius: f64,
        cell: f64,
        noise: &CoronaNoise,
    ) -> LaneFabric {
        let pi = std::f64::consts::PI;
        let tau = 2.0 * pi;
        let s = CORONA_OUTLINE_SHARE;
        let (rho, unit) = self.position(d, body_radius, noise);
        let inward = (s - rho) * unit;
        let spacing = self.faults_on(cell);
        let damage = GROOVE_DAMAGE_SPACINGS * spacing;
        if inward <= -damage.max(0.5 * self.border(cell)) {
            return LaneFabric::default();
        }
        let cover = self.strength * smoothstep(inward / self.border(cell) + 0.5);
        let amplitude = relief_share() * self.spacing * self.strength;
        // The lid's thickness across the cell, as a share of its mean: thinnest over the
        // upwelling, where the flow brings up the most heat.
        let thin = |r: f64| 1.0 / (1.0 + CORONA_FLUX_CONTRAST * (pi * r.clamp(0.0, 1.0)).cos());
        let g = &noise.grooves;
        let set_frequency = body_radius / (GROOVE_SET_SPACINGS * spacing);
        let sets = g.fbm(d, set_frequency, 2, 2.0, 0.5);
        let tones = g.fbm(
            [d[0] - 2.3, d[1] + 5.9, d[2] + 0.7],
            set_frequency,
            2,
            2.0,
            0.5,
        );
        let strength = GROOVE_SLACK + (1.0 - GROOVE_SLACK) * (sets.abs() / 0.3).min(1.0);
        // The belt's faults run parallel to the outline, their spacing as the lid's: the
        // distance in from the outline in units of the mean spacing is the integral of
        // the spacing's share across the cell.
        let across = unit
            * ((s - rho)
                + CORONA_FLUX_CONTRAST / pi * ((pi * s).sin() - (pi * rho.clamp(0.0, 1.0)).sin()));
        let phi = azimuth_from(self.center, self.frame.0, self.frame.1, d);
        let period = tau * self.outline_radius.max(1.0);
        let along = (phi + pi) / tau * period;
        let out = (-inward / damage).clamp(0.0, 1.0);
        let core_edge = s * (1.0 - CORONA_BELT_SHARE);
        let belt = smoothstep((rho - core_edge) / (CORONA_CORE_BLEND * s) + 0.5);
        let level =
            (GROOVE_SEGMENT_LEVEL - (GROOVE_SEGMENT_LEVEL + GROOVE_DAMAGE_LEVEL_DROP) * out) * belt;
        let base = fault_set(
            across,
            along,
            d,
            body_radius,
            spacing,
            sets,
            level,
            0,
            0.0,
            g,
            &noise.draws,
            Some(period),
            GROOVE_FINE_SKEW,
        );
        let mut shape = base.shape;
        let mut lag = base.gate * base.trough * base.trough * base.depth.min(2.0);
        // The core's two sets, straight across it, where the belt gives way: each fault
        // cuts what stood there before across its own trough, as a later pulse's do in a
        // lane (F63).
        if belt < 1.0 {
            let core_spacing = spacing * thin(0.5 * core_edge);
            let cz = dot(d, self.center);
            let u = dot(d, self.frame.0).atan2(cz) * body_radius;
            let v = dot(d, self.frame.1).atan2(cz) * body_radius;
            let level = CORONA_CORE_LEVEL * (1.0 - belt);
            // The sets turn across the core, so that they bend.
            let bend = CORONA_CORE_TURN
                * g.fbm(
                    [d[0] + 4.1, d[1] - 2.6, d[2] + 8.3],
                    body_radius / (CORONA_CORE_TURN_SCALE * self.outline_radius).max(1.0),
                    2,
                    2.0,
                    0.5,
                )
                / 0.5;
            let strike = self.strike + bend;
            for (pulse, turn) in [(1usize, strike), (2, strike + CORONA_CROSS_ANGLE)] {
                let (sin, cos) = turn.sin_cos();
                let set = fault_set(
                    -u * sin + v * cos,
                    u * cos + v * sin,
                    d,
                    body_radius,
                    core_spacing,
                    sets,
                    level,
                    pulse,
                    0.0,
                    g,
                    &noise.draws,
                    None,
                    0.0,
                );
                let own = set.gate * set.trough;
                shape += (set.shape - shape) * own;
                lag += (set.trough * set.trough * set.depth.min(2.0) - lag) * own;
            }
        }
        // The faults die out over the last stretch in to the outline, leaving the plateau
        // they cut, as a lane's before its bounding faults (F132).
        let margin = smoothstep(inward / (GROOVE_MARGIN_SPACINGS * spacing));
        shape = base.plateau + margin * (shape - base.plateau);
        let rough_octaves = ((GROOVE_ROUGH_SPACINGS * spacing) / (2.0 * cell.max(1.0)))
            .log2()
            .ceil()
            .clamp(1.0, 10.0) as u32;
        let rough = GROOVE_ROUGH_SHARE
            * amplitude
            * noise.rough.fbm(
                d,
                body_radius / (GROOVE_ROUGH_SPACINGS * spacing),
                rough_octaves,
                2.0,
                2f64.powf(-ROUGHNESS_HURST),
            )
            / 0.4;
        // The bounding faults face outward (Pappalardo et al. 1997): the corona's side,
        // the footwall, up; the ground beyond, the hanging wall, down; the throw coming and
        // going along the outline as the faults' segments do.
        let throw = {
            let run = noise.rough.fbm(
                [d[0] + 7.7, d[1] + 1.3, d[2] - 0.59],
                body_radius / (GROOVE_EDGE_SEGMENT_SPACINGS * spacing),
                2,
                2.0,
                0.5,
            );
            let x = (GROOVE_BOUNDARY_LEVEL + run).clamp(0.0, 1.0);
            GROOVE_BOUNDARY_THROW * amplitude * x * x * (3.0 - 2.0 * x)
        };
        let up = 0.5 * throw * (-inward.abs() / (GROOVE_FLANK_SPACINGS * spacing)).exp();
        let down = -0.5 * throw * (-inward.abs() / spacing).exp();
        // Beyond the outline the belt's faults run on as open fractures.
        let fracture_depth = GROOVE_FRACTURE_SHARE * amplitude;
        let open = base.gate * base.trough.powi(3);
        let grooves = strength * amplitude * thin(rho) * shape + rough;
        let set_tone = 1.0 - GROOVE_TONE_SETS * (0.5 + 0.5 * (tones / 0.5).clamp(-1.0, 1.0));
        let lag = GROOVE_TONE_TROUGH * strength * lag.clamp(0.0, 1.0);
        LaneFabric {
            cover,
            relief: grooves + up,
            outer: (1.0 - out) * (down - fracture_depth * open),
            tone: self.clean * set_tone * (1.0 - lag),
            fracture: GROOVE_DAMAGE_CLEAN * self.clean * self.strength * (1.0 - out) * open,
        }
    }
}

/// The upwellings' places: each the farthest from those already placed of a few tries,
/// so that they spread over the sphere as convection's cells do, without a lattice.
fn place_upwellings(count: usize, rng: &Rng) -> Vec<[f64; 3]> {
    let mut out: Vec<[f64; 3]> = Vec::with_capacity(count);
    for i in 0..count as u64 {
        let mut best = (f64::NEG_INFINITY, [0.0, 0.0, 1.0]);
        for t in 0..CORONA_PLACEMENT_TRIES {
            let d = rng.unit_direction(i, 100 + t);
            let gap = out
                .iter()
                .map(|&c| angle_between(c, d))
                .fold(std::f64::consts::PI, f64::min);
            if gap > best.0 {
                best = (gap, d);
            }
        }
        out.push(best.1);
    }
    out
}

/// The coronae an ice mantle's convection breaks in the lid over it in a tidal heating
/// episode (F133): one over each upwelling, its time drawn within the episode, its
/// faults spaced by the lid of that time, and its strength the episode's drive times the
/// share of the upwelling's stretching the lid of that time takes in its own faults over
/// that upwelling's own width.
pub fn ice_coronae(graph: &ParamGraph, thermal: &ThermalHistory) -> Vec<Corona> {
    let count = graph.f64("volcanism.corona_count").round().max(0.0) as usize;
    if count == 0 || graph.str("interior.tectonic_regime") != "ice_shell_rock" {
        return Vec::new();
    }
    let body_radius = graph.f64("body.radius");
    let gravity = graph.f64("body.gravity");
    let ice_layer = graph.f64("body.ice_layer_thickness");
    let Some(episode) = IceEpisode::new(thermal, body_radius, ice_layer, gravity) else {
        return Vec::new();
    };
    let age = graph.f64("body.age");
    let mid = (age - graph.f64("volcanism.corona_age")).clamp(0.0, age);
    let width = graph.f64("volcanism.corona_duration").max(0.0);
    let rng = Rng::new(graph.f64("seed.tectonics") as u64 ^ 0xC020_1CE0);
    let centres = place_upwellings(count, &rng);
    // The lid at each corona's own time, from the thermal relations at that time rather
    // than the history's samples, which a short episode falls between.
    let inputs = crate::thermal::ThermalInputs::from_graph(graph);
    let (c_lo, c_hi) = GROOVE_CLEAN_RANGE;
    (0..count)
        .map(|i| {
            let k = i as u64;
            let others: Vec<[f64; 3]> = if count == 1 {
                let c = centres[0];
                vec![[-c[0], -c[1], -c[2]]]
            } else {
                (0..count).filter(|&j| j != i).map(|j| centres[j]).collect()
            };
            let mut corona = Corona::new(k, centres[i], &others, body_radius);
            let time = (mid + (rng.uniform(k, 10) - 0.5) * width).clamp(0.0, age);
            let state = crate::thermal::state_at(&inputs, time);
            let elastic = state.elastic.max(1e3);
            let convecting = Convection::ice(body_radius, ice_layer, state.lithosphere).is_some();
            corona.time_gyr = time;
            corona.spacing = GROOVE_NECKING_WAVELENGTH * GROOVE_RIFT_LID_SHARE * elastic;
            corona.strength = if convecting {
                episode.drive * localized(flexural_parameter(elastic, gravity), corona.width)
            } else {
                0.0
            };
            corona.clean = c_lo + (c_hi - c_lo) * rng.uniform(k, 11);
            corona.strike = 2.0 * std::f64::consts::PI * rng.uniform(k, 12);
            corona
        })
        .filter(|c| c.strength > CORONA_STRENGTH_FLOOR)
        .collect()
}

/// The coronae on the resolved grid, before the finish stamps anything: each corona's
/// hold of the cells about it, for the renewal to take in the order the terrains formed
/// (F60, F61).
pub fn corona_hits(
    grid: &CubeSphere,
    body_radius: f64,
    cell: f64,
    coronae: &[Corona],
) -> Vec<Vec<LaneHit>> {
    coronae
        .iter()
        .map(|c| {
            let noise = CoronaNoise::new(c);
            let cells: Vec<usize> = with_stamper(grid, |stamper| {
                stamper.gather(grid, body_radius, c.center, c.reach(body_radius, cell));
                stamper.footprint.iter().map(|f| f.0).collect()
            });
            let mut hits: Vec<LaneHit> = cells
                .into_par_iter()
                .filter_map(|k| {
                    let f = c.fabric_at(grid.center(grid.cell(k)), body_radius, cell, &noise);
                    (f.cover > 0.0 || f.fracture > 0.004).then_some(LaneHit {
                        k: k as u32,
                        cover: f.cover as f32,
                        tone: f.tone as f32,
                        fracture: f.fracture as f32,
                    })
                })
                .collect();
            hits.sort_by_key(|h| h.k);
            hits
        })
        .collect()
}

/// A corona's fabric on the resolved grid at its time (F133), as a lane's (F61): inside
/// its outline it renews the ground, whatever fine relief stood there giving way by its
/// cover to the floor (`base`, the relief before anything was stamped) with its faults
/// on it; beyond the outline the hanging walls' drop and the fractures are added to the
/// ground as it stands. A cell something younger renewed is left alone.
#[allow(clippy::too_many_arguments)]
pub fn corona_compute(
    layers: &LayerSet,
    grid: &CubeSphere,
    body_radius: f64,
    corona: &Corona,
    cell: f64,
    renewed: &[f32],
    base: &[f32],
    stamper: &mut Stamper,
) -> StampOutput {
    let bedrock = layers.bedrock.as_slice();
    let noise = CoronaNoise::new(corona);
    stamper.gather(
        grid,
        body_radius,
        corona.center,
        corona.reach(body_radius, cell),
    );
    let footprint: Vec<usize> = stamper.footprint.iter().map(|f| f.0).collect();
    let mut cells: Vec<CellWrite> = footprint
        .into_par_iter()
        .filter_map(|k| {
            if renewed[k] as f64 > corona.time_gyr + 1e-3 {
                return None;
            }
            let f = corona.fabric_at(grid.center(grid.cell(k)), body_radius, cell, &noise);
            if f.cover <= 0.0 && f.outer.abs() < 0.5 {
                return None;
            }
            let h = bedrock[k] as f64;
            Some(CellWrite {
                idx: k as u32,
                bedrock: (h + (base[k] as f64 + f.relief - h) * f.cover + f.outer * (1.0 - f.cover))
                    as f32,
                surface: false,
                age: 0.0,
                unit: 0,
                deep: 0.0,
                melt: 0.0,
                cover: 0.0,
            })
        })
        .collect();
    cells.sort_by_key(|c| c.idx);
    if std::env::var_os("EXOGEN_DEBUG").is_some() {
        let (lo, hi) = cells.iter().fold((f32::INFINITY, f32::NEG_INFINITY), |(lo, hi), w| {
            let dh = w.bedrock - bedrock[w.idx as usize];
            (lo.min(dh), hi.max(dh))
        });
        eprintln!(
            "corona {} at {:.3} Gyr: {} cells written, change {lo:.0} to {hi:.0} m",
            corona.id,
            corona.time_gyr,
            cells.len()
        );
    }
    StampOutput {
        cells,
        composition: None,
        stain: None,
        steep: None,
        ejecta_age: None,
        compositions: None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tetrahedron() -> Vec<[f64; 3]> {
        let a = 1.0 / 3f64.sqrt();
        vec![[a, a, a], [a, -a, -a], [-a, a, -a], [-a, -a, a]]
    }

    /// The cell's edge is the bisector with each neighbour: a point halfway to a
    /// neighbour stands at one, the upwelling at nought, and the outline at its share of
    /// the way.
    #[test]
    fn a_cell_is_bounded_by_the_bisectors() {
        let pts = tetrahedron();
        let r = 2.358e5;
        let c = Corona::new(0, pts[0], &pts[1..], r);
        assert_eq!(c.bisectors.len(), 3);
        let half = normalize([
            pts[0][0] + pts[1][0],
            pts[0][1] + pts[1][1],
            pts[0][2] + pts[1][2],
        ]);
        let sharp = |d: [f64; 3]| {
            c.bisectors
                .iter()
                .map(|&b| way(d, b))
                .fold(f64::NEG_INFINITY, f64::max)
        };
        assert!((sharp(half) - 1.0).abs() < 1e-9, "{}", sharp(half));
        assert!(sharp(pts[0]).abs() < 1e-9);
        // A tetrahedron's cell reaches 54.7 degrees to its faces and 70.5 to its corners.
        let face = 0.5 * angle_between(pts[0], pts[1]);
        assert!((face.to_degrees() - 54.74).abs() < 0.1);
        let w = c.width / r;
        assert!(w > face && w < 70.6f64.to_radians(), "{}", w.to_degrees());
        let outline = c.outline_radius / r;
        assert!(
            outline > CORONA_OUTLINE_SHARE * face && outline < CORONA_OUTLINE_SHARE * 1.3,
            "{}",
            outline.to_degrees()
        );
    }

    /// Miranda's resonance heats its ice to convect under a lid a few tens of kilometres
    /// thick; the upwellings are few, as its coronae are, and the lid breaks over them.
    #[test]
    fn mirandas_episode_breaks_its_lid_over_a_few_upwellings() {
        let (g, _) = crate::presets::preset("Miranda").unwrap().to_graph();
        let n = g.f64("volcanism.corona_count");
        assert!((3.0..=8.0).contains(&n), "coronae {n}");
        let strength = g.f64("volcanism.corona_strength");
        assert!(strength > 0.8, "strength {strength}");
        let age = g.f64("volcanism.corona_age");
        assert!((0.8..=1.2).contains(&age), "age {age}");
        let thermal = ThermalHistory::from_graph(&g);
        let coronae = ice_coronae(&g, &thermal);
        assert_eq!(coronae.len() as f64, n);
        for c in &coronae {
            // Each 200 to 400 km across, its faults 3 to 15 km apart: Arden's about 5,
            // Elsinore's wider, as a corona broken while the heating waned has them.
            let across = 2.0 * c.outline_radius / 1e3;
            eprintln!(
                "corona {across:.0} km across, faults {:.1} km apart, strength {:.2}, {:.3} Gyr",
                c.spacing / 1e3,
                c.strength,
                c.time_gyr
            );
            assert!((150.0..=450.0).contains(&across), "{across} km");
            assert!((3e3..=15e3).contains(&c.spacing), "{}", c.spacing);
        }
    }

    /// Bodies with no tidal heating episode, or with an ocean under the ice, have none.
    #[test]
    fn no_coronae_without_an_episode_or_over_an_ocean() {
        for name in ["Ganymede", "Europa", "Enceladus", "Callisto", "Ceres", "Pluto"] {
            let (g, _) = crate::presets::preset(name).unwrap().to_graph();
            assert_eq!(g.f64("volcanism.corona_count"), 0.0, "{name}");
        }
    }

    /// A corona's belt holds faults concentric with its outline, in segments that close
    /// round it without a seam; its footwall stands up at the outline and the ground
    /// beyond is let down.
    #[test]
    fn a_coronas_belt_follows_its_outline() {
        let pts = tetrahedron();
        let r = 2.358e5;
        let mut c = Corona::new(0, pts[0], &pts[1..], r);
        c.strength = 1.0;
        c.spacing = 5e3;
        c.time_gyr = 3.5;
        let noise = CoronaNoise::new(&c);
        let cell = 500.0;
        // Round the belt, at a fixed share of the way in, the cover holds and the faults
        // cut troughs: the relief varies by a good part of the amplitude.
        let (east, north) = c.frame;
        let mut reliefs = Vec::new();
        let mut covers = Vec::new();
        for i in 0..720 {
            let phi = 2.0 * std::f64::consts::PI * i as f64 / 720.0;
            let theta = c.ray_to(phi, 0.7 * CORONA_OUTLINE_SHARE);
            let d = offset(c.center, east, north, theta, phi);
            let f = c.fabric_at(d, r, cell, &noise);
            reliefs.push(f.relief);
            covers.push(f.cover);
        }
        assert!(covers.iter().all(|&v| v > 0.9));
        let amplitude = relief_share() * c.spacing;
        let lo = reliefs.iter().cloned().fold(f64::INFINITY, f64::min);
        let hi = reliefs.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        assert!(hi - lo > amplitude, "{lo}..{hi}");
        // Far outside, nothing.
        let far = offset(c.center, east, north, 1.5 * c.extent, 0.3);
        let f = c.fabric_at(far, r, cell, &noise);
        assert_eq!(f.cover, 0.0);
        assert_eq!(f.relief, 0.0);
    }
}
