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
use crate::flexure::{Plate, Radial, deflect, surface_tension};
use crate::geom::{angle_between, azimuth_from, dot, local_basis, normalize, offset};
use crate::grid::CubeSphere;
use crate::mantle::MantleHistory;
use crate::noise::GradientNoise;
use crate::params::{ParamGraph, Reads};
use crate::rng::Rng;
use crate::stages::ice_shell::{LaneFabric, LaneHit, fault_set, flexural_parameter};
use crate::stages::plume::Rise;
use crate::stages::stress::Limbs;
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
    /// The share of the tilt to lock-up its blocks reach (F145): all of it over an ice
    /// mantle's upwelling, which stretches the lid far longer than its faults take to
    /// rotate; for a rock corona the strain its lid's bending gave.
    pub tilt: f64,
    /// A rock corona's diapir, the melt it made and the profile its loads give (F145).
    pub head: Option<Head>,
    /// How far its outline wanders, as a share of its radius.
    wander: f64,
}

/// A rock corona's diapir and what it left (F145, F146).
#[derive(Clone, Debug)]
pub struct Head {
    /// The flattened head's radius and thickness, m, and its excess temperature over the
    /// mantle about it at its axis, K.
    pub radius: f64,
    pub thickness: f64,
    pub excess: f64,
    /// How long its conduit feeds it, Gyr, and the e-folding time of its buoyancy after,
    /// Gyr: the time its excess takes to diffuse across its half-thickness.
    pub feeding: f64,
    pub cooling: f64,
    /// The melt it made at its pocket's centre, as a column, m.
    pub magma: f64,
    /// The lid's elastic thickness over it at its time, m.
    pub elastic: f64,
    /// The profiles its loads give the surface against the distance from its centre, m.
    /// What lasts: its floor, the melt's cumulates under its pocket and the lava that
    /// floods it, up to the rim's crest and over it as an apron; and its rim, the crust
    /// the cumulates pushed out to its edge. What goes with its buoyancy: its own lift,
    /// and the drip of the lid that founders at its edge. The rim and the drip vary round
    /// the edge.
    pub floor: Radial,
    pub rim: Radial,
    pub lift: Radial,
    pub drip: Radial,
    /// The greatest tension its loads put in its lid's surface about its edge, Pa.
    pub tension: f64,
}

impl Head {
    /// The share of its buoyancy left at `t`, Gyr, for a corona formed at `t0`: all of it
    /// while its conduit feeds it, then falling as the lid conducts its heat away.
    pub fn buoyancy(&self, t0: f64, t: f64) -> f64 {
        if t < t0 {
            0.0
        } else {
            (-(t - t0 - self.feeding).max(0.0) / self.cooling.max(1e-6)).exp()
        }
    }

    /// The surface it gives `r` m from its centre at `t`, Gyr, for a corona formed at `t0`,
    /// where its rim's thickening and its drip stand at these shares of their means.
    pub fn profile_round(&self, r: f64, rim: f64, drip: f64, t0: f64, t: f64) -> f64 {
        if t < t0 {
            return 0.0;
        }
        self.floor.at(r)
            + rim * self.rim.at(r)
            + self.buoyancy(t0, t) * (self.lift.at(r) + drip * self.drip.at(r))
    }

    /// Its mean profile `r` m from its centre at `t`, Gyr.
    pub fn profile(&self, r: f64, t0: f64, t: f64) -> f64 {
        self.profile_round(r, 1.0, 1.0, t0, t)
    }

    /// How far out its profile reaches, m.
    pub fn reach(&self) -> f64 {
        self.floor
            .reach()
            .max(self.rim.reach())
            .max(self.lift.reach())
            .max(self.drip.reach())
    }
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
    (CORONA_FAULT_DIP_DEG - CORONA_LOCKUP_DIP_DEG)
        .to_radians()
        .sin()
        / (1.0 + GROOVE_DEPTH_MEAN)
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
            tilt: 1.0,
            head: None,
            wander: CORONA_OUTLINE_WANDER,
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
        let a = self.wander * self.outline_radius / body_radius;
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

    /// Let its outline wander by this share of its radius, its reach growing with it.
    pub fn set_wander(&mut self, wander: f64) {
        self.extent = (self.extent / (1.0 + 3.0 * self.wander) * (1.0 + 3.0 * wander))
            .min(std::f64::consts::PI);
        self.wander = wander;
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
    pub fn faults_on(&self, cell: f64) -> f64 {
        self.spacing.max(GROOVE_MIN_CELLS * cell)
    }

    /// How far from the upwelling the corona reaches on a grid of a cell, m: its outline
    /// and the damage zone beyond it, and a rock corona's profile.
    pub fn reach(&self, body_radius: f64, cell: f64) -> f64 {
        let fabric = self.extent * body_radius
            + (GROOVE_DAMAGE_SPACINGS * self.faults_on(cell)).max(0.5 * self.border(cell))
            + 2.0 * cell;
        match &self.head {
            Some(h) => fabric.max(h.reach() + 2.0 * cell),
            None => fabric,
        }
    }

    /// How much of its cover a corona renews at a point (F145): all of it over an ice
    /// mantle's upwelling; a rock corona's core, which its lava floods, while its belt
    /// keeps the old ground its faults break.
    pub fn renewing(&self, d: [f64; 3], body_radius: f64, noise: &CoronaNoise) -> f64 {
        if self.head.is_none() {
            return 1.0;
        }
        let (rho, _) = self.position(d, body_radius, noise);
        self.core_share(rho)
    }

    /// A rock corona's profile at a point at `t`, Gyr (F146): its head's, read at the
    /// distance that puts its pocket's edge on its outline, so that the rim follows the
    /// outline's wander; 0 for an icy corona.
    pub fn profile_at(&self, d: [f64; 3], body_radius: f64, noise: &CoronaNoise, t: f64) -> f64 {
        match &self.head {
            Some(h) => {
                let (rho, _) = self.position(d, body_radius, noise);
                let r = rho / CORONA_OUTLINE_SHARE * self.outline_radius;
                // Round the edge, the rim's thickening and the drip's arcs.
                let phi = azimuth_from(self.center, self.frame.0, self.frame.1, d);
                let at = |n: &GradientNoise, f: f64| {
                    n.fbm([phi.cos(), phi.sin(), 0.37], f, 2, 2.0, 0.5) / 0.5
                };
                let rim = (1.0
                    + CORONA_RIM_VARIATION * at(&noise.outline[0], 1.2).clamp(-1.0, 1.0))
                .max(0.0);
                let drip = 2.0 * smoothstep(0.5 + at(&noise.outline[1], 0.8));
                // The variation is the edge's: none toward the centre.
                let edge = smoothstep((r / self.outline_radius - 0.4) / 0.4);
                h.profile_round(
                    r,
                    1.0 + (rim - 1.0) * edge,
                    1.0 + (drip - 1.0) * edge,
                    self.time_gyr,
                    t,
                )
            }
            None => 0.0,
        }
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
    ) -> Vec<(usize, f64, f64)> {
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
                let d = grid.center(grid.cell(idx));
                let c = self.cover_with(d, body_radius, cell, &noise);
                (c > 0.0).then(|| {
                    let (rho, _) = self.position(d, body_radius, &noise);
                    (idx, c, self.core_share(rho))
                })
            })
            .collect()
    }

    /// How much of a cell at a share `rho` of the way out to the cell's edge is the
    /// corona's core rather than its belt (F143): the core within the belt's inner edge,
    /// blending across it.
    pub fn core_share(&self, rho: f64) -> f64 {
        let s = CORONA_OUTLINE_SHARE;
        let core_edge = s * (1.0 - CORONA_BELT_SHARE);
        1.0 - smoothstep((rho - core_edge) / (CORONA_CORE_BLEND * s) + 0.5)
    }

    /// What the corona makes of a cell (F133), as a lane's fabric: its cover; inside its
    /// outline the relief of its faults, an outer belt of segments concentric with the
    /// outline and a core crossed by two sets, with the footwalls of its outward-facing
    /// bounding faults up; beyond the outline their hanging walls down and the belt's
    /// faults running on as open fractures, fewer and shorter the farther out; the tone of
    /// its ice, the core's clean and darker in the troughs, the belt's the old crust's
    /// with clean ice on its scarps (F143), and the clean ice the fractures show.
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
        let amplitude = relief_share() * self.spacing * self.strength * self.tilt;
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
        let belt = 1.0 - self.core_share(rho);
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
        // The core is ice renewed from below, clean; the belt is the old crust broken into
        // tilt blocks, which keeps its tone and shows clean ice only on its faults' scarps,
        // as a fracture's walls do (F143).
        let core_tone = self.clean * set_tone * (1.0 - lag);
        let scarps = GROOVE_DAMAGE_CLEAN
            * self.clean
            * self.strength
            * base.gate
            * base.trough.clamp(0.0, 1.0);
        LaneFabric {
            cover,
            relief: grooves + up,
            outer: (1.0 - out) * (down - fracture_depth * open),
            tone: core_tone * (1.0 - belt) + scarps * belt,
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
    // Each corona forms at a time drawn from the tidal heat its lid carries above the floor
    // that renews an icy surface (F44): the lesser of the tides and the heat the lid
    // conducts, so that none falls while the ice is frozen through, or melting from below
    // under a lid still thick, or freezing from a warm start, or after it froze again
    // (F138). The heat's cumulative share against time.
    let mut cumulative = Vec::with_capacity(thermal.samples.len());
    let mut total = 0.0;
    for w in thermal.samples.windows(2) {
        let excess = |s: &crate::thermal::ThermalState| {
            (s.tidal_flux.min(s.heat_flux) - ICE_RESURFACING_FLUX_FLOOR).max(0.0)
        };
        total += 0.5 * (excess(&w[0]) + excess(&w[1])) * (w[1].time_gyr - w[0].time_gyr);
        cumulative.push((w[1].time_gyr, total));
    }
    let draw = |u: f64| -> f64 {
        if total <= 0.0 {
            return (mid + (u - 0.5) * width).clamp(0.0, age);
        }
        let target = u * total;
        let k = cumulative
            .partition_point(|p| p.1 < target)
            .min(cumulative.len() - 1);
        let (t1, c1) = cumulative[k];
        let (t0, c0) = if k > 0 {
            cumulative[k - 1]
        } else {
            (thermal.samples[0].time_gyr, 0.0)
        };
        if c1 > c0 {
            t0 + (target - c0) / (c1 - c0) * (t1 - t0)
        } else {
            t1
        }
    };
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
            let time = draw(rng.uniform(k, 10)).clamp(0.0, age);
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

/// A rock lid at a place and time (F145): the mantle under it, as hot as the mantle there
/// and the plumes' hot mantle about it make it, and the lid that mantle holds.
struct RockLid<'a> {
    mantle: &'a MantleHistory,
    rises: &'a [Rise],
    /// The upwelling limbs between the plumes, along which the diapirs are drawn as toward
    /// the plumes (F151).
    limbs: &'a Limbs,
    body_radius: f64,
    gravity: f64,
    /// The mantle's depth, m, the width of a large upwelling's convection cell.
    depth: f64,
}

/// The lid's state over a place at a time.
#[derive(Clone, Copy, Debug)]
struct LidState {
    /// The mantle's temperature under it, K, and the depth to the lid's stiff base and its
    /// elastic thickness, m.
    mantle: f64,
    stiff: f64,
    elastic: f64,
}

impl RockLid<'_> {
    /// How much hotter than the mantle at large the plumes' mantle makes it at `d` and
    /// `t`, K: a plume's excess under its rise, falling off over the rise's radius, from
    /// when its rise begins to grow.
    fn upwelling(&self, d: [f64; 3], t: f64) -> f64 {
        self.rises
            .iter()
            .map(|r| {
                let x = angle_between(d, r.center) * self.body_radius
                    / (PLUME_HEAD_SHARE * r.radius).max(1.0);
                PLUME_EXCESS_K * ramp(r, t) * (-x * x).exp()
            })
            .fold(0.0, f64::max)
    }

    fn state(&self, d: [f64; 3], t: f64) -> LidState {
        let m = self.mantle;
        let ts = m.surface_temperature;
        let t_m = m.temperature_at(t);
        let hot = t_m + self.upwelling(d, t);
        let lid = m.lid_at(t) * crate::mantle::hotter_lid_share(t_m, hot - t_m, ts, self.gravity);
        let rheological = GAS_CONSTANT * hot * hot / MANTLE_ACTIVATION;
        let span = (hot - ts).max(1.0);
        LidState {
            mantle: hot,
            stiff: lid * (hot - LID_BASE_SCALES * rheological - ts).max(0.0) / span,
            elastic: lid * ((ELASTIC_ISOTHERM_ROCK - ts) / span).clamp(0.0, 1.0),
        }
    }

    /// Whether a place lies where the rifts run at `t` (F183): on a limb between rises,
    /// where its pull stands past CORONA_FOCUS_PULL, along which the rifts join them
    /// (F150).
    fn on_rift(&self, d: [f64; 3], t: f64) -> bool {
        self.limbs.pull(self.rises, d, t, self.body_radius).0 >= CORONA_FOCUS_PULL
    }

    /// The weight a place gets as a diapir's birthplace at `t`: more within the convection
    /// cells of the large upwellings and along the limbs between them, toward which the
    /// diapirs are drawn (F151).
    fn draw_weight(&self, d: [f64; 3], t: f64) -> f64 {
        let limb = self.limbs.pull(self.rises, d, t, self.body_radius).0;
        1.0 + CORONA_UPWELLING_DRAW
            * (limb
                + self
                    .rises
                    .iter()
                    .map(|r| {
                        let x = angle_between(d, r.center) * self.body_radius / self.depth;
                        ramp(r, t) * (-x * x).exp()
                    })
                    .sum::<f64>())
    }
}

/// How far a rise's plume has come up at `t`, 0 to 1: from when its rise starts to grow
/// until it has.
fn ramp(r: &Rise, t: f64) -> f64 {
    ((t - r.start_gyr) / (r.end_gyr - r.start_gyr).max(1e-6)).clamp(0.0, 1.0)
}

/// The coronae a rock mantle's diapirs make in its lid since the lid last formed (F145).
/// Besides its plumes a mantle sends up transient diapirs off its hot boundary layer, as
/// many as the share of the heat it convects into its lid that they carry allows, each born
/// where the heat comes up and at a time drawn from it. A diapir's head flattens under the
/// lid and melts where it stands above the solidus at the base of the lid over it, which
/// a plume's hot mantle thins and heats; its pocket of melt is the corona: the melt it
/// intrudes and erupts loads the lid, which bends under that and under the head's own lift
/// while the head stays warm, and the corona stands where the bending breaks the lid in an
/// annulus.
pub fn rock_coronae<P: Reads>(
    graph: &P,
    mantle: &MantleHistory,
    rises: &[Rise],
    limbs: &Limbs,
) -> Vec<Corona> {
    let regime = graph.text("interior.tectonic_regime");
    if !matches!(regime, "stagnant_lid" | "episodic_lid") || mantle.times.is_empty() {
        return Vec::new();
    }
    let body_radius = graph.num("body.radius");
    let gravity = graph.num("body.gravity");
    let age = graph.num("body.age");
    let area = 4.0 * std::f64::consts::PI * body_radius * body_radius;
    // Since the lid last formed and the surface has kept a record.
    let start = mantle
        .last_overturn()
        .unwrap_or(0.0)
        .max(graph.num("cratering.record_start"))
        .min(age);
    if age - start <= 1e-6 {
        return Vec::new();
    }
    // The heat the mantle convects into its lid over the window, J, as a cumulative share.
    let steps = 400;
    let mut cumulative = Vec::with_capacity(steps + 1);
    cumulative.push((start, 0.0));
    let mut total = 0.0;
    for i in 0..steps {
        let t0 = start + (age - start) * i as f64 / steps as f64;
        let t1 = start + (age - start) * (i + 1) as f64 / steps as f64;
        let q = MANTLE_FLUX_FRACTION * mantle.heat_flux_at(0.5 * (t0 + t1));
        total += q * area * (t1 - t0) * GYR;
        cumulative.push((t1, total));
    }
    let budget = CORONA_HEAT_SHARE * total;
    let time_of = |u: f64| -> f64 {
        let target = u * total;
        let k = cumulative
            .partition_point(|p| p.1 < target)
            .clamp(1, cumulative.len() - 1);
        let (t0, c0) = cumulative[k - 1];
        let (t1, c1) = cumulative[k];
        if c1 > c0 {
            t0 + (target - c0) / (c1 - c0) * (t1 - t0)
        } else {
            t1
        }
    };
    let lid = RockLid {
        mantle,
        rises,
        limbs,
        body_radius,
        gravity,
        depth: crate::mantle::mantle_depth(
            body_radius,
            graph.num("body.mass"),
            graph.num("body.core_mass_fraction"),
        ),
    };
    let rng = Rng::new(graph.num("seed.tectonics") as u64 ^ 0xC020_0D1A);
    let (e_lo, e_hi) = CORONA_HEAD_EXCESS;
    // The diapirs in order, each taking its heat of the share, then what each makes.
    let mut spent = 0.0;
    let mut diapirs = Vec::new();
    let mut k = 0u64;
    while spent < budget && k < CORONA_MAX_DIAPIRS {
        let radius = CORONA_HEAD_RADIUS * (CORONA_HEAD_SPREAD * rng.normal(k, 1)).exp();
        let thickness = CORONA_HEAD_ASPECT * radius;
        let excess = PLUME_EXCESS_K * (e_lo + (e_hi - e_lo) * rng.uniform(k, 2));
        // The heat its conduit brings while it feeds the head.
        spent +=
            CORONA_BUOYANCY_FLUX * std::f64::consts::PI * radius * radius * MANTLE_HEAT_CAPACITY
                / MANTLE_EXPANSIVITY
                * CORONA_FEED_GYR
                * GYR;
        let time = time_of(rng.uniform(k, 3));
        // Born where the heat comes up: drawn against the greatest weight, the plumes'
        // centres.
        let top = rises
            .iter()
            .map(|r| lid.draw_weight(r.center, time))
            .fold(1.0, f64::max);
        let mut center = rng.unit_direction(k, 4);
        for attempt in 0..64u64 {
            if rng.uniform(k, 100 + attempt) * top <= lid.draw_weight(center, time) {
                break;
            }
            center = rng.unit_direction(k, 200 + attempt);
        }
        diapirs.push((k, center, time, radius, thickness, excess));
        k += 1;
    }
    let made: Vec<Result<Corona, NoCorona>> = diapirs
        .par_iter()
        .map(|&(id, center, time, radius, thickness, excess)| {
            rock_corona(id, center, time, radius, thickness, excess, &lid, &rng)
        })
        .collect();
    let mut out = Vec::new();
    let mut tally = [0usize; 3];
    for m in made {
        match m {
            Ok(c) => out.push(c),
            Err(why) => tally[why as usize] += 1,
        }
    }
    if std::env::var_os("EXOGEN_DEBUG").is_some() {
        eprintln!(
            "rock coronae: {k} diapirs, {} coronae; no pocket under regional melting {}, too cool {}, lid unbroken {}",
            out.len(),
            tally[0],
            tally[1],
            tally[2]
        );
    }
    out
}

/// Why a diapir made no corona.
#[derive(Clone, Copy, Debug)]
enum NoCorona {
    Regional,
    Cool,
    Unbroken,
}

/// The corona a diapir makes, if its head melts under the lid and its loads break the lid
/// in an annulus (F145, F146).
#[allow(clippy::too_many_arguments)]
fn rock_corona(
    id: u64,
    center: [f64; 3],
    time: f64,
    radius: f64,
    thickness: f64,
    excess: f64,
    lid: &RockLid,
    rng: &Rng,
) -> Result<Corona, NoCorona> {
    let gravity = lid.gravity;
    let body_radius = lid.body_radius;
    let state = lid.state(center, time);
    let rho_g = DENSITY_ROCK * gravity;
    let pressure = rho_g * state.stiff;
    let solidus = crate::mantle::solidus(pressure);
    // The head melts where it stands above the solidus at the lid's stiff base, out to where
    // its Gaussian excess falls to what the solidus asks of the mantle about it. Where the
    // mantle about it melts too, the melting is the region's, its floods' and its plains',
    // and no pocket stands out of it: Mars's, Mercury's and the Moon's mantles through most
    // of their histories, and Venus's for a while after its lid founders.
    // Where the region melts, a rift draws the head's melt together, its stretched lid the
    // way up (F183; two thirds of Venus's coronae lie along its chasmata, Stofan et al.
    // 1992): on a limb, where the rifts run, the head's pocket stands as far out
    // as it melts at least twice what the region's mantle does; elsewhere the region's melt
    // spreads with it.
    let need = solidus - state.mantle;
    let threshold = if need > 0.0 {
        need
    } else if lid.on_rift(center, time) {
        CORONA_FOCUS_EXCESS * -need
    } else {
        return Err(NoCorona::Regional);
    };
    if excess <= threshold {
        return Err(NoCorona::Cool);
    }
    let pocket = (radius * (excess / threshold).ln().sqrt()).min(1.5 * radius);
    // The melt column: the head's mantle rose along its adiabat from where it crossed the
    // solidus to the lid's base, melting the more the farther above it, so that its melt
    // fraction stands highest at its top and falls with depth through its thickness, or to
    // nothing where the column it melted over is shorter.
    let p_gpa = pressure / 1e9;
    let adiabat = MANTLE_EXPANSIVITY * state.mantle / (DENSITY_ROCK * MANTLE_HEAT_CAPACITY);
    let slope = ((SOLIDUS_C1 + 2.0 * SOLIDUS_C2 * p_gpa) / 1e9 - adiabat).max(1e-9);
    let head = rho_g * thickness;
    let column = |r: f64| -> f64 {
        let e = (state.mantle + excess * (-(r / radius).powi(2)).exp() - solidus).max(0.0);
        let dp = e / slope;
        let integral = if dp <= head {
            0.5 * dp * dp
        } else {
            head * (dp - 0.5 * head)
        };
        MELT_PRODUCTIVITY * integral / rho_g
    };
    // While its conduit feeds it the head heats the lid over it, which thins toward the
    // stagnant lid the head's mantle holds over the time heat takes to cross it: the lid
    // its loads bend (Smrekar et al. 2023: 11 ± 7 km at Venus's coronae).
    let crossing = state.stiff * state.stiff / (std::f64::consts::PI.powi(2) * MANTLE_DIFFUSIVITY);
    let settled = 1.0 - (-CORONA_FEED_GYR * GYR / crossing.max(1.0)).exp();
    let share = crate::mantle::hotter_lid_share(
        state.mantle,
        excess,
        lid.mantle.surface_temperature,
        gravity,
    );
    let elastic = state.elastic * (1.0 + (share - 1.0) * settled);
    let plate = Plate::rock(elastic.max(1e3), DENSITY_ROCK, gravity);
    let alpha = plate.alpha();
    let step = (pocket.min(alpha) / 12.0).max(500.0);
    let reach = (2.5 * radius).max(pocket + 8.0 * alpha);
    let n = (reach / step).ceil() as usize + 1;
    let rs: Vec<f64> = (0..n).map(|i| i as f64 * step).collect();
    let ring = |r: f64, centre: f64| (-((r - centre) / (CORONA_RING_SHARE * pocket)).powi(2)).exp();
    let inside = |r: f64| if r < pocket { 1.0 } else { 0.0 };
    // The head's own melt, past what the region's mantle makes about it.
    let region = column(f64::INFINITY);
    let intruded: Vec<f64> = rs
        .iter()
        .map(|&r| CORONA_INTRUSIVE_SHARE * (column(r) - region).max(0.0) * inside(r))
        .collect();
    // The crust the cumulates displace thickens the pocket's edge, as much as they take,
    // and the lava that erupts through the annulus builds on it.
    let shell = |i: usize| 2.0 * std::f64::consts::PI * rs[i].max(0.25 * step) * step;
    let volume: f64 = (0..n).map(|i| intruded[i] * shell(i)).sum();
    let ring_area: f64 = (0..n).map(|i| ring(rs[i], pocket) * shell(i)).sum();
    let thickened = volume / ring_area.max(1.0);
    let stack_load: Vec<f64> = intruded
        .iter()
        .map(|m| -CORONA_STACK_DENSITY * gravity * m)
        .collect();
    // The lava erupts through the annulus: a share piles at its vents into the rim, and the
    // rest runs down into the pocket, flooding its floor up to the rim's crest, what more
    // there is spilling over the rim as an apron.
    let erupts = (1.0 - CORONA_INTRUSIVE_SHARE) / CORONA_INTRUSIVE_SHARE * volume;
    let vents = CORONA_VENT_SHARE * erupts / ring_area.max(1.0);
    let rim_load: Vec<f64> = rs
        .iter()
        .map(|&r| {
            ((DENSITY_ROCK - DENSITY_CRUST_ROCK) * thickened - DENSITY_CRUST_ROCK * vents)
                * gravity
                * ring(r, pocket)
        })
        .collect();
    // The plate's bending under the rim, and the surface: the vents' lava on it.
    let rim_plate = deflect(&plate, step, &rim_load);
    let rim_w: Vec<f64> = rim_plate
        .iter()
        .zip(&rs)
        .map(|(w, &r)| w + vents * ring(r, pocket))
        .collect();
    let lava = (1.0 - CORONA_VENT_SHARE) * erupts;
    let ground: Vec<f64> = deflect(&plate, step, &stack_load)
        .iter()
        .zip(&rim_w)
        .map(|(a, b)| a + b)
        .collect();
    let crest = (0..n)
        .filter(|&i| (0.6 * pocket..1.4 * pocket).contains(&rs[i]))
        .max_by(|&a, &b| ground[a].total_cmp(&ground[b]))
        .unwrap_or(0);
    let room: f64 = (0..crest)
        .map(|i| (ground[crest] - ground[i]).max(0.0) * shell(i))
        .sum();
    let mut erupted = vec![0.0f64; n];
    if lava <= room {
        let (mut lo, mut hi) = (
            ground[..crest.max(1)]
                .iter()
                .cloned()
                .fold(f64::INFINITY, f64::min),
            ground[crest],
        );
        for _ in 0..50 {
            let level = 0.5 * (lo + hi);
            let held: f64 = (0..crest)
                .map(|i| (level - ground[i]).max(0.0) * shell(i))
                .sum();
            if held > lava {
                hi = level;
            } else {
                lo = level;
            }
        }
        for i in 0..crest {
            erupted[i] = (lo - ground[i]).max(0.0);
        }
    } else {
        for i in 0..crest {
            erupted[i] = (ground[crest] - ground[i]).max(0.0);
        }
        let apron = |r: f64| (-((r - rs[crest]) / (0.5 * pocket)).powi(2)).exp();
        let spread: f64 = (crest..n).map(|i| apron(rs[i]) * shell(i)).sum();
        for i in crest..n {
            erupted[i] = (lava - room) * apron(rs[i]) / spread.max(1.0);
        }
    }
    let floor_load: Vec<f64> = (0..n)
        .map(|i| stack_load[i] - DENSITY_CRUST_ROCK * gravity * erupted[i])
        .collect();
    let floor: Vec<f64> = deflect(&plate, step, &floor_load)
        .iter()
        .zip(&erupted)
        .map(|(w, t)| w + t)
        .collect();
    // While it stays warm the head lifts the lid, and the lid's lower part founders at the
    // pocket's edge, as cold as the lid's mean.
    let lift = DENSITY_ROCK * MANTLE_EXPANSIVITY * gravity * thickness * excess;
    let drip = DENSITY_ROCK
        * MANTLE_EXPANSIVITY
        * gravity
        * 0.5
        * (state.mantle - lid.mantle.surface_temperature)
        * (state.stiff - state.elastic).max(0.0);
    let lift_load: Vec<f64> = rs
        .iter()
        .map(|&r| lift * (-(r / radius).powi(2)).exp())
        .collect();
    let drip_load: Vec<f64> = rs
        .iter()
        .map(|&r| -drip * ring(r, pocket * (1.0 + CORONA_RING_SHARE)))
        .collect();
    let lift_w = deflect(&plate, step, &lift_load);
    let drip_w = deflect(&plate, step, &drip_load);
    let both: Vec<f64> = (0..n)
        .map(|i| floor[i] - erupted[i] + rim_plate[i] + lift_w[i] + drip_w[i])
        .collect();
    let tension = surface_tension(&plate, step, &both)
        .iter()
        .zip(&rs)
        .filter(|(_, r)| (0.5 * pocket..1.5 * pocket).contains(*r))
        .map(|(s, _)| *s)
        .fold(0.0, f64::max);
    // How fully the lid breaks: by the tension its loads put in its surface, and, as over
    // an ice mantle's upwelling (F133), by the share of the pocket's loads it takes in an
    // annulus of its own rather than spreading them as a swell over its flexural length.
    let strength = ((tension - CORONA_FAULT_STRENGTH) / CORONA_FAULT_STRENGTH).clamp(0.0, 1.0)
        * localized(alpha, pocket);
    if strength <= CORONA_STRENGTH_FLOOR {
        return Err(NoCorona::Unbroken);
    }
    // Trim the profiles where they die out.
    let last = (0..n)
        .rev()
        .find(|&i| {
            floor[i].abs() > 0.5
                || rim_w[i].abs() > 0.5
                || lift_w[i].abs() > 0.5
                || drip_w[i].abs() > 0.5
        })
        .map_or(1, |i| (i + 2).min(n));
    let radial = |v: &[f64]| Radial {
        step,
        values: v[..last].to_vec(),
    };
    // The cell the corona's fabric reads: a ring of neighbours twice its reach away, so
    // that its outline is its pocket's edge.
    let (east, north) = local_basis(center);
    let cell_reach = pocket / CORONA_OUTLINE_SHARE / body_radius;
    let neighbours: Vec<[f64; 3]> = (0..16)
        .map(|j| {
            offset(
                center,
                east,
                north,
                2.0 * cell_reach,
                std::f64::consts::TAU * j as f64 / 16.0,
            )
        })
        .collect();
    let mut corona = Corona::new(id, center, &neighbours, body_radius);
    corona.set_wander(CORONA_HEAD_WANDER);
    corona.time_gyr = time;
    corona.strength = strength;
    corona.spacing = GROOVE_NECKING_WAVELENGTH * GROOVE_RIFT_LID_SHARE * elastic.max(1e3);
    corona.clean = 0.0;
    corona.strike = std::f64::consts::TAU * rng.uniform(id, 12);
    // The strain the bending gave the surface, against the strain at which blocks that broke
    // at the faults' dip lock up.
    let lockup =
        CORONA_FAULT_DIP_DEG.to_radians().sin() / CORONA_LOCKUP_DIP_DEG.to_radians().sin() - 1.0;
    let strain = tension * (1.0 - ROCK_POISSON * ROCK_POISSON) / ROCK_YOUNG_MODULUS;
    corona.tilt = (strain / lockup).clamp(0.0, 1.0);
    let cooling = thickness * thickness / (4.0 * MANTLE_DIFFUSIVITY) / GYR;
    corona.head = Some(Head {
        radius,
        thickness,
        excess,
        feeding: CORONA_FEED_GYR,
        cooling,
        magma: column(0.0),
        elastic,
        floor: radial(&floor),
        rim: radial(&rim_w),
        lift: radial(&lift_w),
        drip: radial(&drip_w),
        tension,
    });
    Ok(corona)
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
                    let d = grid.center(grid.cell(k));
                    let f = c.fabric_at(d, body_radius, cell, &noise);
                    let cover = f.cover * c.renewing(d, body_radius, &noise);
                    (cover > 0.0 || f.fracture > 0.004).then_some(LaneHit {
                        k: k as u32,
                        cover: cover as f32,
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
/// ground as it stands. A cell something younger renewed is left alone. A rock corona
/// (F145, F146) renews only its core, which its lava floods, breaks the old ground of its
/// belt with its faults, and stands on the profile its loads give the lid at `present`.
#[allow(clippy::too_many_arguments)]
pub fn corona_compute(
    layers: &LayerSet,
    grid: &CubeSphere,
    body_radius: f64,
    corona: &Corona,
    cell: f64,
    present: f64,
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
            let d = grid.center(grid.cell(k));
            let f = corona.fabric_at(d, body_radius, cell, &noise);
            let profile = corona.profile_at(d, body_radius, &noise, present);
            if f.cover <= 0.0 && f.outer.abs() < 0.5 && profile.abs() < 0.5 {
                return None;
            }
            let h = bedrock[k] as f64;
            let new = if corona.head.is_some() {
                let renew = f.cover * corona.renewing(d, body_radius, &noise);
                h + (base[k] as f64 - h) * renew
                    + f.relief * f.cover
                    + f.outer * (1.0 - f.cover)
                    + profile
            } else {
                h + (base[k] as f64 + f.relief - h) * f.cover + f.outer * (1.0 - f.cover)
            };
            Some(CellWrite {
                idx: k as u32,
                bedrock: new as f32,
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
        let (lo, hi) = cells
            .iter()
            .fold((f32::INFINITY, f32::NEG_INFINITY), |(lo, hi), w| {
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

    /// Prints a rock body's coronae against Venus's (F145).
    #[test]
    #[ignore]
    fn print_rock_coronae() {
        for name in ["Venus", "Mars", "Mercury", "Moon"] {
            let (g, _) = crate::presets::preset(name).unwrap().to_graph();
            let mantle = crate::mantle::evolve(&crate::mantle::MantleInputs::from_graph(&g));
            let thermal = crate::thermal::ThermalHistory::from_graph(&g);
            let rises = crate::stages::plume::place_rises(&g, &thermal)
                .map(|(r, _)| r)
                .unwrap_or_default();
            let t0 = std::time::Instant::now();
            let coronae = rock_coronae(&g, &mantle, &rises, &Limbs::default());
            let r = g.f64("body.radius");
            let age = g.f64("body.age");
            println!(
                "{name}: {} coronae in {:.2}s, window from {:.3} Gyr, rises {}",
                coronae.len(),
                t0.elapsed().as_secs_f64(),
                mantle
                    .last_overturn()
                    .unwrap_or(0.0)
                    .max(g.f64("cratering.record_start")),
                rises.len()
            );
            if coronae.is_empty() {
                continue;
            }
            let mut d: Vec<f64> = coronae
                .iter()
                .map(|c| 2.0 * c.outline_radius() / 1e3)
                .collect();
            d.sort_by(f64::total_cmp);
            let gm = (d.iter().map(|x| x.ln()).sum::<f64>() / d.len() as f64).exp();
            let sl =
                (d.iter().map(|x| (x.ln() - gm.ln()).powi(2)).sum::<f64>() / d.len() as f64).sqrt();
            let pct = |q: f64| d[((d.len() - 1) as f64 * q) as usize];
            println!(
                "  diameters km: geo mean {gm:.0} sigma_ln {sl:.2}; 5% {:.0} 50% {:.0} 95% {:.0} max {:.0}",
                pct(0.05),
                pct(0.5),
                pct(0.95),
                d[d.len() - 1]
            );
            let near = coronae
                .iter()
                .filter(|c| {
                    rises
                        .iter()
                        .any(|ri| angle_between(c.center, ri.center) * r < 1.5 * ri.radius)
                })
                .count();
            let depth = crate::mantle::MantleInputs::from_graph(&g).depth();
            let cell = coronae
                .iter()
                .filter(|c| {
                    rises
                        .iter()
                        .any(|ri| angle_between(c.center, ri.center) * r < depth)
                })
                .count();
            println!("  within 1.5 rise radii of a rise: {near}; within a mantle depth: {cell}");
            let mut active = 0;
            let mut domes = 0;
            let mut rimmed = 0;
            for c in &coronae {
                let h = c.head.as_ref().unwrap();
                let b = h.buoyancy(c.time_gyr, age);
                if b > 0.3 {
                    active += 1;
                }
                let a = c.outline_radius();
                let centre = h.profile(0.0, c.time_gyr, age);
                let rim = (0..20)
                    .map(|i| h.profile(a * (0.8 + 0.02 * i as f64), c.time_gyr, age))
                    .fold(f64::NEG_INFINITY, f64::max);
                if centre > 0.5 * rim.max(0.0) && centre > 50.0 {
                    domes += 1;
                } else if rim > centre + 100.0 {
                    rimmed += 1;
                }
            }
            println!(
                "  active (buoyancy > 0.3) {active}; raised interiors {domes}; rimmed lows {rimmed}"
            );
            for c in coronae.iter().step_by((coronae.len() / 8).max(1)) {
                let h = c.head.as_ref().unwrap();
                let a = c.outline_radius();
                println!(
                    "  id {} t {:.3} D {:.0} km head {:.0} km dT {:.0} K Te {:.1} km magma {:.2} km tension {:.0} MPa strength {:.2} tilt {:.4} cool {:.0} Myr | centre {:.0} rim {:.0} moat {:.0} m",
                    c.id,
                    c.time_gyr,
                    2.0 * a / 1e3,
                    h.radius / 1e3,
                    h.excess,
                    h.elastic / 1e3,
                    h.magma / 1e3,
                    h.tension / 1e6,
                    c.strength,
                    c.tilt,
                    h.cooling * 1e3,
                    h.profile(0.0, c.time_gyr, age),
                    (0..20)
                        .map(|i| h.profile(a * (0.8 + 0.02 * i as f64), c.time_gyr, age))
                        .fold(f64::NEG_INFINITY, f64::max),
                    (0..40)
                        .map(|i| h.profile(a * (1.0 + 0.05 * i as f64), c.time_gyr, age))
                        .fold(f64::INFINITY, f64::min),
                );
            }
        }
    }

    /// F145: Venus's diapirs make some hundreds of coronae, a couple of hundred kilometres
    /// across on average, some of them still over a buoyant head, on a lid the heads have
    /// thinned to some ten kilometres; the thick lids of Mars and the Moon break under none.
    #[test]
    fn venus_s_diapirs_make_its_coronae_and_thick_lids_make_none() {
        let coronae = |name: &str| {
            let (g, _) = crate::presets::preset(name).unwrap().to_graph();
            let mantle = crate::mantle::evolve(&crate::mantle::MantleInputs::from_graph(&g));
            let thermal = crate::thermal::ThermalHistory::from_graph(&g);
            let rises = crate::stages::plume::place_rises(&g, &thermal)
                .map(|(r, _)| r)
                .unwrap_or_default();
            (
                rock_coronae(&g, &mantle, &rises, &Limbs::default()),
                g.f64("body.age"),
            )
        };
        let (venus, age) = coronae("Venus");
        assert!((400..1200).contains(&venus.len()), "{}", venus.len());
        let d: Vec<f64> = venus.iter().map(|c| 2.0 * c.outline_radius()).collect();
        let mean = (d.iter().map(|x| x.ln()).sum::<f64>() / d.len() as f64).exp();
        assert!((120e3..300e3).contains(&mean), "{mean}");
        let active = venus
            .iter()
            .filter(|c| c.head.as_ref().unwrap().buoyancy(c.time_gyr, age) > 0.3)
            .count() as f64
            / venus.len() as f64;
        assert!((0.1..0.7).contains(&active), "{active}");
        let te: Vec<f64> = venus
            .iter()
            .map(|c| c.head.as_ref().unwrap().elastic)
            .collect();
        let te_mean = te.iter().sum::<f64>() / te.len() as f64;
        assert!((4e3..20e3).contains(&te_mean), "{te_mean}");
        for name in ["Mars", "Moon"] {
            assert!(coronae(name).0.is_empty(), "{name}");
        }
    }

    /// F146: a corona's head lifts its pocket while it is fed and the lift fades after; the
    /// lava it erupted floods its floor no higher than its rim's crest.
    /// Venus's diapirs are drawn along the upwelling limbs between its plumes as toward the
    /// plumes (F151), so a good share of its coronae lie on the limbs, where its rifts run;
    /// their number stays Venus's.
    #[test]
    fn venus_s_coronae_gather_along_its_limbs() {
        let (g, _) = crate::presets::preset("Venus").unwrap().to_graph();
        let thermal = ThermalHistory::from_graph(&g);
        let r = g.f64("body.radius");
        let age = g.f64("body.age");
        let p = crate::stages::plume::PlumeProcess::new(&g, CubeSphere::new(256), &thermal, None)
            .unwrap();
        let mantle = crate::mantle::evolve(&crate::mantle::MantleInputs::from_graph(&g));
        let coronae = rock_coronae(&g, &mantle, &p.rises, &p.limbs);
        assert!(
            coronae.len() >= 400 && coronae.len() <= 1100,
            "{}",
            coronae.len()
        );
        let on_limb = coronae
            .iter()
            .filter(|c| p.limbs.pull(&p.rises, c.center, age, r).0 > 0.5)
            .count();
        assert!(
            on_limb as f64 >= 0.25 * coronae.len() as f64,
            "{on_limb} of {} on the limbs",
            coronae.len()
        );
    }

    #[test]
    fn a_rock_corona_stands_on_its_loads() {
        let (g, _) = crate::presets::preset("Venus").unwrap().to_graph();
        let mantle = crate::mantle::evolve(&crate::mantle::MantleInputs::from_graph(&g));
        let thermal = crate::thermal::ThermalHistory::from_graph(&g);
        let rises = crate::stages::plume::place_rises(&g, &thermal)
            .map(|(r, _)| r)
            .unwrap_or_default();
        let coronae = rock_coronae(&g, &mantle, &rises, &Limbs::default());
        let mut lifted = 0;
        for c in &coronae {
            let h = c.head.as_ref().unwrap();
            let t0 = c.time_gyr;
            assert!((h.buoyancy(t0, t0 + 0.5 * h.feeding) - 1.0).abs() < 1e-9);
            assert!(h.buoyancy(t0, t0 + h.feeding + 3.0 * h.cooling) < 0.06);
            if h.lift.at(0.0) > 0.0 {
                lifted += 1;
            }
            // The floor stands no higher than the mean rim's crest, the lava filling it.
            let a = c.outline_radius();
            let crest = (0..=40)
                .map(|i| {
                    h.floor.at(a * (0.6 + 0.02 * i as f64)) + h.rim.at(a * (0.6 + 0.02 * i as f64))
                })
                .fold(f64::NEG_INFINITY, f64::max);
            let centre = h.floor.at(0.0) + h.rim.at(0.0);
            assert!(centre <= crest + 1.0, "{} {centre} {crest}", c.id);
        }
        assert_eq!(lifted, coronae.len());
    }

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
        for name in [
            "Ganymede",
            "Europa",
            "Enceladus",
            "Callisto",
            "Ceres",
            "Pluto",
        ] {
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
