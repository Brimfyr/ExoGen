// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! The lid's strain from the body's change of radius, and the faults it opens (F147).
//!
//! A body whose radius changes strains its lid the same way in every direction and through
//! its whole thickness. An icy body grows as its ocean freezes, water expanding by nine
//! percent as it turns to ice (Charon's lid stretched by some 1 percent of its area, Beyer
//! et al. 2017), and shrinks where its ocean lies deep enough to freeze as a denser ice;
//! every body shrinks a little as its interior cools (F16). The lid
//! bends elastically until its stress breaks its faulted layer, in tension at under a third
//! of the stress compression needs, and its faults take up whatever strain lies beyond, so
//! that a shell melted and refrozen to where it was breaks nothing new. The stretch since
//! its surface began to keep a record opens normal faults that cut the lid's brittle layer:
//! graben as wide as their faults, dipping at 60°, stand apart where they meet at that
//! layer's base, kilometres deep for their length in ice, which yields some ten times as
//! far as rock for the same stress. With no direction to follow, each runs its own way
//! until it meets an older one at least as long, which takes up the stress it would have
//! run on (a longer one cuts across it), so that the lid breaks into blocks, and half of
//! them run beside one drawn before, so that they gather in belts.

use std::cell::RefCell;
use std::collections::HashMap;

use crate::constants::*;
use crate::geom::{angle_between, local_basis};
use crate::grid::CubeSphere;
use crate::params::ParamGraph;
use crate::rng::Rng;
use crate::stages::plume::{FaultLayer, FaultRules, Trough, populate, trace_in};
use crate::thermal::ThermalHistory;

/// What a lid is made of: its density, kg/m³, Young's modulus, Pa, Poisson's ratio and its
/// faults' friction.
#[derive(Clone, Copy, Debug)]
pub struct LidMaterial {
    pub density: f64,
    pub modulus: f64,
    pub poisson: f64,
    pub friction: f64,
}

impl LidMaterial {
    /// Ice over an icy body's mantle, rock over a rocky one's.
    pub fn of(graph: &ParamGraph) -> Self {
        if graph.f64("body.ice_layer_thickness") > 0.0 {
            LidMaterial {
                density: DENSITY_ICE,
                modulus: ICE_YOUNG_MODULUS,
                poisson: ICE_POISSON,
                friction: ICE_FRICTION,
            }
        } else {
            LidMaterial {
                density: DENSITY_CRUST_ROCK,
                modulus: ROCK_YOUNG_MODULUS,
                poisson: ROCK_POISSON,
                friction: ROCK_FRICTION,
            }
        }
    }

    /// The stress differences that break a layer down to where it bears this pressure, Pa,
    /// in tension and in compression: Anderson's normal and thrust faults on planes at the
    /// friction's best angle, `K = (sqrt(1 + mu²) + mu)²` (Sibson 1974).
    pub fn breaking(&self, pressure: f64) -> (f64, f64) {
        let k = ((1.0 + self.friction.powi(2)).sqrt() + self.friction).powi(2);
        ((1.0 - 1.0 / k) * pressure, (k - 1.0) * pressure)
    }

    /// The stiffness of a membrane of it stretched alike every way, Pa.
    pub fn biaxial(&self) -> f64 {
        self.modulus / (1.0 - self.poisson)
    }
}

/// The depth at which the lid's faults close, m: where the overburden reaches the pressure
/// at which its stress no longer breaks it (F147).
pub fn closure_depth(graph: &ParamGraph) -> f64 {
    let material = LidMaterial::of(graph);
    FAULT_LAYER_PRESSURE / (material.density * graph.f64("body.gravity").max(1e-3))
}

/// The deepest the lid's faults reach at a time, m: within its elastic part, above where it
/// creeps rather than breaks (a shell refreezing over a remelted ocean is thin), within its
/// ice, and within a thin shell's share of the radius.
pub fn brittle_bound(graph: &ParamGraph, thermal: &ThermalHistory, t: f64) -> f64 {
    thermal.at(t).elastic.min(shell_bound(graph)).max(1e3)
}

/// The deepest any of the lid's faults reach, m: within its ice and a thin shell's share of
/// the radius.
pub fn shell_bound(graph: &ParamGraph) -> f64 {
    let ice = graph.f64("body.ice_layer_thickness");
    let within = if ice > 0.0 { ice } else { f64::INFINITY };
    within.min(LID_LAYER_RADIUS_SHARE * graph.f64("body.radius"))
}

/// The lid's strain from the body's change of radius (F147): the extension and the
/// shortening its faults take up, as shares of its area, cumulative from the record's start
/// at times through the record. The record starts as the craters' does: when the surface
/// began to keep one, or as long before the present as a renewed surface keeps its craters
/// (Europa's 100 Myr), the faults before then gone with them.
#[derive(Clone, Debug, Default)]
pub struct LidStrain {
    pub times: Vec<f64>,
    pub extension: Vec<f64>,
    pub shortening: Vec<f64>,
    /// The body's change in radius since it formed at those times, m.
    pub radius_change: Vec<f64>,
}

impl LidStrain {
    pub fn from_graph(graph: &ParamGraph, thermal: &ThermalHistory) -> Self {
        let r = graph.f64("body.radius");
        let ice = graph.f64("body.ice_layer_thickness");
        let age = graph.f64("body.age");
        let start = graph
            .f64("cratering.record_start")
            .max(age - graph.f64("cratering.retention_age"))
            .clamp(0.0, age);
        if age <= 1e-6 {
            return Self::default();
        }
        let material = LidMaterial::of(graph);
        let gravity = graph.f64("body.gravity").max(1e-3);
        let deepest = closure_depth(graph);
        // The ice frozen at a time: the solid shell over the ocean, the whole layer once the
        // ocean is gone; as ice I down to where the pressure reaches its limit, and as the
        // denser high-pressure ices below.
        let shell = |h: f64| 4.0 / 3.0 * std::f64::consts::PI * (r.powi(3) - (r - h).powi(3));
        let ice_one = ICE_I_PRESSURE_LIMIT / (DENSITY_ICE * gravity);
        let frozen = |t: f64| {
            if ice <= 0.0 {
                return (0.0, 0.0);
            }
            let h = thermal.at(t).lithosphere.clamp(0.0, ice);
            let upper = shell(h.min(ice_one));
            (upper, shell(h) - upper)
        };
        let (v0, w0) = frozen(0.0);
        let area = 4.0 * std::f64::consts::PI * r * r;
        let stiffness = material.biaxial();
        let mut out = Self::default();
        // The strain at which the lid would be unstressed, which its faults move.
        let mut relaxed = 0.0;
        let (mut extension, mut shortening) = (0.0, 0.0);
        let steps = 900;
        for i in 0..=steps {
            let t = age * i as f64 / steps as f64;
            let (v, w) = frozen(t);
            let change = (ICE_FREEZING_EXPANSION * (v - v0) + ICE_HP_FREEZING_CHANGE * (w - w0))
                / area
                - thermal.at(t).contraction;
            let strain = change / r;
            let layer = deepest.min(brittle_bound(graph, thermal, t));
            let (tension, compression) = material.breaking(material.density * gravity * layer);
            let recorded = t >= start;
            let stretch = strain - relaxed - tension / stiffness;
            if stretch > 0.0 {
                relaxed += stretch;
                if recorded {
                    extension += 2.0 * stretch;
                }
            }
            let squeeze = relaxed - strain - compression / stiffness;
            if squeeze > 0.0 {
                relaxed -= squeeze;
                if recorded {
                    shortening += 2.0 * squeeze;
                }
            }
            if recorded {
                out.times.push(t);
                out.extension.push(extension);
                out.shortening.push(shortening);
                out.radius_change.push(change);
            }
        }
        out
    }

    /// The time by which the lid had taken up this share of its extension since the record
    /// began, Gyr.
    pub fn extension_time(&self, share: f64) -> f64 {
        let total = self.extension.last().copied().unwrap_or(0.0);
        let n = self.times.len();
        if n < 2 || total <= 0.0 {
            return self.times.first().copied().unwrap_or(0.0);
        }
        let target = share.clamp(0.0, 1.0) * total;
        let k = self
            .extension
            .partition_point(|&e| e < target)
            .clamp(1, n - 1);
        let (e0, e1) = (self.extension[k - 1], self.extension[k]);
        let (t0, t1) = (self.times[k - 1], self.times[k]);
        if e1 > e0 {
            t0 + (target - e0) / (e1 - e0) * (t1 - t0)
        } else {
            t1
        }
    }
}

/// The faults the lid's growth opens since the surface began to keep a record (F147).
pub fn lid_faults(graph: &ParamGraph, grid: CubeSphere, thermal: &ThermalHistory) -> Vec<Trough> {
    if !crate::stages::ice::is_ice_surface(graph) {
        return Vec::new();
    }
    let strain = LidStrain::from_graph(graph, thermal);
    let total = strain.extension.last().copied().unwrap_or(0.0);
    if total <= 0.0 {
        return Vec::new();
    }
    let r = graph.f64("body.radius");
    let area = 4.0 * std::f64::consts::PI * r * r;
    let material = LidMaterial::of(graph);
    // The brittle layer the faults cut: down to where its faults close, within the lid's
    // brittle part at the fault's time.
    let deepest = closure_depth(graph);
    let bound = |t: f64| brittle_bound(graph, thermal, t);
    let layer = deepest.min(shell_bound(graph)).max(1e3);
    let width = 2.0 * layer / GRABEN_DIP_DEG.to_radians().tan();
    // Its faults slip per length as density over stiffness against rock's (a stated choice
    // between stiffness alone, which would cut Charon's chasmata a tenth of their length
    // deep, and the frictional stress at a depth, which would leave them shallower than
    // Mars's graben for their length; Charon's 0.016 of their length at most, F147).
    let compliance =
        (material.density / material.modulus) / (DENSITY_CRUST_ROCK / ROCK_YOUNG_MODULUS);
    let cell = grid.delta() * r;
    let l_min = (LID_FAULT_MIN_WIDTHS * width).max(4.0 * cell);
    let rules = FaultRules {
        l_min,
        l_max: (LID_FAULT_MAX_LENGTH * r).max(2.0 * l_min),
        depth_per_length: RIFT_DEPTH_PER_LENGTH * compliance,
        // A graben throws no more than its width allows (F111).
        depth_cap: GRABEN_DEPTH_PER_WIDTH * width,
        compliance,
        layer: FaultLayer::Lid(deepest),
        cell,
        body_radius: r,
    };
    // Where the faults drawn so far run, by the cells of a grid about a graben wide: each
    // point with the fault's half-width and length.
    let face = ((std::f64::consts::TAU * r / (4.0 * width)).round() as u32).clamp(8, 256);
    let coarse = CubeSphere::new(face);
    type Drawn = Vec<([f64; 3], f64, f64)>;
    let index: RefCell<HashMap<usize, Drawn>> = RefCell::new(HashMap::new());
    let meets = |p: [f64; 3], half: f64, length: f64| -> bool {
        let index = index.borrow();
        let (c, _, _) = coarse.locate(p);
        for di in -1..=1 {
            for dj in -1..=1 {
                let k = coarse.index(coarse.neighbor(c, di, dj));
                if let Some(points) = index.get(&k)
                    && points
                        .iter()
                        .any(|&(q, h, l)| l >= length && angle_between(p, q) * r < half + h)
                {
                    return true;
                }
            }
        }
        false
    };
    let rng = Rng::new(graph.f64("seed.tectonics") as u64 ^ 0x5714_1D00);
    // With the lid stretched alike every way, a fault keeps its own heading.
    let isotropic = |p: [f64; 3]| {
        let (east, north) = local_basis(p);
        (1.0, 1.0, east, east, north)
    };
    let fresh = |draw: u64, length: f64, spacing: f64, width: f64| {
        let start = rng.unit_direction(draw, 20);
        let stop = |p: [f64; 3]| meets(p, 0.5 * width, length);
        trace_in(
            &isotropic,
            start,
            length,
            spacing,
            r,
            &rng,
            draw,
            width,
            None,
            0.0,
            Some(&stop),
        )
    };
    // A fault breaks when the lid has stretched by as much as the faults before it took up,
    // so that the older ones are drawn first and stop the younger.
    let timing = |_u: f64, share: f64| strain.extension_time(share);
    let window = |t: f64, _carried: bool| (t, t);
    let drawn = |t: &Trough| {
        let mut index = index.borrow_mut();
        for &p in &t.points {
            let k = coarse.index(coarse.locate(p).0);
            index
                .entry(k)
                .or_default()
                .push((p, 0.5 * t.width, t.length));
        }
    };
    populate(
        &rules,
        total * area,
        &rng,
        (0x5714_0000_0000, 50_000_000_000),
        &fresh,
        &timing,
        &window,
        &bound,
        Some(&drawn),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::stages::plume::TroughKind;

    /// Charon's lid stretches as its ocean freezes, by at least 1 percent of its area
    /// (Beyer et al. 2017), and breaks into chasmata some 30 to 50 km wide and kilometres
    /// deep, Mandjet's 450 km among them, before its southern plain's flooding.
    #[test]
    fn charon_s_lid_stretches_as_its_ocean_freezes_and_breaks_into_chasmata() {
        let (g, _) = crate::presets::preset("Charon").unwrap().to_graph();
        let thermal = ThermalHistory::from_graph(&g);
        let strain = LidStrain::from_graph(&g, &thermal);
        let extension = strain.extension.last().copied().unwrap_or(0.0);
        assert!(extension > 0.01 && extension < 0.03, "{extension}");
        assert_eq!(strain.shortening.last().copied().unwrap_or(0.0), 0.0);
        let faults = lid_faults(&g, CubeSphere::new(128), &thermal);
        assert!(
            faults.len() >= 20 && faults.len() <= 150,
            "{}",
            faults.len()
        );
        let median = |mut v: Vec<f64>| {
            v.sort_by(f64::total_cmp);
            v[v.len() / 2]
        };
        let width = median(faults.iter().map(|f| f.width).collect());
        let depth = median(faults.iter().map(|f| f.depth).collect());
        let longest = faults.iter().map(|f| f.length).fold(0.0, f64::max);
        assert!(width > 30e3 && width < 60e3, "{width}");
        assert!(depth > 2e3 && depth < 9e3, "{depth}");
        assert!(longest > 400e3, "{longest}");
        assert!(
            faults
                .iter()
                .all(|f| f.kind == TroughKind::Graben && f.start_gyr < 1.0),
            "a rift or a late fault"
        );
        // The older faults are drawn first.
        assert!(
            faults
                .windows(2)
                .all(|w| w[0].start_gyr <= w[1].start_gyr + 1e-12)
        );
    }

    /// A shell melted a little and refrozen to where it was breaks nothing; one melted
    /// beyond what its lid bears in compression takes up the shortening, and refrozen
    /// stretches back by as much less what it bears in tension.
    #[test]
    fn a_shell_melted_and_refrozen_to_where_it_was_breaks_nothing_new() {
        let (g, _) = crate::presets::preset("Charon").unwrap().to_graph();
        let base = ThermalHistory::from_graph(&g);
        let melted = |dip: f64| {
            let mut h = base.clone();
            for s in &mut h.samples {
                let x = ((s.time_gyr - 1.0) / 0.5).clamp(0.0, 2.0);
                let bump = (std::f64::consts::PI * 0.5 * x).sin().max(0.0);
                s.lithosphere = 100e3 - dip * bump;
                s.elastic = 60e3;
                s.contraction = 0.0;
            }
            LidStrain::from_graph(&g, &h)
        };
        let small = melted(1e3);
        assert_eq!(small.extension.last().copied().unwrap(), 0.0);
        assert_eq!(small.shortening.last().copied().unwrap(), 0.0);
        let large = melted(60e3);
        let (e, sh) = (
            large.extension.last().copied().unwrap(),
            large.shortening.last().copied().unwrap(),
        );
        let material = LidMaterial::of(&g);
        let layer = closure_depth(&g).min(60e3);
        let (tension, _) = material.breaking(material.density * g.f64("body.gravity") * layer);
        let expect = sh - 2.0 * tension / material.biaxial();
        assert!(
            sh > 0.0 && (e - expect).abs() < 0.02 * expect,
            "{e} {sh} {expect}"
        );
    }

    /// Rock lids and Io's open no faults by this process.
    #[test]
    fn rock_lids_open_no_lid_faults() {
        for name in ["Io", "Mars", "Mercury"] {
            let (g, _) = crate::presets::preset(name).unwrap().to_graph();
            let thermal = ThermalHistory::from_graph(&g);
            assert!(
                lid_faults(&g, CubeSphere::new(64), &thermal).is_empty(),
                "{name}"
            );
        }
    }

    /// Prints a body's lid through its last 150 Myr (EXOGEN_BODY, Mimas by default).
    #[test]
    #[ignore = "diagnostic: a body's lid through its last 150 Myr"]
    fn print_late_lid() {
        let name = std::env::var("EXOGEN_BODY").unwrap_or_else(|_| "Mimas".into());
        let (g, _) = crate::presets::preset(&name).unwrap().to_graph();
        let thermal = ThermalHistory::from_graph(&g);
        let strain = LidStrain::from_graph(&g, &thermal);
        let age = g.f64("body.age");
        for k in 0..=15 {
            let t = age - 0.15 + 0.01 * k as f64;
            let s = thermal.at(t);
            let i = strain
                .times
                .partition_point(|&x| x < t)
                .min(strain.times.len() - 1);
            println!(
                "{t:.3} Gyr: lid {:.1} km elastic {:.1} km; radius {:+.3} km; extension {:.4}% shortening {:.4}%",
                s.lithosphere / 1e3,
                s.elastic / 1e3,
                strain.radius_change[i] / 1e3,
                100.0 * strain.extension[i],
                100.0 * strain.shortening[i]
            );
        }
    }

    /// Prints the lid faults whose traces jump between points (EXOGEN_BODY, Charon by default).
    #[test]
    #[ignore = "diagnostic: lid faults whose traces jump between points"]
    fn print_trace_gaps() {
        let name = std::env::var("EXOGEN_BODY").unwrap_or_else(|_| "Charon".into());
        let (g, _) = crate::presets::preset(&name).unwrap().to_graph();
        let thermal = ThermalHistory::from_graph(&g);
        let r = g.f64("body.radius");
        let faults = lid_faults(&g, CubeSphere::new(256), &thermal);
        let mut jumping = 0;
        for (i, f) in faults.iter().enumerate() {
            let gaps: Vec<f64> = f
                .points
                .windows(2)
                .map(|w| angle_between(w[0], w[1]) * r / f.spacing)
                .collect();
            let worst = gaps.iter().copied().fold(0.0, f64::max);
            let least = gaps.iter().copied().fold(f64::INFINITY, f64::min);
            let over = gaps.iter().filter(|&&x| x > 1.5).count();
            if worst > 1.5 || least < 0.5 {
                jumping += 1;
                println!(
                    "fault {i}: {} points, step {:.2} to {:.2} spacings, {over} jumps over 1.5",
                    f.points.len(),
                    least,
                    worst
                );
            }
        }
        println!("{name}: {} faults, {jumping} with jumps", faults.len());
    }

    /// Prints each icy body's strain since its record began and the faults it opens.
    #[test]
    #[ignore]
    fn print_lid_faults() {
        for name in crate::presets::PRESET_NAMES {
            let (g, _) = crate::presets::preset(name).unwrap().to_graph();
            if g.f64("body.ice_layer_thickness") <= 0.0 {
                continue;
            }
            let thermal = ThermalHistory::from_graph(&g);
            let strain = LidStrain::from_graph(&g, &thermal);
            let grid = CubeSphere::new(128);
            let t0 = std::time::Instant::now();
            let faults = lid_faults(&g, grid, &thermal);
            let mut l: Vec<f64> = faults.iter().map(|f| f.length / 1e3).collect();
            l.sort_by(f64::total_cmp);
            let deep = faults.iter().map(|f| f.depth).fold(0.0, f64::max);
            let wide = faults.iter().map(|f| f.width).fold(0.0, f64::max);
            let mut widths: Vec<f64> = faults.iter().map(|f| f.width / 1e3).collect();
            widths.sort_by(f64::total_cmp);
            let mut times: Vec<f64> = faults.iter().map(|f| f.start_gyr).collect();
            times.sort_by(f64::total_cmp);
            let mut depths: Vec<f64> = faults.iter().map(|f| f.depth / 1e3).collect();
            depths.sort_by(f64::total_cmp);
            let q = |v: &Vec<f64>, f: f64| {
                v.get(((v.len() as f64 * f) as usize).min(v.len().saturating_sub(1)))
                    .copied()
                    .unwrap_or(0.0)
            };
            let first = strain.radius_change.first().copied().unwrap_or(0.0);
            let last = strain.radius_change.last().copied().unwrap_or(0.0);
            println!(
                "{name}: radius {:+.2} km through the record; extension {:.3}% shortening {:.3}% of the area; {} faults in {:.2}s; length median {:.0} max {:.0} km; width median {:.0} max {:.0} km; depth median {:.1} deepest {:.1} km; times q10 {:.2} q50 {:.2} q90 {:.2}",
                (last - first) / 1e3,
                100.0 * strain.extension.last().copied().unwrap_or(0.0),
                100.0 * strain.shortening.last().copied().unwrap_or(0.0),
                faults.len(),
                t0.elapsed().as_secs_f64(),
                q(&l, 0.5),
                l.last().copied().unwrap_or(0.0),
                q(&widths, 0.5),
                wide / 1e3,
                q(&depths, 0.5),
                deep / 1e3,
                q(&times, 0.1),
                q(&times, 0.5),
                q(&times, 0.9),
            );
        }
    }
}
