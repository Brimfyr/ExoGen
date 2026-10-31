// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! The lid's strain from the body's change of radius, and the faults it opens (F147).
//!
//! A body whose radius changes strains its lid the same way in every direction and through
//! its whole thickness. An icy body grows as its ocean freezes, water expanding by nine
//! percent as it turns to ice (Charon's lid stretched by some 1 percent of its area, Beyer
//! et al. 2017), and shrinks where its water freezes as a denser ice: below ice I's limit
//! under a thick shell, and on the floor of an ocean deep enough that its water passes the
//! high-pressure ices' melting point (F232), so that Ganymede's, Callisto's and Titan's
//! lids hardly stretch as they cool; every body shrinks a little as its interior cools
//! (F16). The lid bends elastically until its stress breaks its faulted layer, in tension
//! at under a third of the stress compression needs, and its faults take up whatever strain
//! lies beyond, so that a shell melted and refrozen to where it was breaks nothing new. The
//! stretch since its surface began to keep a record opens normal faults in the lid's
//! brittle layer, each as deep as its length allows until it cuts the layer through, its
//! graben as wide as its faults, dipping at 60°, meet at their base, kilometres deep for
//! their length in ice, which yields some ten times as far as rock for the same stress.
//! They grow as the lid stretches (F233): the faults standing take up the stretch within
//! their shadows, lengthening, and link where they meet, before new ones open outside
//! them; with no direction to follow, each runs its own way until it meets one at least
//! half as wide, and half the new ones run beside one drawn before, so that they gather in
//! belts.

use crate::constants::*;
use crate::grid::CubeSphere;
use crate::params::ParamGraph;
use crate::stages::plume::{FaultLayer, FaultRules, Trough};
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
        Self::icy(graph.f64("body.ice_layer_thickness") > 0.0)
    }

    /// Ice or rock.
    pub fn icy(icy: bool) -> Self {
        if icy {
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
    within(
        graph.f64("body.ice_layer_thickness"),
        graph.f64("body.radius"),
    )
}

/// Within an ice layer this thick, m, and a thin shell's share of a radius, m.
fn within(ice: f64, radius: f64) -> f64 {
    let layer = if ice > 0.0 { ice } else { f64::INFINITY };
    layer.min(LID_LAYER_RADIUS_SHARE * radius)
}

/// The melting point of ice Ih at a pressure, K (IAPWS R14-08), by bisection on its melting
/// curve; that of the triple point with ice III beyond it.
fn ice_one_melting(p: f64) -> f64 {
    let (tt, pt) = ICE_IH_TRIPLE;
    let curve = |t: f64| {
        let theta = t / tt;
        pt * (1.0
            + ICE_IH_MELTING
                .iter()
                .map(|&(a, b)| a * (1.0 - theta.powf(b)))
                .sum::<f64>())
    };
    let (mut lo, mut hi) = (ICE_HP_MELTING[0].0, tt);
    if p >= curve(lo) {
        return lo;
    }
    for _ in 0..60 {
        let mid = 0.5 * (lo + hi);
        if curve(mid) > p {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

/// The pressure at which water at a temperature freezes as a high-pressure ice, Pa: on the
/// melting curve of ice III, V or VI that spans the temperature (IAPWS R14-08), the triple
/// point of ice Ih, ice III and the melt below them.
fn high_pressure_melting(t: f64) -> f64 {
    let (t_low, p_low, ..) = ICE_HP_MELTING[0];
    if t <= t_low {
        return p_low;
    }
    let &(t_star, p_star, a, b, _) = ICE_HP_MELTING
        .iter()
        .find(|&&(_, _, _, _, top)| t <= top)
        .unwrap_or(&ICE_HP_MELTING[2]);
    p_star * (1.0 - a * (1.0 - (t / t_star).powf(b)))
}

/// F232: the depth of an ocean's floor under a shell `h` thick, m: an ocean stands at the
/// melting point of the ice at its shell's base, and below it its water freezes as the
/// high-pressure ices where its pressure passes their melting point at that temperature
/// (Squyres 1980; Mueller and McKinnon 1988: the oceans of Ganymede, Callisto and Titan lie
/// between ice I and the high-pressure ices, which grow from below as they cool). The
/// shell's base reaching ice I's limit leaves no ocean.
fn ocean_floor(h: f64, gravity: f64) -> f64 {
    let base = DENSITY_ICE * gravity * h;
    let melt = high_pressure_melting(ice_one_melting(base));
    h + (melt - base).max(0.0) / (DENSITY_OCEAN * gravity)
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

/// What the lid's strain reads of the body: its radius, m, its ice layer, m, its age and
/// the start of its record, Gyr, and its gravity, m/s², from the parameter graph or, in a
/// derivation, from the nodes it declares.
#[derive(Clone, Copy, Debug)]
pub struct LidInputs {
    pub radius: f64,
    pub ice: f64,
    pub age: f64,
    pub start: f64,
    pub gravity: f64,
}

impl LidInputs {
    pub fn from_graph(graph: &ParamGraph) -> Self {
        let age = graph.f64("body.age");
        Self {
            radius: graph.f64("body.radius"),
            ice: graph.f64("body.ice_layer_thickness"),
            age,
            start: graph
                .f64("cratering.record_start")
                .max(age - graph.f64("cratering.retention_age"))
                .clamp(0.0, age),
            gravity: graph.f64("body.gravity"),
        }
    }
}

impl LidStrain {
    pub fn from_graph(graph: &ParamGraph, thermal: &ThermalHistory) -> Self {
        Self::from_inputs(&LidInputs::from_graph(graph), thermal)
    }

    pub fn from_inputs(inputs: &LidInputs, thermal: &ThermalHistory) -> Self {
        let r = inputs.radius;
        let ice = inputs.ice;
        let age = inputs.age;
        let start = inputs.start;
        if age <= 1e-6 {
            return Self::default();
        }
        let material = LidMaterial::icy(ice > 0.0);
        let gravity = inputs.gravity.max(1e-3);
        let deepest = FAULT_LAYER_PRESSURE / (material.density * gravity);
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
            // The ocean freezes on its floor as well, as the high-pressure ices, below where
            // its water passes their melting point at the ocean's temperature (F232).
            let floor = ocean_floor(h, gravity).clamp(h, ice);
            (upper, shell(h) - upper + shell(ice) - shell(floor))
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
            let layer = deepest.min(thermal.at(t).elastic.min(within(ice, r)).max(1e3));
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
    let layer = deepest.min(shell_bound(graph)).max(1e3);
    // The widest graben, whose faults cut the whole layer.
    let widest = 2.0 * layer / GRABEN_DIP_DEG.to_radians().tan();
    // Its faults slip per length as density over stiffness against rock's (a stated choice
    // between stiffness alone, which would cut Charon's chasmata a tenth of their length
    // deep, and the frictional stress at a depth, which would leave them shallower than
    // Mars's graben for their length; Charon's 0.016 of their length at most, F147).
    let compliance =
        (material.density / material.modulus) / (DENSITY_CRUST_ROCK / ROCK_YOUNG_MODULUS);
    let cell = grid.delta() * r;
    let l_min = LID_FAULT_MIN_CELLS * cell;
    let rules = FaultRules {
        l_min,
        l_max: (LID_FAULT_MAX_LENGTH * r).max(2.0 * l_min),
        depth_per_length: RIFT_DEPTH_PER_LENGTH * compliance,
        // A graben throws no more than its width allows (F111).
        depth_cap: GRABEN_DEPTH_PER_WIDTH * widest,
        compliance,
        layer: FaultLayer::Lid(layer),
        cell,
        body_radius: r,
    };
    // The faults' shadows read on a grid about a widest graben's cells.
    let face = ((std::f64::consts::TAU * r / (4.0 * widest)).round() as u32).clamp(8, 256);
    let mut rifts = crate::stages::rifting::RiftProcess::for_lid(graph, thermal, rules, face);
    // The extension spent in increments of its whole as the lid took it up, the faults
    // standing taking up what falls within their shadows before new ones open (F233), so
    // that the lid's faults grow and link as it stretches.
    let n = RIFT_GROWTH_INCREMENTS as usize;
    for i in 0..n {
        let span = (
            strain.extension_time(i as f64 / n as f64),
            strain.extension_time((i + 1) as f64 / n as f64),
        );
        rifts.stretch(total * area / n as f64, span);
    }
    let mut faults = rifts.finish(false);
    // The older ones are drawn first.
    faults.sort_by(|a, b| a.start_gyr.total_cmp(&b.start_gyr));
    faults
}

/// The lid's faults through the coarse loop (F153): the wide ones, which the loop's grid
/// resolves, carved as they form so that the floods pond in them and the relaxation sees
/// them, as a rise's canyons are (F105); the narrow ones left to the finish.
pub struct StrainProcess {
    pub canyons: Vec<Trough>,
    pub graben: Vec<Trough>,
    carved: Vec<f32>,
    stamper: crate::crater::stamp::Stamper,
}

impl StrainProcess {
    /// None where the lid opens no fault.
    pub fn new(graph: &ParamGraph, grid: CubeSphere, thermal: &ThermalHistory) -> Option<Self> {
        let faults = lid_faults(graph, grid, thermal);
        if faults.is_empty() {
            return None;
        }
        let cell = grid.delta() * graph.f64("body.radius");
        let (canyons, graben): (Vec<Trough>, Vec<Trough>) = faults
            .into_iter()
            .partition(|f| f.width >= RIFT_LOOP_WIDTH_CELLS * cell);
        Some(Self {
            canyons,
            graben,
            carved: vec![0.0; grid.cell_count()],
            stamper: crate::crater::stamp::Stamper::new(&grid),
        })
    }

    /// Carve the wide faults that form in the step.
    pub fn step(
        &mut self,
        layers: &mut crate::fields::LayerSet,
        grid: CubeSphere,
        body_radius: f64,
        t0: f64,
        t1: f64,
    ) {
        let cell = grid.delta() * body_radius;
        let t_mid = 0.5 * (t0 + t1);
        for fault in &self.canyons {
            let p1 = crate::stages::plume::trough_progress(fault, t1);
            let p0 = crate::stages::plume::trough_progress(fault, t0);
            if p1 <= p0 {
                continue;
            }
            let target = crate::stages::plume::trough_depths(
                &grid,
                body_radius,
                fault,
                p1,
                cell,
                0.0,
                &mut self.stamper,
            );
            let bedrock = layers.bedrock.as_mut_slice();
            let ages = layers.surface_age.as_mut_slice();
            for (idx, depth, face) in target {
                let more = depth - self.carved[idx] as f64;
                if more <= 0.0 {
                    continue;
                }
                bedrock[idx] -= more as f32;
                self.carved[idx] += more as f32;
                if face > 0.5 && more > CANYON_RENEWAL_THICKNESS {
                    ages[idx] = t_mid as f32;
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::angle_between;

    use crate::stages::plume::TroughKind;

    /// Charon's lid stretches as its ocean freezes, by at least 1 percent of its area
    /// (Beyer et al. 2017), and breaks into chasmata 10 to 50 km wide, the long ones that cut
    /// its layer kilometres deep (F233), Serenity the widest, Mandjet's 450 km among them,
    /// before its southern plain's flooding.
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
            faults.len() >= 50 && faults.len() <= 500,
            "{}",
            faults.len()
        );
        let median = |mut v: Vec<f64>| {
            v.sort_by(f64::total_cmp);
            v[v.len() / 2]
        };
        let width = median(faults.iter().map(|f| f.width).collect());
        let long = median(
            faults
                .iter()
                .filter(|f| f.length > 200e3)
                .map(|f| f.depth)
                .collect(),
        );
        let longest = faults.iter().map(|f| f.length).fold(0.0, f64::max);
        let widest = faults.iter().map(|f| f.width).fold(0.0, f64::max);
        assert!(width > 15e3 && width < 35e3, "{width}");
        assert!(widest > 35e3 && widest < 60e3, "{widest}");
        assert!(long > 1.5e3 && long < 6e3, "{long}");
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

    /// F232: an ocean freezes on its floor as the high-pressure ices where its water passes
    /// their melting point at its temperature, which the shell's base sets. Ganymede's,
    /// Callisto's and Titan's oceans lie between ice I and the high-pressure ices, so that as
    /// they cool their floors rise as their shells thicken and their lids hardly stretch
    /// (Squyres 1980: Callisto's lack of grooves); Pluto's and Charon's reach the rock, and
    /// their lids stretch as their oceans freeze (Charon's some 1 percent of its area, Beyer
    /// et al. 2017; Pluto's extension, Bierson et al. 2020).
    #[test]
    fn an_ocean_freezes_on_its_floor_as_the_high_pressure_ices() {
        // The melting curves: ice Ih at 100 MPa (IAPWS R14-08: 264.2 K), ice V at 260 K.
        assert!((ice_one_melting(100e6) - 264.2).abs() < 0.1);
        assert!((high_pressure_melting(260.0) - 403e6).abs() < 2e6);
        for (name, floored) in [
            ("Ganymede", true),
            ("Callisto", true),
            ("Titan", true),
            ("Pluto", false),
            ("Charon", false),
        ] {
            let (g, _) = crate::presets::preset(name).unwrap().to_graph();
            let thermal = ThermalHistory::from_graph(&g);
            let age = g.f64("body.age");
            let (ice, gravity) = (g.f64("body.ice_layer_thickness"), g.f64("body.gravity"));
            let (h0, h1) = (thermal.at(0.0).lithosphere, thermal.at(age).lithosphere);
            let (f0, f1) = (ocean_floor(h0, gravity), ocean_floor(h1, gravity));
            let extension = LidStrain::from_graph(&g, &thermal)
                .extension
                .last()
                .copied()
                .unwrap_or(0.0);
            eprintln!(
                "{name}: shell {:.0} to {:.0} km, floor {:.0} to {:.0} km in {:.0}, extension {:.4}",
                h0 / 1e3,
                h1 / 1e3,
                f0 / 1e3,
                f1 / 1e3,
                ice / 1e3,
                extension
            );
            if floored {
                assert!(f1 < f0 && f1 < ice && f1 > h1, "{name}");
                assert!(extension < 1e-3, "{name} {extension}");
            } else {
                assert!(f0 > ice, "{name}");
                assert!(extension > 0.01, "{name} {extension}");
            }
        }
    }

    /// F233: a lid's faults grow and link as it stretches, each as deep as its length allows
    /// until it cuts the lid's brittle layer: Pluto's long faults stand as wide as its layer
    /// lets them and deeper than its short ones, and a lid stretched little opens few faults
    /// rather than one each time it stretches.
    #[test]
    fn a_lid_s_faults_grow_deeper_with_their_length() {
        let (g, _) = crate::presets::preset("Pluto").unwrap().to_graph();
        let thermal = ThermalHistory::from_graph(&g);
        let faults = lid_faults(&g, CubeSphere::new(128), &thermal);
        let widest =
            2.0 * closure_depth(&g).min(shell_bound(&g)) / GRABEN_DIP_DEG.to_radians().tan();
        let mut by_length: Vec<&Trough> = faults.iter().collect();
        by_length.sort_by(|a, b| a.length.total_cmp(&b.length));
        let tenth = by_length.len() / 10;
        let mean = |set: &[&Trough], f: fn(&Trough) -> f64| {
            set.iter().map(|t| f(t)).sum::<f64>() / set.len().max(1) as f64
        };
        let (short, long) = (&by_length[..tenth], &by_length[by_length.len() - tenth..]);
        let (short_depth, long_depth) = (mean(short, |t| t.depth), mean(long, |t| t.depth));
        let long_width = mean(long, |t| t.width);
        eprintln!(
            "{} faults; short tenth {:.0} km long {short_depth:.0} m deep; long tenth {:.0} km long {long_depth:.0} m deep {long_width:.0} m wide of {widest:.0}",
            faults.len(),
            mean(short, |t| t.length) / 1e3,
            mean(long, |t| t.length) / 1e3,
        );
        assert!(long_depth > 1.5 * short_depth);
        assert!(long_width > 0.8 * widest && long_width <= widest * (1.0 + 1e-9));
        // Europa's lid barely stretches: a few faults, not one each increment.
        let (europa, _) = crate::presets::preset("Europa").unwrap().to_graph();
        let thermal = ThermalHistory::from_graph(&europa);
        let few = lid_faults(&europa, CubeSphere::new(128), &thermal).len();
        assert!(few < 16, "{few}");
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

    /// Charon's wide chasmata are carved through the coarse loop (F153): the strain process
    /// carries the faults the loop's grid resolves, and a step through a fault's time drops
    /// its floor and renews its faces.
    #[test]
    fn charon_s_wide_chasmata_are_carved_through_the_loop() {
        let (g, _) = crate::presets::preset("Charon").unwrap().to_graph();
        let thermal = ThermalHistory::from_graph(&g);
        let grid = CubeSphere::new(128);
        let r = g.f64("body.radius");
        let mut process = StrainProcess::new(&g, grid, &thermal).expect("Charon's lid faults");
        let cell = grid.delta() * r;
        assert!(!process.canyons.is_empty() && !process.graben.is_empty());
        assert!(
            process
                .canyons
                .iter()
                .all(|c| c.width >= RIFT_LOOP_WIDTH_CELLS * cell)
        );
        assert!(process.canyons.iter().all(|c| c.end_gyr > c.start_gyr));
        let widest = process
            .canyons
            .iter()
            .max_by(|a, b| a.width.total_cmp(&b.width))
            .unwrap()
            .clone();
        let mut layers = crate::fields::LayerSet::new(grid);
        let before = layers.bedrock.as_slice().to_vec();
        let (t0, t1) = (widest.start_gyr - 0.01, widest.start_gyr + 0.01);
        process.step(&mut layers, grid, r, t0, t1);
        let centre = widest.points[widest.points.len() / 2];
        let i = grid.index(grid.locate(centre).0);
        let dropped = (before[i] - layers.bedrock.as_slice()[i]) as f64;
        assert!(
            dropped > 0.5 * widest.depth,
            "{dropped} of {}",
            widest.depth
        );
        // A step before its time carves it no more (other faults may form in that step).
        let after = layers.bedrock.as_slice()[i];
        process.step(&mut layers, grid, r, t0 - 0.1, t0);
        assert_eq!(layers.bedrock.as_slice()[i], after);
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
            let total_km: f64 = l.iter().sum();
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
                "{name}: radius {:+.2} km through the record; extension {:.3}% shortening {:.3}% of the area; {} faults ({total_km:.0} km) in {:.2}s; length median {:.0} max {:.0} km; width median {:.0} max {:.0} km; depth median {:.1} deepest {:.1} km; times q10 {:.2} q50 {:.2} q90 {:.2}",
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
