// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! The spin's history (PLAN.md 5A, F238). The present rotation is measured, an input, and
//! its history is derived from it by the tides that change a spin, each a torque of constant
//! Q (Gladman et al. 1996; Murray and Dermott 1999): the tide its primary raises on it, and
//! the tide of its largest satellite, whose orbit takes up the angular momentum the spin
//! loses. A spin the tides hold, locked to its primary or in a resonance with its orbit
//! (Mercury's 3:2), was brought there from the spin it was born with; a free spin is run back
//! from today until the body formed, or until its satellite stands at its Roche limit, where
//! the satellite formed. The history tells when a spin was captured, what it was born with
//! and when a satellite formed; the plausibility rules flag a present state it cannot reach.

use std::f64::consts::TAU;

use crate::constants::*;
use crate::params::Ctx;
use crate::params::value::Curve;

/// What a spin's history reads.
#[derive(Clone, Copy, Debug)]
pub struct SpinInputs {
    /// The body's mass, kg, radius, m, and moment of inertia, kg m².
    pub mass: f64,
    pub radius: f64,
    pub moment: f64,
    /// The body's tidal Love number and its Q over its history.
    pub k2: f64,
    pub q: f64,
    /// The primary's mass, kg, the orbit's semi-major axis, m, and mean motion, rad/s.
    pub primary_mass: f64,
    pub distance: f64,
    pub mean_motion: f64,
    /// The largest satellite's mass, kg, and present distance, m; none at zero mass.
    pub satellite_mass: f64,
    pub satellite_distance: f64,
    /// The present spin, rad/s, negative when it turns against the orbit.
    pub spin: f64,
    /// The spin rate the tides hold it at, rad/s: a resonance with its primary's orbit, or
    /// its satellite's orbital rate when its spin is locked to the satellite; none free.
    pub held: Option<f64>,
    /// The body's age, Gyr.
    pub age: f64,
}

/// A spin's history.
#[derive(Clone, Debug)]
pub struct SpinHistory {
    /// The sidereal period through the body's life, Gyr after formation against hours.
    pub period: Curve,
    /// The satellite's distance through it, Gyr against m, when it has one.
    pub satellite: Option<Curve>,
    /// When the tides brought a held spin to its present state, Gyr after formation; none
    /// when they could not within the body's age, or the spin is free.
    pub capture: Option<f64>,
    /// The spin it was born with, rad/s.
    pub primordial: f64,
    /// When the satellite formed at its Roche limit, Gyr after formation, if the history
    /// reaches it.
    pub satellite_formed: Option<f64>,
}

impl SpinInputs {
    /// The body, its orbit and its satellite, without its present spin: what its locking
    /// time reads, which the default period reads in turn.
    pub fn body(c: &Ctx) -> Self {
        let mass = c.f64("body.mass");
        let radius = c.f64("body.radius");
        // An ocean's tides dissipate today some three times what they did on average.
        let q = if c.str("volatiles.working_liquid") == "water" {
            OCEAN_TIDE_Q_MEAN
        } else {
            c.f64("interior.dissipation_q")
        };
        Self {
            mass,
            radius,
            moment: c.f64("body.moment_of_inertia_factor") * mass * radius * radius,
            k2: c.f64("interior.love_number_k2"),
            q,
            primary_mass: c.f64("orbit.primary_mass"),
            distance: c.f64("orbit.semi_major_axis"),
            mean_motion: c.f64("orbit.mean_motion"),
            satellite_mass: c.f64("satellite.mass").max(0.0),
            satellite_distance: c.f64("satellite.semi_major_axis").max(0.0),
            spin: 0.0,
            held: None,
            age: c.f64("body.age"),
        }
    }

    /// The body with its present spin: held where it is in a resonance with its primary's
    /// orbit, or locked to its satellite's (Pluto to Charon), which a history run back from
    /// today would leave, as a lock is left going back.
    pub fn from_ctx(c: &Ctx) -> Self {
        let body = Self::body(c);
        let period = c.f64("rotation.period");
        let sense = if c.f64("rotation.obliquity") > 90.0 {
            -1.0
        } else {
            1.0
        };
        let spin = sense * TAU / period;
        let satellite = (body.satellite_mass > 0.0 && body.satellite_distance > 0.0)
            .then(|| body.satellite_motion(body.satellite_distance));
        // Held at the present rate itself, which lies within the resonance's tolerance.
        let resonant = crate::insolation::spin_orbit_resonance(c.f64("orbit.period") / period);
        let locked = satellite.is_some_and(|n| {
            crate::insolation::spin_orbit_resonance(n / spin.abs()) == Some((1, 1))
        });
        let held = (resonant.is_some() || locked).then_some(spin.abs());
        Self { spin, held, ..body }
    }

    /// The torque the tide a body of `mass` at `distance` raises on this one exerts, N m,
    /// slowing a spin faster than the body's orbit (Murray and Dermott 1999, eq. 4.159).
    pub fn tide(&self, mass: f64, distance: f64) -> f64 {
        if mass <= 0.0 || distance <= 0.0 {
            return 0.0;
        }
        1.5 * self.k2 / self.q * G * mass * mass * self.radius.powi(5) / distance.powi(6)
    }

    /// The rate at which the body would break up, rad/s.
    pub fn breakup(&self) -> f64 {
        (G * self.mass / self.radius.powi(3)).sqrt()
    }

    /// The spin the body was born with: a share of its breakup rate.
    pub fn primordial(&self) -> f64 {
        PRIMORDIAL_SPIN_SHARE * self.breakup()
    }

    /// The time its primary's tide takes to bring the spin it was born with to its orbit's
    /// rate, Gyr.
    pub fn locking_time(&self) -> f64 {
        let torque = self.tide(self.primary_mass, self.distance);
        if torque <= 0.0 {
            return f64::INFINITY;
        }
        (self.primordial() - self.mean_motion).abs() * self.moment / torque / GYR
    }

    /// The time its satellite's tide takes to bring the spin it was born with to the
    /// satellite's present orbital rate, Gyr; infinite without a satellite.
    pub fn satellite_locking_time(&self) -> f64 {
        let torque = self.tide(self.satellite_mass, self.satellite_distance);
        if torque <= 0.0 {
            return f64::INFINITY;
        }
        (self.primordial() - self.satellite_motion(self.satellite_distance)).abs() * self.moment
            / torque
            / GYR
    }

    /// The satellite's mean motion at a distance, rad/s.
    pub fn satellite_motion(&self, distance: f64) -> f64 {
        (G * (self.mass + self.satellite_mass) / distance.powi(3)).sqrt()
    }

    /// The satellite's fluid Roche limit, m: where it formed.
    fn roche(&self) -> f64 {
        let density = self.mass / (4.0 / 3.0 * std::f64::consts::PI * self.radius.powi(3));
        ROCHE_FLUID * self.radius * (density / SATELLITE_DENSITY).cbrt()
    }

    /// The spin's history.
    pub fn history(&self) -> SpinHistory {
        let hours = |spin: f64| TAU / spin.abs().max(1e-30) / 3600.0;
        if let Some(target) = self.held {
            // Brought from the spin it was born with to its present state by its tides, which
            // slow it while it turns faster than the rate they hold it at.
            let born = self.primordial();
            let torque = self.tide(self.primary_mass, self.distance)
                + self.tide(self.satellite_mass, self.satellite_distance);
            let rate = torque / self.moment;
            let time = if rate > 0.0 {
                (born - target).abs() / rate / GYR
            } else {
                f64::INFINITY
            };
            let reached = time <= self.age;
            let end = if reached {
                target
            } else {
                born - (born - target).signum() * rate * self.age * GYR
            };
            let mut points = vec![(0.0, hours(born))];
            if reached {
                points.push((time, hours(target)));
            }
            points.push((self.age, hours(end)));
            points.dedup_by(|b, a| b.0 <= a.0);
            return SpinHistory {
                period: Curve::new(points),
                satellite: None,
                capture: reached.then_some(time),
                primordial: born,
                satellite_formed: None,
            };
        }
        // A free spin run back from today, its satellite's orbit with it.
        let reduced = self.mass * self.satellite_mass / (self.mass + self.satellite_mass);
        let roche = self.roche();
        let (mut t, mut spin, mut distance) = (self.age, self.spin, self.satellite_distance);
        let mut spins = vec![(t, spin)];
        let mut distances = vec![(t, distance)];
        let mut formed = None;
        let sign = |x: f64, scale: f64| {
            if x.abs() <= 1e-9 * scale {
                0.0
            } else {
                x.signum()
            }
        };
        while t > 0.0 {
            let primary = sign(spin - self.mean_motion, spin.abs())
                * self.tide(self.primary_mass, self.distance);
            let (satellite, recession) = if self.satellite_mass > 0.0 && distance > 0.0 {
                let torque = sign(spin - self.satellite_motion(distance), spin.abs())
                    * self.tide(self.satellite_mass, distance);
                let gain =
                    reduced * 0.5 * (G * (self.mass + self.satellite_mass) / distance).sqrt();
                (torque, torque / gain)
            } else {
                (0.0, 0.0)
            };
            // Forward in time the spin falls at this rate and the satellite recedes at that.
            let fall = (primary + satellite) / self.moment;
            let mut dt = SPIN_STEP_GYR * GYR;
            if fall != 0.0 {
                dt = dt.min(SPIN_STEP_SHARE * spin.abs() / fall.abs());
            }
            if recession != 0.0 {
                dt = dt.min(SPIN_STEP_SHARE * distance / recession.abs());
            }
            dt = dt.min(t * GYR).max(1.0);
            spin += fall * dt;
            distance -= recession * dt;
            t = (t - dt / GYR).max(0.0);
            if self.satellite_mass > 0.0 && distance <= roche {
                formed = Some(t);
                spins.push((t, spin));
                distances.push((t, roche));
                break;
            }
            spins.push((t, spin));
            distances.push((t, distance));
        }
        let primordial = spins.last().map_or(self.spin, |s| s.1);
        let sample = |values: &[(f64, f64)], f: &dyn Fn(f64) -> f64| -> Curve {
            let start = values.last().map_or(0.0, |v| v.0);
            let mut ordered: Vec<(f64, f64)> = values.iter().rev().copied().collect();
            ordered.dedup_by(|b, a| b.0 <= a.0);
            let n = SPIN_HISTORY_SAMPLES;
            let points = (0..n)
                .map(|i| {
                    let t = start + (self.age - start) * i as f64 / (n - 1) as f64;
                    let j = ordered.partition_point(|p| p.0 < t).min(ordered.len() - 1);
                    let value = if j == 0 {
                        ordered[0].1
                    } else {
                        let (a, b) = (ordered[j - 1], ordered[j]);
                        let w = if b.0 > a.0 {
                            (t - a.0) / (b.0 - a.0)
                        } else {
                            1.0
                        };
                        a.1 + (b.1 - a.1) * w
                    };
                    (t, f(value))
                })
                .collect::<Vec<_>>();
            let mut points = points;
            points.dedup_by(|b, a| b.0 <= a.0);
            Curve::new(points)
        };
        SpinHistory {
            period: sample(&spins, &hours),
            satellite: (self.satellite_mass > 0.0).then(|| sample(&distances, &|d| d)),
            capture: None,
            primordial,
            satellite_formed: formed,
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::params::ParamGraph;

    fn graph(name: &str) -> ParamGraph {
        crate::presets::preset(name).unwrap().to_graph().0
    }

    /// Prints each preset's spin and its history: its state, its locking times, when it was
    /// captured, what it was born with, and its period at a few ages before the present.
    #[test]
    #[ignore = "diagnostic: the presets' spin histories"]
    fn print_spins() {
        for name in crate::presets::PRESET_NAMES {
            let g = graph(name);
            let age = g.f64("body.age");
            let at = |ago: f64| match g.get("rotation.period_history") {
                crate::params::Value::Curve(c) => c.eval((age - ago).max(0.0)),
                _ => f64::NAN,
            };
            println!(
                "{name:10} {:9} {:9.2} h  lock {:9.3} Gyr  satellite lock {:9.3}  captured {:6.3}  born {:8.2} h  at 0.62/2/4 Gyr ago {:6.2} {:6.2} {:6.2} h  satellite formed {:6.3}",
                g.str("rotation.mode"),
                g.f64("rotation.period") / 3600.0,
                g.f64("rotation.locking_timescale"),
                g.f64("satellite.locking_timescale"),
                g.f64("rotation.capture_time"),
                g.f64("rotation.primordial_period") / 3600.0,
                at(0.62),
                at(2.0),
                at(4.0),
                g.f64("satellite.formation_time"),
            );
        }
    }

    /// F238: a spin locked to its satellite stays locked going back to its capture: Pluto's
    /// to Charon's orbit, its period through its history its present one.
    #[test]
    fn pluto_has_been_locked_to_charon_since_its_capture() {
        let g = graph("Pluto");
        let capture = g.f64("rotation.capture_time");
        assert!((0.0..0.1).contains(&capture), "capture {capture}");
        let today = g.f64("rotation.period") / 3600.0;
        match g.get("rotation.period_history") {
            crate::params::Value::Curve(c) => {
                assert!(
                    (c.eval(2.0) / today - 1.0).abs() < 1e-6,
                    "{} against {today}",
                    c.eval(2.0)
                )
            }
            _ => panic!("no curve"),
        }
    }

    /// F238: Mercury's 3:2, read from its period, was reached from the spin it was born with
    /// within its age, and a spin its tides hold needs no pin.
    #[test]
    fn mercury_s_resonance_was_reached_within_its_age() {
        let g = graph("Mercury");
        assert_eq!(g.str("rotation.mode"), "resonant");
        let capture = g.f64("rotation.capture_time");
        assert!(
            capture > 0.0 && capture < g.f64("body.age"),
            "capture {capture}"
        );
        let born = g.f64("rotation.primordial_period") / 3600.0;
        assert!(born < 24.0, "born at {born} h");
    }

    /// F238: going back, the Moon comes nearer and Earth's day shorter: some 22 h when the
    /// tidal rhythmites were laid 620 Myr ago (Williams 2000), a few hours where the Moon
    /// stood at its Roche limit near Earth's age.
    #[test]
    fn the_moon_recedes_as_the_earth_slows() {
        let g = graph("Earth");
        let age = g.f64("body.age");
        let day = |ago: f64| match g.get("rotation.period_history") {
            crate::params::Value::Curve(c) => c.eval(age - ago),
            _ => panic!("no curve"),
        };
        assert!((day(0.0) - 23.93).abs() < 0.05, "today {}", day(0.0));
        assert!(
            (19.5..23.5).contains(&day(0.62)),
            "620 Myr ago {}",
            day(0.62)
        );
        let formed = g.f64("satellite.formation_time");
        assert!(
            (0.0..1.0).contains(&formed),
            "the Moon formed {formed} Gyr after Earth"
        );
        let born = g.f64("rotation.primordial_period") / 3600.0;
        assert!((2.0..8.0).contains(&born), "Earth's first day {born} h");
    }

    /// F238: a moon of a planet defaults to locked, its period its orbit's; the present state
    /// of every preset is one its history reaches, without a pin.
    #[test]
    fn the_presets_reach_their_present_spins() {
        for name in ["Moon", "Europa", "Titan", "Charon"] {
            let g = graph(name);
            assert_eq!(g.str("rotation.mode"), "locked", "{name}");
            assert!((g.f64("rotation.period") / g.f64("orbit.period") - 1.0).abs() < 1e-9);
        }
        for name in crate::presets::PRESET_NAMES {
            let g = graph(name);
            let flagged: Vec<String> = crate::params::rules::evaluate(&g)
                .into_iter()
                .filter(|w| w.node.starts_with("rotation") || w.node.starts_with("satellite"))
                .filter(|w| w.severity == crate::params::rules::Severity::Warning)
                .map(|w| w.text)
                .collect();
            assert!(flagged.is_empty(), "{name}: {flagged:?}");
        }
    }

    /// F238: what the history cannot reach is flagged: Mercury turning in a day, which its
    /// tides would have stopped long ago; an Earth with a twin at 30 Earth radii, which would
    /// have locked it to its orbit within a million years.
    #[test]
    fn a_spin_the_tides_would_have_captured_is_flagged() {
        let warned = |g: &ParamGraph| {
            crate::params::rules::evaluate(g)
                .into_iter()
                .any(|w| w.severity == crate::params::rules::Severity::Warning)
        };
        let mut mercury = graph("Mercury");
        mercury.set("rotation.period", 86_164.0).unwrap();
        assert_eq!(mercury.str("rotation.mode"), "free");
        assert!(warned(&mercury));
        let mut twins = graph("Earth");
        twins.set("satellite.mass", 5.972e24).unwrap();
        twins
            .set("satellite.semi_major_axis", 30.0 * 6.371e6)
            .unwrap();
        assert!(twins.f64("satellite.locking_timescale") < 0.01);
        assert!(warned(&twins));
    }
}
