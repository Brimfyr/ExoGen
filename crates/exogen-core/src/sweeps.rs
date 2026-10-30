// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! Interpolated bodies (principle 15, audit V1): sweeps of a continuous input between two
//! reference bodies, read on the parameter graph, so that a class of feature can be seen to
//! come in by degrees where a named physical threshold is crossed, or to switch on at a cut.

#[cfg(test)]
mod tests {
    use crate::params::ParamGraph;
    use crate::presets::preset;

    /// A preset's graph with one input set.
    fn swept(name: &str, key: &str, value: f64) -> ParamGraph {
        let (mut g, _) = preset(name).unwrap().to_graph();
        g.set(key, value).unwrap();
        g
    }

    /// Between Ganymede's thick shell and Europa's thin one a shell renews itself by degrees
    /// (audit S4, S6): as Europa's tides weaken and its shell thickens against its radius,
    /// chaos gives way to lanes and the craters' retention grows, with no step between
    /// neighbouring bodies of the sweep.
    #[test]
    fn a_shell_renews_itself_by_degrees_between_europa_and_ganymede() {
        let (europa, _) = preset("Europa").unwrap().to_graph();
        let e0 = europa.f64("orbit.eccentricity");
        let mut last: Option<(f64, f64, f64)> = None;
        for k in 0..=24 {
            let f = 0.3 + 0.05 * k as f64;
            let g = swept("Europa", "orbit.eccentricity", e0 * f);
            let (chaos, lanes, retention) = (
                g.f64("tectonics.chaos_fraction"),
                g.f64("tectonics.ice_resurfaced_fraction"),
                g.f64("cratering.retention_age"),
            );
            eprintln!(
                "{f:.2}: shell {:.4} of the radius, tide {:.1} mW/m2, chaos {chaos:.3}, lanes {lanes:.3}, retention {retention:.3} Gyr",
                g.f64("interior.ice_shell_thickness") / g.f64("body.radius"),
                1e3 * g.f64("interior.tidal_flux"),
            );
            if let Some((c0, l0, r0)) = last {
                assert!((chaos - c0).abs() < 0.12, "chaos {c0} to {chaos} at {f}");
                assert!((lanes - l0).abs() < 0.3, "lanes {l0} to {lanes} at {f}");
                assert!(
                    retention / r0 < 3.0 && r0 / retention < 3.0,
                    "retention {r0} to {retention} at {f}"
                );
            }
            last = Some((chaos, lanes, retention));
        }
        let europa_now = swept("Europa", "orbit.eccentricity", e0);
        assert!((europa_now.f64("tectonics.chaos_fraction") - 0.3).abs() < 1e-9);
        assert!((europa_now.f64("cratering.retention_age") - 0.1).abs() < 1e-9);
    }

    /// The focusing's episodes and the plumes begin where the thinned shell lets the ocean
    /// reach the surface, and the plumes' fallout grows from there (audit S5): as
    /// Enceladus's tides are turned up from a fifth of its own, the cracks' reach passes one
    /// and the fallout comes in by degrees.
    #[test]
    fn plumes_begin_where_the_ocean_reaches_the_surface() {
        let (enceladus, _) = preset("Enceladus").unwrap().to_graph();
        let e0 = enceladus.f64("orbit.eccentricity");
        let mut last_fallout = 0.0;
        let mut began = false;
        for k in 0..=16 {
            let f = 0.2 + 0.05 * k as f64;
            let g = swept("Enceladus", "orbit.eccentricity", e0 * f);
            let reach = g.f64("tectonics.crack_reach");
            let fallout = g.f64("volcanism.plume_fallout_rate");
            let period = g.f64("tectonics.focusing_period");
            eprintln!(
                "{f:.2}: thinning {:.3}, reach {reach:.3}, fallout {fallout:.0} m/Gyr, focusing every {period:.3} Gyr",
                g.f64("tectonics.polar_thinning")
            );
            assert_eq!(
                fallout > 0.0,
                reach > 1.0 && g.f64("tectonics.ice_crack_density") > 0.0
            );
            assert_eq!(period > 0.0, reach >= 1.0);
            began |= reach > 1.0;
            if fallout > 0.0 && last_fallout > 0.0 {
                assert!(
                    fallout <= 2.5 * last_fallout + 1.0,
                    "fallout {last_fallout} to {fallout} at {f}"
                );
            }
            last_fallout = fallout;
        }
        assert!(began, "the sweep never let the ocean through");
    }

    /// A crater's target goes from ice to rock by degrees as its crust's share of rock grows
    /// (audit S9): along Callisto made denser, its ice degree falls from one to none, and with
    /// it the transition's ice factor, the largest crater it holds and, in a crater at one
    /// and a half times its peak-ring onset, the central pit gives way to a peak ring, none of
    /// them by a step.
    #[test]
    fn craters_go_from_ice_to_rock_by_degrees() {
        use crate::crater::{Scaling, dimensions, hold_diameter};
        let (callisto, _) = preset("Callisto").unwrap().to_graph();
        let rho0 = callisto.f64("body.density");
        let mut last: Option<[f64; 5]> = None;
        for k in 0..=32 {
            let g = swept("Callisto", "body.density", rho0 + 50.0 * k as f64);
            let ice = g.f64("cratering.target_ice");
            let scaling = Scaling {
                simple_complex: g.f64("cratering.simple_complex_transition"),
                peak_ring: g.f64("cratering.peak_ring_onset"),
                multi_ring: g.f64("cratering.multiring_onset"),
                gravity: g.f64("body.gravity"),
                ice,
            };
            // The transition with the body's gravity taken out, and the largest crater
            // against the radius.
            let strength = scaling.simple_complex * scaling.gravity;
            let hold = hold_diameter(g.f64("body.radius"), ice) / g.f64("body.radius");
            let dim = dimensions(1.5 * scaling.peak_ring, &scaling);
            let ring = if dim.ring_count > 0 {
                dim.rings[0].1 / dim.depth
            } else {
                0.0
            };
            eprintln!(
                "{:.0} kg/m3: {} {}, rock {:.3}, ice degree {ice:.3}, transition {:.2} km (x g {strength:.0}), holds {hold:.3} R, {}: pit {:.3} of it, ring {ring:.3} of the depth",
                g.f64("body.density"),
                g.str("body.composition_class"),
                g.str("materials.bedrock"),
                g.f64("body.crust_impurity"),
                scaling.simple_complex / 1e3,
                dim.class.name(),
                dim.pit_radius / dim.radius,
            );
            let pit = dim.pit_radius / dim.radius;
            // Before F229 the transition stepped by 6.7 times and the largest crater by pi.
            if let Some([ice0, strength0, hold0, pit0, ring0]) = last {
                assert!(ice <= ice0 + 1e-12 && ice0 - ice < 0.15, "{ice0} to {ice}");
                let step = |a: f64, b: f64| (a / b).max(b / a);
                assert!(step(strength, strength0) < 1.3, "{strength0} to {strength}");
                assert!(step(hold, hold0) < 1.2, "{hold0} to {hold}");
                assert!(pit <= pit0 + 1e-12 && pit0 - pit < 0.05, "{pit0} to {pit}");
                assert!(
                    ring >= ring0 - 1e-12 && ring - ring0 < 0.05,
                    "{ring0} to {ring}"
                );
            }
            last = Some([ice, strength, hold, pit, ring]);
        }
        let first = swept("Callisto", "body.density", rho0);
        assert_eq!(first.f64("cratering.target_ice"), 1.0);
        assert_eq!(last.unwrap()[0], 0.0);
    }

    /// What an overturn does follows the overturns, and the episodic lid's constants its
    /// degree (audit S1): along Venus's radius from 0.6 of its own, the lid breaks from 0.7
    /// and has broken ever later, its melt's efficiency, its rises' width and its domes
    /// coming in by degrees as its last overturn nears the present, none by the cut on half
    /// its life that stepped its magma fivefold.
    #[test]
    fn a_lid_turns_episodic_by_degrees() {
        let (venus, _) = preset("Venus").unwrap().to_graph();
        let r0 = venus.f64("body.radius");
        assert_eq!(venus.f64("interior.episodic_degree"), 1.0);
        let mut last: Option<[f64; 4]> = None;
        for k in 0..=20 {
            let f = 0.6 + 0.025 * k as f64;
            let g = swept("Venus", "body.radius", r0 * f);
            let degree = g.f64("interior.episodic_degree");
            let magma = g.f64("interior.magma_production_rate");
            let domes = g.f64("volcanism.dome_count");
            let rise = g.f64("tectonics.plume_rise_diameter") / g.f64("body.radius");
            let broken = g.str("interior.tectonic_regime") == "episodic_lid";
            let since = g.f64("tectonics.resurfacing_event_age");
            eprintln!(
                "{f:.3}: {}, last overturn {since:.2} Gyr ago, degree {degree:.3}, magma {magma:.1} km3/yr, domes {domes}, rises {rise:.2} R, tesserae {:.2}, retention {:.2} Gyr",
                g.str("interior.tectonic_regime"),
                g.f64("tectonics.tessera_fraction"),
                g.f64("cratering.retention_age")
            );
            // What an overturn does follows the overturns themselves.
            assert_eq!(g.f64("tectonics.tessera_fraction") > 0.0, broken);
            if broken {
                assert!((g.f64("cratering.retention_age") - since.min(3.5)).abs() < 1e-9);
            }
            if let Some([degree0, magma0, domes0, rise0]) = last {
                let step = |a: f64, b: f64| (a / b).max(b / a);
                assert!((degree - degree0).abs() < 0.5, "{degree0} to {degree}");
                assert!(step(magma, magma0) < 2.5, "{magma0} to {magma}");
                assert!((domes - domes0).abs() < 80.0, "{domes0} to {domes}");
                assert!(step(rise, rise0) < 1.8, "{rise0} to {rise}");
            }
            last = Some([degree, magma, domes, rise]);
        }
    }

    /// Prints the derived nodes the switches read along each sweep.
    #[test]
    #[ignore]
    fn print_switch_sweeps() {
        // Moon to Io (S3, S10): Io's tidal heating from none to half again its own, by its
        // eccentricity.
        let (io, _) = preset("Io").unwrap().to_graph();
        let e0 = io.f64("orbit.eccentricity");
        eprintln!(
            "Io by eccentricity: factor, tidal W/m2, regime, paterae, mountains, magma m3/s, retention Gyr, bedrock"
        );
        for f in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.5] {
            let g = swept("Io", "orbit.eccentricity", e0 * f);
            eprintln!(
                "  {f:.1}: {:.3} {} {:.0} {:.0} {:.1} {:.3} {}",
                g.f64("interior.tidal_flux"),
                g.str("interior.tectonic_regime"),
                g.f64("volcanism.paterae_count"),
                g.f64("tectonics.mountain_count"),
                g.f64("interior.magma_production_rate"),
                g.f64("cratering.retention_age"),
                g.str("materials.bedrock"),
            );
        }
        // Ganymede to Europa (S4 to S6): Europa's tidal heating from a fifth of its own, by
        // its eccentricity, its shell thickening against its radius.
        let (europa, _) = preset("Europa").unwrap().to_graph();
        let e0 = europa.f64("orbit.eccentricity");
        eprintln!(
            "Europa by eccentricity: factor, tidal mW/m2, shell/R, regime, chaos, lanes, thinning, focusing Gyr, fallout, cracks, retention Gyr"
        );
        for f in [0.2, 0.4, 0.6, 0.8, 0.9, 1.0, 1.2, 1.5] {
            let g = swept("Europa", "orbit.eccentricity", e0 * f);
            eprintln!(
                "  {f:.1}: {:.1} {:.4} {} {:.3} {:.3} {:.3} {:.3} {:.2e} {:.2e} {:.3}",
                1e3 * g.f64("interior.tidal_flux"),
                g.f64("interior.ice_shell_thickness") / g.f64("body.radius"),
                g.str("interior.tectonic_regime"),
                g.f64("tectonics.chaos_fraction"),
                g.f64("tectonics.ice_resurfaced_fraction"),
                g.f64("tectonics.polar_thinning"),
                g.f64("tectonics.focusing_period"),
                g.f64("volcanism.plume_fallout_rate"),
                g.f64("tectonics.ice_crack_density"),
                g.f64("cratering.retention_age"),
            );
        }
        // Mars to Venus (S1): Mars made larger, its mass with its volume.
        let (mars, _) = preset("Mars").unwrap().to_graph();
        let r0 = mars.f64("body.radius");
        eprintln!(
            "Mars by radius: factor, regime, retention Gyr, event age Gyr, tesserae, coronae, domes, rises, magma m3/s"
        );
        for f in [1.0, 1.2, 1.4, 1.5, 1.6, 1.7, 1.8] {
            let g = swept("Mars", "body.radius", r0 * f);
            eprintln!(
                "  {f:.1}: {} {:.3} {:.3} {:.3} {:.0} {:.0} {:.0} {:.1}",
                g.str("interior.tectonic_regime"),
                g.f64("cratering.retention_age"),
                g.f64("tectonics.resurfacing_event_age"),
                g.f64("tectonics.tessera_fraction"),
                g.f64("volcanism.corona_count"),
                g.f64("volcanism.dome_count"),
                g.f64("tectonics.plume_rise_count"),
                g.f64("interior.magma_production_rate"),
            );
        }
        // Rock to ice (S2, S9): Callisto's density from its own to Io's.
        let (callisto, _) = preset("Callisto").unwrap().to_graph();
        let d0 = callisto.f64("body.density");
        eprintln!(
            "Callisto by density: kg/m3, class, ice fraction, moment of inertia, regime, bedrock, shell km"
        );
        for d in [
            d0, 2000.0, 2200.0, 2400.0, 2600.0, 2800.0, 3000.0, 3300.0, 3530.0,
        ] {
            let g = swept("Callisto", "body.density", d);
            eprintln!(
                "  {d:.0}: {} {:.3} {:.3} {} {} {:.0}",
                g.str("body.composition_class"),
                g.f64("body.ice_mass_fraction"),
                g.f64("body.moment_of_inertia_factor"),
                g.str("interior.tectonic_regime"),
                g.str("materials.bedrock"),
                g.f64("interior.ice_shell_thickness") / 1e3,
            );
        }
    }
}
