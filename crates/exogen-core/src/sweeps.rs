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
