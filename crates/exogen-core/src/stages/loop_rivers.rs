// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! The loop's rivers (F437): the rivers' erosion through the body's history, each of the coarse
//! loop's steps eroding as the finish's span does (F278), against the ground the loop's own
//! processes raise, under the year's water the body's climate gives on the loop's first ground
//! and the liquid its history keeps (F411), its inputs held at today's (the user's choice,
//! 2026-11-26: one year, at the loop's start). A mobile lid's ground the loop raises none of
//! until Phase 6's plates; its rivers wait for them (the user's choice, 2026-11-26).

use crate::constants::*;
use crate::fields::{Field, LayerSet};
use crate::fluvial::Erosion;
use crate::grid::CubeSphere;
use crate::params::ParamGraph;

/// Bumped when the loop's rivers change.
pub const VERSION: u64 = 1;

/// The year's water on the loop's ground (F437): what the land sheds and what an open water's
/// surface gains at each cell, kg/m²/s, none shed under ice; the liquid's density, kg/m³; and
/// the waves the year's wind raises on the liquid, if any.
pub struct Water {
    pub shed: Vec<f32>,
    pub gain: Vec<f32>,
    pub density: f64,
    pub waves: Option<crate::waves::Waves>,
}

/// What the loop's rivers moved through the history (F437): the ground as the loop began, m,
/// and its liquid's level then, m; the water they ran under, as `Water`'s; the ground they cut
/// at each cell of land less what they laid there, m; their totals, m³; and their time, s.
#[derive(Clone, Debug, Default)]
pub struct Worn {
    pub start: Vec<f32>,
    pub level: f64,
    pub shed: Vec<f32>,
    pub gain: Vec<f32>,
    pub density: f64,
    pub lowered: Vec<f32>,
    pub cut: f64,
    pub laid: f64,
    pub to_sea: f64,
    pub floor: f64,
    pub seconds: f64,
    /// What full seas held back at their mouths for the rivers of the step after, m³ at each
    /// cell: at the loop's end, what its rivers carried that no sea had room for.
    pub held: Vec<f64>,
}

/// Whether the loop raises the ground its rivers would wear (F437): every lid but a mobile one,
/// whose plates wait for Phase 6.
pub fn runs(graph: &ParamGraph) -> bool {
    graph.str("interior.tectonic_regime") != "mobile_lid"
}

/// Each of the cycle's places' ground over its liquid, m: `above`'s mean over the climate's
/// cell each place lies in.
pub(crate) fn place_heights(
    cycle: &crate::climate::Cycle,
    above: &Field<f32>,
    climate_grid: CubeSphere,
) -> Vec<f64> {
    let coarse = crate::climate::coarsen(above, climate_grid);
    cycle
        .sites
        .iter()
        .map(|s| {
            let d = crate::grid::direction_from_lon_lat(s.lon, s.lat);
            coarse.get(coarse.grid().locate(d).0) as f64
        })
        .collect()
}

/// The ground the year's runoff runs down (F434), and soaks into and flows through to the sea
/// at `level`, m (F435): `layers`' terrain and liquid on `grid` brought to YEAR_DRAINAGE_FACE,
/// each cell's cover of the rivers' sediments over the body's rock, frozen through where its
/// place's land stands below its liquid's melting point through `cycle`'s year, its liquid's
/// viscosity that land's, on a body of `radius` m.
pub(crate) fn drains_of(
    graph: &ParamGraph,
    (grid, radius): (CubeSphere, f64),
    air: &crate::moisture::Air,
    layers: &LayerSet,
    cycle: &crate::climate::Cycle,
    level: f64,
) -> crate::moisture::Drains {
    let coarse = CubeSphere::new(YEAR_DRAINAGE_FACE.min(grid.face_size()));
    let height = crate::climate::coarsen(&layers.height(), coarse);
    let wet = Field::from_fn(
        grid,
        |c, _| {
            if layers.liquid.get(c) > 0.0 { 1.0 } else { 0.0 }
        },
    );
    let wet = crate::climate::coarsen(&wet, coarse);
    let cover = crate::climate::coarsen(&layers.regolith, coarse);
    let drains = crate::moisture::Drains::new(
        coarse,
        (
            height.as_slice().to_vec(),
            wet.as_slice().iter().map(|&w| w >= 0.5).collect(),
        ),
        cycle.layout(),
        (radius, air.vapour.liquid_density),
    );
    let places = cycle.sites.len();
    let steps = cycle.surfaces[0].len().checked_div(places).unwrap_or(0);
    let land: Vec<f64> = (0..places)
        .map(|p| {
            let year = &cycle.surfaces[0][p * steps..(p + 1) * steps];
            year.iter().map(|&t| t as f64).sum::<f64>() / steps.max(1) as f64
        })
        .collect();
    let profile = crate::ground::profiles(
        cover.as_slice(),
        &|k| land.get(drains.place_of(k)).copied().unwrap_or(f64::NAN),
        graph.str("materials.bedrock"),
        (air.vapour, air.gravity),
    );
    drains.with_ground(profile, level)
}

/// The year's water on `layers`' ground on `grid` (F437), its liquid standing at `level`, m, as
/// the finish's year gives it: the climate's year, again at its own wind and moisture (F269),
/// its moisture run down the ground and through it (F434, F435) over the ground's shape (F422);
/// not under its own clouds, nor over caps the loop does not lay. None where the body has no
/// working liquid in its air, or its land sheds none of it.
pub fn year(graph: &ParamGraph, grid: CubeSphere, layers: &LayerSet, level: f64) -> Option<Water> {
    let started = std::time::Instant::now();
    let radius = graph.f64("body.radius");
    let balance = crate::climate::EnergyBalance::from_graph(graph);
    let ocean = crate::climate::Ocean::from_graph(graph);
    let share = if ocean.is_some() {
        OCEAN_TRANSPORT_SHARE
    } else {
        0.0
    };
    let air = crate::moisture::Air::from_graph(graph, balance.transport, share)?;
    let climate = balance.solve(CubeSphere::new(CLIMATE_FACE_SIZE));
    let ground = crate::climate::Ground::from_graph(graph);
    let climate_grid = climate.temperature.grid();
    let share = (ocean.is_some() && level.is_finite()).then(|| {
        let wet = Field::from_fn(
            grid,
            |c, _| {
                if layers.liquid.get(c) > 0.0 { 1.0 } else { 0.0 }
            },
        );
        crate::climate::coarsen(&wet, CubeSphere::new(CLIMATE_FACE_SIZE))
    });
    let sea = ocean.as_ref().zip(share.as_ref());
    let height = layers.height();
    let above = Field::from_fn(grid, |c, _| {
        let h = height.get(c);
        if level.is_finite() {
            (h - level as f32).max(0.0)
        } else {
            h
        }
    });
    let cycle = balance.cycle(&climate, &ground, sea, None);
    let drains = drains_of(graph, (grid, radius), &air, layers, &cycle, level);
    let heights = place_heights(&cycle, &above, climate_grid);
    let first = cycle.moisture_over(&air, &heights, Some(&drains));
    let cycle = match cycle.coupling(first.as_ref()) {
        Some(c) => balance.cycle_with(&climate, &ground, sea, None, Some(&c)),
        None => cycle,
    };
    let drains = drains_of(graph, (grid, radius), &air, layers, &cycle, level);
    let heights = place_heights(&cycle, &above, climate_grid);
    let moisture = cycle.moisture_over(&air, &heights, Some(&drains))?;
    let shape = crate::shares::over_cells(&grid, above.as_slice());
    let shape = Field::from_fn(grid, |c, _| shape[grid.index(c)]);
    let fields = moisture.fields(&shape, radius);
    let [mut shed, gain] = moisture.land_water(&fields[0]);
    for (s, &i) in shed.as_mut_slice().iter_mut().zip(layers.ice.as_slice()) {
        if i > 0.0 {
            *s = 0.0;
        }
    }
    // A land that sheds nothing has no rivers to wear it.
    if !shed.as_slice().iter().any(|&s| s > 0.0) {
        return None;
    }
    let waves = cycle.wind.as_ref().and_then(|w| {
        let sea: Vec<bool> = layers.liquid.as_slice().iter().map(|&l| l > 0.0).collect();
        crate::waves::raise(
            grid,
            &sea,
            (
                radius,
                graph.f64("body.gravity"),
                graph.f64("rotation.angular_velocity").abs(),
            ),
            (w.layout.0, w.layout.1, w.steps),
            (&w.east, &w.north, &w.speed),
            &cycle.masses[1],
        )
    });
    if std::env::var_os("EXOGEN_DEBUG").is_some() {
        let s = moisture.summary();
        eprintln!(
            "loop's year: {} falls at {:.0} kg/m² a year over the body; over the land {:.0} falls and {:.0} runs off; {:.1}% of the land under lakes; {:.1} s",
            air.vapour.name,
            s.precipitation,
            s.land_precipitation,
            s.runoff,
            100.0 * s.lakes,
            started.elapsed().as_secs_f64()
        );
    }
    Some(Water {
        shed: shed.as_slice().to_vec(),
        gain: gain.as_slice().to_vec(),
        density: air.vapour.liquid_density,
        waves,
    })
}

/// The rivers' erosion of `layers` on `grid` through `span` years of the loop's history (F437),
/// under `water`, the liquid standing at `level`, m, the ground rebounding as the finish's does
/// (F278), its relief held over the lithosphere's flexural wavelength; what it cut and laid added
/// to `worn`.
pub fn erode(
    graph: &ParamGraph,
    grid: CubeSphere,
    layers: &mut LayerSet,
    water: &Water,
    level: f64,
    span: f64,
    worn: &mut Worn,
) -> Erosion {
    let hold = crate::stages::stress::StressField::flexural_wavelength(
        graph.f64("interior.elastic_thickness"),
        graph.f64("body.gravity"),
    );
    let (crust, mantle) = if crate::stages::ice::is_ice_surface(graph) {
        (DENSITY_CRUST_ICE, DENSITY_ICE)
    } else {
        (DENSITY_CRUST_ROCK, DENSITY_ROCK)
    };
    let started = std::time::Instant::now();
    let sea: Vec<bool> = layers.liquid.as_slice().iter().map(|&l| l > 0.0).collect();
    let before = layers.height();
    let carried = water.waves.as_ref().map(|w| (w, graph.f64("body.gravity")));
    let moved = crate::fluvial::erode_over(
        grid,
        layers,
        &sea,
        (&water.shed, &water.gain),
        (graph.f64("body.radius"), water.density, level),
        (hold, crust, mantle),
        carried,
        (span, (span / LOOP_FLUVIAL_STEP_YR).ceil().max(1.0) as usize),
        Some(&mut worn.held),
    );
    let after = layers.height();
    if worn.lowered.len() != grid.cell_count() {
        worn.lowered = vec![0.0; grid.cell_count()];
    }
    for (k, w) in worn.lowered.iter_mut().enumerate() {
        if !sea[k] {
            *w += before.as_slice()[k] - after.as_slice()[k];
        }
    }
    worn.cut += moved.cut;
    worn.laid += moved.laid;
    worn.to_sea += moved.to_sea;
    worn.floor += moved.floor;
    worn.seconds += started.elapsed().as_secs_f64();
    moved
}

#[cfg(test)]
mod tests {
    use super::*;

    /// F437: the loop's rivers run on every lid but a mobile one, whose ground the loop raises
    /// none of until Phase 6's plates: on Titan's and Mars's, not on the Earth-like's.
    #[test]
    fn the_rivers_wait_for_a_mobile_lid_s_plates() {
        let graph = |name: &str| crate::presets::preset(name).unwrap().to_graph().0;
        assert!(runs(&graph("Titan")));
        assert!(runs(&graph("Mars")));
        assert!(!runs(&graph("Earth")));
    }
}
