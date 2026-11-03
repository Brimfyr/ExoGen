// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! The seas' currents (F272).
//!
//! The wind's stress on the sea, rho_a C_D |V| V, drives the seas' depth-integrated
//! circulation, its streamfunction ψ in the balance Stommel (1948) found for a basin: the
//! stress's curl against the bottom's friction and the change of the planet's spin with
//! latitude, r ∇²ψ + β ∂ψ/∂x = curl τ / ρ, the gyres' return flow pressed against their basins'
//! western coasts within some r / β of them; ψ held at none on land by a friction that grows
//! as the sea's share of a place falls, so that the coasts need no line drawn. On a body
//! spinning slowly, β small, the gyres turn the more evenly, by degrees. The flow across each
//! boundary between places, the difference of ψ along it, carries the seas' heat.

use crate::constants::*;
use crate::wind::Site;
use std::f64::consts::{PI, TAU};

/// The seas' depth-integrated circulation through the year (F272), on the cycle's layout: its
/// streamfunction at each place, m³/s, and the water each place's east and north boundaries
/// pass, m³/s, eastward and northward.
#[derive(Clone, Debug)]
pub struct Currents {
    pub layout: (usize, usize),
    pub radius: f64,
    pub psi: Vec<f64>,
    pub east: Vec<f64>,
    pub north: Vec<f64>,
}

impl Currents {
    /// The largest of the gyres' transports, m³/s: the streamfunction's extremes.
    pub fn strongest(&self) -> f64 {
        self.psi.iter().fold(0.0f64, |a, &v| a.max(v.abs()))
    }
}

/// The circulation that the wind's stress `stress`, N/m² east and north at each place, drives
/// in the seas of `sites` on a layout of `n_lat` latitudes and `n_lon` longitudes, of a body
/// of `radius` m spinning at `omega` rad/s, its sea of `density` kg/m³ (Stommel 1948); `links`
/// each place's neighbours and its conductance to each per unit of its area on the unit
/// sphere. None on a layout by latitude alone, which has no basins, or without a sea.
pub fn circulate(
    sites: &[Site],
    (n_lat, n_lon): (usize, usize),
    stress: &[[f64; 2]],
    (omega, radius, density): (f64, f64, f64),
    links: &[Vec<(usize, f64)>],
) -> Option<Currents> {
    let n = n_lat * n_lon;
    if n_lon < 3 || n_lat < 3 || sites.len() != n || stress.len() != n || links.len() != n {
        return None;
    }
    if sites.iter().all(|s| s.sea <= 0.0) {
        return None;
    }
    let a = radius.max(1.0);
    let (dphi, dlon) = (PI / (n_lat - 1) as f64, TAU / n_lon as f64);
    let lat = |i: usize| -PI / 2.0 + dphi * i as f64;
    let at = |i: usize, j: usize| i * n_lon + (j + n_lon) % n_lon;
    // The stress's curl at each place, by its neighbours east and west and north and south.
    let curl: Vec<f64> = (0..n)
        .map(|p| {
            let (i, j) = (p / n_lon, p % n_lon);
            if i == 0 || i + 1 == n_lat {
                return 0.0;
            }
            let c = lat(i).cos();
            let along = (stress[at(i, j + 1)][1] - stress[at(i, j + n_lon - 1)][1]) / (2.0 * dlon);
            let across = (stress[at(i + 1, j)][0] * lat(i + 1).cos()
                - stress[at(i - 1, j)][0] * lat(i - 1).cos())
                / (2.0 * dphi);
            (along - across) / (a * c.max(1e-6))
        })
        .collect();
    // The operator: the bottom's friction on the streamfunction's Laplacian, the spin's change
    // with latitude upwind from the east, where the long waves that set the interior's flow
    // come from, and the friction of the land, (1 - s) / s times the friction's own across a
    // place, as the sea's share s of it falls.
    let r = OCEAN_BOTTOM_FRICTION;
    let beta = 2.0 * omega / (a * a * dlon);
    let most = links
        .iter()
        .flat_map(|l| l.iter().map(|&(_, g)| g))
        .fold(0.0f64, f64::max)
        / (a * a);
    let land: Vec<f64> = sites
        .iter()
        .map(|s| {
            let sea = s.sea.clamp(0.0, 1.0);
            OCEAN_COAST_FRICTION * r * most * (1.0 - sea) / sea.max(1e-6)
        })
        .collect();
    let apply = |x: &[f64], y: &mut [f64]| {
        for p in 0..n {
            let (i, j) = (p / n_lon, p % n_lon);
            let mut v = -land[p] * x[p];
            for &(q, g) in &links[p] {
                v += r * g / (a * a) * (x[q] - x[p]);
            }
            if i > 0 && i + 1 < n_lat {
                v += beta * (x[at(i, j + 1)] - x[p]);
            }
            y[p] = v;
        }
    };
    let rhs: Vec<f64> = (0..n)
        .map(|p| sites[p].sea.clamp(0.0, 1.0) * curl[p] / density)
        .collect();
    let psi = bicgstab(&apply, &rhs, n)?;
    // The water across each boundary, the streamfunction's difference along it: its corners'
    // the mean of the four places about them, the poles' the mean of their rows'.
    let pole = |i: usize| (0..n_lon).map(|j| psi[at(i, j)]).sum::<f64>() / n_lon as f64;
    let (south, north_pole) = (pole(0), pole(n_lat - 1));
    let corner = |i: isize, j: usize| -> f64 {
        match i {
            i if i < 0 => south,
            i if i as usize + 1 >= n_lat => north_pole,
            i => {
                let i = i as usize;
                let (a, b) = (psi[at(i, j)], psi[at(i, j + 1)]);
                0.25 * (a + b + psi[at(i + 1, j)] + psi[at(i + 1, j + 1)])
            }
        }
    };
    let mut east = vec![0.0; n];
    let mut north = vec![0.0; n];
    for p in 0..n {
        let (i, j) = ((p / n_lon) as isize, p % n_lon);
        if (i as usize) + 1 < n_lat {
            north[p] = corner(i, j) - corner(i, j + n_lon - 1);
        }
        east[p] = corner(i - 1, j) - corner(i, j);
    }
    Some(Currents {
        layout: (n_lat, n_lon),
        radius: a,
        psi,
        east,
        north,
    })
}

/// Solves `apply`(x) = `rhs` by biconjugate gradients, stabilised (van der Vorst 1992), from
/// none, to OCEAN_SOLVE_TOLERANCE of the right-hand side.
fn bicgstab(apply: &dyn Fn(&[f64], &mut [f64]), rhs: &[f64], n: usize) -> Option<Vec<f64>> {
    let dot = |a: &[f64], b: &[f64]| a.iter().zip(b).map(|(x, y)| x * y).sum::<f64>();
    let norm = dot(rhs, rhs).sqrt();
    let mut x = vec![0.0; n];
    if norm <= 0.0 {
        return Some(x);
    }
    let mut r = rhs.to_vec();
    let start = r.clone();
    let (mut rho, mut alpha, mut omega) = (1.0, 1.0, 1.0);
    let (mut v, mut p) = (vec![0.0; n], vec![0.0; n]);
    let (mut s, mut t) = (vec![0.0; n], vec![0.0; n]);
    for _ in 0..OCEAN_SOLVE_STEPS {
        let next = dot(&start, &r);
        if next.abs() < 1e-300 {
            break;
        }
        let b = (next / rho) * (alpha / omega);
        rho = next;
        for k in 0..n {
            p[k] = r[k] + b * (p[k] - omega * v[k]);
        }
        apply(&p, &mut v);
        let sv = dot(&start, &v);
        if sv.abs() < 1e-300 {
            break;
        }
        alpha = rho / sv;
        for k in 0..n {
            s[k] = r[k] - alpha * v[k];
        }
        if dot(&s, &s).sqrt() <= OCEAN_SOLVE_TOLERANCE * norm {
            for k in 0..n {
                x[k] += alpha * p[k];
            }
            return Some(x);
        }
        apply(&s, &mut t);
        let tt = dot(&t, &t);
        omega = if tt > 0.0 { dot(&t, &s) / tt } else { 0.0 };
        for k in 0..n {
            x[k] += alpha * p[k] + omega * s[k];
            r[k] = s[k] - omega * t[k];
        }
        if dot(&r, &r).sqrt() <= OCEAN_SOLVE_TOLERANCE * norm {
            return Some(x);
        }
        if omega == 0.0 {
            break;
        }
    }
    x.iter().all(|v| v.is_finite()).then_some(x)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// F272: trades and westerlies over a basin between two meridional coasts turn a gyre
    /// whose return flow runs against its western coast: the water there flows poleward,
    /// faster than anywhere in the interior flows the other way, and the gyre's transport
    /// stands near Sverdrup's, curl τ over ρ β across the basin.
    #[test]
    fn a_basin_s_gyre_presses_against_its_western_coast() {
        let (n_lat, n_lon) = (46, 90);
        let (dphi, dlon) = (PI / (n_lat - 1) as f64, TAU / n_lon as f64);
        let lat = |i: usize| -PI / 2.0 + dphi * i as f64;
        let lon = |j: usize| -PI + dlon * j as f64;
        // A sea from 60° west to 0° east, 10° to 60° north; land elsewhere.
        let inside = |i: usize, j: usize| {
            let (la, lo) = (lat(i).to_degrees(), lon(j).to_degrees());
            (10.0..=60.0).contains(&la) && (-60.0..=0.0).contains(&lo)
        };
        let sites: Vec<Site> = (0..n_lat * n_lon)
            .map(|p| Site {
                lat: lat(p / n_lon),
                lon: lon(p % n_lon),
                weight: 0.0,
                sea: if inside(p / n_lon, p % n_lon) {
                    1.0
                } else {
                    0.0
                },
                ice: 0.0,
            })
            .collect();
        // Easterlies at 10°, westerlies at 45°: τ_x = -0.1 cos(π (φ - 10°) / 35°) N/m².
        let stress: Vec<[f64; 2]> = (0..n_lat * n_lon)
            .map(|p| {
                let la = lat(p / n_lon).to_degrees();
                [-0.1 * (PI * (la - 10.0) / 35.0).cos(), 0.0]
            })
            .collect();
        let links = crate::climate::place_links((n_lat, n_lon));
        let earth = (7.292e-5, 6.371e6, 1_025.0);
        let c = circulate(&sites, (n_lat, n_lon), &stress, earth, &links).expect("currents");
        // Along the row at 30° north: the flow north at each place of the basin.
        let i = (0..n_lat).min_by_key(|&i| ((lat(i).to_degrees() - 30.0).abs() * 100.0) as i64);
        let i = i.unwrap();
        let row: Vec<(f64, f64)> = (0..n_lon)
            .filter(|&j| inside(i, j))
            .map(|j| (lon(j).to_degrees(), c.north[i * n_lon + j]))
            .collect();
        let west = row.first().unwrap().1;
        let interior = row[row.len() / 2].1;
        println!(
            "30°N: the western coast's flow {:.2e} m³/s, the interior's {:.2e}; the gyre {:.1} Sv",
            west,
            interior,
            c.strongest() / 1e6
        );
        assert!(
            west > 0.0 && interior < 0.0 && west > 3.0 * interior.abs(),
            "{row:?}"
        );
        // Sverdrup's transport across the basin at 30°: curl τ / (ρ β) times its width.
        let (f30, a) = (30f64.to_radians(), 6.371e6);
        let beta = 2.0 * 7.292e-5 * f30.cos() / a;
        let curl = 0.1 * PI / 35f64.to_radians() * (PI * 20.0 / 35.0).sin() / a;
        let width = 60f64.to_radians() * a * f30.cos();
        let sverdrup = curl / (1_025.0 * beta) * width;
        println!("Sverdrup's {:.1} Sv", sverdrup / 1e6);
        let strongest = c.strongest();
        assert!(
            strongest > 0.3 * sverdrup && strongest < 3.0 * sverdrup,
            "{strongest} {sverdrup}"
        );
    }
}
