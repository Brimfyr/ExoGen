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
//!
//! The stress turns the mixed layer's own flow too (M77): its Ekman transport, (eps tau + f k
//! x tau) / (eps² + f²), at right angles to the stress where the spin's f outruns the layer's
//! friction eps and down the wind where it does not, near the equator (Codron 2012). Where it
//! diverges the mixed layer takes water from beneath it, where it converges it gives water
//! down, and beneath it the water returns the other way (Held 2001): the trades' flow away from
//! the equator, their upwelling on it and their sinking in the subtropics, the subtropical
//! cells. The water beneath stands at the temperature of the surface where it went down,
//! carried back along the return flow, each place's the mean of what flows into it by its
//! weight (Klinger and Marotzke 2000): no constant. The overturning runs only where its
//! surface water flows from a warmer sea to a colder: a sea of one salt stands by its warmth
//! alone, and its water goes down only where it is the colder, to rise colder than the sea it
//! rises under. Where the wind drives its surface water toward the warmer sea, the trades'
//! last degrees into the equator or a locked world's sea toward its star, the warm water
//! gathers and stands, its surface's slope against the wind, and turns over none.

use crate::constants::*;
use crate::wind::Site;
use std::f64::consts::{PI, TAU};

/// The seas' depth-integrated circulation through the year (F272), on the cycle's layout: its
/// streamfunction at each place, m³/s, and the water each place's east and north boundaries
/// pass, m³/s, eastward and northward; the wind's stress on the sea that drives it, N/m² east
/// and north at each place (M76); and the water the mixed layer's Ekman flow passes through
/// each place's east and north boundaries, m³/s, eastward and northward (M77).
#[derive(Clone, Debug)]
pub struct Currents {
    pub layout: (usize, usize),
    pub radius: f64,
    pub psi: Vec<f64>,
    pub east: Vec<f64>,
    pub north: Vec<f64>,
    pub stress: Vec<[f64; 2]>,
    pub ekman_east: Vec<f64>,
    pub ekman_north: Vec<f64>,
}

/// The sea's overturning under its winds (M77): the water each boundary between two places
/// passes in the mixed layer's Ekman flow, m³/s, from its first place to its second, the
/// return flow beneath the same the other way; and the water each place's mixed layer takes
/// from beneath it, where its flow diverges, and gives down, where it converges, m³/s.
#[derive(Clone, Debug, Default)]
pub struct Overturning {
    pub faces: Vec<(usize, usize, f64)>,
    pub up: Vec<f64>,
    pub down: Vec<f64>,
}

impl Currents {
    /// The largest of the gyres' transports, m³/s: the streamfunction's extremes.
    pub fn strongest(&self) -> f64 {
        self.psi.iter().fold(0.0f64, |a, &v| a.max(v.abs()))
    }

    /// The water each boundary between two places `keep` holds passes, m³/s, from its first
    /// place to its second (F272, M77): the gyres' and the mixed layer's Ekman flow's, each
    /// carried upwind on its own; and the overturning, over seas at `surface`, K.
    pub fn flows(
        &self,
        keep: &dyn Fn(usize) -> bool,
        surface: &[f64],
    ) -> (Vec<(usize, usize, f64)>, Overturning) {
        let (n_lat, n_lon) = self.layout;
        let n = n_lat * n_lon;
        let mut faces = Vec::with_capacity(4 * n);
        if self.east.len() == n && self.north.len() == n && n_lon > 0 {
            for p in 0..n {
                let (i, j) = (p / n_lon, p % n_lon);
                let east = i * n_lon + (j + 1) % n_lon;
                for (q, f) in [(east, self.east[p]), (p + n_lon, self.north[p])] {
                    if q < n && q != p && keep(p) && keep(q) && f != 0.0 {
                        faces.push(if f > 0.0 { (p, q, f) } else { (q, p, -f) });
                    }
                }
            }
        }
        let over = self.overturning(keep, surface);
        faces.extend_from_slice(&over.faces);
        (faces, over)
    }

    /// The overturning of the seas' mixed layer between the places `keep` holds (M77), over
    /// seas at `surface`, K: the Ekman flow's boundaries between two of them where its water
    /// flows from the warmer sea to the colder, and what each place's mixed layer takes from
    /// beneath it and gives down.
    pub fn overturning(&self, keep: &dyn Fn(usize) -> bool, surface: &[f64]) -> Overturning {
        let (n_lat, n_lon) = self.layout;
        let n = n_lat * n_lon;
        let mut faces = Vec::with_capacity(2 * n);
        let mut net = vec![0.0; n];
        if self.ekman_east.len() != n
            || self.ekman_north.len() != n
            || surface.len() != n
            || n_lon == 0
        {
            return Overturning {
                faces,
                up: net.clone(),
                down: net,
            };
        }
        for p in 0..n {
            let (i, j) = (p / n_lon, p % n_lon);
            let east = i * n_lon + (j + 1) % n_lon;
            for (q, f) in [(east, self.ekman_east[p]), (p + n_lon, self.ekman_north[p])] {
                if q < n && q != p && f != 0.0 && keep(p) && keep(q) {
                    let (a, b, f) = if f > 0.0 { (p, q, f) } else { (q, p, -f) };
                    if surface[a] <= surface[b] {
                        continue;
                    }
                    faces.push((a, b, f));
                    net[a] += f;
                    net[b] -= f;
                }
            }
        }
        Overturning {
            faces,
            up: net.iter().map(|&x| x.max(0.0)).collect(),
            down: net.iter().map(|&x| (-x).max(0.0)).collect(),
        }
    }
}

/// The temperature of the water beneath each place's mixed layer, K (M77): the steady balance
/// of the layer beneath, each place's the mean of what flows into it by its weight, the water
/// its own mixed layer gives down at `surface`, its temperature, and the return flow's from
/// the places downstream of its mixed layer's flow, so that it stands at the surface's
/// temperature where its water went down, carried back along the return flow (Held 2001;
/// Klinger and Marotzke 2000). Swept from `start`, the last step's, until no place moves by
/// more than OCEAN_RETURN_TOLERANCE; a place nothing flows into keeps its start. Returns the
/// temperatures and whether they settled.
pub fn returning(over: &Overturning, surface: &[f64], start: &[f64]) -> (Vec<f64>, bool) {
    let n = surface.len();
    let mut t = start.to_vec();
    if over.down.len() != n || t.len() != n {
        return (surface.to_vec(), false);
    }
    // The return flow into each place, from the places its mixed layer's flow goes to.
    let mut from: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
    for &(a, b, f) in &over.faces {
        from[a].push((b, f));
    }
    let weight: Vec<f64> = (0..n)
        .map(|p| over.down[p] + from[p].iter().map(|&(_, f)| f).sum::<f64>())
        .collect();
    let mut settled = false;
    for sweep in 0..OCEAN_RETURN_SWEEPS {
        let mut most = 0.0f64;
        let mut visit = |p: usize, t: &mut Vec<f64>| {
            if weight[p] > 0.0 {
                let held =
                    over.down[p] * surface[p] + from[p].iter().map(|&(q, f)| f * t[q]).sum::<f64>();
                let v = held / weight[p];
                most = most.max((v - t[p]).abs());
                t[p] = v;
            }
        };
        if sweep % 2 == 0 {
            (0..n).for_each(|p| visit(p, &mut t));
        } else {
            (0..n).rev().for_each(|p| visit(p, &mut t));
        }
        if most <= OCEAN_RETURN_TOLERANCE {
            settled = true;
            break;
        }
    }
    (t, settled)
}

/// The seas' heat carried through a step of `dt` s (F272, M77): the water each boundary in
/// `faces` passes from its first place to its second, m³/s, the gyres' and the Ekman flow's,
/// brings its upwind place's mixed layer's heat to the place downwind; the water each place
/// takes from beneath it, `up`, m³/s, its return water's at `deep`, K; in steps that move no
/// more than OCEAN_CFL of a place's `water`, m³, each place's changing by what flows into it
/// over what it holds. A place `held` keeps its temperature, a sea under its ice at its
/// freezing point, and gathers the heat brought it instead, m³ K. Returns the places' new
/// temperatures and the heat the held places gathered.
pub fn advect(
    faces: &[(usize, usize, f64)],
    (up, deep): (&[f64], &[f64]),
    water: &[f64],
    held: &[bool],
    temps: &[f64],
    dt: f64,
) -> (Vec<f64>, Vec<f64>) {
    let n = temps.len();
    let mut temps = temps.to_vec();
    let mut brought = vec![0.0; n];
    let mut taken = up.to_vec();
    let mut given = vec![0.0; n];
    for &(from, to, f) in faces {
        given[from] += f;
        taken[to] += f;
    }
    let fastest = (0..n)
        .filter(|&p| water[p] > 0.0)
        .map(|p| given[p].max(taken[p]) / water[p])
        .fold(0.0, f64::max);
    let subs = ((dt * fastest / OCEAN_CFL).ceil() as usize).clamp(1, OCEAN_MOST_SUBSTEPS);
    let h = dt / subs as f64;
    let mut change = vec![0.0; n];
    for _ in 0..subs {
        change.fill(0.0);
        for &(from, to, f) in faces {
            change[to] += f * h * (temps[from] - temps[to]);
        }
        for p in 0..n {
            change[p] += up[p] * h * (deep[p] - temps[p]);
        }
        for p in 0..n {
            if held[p] {
                brought[p] += change[p];
            } else if water[p] > 0.0 {
                temps[p] += change[p] / water[p];
            }
        }
    }
    (temps, brought)
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
    let (ekman_east, ekman_north) = ekman(sites, (n_lat, n_lon), stress, (omega, a, density));
    Some(Currents {
        layout: (n_lat, n_lon),
        radius: a,
        psi,
        east,
        north,
        stress: stress.to_vec(),
        ekman_east,
        ekman_north,
    })
}

/// The water the mixed layer's Ekman flow passes through each place's east and north
/// boundaries, m³/s, eastward and northward (M77; Codron 2012): each place's, kg/m/s, (eps
/// tau_x + f tau_y, eps tau_y - f tau_x) / (eps² + f²) under the stress `stress`, N/m², eps
/// OCEAN_EKMAN_FRICTION, the two places' mean across their boundary's length over the sea's
/// `density`, by the lesser of their seas' shares.
pub fn ekman(
    sites: &[Site],
    (n_lat, n_lon): (usize, usize),
    stress: &[[f64; 2]],
    (omega, a, density): (f64, f64, f64),
) -> (Vec<f64>, Vec<f64>) {
    let n = n_lat * n_lon;
    let (mut east, mut north) = (vec![0.0; n], vec![0.0; n]);
    if n_lat < 2 || n_lon == 0 || sites.len() != n || stress.len() != n || density <= 0.0 {
        return (east, north);
    }
    let (dphi, dlon) = (PI / (n_lat - 1) as f64, TAU / n_lon as f64);
    let lat = |i: usize| -PI / 2.0 + dphi * i as f64;
    let eps = OCEAN_EKMAN_FRICTION;
    let flow: Vec<[f64; 2]> = (0..n)
        .map(|p| {
            let f = 2.0 * omega * lat(p / n_lon).sin();
            let [tx, ty] = stress[p];
            let den = eps * eps + f * f;
            [(eps * tx + f * ty) / den, (eps * ty - f * tx) / den]
        })
        .collect();
    let share = |p: usize, q: usize| {
        sites[p]
            .sea
            .clamp(0.0, 1.0)
            .min(sites[q].sea.clamp(0.0, 1.0))
    };
    for p in 0..n {
        let (i, j) = (p / n_lon, p % n_lon);
        let (lo, hi) = (
            (lat(i) - dphi / 2.0).max(-PI / 2.0),
            (lat(i) + dphi / 2.0).min(PI / 2.0),
        );
        let e = i * n_lon + (j + 1) % n_lon;
        if e != p {
            east[p] = 0.5 * (flow[p][0] + flow[e][0]) * a * (hi - lo) * share(p, e) / density;
        }
        if i + 1 < n_lat {
            let q = p + n_lon;
            north[p] =
                0.5 * (flow[p][1] + flow[q][1]) * a * hi.cos() * dlon * share(p, q) / density;
        }
    }
    (east, north)
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

    /// M77: under trades blowing from the east at 0.05 N/m² everywhere, Earth's spin turns the
    /// mixed layer's flow to the right of the wind in the north and to the left in the south,
    /// away from the equator on both sides, tau / f where f outruns the layer's friction, and
    /// on the equator it runs down the wind, tau / eps: by hand at every boundary. It diverges
    /// on the equator, the mixed layer taking water from beneath it there and giving it down
    /// nowhere near it; and the water it gives down is the water it takes, over the body.
    #[test]
    fn the_sea_s_ekman_flow_turns_from_its_wind() {
        let (n_lat, n_lon) = (45, 8);
        let (dphi, dlon) = (PI / (n_lat - 1) as f64, TAU / n_lon as f64);
        let lat = |i: usize| -PI / 2.0 + dphi * i as f64;
        let sites: Vec<Site> = (0..n_lat * n_lon)
            .map(|k| Site {
                lat: lat(k / n_lon),
                lon: dlon * (k % n_lon) as f64,
                weight: 1.0,
                sea: 1.0,
                ice: 0.0,
            })
            .collect();
        let (omega, a, density) = (7.292e-5, 6.371e6, 1_025.0);
        let tau = -0.05;
        let stress = vec![[tau, 0.0]; n_lat * n_lon];
        let (east, north) = ekman(&sites, (n_lat, n_lon), &stress, (omega, a, density));
        let eps = OCEAN_EKMAN_FRICTION;
        for i in 0..n_lat - 1 {
            let p = i * n_lon;
            let my = |i: usize| {
                let f = 2.0 * omega * lat(i).sin();
                -f * tau / (eps * eps + f * f)
            };
            let hi = lat(i) + dphi / 2.0;
            let by_hand = 0.5 * (my(i) + my(i + 1)) * a * hi.cos() * dlon / density;
            assert!(
                (north[p] - by_hand).abs() <= 1e-9 * by_hand.abs().max(1.0),
                "{i}"
            );
            // Away from the equator on both sides of it.
            let away = north[p] * hi.signum();
            assert!(hi.abs() < 1e-9 || away > 0.0, "{i}: {}", north[p]);
        }
        // On the equator, down the wind.
        let e = (n_lat / 2) * n_lon;
        let width = dphi;
        let by_hand = tau / eps * a * width / density;
        assert!(
            (east[e] / by_hand - 1.0).abs() < 1e-9,
            "{} {by_hand}",
            east[e]
        );
        let at_10 = north[((10f64.to_radians() + PI / 2.0) / dphi).floor() as usize * n_lon];
        println!(
            "the mixed layer carries {:.1} Sv across 10 degrees north, {:.1} on the equator westward",
            at_10 * n_lon as f64 / 1e6,
            east[e] * n_lon as f64 / 1e6
        );
        let currents = Currents {
            layout: (n_lat, n_lon),
            radius: a,
            psi: vec![0.0; n_lat * n_lon],
            east: vec![0.0; n_lat * n_lon],
            north: vec![0.0; n_lat * n_lon],
            stress,
            ekman_east: east,
            ekman_north: north,
        };
        // Over a sea warm at its equator, its water flowing away from it to the colder.
        let warm: Vec<f64> = (0..n_lat * n_lon)
            .map(|k| 300.0 - 40.0 * lat(k / n_lon).sin().powi(2))
            .collect();
        let over = currents.overturning(&|_| true, &warm);
        let row = |i: usize, v: &[f64]| (0..n_lon).map(|j| v[i * n_lon + j]).sum::<f64>();
        let (up, down): (f64, f64) = (over.up.iter().sum(), over.down.iter().sum());
        println!(
            "it takes {:.1} Sv from beneath on its equator, {up:.3e} m³/s over the body",
            row(n_lat / 2, &over.up) / 1e6
        );
        assert!(row(n_lat / 2, &over.up) > 0.0 && row(n_lat / 2, &over.down) == 0.0);
        assert!((up - down).abs() <= 1e-9 * up, "{up} {down}");
        // Its return water rises colder than the sea it rises under, everywhere.
        let (deep, settled) = returning(&over, &warm, &warm);
        assert!(settled);
        for p in 0..n_lat * n_lon {
            assert!(
                over.up[p] == 0.0 || deep[p] < warm[p],
                "{p}: {} {}",
                deep[p],
                warm[p]
            );
        }
        // Over a sea warm at its poles the same wind's water flows toward the warmer sea, and
        // turns over none.
        let cold: Vec<f64> = warm.iter().map(|t| 600.0 - t).collect();
        let none = currents.overturning(&|_| true, &cold);
        assert!(none.faces.iter().all(|&(a, b, _)| cold[a] > cold[b]));
        let equator = (n_lat / 2) * n_lon;
        assert_eq!(none.up[equator], 0.0);
    }

    /// M77: the seas' flows carry the gyres' water and the mixed layer's Ekman flow's, each
    /// upwind on its own, between the places with seas and nowhere else; a sea of one
    /// temperature, its return water at it, stays at it through a step.
    #[test]
    fn the_seas_flows_carry_the_gyres_and_the_overturning() {
        let (n_lat, n_lon) = (5, 4);
        let n = n_lat * n_lon;
        let currents = Currents {
            layout: (n_lat, n_lon),
            radius: 1.0,
            psi: vec![0.0; n],
            east: (0..n)
                .map(|p| if p % 2 == 0 { 3.0 } else { -2.0 })
                .collect(),
            north: (0..n)
                .map(|p| if p < n - n_lon { 1.5 } else { 0.0 })
                .collect(),
            stress: vec![[0.0; 2]; n],
            ekman_east: vec![0.5; n],
            ekman_north: (0..n)
                .map(|p| if p < n - n_lon { -0.75 } else { 0.0 })
                .collect(),
        };
        // Every place a sea but one.
        let keep = |p: usize| p != 7;
        let surface: Vec<f64> = (0..n).map(|p| 300.0 - p as f64).collect();
        let (faces, over) = currents.flows(&keep, &surface);
        let gyres = (0..n)
            .flat_map(|p| {
                let (i, j) = (p / n_lon, p % n_lon);
                [
                    (i * n_lon + (j + 1) % n_lon, currents.east[p]),
                    (p + n_lon, currents.north[p]),
                ]
                .into_iter()
                .filter(move |&(q, f)| q < n && q != p && keep(p) && keep(q) && f != 0.0)
                .map(move |(q, f)| (p, q, f))
            })
            .count();
        assert!(!over.faces.is_empty());
        assert_eq!(faces.len(), gyres + over.faces.len());
        for &(a, b, f) in &faces {
            assert!(f > 0.0 && a != 7 && b != 7, "{a} {b} {f}");
        }
        for face in &over.faces {
            assert!(faces.contains(face), "{face:?}");
        }
        let water = vec![1e3; n];
        let held = vec![false; n];
        let temps = vec![290.0; n];
        let (after, _) = advect(&faces, (&over.up, &temps), &water, &held, &temps, 50.0);
        assert!(after.iter().all(|&t| (t - 290.0).abs() < 1e-9), "{after:?}");
    }

    /// M77: the water beneath stands at the surface's temperature where it went down, carried
    /// back along the return flow (Held 2001): along a strip whose mixed layer flows from its
    /// first place to its last, giving down 4, 4 and 2 m³/s at 296, 294 and 292 K, the return
    /// water that rises at the first place stands at their mean by their weights, 294.4 K; each
    /// place's by hand. The heat its flows carry in a short step is the heat they take from the
    /// water beneath and give it, nothing made or lost: where the warm surface water goes down
    /// and the cooler comes up, the first place cools and the rest warm by the same heat.
    #[test]
    fn the_return_water_is_the_water_that_went_down() {
        let over = Overturning {
            faces: vec![(0, 1, 10.0), (1, 2, 10.0), (2, 3, 6.0), (3, 4, 2.0)],
            up: vec![10.0, 0.0, 0.0, 0.0, 0.0],
            down: vec![0.0, 0.0, 4.0, 4.0, 2.0],
        };
        let surface = [300.0, 298.0, 296.0, 294.0, 292.0];
        let (deep, settled) = returning(&over, &surface, &surface);
        assert!(settled);
        let d3 = (4.0 * 294.0 + 2.0 * 292.0) / 6.0;
        let d2 = (4.0 * 296.0 + 6.0 * d3) / 10.0;
        let by_hand = [d2, d2, d2, d3, 292.0];
        for p in 0..5 {
            assert!(
                (deep[p] - by_hand[p]).abs() < 1e-5,
                "{p}: {} {}",
                deep[p],
                by_hand[p]
            );
        }
        assert!((deep[0] - 294.4).abs() < 1e-5, "{}", deep[0]);
        // Its heat through a short step: what the mixed layers gain, the water beneath gives.
        let water = [1e6; 5];
        let held = [false; 5];
        let dt = 100.0;
        let (after, brought) = advect(&over.faces, (&over.up, &deep), &water, &held, &surface, dt);
        let gained: f64 = (0..5).map(|p| water[p] * (after[p] - surface[p])).sum();
        let beneath: f64 = (0..5)
            .map(|p| over.up[p] * deep[p] - over.down[p] * surface[p])
            .sum::<f64>()
            * dt;
        println!("its mixed layers gain {gained:.4e} m³ K, the water beneath gives {beneath:.4e}");
        assert!(beneath.abs() < 1e-6 * dt, "{beneath}");
        assert!(gained.abs() < 1e-3 * dt * 10.0, "{gained}");
        assert!(brought.iter().all(|&b| b == 0.0));
        assert!(after[0] < surface[0] && after[1] > surface[1], "{after:?}");
        // A place held at its freezing point gathers the heat brought it instead.
        let held = [true, false, false, false, false];
        let (after, brought) = advect(&over.faces, (&over.up, &deep), &water, &held, &surface, dt);
        assert_eq!(after[0], surface[0]);
        let by_hand = over.up[0] * dt * (deep[0] - surface[0]);
        assert!(
            (brought[0] - by_hand).abs() < 1e-9 * by_hand.abs(),
            "{} {by_hand}",
            brought[0]
        );
    }

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
