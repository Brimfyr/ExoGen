// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! Axisymmetric flexure of a thin elastic plate over a fluid (F146): the deflection that a
//! load depending only on the distance from a centre gives a lithosphere, and the bending
//! stress at its surface. The plate equation D ∇⁴w + Δρ g w = p in the radius is solved on a
//! uniform grid as the pair ∇²w = u, D ∇²u + Δρ g w = p, a block-tridiagonal system of 2 by
//! 2 blocks that the block Thomas algorithm solves in one sweep each way. The profile is
//! flat at the centre and the plate undisturbed at the grid's end, which lies many flexural
//! parameters out. Against the Hankel transform's solution it agrees to a few parts in ten
//! thousand on a grid a twentieth of the flexural parameter fine.

use crate::constants::{ROCK_POISSON, ROCK_YOUNG_MODULUS};

/// A thin elastic plate over a fluid.
#[derive(Clone, Copy, Debug)]
pub struct Plate {
    /// The flexural rigidity, N m, and the elastic thickness, m.
    pub rigidity: f64,
    pub elastic: f64,
    /// The restoring force per unit deflection, Pa/m: the density of what fills the
    /// deflection times gravity.
    pub restoring: f64,
}

impl Plate {
    /// A rock plate `elastic` m thick over a fluid of density `density`, kg/m³.
    pub fn rock(elastic: f64, density: f64, gravity: f64) -> Self {
        let elastic = elastic.max(1.0);
        Self {
            rigidity: ROCK_YOUNG_MODULUS * elastic.powi(3)
                / (12.0 * (1.0 - ROCK_POISSON * ROCK_POISSON)),
            elastic,
            restoring: density * gravity,
        }
    }

    /// The flexural parameter, m: (4 D / Δρ g)^¼, over which the plate bends.
    pub fn alpha(&self) -> f64 {
        (4.0 * self.rigidity / self.restoring.max(1e-9)).powf(0.25)
    }
}

/// A profile against the distance from a centre, sampled every `step` m from the centre.
#[derive(Clone, Debug, Default)]
pub struct Radial {
    pub step: f64,
    pub values: Vec<f64>,
}

impl Radial {
    /// The profile at `r` m, linear between the samples; the last sample beyond the end.
    pub fn at(&self, r: f64) -> f64 {
        let n = self.values.len();
        if n == 0 || self.step <= 0.0 {
            return 0.0;
        }
        let x = (r.max(0.0) / self.step).min((n - 1) as f64);
        let i = (x.floor() as usize).min(n - 1);
        let j = (i + 1).min(n - 1);
        let f = x - i as f64;
        self.values[i] * (1.0 - f) + self.values[j] * f
    }

    /// How far out the profile reaches, m.
    pub fn reach(&self) -> f64 {
        self.step * self.values.len().saturating_sub(1) as f64
    }
}

/// The coefficients of the radial Laplacian at sample i: (f_{i-1}, f_i, f_{i+1}).
fn laplacian(i: usize, step: f64) -> (f64, f64, f64) {
    let h2 = step * step;
    if i == 0 {
        // Flat at the centre: ∇²f = 2 f'' there, with f_{-1} = f_1.
        (0.0, -4.0 / h2, 4.0 / h2)
    } else {
        let half = 0.5 / i as f64;
        ((1.0 - half) / h2, -2.0 / h2, (1.0 + half) / h2)
    }
}

type Block = [[f64; 2]; 2];

fn mul(a: &Block, b: &Block) -> Block {
    [
        [
            a[0][0] * b[0][0] + a[0][1] * b[1][0],
            a[0][0] * b[0][1] + a[0][1] * b[1][1],
        ],
        [
            a[1][0] * b[0][0] + a[1][1] * b[1][0],
            a[1][0] * b[0][1] + a[1][1] * b[1][1],
        ],
    ]
}

fn apply(a: &Block, v: [f64; 2]) -> [f64; 2] {
    [
        a[0][0] * v[0] + a[0][1] * v[1],
        a[1][0] * v[0] + a[1][1] * v[1],
    ]
}

fn inverse(a: &Block) -> Block {
    let det = a[0][0] * a[1][1] - a[0][1] * a[1][0];
    let det = if det.abs() < 1e-300 { 1e-300 } else { det };
    [
        [a[1][1] / det, -a[0][1] / det],
        [-a[1][0] / det, a[0][0] / det],
    ]
}

/// The deflection, m, up positive, of the plate under a load, Pa, up positive, both
/// sampled every `step` m from the centre; the plate is undisturbed one step past the last
/// sample.
pub fn deflect(plate: &Plate, step: f64, load: &[f64]) -> Vec<f64> {
    let n = load.len();
    if n == 0 {
        return Vec::new();
    }
    let d = plate.rigidity;
    let k = plate.restoring;
    // Unknowns x_i = (w_i, u_i): u_i - L(w)_i = 0 and D L(u)_i + k w_i = p_i.
    let mut cp: Vec<Block> = vec![[[0.0; 2]; 2]; n];
    let mut dp: Vec<[f64; 2]> = vec![[0.0; 2]; n];
    for i in 0..n {
        let (a, b, c) = laplacian(i, step);
        let lower: Block = [[-a, 0.0], [0.0, d * a]];
        let diag: Block = [[-b, 1.0], [k, d * b]];
        let upper: Block = [[-c, 0.0], [0.0, d * c]];
        let rhs = [0.0, load[i]];
        let (m, r) = if i == 0 {
            (diag, rhs)
        } else {
            let lc = mul(&lower, &cp[i - 1]);
            let ld = apply(&lower, dp[i - 1]);
            (
                [
                    [diag[0][0] - lc[0][0], diag[0][1] - lc[0][1]],
                    [diag[1][0] - lc[1][0], diag[1][1] - lc[1][1]],
                ],
                [rhs[0] - ld[0], rhs[1] - ld[1]],
            )
        };
        let mi = inverse(&m);
        cp[i] = mul(&mi, &upper);
        dp[i] = apply(&mi, r);
    }
    let mut x = vec![[0.0f64; 2]; n];
    x[n - 1] = dp[n - 1];
    for i in (0..n - 1).rev() {
        let next = apply(&cp[i], x[i + 1]);
        x[i] = [dp[i][0] - next[0], dp[i][1] - next[1]];
    }
    x.into_iter().map(|v| v[0]).collect()
}

/// The greatest tension at the plate's upper surface along a deflection, Pa, at each sample:
/// E h / (2 (1 - ν²)) times minus the curvature, radial (w'' + ν w'/r) or round the centre
/// (w'/r + ν w''), whichever pulls harder; where the surface bends convex up it is pulled
/// apart.
pub fn surface_tension(plate: &Plate, step: f64, w: &[f64]) -> Vec<f64> {
    let n = w.len();
    let at = |i: isize| -> f64 {
        if i < 0 {
            w[(-i) as usize % n.max(1)]
        } else if (i as usize) < n {
            w[i as usize]
        } else {
            0.0
        }
    };
    let scale = ROCK_YOUNG_MODULUS * plate.elastic / (2.0 * (1.0 - ROCK_POISSON * ROCK_POISSON));
    (0..n)
        .map(|i| {
            let j = i as isize;
            let second = (at(j + 1) - 2.0 * at(j) + at(j - 1)) / (step * step);
            let slope_over_r = if i == 0 {
                second
            } else {
                (at(j + 1) - at(j - 1)) / (2.0 * step) / (i as f64 * step)
            };
            let radial = -(second + ROCK_POISSON * slope_over_r);
            let hoop = -(slope_over_r + ROCK_POISSON * second);
            scale * radial.max(hoop)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A load far wider than the flexural parameter is held up as Airy isostasy holds it,
    /// and one far narrower is spread over the flexural parameter, with a moat and a bulge
    /// beyond it.
    #[test]
    fn a_plate_holds_wide_loads_isostatically_and_spreads_narrow_ones() {
        let plate = Plate::rock(20e3, 3300.0, 8.87);
        let alpha = plate.alpha();
        assert!((40e3..80e3).contains(&alpha), "{alpha}");
        let step = alpha / 20.0;
        let n = 2000;
        let p = 3300.0 * 8.87 * 1000.0;
        // Wide: a disc ten flexural parameters across stands 1 km at its centre.
        let wide: Vec<f64> = (0..n)
            .map(|i| {
                if (i as f64) * step < 10.0 * alpha {
                    p
                } else {
                    0.0
                }
            })
            .collect();
        let w = deflect(&plate, step, &wide);
        assert!((w[0] - 1000.0).abs() < 10.0, "{}", w[0]);
        // Narrow: a disc a tenth of the flexural parameter across sinks the plate far less
        // than isostasy would, and beyond it a moat gives way to a bulge.
        let narrow: Vec<f64> = (0..n)
            .map(|i| {
                if (i as f64) * step < 0.1 * alpha {
                    p
                } else {
                    0.0
                }
            })
            .collect();
        let w = deflect(&plate, step, &narrow);
        assert!(w[0] > 0.0 && w[0] < 50.0, "{}", w[0]);
        let low = w.iter().cloned().fold(f64::INFINITY, f64::min);
        assert!(low < 0.0, "{low}");
        // The bending pulls the surface apart over the load's centre.
        let s = surface_tension(&plate, step, &w);
        assert!(s.iter().cloned().fold(0.0, f64::max) > 0.0);
    }

    /// A Gaussian load on a plate far weaker than it is wide is held up isostatically all
    /// along its profile.
    #[test]
    fn a_weak_plate_follows_its_load() {
        let plate = Plate::rock(1e3, 3300.0, 8.87);
        let step = 2e3;
        let width = 200e3;
        let load: Vec<f64> = (0..1500)
            .map(|i| {
                let r = i as f64 * step;
                3300.0 * 8.87 * 500.0 * (-(r / width).powi(2)).exp()
            })
            .collect();
        let w = deflect(&plate, step, &load);
        for (i, &x) in w.iter().enumerate().step_by(25).take(10) {
            let r = i as f64 * step;
            let airy = 500.0 * (-(r / width).powi(2)).exp();
            assert!((x - airy).abs() < 5.0, "{r}: {x} vs {airy}");
        }
    }
}
