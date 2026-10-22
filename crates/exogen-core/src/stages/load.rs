// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! The volcanic loads on the lid (F157). The lava a flood lays weighs on the lid, which
//! bends under it as a thin elastic plate over the mantle, D ∇⁴w + ρ_m g w = q (Turcotte
//! and Schubert 2002), of the flexural rigidity D = E Te³ / 12 (1 − ν²) the lid's elastic
//! part had when the lava came. Under a load wider than the plate's flexural parameter the
//! lid sinks as Airy isostasy would have it and stands unbent; about the load's margin it
//! bends, its top squeezed just inside and stretched just outside, where a bulge rises
//! beyond the moat: the lunar maria's wrinkle ridges within and the arcuate graben about
//! their margins (Solomon and Head 1980; Freed, Melosh and Solomon 2001). The lid is a shell: a
//! depression shortens it in its own plane as well as bending it, which stiffens it against
//! the load and squeezes the depressed ground in its plane while the ground about it is
//! squeezed outward and stretched round, so that the squeeze within a filled basin is
//! greatest outward and its thrusts run round it, and beyond its margin the lid is pulled
//! both ways and opens graben round it (the shallow shell's equations: D ∇⁴w + E h w / R² +
//! ρ_m g w = q, with the Airy stress function ∇²Φ = −E h w / R; Turcotte et al. 1981;
//! Freed, Melosh and Solomon 2001, whose spherical shells resolve the flat plate's
//! strike-slip paradox). The plate is solved on a cube-sphere grid fine enough for its
//! flexural parameter, by conjugate gradients, and the stresses at its top are read from
//! the curvatures of the deflection and of the stress function in each cell's own tangent
//! plane. Each load bends the lid of its
//! time; the lid that thickens beneath it later bears none of its stress, so the stresses
//! of the loads laid at different times add.

use rayon::prelude::*;

use crate::constants::*;
use crate::fields::det_sum;
use crate::geom::{dot, local_basis};
use crate::grid::CubeSphere;
use crate::stages::stress::FieldHistory;

/// A cell's tangent frame, and the weights that read its curvature (w_uu, w_uv, w_vv) from
/// its eight neighbours' differences from it.
type Frame = ([f64; 3], [f64; 3]);
type Curvature = [[f32; 8]; 3];

/// A thin elastic plate over a fluid on a cube-sphere grid: each cell's area, and the
/// stencils that read the deflection's curvature in each cell's tangent frame, by least
/// squares on a quadratic over the cell and its eight neighbours, exact for a quadratic on
/// every cell, across the cube's edges too. Their trace is the Laplacian G; the plate's
/// bending energy, ½ D ∫ (∇²w)² dA on a surface without edges, gives its operator as
/// D Gᵀ M G + k M with M the areas, symmetric and positive definite.
pub struct SpherePlate {
    pub grid: CubeSphere,
    radius: f64,
    /// Each cell's area, m².
    area: Vec<f64>,
    /// Each cell's tangent frame, its eight neighbours, and the weights that give the
    /// curvature (w_uu, w_uv, w_vv) in that frame from the neighbours' differences from the
    /// cell.
    frame: Vec<Frame>,
    stencil: Vec<[u32; 8]>,
    curvature: Vec<Curvature>,
    /// G's columns, by cell: the rows that read the cell and with what weight, for Gᵀ.
    column_start: Vec<usize>,
    column: Vec<(u32, f64)>,
}

impl SpherePlate {
    pub fn new(face: u32, radius: f64) -> Self {
        let grid = CubeSphere::new(face);
        let n = grid.cell_count();
        let area: Vec<f64> = (0..n)
            .map(|i| grid.cell_area(grid.cell(i)) * radius * radius)
            .collect();
        let centers: Vec<[f64; 3]> = (0..n).map(|i| grid.center(grid.cell(i))).collect();
        let offsets = [
            (-1, -1),
            (0, -1),
            (1, -1),
            (-1, 0),
            (1, 0),
            (-1, 1),
            (0, 1),
            (1, 1),
        ];
        let built: Vec<(Frame, [u32; 8], Curvature)> = (0..n)
            .into_par_iter()
            .map(|i| {
                let c = grid.cell(i);
                let p = centers[i];
                let (u, v) = local_basis(p);
                let mut stencil = [0u32; 8];
                let mut rows = [[0.0f64; 5]; 8];
                for (k, &(di, dj)) in offsets.iter().enumerate() {
                    let j = grid.index(grid.neighbor(c, di, dj));
                    stencil[k] = j as u32;
                    let q = centers[j];
                    let along = dot(q, p).max(1e-6);
                    let x = radius * dot(q, u) / along;
                    let y = radius * dot(q, v) / along;
                    rows[k] = [x, y, 0.5 * x * x, x * y, 0.5 * y * y];
                }
                ((u, v), stencil, curvature_weights(&rows))
            })
            .collect();
        let mut frame = Vec::with_capacity(n);
        let mut stencil = Vec::with_capacity(n);
        let mut curvature = Vec::with_capacity(n);
        for (f, s, w) in built {
            frame.push(f);
            stencil.push(s);
            curvature.push(w);
        }
        // G's rows are the trace's weights on the neighbours and minus their sum on the
        // cell; its columns, gathered once, apply Gᵀ.
        let mut entries: Vec<(u32, u32, f64)> = Vec::with_capacity(9 * n);
        for i in 0..n {
            let mut own = 0.0;
            for k in 0..8 {
                let g = (curvature[i][0][k] + curvature[i][2][k]) as f64;
                entries.push((stencil[i][k], i as u32, g));
                own -= g;
            }
            entries.push((i as u32, i as u32, own));
        }
        entries.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
        let mut column_start = vec![0usize; n + 1];
        for e in &entries {
            column_start[e.0 as usize + 1] += 1;
        }
        for i in 0..n {
            column_start[i + 1] += column_start[i];
        }
        let column = entries.into_iter().map(|(_, r, g)| (r, g)).collect();
        Self {
            grid,
            radius,
            area,
            frame,
            stencil,
            curvature,
            column_start,
            column,
        }
    }

    /// The membrane's stiffness against a deflection, Pa/m, of a shell `elastic` m thick on
    /// this grid's sphere: E h / R².
    pub fn membrane_stiffness(&self, elastic: f64) -> f64 {
        ROCK_YOUNG_MODULUS * elastic / (self.radius * self.radius)
    }

    /// The Laplacian's inverse on the sphere for a field of no mean: x with G x = f, by
    /// BiCGSTAB on G shifted by a hundredth of its least nonzero eigenvalue, which settles
    /// the constant the Laplacian cannot see.
    fn inverse_laplace(&self, f: &[f64]) -> Vec<f64> {
        let n = f.len();
        let shift = 0.02 / (self.radius * self.radius);
        let apply = |x: &[f64], out: &mut [f64]| {
            self.laplace(x, out);
            out.par_iter_mut()
                .enumerate()
                .for_each(|(i, o)| *o -= shift * x[i]);
        };
        let dot = |a: &[f64], b: &[f64]| det_sum(n, |i| a[i] * b[i]);
        let mut x = vec![0.0; n];
        let mut r = f.to_vec();
        let goal = PLATE_SOLVE_TOLERANCE * dot(f, f).sqrt();
        if goal <= 0.0 {
            return x;
        }
        let r_hat = r.clone();
        let (mut rho, mut alpha, mut omega) = (1.0, 1.0, 1.0);
        let mut v = vec![0.0; n];
        let mut p = vec![0.0; n];
        let mut s = vec![0.0; n];
        let mut t = vec![0.0; n];
        for _ in 0..PLATE_SOLVE_MAX_ITERATIONS {
            let rho_next = dot(&r_hat, &r);
            if rho_next.abs() < 1e-300 {
                break;
            }
            let beta = (rho_next / rho) * (alpha / omega);
            p.par_iter_mut()
                .enumerate()
                .for_each(|(i, pi)| *pi = r[i] + beta * (*pi - omega * v[i]));
            apply(&p, &mut v);
            alpha = rho_next / dot(&r_hat, &v);
            s.par_iter_mut()
                .enumerate()
                .for_each(|(i, si)| *si = r[i] - alpha * v[i]);
            if dot(&s, &s).sqrt() <= goal {
                x.par_iter_mut()
                    .enumerate()
                    .for_each(|(i, xi)| *xi += alpha * p[i]);
                break;
            }
            apply(&s, &mut t);
            let tt = dot(&t, &t);
            omega = if tt > 0.0 { dot(&t, &s) / tt } else { 0.0 };
            x.par_iter_mut()
                .enumerate()
                .for_each(|(i, xi)| *xi += alpha * p[i] + omega * s[i]);
            r.par_iter_mut()
                .enumerate()
                .for_each(|(i, ri)| *ri = s[i] - omega * t[i]);
            rho = rho_next;
            if dot(&r, &r).sqrt() <= goal || omega == 0.0 {
                break;
            }
        }
        x
    }

    /// The face of the grid that resolves a plate of this flexural parameter, m, with
    /// PLATE_CELLS_PER_PARAMETER cells across it, within PLATE_FACE_RANGE.
    pub fn face_for(flexural: f64, radius: f64) -> u32 {
        let cell = flexural / PLATE_CELLS_PER_PARAMETER;
        let face = (std::f64::consts::FRAC_PI_2 * radius / cell.max(1.0)).ceil();
        face.clamp(PLATE_FACE_RANGE.0 as f64, PLATE_FACE_RANGE.1 as f64) as u32
    }

    /// The Laplacian of a field at each cell, G x: the trace of its fitted curvature.
    fn laplace(&self, x: &[f64], out: &mut [f64]) {
        out.par_iter_mut().enumerate().for_each(|(i, o)| {
            let xi = x[i];
            let mut sum = 0.0;
            for k in 0..8 {
                let g = (self.curvature[i][0][k] + self.curvature[i][2][k]) as f64;
                sum += g * (x[self.stencil[i][k] as usize] - xi);
            }
            *o = sum;
        });
    }

    /// Gᵀ y at each cell, gathered over G's column.
    fn laplace_transposed(&self, y: &[f64], out: &mut [f64]) {
        out.par_iter_mut().enumerate().for_each(|(j, o)| {
            *o = self.column[self.column_start[j]..self.column_start[j + 1]]
                .iter()
                .map(|&(r, g)| g * y[r as usize])
                .sum();
        });
    }

    /// The plate's operator: D Gᵀ M G x + k M x, its restoring force per cell.
    fn apply(&self, rigidity: f64, restoring: &[f64], x: &[f64], tmp: &mut [f64], out: &mut [f64]) {
        self.laplace(x, tmp);
        tmp.par_iter_mut()
            .enumerate()
            .for_each(|(i, t)| *t *= self.area[i]);
        self.laplace_transposed(tmp, out);
        out.par_iter_mut()
            .enumerate()
            .for_each(|(i, o)| *o = rigidity * *o + restoring[i] * self.area[i] * x[i]);
    }

    /// The deflection, m, down positive, of a plate of this rigidity, N m, over a fluid
    /// whose restoring force is `restoring`, Pa/m, under a load, Pa, down positive, per
    /// cell: D ∇⁴w + k w = q, by conjugate gradients preconditioned by the operator's
    /// diagonal, from the Airy deflection q / k.
    pub fn deflect(&self, rigidity: f64, restoring: f64, load: &[f64]) -> Vec<f64> {
        self.deflect_varying(rigidity, &vec![restoring; self.area.len()], load)
    }

    /// The deflection under a restoring force that varies from cell to cell, Pa/m (F166):
    /// weaker where the sag fills with something heavier than air.
    pub fn deflect_varying(&self, rigidity: f64, restoring: &[f64], load: &[f64]) -> Vec<f64> {
        let n = self.area.len();
        let b: Vec<f64> = (0..n).map(|i| self.area[i] * load[i]).collect();
        let diagonal: Vec<f64> = (0..n)
            .map(|j| {
                let bending: f64 = self.column[self.column_start[j]..self.column_start[j + 1]]
                    .iter()
                    .map(|&(r, g)| g * g * self.area[r as usize])
                    .sum();
                rigidity * bending + restoring[j] * self.area[j]
            })
            .collect();
        let mut x: Vec<f64> = load.iter().zip(restoring).map(|(q, k)| q / k).collect();
        let mut tmp = vec![0.0; n];
        let mut ap = vec![0.0; n];
        self.apply(rigidity, restoring, &x, &mut tmp, &mut ap);
        let mut r: Vec<f64> = (0..n).map(|i| b[i] - ap[i]).collect();
        let mut z: Vec<f64> = (0..n).map(|i| r[i] / diagonal[i]).collect();
        let mut p = z.clone();
        let mut rz = det_sum(n, |i| r[i] * z[i]);
        let goal = PLATE_SOLVE_TOLERANCE * det_sum(n, |i| b[i] * b[i]).sqrt();
        for _ in 0..PLATE_SOLVE_MAX_ITERATIONS {
            if det_sum(n, |i| r[i] * r[i]).sqrt() <= goal {
                break;
            }
            self.apply(rigidity, restoring, &p, &mut tmp, &mut ap);
            let pap = det_sum(n, |i| p[i] * ap[i]);
            if pap <= 0.0 {
                break;
            }
            let alpha = rz / pap;
            x.par_iter_mut()
                .zip(r.par_iter_mut())
                .enumerate()
                .for_each(|(i, (xi, ri))| {
                    *xi += alpha * p[i];
                    *ri -= alpha * ap[i];
                });
            z.par_iter_mut()
                .enumerate()
                .for_each(|(i, zi)| *zi = r[i] / diagonal[i]);
            let rz_next = det_sum(n, |i| r[i] * z[i]);
            let beta = rz_next / rz;
            rz = rz_next;
            p.par_iter_mut()
                .enumerate()
                .for_each(|(i, pi)| *pi = z[i] + beta * *pi);
        }
        x
    }

    /// A field's curvature (f_uu, f_uv, f_vv) at a cell, in its tangent frame.
    fn curvature_at(&self, f: &[f64], i: usize) -> [f64; 3] {
        let mut c = [0.0f64; 3];
        for (k, &j) in self.stencil[i].iter().enumerate() {
            let d = f[j as usize] - f[i];
            for (m, cm) in c.iter_mut().enumerate() {
                *cm += self.curvature[i][m][k] as f64 * d;
            }
        }
        c
    }

    /// The stress at the top of a plate `elastic` m thick bent by a deflection, Pa, tension
    /// positive: E Te / 2 (1 − ν²) times the curvature, (w_uu + ν w_vv) across one way,
    /// (w_vv + ν w_uu) the other and (1 − ν) w_uv between, where the deflection is down
    /// positive, so that the top of a bowl is squeezed and of a bulge stretched; and, for a
    /// shell (`membrane`), its membrane's stress through its thickness, from the Airy stress
    /// function Φ with ∇²Φ = −E h w / R, σ_uu = Φ_vv / h, σ_vv = Φ_uu / h, σ_uv = −Φ_uv / h,
    /// the deflection's mean left out, which an incompressible interior lets no lid make; as
    /// tensors in space.
    pub fn surface_stress(&self, w: &[f64], elastic: f64, membrane: bool) -> LoadStress {
        let s = ROCK_YOUNG_MODULUS * elastic / (2.0 * (1.0 - ROCK_POISSON * ROCK_POISSON));
        let nu = ROCK_POISSON;
        let n = self.area.len();
        let phi = if membrane {
            let total: f64 = det_sum(n, |i| self.area[i]);
            let mean = det_sum(n, |i| self.area[i] * w[i]) / total;
            let scale = -ROCK_YOUNG_MODULUS * elastic / self.radius;
            let f: Vec<f64> = (0..n).map(|i| scale * (w[i] - mean)).collect();
            Some(self.inverse_laplace(&f))
        } else {
            None
        };
        let tensor = (0..n)
            .into_par_iter()
            .map(|i| {
                let c = self.curvature_at(w, i);
                let (mut suu, mut svv, mut suv) = (
                    s * (c[0] + nu * c[2]),
                    s * (c[2] + nu * c[0]),
                    s * (1.0 - nu) * c[1],
                );
                if let Some(phi) = &phi {
                    let m = self.curvature_at(phi, i);
                    suu += m[2] / elastic;
                    svv += m[0] / elastic;
                    suv -= m[1] / elastic;
                }
                let (u, v) = self.frame[i];
                let t = |a: usize, b: usize| {
                    suu * u[a] * u[b] + svv * v[a] * v[b] + suv * (u[a] * v[b] + v[a] * u[b])
                };
                [
                    t(0, 0) as f32,
                    t(1, 1) as f32,
                    t(2, 2) as f32,
                    t(0, 1) as f32,
                    t(0, 2) as f32,
                    t(1, 2) as f32,
                ]
            })
            .collect();
        LoadStress {
            grid: self.grid,
            tensor,
        }
    }
}

/// The weights that read (w_uu, w_uv, w_vv) from eight neighbours' differences from the
/// cell, by least squares on w_j − w_i = a x + b y + ½ A x² + B x y + ½ C y².
fn curvature_weights(rows: &[[f64; 5]; 8]) -> Curvature {
    let mut normal = [[0.0f64; 5]; 5];
    for r in rows {
        for a in 0..5 {
            for b in 0..5 {
                normal[a][b] += r[a] * r[b];
            }
        }
    }
    let inverse = invert5(normal);
    let mut out = [[0.0f32; 8]; 3];
    for (k, r) in rows.iter().enumerate() {
        for m in 0..3 {
            let row = &inverse[m + 2];
            let w: f64 = (0..5).map(|a| row[a] * r[a]).sum();
            out[m][k] = w as f32;
        }
    }
    out
}

/// The inverse of a symmetric positive definite 5×5 matrix by Gauss–Jordan elimination
/// with partial pivoting.
fn invert5(m: [[f64; 5]; 5]) -> [[f64; 5]; 5] {
    let mut a = m;
    let mut inv = [[0.0f64; 5]; 5];
    for (i, row) in inv.iter_mut().enumerate() {
        row[i] = 1.0;
    }
    for col in 0..5 {
        let pivot = (col..5)
            .max_by(|&x, &y| a[x][col].abs().total_cmp(&a[y][col].abs()))
            .unwrap_or(col);
        a.swap(col, pivot);
        inv.swap(col, pivot);
        let d = a[col][col];
        let d = if d.abs() < 1e-300 { 1e-300 } else { d };
        for k in 0..5 {
            a[col][k] /= d;
            inv[col][k] /= d;
        }
        for r in 0..5 {
            if r != col {
                let f = a[r][col];
                if f != 0.0 {
                    for k in 0..5 {
                        a[r][k] -= f * a[col][k];
                        inv[r][k] -= f * inv[col][k];
                    }
                }
            }
        }
    }
    inv
}

/// The stress a load put on the lid when it was laid (F157): at each cell of the plate's
/// grid, at the lid's top, Pa, tension positive, as a tensor in space, xx, yy, zz, xy, xz,
/// yz.
#[derive(Clone, Debug)]
pub struct LoadStress {
    pub grid: CubeSphere,
    pub tensor: Vec<[f32; 6]>,
}

impl LoadStress {
    /// Another load's stress added, on the same grid (F174).
    pub fn add(&mut self, other: &LoadStress) {
        assert_eq!(self.tensor.len(), other.tensor.len());
        for (a, b) in self.tensor.iter_mut().zip(&other.tensor) {
            for (x, y) in a.iter_mut().zip(b) {
                *x += y;
            }
        }
    }

    /// The stress at a point in the frame (u, v) of its tangent plane: (σ_uu, σ_vv, σ_uv),
    /// from the cell that holds it.
    pub fn at(&self, p: [f64; 3], u: [f64; 3], v: [f64; 3]) -> (f64, f64, f64) {
        let t = self.tensor[self.grid.index(self.grid.locate(p).0)];
        let t = t.map(|x| x as f64);
        let apply = |a: [f64; 3]| {
            [
                t[0] * a[0] + t[3] * a[1] + t[4] * a[2],
                t[3] * a[0] + t[1] * a[1] + t[5] * a[2],
                t[4] * a[0] + t[5] * a[1] + t[2] * a[2],
            ]
        };
        let (tu, tv) = (apply(u), apply(v));
        (dot(u, tu), dot(v, tv), dot(u, tv))
    }
}

/// The flexural parameter of a rock lid `elastic` m thick over the mantle, m:
/// (D / ρ_m g)^¼.
pub fn flexural_parameter(elastic: f64, gravity: f64) -> f64 {
    let rigidity = ROCK_YOUNG_MODULUS * elastic.max(1.0).powi(3)
        / (12.0 * (1.0 - ROCK_POISSON * ROCK_POISSON));
    (rigidity / (DENSITY_ROCK * gravity.max(1e-3))).powf(0.25)
}

/// A load laid on the lid (F157, F162): the stress it put on the lid and the lid's
/// deflection under it, m, down positive, on the grid it was solved on.
pub struct Laid {
    pub stress: LoadStress,
    pub grid: CubeSphere,
    pub deflection: Vec<f32>,
}

impl Laid {
    /// The deflection on another grid, bilinear.
    pub fn deflection_on(&self, target: CubeSphere) -> Vec<f32> {
        let mut field = crate::fields::Field::new(self.grid);
        field.as_mut_slice().copy_from_slice(&self.deflection);
        field.resample(target).into_vec()
    }
}

/// How far the ground sinks under a deflection, m (F162): as far as the lid deflects, down
/// positive, but beneath an edifice its lava fills the sag as far as it is thick (`held`,
/// m), so that its surface stays as it built it while the ground about it sinks into a moat;
/// the bulge beyond rises.
pub fn sink(deflection: &[f32], held: Option<&[f32]>) -> Vec<f32> {
    deflection
        .iter()
        .enumerate()
        .map(|(k, &w)| {
            let h = held.map_or(0.0, |s| s[k]);
            if w > 0.0 { w - w.min(h) } else { w }
        })
        .collect()
}

/// The stress the lava laid on a grid puts on a rock lid `elastic` m thick (F157), and the
/// lid's deflection under it (F162): `thickness`, m per cell of `grid`, of lava of
/// `density`, kg/m³, averaged onto a plate grid that resolves the lid's flexural
/// parameter, bent over the mantle. Beneath the edifices standing `held` m high the lid
/// sinks into their own lava, which fills the sag as far as they are thick, so that the
/// mantle pushes back there by the difference between the rock's density and the lava's
/// (F166). None where no lava was laid.
pub fn lay(
    grid: CubeSphere,
    thickness: &[f32],
    body_radius: f64,
    gravity: f64,
    elastic: f64,
    density: f64,
    held: Option<&[f32]>,
) -> Option<Laid> {
    if !thickness.iter().any(|&t| t > 0.0) {
        return None;
    }
    let elastic = elastic.max(1e3);
    let flexural = flexural_parameter(elastic, gravity);
    let plate = SpherePlate::new(SpherePlate::face_for(flexural, body_radius), body_radius);
    let target = plate.grid;
    // The load on the plate's grid: each source cell's weight, by area, into the plate cell
    // that holds its centre.
    let mut weight = vec![0.0f64; target.cell_count()];
    let mut covered = vec![0.0f64; target.cell_count()];
    let mut standing = vec![0.0f64; target.cell_count()];
    for (k, &t) in thickness.iter().enumerate() {
        let c = grid.cell(k);
        let a = grid.cell_area(c);
        let j = target.index(target.locate(grid.center(c)).0);
        weight[j] += a * t as f64;
        covered[j] += a;
        standing[j] += a * held.map_or(0.0, |h| h[k] as f64);
    }
    let edifice: Vec<f64> = standing
        .iter()
        .zip(&covered)
        .map(|(&s, &a)| if a > 0.0 { s / a } else { 0.0 })
        .collect();
    let load: Vec<f64> = weight
        .iter()
        .zip(&covered)
        .map(|(&w, &a)| {
            if a > 0.0 {
                density * gravity * w / a
            } else {
                0.0
            }
        })
        .collect();
    let rigidity =
        ROCK_YOUNG_MODULUS * elastic.powi(3) / (12.0 * (1.0 - ROCK_POISSON * ROCK_POISSON));
    // The lid is a shell: its membrane stiffens it against the load.
    let membrane = plate.membrane_stiffness(elastic);
    let open = DENSITY_ROCK * gravity + membrane;
    // Beneath an edifice its lava fills the sag as far as it is thick (F166): the share of
    // the sag it fills, read from the deflection it is filling, once more from the first.
    let filled = |w: &[f64]| -> Vec<f64> {
        (0..w.len())
            .map(|j| {
                let fill = if edifice[j] <= 0.0 {
                    0.0
                } else if w[j] > 0.0 {
                    (edifice[j] / w[j]).min(1.0)
                } else {
                    1.0
                };
                ((DENSITY_ROCK - density * fill).max(ROOT_MIN_CONTRAST)) * gravity + membrane
            })
            .collect()
    };
    let w = if edifice.iter().any(|&h| h > 0.0) {
        let first = plate.deflect_varying(rigidity, &filled(&vec![0.0; load.len()]), &load);
        plate.deflect_varying(rigidity, &filled(&first), &load)
    } else {
        plate.deflect(rigidity, open, &load)
    };
    let stress = plate.surface_stress(&w, elastic, true);
    if std::env::var_os("EXOGEN_DEBUG").is_some() {
        let most = load.iter().fold(0.0f64, |m, &q| m.max(q));
        let (lo, hi) = stress.tensor.iter().fold((0.0f32, 0.0f32), |(lo, hi), t| {
            let trace = 0.5 * (t[0] + t[1] + t[2]);
            (lo.min(trace), hi.max(trace))
        });
        let (wlo, whi) = w
            .iter()
            .fold((0.0f64, 0.0f64), |(a, b), &v| (a.min(v), b.max(v)));
        eprintln!(
            "load: lid {:.0} km, flexural parameter {:.0} km, plate face {}, load up to {:.1} MPa, mean stress {:.0} to {:.0} MPa, deflection {:.0} to {:.0} m",
            elastic / 1e3,
            flexural / 1e3,
            plate.grid.face_size(),
            most / 1e6,
            lo / 1e6,
            hi / 1e6,
            wlo,
            whi
        );
    }
    Some(Laid {
        stress,
        grid: target,
        deflection: w.iter().map(|&v| v as f32).collect(),
    })
}

/// The graben the loads open (F158): where the floods' loads stretch the lid's top past
/// what breaks its faulted layer, outward beyond the margins of the filled basins, normal
/// faults run round them as the arcuate graben about the lunar maria do (Rimae Hippalus
/// about Mare Humorum; Solomon and Head 1980). They spend the extension the loads add past
/// the yield, the field's tension over the lid's stiffness beyond what the field holds
/// without them, on graben within the faulted layer, as the lid's own faults are (F147),
/// from the faults' generator (F111, F112): each starts where the loads stretch the lid
/// past the yield, runs across the most tensile stress as the field leads it while the
/// loads still stretch the lid there, and forms as the stretch at its start grew through
/// the epochs; those formed before the surface began to keep a record are gone with it.
pub fn graben(
    history: &FieldHistory,
    grid: CubeSphere,
    seed: u64,
    age: f64,
    record_start: f64,
) -> Vec<crate::stages::plume::Trough> {
    use crate::stages::plume::{FaultLayer, FaultRules, populate, trace_in};
    if !history.loaded() {
        return Vec::new();
    }
    let body_radius = history.body_radius;
    let present = history.field(age);
    let yield_tension = present.yield_tension;
    let stiffness = ROCK_YOUNG_MODULUS / (1.0 - ROCK_POISSON);
    // What the loads add to the stretch past the yield, on the coarse grid the loop reads
    // the field on: at the present, and at each epoch that laid a load.
    let coarse = CubeSphere::new(COMPRESSION_FACE);
    let added_at =
        |field: &crate::stages::stress::StressField,
         bare: &crate::stages::stress::StressField,
         p: [f64; 3],
         t: f64| { (field.excess(p, t) - bare.excess(p, t)).max(0.0) };
    let bare = history.field_unloaded(age);
    let added: Vec<f64> = (0..coarse.cell_count())
        .into_par_iter()
        .map(|i| added_at(&present, &bare, coarse.center(coarse.cell(i)), age))
        .collect();
    let stretched: Vec<usize> = (0..added.len()).filter(|&i| added[i] > 0.0).collect();
    if stretched.is_empty() {
        return Vec::new();
    }
    let loads: Vec<f64> = history
        .epochs
        .iter()
        .filter(|e| e.load.is_some())
        .map(|e| e.time)
        .collect();
    // Each stretched cell's added stretch at each loading epoch.
    let mut grown: Vec<Vec<f64>> = vec![Vec::new(); added.len()];
    for &t in &loads {
        let field = history.field(t);
        let bare = history.field_unloaded(t);
        let at: Vec<f64> = stretched
            .par_iter()
            .map(|&i| added_at(&field, &bare, coarse.center(coarse.cell(i)), t))
            .collect();
        for (&i, v) in stretched.iter().zip(at) {
            grown[i].push(v);
        }
    }
    let area = |i: usize| coarse.cell_area(coarse.cell(i)) * body_radius * body_radius;
    let mut cumulative = Vec::with_capacity(added.len() + 1);
    cumulative.push(0.0);
    let mut sum = 0.0;
    for (i, &a) in added.iter().enumerate() {
        sum += a * area(i);
        cumulative.push(sum);
    }
    let budget = sum / stiffness;
    let rng = crate::rng::Rng::new(seed ^ 0x10AD_6A2B);
    let n = added.len();
    let start_of = |draw: u64| -> [f64; 3] {
        let target = rng.uniform(draw, 20) * sum;
        let i = cumulative
            .partition_point(|&c| c <= target)
            .saturating_sub(1)
            .min(n - 1);
        let p = coarse.center(coarse.cell(i));
        let (east, north) = local_basis(p);
        let (u, v) = (rng.uniform(draw, 22) - 0.5, rng.uniform(draw, 23) - 0.5);
        let step = coarse.delta();
        crate::geom::normalize([
            p[0] + step * (u * east[0] + v * north[0]),
            p[1] + step * (u * east[1] + v * north[1]),
            p[2] + step * (u * east[2] + v * north[2]),
        ])
    };
    let cell_of = |p: [f64; 3]| coarse.index(coarse.locate(p).0);
    let rules = FaultRules {
        l_min: RIFT_MIN_LENGTH,
        l_max: 0.5 * std::f64::consts::PI * body_radius,
        depth_per_length: RIFT_DEPTH_PER_LENGTH,
        // Graben within the faulted layer, as the lid's own faults: the loads' stretch is
        // small, and none of its faults would rift (a stated choice).
        depth_cap: FAULT_LAYER_RANGE.0,
        compliance: 1.0,
        layer: FaultLayer::Surface,
        cell: grid.delta() * body_radius,
        body_radius,
    };
    let unstretched = |q: [f64; 3]| added[cell_of(q)] <= 0.0;
    let fresh = |draw: u64, length: f64, spacing: f64, width: f64| {
        trace_in(
            &|q| present.at(q, age),
            start_of(draw),
            length,
            spacing,
            body_radius,
            &rng,
            draw,
            width,
            None,
            yield_tension,
            Some(&unstretched),
        )
    };
    // A graben forms as the stretch the loads added at its start grew: when it reached a
    // uniform share of what it came to, within that epoch.
    let timing = |u: f64, _spent: f64, draw: u64| -> f64 {
        let i = cell_of(start_of(draw));
        let target = u * added[i];
        let mut before = 0.0;
        for (k, &t) in loads.iter().enumerate() {
            let reached = grown[i].get(k).copied().unwrap_or(0.0);
            if reached >= target && reached > 0.0 {
                let previous = if k > 0 { grown[i][k - 1] } else { 0.0 };
                let share = ((target - previous) / (reached - previous).max(1e-9)).clamp(0.0, 1.0);
                return before + (t - before) * share;
            }
            before = t;
        }
        age
    };
    let window = |time: f64, _carried: bool| (time, time);
    let mut faults = populate(
        &rules,
        budget,
        &rng,
        (0x10AD_0000_0000, 60_000_000_000),
        &fresh,
        &timing,
        &window,
        &|_| 0.0,
        None,
    );
    if std::env::var_os("EXOGEN_DEBUG").is_some() {
        let spent: f64 = faults
            .iter()
            .map(|f| f.length * 2.0 * f.depth / GRABEN_DIP_DEG.to_radians().tan())
            .sum();
        eprintln!(
            "load graben: budget {:.3e} m2 over {} stretched cells, {} drawn, heave x length about {:.3e}",
            budget,
            stretched.len(),
            faults.len(),
            spent
        );
    }
    faults.retain(|f| f.start_gyr >= record_start);
    faults
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::flexure::{Plate, deflect};
    use crate::geom::angle_between;

    /// The plate's Laplacian of a spherical harmonic is −l (l + 1) / R² times it, on every
    /// cell, across the cube's edges and corners too.
    #[test]
    fn the_plate_s_laplacian_finds_a_harmonic_s_degree() {
        let r = 1.0e6;
        let plate = SpherePlate::new(48, r);
        let n = plate.grid.cell_count();
        // Y_2^0 ∝ 3 z² − 1, and a tilted Y_2^1 ∝ x z.
        for (shape, degree) in [(0, 2.0), (1, 2.0)] {
            let f: Vec<f64> = (0..n)
                .map(|i| {
                    let p = plate.grid.center(plate.grid.cell(i));
                    if shape == 0 {
                        3.0 * p[2] * p[2] - 1.0
                    } else {
                        p[0] * p[2]
                    }
                })
                .collect();
            let peak = f.iter().fold(0.0f64, |m, v| m.max(v.abs()));
            let mut out = vec![0.0; n];
            plate.laplace(&f, &mut out);
            let expected: f64 = -degree * (degree + 1.0) / (r * r);
            let worst = (0..n)
                .map(|i| (out[i] - expected * f[i]).abs())
                .fold(0.0f64, f64::max);
            assert!(worst < 0.02 * expected.abs() * peak, "{shape}: {worst}");
        }
    }

    /// A disc of lava centred on a cube's corner, where three faces meet, bends the plate
    /// alike every way round it: its top squeezed within the margin and stretched beyond it
    /// by as much across the faces' edges as along their middles.
    #[test]
    fn a_load_on_a_cube_corner_bends_the_plate_alike_every_way() {
        let (r, g, te) = (1.737e6, 1.62, 40e3);
        let flexural = flexural_parameter(te, g);
        let plate = SpherePlate::new(SpherePlate::face_for(flexural, r), r);
        let n = plate.grid.cell_count();
        let s3 = 1.0 / 3f64.sqrt();
        let centre = [s3, s3, s3];
        let (radius, q) = (300e3, 3000.0 * g * 2000.0);
        let load: Vec<f64> = (0..n)
            .map(|i| {
                let p = plate.grid.center(plate.grid.cell(i));
                if angle_between(centre, p) * r < radius {
                    q
                } else {
                    0.0
                }
            })
            .collect();
        let rigidity = ROCK_YOUNG_MODULUS * te.powi(3) / (12.0 * (1.0 - ROCK_POISSON.powi(2)));
        let w = plate.deflect(rigidity, DENSITY_ROCK * g, &load);
        let stress = plate.surface_stress(&w, te, false);
        let (east, north) = local_basis(centre);
        // The squeeze within and the stretch beyond, out from the centre, round it in
        // thirty-six directions, each along a great circle from the centre.
        let mut squeezes = Vec::new();
        let mut stretches = Vec::new();
        for k in 0..36 {
            let phi = k as f64 * std::f64::consts::TAU / 36.0;
            let dir = [
                east[0] * phi.cos() + north[0] * phi.sin(),
                east[1] * phi.cos() + north[1] * phi.sin(),
                east[2] * phi.cos() + north[2] * phi.sin(),
            ];
            let point = |d: f64| {
                let a = d / r;
                [
                    centre[0] * a.cos() + dir[0] * a.sin(),
                    centre[1] * a.cos() + dir[1] * a.sin(),
                    centre[2] * a.cos() + dir[2] * a.sin(),
                ]
            };
            let read = |p: [f64; 3]| {
                // Out from the centre at p, and round.
                let along = dot(centre, p);
                let out = crate::geom::normalize([
                    p[0] * along - centre[0],
                    p[1] * along - centre[1],
                    p[2] * along - centre[2],
                ]);
                let round = crate::geom::normalize(crate::geom::cross(p, out));
                stress.at(p, out, round)
            };
            squeezes.push(read(point(radius - 0.7 * flexural)).0);
            stretches.push(read(point(radius + 0.7 * flexural)).0);
        }
        let spread = |v: &[f64]| {
            let mean = v.iter().sum::<f64>() / v.len() as f64;
            let worst = v.iter().map(|x| (x - mean).abs()).fold(0.0, f64::max);
            (mean, worst / mean.abs())
        };
        let (squeeze, s_spread) = spread(&squeezes);
        let (stretch, t_spread) = spread(&stretches);
        assert!(
            squeeze < -10e6 && s_spread < 0.35,
            "{squeeze} {s_spread} {squeezes:?}"
        );
        assert!(
            stretch > 10e6 && t_spread < 0.35,
            "{stretch} {t_spread} {stretches:?}"
        );
    }

    /// A disc of lava on the Moon's shell bends it as the shallow shell does (the Hankel
    /// transform's solution of D ∇⁴w + E h w / R² + ρ_m g w = q, ∇²Φ = −E h w / R): its
    /// membrane stiffens it, the squeeze within is greatest outward, so that its thrusts run
    /// round it, and beyond its margin the lid is pulled both ways, outward most.
    #[test]
    fn a_disc_load_bends_the_shell_as_the_shallow_shell_does() {
        let (r, g, te) = (1.737e6, 1.62, 60e3);
        let flexural = flexural_parameter(te, g);
        let plate = SpherePlate::new(SpherePlate::face_for(flexural, r), r);
        let n = plate.grid.cell_count();
        let centre = [1.0, 0.0, 0.0];
        let (radius, q) = (434e3, 3000.0 * g * 2000.0);
        let load: Vec<f64> = (0..n)
            .map(|i| {
                let p = plate.grid.center(plate.grid.cell(i));
                if angle_between(centre, p) * r < radius {
                    q
                } else {
                    0.0
                }
            })
            .collect();
        let rigidity = ROCK_YOUNG_MODULUS * te.powi(3) / (12.0 * (1.0 - ROCK_POISSON.powi(2)));
        let restoring = DENSITY_ROCK * g + plate.membrane_stiffness(te);
        let w = plate.deflect(rigidity, restoring, &load);
        let stress = plate.surface_stress(&w, te, true);
        // Along the equator the direction out from the centre is east.
        let read = |d: f64| {
            let p = [(d / r).cos(), (d / r).sin(), 0.0];
            let (u, v) = local_basis(p);
            (
                w[plate.grid.index(plate.grid.locate(p).0)],
                stress.at(p, u, v),
            )
        };
        // The Hankel solution: 1,471 m at the middle, squeezed by 74 MPa; at 300 km 78 MPa
        // outward and 64 round; at 600 km pulled 16 MPa outward and 12 round.
        let (w0, (o0, r0, _)) = read(0.0);
        assert!((w0 / 1471.0 - 1.0).abs() < 0.1, "{w0}");
        assert!(
            (o0 / -74.3e6 - 1.0).abs() < 0.2 && (r0 / -74.3e6 - 1.0).abs() < 0.2,
            "{o0} {r0}"
        );
        let (_, (o3, r3, _)) = read(300e3);
        assert!(o3 < r3 - 5e6, "{o3} {r3}");
        let (_, (o6, r6, _)) = read(600e3);
        assert!(o6 > 5e6 && r6 > 0.0, "{o6} {r6}");
    }

    /// The ground sinks under its load (F162), but an edifice fills the sag beneath it as
    /// far as it is thick: a shield on a plate stands where it built itself, the ground
    /// about it sunk into a moat, deepest at its foot, and raised in a bulge beyond; a flood
    /// with no edifice sinks with the ground.
    #[test]
    fn a_shield_fills_the_sag_beneath_it_and_stands_in_a_moat() {
        let (r, g, te) = (3.39e6, 3.71, 70e3);
        let grid = CubeSphere::new(128);
        let n = grid.cell_count();
        let centre = [1.0, 0.0, 0.0];
        let (radius, height) = (300e3, 20e3);
        let shield: Vec<f32> = (0..n)
            .map(|k| {
                let x = angle_between(centre, grid.center(grid.cell(k))) * r / radius;
                if x < 1.0 {
                    (height * (1.0 - x * x)) as f32
                } else {
                    0.0
                }
            })
            .collect();
        let laid = lay(grid, &shield, r, g, te, 2_900.0, Some(&shield)).expect("a load");
        let w = laid.deflection_on(grid);
        let sunk = sink(&w, Some(&shield));
        let at = |d: f64| grid.index(grid.locate([(d / r).cos(), (d / r).sin(), 0.0]).0);
        // The summit: the lid sags kilometres beneath it, and the edifice fills it all.
        let summit = at(0.0);
        assert!(w[summit] > 2_000.0, "{}", w[summit]);
        assert!(sunk[summit].abs() < 1.0, "{}", sunk[summit]);
        // Just past the foot the ground sinks into a moat, and farther off rises a little.
        let moat = sunk[at(1.1 * radius)];
        assert!(moat > 300.0, "{moat}");
        let far: f32 = (15..40)
            .map(|k| sunk[at(radius + k as f64 * 25e3)])
            .fold(f32::MAX, f32::min);
        assert!(far < 0.0, "{far}");
        // Without an edifice the ground sinks with the lid.
        let bare = sink(&w, None);
        assert_eq!(bare[summit], w[summit]);
    }

    /// An edifice's root weighs on the lid (F166): its lava filling the sag beneath it, the
    /// mantle pushes back there by the difference between the rock's density and the lava's,
    /// so that the lid sinks further beneath it and its moat deepens with it, by some half
    /// again on a thick lid, while the lid's strength keeps it within a few times the sag of
    /// the edifice alone; where no edifice stands the lid bends as before.
    #[test]
    fn an_edifice_s_root_deepens_its_sag_and_its_moat() {
        let (r, g, te) = (3.39e6, 3.71, 70e3);
        let grid = CubeSphere::new(128);
        let n = grid.cell_count();
        let centre = [1.0, 0.0, 0.0];
        let (radius, height) = (300e3, 20e3);
        let shield: Vec<f32> = (0..n)
            .map(|k| {
                let x = angle_between(centre, grid.center(grid.cell(k))) * r / radius;
                if x < 1.0 {
                    (height * (1.0 - x * x)) as f32
                } else {
                    0.0
                }
            })
            .collect();
        let at = |d: f64| grid.index(grid.locate([(d / r).cos(), (d / r).sin(), 0.0]).0);
        let bare = lay(grid, &shield, r, g, te, 2_900.0, None)
            .unwrap()
            .deflection_on(grid);
        let root = lay(grid, &shield, r, g, te, 2_900.0, Some(&shield))
            .unwrap()
            .deflection_on(grid);
        let (summit, moat) = (at(0.0), at(1.1 * radius));
        assert!(
            root[summit] > 1.3 * bare[summit],
            "{} {}",
            root[summit],
            bare[summit]
        );
        assert!(
            root[summit] < 3.0 * bare[summit],
            "{} {}",
            root[summit],
            bare[summit]
        );
        assert!(
            root[moat] > 1.3 * bare[moat],
            "{} {}",
            root[moat],
            bare[moat]
        );
        // The edifice still fills its own sag, standing where it built itself.
        assert!(sink(&root, Some(&shield))[summit].abs() < 1.0);
        // A flood with no edifice bends the lid as an empty sag does.
        let none = vec![0.0f32; n];
        let open = lay(grid, &shield, r, g, te, 2_900.0, Some(&none))
            .unwrap()
            .deflection_on(grid);
        assert_eq!(open[summit], bare[summit]);
    }

    /// A basin filled with lava opens graben beyond its margin (F158), where its load
    /// stretches the lid outward, running round it as the arcuate graben about the maria
    /// do, at the time it was filled, and none before the record began.
    #[test]
    fn a_filled_basin_opens_graben_round_its_margin() {
        use crate::stages::stress::{FieldEpoch, Relief};
        let grid = CubeSphere::new(128);
        let (r, g, te) = (1_737e3, 1.62, 60e3);
        let n = grid.cell_count();
        let centre = crate::geom::normalize([1.0, 0.3, 0.2]);
        let radius = 0.25;
        let thickness: Vec<f32> = (0..n)
            .map(|k| {
                if angle_between(centre, grid.center(grid.cell(k))) < radius {
                    2_000.0
                } else {
                    0.0
                }
            })
            .collect();
        let flat = Relief {
            grid,
            height: vec![0.0; n],
        };
        let history = FieldHistory {
            rises: Vec::new(),
            limbs: Default::default(),
            body_radius: r,
            gravity: g,
            epochs: vec![
                FieldEpoch {
                    time: 0.0,
                    elastic: te,
                    relief: flat.clone(),
                    load: None,
                },
                FieldEpoch {
                    time: 3.0,
                    elastic: te,
                    relief: flat,
                    load: lay(grid, &thickness, r, g, te, 3_000.0, None).map(|l| l.stress),
                },
            ],
        };
        let faults = graben(&history, grid, 7, 4.5, 0.0);
        assert!(faults.len() >= 5, "{}", faults.len());
        let flexural = flexural_parameter(te, g) / r;
        // Near the margin, within a flexural parameter of it, and beyond.
        let (mut steps, mut round) = ([0usize; 2], [0usize; 2]);
        for f in &faults {
            assert!(
                f.start_gyr > 0.0 && f.start_gyr <= 3.0 + 1e-9,
                "{}",
                f.start_gyr
            );
            let m = f.points[f.points.len() / 2];
            let d = angle_between(centre, m);
            assert!(d > radius - flexural && d < radius + 5.0 * flexural, "{d}");
            for w in f.points.windows(2) {
                let band = usize::from(angle_between(centre, w[0]) > radius + flexural);
                let along = dot(centre, w[0]);
                let out = crate::geom::normalize([
                    w[0][0] * along - centre[0],
                    w[0][1] * along - centre[1],
                    w[0][2] * along - centre[2],
                ]);
                let step = crate::geom::normalize([
                    w[1][0] - w[0][0],
                    w[1][1] - w[0][1],
                    w[1][2] - w[0][2],
                ]);
                steps[band] += 1;
                if dot(out, step).abs() < 0.5 {
                    round[band] += 1;
                }
            }
        }
        eprintln!(
            "near: {} of {} round; beyond: {} of {} round",
            round[0], steps[0], round[1], steps[1]
        );
        assert!(round[0] * 10 > steps[0] * 7, "{round:?} of {steps:?}");
        // Formed before the record began, they are gone with it.
        assert!(graben(&history, grid, 7, 4.5, 3.1).is_empty());
    }

    /// A disc of lava on the sphere bends the plate as the axisymmetric plate does: the
    /// same deflection at its centre and a bulge beyond its margin; its top squeezed within
    /// and stretched round its margin.
    #[test]
    fn a_disc_load_bends_the_plate_as_the_axisymmetric_solution() {
        let (r, g, te) = (1.737e6, 1.62, 40e3);
        let flexural = flexural_parameter(te, g);
        let plate = SpherePlate::new(SpherePlate::face_for(flexural, r), r);
        let n = plate.grid.cell_count();
        let centre = [1.0, 0.0, 0.0];
        let (radius, q) = (300e3, 3000.0 * g * 2000.0);
        let load: Vec<f64> = (0..n)
            .map(|i| {
                let p = plate.grid.center(plate.grid.cell(i));
                if angle_between(centre, p) * r < radius {
                    q
                } else {
                    0.0
                }
            })
            .collect();
        let rigidity = ROCK_YOUNG_MODULUS * te.powi(3) / (12.0 * (1.0 - ROCK_POISSON.powi(2)));
        let w = plate.deflect(rigidity, DENSITY_ROCK * g, &load);
        // Against the axisymmetric plate (flexure.rs), up positive there.
        let axis = Plate {
            rigidity,
            elastic: te,
            restoring: DENSITY_ROCK * g,
        };
        let step = 5e3;
        let profile: Vec<f64> = (0..400)
            .map(|k| if (k as f64) * step < radius { -q } else { 0.0 })
            .collect();
        let up = deflect(&axis, step, &profile);
        let at = |d: f64| {
            let p = [(d / r).cos(), (d / r).sin(), 0.0];
            w[plate.grid.index(plate.grid.locate(p).0)]
        };
        let centre_w = at(0.0);
        assert!(
            (centre_w / -up[0] - 1.0).abs() < 0.1,
            "{centre_w} vs {}",
            -up[0]
        );
        // The bulge: the axisymmetric plate's highest point beyond the margin.
        let (k_bulge, bulge) = up
            .iter()
            .enumerate()
            .skip((radius / step) as usize)
            .max_by(|a, b| a.1.total_cmp(b.1))
            .unwrap();
        assert!(*bulge > 0.0);
        assert!(
            at(k_bulge as f64 * step) < 0.0,
            "{}",
            at(k_bulge as f64 * step)
        );
        // The stress: squeezed just within the margin, stretched round it just without.
        let stress = plate.surface_stress(&w, te, false);
        let read = |d: f64| {
            let p = [(d / r).cos(), (d / r).sin(), 0.0];
            let (u, v) = local_basis(p);
            stress.at(p, u, v)
        };
        // Along the equator the direction out from the centre is east: squeezed out and
        // round within the margin, stretched outward beyond it.
        let inside = read(radius - 0.7 * flexural);
        let outside = read(radius + 0.7 * flexural);
        assert!(inside.0 < -10e6 && inside.1 < 0.0, "{inside:?}");
        assert!(outside.0 > 10e6, "{outside:?}");
    }
}
