// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! Equiangular cubesphere.
//!
//! Six faces of `f x f` cells. Cell `(face, i, j)` has centre angles
//! `alpha = -pi/4 + (i + 1/2) delta` and `beta = -pi/4 + (j + 1/2) delta` with
//! `delta = pi / (2 f)`, and the unit direction `normalize(n + tan(alpha) r + tan(beta) t)`
//! for the face basis `(n, r, t)`. The equiangular mapping keeps the cell-area ratio between
//! face centre and corner near 1.3.
//!
//! Neighbours that cross a face edge are resolved geometrically: step one cell in angle
//! beyond the edge and locate the direction. This picks the cell on the adjacent face that
//! contains the stepped direction, which is exact along the edge midline and within one
//! cell toward the corners. Stencil operators that need higher accuracy at edges interpolate
//! from ghost cells; that belongs to the field operators, not the grid.

use std::f64::consts::{FRAC_PI_2, FRAC_PI_4};

use super::{dot, normalize};

/// Face bases: normal, right (alpha axis), up (beta axis). Right-handed, `r x t = n`.
/// Faces 0..4 are +X, -X, +Y, -Y with +Z as up; faces 4 and 5 are +Z and -Z.
const FACE_BASIS: [[[f64; 3]; 3]; 6] = [
    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]],
    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
    [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]],
    [[0.0, 0.0, -1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]],
];

/// A cell address: face index 0..6 and in-face indices `0..f`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Cell {
    pub face: u8,
    pub i: u32,
    pub j: u32,
}

/// An equiangular cubesphere with `f` cells along each face edge.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CubeSphere {
    f: u32,
}

impl CubeSphere {
    pub fn new(face_size: u32) -> Self {
        assert!(face_size >= 2, "face size must be at least 2");
        Self { f: face_size }
    }

    pub fn face_size(&self) -> u32 {
        self.f
    }

    pub fn cell_count(&self) -> usize {
        6 * self.f as usize * self.f as usize
    }

    /// Angular width of a cell, `pi / (2 f)`.
    pub fn delta(&self) -> f64 {
        FRAC_PI_2 / self.f as f64
    }

    /// Linear index of a cell: faces are contiguous, rows (`j`) within a face, `i` fastest.
    pub fn index(&self, c: Cell) -> usize {
        let f = self.f as usize;
        (c.face as usize * f + c.j as usize) * f + c.i as usize
    }

    pub fn cell(&self, index: usize) -> Cell {
        let f = self.f as usize;
        let per_face = f * f;
        let face = index / per_face;
        let r = index % per_face;
        Cell {
            face: face as u8,
            i: (r % f) as u32,
            j: (r / f) as u32,
        }
    }

    /// Cell-centre face angles `(alpha, beta)`.
    pub fn angles(&self, c: Cell) -> (f64, f64) {
        let d = self.delta();
        (
            -FRAC_PI_4 + (c.i as f64 + 0.5) * d,
            -FRAC_PI_4 + (c.j as f64 + 0.5) * d,
        )
    }

    /// Unit direction for face angles. Angles may lie outside `[-pi/4, pi/4]`, in which case
    /// the direction lies on a neighbouring face.
    pub fn direction(face: u8, alpha: f64, beta: f64) -> [f64; 3] {
        let [n, r, t] = FACE_BASIS[face as usize];
        let (u, v) = (alpha.tan(), beta.tan());
        normalize([
            n[0] + u * r[0] + v * t[0],
            n[1] + u * r[1] + v * t[1],
            n[2] + u * r[2] + v * t[2],
        ])
    }

    /// Unit direction of a cell centre.
    pub fn center(&self, c: Cell) -> [f64; 3] {
        let (a, b) = self.angles(c);
        Self::direction(c.face, a, b)
    }

    /// The face containing a direction and the face angles of that direction.
    pub fn face_angles(d: [f64; 3]) -> (u8, f64, f64) {
        let (ax, ay, az) = (d[0].abs(), d[1].abs(), d[2].abs());
        let face = if ax >= ay && ax >= az {
            if d[0] >= 0.0 { 0 } else { 1 }
        } else if ay >= az {
            if d[1] >= 0.0 { 2 } else { 3 }
        } else if d[2] >= 0.0 {
            4
        } else {
            5
        };
        let [n, r, t] = FACE_BASIS[face];
        let dn = dot(d, n);
        (face as u8, (dot(d, r) / dn).atan(), (dot(d, t) / dn).atan())
    }

    /// The cell containing a direction, with continuous in-face coordinates `(x, y)` in
    /// `[0, f)`; the cell centre is at `(i + 1/2, j + 1/2)`.
    pub fn locate(&self, d: [f64; 3]) -> (Cell, f64, f64) {
        let (face, a, b) = Self::face_angles(d);
        let delta = self.delta();
        let x = (a + FRAC_PI_4) / delta;
        let y = (b + FRAC_PI_4) / delta;
        let max = (self.f - 1) as i64;
        let i = (x.floor() as i64).clamp(0, max) as u32;
        let j = (y.floor() as i64).clamp(0, max) as u32;
        (Cell { face, i, j }, x, y)
    }

    /// Neighbour reached by stepping `(di, dj)` cells. Steps that leave the face are
    /// resolved geometrically (see the module documentation).
    pub fn neighbor(&self, c: Cell, di: i32, dj: i32) -> Cell {
        let f = self.f as i64;
        let ni = c.i as i64 + di as i64;
        let nj = c.j as i64 + dj as i64;
        if (0..f).contains(&ni) && (0..f).contains(&nj) {
            return Cell {
                face: c.face,
                i: ni as u32,
                j: nj as u32,
            };
        }
        let (a, b) = self.angles(c);
        let delta = self.delta();
        let d = Self::direction(c.face, a + di as f64 * delta, b + dj as f64 * delta);
        self.locate(d).0
    }

    /// Exact solid angle of a cell in steradians. A cell is a rectangle `[u0, u1] x [v0, v1]`
    /// on the face plane at unit distance, whose solid angle is
    /// `sum over corners of (+/-) atan(u v / sqrt(1 + u^2 + v^2))`.
    pub fn cell_area(&self, c: Cell) -> f64 {
        let delta = self.delta();
        let a0 = -FRAC_PI_4 + c.i as f64 * delta;
        let b0 = -FRAC_PI_4 + c.j as f64 * delta;
        rect_solid_angle(a0.tan(), (a0 + delta).tan(), b0.tan(), (b0 + delta).tan())
    }

    /// Iterate over all cells in index order.
    pub fn cells(&self) -> impl Iterator<Item = Cell> + '_ {
        (0..self.cell_count()).map(move |k| self.cell(k))
    }
}

/// Solid angle subtended at the origin by the rectangle `[u0, u1] x [v0, v1]` on the plane
/// at unit distance along the face normal.
fn rect_solid_angle(u0: f64, u1: f64, v0: f64, v1: f64) -> f64 {
    let term = |u: f64, v: f64| (u * v / (1.0 + u * u + v * v).sqrt()).atan();
    term(u1, v1) - term(u0, v1) - term(u1, v0) + term(u0, v0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn index_round_trip() {
        let g = CubeSphere::new(7);
        for k in 0..g.cell_count() {
            assert_eq!(g.index(g.cell(k)), k);
        }
    }

    #[test]
    fn face_bases_are_right_handed_unit_frames() {
        for basis in FACE_BASIS {
            let [n, r, t] = basis;
            let cross = [
                r[1] * t[2] - r[2] * t[1],
                r[2] * t[0] - r[0] * t[2],
                r[0] * t[1] - r[1] * t[0],
            ];
            for k in 0..3 {
                assert!((cross[k] - n[k]).abs() < 1e-12);
            }
            assert!(dot(n, r).abs() < 1e-12 && dot(n, t).abs() < 1e-12 && dot(r, t).abs() < 1e-12);
        }
    }

    #[test]
    fn locate_recovers_cell_centres() {
        let g = CubeSphere::new(16);
        for c in g.cells() {
            let (found, x, y) = g.locate(g.center(c));
            assert_eq!(found, c);
            assert!((x - (c.i as f64 + 0.5)).abs() < 1e-9);
            assert!((y - (c.j as f64 + 0.5)).abs() < 1e-9);
        }
    }

    #[test]
    fn areas_sum_to_the_sphere_exactly() {
        for f in [2u32, 3, 16, 64] {
            let g = CubeSphere::new(f);
            let total: f64 = g.cells().map(|c| g.cell_area(c)).sum();
            let rel = (total - 4.0 * PI).abs() / (4.0 * PI);
            assert!(rel < 1e-12, "f={f} rel error {rel}");
        }
        // A whole face is one sixth of the sphere.
        assert!((rect_solid_angle(-1.0, 1.0, -1.0, 1.0) - 4.0 * PI / 6.0).abs() < 1e-12);
    }

    #[test]
    fn area_ratio_is_mild() {
        let g = CubeSphere::new(64);
        let centre = g.cell_area(Cell {
            face: 0,
            i: 32,
            j: 32,
        });
        let corner = g.cell_area(Cell {
            face: 0,
            i: 0,
            j: 0,
        });
        let ratio = centre / corner;
        assert!(ratio > 1.2 && ratio < 1.4, "ratio {ratio}");
    }

    #[test]
    fn interior_neighbours_are_index_steps() {
        let g = CubeSphere::new(8);
        let c = Cell {
            face: 2,
            i: 3,
            j: 4,
        };
        assert_eq!(
            g.neighbor(c, 1, 0),
            Cell {
                face: 2,
                i: 4,
                j: 4
            }
        );
        assert_eq!(
            g.neighbor(c, 0, -1),
            Cell {
                face: 2,
                i: 3,
                j: 3
            }
        );
    }

    #[test]
    fn edge_neighbours_are_adjacent_on_the_sphere() {
        // Every edge-crossing neighbour must lie on a different face and be about one cell
        // away in angle, and stepping back must return to within one cell of the start.
        let g = CubeSphere::new(12);
        let delta = g.delta();
        let mut crossings = 0;
        for c in g.cells() {
            for (di, dj) in [(1, 0), (-1, 0), (0, 1), (0, -1)] {
                let n = g.neighbor(c, di, dj);
                let same_face_step = (c.i as i64 + di as i64) >= 0
                    && (c.i as i64 + di as i64) < g.face_size() as i64
                    && (c.j as i64 + dj as i64) >= 0
                    && (c.j as i64 + dj as i64) < g.face_size() as i64;
                if same_face_step {
                    continue;
                }
                crossings += 1;
                assert_ne!(
                    n.face, c.face,
                    "edge step stayed on the face: {c:?} -> {n:?}"
                );
                let ang = dot(g.center(c), g.center(n)).clamp(-1.0, 1.0).acos();
                assert!(
                    ang < 1.6 * delta && ang > 0.4 * delta,
                    "{c:?} -> {n:?} angle {ang} delta {delta}"
                );
                // Step back: the reverse direction on the other face is not simply (-di, -dj),
                // so test by geometry: some 4-neighbour of n is within one cell of c.
                let back = [(1, 0), (-1, 0), (0, 1), (0, -1)]
                    .iter()
                    .map(|&(a, b)| g.neighbor(n, a, b))
                    .map(|m| dot(g.center(m), g.center(c)).clamp(-1.0, 1.0).acos())
                    .fold(f64::INFINITY, f64::min);
                assert!(
                    back < 1.1 * delta,
                    "no way back from {n:?} to {c:?}: {back}"
                );
            }
        }
        assert_eq!(crossings, 6 * 4 * 12);
    }
}
