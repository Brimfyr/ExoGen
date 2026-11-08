// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! Small spherical geometry helpers shared by the stages.

/// Rotate a vector about a unit axis by an angle in radians (Rodrigues).
pub fn rotate(d: [f64; 3], axis: [f64; 3], angle: f64) -> [f64; 3] {
    let (s, c) = angle.sin_cos();
    let dot = axis[0] * d[0] + axis[1] * d[1] + axis[2] * d[2];
    let cross = cross(axis, d);
    [
        d[0] * c + cross[0] * s + axis[0] * dot * (1.0 - c),
        d[1] * c + cross[1] * s + axis[1] * dot * (1.0 - c),
        d[2] * c + cross[2] * s + axis[2] * dot * (1.0 - c),
    ]
}

/// A turn of the body under its spin axis (F104): the rotation carrying the body-fixed
/// directions of before into the frame of the present spin axis; none by default.
#[derive(Clone, Copy, Debug, Default)]
pub struct Turn {
    pub axis: [f64; 3],
    pub angle: f64,
}

impl Turn {
    /// Where a direction of before stands now.
    pub fn apply(&self, d: [f64; 3]) -> [f64; 3] {
        if self.angle == 0.0 {
            d
        } else {
            normalize(rotate(d, self.axis, self.angle))
        }
    }

    /// Where a direction of now stood before.
    pub fn undo(&self, d: [f64; 3]) -> [f64; 3] {
        if self.angle == 0.0 {
            d
        } else {
            normalize(rotate(d, self.axis, -self.angle))
        }
    }
}

#[inline]
pub fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
pub fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
pub fn normalize(v: [f64; 3]) -> [f64; 3] {
    let n = dot(v, v).sqrt().max(1e-300);
    [v[0] / n, v[1] / n, v[2] / n]
}

/// Angle in radians between two unit vectors, accurate for small angles.
#[inline]
pub fn angle_between(a: [f64; 3], b: [f64; 3]) -> f64 {
    let c = cross(a, b);
    dot(c, c).sqrt().atan2(dot(a, b))
}

/// Local east and north unit vectors at a direction on the sphere (+Z is the pole). At the
/// poles an arbitrary consistent frame is returned.
pub fn local_basis(center: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    let up = [0.0, 0.0, 1.0];
    let mut east = cross(up, center);
    if dot(east, east) < 1e-12 {
        east = [1.0, 0.0, 0.0];
    }
    let east = normalize(east);
    let north = normalize(cross(center, east));
    (east, north)
}

/// Azimuth of direction `d` as seen from `center`, radians from east toward north.
pub fn azimuth_from(center: [f64; 3], east: [f64; 3], north: [f64; 3], d: [f64; 3]) -> f64 {
    let rel = [d[0] - center[0], d[1] - center[1], d[2] - center[2]];
    dot(rel, north).atan2(dot(rel, east))
}

/// Point on the unit sphere at an angular distance `theta` from `center` along azimuth `phi`.
pub fn offset(center: [f64; 3], east: [f64; 3], north: [f64; 3], theta: f64, phi: f64) -> [f64; 3] {
    let (st, ct) = theta.sin_cos();
    let (sp, cp) = phi.sin_cos();
    normalize([
        center[0] * ct + (east[0] * cp + north[0] * sp) * st,
        center[1] * ct + (east[1] * cp + north[1] * sp) * st,
        center[2] * ct + (east[2] * cp + north[2] * sp) * st,
    ])
}

/// The eigenvalues of a symmetric 3x3 matrix, largest first, and their unit eigenvectors
/// in the same order (cyclic Jacobi rotations).
pub fn symmetric_eigen(m: [[f64; 3]; 3]) -> ([f64; 3], [[f64; 3]; 3]) {
    let mut a = m;
    // Column k of v is the k-th eigenvector.
    let mut v = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    for _ in 0..64 {
        let scale = a[0][0].abs() + a[1][1].abs() + a[2][2].abs();
        let off = a[0][1].abs() + a[0][2].abs() + a[1][2].abs();
        if off <= 1e-15 * scale.max(f64::MIN_POSITIVE) {
            break;
        }
        for (p, q) in [(0, 1), (0, 2), (1, 2)] {
            if a[p][q] == 0.0 {
                continue;
            }
            let theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
            let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
            let c = 1.0 / (t * t + 1.0).sqrt();
            let s = t * c;
            for row in a.iter_mut() {
                let (x, y) = (row[p], row[q]);
                row[p] = c * x - s * y;
                row[q] = s * x + c * y;
            }
            let (mut row_p, mut row_q) = (a[p], a[q]);
            for (x, y) in row_p.iter_mut().zip(row_q.iter_mut()) {
                let (u, w) = (*x, *y);
                *x = c * u - s * w;
                *y = s * u + c * w;
            }
            a[p] = row_p;
            a[q] = row_q;
            for row in v.iter_mut() {
                let (x, y) = (row[p], row[q]);
                row[p] = c * x - s * y;
                row[q] = s * x + c * y;
            }
        }
    }
    let mut order = [0usize, 1, 2];
    order.sort_by(|&i, &j| a[j][j].total_cmp(&a[i][i]));
    let values = order.map(|i| a[i][i]);
    let vectors = order.map(|i| normalize([v[0][i], v[1][i], v[2][i]]));
    (values, vectors)
}

/// Steffen's (1990) cubic at `t`, 0 to 1 between the second and third of four values a unit
/// apart, `v`: Hermite's between them, each one's slope the gentlest of twice the steps either
/// side and their mean, none where the values turn, so that it never passes the two values it
/// lies between and makes no new extreme.
pub(crate) fn monotone(v: [f64; 4], t: f64) -> f64 {
    let slope = |a: f64, b: f64, c: f64| {
        let (s0, s1) = (b - a, c - b);
        if s0 * s1 <= 0.0 {
            0.0
        } else {
            s0.signum()
                * (2.0 * s0.abs())
                    .min(2.0 * s1.abs())
                    .min(0.5 * (s0 + s1).abs())
        }
    };
    let (m1, m2) = (slope(v[0], v[1], v[2]), slope(v[1], v[2], v[3]));
    let (t2, t3) = (t * t, t * t * t);
    (2.0 * t3 - 3.0 * t2 + 1.0) * v[1]
        + (t3 - 2.0 * t2 + t) * m1
        + (-2.0 * t3 + 3.0 * t2) * v[2]
        + (t3 - t2) * m2
}

/// The standard normal distribution's quantile at `p`, 0 to 1: Acklam's rational
/// approximation, its relative error under 1.2e-9.
pub fn normal_quantile(p: f64) -> f64 {
    const A: [f64; 6] = [
        -3.969683028665376e1,
        2.209460984245205e2,
        -2.759285104469687e2,
        1.38357751867269e2,
        -3.066479806614716e1,
        2.506628277459239,
    ];
    const B: [f64; 5] = [
        -5.447609879822406e1,
        1.615858368580409e2,
        -1.556989798598866e2,
        6.680131188771972e1,
        -1.328068155288572e1,
    ];
    const C: [f64; 6] = [
        -7.784894002430293e-3,
        -3.223964580411365e-1,
        -2.400758277161838,
        -2.549732539343734,
        4.374664141464968,
        2.938163982698783,
    ];
    const D: [f64; 4] = [
        7.784695709041462e-3,
        3.224671290700398e-1,
        2.445134137142996,
        3.754408661907416,
    ];
    // The tails by their own approximation in sqrt(-2 ln p), the middle by its own in p.
    let tail = |q: f64| {
        (((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5])
            / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1.0)
    };
    let p = p.clamp(1e-300, 1.0 - 1e-16);
    if p < 0.02425 {
        tail((-2.0 * p.ln()).sqrt())
    } else if p > 1.0 - 0.02425 {
        -tail((-2.0 * (1.0 - p).ln()).sqrt())
    } else {
        let q = p - 0.5;
        let r = q * q;
        (((((A[0] * r + A[1]) * r + A[2]) * r + A[3]) * r + A[4]) * r + A[5]) * q
            / (((((B[0] * r + B[1]) * r + B[2]) * r + B[3]) * r + B[4]) * r + 1.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The normal quantile at the measured points of the distribution, symmetric about its
    /// middle and rising through it.
    #[test]
    fn the_normal_quantile_stands_where_the_distribution_does() {
        for (p, z) in [
            (0.5, 0.0),
            (0.841_344_746_068_543, 1.0),
            (0.975, 1.959_963_984_540_054),
            (0.999, 3.090_232_306_167_813),
            (1e-6, -4.753_424_308_822_899),
        ] {
            let q = normal_quantile(p);
            assert!(
                (q - z).abs() < 1e-8 * z.abs().max(1.0),
                "{p}: {q} against {z}"
            );
            assert!((normal_quantile(1.0 - p) + q).abs() < 1e-7, "{p}");
        }
        let mut last = f64::NEG_INFINITY;
        for k in 1..1000 {
            let q = normal_quantile(k as f64 / 1000.0);
            assert!(q > last, "{k}");
            last = q;
        }
    }

    #[test]
    fn rotation_preserves_length_and_rotates_by_the_angle() {
        let d = normalize([1.0, 2.0, 0.5]);
        let axis = normalize([0.3, -0.2, 1.0]);
        let r = rotate(d, axis, 0.7);
        assert!((dot(r, r) - 1.0).abs() < 1e-12);
        let perp_d = normalize([
            d[0] - axis[0] * dot(d, axis),
            d[1] - axis[1] * dot(d, axis),
            d[2] - axis[2] * dot(d, axis),
        ]);
        let perp_r = normalize([
            r[0] - axis[0] * dot(r, axis),
            r[1] - axis[1] * dot(r, axis),
            r[2] - axis[2] * dot(r, axis),
        ]);
        assert!((angle_between(perp_d, perp_r) - 0.7).abs() < 1e-9);
    }

    #[test]
    fn offset_and_azimuth_round_trip() {
        let c = normalize([0.2, 0.5, 0.8]);
        let (e, n) = local_basis(c);
        for &(theta, phi) in &[(0.1, 0.3), (0.5, -2.0), (0.01, 3.0)] {
            let p = offset(c, e, n, theta, phi);
            assert!((angle_between(c, p) - theta).abs() < 1e-9);
            let back = azimuth_from(c, e, n, p);
            let diff = (back - phi + std::f64::consts::PI).rem_euclid(2.0 * std::f64::consts::PI)
                - std::f64::consts::PI;
            assert!(diff.abs() < 1e-6, "{back} vs {phi}");
        }
    }
}
