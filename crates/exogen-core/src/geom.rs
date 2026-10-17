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

#[cfg(test)]
mod tests {
    use super::*;

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
