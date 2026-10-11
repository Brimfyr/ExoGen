// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! Gradient noise on the sphere. Evaluated at a 3D point (a unit direction scaled by a
//! frequency), so there are no seams and no poles. Seeded from the counter-based generator
//! (ADR 0004); the same seed gives the same field at any resolution.

use crate::rng::Rng;

/// Classic 3D gradient noise with a seeded permutation.
#[derive(Clone, Debug)]
pub struct GradientNoise {
    perm: [u16; 512],
}

const GRADIENTS: [[f64; 3]; 12] = [
    [1.0, 1.0, 0.0],
    [-1.0, 1.0, 0.0],
    [1.0, -1.0, 0.0],
    [-1.0, -1.0, 0.0],
    [1.0, 0.0, 1.0],
    [-1.0, 0.0, 1.0],
    [1.0, 0.0, -1.0],
    [-1.0, 0.0, -1.0],
    [0.0, 1.0, 1.0],
    [0.0, -1.0, 1.0],
    [0.0, 1.0, -1.0],
    [0.0, -1.0, -1.0],
];

impl GradientNoise {
    pub fn new(seed: u64) -> Self {
        let rng = Rng::new(seed);
        let mut p: Vec<u16> = (0..256).collect();
        // Fisher–Yates with keyed draws.
        for i in (1..256usize).rev() {
            let j = (rng.uniform(i as u64, 0) * (i + 1) as f64) as usize;
            p.swap(i, j.min(i));
        }
        let mut perm = [0u16; 512];
        for i in 0..512 {
            perm[i] = p[i & 255];
        }
        Self { perm }
    }

    #[inline]
    fn fade(t: f64) -> f64 {
        t * t * t * (t * (t * 6.0 - 15.0) + 10.0)
    }

    #[inline]
    fn grad(&self, hash: u16, x: f64, y: f64, z: f64) -> f64 {
        let g = GRADIENTS[(hash % 12) as usize];
        g[0] * x + g[1] * y + g[2] * z
    }

    /// Noise in about [-1, 1] at a 3D point.
    pub fn sample(&self, p: [f64; 3]) -> f64 {
        let (xf, yf, zf) = (p[0].floor(), p[1].floor(), p[2].floor());
        let (xi, yi, zi) = (
            (xf as i64 & 255) as usize,
            (yf as i64 & 255) as usize,
            (zf as i64 & 255) as usize,
        );
        let (x, y, z) = (p[0] - xf, p[1] - yf, p[2] - zf);
        let (u, v, w) = (Self::fade(x), Self::fade(y), Self::fade(z));
        let pm = &self.perm;
        let a = pm[xi] as usize + yi;
        let aa = pm[a] as usize + zi;
        let ab = pm[a + 1] as usize + zi;
        let b = pm[xi + 1] as usize + yi;
        let ba = pm[b] as usize + zi;
        let bb = pm[b + 1] as usize + zi;
        let lerp = |t: f64, a: f64, b: f64| a + t * (b - a);
        lerp(
            w,
            lerp(
                v,
                lerp(
                    u,
                    self.grad(pm[aa], x, y, z),
                    self.grad(pm[ba], x - 1.0, y, z),
                ),
                lerp(
                    u,
                    self.grad(pm[ab], x, y - 1.0, z),
                    self.grad(pm[bb], x - 1.0, y - 1.0, z),
                ),
            ),
            lerp(
                v,
                lerp(
                    u,
                    self.grad(pm[aa + 1], x, y, z - 1.0),
                    self.grad(pm[ba + 1], x - 1.0, y, z - 1.0),
                ),
                lerp(
                    u,
                    self.grad(pm[ab + 1], x, y - 1.0, z - 1.0),
                    self.grad(pm[bb + 1], x - 1.0, y - 1.0, z - 1.0),
                ),
            ),
        )
    }

    /// Fractional Brownian motion on the unit sphere: `octaves` layers starting at
    /// `frequency` in lattice cells per radian along the sphere (a feature of length L on
    /// a body of radius R wants R / L; one cycle round the body wants 1 / 2π), each
    /// `lacunarity` times finer and `gain` times weaker. Normalised so the output amplitude is about 1.
    pub fn fbm(
        &self,
        d: [f64; 3],
        frequency: f64,
        octaves: u32,
        lacunarity: f64,
        gain: f64,
    ) -> f64 {
        let mut sum = 0.0;
        let mut amp = 1.0;
        let mut norm = 0.0;
        let mut freq = frequency;
        // Offset each octave so that octaves do not share zero crossings.
        for o in 0..octaves {
            let off = 17.0 * o as f64;
            sum += amp
                * self.sample([
                    d[0] * freq + off,
                    d[1] * freq + 0.5 * off,
                    d[2] * freq + 0.25 * off,
                ]);
            norm += amp;
            amp *= gain;
            freq *= lacunarity;
        }
        sum / norm
    }

    /// Ridged multifractal variant, in [0, 1], with sharp crests.
    pub fn ridged(
        &self,
        d: [f64; 3],
        frequency: f64,
        octaves: u32,
        lacunarity: f64,
        gain: f64,
    ) -> f64 {
        let mut sum = 0.0;
        let mut amp = 1.0;
        let mut norm = 0.0;
        let mut freq = frequency;
        for o in 0..octaves {
            let off = 31.0 * o as f64;
            let n = 1.0
                - self
                    .sample([
                        d[0] * freq + off,
                        d[1] * freq + 0.3 * off,
                        d[2] * freq + 0.7 * off,
                    ])
                    .abs();
            sum += amp * n * n;
            norm += amp;
            amp *= gain;
            freq *= lacunarity;
        }
        sum / norm
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn noise_is_deterministic_and_bounded() {
        let a = GradientNoise::new(7);
        let b = GradientNoise::new(7);
        let c = GradientNoise::new(8);
        let mut differs = false;
        for i in 0..1000 {
            let t = i as f64 * 0.37;
            let p = [t.sin() * 3.0, t.cos() * 2.0, t * 0.1];
            let va = a.sample(p);
            assert_eq!(va, b.sample(p));
            assert!(va.abs() <= 1.0 + 1e-9);
            if (va - c.sample(p)).abs() > 1e-9 {
                differs = true;
            }
        }
        assert!(differs);
    }

    #[test]
    fn fbm_has_zero_mean_and_unit_scale() {
        let n = GradientNoise::new(3);
        let rng = Rng::new(1);
        let (mut sum, mut sumsq) = (0.0, 0.0);
        let count = 20_000u64;
        for i in 0..count {
            let d = rng.unit_direction(i, 0);
            let v = n.fbm(d, 4.0, 5, 2.0, 0.5);
            sum += v;
            sumsq += v * v;
        }
        let mean = sum / count as f64;
        let rms = (sumsq / count as f64).sqrt();
        assert!(mean.abs() < 0.05, "mean {mean}");
        assert!(rms > 0.1 && rms < 0.6, "rms {rms}");
    }
}
