// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! Counter-based random numbers (ADR 0004).
//!
//! Every draw is a pure function of a key: the master seed, a stream (a stage or a
//! feature population), an index (a cell or a feature) and a draw counter. No sequential
//! state exists, so results do not depend on thread scheduling or evaluation order, and any
//! tile or feature can be regenerated in isolation.
//!
//! The mixer is the 64-bit finaliser of SplitMix64 applied to a folded key. It is not a
//! cryptographic generator and does not need to be; it passes the statistical checks in
//! the tests below and is stable by construction.

/// 64-bit mix (SplitMix64 finaliser; Steele, Lea and Flood 2014).
#[inline]
pub fn mix64(mut z: u64) -> u64 {
    z = z.wrapping_add(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// Fold several words into one 64-bit hash.
#[inline]
pub fn hash_words(words: &[u64]) -> u64 {
    let mut h = 0x243F_6A88_85A3_08D3u64;
    for &w in words {
        h = mix64(h ^ w);
    }
    mix64(h ^ words.len() as u64)
}

/// Hash a string to 64 bits (used for stage ids and node ids).
pub fn hash_str(s: &str) -> u64 {
    let mut h = 0xB7E1_5162_8AED_2A6Bu64;
    for chunk in s.as_bytes().chunks(8) {
        let mut w = 0u64;
        for (i, &b) in chunk.iter().enumerate() {
            w |= (b as u64) << (8 * i);
        }
        h = mix64(h ^ w);
    }
    mix64(h ^ s.len() as u64)
}

/// Sub-seed for a named stream under a master seed.
pub fn stream_seed(master: u64, stream: &str) -> u64 {
    hash_words(&[master, hash_str(stream)])
}

/// A keyed source of random draws.
/// An index offset for a rejection loop's acceptance draws, so that the test never shares
/// a stream with the candidate it tests (a candidate's longitude drawn from the stream its
/// acceptance reads makes the accepted ones crowd toward small longitudes).
pub const ACCEPT_STREAM: u64 = 0xACCE_7000_0000_0000;

#[derive(Clone, Copy, Debug)]
pub struct Rng {
    seed: u64,
}

impl Rng {
    pub fn new(seed: u64) -> Self {
        Self { seed }
    }

    /// Raw 64 bits for (index, draw).
    #[inline]
    pub fn bits(&self, index: u64, draw: u64) -> u64 {
        hash_words(&[self.seed, index, draw])
    }

    /// Uniform in [0, 1) with 53 bits of resolution.
    #[inline]
    pub fn uniform(&self, index: u64, draw: u64) -> f64 {
        (self.bits(index, draw) >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
    }

    /// Uniform in [lo, hi).
    #[inline]
    pub fn range(&self, index: u64, draw: u64, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.uniform(index, draw)
    }

    /// Standard normal by Box–Muller using draws `draw` and `draw + 1`.
    pub fn normal(&self, index: u64, draw: u64) -> f64 {
        let u1 = self.uniform(index, draw).max(1e-300);
        let u2 = self.uniform(index, draw + 1);
        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }

    /// Uniformly distributed unit direction using draws `draw` and `draw + 1`: `draw + 1`
    /// is the longitude, so no other quantity may take it.
    pub fn unit_direction(&self, index: u64, draw: u64) -> [f64; 3] {
        let z = self.range(index, draw, -1.0, 1.0);
        let phi = self.range(index, draw + 1, 0.0, 2.0 * std::f64::consts::PI);
        let r = (1.0 - z * z).max(0.0).sqrt();
        [r * phi.cos(), r * phi.sin(), z]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn draws_are_pure_functions_of_the_key() {
        let a = Rng::new(42);
        let b = Rng::new(42);
        assert_eq!(a.bits(7, 3), b.bits(7, 3));
        assert_ne!(a.bits(7, 3), a.bits(7, 4));
        assert_ne!(a.bits(7, 3), a.bits(8, 3));
        assert_ne!(Rng::new(43).bits(7, 3), a.bits(7, 3));
    }

    #[test]
    fn uniform_moments_are_sane() {
        let r = Rng::new(1);
        let n = 200_000u64;
        let (mut sum, mut sumsq) = (0.0, 0.0);
        let mut bins = [0usize; 10];
        for i in 0..n {
            let u = r.uniform(i, 0);
            assert!((0.0..1.0).contains(&u));
            sum += u;
            sumsq += u * u;
            bins[(u * 10.0) as usize] += 1;
        }
        let mean = sum / n as f64;
        let var = sumsq / n as f64 - mean * mean;
        assert!((mean - 0.5).abs() < 0.003, "mean {mean}");
        assert!((var - 1.0 / 12.0).abs() < 0.002, "variance {var}");
        for b in bins {
            assert!((b as f64 / n as f64 - 0.1).abs() < 0.005);
        }
    }

    #[test]
    fn normal_moments_are_sane() {
        let r = Rng::new(9);
        let n = 100_000u64;
        let (mut sum, mut sumsq) = (0.0, 0.0);
        for i in 0..n {
            let x = r.normal(i, 0);
            sum += x;
            sumsq += x * x;
        }
        let mean = sum / n as f64;
        let var = sumsq / n as f64 - mean * mean;
        assert!(mean.abs() < 0.01, "mean {mean}");
        assert!((var - 1.0).abs() < 0.02, "variance {var}");
    }

    #[test]
    fn unit_directions_are_isotropic() {
        let r = Rng::new(5);
        let n = 100_000u64;
        let mut m = [0.0f64; 3];
        for i in 0..n {
            let d = r.unit_direction(i, 0);
            let len = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
            assert!((len - 1.0).abs() < 1e-12);
            for k in 0..3 {
                m[k] += d[k];
            }
        }
        for mk in m {
            assert!((mk / n as f64).abs() < 0.01);
        }
    }

    #[test]
    fn stream_seeds_are_stable() {
        assert_eq!(stream_seed(1, "cratering"), stream_seed(1, "cratering"));
        assert_ne!(stream_seed(1, "cratering"), stream_seed(1, "tectonics"));
        assert_ne!(stream_seed(1, "cratering"), stream_seed(2, "cratering"));
    }
}
