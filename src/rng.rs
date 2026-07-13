//! A tiny, seedable, dependency-free pseudo-random number generator.
//!
//! This is an `xorshift64` generator, which is more than adequate for
//! sampling measurement outcomes in a simulator (it is *not* cryptographic).
//! Seeding is deterministic: the same seed always produces the same stream,
//! which is exactly what reproducible tests and benchmarks need.

/// A seedable `xorshift64` pseudo-random number generator.
pub struct Rng {
    state: u64,
}

impl Rng {
    /// Create a generator from `seed`.
    ///
    /// The seed is first run through the SplitMix64 finalizer so that even
    /// structured seeds (including `0`) yield a well-distributed, non-zero
    /// internal state. A zero state is fatal for xorshift — it would stay
    /// stuck at zero forever — so we also force the low bit set.
    pub fn new(seed: u64) -> Self {
        let mut z = seed.wrapping_add(0x9E37_79B9_7F4A_7C15);
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        let state = (z ^ (z >> 31)) | 1;
        Rng { state }
    }

    /// Draw the next raw 64-bit value from the stream.
    fn next_u64(&mut self) -> u64 {
        let mut x = self.state;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.state = x;
        x
    }

    /// Draw a uniformly distributed `f64` in the half-open interval `[0, 1)`.
    ///
    /// Uses the top 53 bits so every representable value is equally likely.
    pub fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / ((1u64 << 53) as f64)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn same_seed_same_stream() {
        let mut a = Rng::new(42);
        let mut b = Rng::new(42);
        for _ in 0..100 {
            assert_eq!(a.next_u64(), b.next_u64());
        }
    }

    #[test]
    fn different_seeds_differ() {
        let mut a = Rng::new(1);
        let mut b = Rng::new(2);
        // Overwhelmingly likely to differ within a few draws.
        assert_ne!(a.next_u64(), b.next_u64());
    }

    #[test]
    fn seed_zero_is_not_degenerate() {
        let mut r = Rng::new(0);
        let first = r.next_u64();
        assert_ne!(first, 0);
        assert_ne!(first, r.next_u64());
    }

    #[test]
    fn f64_in_unit_interval() {
        let mut r = Rng::new(7);
        for _ in 0..10_000 {
            let x = r.next_f64();
            assert!((0.0..1.0).contains(&x));
        }
    }
}
