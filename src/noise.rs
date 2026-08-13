//! Single-qubit noise channels, as sets of Kraus operators.
//!
//! A channel maps a density matrix to `sum_i K_i rho K_i^dagger`, and is
//! physical when it is trace preserving: `sum_i K_i^dagger K_i == I`. This
//! engine stores a state vector rather than a density matrix, so a channel is
//! simulated by the quantum-trajectory method — each shot samples one `K_i`
//! with probability `<psi|K_i^dagger K_i|psi>` and renormalizes
//! ([`State::apply_kraus`](crate::State::apply_kraus)).
//!
//! That means a single run is one sample of a random process, not the average.
//! Averages emerge from [`Circuit::sample`](crate::Circuit::sample), which
//! re-runs the circuit per shot, so noisy results need enough shots to be
//! meaningful — the sampling error falls as `1/sqrt(shots)`.
//!
//! The channels here are the standard single-qubit ones. Each takes a
//! parameter in `0..=1` and returns operators in a fixed order, the identity
//! branch first, so a caller can tell which branch was taken.

use num_complex::Complex64;

use crate::gates;

/// A 2x2 operator. Kraus operators are generally not unitary, so this is the
/// same shape as a gate without the same guarantee.
pub type Kraus = [[Complex64; 2]; 2];

/// Scale every entry of `g` by the real factor `s`.
fn scaled(g: Kraus, s: f64) -> Kraus {
    let f = Complex64::new(s, 0.0);
    [[f * g[0][0], f * g[0][1]], [f * g[1][0], f * g[1][1]]]
}

/// Whether `ops` is trace preserving: `sum_i K_i^dagger K_i == I` within
/// `1e-12`. A channel that fails this is not physical.
pub fn is_trace_preserving(ops: &[Kraus]) -> bool {
    let mut sum = [[Complex64::new(0.0, 0.0); 2]; 2];
    for k in ops {
        for i in 0..2 {
            for j in 0..2 {
                sum[i][j] += k[0][i].conj() * k[0][j] + k[1][i].conj() * k[1][j];
            }
        }
    }
    let identity = [
        [Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0)],
        [Complex64::new(0.0, 0.0), Complex64::new(1.0, 0.0)],
    ];
    (0..2).all(|i| (0..2).all(|j| (sum[i][j] - identity[i][j]).norm() < 1e-12))
}

/// Bit flip: applies X with probability `p`.
pub fn bit_flip(p: f64) -> Vec<Kraus> {
    flip_channel(p, gates::x())
}

/// Phase flip: applies Z with probability `p`.
pub fn phase_flip(p: f64) -> Vec<Kraus> {
    flip_channel(p, gates::z())
}

/// A channel that applies `g` with probability `p` and nothing otherwise.
fn flip_channel(p: f64, g: Kraus) -> Vec<Kraus> {
    let p = p.clamp(0.0, 1.0);
    vec![scaled(gates::id(), (1.0 - p).sqrt()), scaled(g, p.sqrt())]
}

/// Depolarizing: with probability `p` the qubit is replaced by the maximally
/// mixed state, which is the same as applying X, Y or Z each with `p/3`.
///
/// Operators are ordered identity, X, Y, Z.
pub fn depolarizing(p: f64) -> Vec<Kraus> {
    let p = p.clamp(0.0, 1.0);
    let each = (p / 3.0).sqrt();
    vec![
        scaled(gates::id(), (1.0 - p).sqrt()),
        scaled(gates::x(), each),
        scaled(gates::y(), each),
        scaled(gates::z(), each),
    ]
}

/// Amplitude damping: decay of |1> towards |0> with probability `gamma`, the
/// channel behind T1 relaxation.
///
/// Unlike the Pauli channels its operators are not proportional to unitaries,
/// so which branch is taken depends on the state — a qubit already in |0>
/// never decays.
pub fn amplitude_damping(gamma: f64) -> Vec<Kraus> {
    let gamma = gamma.clamp(0.0, 1.0);
    let c = |re: f64| Complex64::new(re, 0.0);
    vec![
        [[c(1.0), c(0.0)], [c(0.0), c((1.0 - gamma).sqrt())]],
        [[c(0.0), c(gamma.sqrt())], [c(0.0), c(0.0)]],
    ]
}

/// Phase damping: loss of coherence without loss of energy, with probability
/// `gamma`. Populations are untouched; only the off-diagonal terms shrink.
pub fn phase_damping(gamma: f64) -> Vec<Kraus> {
    let gamma = gamma.clamp(0.0, 1.0);
    let c = |re: f64| Complex64::new(re, 0.0);
    vec![
        [[c(1.0), c(0.0)], [c(0.0), c((1.0 - gamma).sqrt())]],
        [[c(0.0), c(0.0)], [c(0.0), c(gamma.sqrt())]],
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Every channel must be trace preserving at any strength, including the
    /// endpoints, or it is not a physical map.
    #[test]
    fn channels_are_trace_preserving() {
        for p in [0.0, 0.05, 0.3, 0.5, 0.87, 1.0] {
            for ops in [
                bit_flip(p),
                phase_flip(p),
                depolarizing(p),
                amplitude_damping(p),
                phase_damping(p),
            ] {
                assert!(is_trace_preserving(&ops), "not trace preserving at p={p}");
            }
        }
    }

    /// At strength zero every channel is the identity alone.
    #[test]
    fn zero_strength_is_the_identity_branch() {
        for ops in [
            bit_flip(0.0),
            phase_flip(0.0),
            depolarizing(0.0),
            amplitude_damping(0.0),
            phase_damping(0.0),
        ] {
            let id = gates::id();
            for i in 0..2 {
                for j in 0..2 {
                    assert!((ops[0][i][j] - id[i][j]).norm() < 1e-12);
                }
            }
        }
    }

    /// Out-of-range strengths are clamped rather than producing NaN from a
    /// negative square root.
    #[test]
    fn strength_is_clamped() {
        for ops in [depolarizing(-0.5), amplitude_damping(2.0), bit_flip(1.5)] {
            assert!(is_trace_preserving(&ops));
            assert!(ops.iter().all(|k| k
                .iter()
                .flatten()
                .all(|e| e.re.is_finite() && e.im.is_finite())));
        }
    }

    /// A malformed set is rejected: scaling one operator breaks completeness.
    #[test]
    fn detects_a_non_physical_channel() {
        let mut ops = depolarizing(0.3);
        ops[1] = scaled(ops[1], 2.0);
        assert!(!is_trace_preserving(&ops));
    }
}
