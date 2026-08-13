//! Integration tests for noise channels.
//!
//! A channel is simulated by sampling one Kraus operator per shot, so the
//! oracle is the analytic ensemble average: each test compares a sampled
//! distribution against the value the channel is defined to produce.

use qsimulator::{noise, Circuit};

/// Fraction of `shots` in which qubit `q` read |1>.
fn fraction_one(c: &Circuit, q: usize, shots: usize, seed: u64) -> f64 {
    let hist = c.sample(shots, seed);
    let ones: usize = hist
        .iter()
        .filter(|(state, _)| *state >> q & 1 == 1)
        .map(|(_, n)| n)
        .sum();
    ones as f64 / shots as f64
}

/// A bit flip with probability `p` on |0> leaves the qubit reading |1> exactly
/// `p` of the time.
#[test]
fn bit_flip_matches_its_probability() {
    for p in [0.0, 0.1, 0.5, 0.9, 1.0] {
        let mut c = Circuit::new(1);
        c.bit_flip(p, 0);
        let measured = fraction_one(&c, 0, 20_000, 7);
        assert!((measured - p).abs() < 0.02, "bit_flip({p}) gave {measured}");
    }
}

/// Depolarizing sends |0> to |1> whenever it picks X or Y, which is `2p/3`.
#[test]
fn depolarizing_matches_its_analytic_rate() {
    for p in [0.0, 0.3, 0.75, 1.0] {
        let mut c = Circuit::new(1);
        c.depolarizing(p, 0);
        let measured = fraction_one(&c, 0, 20_000, 11);
        let expected = 2.0 * p / 3.0;
        assert!(
            (measured - expected).abs() < 0.02,
            "depolarizing({p}) gave {measured}, expected {expected}"
        );
    }
}

/// Amplitude damping decays |1> to |0> with probability gamma, and leaves |0>
/// alone — the asymmetry that makes it more than a Pauli channel.
#[test]
fn amplitude_damping_only_decays_the_excited_state() {
    for gamma in [0.0, 0.25, 0.8, 1.0] {
        let mut excited = Circuit::new(1);
        excited.x(0).amplitude_damping(gamma, 0);
        let survived = fraction_one(&excited, 0, 20_000, 3);
        assert!(
            (survived - (1.0 - gamma)).abs() < 0.02,
            "from |1> with gamma={gamma}: {survived} survived"
        );

        // From the ground state there is nothing to lose.
        let mut ground = Circuit::new(1);
        ground.amplitude_damping(gamma, 0);
        assert_eq!(fraction_one(&ground, 0, 2_000, 3), 0.0);
    }
}

/// Phase damping destroys coherence without moving population: the qubit still
/// reads |1> half the time after H, but the interference that H;H would undo
/// is gone.
#[test]
fn phase_damping_kills_coherence_not_population() {
    let mut populations = Circuit::new(1);
    populations.h(0).phase_damping(0.5, 0);
    let ones = fraction_one(&populations, 0, 20_000, 5);
    assert!((ones - 0.5).abs() < 0.02, "population moved: {ones}");

    // H, full dephasing, H would return to |0> with no noise; complete
    // dephasing leaves it a coin flip instead.
    let mut interference = Circuit::new(1);
    interference.h(0).phase_damping(1.0, 0).h(0);
    let ones = fraction_one(&interference, 0, 20_000, 5);
    assert!((ones - 0.5).abs() < 0.02, "coherence survived: {ones}");

    // With no dephasing the two Hadamards cancel exactly.
    let mut clean = Circuit::new(1);
    clean.h(0).phase_damping(0.0, 0).h(0);
    assert_eq!(fraction_one(&clean, 0, 2_000, 5), 0.0);
}

/// Noise on one half of a Bell pair breaks the correlation at the rate the
/// channel prescribes, which a single-qubit test could not show.
#[test]
fn noise_breaks_entangled_correlation() {
    let mut c = Circuit::new(2);
    c.h(0).cnot(0, 1).bit_flip(0.25, 1);

    let hist = c.sample(20_000, 9);
    let disagreeing: usize = hist
        .iter()
        .filter(|(state, _)| (*state & 1) != (*state >> 1 & 1))
        .map(|(_, n)| n)
        .sum();
    let rate = disagreeing as f64 / 20_000.0;
    assert!((rate - 0.25).abs() < 0.02, "disagreement rate {rate}");
}

/// A custom channel goes through the same path as the named ones.
#[test]
fn a_custom_channel_runs() {
    let mut c = Circuit::new(1);
    c.channel(noise::bit_flip(0.5), 0);
    let measured = fraction_one(&c, 0, 20_000, 2);
    assert!((measured - 0.5).abs() < 0.02, "{measured}");
}

/// A channel that is not trace preserving is rejected rather than quietly
/// changing the total probability.
#[test]
#[should_panic(expected = "not trace preserving")]
fn non_physical_channel_is_rejected() {
    let mut ops = noise::depolarizing(0.3);
    ops.truncate(2);
    Circuit::new(1).channel(ops, 0);
}
