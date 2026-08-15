//! Integration tests for the density-matrix backend.
//!
//! The oracle for unitary evolution is the state vector: the two must agree on
//! every probability. For noise the oracle is analytic, since a density matrix
//! is exact where trajectories only converge.

use qsimulator::density::DensityMatrix;
use qsimulator::{gates, noise, Circuit, DensityError, State};

fn assert_probs_match(rho: &DensityMatrix, state: &State, what: &str) {
    for i in 0..state.amplitudes().len() {
        assert!(
            (rho.probability(i) - state.probability(i)).abs() < 1e-12,
            "{what}: basis {i} differs, {} vs {}",
            rho.probability(i),
            state.probability(i)
        );
    }
}

/// Unitary evolution must match the state vector exactly, for single-qubit,
/// controlled, multi-controlled gates and swap alike.
#[test]
fn unitary_evolution_matches_the_state_vector() {
    let n = 4;
    let mut state = State::new(n);
    let mut rho = DensityMatrix::new(n);

    let steps: Vec<(&str, usize, usize)> = vec![
        ("h", 0, 0),
        ("u3", 1, 0),
        ("cx", 0, 1),
        ("h", 2, 0),
        ("ccx", 0, 3),
        ("swap", 1, 3),
        ("u3", 2, 0),
        ("cx", 2, 3),
    ];
    for (kind, a, b) in steps {
        match kind {
            "h" => {
                state.apply_1q(&gates::h(), a);
                rho.apply_1q(&gates::h(), a);
            }
            "u3" => {
                let g = gates::u3(0.7 + a as f64, -0.3, 1.1);
                state.apply_1q(&g, a);
                rho.apply_1q(&g, a);
            }
            "cx" => {
                state.apply_controlled_1q(&gates::x(), a, b);
                rho.apply_controlled_1q(&gates::x(), a, b);
            }
            "ccx" => {
                state.apply_multi_controlled_1q(&gates::x(), &[a, 1], b);
                rho.apply_multi_controlled_1q(&gates::x(), &[a, 1], b);
            }
            _ => {
                state.swap_qubits(a, b);
                rho.swap_qubits(a, b);
            }
        }
        assert_probs_match(&rho, &state, kind);
    }

    // A pure state stays pure under unitaries, and normalized.
    assert!((rho.trace() - 1.0).abs() < 1e-12, "trace {}", rho.trace());
    assert!(
        (rho.purity() - 1.0).abs() < 1e-12,
        "purity {}",
        rho.purity()
    );
}

/// `from_state` reproduces the state it came from.
#[test]
fn from_state_matches_its_source() {
    let mut c = Circuit::new(3);
    c.h(0).cnot(0, 1).u3(0.4, 0.2, -0.9, 2);
    let state = c.run();
    let rho = DensityMatrix::from_state(&state);
    assert_probs_match(&rho, &state, "from_state");
    assert!((rho.purity() - 1.0).abs() < 1e-12);
}

/// Noise is applied exactly, so depolarizing lands on 2p/3 with no sampling
/// error at all — the reason to pay for the memory.
#[test]
fn noise_is_exact_not_sampled() {
    for p in [0.0, 0.1, 0.37, 0.5, 1.0] {
        let mut rho = DensityMatrix::new(1);
        rho.apply_kraus(&noise::depolarizing(p), 0);
        let expected = 2.0 * p / 3.0;
        assert!(
            (rho.prob_qubit_one(0) - expected).abs() < 1e-12,
            "depolarizing({p}) gave {} not {expected}",
            rho.prob_qubit_one(0)
        );
        assert!((rho.trace() - 1.0).abs() < 1e-12);
    }
}

/// Amplitude damping decays only the excited state, exactly.
#[test]
fn amplitude_damping_is_exact() {
    for gamma in [0.0, 0.25, 0.8, 1.0] {
        let mut excited = DensityMatrix::new(1);
        excited.apply_1q(&gates::x(), 0);
        excited.apply_kraus(&noise::amplitude_damping(gamma), 0);
        assert!(
            (excited.prob_qubit_one(0) - (1.0 - gamma)).abs() < 1e-12,
            "gamma={gamma}: {}",
            excited.prob_qubit_one(0)
        );

        let mut ground = DensityMatrix::new(1);
        ground.apply_kraus(&noise::amplitude_damping(gamma), 0);
        assert!(ground.prob_qubit_one(0).abs() < 1e-12);
    }
}

/// Noise mixes the state, which purity measures directly: full depolarizing
/// takes one qubit to the maximally mixed state, purity 1/2.
#[test]
fn purity_tracks_mixing() {
    let mut pure = DensityMatrix::new(1);
    pure.apply_1q(&gates::h(), 0);
    assert!((pure.purity() - 1.0).abs() < 1e-12);

    // Full depolarization is p = 3/4 in this parameterization, not p = 1: at
    // 3/4 all four Paulis including the identity are equally likely, which is
    // what averages the qubit away to purity 1/2.
    let mut mixed = DensityMatrix::new(1);
    mixed.apply_1q(&gates::h(), 0);
    mixed.apply_kraus(&noise::depolarizing(0.75), 0);
    assert!(
        (mixed.purity() - 0.5).abs() < 1e-12,
        "fully depolarized purity {}",
        mixed.purity()
    );

    // At p = 1 the channel is rho -> (2I - rho)/3, whose eigenvalues are 1/3
    // and 2/3 — mixed, but not maximally so.
    let mut overshot = DensityMatrix::new(1);
    overshot.apply_1q(&gates::h(), 0);
    overshot.apply_kraus(&noise::depolarizing(1.0), 0);
    assert!(
        (overshot.purity() - 5.0 / 9.0).abs() < 1e-12,
        "{}",
        overshot.purity()
    );

    // Partial noise sits strictly between pure and maximally mixed.
    let mut partial = DensityMatrix::new(1);
    partial.apply_1q(&gates::h(), 0);
    partial.apply_kraus(&noise::depolarizing(0.4), 0);
    assert!(partial.purity() > 0.5 && partial.purity() < 1.0);
}

/// An unread measurement is exactly dephasing: populations stay, coherence
/// goes — so the interference that makes H twice an identity is destroyed.
#[test]
fn measure_dephase_erases_only_coherence() {
    let mut rho = DensityMatrix::new(1);
    rho.apply_1q(&gates::h(), 0);
    let before = rho.prob_qubit_one(0);
    rho.measure_dephase(0);
    assert!(
        (rho.prob_qubit_one(0) - before).abs() < 1e-12,
        "population moved"
    );
    assert!((rho.entry(0, 1).norm()) < 1e-12, "coherence survived");

    // H, dephase, H is now a coin flip rather than a return to |0>.
    rho.apply_1q(&gates::h(), 0);
    assert!((rho.prob_qubit_one(0) - 0.5).abs() < 1e-12);
}

/// Reset forces |0> exactly, from any state.
#[test]
fn reset_is_exact() {
    let mut rho = DensityMatrix::new(2);
    rho.apply_1q(&gates::h(), 0);
    rho.apply_controlled_1q(&gates::x(), 0, 1);
    rho.reset(0);
    assert!(rho.prob_qubit_one(0).abs() < 1e-12);
    assert!((rho.trace() - 1.0).abs() < 1e-12);
    // Qubit 1 keeps the population it had, now classically correlated.
    assert!((rho.prob_qubit_one(1) - 0.5).abs() < 1e-12);
}

/// The trajectory backend must converge to what the density matrix says
/// exactly — the two representations of the same physics.
#[test]
fn trajectories_converge_to_the_density_matrix() {
    let mut c = Circuit::new(2);
    c.h(0)
        .cnot(0, 1)
        .depolarizing(0.3, 0)
        .amplitude_damping(0.2, 1);

    let mut rho = DensityMatrix::new(2);
    rho.apply_1q(&gates::h(), 0);
    rho.apply_controlled_1q(&gates::x(), 0, 1);
    rho.apply_kraus(&noise::depolarizing(0.3), 0);
    rho.apply_kraus(&noise::amplitude_damping(0.2), 1);

    let shots = 40_000;
    let hist = c.sample(shots, 17);
    for i in 0..4 {
        let sampled = hist.get(&i).copied().unwrap_or(0) as f64 / shots as f64;
        let exact = rho.probability(i);
        assert!(
            (sampled - exact).abs() < 0.02,
            "basis {i}: sampled {sampled:.4} vs exact {exact:.4}"
        );
    }
}

/// The register ceiling is enforced, since 4^n grows fast enough to matter.
#[test]
#[should_panic(expected = "exceeds the maximum")]
fn oversized_register_is_rejected() {
    DensityMatrix::new(qsimulator::density::MAX_DENSITY_QUBITS + 1);
}

/// `run_density` agrees with sampling, but exactly: the same noisy circuit run
/// both ways must land on the same distribution, one with sampling error and
/// one without.
#[test]
fn run_density_agrees_with_sampling() {
    let mut c = Circuit::new(3);
    c.h(0)
        .cnot(0, 1)
        .depolarizing(0.25, 0)
        .u3(0.6, -0.2, 0.9, 2)
        .cnot(1, 2)
        .amplitude_damping(0.3, 1)
        .phase_damping(0.4, 2);

    let rho = c.run_density().expect("no conditional, so it runs");
    assert!((rho.trace() - 1.0).abs() < 1e-12, "trace {}", rho.trace());
    // Noise mixed the state, so it is no longer pure.
    assert!(rho.purity() < 0.9, "purity {}", rho.purity());

    let shots = 60_000;
    let hist = c.sample(shots, 23);
    for i in 0..8 {
        let sampled = hist.get(&i).copied().unwrap_or(0) as f64 / shots as f64;
        assert!(
            (sampled - rho.probability(i)).abs() < 0.02,
            "basis {i}: sampled {sampled:.4} vs exact {:.4}",
            rho.probability(i)
        );
    }
}

/// A noiseless circuit gives the same probabilities either way, so the two
/// backends agree on the unitary case too.
#[test]
fn run_density_matches_run_without_noise() {
    let mut c = Circuit::new(3);
    c.h(0).cnot(0, 1).u3(0.4, 0.2, -0.9, 2).toffoli(0, 1, 2);

    let rho = c.run_density().unwrap();
    let state = c.run();
    for i in 0..8 {
        assert!((rho.probability(i) - state.probability(i)).abs() < 1e-12);
    }
    assert!((rho.purity() - 1.0).abs() < 1e-12, "should stay pure");
}

/// Feed-forward needs a classical outcome the matrix does not carry, so it is
/// refused rather than approximated.
#[test]
fn feed_forward_is_refused() {
    let mut c = Circuit::new(2);
    c.h(0).measure(0);
    c.if_classical_eq(1, |b| {
        b.x(1);
    });
    assert_eq!(c.run_density(), Err(DensityError::ClassicalFeedForward));
    assert!(c
        .run_density()
        .unwrap_err()
        .to_string()
        .contains("sample it instead"));
}

/// Too large a register is refused with the numbers, rather than attempting a
/// multi-gigabyte allocation.
#[test]
fn oversized_circuit_is_refused() {
    let c = Circuit::new(qsimulator::density::MAX_DENSITY_QUBITS + 2);
    match c.run_density() {
        Err(DensityError::TooManyQubits { qubits, max }) => {
            assert_eq!(qubits, qsimulator::density::MAX_DENSITY_QUBITS + 2);
            assert_eq!(max, qsimulator::density::MAX_DENSITY_QUBITS);
        }
        other => panic!("expected TooManyQubits, got {other:?}"),
    }
}
