//! Integration tests for measurement, collapse, and sampling.

use approx::assert_relative_eq;
use qsimulator::{Circuit, Rng, State};

/// Measuring a fresh |0...0> register is deterministic: every qubit reads 0.
#[test]
fn basis_state_zero_is_deterministic() {
    let mut rng = Rng::new(123);
    let mut state = State::new(3);
    assert_eq!(state.measure_all(&mut rng), 0);

    // measure_qubit on |0> must also always return false.
    for q in 0..3 {
        let mut s = State::new(3);
        assert!(!s.measure_qubit(q, &mut rng));
    }
}

/// Applying X then measuring yields |1> with certainty.
#[test]
fn x_then_measure_is_one() {
    let mut rng = Rng::new(7);
    let mut circuit = Circuit::new(1);
    circuit.x(0);
    let mut state = circuit.run();
    assert!(state.measure_qubit(0, &mut rng));
    assert_eq!(state.norm(), 1.0);
}

/// X on qubit 1 of a 3-qubit register lands on basis index 0b010 every time.
#[test]
fn x_on_middle_qubit_measures_to_that_bit() {
    let mut rng = Rng::new(99);
    let mut circuit = Circuit::new(3);
    circuit.x(1);
    for _ in 0..50 {
        let outcome = circuit.run().measure_all(&mut rng);
        assert_eq!(outcome, 0b010);
    }
}

/// A Bell state only ever collapses to 00 or 11 — never 01 or 10.
#[test]
fn bell_sampling_only_hits_00_and_11() {
    let mut circuit = Circuit::new(2);
    circuit.h(0).cnot(0, 1);
    let histogram = circuit.sample(2000, 0xABCD);

    assert_eq!(histogram.get(&0b01).copied().unwrap_or(0), 0);
    assert_eq!(histogram.get(&0b10).copied().unwrap_or(0), 0);

    let n00 = histogram.get(&0b00).copied().unwrap_or(0);
    let n11 = histogram.get(&0b11).copied().unwrap_or(0);
    assert_eq!(n00 + n11, 2000);
    // Roughly balanced; generous bounds keep this from being flaky.
    assert!(n00 > 800 && n00 < 1200, "n00 = {n00}");
    assert!(n11 > 800 && n11 < 1200, "n11 = {n11}");
}

/// Measuring one half of a Bell pair determines the other half.
#[test]
fn single_qubit_collapse_correlates_the_pair() {
    let mut rng = Rng::new(2024);
    for _ in 0..100 {
        let mut circuit = Circuit::new(2);
        circuit.h(0).cnot(0, 1);
        let mut state = circuit.run();

        let first = state.measure_qubit(0, &mut rng);
        let second = state.measure_qubit(1, &mut rng);
        assert_eq!(first, second, "Bell pair qubits must agree after collapse");
        assert!((state.norm() - 1.0).abs() < 1e-12);
    }
}

/// Identical seeds reproduce identical histograms; a different seed generally
/// does not.
#[test]
fn same_seed_is_reproducible() {
    let mut circuit = Circuit::new(3);
    circuit.h(0).h(1).h(2);

    let a = circuit.sample(500, 42);
    let b = circuit.sample(500, 42);
    assert_eq!(a, b);

    let c = circuit.sample(500, 43);
    assert_ne!(a, c);
}

/// Regression: a mid-circuit measurement used to be ignored on import, so
/// `h; measure; h` collapsed into `h; h` — the identity — and reported |0>
/// with certainty. Measuring first makes the second H act on a definite state,
/// so the real answer is a coin flip.
#[test]
fn mid_circuit_measurement_is_not_a_no_op() {
    let mut c = Circuit::new(1);
    c.h(0).measure(0).h(0);

    let hist = c.sample(4000, 7);
    let zeros = hist.get(&0).copied().unwrap_or(0);
    let ones = hist.get(&1).copied().unwrap_or(0);
    assert_eq!(zeros + ones, 4000);
    // Both outcomes must be well represented; without the collapse this was
    // 4000/0. Wide bounds keep it robust — the point is "not certain".
    assert!(
        (1700..=2300).contains(&zeros),
        "expected ~2000 zeros, got {zeros}"
    );

    // Without the measurement the two H gates really do cancel.
    let mut plain = Circuit::new(1);
    plain.h(0).h(0);
    assert_eq!(plain.sample(4000, 7).get(&0).copied().unwrap_or(0), 4000);
}

/// Measuring collapses the register: a gate after it sees a definite state,
/// not a superposition. The trailing `id` is what makes the measurement
/// mid-circuit rather than readout.
#[test]
fn measurement_collapses_the_register() {
    let mut c = Circuit::new(1);
    c.h(0).measure(0).id(0);
    let p1 = c.run().prob_qubit_one(0);
    assert!(
        !(1e-12..=1.0 - 1e-12).contains(&p1),
        "qubit should be collapsed, got p(1) = {p1}"
    );
}

/// A collapse propagates: measuring one half of what would become a Bell pair
/// still leaves the two qubits perfectly correlated, never anti-correlated.
#[test]
fn collapse_propagates_to_later_gates() {
    let mut c = Circuit::new(2);
    c.h(0).measure(0).cnot(0, 1);
    let hist = c.sample(2000, 3);
    assert_eq!(hist.get(&0b01).copied().unwrap_or(0), 0);
    assert_eq!(hist.get(&0b10).copied().unwrap_or(0), 0);
    assert_eq!(
        hist.get(&0b00).copied().unwrap_or(0) + hist.get(&0b11).copied().unwrap_or(0),
        2000
    );
}

/// Sampling stays a pure function of the seed even when the circuit branches,
/// and different seeds really do explore different branches.
#[test]
fn stochastic_sampling_is_reproducible() {
    let mut c = Circuit::new(1);
    c.h(0).measure(0).h(0);
    assert_eq!(c.sample(500, 11), c.sample(500, 11));
    assert_ne!(c.sample(500, 11), c.sample(500, 12));
}

/// `run_seeded` picks the collapse stream; with no measurement the seed is
/// irrelevant and every seed gives the same state.
#[test]
fn run_seeded_only_matters_when_measuring() {
    let mut stochastic = Circuit::new(1);
    stochastic.h(0).measure(0).id(0);
    let outcomes: Vec<f64> = (0..16)
        .map(|s| stochastic.run_seeded(s).prob_qubit_one(0))
        .collect();
    assert!(
        outcomes.iter().any(|p| *p > 0.5) && outcomes.iter().any(|p| *p < 0.5),
        "seeds should reach both branches, got {outcomes:?}"
    );

    let mut unitary = Circuit::new(2);
    unitary.h(0).cnot(0, 1);
    for s in 0..4 {
        assert_relative_eq!(
            unitary.run_seeded(s).probability(0b11),
            unitary.run().probability(0b11),
            epsilon = 1e-12
        );
    }
}

/// A measurement at the end of a circuit is readout: nothing follows it, so
/// collapsing would only discard the prepared state and report one arbitrary
/// branch. `run` therefore reports the superposition, and stays deterministic
/// — which matters because nearly every written-out program ends this way.
#[test]
fn trailing_measurement_is_readout_not_collapse() {
    let mut c = Circuit::new(2);
    c.h(0).cnot(0, 1).measure(0).measure(1);

    let state = c.run();
    assert_relative_eq!(state.probability(0b00), 0.5, epsilon = 1e-12);
    assert_relative_eq!(state.probability(0b11), 0.5, epsilon = 1e-12);

    // Deterministic across seeds, unlike a genuinely mid-circuit measurement.
    for seed in 0..4 {
        assert_relative_eq!(c.run_seeded(seed).probability(0b00), 0.5, epsilon = 1e-12);
    }

    // And sampling it matches the same circuit without the readout.
    let mut plain = Circuit::new(2);
    plain.h(0).cnot(0, 1);
    assert_eq!(c.sample(500, 4), plain.sample(500, 4));
}

/// `reset` forces a qubit to |0> whatever state it was in — including a
/// superposition, which it collapses first.
#[test]
fn reset_forces_zero() {
    let mut c = Circuit::new(1);
    c.h(0).reset(0);
    assert_relative_eq!(c.run().probability(0), 1.0, epsilon = 1e-12);

    let mut from_one = Circuit::new(1);
    from_one.x(0).reset(0);
    assert_relative_eq!(from_one.run().probability(0), 1.0, epsilon = 1e-12);
}

/// Resetting half of a Bell pair collapses the pair first, so the partner is
/// left on whichever branch came up — |00> or |01>, never a superposition of
/// the two and never |11>.
#[test]
fn reset_leaves_the_partner_on_the_collapsed_branch() {
    let mut c = Circuit::new(2);
    c.h(0).cnot(0, 1).reset(0);

    let hist = c.sample(2000, 5);
    // Index bit 0 is qubit 0, which the reset drove to |0>.
    assert_eq!(
        hist.get(&0b01).copied().unwrap_or(0),
        0,
        "qubit 0 must be 0"
    );
    assert_eq!(
        hist.get(&0b11).copied().unwrap_or(0),
        0,
        "qubit 0 must be 0"
    );
    let zero = hist.get(&0b00).copied().unwrap_or(0);
    let two = hist.get(&0b10).copied().unwrap_or(0);
    assert_eq!(zero + two, 2000);
    assert!((800..=1200).contains(&zero), "expected ~1000, got {zero}");
}

/// Unlike a trailing measurement, a trailing reset is not readout: it changes
/// the state that gets sampled, so it always applies.
#[test]
fn trailing_reset_still_applies() {
    let mut c = Circuit::new(1);
    c.h(0).reset(0);
    // If this were treated as readout the state would still be a superposition.
    assert_relative_eq!(c.run().probability(0), 1.0, epsilon = 1e-12);
    assert_eq!(c.sample(300, 1).get(&0).copied().unwrap_or(0), 300);
}

/// A conditional runs only when the classical register matches. Measuring
/// qubit 0 into bit 0 and guarding on `c == 1` correlates the two qubits.
#[test]
fn conditional_fires_on_the_measured_value() {
    let mut c = Circuit::new(2);
    c.h(0).measure(0);
    c.if_classical_eq(1, |b| {
        b.x(1);
    });
    c.measure(1);

    let hist = c.sample(2000, 4);
    assert_eq!(hist.get(&0b01).copied().unwrap_or(0), 0);
    assert_eq!(hist.get(&0b10).copied().unwrap_or(0), 0);
    let both = hist.get(&0b00).copied().unwrap_or(0) + hist.get(&0b11).copied().unwrap_or(0);
    assert_eq!(both, 2000);
}

/// Bits that were never measured read 0, so a guard on 0 fires and any other
/// guard does not.
#[test]
fn unmeasured_classical_bits_read_zero() {
    let mut fires = Circuit::new(1);
    fires.if_classical_eq(0, |b| {
        b.x(0);
    });
    assert_relative_eq!(fires.run().probability(1), 1.0, epsilon = 1e-12);

    let mut quiet = Circuit::new(1);
    quiet.if_classical_eq(1, |b| {
        b.x(0);
    });
    assert_relative_eq!(quiet.run().probability(0), 1.0, epsilon = 1e-12);
}

/// The whole register is compared, so a guard can depend on several bits.
#[test]
fn conditional_compares_the_whole_register() {
    let mut c = Circuit::new(3);
    c.x(0).x(1).measure(0).measure(1);
    c.if_classical_eq(0b11, |b| {
        b.x(2);
    });
    assert_relative_eq!(c.run().probability(0b111), 1.0, epsilon = 1e-12);
}

/// `measure_into` sends the outcome to a chosen bit, which is what an imported
/// program needs when it compacts measurements into low bits.
#[test]
fn measure_into_targets_a_chosen_bit() {
    let mut c = Circuit::new(2);
    c.x(1).measure_into(1, 0); // qubit 1 -> classical bit 0
    c.if_classical_eq(1, |b| {
        b.x(0);
    });
    assert_relative_eq!(c.run().probability(0b11), 1.0, epsilon = 1e-12);
}

/// A conditional block may only hold gates: anything else would change the
/// value being tested part-way through, or have no OpenQASM form. Dropping it
/// silently would simulate something other than what was written.
#[test]
#[should_panic(expected = "may only contain gates")]
fn conditional_block_rejects_a_measurement() {
    let mut c = Circuit::new(1);
    c.if_classical_eq(0, |b| {
        b.measure(0);
    });
}

#[test]
#[should_panic(expected = "nested conditional")]
fn conditional_block_rejects_nesting() {
    let mut c = Circuit::new(1);
    c.if_classical_eq(0, |b| {
        b.if_classical_eq(0, |inner| {
            inner.x(0);
        });
    });
}
