//! Integration tests for measurement, collapse, and sampling.

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
