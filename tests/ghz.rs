//! GHZ state preparation and sampling.

use approx::assert_relative_eq;
use qsimulator::Circuit;

/// Prepare a 3-qubit GHZ state: (|000> + |111>) / sqrt(2).
/// H(0), CNOT(0,1), CNOT(0,2).
#[test]
fn ghz_3_qubit_probabilities() {
    let mut circuit = Circuit::new(3);
    circuit.h(0).cnot(0, 1).cnot(0, 2);
    let state = circuit.run();

    assert_relative_eq!(state.probability(0b000), 0.5, epsilon = 1e-12);
    assert_relative_eq!(state.probability(0b111), 0.5, epsilon = 1e-12);
    for i in 1..7 {
        assert_relative_eq!(state.probability(i), 0.0, epsilon = 1e-12);
    }
    assert_relative_eq!(state.norm(), 1.0, epsilon = 1e-12);
}

/// Sampling a GHZ state only ever produces 000 or 111.
#[test]
fn ghz_sampling_only_000_and_111() {
    let mut circuit = Circuit::new(3);
    circuit.h(0).cnot(0, 1).cnot(0, 2);
    let histogram = circuit.sample(2000, 0xDEAD);

    let n000 = histogram.get(&0b000).copied().unwrap_or(0);
    let n111 = histogram.get(&0b111).copied().unwrap_or(0);
    assert_eq!(n000 + n111, 2000);

    // No other outcomes.
    for outcome in 1..7 {
        assert_eq!(histogram.get(&outcome).copied().unwrap_or(0), 0);
    }

    assert!(n000 > 800 && n000 < 1200, "n000 = {n000}");
}

/// 4-qubit GHZ: (|0000> + |1111>) / sqrt(2).
#[test]
fn ghz_4_qubit() {
    let mut circuit = Circuit::new(4);
    circuit.h(0).cnot(0, 1).cnot(0, 2).cnot(0, 3);
    let state = circuit.run();

    assert_relative_eq!(state.probability(0b0000), 0.5, epsilon = 1e-12);
    assert_relative_eq!(state.probability(0b1111), 0.5, epsilon = 1e-12);
    for i in 1..15 {
        assert_relative_eq!(state.probability(i), 0.0, epsilon = 1e-12);
    }
}
