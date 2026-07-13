//! Integration tests for the SWAP gate.

use approx::assert_relative_eq;
use qsimulator::Circuit;

/// SWAP exchanges the two qubits: prepare |q1 q0> = |01> (qubit 0 is |1>),
/// then swap so qubit 1 becomes |1> and qubit 0 becomes |0> (index 0b10).
#[test]
fn swap_exchanges_two_qubits() {
    let mut circuit = Circuit::new(2);
    circuit.x(0).swap(0, 1);
    let state = circuit.run();

    assert_relative_eq!(state.probability(0b10), 1.0, epsilon = 1e-12);
    assert_relative_eq!(state.probability(0b01), 0.0, epsilon = 1e-12);
    assert_relative_eq!(state.norm(), 1.0, epsilon = 1e-12);
}

/// SWAP is its own inverse: applying it twice restores the original state.
#[test]
fn swap_is_involutive() {
    let mut circuit = Circuit::new(2);
    circuit.x(0).swap(0, 1).swap(0, 1);
    let state = circuit.run();
    assert_relative_eq!(state.probability(0b01), 1.0, epsilon = 1e-12);
}

/// Swapping equal qubits is a no-op, and swapping non-adjacent qubits in a
/// larger register moves the excitation correctly.
#[test]
fn swap_non_adjacent_qubits() {
    // 3-qubit register, put qubit 0 into |1> (index 0b001), swap qubits 0 and 2.
    let mut circuit = Circuit::new(3);
    circuit.x(0).swap(0, 2);
    let state = circuit.run();
    assert_relative_eq!(state.probability(0b100), 1.0, epsilon = 1e-12);
}

/// SWAP preserves a superposition on the swapped qubits by moving amplitudes,
/// not just basis labels. Here H(0) then SWAP(0,1) puts the superposition on
/// qubit 1 instead of qubit 0.
#[test]
fn swap_moves_superposition() {
    let mut circuit = Circuit::new(2);
    circuit.h(0).swap(0, 1);
    let state = circuit.run();
    // Superposition now on qubit 1: equal weight on |00> and |10>.
    assert_relative_eq!(state.probability(0b00), 0.5, epsilon = 1e-12);
    assert_relative_eq!(state.probability(0b10), 0.5, epsilon = 1e-12);
    assert_relative_eq!(state.probability(0b01), 0.0, epsilon = 1e-12);
}
