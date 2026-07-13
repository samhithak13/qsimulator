//! Integration tests for the Rx/Ry/Rz rotation-gate builder methods.

use approx::assert_relative_eq;
use qsimulator::Circuit;
use std::f64::consts::PI;

/// Rx(π) acts like X on the state: |0> -> (phase)|1>, so p(1) = 1.
#[test]
fn rx_pi_flips_qubit() {
    let mut circuit = Circuit::new(1);
    circuit.rx(PI, 0);
    let state = circuit.run();
    assert_relative_eq!(state.probability(1), 1.0, epsilon = 1e-12);
    assert_relative_eq!(state.norm(), 1.0, epsilon = 1e-12);
}

/// Ry(π/2) on |0> produces an equal superposition (the |+>-like state).
#[test]
fn ry_half_pi_makes_equal_superposition() {
    let mut circuit = Circuit::new(1);
    circuit.ry(PI / 2.0, 0);
    let state = circuit.run();
    assert_relative_eq!(state.probability(0), 0.5, epsilon = 1e-12);
    assert_relative_eq!(state.probability(1), 0.5, epsilon = 1e-12);
}

/// Rz only changes phases, so it leaves computational-basis probabilities
/// untouched.
#[test]
fn rz_preserves_basis_probabilities() {
    let mut circuit = Circuit::new(1);
    circuit.x(0).rz(0.9, 0);
    let state = circuit.run();
    assert_relative_eq!(state.probability(0), 0.0, epsilon = 1e-12);
    assert_relative_eq!(state.probability(1), 1.0, epsilon = 1e-12);
}

/// Two half-rotations compose into a full rotation: Rx(π/2)·Rx(π/2) = Rx(π).
#[test]
fn rx_half_rotations_compose() {
    let mut circuit = Circuit::new(1);
    circuit.rx(PI / 2.0, 0).rx(PI / 2.0, 0);
    let state = circuit.run();
    assert_relative_eq!(state.probability(1), 1.0, epsilon = 1e-12);
}
