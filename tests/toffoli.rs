//! Integration tests for the Toffoli (CCNOT) gate.

use approx::assert_relative_eq;
use qsimulator::Circuit;

/// Full truth table: with controls on qubits 0 and 1 and target on qubit 2,
/// the target flips iff both controls are 1. All other bits are preserved.
#[test]
fn toffoli_truth_table() {
    for c1 in 0..2 {
        for c2 in 0..2 {
            for t in 0..2 {
                let mut circuit = Circuit::new(3);
                if c1 == 1 {
                    circuit.x(0);
                }
                if c2 == 1 {
                    circuit.x(1);
                }
                if t == 1 {
                    circuit.x(2);
                }
                circuit.toffoli(0, 1, 2);
                let state = circuit.run();

                let expected_t = t ^ (c1 & c2);
                let expected_index = c1 | (c2 << 1) | (expected_t << 2);

                assert_relative_eq!(state.probability(expected_index), 1.0, epsilon = 1e-12);
                assert_relative_eq!(state.norm(), 1.0, epsilon = 1e-12);
            }
        }
    }
}

/// Toffoli with only one control satisfied leaves the target unchanged.
#[test]
fn toffoli_needs_both_controls() {
    // Only control qubit 0 is set; target must stay |0>.
    let mut circuit = Circuit::new(3);
    circuit.x(0).toffoli(0, 1, 2);
    let state = circuit.run();
    assert_relative_eq!(state.probability(0b001), 1.0, epsilon = 1e-12);
}

/// Toffoli can create a GHZ-like entangled state from a superposition:
/// H(0), H(1), then CCX(0,1,2) flips qubit 2 only on the |11> branch.
#[test]
fn toffoli_on_superposition_entangles() {
    let mut circuit = Circuit::new(3);
    circuit.h(0).h(1).toffoli(0, 1, 2);
    let state = circuit.run();

    // Four equally weighted branches; the |11> control branch has qubit 2 set.
    assert_relative_eq!(state.probability(0b000), 0.25, epsilon = 1e-12);
    assert_relative_eq!(state.probability(0b001), 0.25, epsilon = 1e-12);
    assert_relative_eq!(state.probability(0b010), 0.25, epsilon = 1e-12);
    assert_relative_eq!(state.probability(0b111), 0.25, epsilon = 1e-12);
    // The un-flipped |011> outcome must have zero weight.
    assert_relative_eq!(state.probability(0b011), 0.0, epsilon = 1e-12);
}
