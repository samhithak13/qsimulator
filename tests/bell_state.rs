use approx::assert_relative_eq;
use qsimulator::Circuit;

#[test]
fn bell_state_probabilities() {
    let mut circuit = Circuit::new(2);
    circuit.h(0).cnot(0, 1);
    let state = circuit.run();

    // Expect (|00> + |11>)/sqrt(2): p(00) = p(11) = 0.5, others 0.
    assert_relative_eq!(state.probability(0b00), 0.5, epsilon = 1e-9);
    assert_relative_eq!(state.probability(0b01), 0.0, epsilon = 1e-9);
    assert_relative_eq!(state.probability(0b10), 0.0, epsilon = 1e-9);
    assert_relative_eq!(state.probability(0b11), 0.5, epsilon = 1e-9);
    assert_relative_eq!(state.norm(), 1.0, epsilon = 1e-9);
}

#[test]
fn x_gate_flips_qubit() {
    let mut circuit = Circuit::new(1);
    circuit.x(0);
    let state = circuit.run();
    assert_relative_eq!(state.probability(1), 1.0, epsilon = 1e-9);
}
