//! Integration tests for the new builder methods: y, s, t, cz, cu, mcx, mcu,
//! and circuit Display.

use approx::assert_relative_eq;
use qsimulator::Circuit;

#[test]
fn y_gate_flips_with_phase() {
    let mut c = Circuit::new(1);
    c.y(0);
    let state = c.run();
    assert_relative_eq!(state.probability(0), 0.0, epsilon = 1e-12);
    assert_relative_eq!(state.probability(1), 1.0, epsilon = 1e-12);
    assert_relative_eq!(state.norm(), 1.0, epsilon = 1e-12);
}

#[test]
fn s_gate_is_phase() {
    // S|+> should give a state where p(0) = p(1) = 0.5 (S only changes phase).
    let mut c = Circuit::new(1);
    c.h(0).s(0);
    let state = c.run();
    assert_relative_eq!(state.probability(0), 0.5, epsilon = 1e-12);
    assert_relative_eq!(state.probability(1), 0.5, epsilon = 1e-12);
}

#[test]
fn t_gate_preserves_probabilities() {
    let mut c = Circuit::new(1);
    c.h(0).t(0);
    let state = c.run();
    assert_relative_eq!(state.probability(0), 0.5, epsilon = 1e-12);
    assert_relative_eq!(state.probability(1), 0.5, epsilon = 1e-12);
}

#[test]
fn s_squared_is_z() {
    // S·S = Z, so H·S·S·H should map |0> to |0> with prob 0 for |1>
    // (H·Z·H = X, so result is |1>).
    let mut c = Circuit::new(1);
    c.h(0).s(0).s(0).h(0);
    let state = c.run();
    assert_relative_eq!(state.probability(1), 1.0, epsilon = 1e-12);
}

#[test]
fn cz_flips_phase_when_both_one() {
    // CZ|11> = -|11>. Since the global phase is unobservable, test via
    // interference: H(0) · CZ(0,1) · H(0) on |11> should flip qubit 0.
    let mut c = Circuit::new(2);
    c.x(0).x(1).h(0).cz(0, 1).h(0);
    let state = c.run();
    // qubit 0 flipped to |0>, qubit 1 stays |1> → |10> = index 0b10.
    assert_relative_eq!(state.probability(0b10), 1.0, epsilon = 1e-12);
}

#[test]
fn cz_no_flip_when_control_zero() {
    let mut c = Circuit::new(2);
    c.x(1).cz(0, 1);
    let state = c.run();
    // Control is |0>, so CZ is identity. State stays |10>.
    assert_relative_eq!(state.probability(0b10), 1.0, epsilon = 1e-12);
}

#[test]
fn cu_with_h_makes_controlled_h() {
    // Controlled-H: if control is |1>, apply H to target.
    // |10> → control set, target |0> → H|0> = |+> on target.
    let mut c = Circuit::new(2);
    c.x(0).cu(qsimulator::gates::h(), 0, 1);
    let state = c.run();
    // qubit 0 = |1>, qubit 1 = |+>: |01> and |11> each with 0.5.
    assert_relative_eq!(state.probability(0b01), 0.5, epsilon = 1e-12);
    assert_relative_eq!(state.probability(0b11), 0.5, epsilon = 1e-12);
}

#[test]
fn mcx_zero_controls_is_unconditional_x() {
    let mut c = Circuit::new(1);
    c.mcx(&[], 0);
    let state = c.run();
    assert_relative_eq!(state.probability(1), 1.0, epsilon = 1e-12);
}

#[test]
fn mcx_one_control_is_cnot() {
    let mut c = Circuit::new(2);
    c.x(0).mcx(&[0], 1);
    let state = c.run();
    assert_relative_eq!(state.probability(0b11), 1.0, epsilon = 1e-12);
}

#[test]
fn mcx_three_controls() {
    // 4-qubit: flip qubit 3 only when qubits 0, 1, 2 are all |1>.
    let mut c = Circuit::new(4);
    c.x(0).x(1).x(2).mcx(&[0, 1, 2], 3);
    let state = c.run();
    assert_relative_eq!(state.probability(0b1111), 1.0, epsilon = 1e-12);

    // Missing one control: target should stay |0>.
    let mut c2 = Circuit::new(4);
    c2.x(0).x(2).mcx(&[0, 1, 2], 3);
    let state2 = c2.run();
    assert_relative_eq!(state2.probability(0b0101), 1.0, epsilon = 1e-12);
}

#[test]
fn mcu_with_z_gate() {
    // MCU with Z and one control = CZ.
    let mut c = Circuit::new(2);
    c.x(0).x(1).h(0).mcu(qsimulator::gates::z(), &[0], 1).h(0);
    let state = c.run();
    // Same as the CZ test: qubit 0 flips.
    assert_relative_eq!(state.probability(0b10), 1.0, epsilon = 1e-12);
}

#[test]
fn circuit_display_smoke() {
    let mut c = Circuit::new(3);
    c.h(0).cnot(0, 1).cnot(0, 2);
    let diagram = format!("{c}");
    assert!(diagram.contains("q0:"));
    assert!(diagram.contains("q1:"));
    assert!(diagram.contains("q2:"));
    assert!(diagram.contains("[H]"));
    assert!(diagram.contains("[X]"));
    assert!(diagram.contains("●"));
}
