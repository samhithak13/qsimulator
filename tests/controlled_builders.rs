//! Integration tests for the controlled builders: `cz`, `cu`, `mcx`, `mcu`,
//! and the controlled rotations `crz` / `cp`.

use approx::assert_relative_eq;
use qsimulator::{gates, Circuit};

/// Controlled-Z flips the sign of the |11> branch only, and is symmetric in
/// its two arguments.
#[test]
fn cz_phases_the_11_branch() {
    // H(0), H(1) makes the uniform superposition; CZ negates |11>.
    let mut circuit = Circuit::new(2);
    circuit.h(0).h(1).cz(0, 1);
    let state = circuit.run();

    let s = 0.5; // amplitude magnitude 1/2 for each of four branches
    assert_relative_eq!(state.amplitudes()[0b00].re, s, epsilon = 1e-12);
    assert_relative_eq!(state.amplitudes()[0b01].re, s, epsilon = 1e-12);
    assert_relative_eq!(state.amplitudes()[0b10].re, s, epsilon = 1e-12);
    assert_relative_eq!(state.amplitudes()[0b11].re, -s, epsilon = 1e-12);
}

/// Controlled phase cp(λ) multiplies the |11> branch by e^{iλ} and leaves the
/// other branches alone; cp(π) matches CZ.
#[test]
fn cp_phases_the_11_branch() {
    use std::f64::consts::PI;
    let mut circuit = Circuit::new(2);
    circuit.x(0).x(1).cp(PI / 2.0, 0, 1);
    let state = circuit.run();
    // e^{i·π/2} = i on the |11> amplitude; probability unchanged.
    assert_relative_eq!(state.probability(0b11), 1.0, epsilon = 1e-12);
    assert_relative_eq!(state.amplitudes()[0b11].re, 0.0, epsilon = 1e-12);
    assert_relative_eq!(state.amplitudes()[0b11].im, 1.0, epsilon = 1e-12);

    // cp(π) == cz on the same input.
    let mut viacz = Circuit::new(2);
    viacz.x(0).x(1).cz(0, 1);
    let mut viacp = Circuit::new(2);
    viacp.x(0).x(1).cp(PI, 0, 1);
    assert_relative_eq!(
        (viacz.run().amplitudes()[0b11] - viacp.run().amplitudes()[0b11]).norm(),
        0.0,
        epsilon = 1e-12
    );
}

/// crz only acts when the control is set: with the control at |0>, crz is a
/// no-op on the target.
#[test]
fn crz_requires_control() {
    let mut circuit = Circuit::new(2);
    // control = qubit 0 = |0>; target = qubit 1 = |1>.
    circuit.x(1).crz(1.3, 0, 1);
    let state = circuit.run();
    // rz not applied, so the |1> amplitude on qubit 1 stays real +1.
    assert_relative_eq!(state.amplitudes()[0b10].re, 1.0, epsilon = 1e-12);
    assert_relative_eq!(state.amplitudes()[0b10].im, 0.0, epsilon = 1e-12);
}

/// CZ is symmetric: cz(0, 1) and cz(1, 0) produce the same state.
#[test]
fn cz_is_symmetric() {
    let mut a = Circuit::new(2);
    a.h(0).h(1).cz(0, 1);
    let mut b = Circuit::new(2);
    b.h(0).h(1).cz(1, 0);

    let sa = a.run();
    let sb = b.run();
    for i in 0..4 {
        assert_relative_eq!(
            (sa.amplitudes()[i] - sb.amplitudes()[i]).norm(),
            0.0,
            epsilon = 1e-12
        );
    }
}

/// Controlled-U with an arbitrary gate: a controlled-Y flips and phases the
/// target when the control is set.
#[test]
fn cu_with_y_matches_manual_control() {
    // Control set -> apply Y to target: Y|0> = i|1>.
    let mut circuit = Circuit::new(2);
    circuit.x(0).cu(gates::y(), 0, 1);
    let state = circuit.run();

    // qubit0 = 1 (control), qubit1 = 1 (target flipped) -> index 0b11.
    assert_relative_eq!(state.probability(0b11), 1.0, epsilon = 1e-12);
    assert_relative_eq!(state.amplitudes()[0b11].im, 1.0, epsilon = 1e-12);
}

/// `mcx` with two controls is exactly a Toffoli.
#[test]
fn mcx_two_controls_equals_toffoli() {
    let mut viamcx = Circuit::new(3);
    viamcx.x(0).x(1).mcx(&[0, 1], 2);
    let mut viatoffoli = Circuit::new(3);
    viatoffoli.x(0).x(1).toffoli(0, 1, 2);

    let a = viamcx.run();
    let b = viatoffoli.run();
    for i in 0..8 {
        assert_relative_eq!(
            (a.amplitudes()[i] - b.amplitudes()[i]).norm(),
            0.0,
            epsilon = 1e-12
        );
    }
    // Both controls set -> target flips -> |111>.
    assert_relative_eq!(a.probability(0b111), 1.0, epsilon = 1e-12);
}

/// `mcx` with three controls (a C3X) only flips when all three are set.
#[test]
fn mcx_three_controls_truth_table() {
    // All three controls set: target flips.
    let mut all_set = Circuit::new(4);
    all_set.x(0).x(1).x(2).mcx(&[0, 1, 2], 3);
    assert_relative_eq!(all_set.run().probability(0b1111), 1.0, epsilon = 1e-12);

    // One control missing: target stays |0>.
    let mut one_missing = Circuit::new(4);
    one_missing.x(0).x(1).mcx(&[0, 1, 2], 3);
    assert_relative_eq!(one_missing.run().probability(0b0011), 1.0, epsilon = 1e-12);
}

/// Zero controls makes `mcu` an unconditional gate.
#[test]
fn mcu_zero_controls_is_unconditional() {
    let mut circuit = Circuit::new(1);
    circuit.mcu(gates::x(), &[], 0);
    assert_relative_eq!(circuit.run().probability(0b1), 1.0, epsilon = 1e-12);
}
