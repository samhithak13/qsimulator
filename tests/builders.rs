//! Integration tests for the extended circuit builders: single-qubit `y`,
//! `s`, `t`, controlled-Z (`cz`), arbitrary controlled-U (`cu`), and the
//! generalized multi-controlled `mcx` / `mcu`.

use approx::assert_relative_eq;
use num_complex::Complex64;
use qsimulator::{gates, Circuit};

/// Y|0> = i|1>: probability lands entirely on |1> and the amplitude is +i.
#[test]
fn y_maps_zero_to_i_one() {
    let mut circuit = Circuit::new(1);
    circuit.y(0);
    let state = circuit.run();
    assert_relative_eq!(state.probability(1), 1.0, epsilon = 1e-12);
    assert_relative_eq!(state.amplitudes()[1].re, 0.0, epsilon = 1e-12);
    assert_relative_eq!(state.amplitudes()[1].im, 1.0, epsilon = 1e-12);
}

/// S = diag(1, i): S applied to |1> multiplies its amplitude by i.
#[test]
fn s_phases_one_by_i() {
    let mut circuit = Circuit::new(1);
    circuit.x(0).s(0);
    let state = circuit.run();
    assert_relative_eq!(state.amplitudes()[1].re, 0.0, epsilon = 1e-12);
    assert_relative_eq!(state.amplitudes()[1].im, 1.0, epsilon = 1e-12);
}

/// T = diag(1, e^{i·π/4}): T on |1> applies a π/4 phase, and T·T = S.
#[test]
fn t_squared_is_s() {
    let mut circuit = Circuit::new(1);
    circuit.x(0).t(0).t(0);
    let state = circuit.run();
    // e^{iπ/4} · e^{iπ/4} = e^{iπ/2} = i, matching S on |1>.
    assert_relative_eq!(state.amplitudes()[1].re, 0.0, epsilon = 1e-12);
    assert_relative_eq!(state.amplitudes()[1].im, 1.0, epsilon = 1e-12);
}

/// CZ flips the phase of |11> only, leaving the other basis states alone.
#[test]
fn cz_phases_only_eleven() {
    // Equal superposition over 2 qubits, then CZ. Amplitude of |11> negates.
    let mut circuit = Circuit::new(2);
    circuit.h(0).h(1).cz(0, 1);
    let state = circuit.run();
    let half = 0.5;
    assert_relative_eq!(state.amplitudes()[0b00].re, half, epsilon = 1e-12);
    assert_relative_eq!(state.amplitudes()[0b01].re, half, epsilon = 1e-12);
    assert_relative_eq!(state.amplitudes()[0b10].re, half, epsilon = 1e-12);
    assert_relative_eq!(state.amplitudes()[0b11].re, -half, epsilon = 1e-12);
}

/// CZ is symmetric in its arguments: cz(0,1) == cz(1,0).
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
            sa.amplitudes()[i].re,
            sb.amplitudes()[i].re,
            epsilon = 1e-12
        );
        assert_relative_eq!(
            sa.amplitudes()[i].im,
            sb.amplitudes()[i].im,
            epsilon = 1e-12
        );
    }
}

/// A controlled-U with U = X reproduces CNOT.
#[test]
fn cu_with_x_is_cnot() {
    let mut cu = Circuit::new(2);
    cu.x(0).cu(gates::x(), 0, 1);
    let mut cnot = Circuit::new(2);
    cnot.x(0).cnot(0, 1);
    let a = cu.run();
    let b = cnot.run();
    for i in 0..4 {
        assert_relative_eq!(a.probability(i), b.probability(i), epsilon = 1e-12);
    }
}

/// mcx with zero controls is an unconditional X.
#[test]
fn mcx_zero_controls_is_x() {
    let mut circuit = Circuit::new(1);
    circuit.mcx(&[], 0);
    let state = circuit.run();
    assert_relative_eq!(state.probability(1), 1.0, epsilon = 1e-12);
}

/// mcx with one control is a CNOT.
#[test]
fn mcx_one_control_is_cnot() {
    let mut circuit = Circuit::new(2);
    circuit.x(0).mcx(&[0], 1);
    let state = circuit.run();
    assert_relative_eq!(state.probability(0b11), 1.0, epsilon = 1e-12);
}

/// mcx with three controls flips the target only when all three are set.
#[test]
fn mcx_three_controls() {
    // Only two of the three controls set: target stays |0>.
    let mut partial = Circuit::new(4);
    partial.x(0).x(1).mcx(&[0, 1, 2], 3);
    let sp = partial.run();
    assert_relative_eq!(sp.probability(0b0011), 1.0, epsilon = 1e-12);

    // All three controls set: target flips.
    let mut full = Circuit::new(4);
    full.x(0).x(1).x(2).mcx(&[0, 1, 2], 3);
    let sf = full.run();
    assert_relative_eq!(sf.probability(0b1111), 1.0, epsilon = 1e-12);
}

/// mcu applies an arbitrary gate under multiple controls: a controlled-Z
/// built from mcu phases |111> by -1 when both controls are set.
#[test]
fn mcu_with_z_phases_target() {
    let mut circuit = Circuit::new(3);
    circuit.x(0).x(1).x(2).mcu(gates::z(), &[0, 1], 2);
    let state = circuit.run();
    let amp: Complex64 = state.amplitudes()[0b111];
    assert_relative_eq!(amp.re, -1.0, epsilon = 1e-12);
    assert_relative_eq!(amp.im, 0.0, epsilon = 1e-12);
}
