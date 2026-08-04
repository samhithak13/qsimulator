//! Integration tests for the `y`, `s`, `t`, `sdg`, `tdg`, and `p` builders.

use approx::assert_relative_eq;
use qsimulator::Circuit;
use std::f64::consts::{FRAC_PI_4, PI};

/// Y|0> = i|1>.
#[test]
fn y_flips_and_phases() {
    let mut circuit = Circuit::new(1);
    circuit.y(0);
    let state = circuit.run();
    assert_relative_eq!(state.probability(0b1), 1.0, epsilon = 1e-12);
    assert_relative_eq!(state.amplitudes()[0b1].im, 1.0, epsilon = 1e-12);
}

/// S = diag(1, i): applied to |1> it multiplies the amplitude by i.
#[test]
fn s_phases_the_one_state() {
    let mut circuit = Circuit::new(1);
    circuit.x(0).s(0);
    let state = circuit.run();
    assert_relative_eq!(state.amplitudes()[0b1].re, 0.0, epsilon = 1e-12);
    assert_relative_eq!(state.amplitudes()[0b1].im, 1.0, epsilon = 1e-12);
    // S = T^2: applying T twice must match S.
    let mut viat = Circuit::new(1);
    viat.x(0).t(0).t(0);
    let s2 = viat.run();
    assert_relative_eq!(
        (state.amplitudes()[0b1] - s2.amplitudes()[0b1]).norm(),
        0.0,
        epsilon = 1e-12
    );
}

/// T = diag(1, e^{i pi/4}): applied to |1> it rotates the phase by pi/4.
#[test]
fn t_rotates_phase_by_pi_over_four() {
    let mut circuit = Circuit::new(1);
    circuit.x(0).t(0);
    let state = circuit.run();
    assert_relative_eq!(state.amplitudes()[0b1].re, FRAC_PI_4.cos(), epsilon = 1e-12);
    assert_relative_eq!(state.amplitudes()[0b1].im, FRAC_PI_4.sin(), epsilon = 1e-12);
}

/// S·S† = identity: on |1>, applying S then S-dagger leaves the phase at +1.
#[test]
fn sdg_undoes_s() {
    let mut circuit = Circuit::new(1);
    circuit.x(0).s(0).sdg(0);
    let state = circuit.run();
    assert_relative_eq!(state.amplitudes()[0b1].re, 1.0, epsilon = 1e-12);
    assert_relative_eq!(state.amplitudes()[0b1].im, 0.0, epsilon = 1e-12);
}

/// T·T† = identity on |1>.
#[test]
fn tdg_undoes_t() {
    let mut circuit = Circuit::new(1);
    circuit.x(0).t(0).tdg(0);
    let state = circuit.run();
    assert_relative_eq!(state.amplitudes()[0b1].re, 1.0, epsilon = 1e-12);
    assert_relative_eq!(state.amplitudes()[0b1].im, 0.0, epsilon = 1e-12);
}

/// Phase gate p(π) acts like Z on |1> (amplitude -> -1), and p only changes
/// phases so basis-state probabilities are untouched.
#[test]
fn p_phase_gate_matches_z_at_pi() {
    let mut circuit = Circuit::new(1);
    circuit.x(0).p(PI, 0);
    let state = circuit.run();
    assert_relative_eq!(state.probability(0b1), 1.0, epsilon = 1e-12);
    assert_relative_eq!(state.amplitudes()[0b1].re, -1.0, epsilon = 1e-12);
    assert_relative_eq!(state.amplitudes()[0b1].im, 0.0, epsilon = 1e-12);
}
