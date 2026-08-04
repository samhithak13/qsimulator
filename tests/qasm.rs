//! Integration tests for the OpenQASM 2.0 subset importer.
//!
//! The oracle is the builder API: an imported circuit must produce the same
//! state as the equivalent hand-built circuit.

use approx::assert_relative_eq;
use qsimulator::{qasm, Circuit, State};

/// Assert two states have the same per-basis-state probabilities.
fn assert_same_probs(a: &State, b: &State) {
    assert_eq!(a.n_qubits(), b.n_qubits());
    for i in 0..a.amplitudes().len() {
        assert_relative_eq!(a.probability(i), b.probability(i), epsilon = 1e-12);
    }
}

#[test]
fn imports_bell_state() {
    let src = "\
OPENQASM 2.0;
include \"qelib1.inc\";
qreg q[2];
creg c[2];
h q[0];
cx q[0],q[1];
measure q -> c;
";
    let imported = qasm::parse(src).expect("should parse").run();

    let mut expected = Circuit::new(2);
    expected.h(0).cnot(0, 1);
    assert_same_probs(&imported, &expected.run());
}

#[test]
fn imports_ghz_and_ignores_barrier() {
    let src = "\
OPENQASM 2.0;
qreg q[3];
h q[0];
barrier q;
cx q[0],q[1];
cx q[1],q[2];
";
    let state = qasm::parse(src).expect("should parse").run();
    assert_relative_eq!(state.probability(0b000), 0.5, epsilon = 1e-12);
    assert_relative_eq!(state.probability(0b111), 0.5, epsilon = 1e-12);
}

#[test]
fn imports_rotation_with_pi_angle() {
    // rx(pi) flips |0> to |1> up to global phase.
    let state = qasm::parse("qreg q[1];\nrx(pi) q[0];\n")
        .expect("should parse")
        .run();
    assert_relative_eq!(state.probability(1), 1.0, epsilon = 1e-12);
}

#[test]
fn imports_ccx_toffoli() {
    // Both controls set -> target flips.
    let src = "qreg q[3];\nx q[0];\nx q[1];\nccx q[0],q[1],q[2];\n";
    let state = qasm::parse(src).expect("should parse").run();
    assert_relative_eq!(state.probability(0b111), 1.0, epsilon = 1e-12);
}

#[test]
fn imports_sdg_tdg_and_u1_phase() {
    // x then t then tdg then sdg then s -> net phase +1 on |1> (identity).
    let src = "qreg q[1];\nx q[0];\nt q[0];\ntdg q[0];\nsdg q[0];\ns q[0];\n";
    let imported = qasm::parse(src).expect("should parse").run();
    assert_relative_eq!(imported.amplitudes()[1].re, 1.0, epsilon = 1e-12);
    assert_relative_eq!(imported.amplitudes()[1].im, 0.0, epsilon = 1e-12);

    // u1(pi) is the phase gate; on |1> it matches Z.
    let mut expected = Circuit::new(1);
    expected.x(0).p(std::f64::consts::PI, 0);
    let via_u1 = qasm::parse("qreg q[1];\nx q[0];\nu1(pi) q[0];\n").unwrap();
    assert_same_probs(&via_u1.run(), &expected.run());
    assert_relative_eq!(via_u1.run().amplitudes()[1].re, -1.0, epsilon = 1e-12);

    // `p(lambda)` is accepted as an alias for `u1(lambda)`.
    let via_p = qasm::parse("qreg q[1];\nx q[0];\np(pi) q[0];\n").unwrap();
    assert_same_probs(&via_p.run(), &expected.run());
}

#[test]
fn imports_u2_and_u3() {
    // u3(pi,0,pi) = X, so on |0> gives |1>.
    let s1 = qasm::parse("qreg q[1];\nu3(pi,0,pi) q[0];\n")
        .expect("should parse")
        .run();
    assert_relative_eq!(s1.probability(1), 1.0, epsilon = 1e-12);

    // u2(0,pi) = H, so on |0> gives an equal superposition.
    let s2 = qasm::parse("qreg q[1];\nu2(0,pi) q[0];\n")
        .expect("should parse")
        .run();
    assert_relative_eq!(s2.probability(0), 0.5, epsilon = 1e-12);
    assert_relative_eq!(s2.probability(1), 0.5, epsilon = 1e-12);
}

#[test]
fn error_u3_wrong_angle_count() {
    let err = qasm::parse("qreg q[1];\nu3(pi,0) q[0];\n").unwrap_err();
    assert!(err.contains("3 angle(s)"), "{err}");
}

#[test]
fn imports_cz_and_swap() {
    let src = "qreg q[2];\nx q[0];\nswap q[0],q[1];\ncz q[0],q[1];\n";
    let imported = qasm::parse(src).expect("should parse").run();

    let mut expected = Circuit::new(2);
    expected.x(0).swap(0, 1).cz(0, 1);
    assert_same_probs(&imported, &expected.run());
}

#[test]
fn multiple_registers_map_to_flat_space() {
    // q[0..2] -> 0,1 ; r[0] -> 2. cx q[0],r[0] entangles qubit 0 and qubit 2.
    let src = "qreg q[2];\nqreg r[1];\nh q[0];\ncx q[0],r[0];\n";
    let imported = qasm::parse(src).expect("should parse").run();

    let mut expected = Circuit::new(3);
    expected.h(0).cnot(0, 2);
    assert_same_probs(&imported, &expected.run());
}

#[test]
fn strips_line_and_block_comments() {
    let src = "\
// leading comment
qreg q[1]; /* inline */ x q[0]; // trailing
/* multi
   line */
";
    let state = qasm::parse(src).expect("should parse").run();
    assert_relative_eq!(state.probability(1), 1.0, epsilon = 1e-12);
}

#[test]
fn error_unsupported_gate() {
    let err = qasm::parse("qreg q[1];\nsx q[0];\n").unwrap_err();
    assert!(err.contains("unsupported gate `sx`"), "{err}");
}

#[test]
fn error_unsupported_feature() {
    let err = qasm::parse("qreg q[1];\nreset q[0];\n").unwrap_err();
    assert!(
        err.contains("unsupported OpenQASM feature `reset`"),
        "{err}"
    );
}

#[test]
fn error_no_qreg() {
    let err = qasm::parse("OPENQASM 2.0;\nh q[0];\n").unwrap_err();
    assert!(
        err.contains("no `qreg`") || err.contains("unknown register"),
        "{err}"
    );
}

#[test]
fn error_qubit_out_of_range() {
    let err = qasm::parse("qreg q[2];\nx q[5];\n").unwrap_err();
    assert!(err.contains("out of range"), "{err}");
}

#[test]
fn error_unknown_register() {
    let err = qasm::parse("qreg q[2];\nx r[0];\n").unwrap_err();
    assert!(err.contains("unknown register `r`"), "{err}");
}

#[test]
fn bundled_bell_qasm_matches_builder() {
    let src = include_str!("../programs/bell.qasm");
    let imported = qasm::parse(src)
        .expect("bundled example should parse")
        .run();

    let mut expected = Circuit::new(2);
    expected.h(0).cnot(0, 1);
    assert_same_probs(&imported, &expected.run());
}
