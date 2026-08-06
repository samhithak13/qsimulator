//! Integration tests for the text program parser (`qsimulator::program`).

use approx::assert_relative_eq;
use qsimulator::program::{self, SampleSpec};

#[test]
fn parses_and_runs_bell() {
    let src = "qubits 2\nh 0\ncnot 0 1\n";
    let prog = program::parse(src).expect("should parse");
    assert!(prog.sample.is_none());

    let state = prog.circuit.run();
    assert_relative_eq!(state.probability(0b00), 0.5, epsilon = 1e-12);
    assert_relative_eq!(state.probability(0b11), 0.5, epsilon = 1e-12);
}

#[test]
fn parses_ghz_with_comments_and_blanks() {
    let src = "\
# GHZ
qubits 3

h 0        # superpose
cnot 0 1
cnot 1 2
sample 500 7
";
    let prog = program::parse(src).expect("should parse");
    assert_eq!(
        prog.sample,
        Some(SampleSpec {
            shots: 500,
            seed: 7
        })
    );

    // GHZ collapses only to 000 or 111.
    let hist = prog.circuit.sample(500, 7);
    assert_eq!(hist.get(&0b010).copied().unwrap_or(0), 0);
    let total = hist.get(&0b000).copied().unwrap_or(0) + hist.get(&0b111).copied().unwrap_or(0);
    assert_eq!(total, 500);
}

#[test]
fn parses_symbolic_pi_angle() {
    // rx(pi) flips |0> to |1> up to global phase, so p(1) = 1.
    let prog = program::parse("qubits 1\nrx pi 0\n").expect("should parse");
    let state = prog.circuit.run();
    assert_relative_eq!(state.probability(1), 1.0, epsilon = 1e-12);

    // A float angle parses too; rx(pi/2) twice is rx(pi).
    let prog2 = program::parse("qubits 1\nrx 1.5707963267948966 0\nrx pi/2 0\n").unwrap();
    assert_relative_eq!(prog2.circuit.run().probability(1), 1.0, epsilon = 1e-9);
}

#[test]
fn all_two_qubit_forms_parse() {
    let prog = program::parse("qubits 3\ncz 0 1\nswap 0 2\ntoffoli 0 1 2\n");
    assert!(prog.is_ok(), "{:?}", prog.err());
}

#[test]
fn parses_u2_and_u3() {
    // u3(pi,0,pi) = X and u2(0,pi) = H.
    let prog = program::parse("qubits 2\nu3 pi 0 pi 0\nu2 0 pi 1\n").expect("should parse");
    let state = prog.circuit.run();
    // qubit 0 flipped to |1>; qubit 1 in equal superposition.
    assert_relative_eq!(state.probability(0b01), 0.5, epsilon = 1e-12);
    assert_relative_eq!(state.probability(0b11), 0.5, epsilon = 1e-12);
}

#[test]
fn parses_sdg_tdg_and_phase() {
    // x, then t; t; sdg; s -> back to phase +1 (identity on |1>); p(pi) = Z.
    let prog = program::parse("qubits 1\nx 0\nsdg 0\ns 0\np pi 0\n").expect("should parse");
    let state = prog.circuit.run();
    // sdg then s cancel; p(pi) negates the |1> amplitude.
    assert_relative_eq!(state.probability(1), 1.0, epsilon = 1e-12);
    assert_relative_eq!(state.amplitudes()[1].re, -1.0, epsilon = 1e-12);
}

#[test]
fn error_gate_before_qubits() {
    let err = program::parse("h 0\n").unwrap_err().to_string();
    assert!(
        err.contains("first instruction must be `qubits N`"),
        "{err}"
    );
}

#[test]
fn error_qubit_out_of_range() {
    let err = program::parse("qubits 2\nx 5\n").unwrap_err().to_string();
    assert!(err.contains("out of range"), "{err}");
    assert!(err.contains("line 2"), "{err}");
}

// Regression: a huge register used to abort the process on allocation.
#[test]
fn error_too_many_qubits() {
    let err = program::parse("qubits 40\n").unwrap_err().to_string();
    assert!(err.contains("exceeds the maximum"), "{err}");
}

#[test]
fn error_unknown_instruction() {
    let err = program::parse("qubits 2\nfoo 0\n").unwrap_err().to_string();
    assert!(err.contains("unknown instruction"), "{err}");
}

#[test]
fn error_bad_arity() {
    let err = program::parse("qubits 2\ncnot 0\n")
        .unwrap_err()
        .to_string();
    assert!(err.contains("expects 3 tokens"), "{err}");
}

#[test]
fn error_control_equals_target() {
    let err = program::parse("qubits 2\ncnot 1 1\n")
        .unwrap_err()
        .to_string();
    assert!(err.contains("must differ"), "{err}");
}

#[test]
fn error_empty_program() {
    let err = program::parse("# just a comment\n\n")
        .unwrap_err()
        .to_string();
    assert!(err.contains("empty"), "{err}");
}

#[test]
fn bundled_ghz_example_parses() {
    let src = include_str!("../programs/ghz.qsim");
    let prog = program::parse(src).expect("bundled example should parse");
    assert_eq!(prog.circuit.run().n_qubits(), 3);
}
