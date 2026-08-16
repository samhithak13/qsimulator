//! Integration tests for the OpenQASM 3 front end.
//!
//! The oracle is the OpenQASM 2 importer: since OpenQASM 3 is normalized into
//! that subset, each test asserts the two spellings of a circuit agree.

use approx::assert_relative_eq;
use qsimulator::{qasm, qasm3};

/// Assert an OpenQASM 3 program and its OpenQASM 2 equivalent give the same
/// state, amplitude by amplitude.
fn assert_same_as_qasm2(three: &str, two: &str) {
    let a = qasm3::parse(three).expect("OpenQASM 3 should parse").run();
    let b = qasm::parse(two).expect("OpenQASM 2 should parse").run();
    assert_eq!(a.n_qubits(), b.n_qubits(), "register sizes differ");
    for (i, (x, y)) in a.amplitudes().iter().zip(b.amplitudes()).enumerate() {
        assert!((x - y).norm() < 1e-12, "amplitude {i}: {x} vs {y}");
    }
}

/// The declaration forms are what differ most between the two languages.
#[test]
fn declarations_translate() {
    assert_same_as_qasm2(
        "OPENQASM 3.0;\ninclude \"stdgates.inc\";\nqubit[2] q;\nh q[0];\ncx q[0], q[1];\n",
        "OPENQASM 2.0;\nqreg q[2];\nh q[0];\ncx q[0],q[1];\n",
    );

    // An unsized declaration is a single qubit.
    assert_same_as_qasm2(
        "OPENQASM 3.0;\nqubit q;\nh q[0];\n",
        "OPENQASM 2.0;\nqreg q[1];\nh q[0];\n",
    );

    // Several registers still share one flat index space, in order.
    assert_same_as_qasm2(
        "OPENQASM 3.0;\nqubit[2] q;\nqubit[1] r;\nh q[0];\ncx q[0], r[0];\n",
        "OPENQASM 2.0;\nqreg q[2];\nqreg r[1];\nh q[0];\ncx q[0],r[0];\n",
    );
}

/// Measurement is an assignment in OpenQASM 3, and the arrow form is still
/// legal; both must reach the same place.
#[test]
fn both_measurement_forms_work() {
    let assigned = "OPENQASM 3.0;\nqubit[2] q;\nbit[2] c;\nh q[0];\nc[0] = measure q[0];\nif (c == 1) {\n  x q[1];\n}\n";
    let arrow = "OPENQASM 3.0;\nqubit[2] q;\nbit[2] c;\nh q[0];\nmeasure q[0] -> c[0];\nif (c == 1) {\n  x q[1];\n}\n";

    for src in [assigned, arrow] {
        let c = qasm3::parse(src).expect("should parse");
        let hist = c.sample(2000, 5);
        // The guarded flip correlates the qubits: only 00 and 11 occur.
        assert_eq!(hist.get(&0b01).copied().unwrap_or(0), 0, "{src}");
        assert_eq!(hist.get(&0b10).copied().unwrap_or(0), 0, "{src}");
    }
}

/// A braced conditional may hold several statements; guarding each one is the
/// same as guarding the block, since a block holds only gates.
#[test]
fn braced_conditional_guards_every_statement() {
    let src = "OPENQASM 3.0;
qubit[3] q;
bit[3] c;
x q[0];
c[0] = measure q[0];
if (c == 1) {
  x q[1];
  x q[2];
}
";
    let state = qasm3::parse(src).expect("should parse").run();
    assert_relative_eq!(state.probability(0b111), 1.0, epsilon = 1e-12);

    // With the condition unmet, neither guarded gate runs.
    let quiet = src.replace("x q[0];\n", "");
    let state = qasm3::parse(&quiet).expect("should parse").run();
    assert_relative_eq!(state.probability(0b000), 1.0, epsilon = 1e-12);
}

/// `gate` declarations read the same in both languages, including the
/// underscore-prefixed argument names Qiskit generates.
#[test]
fn gate_declarations_work() {
    assert_same_as_qasm2(
        "OPENQASM 3.0;
gate bell _gate_q_0, _gate_q_1 {
  h _gate_q_0;
  cx _gate_q_0, _gate_q_1;
}
qubit[2] q;
bell q[0], q[1];
",
        "OPENQASM 2.0;\nqreg q[2];\nh q[0];\ncx q[0],q[1];\n",
    );
}

/// `stdgates.inc` names map onto the same gates, including the ones spelled
/// differently from qelib1.
#[test]
fn stdgates_names_work() {
    assert_same_as_qasm2(
        "OPENQASM 3.0;\nqubit[2] q;\np(0.3) q[0];\ncp(0.4) q[0], q[1];\nsx q[1];\n",
        "OPENQASM 2.0;\nqreg q[2];\nu1(0.3) q[0];\ncu1(0.4) q[0],q[1];\nsx q[1];\n",
    );
}

/// The version header picks the front end, and only a 3 means OpenQASM 3.
#[test]
fn version_detection() {
    assert!(qasm3::is_openqasm3("OPENQASM 3.0;\nqubit[1] q;\n"));
    assert!(qasm3::is_openqasm3("// a comment\nOPENQASM 3;\n"));
    assert!(!qasm3::is_openqasm3("OPENQASM 2.0;\nqreg q[1];\n"));
    assert!(!qasm3::is_openqasm3("qreg q[1];\nh q[0];\n"));
}

/// Constructs with no OpenQASM 2 counterpart are named, not ignored.
#[test]
fn unsupported_constructs_are_reported() {
    for (src, needle) in [
        (
            "OPENQASM 3.0;\nqubit[1] q;\nfor int i in [0:3] { h q[0]; }\n",
            "`for`",
        ),
        (
            "OPENQASM 3.0;\nqubit[1] q;\nwhile (true) { h q[0]; }\n",
            "`while`",
        ),
        (
            "OPENQASM 3.0;\nqubit[1] q;\ndef thing() { h q[0]; }\n",
            "`def`",
        ),
        (
            "OPENQASM 3.0;\nqubit[1] q;\nbit[1] c;\nif (c == 0) { x q[0]; } else { h q[0]; }\n",
            "`else`",
        ),
    ] {
        let err = qasm3::parse(src).unwrap_err().to_string();
        assert!(err.contains(needle), "for `{src}`: {err}");
    }
}

/// A malformed declaration is reported rather than half-understood.
#[test]
fn malformed_declarations_are_reported() {
    for (src, needle) in [
        ("OPENQASM 3.0;\nqubit[2 q;\n", "needs `]`"),
        ("OPENQASM 3.0;\nqubit[2] ;\n", "needs a name"),
        (
            "OPENQASM 3.0;\nqubit[1] q;\nbit[1] c;\nif (c == 0) x q[0];\n",
            "needs a `{",
        ),
    ] {
        let err = qasm3::parse(src).unwrap_err().to_string();
        assert!(err.contains(needle), "for `{src}`: {err}");
    }
}

/// Regression: version detection used to scan raw lines, so a comment whose
/// line began with `OPENQASM 3` sent a perfectly good OpenQASM 2 file to the
/// wrong front end — where its single-statement `if` was then rejected for
/// having no braces.
#[test]
fn a_comment_does_not_pick_the_parser() {
    let src = "\
/*
OPENQASM 3 was the original source language.
This file is the OpenQASM 2 translation.
*/
OPENQASM 2.0;
qreg q[1];
creg c[1];
h q[0];
measure q[0] -> c[0];
if(c==1) x q[0];
";
    assert!(
        !qasm3::is_openqasm3(src),
        "comment text must not select OpenQASM 3"
    );
    qasm::parse(src).expect("a valid OpenQASM 2 file must still parse");

    // And a real OpenQASM 3 header is still found when a comment precedes it.
    let three = "// a leading note\n/* and a block */\nOPENQASM 3.0;\nqubit[1] q;\n";
    assert!(qasm3::is_openqasm3(three));
}
