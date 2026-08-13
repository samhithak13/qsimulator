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
fn parses_cu3_as_controlled_x() {
    // cu3(pi, 0, pi) is a controlled-X, so |01> -> |11>.
    let prog = program::parse("qubits 2\nx 0\ncu3 pi 0 pi 0 1\n").expect("should parse");
    assert_relative_eq!(prog.circuit.run().probability(0b11), 1.0, epsilon = 1e-12);
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

/// Every native instruction parses (exercises each builder arm).
#[test]
fn every_native_instruction_parses() {
    let src = "\
qubits 3
id 0
h 0
x 0
y 1
z 2
s 0
t 1
sdg 2
tdg 0
rx pi/2 0
ry -pi/4 1
rz 0.3 2
p pi 0
u2 0.1 0.2 1
u3 0.4 0.5 0.6 2
cnot 0 1
cy 1 2
cz 0 2
ch 0 1
crz 0.7 1 2
cp pi/2 0 1
cu3 0.1 0.2 0.3 0 2
swap 0 2
toffoli 0 1 2
cswap 0 1 2
measure 0
reset 1
if 1 x 2
mcx 0 1 2
mcu3 0.1 0.2 0.3 0 1 2
sample 100 42
";
    let prog = program::parse(src).expect("should parse");
    assert!(prog.sample.is_some());
    assert!((prog.circuit.run().norm() - 1.0).abs() < 1e-9);
}

/// `mcx` takes any number of controls: with all three controls set it flips
/// the target, and it leaves every other basis state alone.
#[test]
fn mcx_flips_only_the_all_controls_set_state() {
    let prog = program::parse("qubits 4\nh 0\nh 1\nh 2\nmcx 0 1 2 3\n").expect("should parse");
    let state = prog.circuit.run();
    // |0111> (controls set, target clear) is the one state that moved, to
    // |1111>; the other seven keep their target at 0.
    assert_relative_eq!(state.probability(0b0111), 0.0, epsilon = 1e-12);
    assert_relative_eq!(state.probability(0b1111), 0.125, epsilon = 1e-12);
    for i in [0b0000, 0b0001, 0b0010, 0b0011, 0b0100, 0b0101, 0b0110] {
        assert_relative_eq!(state.probability(i), 0.125, epsilon = 1e-12);
    }
}

/// `mcu3` with theta = 0 is a multi-controlled phase, so it flips the sign of
/// the all-ones amplitude and leaves the probabilities untouched.
#[test]
fn mcu3_applies_a_multi_controlled_phase() {
    let prog = program::parse("qubits 3\nh 0\nh 1\nh 2\nmcu3 0 0 pi 0 1 2\n").expect("parse");
    let amps = prog.circuit.run().amplitudes().to_vec();
    for (i, a) in amps.iter().enumerate() {
        let want = if i == 0b111 {
            -0.125f64.sqrt()
        } else {
            0.125f64.sqrt()
        };
        assert_relative_eq!(a.re, want, epsilon = 1e-12);
        assert_relative_eq!(a.im, 0.0, epsilon = 1e-12);
    }
}

/// `if VALUE INSTRUCTION` guards an instruction on the classical register, the
/// native spelling of OpenQASM's `if (c == VALUE)`.
#[test]
fn if_guards_an_instruction() {
    let prog = program::parse("qubits 2\nh 0\nmeasure 0\nif 1 x 1\nmeasure 1\n").expect("parse");
    let hist = prog.circuit.sample(2000, 4);
    // The guarded flip correlates the qubits: only 00 and 11 occur.
    assert_eq!(hist.get(&0b01).copied().unwrap_or(0), 0);
    assert_eq!(hist.get(&0b10).copied().unwrap_or(0), 0);
    assert_eq!(
        hist.get(&0b00).copied().unwrap_or(0) + hist.get(&0b11).copied().unwrap_or(0),
        2000
    );

    // A guard that never matches leaves the circuit alone.
    let quiet = program::parse("qubits 1\nif 1 x 0\n").expect("parse");
    assert_relative_eq!(quiet.circuit.run().probability(0), 1.0, epsilon = 1e-12);
}

#[test]
fn native_error_paths() {
    let cases = [
        ("qubits 1\nqubits 2\n", "may only appear once"),
        ("qubits 0\n", "must be >= 1"),
        ("qubits 2\ncrz 0.5 1 1\n", "must differ"),
        ("qubits 3\ncswap 0 1 1\n", "must be distinct"),
        ("qubits 3\ntoffoli 0 1 1\n", "must differ"),
        ("qubits 3\nmcx 0 1 1\n", "is repeated"),
        ("qubits 3\nmcx 2\n", "at least 3 tokens"),
        ("qubits 3\nmcu3 0.1 0.2 0.3 2\n", "at least 6 tokens"),
        ("qubits 3\nmcu3 nope 0.2 0.3 0 2\n", "invalid angle"),
        ("qubits 3\nmcx 0 9\n", "out of range"),
        ("qubits 1\nif 1\n", "expects a value then an instruction"),
        ("qubits 1\nif nope x 0\n", "invalid conditional value"),
        (
            "qubits 1\nif 1 measure 0\n",
            "only a gate may be conditional",
        ),
        (
            "qubits 1\nif 1 if 1 x 0\n",
            "only a gate may be conditional",
        ),
        ("qubits 1\nif 1 bogus 0\n", "unknown instruction"),
        (
            "qubits 2\nh 0\nsample 1 1\nsample 2 2\n",
            "may only appear once",
        ),
    ];
    for (src, needle) in cases {
        let err = program::parse(src).unwrap_err().to_string();
        assert!(err.contains(needle), "for `{src}`: {err}");
    }
}

#[test]
fn angle_forms_and_errors() {
    // Negative pi form and a coefficient form both parse.
    let p = program::parse("qubits 1\nrz -pi/2 0\nrz 2pi 0\n").expect("should parse");
    assert!((p.circuit.run().norm() - 1.0).abs() < 1e-12);

    // Angles are full arithmetic expressions, so a sum, a parenthesised group,
    // and a function call all work and agree with the equivalent literal. The
    // native format splits on whitespace, so an angle is one unspaced token.
    let expr = program::parse("qubits 1\nrx (pi/4+pi/4)*2 0\n").expect("should parse");
    let literal = program::parse("qubits 1\nrx pi 0\n").unwrap();
    assert_relative_eq!(
        expr.circuit.run().probability(1),
        literal.circuit.run().probability(1),
        epsilon = 1e-12
    );
    assert!(program::parse("qubits 1\nrz sqrt(4)*pi 0\n").is_ok());

    for (src, needle) in [
        ("qubits 1\nrx bogus 0\n", "invalid angle"),
        ("qubits 1\nrx pi/0 0\n", "division by zero"),
    ] {
        let err = program::parse(src).unwrap_err().to_string();
        assert!(err.contains(needle), "for `{src}`: {err}");
    }
}
