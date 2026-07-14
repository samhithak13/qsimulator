//! Integration tests for the text program parser (`qsimulator::program`).

use qsimulator::program::parse;

/// A GHZ program produces the (|000> + |111>)/sqrt(2) state.
#[test]
fn ghz_program_builds_ghz_state() {
    let src = "\
qubits 3
h 0
cnot 0 1
cnot 0 2
sample 500 1
";
    let prog = parse(src).expect("program should parse");
    assert_eq!(prog.shots, Some(500));
    assert_eq!(prog.seed, 1);

    let state = prog.circuit.run();
    // Only |000> and |111> carry probability, each 1/2.
    assert!((state.probability(0) - 0.5).abs() < 1e-12);
    assert!((state.probability(7) - 0.5).abs() < 1e-12);
    for i in 1..7 {
        assert!(state.probability(i) < 1e-12, "state {i} should be empty");
    }
}

/// Sampling a parsed program is reproducible for a fixed seed, and every
/// shot lands on one of the two GHZ outcomes.
#[test]
fn ghz_sampling_is_reproducible_and_valid() {
    let src = "qubits 3\nh 0\ncnot 0 1\ncnot 0 2\nsample 400 12345\n";
    let prog = parse(src).unwrap();
    let shots = prog.shots.unwrap();

    let a = prog.circuit.sample(shots, prog.seed);
    let b = prog.circuit.sample(shots, prog.seed);
    assert_eq!(a, b, "same seed must give the same histogram");

    let total: usize = a.values().sum();
    assert_eq!(total, shots);
    // All mass is on |000> and |111>.
    assert_eq!(
        a.get(&0).copied().unwrap_or(0) + a.get(&7).copied().unwrap_or(0),
        shots
    );
}

/// The example file shipped in `examples/` parses and runs.
#[test]
fn shipped_example_parses() {
    let src = include_str!("../examples/ghz.qsim");
    let prog = parse(src).expect("examples/ghz.qsim should parse");
    let state = prog.circuit.run();
    assert!((state.probability(0) - 0.5).abs() < 1e-12);
    assert!((state.probability(7) - 0.5).abs() < 1e-12);
}

/// Parse errors carry a 1-based line number.
#[test]
fn parse_errors_point_at_the_line() {
    let err = parse("qubits 2\nh 0\nbogus 0\n").err().unwrap();
    assert!(err.starts_with("line 3"), "got: {err}");
}
