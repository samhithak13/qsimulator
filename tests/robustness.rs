//! Robustness: the parsers must never panic on arbitrary input — only ever
//! return `Ok` or `Err`. This is a stable-toolchain complement to the
//! `cargo fuzz` targets under `fuzz/`, and it runs in normal CI.
//!
//! Only `parse` is exercised (never `run`), so no large state vector is ever
//! allocated: parsing is `O(input length)` regardless of the declared qubit
//! count.

use qsimulator::{program, qasm, Rng};

/// Build a pseudo-random string from a small alphabet biased toward the tokens
/// the parsers recognize, so the fuzzing reaches deep into the grammar rather
/// than bouncing off the first byte.
fn random_input(rng: &mut Rng, tokens: usize) -> String {
    const ALPHABET: &[&str] = &[
        "OPENQASM", "2.0", "qreg", "creg", "qubits", "q", "r", "c", "h", "x", "y", "z", "s", "t",
        "sdg", "tdg", "rx", "ry", "rz", "p", "u1", "u2", "u3", "cx", "cz", "crz", "cu1", "cp",
        "swap", "ccx", "cnot", "toffoli", "measure", "barrier", "sample", "include", "pi", "->",
        "[", "]", "(", ")", ",", ";", "0", "1", "2", "40", "-", "/", ".", "#", "//", "/*", "*/",
        "\"", " ", "\n",
    ];
    let mut s = String::new();
    for _ in 0..tokens {
        let idx = ((rng.next_f64() * ALPHABET.len() as f64) as usize).min(ALPHABET.len() - 1);
        s.push_str(ALPHABET[idx]);
        s.push(' ');
    }
    s
}

#[test]
fn parsers_never_panic_on_random_input() {
    let mut rng = Rng::new(0xF00D_CAFE);
    for _ in 0..50_000 {
        let tokens = 1 + ((rng.next_f64() * 40.0) as usize);
        let input = random_input(&mut rng, tokens);
        // The only contract under test: these return, they do not panic.
        let _ = program::parse(&input);
        let _ = qasm::parse(&input);
    }
}

#[test]
fn parsers_handle_pathological_bytes() {
    // A few hand-picked shapes that have caused trouble in parsers.
    for input in [
        "",
        ";",
        ";;;;",
        "[]()",
        "q]0[",
        "qreg q]9[;",
        "OPENQASM",
        "qreg q[;",
        "qubits",
        "qubits -1",
        "rx(",
        "u3(,,)",
        "cx q[0],",
        "\n\n\n",
        "// only a comment",
        "/* unclosed",
        "qreg\0q[2];",
    ] {
        let _ = program::parse(input);
        let _ = qasm::parse(input);
    }
}
