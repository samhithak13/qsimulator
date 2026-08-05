//! qsimulator CLI.
//!
//! With no arguments, runs a built-in Bell-state demo. Given a program file
//! (or `-` for stdin), parses and runs it — see `qsimulator --help` or the
//! `qsimulator::program` module for the format.

use qsimulator::program::{self, SampleSpec};
use qsimulator::{qasm, Circuit};
use std::io::Read;
use std::process::ExitCode;

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    match args.as_slice() {
        [] => {
            run_bell_demo();
            ExitCode::SUCCESS
        }
        [flag] if flag == "-h" || flag == "--help" => {
            print_help();
            ExitCode::SUCCESS
        }
        // `--emit-qasm FILE`: parse the input and print it back as OpenQASM.
        [flag, path] if flag == "--emit-qasm" => match load_circuit(path) {
            Ok((circuit, _sample)) => match circuit.to_qasm() {
                Ok(qasm) => {
                    print!("{qasm}");
                    ExitCode::SUCCESS
                }
                Err(e) => {
                    eprintln!("error: {e}");
                    ExitCode::FAILURE
                }
            },
            Err(e) => {
                eprintln!("error: {e}");
                ExitCode::FAILURE
            }
        },
        // `--statevector FILE`: print the final amplitudes as JSON, for
        // machine consumption (e.g. the Qiskit cross-validation harness).
        [flag, path] if flag == "--statevector" => match load_circuit(path) {
            Ok((circuit, _sample)) => {
                print_statevector(&circuit.run());
                ExitCode::SUCCESS
            }
            Err(e) => {
                eprintln!("error: {e}");
                ExitCode::FAILURE
            }
        },
        [path] => match load_circuit(path) {
            Ok((circuit, sample)) => {
                run_circuit(&circuit, sample);
                ExitCode::SUCCESS
            }
            Err(e) => {
                eprintln!("error: {e}");
                ExitCode::FAILURE
            }
        },
        _ => {
            eprintln!("error: unexpected arguments; expected a program file, `-`, or none");
            eprintln!("try `qsimulator --help`");
            ExitCode::from(2)
        }
    }
}

/// Read a program from `path` (or stdin when `path` is `-`) and parse it into a
/// circuit, choosing the OpenQASM importer by `.qasm` extension or `OPENQASM`
/// header, otherwise the native text program format.
fn load_circuit(path: &str) -> Result<(Circuit, Option<SampleSpec>), String> {
    let src = read_source(path).map_err(|e| format!("cannot read {path}: {e}"))?;
    let is_qasm = path.ends_with(".qasm") || src.trim_start().starts_with("OPENQASM");
    if is_qasm {
        qasm::parse(&src).map(|circuit| (circuit, None))
    } else {
        program::parse(&src).map(|prog| (prog.circuit, prog.sample))
    }
}

/// Read a program from a file path, or from stdin when `path` is `-`.
fn read_source(path: &str) -> std::io::Result<String> {
    if path == "-" {
        let mut buf = String::new();
        std::io::stdin().read_to_string(&mut buf)?;
        Ok(buf)
    } else {
        std::fs::read_to_string(path)
    }
}

/// Print the final state as a JSON array of `[re, im]` amplitude pairs, in
/// little-endian basis-state order. Full `f64` precision, one flat line.
fn print_statevector(state: &qsimulator::State) {
    let mut out = String::from("[");
    for (i, a) in state.amplitudes().iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        out.push_str(&format!("[{},{}]", a.re, a.im));
    }
    out.push(']');
    println!("{out}");
}

/// Print the circuit diagram, the final-state probabilities, and (if a
/// `sample` directive was present) the sampled histogram.
fn run_circuit(circuit: &Circuit, sample: Option<SampleSpec>) {
    println!("Circuit:");
    println!("{circuit}\n");

    let state = circuit.run();
    println!("Final-state probabilities:");
    print_probabilities(&state);

    if let Some(spec) = sample {
        let histogram = circuit.sample(spec.shots, spec.seed);
        println!("\nSampling {} shots (seed = {}):", spec.shots, spec.seed);
        let width = state.n_qubits();
        for outcome in 0..(1usize << width) {
            let count = histogram.get(&outcome).copied().unwrap_or(0);
            if count > 0 {
                println!("  |{outcome:0width$b}>: {count} shots", width = width);
            }
        }
    }
}

/// Print every basis state with nonzero probability.
fn print_probabilities(state: &qsimulator::State) {
    let width = state.n_qubits();
    for i in 0..state.amplitudes().len() {
        let p = state.probability(i);
        if p > 1e-12 {
            println!("  |{i:0width$b}>: p = {p:.3}", width = width);
        }
    }
    println!("  total probability = {:.6}", state.norm());
}

/// The default demo: a Bell state, printed and sampled.
fn run_bell_demo() {
    let mut circuit = Circuit::new(2);
    circuit.h(0).cnot(0, 1);

    println!("Bell-state demo (run `qsimulator --help` for the program format)\n");
    println!("Circuit:");
    println!("{circuit}\n");

    let state = circuit.run();
    println!("Final-state probabilities:");
    print_probabilities(&state);

    let shots = 1000;
    let seed = 0xC0FF_EE00;
    let histogram = circuit.sample(shots, seed);
    println!("\nSampling {shots} shots (seed = {seed:#x}):");
    for outcome in 0..(1usize << state.n_qubits()) {
        let count = histogram.get(&outcome).copied().unwrap_or(0);
        println!("  |{outcome:02b}>: {count} shots");
    }
}

fn print_help() {
    println!(
        "qsimulator — a state-vector quantum circuit simulator

USAGE:
    qsimulator                 Run the built-in Bell-state demo
    qsimulator <FILE>          Parse and run a program file (.qasm = OpenQASM)
    qsimulator -               Read a program from stdin
    qsimulator --emit-qasm <FILE>   Print the circuit as OpenQASM 2.0
    qsimulator --statevector <FILE> Print final amplitudes as JSON
    qsimulator --help          Show this help

An input is treated as OpenQASM 2.0 if it ends in `.qasm` or begins with an
`OPENQASM` header (see programs/bell.qasm); otherwise the native format below.

PROGRAM FORMAT (one instruction per line; `#` starts a comment):
    qubits N                   Declare the register size (must come first)
    h|x|y|z|s|t|sdg|tdg Q      Single-qubit gate on qubit Q
    rx|ry|rz|p THETA Q         Rotation/phase by THETA (float, or pi, pi/2, ...)
    u2 PHI LAMBDA Q            General single-qubit gate U2
    u3 THETA PHI LAMBDA Q      General single-qubit gate U3
    cnot|cz C T                Two-qubit controlled gate (control C, target T)
    crz|cp THETA C T           Controlled-Rz / controlled-phase by THETA
    swap A B                   Exchange qubits A and B
    toffoli C1 C2 T            CCNOT
    sample SHOTS SEED          Sample the final state (optional, once)

EXAMPLE (GHZ state):
    qubits 3
    h 0
    cnot 0 1
    cnot 1 2
    sample 1000 42"
    );
}
