//! qsimulator CLI.
//!
//! With no arguments, runs a built-in Bell-state demo. Given a program file
//! (or `-` for stdin), parses and runs it — see `qsimulator --help` or the
//! `qsimulator::program` module for the format.

use qsimulator::program::{self, Program};
use qsimulator::Circuit;
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
        [path] => {
            let src = match read_source(path) {
                Ok(s) => s,
                Err(e) => {
                    eprintln!("error: cannot read {path}: {e}");
                    return ExitCode::FAILURE;
                }
            };
            match program::parse(&src) {
                Ok(prog) => {
                    run_program(&prog);
                    ExitCode::SUCCESS
                }
                Err(e) => {
                    eprintln!("error: {e}");
                    ExitCode::FAILURE
                }
            }
        }
        _ => {
            eprintln!("error: too many arguments; expected a program file, `-`, or none");
            eprintln!("try `qsimulator --help`");
            ExitCode::from(2)
        }
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

/// Print the circuit diagram, the final-state probabilities, and (if the
/// program had a `sample` directive) the sampled histogram.
fn run_program(prog: &Program) {
    println!("Circuit:");
    println!("{}\n", prog.circuit);

    let state = prog.circuit.run();
    println!("Final-state probabilities:");
    print_probabilities(&state);

    if let Some(spec) = prog.sample {
        let histogram = prog.circuit.sample(spec.shots, spec.seed);
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
    qsimulator <FILE>          Parse and run a program file
    qsimulator -               Read a program from stdin
    qsimulator --help          Show this help

PROGRAM FORMAT (one instruction per line; `#` starts a comment):
    qubits N                   Declare the register size (must come first)
    h|x|y|z|s|t Q              Single-qubit gate on qubit Q
    rx|ry|rz THETA Q           Rotation by THETA (float, or pi, pi/2, -pi/4, 2pi)
    cnot|cz C T                Two-qubit controlled gate (control C, target T)
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
