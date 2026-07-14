//! Command-line front end for `qsimulator`.
//!
//! Usage:
//!
//! ```text
//! qsimulator             # run the built-in Bell-state demo
//! qsimulator <file>      # run a program file (see src/program.rs)
//! qsimulator -           # read a program from standard input
//! ```
//!
//! A program is a sequence of gate instructions; an optional `sample`
//! directive draws measurement shots. See `src/program.rs` for the grammar
//! and `examples/` for sample programs.

use qsimulator::program::{self, Program};
use qsimulator::state::State;
use qsimulator::Circuit;
use std::io::Read;
use std::process::ExitCode;

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().collect();
    match args.get(1).map(String::as_str) {
        None => {
            run_demo();
            ExitCode::SUCCESS
        }
        Some("-h" | "--help") => {
            print_usage();
            ExitCode::SUCCESS
        }
        Some(source) => match load(source) {
            Ok(src) => run_source(&src),
            Err(e) => {
                eprintln!("error reading {source}: {e}");
                ExitCode::FAILURE
            }
        },
    }
}

fn print_usage() {
    eprintln!(
        "qsimulator — a state-vector quantum circuit simulator\n\n\
         Usage:\n  \
         qsimulator             run the built-in Bell-state demo\n  \
         qsimulator <file>      run a program file\n  \
         qsimulator -           read a program from standard input\n\n\
         Program directives (one per line, `#` starts a comment):\n  \
         qubits N               declare an N-qubit register (required first)\n  \
         h|x|z <t>              single-qubit gate on qubit t\n  \
         rx|ry|rz <angle> <t>   rotation; angle is a float or a pi-expression\n  \
         cnot <c> <t>           controlled-NOT\n  \
         swap <a> <b>           swap two qubits\n  \
         toffoli <c1> <c2> <t>  doubly-controlled NOT\n  \
         sample <shots> [seed]  measure `shots` times (default seed 0)"
    );
}

/// Read a program from a path, or from stdin when `source` is `-`.
fn load(source: &str) -> std::io::Result<String> {
    if source == "-" {
        let mut buf = String::new();
        std::io::stdin().read_to_string(&mut buf)?;
        Ok(buf)
    } else {
        std::fs::read_to_string(source)
    }
}

/// Parse and run a program's source, printing the outcome.
fn run_source(src: &str) -> ExitCode {
    match program::parse(src) {
        Ok(prog) => {
            run_program(&prog);
            ExitCode::SUCCESS
        }
        Err(e) => {
            eprintln!("parse error: {e}");
            ExitCode::FAILURE
        }
    }
}

/// Execute a parsed program: print the final amplitudes, and — if the
/// program asked for it — a sampling histogram.
fn run_program(prog: &Program) {
    let state = prog.circuit.run();
    print_amplitudes(&state);
    if let Some(shots) = prog.shots {
        let histogram = prog.circuit.sample(shots, prog.seed);
        print_histogram(&state, shots, prog.seed, &histogram);
    }
}

/// Print each basis state's amplitude and probability, then the total norm.
fn print_amplitudes(state: &State) {
    let width = state.n_qubits();
    println!("Final amplitudes:");
    for (i, amp) in state.amplitudes().iter().enumerate() {
        if state.probability(i) < 1e-12 {
            continue;
        }
        println!(
            "  |{i:0width$b}>: amp = {:+.3}{:+.3}i, p = {:.3}",
            amp.re,
            amp.im,
            state.probability(i)
        );
    }
    println!("total probability = {:.6}", state.norm());
}

/// Print a sampling histogram over every basis state.
fn print_histogram(
    state: &State,
    shots: usize,
    seed: u64,
    histogram: &std::collections::HashMap<usize, usize>,
) {
    let width = state.n_qubits();
    println!("\nSampling {shots} shots (seed = {seed}):");
    for outcome in 0..(1usize << width) {
        let count = histogram.get(&outcome).copied().unwrap_or(0);
        if count == 0 {
            continue;
        }
        println!("  |{outcome:0width$b}>: {count} shots");
    }
}

/// The built-in demo, run when no program is given: a Bell state, sampled.
fn run_demo() {
    let mut circuit = Circuit::new(2);
    circuit.h(0).cnot(0, 1);
    let prog = Program {
        circuit,
        shots: Some(1000),
        seed: 0xC0FF_EE00,
    };
    println!("(no program given — running the built-in Bell-state demo)\n");
    run_program(&prog);
}
