//! qsimulator CLI.
//!
//! With no arguments, runs a built-in Bell-state demo. Given a program file
//! (or `-` for stdin), parses and runs it — see `qsimulator --help` or the
//! `qsimulator::program` module for the format.

use qsimulator::program::{self, SampleSpec};
use qsimulator::{qasm, Circuit};
use std::io::Read;
use std::process::ExitCode;

/// What the CLI was asked to produce.
enum Mode {
    /// No arguments: the built-in demo.
    Demo,
    Help,
    /// Run the program and print the diagram, probabilities, and any sampling.
    Run,
    /// Print the circuit as OpenQASM 2.0.
    EmitQasm,
    /// Print the final amplitudes as JSON.
    Statevector,
}

/// The parsed command line.
struct Cli {
    mode: Mode,
    path: Option<String>,
    shots: Option<usize>,
    seed: Option<u64>,
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let cli = match parse_args(&args) {
        Ok(cli) => cli,
        Err(e) => {
            eprintln!("error: {e}");
            eprintln!("try `qsimulator --help`");
            return ExitCode::from(2);
        }
    };

    match cli.mode {
        Mode::Help => {
            print_help();
            return ExitCode::SUCCESS;
        }
        Mode::Demo => {
            run_bell_demo();
            return ExitCode::SUCCESS;
        }
        _ => {}
    }

    let path = cli.path.as_deref().expect("a path outside Demo/Help mode");
    let (circuit, spec) = match load_circuit(path) {
        Ok(loaded) => loaded,
        Err(e) => {
            eprintln!("error: {e}");
            return ExitCode::FAILURE;
        }
    };

    match cli.mode {
        Mode::EmitQasm => match circuit.to_qasm() {
            Ok(qasm) => print!("{qasm}"),
            Err(e) => {
                eprintln!("error: {e}");
                return ExitCode::FAILURE;
            }
        },
        Mode::Statevector => print_statevector(&circuit.run()),
        _ => {
            let sample = match effective_sample(spec, cli.shots, cli.seed) {
                Ok(sample) => sample,
                Err(e) => {
                    eprintln!("error: {e}");
                    return ExitCode::from(2);
                }
            };
            run_circuit(&circuit, sample);
        }
    }
    ExitCode::SUCCESS
}

/// Parse the command line by hand — the surface is small enough that an
/// argument-parsing dependency would cost more than it saves.
fn parse_args(args: &[String]) -> Result<Cli, String> {
    let mut cli = Cli {
        mode: Mode::Run,
        path: None,
        shots: None,
        seed: None,
    };
    let mut saw_mode_flag = false;
    let mut rest = args.iter();

    while let Some(arg) = rest.next() {
        // A flag that takes a value, consuming the next argument.
        let mut value = |flag: &str| -> Result<String, String> {
            rest.next()
                .cloned()
                .ok_or_else(|| format!("`{flag}` needs a value"))
        };
        match arg.as_str() {
            "-h" | "--help" => {
                return Ok(Cli {
                    mode: Mode::Help,
                    ..cli
                })
            }
            "--emit-qasm" | "--statevector" => {
                if saw_mode_flag {
                    return Err("`--emit-qasm` and `--statevector` are exclusive".into());
                }
                saw_mode_flag = true;
                cli.mode = if arg == "--emit-qasm" {
                    Mode::EmitQasm
                } else {
                    Mode::Statevector
                };
            }
            "--shots" => {
                let v = value("--shots")?;
                cli.shots = Some(v.parse().map_err(|_| format!("invalid shot count `{v}`"))?);
            }
            "--seed" => {
                let v = value("--seed")?;
                cli.seed = Some(v.parse().map_err(|_| format!("invalid seed `{v}`"))?);
            }
            // `-` is stdin, so only a longer leading dash is an unknown flag.
            other if other.starts_with('-') && other != "-" => {
                return Err(format!("unknown flag `{other}`"));
            }
            other => {
                if cli.path.is_some() {
                    return Err(format!("unexpected extra argument `{other}`"));
                }
                cli.path = Some(other.to_string());
            }
        }
    }

    if cli.path.is_none() {
        if saw_mode_flag || cli.shots.is_some() || cli.seed.is_some() {
            return Err("expected a program file, or `-` to read stdin".into());
        }
        cli.mode = Mode::Demo;
    }
    if saw_mode_flag && (cli.shots.is_some() || cli.seed.is_some()) {
        return Err("`--shots`/`--seed` do not apply to `--emit-qasm` or `--statevector`".into());
    }
    Ok(cli)
}

/// Combine the program's own `sample` directive with the `--shots`/`--seed`
/// flags, which override it. A seed with nothing to sample is an error rather
/// than a silently ignored flag.
fn effective_sample(
    spec: Option<SampleSpec>,
    shots: Option<usize>,
    seed: Option<u64>,
) -> Result<Option<SampleSpec>, String> {
    match (shots, seed) {
        (None, None) => Ok(spec),
        (Some(shots), seed) => Ok(Some(SampleSpec {
            shots,
            seed: seed.or(spec.map(|s| s.seed)).unwrap_or(0),
        })),
        (None, Some(seed)) => match spec {
            Some(s) => Ok(Some(SampleSpec {
                shots: s.shots,
                seed,
            })),
            None => Err("`--seed` needs `--shots`, or a `sample` directive in the program".into()),
        },
    }
}

/// Read a program from `path` (or stdin when `path` is `-`) and parse it into a
/// circuit, choosing the OpenQASM importer by `.qasm` extension or `OPENQASM`
/// header, otherwise the native text program format.
fn load_circuit(path: &str) -> Result<(Circuit, Option<SampleSpec>), Box<dyn std::error::Error>> {
    let src = read_source(path).map_err(|e| format!("cannot read {path}: {e}"))?;
    let is_qasm = path.ends_with(".qasm") || src.trim_start().starts_with("OPENQASM");
    if is_qasm {
        Ok((qasm::parse(&src)?, None))
    } else {
        let prog = program::parse(&src)?;
        Ok((prog.circuit, prog.sample))
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

OPTIONS:
    --shots N                  Sample the final state N times, whatever the
                               program says -- so an OpenQASM file, or one
                               with no `sample` directive, can be sampled too
    --seed S                   Seed the sampling RNG (default 0). Sampling is
                               a pure function of the seed, so a run repeats
                               exactly. On its own, re-seeds the program's own
                               `sample` directive.

An input is treated as OpenQASM 2.0 if it ends in `.qasm` or begins with an
`OPENQASM` header (see programs/bell.qasm); otherwise the native format below.

PROGRAM FORMAT (one instruction per line; `#` starts a comment):
    qubits N                   Declare the register size (must come first)
    id|h|x|y|z|s|t|sdg|tdg|sx|sxdg Q
                               Single-qubit gate (`id` is a no-op)
    rx|ry|rz|p THETA Q         Rotation/phase by THETA (float, or pi, pi/2, ...)
    u2 PHI LAMBDA Q            General single-qubit gate U2
    u3 THETA PHI LAMBDA Q      General single-qubit gate U3
    cnot|cy|cz|ch C T          Two-qubit controlled gate (control C, target T)
    crz|cp THETA C T           Controlled-Rz / controlled-phase by THETA
    cu3 THETA PHI LAMBDA C T   Controlled-U3
    swap A B                   Exchange qubits A and B
    cswap C A B                Controlled swap (Fredkin)
    toffoli C1 C2 T            CCNOT
    mcx C... T                 Multi-controlled X, any number of controls
    mcu3 THETA PHI LAMBDA C... T   Multi-controlled U3
    measure Q                  Collapse qubit Q (readout if nothing follows)
    sample SHOTS SEED          Sample the final state (optional, once)

EXAMPLE (GHZ state):
    qubits 3
    h 0
    cnot 0 1
    cnot 1 2
    sample 1000 42"
    );
}
