//! End-to-end tests of the `qsimulator` CLI binary. Cargo builds the binary
//! and exposes its path via `CARGO_BIN_EXE_qsimulator`.

use std::io::Write;
use std::process::{Command, Stdio};

/// Run the CLI with `args`, optionally piping `stdin`. Returns
/// `(stdout, stderr, exit_code)`.
fn run(args: &[&str], stdin: Option<&str>) -> (String, String, i32) {
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_qsimulator"));
    cmd.args(args).stdout(Stdio::piped()).stderr(Stdio::piped());
    if stdin.is_some() {
        cmd.stdin(Stdio::piped());
    }
    let mut child = cmd.spawn().expect("failed to spawn qsimulator");
    if let Some(input) = stdin {
        child
            .stdin
            .take()
            .unwrap()
            .write_all(input.as_bytes())
            .unwrap();
    }
    let out = child.wait_with_output().unwrap();
    (
        String::from_utf8_lossy(&out.stdout).into_owned(),
        String::from_utf8_lossy(&out.stderr).into_owned(),
        out.status.code().unwrap_or(-1),
    )
}

#[test]
fn no_args_runs_the_demo() {
    let (stdout, _err, code) = run(&[], None);
    assert_eq!(code, 0);
    assert!(stdout.contains("Bell-state demo"), "{stdout}");
    assert!(stdout.contains("|00>"), "{stdout}");
}

#[test]
fn help_lists_usage() {
    for flag in ["--help", "-h"] {
        let (stdout, _err, code) = run(&[flag], None);
        assert_eq!(code, 0);
        assert!(stdout.contains("USAGE"), "{stdout}");
        assert!(stdout.contains("PROGRAM FORMAT"), "{stdout}");
    }
}

#[test]
fn runs_a_program_from_stdin() {
    let (stdout, _err, code) = run(&["-"], Some("qubits 1\nx 0\n"));
    assert_eq!(code, 0);
    assert!(stdout.contains("|1>: p = 1.000"), "{stdout}");
}

#[test]
fn runs_a_program_with_sampling() {
    let (stdout, _err, code) = run(&["-"], Some("qubits 2\nh 0\ncnot 0 1\nsample 100 7\n"));
    assert_eq!(code, 0);
    assert!(stdout.contains("Sampling 100 shots"), "{stdout}");
}

#[test]
fn runs_a_qasm_file() {
    let (stdout, _err, code) = run(&["programs/bell.qasm"], None);
    assert_eq!(code, 0);
    assert!(stdout.contains("|00>: p = 0.500"), "{stdout}");
    assert!(stdout.contains("|11>: p = 0.500"), "{stdout}");
}

#[test]
fn emit_qasm_prints_openqasm() {
    let (stdout, _err, code) = run(&["--emit-qasm", "-"], Some("qubits 2\nh 0\ncnot 0 1\n"));
    assert_eq!(code, 0);
    assert!(stdout.starts_with("OPENQASM 2.0;"), "{stdout}");
    assert!(stdout.contains("cx q[0],q[1];"), "{stdout}");
}

#[test]
fn statevector_prints_json() {
    let (stdout, _err, code) = run(&["--statevector", "-"], Some("qubits 1\nx 0\n"));
    assert_eq!(code, 0);
    // X|0> = |1>: amplitude vector [[0,0],[1,0]].
    assert_eq!(stdout.trim(), "[[0,0],[1,0]]");
}

#[test]
fn parse_error_exits_nonzero() {
    let (_out, err, code) = run(&["-"], Some("qubits 2\nbogus 0\n"));
    assert_eq!(code, 1);
    assert!(err.contains("error:"), "{err}");
    assert!(err.contains("unknown instruction"), "{err}");
}

#[test]
fn missing_file_exits_nonzero() {
    let (_out, err, code) = run(&["/no/such/file.qsim"], None);
    assert_eq!(code, 1);
    assert!(err.contains("cannot read"), "{err}");
}

#[test]
fn too_many_args_exits_with_usage_code() {
    let (_out, err, code) = run(&["a", "b", "c"], None);
    assert_eq!(code, 2);
    // The message names the first argument that had nowhere to go.
    assert!(err.contains("unexpected extra argument `b`"), "{err}");
}

/// `--shots` samples a program that never asked to be sampled — including an
/// OpenQASM file, which has no way to carry a `sample` directive at all.
#[test]
fn shots_flag_samples_a_program_without_a_directive() {
    let (stdout, _err, code) = run(&["--shots", "16", "-"], Some("qubits 2\nh 0\ncnot 0 1\n"));
    assert_eq!(code, 0);
    assert!(stdout.contains("Sampling 16 shots"), "{stdout}");

    let (qasm_out, _err, code) = run(&["--shots", "16", "programs/bell.qasm"], None);
    assert_eq!(code, 0);
    assert!(qasm_out.contains("Sampling 16 shots"), "{qasm_out}");
}

/// Sampling is a pure function of the seed: same seed, same histogram;
/// different seed, different histogram.
#[test]
fn seed_flag_makes_sampling_reproducible() {
    let args = ["--shots", "200", "--seed", "7", "programs/bell.qasm"];
    let (a, _, _) = run(&args, None);
    let (b, _, _) = run(&args, None);
    assert_eq!(a, b);

    let (c, _, _) = run(
        &["--shots", "200", "--seed", "8", "programs/bell.qasm"],
        None,
    );
    assert_ne!(a, c);
}

/// `--shots` overrides the program's own directive; `--seed` alone re-seeds it
/// while keeping the shot count.
#[test]
fn flags_override_the_sample_directive() {
    let program = "qubits 2\nh 0\ncnot 0 1\nsample 50 1\n";
    let (stdout, _err, code) = run(&["--shots", "9", "-"], Some(program));
    assert_eq!(code, 0);
    assert!(stdout.contains("Sampling 9 shots (seed = 1)"), "{stdout}");

    let (stdout, _err, code) = run(&["--seed", "42", "-"], Some(program));
    assert_eq!(code, 0);
    assert!(stdout.contains("Sampling 50 shots (seed = 42)"), "{stdout}");
}

/// Flag misuse is a usage error (exit 2), not a silently ignored flag.
#[test]
fn sampling_flag_misuse_exits_with_usage_code() {
    let cases: [&[&str]; 5] = [
        &["--seed", "3", "programs/bell.qasm"],
        &["--shots", "4", "--statevector", "programs/bell.qasm"],
        &["--shots", "abc", "programs/bell.qasm"],
        &["--shots"],
        &["--bogus", "programs/bell.qasm"],
    ];
    for args in cases {
        let (_out, err, code) = run(args, None);
        assert_eq!(code, 2, "for {args:?}: {err}");
        assert!(err.starts_with("error: "), "for {args:?}: {err}");
    }
}
