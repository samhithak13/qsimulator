//! A tiny line-based program format for describing circuits as text, so the
//! CLI can run a circuit without recompiling.
//!
//! One instruction per line. `#` starts a comment; blank lines are ignored.
//! The first instruction must be `qubits N`. Example:
//!
//! ```text
//! qubits 3
//! h 0
//! cnot 0 1
//! cnot 1 2      # GHZ state
//! sample 1000 42
//! ```
//!
//! Supported instructions: `qubits N`; single-qubit `id/h/x/y/z/s/t/sdg/tdg/sx/sxdg Q`;
//! rotations `rx/ry/rz THETA Q`, phase `p THETA Q`, and the general
//! `u2 PHI LAMBDA Q` / `u3 THETA PHI LAMBDA Q`; two-qubit `cnot/cy/cz/ch C T`,
//! `crz/cp THETA C T`, `cu3 THETA PHI LAMBDA C T`, and `swap A B`; three-qubit
//! `toffoli C1 C2 T` and `cswap C A B`; the open-ended multi-controlled
//! `mcx C... T` and `mcu3 THETA PHI LAMBDA C... T`; `measure Q`, which
//! collapses a qubit mid-circuit; and a terminal
//! `sample SHOTS SEED`. An angle is an arithmetic expression over numbers and
//! `pi` — `0.7`, `pi`, `pi/2`, `-pi/4`, `2pi`, or `(pi/4 + 0.1)*2` — with the
//! operators `+ - * / ^`, parentheses, and `sin`/`cos`/`tan`/`exp`/`ln`/`sqrt`.

use crate::error::ParseError;
use crate::Circuit;

/// Largest register a program may declare, matching the OpenQASM importer.
/// Guards against a malformed huge `qubits N` aborting on allocation.
const MAX_QUBITS: usize = 30;

/// A parsed program: the built circuit plus an optional `sample` directive.
#[derive(Debug, Clone)]
pub struct Program {
    /// The circuit built from the program's gate instructions.
    pub circuit: Circuit,
    /// The `sample SHOTS SEED` directive, if the program had one.
    pub sample: Option<SampleSpec>,
}

/// Parameters from a `sample SHOTS SEED` directive.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SampleSpec {
    /// Number of shots to sample.
    pub shots: usize,
    /// RNG seed, for reproducible sampling.
    pub seed: u64,
}

/// Parse a program from its textual form, or return a [`ParseError`] tagged
/// with the offending line number.
pub fn parse(src: &str) -> Result<Program, ParseError> {
    let mut circuit: Option<Circuit> = None;
    let mut n = 0usize;
    let mut sample: Option<SampleSpec> = None;

    for (lineno, raw) in src.lines().enumerate() {
        // Strip comments and surrounding whitespace.
        let line = raw.split('#').next().unwrap_or("").trim();
        if line.is_empty() {
            continue;
        }
        let toks: Vec<&str> = line.split_whitespace().collect();
        let cmd = toks[0];
        let at = |msg: String| ParseError::at_line(lineno + 1, msg);

        if cmd == "qubits" {
            if circuit.is_some() {
                return Err(at("`qubits` may only appear once".into()));
            }
            expect_arity(&toks, 2).map_err(&at)?;
            n = toks[1]
                .parse()
                .map_err(|_| at(format!("invalid qubit count `{}`", toks[1])))?;
            if n == 0 {
                return Err(at("qubit count must be >= 1".into()));
            }
            if n > MAX_QUBITS {
                return Err(at(format!(
                    "qubit count {n} exceeds the maximum of {MAX_QUBITS}"
                )));
            }
            circuit = Some(Circuit::new(n));
            continue;
        }

        let c = circuit
            .as_mut()
            .ok_or_else(|| at("first instruction must be `qubits N`".into()))?;

        match cmd {
            "id" | "h" | "x" | "y" | "z" | "s" | "t" | "sdg" | "tdg" | "sx" | "sxdg" => {
                expect_arity(&toks, 2).map_err(&at)?;
                let q = parse_qubit(&toks, 1, n).map_err(&at)?;
                match cmd {
                    "id" => c.id(q),
                    "h" => c.h(q),
                    "x" => c.x(q),
                    "y" => c.y(q),
                    "z" => c.z(q),
                    "s" => c.s(q),
                    "t" => c.t(q),
                    "sdg" => c.sdg(q),
                    "sx" => c.sx(q),
                    "sxdg" => c.sxdg(q),
                    "tdg" => c.tdg(q),
                    _ => unreachable!(),
                };
            }
            "rx" | "ry" | "rz" | "p" => {
                expect_arity(&toks, 3).map_err(&at)?;
                let theta = parse_angle(toks[1]).map_err(&at)?;
                let q = parse_qubit(&toks, 2, n).map_err(&at)?;
                match cmd {
                    "rx" => c.rx(theta, q),
                    "ry" => c.ry(theta, q),
                    "rz" => c.rz(theta, q),
                    "p" => c.p(theta, q),
                    _ => unreachable!(),
                };
            }
            "u2" => {
                expect_arity(&toks, 4).map_err(&at)?;
                let phi = parse_angle(toks[1]).map_err(&at)?;
                let lambda = parse_angle(toks[2]).map_err(&at)?;
                let q = parse_qubit(&toks, 3, n).map_err(&at)?;
                c.u2(phi, lambda, q);
            }
            "u3" => {
                expect_arity(&toks, 5).map_err(&at)?;
                let theta = parse_angle(toks[1]).map_err(&at)?;
                let phi = parse_angle(toks[2]).map_err(&at)?;
                let lambda = parse_angle(toks[3]).map_err(&at)?;
                let q = parse_qubit(&toks, 4, n).map_err(&at)?;
                c.u3(theta, phi, lambda, q);
            }
            "cnot" | "cz" | "cy" | "ch" => {
                expect_arity(&toks, 3).map_err(&at)?;
                let ctrl = parse_qubit(&toks, 1, n).map_err(&at)?;
                let tgt = parse_qubit(&toks, 2, n).map_err(&at)?;
                if ctrl == tgt {
                    return Err(at("control and target must differ".into()));
                }
                match cmd {
                    "cnot" => c.cnot(ctrl, tgt),
                    "cy" => c.cy(ctrl, tgt),
                    "ch" => c.ch(ctrl, tgt),
                    _ => c.cz(ctrl, tgt),
                };
            }
            "cswap" => {
                expect_arity(&toks, 4).map_err(&at)?;
                let ctrl = parse_qubit(&toks, 1, n).map_err(&at)?;
                let a = parse_qubit(&toks, 2, n).map_err(&at)?;
                let b = parse_qubit(&toks, 3, n).map_err(&at)?;
                if ctrl == a || ctrl == b || a == b {
                    return Err(at("cswap qubits must be distinct".into()));
                }
                c.cswap(ctrl, a, b);
            }
            "crz" | "cp" => {
                expect_arity(&toks, 4).map_err(&at)?;
                let theta = parse_angle(toks[1]).map_err(&at)?;
                let ctrl = parse_qubit(&toks, 2, n).map_err(&at)?;
                let tgt = parse_qubit(&toks, 3, n).map_err(&at)?;
                if ctrl == tgt {
                    return Err(at("control and target must differ".into()));
                }
                match cmd {
                    "crz" => c.crz(theta, ctrl, tgt),
                    _ => c.cp(theta, ctrl, tgt),
                };
            }
            "cu3" => {
                expect_arity(&toks, 6).map_err(&at)?;
                let theta = parse_angle(toks[1]).map_err(&at)?;
                let phi = parse_angle(toks[2]).map_err(&at)?;
                let lambda = parse_angle(toks[3]).map_err(&at)?;
                let ctrl = parse_qubit(&toks, 4, n).map_err(&at)?;
                let tgt = parse_qubit(&toks, 5, n).map_err(&at)?;
                if ctrl == tgt {
                    return Err(at("control and target must differ".into()));
                }
                c.cu3(theta, phi, lambda, ctrl, tgt);
            }
            "swap" => {
                expect_arity(&toks, 3).map_err(&at)?;
                let a = parse_qubit(&toks, 1, n).map_err(&at)?;
                let b = parse_qubit(&toks, 2, n).map_err(&at)?;
                c.swap(a, b);
            }
            "toffoli" => {
                expect_arity(&toks, 4).map_err(&at)?;
                let c1 = parse_qubit(&toks, 1, n).map_err(&at)?;
                let c2 = parse_qubit(&toks, 2, n).map_err(&at)?;
                let tgt = parse_qubit(&toks, 3, n).map_err(&at)?;
                if c1 == tgt || c2 == tgt {
                    return Err(at("control and target must differ".into()));
                }
                c.toffoli(c1, c2, tgt);
            }
            // Multi-controlled gates: every token but the last is a control,
            // the last is the target, so the arity is open-ended.
            "mcx" => {
                let (controls, tgt) = parse_control_list(&toks, 1, n).map_err(&at)?;
                c.mcx(&controls, tgt);
            }
            "mcu3" => {
                let theta = parse_angle(toks.get(1).copied().unwrap_or("")).map_err(&at)?;
                let phi = parse_angle(toks.get(2).copied().unwrap_or("")).map_err(&at)?;
                let lambda = parse_angle(toks.get(3).copied().unwrap_or("")).map_err(&at)?;
                let (controls, tgt) = parse_control_list(&toks, 4, n).map_err(&at)?;
                c.mcu(crate::gates::u3(theta, phi, lambda), &controls, tgt);
            }
            "measure" => {
                expect_arity(&toks, 2).map_err(&at)?;
                let q = parse_qubit(&toks, 1, n).map_err(&at)?;
                c.measure(q);
            }
            "sample" => {
                if sample.is_some() {
                    return Err(at("`sample` may only appear once".into()));
                }
                expect_arity(&toks, 3).map_err(&at)?;
                let shots = toks[1]
                    .parse()
                    .map_err(|_| at(format!("invalid shot count `{}`", toks[1])))?;
                let seed = toks[2]
                    .parse()
                    .map_err(|_| at(format!("invalid seed `{}`", toks[2])))?;
                sample = Some(SampleSpec { shots, seed });
            }
            other => return Err(at(format!("unknown instruction `{other}`"))),
        }
    }

    let circuit =
        circuit.ok_or_else(|| ParseError::new("program is empty (expected `qubits N`)"))?;
    Ok(Program { circuit, sample })
}

/// Require a line to have exactly `k` whitespace tokens.
fn expect_arity(toks: &[&str], k: usize) -> Result<(), String> {
    if toks.len() == k {
        Ok(())
    } else {
        Err(format!(
            "`{}` expects {} tokens, got {}",
            toks[0],
            k,
            toks.len()
        ))
    }
}

/// Parse the qubit list of a multi-controlled instruction: tokens from `start`
/// to the end, where the last is the target and the rest are controls.
///
/// Requires at least one control (use the plain gate for none) and rejects a
/// repeated qubit, which would mean controlling a gate on its own target.
fn parse_control_list(
    toks: &[&str],
    start: usize,
    n: usize,
) -> Result<(Vec<usize>, usize), String> {
    if toks.len() < start + 2 {
        return Err(format!(
            "`{}` expects at least {} tokens (one or more controls, then a target), got {}",
            toks[0],
            start + 2,
            toks.len()
        ));
    }
    let qubits: Vec<usize> = (start..toks.len())
        .map(|i| parse_qubit(toks, i, n))
        .collect::<Result<_, _>>()?;
    for (i, q) in qubits.iter().enumerate() {
        if qubits[..i].contains(q) {
            return Err(format!("qubit {q} is repeated; all must be distinct"));
        }
    }
    let (&target, controls) = qubits.split_last().expect("checked non-empty above");
    Ok((controls.to_vec(), target))
}

/// Parse token `idx` as a qubit index and bounds-check it against `n`.
fn parse_qubit(toks: &[&str], idx: usize, n: usize) -> Result<usize, String> {
    let s = toks
        .get(idx)
        .ok_or_else(|| "missing qubit argument".to_string())?;
    let q: usize = s.parse().map_err(|_| format!("invalid qubit `{s}`"))?;
    if q >= n {
        return Err(format!("qubit {q} out of range (0..{n})"));
    }
    Ok(q)
}

/// Parse an angle: an arithmetic expression over numbers and `pi`, such as
/// `0.7`, `pi`, `pi/2`, `-pi/4`, `2pi`, `2*pi`, or `pi/2 + 0.3`.
///
/// Shared with the OpenQASM importer, which uses the same angle syntax — see
/// [`crate::expr`] for the full grammar.
pub(crate) fn parse_angle(s: &str) -> Result<f64, String> {
    crate::expr::eval(s, &std::collections::HashMap::new())
}
