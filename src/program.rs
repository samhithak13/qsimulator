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
//! Supported instructions: `qubits N`; single-qubit `h/x/y/z/s/t Q`;
//! rotations `rx/ry/rz THETA Q`; two-qubit `cnot/cz C T` and `swap A B`;
//! `toffoli C1 C2 T`; and a terminal `sample SHOTS SEED`. Angles are a plain
//! float or a symbolic multiple of pi such as `pi`, `pi/2`, `-pi/4`, or `2pi`.

use crate::Circuit;
use std::f64::consts::PI;

/// A parsed program: the built circuit plus an optional `sample` directive.
#[derive(Debug, Clone)]
pub struct Program {
    pub circuit: Circuit,
    pub sample: Option<SampleSpec>,
}

/// Parameters from a `sample SHOTS SEED` directive.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SampleSpec {
    pub shots: usize,
    pub seed: u64,
}

/// Parse a program from its textual form, or return a human-readable error
/// message tagged with the offending line number.
pub fn parse(src: &str) -> Result<Program, String> {
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
        let at = |msg: String| format!("line {}: {}", lineno + 1, msg);

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
            circuit = Some(Circuit::new(n));
            continue;
        }

        let c = circuit
            .as_mut()
            .ok_or_else(|| at("first instruction must be `qubits N`".into()))?;

        match cmd {
            "h" | "x" | "y" | "z" | "s" | "t" => {
                expect_arity(&toks, 2).map_err(&at)?;
                let q = parse_qubit(&toks, 1, n).map_err(&at)?;
                match cmd {
                    "h" => c.h(q),
                    "x" => c.x(q),
                    "y" => c.y(q),
                    "z" => c.z(q),
                    "s" => c.s(q),
                    "t" => c.t(q),
                    _ => unreachable!(),
                };
            }
            "rx" | "ry" | "rz" => {
                expect_arity(&toks, 3).map_err(&at)?;
                let theta = parse_angle(toks[1]).map_err(&at)?;
                let q = parse_qubit(&toks, 2, n).map_err(&at)?;
                match cmd {
                    "rx" => c.rx(theta, q),
                    "ry" => c.ry(theta, q),
                    "rz" => c.rz(theta, q),
                    _ => unreachable!(),
                };
            }
            "cnot" | "cz" => {
                expect_arity(&toks, 3).map_err(&at)?;
                let ctrl = parse_qubit(&toks, 1, n).map_err(&at)?;
                let tgt = parse_qubit(&toks, 2, n).map_err(&at)?;
                if ctrl == tgt {
                    return Err(at("control and target must differ".into()));
                }
                match cmd {
                    "cnot" => c.cnot(ctrl, tgt),
                    _ => c.cz(ctrl, tgt),
                };
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

    let circuit = circuit.ok_or_else(|| "program is empty (expected `qubits N`)".to_string())?;
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

/// Parse an angle: a plain float, or a symbolic multiple of pi such as `pi`,
/// `pi/2`, `-pi/4`, `2pi`, `2*pi`, or `0.5*pi`.
///
/// Shared with the OpenQASM importer, which uses the same angle syntax.
pub(crate) fn parse_angle(s: &str) -> Result<f64, String> {
    if let Ok(v) = s.parse::<f64>() {
        return Ok(v);
    }
    let (sign, rest) = match s.strip_prefix('-') {
        Some(r) => (-1.0, r),
        None => (1.0, s),
    };
    let body = rest.replace('*', "");
    let (coeff_str, denom_str) = body
        .find("pi")
        .map(|i| (&body[..i], &body[i + 2..]))
        .ok_or_else(|| format!("invalid angle `{s}`"))?;

    let coeff: f64 = if coeff_str.is_empty() {
        1.0
    } else {
        coeff_str
            .parse()
            .map_err(|_| format!("invalid angle `{s}`"))?
    };
    let denom: f64 = if denom_str.is_empty() {
        1.0
    } else if let Some(d) = denom_str.strip_prefix('/') {
        d.parse().map_err(|_| format!("invalid angle `{s}`"))?
    } else {
        return Err(format!("invalid angle `{s}`"));
    };
    if denom == 0.0 {
        return Err(format!("invalid angle `{s}` (division by zero)"));
    }
    Ok(sign * coeff * PI / denom)
}
