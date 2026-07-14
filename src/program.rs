//! A tiny text program format for describing circuits on the command line.
//!
//! A program is a sequence of newline-separated directives. Blank lines and
//! lines beginning with `#` are ignored. The first meaningful line must be a
//! `qubits N` declaration; the remaining lines are gate instructions and an
//! optional `sample` directive. Example:
//!
//! ```text
//! # A 3-qubit GHZ state, sampled 1000 times.
//! qubits 3
//! h 0
//! cnot 0 1
//! cnot 0 2
//! sample 1000 42
//! ```
//!
//! Angles for the rotation gates accept a plain float (`1.5708`) or a
//! `pi`-expression: `pi`, `-pi`, `pi/2`, `2pi`, `2*pi`, `3pi/4`.

use crate::circuit::Circuit;
use std::f64::consts::PI;

/// A parsed program: the built circuit plus an optional sampling request.
pub struct Program {
    /// The circuit assembled from the gate instructions.
    pub circuit: Circuit,
    /// If a `sample` directive was present, the number of shots to draw.
    pub shots: Option<usize>,
    /// The seed for sampling (defaults to `0` when a `sample` directive
    /// omits it).
    pub seed: u64,
}

/// Parse a program from its textual source.
///
/// Returns a human-readable error (prefixed with the 1-based line number)
/// on the first malformed directive.
pub fn parse(src: &str) -> Result<Program, String> {
    let mut circuit: Option<Circuit> = None;
    let mut n_qubits = 0usize;
    let mut shots = None;
    let mut seed = 0u64;

    for (lineno, raw) in src.lines().enumerate() {
        let line = strip_comment(raw).trim();
        if line.is_empty() {
            continue;
        }
        let at = |msg: String| format!("line {}: {msg}", lineno + 1);

        let mut tok = line.split_whitespace();
        let op = tok.next().unwrap(); // non-empty, so at least one token
        let rest: Vec<&str> = tok.collect();

        if op == "qubits" {
            if circuit.is_some() {
                return Err(at("duplicate `qubits` declaration".into()));
            }
            let [n] = expect_args(&rest, op).map_err(at)?;
            n_qubits = parse_index(n).map_err(at)?;
            if n_qubits == 0 {
                return Err(at("`qubits` must be at least 1".into()));
            }
            circuit = Some(Circuit::new(n_qubits));
            continue;
        }

        let c = circuit
            .as_mut()
            .ok_or_else(|| at("first directive must be `qubits N`".into()))?;
        let check = |q: usize| -> Result<usize, String> {
            if q >= n_qubits {
                Err(at(format!("qubit {q} out of range (have {n_qubits})")))
            } else {
                Ok(q)
            }
        };

        match op {
            "h" | "x" | "z" => {
                let [t] = expect_args(&rest, op).map_err(at)?;
                let t = check(parse_index(t).map_err(at)?)?;
                match op {
                    "h" => c.h(t),
                    "x" => c.x(t),
                    _ => c.z(t),
                };
            }
            "rx" | "ry" | "rz" => {
                let [theta, t] = expect_args(&rest, op).map_err(at)?;
                let theta = parse_angle(theta).map_err(at)?;
                let t = check(parse_index(t).map_err(at)?)?;
                match op {
                    "rx" => c.rx(theta, t),
                    "ry" => c.ry(theta, t),
                    _ => c.rz(theta, t),
                };
            }
            "cnot" => {
                let [ctrl, t] = expect_args(&rest, op).map_err(at)?;
                let ctrl = check(parse_index(ctrl).map_err(at)?)?;
                let t = check(parse_index(t).map_err(at)?)?;
                if ctrl == t {
                    return Err(at("`cnot` control and target must differ".into()));
                }
                c.cnot(ctrl, t);
            }
            "swap" => {
                let [a, b] = expect_args(&rest, op).map_err(at)?;
                let a = check(parse_index(a).map_err(at)?)?;
                let b = check(parse_index(b).map_err(at)?)?;
                c.swap(a, b);
            }
            "toffoli" => {
                let [c1, c2, t] = expect_args(&rest, op).map_err(at)?;
                let c1 = check(parse_index(c1).map_err(at)?)?;
                let c2 = check(parse_index(c2).map_err(at)?)?;
                let t = check(parse_index(t).map_err(at)?)?;
                if c1 == t || c2 == t || c1 == c2 {
                    return Err(at("`toffoli` controls and target must be distinct".into()));
                }
                c.toffoli(c1, c2, t);
            }
            "sample" => {
                if shots.is_some() {
                    return Err(at("duplicate `sample` directive".into()));
                }
                match rest.as_slice() {
                    [s] => shots = Some(parse_index(s).map_err(at)?),
                    [s, sd] => {
                        shots = Some(parse_index(s).map_err(at)?);
                        seed = sd
                            .parse::<u64>()
                            .map_err(|_| at(format!("invalid seed `{sd}`")))?;
                    }
                    _ => return Err(at("`sample` takes <shots> [seed]".into())),
                }
            }
            other => return Err(at(format!("unknown instruction `{other}`"))),
        }
    }

    let circuit = circuit.ok_or_else(|| "empty program: missing `qubits N`".to_string())?;
    Ok(Program {
        circuit,
        shots,
        seed,
    })
}

/// Drop everything from the first `#` onward.
fn strip_comment(line: &str) -> &str {
    match line.find('#') {
        Some(i) => &line[..i],
        None => line,
    }
}

/// Collect exactly `N` arguments, erroring if the count is wrong.
fn expect_args<'a, const N: usize>(rest: &[&'a str], op: &str) -> Result<[&'a str; N], String> {
    if rest.len() != N {
        return Err(format!(
            "`{op}` expects {N} argument(s), got {}",
            rest.len()
        ));
    }
    Ok(std::array::from_fn(|i| rest[i]))
}

/// Parse a non-negative qubit index or count.
fn parse_index(s: &str) -> Result<usize, String> {
    s.parse::<usize>()
        .map_err(|_| format!("expected a non-negative integer, got `{s}`"))
}

/// Parse a rotation angle: a plain float, or a `pi`-expression such as
/// `pi`, `-pi`, `pi/2`, `2pi`, `2*pi`, or `3pi/4`.
fn parse_angle(tok: &str) -> Result<f64, String> {
    let t = tok.trim();
    if let Ok(v) = t.parse::<f64>() {
        return Ok(v);
    }
    let (sign, body) = match t.strip_prefix('-') {
        Some(r) => (-1.0, r),
        None => (1.0, t.strip_prefix('+').unwrap_or(t)),
    };
    let (num_part, div) = match body.split_once('/') {
        Some((n, d)) => (
            n,
            d.parse::<f64>()
                .map_err(|_| format!("bad divisor in angle `{tok}`"))?,
        ),
        None => (body, 1.0),
    };
    let np = num_part.replace('*', "");
    let coeff = if np == "pi" {
        1.0
    } else if let Some(c) = np.strip_suffix("pi") {
        c.parse::<f64>()
            .map_err(|_| format!("bad coefficient in angle `{tok}`"))?
    } else {
        return Err(format!("unrecognized angle `{tok}`"));
    };
    Ok(sign * coeff * PI / div)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_bell_program() {
        let prog = parse("qubits 2\nh 0\ncnot 0 1\nsample 100 7").unwrap();
        assert_eq!(prog.shots, Some(100));
        assert_eq!(prog.seed, 7);
        let state = prog.circuit.run();
        // Bell state: |00> and |11> each with probability 1/2.
        assert!((state.probability(0) - 0.5).abs() < 1e-12);
        assert!((state.probability(3) - 0.5).abs() < 1e-12);
    }

    #[test]
    fn comments_and_blank_lines_ignored() {
        let src = "# header\n\nqubits 1  # inline comment\nx 0\n";
        let prog = parse(src).unwrap();
        assert!(prog.shots.is_none());
        assert_eq!(prog.seed, 0);
        let state = prog.circuit.run();
        assert!((state.probability(1) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn angle_forms() {
        assert!((parse_angle("pi").unwrap() - PI).abs() < 1e-12);
        assert!((parse_angle("-pi").unwrap() + PI).abs() < 1e-12);
        assert!((parse_angle("pi/2").unwrap() - PI / 2.0).abs() < 1e-12);
        assert!((parse_angle("2pi").unwrap() - 2.0 * PI).abs() < 1e-12);
        assert!((parse_angle("2*pi").unwrap() - 2.0 * PI).abs() < 1e-12);
        assert!((parse_angle("3pi/4").unwrap() - 3.0 * PI / 4.0).abs() < 1e-12);
        assert!((parse_angle("0.75").unwrap() - 0.75).abs() < 1e-12);
    }

    #[test]
    fn rx_pi_flips_qubit() {
        // Rx(pi) is X up to global phase, so it fully flips |0> to |1>.
        let prog = parse("qubits 1\nrx pi 0").unwrap();
        let state = prog.circuit.run();
        assert!((state.probability(1) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn missing_qubits_is_error() {
        assert!(parse("h 0").is_err());
        assert!(parse("").is_err());
    }

    #[test]
    fn out_of_range_qubit_is_error() {
        let err = parse("qubits 2\nh 5").err().unwrap();
        assert!(err.contains("out of range"), "got: {err}");
    }

    #[test]
    fn bad_instruction_reports_line() {
        let err = parse("qubits 2\nh 0\nfoo 1").err().unwrap();
        assert!(err.starts_with("line 3"), "got: {err}");
    }

    #[test]
    fn wrong_arity_is_error() {
        assert!(parse("qubits 2\ncnot 0").is_err());
        assert!(parse("qubits 2\nh 0 1").is_err());
    }

    #[test]
    fn duplicate_qubits_is_error() {
        assert!(parse("qubits 2\nqubits 3").is_err());
    }

    #[test]
    fn control_equals_target_is_error() {
        assert!(parse("qubits 2\ncnot 1 1").is_err());
        assert!(parse("qubits 3\ntoffoli 0 1 1").is_err());
    }
}
