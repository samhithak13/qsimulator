//! A small OpenQASM 2.0 subset importer.
//!
//! This parses the common, hand-written subset of OpenQASM 2.0 into a
//! [`Circuit`], enough to import textbook circuits and cross-check against
//! other tools. It is deliberately *not* a full implementation.
//!
//! Supported:
//! - `OPENQASM 2.0;` header and `include "...";` (both accepted and ignored).
//! - `qreg name[N];` — one or more; qubits are mapped into a single flat index
//!   space in declaration order (first declared register gets the lowest
//!   indices).
//! - `creg name[N];`, `barrier ...;`, and `measure ... -> ...;` are accepted
//!   and ignored (sampling is done separately via [`Circuit::sample`]).
//! - Gates: `x y z h s t sdg tdg` (1 qubit), `rx ry rz(theta)`, the phase gate
//!   `u1(lambda)` / `p(lambda)`, and the general `u2(phi,lambda)` /
//!   `u3(theta,phi,lambda)` (1 qubit); `cx cz swap` (2 qubits); `ccx`
//!   (3 qubits); `cy cz ch` (2 qubits); `cswap` (3 qubits); the controlled
//!   rotations `crz(theta)`, controlled phase `cu1(lambda)` / `cp(lambda)`, and
//!   controlled-U3 `cu3(theta,phi,lambda)` (2 qubits). Angles use the same
//!   syntax as the text program format (`pi`, `pi/2`, `-pi/4`, `2*pi`, or a
//!   float).
//! - `//` line comments and `/* ... */` block comments.
//!
//! Anything else — custom `gate` definitions, `if`, `reset`, etc. — is
//! reported as an unsupported-feature error rather than silently
//! mis-simulated.

use crate::error::ParseError;
use crate::program::parse_angle;
use crate::Circuit;
use std::collections::HashMap;

/// Largest register the importer will build. A dense state vector needs
/// `16·2^n` bytes, so this caps memory and stops a malformed huge `qreg` from
/// aborting the process on allocation.
const MAX_QUBITS: usize = 30;

/// Parse an OpenQASM 2.0 subset program into a [`Circuit`], or return a
/// [`ParseError`] naming the offending statement.
pub fn parse(src: &str) -> Result<Circuit, ParseError> {
    parse_inner(src).map_err(ParseError::new)
}

fn parse_inner(src: &str) -> Result<Circuit, String> {
    let clean = strip_comments(src);
    let statements: Vec<&str> = clean
        .split(';')
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .collect();

    // Pass 1: collect quantum registers and assign flat index offsets.
    let mut regs: HashMap<String, Reg> = HashMap::new();
    let mut total = 0usize;
    for stmt in &statements {
        if keyword(stmt) == "qreg" {
            let (name, size) = parse_reg_decl(stmt)?;
            if regs.contains_key(&name) {
                return Err(format!("duplicate register `{name}`"));
            }
            // Bound each register and the running total *before* accumulating,
            // so a huge size cannot overflow `total`.
            if size > MAX_QUBITS || total + size > MAX_QUBITS {
                return Err(format!(
                    "register total exceeds the maximum of {MAX_QUBITS} qubits"
                ));
            }
            regs.insert(
                name,
                Reg {
                    offset: total,
                    size,
                },
            );
            total += size;
        }
    }
    if total == 0 {
        return Err("no `qreg` declared".to_string());
    }
    let mut circuit = Circuit::new(total);

    // Pass 2: apply gates.
    for stmt in &statements {
        match keyword(stmt) {
            // Declarations and no-ops we accept and skip.
            "OPENQASM" | "include" | "qreg" | "creg" | "barrier" | "measure" => continue,
            // Features we deliberately reject rather than mis-simulate.
            "gate" | "opaque" | "if" | "reset" => {
                return Err(format!("unsupported OpenQASM feature `{}`", keyword(stmt)));
            }
            _ => apply_gate(&mut circuit, stmt, &regs)?,
        }
    }
    Ok(circuit)
}

/// A quantum register's placement in the flat qubit index space.
struct Reg {
    offset: usize,
    size: usize,
}

/// The leading keyword of a statement (up to the first whitespace, `(`, or `[`).
fn keyword(stmt: &str) -> &str {
    stmt.split(|c: char| c.is_whitespace() || c == '(' || c == '[')
        .next()
        .unwrap_or("")
}

/// Parse `qreg name[N]` (the trailing `;` already stripped).
fn parse_reg_decl(stmt: &str) -> Result<(String, usize), String> {
    let rest = stmt
        .strip_prefix("qreg")
        .ok_or_else(|| format!("malformed register declaration `{stmt}`"))?
        .trim();
    let open = rest
        .find('[')
        .ok_or_else(|| format!("register declaration needs a size: `{stmt}`"))?;
    let close = rest
        .find(']')
        .ok_or_else(|| format!("register declaration needs `]`: `{stmt}`"))?;
    if close < open {
        return Err(format!("malformed register declaration `{stmt}`"));
    }
    let name = rest[..open].trim().to_string();
    if name.is_empty() {
        return Err(format!("register needs a name: `{stmt}`"));
    }
    let size: usize = rest[open + 1..close]
        .trim()
        .parse()
        .map_err(|_| format!("invalid register size in `{stmt}`"))?;
    if size == 0 {
        return Err(format!("register size must be >= 1: `{stmt}`"));
    }
    Ok((name, size))
}

/// Parse and apply a single gate statement to `circuit`.
fn apply_gate(
    circuit: &mut Circuit,
    stmt: &str,
    regs: &HashMap<String, Reg>,
) -> Result<(), String> {
    let (head, operands_str) = split_head(stmt);
    let (name, params) = parse_head(head);
    let q = parse_operands(operands_str, regs)?;

    // Parse the parenthesized angle parameters, if any.
    let angles: Vec<f64> = match params {
        Some(p) if !p.trim().is_empty() => p
            .split(',')
            .map(|a| parse_angle(a.trim()))
            .collect::<Result<_, _>>()
            .map_err(|e| format!("{e} in `{stmt}`"))?,
        _ => Vec::new(),
    };

    // Validate the operand and angle counts for the named gate.
    let want = |nq: usize, na: usize| -> Result<(), String> {
        if q.len() != nq {
            return Err(format!(
                "`{name}` takes {nq} qubit(s), got {} in `{stmt}`",
                q.len()
            ));
        }
        if angles.len() != na {
            return Err(format!(
                "`{name}` takes {na} angle(s), got {} in `{stmt}`",
                angles.len()
            ));
        }
        Ok(())
    };

    match name {
        "x" => {
            want(1, 0)?;
            circuit.x(q[0]);
        }
        "y" => {
            want(1, 0)?;
            circuit.y(q[0]);
        }
        "z" => {
            want(1, 0)?;
            circuit.z(q[0]);
        }
        "h" => {
            want(1, 0)?;
            circuit.h(q[0]);
        }
        "s" => {
            want(1, 0)?;
            circuit.s(q[0]);
        }
        "t" => {
            want(1, 0)?;
            circuit.t(q[0]);
        }
        "sdg" => {
            want(1, 0)?;
            circuit.sdg(q[0]);
        }
        "tdg" => {
            want(1, 0)?;
            circuit.tdg(q[0]);
        }
        // `u1(lambda)` is the OpenQASM 2 phase gate; `p` is its OpenQASM 3 name.
        "u1" | "p" => {
            want(1, 1)?;
            circuit.p(angles[0], q[0]);
        }
        "u2" => {
            want(1, 2)?;
            circuit.u2(angles[0], angles[1], q[0]);
        }
        "u3" => {
            want(1, 3)?;
            circuit.u3(angles[0], angles[1], angles[2], q[0]);
        }
        "rx" => {
            want(1, 1)?;
            circuit.rx(angles[0], q[0]);
        }
        "ry" => {
            want(1, 1)?;
            circuit.ry(angles[0], q[0]);
        }
        "rz" => {
            want(1, 1)?;
            circuit.rz(angles[0], q[0]);
        }
        "cx" => {
            want(2, 0)?;
            require_distinct(&q, stmt)?;
            circuit.cnot(q[0], q[1]);
        }
        "cy" => {
            want(2, 0)?;
            require_distinct(&q, stmt)?;
            circuit.cy(q[0], q[1]);
        }
        "cz" => {
            want(2, 0)?;
            require_distinct(&q, stmt)?;
            circuit.cz(q[0], q[1]);
        }
        "ch" => {
            want(2, 0)?;
            require_distinct(&q, stmt)?;
            circuit.ch(q[0], q[1]);
        }
        "crz" => {
            want(2, 1)?;
            require_distinct(&q, stmt)?;
            circuit.crz(angles[0], q[0], q[1]);
        }
        // `cu1(lambda)` is the OpenQASM 2 controlled phase; `cp` is its
        // OpenQASM 3 name.
        "cu1" | "cp" => {
            want(2, 1)?;
            require_distinct(&q, stmt)?;
            circuit.cp(angles[0], q[0], q[1]);
        }
        "cu3" => {
            want(2, 3)?;
            require_distinct(&q, stmt)?;
            circuit.cu3(angles[0], angles[1], angles[2], q[0], q[1]);
        }
        "swap" => {
            want(2, 0)?;
            require_distinct(&q, stmt)?;
            circuit.swap(q[0], q[1]);
        }
        "cswap" => {
            want(3, 0)?;
            require_distinct(&q, stmt)?;
            circuit.cswap(q[0], q[1], q[2]);
        }
        "ccx" => {
            want(3, 0)?;
            require_distinct(&q, stmt)?;
            circuit.toffoli(q[0], q[1], q[2]);
        }
        other => return Err(format!("unsupported gate `{other}` in `{stmt}`")),
    }
    Ok(())
}

/// Reject repeated qubit operands (e.g. `cx q[0],q[0]`).
fn require_distinct(q: &[usize], stmt: &str) -> Result<(), String> {
    for i in 0..q.len() {
        for j in (i + 1)..q.len() {
            if q[i] == q[j] {
                return Err(format!("operands must be distinct qubits: `{stmt}`"));
            }
        }
    }
    Ok(())
}

/// Split a gate statement into its head (`gate` + optional `(params)`) and its
/// operand list, at the first top-level whitespace.
fn split_head(stmt: &str) -> (&str, &str) {
    let mut depth = 0i32;
    for (i, c) in stmt.char_indices() {
        match c {
            '(' => depth += 1,
            ')' => depth -= 1,
            _ if c.is_whitespace() && depth == 0 => {
                return (stmt[..i].trim(), stmt[i..].trim());
            }
            _ => {}
        }
    }
    (stmt.trim(), "")
}

/// Split a head into `(gate_name, Some(params))` or `(gate_name, None)`.
fn parse_head(head: &str) -> (&str, Option<&str>) {
    match (head.find('('), head.rfind(')')) {
        (Some(open), Some(close)) if close > open => {
            (head[..open].trim(), Some(head[open + 1..close].trim()))
        }
        _ => (head.trim(), None),
    }
}

/// Parse a comma-separated operand list like `q[0],q[1]` into flat qubit
/// indices using the register map.
fn parse_operands(operands: &str, regs: &HashMap<String, Reg>) -> Result<Vec<usize>, String> {
    if operands.trim().is_empty() {
        return Ok(Vec::new());
    }
    operands
        .split(',')
        .map(|op| parse_qubit_ref(op.trim(), regs))
        .collect()
}

/// Resolve a single `name[index]` reference to a flat qubit index.
fn parse_qubit_ref(op: &str, regs: &HashMap<String, Reg>) -> Result<usize, String> {
    let open = op
        .find('[')
        .ok_or_else(|| format!("expected a qubit reference like `q[0]`, got `{op}`"))?;
    let close = op
        .find(']')
        .ok_or_else(|| format!("qubit reference missing `]`: `{op}`"))?;
    if close < open {
        return Err(format!("malformed qubit reference `{op}`"));
    }
    let name = op[..open].trim();
    let index: usize = op[open + 1..close]
        .trim()
        .parse()
        .map_err(|_| format!("invalid qubit index in `{op}`"))?;
    let reg = regs
        .get(name)
        .ok_or_else(|| format!("unknown register `{name}`"))?;
    if index >= reg.size {
        return Err(format!(
            "qubit {name}[{index}] out of range (register size {})",
            reg.size
        ));
    }
    Ok(reg.offset + index)
}

/// Remove `//` line comments and `/* ... */` block comments (UTF-8 safe).
fn strip_comments(src: &str) -> String {
    let mut out = String::with_capacity(src.len());
    let mut chars = src.chars().peekable();
    let mut in_block = false;
    while let Some(c) = chars.next() {
        if in_block {
            if c == '*' && chars.peek() == Some(&'/') {
                chars.next();
                in_block = false;
            }
            continue;
        }
        if c == '/' {
            match chars.peek() {
                Some('*') => {
                    chars.next();
                    in_block = true;
                    continue;
                }
                Some('/') => {
                    // Consume the rest of the line, keeping the newline.
                    for nc in chars.by_ref() {
                        if nc == '\n' {
                            out.push('\n');
                            break;
                        }
                    }
                    continue;
                }
                _ => {}
            }
        }
        out.push(c);
    }
    out
}
