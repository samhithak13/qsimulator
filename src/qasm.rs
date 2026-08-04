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
//! - Gates: `x y z h s t` (1 qubit), `rx ry rz(theta)` (1 qubit),
//!   `cx cz swap` (2 qubits), `ccx` (3 qubits). Angles use the same syntax as
//!   the text program format (`pi`, `pi/2`, `-pi/4`, `2*pi`, or a float).
//! - `//` line comments and `/* ... */` block comments.
//!
//! Anything else — custom `gate` definitions, `if`, `reset`, `u1/u2/u3`,
//! `sdg/tdg`, etc. — is reported as an unsupported-feature error rather than
//! silently mis-simulated.

use crate::program::parse_angle;
use crate::Circuit;
use std::collections::HashMap;

/// Parse an OpenQASM 2.0 subset program into a [`Circuit`], or return a
/// human-readable error naming the offending statement.
pub fn parse(src: &str) -> Result<Circuit, String> {
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

    // Validate the operand count for the named gate.
    let want = |k: usize| -> Result<(), String> {
        if q.len() == k {
            Ok(())
        } else {
            Err(format!(
                "`{name}` takes {k} qubit(s), got {} in `{stmt}`",
                q.len()
            ))
        }
    };
    let angle = || -> Result<f64, String> {
        let p = params.ok_or_else(|| format!("`{name}` needs an angle in `{stmt}`"))?;
        parse_angle(p.trim()).map_err(|e| format!("{e} in `{stmt}`"))
    };

    match name {
        "x" => {
            want(1)?;
            circuit.x(q[0]);
        }
        "y" => {
            want(1)?;
            circuit.y(q[0]);
        }
        "z" => {
            want(1)?;
            circuit.z(q[0]);
        }
        "h" => {
            want(1)?;
            circuit.h(q[0]);
        }
        "s" => {
            want(1)?;
            circuit.s(q[0]);
        }
        "t" => {
            want(1)?;
            circuit.t(q[0]);
        }
        "rx" => {
            want(1)?;
            circuit.rx(angle()?, q[0]);
        }
        "ry" => {
            want(1)?;
            circuit.ry(angle()?, q[0]);
        }
        "rz" => {
            want(1)?;
            circuit.rz(angle()?, q[0]);
        }
        "cx" => {
            want(2)?;
            require_distinct(&q, stmt)?;
            circuit.cnot(q[0], q[1]);
        }
        "cz" => {
            want(2)?;
            require_distinct(&q, stmt)?;
            circuit.cz(q[0], q[1]);
        }
        "swap" => {
            want(2)?;
            require_distinct(&q, stmt)?;
            circuit.swap(q[0], q[1]);
        }
        "ccx" => {
            want(3)?;
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
