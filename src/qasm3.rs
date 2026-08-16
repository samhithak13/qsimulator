//! OpenQASM 3 support, by normalizing into the OpenQASM 2 subset.
//!
//! The two languages differ in how they *declare* things far more than in what
//! they *do*. Gate calls, `gate` declarations, angle expressions, `reset` and
//! conditionals are the same or nearly so, and the standard library
//! `stdgates.inc` overlaps `qelib1.inc` almost exactly — `p` and `cp` are the
//! OpenQASM 3 spellings of `u1` and `cu1`, which this engine already accepts.
//!
//! So rather than duplicate the gate table, the operand parsing and the
//! expression evaluator, this rewrites the handful of constructs that differ
//! and hands the result to [`crate::qasm::parse`]:
//!
//! | OpenQASM 3            | rewritten to          |
//! |-----------------------|-----------------------|
//! | `qubit[n] q;`         | `qreg q[n];`          |
//! | `qubit q;`            | `qreg q[1];`          |
//! | `bit[n] c;`           | `creg c[n];`          |
//! | `c[i] = measure q[j];`| `measure q[j] -> c[i];` |
//! | `if (c == v) { a; b; }` | `if(c==v) a; if(c==v) b;` |
//!
//! The last is exact rather than approximate: a conditional block holds only
//! gates, so nothing inside it can change the value being tested, and guarding
//! each statement is the same as guarding the block.
//!
//! Not supported: `for`, `while`, `def`, `let`, classical arithmetic on bits,
//! physical qubits (`$0`), and timing. These are the parts of OpenQASM 3 that
//! genuinely have no OpenQASM 2 counterpart, and each is reported rather than
//! ignored.

use crate::error::ParseError;
use crate::qasm;
use crate::Circuit;

/// Parse an OpenQASM 3 program into a [`Circuit`].
pub fn parse(src: &str) -> Result<Circuit, ParseError> {
    let normalized = normalize(src).map_err(ParseError::new)?;
    qasm::parse(&normalized)
}

/// Whether `src` looks like OpenQASM 3 rather than 2, by its version header.
///
/// Comments are stripped first: a line *inside* a comment can begin with
/// `OPENQASM`, and reading it would route a perfectly good OpenQASM 2 file to
/// the wrong front end.
pub fn is_openqasm3(src: &str) -> bool {
    let clean = qasm::strip_comments(src);
    for line in clean.lines() {
        if let Some(rest) = line.trim().strip_prefix("OPENQASM") {
            return rest.trim_start().starts_with('3');
        }
    }
    false
}

/// Rewrite OpenQASM 3 declarations and statements into the OpenQASM 2 subset.
fn normalize(src: &str) -> Result<String, String> {
    let clean = qasm::strip_comments(src);
    let statements = qasm::split_statements(&clean)?;
    let mut out = String::from("OPENQASM 2.0;\n");

    for stmt in statements {
        // The version header is replaced by the one written above.
        if stmt.starts_with("OPENQASM") {
            continue;
        }
        for rewritten in rewrite(stmt)? {
            out.push_str(&rewritten);
            out.push_str(";\n");
        }
    }
    Ok(out)
}

/// Rewrite one statement, which a conditional block may expand into several.
fn rewrite(stmt: &str) -> Result<Vec<String>, String> {
    let head = stmt
        .split(|c: char| c.is_whitespace() || c == '(' || c == '[')
        .next()
        .unwrap_or("");

    match head {
        // Constructs with no OpenQASM 2 counterpart, named rather than ignored.
        // A statement ends at its closing brace, so an `else` clause arrives
        // here as a statement of its own rather than as part of the `if`.
        "for" | "while" | "def" | "defcal" | "let" | "cal" | "extern" | "duration" | "stretch"
        | "delay" | "box" | "else" => Err(format!(
            "unsupported OpenQASM 3 construct `{head}`; this engine covers the \
             subset that maps onto OpenQASM 2"
        )),
        "qubit" => Ok(vec![declaration(stmt, "qubit", "qreg")?]),
        "bit" => Ok(vec![declaration(stmt, "bit", "creg")?]),
        "if" => conditional(stmt),
        // `gate`, `reset`, `barrier`, `include` and every gate call already read
        // the same in both languages.
        _ => match measurement_assignment(stmt) {
            Some(rewritten) => Ok(vec![rewritten]),
            None => Ok(vec![stmt.to_string()]),
        },
    }
}

/// `qubit[n] q` -> `qreg q[n]`, and the unsized `qubit q` -> `qreg q[1]`.
fn declaration(stmt: &str, from: &str, to: &str) -> Result<String, String> {
    let rest = stmt
        .strip_prefix(from)
        .ok_or_else(|| format!("malformed declaration `{stmt}`"))?
        .trim();
    let (size, name) = match rest.strip_prefix('[') {
        Some(after) => {
            let close = after
                .find(']')
                .ok_or_else(|| format!("declaration needs `]`: `{stmt}`"))?;
            (after[..close].trim().to_string(), after[close + 1..].trim())
        }
        // An unsized declaration is a single qubit or bit.
        None => ("1".to_string(), rest),
    };
    if name.is_empty() {
        return Err(format!("declaration needs a name: `{stmt}`"));
    }
    Ok(format!("{to} {name}[{size}]"))
}

/// `c[i] = measure q[j]` -> `measure q[j] -> c[i]`, if that is what this is.
fn measurement_assignment(stmt: &str) -> Option<String> {
    let (target, rest) = stmt.split_once('=')?;
    // Do not mistake `==` inside some other statement for an assignment.
    let source = rest.strip_prefix(' ').unwrap_or(rest).trim();
    let source = source.strip_prefix("measure")?;
    Some(format!("measure {} -> {}", source.trim(), target.trim()))
}

/// `if (cond) { a; b; }` -> one guarded statement per statement in the block.
fn conditional(stmt: &str) -> Result<Vec<String>, String> {
    let open = stmt
        .find('{')
        .ok_or_else(|| format!("conditional needs a `{{ ... }}` block: `{stmt}`"))?;
    let close = stmt
        .rfind('}')
        .ok_or_else(|| format!("conditional needs `}}`: `{stmt}`"))?;
    if close < open {
        return Err(format!("malformed conditional `{stmt}`"));
    }
    let condition = stmt[..open].trim();

    let mut guarded = Vec::new();
    for inner in stmt[open + 1..close].split(';') {
        let inner = inner.trim();
        if inner.is_empty() {
            continue;
        }
        // Guarding each statement matches guarding the block, since a block
        // holds only gates and so cannot change the value being tested.
        guarded.push(format!("{condition} {inner}"));
    }
    if guarded.is_empty() {
        return Err(format!("conditional has no statements: `{stmt}`"));
    }
    Ok(guarded)
}
