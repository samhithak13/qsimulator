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
//! - `creg name[N];` and `barrier ...;` are accepted and ignored.
//! - `measure q[i] -> c[j];` collapses qubit `i` onto a sampled outcome. The
//!   classical bit is not modelled — nothing in this subset can read one back —
//!   but the collapse is what makes a following gate see a definite state. A
//!   circuit containing one is stochastic: see [`Circuit::run_seeded`].
//! - Gates: `id x y z h s t sdg tdg` (1 qubit), `rx ry rz(theta)`, the phase gate
//!   `u1(lambda)` / `p(lambda)`, and the general `u2(phi,lambda)` /
//!   `u3(theta,phi,lambda)` (1 qubit); `cx cz swap` (2 qubits); `ccx`
//!   (3 qubits); `cy cz ch` (2 qubits); `cswap` (3 qubits); the controlled
//!   rotations `crz(theta)`, controlled phase `cu1(lambda)` / `cp(lambda)`, and
//!   controlled-U3 `cu3(theta,phi,lambda)` (2 qubits). Angles are arithmetic
//!   expressions, shared with the text program format: numbers, `pi`, the
//!   operators `+ - * / ^`, parentheses, and
//!   `sin`/`cos`/`tan`/`exp`/`ln`/`sqrt`.
//! - `U(theta,phi,lambda)` and `CX`, the two OpenQASM 2 primitives.
//! - `gate name(params) qargs { ... }` declarations, expanded at each call
//!   site with the actual angles and qubits substituted in. Bodies may call
//!   other declarations, and their angles are expressions over the formal
//!   parameters. A file's own declaration wins over the built-in of the same
//!   name, so a program written against the primitives behaves as written.
//! - `//` line comments and `/* ... */` block comments.
//!
//! - The rest of `qelib1.inc`: `u`, `u0`, `sx`, `sxdg`, `crx`, `cry`, `csx`,
//!   `cu`, `rxx`, `rzz`, the relative-phase Toffolis `rccx` and `rc3x`, and
//!   `c3x`, `c3sqrtx`, `c4x`. `include` is *ignored* rather than honoured, so
//!   these are implemented here rather than read from the header. Each is
//!   checked against Qiskit's unitary with its phase; `rxx` and `rzz` follow
//!   qelib1's decomposition and so differ from Qiskit's gate object by a global
//!   phase of `theta/2`, which is unobservable and which neither OpenQASM 2 nor
//!   this engine's circuit representation can express.
//!
//! - `reset q[i];` and `reset q;` collapse a qubit and force it to |0>.
//!
//! Anything else — `if`, `opaque` — is reported as an unsupported-feature
//! error rather than silently mis-simulated.

use crate::error::ParseError;
use crate::gates;
use crate::Circuit;
use std::collections::HashMap;
use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI};

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
    let statements = split_statements(&clean)?;

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

    // Pass 2: collect `gate` declarations, so a call may precede its
    // declaration textually even though OpenQASM does not require that.
    let mut defs: HashMap<String, GateDef> = HashMap::new();
    for stmt in &statements {
        if keyword(stmt) == "gate" {
            let (name, def) = parse_gate_def(stmt)?;
            if defs.contains_key(&name) {
                return Err(format!("duplicate `gate {name}` declaration"));
            }
            defs.insert(name, def);
        }
    }

    // Pass 3: apply gates.
    let mut budget = MAX_EXPANDED_OPS;
    for stmt in &statements {
        match keyword(stmt) {
            // Declarations and no-ops we accept and skip.
            "OPENQASM" | "include" | "qreg" | "creg" | "barrier" | "gate" => continue,
            // A measurement collapses the register, so ignoring one would
            // silently give the wrong answer for anything that follows it.
            "measure" => apply_measure(&mut circuit, stmt, &regs)?,
            "reset" => apply_reset(&mut circuit, stmt, &regs)?,
            // Features we deliberately reject rather than mis-simulate.
            "opaque" | "if" => {
                return Err(format!("unsupported OpenQASM feature `{}`", keyword(stmt)));
            }
            _ => apply_gate(&mut circuit, stmt, &regs, &defs, &mut budget)?,
        }
    }
    Ok(circuit)
}

/// Split a program into statements.
///
/// A statement normally ends at `;`, but a `gate` declaration's `{ ... }` body
/// holds semicolons of its own and ends at the closing brace with no `;` after
/// it. Splitting naively on `;` therefore glues that trailing `}` onto whatever
/// follows — and since `gate` blocks conventionally precede `qreg`, the register
/// declaration would vanish into a fragment starting with `}`.
fn split_statements(src: &str) -> Result<Vec<&str>, String> {
    let mut statements = Vec::new();
    let mut depth = 0usize;
    let mut start = 0usize;
    for (i, c) in src.char_indices() {
        match c {
            '{' => depth += 1,
            '}' => {
                depth = depth
                    .checked_sub(1)
                    .ok_or_else(|| "unexpected `}`".to_string())?;
                // The brace closes the statement it belongs to.
                if depth == 0 {
                    statements.push(&src[start..=i]);
                    start = i + c.len_utf8();
                }
            }
            ';' if depth == 0 => {
                statements.push(&src[start..i]);
                start = i + c.len_utf8();
            }
            _ => {}
        }
    }
    if depth != 0 {
        return Err("unterminated `{` (missing `}`)".to_string());
    }
    statements.push(&src[start..]);
    Ok(statements
        .into_iter()
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .collect())
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

/// A user `gate` declaration: its formal angle parameters, its formal qubit
/// arguments, and the statements of its body, all unresolved.
struct GateDef {
    params: Vec<String>,
    qargs: Vec<String>,
    body: Vec<String>,
}

/// How deeply user gates may nest. OpenQASM 2 forbids recursion (a body may
/// only call gates already declared), but a malformed file can still describe a
/// cycle, and this bounds it.
const MAX_GATE_DEPTH: usize = 64;

/// How many primitive gates a program may expand to. A short file can describe
/// an exponential expansion — each of `n` definitions calling the previous one
/// twice is `2^n` gates — so the budget, not the input size, is what bounds the
/// work and the memory it would take. Far above any real circuit this engine
/// could simulate anyway: 30 qubits is the register ceiling.
const MAX_EXPANDED_OPS: usize = 100_000;

/// Parse `gate name(params) qargs { body }` into its declaration.
fn parse_gate_def(stmt: &str) -> Result<(String, GateDef), String> {
    let rest = stmt
        .strip_prefix("gate")
        .ok_or_else(|| format!("malformed gate declaration `{stmt}`"))?;
    let open = rest
        .find('{')
        .ok_or_else(|| format!("gate declaration needs a `{{ ... }}` body: `{stmt}`"))?;
    let close = rest
        .rfind('}')
        .ok_or_else(|| format!("gate declaration needs `}}`: `{stmt}`"))?;
    if close < open {
        return Err(format!("malformed gate declaration `{stmt}`"));
    }

    let (head, qargs_str) = split_head(rest[..open].trim());
    let (name, params) = parse_head(head);
    if name.is_empty() {
        return Err(format!("gate declaration needs a name: `{stmt}`"));
    }
    let params = parse_identifiers(params.unwrap_or(""), "parameter", stmt)?;
    let qargs = parse_identifiers(qargs_str, "qubit argument", stmt)?;
    if qargs.is_empty() {
        return Err(format!("`gate {name}` needs at least one qubit argument"));
    }

    let body = rest[open + 1..close]
        .split(';')
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .map(str::to_string)
        .collect();
    Ok((
        name.to_string(),
        GateDef {
            params,
            qargs,
            body,
        },
    ))
}

/// Parse `measure q[i] -> c[j]` and apply the collapse. The classical target
/// is checked for shape but not modelled.
fn apply_measure(
    circuit: &mut Circuit,
    stmt: &str,
    regs: &HashMap<String, Reg>,
) -> Result<(), String> {
    let rest = stmt
        .strip_prefix("measure")
        .ok_or_else(|| format!("malformed measurement `{stmt}`"))?;
    let (source, dest) = rest
        .split_once("->")
        .ok_or_else(|| format!("measurement needs `-> target`: `{stmt}`"))?;
    if dest.trim().is_empty() {
        return Err(format!("measurement needs a classical target: `{stmt}`"));
    }
    // `measure q[i] -> c[j]` measures one qubit; `measure q -> c` measures the
    // whole register, which is how most hand-written programs end.
    for qubit in resolve_qubits(source.trim(), regs, stmt)? {
        circuit.measure(qubit);
    }
    Ok(())
}

/// Parse `reset q[i]` or `reset q` and apply the collapse-to-|0>.
fn apply_reset(
    circuit: &mut Circuit,
    stmt: &str,
    regs: &HashMap<String, Reg>,
) -> Result<(), String> {
    let target = stmt
        .strip_prefix("reset")
        .ok_or_else(|| format!("malformed reset `{stmt}`"))?
        .trim();
    for qubit in resolve_qubits(target, regs, stmt)? {
        circuit.reset(qubit);
    }
    Ok(())
}

/// Resolve `q[i]` to that one qubit, or a bare `q` to the whole register.
fn resolve_qubits(
    operand: &str,
    regs: &HashMap<String, Reg>,
    stmt: &str,
) -> Result<Vec<usize>, String> {
    if operand.is_empty() {
        return Err(format!("expected a qubit in `{stmt}`"));
    }
    if operand.contains('[') {
        let qubits = parse_operands(operand, regs)?;
        if qubits.len() != 1 {
            return Err(format!(
                "expected one qubit, got {} in `{stmt}`",
                qubits.len()
            ));
        }
        return Ok(qubits);
    }
    let reg = regs
        .get(operand)
        .ok_or_else(|| format!("unknown register `{operand}` in `{stmt}`"))?;
    Ok((0..reg.size).map(|i| reg.offset + i).collect())
}

/// Parse a comma-separated list of formal names, rejecting duplicates and
/// anything that is not an identifier.
fn parse_identifiers(list: &str, kind: &str, stmt: &str) -> Result<Vec<String>, String> {
    let list = list.trim();
    if list.is_empty() {
        return Ok(Vec::new());
    }
    let mut names = Vec::new();
    for raw in list.split(',') {
        let name = raw.trim();
        let is_identifier = !name.is_empty()
            && name.starts_with(|c: char| c.is_ascii_alphabetic() || c == '_')
            && name.chars().all(|c| c.is_ascii_alphanumeric() || c == '_');
        if !is_identifier {
            return Err(format!("invalid {kind} `{name}` in `{stmt}`"));
        }
        if names.iter().any(|n| n == name) {
            return Err(format!("duplicate {kind} `{name}` in `{stmt}`"));
        }
        names.push(name.to_string());
    }
    Ok(names)
}

/// Evaluate a gate's parenthesised angle list against `vars`.
fn parse_angle_list(
    params: Option<&str>,
    vars: &HashMap<String, f64>,
    stmt: &str,
) -> Result<Vec<f64>, String> {
    match params {
        Some(p) if !p.trim().is_empty() => p
            .split(',')
            .map(|a| crate::expr::eval(a.trim(), vars))
            .collect::<Result<_, _>>()
            .map_err(|e| format!("{e} in `{stmt}`")),
        _ => Ok(Vec::new()),
    }
}

/// Parse and apply a single gate statement to `circuit`, resolving its operands
/// against the register map.
fn apply_gate(
    circuit: &mut Circuit,
    stmt: &str,
    regs: &HashMap<String, Reg>,
    defs: &HashMap<String, GateDef>,
    budget: &mut usize,
) -> Result<(), String> {
    let (head, operands_str) = split_head(stmt);
    let (name, params) = parse_head(head);
    let q = parse_operands(operands_str, regs)?;
    let angles = parse_angle_list(params, &HashMap::new(), stmt)?;
    apply_resolved(circuit, name, &angles, &q, defs, budget, 0, stmt)
}

/// Apply a gate whose angles and qubits are already resolved to values —
/// expanding it if it is a user declaration, otherwise applying the builtin.
#[allow(clippy::too_many_arguments)]
fn apply_resolved(
    circuit: &mut Circuit,
    name: &str,
    angles: &[f64],
    q: &[usize],
    defs: &HashMap<String, GateDef>,
    budget: &mut usize,
    depth: usize,
    stmt: &str,
) -> Result<(), String> {
    // A file's own declarations win, so a program that defines `cx` in terms of
    // the OpenQASM primitives (the way qelib1 itself does) behaves as written.
    if let Some(def) = defs.get(name) {
        return expand_user_gate(circuit, name, def, angles, q, defs, budget, depth, stmt);
    }
    *budget = budget.checked_sub(1).ok_or_else(|| {
        format!("program expands to more than {MAX_EXPANDED_OPS} gates; is a `gate` recursive?")
    })?;
    apply_builtin(circuit, name, angles, q, stmt)
}

/// Substitute a call's actual angles and qubits into a declaration's body and
/// apply each statement.
#[allow(clippy::too_many_arguments)]
fn expand_user_gate(
    circuit: &mut Circuit,
    name: &str,
    def: &GateDef,
    angles: &[f64],
    q: &[usize],
    defs: &HashMap<String, GateDef>,
    budget: &mut usize,
    depth: usize,
    stmt: &str,
) -> Result<(), String> {
    if depth >= MAX_GATE_DEPTH {
        return Err(format!(
            "`gate` calls nested more than {MAX_GATE_DEPTH} deep at `{name}`; is it recursive?"
        ));
    }
    if def.params.len() != angles.len() || def.qargs.len() != q.len() {
        return Err(format!(
            "`{name}` takes {} qubit(s) and {} angle(s), got {} and {} in `{stmt}`",
            def.qargs.len(),
            def.params.len(),
            q.len(),
            angles.len()
        ));
    }

    let vars: HashMap<String, f64> = def
        .params
        .iter()
        .cloned()
        .zip(angles.iter().copied())
        .collect();
    let qubits: HashMap<&str, usize> = def
        .qargs
        .iter()
        .map(String::as_str)
        .zip(q.iter().copied())
        .collect();

    for body_stmt in &def.body {
        let (head, operands_str) = split_head(body_stmt);
        let (body_name, body_params) = parse_head(head);
        // `barrier` is scheduling only, and carries no operands we must resolve.
        if body_name == "barrier" {
            continue;
        }
        let body_angles = parse_angle_list(body_params, &vars, body_stmt)?;
        let body_qubits = operands_str
            .split(',')
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .map(|operand| {
                qubits.get(operand).copied().ok_or_else(|| {
                    format!("`{operand}` is not a qubit argument of `{name}` in `{body_stmt}`")
                })
            })
            .collect::<Result<Vec<_>, _>>()?;
        apply_resolved(
            circuit,
            body_name,
            &body_angles,
            &body_qubits,
            defs,
            budget,
            depth + 1,
            body_stmt,
        )?;
    }
    Ok(())
}

/// Apply one of the built-in gates by name.
fn apply_builtin(
    circuit: &mut Circuit,
    name: &str,
    angles: &[f64],
    q: &[usize],
    stmt: &str,
) -> Result<(), String> {
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
        "id" => {
            want(1, 0)?;
            circuit.id(q[0]);
        }
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
        // `U` and `CX` are the two OpenQASM 2 primitives that qelib1 itself is
        // defined in terms of, so a program that pastes those definitions in
        // rather than relying on the include still works.
        "u3" | "U" => {
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
        "cx" | "CX" => {
            want(2, 0)?;
            require_distinct(q, stmt)?;
            circuit.cnot(q[0], q[1]);
        }
        "cy" => {
            want(2, 0)?;
            require_distinct(q, stmt)?;
            circuit.cy(q[0], q[1]);
        }
        "cz" => {
            want(2, 0)?;
            require_distinct(q, stmt)?;
            circuit.cz(q[0], q[1]);
        }
        "ch" => {
            want(2, 0)?;
            require_distinct(q, stmt)?;
            circuit.ch(q[0], q[1]);
        }
        "crz" => {
            want(2, 1)?;
            require_distinct(q, stmt)?;
            circuit.crz(angles[0], q[0], q[1]);
        }
        // `cu1(lambda)` is the OpenQASM 2 controlled phase; `cp` is its
        // OpenQASM 3 name.
        "cu1" | "cp" => {
            want(2, 1)?;
            require_distinct(q, stmt)?;
            circuit.cp(angles[0], q[0], q[1]);
        }
        "cu3" => {
            want(2, 3)?;
            require_distinct(q, stmt)?;
            circuit.cu3(angles[0], angles[1], angles[2], q[0], q[1]);
        }
        "swap" => {
            want(2, 0)?;
            require_distinct(q, stmt)?;
            circuit.swap(q[0], q[1]);
        }
        "cswap" => {
            want(3, 0)?;
            require_distinct(q, stmt)?;
            circuit.cswap(q[0], q[1], q[2]);
        }
        "ccx" => {
            want(3, 0)?;
            require_distinct(q, stmt)?;
            circuit.toffoli(q[0], q[1], q[2]);
        }

        // --- The rest of qelib1 ---
        //
        // `include` is ignored rather than honoured, so these are implemented
        // here instead of being read from the header. Where a gate is exactly
        // one of this engine's primitives it is applied as one — `crx` really
        // is a controlled Rx — which keeps diagrams and re-export readable.
        // The rest follow qelib1's own decomposition, and every one of them is
        // checked against Qiskit's unitary, phase included.
        "u" => {
            want(1, 3)?;
            circuit.u3(angles[0], angles[1], angles[2], q[0]);
        }
        // An idle of duration `gamma`: no rotation, only a delay this engine
        // has no notion of.
        "u0" => {
            want(1, 1)?;
            circuit.id(q[0]);
        }
        "sx" => {
            want(1, 0)?;
            circuit.sx(q[0]);
        }
        "sxdg" => {
            want(1, 0)?;
            circuit.sxdg(q[0]);
        }
        "crx" => {
            want(2, 1)?;
            require_distinct(q, stmt)?;
            circuit.cu(gates::rx(angles[0]), q[0], q[1]);
        }
        "cry" => {
            want(2, 1)?;
            require_distinct(q, stmt)?;
            circuit.cu(gates::ry(angles[0]), q[0], q[1]);
        }
        "csx" => {
            want(2, 0)?;
            require_distinct(q, stmt)?;
            circuit.cu(gates::sx(), q[0], q[1]);
        }
        // Controlled-U3 with an extra phase `gamma` on the control — the same
        // shape `to_qasm` emits for an arbitrary controlled-U.
        "cu" => {
            want(2, 4)?;
            require_distinct(q, stmt)?;
            circuit.p(angles[3], q[0]);
            circuit.cu3(angles[0], angles[1], angles[2], q[0], q[1]);
        }
        "rxx" => {
            want(2, 1)?;
            require_distinct(q, stmt)?;
            let theta = angles[0];
            circuit.u3(FRAC_PI_2, theta, 0.0, q[0]).h(q[1]);
            circuit.cnot(q[0], q[1]).p(-theta, q[1]).cnot(q[0], q[1]);
            circuit.h(q[1]).u2(-PI, PI - theta, q[0]);
        }
        "rzz" => {
            want(2, 1)?;
            require_distinct(q, stmt)?;
            circuit.cnot(q[0], q[1]).p(angles[0], q[1]).cnot(q[0], q[1]);
        }
        // A *relative-phase* Toffoli: it differs from `ccx` by phases on some
        // basis states, so it is not interchangeable with one.
        "rccx" => {
            want(3, 0)?;
            require_distinct(q, stmt)?;
            let (a, b, c) = (q[0], q[1], q[2]);
            circuit.u2(0.0, PI, c).p(FRAC_PI_4, c);
            circuit.cnot(b, c).p(-FRAC_PI_4, c);
            circuit.cnot(a, c).p(FRAC_PI_4, c);
            circuit.cnot(b, c).p(-FRAC_PI_4, c);
            circuit.u2(0.0, PI, c);
        }
        // The relative-phase three-control X, likewise not a plain `c3x`.
        "rc3x" => {
            want(4, 0)?;
            require_distinct(q, stmt)?;
            let (a, b, c, d) = (q[0], q[1], q[2], q[3]);
            circuit.u2(0.0, PI, d).p(FRAC_PI_4, d);
            circuit.cnot(c, d).p(-FRAC_PI_4, d).u2(0.0, PI, d);
            circuit.cnot(a, d).p(FRAC_PI_4, d);
            circuit.cnot(b, d).p(-FRAC_PI_4, d);
            circuit.cnot(a, d).p(FRAC_PI_4, d);
            circuit.cnot(b, d).p(-FRAC_PI_4, d);
            circuit.u2(0.0, PI, d).p(FRAC_PI_4, d);
            circuit.cnot(c, d).p(-FRAC_PI_4, d).u2(0.0, PI, d);
        }
        "c3x" => {
            want(4, 0)?;
            require_distinct(q, stmt)?;
            circuit.mcx(&q[..3], q[3]);
        }
        "c3sqrtx" => {
            want(4, 0)?;
            require_distinct(q, stmt)?;
            circuit.mcu(gates::sx(), &q[..3], q[3]);
        }
        "c4x" => {
            want(5, 0)?;
            require_distinct(q, stmt)?;
            circuit.mcx(&q[..4], q[4]);
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
