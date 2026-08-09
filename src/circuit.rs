//! Circuit builder and executor.

use crate::gates;
use crate::rng::Rng;
use crate::state::State;
use num_complex::Complex64;
use std::collections::HashMap;
use std::fmt;

type Gate = [[Complex64; 2]; 2];

/// Reason a circuit could not be exported to the OpenQASM 2.0 subset by
/// [`Circuit::to_qasm`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ExportError {
    /// An arbitrary single-qubit gate with no OpenQASM name. Built-in gates
    /// never produce this; it can only arise from a manually built `Op`.
    SingleGate {
        /// The gate's diagram label.
        label: &'static str,
    },
}

impl fmt::Display for ExportError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ExportError::SingleGate { label } => {
                write!(f, "cannot export single-qubit gate `{label}` to OpenQASM 2")
            }
        }
    }
}

impl std::error::Error for ExportError {}

/// A single instruction in a circuit.
///
/// Each gate-bearing variant carries a short `label` (e.g. `"H"`, `"RX"`)
/// used for diagram rendering and QASM export; it never affects execution.
/// `Single` additionally records the gate's angle `params` (empty for
/// non-parametric gates, one for rotations/phase, up to three for `u3`) so
/// that export is lossless.
#[derive(Debug, Clone)]
enum Op {
    Single {
        gate: Gate,
        target: usize,
        label: &'static str,
        params: Vec<f64>,
    },
    Controlled {
        gate: Gate,
        control: usize,
        target: usize,
        label: &'static str,
        params: Vec<f64>,
    },
    Swap {
        a: usize,
        b: usize,
    },
    MultiControlled {
        gate: Gate,
        controls: Vec<usize>,
        target: usize,
        label: &'static str,
    },
}

/// A quantum circuit: an ordered list of gate operations on `n_qubits`.
#[derive(Debug, Clone)]
pub struct Circuit {
    n_qubits: usize,
    ops: Vec<Op>,
}

impl Circuit {
    /// Create an empty circuit over `n_qubits` qubits.
    pub fn new(n_qubits: usize) -> Self {
        Circuit {
            n_qubits,
            ops: Vec::new(),
        }
    }

    /// Hadamard gate on `target`.
    pub fn h(&mut self, target: usize) -> &mut Self {
        self.ops.push(Op::Single {
            gate: gates::h(),
            target,
            label: "H",
            params: Vec::new(),
        });
        self
    }

    /// Pauli-X (NOT) gate on `target`.
    pub fn x(&mut self, target: usize) -> &mut Self {
        self.ops.push(Op::Single {
            gate: gates::x(),
            target,
            label: "X",
            params: Vec::new(),
        });
        self
    }

    /// Pauli-Z gate on `target`.
    pub fn z(&mut self, target: usize) -> &mut Self {
        self.ops.push(Op::Single {
            gate: gates::z(),
            target,
            label: "Z",
            params: Vec::new(),
        });
        self
    }

    /// Identity gate on `target`: a no-op, kept so a circuit can record an
    /// explicit idle step (OpenQASM's `id`).
    pub fn id(&mut self, target: usize) -> &mut Self {
        self.ops.push(Op::Single {
            gate: gates::id(),
            target,
            label: "I",
            params: Vec::new(),
        });
        self
    }

    /// Square root of X on `target` (OpenQASM `sx`): `sx` twice is X.
    pub fn sx(&mut self, target: usize) -> &mut Self {
        self.ops.push(Op::Single {
            gate: gates::sx(),
            target,
            label: "SX",
            params: Vec::new(),
        });
        self
    }

    /// Inverse square root of X on `target` (OpenQASM `sxdg`).
    pub fn sxdg(&mut self, target: usize) -> &mut Self {
        self.ops.push(Op::Single {
            gate: gates::sxdg(),
            target,
            label: "SXDG",
            params: Vec::new(),
        });
        self
    }

    /// Pauli-Y gate on `target`.
    pub fn y(&mut self, target: usize) -> &mut Self {
        self.ops.push(Op::Single {
            gate: gates::y(),
            target,
            label: "Y",
            params: Vec::new(),
        });
        self
    }

    /// Phase gate S = diag(1, i).
    pub fn s(&mut self, target: usize) -> &mut Self {
        self.ops.push(Op::Single {
            gate: gates::s(),
            target,
            label: "S",
            params: Vec::new(),
        });
        self
    }

    /// T gate = diag(1, e^{i pi/4}).
    pub fn t(&mut self, target: usize) -> &mut Self {
        self.ops.push(Op::Single {
            gate: gates::t(),
            target,
            label: "T",
            params: Vec::new(),
        });
        self
    }

    /// S-dagger gate = S† = diag(1, -i).
    pub fn sdg(&mut self, target: usize) -> &mut Self {
        self.ops.push(Op::Single {
            gate: gates::sdg(),
            target,
            label: "SDG",
            params: Vec::new(),
        });
        self
    }

    /// T-dagger gate = T† = diag(1, e^{-i pi/4}).
    pub fn tdg(&mut self, target: usize) -> &mut Self {
        self.ops.push(Op::Single {
            gate: gates::tdg(),
            target,
            label: "TDG",
            params: Vec::new(),
        });
        self
    }

    /// Phase gate P(lambda) = diag(1, e^{i·lambda}) on `target`.
    ///
    /// Generalizes Z/S/T (`p(π)` = Z, `p(π/2)` = S, `p(π/4)` = T) and maps to
    /// OpenQASM's `u1(lambda)`.
    pub fn p(&mut self, lambda: f64, target: usize) -> &mut Self {
        self.ops.push(Op::Single {
            gate: gates::p(lambda),
            target,
            label: "P",
            params: vec![lambda],
        });
        self
    }

    /// General single-qubit gate U3(theta, phi, lambda) on `target`
    /// (the OpenQASM `u3`).
    pub fn u3(&mut self, theta: f64, phi: f64, lambda: f64, target: usize) -> &mut Self {
        self.ops.push(Op::Single {
            gate: gates::u3(theta, phi, lambda),
            target,
            label: "U3",
            params: vec![theta, phi, lambda],
        });
        self
    }

    /// Single-qubit gate U2(phi, lambda) = U3(pi/2, phi, lambda) on `target`
    /// (the OpenQASM `u2`).
    pub fn u2(&mut self, phi: f64, lambda: f64, target: usize) -> &mut Self {
        self.ops.push(Op::Single {
            gate: gates::u2(phi, lambda),
            target,
            label: "U2",
            params: vec![phi, lambda],
        });
        self
    }

    /// Rotation about the X axis by `theta` on `target`.
    pub fn rx(&mut self, theta: f64, target: usize) -> &mut Self {
        self.ops.push(Op::Single {
            gate: gates::rx(theta),
            target,
            label: "RX",
            params: vec![theta],
        });
        self
    }

    /// Rotation about the Y axis by `theta` on `target`.
    pub fn ry(&mut self, theta: f64, target: usize) -> &mut Self {
        self.ops.push(Op::Single {
            gate: gates::ry(theta),
            target,
            label: "RY",
            params: vec![theta],
        });
        self
    }

    /// Rotation about the Z axis by `theta` on `target`.
    pub fn rz(&mut self, theta: f64, target: usize) -> &mut Self {
        self.ops.push(Op::Single {
            gate: gates::rz(theta),
            target,
            label: "RZ",
            params: vec![theta],
        });
        self
    }

    /// Controlled-NOT: flips `target` when `control` is |1>.
    pub fn cnot(&mut self, control: usize, target: usize) -> &mut Self {
        self.ops.push(Op::Controlled {
            gate: gates::x(),
            control,
            target,
            label: "X",
            params: Vec::new(),
        });
        self
    }

    /// Controlled-U: apply the arbitrary 2x2 unitary `gate` to `target` only
    /// on basis states where `control` is |1>.
    pub fn cu(&mut self, gate: Gate, control: usize, target: usize) -> &mut Self {
        self.ops.push(Op::Controlled {
            gate,
            control,
            target,
            label: "U",
            params: Vec::new(),
        });
        self
    }

    /// Controlled-Rz: apply `rz(theta)` to `target` when `control` is |1>
    /// (the OpenQASM `crz`).
    pub fn crz(&mut self, theta: f64, control: usize, target: usize) -> &mut Self {
        self.ops.push(Op::Controlled {
            gate: gates::rz(theta),
            control,
            target,
            label: "CRZ",
            params: vec![theta],
        });
        self
    }

    /// Controlled phase: apply a phase e^{i·lambda} to the |11> component of
    /// `control` and `target` (the OpenQASM `cu1`). Symmetric in its two
    /// arguments.
    pub fn cp(&mut self, lambda: f64, control: usize, target: usize) -> &mut Self {
        self.ops.push(Op::Controlled {
            gate: gates::p(lambda),
            control,
            target,
            label: "CP",
            params: vec![lambda],
        });
        self
    }

    /// Controlled-U3: apply `u3(theta, phi, lambda)` to `target` when `control`
    /// is |1> (the OpenQASM `cu3`).
    pub fn cu3(
        &mut self,
        theta: f64,
        phi: f64,
        lambda: f64,
        control: usize,
        target: usize,
    ) -> &mut Self {
        self.ops.push(Op::Controlled {
            gate: gates::u3(theta, phi, lambda),
            control,
            target,
            label: "CU3",
            params: vec![theta, phi, lambda],
        });
        self
    }

    /// Controlled-Z: apply a phase of -1 to the |11> component of `control`
    /// and `target`. Symmetric in its two arguments.
    pub fn cz(&mut self, control: usize, target: usize) -> &mut Self {
        self.ops.push(Op::Controlled {
            gate: gates::z(),
            control,
            target,
            label: "Z",
            params: Vec::new(),
        });
        self
    }

    /// Controlled-Y: apply Y to `target` when `control` is |1> (OpenQASM `cy`).
    pub fn cy(&mut self, control: usize, target: usize) -> &mut Self {
        self.ops.push(Op::Controlled {
            gate: gates::y(),
            control,
            target,
            label: "Y",
            params: Vec::new(),
        });
        self
    }

    /// Controlled-Hadamard: apply H to `target` when `control` is |1>
    /// (OpenQASM `ch`).
    pub fn ch(&mut self, control: usize, target: usize) -> &mut Self {
        self.ops.push(Op::Controlled {
            gate: gates::h(),
            control,
            target,
            label: "H",
            params: Vec::new(),
        });
        self
    }

    /// SWAP: exchange the states of qubits `a` and `b`.
    pub fn swap(&mut self, a: usize, b: usize) -> &mut Self {
        self.ops.push(Op::Swap { a, b });
        self
    }

    /// Controlled-SWAP (Fredkin): exchange qubits `a` and `b` when `control`
    /// is |1> (OpenQASM `cswap`). Implemented via the standard identity
    /// `CSWAP(c,a,b) = CNOT(b,a)·CCX(c,a,b)·CNOT(b,a)`.
    pub fn cswap(&mut self, control: usize, a: usize, b: usize) -> &mut Self {
        self.cnot(b, a);
        self.toffoli(control, a, b);
        self.cnot(b, a);
        self
    }

    /// Multi-controlled-U: apply the arbitrary 2x2 unitary `gate` to `target`
    /// only on basis states where *every* qubit in `controls` is |1>.
    ///
    /// Zero controls is an unconditional gate, one control matches [`cu`], and
    /// two controls with X is a Toffoli.
    ///
    /// [`cu`]: Circuit::cu
    pub fn mcu(&mut self, gate: Gate, controls: &[usize], target: usize) -> &mut Self {
        self.ops.push(Op::MultiControlled {
            gate,
            controls: controls.to_vec(),
            target,
            label: "U",
        });
        self
    }

    /// Multi-controlled-X: flip `target` only when every qubit in `controls`
    /// is |1>. The generalization of [`cnot`] and [`toffoli`] to any number of
    /// controls.
    ///
    /// [`cnot`]: Circuit::cnot
    /// [`toffoli`]: Circuit::toffoli
    pub fn mcx(&mut self, controls: &[usize], target: usize) -> &mut Self {
        self.ops.push(Op::MultiControlled {
            gate: gates::x(),
            controls: controls.to_vec(),
            target,
            label: "X",
        });
        self
    }

    /// Toffoli (CCNOT): flip `target` only when both `control1` and
    /// `control2` are |1>.
    pub fn toffoli(&mut self, control1: usize, control2: usize, target: usize) -> &mut Self {
        self.mcx(&[control1, control2], target)
    }

    /// Run the circuit starting from |0...0> and return the final state.
    pub fn run(&self) -> State {
        let mut state = State::new(self.n_qubits);
        for op in &self.ops {
            match op {
                Op::Single { gate, target, .. } => state.apply_1q(gate, *target),
                Op::Controlled {
                    gate,
                    control,
                    target,
                    ..
                } => state.apply_controlled_1q(gate, *control, *target),
                Op::Swap { a, b } => state.swap_qubits(*a, *b),
                Op::MultiControlled {
                    gate,
                    controls,
                    target,
                    ..
                } => state.apply_multi_controlled_1q(gate, controls, *target),
            }
        }
        state
    }

    /// Run the circuit `shots` times and return a histogram of measured
    /// basis-state outcomes.
    ///
    /// The circuit is executed once, then each shot measures a fresh clone of
    /// the resulting state so the shots are independent. `seed` makes the
    /// whole sampling run deterministic and reproducible. Keys of the returned
    /// map are little-endian basis-state indices; values are counts.
    pub fn sample(&self, shots: usize, seed: u64) -> HashMap<usize, usize> {
        let final_state = self.run();
        let mut rng = Rng::new(seed);
        let mut histogram = HashMap::new();
        for _ in 0..shots {
            let outcome = final_state.clone().measure_all(&mut rng);
            *histogram.entry(outcome).or_insert(0) += 1;
        }
        histogram
    }

    /// Export the circuit as an OpenQASM 2.0 program string.
    ///
    /// Emits the standard header, a single `qreg q[n]`, and one gate per line.
    /// The output round-trips through [`qasm::parse`](crate::qasm::parse):
    /// re-importing it yields an equivalent circuit. Angles are written at full
    /// `f64` precision so the round trip is exact.
    ///
    /// Gates without a direct OpenQASM 2 equivalent are decomposed into ones
    /// that have one: an arbitrary controlled-U ([`cu`](Circuit::cu)) into a
    /// control phase (`u1`) plus a `cu3`, and a multi-controlled gate
    /// ([`mcx`](Circuit::mcx), [`mcu`](Circuit::mcu)) into Toffolis and
    /// single-qubit rotations. Only an unconditional [`mcu`](Circuit::mcu) —
    /// one with no controls at all — loses anything: its global phase, which
    /// OpenQASM 2 cannot express and no measurement can observe.
    ///
    /// Every built-in gate exports, so the error case is unreachable for
    /// circuits built through this API.
    ///
    /// # Example
    ///
    /// ```
    /// use qsimulator::Circuit;
    /// let mut c = Circuit::new(2);
    /// c.h(0).cnot(0, 1);
    /// let qasm = c.to_qasm().unwrap();
    /// assert!(qasm.starts_with("OPENQASM 2.0;"));
    /// assert!(qasm.contains("cx q[0],q[1];"));
    /// ```
    pub fn to_qasm(&self) -> Result<String, ExportError> {
        let mut out = String::from("OPENQASM 2.0;\ninclude \"qelib1.inc\";\n");
        out.push_str(&format!("qreg q[{}];\n", self.n_qubits));

        for op in &self.ops {
            match op {
                Op::Single {
                    target,
                    label,
                    params,
                    ..
                } => {
                    let name = match *label {
                        "I" => "id",
                        "H" => "h",
                        "X" => "x",
                        "Y" => "y",
                        "Z" => "z",
                        "S" => "s",
                        "T" => "t",
                        "SDG" => "sdg",
                        "SX" => "sx",
                        "SXDG" => "sxdg",
                        "TDG" => "tdg",
                        "P" => "u1",
                        "U2" => "u2",
                        "U3" => "u3",
                        "RX" => "rx",
                        "RY" => "ry",
                        "RZ" => "rz",
                        other => return Err(ExportError::SingleGate { label: other }),
                    };
                    out.push_str(&format!("{name}{} q[{target}];\n", format_params(params)));
                }
                Op::Controlled {
                    gate,
                    control,
                    target,
                    label,
                    params,
                } => {
                    match *label {
                        "X" => out.push_str(&format!("cx q[{control}],q[{target}];\n")),
                        "Y" => out.push_str(&format!("cy q[{control}],q[{target}];\n")),
                        "Z" => out.push_str(&format!("cz q[{control}],q[{target}];\n")),
                        "H" => out.push_str(&format!("ch q[{control}],q[{target}];\n")),
                        "CRZ" | "CP" | "CU3" => {
                            let name = match *label {
                                "CRZ" => "crz",
                                "CP" => "cu1",
                                _ => "cu3",
                            };
                            out.push_str(&format!(
                                "{name}{} q[{control}],q[{target}];\n",
                                format_params(params)
                            ));
                        }
                        // Any other controlled single-qubit unitary (e.g. `cu`
                        // with an arbitrary matrix): decompose into a phase on
                        // the control followed by a controlled-U3.
                        _ => emit_controlled_u(&mut out, gate, *control, *target),
                    }
                }
                Op::Swap { a, b } => {
                    out.push_str(&format!("swap q[{a}],q[{b}];\n"));
                }
                Op::MultiControlled {
                    gate,
                    controls,
                    target,
                    label,
                } => {
                    if *label == "X" {
                        emit_mcx(&mut out, self.n_qubits, controls, *target);
                    } else {
                        emit_mcu(&mut out, self.n_qubits, gate, controls, *target);
                    }
                }
            }
        }
        Ok(out)
    }

    /// Render the circuit as an ASCII diagram.
    ///
    /// One column per operation, time flowing left to right, with qubit `q0`
    /// on the top row. Controls are drawn as `*`, targets as their gate label
    /// (CNOT/Toffoli targets as `X`), SWAP endpoints as `x`, and `|` connects
    /// the qubits an operation spans. The diagram is presentational only —
    /// each operation occupies its own column, so it shows program order, not
    /// a parallel-scheduled timeline.
    ///
    /// # Example
    ///
    /// ```
    /// use qsimulator::Circuit;
    /// let mut c = Circuit::new(2);
    /// c.h(0).cnot(0, 1);
    /// assert_eq!(c.diagram(), "q0: -H-*-\nq1: ---X-");
    /// ```
    pub fn diagram(&self) -> String {
        let n = self.n_qubits;

        // For each op build one column: per qubit, the token to place, or
        // `None` for a plain wire.
        let mut columns: Vec<Vec<Option<&'static str>>> = Vec::with_capacity(self.ops.len());
        for op in &self.ops {
            let mut col: Vec<Option<&'static str>> = vec![None; n];
            match op {
                Op::Single { target, label, .. } => {
                    col[*target] = Some(label);
                }
                Op::Controlled {
                    control,
                    target,
                    label,
                    ..
                } => {
                    col[*control] = Some("*");
                    col[*target] = Some(label);
                    fill_connector(&mut col, &[*control, *target]);
                }
                Op::Swap { a, b } => {
                    col[*a] = Some("x");
                    col[*b] = Some("x");
                    fill_connector(&mut col, &[*a, *b]);
                }
                Op::MultiControlled {
                    controls,
                    target,
                    label,
                    ..
                } => {
                    for &c in controls {
                        col[c] = Some("*");
                    }
                    col[*target] = Some(label);
                    let mut involved = controls.clone();
                    involved.push(*target);
                    fill_connector(&mut col, &involved);
                }
            }
            columns.push(col);
        }

        // Cell width = widest token present (at least 1).
        let cell_w = columns
            .iter()
            .flatten()
            .flatten()
            .map(|s| s.len())
            .max()
            .unwrap_or(1)
            .max(1);

        // Width of the "qN:" label gutter, sized for the largest qubit index.
        let gutter = format!("q{}:", n.saturating_sub(1)).len();

        let mut lines = Vec::with_capacity(n);
        for r in 0..n {
            let mut line = format!("{:<width$} -", format!("q{r}:"), width = gutter);
            if columns.is_empty() {
                line.push_str(&"-".repeat(cell_w));
            }
            for col in &columns {
                line.push_str(&center(col[r].unwrap_or(""), cell_w));
                line.push('-');
            }
            lines.push(line);
        }
        lines.join("\n")
    }
}

/// Angles smaller than this are treated as zero, and their gate is left out of
/// the exported program.
const ANGLE_EPS: f64 = 1e-12;

/// Emit `name(angle) q[target];`, unless `angle` is negligible.
fn emit_rotation(out: &mut String, name: &str, angle: f64, target: usize) {
    if angle.abs() > ANGLE_EPS {
        out.push_str(&format!("{name}({angle}) q[{target}];\n"));
    }
}

/// Emit a controlled arbitrary 2x2 unitary as a phase on the control followed
/// by a `cu3`, using the Euler decomposition `gate == e^{iγ}·u3(θ, φ, λ)`.
fn emit_controlled_u(out: &mut String, gate: &Gate, control: usize, target: usize) {
    let (gamma, theta, phi, lambda) = gates::u3_decompose(gate);
    emit_rotation(out, "u1", gamma, control);
    out.push_str(&format!(
        "cu3({theta},{phi},{lambda}) q[{control}],q[{target}];\n"
    ));
}

/// The lowest-numbered qubit of the register that an operation on
/// `controls`/`target` leaves untouched, if any.
fn free_qubit(n_qubits: usize, controls: &[usize], target: usize) -> Option<usize> {
    (0..n_qubits).find(|q| *q != target && !controls.contains(q))
}

/// Emit a multi-controlled X (`target ^= AND(controls)`) as OpenQASM 2 gates.
///
/// Up to two controls this is `x`/`cx`/`ccx` directly. Beyond that it is a
/// Barenco-style decomposition into Toffolis, in one of two ways:
///
/// - **With a spare qubit** — any register qubit the operation does not touch
///   can be *borrowed* even though its state is unknown ("dirty"). The
///   controls are split into two halves; the borrowed qubit is toggled by the
///   first half and used as an extra control for the second, and the whole
///   pair is run twice so the qubit's unknown initial value cancels out of the
///   target and the qubit itself is left as it was found. Each half is a
///   strictly smaller multi-controlled X, so this recurses to Toffolis in
///   `O(controls²)` gates.
/// - **With no spare qubit** (the controls plus the target are the whole
///   register) — the square-root recursion, with `V·V = X`:
///   `C^m(X) = C(V)[c_m,t] · C^{m-1}(X)[c_<m→c_m] · C(V†)[c_m,t] ·
///   C^{m-1}(X)[c_<m→c_m] · C^{m-1}(V)[c_<m→t]`. Every inner operation leaves
///   a qubit untouched, so the recursion lands in the borrowed-qubit case.
fn emit_mcx(out: &mut String, n_qubits: usize, controls: &[usize], target: usize) {
    match controls {
        [] => out.push_str(&format!("x q[{target}];\n")),
        [c] => out.push_str(&format!("cx q[{c}],q[{target}];\n")),
        [a, b] => out.push_str(&format!("ccx q[{a}],q[{b}],q[{target}];\n")),
        _ => match free_qubit(n_qubits, controls, target) {
            Some(borrowed) => {
                // Both halves are shorter than `controls` (the split point is
                // between 2 and controls.len() - 1), so this terminates.
                let (lo, hi) = controls.split_at(controls.len().div_ceil(2));
                let mut hi_borrowed = hi.to_vec();
                hi_borrowed.push(borrowed);
                for _ in 0..2 {
                    emit_mcx(out, n_qubits, lo, borrowed);
                    emit_mcx(out, n_qubits, &hi_borrowed, target);
                }
            }
            None => {
                let (&last, rest) = controls.split_last().expect("at least three controls");
                emit_controlled_u(out, &gates::sx(), last, target);
                emit_mcx(out, n_qubits, rest, last);
                emit_controlled_u(out, &gates::sxdg(), last, target);
                emit_mcx(out, n_qubits, rest, last);
                emit_mcu(out, n_qubits, &gates::sx(), rest, target);
            }
        },
    }
}

/// Emit a multi-controlled arbitrary 2x2 unitary as OpenQASM 2 gates.
///
/// Writes `gate` as `e^{iγ''}·A·X·B·X·C` with `A·B·C = I`, where `A`, `B` and
/// `C` are products of `rz`/`ry` built from the Euler angles of
/// `gates::u3_decompose`. Conditioning only the two `X`s (on the full control
/// set, via [`emit_mcx`]) and the phase `γ''` (via [`emit_mcphase`]) then gives
/// the controlled gate: when the controls are not all set, `A·B·C` collapses to
/// the identity and no phase is applied.
///
/// A diagonal `gate` — a multi-controlled Z, S, T or phase — skips the `X`s
/// entirely and exports as two phase terms.
///
/// With no controls the gate is unconditional, so its global phase is
/// unobservable and is dropped (OpenQASM 2 cannot express one).
fn emit_mcu(out: &mut String, n_qubits: usize, gate: &Gate, controls: &[usize], target: usize) {
    let (gamma, theta, phi, lambda) = gates::u3_decompose(gate);

    if controls.is_empty() {
        out.push_str(&format!("u3({theta},{phi},{lambda}) q[{target}];\n"));
        return;
    }

    // Diagonal gate: diag(e^{iγ}, e^{i(γ+λ)}). The first factor is a phase on
    // the controls alone, the second a phase on the controls and the target.
    // Cheaper than the general path at every width, so it comes first.
    if theta.abs() < ANGLE_EPS {
        emit_mcphase(out, n_qubits, gamma, controls);
        let mut controls_and_target = controls.to_vec();
        controls_and_target.push(target);
        emit_mcphase(out, n_qubits, lambda, &controls_and_target);
        return;
    }

    if controls.len() == 1 {
        emit_controlled_u(out, gate, controls[0], target);
        return;
    }

    // C = Rz((λ-φ)/2), B = Ry(-θ/2)·Rz(-(λ+φ)/2), A = Rz(φ)·Ry(θ/2). Then
    // A·B·C = I and A·X·B·X·C = Rz(φ)·Ry(θ)·Rz(λ), which is u3(θ,φ,λ) up to
    // the phase e^{i(φ+λ)/2} — folded into γ'' below. Gates are emitted in
    // time order, so each matrix product is written back to front.
    emit_rotation(out, "rz", (lambda - phi) / 2.0, target);
    emit_mcx(out, n_qubits, controls, target);
    emit_rotation(out, "rz", -(lambda + phi) / 2.0, target);
    emit_rotation(out, "ry", -theta / 2.0, target);
    emit_mcx(out, n_qubits, controls, target);
    emit_rotation(out, "ry", theta / 2.0, target);
    emit_rotation(out, "rz", phi, target);
    emit_mcphase(out, n_qubits, gamma + (phi + lambda) / 2.0, controls);
}

/// Emit a phase of `e^{i·lambda}` applied only when every qubit in `qubits` is
/// |1>, as OpenQASM 2 gates.
///
/// One qubit is `u1` and two are `cu1`. Beyond that it uses the same identity
/// that decomposes `cu1` — half the phase on the leading qubits, then the
/// remaining qubit conjugated by a multi-controlled X — which peels off one
/// qubit per step. Since a phase is diagonal, the emitted gates commute and the
/// order among them does not matter.
fn emit_mcphase(out: &mut String, n_qubits: usize, lambda: f64, qubits: &[usize]) {
    if lambda.abs() < ANGLE_EPS {
        return;
    }
    match qubits {
        // A global phase, which OpenQASM 2 cannot express and no measurement
        // can see.
        [] => {}
        [q] => out.push_str(&format!("u1({lambda}) q[{q}];\n")),
        [a, b] => out.push_str(&format!("cu1({lambda}) q[{a}],q[{b}];\n")),
        _ => {
            let (&last, rest) = qubits.split_last().expect("at least three qubits");
            let half = lambda / 2.0;
            emit_mcphase(out, n_qubits, half, rest);
            emit_mcx(out, n_qubits, rest, last);
            emit_rotation(out, "u1", -half, last);
            emit_mcx(out, n_qubits, rest, last);
            emit_rotation(out, "u1", half, last);
        }
    }
}

/// Format gate angle parameters as an OpenQASM parameter list: `""` for none,
/// `"(θ)"` for one, `"(θ,φ,λ)"` for several. Angles use full `f64` precision
/// so exported programs round-trip exactly.
fn format_params(params: &[f64]) -> String {
    if params.is_empty() {
        return String::new();
    }
    let joined = params
        .iter()
        .map(|p| p.to_string())
        .collect::<Vec<_>>()
        .join(",");
    format!("({joined})")
}

/// Mark the rows strictly between the outermost involved qubits (that are not
/// themselves involved) with a `|` connector.
fn fill_connector(col: &mut [Option<&'static str>], involved: &[usize]) {
    let min = *involved.iter().min().unwrap();
    let max = *involved.iter().max().unwrap();
    for (r, cell) in col.iter_mut().enumerate() {
        if r > min && r < max && cell.is_none() {
            *cell = Some("|");
        }
    }
}

/// Center `token` in a field of `width`, padded with `-` (the wire character).
fn center(token: &str, width: usize) -> String {
    if token.len() >= width {
        return token.to_string();
    }
    let pad = width - token.len();
    let left = pad / 2;
    let right = pad - left;
    format!("{}{}{}", "-".repeat(left), token, "-".repeat(right))
}

impl fmt::Display for Circuit {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.diagram())
    }
}
