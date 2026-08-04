//! Circuit builder and executor.

use crate::gates;
use crate::rng::Rng;
use crate::state::State;
use num_complex::Complex64;
use std::collections::HashMap;
use std::fmt;

type Gate = [[Complex64; 2]; 2];

/// A single instruction in a circuit.
///
/// Each gate-bearing variant carries a short `label` (e.g. `"H"`, `"RX"`)
/// used for diagram rendering and QASM export; it never affects execution.
/// `Single` additionally records a `param` (the angle) for rotation gates so
/// that export is lossless; it is `None` for non-parametric gates.
#[derive(Debug, Clone)]
enum Op {
    Single {
        gate: Gate,
        target: usize,
        label: &'static str,
        param: Option<f64>,
    },
    Controlled {
        gate: Gate,
        control: usize,
        target: usize,
        label: &'static str,
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

    pub fn h(&mut self, target: usize) -> &mut Self {
        self.ops.push(Op::Single {
            gate: gates::h(),
            target,
            label: "H",
            param: None,
        });
        self
    }

    pub fn x(&mut self, target: usize) -> &mut Self {
        self.ops.push(Op::Single {
            gate: gates::x(),
            target,
            label: "X",
            param: None,
        });
        self
    }

    pub fn z(&mut self, target: usize) -> &mut Self {
        self.ops.push(Op::Single {
            gate: gates::z(),
            target,
            label: "Z",
            param: None,
        });
        self
    }

    pub fn y(&mut self, target: usize) -> &mut Self {
        self.ops.push(Op::Single {
            gate: gates::y(),
            target,
            label: "Y",
            param: None,
        });
        self
    }

    /// Phase gate S = diag(1, i).
    pub fn s(&mut self, target: usize) -> &mut Self {
        self.ops.push(Op::Single {
            gate: gates::s(),
            target,
            label: "S",
            param: None,
        });
        self
    }

    /// T gate = diag(1, e^{i pi/4}).
    pub fn t(&mut self, target: usize) -> &mut Self {
        self.ops.push(Op::Single {
            gate: gates::t(),
            target,
            label: "T",
            param: None,
        });
        self
    }

    /// Rotation about the X axis by `theta` on `target`.
    pub fn rx(&mut self, theta: f64, target: usize) -> &mut Self {
        self.ops.push(Op::Single {
            gate: gates::rx(theta),
            target,
            label: "RX",
            param: Some(theta),
        });
        self
    }

    /// Rotation about the Y axis by `theta` on `target`.
    pub fn ry(&mut self, theta: f64, target: usize) -> &mut Self {
        self.ops.push(Op::Single {
            gate: gates::ry(theta),
            target,
            label: "RY",
            param: Some(theta),
        });
        self
    }

    /// Rotation about the Z axis by `theta` on `target`.
    pub fn rz(&mut self, theta: f64, target: usize) -> &mut Self {
        self.ops.push(Op::Single {
            gate: gates::rz(theta),
            target,
            label: "RZ",
            param: Some(theta),
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
        });
        self
    }

    /// SWAP: exchange the states of qubits `a` and `b`.
    pub fn swap(&mut self, a: usize, b: usize) -> &mut Self {
        self.ops.push(Op::Swap { a, b });
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
    /// re-importing it yields an equivalent circuit. Rotation angles are
    /// written at full `f64` precision so the round trip is exact.
    ///
    /// Returns an error for gates the OpenQASM subset cannot express directly:
    /// an arbitrary controlled-U ([`cu`](Circuit::cu)/[`mcu`](Circuit::mcu)),
    /// or a multi-controlled-X with a control count other than two (only
    /// [`toffoli`](Circuit::toffoli)/`ccx` is representable).
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
    pub fn to_qasm(&self) -> Result<String, String> {
        let mut out = String::from("OPENQASM 2.0;\ninclude \"qelib1.inc\";\n");
        out.push_str(&format!("qreg q[{}];\n", self.n_qubits));

        for op in &self.ops {
            match op {
                Op::Single {
                    target,
                    label,
                    param,
                    ..
                } => {
                    let name = match *label {
                        "H" => "h",
                        "X" => "x",
                        "Y" => "y",
                        "Z" => "z",
                        "S" => "s",
                        "T" => "t",
                        "RX" => "rx",
                        "RY" => "ry",
                        "RZ" => "rz",
                        other => return Err(format!("cannot export single-qubit gate `{other}`")),
                    };
                    match param {
                        Some(theta) => out.push_str(&format!("{name}({theta}) q[{target}];\n")),
                        None => out.push_str(&format!("{name} q[{target}];\n")),
                    }
                }
                Op::Controlled {
                    control,
                    target,
                    label,
                    ..
                } => {
                    let name = match *label {
                        "X" => "cx",
                        "Z" => "cz",
                        _ => {
                            return Err(
                                "cannot export an arbitrary controlled-U (cu) to OpenQASM".into()
                            )
                        }
                    };
                    out.push_str(&format!("{name} q[{control}],q[{target}];\n"));
                }
                Op::Swap { a, b } => {
                    out.push_str(&format!("swap q[{a}],q[{b}];\n"));
                }
                Op::MultiControlled {
                    controls,
                    target,
                    label,
                    ..
                } => {
                    if *label != "X" {
                        return Err("cannot export a multi-controlled-U (mcu) to OpenQASM".into());
                    }
                    if controls.len() != 2 {
                        return Err(format!(
                            "cannot export a multi-controlled-X with {} controls (only 2 = ccx)",
                            controls.len()
                        ));
                    }
                    out.push_str(&format!(
                        "ccx q[{}],q[{}],q[{target}];\n",
                        controls[0], controls[1]
                    ));
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
