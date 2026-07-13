//! Circuit builder and executor.

use crate::gates;
use crate::rng::Rng;
use crate::state::State;
use num_complex::Complex64;
use std::collections::HashMap;
use std::fmt;

type Gate = [[Complex64; 2]; 2];

/// A single instruction in a circuit.
enum Op {
    Single {
        gate: Gate,
        target: usize,
        name: &'static str,
    },
    Controlled {
        gate: Gate,
        control: usize,
        target: usize,
        name: &'static str,
    },
    Swap {
        a: usize,
        b: usize,
    },
    MultiControlled {
        gate: Gate,
        controls: Vec<usize>,
        target: usize,
        name: &'static str,
    },
}

/// A quantum circuit: an ordered list of gate operations on `n_qubits`.
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
            name: "H",
        });
        self
    }

    pub fn x(&mut self, target: usize) -> &mut Self {
        self.ops.push(Op::Single {
            gate: gates::x(),
            target,
            name: "X",
        });
        self
    }

    pub fn z(&mut self, target: usize) -> &mut Self {
        self.ops.push(Op::Single {
            gate: gates::z(),
            target,
            name: "Z",
        });
        self
    }

    /// Rotation about the X axis by `theta` on `target`.
    pub fn rx(&mut self, theta: f64, target: usize) -> &mut Self {
        self.ops.push(Op::Single {
            gate: gates::rx(theta),
            target,
            name: "Rx",
        });
        self
    }

    /// Rotation about the Y axis by `theta` on `target`.
    pub fn ry(&mut self, theta: f64, target: usize) -> &mut Self {
        self.ops.push(Op::Single {
            gate: gates::ry(theta),
            target,
            name: "Ry",
        });
        self
    }

    /// Rotation about the Z axis by `theta` on `target`.
    pub fn rz(&mut self, theta: f64, target: usize) -> &mut Self {
        self.ops.push(Op::Single {
            gate: gates::rz(theta),
            target,
            name: "Rz",
        });
        self
    }

    pub fn y(&mut self, target: usize) -> &mut Self {
        self.ops.push(Op::Single {
            gate: gates::y(),
            target,
            name: "Y",
        });
        self
    }

    pub fn s(&mut self, target: usize) -> &mut Self {
        self.ops.push(Op::Single {
            gate: gates::s(),
            target,
            name: "S",
        });
        self
    }

    pub fn t(&mut self, target: usize) -> &mut Self {
        self.ops.push(Op::Single {
            gate: gates::t(),
            target,
            name: "T",
        });
        self
    }

    /// Controlled-Z: applies Z to `target` when `control` is |1>.
    pub fn cz(&mut self, control: usize, target: usize) -> &mut Self {
        self.ops.push(Op::Controlled {
            gate: gates::z(),
            control,
            target,
            name: "Z",
        });
        self
    }

    /// Controlled-U: applies arbitrary single-qubit `gate` to `target`
    /// when `control` is |1>.
    pub fn cu(&mut self, gate: Gate, control: usize, target: usize) -> &mut Self {
        self.ops.push(Op::Controlled {
            gate,
            control,
            target,
            name: "U",
        });
        self
    }

    /// Controlled-NOT: flips `target` when `control` is |1>.
    pub fn cnot(&mut self, control: usize, target: usize) -> &mut Self {
        self.ops.push(Op::Controlled {
            gate: gates::x(),
            control,
            target,
            name: "X",
        });
        self
    }

    /// SWAP: exchange the states of qubits `a` and `b`.
    pub fn swap(&mut self, a: usize, b: usize) -> &mut Self {
        self.ops.push(Op::Swap { a, b });
        self
    }

    /// Toffoli (CCNOT): flip `target` only when both `control1` and
    /// `control2` are |1>.
    pub fn toffoli(&mut self, control1: usize, control2: usize, target: usize) -> &mut Self {
        self.ops.push(Op::MultiControlled {
            gate: gates::x(),
            controls: vec![control1, control2],
            target,
            name: "X",
        });
        self
    }

    /// Multi-controlled X: flip `target` when all `controls` are |1>.
    pub fn mcx(&mut self, controls: &[usize], target: usize) -> &mut Self {
        self.ops.push(Op::MultiControlled {
            gate: gates::x(),
            controls: controls.to_vec(),
            target,
            name: "X",
        });
        self
    }

    /// Multi-controlled arbitrary gate: apply `gate` to `target` when all
    /// `controls` are |1>.
    pub fn mcu(&mut self, gate: Gate, controls: &[usize], target: usize) -> &mut Self {
        self.ops.push(Op::MultiControlled {
            gate,
            controls: controls.to_vec(),
            target,
            name: "U",
        });
        self
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
}

impl fmt::Display for Circuit {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let n = self.n_qubits;
        if n == 0 || self.ops.is_empty() {
            return write!(f, "(empty circuit)");
        }

        // Build a column per op; each column is n_qubits cells.
        // Cell content: gate label, control dot "●", swap "×", or wire "─".
        let mut columns: Vec<Vec<String>> = Vec::new();

        for op in &self.ops {
            let mut col = vec![String::new(); n];
            match op {
                Op::Single { target, name, .. } => {
                    col[*target] = format!("[{name}]");
                }
                Op::Controlled {
                    control,
                    target,
                    name,
                    ..
                } => {
                    col[*control] = "●".to_string();
                    col[*target] = format!("[{name}]");
                    let lo = (*control).min(*target);
                    let hi = (*control).max(*target);
                    for cell in col.iter_mut().take(hi).skip(lo + 1) {
                        if cell.is_empty() {
                            *cell = "│".to_string();
                        }
                    }
                }
                Op::Swap { a, b } => {
                    col[*a] = "×".to_string();
                    col[*b] = "×".to_string();
                    let lo = (*a).min(*b);
                    let hi = (*a).max(*b);
                    for cell in col.iter_mut().take(hi).skip(lo + 1) {
                        if cell.is_empty() {
                            *cell = "│".to_string();
                        }
                    }
                }
                Op::MultiControlled {
                    controls,
                    target,
                    name,
                    ..
                } => {
                    for &c in controls {
                        col[c] = "●".to_string();
                    }
                    col[*target] = format!("[{name}]");
                    let all: Vec<usize> = controls
                        .iter()
                        .copied()
                        .chain(std::iter::once(*target))
                        .collect();
                    let lo = *all.iter().min().unwrap();
                    let hi = *all.iter().max().unwrap();
                    for cell in col.iter_mut().take(hi).skip(lo + 1) {
                        if cell.is_empty() {
                            *cell = "│".to_string();
                        }
                    }
                }
            }
            columns.push(col);
        }

        // Determine the display width of each column.
        let col_widths: Vec<usize> = columns
            .iter()
            .map(|col| col.iter().map(|s| s.len()).max().unwrap_or(1).max(1))
            .collect();

        for q in 0..n {
            write!(f, "q{q}: ")?;
            for (ci, col) in columns.iter().enumerate() {
                let w = col_widths[ci];
                let cell = &col[q];
                if cell.is_empty() {
                    // Plain wire
                    for _ in 0..w {
                        write!(f, "─")?;
                    }
                } else {
                    let pad = w.saturating_sub(cell.len());
                    let left = pad / 2;
                    let right = pad - left;
                    for _ in 0..left {
                        write!(f, "─")?;
                    }
                    write!(f, "{cell}")?;
                    for _ in 0..right {
                        write!(f, "─")?;
                    }
                }
                if ci + 1 < columns.len() {
                    write!(f, "─")?;
                }
            }
            writeln!(f)?;
        }
        Ok(())
    }
}
