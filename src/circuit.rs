//! Circuit builder and executor.

use crate::gates;
use crate::rng::Rng;
use crate::state::State;
use num_complex::Complex64;
use std::collections::HashMap;

type Gate = [[Complex64; 2]; 2];

/// A single instruction in a circuit.
enum Op {
    Single {
        gate: Gate,
        target: usize,
    },
    Controlled {
        gate: Gate,
        control: usize,
        target: usize,
    },
    Swap {
        a: usize,
        b: usize,
    },
    MultiControlled {
        gate: Gate,
        controls: Vec<usize>,
        target: usize,
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
        });
        self
    }

    pub fn x(&mut self, target: usize) -> &mut Self {
        self.ops.push(Op::Single {
            gate: gates::x(),
            target,
        });
        self
    }

    pub fn z(&mut self, target: usize) -> &mut Self {
        self.ops.push(Op::Single {
            gate: gates::z(),
            target,
        });
        self
    }

    pub fn y(&mut self, target: usize) -> &mut Self {
        self.ops.push(Op::Single {
            gate: gates::y(),
            target,
        });
        self
    }

    /// Phase gate S = diag(1, i).
    pub fn s(&mut self, target: usize) -> &mut Self {
        self.ops.push(Op::Single {
            gate: gates::s(),
            target,
        });
        self
    }

    /// T gate = diag(1, e^{i pi/4}).
    pub fn t(&mut self, target: usize) -> &mut Self {
        self.ops.push(Op::Single {
            gate: gates::t(),
            target,
        });
        self
    }

    /// Rotation about the X axis by `theta` on `target`.
    pub fn rx(&mut self, theta: f64, target: usize) -> &mut Self {
        self.ops.push(Op::Single {
            gate: gates::rx(theta),
            target,
        });
        self
    }

    /// Rotation about the Y axis by `theta` on `target`.
    pub fn ry(&mut self, theta: f64, target: usize) -> &mut Self {
        self.ops.push(Op::Single {
            gate: gates::ry(theta),
            target,
        });
        self
    }

    /// Rotation about the Z axis by `theta` on `target`.
    pub fn rz(&mut self, theta: f64, target: usize) -> &mut Self {
        self.ops.push(Op::Single {
            gate: gates::rz(theta),
            target,
        });
        self
    }

    /// Controlled-NOT: flips `target` when `control` is |1>.
    pub fn cnot(&mut self, control: usize, target: usize) -> &mut Self {
        self.ops.push(Op::Controlled {
            gate: gates::x(),
            control,
            target,
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
        });
        self
    }

    /// Controlled-Z: apply a phase of -1 to the |11> component of `control`
    /// and `target`. Symmetric in its two arguments.
    pub fn cz(&mut self, control: usize, target: usize) -> &mut Self {
        self.cu(gates::z(), control, target)
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
        self.mcu(gates::x(), controls, target)
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
                Op::Single { gate, target } => state.apply_1q(gate, *target),
                Op::Controlled {
                    gate,
                    control,
                    target,
                } => state.apply_controlled_1q(gate, *control, *target),
                Op::Swap { a, b } => state.swap_qubits(*a, *b),
                Op::MultiControlled {
                    gate,
                    controls,
                    target,
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
