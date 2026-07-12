//! Circuit builder and executor.

use crate::gates;
use crate::state::State;
use num_complex::Complex64;

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

    /// Controlled-NOT: flips `target` when `control` is |1>.
    pub fn cnot(&mut self, control: usize, target: usize) -> &mut Self {
        self.ops.push(Op::Controlled {
            gate: gates::x(),
            control,
            target,
        });
        self
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
            }
        }
        state
    }
}
