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
#[non_exhaustive]
pub enum ExportError {
    /// An arbitrary single-qubit gate with no OpenQASM name. Built-in gates
    /// never produce this; it can only arise from a manually built `Op`.
    SingleGate {
        /// The gate's diagram label.
        label: &'static str,
    },
    /// A noise channel. OpenQASM 2 describes unitary circuits and measurement;
    /// it has no syntax for a quantum channel, so a noisy circuit cannot be
    /// written out without silently dropping the noise.
    Noise,
}

impl fmt::Display for ExportError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ExportError::SingleGate { label } => {
                write!(f, "cannot export single-qubit gate `{label}` to OpenQASM 2")
            }
            ExportError::Noise => f.write_str(
                "cannot export a noise channel to OpenQASM 2, which has no syntax for one",
            ),
        }
    }
}

impl std::error::Error for ExportError {}

/// Reason a circuit could not be run as a density matrix by
/// [`Circuit::run_density`].
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum DensityError {
    /// The register is too large for a `2^n x 2^n` matrix. The state-vector
    /// backend reaches roughly twice as many qubits.
    TooManyQubits {
        /// The circuit's qubit count.
        qubits: usize,
        /// The largest a density matrix will allocate.
        max: usize,
    },
    /// Too many distinct classical outcomes to track. Feed-forward is handled
    /// by carrying one density matrix per reachable classical register value,
    /// so a circuit measuring many bits before branching on them multiplies
    /// what has to be held at once.
    TooManyBranches {
        /// The cap that was exceeded.
        max: usize,
    },
}

impl fmt::Display for DensityError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            DensityError::TooManyQubits { qubits, max } => write!(
                f,
                "a density matrix for {qubits} qubits is 4^{qubits} entries, over the \
                 maximum of {max}"
            ),
            DensityError::TooManyBranches { max } => write!(
                f,
                "classical feed-forward needs one density matrix per reachable register \
                 value, and this circuit exceeds the limit of {max}; sample it instead"
            ),
        }
    }
}

impl std::error::Error for DensityError {}

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
    /// A mid-circuit measurement: collapses `qubit` onto the sampled outcome.
    /// The classical bit it would be written to is not modelled, since nothing
    /// in the supported subset can read one back.
    Measure {
        qubit: usize,
        clbit: usize,
    },
    /// Collapse `qubit` and force it to |0>, leaving the rest of the register
    /// on whichever branch the collapse took.
    Reset {
        qubit: usize,
    },
    MultiControlled {
        gate: Gate,
        controls: Vec<usize>,
        target: usize,
        label: &'static str,
    },
    /// A noise channel on `qubit`: one Kraus operator is sampled per shot.
    Kraus {
        ops: Vec<[[Complex64; 2]; 2]>,
        qubit: usize,
        label: &'static str,
    },
    /// Gates that run only when the whole classical register equals `value`.
    /// The block holds no measurement or reset, so `value` cannot change
    /// while it executes.
    Conditional {
        value: u64,
        ops: Vec<Op>,
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

    /// How many qubits the circuit acts on.
    pub fn n_qubits(&self) -> usize {
        self.n_qubits
    }

    /// Append another circuit's operations to this one.
    ///
    /// Both must be over the same number of qubits. Used by the OpenQASM
    /// importer to build a guarded statement in isolation — so its operands
    /// resolve and its errors surface — before folding it into a conditional.
    pub(crate) fn extend_from(&mut self, other: &Circuit) {
        debug_assert_eq!(self.n_qubits, other.n_qubits);
        self.ops.extend(other.ops.iter().cloned());
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

    /// Apply a noise channel to `qubit`, given as Kraus operators.
    ///
    /// The channel is simulated by sampling one operator per shot rather than
    /// by carrying a density matrix, so a single [`run`](Circuit::run) is one
    /// trajectory and averages come from [`sample`](Circuit::sample). See the
    /// [`noise`](crate::noise) module for the standard channels and for what
    /// that costs in shots.
    ///
    /// # Panics
    ///
    /// If `ops` is empty or not trace preserving, since a channel that is not
    /// would quietly change the total probability.
    pub fn channel(&mut self, ops: Vec<[[Complex64; 2]; 2]>, qubit: usize) -> &mut Self {
        self.named_channel(ops, qubit, "N")
    }

    fn named_channel(
        &mut self,
        ops: Vec<[[Complex64; 2]; 2]>,
        qubit: usize,
        label: &'static str,
    ) -> &mut Self {
        assert!(
            !ops.is_empty(),
            "a channel needs at least one Kraus operator"
        );
        assert!(
            crate::noise::is_trace_preserving(&ops),
            "channel is not trace preserving (sum K^dagger K must be the identity)"
        );
        self.ops.push(Op::Kraus { ops, qubit, label });
        self
    }

    /// Depolarizing noise on `qubit` with probability `p`.
    pub fn depolarizing(&mut self, p: f64, qubit: usize) -> &mut Self {
        self.named_channel(crate::noise::depolarizing(p), qubit, "DEP")
    }

    /// Bit-flip noise on `qubit`: X with probability `p`.
    pub fn bit_flip(&mut self, p: f64, qubit: usize) -> &mut Self {
        self.named_channel(crate::noise::bit_flip(p), qubit, "BF")
    }

    /// Phase-flip noise on `qubit`: Z with probability `p`.
    pub fn phase_flip(&mut self, p: f64, qubit: usize) -> &mut Self {
        self.named_channel(crate::noise::phase_flip(p), qubit, "PF")
    }

    /// Amplitude damping on `qubit` with probability `gamma` (T1 relaxation).
    pub fn amplitude_damping(&mut self, gamma: f64, qubit: usize) -> &mut Self {
        self.named_channel(crate::noise::amplitude_damping(gamma), qubit, "AD")
    }

    /// Phase damping on `qubit` with probability `gamma`: coherence decays
    /// while populations stay put.
    pub fn phase_damping(&mut self, gamma: f64, qubit: usize) -> &mut Self {
        self.named_channel(crate::noise::phase_damping(gamma), qubit, "PD")
    }

    /// Reset `qubit` to |0>, whatever state it was in.
    ///
    /// This is a collapse followed by a flip if the outcome was |1>, so it is
    /// not a unitary: it destroys any superposition on `qubit`, and where that
    /// qubit was entangled it leaves the rest of the register on whichever
    /// branch the collapse took. Unlike a trailing
    /// [`measure`](Circuit::measure), a reset always applies — it changes the
    /// state that gets sampled even as the last operation.
    pub fn reset(&mut self, qubit: usize) -> &mut Self {
        self.ops.push(Op::Reset { qubit });
        self
    }

    /// Measure `qubit` in the computational basis, collapsing the register onto
    /// the observed outcome.
    ///
    /// This makes the circuit *stochastic*: everything after it depends on
    /// which outcome came up, so [`run`](Circuit::run) is no longer a pure
    /// function of the circuit alone and [`run_seeded`](Circuit::run_seeded)
    /// is how you choose the stream. The outcome itself is discarded — the
    /// supported subset has no classical control to read it back — but the
    /// collapse is what makes a later gate see a definite state rather than a
    /// superposition.
    pub fn measure(&mut self, qubit: usize) -> &mut Self {
        self.measure_into(qubit, qubit)
    }

    /// Measure `qubit`, recording the outcome in classical bit `clbit`.
    ///
    /// The circuit has one classical register, as wide as the quantum one, and
    /// [`measure`](Circuit::measure) writes qubit `i` to bit `i`. Use this when
    /// the destination differs — an imported program may compact several
    /// measurements into low bits.
    pub fn measure_into(&mut self, qubit: usize, clbit: usize) -> &mut Self {
        self.ops.push(Op::Measure { qubit, clbit });
        self
    }

    /// Run `build`'s gates only when the classical register equals `value`.
    ///
    /// This is OpenQASM's `if (c == value)`. The register is the one written by
    /// [`measure`](Circuit::measure), compared as a whole: bit `i` is the last
    /// outcome recorded for classical bit `i`, and unmeasured bits read 0.
    ///
    /// # Panics
    ///
    /// The block may only contain gates. A measurement or reset inside it would
    /// change the value being tested part-way through, so guarding each
    /// statement — which is all OpenQASM's single-statement `if` can express —
    /// would stop matching the block as a whole; and a nested conditional has
    /// no OpenQASM form at all. Rather than drop such an operation and quietly
    /// simulate something else, this panics.
    pub fn if_classical_eq(&mut self, value: u64, build: impl FnOnce(&mut Circuit)) -> &mut Self {
        let mut inner = Circuit::new(self.n_qubits);
        build(&mut inner);
        for op in &inner.ops {
            let kind = match op {
                Op::Measure { .. } => "a measurement",
                Op::Reset { .. } => "a reset",
                Op::Conditional { .. } => "a nested conditional",
                _ => continue,
            };
            panic!("a conditional block may only contain gates, but it contains {kind}");
        }
        if !inner.ops.is_empty() {
            self.ops.push(Op::Conditional {
                value,
                ops: inner.ops,
            });
        }
        self
    }

    /// Whether the circuit measures at all.
    fn has_measurement(&self) -> bool {
        self.ops.iter().any(|op| matches!(op, Op::Measure { .. }))
    }

    /// Whether the circuit reads the classical register.
    fn has_conditional(&self) -> bool {
        self.ops
            .iter()
            .any(|op| matches!(op, Op::Conditional { .. }))
    }

    /// Whether the circuit contains a reset or a noise channel. Both always
    /// branch — even as the final operation, since both change the state that
    /// is sampled.
    fn has_reset(&self) -> bool {
        self.ops
            .iter()
            .any(|op| matches!(op, Op::Reset { .. } | Op::Kraus { .. }))
    }

    /// How many of the final operations are measurements — the circuit's
    /// readout, as opposed to a measurement something else depends on.
    fn trailing_measurements(&self) -> usize {
        self.ops
            .iter()
            .rev()
            .take_while(|op| matches!(op, Op::Measure { .. }))
            .count()
    }

    /// Whether a collapse can affect anything that follows it.
    ///
    /// Only a measurement with an operation after it can: it changes what that
    /// operation sees. A measurement at the very end cannot, and it is that
    /// case alone that forces [`sample`](Circuit::sample) to re-run the circuit
    /// per shot.
    fn measures_mid_circuit(&self) -> bool {
        self.has_reset()
            || self.ops[..self.ops.len() - self.trailing_measurements()]
                .iter()
                .any(|op| matches!(op, Op::Measure { .. }))
    }

    /// Run the circuit starting from |0...0> and return the state it prepares.
    ///
    /// Measurements at the *end* of the circuit are readout and are not
    /// applied: collapsing them would discard the prepared state and report one
    /// arbitrary branch, and [`sample`](Circuit::sample) draws the same
    /// distribution either way. This keeps the result deterministic for the
    /// usual written-out program, which ends in a measurement.
    ///
    /// A measurement with gates *after* it does collapse, since those gates
    /// must see a definite state. That makes the result depend on the collapse
    /// outcomes, and this uses seed 0 — see
    /// [`run_seeded`](Circuit::run_seeded) to choose the stream.
    pub fn run(&self) -> State {
        self.run_seeded(0)
    }

    /// Run the circuit with an explicit seed for any mid-circuit measurement.
    ///
    /// Without a mid-circuit [`measure`](Circuit::measure) the seed changes
    /// nothing and this is exactly [`run`](Circuit::run).
    pub fn run_seeded(&self, seed: u64) -> State {
        let mut rng = Rng::new(seed);
        self.run_with(&mut rng)
    }

    /// Execute against a caller-supplied RNG, so repeated shots of a
    /// stochastic circuit draw from one stream rather than restarting it.
    fn run_with(&self, rng: &mut Rng) -> State {
        let mut state = State::new(self.n_qubits);
        // The circuit's one classical register: bit `i` is the last outcome
        // recorded for classical bit `i`, and unmeasured bits read 0.
        let mut clbits: u64 = 0;
        // Trailing measurements are readout: nothing follows them, so
        // collapsing would only discard the prepared state and make the result
        // depend on the seed. `sample` draws the same distribution either way.
        let executable = self.ops.len() - self.trailing_measurements();
        for op in &self.ops[..executable] {
            apply_op(&mut state, op, rng, &mut clbits);
        }
        state
    }
}

impl Op {
    /// The qubit a diagram should mark for this operation — its target, or the
    /// first of a pair. Used only to draw a conditional block in one column.
    fn primary_qubit(&self) -> Option<usize> {
        match self {
            Op::Single { target, .. } => Some(*target),
            Op::Controlled { target, .. } => Some(*target),
            Op::MultiControlled { target, .. } => Some(*target),
            Op::Swap { a, .. } => Some(*a),
            Op::Measure { qubit, .. } | Op::Reset { qubit } | Op::Kraus { qubit, .. } => {
                Some(*qubit)
            }
            Op::Conditional { .. } => None,
        }
    }
}

/// Largest number of distinct classical register values a density-matrix run
/// will carry at once. Each is a full `4^n` matrix, so this bounds the memory
/// a circuit that measures many bits before branching can demand.
const MAX_DENSITY_BRANCHES: usize = 64;

/// Apply one operation to a density matrix. Measurement and conditionals are
/// handled by the caller, which owns the classical mixture.
fn apply_density_op(rho: &mut crate::density::DensityMatrix, op: &Op) {
    match op {
        Op::Single { gate, target, .. } => rho.apply_1q(gate, *target),
        Op::Controlled {
            gate,
            control,
            target,
            ..
        } => rho.apply_controlled_1q(gate, *control, *target),
        Op::Swap { a, b } => rho.swap_qubits(*a, *b),
        Op::MultiControlled {
            gate,
            controls,
            target,
            ..
        } => rho.apply_multi_controlled_1q(gate, controls, *target),
        Op::Reset { qubit } => rho.reset(*qubit),
        Op::Kraus { ops, qubit, .. } => rho.apply_kraus(ops, *qubit),
        // A conditional block holds only gates, so nesting cannot occur; a
        // measurement needs the mixture the caller owns.
        Op::Measure { qubit, .. } => rho.measure_dephase(*qubit),
        Op::Conditional { .. } => unreachable!("conditionals are handled by the caller"),
    }
}

/// Apply one operation, updating the state and the classical register.
fn apply_op(state: &mut State, op: &Op, rng: &mut Rng, clbits: &mut u64) {
    match op {
        Op::Measure { qubit, clbit } => {
            assert!(
                *clbit < state.n_qubits(),
                "classical bit {clbit} out of range"
            );
            let outcome = state.measure_qubit(*qubit, rng);
            let mask = 1u64 << clbit;
            if outcome {
                *clbits |= mask;
            } else {
                *clbits &= !mask;
            }
        }
        // Collapse, then flip the |1> branch back down to |0>.
        Op::Reset { qubit } => {
            if state.measure_qubit(*qubit, rng) {
                state.apply_1q(&gates::x(), *qubit);
            }
        }
        Op::Kraus { ops, qubit, .. } => {
            state.apply_kraus(ops, *qubit, rng);
        }
        Op::Conditional { value, ops } => {
            if *clbits == *value {
                for inner in ops {
                    apply_op(state, inner, rng, clbits);
                }
            }
        }
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

impl Circuit {
    /// Run the circuit as a density matrix, exactly.
    ///
    /// Where [`run`](Circuit::run) samples one trajectory through any noise or
    /// collapse, this carries the whole mixture, so noise channels, measurement
    /// and reset are applied exactly and the result is deterministic — no
    /// shots, no sampling error. An unread measurement is precisely the channel
    /// that erases coherence between its outcomes.
    ///
    /// The cost is `4^n` entries against the state vector's `2^n`, so the
    /// register ceiling is about half; see
    /// [`MAX_DENSITY_QUBITS`](crate::density::MAX_DENSITY_QUBITS).
    ///
    /// Classical feed-forward works too, by carrying one matrix per reachable
    /// classical register value — the quantum state is correlated with the
    /// outcome, so a single matrix cannot express it. That costs memory per
    /// branch, so a circuit measuring many bits before branching on them can
    /// exceed [`DensityError::TooManyBranches`].
    pub fn run_density(&self) -> Result<crate::density::DensityMatrix, DensityError> {
        if self.n_qubits > crate::density::MAX_DENSITY_QUBITS {
            return Err(DensityError::TooManyQubits {
                qubits: self.n_qubits,
                max: crate::density::MAX_DENSITY_QUBITS,
            });
        }
        // Without feed-forward one matrix suffices: an unread measurement is
        // just dephasing. With it, the quantum state is correlated with the
        // classical register, so carry one matrix per reachable register value
        // — unnormalized, its trace being that branch's probability.
        let mut branches: HashMap<u64, crate::density::DensityMatrix> = HashMap::new();
        branches.insert(0, crate::density::DensityMatrix::new(self.n_qubits));
        let tracks_outcomes = self.has_conditional();

        for op in &self.ops {
            match op {
                Op::Conditional { value, ops } => {
                    if let Some(rho) = branches.get_mut(value) {
                        for inner in ops {
                            apply_density_op(rho, inner);
                        }
                    }
                }
                Op::Measure { qubit, clbit } => {
                    if !tracks_outcomes {
                        // Nothing can read the outcome, so the split would only
                        // be summed back together: dephase in place instead.
                        for rho in branches.values_mut() {
                            rho.measure_dephase(*qubit);
                        }
                        continue;
                    }
                    let mut split: HashMap<u64, crate::density::DensityMatrix> = HashMap::new();
                    for (value, rho) in branches.drain() {
                        for outcome in [false, true] {
                            let mut branch = rho.clone();
                            branch.project(*qubit, outcome);
                            // A branch with no weight is not reachable, and
                            // keeping it would double the count for nothing.
                            if branch.trace() < 1e-15 {
                                continue;
                            }
                            let mut key = value & !(1u64 << clbit);
                            if outcome {
                                key |= 1u64 << clbit;
                            }
                            split
                                .entry(key)
                                .and_modify(|existing| existing.add_assign(&branch))
                                .or_insert(branch);
                        }
                    }
                    branches = split;
                    if branches.len() > MAX_DENSITY_BRANCHES {
                        return Err(DensityError::TooManyBranches {
                            max: MAX_DENSITY_BRANCHES,
                        });
                    }
                }
                other => {
                    for rho in branches.values_mut() {
                        apply_density_op(rho, other);
                    }
                }
            }
        }

        // The state of the register, ignoring which classical value came up,
        // is the sum over branches — each already weighted by its probability.
        let mut total = crate::density::DensityMatrix::new(self.n_qubits);
        let mut drained = branches.into_values();
        let mut result = drained.next().expect("at least one branch");
        for rest in drained {
            result.add_assign(&rest);
        }
        std::mem::swap(&mut total, &mut result);
        Ok(total)
    }

    /// Run the circuit `shots` times and return a histogram of measured
    /// basis-state outcomes.
    ///
    /// The circuit is executed once, then each shot measures a fresh clone of
    /// the resulting state so the shots are independent. `seed` makes the
    /// whole sampling run deterministic and reproducible. Keys of the returned
    /// map are little-endian basis-state indices; values are counts.
    pub fn sample(&self, shots: usize, seed: u64) -> HashMap<usize, usize> {
        let mut rng = Rng::new(seed);
        let mut histogram = HashMap::new();

        // A mid-circuit measurement branches the state, so every shot has to
        // re-run the circuit; there is no single final state to draw from.
        // Otherwise — including when the circuit ends with measurements, as
        // most written-out programs do — running once and measuring clones
        // draws from the same distribution, and executes the circuit once
        // instead of `shots` times.
        if self.measures_mid_circuit() {
            for _ in 0..shots {
                let outcome = self.run_with(&mut rng).measure_all(&mut rng);
                *histogram.entry(outcome).or_insert(0) += 1;
            }
        } else {
            let final_state = self.run();
            for _ in 0..shots {
                let outcome = final_state.clone().measure_all(&mut rng);
                *histogram.entry(outcome).or_insert(0) += 1;
            }
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
        if self.has_measurement() || self.has_conditional() {
            out.push_str(&format!("creg c[{}];\n", self.n_qubits));
        }

        for op in &self.ops {
            emit_op(&mut out, self.n_qubits, op)?;
        }
        Ok(out)
    }
}

/// Emit one operation's OpenQASM statements.
fn emit_op(out: &mut String, n_qubits: usize, op: &Op) -> Result<(), ExportError> {
    {
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
                    _ => emit_controlled_u(out, gate, *control, *target),
                }
            }
            Op::Swap { a, b } => {
                out.push_str(&format!("swap q[{a}],q[{b}];\n"));
            }
            // The classical bit is not modelled, so each measurement gets
            // its own bit, named after the qubit it came from.
            Op::Measure { qubit, clbit } => {
                out.push_str(&format!("measure q[{qubit}] -> c[{clbit}];\n"));
            }
            Op::Reset { qubit } => {
                out.push_str(&format!("reset q[{qubit}];\n"));
            }
            // OpenQASM 2 has no syntax for a channel; writing the circuit out
            // without it would silently export a different circuit.
            Op::Kraus { .. } => return Err(ExportError::Noise),
            Op::MultiControlled {
                gate,
                controls,
                target,
                label,
            } => {
                if *label == "X" {
                    emit_mcx(out, n_qubits, controls, *target);
                } else {
                    emit_mcu(out, n_qubits, gate, controls, *target);
                }
            }
            // OpenQASM guards a single statement, so render the block and
            // prefix every line it produced. The condition cannot change
            // inside the block, so guarding each line is the same as
            // guarding the whole.
            Op::Conditional { value, ops } => {
                let mut block = String::new();
                for inner in ops {
                    emit_op(&mut block, n_qubits, inner)?;
                }
                for line in block.lines() {
                    out.push_str(&format!("if(c=={value}) {line}\n"));
                }
            }
        }
    }
    Ok(())
}

impl Circuit {
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
                Op::Measure { qubit, .. } => {
                    col[*qubit] = Some("M");
                }
                // The guarded gates share one column, marked to show they are
                // conditional rather than unconditional.
                Op::Conditional { ops, .. } => {
                    for inner in ops {
                        if let Some(q) = inner.primary_qubit() {
                            col[q] = Some("?");
                        }
                    }
                }
                Op::Reset { qubit } => {
                    col[*qubit] = Some("R");
                }
                Op::Kraus { qubit, label, .. } => {
                    col[*qubit] = Some(label);
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
