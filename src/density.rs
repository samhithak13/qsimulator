//! Density-matrix representation, for exact noise.
//!
//! A state vector cannot represent a mixed state, so [`State`](crate::State)
//! simulates a noise channel by sampling one Kraus operator per shot and
//! averaging over many runs. A density matrix `rho` represents the mixture
//! directly, so the same channel is applied exactly — `rho` becomes
//! `sum_i K_i rho K_i^dagger` in one step, with no sampling error at all.
//!
//! The cost is memory. An `n`-qubit density matrix is `2^n x 2^n`, so it holds
//! `4^n` amplitudes where the state vector holds `2^n`: 12 qubits here is the
//! same footprint as 24 there. That is the whole trade — exactness and mixed
//! states, against half the reach.
//!
//! Measurement and reset are exact too. An outcome nobody reads is the channel
//! `rho -> P_0 rho P_0 + P_1 rho P_1`, which is deterministic. Classical
//! feed-forward is exact as well, but needs more than one matrix: the quantum
//! state becomes correlated with the classical outcome, so
//! [`Circuit::run_density`](crate::Circuit::run_density) carries one matrix per
//! reachable register value and sums them at the end.

use num_complex::Complex64;

use crate::gates;

/// Largest register a density matrix will allocate. `4^n` entries at 16 bytes
/// each puts 12 qubits at ~268 MB; 13 would be over a gigabyte.
pub const MAX_DENSITY_QUBITS: usize = 12;

/// The density matrix of an `n`-qubit register, row-major over `2^n x 2^n`.
///
/// Rows and columns use the same little-endian basis-state indexing as
/// [`State`](crate::State), so entry `(i, j)` is `<i|rho|j>`.
#[derive(Debug, Clone, PartialEq)]
pub struct DensityMatrix {
    n_qubits: usize,
    dim: usize,
    entries: Vec<Complex64>,
}

impl DensityMatrix {
    /// The pure state |0...0><0...0| for `n_qubits` qubits.
    ///
    /// # Panics
    ///
    /// Above [`MAX_DENSITY_QUBITS`], where the allocation would be measured in
    /// gigabytes.
    pub fn new(n_qubits: usize) -> Self {
        assert!(
            n_qubits <= MAX_DENSITY_QUBITS,
            "a density matrix is 4^n entries; {n_qubits} qubits exceeds the \
             maximum of {MAX_DENSITY_QUBITS}"
        );
        let dim = 1usize << n_qubits;
        let mut entries = vec![Complex64::new(0.0, 0.0); dim * dim];
        entries[0] = Complex64::new(1.0, 0.0);
        DensityMatrix {
            n_qubits,
            dim,
            entries,
        }
    }

    /// A matrix of zeros — not a physical state, but the identity for summing
    /// the branches of a classical mixture.
    pub(crate) fn zeros(n_qubits: usize) -> Self {
        let mut rho = DensityMatrix::new(n_qubits);
        rho.entries[0] = Complex64::new(0.0, 0.0);
        rho
    }

    /// The pure state `|psi><psi|` for a state vector.
    pub fn from_state(state: &crate::State) -> Self {
        let n_qubits = state.n_qubits();
        let mut rho = DensityMatrix::new(n_qubits);
        let amps = state.amplitudes();
        for (i, a) in amps.iter().enumerate() {
            for (j, b) in amps.iter().enumerate() {
                rho.entries[i * rho.dim + j] = a * b.conj();
            }
        }
        rho
    }

    /// Number of qubits in the register.
    pub fn n_qubits(&self) -> usize {
        self.n_qubits
    }

    /// Read-only access to the row-major entries.
    pub fn entries(&self) -> &[Complex64] {
        &self.entries
    }

    /// The entry `<row|rho|col>`.
    pub fn entry(&self, row: usize, col: usize) -> Complex64 {
        self.entries[row * self.dim + col]
    }

    /// Probability of measuring basis state `index`, the diagonal entry.
    pub fn probability(&self, index: usize) -> f64 {
        self.entry(index, index).re
    }

    /// Trace of the matrix, which a physical state keeps at 1.
    pub fn trace(&self) -> f64 {
        (0..self.dim).map(|i| self.probability(i)).sum()
    }

    /// Purity `Tr(rho^2)`: 1 for a pure state, down to `1/2^n` for the
    /// maximally mixed one. The direct measure of how much noise has cost.
    pub fn purity(&self) -> f64 {
        // Tr(rho^2) = sum_{i,j} rho[i][j] * rho[j][i], and rho is Hermitian, so
        // that is the sum of |rho[i][j]|^2.
        self.entries.iter().map(|e| e.norm_sqr()).sum()
    }

    /// Probability that qubit `q` reads |1>.
    pub fn prob_qubit_one(&self, q: usize) -> f64 {
        assert!(q < self.n_qubits, "qubit out of range");
        let mask = 1usize << q;
        (0..self.dim)
            .filter(|i| i & mask != 0)
            .map(|i| self.probability(i))
            .sum()
    }

    /// Apply a single-qubit unitary: `rho -> U rho U^dagger`.
    pub fn apply_1q(&mut self, gate: &[[Complex64; 2]; 2], target: usize) {
        assert!(target < self.n_qubits, "target qubit out of range");
        self.apply_masked(gate, 0, target);
    }

    /// Apply a controlled single-qubit unitary.
    pub fn apply_controlled_1q(
        &mut self,
        gate: &[[Complex64; 2]; 2],
        control: usize,
        target: usize,
    ) {
        assert_ne!(control, target, "control and target must differ");
        self.apply_multi_controlled_1q(gate, &[control], target);
    }

    /// Apply a single-qubit unitary where every qubit in `controls` is |1>.
    pub fn apply_multi_controlled_1q(
        &mut self,
        gate: &[[Complex64; 2]; 2],
        controls: &[usize],
        target: usize,
    ) {
        assert!(target < self.n_qubits, "target qubit out of range");
        let mut cmask = 0usize;
        for &c in controls {
            assert!(c < self.n_qubits, "control qubit out of range");
            assert_ne!(c, target, "control and target must differ");
            cmask |= 1usize << c;
        }
        self.apply_masked(gate, cmask, target);
    }

    /// `rho -> U rho U^dagger` restricted to basis states matching `cmask`.
    ///
    /// Conjugation is two sweeps: `U` acting on the row index, then `U`
    /// conjugated acting on the column index. Writing it that way means the
    /// same butterfly serves both, rather than materializing `U rho` and
    /// multiplying again.
    fn apply_masked(&mut self, gate: &[[Complex64; 2]; 2], cmask: usize, target: usize) {
        let mask = 1usize << target;
        let dim = self.dim;

        // Left: rows carry U.
        for row in 0..dim {
            if row & mask != 0 || (row & cmask) != cmask {
                continue;
            }
            for col in 0..dim {
                let (a, b) = (
                    self.entries[row * dim + col],
                    self.entries[(row | mask) * dim + col],
                );
                self.entries[row * dim + col] = gate[0][0] * a + gate[0][1] * b;
                self.entries[(row | mask) * dim + col] = gate[1][0] * a + gate[1][1] * b;
            }
        }

        // Right: columns carry U^dagger, which is `conj(U)` applied along the
        // column index.
        for row in 0..dim {
            for col in 0..dim {
                if col & mask != 0 || (col & cmask) != cmask {
                    continue;
                }
                let (a, b) = (
                    self.entries[row * dim + col],
                    self.entries[row * dim + (col | mask)],
                );
                self.entries[row * dim + col] = gate[0][0].conj() * a + gate[0][1].conj() * b;
                self.entries[row * dim + (col | mask)] =
                    gate[1][0].conj() * a + gate[1][1].conj() * b;
            }
        }
    }

    /// Exchange qubits `a` and `b`.
    pub fn swap_qubits(&mut self, a: usize, b: usize) {
        if a == b {
            return;
        }
        // SWAP is three CNOTs, which keeps this to the conjugation already
        // implemented rather than a second index-permutation routine.
        self.apply_controlled_1q(&gates::x(), a, b);
        self.apply_controlled_1q(&gates::x(), b, a);
        self.apply_controlled_1q(&gates::x(), a, b);
    }

    /// Apply a quantum channel exactly: `rho -> sum_i K_i rho K_i^dagger`.
    ///
    /// This is where a density matrix earns its memory. The state vector has to
    /// sample one operator per shot and average; here the whole mixture is
    /// carried, so one call gives the exact result.
    pub fn apply_kraus(&mut self, ops: &[[[Complex64; 2]; 2]], target: usize) {
        assert!(target < self.n_qubits, "target qubit out of range");
        assert!(
            !ops.is_empty(),
            "a channel needs at least one Kraus operator"
        );

        // One scratch matrix, refilled from `self` each time, rather than a
        // fresh clone per operator: depolarizing has four, and at the register
        // ceiling each clone is hundreds of megabytes.
        let mut total = vec![Complex64::new(0.0, 0.0); self.entries.len()];
        let mut scratch = self.clone();
        for k in ops {
            scratch.entries.copy_from_slice(&self.entries);
            scratch.apply_masked(k, 0, target);
            for (acc, term) in total.iter_mut().zip(scratch.entries.iter()) {
                *acc += term;
            }
        }
        self.entries = total;
    }

    /// Project onto qubit `q` reading `outcome`, without renormalizing:
    /// `rho -> P rho P`.
    ///
    /// The trace of the result is the probability of that outcome, which is
    /// what lets a mixture over classical outcomes carry its own weights.
    pub fn project(&mut self, q: usize, outcome: bool) {
        assert!(q < self.n_qubits, "qubit out of range");
        let mask = 1usize << q;
        let want = if outcome { mask } else { 0 };
        let dim = self.dim;
        for row in 0..dim {
            for col in 0..dim {
                if (row & mask) != want || (col & mask) != want {
                    self.entries[row * dim + col] = Complex64::new(0.0, 0.0);
                }
            }
        }
    }

    /// Add another matrix into this one, entry by entry.
    ///
    /// Used to recombine the branches of a classical mixture; the operands are
    /// unnormalized, and their traces are the branch weights.
    pub fn add_assign(&mut self, other: &DensityMatrix) {
        assert_eq!(self.n_qubits, other.n_qubits, "registers must match");
        for (a, b) in self.entries.iter_mut().zip(other.entries.iter()) {
            *a += b;
        }
    }

    /// Measure `q` without reading the outcome: `rho -> P_0 rho P_0 + P_1 rho
    /// P_1`.
    ///
    /// Exact and deterministic — an unread measurement is just the channel that
    /// erases coherence between the two outcomes. Reading the outcome, and
    /// acting on it, is what a density matrix alone cannot express.
    pub fn measure_dephase(&mut self, q: usize) {
        assert!(q < self.n_qubits, "qubit out of range");
        let mask = 1usize << q;
        let dim = self.dim;
        for row in 0..dim {
            for col in 0..dim {
                // The projectors keep only entries whose row and column agree
                // on qubit `q`; the rest are coherences between outcomes.
                if (row & mask) != (col & mask) {
                    self.entries[row * dim + col] = Complex64::new(0.0, 0.0);
                }
            }
        }
    }

    /// Reset `q` to |0>, exactly: dephase, then move the |1> population down.
    pub fn reset(&mut self, q: usize) {
        assert!(q < self.n_qubits, "qubit out of range");
        let c = |re: f64| Complex64::new(re, 0.0);
        // K0 = |0><0|, K1 = |0><1| is the reset channel.
        let k0 = [[c(1.0), c(0.0)], [c(0.0), c(0.0)]];
        let k1 = [[c(0.0), c(1.0)], [c(0.0), c(0.0)]];
        self.apply_kraus(&[k0, k1], q);
    }
}
