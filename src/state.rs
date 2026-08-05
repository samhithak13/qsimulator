//! State-vector representation of an n-qubit register.

use crate::rng::Rng;
use num_complex::Complex64;

/// The full state vector of an `n`-qubit quantum register.
///
/// Amplitudes are stored in little-endian order: index `i` corresponds to
/// the computational basis state whose bit `q` equals `(i >> q) & 1`.
#[derive(Debug, Clone, PartialEq)]
pub struct State {
    n_qubits: usize,
    amps: Vec<Complex64>,
}

impl State {
    /// Create the |0...0> state for `n_qubits` qubits.
    pub fn new(n_qubits: usize) -> Self {
        let mut amps = vec![Complex64::new(0.0, 0.0); 1usize << n_qubits];
        amps[0] = Complex64::new(1.0, 0.0);
        State { n_qubits, amps }
    }

    /// Number of qubits in the register.
    pub fn n_qubits(&self) -> usize {
        self.n_qubits
    }

    /// Read-only access to the amplitude vector.
    pub fn amplitudes(&self) -> &[Complex64] {
        &self.amps
    }

    /// Probability of measuring basis state `index`.
    pub fn probability(&self, index: usize) -> f64 {
        self.amps[index].norm_sqr()
    }

    /// Apply a single-qubit gate (2x2 unitary) to qubit `target`.
    ///
    /// The amplitude vector splits into contiguous blocks of `2·2^target`,
    /// and within each block the low and high halves are the `|0>`/`|1>`
    /// partners of the target bit. Iterating over those halves with slice
    /// iterators keeps the inner loop bounds-check-free.
    pub fn apply_1q(&mut self, gate: &[[Complex64; 2]; 2], target: usize) {
        assert!(target < self.n_qubits, "target qubit out of range");
        // Hoist the matrix entries so the inner loop does not reload them.
        let (g00, g01) = (gate[0][0], gate[0][1]);
        let (g10, g11) = (gate[1][0], gate[1][1]);
        let step = 1usize << target;
        for block in self.amps.chunks_exact_mut(step << 1) {
            let (low, high) = block.split_at_mut(step);
            for (a, b) in low.iter_mut().zip(high.iter_mut()) {
                let (x, y) = (*a, *b);
                *a = g00 * x + g01 * y;
                *b = g10 * x + g11 * y;
            }
        }
    }

    /// Apply a controlled single-qubit gate: `gate` acts on `target`
    /// only when `control` is |1>.
    pub fn apply_controlled_1q(
        &mut self,
        gate: &[[Complex64; 2]; 2],
        control: usize,
        target: usize,
    ) {
        assert!(control < self.n_qubits, "control qubit out of range");
        assert!(target < self.n_qubits, "target qubit out of range");
        assert_ne!(control, target, "control and target must differ");
        let (g00, g01) = (gate[0][0], gate[0][1]);
        let (g10, g11) = (gate[1][0], gate[1][1]);
        let cmask = 1usize << control;
        let tstep = 1usize << target;
        for j in 0..self.amps.len() {
            // Act on the pair (j, j+tstep) where the target bit of j is 0
            // and the control bit is 1.
            if (j & tstep) == 0 && (j & cmask) != 0 {
                let a = self.amps[j];
                let b = self.amps[j + tstep];
                self.amps[j] = g00 * a + g01 * b;
                self.amps[j + tstep] = g10 * a + g11 * b;
            }
        }
    }

    /// Apply a single-qubit `gate` to `target`, but only on basis states
    /// where **every** qubit in `controls` is |1>.
    ///
    /// This generalizes [`apply_controlled_1q`](Self::apply_controlled_1q) to
    /// any number of controls: zero controls degenerate to an unconditional
    /// gate, one control reproduces `apply_controlled_1q`, and two controls
    /// with an X gate give a Toffoli. Controls must be distinct from `target`.
    pub fn apply_multi_controlled_1q(
        &mut self,
        gate: &[[Complex64; 2]; 2],
        controls: &[usize],
        target: usize,
    ) {
        assert!(target < self.n_qubits, "target qubit out of range");
        let (g00, g01) = (gate[0][0], gate[0][1]);
        let (g10, g11) = (gate[1][0], gate[1][1]);
        let mut cmask = 0usize;
        for &ctrl in controls {
            assert!(ctrl < self.n_qubits, "control qubit out of range");
            assert_ne!(ctrl, target, "control and target must differ");
            cmask |= 1usize << ctrl;
        }
        let tstep = 1usize << target;
        for j in 0..self.amps.len() {
            // Act once per pair (target bit 0) and only when all controls set.
            if (j & tstep) == 0 && (j & cmask) == cmask {
                let a = self.amps[j];
                let b = self.amps[j + tstep];
                self.amps[j] = g00 * a + g01 * b;
                self.amps[j + tstep] = g10 * a + g11 * b;
            }
        }
    }

    /// Exchange the states of qubits `a` and `b`.
    ///
    /// Implemented by swapping the amplitudes of every basis-state pair that
    /// differs only in these two bits (i.e. |..1..0..> <-> |..0..1..>). Each
    /// pair is touched exactly once.
    pub fn swap_qubits(&mut self, a: usize, b: usize) {
        assert!(a < self.n_qubits, "qubit a out of range");
        assert!(b < self.n_qubits, "qubit b out of range");
        if a == b {
            return;
        }
        let amask = 1usize << a;
        let bmask = 1usize << b;
        for i in 0..self.amps.len() {
            // Act once per pair: pick the index with bit a set and bit b clear,
            // then swap it with its partner (bit a clear, bit b set).
            if (i & amask) != 0 && (i & bmask) == 0 {
                let j = (i & !amask) | bmask;
                self.amps.swap(i, j);
            }
        }
    }

    /// Total probability (should stay ~1.0). Useful for sanity checks/tests.
    pub fn norm(&self) -> f64 {
        self.amps.iter().map(|c| c.norm_sqr()).sum()
    }

    /// Probability that qubit `q` is measured in the |1> state.
    ///
    /// This is the sum of `|amplitude|^2` over every basis state whose
    /// bit `q` is set. Does not modify the state.
    pub fn prob_qubit_one(&self, q: usize) -> f64 {
        assert!(q < self.n_qubits, "qubit out of range");
        let mask = 1usize << q;
        self.amps
            .iter()
            .enumerate()
            .filter(|(i, _)| (i & mask) != 0)
            .map(|(_, a)| a.norm_sqr())
            .sum()
    }

    /// Measure qubit `q` in the computational basis.
    ///
    /// Samples an outcome (`true` = |1>, `false` = |0>) via the Born rule,
    /// then collapses the state onto the measured subspace and renormalizes
    /// so the surviving amplitudes again sum to probability 1.
    pub fn measure_qubit(&mut self, q: usize, rng: &mut Rng) -> bool {
        assert!(q < self.n_qubits, "qubit out of range");
        let mask = 1usize << q;
        let p1 = self.prob_qubit_one(q);
        let outcome = rng.next_f64() < p1;

        // Probability mass of the branch we are keeping.
        let branch_prob = if outcome { p1 } else { 1.0 - p1 };
        // Guard against dividing by zero if the branch has (numerically) no
        // amplitude; in that case the surviving amplitudes are all zero anyway.
        let inv_norm = if branch_prob > 0.0 {
            1.0 / branch_prob.sqrt()
        } else {
            0.0
        };

        for (i, a) in self.amps.iter_mut().enumerate() {
            let bit_set = (i & mask) != 0;
            if bit_set == outcome {
                *a *= inv_norm;
            } else {
                *a = Complex64::new(0.0, 0.0);
            }
        }
        outcome
    }

    /// Measure every qubit at once, returning the sampled basis-state index.
    ///
    /// The outcome is drawn from the full Born-rule distribution
    /// `p(i) = |amplitude(i)|^2` and the state collapses onto that single
    /// basis state. The returned `usize` is little-endian: bit `q` is the
    /// measured value of qubit `q`.
    pub fn measure_all(&mut self, rng: &mut Rng) -> usize {
        let r = rng.next_f64();
        let mut cumulative = 0.0;
        // Default to the last basis state carrying nonzero amplitude. This
        // guards the edge case where floating-point roundoff leaves the
        // cumulative sum just shy of `r`.
        let mut outcome = 0;
        for (i, a) in self.amps.iter().enumerate() {
            let p = a.norm_sqr();
            if p > 0.0 {
                outcome = i;
            }
            cumulative += p;
            if r < cumulative {
                outcome = i;
                break;
            }
        }

        for (i, a) in self.amps.iter_mut().enumerate() {
            *a = if i == outcome {
                Complex64::new(1.0, 0.0)
            } else {
                Complex64::new(0.0, 0.0)
            };
        }
        outcome
    }
}
