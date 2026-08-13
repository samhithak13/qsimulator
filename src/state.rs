//! State-vector representation of an n-qubit register.

use crate::rng::Rng;
use num_complex::Complex64;
#[cfg(feature = "parallel")]
use rayon::prelude::*;

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
    ///
    /// With the `parallel` feature this runs across threads via rayon,
    /// parallelizing over blocks for a low target qubit and within a block for
    /// a high one — whichever axis carries more work.
    pub fn apply_1q(&mut self, gate: &[[Complex64; 2]; 2], target: usize) {
        assert!(target < self.n_qubits, "target qubit out of range");
        // Hoist the matrix entries so the inner loop does not reload them.
        let (g00, g01) = (gate[0][0], gate[0][1]);
        let (g10, g11) = (gate[1][0], gate[1][1]);
        let butterfly = move |a: &mut Complex64, b: &mut Complex64| {
            let (x, y) = (*a, *b);
            *a = g00 * x + g01 * y;
            *b = g10 * x + g11 * y;
        };
        let step = 1usize << target;
        let block = step << 1;

        #[cfg(feature = "parallel")]
        {
            // A coarse tile so parallel tasks are large enough to amortize
            // scheduling overhead (a block can be as small as two amplitudes).
            const TILE: usize = 1 << 14;
            if block <= TILE {
                // Low target: parallelize over tiles, each holding many blocks.
                self.amps.par_chunks_mut(TILE).for_each(|region| {
                    for blk in region.chunks_exact_mut(block) {
                        let (low, high) = blk.split_at_mut(step);
                        low.iter_mut()
                            .zip(high.iter_mut())
                            .for_each(|(a, b)| butterfly(a, b));
                    }
                });
            } else {
                // High target: few large blocks, so parallelize their halves.
                for blk in self.amps.chunks_exact_mut(block) {
                    let (low, high) = blk.split_at_mut(step);
                    low.par_iter_mut()
                        .zip(high.par_iter_mut())
                        .for_each(|(a, b)| butterfly(a, b));
                }
            }
        }
        #[cfg(not(feature = "parallel"))]
        {
            for blk in self.amps.chunks_exact_mut(block) {
                let (low, high) = blk.split_at_mut(step);
                low.iter_mut()
                    .zip(high.iter_mut())
                    .for_each(|(a, b)| butterfly(a, b));
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
        // With a single control bit, `(j & cmask) != 0` and `== cmask` agree.
        self.apply_masked(gate, 1usize << control, target);
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
        let mut cmask = 0usize;
        for &ctrl in controls {
            assert!(ctrl < self.n_qubits, "control qubit out of range");
            assert_ne!(ctrl, target, "control and target must differ");
            cmask |= 1usize << ctrl;
        }
        self.apply_masked(gate, cmask, target);
    }

    /// Apply `gate`'s butterfly to the `target` bit's `|0>`/`|1>` pairs, but
    /// only for basis states where every bit in `cmask` is set. Shared by the
    /// controlled and multi-controlled kernels (a single control is just a
    /// one-bit `cmask`).
    ///
    /// Structurally identical to [`apply_1q`](Self::apply_1q) — a walk over the
    /// target bit's halves — with a per-pair `cmask` test on the absolute
    /// index. With the `parallel` feature it threads the same way.
    fn apply_masked(&mut self, gate: &[[Complex64; 2]; 2], cmask: usize, target: usize) {
        let (g00, g01) = (gate[0][0], gate[0][1]);
        let (g10, g11) = (gate[1][0], gate[1][1]);
        // `base` is the absolute index of `low[0]` in its block; the pair for
        // `low[k]` sits at index `base + k`, whose control bits we test.
        let butterfly = move |base: usize, k: usize, a: &mut Complex64, b: &mut Complex64| {
            if ((base + k) & cmask) == cmask {
                let (x, y) = (*a, *b);
                *a = g00 * x + g01 * y;
                *b = g10 * x + g11 * y;
            }
        };
        let step = 1usize << target;
        let block = step << 1;

        #[cfg(feature = "parallel")]
        {
            const TILE: usize = 1 << 14;
            if block <= TILE {
                self.amps
                    .par_chunks_mut(TILE)
                    .enumerate()
                    .for_each(|(t, region)| {
                        let region_base = t * TILE;
                        for (bi, blk) in region.chunks_exact_mut(block).enumerate() {
                            let base = region_base + bi * block;
                            let (low, high) = blk.split_at_mut(step);
                            for (k, (a, b)) in low.iter_mut().zip(high.iter_mut()).enumerate() {
                                butterfly(base, k, a, b);
                            }
                        }
                    });
            } else {
                for (bi, blk) in self.amps.chunks_exact_mut(block).enumerate() {
                    let base = bi * block;
                    let (low, high) = blk.split_at_mut(step);
                    low.par_iter_mut()
                        .zip(high.par_iter_mut())
                        .enumerate()
                        .for_each(|(k, (a, b))| butterfly(base, k, a, b));
                }
            }
        }
        #[cfg(not(feature = "parallel"))]
        {
            for (bi, blk) in self.amps.chunks_exact_mut(block).enumerate() {
                let base = bi * block;
                let (low, high) = blk.split_at_mut(step);
                for (k, (a, b)) in low.iter_mut().zip(high.iter_mut()).enumerate() {
                    butterfly(base, k, a, b);
                }
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

    /// Expectation `<psi|M|psi>` of a single-qubit Hermitian `m` on qubit `q`.
    ///
    /// One read-only pass over the amplitudes, so a Kraus branch probability
    /// costs no allocation and no clone of the state vector.
    fn expectation_1q(&self, m: &[[Complex64; 2]; 2], q: usize) -> f64 {
        let mask = 1usize << q;
        let mut total = 0.0;
        for (i, a0) in self.amps.iter().enumerate() {
            if i & mask != 0 {
                continue;
            }
            let a1 = self.amps[i | mask];
            let v0 = m[0][0] * a0 + m[0][1] * a1;
            let v1 = m[1][0] * a0 + m[1][1] * a1;
            // `m` is Hermitian, so this is real; the imaginary part is rounding.
            total += (a0.conj() * v0 + a1.conj() * v1).re;
        }
        total
    }

    /// Apply a quantum channel to `q` by sampling one Kraus operator, returning
    /// which one was taken.
    ///
    /// This is the quantum-trajectory (Monte Carlo wave function) method: a
    /// channel that would otherwise need a density matrix is simulated by
    /// choosing operator `K_i` with probability `<psi|K_i^dagger K_i|psi>` and
    /// renormalizing. Averaged over shots the ensemble reproduces the density
    /// matrix, so a state vector suffices — at the cost of needing many shots.
    ///
    /// `ops` should be trace preserving; see
    /// [`noise::is_trace_preserving`](crate::noise::is_trace_preserving). The
    /// probabilities are normalized before sampling regardless, so a set that
    /// is slightly off from rounding still behaves.
    pub fn apply_kraus(&mut self, ops: &[[[Complex64; 2]; 2]], q: usize, rng: &mut Rng) -> usize {
        assert!(q < self.n_qubits, "qubit out of range");
        assert!(
            !ops.is_empty(),
            "a channel needs at least one Kraus operator"
        );

        let probs: Vec<f64> = ops
            .iter()
            .map(|k| self.expectation_1q(&dagger_product(k), q).max(0.0))
            .collect();
        let total: f64 = probs.iter().sum();
        assert!(total > 0.0, "channel has zero probability on this state");

        // Walk the cumulative distribution to pick a branch.
        let mut r = rng.next_f64() * total;
        let mut chosen = probs.len() - 1;
        for (i, p) in probs.iter().enumerate() {
            if r < *p {
                chosen = i;
                break;
            }
            r -= p;
        }

        self.apply_1q(&ops[chosen], q);
        // The chosen branch carried probability `probs[chosen]`; scale back to
        // a unit vector.
        let scale = 1.0 / probs[chosen].sqrt();
        for a in self.amps.iter_mut() {
            *a *= scale;
        }
        chosen
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

/// `K^dagger K` for a 2x2 operator — the positive matrix whose expectation on
/// the state is the probability of taking that Kraus branch.
fn dagger_product(k: &[[Complex64; 2]; 2]) -> [[Complex64; 2]; 2] {
    let entry = |i: usize, j: usize| k[0][i].conj() * k[0][j] + k[1][i].conj() * k[1][j];
    [[entry(0, 0), entry(0, 1)], [entry(1, 0), entry(1, 1)]]
}
