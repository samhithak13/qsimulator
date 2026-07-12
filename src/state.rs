//! State-vector representation of an n-qubit register.

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
    pub fn apply_1q(&mut self, gate: &[[Complex64; 2]; 2], target: usize) {
        assert!(target < self.n_qubits, "target qubit out of range");
        let step = 1usize << target;
        let mut i = 0;
        while i < self.amps.len() {
            for j in i..i + step {
                let a = self.amps[j];
                let b = self.amps[j + step];
                self.amps[j] = gate[0][0] * a + gate[0][1] * b;
                self.amps[j + step] = gate[1][0] * a + gate[1][1] * b;
            }
            i += step << 1;
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
        let cmask = 1usize << control;
        let tstep = 1usize << target;
        for j in 0..self.amps.len() {
            // Act on the pair (j, j+tstep) where the target bit of j is 0
            // and the control bit is 1.
            if (j & tstep) == 0 && (j & cmask) != 0 {
                let a = self.amps[j];
                let b = self.amps[j + tstep];
                self.amps[j] = gate[0][0] * a + gate[0][1] * b;
                self.amps[j + tstep] = gate[1][0] * a + gate[1][1] * b;
            }
        }
    }

    /// Total probability (should stay ~1.0). Useful for sanity checks/tests.
    pub fn norm(&self) -> f64 {
        self.amps.iter().map(|c| c.norm_sqr()).sum()
    }
}
