//! Pauli strings and their expectation values.
//!
//! Sampling answers "what outcomes come up"; an expectation value answers
//! "what is `<P>`", which is what a variational algorithm actually optimizes.
//! Getting `<Z_0 Z_1>` by sampling costs shots and converges as `1/sqrt(shots)`;
//! reading it off the state costs one pass and is exact.
//!
//! A [`PauliString`] is a product of single-qubit Paulis over named qubits,
//! with the identity everywhere else — so `Z_0 Z_1` on a five-qubit register is
//! two terms, not five. Both backends evaluate it:
//! [`State::expectation`](crate::State::expectation) and
//! [`DensityMatrix::expectation`](crate::density::DensityMatrix::expectation).

/// A single-qubit Pauli operator.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Pauli {
    /// The identity, which contributes nothing.
    I,
    /// Pauli-X.
    X,
    /// Pauli-Y.
    Y,
    /// Pauli-Z.
    Z,
}

/// A product of single-qubit Paulis, one per named qubit.
///
/// Qubits not named carry the identity. Building one validates that no qubit is
/// named twice, since `Z_0 Z_0` is the identity rather than a Pauli string and
/// almost always means a mistake at the call site.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct PauliString {
    terms: Vec<(usize, Pauli)>,
}

impl PauliString {
    /// The identity: `<I> == 1` for any state.
    pub fn identity() -> Self {
        PauliString { terms: Vec::new() }
    }

    /// Build from `(qubit, pauli)` pairs, dropping identities.
    ///
    /// # Panics
    ///
    /// If a qubit is named more than once.
    pub fn new(terms: &[(usize, Pauli)]) -> Self {
        let mut kept: Vec<(usize, Pauli)> = Vec::with_capacity(terms.len());
        for &(qubit, pauli) in terms {
            assert!(
                !kept.iter().any(|(q, _)| *q == qubit),
                "qubit {qubit} appears twice in a Pauli string"
            );
            // The identity multiplies to nothing, so carrying it would only
            // cost work per amplitude.
            if pauli != Pauli::I {
                kept.push((qubit, pauli));
            }
        }
        PauliString { terms: kept }
    }

    /// The non-identity terms, as `(qubit, pauli)` pairs.
    pub fn terms(&self) -> &[(usize, Pauli)] {
        &self.terms
    }

    /// Highest qubit index named, or `None` for the identity. A caller can use
    /// this to check the string fits the register it is about to be applied to.
    pub fn max_qubit(&self) -> Option<usize> {
        self.terms.iter().map(|(q, _)| *q).max()
    }

    /// The bits this string flips: the qubits carrying X or Y.
    ///
    /// `P|j>` is always a single basis state `c|j ^ flip_mask>`, which is what
    /// makes an expectation value one pass rather than a matrix product.
    pub(crate) fn flip_mask(&self) -> usize {
        self.terms
            .iter()
            .filter(|(_, p)| matches!(p, Pauli::X | Pauli::Y))
            .fold(0, |mask, (q, _)| mask | 1usize << q)
    }

    /// The coefficient `c` in `P|basis> = c |basis ^ flip_mask>`.
    ///
    /// Z contributes `-1` when its qubit is set; Y contributes `i` when its
    /// qubit is clear and `-i` when set (`Y|0> = i|1>`, `Y|1> = -i|0>`); X
    /// contributes nothing beyond the flip.
    pub(crate) fn coefficient(&self, basis: usize) -> num_complex::Complex64 {
        let mut real = 1.0f64;
        // Track the power of `i` separately, so the result stays exact rather
        // than accumulating rounding through repeated complex multiplies.
        let mut i_power = 0u32;
        for &(qubit, pauli) in &self.terms {
            let set = basis >> qubit & 1 == 1;
            match pauli {
                Pauli::I | Pauli::X => {}
                Pauli::Z => {
                    if set {
                        real = -real;
                    }
                }
                Pauli::Y => {
                    i_power += 1;
                    if set {
                        real = -real;
                    }
                }
            }
        }
        let unit = match i_power % 4 {
            0 => num_complex::Complex64::new(1.0, 0.0),
            1 => num_complex::Complex64::new(0.0, 1.0),
            2 => num_complex::Complex64::new(-1.0, 0.0),
            _ => num_complex::Complex64::new(0.0, -1.0),
        };
        unit * real
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identities_are_dropped_but_the_string_still_works() {
        let p = PauliString::new(&[(0, Pauli::I), (1, Pauli::Z), (2, Pauli::I)]);
        assert_eq!(p.terms(), &[(1, Pauli::Z)]);
        assert_eq!(p.max_qubit(), Some(1));
        assert_eq!(PauliString::identity().max_qubit(), None);
    }

    #[test]
    #[should_panic(expected = "appears twice")]
    fn a_repeated_qubit_is_rejected() {
        PauliString::new(&[(0, Pauli::X), (0, Pauli::Z)]);
    }

    #[test]
    fn flip_mask_covers_x_and_y_only() {
        let p = PauliString::new(&[(0, Pauli::X), (1, Pauli::Y), (2, Pauli::Z)]);
        assert_eq!(p.flip_mask(), 0b011);
    }

    /// The coefficients are the Pauli matrices' own entries, so check them
    /// against the definitions rather than against a reimplementation.
    #[test]
    fn coefficients_match_the_matrices() {
        let c = |re: f64, im: f64| num_complex::Complex64::new(re, im);

        // Z|0> = |0>, Z|1> = -|1>.
        let z = PauliString::new(&[(0, Pauli::Z)]);
        assert_eq!(z.coefficient(0b0), c(1.0, 0.0));
        assert_eq!(z.coefficient(0b1), c(-1.0, 0.0));

        // Y|0> = i|1>, Y|1> = -i|0>.
        let y = PauliString::new(&[(0, Pauli::Y)]);
        assert_eq!(y.coefficient(0b0), c(0.0, 1.0));
        assert_eq!(y.coefficient(0b1), c(0.0, -1.0));

        // X carries no phase.
        let x = PauliString::new(&[(0, Pauli::X)]);
        assert_eq!(x.coefficient(0b0), c(1.0, 0.0));
        assert_eq!(x.coefficient(0b1), c(1.0, 0.0));

        // Two Ys give i^2 = -1 when both qubits are clear.
        let yy = PauliString::new(&[(0, Pauli::Y), (1, Pauli::Y)]);
        assert_eq!(yy.coefficient(0b00), c(-1.0, 0.0));
        // And both set flips the sign twice more, back to -1.
        assert_eq!(yy.coefficient(0b11), c(-1.0, 0.0));

        assert_eq!(PauliString::identity().coefficient(0b101), c(1.0, 0.0));
    }
}
