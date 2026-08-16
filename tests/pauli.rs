//! Integration tests for Pauli expectation values.
//!
//! The oracles are analytic: `<Z>` on a known state, `<ZZ>` on a Bell pair, and
//! the decay of `<Z>` under a channel whose rate is known. Where a value is
//! also reachable by sampling, both are checked, since they are independent
//! routes to the same number.

use approx::assert_relative_eq;
use qsimulator::{Circuit, Pauli, PauliString};

fn z(qubit: usize) -> PauliString {
    PauliString::new(&[(qubit, Pauli::Z)])
}

/// `<Z>` is `p(0) - p(1)`, which is the definition and also a second way to
/// compute it from the probabilities.
#[test]
fn single_z_matches_the_probabilities() {
    for theta in [0.0, 0.3, 1.1, std::f64::consts::PI] {
        let mut c = Circuit::new(1);
        c.ry(theta, 0);
        let state = c.run();
        let expected = state.probability(0) - state.probability(1);
        assert_relative_eq!(state.expectation(&z(0)), expected, epsilon = 1e-12);
        // Ry(theta)|0> has <Z> = cos(theta).
        assert_relative_eq!(state.expectation(&z(0)), theta.cos(), epsilon = 1e-12);
    }
}

/// X and Y read the axes a computational-basis histogram cannot see: |+> has
/// `<X> = 1` but `<Z> = 0`, and both look identical when sampled.
#[test]
fn x_and_y_see_what_sampling_cannot() {
    let mut plus = Circuit::new(1);
    plus.h(0);
    let state = plus.run();
    assert_relative_eq!(
        state.expectation(&PauliString::new(&[(0, Pauli::X)])),
        1.0,
        epsilon = 1e-12
    );
    assert_relative_eq!(state.expectation(&z(0)), 0.0, epsilon = 1e-12);

    // S|+> = |i>, which points along +Y.
    let mut plus_i = Circuit::new(1);
    plus_i.h(0).s(0);
    let state = plus_i.run();
    assert_relative_eq!(
        state.expectation(&PauliString::new(&[(0, Pauli::Y)])),
        1.0,
        epsilon = 1e-12
    );
    assert_relative_eq!(
        state.expectation(&PauliString::new(&[(0, Pauli::X)])),
        0.0,
        epsilon = 1e-12
    );
}

/// A Bell pair is perfectly correlated: `<ZZ> = 1` and `<XX> = 1`, while each
/// qubit alone is unbiased.
#[test]
fn bell_pair_correlations() {
    let mut c = Circuit::new(2);
    c.h(0).cnot(0, 1);
    let state = c.run();

    let zz = PauliString::new(&[(0, Pauli::Z), (1, Pauli::Z)]);
    let xx = PauliString::new(&[(0, Pauli::X), (1, Pauli::X)]);
    let yy = PauliString::new(&[(0, Pauli::Y), (1, Pauli::Y)]);
    assert_relative_eq!(state.expectation(&zz), 1.0, epsilon = 1e-12);
    assert_relative_eq!(state.expectation(&xx), 1.0, epsilon = 1e-12);
    // <YY> = -1 for this Bell state, the sign that distinguishes it.
    assert_relative_eq!(state.expectation(&yy), -1.0, epsilon = 1e-12);
    assert_relative_eq!(state.expectation(&z(0)), 0.0, epsilon = 1e-12);
    assert_relative_eq!(state.expectation(&z(1)), 0.0, epsilon = 1e-12);
}

/// The identity is 1 on any state, normalized or not.
#[test]
fn identity_is_one() {
    let mut c = Circuit::new(3);
    c.h(0).cnot(0, 1).u3(0.4, 0.2, -0.9, 2);
    assert_relative_eq!(
        c.run().expectation(&PauliString::identity()),
        1.0,
        epsilon = 1e-12
    );
}

/// Sampling estimates the same number, which is the point of having the exact
/// one: agreement here means the exact route is not measuring something else.
#[test]
fn sampling_estimates_the_same_value() {
    let mut c = Circuit::new(2);
    c.ry(0.9, 0).cnot(0, 1);
    let exact = c
        .run()
        .expectation(&PauliString::new(&[(0, Pauli::Z), (1, Pauli::Z)]));

    let shots = 40_000;
    let hist = c.sample(shots, 13);
    // <ZZ> is (+1) for agreeing outcomes and (-1) for disagreeing ones.
    let estimate: f64 = hist
        .iter()
        .map(|(state, count)| {
            let parity = (state & 1) ^ (state >> 1 & 1);
            let sign = if parity == 0 { 1.0 } else { -1.0 };
            sign * *count as f64
        })
        .sum::<f64>()
        / shots as f64;
    assert!(
        (exact - estimate).abs() < 0.02,
        "exact {exact:.4} vs sampled {estimate:.4}"
    );
}

/// The density backend agrees with the state vector on a pure state, and keeps
/// working on a mixed one — where `<Z>` decaying under damping is the whole
/// reason to want an expectation value under noise.
#[test]
fn density_agrees_and_handles_mixed_states() {
    let mut c = Circuit::new(2);
    c.h(0).cnot(0, 1);
    let zz = PauliString::new(&[(0, Pauli::Z), (1, Pauli::Z)]);
    assert_relative_eq!(
        c.run_density().unwrap().expectation(&zz),
        c.run().expectation(&zz),
        epsilon = 1e-12
    );

    // Amplitude damping on |1> decays <Z> from -1 towards +1 as 2*gamma - 1.
    for gamma in [0.0, 0.25, 0.5, 1.0] {
        let mut damped = Circuit::new(1);
        damped.x(0).amplitude_damping(gamma, 0);
        let rho = damped.run_density().unwrap();
        assert_relative_eq!(rho.expectation(&z(0)), 2.0 * gamma - 1.0, epsilon = 1e-12);
    }

    // Full dephasing leaves <Z> untouched but destroys <X>.
    let mut dephased = Circuit::new(1);
    dephased.h(0).phase_damping(1.0, 0);
    let rho = dephased.run_density().unwrap();
    assert_relative_eq!(rho.expectation(&z(0)), 0.0, epsilon = 1e-12);
    assert_relative_eq!(
        rho.expectation(&PauliString::new(&[(0, Pauli::X)])),
        0.0,
        epsilon = 1e-12
    );
}

/// A string naming a qubit the register does not have is rejected, rather than
/// reading a neighbouring bit.
#[test]
#[should_panic(expected = "outside a 2-qubit register")]
fn a_qubit_outside_the_register_is_rejected() {
    Circuit::new(2).run().expectation(&z(5));
}

/// The one-pass expectation is checked against an explicit matrix product:
/// build the full `2^n x 2^n` Pauli operator by Kronecker product and compute
/// `<psi|P|psi>` directly. That is a different algorithm, so it catches the
/// bugs the fast path could hide — a wrong flip mask, or a Y phase with the
/// sign or the qubit-parity backwards — which analytic spot checks on a few
/// states might not.
#[test]
fn one_pass_expectation_matches_an_explicit_matrix() {
    use num_complex::Complex64;

    fn matrix_of(p: Pauli) -> [[Complex64; 2]; 2] {
        match p {
            Pauli::I => qsimulator::gates::id(),
            Pauli::X => qsimulator::gates::x(),
            Pauli::Y => qsimulator::gates::y(),
            Pauli::Z => qsimulator::gates::z(),
        }
    }

    /// The full operator for `paulis[q]` acting on qubit `q`, little-endian:
    /// entry (r, c) is the product of each qubit's 2x2 entry at that qubit's
    /// bits of r and c.
    fn full_operator(paulis: &[Pauli]) -> Vec<Vec<Complex64>> {
        let dim = 1usize << paulis.len();
        let mut m = vec![vec![Complex64::new(0.0, 0.0); dim]; dim];
        for (r, row) in m.iter_mut().enumerate() {
            for (c, cell) in row.iter_mut().enumerate() {
                let mut product = Complex64::new(1.0, 0.0);
                for (q, p) in paulis.iter().enumerate() {
                    product *= matrix_of(*p)[r >> q & 1][c >> q & 1];
                }
                *cell = product;
            }
        }
        m
    }

    let all = [Pauli::I, Pauli::X, Pauli::Y, Pauli::Z];
    for n in 1..=3usize {
        // A state with distinct amplitudes and real phases, so no term cancels
        // by accident.
        let mut c = Circuit::new(n);
        for q in 0..n {
            let k = q as f64;
            c.u3(0.6 + 0.4 * k, 0.3 - 0.2 * k, 1.2 + 0.5 * k, q);
        }
        if n >= 2 {
            c.cnot(0, n - 1);
        }
        let state = c.run();
        let amps = state.amplitudes();

        // Every Pauli string over `n` qubits.
        let mut indices = vec![0usize; n];
        loop {
            let paulis: Vec<Pauli> = indices.iter().map(|i| all[*i]).collect();
            let terms: Vec<(usize, Pauli)> = paulis.iter().copied().enumerate().collect();
            let fast = state.expectation(&PauliString::new(&terms));

            // <psi|P|psi>, summed straight from the matrix.
            let m = full_operator(&paulis);
            let mut slow = Complex64::new(0.0, 0.0);
            for (r, row) in m.iter().enumerate() {
                for (col, cell) in row.iter().enumerate() {
                    slow += amps[r].conj() * cell * amps[col];
                }
            }
            assert_relative_eq!(fast, slow.re, epsilon = 1e-12);
            assert!(slow.im.abs() < 1e-12, "a Pauli string must be Hermitian");

            // Odometer over all 4^n strings.
            let mut carry = 0;
            while carry < n {
                indices[carry] += 1;
                if indices[carry] < 4 {
                    break;
                }
                indices[carry] = 0;
                carry += 1;
            }
            if carry == n {
                break;
            }
        }
    }
}
