//! Property-based tests of invariants that must hold for every circuit,
//! independent of the specific gates: unitarity (the norm is preserved) and
//! reversibility (a gate followed by its inverse is the identity).

use qsimulator::{Circuit, Rng};
use std::f64::consts::PI;

/// Append a random gate drawn from the full gate set to `c`.
fn random_gate(c: &mut Circuit, n: usize, rng: &mut Rng) {
    let pick = (rng.next_f64() * 16.0) as u32;
    let q = |rng: &mut Rng| (rng.next_f64() * n as f64) as usize % n;
    let angle = |rng: &mut Rng| (rng.next_f64() - 0.5) * 4.0 * PI;
    // Two distinct qubits (n >= 2 guaranteed by the caller for these arms).
    let two = |rng: &mut Rng| {
        let a = q(rng);
        let mut b = q(rng);
        while b == a {
            b = q(rng);
        }
        (a, b)
    };
    match pick {
        0 => c.h(q(rng)),
        1 => c.x(q(rng)),
        2 => c.y(q(rng)),
        3 => c.z(q(rng)),
        4 => c.s(q(rng)),
        5 => c.t(q(rng)),
        6 => c.rx(angle(rng), q(rng)),
        7 => c.ry(angle(rng), q(rng)),
        8 => c.rz(angle(rng), q(rng)),
        9 => c.p(angle(rng), q(rng)),
        10 => c.u3(angle(rng), angle(rng), angle(rng), q(rng)),
        11 => {
            let (a, b) = two(rng);
            c.cnot(a, b)
        }
        12 => {
            let (a, b) = two(rng);
            c.cz(a, b)
        }
        13 => {
            let (a, b) = two(rng);
            c.swap(a, b)
        }
        14 => {
            let (a, b) = two(rng);
            c.cp(angle(rng), a, b)
        }
        _ => {
            let (a, b) = two(rng);
            c.cu3(angle(rng), angle(rng), angle(rng), a, b)
        }
    };
}

/// A random circuit of unitary gates keeps the state normalized.
#[test]
fn random_circuits_preserve_norm() {
    let mut rng = Rng::new(0xB0BA_CAFE);
    for _ in 0..2000 {
        let n = 2 + (rng.next_f64() * 4.0) as usize; // 2..=5 qubits
        let depth = 1 + (rng.next_f64() * 30.0) as usize;
        let mut c = Circuit::new(n);
        for _ in 0..depth {
            random_gate(&mut c, n, &mut rng);
        }
        let norm = c.run().norm();
        assert!((norm - 1.0).abs() < 1e-9, "norm drifted to {norm}");
    }
}

/// Self-inverse gates applied twice return to the start.
#[test]
fn self_inverse_gates_cancel() {
    let mut rng = Rng::new(7);
    for _ in 0..500 {
        let mut c = Circuit::new(3);
        // Scramble first so the check is not just on |000>.
        for _ in 0..6 {
            random_gate(&mut c, 3, &mut rng);
        }
        let before = c.clone().run();

        // Each of these is its own inverse; applying twice is the identity.
        c.h(0).h(0);
        c.x(1).x(1);
        c.cnot(0, 2).cnot(0, 2);
        c.swap(1, 2).swap(1, 2);

        let after = c.run();
        for i in 0..(1usize << 3) {
            let d = (before.amplitudes()[i] - after.amplitudes()[i]).norm();
            assert!(d < 1e-12, "state changed at index {i} by {d}");
        }
    }
}

/// A rotation followed by its negation is the identity.
#[test]
fn inverse_rotations_cancel() {
    let mut c = Circuit::new(1);
    c.h(0).rz(0.7, 0).rz(-0.7, 0).rx(1.3, 0).rx(-1.3, 0);
    let state = c.run();
    // Back to H|0> = |+>: equal, real, positive amplitudes.
    let inv_sqrt2 = std::f64::consts::FRAC_1_SQRT_2;
    assert!((state.amplitudes()[0].re - inv_sqrt2).abs() < 1e-12);
    assert!((state.amplitudes()[1].re - inv_sqrt2).abs() < 1e-12);
    assert!(state.amplitudes()[0].im.abs() < 1e-12);
    assert!(state.amplitudes()[1].im.abs() < 1e-12);
}
